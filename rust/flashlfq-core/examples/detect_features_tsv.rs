//! End-to-end untargeted feature-detection runner.
//!
//! Reads a spectra file (mzML or Thermo `.raw`), runs the full untargeted pipeline
//! (index → `detect_features` → `refine_feature` → `resolve_charge_state_consensus`), writes the
//! resolved peptide-level features to a human-readable TSV, and — if given a base-FlashLFQ
//! `AllQuantifiedPeaks.tsv` — reports how many of those PSM-based peaks the untargeted run
//! independently rediscovered (mass + RT match).
//!
//! Usage:
//!   cargo run --release --example detect_features_tsv -- <spectra_file> <out.tsv> [reference.tsv]

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

use flashlfq_core::deconvolution::{ClassicDeconvolutionParameters, Polarity};
use flashlfq_core::feature_refinement::{
    four_way_decon, four_way_decon_detector_anchor, four_way_decon_gated, refine_feature,
    refine_feature_censored, refine_feature_shift, refine_feature_shift_neighbor,
    resolve_charge_state_consensus, DeconView, FourWayDecon, NeighborIndex, RefinedFeature,
    ResolvedFeature,
};
use flashlfq_core::isotopic_envelope::{mass_to_mz_f64, C13_MINUS_C12};
use flashlfq_core::peak_indexing::{read_ms1_scans, PeakIndexingEngine, PeakKey};
use flashlfq_core::trace_kernel::{
    detect_features, estimate_noise_floor, median_ms1_scan_spacing_minutes, CombWeightModel,
    DetectedFeature, ScoreModel, TraceKernelParameters, FWHM_TO_SIGMA,
};

/// Neighbour isotope-m/z positions in feature `f`'s window that are **not** on its own grid — the peaks
/// to mask from its shift fit (`NEIGHBOR_REFINE`). Own grid = `mono_mz + k·spacing`, `k ∈ [-1, kmax+2]`.
/// Only neighbours at least `min_ratio ×` this feature's intensity contribute (defer to stronger
/// species only). `idx` may be built from detected OR refined features; querying by `f.apex_rt` avoids
/// any self-index alignment, and `min_ratio > 1` guarantees a feature never masks its own peaks.
fn neighbor_mask_for(idx: &NeighborIndex, f: &DetectedFeature, min_ratio: f64) -> Vec<f64> {
    const GRID_PPM: f64 = 15.0;
    let spacing = C13_MINUS_C12 / f.charge.max(1) as f64;
    let win_min = (f.mono_mz - 1.5).max(0.0);
    let win_max = f.mono_mz + (f.num_isotopes_observed as f64 + 3.0) * spacing + 1.0;
    let kmax = f.num_isotopes_observed as i32 + 2;
    let min_intensity = f.summed_intensity * min_ratio;
    idx.forbidden_positions_query(f.apex_rt, win_min, win_max, min_intensity)
        .into_iter()
        .filter(|&p| {
            !(-1..=kmax).any(|k| {
                let own = f.mono_mz + k as f64 * spacing;
                (p - own).abs() / p * 1e6 <= GRID_PPM
            })
        })
        .collect()
}

/// A growing, RT-bucketed set of *already-refined* feature grids, used by the iterative
/// (`NEIGHBOR_REFINE=iterative`) pass: features are refined in descending score order, and each locks
/// its corrected isotope grid here so subsequent (lower-scoring) features can mask its peaks.
struct LockedGrids {
    /// RT bin (`floor(apex_rt / rt_tol)`) → `(apex_rt, mono_mz, spacing, kmax)`.
    buckets: std::collections::HashMap<i64, Vec<(f64, f64, f64, i32)>>,
    rt_tol: f64,
}

impl LockedGrids {
    fn new(rt_tol: f64) -> Self {
        LockedGrids { buckets: std::collections::HashMap::new(), rt_tol }
    }
    fn bin(&self, rt: f64) -> i64 {
        (rt / self.rt_tol).floor() as i64
    }
    fn add(&mut self, apex_rt: f64, mono_mz: f64, charge: i32, kmax: i32) {
        let spacing = C13_MINUS_C12 / charge.max(1) as f64;
        let b = self.bin(apex_rt);
        self.buckets.entry(b).or_default().push((apex_rt, mono_mz, spacing, kmax));
    }
    /// Ascending isotope m/z of locked features co-eluting with `apex_rt` and in `[win_min, win_max]`.
    fn positions(&self, apex_rt: f64, win_min: f64, win_max: f64) -> Vec<f64> {
        let b = self.bin(apex_rt);
        let mut out = Vec::new();
        for bb in (b - 1)..=(b + 1) {
            let Some(v) = self.buckets.get(&bb) else { continue };
            for &(rt, mono_mz, spacing, kmax) in v {
                if (rt - apex_rt).abs() > self.rt_tol {
                    continue;
                }
                for k in 0..=kmax {
                    let m = mono_mz + k as f64 * spacing;
                    if m < win_min {
                        continue;
                    }
                    if m > win_max {
                        break;
                    }
                    out.push(m);
                }
            }
        }
        out.sort_by(f64::total_cmp);
        out
    }
}

/// Off-own-grid filter shared by the neighbour-mask paths: keep only positions not within `GRID_PPM`
/// of `f`'s own isotope grid (`mono_mz + k·spacing`, `k ∈ [-1, kmax+2]`).
fn filter_off_own_grid(raw: Vec<f64>, f: &DetectedFeature) -> Vec<f64> {
    const GRID_PPM: f64 = 15.0;
    let spacing = C13_MINUS_C12 / f.charge.max(1) as f64;
    let kmax = f.num_isotopes_observed as i32 + 2;
    raw.into_iter()
        .filter(|&p| {
            !(-1..=kmax).any(|k| {
                let own = f.mono_mz + k as f64 * spacing;
                (p - own).abs() / p * 1e6 <= GRID_PPM
            })
        })
        .collect()
}

/// Derives a sibling output path from the final path: `out.tsv` + tag `detected` -> `out.detected.tsv`.
fn sibling(out: &str, tag: &str) -> String {
    match out.strip_suffix(".tsv") {
        Some(stem) => format!("{stem}.{tag}.tsv"),
        None => format!("{out}.{tag}.tsv"),
    }
}

/// Opens an output writer, tolerating a locked target (e.g. the file is open in Excel for manual
/// validation): on failure it falls back to `<path>.new` and warns, rather than panicking and
/// discarding the whole run. Returns `None` only if even the fallback cannot be created.
fn open_out(path: &str) -> Option<BufWriter<File>> {
    match File::create(path) {
        Ok(f) => Some(BufWriter::new(f)),
        Err(e) => {
            let alt = format!("{path}.new");
            eprintln!("  WARN: could not write {path} ({e}) — is it open? writing {alt} instead");
            match File::create(&alt) {
                Ok(f) => Some(BufWriter::new(f)),
                Err(e2) => {
                    eprintln!("  WARN: fallback {alt} also failed ({e2}); skipping this file");
                    None
                }
            }
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: detect_features_tsv <spectra_file> <out.tsv> [reference.tsv]");
        std::process::exit(2);
    }
    let spectra_path = &args[1];
    let out_path = &args[2];
    let reference_path = args.get(3);

    let detected_path = sibling(out_path, "detected");
    let refined_path = sibling(out_path, "refined");
    let log_path = match out_path.strip_suffix(".tsv") {
        Some(stem) => format!("{stem}.log"),
        None => format!("{out_path}.log"),
    };
    eprintln!("output files:");
    eprintln!("  detected (pre-refinement): {detected_path}");
    eprintln!("  refined (post-decon):      {refined_path}");
    eprintln!("  resolved (final):          {out_path}");
    eprintln!("  timing log:                {log_path}");

    // Per-step wall-clock timings, reported to the chat and written to the log file at the end.
    let run_start = Instant::now();
    let mut timings: Vec<(String, f64)> = Vec::new();

    // --- read + index --------------------------------------------------------------------------
    let t0 = Instant::now();
    eprintln!("reading MS1 scans from {spectra_path} ...");
    let scans = read_ms1_scans(spectra_path).expect("failed to read spectra file");
    let engine = PeakIndexingEngine::index_peaks(&scans).expect("no indexable MS1 peaks");
    let n_peaks: usize = scans.iter().map(|s| s.mz.len()).sum();
    let total_intensity: f64 = scans.iter().flat_map(|s| s.intensity.iter()).sum();
    let read_dur = t0.elapsed();
    timings.push(("read + index".into(), read_dur.as_secs_f64()));
    eprintln!(
        "  {} MS1 scans, {} peaks, ΣTIC {:.3e}, median scan spacing {:.4} min  ({:.1?})",
        scans.len(),
        n_peaks,
        total_intensity,
        median_ms1_scan_spacing_minutes(engine.scan_info()),
        read_dur
    );

    // --- detect --------------------------------------------------------------------------------
    // Comb-weight model selectable via COMB_MODEL=averagine|poisson (default averagine) for benchmarking.
    let weight_model = match std::env::var("COMB_MODEL").as_deref() {
        Ok("poisson") => CombWeightModel::Poisson,
        _ => CombWeightModel::Averagine,
    };
    // Coverage target (fraction of ΣTIC to explain) selectable via COVERAGE_TARGET=0.80|0.90|0.99…
    let coverage_target = std::env::var("COVERAGE_TARGET")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.90);
    // Assumed chromatographic FWHM (seconds) that sets the RT Gaussian σ and matched-filter window,
    // selectable via ASSUMED_FWHM_SEC (default 36). Real CA/Lumos peaks are ~3-20 s wide, so smaller
    // values narrow the window to the data and reduce broad-hypothesis interference.
    let assumed_fwhm_sec = std::env::var("ASSUMED_FWHM_SEC")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(36.0);
    // Claim-extent trace knobs (Change A): missed-scan tolerance and the RT half-width guard for the
    // XIC that follows a real elution beyond the ~2σ scoring window. TRACE_MAX_HALF_WIDTH_SEC is in
    // seconds; default 30 s (0.5 min).
    let trace_missed = std::env::var("TRACE_MISSED_SCANS")
        .ok()
        .and_then(|s| s.parse::<i32>().ok())
        .unwrap_or(1);
    let trace_half_width_min = std::env::var("TRACE_MAX_HALF_WIDTH_SEC")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .map(|s| s / 60.0)
        .unwrap_or(0.5);
    // Chromatographic-persistence gate (Change A): reject features whose traced extent spans fewer
    // than this many distinct scans. Default 2 drops single-scan noise doublets; MIN_FEATURE_SCANS=1
    // disables it.
    let min_feature_scans = std::env::var("MIN_FEATURE_SCANS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(2);
    // Seed intensity floor (bounds detection cost). Lower it to reach higher coverage (the default
    // 1000 exhausts seeds ~90% ΣTIC on CA/Lumos before the coverage target is hit).
    let min_seed_intensity = std::env::var("MIN_SEED_INTENSITY")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1000.0);
    // Score model (Change B): SCORE_MODEL=normalized selects the noise-floor-truncated normalised
    // correlation (default raw sum). NOISE_PCT is the percentile of peak intensity used as η (default 5).
    let score_model = match std::env::var("SCORE_MODEL").as_deref() {
        Ok("normalized") | Ok("normalised") => ScoreModel::NormalizedNoiseFloor,
        _ => ScoreModel::RawSum,
    };
    let noise_pct = std::env::var("NOISE_PCT")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(5.0);
    let noise_floor = if matches!(score_model, ScoreModel::NormalizedNoiseFloor) {
        estimate_noise_floor(&engine, noise_pct)
    } else {
        0.0
    };
    // Change B variant knobs: AMP_SEED=1 uses the seed intensity as the apex amplitude A (vs the
    // least-squares fit); SCORE_COSINE=1 divides additionally by ‖I‖ for a bounded cosine shape-fit.
    let score_use_seed_amplitude = matches!(std::env::var("AMP_SEED").as_deref(), Ok("1") | Ok("true"));
    let score_cosine = matches!(std::env::var("SCORE_COSINE").as_deref(), Ok("1") | Ok("true"));
    let base = TraceKernelParameters {
        ppm_tolerance: 10.0,
        min_seed_intensity,
        coverage_target,
        weight_model,
        score_model,
        noise_floor,
        score_use_seed_amplitude,
        score_cosine,
        trace_missed_scans_allowed: trace_missed,
        trace_max_half_width_minutes: trace_half_width_min,
        min_feature_scans,
        ..TraceKernelParameters::default()
    };
    eprintln!(
        "score model: {:?}  (η = {:.1} @ p{:.0}, A = {}, {})",
        score_model,
        noise_floor,
        noise_pct,
        if score_use_seed_amplitude { "seed" } else { "least-squares" },
        if score_cosine { "cosine" } else { "template-norm" }
    );
    // σ from the data by default (Change A: measure FWHM via XIC half-max, ASSUMED_FWHM_SEC is the
    // fallback). FIXED_SIGMA=1 forces the assumed-FWHM path for A/B comparison.
    let data_driven_sigma = std::env::var("FIXED_SIGMA").is_err();
    let params = if data_driven_sigma {
        base.with_rt_from_index(&engine, assumed_fwhm_sec)
    } else {
        base.with_rt_from_scans(engine.scan_info(), assumed_fwhm_sec)
    };
    if data_driven_sigma {
        eprintln!(
            "σ source: data-driven (measured FWHM via XIC half-max; fallback {assumed_fwhm_sec} s) \
             → {:.2} s FWHM",
            params.rt_sigma_minutes * FWHM_TO_SIGMA * 60.0
        );
    } else {
        eprintln!("σ source: assumed FWHM {assumed_fwhm_sec} s (FIXED_SIGMA)");
    }
    eprintln!("comb weight model: {weight_model:?}");
    eprintln!(
        "detecting (charge {}..={}, {} ppm, σ_rt {:.4} min, ±{} scans, trace: {} missed / {:.0} s half-width, min {} scans, seed floor {:.0}, coverage {:.0}%) ...",
        params.min_charge,
        params.max_charge,
        params.ppm_tolerance,
        params.rt_sigma_minutes,
        params.half_window_scans,
        params.trace_missed_scans_allowed,
        params.trace_max_half_width_minutes * 60.0,
        params.min_feature_scans,
        params.min_seed_intensity,
        params.coverage_target * 100.0
    );
    let t1 = Instant::now();
    let detected = detect_features(&engine, &params);
    let detect_dur = t1.elapsed();
    timings.push(("detect".into(), detect_dur.as_secs_f64()));
    let detected_intensity: f64 = detected.iter().map(|f| f.summed_intensity).sum();
    eprintln!(
        "  {} features detected, explained {:.1}% of ΣTIC  ({:.1?})",
        detected.len(),
        100.0 * detected_intensity / total_intensity,
        detect_dur
    );
    write_detected_tsv(&detected_path, &detected);
    eprintln!("  wrote {} detected features -> {detected_path}", detected.len());

    // Escape hatch for diagnostics: skip the (potentially intractable at huge feature counts)
    // refine + O(n^2) charge-consensus and stop after detection.
    if std::env::var("DETECT_ONLY").is_ok() {
        eprintln!("DETECT_ONLY set — skipping refine/resolve/compare after {} detected features.", detected.len());
        return;
    }

    // --- refine --------------------------------------------------------------------------------
    let avg = flashlfq_core::spectral_averaging::SpectralAveragingParameters::default();
    let decon = ClassicDeconvolutionParameters::new(
        params.min_charge,
        params.max_charge,
        10.0,
        3.0,
        Polarity::Positive,
    );
    // CENSOR_CLAIMED=1 removes peaks claimed by OTHER (stronger, already-assigned) features from
    // each feature's composite before deconvolution — a subtractive-decon experiment. Since the
    // detector claims greedily tallest-first, this hands weaker features a window with co-eluting
    // interferents removed.
    let censor_claimed = std::env::var("CENSOR_CLAIMED").is_ok();
    let all_claimed: HashSet<PeakKey> = if censor_claimed {
        detected
            .iter()
            .flat_map(|f| f.peaks.iter().map(|p| p.key()))
            .collect()
    } else {
        HashSet::new()
    };
    if censor_claimed {
        eprintln!(
            "  CENSOR_CLAIMED: subtracting {} claimed peaks from other features' decon windows",
            all_claimed.len()
        );
    }

    // REFINE_METHOD selects the deconvolution used to place the refined monoisotope:
    //   classic (default) | shift_composite | shift_apex  (detector-anchored FlashLFQ-style shift).
    // Default: the detector-anchored shift decon on the apex scan (REFINE_METHOD=classic to opt out
    // back to the parity-locked classic deconvolution; shift_composite selects the averaged composite).
    let refine_method = std::env::var("REFINE_METHOD").unwrap_or_else(|_| "shift_apex".to_string());
    let use_shift_apex = refine_method == "shift_apex";
    let use_shift = use_shift_apex || refine_method == "shift_composite";
    // Charge re-selection by the envelope-fit cosine (fit + explained + completeness) defaults ON for
    // the shift methods; disable with RECHARGE=0. Recovers charge-halved features (the light-z2 class).
    let recharge = std::env::var("RECHARGE")
        .map(|v| v != "0" && !v.eq_ignore_ascii_case("false"))
        .unwrap_or(true);
    eprintln!("  refine method: {refine_method}{}", if recharge && use_shift { " + recharge (envelope-fit cosine)" } else { "" });

    // NEIGHBOR_REFINE=1: mask co-eluting neighbours' peaks (off this feature's own grid) from the
    // charge-selection fit, walk-back and score — so a low-scoring feature in a crowded window is judged
    // on the signal plausibly its own. Experiment path; builds a NeighborIndex over the detections.
    // NEIGHBOR_REFINE: mask co-eluting neighbours' peaks from the shift fit. Value picks the neighbour
    // context: "detected"/"1" = raw detections; "refined" = a two-pass build (refine once, then mask
    // against those corrected placements). NEIGHBOR_MIN_RATIO (default 5.0) = only defer to neighbours
    // at least that many times more intense — a weak feature in a strong neighbour's shadow.
    let neighbor_mode = std::env::var("NEIGHBOR_REFINE").ok();
    // "iterative": refine in descending-score order, each feature locking its corrected grid so later
    // (lower-scoring) features mask its peaks. The confident features claim their signal first.
    let iterative = use_shift && neighbor_mode.as_deref() == Some("iterative");
    let neighbor_refine = use_shift && !iterative && neighbor_mode.is_some();
    let neighbor_refined_context = neighbor_mode.as_deref() == Some("refined");
    let neighbor_min_ratio = std::env::var("NEIGHBOR_MIN_RATIO")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(5.0);
    let neighbor_idx = if neighbor_refine {
        let src = if neighbor_refined_context {
            // Pass 1: plain refine, then index the corrected placements for the masked pass below.
            let pass1: Vec<RefinedFeature> = detected
                .iter()
                .filter_map(|f| refine_feature_shift(f, &scans, &avg, 20.0, use_shift_apex, recharge))
                .collect();
            eprintln!(
                "  NEIGHBOR_REFINE=refined: two-pass, masking neighbours >= {neighbor_min_ratio}x (from {} refined)",
                pass1.len()
            );
            NeighborIndex::build_from_refined(&pass1, 0.05)
        } else {
            eprintln!(
                "  NEIGHBOR_REFINE: masking peaks of co-eluting detected neighbours >= {neighbor_min_ratio}x this feature's intensity"
            );
            NeighborIndex::build(&detected, 0.05)
        };
        Some(src)
    } else {
        None
    };

    let t2 = Instant::now();
    let mut refined: Vec<RefinedFeature> = Vec::with_capacity(detected.len());
    let progress_every = 1000usize;
    if iterative {
        // Only features whose fit clears this bar lock their grid into the mask — "higher-scoring"
        // is not enough; a mediocre-but-higher feature is still uncertain and would add collateral.
        let lock_min_score = std::env::var("NEIGHBOR_LOCK_MIN_SCORE")
            .ok()
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(0.9);
        eprintln!("  NEIGHBOR_REFINE=iterative: score-ordered; only features with fit >= {lock_min_score} lock their grids");
        // Pass 1: plain refine to get each feature's initial score.
        let init: Vec<Option<RefinedFeature>> = detected
            .iter()
            .map(|f| refine_feature_shift(f, &scans, &avg, 20.0, use_shift_apex, recharge))
            .collect();
        // Process indices in descending initial score.
        let mut order: Vec<usize> = (0..detected.len()).filter(|&i| init[i].is_some()).collect();
        order.sort_by(|&a, &b| {
            init[b].as_ref().unwrap().decon_score.total_cmp(&init[a].as_ref().unwrap().decon_score)
        });
        let mut locked = LockedGrids::new(0.05);
        let mut out: Vec<Option<RefinedFeature>> = vec![None; detected.len()];
        let mut n_locked = 0usize;
        for (n, &i) in order.iter().enumerate() {
            let f = &detected[i];
            let spacing = C13_MINUS_C12 / f.charge.max(1) as f64;
            let win_min = (f.mono_mz - 1.5).max(0.0);
            let win_max = f.mono_mz + (f.num_isotopes_observed as f64 + 3.0) * spacing + 1.0;
            let mask = filter_off_own_grid(locked.positions(f.apex_rt, win_min, win_max), f);
            let r = refine_feature_shift_neighbor(f, &scans, &avg, 20.0, use_shift_apex, recharge, &mask)
                .or_else(|| init[i].clone());
            if let Some(rr) = &r {
                // Only confident features become mask sources for the lower-scoring ones that follow.
                if rr.decon_score >= lock_min_score {
                    let mono_mz = mass_to_mz_f64(rr.refined_monoisotopic_mass, rr.refined_charge);
                    let kmax = f.num_isotopes_observed as i32 + 2;
                    locked.add(rr.detected.apex_rt, mono_mz, rr.refined_charge, kmax);
                    n_locked += 1;
                }
            }
            out[i] = r;
            if (n + 1) % progress_every == 0 || n + 1 == order.len() {
                eprintln!("    iterative refined {}/{} ({n_locked} locked)  [{:?}]", n + 1, order.len(), t2.elapsed());
            }
        }
        refined = out.into_iter().flatten().collect();
    } else {
        for (i, f) in detected.iter().enumerate() {
            let r = if let Some(idx) = &neighbor_idx {
                let mask = neighbor_mask_for(idx, f, neighbor_min_ratio);
                refine_feature_shift_neighbor(f, &scans, &avg, 20.0, use_shift_apex, recharge, &mask)
            } else if use_shift {
                refine_feature_shift(f, &scans, &avg, 20.0, use_shift_apex, recharge)
            } else if censor_claimed {
                refine_feature_censored(f, &scans, &avg, &decon, &all_claimed)
            } else {
                refine_feature(f, &scans, &avg, &decon)
            };
            if let Some(r) = r {
                refined.push(r);
            }
            if (i + 1) % progress_every == 0 || i + 1 == detected.len() {
                eprintln!(
                    "    refined {}/{} ({} kept)  [{:?}]",
                    i + 1,
                    detected.len(),
                    refined.len(),
                    t2.elapsed()
                );
            }
        }
    }
    let refine_dur = t2.elapsed();
    timings.push(("refine".into(), refine_dur.as_secs_f64()));
    eprintln!(
        "  {} / {} features refined against averaged composites  ({:.1?})",
        refined.len(),
        detected.len(),
        refine_dur
    );
    write_refined_tsv(&refined_path, &refined);
    eprintln!("  wrote {} refined features -> {refined_path}", refined.len());

    // --- four-way decon comparator (opt-in diagnostic) -----------------------------------------
    // FOUR_WAY_DECON=1 runs the classic×shift on composite×apex comparator over every detected
    // feature, reporting how often the four monoisotope views disagree (→ candidates for advanced
    // multi-envelope decon) and writing a per-feature disagreements TSV. Gated because it roughly
    // doubles the decon cost; off by default so normal runs are unaffected.
    if std::env::var("FOUR_WAY_DECON").is_ok() {
        let t_fw = Instant::now();
        let fw_path = sibling(out_path, "disagreements");
        run_four_way(&detected, &scans, &avg, &decon, &fw_path);
        timings.push(("four-way decon".into(), t_fw.elapsed().as_secs_f64()));
    }

    // DIFF_REFINE=1 exports the features CLASSIC refine drops (no envelope) but the detector-anchored
    // SHIFT refine keeps — the completeness gain behind the +2.4% recall. Feeds make_gain_targets.py.
    if std::env::var("DIFF_REFINE").is_ok() {
        let dpath = sibling(out_path, "refinediff");
        run_refine_diff(&detected, &scans, &avg, &decon, &dpath);
    }

    // --- resolve charge-state consensus --------------------------------------------------------
    let t3 = Instant::now();
    let resolved = resolve_charge_state_consensus(&refined, 10.0, 0.1);
    let consensus_dur = t3.elapsed();
    timings.push(("charge-state consensus".into(), consensus_dur.as_secs_f64()));
    eprintln!(
        "  {} peptide-level features after charge-state consensus  ({:.1?})",
        resolved.len(),
        consensus_dur
    );

    let t4 = Instant::now();
    write_tsv(out_path, &resolved);
    let write_dur = t4.elapsed();
    timings.push(("write resolved".into(), write_dur.as_secs_f64()));
    eprintln!("wrote {} -> {}", resolved.len(), out_path);

    if let Some(ref_path) = reference_path {
        let t5 = Instant::now();
        compare_to_reference(ref_path, &resolved);
        timings.push(("compare".into(), t5.elapsed().as_secs_f64()));
    }

    // --- timing summary (chat + log file) ------------------------------------------------------
    let total = run_start.elapsed().as_secs_f64();
    let mut lines: Vec<String> = Vec::new();
    lines.push("=== timing summary ===".into());
    lines.push(format!("spectra file: {spectra_path}"));
    lines.push(format!(
        "{} MS1 scans, {} peaks, {} detected, {} refined, {} resolved features",
        scans.len(),
        n_peaks,
        detected.len(),
        refined.len(),
        resolved.len()
    ));
    for (name, secs) in &timings {
        lines.push(format!("  {name:<24} {secs:8.2} s  ({:4.1}%)", 100.0 * secs / total));
    }
    lines.push(format!("  {:<24} {total:8.2} s", "TOTAL"));

    for l in &lines {
        eprintln!("{l}");
    }
    // Append to the log file so repeated runs accumulate a history (never blocks the run).
    match std::fs::OpenOptions::new().create(true).append(true).open(&log_path) {
        Ok(mut f) => {
            for l in &lines {
                let _ = writeln!(f, "{l}");
            }
            let _ = writeln!(f);
            eprintln!("timing log appended to {log_path}");
        }
        Err(e) => eprintln!("  WARN: could not write timing log {log_path} ({e})"),
    }
}

/// Exports features that classic `refine_feature` drops (returns `None`) but detector-anchored
/// `refine_feature_shift` (apex) keeps — the completeness gain. Columns feed `make_gain_targets.py`
/// (which joins to the reference): shift-refined mono, apex RT, charge, most-abundant peak m/z.
fn run_refine_diff(
    detected: &[DetectedFeature],
    scans: &[flashlfq_core::peak_indexing::Scan],
    avg: &flashlfq_core::spectral_averaging::SpectralAveragingParameters,
    decon: &ClassicDeconvolutionParameters,
    path: &str,
) {
    let mut w = match open_out(path) {
        Some(w) => w,
        None => return,
    };
    writeln!(w, "Mono\tRT\tCharge\tPk MZ\tSummed Intensity").unwrap();
    let mut n_gain = 0usize;
    for f in detected {
        // Classic kept it → not a gain.
        if refine_feature(f, scans, avg, decon).is_some() {
            continue;
        }
        if let Some(r) = refine_feature_shift(f, scans, avg, 20.0, true, false) {
            n_gain += 1;
            let pk = f
                .peaks
                .iter()
                .max_by(|a, b| a.intensity.total_cmp(&b.intensity))
                .map(|p| p.m() as f64)
                .unwrap_or(f.mono_mz);
            writeln!(
                w,
                "{:.5}\t{:.4}\t{}\t{:.5}\t{:.4e}",
                r.refined_monoisotopic_mass, f.apex_rt, f.charge, pk, f.summed_intensity
            )
            .unwrap();
        }
    }
    let _ = w.flush();
    eprintln!(
        "  DIFF_REFINE: {n_gain} features shift-kept but classic-dropped -> {path}"
    );
}

/// Runs the four-way decon comparator over every detected feature, reports the disagreement rate to
/// the chat, and writes a per-disagreeing-feature TSV (the candidates for advanced multi-envelope
/// decon). `decon` uses the pipeline's classic parameters; the shift views use a 20 ppm match tol.
fn run_four_way(
    detected: &[DetectedFeature],
    scans: &[flashlfq_core::peak_indexing::Scan],
    avg: &flashlfq_core::spectral_averaging::SpectralAveragingParameters,
    decon: &ClassicDeconvolutionParameters,
    path: &str,
) {
    let short = |v: DeconView| match v {
        DeconView::ClassicComposite => "cc",
        DeconView::ClassicApex => "ca",
        DeconView::ShiftComposite => "sc",
        DeconView::ShiftApex => "sa",
    };

    let mut w = open_out(path);
    if let Some(w) = w.as_mut() {
        writeln!(
            w,
            "Detector Mono\tCharge\tApex RT\tConsensus k\tNum Verdicts\tViews (view:k:conf)"
        )
        .unwrap();
    }

    // Full per-feature verdict table: each view's monoisotopic mass (blank if the view produced
    // nothing). Lets the within-method (composite-vs-apex) and ground-truth analysis run in Python.
    let vpath = path.replace(".disagreements.tsv", ".verdicts.tsv");
    let mut vw = open_out(&vpath);
    if let Some(vw) = vw.as_mut() {
        writeln!(vw, "Detector Mono\tCharge\tApex RT\tcc_mono\tca_mono\tsc_mono\tsa_mono").unwrap();
    }
    let mono_of = |fw: &FourWayDecon, view: DeconView| -> Option<f64> {
        fw.verdicts
            .iter()
            .find(|v| v.view == view)
            .map(|v| v.monoisotopic_mass)
    };
    let fmt = |m: Option<f64>| m.map(|x| format!("{x:.5}")).unwrap_or_default();
    // Within-method composite-vs-apex agreement (same integer ¹³C offset), counted where both views exist.
    let mut classic_agree = 0usize;
    let mut classic_disagree = 0usize;
    let mut shift_agree = 0usize;
    let mut shift_disagree = 0usize;
    let same_k = |a: f64, b: f64| ((a - b) / C13_MINUS_C12).round() == 0.0;

    // NEIGHBOR_AWARE=1 gates the shift-decon anchor: peaks belonging to a co-eluting already-detected
    // feature (its predicted isotope grid) are excluded as anchor candidates, so shift-decon can't
    // "distant grab" a stronger neighbour several isotopes away. Builds an RT-bucketed neighbour index.
    // DETECTOR_ANCHOR=1 takes precedence: anchor the shift views on the detector's own most-abundant
    // claimed peak (cannot grab any foreign peak). NEIGHBOR_AWARE=1 gates against co-eluting neighbours.
    let detector_anchor = std::env::var("DETECTOR_ANCHOR").is_ok();
    let neighbor_aware = !detector_anchor && std::env::var("NEIGHBOR_AWARE").is_ok();
    let neighbors = if neighbor_aware {
        eprintln!("  NEIGHBOR_AWARE: gating shift anchors against co-eluting detected features");
        Some(NeighborIndex::build(detected, 0.05))
    } else {
        None
    };
    if detector_anchor {
        eprintln!("  DETECTOR_ANCHOR: shift views anchor on the detector's most-abundant claimed peak");
    }
    const GRID_PPM: f64 = 15.0;

    let mut total = 0usize; // features that produced >=1 verdict
    let mut with_verdicts_hist = [0usize; 5]; // count by number of verdicts (0..=4)
    let mut unanimous = 0usize;
    let mut disagreed = 0usize;
    for (i, f) in detected.iter().enumerate() {
        let fw: FourWayDecon = if detector_anchor {
            four_way_decon_detector_anchor(f, scans, avg, decon, 20.0)
        } else {
            match &neighbors {
                Some(idx) => {
                    // Same shift window the comparator uses, to gather the neighbours' forbidden teeth.
                    let spacing = C13_MINUS_C12 / f.charge.max(1) as f64;
                    let win_min = (f.mono_mz - 1.5).max(0.0);
                    let win_max = f.mono_mz + (f.num_isotopes_observed as f64 + 3.0) * spacing + 1.0;
                    let forbidden = idx.forbidden_positions(i, win_min, win_max);
                    four_way_decon_gated(f, scans, avg, decon, 20.0, &forbidden, GRID_PPM)
                }
                None => four_way_decon(f, scans, avg, decon, 20.0),
            }
        };
        with_verdicts_hist[fw.verdicts.len().min(4)] += 1;
        if fw.verdicts.is_empty() {
            continue;
        }
        total += 1;

        // Per-feature verdict row + within-method (composite vs apex) agreement.
        let (cc, ca, sc, sa) = (
            mono_of(&fw, DeconView::ClassicComposite),
            mono_of(&fw, DeconView::ClassicApex),
            mono_of(&fw, DeconView::ShiftComposite),
            mono_of(&fw, DeconView::ShiftApex),
        );
        if let Some(vw) = vw.as_mut() {
            writeln!(
                vw,
                "{:.5}\t{}\t{:.4}\t{}\t{}\t{}\t{}",
                fw.detector_mono,
                fw.charge,
                f.apex_rt,
                fmt(cc),
                fmt(ca),
                fmt(sc),
                fmt(sa)
            )
            .unwrap();
        }
        if let (Some(a), Some(b)) = (cc, ca) {
            if same_k(a, b) {
                classic_agree += 1;
            } else {
                classic_disagree += 1;
            }
        }
        if let (Some(a), Some(b)) = (sc, sa) {
            if same_k(a, b) {
                shift_agree += 1;
            } else {
                shift_disagree += 1;
            }
        }
        if fw.needs_advanced {
            disagreed += 1;
            if let Some(w) = w.as_mut() {
                let views: Vec<String> = fw
                    .verdicts
                    .iter()
                    .map(|v| format!("{}:{:+}:{}", short(v.view), v.offset_k, if v.confident { 1 } else { 0 }))
                    .collect();
                writeln!(
                    w,
                    "{:.5}\t{}\t{:.4}\t{}\t{}\t{}",
                    fw.detector_mono,
                    fw.charge,
                    f.apex_rt,
                    fw.consensus_k.map(|k| k.to_string()).unwrap_or_default(),
                    fw.verdicts.len(),
                    views.join(" ")
                )
                .unwrap();
            }
        } else {
            unanimous += 1;
        }
    }
    if let Some(mut w) = w {
        let _ = w.flush();
    }
    if let Some(mut vw) = vw {
        let _ = vw.flush();
        eprintln!("  wrote per-feature verdicts -> {vpath}");
    }

    let cpair = classic_agree + classic_disagree;
    let spair = shift_agree + shift_disagree;
    eprintln!("\n=== within-method composite-vs-apex agreement ===");
    eprintln!(
        "  classic (cc vs ca): {} agree / {} disagree  ({:.1}% agree of {} with both)",
        classic_agree, classic_disagree, 100.0 * classic_agree as f64 / cpair.max(1) as f64, cpair
    );
    eprintln!(
        "  shift   (sc vs sa): {} agree / {} disagree  ({:.1}% agree of {} with both)",
        shift_agree, shift_disagree, 100.0 * shift_agree as f64 / spair.max(1) as f64, spair
    );

    eprintln!("\n=== four-way decon comparator ===");
    eprintln!("  detected features: {}", detected.len());
    eprintln!(
        "  produced >=1 verdict: {} (verdict-count histogram [0..4]: {:?})",
        total, with_verdicts_hist
    );
    eprintln!(
        "  unanimous (confident placement): {} ({:.1}%)",
        unanimous,
        100.0 * unanimous as f64 / total.max(1) as f64
    );
    eprintln!(
        "  disagreed (→ advanced multi-envelope): {} ({:.1}%)  → {path}",
        disagreed,
        100.0 * disagreed as f64 / total.max(1) as f64
    );
}

/// Writes the raw detected features (pre-refinement, straight from the trace kernel) to a TSV,
/// sorted by summed intensity descending. Lets you inspect what the detector alone produced.
fn write_detected_tsv(path: &str, detected: &[DetectedFeature]) {
    let mut rows: Vec<&DetectedFeature> = detected.iter().collect();
    rows.sort_by(|a, b| b.summed_intensity.total_cmp(&a.summed_intensity));
    let mut w = match open_out(path) {
        Some(w) => w,
        None => return,
    };
    writeln!(
        w,
        "Monoisotopic Mass\tCharge\tMono m/z\tApex RT\tRT Start\tRT End\tSummed Intensity\t\
         Detector Score\tNum Isotopes\tNum Peaks"
    )
    .unwrap();
    for d in rows {
        writeln!(
            w,
            "{:.5}\t{}\t{:.5}\t{:.4}\t{:.4}\t{:.4}\t{:.4e}\t{:.4e}\t{}\t{}",
            d.monoisotopic_mass,
            d.charge,
            d.mono_mz,
            d.apex_rt,
            d.start_rt,
            d.end_rt,
            d.summed_intensity,
            d.score,
            d.num_isotopes_observed,
            d.peaks.len()
        )
        .unwrap();
    }
    w.flush().unwrap();
}

/// Writes the refined features (post composite-deconvolution, pre charge-consensus) to a TSV,
/// sorted by summed intensity descending.
fn write_refined_tsv(path: &str, refined: &[RefinedFeature]) {
    let mut rows: Vec<&RefinedFeature> = refined.iter().collect();
    rows.sort_by(|a, b| {
        b.detected
            .summed_intensity
            .total_cmp(&a.detected.summed_intensity)
    });
    let mut w = match open_out(path) {
        Some(w) => w,
        None => return,
    };
    writeln!(
        w,
        "Refined Monoisotopic Mass\tCharge\tApex RT\tSummed Intensity\tDecon Score\t\
         Num Candidate Masses\tDetector Mono Mass"
    )
    .unwrap();
    for r in rows {
        writeln!(
            w,
            "{:.5}\t{}\t{:.4}\t{:.4e}\t{:.4e}\t{}\t{:.5}",
            r.refined_monoisotopic_mass,
            r.refined_charge,
            r.detected.apex_rt,
            r.detected.summed_intensity,
            r.decon_score,
            r.candidate_masses.len(),
            r.detected.monoisotopic_mass
        )
        .unwrap();
    }
    w.flush().unwrap();
}

/// Writes the resolved features to a human-readable TSV, sorted by summed intensity descending.
fn write_tsv(path: &str, resolved: &[ResolvedFeature]) {
    let mut rows: Vec<&ResolvedFeature> = resolved.iter().collect();
    rows.sort_by(|a, b| b.summed_intensity.total_cmp(&a.summed_intensity));

    let mut w = match open_out(path) {
        Some(w) => w,
        None => return,
    };
    // Columns lead with the observed-feature answer — the detected m/z, the RT extent, and every
    // charge state seen — then follow with the derived mass/intensity fields.
    //   `Detected m/z (primary)` is the tallest observed isotope-peak m/z of the tallest charge
    //     (the detector seed) — the direct analogue of base FlashLFQ's observed `Peak MZ`, which for
    //     heavier peptides sits ~1 ¹³C step above the monoisotope.
    //   `Per-Charge Detected m/z` lists that same observed m/z for EVERY detected charge, as
    //     `z<charge>:<m/z>` pairs (ascending by charge) — so a peptide seen at z2 and z3 shows both.
    //   `Mono m/z (primary)` is the monoisotopic-peak m/z at the primary charge (derived from the
    //     consensus neutral mass), for reference alongside the observed detected m/z.
    writeln!(
        w,
        "Detected m/z (primary)\tRT Start\tRT Apex\tRT End\tCharge States\tPer-Charge Detected m/z\t\
         Num Charge States\tPrimary Charge\tMonoisotopic Mass\tMono m/z (primary)\tSummed Intensity\t\
         Cross-Charge Support\tNum Members"
    )
    .unwrap();
    // Most-abundant observed isotope-peak m/z of one member's detection (its tallest claimed peak).
    let member_detected_mz = |m: &RefinedFeature| -> Option<f64> {
        m.detected
            .peaks
            .iter()
            .max_by(|a, b| a.intensity.total_cmp(&b.intensity))
            .map(|p| p.m() as f64)
    };
    for r in rows {
        // Primary member = the tallest member; its charge and its seed (tallest) peak drive the m/z.
        let primary_member = r.members.iter().max_by(|a, b| {
            a.detected
                .summed_intensity
                .total_cmp(&b.detected.summed_intensity)
        });
        let primary_charge = primary_member.map(|m| m.refined_charge).unwrap_or(0);
        let mono_mz = if primary_charge != 0 {
            mass_to_mz_f64(r.monoisotopic_mass, primary_charge)
        } else {
            0.0
        };
        let detected_mz = primary_member.and_then(member_detected_mz).unwrap_or(0.0);
        // Observed detected m/z per charge state: for each detected charge, the tallest member of
        // that charge and its most-abundant peak m/z, ascending by charge (`z2:497.2584;z3:331.8416`).
        let per_charge_mz: Vec<String> = r
            .charge_states
            .iter()
            .map(|&z| {
                let mz = r
                    .members
                    .iter()
                    .filter(|m| m.refined_charge == z)
                    .max_by(|a, b| {
                        a.detected
                            .summed_intensity
                            .total_cmp(&b.detected.summed_intensity)
                    })
                    .and_then(member_detected_mz)
                    .unwrap_or(0.0);
                format!("z{z}:{mz:.4}")
            })
            .collect();
        let charges: Vec<String> = r.charge_states.iter().map(|c| c.to_string()).collect();
        writeln!(
            w,
            "{:.5}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{:.5}\t{:.5}\t{:.4e}\t{}\t{}",
            detected_mz,
            r.start_rt,
            r.apex_rt,
            r.end_rt,
            charges.join(";"),
            per_charge_mz.join(";"),
            r.charge_states.len(),
            primary_charge,
            r.monoisotopic_mass,
            mono_mz,
            r.summed_intensity,
            r.cross_charge_support,
            r.members.len()
        )
        .unwrap();
    }
    w.flush().unwrap();
}

/// Loads the base-FlashLFQ `AllQuantifiedPeaks.tsv` and reports how many of its peaks a resolved
/// feature independently matches by monoisotopic mass (±20 ppm) and apex RT (±0.3 min). The
/// reference is the calibrated file, so a small mass/RT offset vs. the raw is expected — tolerances
/// are generous accordingly.
fn compare_to_reference(path: &str, resolved: &[ResolvedFeature]) {
    let text = match std::fs::read_to_string(path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("could not read reference {path}: {e}");
            return;
        }
    };
    let mut lines = text.lines();
    let header = match lines.next() {
        Some(h) => h,
        None => return,
    };
    let cols: HashMap<&str, usize> = header.split('\t').enumerate().map(|(i, c)| (c, i)).collect();
    let mass_i = cols["Peptide Monoisotopic Mass"];
    let rt_i = cols["Peak RT Apex"];
    let charge_i = cols["Peak Charge"];
    let seq_i = cols["Full Sequence"];

    // Reference rows: (mass, apex_rt, charge, sequence). Skip malformed / empty-mass rows.
    struct RefPeak {
        mass: f64,
        rt: f64,
        charge: i32,
        seq: String,
    }
    let mut refs: Vec<RefPeak> = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let (m, rt) = match (
            fields.get(mass_i).and_then(|s| s.parse::<f64>().ok()),
            fields.get(rt_i).and_then(|s| s.parse::<f64>().ok()),
        ) {
            (Some(m), Some(rt)) => (m, rt),
            _ => continue,
        };
        let charge = fields
            .get(charge_i)
            .and_then(|s| s.parse::<i32>().ok())
            .unwrap_or(0);
        let seq = fields.get(seq_i).map(|s| s.to_string()).unwrap_or_default();
        refs.push(RefPeak { mass: m, rt, charge, seq });
    }

    const MASS_PPM: f64 = 20.0;
    const RT_MIN: f64 = 0.3;

    let mut matched = 0usize;
    let mut matched_with_charge = 0usize;
    let mut unmatched_examples: Vec<String> = Vec::new();
    for rp in &refs {
        let hit = resolved.iter().find(|f| {
            (f.monoisotopic_mass - rp.mass).abs() / rp.mass * 1e6 <= MASS_PPM
                && (f.apex_rt - rp.rt).abs() <= RT_MIN
        });
        match hit {
            Some(f) => {
                matched += 1;
                if f.charge_states.contains(&rp.charge) {
                    matched_with_charge += 1;
                }
            }
            None => {
                if unmatched_examples.len() < 15 {
                    unmatched_examples.push(format!(
                        "{:.4} Da  RT {:.3}  z{}  {}",
                        rp.mass, rp.rt, rp.charge, rp.seq
                    ));
                }
            }
        }
    }

    eprintln!("\n=== comparison to base FlashLFQ ({}) ===", path);
    eprintln!("  reference PSM-based peaks: {}", refs.len());
    eprintln!("  resolved untargeted features: {}", resolved.len());
    eprintln!(
        "  reference peaks rediscovered (±{:.0} ppm mass, ±{:.1} min RT): {} / {}  ({:.1}%)",
        MASS_PPM,
        RT_MIN,
        matched,
        refs.len(),
        100.0 * matched as f64 / refs.len().max(1) as f64
    );
    eprintln!(
        "    ...of which the charge state also matched: {} ({:.1}%)",
        matched_with_charge,
        100.0 * matched_with_charge as f64 / refs.len().max(1) as f64
    );
    if !unmatched_examples.is_empty() {
        eprintln!("  examples of reference peaks NOT rediscovered:");
        for e in &unmatched_examples {
            eprintln!("    - {e}");
        }
    }
}
