//! Off-by-one envelope probe — a diagnostic for the monoisotope off-by-one misses.
//!
//! For a handful of reference peaks the detector got off by one ¹³C unit, this pulls out:
//!   1. the **individual per-scan envelopes** in the apex ± 3 scan window (the exact averaging
//!      window `refine_feature` uses), and
//!   2. the **averaged composite envelope** (`average_spectra`, the same call `refine_feature` makes),
//! then runs the parity-gated `classic_deconvolute` on every one of them and reports the monoisotopic
//! mass / charge each yields, labelled by its ¹³C offset from the reference (theoretical) mass.
//!
//! The point: see WHERE the off-by-one enters. Does classic decon misplace the mono even on a clean
//! single scan? On the average? Or does it get it right (implying the miss is a seeding/scoring issue
//! upstream, not a deconvolution failure)?
//!
//! Peaks arrive vendor-centroided from `read_ms1_scans`. Masses run ~+10 ppm high vs the theoretical
//! reference (uncalibrated raw — see the `lumos-test-data-paths` memory), so matches are reported by
//! ¹³C offset with the residual ppm after removing that integer offset.
//!
//! Usage:
//!   cargo run --release --example obo_envelope_probe -- [spectra_file]
//! (spectra_file defaults to the CA/Lumos raw.)

use std::fs::File;
use std::io::{BufWriter, Write};

use flashlfq_core::deconvolution::{
    averagine_intensities_from_mono, classic_deconvolute, ClassicDeconvolutionParameters,
    DeconEnvelope, Polarity,
};
use flashlfq_core::isotopic_envelope::{C13_MINUS_C12, PROTON_MASS};
use flashlfq_core::peak_indexing::{read_ms1_scans, PeakIndexingEngine};
use flashlfq_core::spectral_averaging::{average_spectra, SpectralAveragingParameters};
use flashlfq_core::tolerance::PpmTolerance;

/// One reference off-by-one case to dig into.
struct Target {
    label: &'static str,
    ref_mono: f64, // theoretical neutral monoisotopic mass from the reference
    rt: f64,       // reference peak RT apex (min)
    z: i32,        // reference peak charge
    pk_mz: f64,    // reference "Peak MZ" — the observed most-abundant isotope m/z
    reported_k: i32, // the ¹³C offset at which the detector placed a feature (from offbyone_diag)
}

/// Neutral mass -> m/z at charge z (positive mode, proton adduct).
fn to_mz(mass: f64, z: i32) -> f64 {
    mass / z as f64 + PROTON_MASS
}
/// m/z -> neutral mass at charge z.
fn to_mass(mz: f64, z: i32) -> f64 {
    mz * z as f64 - z as f64 * PROTON_MASS
}

/// Isotope index of an m/z relative to the reference monoisotope, at charge z.
fn iso_index(mz: f64, ref_mono: f64, z: i32) -> f64 {
    (to_mass(mz, z) - ref_mono) / C13_MINUS_C12
}

/// Summarize one deconvoluted envelope against the reference mass: the integer ¹³C offset of its
/// mono from `ref_mono`, and the residual ppm once that offset is removed.
fn describe_env(e: &DeconEnvelope, ref_mono: f64) -> String {
    let k = ((e.monoisotopic_mass - ref_mono) / C13_MINUS_C12).round();
    let corrected_ref = ref_mono + k * C13_MINUS_C12;
    let ppm = (e.monoisotopic_mass - corrected_ref) / corrected_ref * 1e6;
    format!(
        "z{} mono {:.4} [k={:+.0}, {:+.1} ppm] peaks={} score={:.2e}",
        e.charge,
        e.monoisotopic_mass,
        k,
        ppm,
        e.peaks.len(),
        e.score
    )
}

/// Run classic decon on a local (mz, intensity) slice and print the envelopes, flagging the one at
/// the target charge whose mono m/z is closest to the reference mono m/z (what `refine_feature` picks).
fn decon_and_report(
    mz: &[f64],
    inten: &[f64],
    t: &Target,
    decon: &ClassicDeconvolutionParameters,
    range_min: f64,
    range_max: f64,
    indent: &str,
) {
    if mz.len() < 2 {
        println!("{indent}(only {} peak(s) — nothing to deconvolute)", mz.len());
        return;
    }
    let envs = classic_deconvolute(mz, inten, range_min, range_max, decon);
    if envs.is_empty() {
        println!("{indent}decon: (no envelopes)");
        return;
    }
    let mono_mz = to_mz(t.ref_mono, t.z);
    // The refine_feature pick: matching charge, mono m/z closest to the (true) feature m/z.
    let pick = envs
        .iter()
        .filter(|e| e.charge == t.z)
        .min_by(|a, b| {
            let da = (to_mz(a.monoisotopic_mass, a.charge) - mono_mz).abs();
            let db = (to_mz(b.monoisotopic_mass, b.charge) - mono_mz).abs();
            da.total_cmp(&db)
        });
    for (i, e) in envs.iter().take(6).enumerate() {
        let is_pick = pick.map_or(false, |p| std::ptr::eq(p, e));
        let flag = if is_pick { " <-- refine pick" } else { "" };
        println!("{indent}[{i}] {}{flag}", describe_env(e, t.ref_mono));
    }
    if let Some(p) = pick {
        let k = ((p.monoisotopic_mass - t.ref_mono) / C13_MINUS_C12).round() as i32;
        let verdict = if k == 0 { "CORRECT" } else { "OFF-BY-ONE" };
        println!("{indent}=> pick at z{} is {verdict} (k={k:+})", t.z);
    } else {
        println!("{indent}=> no envelope at target charge z{}", t.z);
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let spectra_path = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| r"D:\SP_Tutorial\Lumos\04-17-23_CA_Tryp_HCD_10min.raw".to_string());

    // Off-by-one cases pulled from offbyone_diag / obo_list. Clean single-charge-carrier peptides
    // spanning both k directions, plus one higher-mass multi-charge case.
    let targets = [
        Target { label: "YENEVALR",       ref_mono: 992.4927,  rt: 11.71, z: 2, pk_mz: 497.2584,  reported_k: 1 },
        Target { label: "VLDELTLTK",      ref_mono: 1030.5910, rt: 13.84, z: 2, pk_mz: 516.3077,  reported_k: -1 },
        Target { label: "TAGWNIPMGLLYSK", ref_mono: 1549.7963, rt: 16.50, z: 2, pk_mz: 775.9128,  reported_k: -1 },
        Target { label: "AVVQDPALKP(z3)", ref_mono: 2198.1947, rt: 15.91, z: 3, pk_mz: 734.0803,  reported_k: -1 },
    ];

    eprintln!("reading MS1 scans from {spectra_path} ...");
    let scans = read_ms1_scans(&spectra_path).expect("failed to read spectra file");
    let engine = PeakIndexingEngine::index_peaks(&scans).expect("no indexable MS1 peaks");
    eprintln!("  {} MS1 scans indexed\n", scans.len());

    // Same params as the pipeline: decon 10 ppm / ratio 3, charge 1..6; default averaging.
    let decon = ClassicDeconvolutionParameters::new(1, 6, 10.0, 3.0, Polarity::Positive);
    let avg = SpectralAveragingParameters::default();
    let ppm20 = PpmTolerance::new(20.0);

    // Optional data export for plotting: EXPORT_DIR=<dir> writes obo_series.tsv (long-form peaks)
    // and obo_meta.tsv (per-case scalars). Python (plot_obo_envelopes.py) turns them into figures.
    let export_dir = std::env::var("EXPORT_DIR").ok();
    let (mut series_w, mut meta_w, mut scans_w) = match &export_dir {
        Some(d) => {
            let mut s = BufWriter::new(File::create(format!("{d}/obo_series.tsv")).unwrap());
            let mut m = BufWriter::new(File::create(format!("{d}/obo_meta.tsv")).unwrap());
            let mut sc = BufWriter::new(File::create(format!("{d}/obo_scans.tsv")).unwrap());
            writeln!(s, "case\tseries\tmz\tintensity").unwrap();
            writeln!(m, "case\tref_mono\tz\trt\tpk_mz\treported_k\tspacing\tmono_mz").unwrap();
            writeln!(sc, "case\tscan_offset\tone_based_scan\trt\tmz\tintensity").unwrap();
            (Some(s), Some(m), Some(sc))
        }
        None => (None, None, None),
    };

    for t in &targets {
        let mono_mz = to_mz(t.ref_mono, t.z);
        let spacing = C13_MINUS_C12 / t.z as f64;
        println!("========================================================================");
        println!(
            "{}  ref_mono={:.4}  z{}  RT {:.2}  pkMZ {:.4}  mono m/z {:.4}  (detector reported k={:+})",
            t.label, t.ref_mono, t.z, t.rt, t.pk_mz, mono_mz, t.reported_k
        );
        println!("  (isotope spacing {:.4} Th; peaks labelled by k = isotope index vs ref mono)", spacing);

        // --- find the apex scan: max intensity at the observed most-abundant m/z near the ref RT ---
        let info = engine.scan_info();
        let mut apex_idx: i32 = -1;
        let mut apex_int = f64::NEG_INFINITY;
        for si in info {
            if (si.retention_time - t.rt).abs() > 0.4 {
                continue;
            }
            if let Some(p) = engine.get_indexed_peak(t.pk_mz, si.zero_based_scan_index, &ppm20) {
                if p.intensity as f64 > apex_int {
                    apex_int = p.intensity as f64;
                    apex_idx = si.zero_based_scan_index;
                }
            }
        }
        if apex_idx < 0 {
            // Fall back to the scan nearest the ref RT.
            apex_idx = info
                .iter()
                .min_by(|a, b| {
                    (a.retention_time - t.rt)
                        .abs()
                        .total_cmp(&(b.retention_time - t.rt).abs())
                })
                .map(|s| s.zero_based_scan_index)
                .unwrap_or(0);
            println!("  (no peak at pkMZ near RT; using nearest-RT scan {apex_idx})");
        }

        // --- window = apex ± 3 (<= 7 scans), exactly as refine_feature ---
        let half = 3i32;
        let lo = (apex_idx - half).max(0) as usize;
        let hi = ((apex_idx + half).max(0) as usize).min(scans.len() - 1);
        println!(
            "  apex scan {} (RT {:.3}, pkMZ int {:.2e}); window scans {}..={}",
            apex_idx, info[apex_idx as usize].retention_time, apex_int, lo, hi
        );

        // Deconvolution range: one isotope below the mono (to catch a k=-1 miss) up through +6
        // isotopes. Kept tight to the peptide envelope so classic decon isn't captured by an intense
        // co-eluting interferent several isotopes away. The averaging slice pads this slightly.
        let n_iso = 6.0;
        let range_min = mono_mz - 1.1 * spacing;
        let range_max = mono_mz + n_iso * spacing;
        let slice_lo = range_min - 0.3;
        let slice_hi = range_max + 0.3;

        let window = &scans[lo..=hi];
        let mut x_arrays: Vec<Vec<f64>> = Vec::with_capacity(window.len());
        let mut y_arrays: Vec<Vec<f64>> = Vec::with_capacity(window.len());

        // --- individual per-scan envelopes ---
        println!("\n  --- individual per-scan envelopes (peaks: k:m/z:intensity) ---");
        for s in window {
            let a = s.mz.partition_point(|&m| m < slice_lo);
            let b = s.mz.partition_point(|&m| m <= slice_hi);
            let smz = &s.mz[a..b];
            let sint = &s.intensity[a..b];
            x_arrays.push(smz.to_vec());
            y_arrays.push(sint.to_vec());

            let peaks: String = smz
                .iter()
                .zip(sint.iter())
                .map(|(&m, &i)| format!("{:+.2}:{:.3}:{:.1e}", iso_index(m, t.ref_mono, t.z), m, i))
                .collect::<Vec<_>>()
                .join("  ");
            println!(
                "  scan {:>4} RT {:.3}: {}",
                s.one_based_scan_number, s.retention_time, peaks
            );
            decon_and_report(smz, sint, t, &decon, range_min, range_max, "      ");
        }

        // --- averaged composite ---
        println!("\n  --- averaged composite envelope (average_spectra, {} scans) ---", window.len());
        let (cmz, cint) = average_spectra(&x_arrays, &y_arrays, &avg);
        let cpeaks: String = cmz
            .iter()
            .zip(cint.iter())
            .map(|(&m, &i)| format!("{:+.2}:{:.3}:{:.1e}", iso_index(m, t.ref_mono, t.z), m, i))
            .collect::<Vec<_>>()
            .join("  ");
        println!("  composite: {cpeaks}");
        decon_and_report(&cmz, &cint, t, &decon, range_min, range_max, "      ");
        println!();

        // --- data export for plotting ---
        if let (Some(sw), Some(mw)) = (series_w.as_mut(), meta_w.as_mut()) {
            writeln!(
                mw,
                "{}\t{:.4}\t{}\t{:.3}\t{:.4}\t{}\t{:.5}\t{:.5}",
                t.label, t.ref_mono, t.z, t.rt, t.pk_mz, t.reported_k, spacing, mono_mz
            )
            .unwrap();

            // Observed composite peaks.
            for (&m, &i) in cmz.iter().zip(cint.iter()) {
                writeln!(sw, "{}\tobserved_composite\t{:.5}\t{:.6e}", t.label, m, i).unwrap();
            }
            // Observed apex-scan peaks (single strongest scan, sliced to the same window).
            let apex = &scans[apex_idx as usize];
            let a = apex.mz.partition_point(|&m| m < slice_lo);
            let b = apex.mz.partition_point(|&m| m <= slice_hi);
            for (&m, &i) in apex.mz[a..b].iter().zip(apex.intensity[a..b].iter()) {
                writeln!(sw, "{}\tobserved_apex\t{:.5}\t{:.6e}", t.label, m, i).unwrap();
            }
            // Every window scan, tagged by offset from the apex, for elution exploration.
            if let Some(scw) = scans_w.as_mut() {
                for idx in lo..=hi {
                    let s = &scans[idx];
                    let a = s.mz.partition_point(|&m| m < slice_lo);
                    let b = s.mz.partition_point(|&m| m <= slice_hi);
                    for (&m, &i) in s.mz[a..b].iter().zip(s.intensity[a..b].iter()) {
                        writeln!(
                            scw,
                            "{}\t{}\t{}\t{:.4}\t{:.5}\t{:.6e}",
                            t.label,
                            idx as i32 - apex_idx,
                            s.one_based_scan_number,
                            s.retention_time,
                            m,
                            i
                        )
                        .unwrap();
                    }
                }
            }
            // Averagine theoretical templates anchored at three candidate monoisotopes: the
            // reference mono (k=0) and its ±1 ¹³C neighbours. Weights are max-normalised to 1.0; the
            // Python side scales them to the observed envelope. Placed at m/z of each isotope tooth.
            for (tag, k) in [("avg_mono", 0i32), ("avg_mono_m1", -1), ("avg_mono_p1", 1)] {
                let cand = t.ref_mono + k as f64 * C13_MINUS_C12;
                let weights = averagine_intensities_from_mono(cand, 1e-3, 14);
                for (j, &w) in weights.iter().enumerate() {
                    let mz = to_mz(cand + j as f64 * C13_MINUS_C12, t.z);
                    writeln!(sw, "{}\t{}\t{:.5}\t{:.6e}", t.label, tag, mz, w).unwrap();
                }
            }
        }
    }

    if let Some(mut sw) = series_w.take() {
        sw.flush().unwrap();
    }
    if let Some(mut mw) = meta_w.take() {
        mw.flush().unwrap();
    }
    if let Some(mut scw) = scans_w.take() {
        scw.flush().unwrap();
    }
    if let Some(d) = &export_dir {
        eprintln!("exported plot data to {d}/obo_series.tsv, obo_meta.tsv, obo_scans.tsv");
    }
}
