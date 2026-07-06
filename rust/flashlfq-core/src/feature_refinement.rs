//! Detect-then-refine and charge-state consensus — the assembly layer above the trace-kernel
//! detector.
//!
//! This is **new algorithm work** (no mzLib counterpart, so no C# golden). It wires the two
//! parity-gated primitives — [`crate::spectral_averaging`] and [`crate::deconvolution`] — onto the
//! [`crate::trace_kernel`] detector's output, following the design in
//! `agent_info/Feature-Detection-Design.md` ("Window + averaging", "Final deconvolution",
//! "Feature assembly — charge-state consensus on candidate masses").
//!
//! ## Two stages
//! 1. **[`refine_feature`] (detect-then-refine).** A [`crate::trace_kernel::DetectedFeature`] gives
//!    a coarse mass/charge and, crucially, the *scan extent* of its elution. We average that extent's
//!    MS1 scans into a single high-SNR composite spectrum (SNR ~√N; the m/z binning gives it its mass
//!    accuracy) and run the parity-gated `classic_deconvolute` on it. The composite envelope matching
//!    the detected charge yields a **precise** monoisotopic mass plus the raw per-peak candidate-mass
//!    list (see [`crate::deconvolution::DeconEnvelope::monoisotopic_mass_predictions`]).
//! 2. **[`resolve_charge_state_consensus`].** Co-eluting charge states of the same peptide are
//!    grouped, then the true monoisotopic mass is *resolved* by pooling every member's candidate
//!    masses and picking the cluster with the greatest **cross-charge support**. The true mass
//!    recurs across charges; a spurious ±1 Da monoisotopic off-by-one does not — so this resolves the
//!    off-by-one rather than merely averaging noise down.
//!
//! ## Design choices (flagged; this is de-novo, not a port)
//! - **Refinement window = the traced elution extent.** The averaging window is exactly the
//!   scan-index span of the feature's claimed peaks (`min..=max`), which the detector already limited
//!   to ~FWHM via its RT-Gaussian window — so co-eluting neighbours are not pulled into the composite.
//! - **Off-by-one-aware grouping.** Two co-eluting features are grouped when their neutral masses
//!   agree within `mass_tolerance_ppm` *after allowing for an integer number of ¹³C units* of
//!   difference (0, ±1, ±2). The mono off-by-one is precisely a ±1 ¹³C shift, so a feature whose
//!   *refined* mass is off by one must still be able to join the group whose consensus will correct
//!   it. The candidate-mass clustering *inside* the group stays at the tight `mass_tolerance_ppm`
//!   (with a 0.01 Da absolute floor), which is what keeps the true-mass and +1 candidates in separate
//!   clusters so cross-charge support can discriminate them.
//! - **Cross-charge support** = number of *distinct charges* contributing a candidate to a cluster.
//!   Ties break by candidate count, then by summed contributing intensity. The strict "present in all
//!   charges" is softened to "the max-support cluster" (design's "≥2 charges, weighted by support").

use crate::deconvolution::{classic_deconvolute, ClassicDeconvolutionParameters};
use crate::isotope_shift_decon::{
    best_charge_by_fit, envelope_fit_cosine_masked, shift_decon, shift_decon_gated,
    shift_decon_in_window, walkback_mono_high_charge, NEIGHBOR_MASK_PPM, RECHARGE_PREFER_MARGIN,

};
use crate::isotopic_envelope::{mass_to_mz_f64, C13_MINUS_C12};
use crate::peak_indexing::{PeakKey, Scan};
use crate::spectral_averaging::{average_spectra, SpectralAveragingParameters};
use crate::trace_kernel::DetectedFeature;
use std::collections::{HashMap, HashSet};

/// Absolute mass-clustering floor (Da) for small masses, where a ppm window would be tighter than
/// real mass precision. Applied as `max(mass · ppm/1e6, MASS_CLUSTER_ABS_FLOOR_DA)`.
pub const MASS_CLUSTER_ABS_FLOOR_DA: f64 = 0.01;

/// Hard cap on how many MS1 scans [`refine_feature`] averages into the composite. The composite is a
/// high-SNR snapshot at the feature apex, so this stays small — just the apex plus one scan on either
/// side; averaging more pulls in co-eluting interference. A `> MAX_SCANS_TO_AVERAGE` window is treated
/// as a bug (asserted), not silently accepted.
pub const MAX_SCANS_TO_AVERAGE: usize = 3;

/// Largest integer ¹³C off-by-one offset tolerated when grouping features by neutral mass. The mono
/// off-by-one is normally ±1; ±2 is allowed for robustness against a doubly-mis-assigned monoisotope.
const MAX_OFFBYONE_UNITS: i32 = 2;

/// RT padding (minutes) added to each side of a feature's elution window when testing tight co-elution
/// for an off-by-one link ([`features_coelute_tightly`]). ~0.6 s — about half an MS1 cycle on this
/// data — so a feature detected in only one or two scans can still match a genuinely co-eluting partner
/// whose apex sits a scan away, without admitting a species several scans distant.
const COELUTION_PAD_MIN: f64 = 0.01;

/// A detected feature refined against an averaged composite spectrum via the parity-gated classic
/// deconvolution. Carries the raw candidate-mass list the charge-state consensus intersects.
#[derive(Debug, Clone)]
pub struct RefinedFeature {
    /// The originating trace-kernel detection (cloned; RT bounds, apex, claimed peaks, coarse mass).
    pub detected: DetectedFeature,
    /// Precise monoisotopic mass from the composite-spectrum deconvolution (median of
    /// [`candidate_masses`](Self::candidate_masses), as `classic_deconvolute` sets it).
    pub refined_monoisotopic_mass: f64,
    /// Charge state of the matched composite envelope (equals `detected.charge`).
    pub refined_charge: i32,
    /// Raw per-peak monoisotopic-mass predictions from the composite envelope — the values the
    /// median in [`refined_monoisotopic_mass`](Self::refined_monoisotopic_mass) collapses. This is
    /// what the cross-charge consensus intersects.
    pub candidate_masses: Vec<f64>,
    /// The composite envelope's classic-deconvolution score.
    pub decon_score: f64,
}

/// A peptide feature resolved across co-eluting charge states. The monoisotopic mass is the
/// cross-charge consensus (off-by-one corrected), not any single charge's refined value.
#[derive(Debug, Clone)]
pub struct ResolvedFeature {
    /// Consensus monoisotopic neutral mass (the intensity/count-weighted mean of the winning
    /// candidate cluster).
    pub monoisotopic_mass: f64,
    /// The distinct charge states observed for this peptide, ascending.
    pub charge_states: Vec<i32>,
    /// Apex retention time (minutes) — taken from the most intense member.
    pub apex_rt: f64,
    /// Earliest start RT across members (minutes).
    pub start_rt: f64,
    /// Latest end RT across members (minutes).
    pub end_rt: f64,
    /// Summed intensity across all member features.
    pub summed_intensity: f64,
    /// Number of distinct charges supporting the resolved mass (the winning cluster's cross-charge
    /// support; 1 for a single-charge/fallback group).
    pub cross_charge_support: usize,
    /// The member refined features that make up this resolved feature.
    pub members: Vec<RefinedFeature>,
}

/// Refines a single detected feature against an averaged composite of its elution window.
///
/// The averaging window is the scan-index range spanned by the feature's claimed peaks
/// (`min..=max` of `feature.peaks[*].zero_based_scan_index`), clamped to `0..scans.len()`. `scans`
/// must be the same MS1 scan array (in the same order) that produced the index the feature was
/// detected on, so that a zero-based scan index addresses `scans[index]`.
///
/// Returns `None` when the window is empty, the composite is empty, or no composite envelope matches
/// the detected charge state.
/// Refines a feature using the **detector-anchored shift decon** instead of classic deconvolution:
/// anchors on the detector's own most-abundant claimed peak, then runs the FlashLFQ-style −1/0/+1
/// placement. `use_apex` selects the single apex scan (more accurate vs FlashLFQ truth in testing)
/// over the averaged composite. Returns `None` on the same degenerate conditions as [`refine_feature`],
/// or when the shift decon finds no envelope.
///
/// The resulting [`RefinedFeature`] carries a single candidate mass (the shift mono); the cross-charge
/// consensus still clusters those across charges. Experiment path (pipeline `REFINE_METHOD`).
///
/// When `recharge` is set, the charge is **re-selected** by the envelope-fit cosine
/// ([`best_charge_by_fit`]) over the harmonic candidates {z/2, z, 2z}: the charge whose placed
/// envelope best fits *and* explains the window wins. This recovers charge-halved features (a real
/// z=2 the detector labelled z=1) and rejects the doubled-mass harmonic, by maximising the unified
/// fit/explained/completeness metric rather than trusting the detector's charge.
pub fn refine_feature_shift(
    feature: &DetectedFeature,
    scans: &[Scan],
    averaging_params: &SpectralAveragingParameters,
    shift_tol_ppm: f64,
    use_apex: bool,
    recharge: bool,
) -> Option<RefinedFeature> {
    refine_feature_shift_inner(
        feature,
        scans,
        averaging_params,
        shift_tol_ppm,
        use_apex,
        recharge,
        &[],
    )
}

/// **Neighbour-aware** [`refine_feature_shift`]: `neighbor_mz` is the ascending isotope m/z of
/// co-eluting *other* features that are **not** on this feature's own grid (from [`NeighborIndex`]).
/// Those peaks are masked from the charge-selection fit, the walk-back and the reported score, so a
/// low-scoring feature in a crowded window is judged on the signal plausibly its own — a competing
/// charge cannot borrow a neighbour's peaks (the z↔2z harmonic), and the walk-back cannot anchor on a
/// neighbour's peak. Experiment path (pipeline `NEIGHBOR_REFINE`). With `neighbor_mz` empty this is
/// exactly [`refine_feature_shift`].
pub fn refine_feature_shift_neighbor(
    feature: &DetectedFeature,
    scans: &[Scan],
    averaging_params: &SpectralAveragingParameters,
    shift_tol_ppm: f64,
    use_apex: bool,
    recharge: bool,
    neighbor_mz: &[f64],
) -> Option<RefinedFeature> {
    refine_feature_shift_inner(
        feature,
        scans,
        averaging_params,
        shift_tol_ppm,
        use_apex,
        recharge,
        neighbor_mz,
    )
}

#[allow(clippy::too_many_arguments)]
fn refine_feature_shift_inner(
    feature: &DetectedFeature,
    scans: &[Scan],
    averaging_params: &SpectralAveragingParameters,
    shift_tol_ppm: f64,
    use_apex: bool,
    recharge: bool,
    neighbor_mz: &[f64],
) -> Option<RefinedFeature> {
    let s = build_feature_slices(feature, scans, averaging_params)?;
    // Detector's own most-abundant claimed peak — the anchor that cannot grab a foreign peak.
    let anchor_mz = feature
        .peaks
        .iter()
        .max_by(|a, b| a.intensity.total_cmp(&b.intensity))
        .map(|p| p.m() as f64)?;
    let (mz, inten) = if use_apex {
        (&s.apex_mz, &s.apex_int)
    } else {
        (&s.comp_mz, &s.comp_int)
    };

    let (refined_charge, refined_mono0) = if recharge {
        let candidates = charge_candidates(feature.charge);
        // Keep the detector's charge unless another candidate fits clearly better — blocks a spurious
        // z<->2z harmonic flip in crowded windows (e.g. a real z=2 re-charged to z=4, mass doubled).
        let (z, mono, _cos) = best_charge_by_fit(
            mz,
            inten,
            anchor_mz,
            &candidates,
            shift_tol_ppm,
            0.0,
            feature.charge,
            RECHARGE_PREFER_MARGIN,
            neighbor_mz,
            NEIGHBOR_MASK_PPM,
        )?;
        (z, mono)
    } else {
        let r = shift_decon(mz, inten, anchor_mz, feature.charge, shift_tol_ppm)?;
        (feature.charge, r.monoisotopic_mass)
    };
    // High-charge double-check: heavy peptides can seed the mono one or two ¹³C too high (the envelope
    // mode sits well above the mono), and the fit window never sees the unexplained peak beneath it.
    // Walk the mono back up to 2 ¹³C and keep the lowest that still fits (no-op for |z| < 4 / good mono).
    let refined_mono = walkback_mono_high_charge(
        mz,
        inten,
        refined_mono0,
        refined_charge,
        shift_tol_ppm,
        0.0,
        neighbor_mz,
        NEIGHBOR_MASK_PPM,
    );
    Some(RefinedFeature {
        detected: feature.clone(),
        refined_monoisotopic_mass: refined_mono,
        refined_charge,
        candidate_masses: vec![refined_mono],
        decon_score: envelope_fit_cosine_masked(
            mz,
            inten,
            mass_to_mz_f64(refined_mono, refined_charge),
            refined_charge,
            shift_tol_ppm,
            0.2,
            0.0,
            neighbor_mz,
            NEIGHBOR_MASK_PPM,
        ),
    })
}

/// Harmonic charge candidates for re-selection: the detector's charge plus its half and double,
/// bounded to `[1, 6]` and deduped. Catches both charge-halving (real z=2 labelled z=1 → include 2z)
/// and the doubled-mass harmonic (real z=1 labelled z=2 → include z/2).
fn charge_candidates(z: i32) -> Vec<i32> {
    let mut c = vec![z];
    if z * 2 <= 6 {
        c.push(z * 2);
    }
    if z % 2 == 0 && z / 2 >= 1 {
        c.push(z / 2);
    }
    c.retain(|&x| (1..=6).contains(&x));
    c.sort_unstable();
    c.dedup();
    c
}

/// Whether `mz` falls on this feature's isotope grid — within [`GRID_PPM`] of `mono_mz + k·spacing`
/// for some isotope index `k ∈ [-1, num_isotopes_observed + 2]`. Grid-gated censoring keeps on-grid
/// peaks (this feature's own envelope teeth, whoever claimed them) and only removes off-grid claimed
/// interferents.
fn on_isotope_grid(mz: f64, feature: &DetectedFeature, spacing: f64) -> bool {
    /// ppm window for treating an observed peak as one of the feature's expected isotope teeth.
    const GRID_PPM: f64 = 15.0;
    if spacing <= 0.0 {
        return false;
    }
    let k = ((mz - feature.mono_mz) / spacing).round();
    let kmax = feature.num_isotopes_observed as f64 + 2.0;
    if k < -1.0 || k > kmax {
        return false;
    }
    let expected = feature.mono_mz + k * spacing;
    (mz - expected).abs() / expected * 1e6 <= GRID_PPM
}

pub fn refine_feature(
    feature: &DetectedFeature,
    scans: &[Scan],
    averaging_params: &SpectralAveragingParameters,
    decon_params: &ClassicDeconvolutionParameters,
) -> Option<RefinedFeature> {
    refine_feature_inner(feature, scans, averaging_params, decon_params, None)
}

/// [`refine_feature`] with **subtractive peak censoring**: peaks claimed by *other* features are
/// removed from this feature's composite before deconvolution.
///
/// `all_claimed` is the union of every detected feature's claimed-peak keys (from
/// [`DetectedFeature::peaks`]). Because the trace-kernel detector claims greedily tallest-first,
/// every peak is owned by exactly one feature; when refining feature X we drop any window peak that
/// is claimed but not owned by X — i.e. the peaks that stronger, already-assigned features took.
/// This gives the deconvolution a cleaner window (co-eluting interferents removed) so it cannot
/// anchor an envelope on a neighbouring species' peak. This feature's *own* claimed peaks and any
/// **unclaimed** peaks (noise / not-yet-explained signal) are kept.
///
/// This is an experiment path (see the pipeline's `CENSOR_CLAIMED` flag); the uncensored
/// [`refine_feature`] remains the default.
pub fn refine_feature_censored(
    feature: &DetectedFeature,
    scans: &[Scan],
    averaging_params: &SpectralAveragingParameters,
    decon_params: &ClassicDeconvolutionParameters,
    all_claimed: &HashSet<PeakKey>,
) -> Option<RefinedFeature> {
    refine_feature_inner(feature, scans, averaging_params, decon_params, Some(all_claimed))
}

fn refine_feature_inner(
    feature: &DetectedFeature,
    scans: &[Scan],
    averaging_params: &SpectralAveragingParameters,
    decon_params: &ClassicDeconvolutionParameters,
    censor: Option<&HashSet<PeakKey>>,
) -> Option<RefinedFeature> {
    if scans.is_empty() || feature.peaks.is_empty() {
        return None;
    }

    // Averaging window = a SMALL number of scans centred on the feature's apex. The composite is a
    // high-SNR snapshot of the envelope at its strongest point, NOT the whole elution — a wide
    // window pulls in co-eluting interference and defeats the purpose (SpectralAveraging's own
    // default is 5 scans). We deliberately do NOT use the feature's full claimed-peak scan extent:
    // the detector's RT window is ~±2σ, which in dense MS1 regions is ~100 scans.
    let apex = feature.apex_scan_index;
    let half = (MAX_SCANS_TO_AVERAGE / 2) as i32; // 1 → up to 3 scans (apex ± 1)
    let lo = (apex - half).max(0) as usize;
    let hi = ((apex + half).max(0) as usize).min(scans.len() - 1);
    if lo > hi {
        return None;
    }
    // Safety invariant: averaging more than a handful of scans means the window logic is wrong.
    let n_avg = hi - lo + 1;
    assert!(
        n_avg <= MAX_SCANS_TO_AVERAGE,
        "refine_feature would average {n_avg} scans (> {MAX_SCANS_TO_AVERAGE}); the averaging \
         window must stay small — something is wrong with the window computation"
    );

    // m/z window around the feature. This is both the deconvolution range AND the slice we average
    // over: averaging only the local neighbourhood — instead of binning every peak of every window
    // scan — is the key cost reduction (a Lumos MS1 scan has thousands of peaks; binning all of them
    // for all ~10^4 features is what made refinement take minutes). The composite only needs the
    // envelope's vicinity. Trade-off: `RelativeToTics` normalization now uses the sliced TIC rather
    // than the full-scan TIC, changing inter-scan weighting slightly — negligible for the charge/mass
    // inference the deconvolution performs.
    let spacing = C13_MINUS_C12 / feature.charge as f64;
    let range_min0 = (feature.mono_mz - 1.5).max(0.0);
    let range_max =
        feature.mono_mz + (feature.num_isotopes_observed as f64 + 3.0) * spacing + 1.0;
    let slice_lo = range_min0 - 0.5;
    let slice_hi = range_max + 0.5;

    // Gather the (mz, intensity) arrays of the window's scans, sliced to the local m/z window. Each
    // scan's m/z is ascending, so binary-search the bounds and copy just that sub-range.
    let window = &scans[lo..=hi];
    let mut x_arrays: Vec<Vec<f64>> = Vec::with_capacity(window.len());
    let mut y_arrays: Vec<Vec<f64>> = Vec::with_capacity(window.len());
    for (wi, s) in window.iter().enumerate() {
        let a = s.mz.partition_point(|&m| m < slice_lo);
        let b = s.mz.partition_point(|&m| m <= slice_hi);
        match censor {
            None => {
                x_arrays.push(s.mz[a..b].to_vec());
                y_arrays.push(s.intensity[a..b].to_vec());
            }
            Some(claimed) => {
                // Absolute scan index == position in `scans` == the peaks' zero_based_scan_index.
                let abs = (lo + wi) as i32;
                let mut mzv = Vec::with_capacity(b - a);
                let mut inv = Vec::with_capacity(b - a);
                for j in a..b {
                    // Reproduce the indexed peak's f32-narrowed key (index_peaks built keys from the
                    // same raw scans via `mz as f32` / `intensity as f32`).
                    let key: PeakKey = (
                        abs,
                        (s.mz[j] as f32).to_bits(),
                        (s.intensity[j] as f32).to_bits(),
                    );
                    // Grid-gated censoring: drop a peak only if another feature claimed it AND it is
                    // NOT on this feature's own isotope grid. This removes off-grid interferents
                    // (the peaks a competing envelope could anchor on) while preserving every tooth
                    // of this feature's envelope — even isotopes a marginally-taller neighbour
                    // claimed first — so a noisy claim partition can no longer delete real signal.
                    if claimed.contains(&key) && !on_isotope_grid(s.mz[j], feature, spacing) {
                        continue;
                    }
                    mzv.push(s.mz[j]);
                    inv.push(s.intensity[j]);
                }
                x_arrays.push(mzv);
                y_arrays.push(inv);
            }
        }
    }
    if x_arrays.iter().all(|x| x.is_empty()) {
        return None;
    }

    let (comp_mz, comp_intensity) = average_spectra(&x_arrays, &y_arrays, averaging_params);
    if comp_mz.is_empty() {
        return None;
    }
    // Never ask the deconvolution for a range below the smallest composite peak.
    let range_min = range_min0.max(comp_mz[0]);

    let envelopes =
        classic_deconvolute(&comp_mz, &comp_intensity, range_min, range_max, decon_params);

    // Pick the envelope matching the detected charge whose mono m/z is closest to the feature's.
    let best = envelopes
        .into_iter()
        .filter(|e| e.charge == feature.charge)
        .min_by(|a, b| {
            let da = (mass_to_mz_f64(a.monoisotopic_mass, a.charge) - feature.mono_mz).abs();
            let db = (mass_to_mz_f64(b.monoisotopic_mass, b.charge) - feature.mono_mz).abs();
            da.total_cmp(&db)
        })?;

    // The deconvolution's monoisotope is trusted here. A per-feature off-by-one corrector was attempted
    // twice — a cosine envelope match and a detector-comb-anchor snap — and BOTH regressed reference
    // recall on chimeric composites (the snap: 88.4% → 80.8%), because the discriminating signal is too
    // weak on real co-eluting data. Off-by-one is deferred to dedicated discriminator work; where
    // co-eluting charges agree, `resolve_charge_state_consensus` already corrects it downstream.
    Some(RefinedFeature {
        detected: feature.clone(),
        refined_monoisotopic_mass: best.monoisotopic_mass,
        refined_charge: best.charge,
        candidate_masses: best.monoisotopic_mass_predictions.clone(),
        decon_score: best.score,
    })
}

// ---------------------------------------------------------------------------
// Four-way decon comparator
// ---------------------------------------------------------------------------
//
// Runs two deconvolution *strategies* — the parity-gated classic `classic_deconvolute` and the
// untargeted FlashLFQ-style [`crate::isotope_shift_decon`] — over two *spectrum views* — the
// averaged apex±1 composite and the single apex scan — giving up to **four** monoisotope verdicts
// for one detected feature. When the four agree (same integer ¹³C offset) the placement is
// confident; when they disagree the feature is flagged for a more advanced multi-envelope decon
// (the disagreement is the signal that a single averagine envelope does not explain the window —
// typically a chimeric region). See the `isotope-shift-decon` validation: classic and shift can
// each be individually fooled, and heavy peptides can be *unanimously* wrong, so this comparator
// deliberately triggers on **disagreement only**; unanimous-but-wrong is left to the downstream
// cross-charge consensus, which is the correct backstop for the heavy-peptide case.

/// An RT-bucketed index of detected features, used to find the co-eluting neighbours of a feature so
/// their predicted isotope m/z can be forbidden as shift-decon anchors ([`four_way_decon_gated`]).
///
/// Built once over all detections. `forbidden_positions` returns, for one feature, the ascending
/// isotope m/z of every *other* feature whose apex RT is within `rt_tol` and whose teeth fall in the
/// query window — i.e. the peaks that belong to a co-eluting neighbour, not to this feature.
pub struct NeighborIndex {
    /// One entry per feature (parallel to the input slice): `(apex_rt, mono_mz, spacing, kmax)`.
    grids: Vec<(f64, f64, f64, i32)>,
    /// Per-feature summed intensity (parallel to `grids`), for stronger-neighbour filtering.
    intensities: Vec<f64>,
    /// RT bin (`floor(apex_rt / rt_tol)`) → feature indices.
    buckets: HashMap<i64, Vec<usize>>,
    rt_tol: f64,
}

impl NeighborIndex {
    /// Builds the index over `features`, bucketing by `rt_tol`-wide RT bins. Each feature's isotope
    /// grid spans `k ∈ [0, num_isotopes_observed + 2]` (a little past the observed envelope).
    pub fn build(features: &[DetectedFeature], rt_tol: f64) -> Self {
        let mut grids = Vec::with_capacity(features.len());
        let mut intensities = Vec::with_capacity(features.len());
        let mut buckets: HashMap<i64, Vec<usize>> = HashMap::new();
        let bin_of = |rt: f64| -> i64 { (rt / rt_tol).floor() as i64 };
        for (i, f) in features.iter().enumerate() {
            let spacing = C13_MINUS_C12 / f.charge.max(1) as f64;
            grids.push((
                f.apex_rt,
                f.mono_mz,
                spacing,
                f.num_isotopes_observed as i32 + 2,
            ));
            intensities.push(f.summed_intensity);
            buckets.entry(bin_of(f.apex_rt)).or_default().push(i);
        }
        NeighborIndex { grids, intensities, buckets, rt_tol }
    }

    /// Builds the index from **refined** features (their corrected mono m/z / charge), for a second
    /// neighbour-aware refinement pass whose masks come from placements already corrected in pass one.
    pub fn build_from_refined(refined: &[RefinedFeature], rt_tol: f64) -> Self {
        let mut grids = Vec::with_capacity(refined.len());
        let mut intensities = Vec::with_capacity(refined.len());
        let mut buckets: HashMap<i64, Vec<usize>> = HashMap::new();
        let bin_of = |rt: f64| -> i64 { (rt / rt_tol).floor() as i64 };
        for (i, r) in refined.iter().enumerate() {
            let spacing = C13_MINUS_C12 / r.refined_charge.max(1) as f64;
            let mono_mz = mass_to_mz_f64(r.refined_monoisotopic_mass, r.refined_charge);
            grids.push((
                r.detected.apex_rt,
                mono_mz,
                spacing,
                r.detected.num_isotopes_observed as i32 + 2,
            ));
            intensities.push(r.detected.summed_intensity);
            buckets.entry(bin_of(r.detected.apex_rt)).or_default().push(i);
        }
        NeighborIndex { grids, intensities, buckets, rt_tol }
    }

    /// Ascending isotope m/z of every feature other than `self_idx` that co-elutes with it
    /// (`|Δapex_rt| ≤ rt_tol`) and whose teeth fall in `[win_min_mz, win_max_mz]`.
    pub fn forbidden_positions(
        &self,
        self_idx: usize,
        win_min_mz: f64,
        win_max_mz: f64,
    ) -> Vec<f64> {
        self.collect_positions(self_idx, win_min_mz, win_max_mz, 0.0)
    }

    /// Ascending isotope m/z of indexed features that co-elute with `query_rt` (`|Δrt| ≤ rt_tol`), are
    /// at least `min_intensity` intense, and whose teeth fall in `[win_min_mz, win_max_mz]`. Unlike
    /// [`forbidden_positions`](Self::forbidden_positions) this takes the query's RT directly (no
    /// `self_idx`), for querying an index built from a *different* feature set (e.g. refined vs the
    /// detected loop). A feature never masks its own peaks here as long as `min_intensity` exceeds its
    /// own intensity — guaranteed when the caller uses a strength ratio > 1.
    pub fn forbidden_positions_query(
        &self,
        query_rt: f64,
        win_min_mz: f64,
        win_max_mz: f64,
        min_intensity: f64,
    ) -> Vec<f64> {
        let bin = (query_rt / self.rt_tol).floor() as i64;
        let mut out: Vec<f64> = Vec::new();
        for b in (bin - 1)..=(bin + 1) {
            let Some(idxs) = self.buckets.get(&b) else { continue };
            for &j in idxs {
                if self.intensities[j] < min_intensity {
                    continue;
                }
                let (rt, mono_mz, spacing, kmax) = self.grids[j];
                if (rt - query_rt).abs() > self.rt_tol {
                    continue;
                }
                for k in 0..=kmax {
                    let m = mono_mz + k as f64 * spacing;
                    if m < win_min_mz {
                        continue;
                    }
                    if m > win_max_mz {
                        break;
                    }
                    out.push(m);
                }
            }
        }
        out.sort_by(f64::total_cmp);
        out
    }

    /// Like [`forbidden_positions`](Self::forbidden_positions) but only from neighbours whose summed
    /// intensity is at least `min_ratio ×` this feature's — i.e. defer only to **stronger** co-eluting
    /// species (a weak feature in a strong neighbour's shadow), never mask in favour of a weaker one.
    pub fn forbidden_positions_stronger(
        &self,
        self_idx: usize,
        win_min_mz: f64,
        win_max_mz: f64,
        min_ratio: f64,
    ) -> Vec<f64> {
        let min_intensity = self.intensities[self_idx] * min_ratio;
        self.collect_positions(self_idx, win_min_mz, win_max_mz, min_intensity)
    }

    fn collect_positions(
        &self,
        self_idx: usize,
        win_min_mz: f64,
        win_max_mz: f64,
        min_intensity: f64,
    ) -> Vec<f64> {
        let (self_rt, ..) = self.grids[self_idx];
        let bin = (self_rt / self.rt_tol).floor() as i64;
        let mut out: Vec<f64> = Vec::new();
        for b in (bin - 1)..=(bin + 1) {
            let Some(idxs) = self.buckets.get(&b) else { continue };
            for &j in idxs {
                if j == self_idx || self.intensities[j] < min_intensity {
                    continue;
                }
                let (rt, mono_mz, spacing, kmax) = self.grids[j];
                if (rt - self_rt).abs() > self.rt_tol {
                    continue;
                }
                for k in 0..=kmax {
                    let m = mono_mz + k as f64 * spacing;
                    if m < win_min_mz {
                        continue;
                    }
                    if m > win_max_mz {
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

/// Which of the four strategy×view combinations produced a verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeconView {
    /// Classic deconvolution on the averaged apex±1 composite (what [`refine_feature`] uses).
    ClassicComposite,
    /// Classic deconvolution on the single apex scan.
    ClassicApex,
    /// FlashLFQ-style shift decon on the averaged composite.
    ShiftComposite,
    /// FlashLFQ-style shift decon on the single apex scan.
    ShiftApex,
}

/// One monoisotope verdict from a single [`DeconView`].
#[derive(Debug, Clone, PartialEq)]
pub struct DeconVerdict {
    /// Which strategy×view produced this.
    pub view: DeconView,
    /// The monoisotopic neutral mass this view placed the feature at.
    pub monoisotopic_mass: f64,
    /// Charge state (equals the detected charge).
    pub charge: i32,
    /// Integer ¹³C offset of this verdict's mono from the detector's coarse mono
    /// (`round((mono − detector_mono) / C13)`); the quantity the four views are compared on.
    pub offset_k: i32,
    /// Confidence for this view: classic views are confident when an envelope at the detected
    /// charge was picked; shift views are confident when the shift-0 acceptance gate passed. A
    /// gate-failing / unpicked view still contributes an `offset_k` but a low-confidence one.
    pub confident: bool,
    /// Classic decon score, or (for shift views) the winning shift's Pearson correlation.
    pub score: f64,
}

/// The four-way decon comparison for one feature.
#[derive(Debug, Clone, PartialEq)]
pub struct FourWayDecon {
    /// The detector's coarse monoisotopic mass (the reference the offsets are measured against).
    pub detector_mono: f64,
    /// The feature's charge state.
    pub charge: i32,
    /// The verdicts that were produced (0–4; a view is omitted when its decon yields nothing).
    pub verdicts: Vec<DeconVerdict>,
    /// Whether every produced verdict shares the same [`offset_k`](DeconVerdict::offset_k).
    pub unanimous: bool,
    /// The agreed offset when unanimous, else the plurality offset (ties broken toward the most
    /// confident, then the smallest offset). `None` only when no verdict was produced.
    pub consensus_k: Option<i32>,
    /// Consensus monoisotopic mass to adopt when *not* routing to advanced decon: the mean mass of
    /// the verdicts at [`consensus_k`](Self::consensus_k). `None` when no verdict was produced.
    pub consensus_mono: Option<f64>,
    /// `true` when the produced verdicts disagree on the offset — route to advanced multi-envelope
    /// decon. `false` for unanimous (or single/empty) results.
    pub needs_advanced: bool,
}

/// Locally sliced composite + apex spectra for one feature, plus the deconvolution m/z range.
struct FeatureSlices {
    comp_mz: Vec<f64>,
    comp_int: Vec<f64>,
    apex_mz: Vec<f64>,
    apex_int: Vec<f64>,
    /// Classic-decon range floor for the composite (`range_min0` clamped to the composite's first peak).
    range_min_comp: f64,
    /// Classic-decon range floor for the apex slice.
    range_min_apex: f64,
    /// Shared upper m/z bound of the envelope window; also the anchor-search upper bound for shift decon.
    range_max: f64,
    /// Lower bound of the anchor-search window for shift decon (`range_min0`, unclamped).
    anchor_min: f64,
}

/// Builds the averaged apex±1 composite and the apex-scan slice for a feature, over the same m/z
/// window [`refine_feature`] uses. Returns `None` on the same degenerate conditions
/// (`refine_feature` would also return `None`): empty scans/peaks, empty window, or empty composite.
fn build_feature_slices(
    feature: &DetectedFeature,
    scans: &[Scan],
    averaging_params: &SpectralAveragingParameters,
) -> Option<FeatureSlices> {
    if scans.is_empty() || feature.peaks.is_empty() {
        return None;
    }
    let apex = feature.apex_scan_index;
    let half = (MAX_SCANS_TO_AVERAGE / 2) as i32;
    let lo = (apex - half).max(0) as usize;
    let hi = ((apex + half).max(0) as usize).min(scans.len() - 1);
    if lo > hi {
        return None;
    }

    let spacing = C13_MINUS_C12 / feature.charge as f64;
    let range_min0 = (feature.mono_mz - 1.5).max(0.0);
    let range_max = feature.mono_mz + (feature.num_isotopes_observed as f64 + 3.0) * spacing + 1.0;
    let slice_lo = range_min0 - 0.5;
    let slice_hi = range_max + 0.5;

    let window = &scans[lo..=hi];
    let mut x_arrays: Vec<Vec<f64>> = Vec::with_capacity(window.len());
    let mut y_arrays: Vec<Vec<f64>> = Vec::with_capacity(window.len());
    for s in window {
        let a = s.mz.partition_point(|&m| m < slice_lo);
        let b = s.mz.partition_point(|&m| m <= slice_hi);
        x_arrays.push(s.mz[a..b].to_vec());
        y_arrays.push(s.intensity[a..b].to_vec());
    }
    if x_arrays.iter().all(|x| x.is_empty()) {
        return None;
    }

    let (comp_mz, comp_int) = average_spectra(&x_arrays, &y_arrays, averaging_params);
    if comp_mz.is_empty() {
        return None;
    }
    let range_min_comp = range_min0.max(comp_mz[0]);

    // Apex-scan slice over the same m/z window.
    let apex_usize = (feature.apex_scan_index.max(0) as usize).min(scans.len() - 1);
    let ascan = &scans[apex_usize];
    let a = ascan.mz.partition_point(|&m| m < slice_lo);
    let b = ascan.mz.partition_point(|&m| m <= slice_hi);
    let apex_mz = ascan.mz[a..b].to_vec();
    let apex_int = ascan.intensity[a..b].to_vec();
    let range_min_apex = range_min0.max(apex_mz.first().copied().unwrap_or(range_min0));

    Some(FeatureSlices {
        comp_mz,
        comp_int,
        apex_mz,
        apex_int,
        range_min_comp,
        range_min_apex,
        range_max,
        anchor_min: range_min0,
    })
}

/// Runs `classic_deconvolute` on a slice and returns the (mono, score) of the envelope at the
/// feature's charge whose mono m/z is closest to the feature's — the same pick [`refine_feature`]
/// makes. `None` when the slice is too small or no envelope matches the charge.
fn pick_classic_mono(
    mz: &[f64],
    inten: &[f64],
    feature: &DetectedFeature,
    range_min: f64,
    range_max: f64,
    decon_params: &ClassicDeconvolutionParameters,
) -> Option<(f64, f64)> {
    if mz.len() < 2 {
        return None;
    }
    classic_deconvolute(mz, inten, range_min, range_max, decon_params)
        .into_iter()
        .filter(|e| e.charge == feature.charge)
        .min_by(|a, b| {
            let da = (mass_to_mz_f64(a.monoisotopic_mass, a.charge) - feature.mono_mz).abs();
            let db = (mass_to_mz_f64(b.monoisotopic_mass, b.charge) - feature.mono_mz).abs();
            da.total_cmp(&db)
        })
        .map(|e| (e.monoisotopic_mass, e.score))
}

/// Builds a [`DeconVerdict`], computing the integer ¹³C offset from the detector's coarse mono.
fn make_verdict(
    view: DeconView,
    mono: f64,
    charge: i32,
    detector_mono: f64,
    confident: bool,
    score: f64,
) -> DeconVerdict {
    let offset_k = ((mono - detector_mono) / C13_MINUS_C12).round() as i32;
    DeconVerdict {
        view,
        monoisotopic_mass: mono,
        charge,
        offset_k,
        confident,
        score,
    }
}

/// Runs the four decon views for one detected feature and compares their monoisotope placements.
///
/// The two composite views reuse the exact averaging window and pick logic of [`refine_feature`];
/// the two apex views run on the single apex scan. `shift_tol_ppm` is the ppm tolerance for the
/// shift-decon peak matching (the classic views use `decon_params`).
/// How the shift-decon views choose their anchor (most-abundant reference) peak.
enum ShiftAnchor<'a> {
    /// Tallest peak in the envelope window (original; prone to grabbing a stronger neighbour).
    TallestInWindow,
    /// Tallest peak in the window that is not attributed to a co-eluting neighbour ([`shift_decon_gated`]).
    NeighborGated { forbidden: &'a [f64], grid_ppm: f64 },
    /// The detector's own most-abundant claimed peak (tallest of `feature.peaks`) — cannot grab any
    /// foreign peak because it uses the detector's own attribution of the feature's signal.
    DetectorPeak,
}

pub fn four_way_decon(
    feature: &DetectedFeature,
    scans: &[Scan],
    averaging_params: &SpectralAveragingParameters,
    decon_params: &ClassicDeconvolutionParameters,
    shift_tol_ppm: f64,
) -> FourWayDecon {
    four_way_decon_inner(
        feature,
        scans,
        averaging_params,
        decon_params,
        shift_tol_ppm,
        ShiftAnchor::TallestInWindow,
    )
}

/// **Detector-anchored** [`four_way_decon`]: the shift views anchor on the detector's own
/// most-abundant claimed peak (the tallest of `feature.peaks`), which can never be a foreign peak.
pub fn four_way_decon_detector_anchor(
    feature: &DetectedFeature,
    scans: &[Scan],
    averaging_params: &SpectralAveragingParameters,
    decon_params: &ClassicDeconvolutionParameters,
    shift_tol_ppm: f64,
) -> FourWayDecon {
    four_way_decon_inner(
        feature,
        scans,
        averaging_params,
        decon_params,
        shift_tol_ppm,
        ShiftAnchor::DetectorPeak,
    )
}

/// **Neighbor-aware** [`four_way_decon`]: the two shift views anchor with [`shift_decon_gated`],
/// excluding peaks that belong to an already-detected neighbour (any peak within `grid_ppm` of a
/// `forbidden_mz` position that is not on this feature's own isotope grid). `forbidden_mz` — the
/// predicted isotope m/z of co-eluting neighbours in this feature's window — must be ascending. The
/// classic views are unchanged (they are already tethered to the feature's mono and rarely drift).
pub fn four_way_decon_gated(
    feature: &DetectedFeature,
    scans: &[Scan],
    averaging_params: &SpectralAveragingParameters,
    decon_params: &ClassicDeconvolutionParameters,
    shift_tol_ppm: f64,
    forbidden_mz: &[f64],
    grid_ppm: f64,
) -> FourWayDecon {
    four_way_decon_inner(
        feature,
        scans,
        averaging_params,
        decon_params,
        shift_tol_ppm,
        ShiftAnchor::NeighborGated { forbidden: forbidden_mz, grid_ppm },
    )
}

fn four_way_decon_inner(
    feature: &DetectedFeature,
    scans: &[Scan],
    averaging_params: &SpectralAveragingParameters,
    decon_params: &ClassicDeconvolutionParameters,
    shift_tol_ppm: f64,
    anchor: ShiftAnchor,
) -> FourWayDecon {
    let detector_mono = feature.monoisotopic_mass;
    let charge = feature.charge;
    let mut verdicts: Vec<DeconVerdict> = Vec::new();
    // Own isotope grid for the gated anchor: peaks at `mono_mz + k·spacing`, k up to n_iso + 2.
    let own_kmax = feature.num_isotopes_observed + 2;
    // The detector's most-abundant claimed peak m/z (tallest of feature.peaks), for DetectorPeak.
    let detector_anchor_mz = feature
        .peaks
        .iter()
        .max_by(|a, b| a.intensity.total_cmp(&b.intensity))
        .map(|p| p.m() as f64);

    if let Some(s) = build_feature_slices(feature, scans, averaging_params) {
        // Classic composite.
        if let Some((mono, score)) = pick_classic_mono(
            &s.comp_mz,
            &s.comp_int,
            feature,
            s.range_min_comp,
            s.range_max,
            decon_params,
        ) {
            verdicts.push(make_verdict(
                DeconView::ClassicComposite,
                mono,
                charge,
                detector_mono,
                true,
                score,
            ));
        }
        // Classic apex.
        if let Some((mono, score)) = pick_classic_mono(
            &s.apex_mz,
            &s.apex_int,
            feature,
            s.range_min_apex,
            s.range_max,
            decon_params,
        ) {
            verdicts.push(make_verdict(
                DeconView::ClassicApex,
                mono,
                charge,
                detector_mono,
                true,
                score,
            ));
        }
        // Shift composite — anchor per the selected strategy.
        let sc = match &anchor {
            ShiftAnchor::TallestInWindow => shift_decon_in_window(
                &s.comp_mz, &s.comp_int, s.anchor_min, s.range_max, charge, shift_tol_ppm,
            ),
            ShiftAnchor::NeighborGated { forbidden, grid_ppm } => shift_decon_gated(
                &s.comp_mz, &s.comp_int, s.anchor_min, s.range_max, charge, shift_tol_ppm,
                feature.mono_mz, own_kmax, forbidden, *grid_ppm,
            ),
            ShiftAnchor::DetectorPeak => detector_anchor_mz
                .and_then(|a| shift_decon(&s.comp_mz, &s.comp_int, a, charge, shift_tol_ppm)),
        };
        if let Some(r) = sc {
            verdicts.push(make_verdict(
                DeconView::ShiftComposite,
                r.monoisotopic_mass,
                charge,
                detector_mono,
                r.shift0_passes_gate,
                r.shift0_correlation(),
            ));
        }
        // Shift apex.
        let sa = match &anchor {
            ShiftAnchor::TallestInWindow => shift_decon_in_window(
                &s.apex_mz, &s.apex_int, s.anchor_min, s.range_max, charge, shift_tol_ppm,
            ),
            ShiftAnchor::NeighborGated { forbidden, grid_ppm } => shift_decon_gated(
                &s.apex_mz, &s.apex_int, s.anchor_min, s.range_max, charge, shift_tol_ppm,
                feature.mono_mz, own_kmax, forbidden, *grid_ppm,
            ),
            ShiftAnchor::DetectorPeak => detector_anchor_mz
                .and_then(|a| shift_decon(&s.apex_mz, &s.apex_int, a, charge, shift_tol_ppm)),
        };
        if let Some(r) = sa {
            verdicts.push(make_verdict(
                DeconView::ShiftApex,
                r.monoisotopic_mass,
                charge,
                detector_mono,
                r.shift0_passes_gate,
                r.shift0_correlation(),
            ));
        }
    }

    analyze_agreement(detector_mono, charge, verdicts)
}

/// Tallies the verdicts' integer offsets, decides unanimity / plurality, and assembles the
/// [`FourWayDecon`]. Unanimous ⇔ all produced verdicts share one offset. When not unanimous the
/// plurality offset wins (ties → most confident verdicts at that offset, then smallest offset), and
/// `needs_advanced` is set. `consensus_mono` is the mean mass of the verdicts at the chosen offset.
fn analyze_agreement(
    detector_mono: f64,
    charge: i32,
    verdicts: Vec<DeconVerdict>,
) -> FourWayDecon {
    if verdicts.is_empty() {
        return FourWayDecon {
            detector_mono,
            charge,
            verdicts,
            unanimous: false,
            consensus_k: None,
            consensus_mono: None,
            needs_advanced: false,
        };
    }

    // Distinct offsets present.
    let first_k = verdicts[0].offset_k;
    let unanimous = verdicts.iter().all(|v| v.offset_k == first_k);

    // Per-offset tally: (count, confident-count) keyed by offset.
    let mut tally: HashMap<i32, (usize, usize)> = HashMap::new();
    for v in &verdicts {
        let e = tally.entry(v.offset_k).or_insert((0, 0));
        e.0 += 1;
        if v.confident {
            e.1 += 1;
        }
    }
    // Plurality: most verdicts, then most confident, then smallest offset.
    let consensus_k = tally
        .iter()
        .max_by(|a, b| {
            // a, b: (&offset, &(count, confident_count))
            a.1 .0
                .cmp(&b.1 .0) // count
                .then(a.1 .1.cmp(&b.1 .1)) // confident count
                .then(b.0.cmp(a.0)) // smaller offset wins (reversed so it ranks as "greater")
        })
        .map(|(k, _)| *k);

    let consensus_mono = consensus_k.map(|k| {
        let masses: Vec<f64> = verdicts
            .iter()
            .filter(|v| v.offset_k == k)
            .map(|v| v.monoisotopic_mass)
            .collect();
        masses.iter().sum::<f64>() / masses.len() as f64
    });

    FourWayDecon {
        detector_mono,
        charge,
        verdicts,
        unanimous,
        consensus_k,
        consensus_mono,
        needs_advanced: !unanimous,
    }
}

/// Resolves refined features into peptide features by grouping co-eluting charge states of the same
/// neutral mass and taking a cross-charge consensus on the monoisotopic mass.
///
/// Grouping links two features when they co-elute (`|apex_rt_i − apex_rt_j| ≤ rt_tolerance_minutes`)
/// AND their neutral masses agree within `mass_tolerance_ppm` after allowing for an integer ¹³C
/// off-by-one (see module docs). Linkage is single-linkage (transitive) via union-find. Within each
/// group the mass is resolved by pooling candidate masses across members and choosing the cluster
/// with the greatest cross-charge support.
pub fn resolve_charge_state_consensus(
    refined: &[RefinedFeature],
    mass_tolerance_ppm: f64,
    rt_tolerance_minutes: f64,
) -> Vec<ResolvedFeature> {
    if refined.is_empty() {
        return Vec::new();
    }

    group_features(refined, mass_tolerance_ppm, rt_tolerance_minutes)
        .into_iter()
        .map(|idxs| {
            let members: Vec<RefinedFeature> = idxs.iter().map(|&i| refined[i].clone()).collect();
            resolve_group(members, mass_tolerance_ppm)
        })
        .collect()
}

/// Groups refined features into single-linkage connected components under [`features_link`], returning
/// each component as its member indices. Component order and within-component order are ascending by
/// first-seen index, so the result depends **only** on the linkage partition — not on the order edges
/// were discovered.
///
/// Spatially pruned to ~O(n log n) instead of the naive O(n²) all-pairs scan. Both arms of
/// `features_link` are *local*: a partner co-elutes (`|Δapex_rt| ≤ rt_tol`) and its mass sits within a
/// bounded ±¹³C-off-by-one band of ours. So we bucket features by RT bin (`floor(apex_rt / rt_tol)`; a
/// partner within `rt_tol` lies in `bin ± 1`) with each bucket sorted by mass, then within those ≤3
/// buckets binary-search a mass window that is a *superset* of every off-by-one hit. The **exact**
/// `features_link` predicate is still applied to each surviving candidate — the window only prunes,
/// it never decides — so the connected components are byte-identical to the all-pairs version.
fn group_features(
    refined: &[RefinedFeature],
    mass_tolerance_ppm: f64,
    rt_tolerance_minutes: f64,
) -> Vec<Vec<usize>> {
    let n = refined.len();

    // Union-find.
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }

    // The RT binning divides by `rt_tol`; a non-positive or non-finite tolerance makes that
    // meaningless, so fall back to the exact all-pairs scan (this config does not occur in practice
    // and keeps exactness the priority).
    if !(rt_tolerance_minutes.is_finite() && rt_tolerance_minutes > 0.0) {
        for i in 0..n {
            for j in (i + 1)..n {
                if features_link(&refined[i], &refined[j], mass_tolerance_ppm, rt_tolerance_minutes) {
                    let ri = find(&mut parent, i);
                    let rj = find(&mut parent, j);
                    if ri != rj {
                        parent[ri] = rj;
                    }
                }
            }
        }
    } else {
        let rt_tol = rt_tolerance_minutes;
        let bin_of = |rt: f64| -> i64 { (rt / rt_tol).floor() as i64 };

        // Bucket indices by RT bin, each bucket sorted ascending by refined mass; carry a parallel
        // mass array so the candidate sub-range is a binary search.
        let mut buckets: HashMap<i64, Vec<usize>> = HashMap::new();
        for (i, feat) in refined.iter().enumerate() {
            buckets
                .entry(bin_of(feat.detected.apex_rt))
                .or_default()
                .push(i);
        }
        let mut bucket_masses: HashMap<i64, Vec<f64>> = HashMap::with_capacity(buckets.len());
        for (bin, idxs) in buckets.iter_mut() {
            idxs.sort_by(|&a, &b| {
                refined[a]
                    .refined_monoisotopic_mass
                    .total_cmp(&refined[b].refined_monoisotopic_mass)
            });
            bucket_masses.insert(
                *bin,
                idxs.iter()
                    .map(|&i| refined[i].refined_monoisotopic_mass)
                    .collect(),
            );
        }

        for i in 0..n {
            let mi = refined[i].refined_monoisotopic_mass;
            // Widest mass reach: the ±MAX_OFFBYONE_UNITS ¹³C band plus the ppm window (the ppm term
            // is relative to the partner mass, which at the band edge equals `mi`), plus a small
            // absolute epsilon to defend the binary-search boundary against float round-off. This is
            // a superset window — false candidates are removed by the exact predicate below.
            let w = MAX_OFFBYONE_UNITS as f64 * C13_MINUS_C12
                + mass_tolerance_ppm * 1e-6 * mi.abs()
                + 1e-6;
            let lo_mass = mi - w;
            let hi_mass = mi + w;
            let bi = bin_of(refined[i].detected.apex_rt);
            for b in (bi - 1)..=(bi + 1) {
                let Some(idxs) = buckets.get(&b) else {
                    continue;
                };
                let masses = &bucket_masses[&b];
                let start = masses.partition_point(|&m| m < lo_mass);
                let end = masses.partition_point(|&m| m <= hi_mass);
                for &j in &idxs[start..end] {
                    if j == i {
                        continue;
                    }
                    if features_link(
                        &refined[i],
                        &refined[j],
                        mass_tolerance_ppm,
                        rt_tolerance_minutes,
                    ) {
                        let ri = find(&mut parent, i);
                        let rj = find(&mut parent, j);
                        if ri != rj {
                            parent[ri] = rj;
                        }
                    }
                }
            }
        }
    }

    // Collect connected components, preserving first-seen order for determinism.
    let mut group_of: Vec<Option<usize>> = vec![None; n];
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for i in 0..n {
        let root = find(&mut parent, i);
        let g = match group_of[root] {
            Some(g) => g,
            None => {
                let g = groups.len();
                group_of[root] = Some(g);
                groups.push(Vec::new());
                g
            }
        };
        groups[g].push(i);
    }
    groups
}

/// Whether two refined features should be grouped: co-elution AND off-by-one-aware neutral-mass
/// agreement.
///
/// Co-elution requires the two **apexes** to fall within `rt_tolerance_minutes` — the signature of one
/// peptide seen at two charge states, which co-elute apex-to-apex. An earlier RT-range-*overlap* arm
/// was removed: since a detected feature spans ~1 min and single-linkage is transitive, overlapping
/// ranges daisy-chained long runs of features (apex A ≈ apex B's tail ≈ apex C's tail …) into one
/// component, inflating the resolved RT extent to many minutes (megagroups of 1000+ members). Charge
/// states of the same peptide share an apex, so apex proximity is the correct, chain-resistant test.
fn features_link(
    a: &RefinedFeature,
    b: &RefinedFeature,
    mass_tolerance_ppm: f64,
    rt_tolerance_minutes: f64,
) -> bool {
    if (a.detected.apex_rt - b.detected.apex_rt).abs() > rt_tolerance_minutes {
        return false;
    }
    let ma = a.refined_monoisotopic_mass;
    let mb = b.refined_monoisotopic_mass;
    for k in -MAX_OFFBYONE_UNITS..=MAX_OFFBYONE_UNITS {
        let shifted = mb + k as f64 * C13_MINUS_C12;
        if ppm_diff(ma, shifted) <= mass_tolerance_ppm {
            // k == 0 is a plain same-mass cross-charge match (the same peptide seen at another charge);
            // apex proximity is enough. k != 0 is an *off-by-one* bridge — two features whose masses
            // differ by 1–2 ¹³C. That is a real mono-placement disagreement only when they are the
            // same peptide, which co-elutes tightly; a different species sitting ~1 Da away that merely
            // drifts inside the loose apex window would otherwise be knitted as a spurious off-by-one
            // (and could then out-vote the true placement in consensus). So gate k != 0 on tight
            // co-elution: each apex must fall inside the other's elution bounds.
            if k == 0 || features_coelute_tightly(a, b) {
                return true;
            }
        }
    }
    false
}

/// Whether `a` and `b` share a chromatographic peak closely enough to be one peptide's off-by-one /
/// cross-charge pair, rather than two different species ~1 Da apart that merely pass the loose apex-RT
/// window. Requires **each** feature's apex to fall within the other's `[start_rt, end_rt]` elution
/// bounds (padded by [`COELUTION_PAD_MIN`] so a narrow, few-scan feature can still match a co-eluting
/// partner). Charge states of the same peptide share an apex and pass; a species eluting several scans
/// away has its apex outside the other's peak and fails.
fn features_coelute_tightly(a: &RefinedFeature, b: &RefinedFeature) -> bool {
    let a = &a.detected;
    let b = &b.detected;
    let within = |apex: f64, s: f64, e: f64| apex >= s - COELUTION_PAD_MIN && apex <= e + COELUTION_PAD_MIN;
    within(a.apex_rt, b.start_rt, b.end_rt) && within(b.apex_rt, a.start_rt, a.end_rt)
}

/// Resolves one group of grouped-and-cloned members into a [`ResolvedFeature`].
fn resolve_group(members: Vec<RefinedFeature>, mass_tolerance_ppm: f64) -> ResolvedFeature {
    // Distinct charges, ascending.
    let mut charge_states: Vec<i32> = members.iter().map(|m| m.refined_charge).collect();
    charge_states.sort_unstable();
    charge_states.dedup();

    // RT / intensity aggregates (from the detected extents).
    let start_rt = members
        .iter()
        .map(|m| m.detected.start_rt)
        .fold(f64::INFINITY, f64::min);
    let end_rt = members
        .iter()
        .map(|m| m.detected.end_rt)
        .fold(f64::NEG_INFINITY, f64::max);
    let summed_intensity: f64 = members.iter().map(|m| m.detected.summed_intensity).sum();
    // Apex from the most intense member.
    let tallest = members
        .iter()
        .max_by(|a, b| {
            a.detected
                .summed_intensity
                .total_cmp(&b.detected.summed_intensity)
        })
        .expect("group is non-empty");
    let apex_rt = tallest.detected.apex_rt;

    // Resolve the monoisotopic mass.
    let (monoisotopic_mass, cross_charge_support) = if charge_states.len() >= 2 {
        resolve_mass_by_cross_charge(&members, mass_tolerance_ppm)
    } else {
        // Single-charge (or singleton) group: fall back to the tallest member's refined mass.
        (tallest.refined_monoisotopic_mass, 1)
    };

    ResolvedFeature {
        monoisotopic_mass,
        charge_states,
        apex_rt,
        start_rt,
        end_rt,
        summed_intensity,
        cross_charge_support,
        members,
    }
}

/// A single pooled candidate mass tagged with its source charge, an intensity weight (for the
/// weighted-mean mass) and the source feature's envelope-fit `decon_score` (for the ranking tiebreak).
struct Candidate {
    mass: f64,
    charge: i32,
    weight: f64,
    fit: f64,
}

/// Pools every member's candidate masses, clusters them within `mass_tolerance_ppm` (0.01 Da floor),
/// and returns the (weighted-mean mass, cross-charge support) of the best cluster, ranked by
/// distinct-charge support, then candidate count, then **best envelope fit**, then summed intensity.
///
/// The fit tiebreak matters when two mass clusters tie on charge support — e.g. one peptide detected at
/// z=1 and z=2 whose z=1 envelope is contaminated by co-eluting interferents (refine mis-places its
/// mono, low `decon_score`) while its z=2 envelope is clean (correct mono, high score). Both are lone
/// singletons (support 1), so the earlier intensity-only tiebreak handed the mass to the *stronger* z=1
/// even though its placement is wrong; ranking the better-fitting placement first picks the clean z=2.
/// Support stays the **primary** key (a genuine multi-charge agreement must still outrank a lone
/// high-fit placement — making summed fit the primary key instead regressed recall).
fn resolve_mass_by_cross_charge(
    members: &[RefinedFeature],
    mass_tolerance_ppm: f64,
) -> (f64, usize) {
    let mut candidates: Vec<Candidate> = Vec::new();
    for m in members {
        // Use the member's summed intensity as the per-candidate weight (fall back to 1.0).
        let w = if m.detected.summed_intensity > 0.0 {
            m.detected.summed_intensity
        } else {
            1.0
        };
        let fit = if m.decon_score.is_finite() {
            m.decon_score.max(0.0)
        } else {
            0.0
        };
        for &mass in &m.candidate_masses {
            candidates.push(Candidate {
                mass,
                charge: m.refined_charge,
                weight: w,
                fit,
            });
        }
    }

    if candidates.is_empty() {
        // Nothing to intersect — fall back to the tallest member's refined mass.
        let tallest = members
            .iter()
            .max_by(|a, b| {
                a.detected
                    .summed_intensity
                    .total_cmp(&b.detected.summed_intensity)
            })
            .expect("group is non-empty");
        return (tallest.refined_monoisotopic_mass, 1);
    }

    candidates.sort_by(|a, b| a.mass.total_cmp(&b.mass));

    // Single-linkage clustering over the sorted candidates: break to a new cluster when the gap to
    // the previous candidate exceeds the local tolerance.
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    let mut current: Vec<usize> = vec![0];
    for i in 1..candidates.len() {
        let prev = candidates[i - 1].mass;
        let cur = candidates[i].mass;
        let tol = (cur.max(prev) * mass_tolerance_ppm / 1e6).max(MASS_CLUSTER_ABS_FLOOR_DA);
        if (cur - prev) <= tol {
            current.push(i);
        } else {
            clusters.push(std::mem::take(&mut current));
            current.push(i);
        }
    }
    clusters.push(current);

    // Score each cluster; pick the best. Tuple: (support, count, max_fit, sum_w, weighted_mass).
    let mut best: Option<(usize, usize, f64, f64, f64)> = None;
    for cluster in &clusters {
        let mut charges: Vec<i32> = cluster.iter().map(|&i| candidates[i].charge).collect();
        charges.sort_unstable();
        charges.dedup();
        let support = charges.len();
        let count = cluster.len();
        let sum_w: f64 = cluster.iter().map(|&i| candidates[i].weight).sum();
        // Best envelope fit in this cluster — the tiebreak when support and count are equal.
        let max_fit: f64 = cluster
            .iter()
            .map(|&i| candidates[i].fit)
            .fold(0.0_f64, f64::max);
        let weighted_mass = if sum_w > 0.0 {
            cluster
                .iter()
                .map(|&i| candidates[i].mass * candidates[i].weight)
                .sum::<f64>()
                / sum_w
        } else {
            cluster.iter().map(|&i| candidates[i].mass).sum::<f64>() / count as f64
        };
        // (support, count, max_fit, sum_w): support primary, then count, then best fit, then intensity.
        let better = match &best {
            None => true,
            Some((bs, bc, bf, bw, _)) => {
                support > *bs
                    || (support == *bs && count > *bc)
                    || (support == *bs && count == *bc && max_fit > *bf + 1e-9)
                    || (support == *bs
                        && count == *bc
                        && (max_fit - *bf).abs() <= 1e-9
                        && sum_w > *bw)
            }
        };
        if better {
            best = Some((support, count, max_fit, sum_w, weighted_mass));
        }
    }

    let (support, _, _, _, mass) = best.expect("at least one cluster exists");
    (mass, support)
}

/// Absolute ppm difference between two masses, relative to `b`.
#[inline]
fn ppm_diff(a: f64, b: f64) -> f64 {
    (a - b).abs() / b * 1e6
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::deconvolution::Polarity;
    use crate::isotopic_envelope::mass_to_mz_f64;
    use crate::peak_indexing::{IndexedMassSpectralPeak, PeakIndexingEngine};
    use crate::trace_kernel::{detect_features, poisson_comb_weights, TraceKernelParameters};

    /// Gaussian `exp(-½(Δ/σ)²)` — local copy (the one in `trace_kernel` is private).
    fn gaussian(delta: f64, sigma: f64) -> f64 {
        let z = delta / sigma;
        (-0.5 * z * z).exp()
    }

    /// Replicates `trace_kernel`'s `synthetic_envelope_scans`: a clean charge-2 isotope envelope
    /// (monoisotopic neutral mass 1000) eluting across 9 scans with a Gaussian RT profile
    /// (apex at scan 4). Returns the scans and the true mono m/z.
    fn synthetic_envelope_scans() -> (Vec<Scan>, f64, f64) {
        let mono_mass = 1000.0;
        let charge = 2;
        let mono_mz = mass_to_mz_f64(mono_mass, charge);
        let spacing = C13_MINUS_C12 / charge as f64;
        let weights = poisson_comb_weights(mono_mass, 1e-4, 8);
        let apex_intensity = 1.0e7;
        let rt_sigma = 0.15;

        let n_scans = 9;
        let apex_scan = 4;
        let mut scans = Vec::new();
        for s in 0..n_scans {
            let rt = 10.0 + s as f64 * 0.1;
            let g = gaussian(rt - (10.0 + apex_scan as f64 * 0.1), rt_sigma);
            let mut mz = Vec::new();
            let mut intensity = Vec::new();
            for (k, &wk) in weights.iter().enumerate() {
                mz.push(mono_mz + k as f64 * spacing);
                intensity.push(apex_intensity * wk * g);
            }
            scans.push(Scan {
                mz,
                intensity,
                one_based_scan_number: s + 1,
                retention_time: rt,
                msn_order: 1,
            });
        }
        (scans, mono_mz, mono_mass)
    }

    fn ppm_of(a: f64, b: f64) -> f64 {
        (a - b).abs() / b * 1e6
    }

    #[test]
    fn refine_recovers_precise_mass_and_candidates() {
        let (scans, _mono_mz, mono_mass) = synthetic_envelope_scans();
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let det_params = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.15,
            half_window_scans: 4,
            ..TraceKernelParameters::default()
        };
        let features = detect_features(&engine, &det_params);
        assert!(!features.is_empty(), "detector should find the envelope");
        let feature = &features[0];
        assert_eq!(feature.charge, 2);

        let avg = SpectralAveragingParameters::default();
        let decon = ClassicDeconvolutionParameters::new(1, 6, 20.0, 3.0, Polarity::Positive);
        let refined = refine_feature(feature, &scans, &avg, &decon)
            .expect("refinement should match the detected charge");

        assert_eq!(refined.refined_charge, 2);
        assert!(
            !refined.candidate_masses.is_empty(),
            "candidate masses should be surfaced"
        );
        assert!(
            ppm_of(refined.refined_monoisotopic_mass, mono_mass) <= 30.0,
            "refined mass {} not within 30 ppm of {} (ppm={})",
            refined.refined_monoisotopic_mass,
            mono_mass,
            ppm_of(refined.refined_monoisotopic_mass, mono_mass)
        );
    }

    /// Builds a `RefinedFeature` directly with controlled fields (bypassing detection/refinement),
    /// for consensus tests. The detected extent carries apex/start/end RT and summed intensity.
    fn make_refined(
        charge: i32,
        refined_mass: f64,
        candidate_masses: Vec<f64>,
        apex_rt: f64,
        intensity: f64,
    ) -> RefinedFeature {
        let mono_mz = mass_to_mz_f64(refined_mass, charge);
        let detected = DetectedFeature {
            monoisotopic_mass: refined_mass,
            charge,
            mono_mz,
            apex_scan_index: 0,
            apex_rt,
            start_rt: apex_rt - 0.2,
            end_rt: apex_rt + 0.2,
            summed_intensity: intensity,
            score: intensity,
            num_isotopes_observed: candidate_masses.len().max(2),
            peaks: Vec::<IndexedMassSpectralPeak>::new(),
        };
        RefinedFeature {
            detected,
            refined_monoisotopic_mass: refined_mass,
            refined_charge: charge,
            candidate_masses,
            decon_score: intensity,
        }
    }

    #[test]
    fn consensus_corrects_mono_off_by_one_across_charges() {
        let true_mass = 1000.0;
        let off = true_mass + C13_MINUS_C12; // +1.00335 Da monoisotopic off-by-one

        // Charge 2: refined correctly; candidates cluster around the true mass.
        let a = make_refined(
            2,
            true_mass,
            vec![true_mass, true_mass + 0.0008, true_mass - 0.0006],
            20.0,
            5.0e6,
        );
        // Charge 3: refined mass is off by +1, but its candidate list still contains the true mass
        // (alongside +1 predictions) — the classic decon's per-peak predictions disagree.
        let b = make_refined(
            3,
            off,
            vec![off, off + 0.0007, true_mass, true_mass + 0.0005],
            20.02,
            3.0e6,
        );

        // Grouping tolerance is tight (off-by-one-aware grouping bridges the 1 Da gap); candidate
        // clustering at the same 15 ppm keeps true-mass and +1 candidates in separate clusters.
        let resolved = resolve_charge_state_consensus(&[a, b], 15.0, 0.1);
        assert_eq!(resolved.len(), 1, "the two charges should form one feature");
        let r = &resolved[0];
        assert_eq!(r.charge_states, vec![2, 3]);
        assert_eq!(r.cross_charge_support, 2, "true mass is supported by both charges");
        assert!(
            ppm_of(r.monoisotopic_mass, true_mass) <= 5.0,
            "consensus mass {} should be the TRUE mass {} (ppm={})",
            r.monoisotopic_mass,
            true_mass,
            ppm_of(r.monoisotopic_mass, true_mass)
        );
    }

    #[test]
    fn offbyone_across_separate_elution_peaks_are_not_knitted() {
        // ECCHGDLLECADDRADLAK scenario: a real z=4 feature at the true mass, and a DIFFERENT species
        // exactly +1 ¹³C away (a co-eluting z=3 peptide that only *looks* like an off-by-one). Their
        // apexes are 0.099 min apart — inside the 0.1 min apex window — but the interloper's apex sits
        // OUTSIDE the true feature's narrow elution bounds, so they must stay two distinct features. If
        // knitted, the stronger interloper's +1 mass would win consensus and the true mono would be
        // lost (the observed off-by-one miss).
        let true_mass = 2246.9533;
        let interloper = true_mass + C13_MINUS_C12;
        let mut a = make_refined(4, true_mass, vec![true_mass], 12.153, 9.0e6);
        a.detected.start_rt = 12.153;
        a.detected.end_rt = 12.169; // narrow, weak, few-scan peak
        let mut b = make_refined(3, interloper, vec![interloper], 12.252, 2.0e8);
        b.detected.start_rt = 12.136;
        b.detected.end_rt = 12.252; // broad, strong, apex outside a's window
        let resolved = resolve_charge_state_consensus(&[a, b], 10.0, 0.1);
        assert_eq!(
            resolved.len(),
            2,
            "separately-eluting off-by-one species must not be knitted (got {} group(s))",
            resolved.len()
        );
    }

    #[test]
    fn coeluting_offbyone_charges_still_knit() {
        // Counterpart: a genuine same-peptide off-by-one across charges, sharing an apex, MUST still
        // knit (the co-elution gate only blocks the separately-eluting case above).
        let true_mass = 1500.0;
        let off = true_mass + C13_MINUS_C12;
        let a = make_refined(2, true_mass, vec![true_mass], 20.0, 5.0e6); // window [19.8, 20.2]
        let b = make_refined(3, off, vec![off, true_mass], 20.03, 3.0e6); // window [19.83, 20.23]
        let resolved = resolve_charge_state_consensus(&[a, b], 15.0, 0.1);
        assert_eq!(resolved.len(), 1, "co-eluting off-by-one charges should knit into one feature");
    }

    #[test]
    fn consensus_tiebreak_prefers_better_fit_over_intensity() {
        // AAVTAFWGK scenario: the peptide is seen at z=1 (envelope contaminated by co-eluting
        // interferents, so refine mis-places its mono +1 -> low fit) and z=2 (clean -> correct mono,
        // high fit). The z=1 is the STRONGER peak. Both are lone singletons (support 1), so the mass
        // vote is a tiebreak: it must take the better-fitting z=2 placement, not the intense z=1 +1.
        let true_mass = 992.5174;
        let z1_wrong = true_mass + C13_MINUS_C12;
        let mut z1 = make_refined(1, z1_wrong, vec![z1_wrong], 16.309, 9.4e7);
        z1.decon_score = 0.39; // contaminated envelope, poor fit
        let mut z2 = make_refined(2, true_mass, vec![true_mass], 16.310, 6.5e7);
        z2.decon_score = 0.99; // clean envelope
        let resolved = resolve_charge_state_consensus(&[z1, z2], 10.0, 0.1);
        assert_eq!(resolved.len(), 1, "same peptide at z1/z2 should knit into one feature");
        let r = &resolved[0];
        assert!(
            ppm_of(r.monoisotopic_mass, true_mass) <= 5.0,
            "consensus should take the better-fitting z=2 mass {} not the more intense z=1 +1 (got {})",
            true_mass,
            r.monoisotopic_mass
        );
    }

    #[test]
    fn four_way_unanimous_on_clean_synthetic_feature() {
        // A clean single-charge envelope: all four views should place the mono at the same offset,
        // so the result is unanimous and does NOT route to advanced decon.
        let (scans, _mono_mz, mono_mass) = synthetic_envelope_scans();
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let det_params = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.15,
            half_window_scans: 4,
            ..TraceKernelParameters::default()
        };
        let features = detect_features(&engine, &det_params);
        let feature = &features[0];

        let avg = SpectralAveragingParameters::default();
        let decon = ClassicDeconvolutionParameters::new(1, 6, 20.0, 3.0, Polarity::Positive);
        let fw = four_way_decon(feature, &scans, &avg, &decon, 20.0);

        assert!(!fw.verdicts.is_empty(), "clean feature should yield verdicts");
        assert!(fw.unanimous, "clean feature: all views agree (verdicts={:?})", fw.verdicts);
        assert!(!fw.needs_advanced, "unanimous ⇒ not routed to advanced");
        let mono = fw.consensus_mono.expect("consensus mono present");
        assert!(
            ppm_of(mono, mono_mass) <= 30.0,
            "consensus mono {} not within 30 ppm of {} (ppm={})",
            mono,
            mono_mass,
            ppm_of(mono, mono_mass)
        );
    }

    /// Builds a `DeconVerdict` directly for agreement-logic tests.
    fn verdict(view: DeconView, offset_k: i32, confident: bool) -> DeconVerdict {
        let detector_mono = 1000.0;
        let mono = detector_mono + offset_k as f64 * C13_MINUS_C12;
        DeconVerdict {
            view,
            monoisotopic_mass: mono,
            charge: 2,
            offset_k,
            confident,
            score: 1.0,
        }
    }

    #[test]
    fn agreement_unanimous_when_all_offsets_match() {
        let vs = vec![
            verdict(DeconView::ClassicComposite, 0, true),
            verdict(DeconView::ClassicApex, 0, true),
            verdict(DeconView::ShiftComposite, 0, true),
            verdict(DeconView::ShiftApex, 0, false),
        ];
        let fw = analyze_agreement(1000.0, 2, vs);
        assert!(fw.unanimous);
        assert!(!fw.needs_advanced);
        assert_eq!(fw.consensus_k, Some(0));
    }

    #[test]
    fn agreement_flags_disagreement_and_takes_plurality() {
        // Three views say k=0, one says k=-1 ⇒ disagreement (advanced) with plurality k=0.
        let vs = vec![
            verdict(DeconView::ClassicComposite, -1, true),
            verdict(DeconView::ClassicApex, 0, true),
            verdict(DeconView::ShiftComposite, 0, true),
            verdict(DeconView::ShiftApex, 0, false),
        ];
        let fw = analyze_agreement(1000.0, 2, vs);
        assert!(!fw.unanimous);
        assert!(fw.needs_advanced, "disagreement must route to advanced");
        assert_eq!(fw.consensus_k, Some(0), "plurality offset is 0");
        // consensus_mono is the mean of the three k=0 masses (all identical here).
        let expected = 1000.0;
        assert!((fw.consensus_mono.unwrap() - expected).abs() < 1e-9);
    }

    #[test]
    fn agreement_tie_breaks_toward_confident_then_smaller_offset() {
        // Two offsets each with one verdict; the +1 one is confident, the -1 one is not.
        // Equal count ⇒ confident count breaks the tie ⇒ +1 wins even though it is the larger offset.
        let vs = vec![
            verdict(DeconView::ClassicComposite, 1, true),
            verdict(DeconView::ShiftComposite, -1, false),
        ];
        let fw = analyze_agreement(1000.0, 2, vs);
        assert!(!fw.unanimous);
        assert_eq!(fw.consensus_k, Some(1), "confident offset wins the tie");
    }

    #[test]
    fn agreement_empty_is_not_advanced() {
        let fw = analyze_agreement(1000.0, 2, Vec::new());
        assert!(!fw.unanimous);
        assert!(!fw.needs_advanced);
        assert_eq!(fw.consensus_k, None);
        assert_eq!(fw.consensus_mono, None);
    }

    #[test]
    fn singleton_group_falls_back_to_refined_mass() {
        let m = make_refined(2, 1234.5678, vec![1234.5678, 1234.5690], 30.0, 1.0e6);
        let resolved = resolve_charge_state_consensus(std::slice::from_ref(&m), 15.0, 0.1);
        assert_eq!(resolved.len(), 1);
        let r = &resolved[0];
        assert_eq!(r.charge_states, vec![2]);
        assert_eq!(r.cross_charge_support, 1);
        assert_eq!(
            r.monoisotopic_mass, 1234.5678,
            "singleton resolves to its own refined mass"
        );
    }

    /// Reference O(n²) all-pairs grouping — the exact semantics the bucketed [`group_features`] must
    /// reproduce. Same union-find and same first-seen component collection, so equality is byte-for-byte.
    fn naive_groups(
        refined: &[RefinedFeature],
        mass_ppm: f64,
        rt_tol: f64,
    ) -> Vec<Vec<usize>> {
        let n = refined.len();
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(parent: &mut [usize], mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }
        for i in 0..n {
            for j in (i + 1)..n {
                if features_link(&refined[i], &refined[j], mass_ppm, rt_tol) {
                    let ri = find(&mut parent, i);
                    let rj = find(&mut parent, j);
                    if ri != rj {
                        parent[ri] = rj;
                    }
                }
            }
        }
        let mut group_of: Vec<Option<usize>> = vec![None; n];
        let mut groups: Vec<Vec<usize>> = Vec::new();
        for i in 0..n {
            let root = find(&mut parent, i);
            let g = match group_of[root] {
                Some(g) => g,
                None => {
                    let g = groups.len();
                    group_of[root] = Some(g);
                    groups.push(Vec::new());
                    g
                }
            };
            groups[g].push(i);
        }
        groups
    }

    #[test]
    fn bucketed_grouping_matches_naive_on_random_features() {
        // Deterministic LCG (numerical-recipes constants) — no rand dependency, reproducible.
        let mut state: u64 = 0x1234_5678_9abc_def0;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as f64 / (1u64 << 31) as f64 // in [0, 1)
        };

        let mass_ppm = 10.0;
        let rt_tol = 0.1;

        // Build ~2.5k features. Cluster masses around a few hundred base values, each drawn near a
        // base ± an integer ¹³C offset (0/±1/±2) plus ppm-scale jitter, and apex RTs bunched into a
        // handful of RT neighbourhoods — this deliberately exercises the off-by-one band, the RT-bin
        // boundaries, and multi-charge co-elution the bucketing must not miss.
        let n = 2500;
        let mut feats: Vec<RefinedFeature> = Vec::with_capacity(n);
        for _ in 0..n {
            let base = 600.0 + (next() * 300.0).floor() * 3.0; // discrete base masses ~600..1500
            let off = ((next() * 5.0).floor() as i32 - 2) as f64 * C13_MINUS_C12; // -2..+2 ¹³C
            let jitter = (next() - 0.5) * 2.0 * (mass_ppm * 1e-6 * base) * 1.5; // straddle the ppm edge
            let mass = base + off + jitter;
            let charge = 1 + (next() * 4.0).floor() as i32; // 1..4
            // RT bunched into ~15 neighbourhoods, each a few multiples of rt_tol wide, so groups form.
            let hub = (next() * 15.0).floor() * (rt_tol * 4.0) + 10.0;
            let apex_rt = hub + (next() - 0.5) * 2.0 * rt_tol * 1.5; // straddle the ±1-bin boundary
            feats.push(make_refined(charge, mass, vec![mass], apex_rt, 1.0e6 * (1.0 + next())));
        }

        let got = group_features(&feats, mass_ppm, rt_tol);
        let want = naive_groups(&feats, mass_ppm, rt_tol);
        assert_eq!(
            got, want,
            "bucketed grouping must equal naive O(n²) grouping (components and order)"
        );

        // Sanity: the fixture actually produced non-trivial structure (some multi-member groups),
        // otherwise the test would pass vacuously.
        assert!(
            want.iter().any(|g| g.len() >= 2),
            "fixture should form at least one multi-member group"
        );
        assert!(want.len() < feats.len(), "fixture should merge at least some features");
    }
}
