//! Trace-kernel feature detection — the untargeted MS1 detector.
//!
//! This is the **new** algorithm (no mzLib counterpart, so no C# golden): it locates
//! "peptide-shaped objects" directly in the raw MS1 data via a sparse 2D matched filter — an
//! **isotope comb** in the m/z dimension × a **Gaussian** in the retention-time dimension — scored
//! against the [`crate::peak_indexing::PeakIndexingEngine`]. See
//! `agent_info/Feature-Detection-Design.md` ("Trace kernel (detection)") for the rationale.
//!
//! ## What it does
//! For each seed peak (tallest first, greedy claim-as-you-go), it scores charge hypotheses
//! `z = min..=max`. Each hypothesis lays an isotope comb spaced `(C13 − C12)/z` in m/z and weighted
//! by the expected isotope envelope, times a Gaussian in RT centred on the seed's scan, and sums the
//! observed intensity the comb lands on. The best-scoring `z` wins (**cross-z non-max
//! suppression**); its peaks are claimed so overlapping harmonics (a real z=2 is a subset of z=4/z=6
//! combs) cannot re-fire. Accepted hypotheses become [`DetectedFeature`] records.
//!
//! ## Design decisions realised here
//! - **The Gaussian template *is* the shape test** — there is no separate data-vs-data correlation
//!   gate at detection time. A matched filter degrades gracefully on tailed (real) peaks; tail-aware
//!   integration bounds are `cut_peak`'s job, downstream.
//! - **The comb runs both directions** from the seed: the seed is the *most intense* peak, which for
//!   heavier masses is not the monoisotopic one, so the mono is placed at `seed − i*·spacing/z` where
//!   `i*` is the most-abundant isotope index of the envelope model.
//! - **Comb weights: the averagine table is the default.** The comb teeth are weighted by the real
//!   averagine isotope envelope ([`crate::deconvolution::averagine_comb_weights`], a table lookup),
//!   which places the monoisotope correctly across the mass range — including near ~1.8 kDa where the
//!   envelope mode shifts off the monoisotope and a single-parameter Poisson `i*` can be off by one
//!   ¹³C unit. The closed-form `Poisson(λ = 0.00048·M)` ([`poisson_comb_weights`]) is retained as a
//!   faster, table-free alternative selectable via [`CombWeightModel`].
//! - **Evaluate sparsely.** Only the comb's expected `(m/z, scan)` points are looked up in the index;
//!   nothing is rasterised.
//!
//! ## Not yet here (follow-ups)
//! - **Detect-then-refine**: averaging the RT window ([`crate::spectral_averaging`]) + a final
//!   [`crate::deconvolution`] pass on the composite. This module is the *detector/assembler*; the
//!   refiner is wired separately.
//! - **Charge-state consensus** across co-eluting z of the same neutral mass (needs the
//!   `DeconEnvelope` candidate-mass list). Grouping here is per-hypothesis, one charge at a time.
//! - **Averagine comb weights** (benchmark alternative to Poisson).

use std::collections::HashSet;

use crate::isotopic_envelope::{C13_MINUS_C12, PROTON_MASS};
use crate::peak_indexing::{IndexedMassSpectralPeak, PeakIndexingEngine, PeakKey, ScanInfo};
use crate::tolerance::PpmTolerance;

/// Poisson rate per dalton for the closed-form comb: `λ ≈ 0.00048·M`. This is (carbons per Da)
/// × (¹³C natural abundance) ≈ `(1 / averagineUnitMass · averageC) · 0.0107`, i.e. how many ¹³C
/// substitutions a peptide of mass `M` carries on average. The Poisson in that count *is* the
/// isotope envelope.
pub const POISSON_LAMBDA_PER_DA: f64 = 0.00048;

/// Full-width-at-half-maximum → Gaussian σ conversion factor: `FWHM = 2·√(2·ln2)·σ ≈ 2.3548·σ`.
pub const FWHM_TO_SIGMA: f64 = 2.354_820_045_030_949;

/// How the isotope comb's per-peak weights are produced.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CombWeightModel {
    /// Closed-form `Poisson(λ = 0.00048·M)` over the ¹³C-substitution count. One parameter, no table
    /// lookup; a faster, table-free alternative to [`Self::Averagine`].
    Poisson,
    /// The real averagine isotope envelope
    /// ([`crate::deconvolution::averagine_comb_weights`]), a table lookup. **The default.** More
    /// accurate than Poisson near ~1.8 kDa where the envelope mode shifts off the monoisotope — the
    /// regime where a Poisson `i*` can misplace the monoisotope by one ¹³C unit (off-by-one).
    Averagine,
}

/// How a charge hypothesis's matched-filter response is scored for cross-z non-max suppression.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoreModel {
    /// Raw inner product `Σ (wₖ·gₛ·I)` over the observed comb slots. **Unnormalised** — it grows with
    /// how many teeth/scans a hypothesis spans, so a broad or higher-charge comb can out-score the
    /// correct one just by covering more of the window. The original detector score; kept as default
    /// until the normalised model is validated to not regress.
    RawSum,
    /// **Noise-floor-truncated normalised correlation** (Change B). `score = Σ_S(wₖ·gₛ·I) /
    /// sqrt(Σ_S (wₖ·gₛ)²)` over the *expected-observable support* `S = { (k,s) : A·wₖ·gₛ ≥ η }`, where
    /// `A` is the matched-filter least-squares apex amplitude and `η` is the run-level noise floor
    /// ([`TraceKernelParameters::noise_floor`]). Slots the model predicts fall **below** the noise
    /// floor are dropped from both numerator and denominator (a faint real peak is not penalised for
    /// teeth the instrument could never record); a slot in `S` with no observed peak stays in the
    /// denominator and correctly penalises (a real miss). `η = 0` degenerates to a full-template norm.
    NormalizedNoiseFloor,
}

/// Parameters governing the trace-kernel detector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TraceKernelParameters {
    /// Lowest charge hypothesis to test (bottom-up default: 1).
    pub min_charge: i32,
    /// Highest charge hypothesis to test (bottom-up default: 6).
    pub max_charge: i32,
    /// m/z match tolerance, ppm, for looking up comb teeth in the index.
    pub ppm_tolerance: f64,
    /// Gaussian σ along retention time, in **minutes** (the same unit as `ScanInfo::retention_time`).
    pub rt_sigma_minutes: f64,
    /// **Legacy / superseded by [`Self::rt_half_window_minutes`].** Formerly bounded the window by a
    /// fixed scan count; retained for API/compat but no longer consumed by the detector (the window
    /// is now a time window). Still set by `with_rt_from_scans` for reference.
    pub half_window_scans: i32,
    /// Which comb-weight model to use.
    pub weight_model: CombWeightModel,
    /// How the hypothesis response is scored (raw sum vs normalised; see [`ScoreModel`]).
    pub score_model: ScoreModel,
    /// Run-level noise floor `η` for [`ScoreModel::NormalizedNoiseFloor`] — comb slots whose
    /// model-predicted intensity `A·wₖ·gₛ` falls below this are treated as unobservable and dropped
    /// from the normalised score. `0.0` disables the truncation (full-template norm). Ignored by
    /// [`ScoreModel::RawSum`]. Set from the data (e.g. [`estimate_noise_floor`]).
    pub noise_floor: f64,
    /// Apex-amplitude `A` estimator for [`ScoreModel::NormalizedNoiseFloor`]'s support set. `false`
    /// (default) = matched-filter least-squares `A = Σ(t·I)/Σ(t²)` over observed slots; `true` = the
    /// seed (most-abundant tooth) intensity. Only affects which slots clear `A·t ≥ η`.
    pub score_use_seed_amplitude: bool,
    /// If `true`, additionally divide the normalised score by `‖I‖` over the support `S`, giving a
    /// bounded `[0, 1]` cosine shape-fit instead of the template-normalised correlation. Only applies
    /// to [`ScoreModel::NormalizedNoiseFloor`]. Default `false`.
    pub score_cosine: bool,
    /// Smallest envelope weight (relative to the tallest = 1) that still contributes a comb tooth.
    pub min_isotope_weight: f64,
    /// Hard cap on the number of comb teeth (isotopes) considered.
    pub max_isotopes: usize,
    /// Minimum number of distinct comb teeth that must be observed for a hypothesis to be accepted.
    /// Two (a doublet) is the floor — a lone peak is not a feature.
    pub min_isotopes_observed: usize,
    /// Minimum number of **distinct scans** the accepted feature's *traced extent* (Change A) must
    /// span — a chromatographic-persistence gate. A real peptide elutes over several scans; a feature
    /// claiming only one scan is a noise doublet (two peaks at comb spacing that happened to co-occur
    /// in a single scan), not an elution. `1` disables the gate (any accepted hypothesis passes).
    pub min_feature_scans: usize,
    /// Seeds with intensity below this are not considered (noise floor). `0.0` disables the floor
    /// and seeds from every peak. Because seeds are visited intensity-descending, this also bounds
    /// runtime: the seed loop stops as soon as it drops below the floor.
    pub min_seed_intensity: f64,
    /// Stop once this fraction of the total MS1 intensity has been explained (claimed by accepted
    /// features). This is the design's primary stopping criterion — "explain most of the big
    /// signal, not every peak". Denominator = Σ of all indexed peak intensities. `1.0` (or more)
    /// disables the cap and detects until seeds are exhausted / fall below `min_seed_intensity`.
    pub coverage_target: f64,
    /// Half-width of the retention-time window (minutes) the matched filter evaluates around the
    /// seed apex. This bounds the window **in time**, not in scan count: DDA interleaves a variable
    /// number of MS2 scans between MS1 scans, so a fixed scan-count window spans wildly different
    /// times — producing over-wide features that over-claim and split one elution into several.
    /// Typically ≈ 2σ (`with_rt_from_scans` sets it there). Supersedes `half_window_scans`.
    pub rt_half_window_minutes: f64,
    /// Consecutive-miss tolerance for the **claim-extent XIC trace** (Change A). When an accepted
    /// hypothesis claims its true elution, the most-abundant tooth is followed outward in RT via
    /// [`crate::peak_indexing::PeakIndexingEngine::get_xic_by_scan_index`]; the walk stops after this
    /// many consecutive scans with no matching peak. Small (1) so genuinely co-eluting neighbours of
    /// the same m/z (separated by a valley) are not merged. Independent of the scoring window.
    pub trace_missed_scans_allowed: i32,
    /// Maximum RT half-width (minutes) the claim-extent trace may reach from the apex — the runaway
    /// guard on the XIC walk. Decoupled from (and much wider than) the ~2σ *scoring* window: the whole
    /// point of Change A is to claim a real elution wider than 2σ, so this bounds only pathological
    /// traces, not real peaks.
    pub trace_max_half_width_minutes: f64,
}

impl Default for TraceKernelParameters {
    /// Bottom-up defaults: charge 1–6, 10 ppm, averagine weights, ≥2 observed isotopes. The RT σ and
    /// window are left at ~6 s / ±3 scans placeholders — callers should set them from the data via
    /// [`TraceKernelParameters::with_rt_from_scans`].
    fn default() -> Self {
        TraceKernelParameters {
            min_charge: 1,
            max_charge: 6,
            ppm_tolerance: 10.0,
            rt_sigma_minutes: 0.1,
            half_window_scans: 3,
            weight_model: CombWeightModel::Averagine,
            score_model: ScoreModel::RawSum,
            noise_floor: 0.0,
            score_use_seed_amplitude: false,
            score_cosine: false,
            min_isotope_weight: 1e-3,
            max_isotopes: 12,
            min_isotopes_observed: 2,
            min_feature_scans: 2,
            min_seed_intensity: 0.0,
            coverage_target: 1.0,
            rt_half_window_minutes: 0.5,
            trace_missed_scans_allowed: 1,
            trace_max_half_width_minutes: 0.5,
        }
    }
}

impl TraceKernelParameters {
    /// Derives the RT σ and scan half-window from the data, given an **assumed** chromatographic peak
    /// width. `assumed_fwhm_seconds` is the design's "~36 s peaks" starting assumption; the σ is
    /// `FWHM / 2.3548` and the half-window spans ±2σ in scans, using the median MS1 scan spacing.
    ///
    /// The FWHM (not the full peak width) drives the averaging window so co-eluting neighbours are
    /// not pulled into the composite; the same σ is reused as the detector's RT Gaussian width.
    ///
    /// Prefer [`Self::with_rt_from_index`] when the built index is available — it *measures* the FWHM
    /// from the data instead of assuming it.
    pub fn with_rt_from_scans(self, scan_info: &[ScanInfo], assumed_fwhm_seconds: f64) -> Self {
        self.apply_fwhm(scan_info, assumed_fwhm_seconds)
    }

    /// Derives the RT σ and window from the run's **measured** chromatographic FWHM (Change A).
    ///
    /// Estimates the true FWHM from XIC half-max over a bounded sample of the tallest clean XICs (see
    /// [`estimate_fwhm_seconds`]), clamps it to a sane `[FWHM_FLOOR_SEC, FWHM_CEIL_SEC]` band so a
    /// pathological run cannot drive σ to a degenerate value, and sets σ / windows from it. Falls back
    /// to `fallback_fwhm_seconds` (the assumed-FWHM path) when too few clean XICs are found.
    ///
    /// Unlike [`Self::with_rt_from_scans`], this needs the **built** [`PeakIndexingEngine`], so it
    /// introduces an ordering dependency: build the index → finalize params with this → detect.
    pub fn with_rt_from_index(self, engine: &PeakIndexingEngine, fallback_fwhm_seconds: f64) -> Self {
        let ppm = PpmTolerance::new(self.ppm_tolerance);
        let fwhm_seconds = estimate_fwhm_seconds(engine, &ppm)
            .unwrap_or(fallback_fwhm_seconds)
            .clamp(FWHM_FLOOR_SEC, FWHM_CEIL_SEC);
        self.apply_fwhm(engine.scan_info(), fwhm_seconds)
    }

    /// Sets σ, the scan half-window, and the RT time-window from a chromatographic FWHM (seconds).
    /// Shared by the assumed-FWHM ([`Self::with_rt_from_scans`]) and measured-FWHM
    /// ([`Self::with_rt_from_index`]) constructors. Leaves the claim-extent trace knobs alone — those
    /// are deliberately independent of the scoring σ.
    fn apply_fwhm(mut self, scan_info: &[ScanInfo], fwhm_seconds: f64) -> Self {
        let sigma_minutes = (fwhm_seconds / 60.0) / FWHM_TO_SIGMA;
        let spacing = median_ms1_scan_spacing_minutes(scan_info).max(f64::MIN_POSITIVE);
        self.rt_sigma_minutes = sigma_minutes;
        self.half_window_scans = ((2.0 * sigma_minutes) / spacing).round().max(1.0) as i32;
        // The matched filter is bounded in *time* (see `rt_half_window_minutes`); ±2σ covers the peak.
        self.rt_half_window_minutes = 2.0 * sigma_minutes;
        self
    }
}

/// Lower clamp (seconds) for the measured-FWHM estimate — below this, σ would be so tight the RT
/// Gaussian is essentially a delta and the window collapses to the apex scan.
pub const FWHM_FLOOR_SEC: f64 = 1.0;
/// Upper clamp (seconds) for the measured-FWHM estimate — above this we distrust the measurement (a
/// pathological / co-eluting-dominated run) and cap it.
pub const FWHM_CEIL_SEC: f64 = 60.0;

/// Target number of clean XIC FWHM measurements to accumulate before taking the median.
const FWHM_PROBE_SAMPLE_TARGET: usize = 500;
/// Hard cap on seeds examined by the probe, so a run of mostly-unmeasurable XICs still returns
/// promptly (bounded startup cost regardless of how many clean XICs exist).
const FWHM_PROBE_MAX_SEEDS: usize = 20_000;
/// Minimum clean measurements required to trust the median; below this the probe returns `None` and
/// the caller falls back to the assumed FWHM.
const FWHM_PROBE_MIN_SAMPLE: usize = 12;

/// Estimates the run's chromatographic FWHM (**seconds**) from XIC half-max, over a bounded sample of
/// the tallest clean XICs. Returns `None` when fewer than [`FWHM_PROBE_MIN_SAMPLE`] clean XICs are
/// measurable (the caller then uses its assumed-FWHM fallback).
///
/// Bounded by design: seeds are visited tallest-first (each surviving XIC's peaks are marked so later
/// seeds skip them, mirroring `get_all_xics`), stopping once [`FWHM_PROBE_SAMPLE_TARGET`] clean
/// measurements are collected or [`FWHM_PROBE_MAX_SEEDS`] seeds have been examined. The median (not
/// the mean) is returned, to resist tails and co-elution.
pub fn estimate_fwhm_seconds(engine: &PeakIndexingEngine, ppm: &PpmTolerance) -> Option<f64> {
    let mut seeds = engine.all_peaks();
    seeds.sort_by(|a, b| b.intensity.total_cmp(&a.intensity));

    let mut claimed: HashSet<PeakKey> = HashSet::new();
    let mut widths: Vec<f64> = Vec::new();
    let mut examined = 0usize;

    for seed in &seeds {
        if widths.len() >= FWHM_PROBE_SAMPLE_TARGET || examined >= FWHM_PROBE_MAX_SEEDS {
            break;
        }
        if claimed.contains(&seed.key()) {
            continue;
        }
        examined += 1;
        // Generous RT cap (2 min) so a real peak is never clipped before its half-max shoulders.
        let xic = engine.get_xic_by_scan_index(
            seed.m() as f64,
            seed.zero_based_scan_index,
            ppm,
            1,
            2.0,
            Some(&claimed),
        );
        for p in &xic {
            claimed.insert(p.key());
        }
        if let Some(w) = xic_fwhm_minutes(&xic) {
            widths.push(w);
        }
    }

    if widths.len() < FWHM_PROBE_MIN_SAMPLE {
        return None;
    }
    widths.sort_by(|a, b| a.total_cmp(b));
    let n = widths.len();
    let median_minutes = if n % 2 == 1 {
        widths[n / 2]
    } else {
        (widths[n / 2 - 1] + widths[n / 2]) / 2.0
    };
    Some(median_minutes * 60.0)
}

/// FWHM (minutes) of one XIC via linear-interpolated half-max crossings. `None` if the trace has
/// fewer than 3 points, the apex sits at an edge (not a real rise-then-fall), or it does not fall
/// below half-max on both sides. Peaks must be RT-ascending (as `get_xic_by_scan_index` returns).
fn xic_fwhm_minutes(xic: &[IndexedMassSpectralPeak]) -> Option<f64> {
    let n = xic.len();
    if n < 3 {
        return None;
    }
    let pts: Vec<(f64, f64)> = xic
        .iter()
        .map(|p| (p.retention_time as f64, p.intensity as f64))
        .collect();
    // Apex = max-intensity sample; require it internal (a genuine rise-then-fall peak).
    let mut ai = 0usize;
    for i in 1..n {
        if pts[i].1 > pts[ai].1 {
            ai = i;
        }
    }
    if ai == 0 || ai == n - 1 {
        return None;
    }
    let half = pts[ai].1 / 2.0;
    if half <= 0.0 {
        return None;
    }
    // Left crossing: nearest sample left of apex at or below half, interpolated to `half`.
    let mut left = None;
    for i in (0..ai).rev() {
        if pts[i].1 <= half {
            let (t0, y0) = pts[i];
            let (t1, y1) = pts[i + 1];
            left = Some(if y1 != y0 {
                t0 + (half - y0) * (t1 - t0) / (y1 - y0)
            } else {
                t0
            });
            break;
        }
    }
    // Right crossing.
    let mut right = None;
    for i in (ai + 1)..n {
        if pts[i].1 <= half {
            let (t0, y0) = pts[i - 1];
            let (t1, y1) = pts[i];
            right = Some(if y1 != y0 {
                t0 + (half - y0) * (t1 - t0) / (y1 - y0)
            } else {
                t1
            });
            break;
        }
    }
    match (left, right) {
        (Some(l), Some(r)) if r > l => Some(r - l),
        _ => None,
    }
}

/// A detected untargeted MS1 feature: one charge state's isotope envelope traced across RT.
#[derive(Debug, Clone)]
pub struct DetectedFeature {
    /// Monoisotopic neutral mass inferred from the mono comb position and charge.
    pub monoisotopic_mass: f64,
    /// Charge state (the winning hypothesis).
    pub charge: i32,
    /// m/z of the monoisotopic comb tooth (`seed − i*·spacing/z`).
    pub mono_mz: f64,
    /// Zero-based scan index of the feature's most intense claimed peak.
    pub apex_scan_index: i32,
    /// Retention time of the apex, minutes.
    pub apex_rt: f64,
    /// Earliest RT among claimed peaks, minutes.
    pub start_rt: f64,
    /// Latest RT among claimed peaks, minutes.
    pub end_rt: f64,
    /// Sum of claimed peak intensities (the feature's explained signal).
    pub summed_intensity: f64,
    /// Matched-filter response — the detector's score for this feature.
    pub score: f64,
    /// Number of distinct isotope teeth that contributed at least one observed peak.
    pub num_isotopes_observed: usize,
    /// The peaks this feature claims (across the RT window and all matched isotopes).
    pub peaks: Vec<IndexedMassSpectralPeak>,
}

/// Neutral mass from an m/z at a given charge: `|z|·mz − z·ProtonMass` (f64). Matches
/// `ClassExtensions.ToMass` / [`crate::isotopic_envelope::mz_to_mass_f32`] but keeps f64 precision.
#[inline]
fn mz_to_mass(mz: f64, charge: i32) -> f64 {
    (charge.abs() as f64) * mz - (charge as f64) * PROTON_MASS
}

/// Closed-form Poisson comb weights for a peptide of neutral mass `neutral_mass`.
///
/// Returns the isotope envelope `w[k] = e^{-λ} λ^k / k!` (`λ = 0.00048·neutral_mass`) as a
/// probability-mass vector, truncated once a weight falls below `min_weight × w_max` (past the mode)
/// or `max_isotopes` teeth are reached. The vector is **normalised so its maximum weight is 1.0**,
/// which makes the seed (the tallest observed peak) align naturally with the tallest template tooth.
pub fn poisson_comb_weights(neutral_mass: f64, min_weight: f64, max_isotopes: usize) -> Vec<f64> {
    let lambda = POISSON_LAMBDA_PER_DA * neutral_mass.max(0.0);
    let mut weights: Vec<f64> = Vec::new();
    // w[0] = e^{-λ}; w[k] = w[k-1]·λ/k. Build up to max_isotopes, tracking the max for normalisation.
    let mut w = (-lambda).exp();
    let mut max_w = w;
    weights.push(w);
    for k in 1..max_isotopes {
        w = w * lambda / (k as f64);
        weights.push(w);
        if w > max_w {
            max_w = w;
        }
        // Stop once we are past the mode (weights descending) and below the relative floor.
        if w < min_weight * max_w && w < weights[k - 1] {
            break;
        }
    }
    if max_w > 0.0 {
        for wk in weights.iter_mut() {
            *wk /= max_w;
        }
    }
    weights
}

/// The isotope comb weights for a neutral mass under the configured [`CombWeightModel`]. Shared by
/// the scorer and the claim-extent tracer so both lay the identical comb.
fn comb_weights(neutral_mass: f64, params: &TraceKernelParameters) -> Vec<f64> {
    match params.weight_model {
        CombWeightModel::Poisson => {
            poisson_comb_weights(neutral_mass, params.min_isotope_weight, params.max_isotopes)
        }
        CombWeightModel::Averagine => crate::deconvolution::averagine_comb_weights(
            neutral_mass,
            params.min_isotope_weight,
            params.max_isotopes,
        ),
    }
}

/// Index of the most-abundant (tallest) tooth in a weight vector. Ties resolve to the lower index.
fn most_abundant_index(weights: &[f64]) -> usize {
    let mut best = 0;
    for (i, &w) in weights.iter().enumerate() {
        if w > weights[best] {
            best = i;
        }
    }
    best
}

/// Gaussian value `exp(-½ (Δ/σ)²)`. A non-positive σ degenerates to a delta function (only the
/// apex, `Δ == 0`, contributes) rather than dividing by zero and poisoning the response with `NaN`.
#[inline]
fn gaussian(delta: f64, sigma: f64) -> f64 {
    if sigma <= 0.0 {
        return if delta == 0.0 { 1.0 } else { 0.0 };
    }
    let z = delta / sigma;
    (-0.5 * z * z).exp()
}

/// Median spacing between consecutive MS1 scan retention times (minutes). Returns 0.0 for < 2 scans.
///
/// `scan_info` is assumed ordered by scan (as the index builds it). Uses the standard median
/// convention (average of the two middle order statistics for an even count).
pub fn median_ms1_scan_spacing_minutes(scan_info: &[ScanInfo]) -> f64 {
    if scan_info.len() < 2 {
        return 0.0;
    }
    let mut diffs: Vec<f64> = scan_info
        .windows(2)
        .map(|w| w[1].retention_time - w[0].retention_time)
        .collect();
    diffs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = diffs.len();
    if n % 2 == 1 {
        diffs[n / 2]
    } else {
        (diffs[n / 2 - 1] + diffs[n / 2]) / 2.0
    }
}

/// The outcome of scoring one `(seed, charge)` hypothesis.
struct HypothesisScore {
    response: f64,
    charge: i32,
    mono_mz: f64,
    peaks: Vec<IndexedMassSpectralPeak>,
    num_isotopes_observed: usize,
}

/// The RT window a seed's matched filter evaluates: `(scan_index, gaussian_weight)` for every scan
/// within `rt_half_window_minutes` of the seed apex. Walks outward from the apex in both directions,
/// stopping as soon as RT leaves the window — scans are RT-ordered, so this is a bounded walk that
/// adapts to the local (uneven) MS1 spacing instead of a fixed scan count.
///
/// This is **charge-independent**, so `detect_features` computes it once per seed and shares it
/// across every charge hypothesis (the Gaussian weight likewise depends only on the seed apex).
fn seed_rt_window(
    engine: &PeakIndexingEngine,
    seed: &IndexedMassSpectralPeak,
    params: &TraceKernelParameters,
) -> Vec<(i32, f64)> {
    let scan_info = engine.scan_info();
    let n_scans = scan_info.len() as i32;
    let apex = seed.zero_based_scan_index;
    let rt_apex = seed.retention_time as f64;
    let rt_win = params.rt_half_window_minutes;
    let sigma = params.rt_sigma_minutes;

    let mut window: Vec<(i32, f64)> = Vec::new();
    let mut s = apex;
    while s >= 0 && (scan_info[s as usize].retention_time - rt_apex).abs() <= rt_win {
        window.push((s, gaussian(scan_info[s as usize].retention_time - rt_apex, sigma)));
        s -= 1;
    }
    let mut s = apex + 1;
    while s < n_scans && (scan_info[s as usize].retention_time - rt_apex).abs() <= rt_win {
        window.push((s, gaussian(scan_info[s as usize].retention_time - rt_apex, sigma)));
        s += 1;
    }
    window
}

/// Scores a single charge hypothesis for a seed peak: lays the isotope comb (anchored so the
/// most-abundant tooth sits on the seed), evaluates the RT Gaussian across the scan window, and
/// sums `weight · gaussian · observed_intensity` over every comb `(m/z, scan)` point. Peaks already
/// in `claimed` are treated as absent (this is what makes cross-feature NMS work).
///
/// `window` is the seed's precomputed [`seed_rt_window`] — `(scan_index, gaussian_weight)` pairs,
/// shared across all charge hypotheses of the same seed.
fn score_hypothesis(
    engine: &PeakIndexingEngine,
    seed: &IndexedMassSpectralPeak,
    charge: i32,
    params: &TraceKernelParameters,
    ppm: &PpmTolerance,
    claimed: &HashSet<PeakKey>,
    window: &[(i32, f64)],
) -> HypothesisScore {
    let seed_mz = seed.m() as f64;
    let seed_mass = mz_to_mass(seed_mz, charge);
    let weights = comb_weights(seed_mass, params);
    // An empty envelope (e.g. a degenerate weight model) has no comb to lay — score it as a miss
    // rather than indexing into an empty vector.
    if weights.is_empty() {
        return HypothesisScore {
            response: 0.0,
            charge,
            mono_mz: seed_mz,
            peaks: Vec::new(),
            num_isotopes_observed: 0,
        };
    }
    let i_star = most_abundant_index(&weights);
    let spacing = C13_MINUS_C12 / charge as f64;
    let mono_mz = seed_mz - (i_star as f64) * spacing;

    let mut peaks: Vec<IndexedMassSpectralPeak> = Vec::new();
    let mut observed_isotopes: HashSet<usize> = HashSet::new();
    // Peaks already used *within this hypothesis*. For higher charges the comb spacing (1.0033/z) is
    // small, so two adjacent isotope slots can resolve to the same physical peak; without this a peak
    // would be double-counted in the response and intensity and would inflate the isotope count.
    let mut used: HashSet<PeakKey> = HashSet::new();
    // Every comb (isotope k, scan s) slot as `(template = wₖ·gₛ, observed intensity)`. A missing or
    // already-claimed/used peak keeps its template weight but contributes zero observed intensity, so
    // the normalised score can penalise a predicted-but-absent tooth. `RawSum` only reads `t · I`.
    let mut slots: Vec<(f64, f64)> = Vec::with_capacity(window.len() * weights.len());

    for &(s, g) in window {
        for (k, &wk) in weights.iter().enumerate() {
            let template = wk * g;
            let expected_mz = mono_mz + (k as f64) * spacing;
            let observed = if let Some(peak) = engine.get_indexed_peak(expected_mz, s, ppm) {
                let key = peak.key();
                if !claimed.contains(&key) && used.insert(key) {
                    peaks.push(*peak);
                    observed_isotopes.insert(k);
                    peak.intensity as f64
                } else {
                    // Claimed by another feature, or already consumed by another slot of this
                    // hypothesis — absent for this slot.
                    0.0
                }
            } else {
                0.0
            };
            slots.push((template, observed));
        }
    }

    HypothesisScore {
        response: hypothesis_response(
            &slots,
            params.score_model,
            params.noise_floor,
            seed.intensity as f64,
            params.score_use_seed_amplitude,
            params.score_cosine,
        ),
        charge,
        mono_mz,
        peaks,
        num_isotopes_observed: observed_isotopes.len(),
    }
}

/// Reduces a hypothesis's comb slots `(template = wₖ·gₛ, observed_intensity)` to the scalar response
/// used for cross-z NMS, under the chosen [`ScoreModel`].
///
/// - [`ScoreModel::RawSum`]: `Σ (t · I)` — identical to the pre-Change-B accumulation.
/// - [`ScoreModel::NormalizedNoiseFloor`]: estimate the apex amplitude `A` (least-squares
///   `Σ(t·I)/Σ(t²)` over observed slots, or `seed_intensity` when `use_seed_amplitude`), form the
///   expected-observable support `S = { slots : A·t ≥ η }`, and return the template-normalised
///   correlation `Σ_S(t·I)/sqrt(Σ_S t²)` (or, when `cosine`, the bounded cosine
///   `Σ_S(t·I)/(sqrt(Σ_S t²)·sqrt(Σ_S I²))`). Slots predicted below `η` are dropped from both sums;
///   predicted-and-present teeth reward, predicted-and-absent teeth (in `S`) penalise.
fn hypothesis_response(
    slots: &[(f64, f64)],
    model: ScoreModel,
    noise_floor: f64,
    seed_intensity: f64,
    use_seed_amplitude: bool,
    cosine: bool,
) -> f64 {
    match model {
        ScoreModel::RawSum => slots.iter().map(|(t, i)| t * i).sum(),
        ScoreModel::NormalizedNoiseFloor => {
            let a = if use_seed_amplitude {
                seed_intensity
            } else {
                let mut num_a = 0.0;
                let mut den_a = 0.0;
                for &(t, i) in slots {
                    if i > 0.0 {
                        num_a += t * i;
                        den_a += t * t;
                    }
                }
                if den_a <= 0.0 {
                    return 0.0;
                }
                num_a / den_a
            };
            let mut num = 0.0;
            let mut den_t = 0.0;
            let mut den_i = 0.0;
            for &(t, i) in slots {
                if a * t >= noise_floor {
                    num += t * i;
                    den_t += t * t;
                    den_i += i * i;
                }
            }
            if den_t <= 0.0 {
                return 0.0;
            }
            let template_norm = num / den_t.sqrt();
            if cosine {
                if den_i <= 0.0 {
                    return 0.0;
                }
                template_norm / den_i.sqrt()
            } else {
                template_norm
            }
        }
    }
}

/// Estimates the run-level MS1 noise floor `η` as a low percentile of the positive peak intensities —
/// a simple global baseline for [`ScoreModel::NormalizedNoiseFloor`]. `percentile` is in `[0, 100]`
/// (e.g. `5.0` for the 5th percentile). Returns `0.0` for an empty index (which disables the
/// noise-floor truncation, i.e. a full-template norm).
pub fn estimate_noise_floor(engine: &PeakIndexingEngine, percentile: f64) -> f64 {
    let mut intensities: Vec<f64> = engine
        .all_peaks()
        .iter()
        .map(|p| p.intensity as f64)
        .filter(|&i| i > 0.0)
        .collect();
    if intensities.is_empty() {
        return 0.0;
    }
    intensities.sort_by(|a, b| a.total_cmp(b));
    let p = percentile.clamp(0.0, 100.0) / 100.0;
    let idx = (((intensities.len() - 1) as f64) * p).round() as usize;
    intensities[idx]
}

/// Traces the accepted hypothesis's **true elution extent** and returns the peaks the feature will
/// claim (Change A). This decouples the *claim* from the narrow ~2σ *scoring* window: a real peak
/// wider than 2σ is claimed whole, so its smaller adjacent seeds are already claimed and never fire —
/// fragmentation never forms, and there is nothing to merge downstream.
///
/// Two steps:
/// 1. **Extent.** Follow the most-abundant tooth (the seed's m/z — highest SNR, most reliable
///    boundary) outward in RT with [`PeakIndexingEngine::get_xic_by_scan_index`], stopping on
///    `trace_missed_scans_allowed` consecutive misses or the `trace_max_half_width_minutes` guard.
///    Its peaks' scan indices give the extent `[s_lo, s_hi]`. Peaks already in `claimed` count as
///    misses (a taller neighbour claimed them first), which is what splits co-eluting same-m/z peaks
///    greedily instead of merging them.
/// 2. **Gather.** Collect every comb tooth's peak at each scan in `[s_lo, s_hi]`, excluding anything
///    already `claimed`. The union (deduped) is the feature's peak set — the single set that backs its
///    RT bounds, summed intensity, coverage contribution, and the NMS claim mask alike.
fn trace_claim_extent(
    engine: &PeakIndexingEngine,
    seed: &IndexedMassSpectralPeak,
    hyp: &HypothesisScore,
    params: &TraceKernelParameters,
    ppm: &PpmTolerance,
    claimed: &HashSet<PeakKey>,
) -> Vec<IndexedMassSpectralPeak> {
    let charge = hyp.charge;
    let seed_mz = seed.m() as f64;
    let weights = comb_weights(mz_to_mass(seed_mz, charge), params);
    if weights.is_empty() {
        // Degenerate comb — nothing to trace; claim the scored peaks (minus any already claimed).
        return hyp
            .peaks
            .iter()
            .filter(|p| !claimed.contains(&p.key()))
            .copied()
            .collect();
    }
    let spacing = C13_MINUS_C12 / charge as f64;
    let mono_mz = hyp.mono_mz;
    let apex_scan = seed.zero_based_scan_index;

    // 1) Extent: follow the most-abundant tooth to fix the scan span.
    let trace = engine.get_xic_by_scan_index(
        seed_mz,
        apex_scan,
        ppm,
        params.trace_missed_scans_allowed,
        params.trace_max_half_width_minutes,
        Some(claimed),
    );
    let mut s_lo = apex_scan;
    let mut s_hi = apex_scan;
    for p in &trace {
        s_lo = s_lo.min(p.zero_based_scan_index);
        s_hi = s_hi.max(p.zero_based_scan_index);
    }

    // 2) Gather every comb tooth's peaks across the extent, skipping already-claimed peaks and
    //    de-duplicating (adjacent teeth of a high charge can resolve to the same physical peak).
    let mut seen: HashSet<PeakKey> = HashSet::new();
    let mut peaks: Vec<IndexedMassSpectralPeak> = Vec::new();
    for k in 0..weights.len() {
        let tooth_mz = mono_mz + (k as f64) * spacing;
        for s in s_lo..=s_hi {
            if let Some(p) = engine.get_indexed_peak(tooth_mz, s, ppm) {
                let key = p.key();
                if claimed.contains(&key) || !seen.insert(key) {
                    continue;
                }
                peaks.push(*p);
            }
        }
    }

    // Defensive: never let an accepted feature end up with an empty peak set (the seed alone should
    // always survive), which would break `build_feature`'s apex/extent derivation.
    if peaks.is_empty() {
        for p in &hyp.peaks {
            let key = p.key();
            if !claimed.contains(&key) && seen.insert(key) {
                peaks.push(*p);
            }
        }
    }
    peaks
}

/// Runs the trace-kernel detector over an indexed run.
///
/// Seeds are every indexed peak taken tallest-first (the busiest RT regions hold the most signal, so
/// this maximises explained-signal-per-feature — the coverage objective). For each unclaimed seed,
/// charge hypotheses `min..=max` are scored and the best-response one wins (cross-z NMS). A
/// hypothesis is accepted when it observes at least `min_isotopes_observed` teeth and has positive
/// response; its peaks are then claimed so overlapping harmonics cannot re-fire. Returns the detected
/// features in acceptance order (tallest-seed first).
pub fn detect_features(
    engine: &PeakIndexingEngine,
    params: &TraceKernelParameters,
) -> Vec<DetectedFeature> {
    let ppm = PpmTolerance::new(params.ppm_tolerance);

    // Seeds: all peaks, tallest first. Stable ordering (intensity desc) mirrors get_all_xics.
    let mut seeds = engine.all_peaks();
    seeds.sort_by(|a, b| b.intensity.total_cmp(&a.intensity));

    // Coverage bookkeeping: denominator = Σ all peak intensities; stop once the claimed fraction
    // reaches `coverage_target`. The sum is only needed when the cap is actually engaged
    // (`coverage_target < 1.0`), so skip the O(peaks) pass for the common "detect everything" case.
    let coverage_stop = if params.coverage_target < 1.0 {
        let total_intensity: f64 = seeds.iter().map(|p| p.intensity as f64).sum();
        if total_intensity > 0.0 {
            params.coverage_target * total_intensity
        } else {
            f64::INFINITY
        }
    } else {
        f64::INFINITY
    };
    let mut explained_intensity = 0.0;

    let mut claimed: HashSet<PeakKey> = HashSet::new();
    let mut features: Vec<DetectedFeature> = Vec::new();

    for seed in &seeds {
        // Seeds are intensity-descending, so once we fall below the floor every remaining seed is
        // too — stop rather than continue.
        if (seed.intensity as f64) < params.min_seed_intensity {
            break;
        }
        if claimed.contains(&seed.key()) {
            continue;
        }

        // The RT window (scan indices + Gaussian weights) is charge-independent — compute it once
        // per seed and share it across every charge hypothesis.
        let window = seed_rt_window(engine, seed, params);

        // Score every charge hypothesis; keep the highest response (cross-z non-max suppression).
        let mut best: Option<HypothesisScore> = None;
        for z in params.min_charge..=params.max_charge {
            if z == 0 {
                continue;
            }
            let score = score_hypothesis(engine, seed, z, params, &ppm, &claimed, &window);
            let better = match &best {
                None => true,
                Some(b) => score.response > b.response,
            };
            if better {
                best = Some(score);
            }
        }

        let best = match best {
            Some(b) => b,
            None => continue,
        };

        if best.num_isotopes_observed < params.min_isotopes_observed || best.response <= 0.0 {
            // Not a feature; retire this seed so we do not reconsider it.
            claimed.insert(seed.key());
            continue;
        }

        // Claim the feature's TRUE traced extent (not just the narrow scored window), so the whole
        // elution is claimed at once and its smaller adjacent seeds cannot re-fire as fragments.
        let traced = trace_claim_extent(engine, seed, &best, params, &ppm, &claimed);

        // Chromatographic-persistence gate: a real elution spans several scans; a feature whose
        // traced extent covers fewer than `min_feature_scans` distinct scans is a single-scan noise
        // doublet, not a peak. Retire the seed (as with the isotope-count gate) without emitting it.
        if params.min_feature_scans > 1 {
            let distinct_scans: HashSet<i32> =
                traced.iter().map(|p| p.zero_based_scan_index).collect();
            if distinct_scans.len() < params.min_feature_scans {
                claimed.insert(seed.key());
                continue;
            }
        }

        for p in &traced {
            claimed.insert(p.key());
        }
        let feature = build_feature(best, traced);
        explained_intensity += feature.summed_intensity;
        features.push(feature);

        // Stop once we have explained the target fraction of the total MS1 signal.
        if explained_intensity >= coverage_stop {
            break;
        }
    }

    features
}

/// Assembles a [`DetectedFeature`] from an accepted hypothesis and its **traced** peak set (Change A).
///
/// The apex, RT bounds, and summed intensity are derived from `peaks` — the true traced extent from
/// [`trace_claim_extent`], the same set that was written into the claim mask — so the reported RT
/// bounds, the coverage/quant intensity, and the NMS claim are all backed by one peak set. The `score`
/// stays the hypothesis's narrow matched-filter `response` (a shape-fit, deliberately not an
/// extent-sum), and `num_isotopes_observed` stays the scored-window tooth count.
fn build_feature(hyp: HypothesisScore, peaks: Vec<IndexedMassSpectralPeak>) -> DetectedFeature {
    let apex = peaks
        .iter()
        .max_by(|a, b| a.intensity.total_cmp(&b.intensity))
        .expect("accepted feature has at least one traced peak");
    let apex_scan_index = apex.zero_based_scan_index;
    let apex_rt = apex.retention_time as f64;
    let start_rt = peaks
        .iter()
        .map(|p| p.retention_time as f64)
        .fold(f64::INFINITY, f64::min);
    let end_rt = peaks
        .iter()
        .map(|p| p.retention_time as f64)
        .fold(f64::NEG_INFINITY, f64::max);
    let summed_intensity = peaks.iter().map(|p| p.intensity as f64).sum();

    DetectedFeature {
        monoisotopic_mass: mz_to_mass(hyp.mono_mz, hyp.charge),
        charge: hyp.charge,
        mono_mz: hyp.mono_mz,
        apex_scan_index,
        apex_rt,
        start_rt,
        end_rt,
        summed_intensity,
        score: hyp.response,
        num_isotopes_observed: hyp.num_isotopes_observed,
        peaks,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::deconvolution::averagine_intensities_from_mono;
    use crate::isotopic_envelope::mass_to_mz_f64;
    use crate::peak_indexing::Scan;

    fn approx(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "expected {b}, got {a} (tol {tol})");
    }

    #[test]
    fn poisson_weights_light_mass_mono_is_tallest() {
        // M = 1000 → λ = 0.48, mode at k=0 (mono tallest), strictly descending.
        let w = poisson_comb_weights(1000.0, 1e-3, 12);
        assert_eq!(most_abundant_index(&w), 0);
        assert!(w[0] >= w[1] && w[1] >= w[2]);
        approx(w[0], 1.0, 1e-12); // normalised so the max is 1
    }

    #[test]
    fn poisson_weights_heavy_mass_mode_shifts_up() {
        // M = 5000 → λ = 2.4, mode at k=2 — the seed (tallest peak) is NOT the monoisotopic peak,
        // which is exactly why the comb must look below the seed.
        let w = poisson_comb_weights(5000.0, 1e-3, 20);
        assert_eq!(most_abundant_index(&w), 2);
    }

    #[test]
    fn median_scan_spacing_is_robust() {
        let scans = synthetic_envelope_scans().0;
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        // RTs are 0.1 min apart.
        approx(median_ms1_scan_spacing_minutes(engine.scan_info()), 0.1, 1e-9);
    }

    /// Builds a clean charge-2 isotope envelope (monoisotopic neutral mass 1000) eluting across
    /// 9 scans with a Gaussian RT profile (apex at scan 4). Each isotope tooth carries its Poisson
    /// weight × the apex intensity × the RT Gaussian. Returns the scans and the true mono m/z.
    fn synthetic_envelope_scans() -> (Vec<Scan>, f64) {
        let mono_mass = 1000.0;
        let charge = 2;
        let mono_mz = mass_to_mz_f64(mono_mass, charge); // ~501.007
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
        (scans, mono_mz)
    }

    #[test]
    fn detects_charge_two_envelope_with_correct_mass() {
        let (scans, mono_mz) = synthetic_envelope_scans();
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let params = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.15,
            half_window_scans: 4,
            ..TraceKernelParameters::default()
        };
        let features = detect_features(&engine, &params);
        assert!(!features.is_empty(), "should detect at least one feature");

        // The top feature (tallest seed) must be the charge-2 envelope with mono mass ~1000.
        let top = &features[0];
        assert_eq!(top.charge, 2, "charge-2 envelope must beat other z hypotheses");
        approx(top.monoisotopic_mass, 1000.0, 0.01);
        approx(top.mono_mz, mono_mz, 1e-4);
        assert!(top.num_isotopes_observed >= 2);
        assert_eq!(top.apex_scan_index, 4, "apex is the max-intensity scan");
        assert!(top.score > 0.0);
    }

    #[test]
    fn detects_charge_two_envelope_with_averagine_weights() {
        // The averagine comb-weight model should detect the same synthetic z=2 envelope with the
        // correct mass/charge, exercising the CombWeightModel::Averagine path end to end.
        let (scans, mono_mz) = synthetic_envelope_scans();
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let params = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.15,
            half_window_scans: 4,
            weight_model: CombWeightModel::Averagine,
            ..TraceKernelParameters::default()
        };
        let features = detect_features(&engine, &params);
        assert!(!features.is_empty(), "averagine model should detect the envelope");
        let top = &features[0];
        assert_eq!(top.charge, 2);
        approx(top.monoisotopic_mass, 1000.0, 0.01);
        approx(top.mono_mz, mono_mz, 1e-4);
    }

    #[test]
    fn charge_two_beats_charge_one_and_four_on_response() {
        // Directly confirm the cross-z ranking on the same seed: score z=1,2,4 for the apex mono
        // peak and check z=2 wins.
        let (scans, _) = synthetic_envelope_scans();
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let params = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.15,
            half_window_scans: 4,
            ..TraceKernelParameters::default()
        };
        let ppm = PpmTolerance::new(params.ppm_tolerance);
        let claimed = HashSet::new();

        // Seed = the tallest peak (apex scan mono).
        let mut seeds = engine.all_peaks();
        seeds.sort_by(|a, b| b.intensity.total_cmp(&a.intensity));
        let seed = seeds[0];

        let window = seed_rt_window(&engine, &seed, &params);
        let s1 = score_hypothesis(&engine, &seed, 1, &params, &ppm, &claimed, &window).response;
        let s2 = score_hypothesis(&engine, &seed, 2, &params, &ppm, &claimed, &window).response;
        let s4 = score_hypothesis(&engine, &seed, 4, &params, &ppm, &claimed, &window).response;
        assert!(s2 > s1, "z=2 ({s2}) should beat z=1 ({s1})");
        assert!(s2 > s4, "z=2 ({s2}) should beat z=4 ({s4})");
    }

    #[test]
    fn raw_sum_response_is_plain_inner_product() {
        // RawSum must equal Σ(template · observed); a missing tooth (observed 0) contributes nothing.
        let slots = [(1.0, 10.0), (0.5, 4.0), (0.25, 0.0)];
        approx(
            hypothesis_response(&slots, ScoreModel::RawSum, 0.0, 0.0, false, false),
            1.0 * 10.0 + 0.5 * 4.0,
            1e-9,
        );
    }

    #[test]
    fn normalized_response_penalizes_missing_predicted_tooth() {
        // Same tall tooth, but one hypothesis is missing a second tooth the template predicts well
        // above the floor. It must rank below the complete one (the absent tooth stays in the norm).
        let eta = 1.0;
        let complete = hypothesis_response(&[(1.0, 100.0), (0.6, 60.0)], ScoreModel::NormalizedNoiseFloor, eta, 0.0, false, false);
        let missing = hypothesis_response(&[(1.0, 100.0), (0.6, 0.0)], ScoreModel::NormalizedNoiseFloor, eta, 0.0, false, false);
        assert!(
            complete > missing,
            "complete ({complete}) should beat a missing predicted tooth ({missing})"
        );
    }

    #[test]
    fn normalized_noise_floor_excludes_below_floor_teeth() {
        // A faint tooth the model predicts BELOW the floor must not penalise: a hypothesis missing
        // only that below-floor tooth scores the same as one where the slot never existed.
        // A = 100 (from the tall tooth); a tooth with template 0.005 → A·t = 0.5 < η = 1.0 → excluded.
        let eta = 1.0;
        let with_faint = hypothesis_response(&[(1.0, 100.0), (0.005, 0.0)], ScoreModel::NormalizedNoiseFloor, eta, 0.0, false, false);
        let without = hypothesis_response(&[(1.0, 100.0)], ScoreModel::NormalizedNoiseFloor, eta, 0.0, false, false);
        approx(with_faint, without, 1e-9);
    }

    #[test]
    fn normalized_score_still_selects_charge_two() {
        // The normalised model must not break basic charge selection on a clean z=2 envelope.
        let (scans, _) = synthetic_envelope_scans();
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let params = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.15,
            half_window_scans: 4,
            score_model: ScoreModel::NormalizedNoiseFloor,
            noise_floor: 0.0,
            ..TraceKernelParameters::default()
        };
        let features = detect_features(&engine, &params);
        assert!(!features.is_empty(), "normalised model should still detect the envelope");
        assert_eq!(features[0].charge, 2, "normalised score must still pick z=2");
    }

    #[test]
    fn cosine_response_is_one_for_a_perfect_fit() {
        // Observed I = A·t exactly for every slot → the cosine variant returns 1.0 (bounded shape-fit).
        let a = 50.0;
        let slots = [(1.0, a * 1.0), (0.5, a * 0.5), (0.25, a * 0.25)];
        let r = hypothesis_response(&slots, ScoreModel::NormalizedNoiseFloor, 0.0, 0.0, false, true);
        approx(r, 1.0, 1e-9);
    }

    #[test]
    fn seed_amplitude_matches_ls_when_amplitudes_agree() {
        // The LS amplitude of these observed slots is (100+36)/(1+0.36) = 100; passing seed_intensity
        // = 100 to the seed-amplitude path must yield the same support S and thus the same score.
        let slots = [(1.0, 100.0), (0.6, 60.0), (0.01, 0.0)];
        let ls = hypothesis_response(&slots, ScoreModel::NormalizedNoiseFloor, 1.0, 0.0, false, false);
        let seed = hypothesis_response(&slots, ScoreModel::NormalizedNoiseFloor, 1.0, 100.0, true, false);
        approx(ls, seed, 1e-9);
    }

    #[test]
    fn claiming_prevents_double_detection() {
        // A single clean envelope should yield exactly one feature — its peaks get claimed, so no
        // second feature is assembled from the same signal.
        let (scans, _) = synthetic_envelope_scans();
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let params = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.15,
            half_window_scans: 4,
            ..TraceKernelParameters::default()
        };
        let features = detect_features(&engine, &params);
        assert_eq!(features.len(), 1, "one envelope → one feature");
    }

    #[test]
    fn with_rt_from_scans_derives_sigma_and_window() {
        let (scans, _) = synthetic_envelope_scans();
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        // 36 s FWHM → σ = 0.6 min / 2.3548 ≈ 0.2548 min; spacing 0.1 min → half-window ≈ round(5.1) = 5.
        let params = TraceKernelParameters::default()
            .with_rt_from_scans(engine.scan_info(), 36.0);
        approx(params.rt_sigma_minutes, (36.0 / 60.0) / FWHM_TO_SIGMA, 1e-9);
        assert_eq!(params.half_window_scans, 5);
    }

    #[test]
    fn default_weight_model_is_averagine() {
        // The detector defaults to the table-driven averagine envelope (Poisson is opt-in now).
        assert_eq!(
            TraceKernelParameters::default().weight_model,
            CombWeightModel::Averagine
        );
    }

    /// Builds a synthetic isotope envelope whose per-tooth intensities come from the **real averagine
    /// distribution** (not Poisson), for a peptide of monoisotopic neutral mass `mono_mass` at
    /// `charge`, eluting across `n_scans` with a Gaussian RT profile (apex at the middle scan).
    /// Returns the scans and the true mono m/z. This exercises the averagine comb model on
    /// averagine-shaped data (the previous averagine test reused Poisson-shaped input).
    fn averagine_envelope_scans(mono_mass: f64, charge: i32, n_scans: i32) -> (Vec<Scan>, f64) {
        let mono_mz = mass_to_mz_f64(mono_mass, charge);
        let spacing = C13_MINUS_C12 / charge as f64;
        // Mono-keyed averagine intensities: index 0 is the monoisotope, then +1 ¹³C, +2, …
        let weights = averagine_intensities_from_mono(mono_mass, 1e-4, 20);
        let apex_intensity = 1.0e7;
        let rt_sigma = 0.15;
        let apex_scan = n_scans / 2;
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
        (scans, mono_mz)
    }

    #[test]
    fn averagine_detects_light_envelope_default_model() {
        // A light peptide with the default (now averagine) params — mono is the tallest tooth here,
        // so this is the easy case that must keep working after the default flip.
        let (scans, mono_mz) = averagine_envelope_scans(1200.0, 2, 9);
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let params = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.15,
            half_window_scans: 4,
            ..TraceKernelParameters::default()
        };
        assert_eq!(params.weight_model, CombWeightModel::Averagine);
        let features = detect_features(&engine, &params);
        assert!(!features.is_empty(), "should detect the light averagine envelope");
        let top = &features[0];
        assert_eq!(top.charge, 2);
        approx(top.monoisotopic_mass, 1200.0, 0.02);
        approx(top.mono_mz, mono_mz, 1e-3);
        assert!(top.num_isotopes_observed >= 2);
    }

    #[test]
    fn averagine_places_mono_below_seed_in_mode_shift_regime() {
        // The regime that motivates averagine: heavy peptides whose envelope mode sits *above* the
        // monoisotope. The seed (tallest peak) is then NOT the mono, so the comb must look below it
        // by exactly `i*` ¹³C units. If averagine misplaces `i*`, the recovered mono mass is off by
        // ~1 Da/charge; a tight tolerance here is what proves the placement is correct.
        for &(mono_mass, charge) in &[(2400.0, 3), (4000.0, 4), (5200.0, 4)] {
            let (scans, mono_mz) = averagine_envelope_scans(mono_mass, charge, 9);
            let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
            let params = TraceKernelParameters {
                ppm_tolerance: 5.0,
                rt_sigma_minutes: 0.15,
                half_window_scans: 4,
                ..TraceKernelParameters::default()
            };
            let features = detect_features(&engine, &params);
            assert!(
                !features.is_empty(),
                "should detect the averagine envelope at mono {mono_mass}, z{charge}"
            );
            let top = &features[0];
            assert_eq!(top.charge, charge, "charge for mono {mono_mass}");
            // Recovered mono within ~1/4 of a ¹³C unit at this charge — far tighter than an
            // off-by-one error (which would be ~1 Da) would allow.
            approx(top.monoisotopic_mass, mono_mass, 0.05);
            approx(top.mono_mz, mono_mz, 1e-3);
        }
    }

    #[test]
    fn averagine_seed_is_above_mono_for_heavy_mass() {
        // Sanity-check the fixture itself: for a heavy peptide the most-intense (seed) tooth really
        // is above the monoisotope, so the mode-shift test above is exercising the intended path.
        let w = averagine_intensities_from_mono(4000.0, 1e-4, 20);
        let mode = most_abundant_index(&w);
        assert!(mode >= 1, "heavy averagine mode should sit above the mono, got {mode}");
    }

    /// Builds a clean averagine-shaped z=`charge` envelope eluting across `n_scans` (0.1-min spacing
    /// from `rt0`), with the RT Gaussian centred at local scan `apex_scan` and width `rt_sigma`. Every
    /// scan in the span carries the full comb, so the elution is contiguous (no missed scans) and the
    /// claim trace can cover the whole extent. `first_scan_number` sets the (cosmetic) one-based
    /// numbering; the zero-based scan index comes from position in the concatenated scan array.
    fn elution_scans(
        mono_mass: f64,
        charge: i32,
        n_scans: i32,
        apex_scan: i32,
        rt_sigma: f64,
        rt0: f64,
        first_scan_number: i32,
    ) -> Vec<Scan> {
        let mono_mz = mass_to_mz_f64(mono_mass, charge);
        let spacing = C13_MINUS_C12 / charge as f64;
        let weights = averagine_intensities_from_mono(mono_mass, 1e-4, 12);
        let apex_intensity = 1.0e7;
        let mut scans = Vec::new();
        for s in 0..n_scans {
            let rt = rt0 + s as f64 * 0.1;
            let g = gaussian(rt - (rt0 + apex_scan as f64 * 0.1), rt_sigma);
            let mut mz = Vec::new();
            let mut intensity = Vec::new();
            for (k, &wk) in weights.iter().enumerate() {
                mz.push(mono_mz + k as f64 * spacing);
                intensity.push(apex_intensity * wk * g);
            }
            scans.push(Scan {
                mz,
                intensity,
                one_based_scan_number: first_scan_number + s,
                retention_time: rt,
                msn_order: 1,
            });
        }
        scans
    }

    #[test]
    fn trace_claims_full_elution_not_just_scoring_window() {
        // A z=2 elution ~1.5 min wide (16 scans). The SCORING window is deliberately narrow
        // (±~0.12 min ≈ 1 scan each side) — the regime that used to split a wide peak into many
        // seeds. Trace-following (Change A) must claim the WHOLE elution → exactly one feature that
        // spans all 16 scans, with its RT bounds and peak set covering the true extent.
        let n = 16;
        let apex = 8;
        let scans = elution_scans(1200.0, 2, n, apex, 0.35, 10.0, 1);
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let params = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.1,
            rt_half_window_minutes: 0.12, // narrow scoring window on purpose
            trace_missed_scans_allowed: 1,
            trace_max_half_width_minutes: 1.5, // wide enough to follow the whole elution
            ..TraceKernelParameters::default()
        };
        let features = detect_features(&engine, &params);
        assert_eq!(
            features.len(),
            1,
            "a wide elution must be ONE feature, not fragmented into many"
        );
        let f = &features[0];
        assert_eq!(f.charge, 2);
        approx(f.start_rt, 10.0, 1e-3);
        approx(f.end_rt, 10.0 + (n - 1) as f64 * 0.1, 1e-3);
        // The claim covers every scan of the elution across the comb teeth.
        let distinct_scans: HashSet<i32> =
            f.peaks.iter().map(|p| p.zero_based_scan_index).collect();
        assert_eq!(
            distinct_scans.len(),
            n as usize,
            "claim should span every scan of the elution, got {} of {n}",
            distinct_scans.len()
        );
    }

    #[test]
    fn trace_stops_at_gap_between_co_eluting_same_mz_peaks() {
        // Two z=2 elutions at the SAME m/z separated by a 3-scan empty gap. With
        // trace_missed_scans_allowed = 1 the claim trace stops in the gap rather than merging the
        // two, so greedy tallest-first detection yields TWO features, each on its own side.
        let mut scans = elution_scans(1200.0, 2, 5, 2, 0.15, 10.0, 1); // A: idx 0..4, apex idx 2
        for i in 0..3 {
            // Empty gap scans — present in scan_info, so the scan-index gap is real.
            scans.push(Scan {
                mz: vec![],
                intensity: vec![],
                one_based_scan_number: 6 + i,
                retention_time: 10.5 + i as f64 * 0.1,
                msn_order: 1,
            });
        }
        // B: same m/z, later RT, scaled down so A (taller) seeds and claims first.
        let mut b = elution_scans(1200.0, 2, 5, 2, 0.15, 10.8, 9); // B: idx 8..12, apex idx 10
        for s in &mut b {
            for y in &mut s.intensity {
                *y *= 0.4;
            }
        }
        scans.extend(b);

        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let params = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.15,
            rt_half_window_minutes: 0.3,
            trace_missed_scans_allowed: 1,
            trace_max_half_width_minutes: 1.5,
            ..TraceKernelParameters::default()
        };
        let features = detect_features(&engine, &params);
        assert_eq!(
            features.len(),
            2,
            "a valley gap must keep the two same-m/z elutions separate"
        );
        for f in &features {
            assert!(
                f.end_rt - f.start_rt < 0.6,
                "feature spanning [{:.3}, {:.3}] merged across the gap",
                f.start_rt,
                f.end_rt
            );
        }
    }

    #[test]
    fn min_feature_scans_rejects_single_scan_noise() {
        // A 2-tooth envelope present in exactly ONE scan (neighbours empty), i.e. a single-scan
        // noise doublet. It clears the isotope-count gate but its traced extent is one scan.
        let mono_mass = 1200.0;
        let charge = 2;
        let mono_mz = mass_to_mz_f64(mono_mass, charge);
        let spacing = C13_MINUS_C12 / charge as f64;
        let weights = averagine_intensities_from_mono(mono_mass, 1e-4, 6);
        let mut scans = Vec::new();
        for s in 0..5 {
            let (mz, intensity) = if s == 2 {
                let mut mz = Vec::new();
                let mut inten = Vec::new();
                for (k, &wk) in weights.iter().enumerate() {
                    mz.push(mono_mz + k as f64 * spacing);
                    inten.push(1.0e7 * wk);
                }
                (mz, inten)
            } else {
                (Vec::new(), Vec::new())
            };
            scans.push(Scan {
                mz,
                intensity,
                one_based_scan_number: s + 1,
                retention_time: 10.0 + s as f64 * 0.1,
                msn_order: 1,
            });
        }
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let base = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.1,
            rt_half_window_minutes: 0.3,
            ..TraceKernelParameters::default()
        };

        // Gate off (=1): the single-scan doublet is detected.
        let allow = TraceKernelParameters { min_feature_scans: 1, ..base };
        assert_eq!(
            detect_features(&engine, &allow).len(),
            1,
            "with the persistence gate off, the single-scan doublet is a feature"
        );

        // Gate on (=2, the default): it is rejected as non-chromatographic.
        let gate = TraceKernelParameters { min_feature_scans: 2, ..base };
        assert!(
            detect_features(&engine, &gate).is_empty(),
            "min_feature_scans=2 must reject a feature whose traced extent is one scan"
        );
    }

    #[test]
    fn with_rt_from_index_measures_sigma_from_data() {
        // Several clean elutions with a KNOWN RT width (σ_data = 0.15 min → FWHM ≈ 0.353 min ≈ 21 s).
        // with_rt_from_index must recover σ ≈ 0.15 from XIC half-max — NOT the (different) fallback.
        let mut scans = elution_scans(1500.0, 2, 15, 7, 0.15, 10.0, 1);
        scans.extend(elution_scans(2000.0, 3, 15, 7, 0.15, 20.0, 100));
        scans.extend(elution_scans(2500.0, 3, 15, 7, 0.15, 30.0, 200));
        scans.extend(elution_scans(3000.0, 3, 15, 7, 0.15, 40.0, 300));
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");

        let ppm = PpmTolerance::new(5.0);
        let measured = estimate_fwhm_seconds(&engine, &ppm).expect("enough clean XICs to measure");
        // ~21 s, well within the sane band.
        approx(measured, 0.15 * FWHM_TO_SIGMA * 60.0, 3.0);

        let params = TraceKernelParameters {
            ppm_tolerance: 5.0,
            ..TraceKernelParameters::default()
        }
        // Fallback 40 s (σ ≈ 0.283) is deliberately far from the true 0.15 so a wrong fallback shows.
        .with_rt_from_index(&engine, 40.0);
        approx(params.rt_sigma_minutes, 0.15, 0.03);
        assert!(params.rt_half_window_minutes > 0.0);
    }

    #[test]
    fn with_rt_from_index_falls_back_when_too_few_clean_xics() {
        // A 2-scan blip has no measurable clean XIC (apex at an edge, < 3 points) → the estimator
        // returns None → params take the fallback FWHM.
        let scans = elution_scans(1000.0, 2, 2, 0, 0.15, 10.0, 1);
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let ppm = PpmTolerance::new(10.0);
        assert!(
            estimate_fwhm_seconds(&engine, &ppm).is_none(),
            "a 2-scan run should not yield a measurable FWHM"
        );

        let fallback = 30.0;
        let params = TraceKernelParameters::default().with_rt_from_index(&engine, fallback);
        approx(
            params.rt_sigma_minutes,
            (fallback / 60.0) / FWHM_TO_SIGMA,
            1e-9,
        );
    }
}
