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
use std::time::{Duration, Instant};

use rayon::prelude::*;

use crate::isotopic_envelope::{C13_MINUS_C12, PROTON_MASS};
use crate::peak_indexing::{
    IndexedMassSpectralPeak, PeakIndexingEngine, PeakKey, PeakSource, ScanInfo,
};
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
    /// **Opt-in, default `false`.** Auto-stop the seed walk at the *knee* of the coverage curve —
    /// where the marginal %ΣTIC claimed per seed collapses. Because seeds are visited tallest-first
    /// the marginal claim is monotonically non-increasing (the curve is concave by construction), so
    /// a rolling 2-point slope suffices. A pure speed/recall knob, not a correctness fix: it trades
    /// the low-abundance tail (the bulk of seeds/wall-clock) for early exit. `false` => never stops
    /// early (byte-identical to the un-instrumented detector).
    pub knee_stop_enabled: bool,
    /// Width (in seeds visited) of the rolling window over which the knee slope is measured.
    pub knee_window_seeds: usize,
    /// Stop when the current window's slope falls below this fraction of the FIRST window's slope.
    pub knee_slope_frac: f64,
    /// Absolute floor on the rolling slope (fraction of ΣTIC per seed) below which the tail is flat
    /// enough to stop regardless of the relative test.
    pub knee_abs_eps: f64,
    /// **Opt-in, default `false` — the shipped detector runs UNCAPPED.** Auto-stop the seed walk once the
    /// rolling fraction of *considered* (scored) seeds that get rejected exceeds [`Self::reject_stop_frac`].
    /// Unlike the TIC knee this keys on the detector's own hit-rate — a self-referential,
    /// machine/sample-portable signal — and fires where the tail turns mostly to noise
    /// (persistence-rejected single-scan spikes). On the parallel paths it is applied **per tile** (see
    /// [`tile_reject_cfg`] / [`detect_bin`]), which is what lets a capped run stay parallel. A speed/recall
    /// knob, not a correctness fix.
    ///
    /// **Why it defaults off (2026-07-08 A/B, `COVERAGE_TARGET=1.0`, 2-D tiling path, PSM-recall vs the
    /// uncapped baseline on 10-min CA / 65-min glyco / 2-hr IonStar).** At the best setting found
    /// (`frac 0.80`, `window_frac 0.05`) the cap cost **−0.8 / −1.9 / −0.2 pp recall** for only a
    /// **1.48× / 1.44× / 1.35× detect** speedup; the fragmented glyco file was always the worst case
    /// because much of its *real* low-abundance signal lives exactly where a tile's local reject rate
    /// crosses the threshold. `frac 0.50` over-cut badly (−9 pp on two of three files). The recall loss
    /// tracks the ΣTIC loss, i.e. the cap removes genuine explained signal, not just noise. **Verdict: not
    /// worth it — ship uncapped.** Kept opt-in for callers that want to trade recall for detect speed on
    /// large files. NB: the *real* bound on an "uncapped" run is [`Self::min_seed_intensity`], not this
    /// cap; a data-dependent seed floor (see `agent_info/TODO.md`) is the better lever and the intended
    /// successor to this experiment.
    pub reject_stop_enabled: bool,
    /// Width (in scored seeds) of the rolling window over which the reject fraction is measured *on the
    /// serial path*. The parallel paths size each tile's window from data instead — see
    /// [`tile_reject_cfg`] and `DETECT_TILE2D_REJECT_WINDOW_FRAC`.
    pub reject_stop_window: usize,
    /// Stop when the window's reject fraction reaches this value. Default `0.80` (the A/B sweet spot; see
    /// [`Self::reject_stop_enabled`]) — `0.50` was too eager and over-cut real low-abundance features.
    pub reject_stop_frac: f64,
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
            knee_stop_enabled: false,
            knee_window_seeds: 20_000,
            knee_slope_frac: 0.02,
            knee_abs_eps: 1e-7,
            reject_stop_enabled: false,
            reject_stop_window: 20_000,
            reject_stop_frac: 0.80,
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
    engine: &impl PeakSource,
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
    engine: &impl PeakSource,
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
/// Two steps, split so step 1 can run **before** charge scoring as a cheap persistence pre-gate:
/// 1. [`trace_seed_extent`] — follow the most-abundant tooth (the seed's m/z) to fix `[s_lo, s_hi]`.
/// 2. [`gather_extent_peaks`] — collect every comb tooth's peaks across that extent.
///
/// **Step 1: extent.** Follow the most-abundant tooth (the seed's m/z — highest SNR, most reliable
/// boundary) outward in RT with [`PeakIndexingEngine::get_xic_by_scan_index`], stopping on
/// `trace_missed_scans_allowed` consecutive misses or the `trace_max_half_width_minutes` guard. Its
/// peaks' scan indices give the extent `[s_lo, s_hi]`. Peaks already in `claimed` count as misses (a
/// taller neighbour claimed them first), which is what splits co-eluting same-m/z peaks greedily
/// instead of merging them. Charge-independent (uses only the seed m/z), so it is safe to run before
/// scoring: a seed whose own XIC spans too few scans is a single-scan noise spike and can be retired
/// without paying for the (six-charge) envelope scoring — that is where most of the tail's wasted
/// compute goes (the persistence gate, not the envelope gate, rejects the bulk of low-abundance seeds).
fn trace_seed_extent(
    engine: &impl PeakSource,
    seed: &IndexedMassSpectralPeak,
    params: &TraceKernelParameters,
    ppm: &PpmTolerance,
    claimed: &HashSet<PeakKey>,
) -> (i32, i32) {
    let seed_mz = seed.m() as f64;
    let apex_scan = seed.zero_based_scan_index;
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
    (s_lo, s_hi)
}

/// **Step 2: gather.** Collect every comb tooth's peak at each scan in `[s_lo, s_hi]`, excluding
/// anything already `claimed`. The union (deduped) is the feature's peak set — the single set that
/// backs its RT bounds, summed intensity, coverage contribution, and the NMS claim mask alike. Needs
/// the accepted hypothesis (charge + monoisotopic m/z), so it runs after scoring; the extent it walks
/// is the one [`trace_seed_extent`] already fixed (no second XIC walk).
fn gather_extent_peaks(
    engine: &impl PeakSource,
    seed: &IndexedMassSpectralPeak,
    hyp: &HypothesisScore,
    s_lo: i32,
    s_hi: i32,
    params: &TraceKernelParameters,
    ppm: &PpmTolerance,
    claimed: &HashSet<PeakKey>,
) -> Vec<IndexedMassSpectralPeak> {
    let charge = hyp.charge;
    let seed_mz = seed.m() as f64;
    let weights = comb_weights(mz_to_mass(seed_mz, charge), params);
    if weights.is_empty() {
        // Degenerate comb — nothing to gather; claim the scored peaks (minus any already claimed).
        return hyp
            .peaks
            .iter()
            .filter(|p| !claimed.contains(&p.key()))
            .copied()
            .collect();
    }
    let spacing = C13_MINUS_C12 / charge as f64;
    let mono_mz = hyp.mono_mz;

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

/// Running tally of what happens to each *scored* seed (one that was unclaimed and reached charge
/// scoring — claimed seeds are skipped before this and never counted). `scored == accepted +
/// rej_env + rej_pers`. Reported in the progress line so the seed-rejection rate (fraction of
/// considered seeds that fail to yield a viable envelope) can be plotted vs progress.
#[derive(Clone, Copy, Default)]
struct SeedTally {
    /// Unclaimed seeds that reached charge scoring (the denominator of "considered").
    scored: u64,
    /// Rejected: no viable envelope (too few isotopes observed, or non-positive response).
    rej_env: u64,
    /// Rejected: had an envelope but its traced extent failed the persistence gate.
    rej_pers: u64,
}

/// Opt-in live progress reporter for the detect loop (env `DETECT_PROGRESS`). Loop-local and
/// non-`Sync`: constructed once before the serial acceptance walk, mutated only from that walk.
///
/// Throttled so the hot path pays only a single integer compare per accepted feature before any
/// `Instant::now()`/formatting: emit only after both `every_seeds` seeds AND `every` wall-time have
/// elapsed since the last line. Output is a single greppable stderr line prefixed
/// `[DETECT_PROGRESS]`; the feature TSV (stdout/files) is never touched.
struct DetectProgress {
    enabled: bool,
    every_seeds: u64,
    every: Duration,
    start_time: Instant,
    last_emit_seeds: u64,
    last_emit_time: Instant,
    total_tic: f64,
    pool_total: usize,
}

impl DetectProgress {
    /// Throttled emit. Ordered cheapest-check-first so the common (no-emit) path is a single
    /// integer subtraction+compare and never calls `Instant::now()` or allocates.
    fn maybe_emit(&mut self, seeds_considered: u64, accepted: usize, explained: f64, seed_intensity: f64, tally: SeedTally) {
        if !self.enabled {
            return;
        }
        if seeds_considered.saturating_sub(self.last_emit_seeds) < self.every_seeds {
            return;
        }
        let now = Instant::now();
        if now.duration_since(self.last_emit_time) < self.every {
            return;
        }
        self.last_emit_seeds = seeds_considered;
        self.last_emit_time = now;
        self.emit(seeds_considered, accepted, explained, seed_intensity, tally);
    }

    /// Unconditional emit (used for the final line at loop exit); still gated on `enabled`.
    fn final_emit(&self, seeds_considered: u64, accepted: usize, explained: f64, seed_intensity: f64, tally: SeedTally) {
        if self.enabled {
            self.emit(seeds_considered, accepted, explained, seed_intensity, tally);
        }
    }

    /// `seed_intensity` is the intensity of the seed being processed at this emit. Because seeds are
    /// visited intensity-descending, it is the effective *seed-intensity floor* reached so far — i.e.
    /// every peak still unclaimed below it is untouched. This is the physical quantity a S/N-based
    /// auto-stop keys on, and lets a coverage %ΣTIC be mapped to the seed intensity that produced it.
    /// `tally` carries the cumulative scored/rejected seed counts for the rejection-rate curve.
    fn emit(&self, seeds_considered: u64, accepted: usize, explained: f64, seed_intensity: f64, tally: SeedTally) {
        let pool_total = self.pool_total as u64;
        let pool_remaining = pool_total.saturating_sub(seeds_considered.min(pool_total));
        let pct = if self.total_tic > 0.0 {
            100.0 * explained / self.total_tic
        } else {
            0.0
        };
        eprintln!(
            "[DETECT_PROGRESS] seeds {}/{} (pool_remaining {}) | accepted {} | \
             seed_int {:.0} | scored {} | rej_env {} | rej_pers {} | \
             TIC {:.2e} ({:.1}% of ΣTIC) | {:.1}s",
            seeds_considered,
            pool_total,
            pool_remaining,
            accepted,
            seed_intensity,
            tally.scored,
            tally.rej_env,
            tally.rej_pers,
            explained,
            pct,
            self.start_time.elapsed().as_secs_f64(),
        );
    }
}

/// Opt-in auto-stop at the *knee* of the coverage curve (params `knee_stop_enabled`; driver env
/// `DETECT_KNEE`). Loop-local, non-`Sync`, constructed once before the serial acceptance walk.
///
/// Since seeds are visited tallest-first, the marginal %ΣTIC claimed per seed is monotonically
/// non-increasing — the explained-vs-seeds curve is concave — so a rolling 2-point slope over a
/// fixed seed window captures the knee with O(1) state (no history buffer). `observe` returns `true`
/// exactly once, when the loop should break.
struct KneeDetector {
    enabled: bool,
    window_seeds: u64,
    slope_frac: f64,
    abs_eps: f64,
    total_tic: f64,
    win_start_seeds: u64,
    win_start_explained: f64,
    early_slope: Option<f64>,
    windows_seen: u32,
}

impl KneeDetector {
    /// Feed the running (seeds_considered, explained) state. Returns `true` when the rolling slope
    /// has collapsed — below `slope_frac × early_slope` (relative) OR below `abs_eps` (absolute,
    /// flat tail) — but only after at least two full windows have elapsed (the first sets the
    /// reference slope). A no-op returning `false` when disabled.
    fn observe(&mut self, seeds_considered: u64, explained: f64) -> bool {
        if !self.enabled {
            return false;
        }
        let dseeds = seeds_considered.saturating_sub(self.win_start_seeds);
        if dseeds < self.window_seeds {
            return false;
        }
        // A full window elapsed: rolling slope = Δexplained/Δseeds, normalised by ΣTIC so it reads
        // as a fraction of total signal claimed per seed.
        let denom = if self.total_tic > 0.0 { self.total_tic } else { 1.0 };
        let rolling_slope = (explained - self.win_start_explained) / dseeds as f64 / denom;
        self.windows_seen += 1;
        // Advance the window origin for the next window.
        self.win_start_seeds = seeds_considered;
        self.win_start_explained = explained;
        // The first completed window only establishes the reference slope; never stop on it.
        let early = match self.early_slope {
            None => {
                self.early_slope = Some(rolling_slope);
                return false;
            }
            Some(s) => s,
        };
        if self.windows_seen < 2 {
            return false;
        }
        rolling_slope < self.slope_frac * early || rolling_slope < self.abs_eps
    }
}

/// Opt-in auto-stop keyed on the detector's *seed-rejection rate* (params `reject_stop_enabled`;
/// driver env `DETECT_REJECT_STOP`). Loop-local, O(1) state. Fires once the rolling fraction of
/// considered (scored) seeds being rejected — no persistent envelope — reaches `frac` over a window
/// of `window` scored seeds. Because seeds are visited tallest-first this fraction rises
/// monotonically, so a single rolling window (no reference/warmup) captures the crossing; unlike the
/// TIC knee it needs no ΣTIC normalisation and is self-referential (portable across files/instruments).
struct RejectStop {
    enabled: bool,
    window: u64,
    frac: f64,
    win_start_scored: u64,
    win_start_rejected: u64,
    /// Reject fraction of the last completed window (for the stop-line log).
    last_rate: f64,
}

impl RejectStop {
    /// Feed the running cumulative (scored, rejected) counts. Returns `true` when the most recent full
    /// window's reject fraction reached `frac`. A no-op returning `false` when disabled.
    fn observe(&mut self, scored: u64, rejected: u64) -> bool {
        if !self.enabled {
            return false;
        }
        let dscored = scored.saturating_sub(self.win_start_scored);
        if dscored < self.window {
            return false;
        }
        let drej = rejected.saturating_sub(self.win_start_rejected);
        let rate = drej as f64 / dscored as f64;
        self.last_rate = rate;
        self.win_start_scored = scored;
        self.win_start_rejected = rejected;
        rate >= self.frac
    }

    /// Build a fresh loop-local instance. `window` is resolved per tile from the tile's own walked-seed
    /// count (see [`detect_bin`]), so every tile scales its rolling window to how many seeds it actually
    /// considers — the windows are independent per tile (they run on separate threads, never shared).
    fn new(enabled: bool, window: u64, frac: f64) -> Self {
        Self {
            enabled,
            window: window.max(1),
            frac,
            win_start_scored: 0,
            win_start_rejected: 0,
            last_rate: 0.0,
        }
    }
}

/// Per-tile configuration for the reject-rate auto-stop on the parallel detectors. `Copy` so each tile's
/// [`detect_bin`] constructs its OWN [`RejectStop`] from it — the windows are independent per tile (that
/// is exactly what makes the cap compose with the parallel schedule: no global ΣTIC, no cross-tile
/// coordination, unlike the coverage-target and knee stops).
///
/// `window_frac` is the rolling-window width as a **fraction of the tile's own walked-seed count**
/// (data-dependent), resolved to an absolute seed count inside [`detect_bin`]; `frac` is the reject-rate
/// threshold at which the tile stops.
#[derive(Clone, Copy)]
struct RejectStopCfg {
    enabled: bool,
    window_frac: f64,
    frac: f64,
}

impl RejectStopCfg {
    /// The cap turned off — the serial reference and every test that does not exercise it pass this.
    const DISABLED: Self = Self {
        enabled: false,
        window_frac: 0.0,
        frac: 1.0,
    };
}

/// What happened to one seed run through [`process_seed`]. The three-way split lets the per-tile
/// reject-rate cap count **only genuine no-signal rejections** — a seed skipped because a neighbour tile
/// (via strip-seed) or a re-anchored apex already claimed its peaks is not evidence the tail has turned to
/// noise, and counting it would make border/late-colour tiles stop prematurely and under-detect.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum SeedOutcome {
    /// A feature was emitted.
    Accepted,
    /// The seed was scored but yielded no persistent envelope (persistence or envelope/isotope gate).
    RejectedNoSignal,
    /// The seed's peaks were already claimed; nothing was scored. Invisible to the reject-rate cap.
    SkippedClaimed,
}

/// Runs the untargeted MS1 detector.
///
/// **Defaults to the intra-file 2-D m/z×RT tiling path** ([`detect_features_tile2d`]): the file is cut
/// into a 4-coloured grid of m/z×RT tiles detected in parallel (rayon), with per-tile seed re-anchoring so
/// an isotope tooth split across an m/z border cannot mis-anchor. Env overrides:
/// - **`DETECT_SERIAL`** — force the single-threaded greedy reference ([`detect_features_serial`]): the
///   tallest-first, claim-as-you-go path the tiling detector is validated bit-against.
/// - **`DETECT_PARALLEL`** — use the 1-D red-black RT-binning path ([`detect_features_parallel`]) instead
///   of the 2-D grid.
///
/// Neither parallel path is bit-identical to serial — a feature straddling a tile/bin boundary can be
/// claimed by the earlier phase regardless of global intensity order — so both are **automatically
/// disabled (→ serial) whenever a *global* stopping heuristic (coverage target < 1, knee stop) is
/// engaged**, since those need a global running ΣTIC that cannot be evaluated inside one independent
/// tile/bin. The **reject-rate** stop is the exception — being self-referential it is applied per
/// tile/bin (see [`detect_bin`] / [`tile_reject_cfg`]), so a reject-capped run stays parallel. A tiny
/// input that cannot be split (too few scans / too narrow an m/z range) also falls back to serial.
pub fn detect_features(
    engine: &PeakIndexingEngine,
    params: &TraceKernelParameters,
) -> Vec<DetectedFeature> {
    // 2-D tiling is the default; DETECT_SERIAL forces the reference path, DETECT_PARALLEL selects the 1-D
    // RT-binning path. The coverage target and knee stops need a *global* running ΣTIC, so they cannot be
    // evaluated inside an independent tile/bin and force serial. The reject-rate stop is deliberately NOT
    // one of these — it is self-referential and applied per tile inside [`detect_bin`] via
    // [`tile_reject_cfg`], so a reject-capped run stays parallel.
    let force_serial = std::env::var("DETECT_SERIAL").is_ok();
    let parallel_1d_requested = std::env::var("DETECT_PARALLEL").is_ok();
    let global_stop_engaged = params.coverage_target < 1.0 || params.knee_stop_enabled;

    if !force_serial && !global_stop_engaged {
        let parallel = if parallel_1d_requested {
            detect_features_parallel(engine, params)
        } else {
            detect_features_tile2d(engine, params)
        };
        if let Some(features) = parallel {
            return features;
        }
    }
    detect_features_serial(engine, params)
}

fn detect_features_serial(
    engine: &PeakIndexingEngine,
    params: &TraceKernelParameters,
) -> Vec<DetectedFeature> {
    let ppm = PpmTolerance::new(params.ppm_tolerance);

    // Opt-in sub-stage profiling (env `DETECT_PROFILE=1`). All timing work is gated behind this
    // bool, read once here, so the default path pays only a predictable-branch check per section —
    // no `Instant::now()` in the hot loop unless profiling is explicitly requested.
    let profile = std::env::var("DETECT_PROFILE").is_ok();
    let prof_t0 = Instant::now();
    let mut t_window = 0.0f64;
    let mut t_score = 0.0f64;
    let mut t_trace = 0.0f64;
    let mut n_seed_considered = 0u64;
    let mut n_score_calls = 0u64;

    // Seeds: all peaks, tallest first. Stable ordering (intensity desc) mirrors get_all_xics.
    let mut seeds = engine.all_peaks();
    seeds.sort_by(|a, b| b.intensity.total_cmp(&a.intensity));
    let t_seed_prep = prof_t0.elapsed().as_secs_f64();

    // Opt-in live progress reporting (env `DETECT_PROGRESS`; throttle override `DETECT_PROGRESS_EVERY`
    // in seeds, default 50_000). Read once here so the hot loop only checks a bool.
    let progress_enabled = std::env::var("DETECT_PROGRESS").is_ok();
    let progress_every: u64 = std::env::var("DETECT_PROGRESS_EVERY")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .filter(|&n| n > 0)
        .unwrap_or(50_000);

    // Coverage bookkeeping: denominator = Σ all peak intensities (ΣTIC). The sum is only needed when
    // a consumer engages — the coverage cap (`coverage_target < 1.0`), live progress reporting, OR
    // the knee auto-stop — so the common "detect everything, no progress" default still skips the
    // O(peaks) pass.
    let need_total_tic = params.coverage_target < 1.0
        || progress_enabled
        || params.knee_stop_enabled
        || params.reject_stop_enabled;
    let total_intensity: f64 = if need_total_tic {
        seeds.iter().map(|p| p.intensity as f64).sum()
    } else {
        0.0
    };
    let coverage_stop = if params.coverage_target < 1.0 && total_intensity > 0.0 {
        params.coverage_target * total_intensity
    } else {
        f64::INFINITY
    };

    let mut progress = DetectProgress {
        enabled: progress_enabled,
        every_seeds: progress_every,
        every: Duration::from_millis(500),
        start_time: prof_t0,
        last_emit_seeds: 0,
        last_emit_time: prof_t0,
        total_tic: total_intensity,
        pool_total: seeds.len(),
    };

    // Opt-in auto-stop at the coverage knee (params `knee_stop_enabled`; default OFF).
    let mut knee = KneeDetector {
        enabled: params.knee_stop_enabled,
        window_seeds: params.knee_window_seeds.max(1) as u64,
        slope_frac: params.knee_slope_frac,
        abs_eps: params.knee_abs_eps,
        total_tic: total_intensity,
        win_start_seeds: 0,
        win_start_explained: 0.0,
        early_slope: None,
        windows_seen: 0,
    };
    let knee_log = progress_enabled || params.knee_stop_enabled;

    // Opt-in auto-stop at the seed-rejection-rate threshold (params `reject_stop_enabled`; default OFF).
    let mut reject_stop = RejectStop {
        enabled: params.reject_stop_enabled,
        window: (params.reject_stop_window.max(1)) as u64,
        frac: params.reject_stop_frac,
        win_start_scored: 0,
        win_start_rejected: 0,
        last_rate: 0.0,
    };
    let reject_log = progress_enabled || params.reject_stop_enabled;

    let mut explained_intensity = 0.0;

    let mut claimed: HashSet<PeakKey> = HashSet::new();
    let mut features: Vec<DetectedFeature> = Vec::new();
    let mut seeds_visited: u64 = 0;
    // Intensity of the most recent above-floor seed; == the effective seed-intensity floor at any
    // loop-exit point (seeds are intensity-descending). Reported by the final progress line.
    let mut last_seed_intensity: f64 = 0.0;
    // Cumulative outcome tally over scored seeds (feeds the rejection-rate curve in the progress line).
    let mut tally = SeedTally::default();

    for (seed_idx, seed) in seeds.iter().enumerate() {
        seeds_visited = seed_idx as u64 + 1;
        // Seeds are intensity-descending, so once we fall below the floor every remaining seed is
        // too — stop rather than continue.
        if (seed.intensity as f64) < params.min_seed_intensity {
            break;
        }
        last_seed_intensity = seed.intensity as f64;
        if claimed.contains(&seed.key()) {
            continue;
        }

        // Auto-stop once the rolling reject rate crosses the threshold (opt-in; checked before the
        // current seed so it fires even during a long run of consecutive rejections). Uses the
        // cumulative tally of already-completed seeds.
        if reject_stop.observe(tally.scored, tally.rej_env + tally.rej_pers) {
            if reject_log {
                let pct = if total_intensity > 0.0 {
                    100.0 * explained_intensity / total_intensity
                } else {
                    0.0
                };
                eprintln!(
                    "[DETECT_REJECT_STOP] stopping at seed {}/{} | scored {} | accepted {} | \
                     window reject {:.1}% | {:.1}% of ΣTIC | {:.1}s",
                    seed_idx as u64 + 1,
                    seeds.len(),
                    tally.scored,
                    features.len(),
                    100.0 * reject_stop.last_rate,
                    pct,
                    prof_t0.elapsed().as_secs_f64(),
                );
            }
            break;
        }

        // Committed to considering this seed: it counts toward the rejection-rate curve.
        tally.scored += 1;
        if profile {
            n_seed_considered += 1;
        }

        // Trace the most-abundant tooth (the seed's own m/z) FIRST, before any charge scoring. Its
        // scan span is a cheap, charge-independent persistence pre-gate: a seed whose XIC spans fewer
        // scans than `min_feature_scans` is a single-scan noise spike that the persistence gate would
        // reject anyway (the gathered extent can span no more scans than this seed XIC), so retire it
        // now and skip the six-charge envelope scoring. This front-loads the tail's dominant rejection
        // (persistence, not envelope) ahead of its most expensive step.
        let ts = if profile { Some(Instant::now()) } else { None };
        let (s_lo, s_hi) = trace_seed_extent(engine, seed, params, &ppm, &claimed);
        if let Some(ts) = ts {
            t_trace += ts.elapsed().as_secs_f64();
        }
        if params.min_feature_scans > 1 && (s_hi - s_lo + 1) < params.min_feature_scans as i32 {
            tally.rej_pers += 1;
            claimed.insert(seed.key());
            continue;
        }

        // The RT window (scan indices + Gaussian weights) is charge-independent — compute it once
        // per seed and share it across every charge hypothesis.
        let ts = if profile { Some(Instant::now()) } else { None };
        let window = seed_rt_window(engine, seed, params);
        if let Some(ts) = ts {
            t_window += ts.elapsed().as_secs_f64();
        }

        // Score every charge hypothesis; keep the highest response (cross-z non-max suppression).
        let ts = if profile { Some(Instant::now()) } else { None };
        let mut best: Option<HypothesisScore> = None;
        for z in params.min_charge..=params.max_charge {
            if z == 0 {
                continue;
            }
            if profile {
                n_score_calls += 1;
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
        if let Some(ts) = ts {
            t_score += ts.elapsed().as_secs_f64();
        }

        let best = match best {
            Some(b) => b,
            None => {
                tally.rej_env += 1;
                continue;
            }
        };

        if best.num_isotopes_observed < params.min_isotopes_observed || best.response <= 0.0 {
            // Not a feature; retire this seed so we do not reconsider it.
            tally.rej_env += 1;
            claimed.insert(seed.key());
            continue;
        }

        // Gather the feature's TRUE traced extent (not just the narrow scored window) over the span
        // fixed above, so the whole elution is claimed at once and its smaller adjacent seeds cannot
        // re-fire as fragments.
        let ts = if profile { Some(Instant::now()) } else { None };
        let traced = gather_extent_peaks(engine, seed, &best, s_lo, s_hi, params, &ppm, &claimed);
        if let Some(ts) = ts {
            t_trace += ts.elapsed().as_secs_f64();
        }

        // Full chromatographic-persistence gate: the pre-gate above bounds the *span*, but the gathered
        // teeth can still cover fewer than `min_feature_scans` *distinct* scans if they are sparse
        // within that span. Retire such a seed (as with the isotope-count gate) without emitting it.
        if params.min_feature_scans > 1 {
            let distinct_scans: HashSet<i32> =
                traced.iter().map(|p| p.zero_based_scan_index).collect();
            if distinct_scans.len() < params.min_feature_scans {
                tally.rej_pers += 1;
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

        // Live progress (throttled, opt-in; no-op on the default path).
        progress.maybe_emit(seed_idx as u64 + 1, features.len(), explained_intensity, seed.intensity as f64, tally);

        // Stop once we have explained the target fraction of the total MS1 signal.
        if explained_intensity >= coverage_stop {
            break;
        }

        // Auto-stop at the coverage knee (opt-in; can fire earlier than the coverage target).
        if knee.observe(seed_idx as u64 + 1, explained_intensity) {
            if knee_log {
                let pct = if total_intensity > 0.0 {
                    100.0 * explained_intensity / total_intensity
                } else {
                    0.0
                };
                eprintln!(
                    "[DETECT_KNEE] stopping at seed {}/{} | accepted {} | {:.1}% of ΣTIC | {:.1}s",
                    seed_idx as u64 + 1,
                    seeds.len(),
                    features.len(),
                    pct,
                    prof_t0.elapsed().as_secs_f64(),
                );
            }
            break;
        }
    }

    // Final progress line at loop exit (regardless of throttle), so the last state is always logged.
    progress.final_emit(seeds_visited, features.len(), explained_intensity, last_seed_intensity, tally);

    if profile {
        let total = prof_t0.elapsed().as_secs_f64();
        eprintln!(
            "  [DETECT_PROFILE] total {:.2}s | seed_prep(all_peaks+sort) {:.2}s | \
             rt_window {:.2}s | score_hypothesis {:.2}s | trace_claim {:.2}s | \
             other {:.2}s",
            total,
            t_seed_prep,
            t_window,
            t_score,
            t_trace,
            (total - t_seed_prep - t_window - t_score - t_trace).max(0.0),
        );
        eprintln!(
            "  [DETECT_PROFILE] seeds_considered(unclaimed) {} | score_hypothesis calls {} | \
             features {}",
            n_seed_considered,
            n_score_calls,
            features.len(),
        );
    }

    features
}

/// Minimum RT width, in minutes, of a detector bin. Chosen well above `2 × reach` (the maximum RT
/// distance any per-seed computation extends — see [`bin_reach_minutes`]) so that two bins processed
/// concurrently in the same red-black phase — always separated by at least one full bin — can never
/// claim the same peak. At the observed data-driven windows (tens of seconds) 4 min leaves a wide
/// margin; do not lower it without re-deriving the separation guarantee.
const DETECT_BIN_MIN_WIDTH_MINUTES: f64 = 4.0;

/// Target m/z column width (daltons) for the 2-D tiling detector (§2). Deliberately **wide** — far
/// above the strip-pad-derived floor — so only a small fraction of seeds sit near an m/z border and the
/// re-anchoring halo scan (§6a) rarely fires; the RT axis supplies the bulk of the parallelism.
/// Overridable via `DETECT_TILE2D_MZ_WIDTH` for tuning (raised to the strip-pad-derived floor if smaller).
const DETECT_TILE_MZ_WIDTH_DALTONS: f64 = 96.0;

/// Default m/z **re-anchoring reach** (thomson) for the 2-D tiling detector — the halo radius and the
/// border-seed threshold (§6a). This is the *envelope span* in m/z: how far a minor isotope tooth can sit
/// from its own envelope apex. Measured on the IonStar 2-hr gradient (1.1 M detections), the per-feature
/// m/z reach `= observed_isotopes × (C13−C12)/z` is bounded in **thomson and shrinks with charge** (isotope
/// count grows with mass but the `1.0033/z` spacing shrinks faster): z=1 is the worst case at ~3 Th median /
/// 7 Th max, everything heavier is tighter. Overall p99 = 4.0 Th, p99.9 = 5.0 Th; only ~2 % of features
/// exceed 4 Th and ~0.008 % exceed 6 Th. So **4 Th covers ~98 % of straddles directly**; the rare tail
/// beyond is a mis-anchored bogus feature backstopped by the collision detector + fall-through (§6a residual
/// edge), not a lost or double-claimed peak. Chosen over 6 Th because the halo scan cost grows ~radius²
/// (scan width × border-seed fraction) — 6→4 Th cut detect time 16 % (10-min) / 6 % (2-hr) with <0.1 %
/// feature change. Deliberately **smaller** than the theoretical claim bound
/// (`max_isotopes × (C13−C12)/min_charge ≈ 12 Da`), which is used only for the strip pad / tile floor (§4).
/// Overridable via `DETECT_TILE2D_REACH_MZ`.
const DETECT_TILE_REACH_MZ_DALTONS: f64 = 4.0;

/// Default rolling-window width for the **per-tile** reject-rate auto-stop on the parallel detectors, as a
/// **fraction of each tile's own walked-seed count** (data-dependent). A fixed count cannot serve tiles
/// that span three orders of magnitude in seed count (short-gradient tiles vs. fragmented-glyco tiles); a
/// fraction sizes each tile's window to how many seeds it actually walks. Overridable via
/// `DETECT_TILE2D_REJECT_WINDOW_FRAC`. Only consulted when the run enables the reject stop
/// (`reject_stop_enabled`, itself off by default); the reject-fraction threshold is shared with the serial
/// path (`reject_stop_frac`). Default `0.05` = the A/B sweet spot (see [`TraceKernelParameters::reject_stop_enabled`]):
/// 1% was faster but roughly doubled the recall cost on the fragmented glyco file.
const DETECT_TILE_REJECT_WINDOW_FRAC: f64 = 0.05;

/// Build the per-tile reject-rate cap config for the parallel paths from the run params. Disabled unless
/// the run requested the reject stop; when enabled, uses the data-dependent window fraction (default
/// [`DETECT_TILE_REJECT_WINDOW_FRAC`], overridable via `DETECT_TILE2D_REJECT_WINDOW_FRAC`) and the run's
/// `reject_stop_frac`. This is what lets a *capped* run stay on the parallel path instead of falling back
/// to serial — the cap is applied independently inside each tile's [`detect_bin`].
fn tile_reject_cfg(params: &TraceKernelParameters) -> RejectStopCfg {
    if !params.reject_stop_enabled {
        return RejectStopCfg::DISABLED;
    }
    let window_frac = std::env::var("DETECT_TILE2D_REJECT_WINDOW_FRAC")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .filter(|&v| v > 0.0 && v <= 1.0)
        .unwrap_or(DETECT_TILE_REJECT_WINDOW_FRAC);
    RejectStopCfg {
        enabled: true,
        window_frac,
        frac: params.reject_stop_frac,
    }
}

/// The maximum RT half-width (minutes) any single seed's scoring or claim-trace reaches from its
/// apex. A feature seeded in one bin can therefore touch peaks at most this far into an adjacent bin,
/// and no further. The scoring window is `rt_half_window_minutes`; the claim-extent trace is capped
/// at `trace_max_half_width_minutes`; the gathered extent never exceeds the trace. The bin width must
/// exceed `2 ×` this for the within-phase no-collision guarantee to hold.
fn bin_reach_minutes(params: &TraceKernelParameters) -> f64 {
    params
        .rt_half_window_minutes
        .max(params.trace_max_half_width_minutes)
}

/// Inclusive `[scan_lo, scan_hi]` zero-based scan indices covering every scan whose RT falls in
/// `[rt_lo, rt_hi]`, plus a one-scan guard each side. `scan_info` is RT-ascending. Used to size a
/// [`PeakIndexView`]'s scan window to a tile's padded RT box: every scan any of the tile's queries can
/// walk to lies within `rt_lo..rt_hi` (the box already carries the ± reach margin), so the guarded
/// window is a strict superset of the queried scans and the view returns identical results.
fn scan_bounds_for_rt(scan_info: &[ScanInfo], rt_lo: f64, rt_hi: f64) -> (i32, i32) {
    let n = scan_info.len();
    if n == 0 {
        return (0, -1);
    }
    let lo = scan_info.partition_point(|s| s.retention_time < rt_lo);
    let hi = scan_info.partition_point(|s| s.retention_time <= rt_hi);
    let scan_lo = (lo as i64 - 1).max(0) as i32;
    let scan_hi = (hi as i64).min(n as i64 - 1) as i32;
    (scan_lo, scan_hi)
}

/// Whether the zero-copy [`PeakIndexView`] narrowing is active. **On by default** — it is a
/// byte-identical detect speedup (≈1.43× on the 65-min 2-D path; the narrowed per-tile slices give the
/// hot lookups shorter, cache-local binary searches). Opt out with `DETECT_INDEX_VIEW=0` (or `false`)
/// to fall back to the full shared engine.
fn index_view_enabled() -> bool {
    !matches!(
        std::env::var("DETECT_INDEX_VIEW").as_deref(),
        Ok("0") | Ok("false")
    )
}

/// Intra-file red-black RT-binning detector (opt-in; see [`detect_features`]). Returns `None` when the
/// run is too short to split into ≥ 2 bins, in which case the caller falls back to the serial path.
///
/// **Binning.** Seeds (all peaks) are partitioned into contiguous RT bins that are simultaneously
/// *equal-work* (≈ equal seed count, the load-balancing target) and *≥ [`DETECT_BIN_MIN_WIDTH_MINUTES`]
/// wide* (the correctness floor). A bin is closed only once both hold, so dense regions produce
/// min-width bins and sparse regions widen until they carry a full share of seeds.
///
/// **Red-black scheduling.** Bins are 2-colored by index parity. Even bins are detected fully in
/// parallel (phase 0), a barrier merges their claims, then odd bins are detected in parallel (phase 1)
/// against those claims. Because every pair of same-color bins is separated by a full (≥ 4 min) bin and
/// a seed reaches at most [`bin_reach_minutes`] (≪ 2 min) from its apex, concurrently-processed bins
/// have disjoint claim sets — so each bin runs the ordinary serial [`detect_bin`] with its own
/// `claimed` set and no locking on the hot path.
///
/// **Parity.** Not bit-identical to serial: at a bin boundary the even phase claims first regardless of
/// which side holds the taller feature, so a small, boundary-localized set of features differ from the
/// global tallest-first result. Everything away from a boundary is identical.
fn detect_features_parallel(
    engine: &PeakIndexingEngine,
    params: &TraceKernelParameters,
) -> Option<Vec<DetectedFeature>> {
    let ppm = PpmTolerance::new(params.ppm_tolerance);
    let profile = std::env::var("DETECT_PROFILE").is_ok();
    let t0 = Instant::now();

    let scan_info = engine.scan_info();
    if scan_info.len() < 2 {
        return None;
    }
    let rt_min = scan_info.first().unwrap().retention_time;
    let rt_max = scan_info.last().unwrap().retention_time;
    let rt_span = rt_max - rt_min;
    if rt_span < 2.0 * DETECT_BIN_MIN_WIDTH_MINUTES {
        return None; // too short to yield ≥ 2 min-width bins; not worth splitting
    }

    // Seeds ordered globally tallest-first, then re-grouped by bin. `all_peaks` mirrors the serial
    // seed pool exactly.
    let seeds = engine.all_peaks();
    if seeds.is_empty() {
        return Some(Vec::new());
    }

    // --- Build bin boundaries: equal-work with a hard min-width floor. -----------------------------
    let threads = rayon::current_num_threads().max(1);
    // Aim for a few bins per thread per phase so rayon can balance; more bins = finer balance, and the
    // min-width floor caps how many actually fit. `2 × threads` per phase → `4 × threads` total target.
    let target_bins = (((rt_span / DETECT_BIN_MIN_WIDTH_MINUTES).floor() as usize).max(2))
        .min((4 * threads).max(2));
    let target_work = seeds.len().div_ceil(target_bins).max(1);

    // Boundaries are ascending RT cut points; a seed with `rt < cut` belongs to the lower bin. Built
    // by walking seeds in RT order and closing a bin once it holds ≥ target_work seeds AND spans
    // ≥ the min width.
    let mut sorted_rts: Vec<f64> = seeds.iter().map(|p| p.retention_time as f64).collect();
    sorted_rts.sort_by(f64::total_cmp);
    let mut boundaries: Vec<f64> = Vec::new();
    let mut count = 0usize;
    let mut bin_start_rt = rt_min;
    for &rt in &sorted_rts {
        count += 1;
        if count >= target_work && (rt - bin_start_rt) >= DETECT_BIN_MIN_WIDTH_MINUTES {
            boundaries.push(rt);
            bin_start_rt = rt;
            count = 0;
        }
    }
    drop(sorted_rts);
    let n_bins = boundaries.len() + 1;
    if n_bins < 2 {
        return None;
    }

    // Bin index of an RT = number of boundaries at or below it.
    let bin_of = |rt: f64| -> usize { boundaries.partition_point(|&b| b <= rt) };

    // --- Group seeds by bin, tallest-first within each bin, without copying peaks. -----------------
    // Sort an index permutation by (bin asc, intensity desc), then materialize one reordered peak
    // Vec so each bin is a contiguous slice. `seeds` itself is dropped afterward, so peak memory
    // matches the serial path (one seed pool resident).
    let seed_bins: Vec<u32> = seeds
        .iter()
        .map(|p| bin_of(p.retention_time as f64) as u32)
        .collect();
    let mut order: Vec<u32> = (0..seeds.len() as u32).collect();
    order.sort_by(|&a, &b| {
        let (a, b) = (a as usize, b as usize);
        seed_bins[a]
            .cmp(&seed_bins[b])
            .then_with(|| seeds[b].intensity.total_cmp(&seeds[a].intensity))
    });
    let ordered: Vec<IndexedMassSpectralPeak> =
        order.iter().map(|&i| seeds[i as usize]).collect();
    let ordered_bins: Vec<u32> = order.iter().map(|&i| seed_bins[i as usize]).collect();
    drop(seeds);
    drop(seed_bins);
    drop(order);

    // Contiguous [start, end) slice range for each bin in `ordered`.
    let mut ranges: Vec<(usize, usize)> = vec![(0, 0); n_bins];
    let mut start = 0usize;
    for (b, range) in ranges.iter_mut().enumerate() {
        let mut end = start;
        while end < ordered.len() && ordered_bins[end] as usize == b {
            end += 1;
        }
        *range = (start, end);
        start = end;
    }
    drop(ordered_bins);

    // Bin edges (n_bins + 1 RT points): [rt_min, boundaries.., rt_max⁺]. Bin b spans
    // [bin_edges[b], bin_edges[b+1]). The final edge is nudged past rt_max so the last seed is inside.
    let mut bin_edges: Vec<f64> = Vec::with_capacity(n_bins + 1);
    bin_edges.push(rt_min);
    bin_edges.extend_from_slice(&boundaries);
    bin_edges.push(rt_max + 1.0);

    let reach = bin_reach_minutes(params);
    // Per-bin reject-rate cap (opt-in): each bin stops its own noise tail independently, so a capped run
    // stays parallel instead of falling back to serial. Disabled cfg when the run didn't request it.
    let reject = tile_reject_cfg(params);
    // Zero-copy narrowing (on by default; DETECT_INDEX_VIEW=0 opts out): each bin runs against a
    // [`PeakIndexView`] restricted to its RT span ± `reach` (full m/z — the 1-D path tiles RT only), so
    // the hot lookups binary-search short scan-window slices instead of whole-run bins. Byte-identical.
    let use_view = index_view_enabled();

    // --- Phase 0: even bins in parallel, each with a fresh (empty) claim set. ----------------------
    let even_indices: Vec<usize> = (0..n_bins).step_by(2).collect();
    let even_results: Vec<(usize, Vec<DetectedFeature>, Vec<IndexedMassSpectralPeak>)> = even_indices
        .par_iter()
        .map(|&b| {
            let (lo, hi) = ranges[b];
            let (features, claimed) = if use_view {
                let (s_lo, s_hi) =
                    scan_bounds_for_rt(engine.scan_info(), bin_edges[b] - reach, bin_edges[b + 1] + reach);
                let view = engine.view_scans(s_lo, s_hi);
                detect_bin(&view, params, &ppm, &ordered[lo..hi], HashSet::new(), None, reject)
            } else {
                detect_bin(engine, params, &ppm, &ordered[lo..hi], HashSet::new(), None, reject)
            };
            (b, features, claimed)
        })
        .collect();

    // Merge even-phase claims, sorted by RT, so each odd bin can be seeded from the thin boundary
    // strip within `reach` of its span (all it can possibly observe from the earlier phase).
    let mut even_claims: Vec<IndexedMassSpectralPeak> = even_results
        .iter()
        .flat_map(|(_, _, claimed)| claimed.iter().copied())
        .collect();
    even_claims.sort_by(|a, b| a.retention_time.total_cmp(&b.retention_time));
    let claim_rts: Vec<f64> = even_claims.iter().map(|p| p.retention_time as f64).collect();

    // --- Phase 1: odd bins in parallel, each pre-seeded with the even claims in its widened span. --
    let odd_indices: Vec<usize> = (1..n_bins).step_by(2).collect();
    let odd_results: Vec<(usize, Vec<DetectedFeature>)> = odd_indices
        .par_iter()
        .map(|&b| {
            let lo_rt = bin_edges[b] - reach;
            let hi_rt = bin_edges[b + 1] + reach;
            let s = claim_rts.partition_point(|&r| r < lo_rt);
            let e = claim_rts.partition_point(|&r| r <= hi_rt);
            let seeded: HashSet<PeakKey> = even_claims[s..e].iter().map(|p| p.key()).collect();
            let (lo, hi) = ranges[b];
            let (features, _) = if use_view {
                let (s_lo, s_hi) = scan_bounds_for_rt(engine.scan_info(), lo_rt, hi_rt);
                let view = engine.view_scans(s_lo, s_hi);
                detect_bin(&view, params, &ppm, &ordered[lo..hi], seeded, None, reject)
            } else {
                detect_bin(engine, params, &ppm, &ordered[lo..hi], seeded, None, reject)
            };
            (b, features)
        })
        .collect();

    // --- Reassemble in bin order (deterministic regardless of thread scheduling). ------------------
    let mut per_bin: Vec<Vec<DetectedFeature>> = (0..n_bins).map(|_| Vec::new()).collect();
    for (b, features, _) in even_results {
        per_bin[b] = features;
    }
    for (b, features) in odd_results {
        per_bin[b] = features;
    }
    let features: Vec<DetectedFeature> = per_bin.into_iter().flatten().collect();

    if profile {
        eprintln!(
            "  [DETECT_PARALLEL] total {:.2}s | {} bins ({} even, {} odd) | {} threads | \
             reach {:.2} min / min-width {:.1} min | index {} | features {}",
            t0.elapsed().as_secs_f64(),
            n_bins,
            even_indices.len(),
            odd_indices.len(),
            threads,
            reach,
            DETECT_BIN_MIN_WIDTH_MINUTES,
            if use_view { "view" } else { "shared" },
            features.len(),
        );
    }

    Some(features)
}

/// The 4 colour of a tile from its grid coordinates: `(i_rt & 1, i_mz & 1)` packed as
/// `((i_rt & 1) << 1) | (i_mz & 1)` ∈ 0..4. Same-colour tiles differ by ≥ 2 on at least one axis
/// (§3), hence are separated by a full margin-respecting tile and can never claim the same peak.
#[inline]
fn tile_color(i_rt: usize, i_mz: usize) -> u8 {
    (((i_rt & 1) << 1) | (i_mz & 1)) as u8
}

/// Intra-file **2-D m/z×RT tiling** detector (the **default** path — see [`detect_features`]; opt out with
/// `DETECT_SERIAL` — and `agent_info/Detector-2D-Tiling-Spec.md`). Generalises the 1-D red-black RT path
/// ([`detect_features_parallel`]) to a 2-D grid 4-coloured so non-adjacent tiles run concurrently.
/// Returns `None` (caller falls back to serial) when the run is too short/narrow to yield ≥ 2
/// margin-respecting tiles on each axis.
///
/// **Why 2-D.** The RT claim reach grows with the trace cap (forcing coarse RT slices), but the m/z
/// claim reach is fixed by the isotope model — so tiling the m/z axis keeps exposing independent
/// regions exactly when RT slices go coarse.
///
/// **Grid.** RT bands are equal-work with a `2·reach_rt` (≥ 4 min) floor; m/z columns are fixed ~96-Da
/// wide (§2). Every tile is validated ≥ `2·reach` on each axis so same-colour tiles are provably
/// collision-free (§4).
///
/// **Schedule.** Colours `(i_rt&1, i_mz&1)` run in series 0→3; within a colour tiles run in parallel,
/// each strip-seeded with the claims of its already-processed (earlier-colour) king-neighbours in the
/// padded box (§5), then run through the ordinary serial [`detect_bin`] with seed re-anchoring (§6a).
/// A per-colour barrier merges each tile's claims into a global set, and that merge *is* the collision
/// detector (§6): a duplicate insert (impossible in nominal margin-respecting mode) is counted and the
/// loser dropped. Features are reassembled in fixed (colour, tile) order — deterministic run-to-run.
fn detect_features_tile2d(
    engine: &PeakIndexingEngine,
    params: &TraceKernelParameters,
) -> Option<Vec<DetectedFeature>> {
    let ppm = PpmTolerance::new(params.ppm_tolerance);
    let profile = std::env::var("DETECT_PROFILE").is_ok();
    let t0 = Instant::now();

    let scan_info = engine.scan_info();
    if scan_info.len() < 2 {
        return None;
    }
    let rt_min = scan_info.first().unwrap().retention_time;
    let rt_max = scan_info.last().unwrap().retention_time;
    let rt_span = rt_max - rt_min;

    // --- Reach on each axis (§1). ------------------------------------------------------------------
    // Two distinct m/z quantities (decoupled for performance — the spec conflates them):
    //  • `reach_mz` — the re-anchoring halo radius / border threshold (§6a). This is the envelope *span*
    //    (how far a minor tooth sits from its apex), ~6 Da; a wider halo just costs more for no gain.
    //  • `claim_reach_mz` — the theoretical bound on how far *any* seed's comb reaches from it:
    //    `max_isotopes × (C13−C12)/min_charge` (~12 Da at z=1). This is what a claim can actually cross a
    //    border by, so the strip pad and tile floor must be built from it (else a neighbour re-claims the
    //    far tail of a long z=1 comb → a collision). RT reach grows with the trace cap; RT re-anchoring
    //    stays on a single scan, so the RT reach is undoubled.
    let reach_rt = bin_reach_minutes(params);
    let claim_reach_mz =
        params.max_isotopes as f64 * (C13_MINUS_C12 / params.min_charge.max(1) as f64);
    let reach_mz = std::env::var("DETECT_TILE2D_REACH_MZ")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .filter(|&v| v > 0.0)
        .unwrap_or(DETECT_TILE_REACH_MZ_DALTONS)
        .min(claim_reach_mz);
    // Worst-case m/z distance a claim reaches past a tile border: a re-anchored border seed finds its
    // apex up to `reach_mz` across the border, and that apex's own comb reaches a further `claim_reach_mz`
    // → `reach_mz + claim_reach_mz`. (An interior seed reaches only `claim_reach_mz`, which is smaller.)
    // The strip pad must cover it so the neighbour strips those peaks; the column floor is `2×` the pad so
    // two same-colour columns, separated by one intervening column, cannot both reach into its middle.
    let strip_pad_mz = reach_mz + claim_reach_mz;
    let min_axis_mz = 2.0 * strip_pad_mz;
    let mz_width = std::env::var("DETECT_TILE2D_MZ_WIDTH")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .filter(|&v| v >= min_axis_mz)
        .unwrap_or(DETECT_TILE_MZ_WIDTH_DALTONS.max(min_axis_mz));

    let seeds = engine.all_peaks();
    if seeds.is_empty() {
        return Some(Vec::new());
    }

    // --- m/z columns: fixed-width, global cut points, thin trailing column merged into its neighbour.
    let (mz_min, mz_max) = seeds.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
        let m = p.m() as f64;
        (lo.min(m), hi.max(m))
    });
    let mut mz_edges: Vec<f64> = vec![mz_min];
    let mut e = mz_min + mz_width;
    while e < mz_max {
        mz_edges.push(e);
        e += mz_width;
    }
    mz_edges.push(mz_max + 1e-6); // nudge the final edge past the max so that seed lands in the last column
    // Coarsen a too-thin trailing column (the remainder) until it clears the `min_axis_mz` floor.
    while mz_edges.len() >= 3 {
        let n = mz_edges.len();
        if mz_edges[n - 1] - mz_edges[n - 2] < min_axis_mz {
            mz_edges.remove(n - 2);
        } else {
            break;
        }
    }
    let n_mz = mz_edges.len() - 1;
    if n_mz < 2 {
        if profile {
            eprintln!(
                "  [DETECT_TILE2D] fallback: m/z range {:.1}–{:.1} Da too narrow for ≥ 2 columns of {:.0} Da",
                mz_min, mz_max, mz_width
            );
        }
        return None; // too narrow for 2 m/z tiles — the 1-D path or serial covers this
    }

    // --- RT bands: equal-work with a `2·reach_rt` (≥ 4 min) floor (mirror the 1-D walk). ------------
    let rt_floor = (2.0 * reach_rt).max(DETECT_BIN_MIN_WIDTH_MINUTES);
    let max_rt_bands = (rt_span / rt_floor).floor() as usize;
    if max_rt_bands < 2 {
        if profile {
            eprintln!(
                "  [DETECT_TILE2D] fallback: RT span {:.1} min too short for ≥ 2 bands of {:.1} min",
                rt_span, rt_floor
            );
        }
        return None;
    }
    // Aim for ~2–4 tiles per thread per colour → ~12·threads total tiles; the m/z axis already fixes
    // `n_mz` columns, so target the RT band count that hits the tile budget, clamped to what the floor
    // allows.
    let threads = rayon::current_num_threads().max(1);
    let target_total = (3 * threads * 4).max(4);
    let target_n_rt = (target_total / n_mz).clamp(2, max_rt_bands);
    let target_work = seeds.len().div_ceil(target_n_rt).max(1);

    let mut sorted_rts: Vec<f64> = seeds.iter().map(|p| p.retention_time as f64).collect();
    sorted_rts.sort_by(f64::total_cmp);
    let mut rt_cuts: Vec<f64> = Vec::new();
    let mut count = 0usize;
    let mut band_start_rt = rt_min;
    for &rt in &sorted_rts {
        count += 1;
        if count >= target_work && (rt - band_start_rt) >= rt_floor {
            rt_cuts.push(rt);
            band_start_rt = rt;
            count = 0;
        }
    }
    drop(sorted_rts);
    let n_rt = rt_cuts.len() + 1;
    if n_rt < 2 {
        return None;
    }

    // RT edges (n_rt + 1 points): [rt_min, cuts.., rt_max⁺]. Band i spans [rt_edges[i], rt_edges[i+1]).
    let mut rt_edges: Vec<f64> = Vec::with_capacity(n_rt + 1);
    rt_edges.push(rt_min);
    rt_edges.extend_from_slice(&rt_cuts);
    rt_edges.push(rt_max + 1.0);

    let n_tiles = n_rt * n_mz;

    // Interior cut points for O(log) tile lookup (partition_point counts cuts at/below the coordinate).
    let rt_band_of = |rt: f64| -> usize { rt_cuts.partition_point(|&b| b <= rt) };
    let mz_col_of = |mz: f64| -> usize { mz_edges[1..n_mz].partition_point(|&b| b <= mz) };

    // --- Group seeds by tile, tallest-first within each tile, without copying peaks. ----------------
    // Sort an index permutation by (tile asc, intensity desc), then materialise one reordered peak Vec
    // so each tile is a contiguous slice (mirrors the 1-D path; one seed pool resident).
    let seed_tiles: Vec<u32> = seeds
        .iter()
        .map(|p| (rt_band_of(p.retention_time as f64) * n_mz + mz_col_of(p.m() as f64)) as u32)
        .collect();
    let mut order: Vec<u32> = (0..seeds.len() as u32).collect();
    order.sort_by(|&a, &b| {
        let (a, b) = (a as usize, b as usize);
        seed_tiles[a]
            .cmp(&seed_tiles[b])
            .then_with(|| seeds[b].intensity.total_cmp(&seeds[a].intensity))
    });
    let ordered: Vec<IndexedMassSpectralPeak> = order.iter().map(|&i| seeds[i as usize]).collect();
    let ordered_tiles: Vec<u32> = order.iter().map(|&i| seed_tiles[i as usize]).collect();
    drop(seeds);
    drop(seed_tiles);
    drop(order);

    // Contiguous [start, end) slice range for each tile in `ordered`.
    let mut ranges: Vec<(usize, usize)> = vec![(0, 0); n_tiles];
    let mut start = 0usize;
    for (t, range) in ranges.iter_mut().enumerate() {
        let mut end = start;
        while end < ordered.len() && ordered_tiles[end] as usize == t {
            end += 1;
        }
        *range = (start, end);
        start = end;
    }
    drop(ordered_tiles);

    // --- 4-colour strip-seeded schedule (§3, §5). --------------------------------------------------
    // Claims accumulate across colours: `per_tile_claims[t]` is the peaks tile `t` newly claimed (used
    // to strip-seed its later-colour king-neighbours); `global_claimed` is the merge/collision detector.
    let mut per_tile_claims: Vec<Vec<IndexedMassSpectralPeak>> =
        (0..n_tiles).map(|_| Vec::new()).collect();
    let mut per_tile_features: Vec<Vec<DetectedFeature>> =
        (0..n_tiles).map(|_| Vec::new()).collect();
    let mut global_claimed: HashSet<PeakKey> = HashSet::new();
    let mut collisions = 0usize;

    // Per-tile reject-rate cap (opt-in): each tile stops its own noise tail independently. This is what
    // keeps a capped run on the 2-D path — the coverage/knee stops can't (they need a global ΣTIC), but
    // the reject rate is self-referential per tile. Disabled cfg when the run didn't request it.
    let reject = tile_reject_cfg(params);
    // Zero-copy narrowing (on by default; DETECT_INDEX_VIEW=0 opts out): each tile runs against a
    // [`PeakIndexView`] restricted to its padded m/z × RT box, so the hot lookups binary-search short
    // slices. The box is exactly the strip-pad / reach margin the tiling already enforces, so the view is
    // a strict superset of every query the tile makes → byte-identical output.
    let use_view = index_view_enabled();

    for color in 0..4u8 {
        let tiles: Vec<usize> = (0..n_tiles)
            .filter(|&t| tile_color(t / n_mz, t % n_mz) == color)
            .collect();
        if tiles.is_empty() {
            continue;
        }
        // Snapshot the accumulated claims so the parallel closure borrows them immutably.
        let claims_ref = &per_tile_claims;
        let results: Vec<(usize, Vec<DetectedFeature>, Vec<IndexedMassSpectralPeak>)> = tiles
            .par_iter()
            .map(|&t| {
                let i_rt = t / n_mz;
                let i_mz = t % n_mz;
                let box_rt_lo = rt_edges[i_rt] - reach_rt;
                let box_rt_hi = rt_edges[i_rt + 1] + reach_rt;
                let box_mz_lo = mz_edges[i_mz] - strip_pad_mz;
                let box_mz_hi = mz_edges[i_mz + 1] + strip_pad_mz;

                // Strip-seed from earlier-colour king-neighbours' claims that fall in the padded box.
                let mut seeded: HashSet<PeakKey> = HashSet::new();
                for di in -1i64..=1 {
                    for dj in -1i64..=1 {
                        if di == 0 && dj == 0 {
                            continue;
                        }
                        let ni = i_rt as i64 + di;
                        let nj = i_mz as i64 + dj;
                        if ni < 0 || ni >= n_rt as i64 || nj < 0 || nj >= n_mz as i64 {
                            continue;
                        }
                        let (ni, nj) = (ni as usize, nj as usize);
                        if tile_color(ni, nj) >= color {
                            continue; // only already-processed (earlier) colours hold claims
                        }
                        for p in &claims_ref[ni * n_mz + nj] {
                            let rt = p.retention_time as f64;
                            let mz = p.m() as f64;
                            if rt >= box_rt_lo
                                && rt <= box_rt_hi
                                && mz >= box_mz_lo
                                && mz <= box_mz_hi
                            {
                                seeded.insert(p.key());
                            }
                        }
                    }
                }

                let reanchor = Reanchor {
                    reach_mz,
                    mz_lo: mz_edges[i_mz],
                    mz_hi: mz_edges[i_mz + 1],
                };
                let (lo, hi) = ranges[t];
                let (features, claimed) = if use_view {
                    let (s_lo, s_hi) = scan_bounds_for_rt(engine.scan_info(), box_rt_lo, box_rt_hi);
                    let view = engine.view_box(box_mz_lo, box_mz_hi, s_lo, s_hi);
                    detect_bin(
                        &view,
                        params,
                        &ppm,
                        &ordered[lo..hi],
                        seeded,
                        Some(&reanchor),
                        reject,
                    )
                } else {
                    detect_bin(
                        engine,
                        params,
                        &ppm,
                        &ordered[lo..hi],
                        seeded,
                        Some(&reanchor),
                        reject,
                    )
                };
                (t, features, claimed)
            })
            .collect();

        // Barrier: merge each tile's claims serially (single thread) — the merge is the detector.
        for (t, features, claimed) in results {
            for p in &claimed {
                if !global_claimed.insert(p.key()) {
                    collisions += 1;
                }
            }
            per_tile_claims[t] = claimed;
            per_tile_features[t] = features;
        }
    }

    // --- Reassemble in fixed (colour, tile) order — deterministic regardless of thread scheduling. --
    let mut features: Vec<DetectedFeature> = Vec::new();
    for color in 0..4u8 {
        for t in 0..n_tiles {
            if tile_color(t / n_mz, t % n_mz) == color {
                features.append(&mut per_tile_features[t]);
            }
        }
    }

    if collisions > 0 {
        // Nominal (margin-respecting, re-anchored) mode must produce zero collisions — a nonzero count
        // is the correctness canary that the width floor was violated or reach was mis-set.
        eprintln!(
            "[DETECT_TILE2D] WARNING: {collisions} same-colour claim collisions (expected 0 in \
             margin-respecting mode) — losers dropped; check reach_mz / tile floor."
        );
    }
    if profile {
        let reject_desc = if reject.enabled {
            format!(
                " | reject: frac {:.2}, window {:.1}% of tile walked-seeds",
                reject.frac,
                100.0 * reject.window_frac
            )
        } else {
            String::new()
        };
        eprintln!(
            "  [DETECT_TILE2D] total {:.2}s | grid {}×{} = {} tiles ({} RT bands, {} m/z cols) | \
             {} threads | reach_rt {:.2} min / reach_mz {:.1} Da (strip pad {:.1} Da) | \
             mz_width {:.0} Da | index {} | collisions {} | features {}{}",
            t0.elapsed().as_secs_f64(),
            n_rt,
            n_mz,
            n_tiles,
            n_rt,
            n_mz,
            threads,
            reach_rt,
            reach_mz,
            strip_pad_mz,
            mz_width,
            if use_view { "view" } else { "shared" },
            collisions,
            features.len(),
            reject_desc,
        );
    }

    Some(features)
}

/// Per-tile seed re-anchoring context for the 2-D tiling detector (§6a). Present only on the 2-D
/// path; the serial reference and the 1-D RT-binning path pass `None`, keeping them byte-identical.
///
/// `mz_lo`/`mz_hi` are the tile's own m/z column bounds; a seed within `reach_mz` of either is a
/// *border* seed whose envelope apex might sit across the border, and only those seeds pay the halo
/// scan. `reach_mz` is the halo radius (also the strip pad and the tile-width floor — see §1).
struct Reanchor {
    reach_mz: f64,
    mz_lo: f64,
    mz_hi: f64,
}

impl Reanchor {
    /// True when `seed_mz` is within `reach_mz` of this tile's m/z border — the only seeds whose
    /// envelope can straddle into a neighbouring tile, so the only ones that need the halo scan.
    #[inline]
    fn is_border(&self, seed_mz: f64) -> bool {
        seed_mz - self.mz_lo < self.reach_mz || self.mz_hi - seed_mz < self.reach_mz
    }
}

/// Processes one seed through the serial accept/reject pipeline against `claimed`, pushing an accepted
/// feature onto `features` and every newly-claimed peak onto both `claimed` and `claimed_peaks`. Factored
/// out of [`detect_bin`] so a re-anchored apex peak (§6a) can be run through the identical logic before
/// its originating border seed. A seed already in `claimed` is a no-op (so a re-anchored apex that also
/// appears later as its own tile's seed is processed exactly once).
#[allow(clippy::too_many_arguments)]
fn process_seed(
    engine: &impl PeakSource,
    params: &TraceKernelParameters,
    ppm: &PpmTolerance,
    seed: &IndexedMassSpectralPeak,
    claimed: &mut HashSet<PeakKey>,
    claimed_peaks: &mut Vec<IndexedMassSpectralPeak>,
    features: &mut Vec<DetectedFeature>,
) -> SeedOutcome {
    if claimed.contains(&seed.key()) {
        return SeedOutcome::SkippedClaimed;
    }

    // Cheap charge-independent persistence pre-gate (same as serial): a seed whose own-m/z XIC
    // spans too few scans is a noise spike; retire it before the six-charge scoring.
    let (s_lo, s_hi) = trace_seed_extent(engine, seed, params, ppm, claimed);
    if params.min_feature_scans > 1 && (s_hi - s_lo + 1) < params.min_feature_scans as i32 {
        if claimed.insert(seed.key()) {
            claimed_peaks.push(*seed);
        }
        return SeedOutcome::RejectedNoSignal;
    }

    let window = seed_rt_window(engine, seed, params);
    let mut best: Option<HypothesisScore> = None;
    for z in params.min_charge..=params.max_charge {
        if z == 0 {
            continue;
        }
        let score = score_hypothesis(engine, seed, z, params, ppm, claimed, &window);
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
        None => return SeedOutcome::RejectedNoSignal,
    };
    if best.num_isotopes_observed < params.min_isotopes_observed || best.response <= 0.0 {
        if claimed.insert(seed.key()) {
            claimed_peaks.push(*seed);
        }
        return SeedOutcome::RejectedNoSignal;
    }

    let traced = gather_extent_peaks(engine, seed, &best, s_lo, s_hi, params, ppm, claimed);
    if params.min_feature_scans > 1 {
        let distinct_scans: HashSet<i32> =
            traced.iter().map(|p| p.zero_based_scan_index).collect();
        if distinct_scans.len() < params.min_feature_scans {
            if claimed.insert(seed.key()) {
                claimed_peaks.push(*seed);
            }
            return SeedOutcome::RejectedNoSignal;
        }
    }

    for p in &traced {
        if claimed.insert(p.key()) {
            claimed_peaks.push(*p);
        }
    }
    features.push(build_feature(best, traced));
    SeedOutcome::Accepted
}

/// Detects features within a single tile (or RT bin): the serial greedy accept/reject loop, run over one
/// tile's `seeds` (already tallest-first) against a tile-local `claimed` set. Pre-seed `claimed` with any
/// prior phase's claims that overlap this tile (odd RT bin, or an earlier-colour neighbour's strip); pass
/// an empty set for a fresh bin.
///
/// `reanchor` is `Some` only on the 2-D tiling path: for a seed within `reach_mz` of an m/z tile border it
/// runs the seed-re-anchoring repair (§6a) — search the m/z halo at the seed's apex scan for the tallest
/// *unclaimed* peak and process THAT first, so a minor tooth whose apex straddled into a neighbour tile
/// cannot anchor a mis-placed comb. `None` (serial reference and 1-D RT binning, where an envelope's teeth
/// share a scan and never split across the tiling axis) keeps the loop byte-identical to before.
///
/// Returns the tile's accepted features and **the peaks it newly claimed** (paired with their RT/m-z via the
/// peak record). The `insert` return value guarantees only genuinely-new keys are reported, so a pre-seeded
/// key is never echoed back. This mirrors the body of [`detect_features_serial`] minus the global progress /
/// coverage / knee bookkeeping, which does not compose with per-tile execution.
///
/// `reject` is the per-tile reject-rate auto-stop (§ parallel-capping): unlike the coverage/knee stops it
/// *does* compose with tiling because it is self-referential (fraction of this tile's own scored seeds that
/// yield no signal) — no global ΣTIC, no cross-tile state. Each tile runs its own rolling window and stops
/// its tail independently. `RejectStopCfg::DISABLED` (serial reference, 1-D-without-cap, tests) keeps the
/// loop running to the seed floor exactly as before. Only [`SeedOutcome::RejectedNoSignal`] feeds it —
/// seeds skipped as already-claimed (a neighbour's strip-seed or a re-anchored apex) must not, or borders
/// and late colours would stop early.
#[allow(clippy::too_many_arguments)]
fn detect_bin(
    engine: &impl PeakSource,
    params: &TraceKernelParameters,
    ppm: &PpmTolerance,
    seeds: &[IndexedMassSpectralPeak],
    mut claimed: HashSet<PeakKey>,
    reanchor: Option<&Reanchor>,
    reject: RejectStopCfg,
) -> (Vec<DetectedFeature>, Vec<IndexedMassSpectralPeak>) {
    let mut features: Vec<DetectedFeature> = Vec::new();
    let mut claimed_peaks: Vec<IndexedMassSpectralPeak> = Vec::new();

    // Per-tile reject-rate auto-stop state. `scored`/`rejected` are this tile's cumulative no-skip counts
    // (accepted + no-signal, and no-signal alone) — exactly the serial tally, scoped to one tile. The
    // rolling window is sized from THIS tile's walked-seed count (seeds at/above the intensity floor — the
    // ones the loop will actually consider), so it is data-dependent and scales with tile size rather than
    // a fixed count. `seeds` is intensity-descending, so the above-floor count is a partition point.
    let n_walked = seeds.partition_point(|s| (s.intensity as f64) >= params.min_seed_intensity);
    let window = (reject.window_frac * n_walked as f64).ceil() as u64;
    let mut reject_stop = RejectStop::new(reject.enabled, window, reject.frac);
    let mut scored: u64 = 0;
    let mut rejected: u64 = 0;

    for seed in seeds {
        // Seeds are intensity-descending within the tile, so once one falls below the floor every
        // remaining seed does too.
        if (seed.intensity as f64) < params.min_seed_intensity {
            break;
        }

        // Once this tile's rolling reject fraction crosses the threshold its tail has turned to noise —
        // stop. Checked before the seed (as in serial) so it fires even mid-run of consecutive rejects.
        // A no-op returning false when the cap is disabled.
        if reject_stop.observe(scored, rejected) {
            break;
        }

        // §6a seed re-anchoring: a border seed may be a minor tooth of an envelope whose apex sits in
        // a neighbouring m/z tile. Enforce "the taller peak goes first" locally — find the tallest
        // unclaimed peak in the m/z halo at this seed's scan and run it first. It is unconditionally
        // safe (if it is this seed's true apex the comb anchors correctly and claims the seed as a
        // tooth; if it is a distinct co-eluting species both are detected). Falls through to the seed
        // below, which `process_seed` no-ops if the re-anchored feature already claimed it.
        if let Some(re) = reanchor {
            let seed_mz = seed.m() as f64;
            if re.is_border(seed_mz) && !claimed.contains(&seed.key()) {
                if let Some(apex) = engine.tallest_unclaimed_in_mz_at_scan(
                    seed_mz - re.reach_mz,
                    seed_mz + re.reach_mz,
                    seed.zero_based_scan_index,
                    &claimed,
                ) {
                    if apex.key() != seed.key() {
                        let outcome = process_seed(
                            engine,
                            params,
                            ppm,
                            &apex,
                            &mut claimed,
                            &mut claimed_peaks,
                            &mut features,
                        );
                        tally_seed_outcome(outcome, &mut scored, &mut rejected);
                    }
                }
            }
        }

        let outcome = process_seed(
            engine,
            params,
            ppm,
            seed,
            &mut claimed,
            &mut claimed_peaks,
            &mut features,
        );
        tally_seed_outcome(outcome, &mut scored, &mut rejected);
    }

    (features, claimed_peaks)
}

/// Fold one [`SeedOutcome`] into a tile's `(scored, rejected)` reject-rate counters. Mirrors the serial
/// tally: a scored seed is one that was not skipped-as-claimed; a rejected seed is a scored one that
/// produced no feature. `SkippedClaimed` touches neither, so it is invisible to the reject-rate cap.
#[inline]
fn tally_seed_outcome(outcome: SeedOutcome, scored: &mut u64, rejected: &mut u64) {
    match outcome {
        SeedOutcome::Accepted => *scored += 1,
        SeedOutcome::RejectedNoSignal => {
            *scored += 1;
            *rejected += 1;
        }
        SeedOutcome::SkippedClaimed => {}
    }
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

    #[test]
    fn knee_detector_stops_on_concave_sequence() {
        // Saturating (concave) coverage curve explained(s) = TOTAL*(1 - exp(-s/tau)). Sampling at
        // window boundaries s = k*W (with W = tau) gives a per-window slope ratio r_k = exp(-(k-1))
        // versus the first window. For slope_frac = 0.02 the first window with r_k < 0.02 is k=5
        // (r_5 = exp(-4) ≈ 0.0183), so the detector must first return true at seed 5*W.
        let total = 1.0_f64;
        let w = 100u64;
        let tau = 100.0_f64;
        let mut knee = KneeDetector {
            enabled: true,
            window_seeds: w,
            slope_frac: 0.02,
            abs_eps: 1e-9,
            total_tic: total,
            win_start_seeds: 0,
            win_start_explained: 0.0,
            early_slope: None,
            windows_seen: 0,
        };
        let explained = |s: f64| total * (1.0 - (-s / tau).exp());
        let mut stopped_at: Option<u64> = None;
        for k in 1..=6u64 {
            let seeds = k * w;
            if knee.observe(seeds, explained(seeds as f64)) {
                stopped_at = Some(k);
                break;
            }
        }
        assert_eq!(
            stopped_at,
            Some(5),
            "knee should first fire at window k=5 for slope_frac=0.02"
        );

        // Absolute-floor (flatline) path: with slope_frac = 0 the relative test can never fire on a
        // positive slope, so a near-flat tail must stop solely via abs_eps.
        let mut flat = KneeDetector {
            enabled: true,
            window_seeds: 10,
            slope_frac: 0.0,
            abs_eps: 1e-3,
            total_tic: 1.0,
            win_start_seeds: 0,
            win_start_explained: 0.0,
            early_slope: None,
            windows_seen: 0,
        };
        assert!(!flat.observe(10, 0.5), "first window only sets the reference slope");
        // Second window slope ≈ 1e-5/seed << abs_eps 1e-3 (and relative test disabled) => stop.
        assert!(flat.observe(20, 0.5001), "flat tail must stop via abs_eps");
    }

    #[test]
    fn knee_detector_disabled_never_stops() {
        // Disabled detector must return false forever, even with thresholds that would trivially
        // fire (tiny window, huge slope_frac/abs_eps) on a flat sequence.
        let mut knee = KneeDetector {
            enabled: false,
            window_seeds: 1,
            slope_frac: 1.0,
            abs_eps: 1.0,
            total_tic: 1.0,
            win_start_seeds: 0,
            win_start_explained: 0.0,
            early_slope: None,
            windows_seen: 0,
        };
        for k in 1..=1000u64 {
            assert!(!knee.observe(k, 0.0), "disabled knee must never stop");
        }
    }

    #[test]
    fn per_tile_reject_cap_counts_only_no_signal_rejects() {
        // The per-tile reject-rate cap must count only genuine no-signal rejections. This locks the one
        // subtle invariant: a seed skipped because a neighbour (strip-seed) or a re-anchored apex already
        // claimed its peaks is NOT evidence the tail turned to noise, and must not advance the window —
        // else border/late-colour tiles would stop early and under-detect. `window` here is the absolute
        // per-tile window that `detect_bin` would resolve from the tile's walked-seed count.
        let (window, frac) = (10u64, 0.5f64);

        // 1) A thousand SkippedClaimed outcomes advance neither counter and never trip the cap.
        let mut rs = RejectStop::new(true, window, frac);
        let (mut scored, mut rejected) = (0u64, 0u64);
        for _ in 0..1000 {
            assert!(
                !rs.observe(scored, rejected),
                "skipped-claimed seeds must not trip the per-tile cap"
            );
            tally_seed_outcome(SeedOutcome::SkippedClaimed, &mut scored, &mut rejected);
        }
        assert_eq!(
            (scored, rejected),
            (0, 0),
            "SkippedClaimed advances neither the scored nor the rejected counter"
        );

        // 2) A full window of no-signal rejects (100% > 50%) trips exactly when the window first fills.
        let mut rs = RejectStop::new(true, window, frac);
        let (mut scored, mut rejected) = (0u64, 0u64);
        let mut tripped_at = None;
        for i in 0..25u64 {
            if rs.observe(scored, rejected) {
                tripped_at = Some(i);
                break;
            }
            tally_seed_outcome(SeedOutcome::RejectedNoSignal, &mut scored, &mut rejected);
        }
        assert_eq!(
            tripped_at,
            Some(10),
            "cap trips at the first full window of all-reject seeds"
        );

        // 3) A window at exactly the fraction (5 rejects / 10 scored = 0.5) trips (rate >= frac).
        let mut rs = RejectStop::new(true, window, frac);
        let (mut scored, mut rejected) = (0u64, 0u64);
        let mut tripped = false;
        for i in 0..40u64 {
            if rs.observe(scored, rejected) {
                tripped = true;
                break;
            }
            let o = if i % 2 == 0 {
                SeedOutcome::RejectedNoSignal
            } else {
                SeedOutcome::Accepted
            };
            tally_seed_outcome(o, &mut scored, &mut rejected);
        }
        assert!(tripped, "a window at exactly frac (0.5) must trip");

        // 4) A DISABLED cap never trips, even on an all-reject stream.
        let mut rs = RejectStop::new(false, 1, 1.0);
        let (mut scored, mut rejected) = (0u64, 0u64);
        for _ in 0..1000 {
            assert!(
                !rs.observe(scored, rejected),
                "disabled per-tile cap must never trip"
            );
            tally_seed_outcome(SeedOutcome::RejectedNoSignal, &mut scored, &mut rejected);
        }
    }

    #[test]
    fn reanchoring_recovers_apex_when_seeded_from_minor_tooth() {
        // §6a: a clean z=2 envelope whose monoisotope is the most-abundant tooth (the apex). Seed the
        // detector from the +1 *minor* tooth alone — exactly what an m/z tile border does when the apex
        // lands in the neighbouring column. WITHOUT re-anchoring the comb anchors its most-abundant
        // tooth on the +1 peak → mono mis-placed by ~1 unit (a bogus feature). WITH re-anchoring the
        // halo scan finds the taller mono and processes the true envelope, so the recovered feature is
        // identical to seeding from the apex — regardless of which side's tile runs first.
        let (scans, mono_mz) = averagine_envelope_scans(1200.0, 2, 9);
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let params = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.15,
            rt_half_window_minutes: 0.3,
            ..TraceKernelParameters::default()
        };
        let ppm = PpmTolerance::new(params.ppm_tolerance);
        let spacing = C13_MINUS_C12 / 2.0;
        let apex_scan = 4;

        let mono_peak = *engine.get_indexed_peak(mono_mz, apex_scan, &ppm).expect("mono tooth");
        let plus1_peak = *engine
            .get_indexed_peak(mono_mz + spacing, apex_scan, &ppm)
            .expect("+1 tooth");
        assert!(
            plus1_peak.intensity < mono_peak.intensity,
            "the +1 tooth must be a minor tooth for this fixture"
        );

        // reach_mz comfortably covers a tooth spacing; the tight mz_lo/mz_hi make the seed a border seed.
        let ctx = |p: &IndexedMassSpectralPeak| Reanchor {
            reach_mz: 6.0,
            mz_lo: p.m() as f64 - 0.1,
            mz_hi: p.m() as f64 + 0.1,
        };

        // Seed from the minor tooth, re-anchoring ON.
        let re_b = ctx(&plus1_peak);
        let (feat_b, _) = detect_bin(
            &engine,
            &params,
            &ppm,
            &[plus1_peak],
            HashSet::new(),
            Some(&re_b),
            RejectStopCfg::DISABLED,
        );
        // Seed from the apex, re-anchoring ON — the order-independence reference.
        let re_a = ctx(&mono_peak);
        let (feat_a, _) = detect_bin(
            &engine,
            &params,
            &ppm,
            &[mono_peak],
            HashSet::new(),
            Some(&re_a),
            RejectStopCfg::DISABLED,
        );

        assert_eq!(feat_b.len(), 1, "re-anchored minor-tooth seed yields exactly one feature");
        assert_eq!(feat_a.len(), 1);
        assert_eq!(feat_b[0].charge, 2, "recovered charge");
        approx(feat_b[0].monoisotopic_mass, 1200.0, 0.02);
        // Same envelope regardless of which tooth seeded it.
        assert_eq!(feat_b[0].charge, feat_a[0].charge);
        approx(feat_b[0].monoisotopic_mass, feat_a[0].monoisotopic_mass, 1e-6);

        // Contrast: WITHOUT re-anchoring the minor-tooth seed mis-places the mono (off by ~1 unit),
        // so it must NOT recover 1200 — proving re-anchoring is what fixes the straddle.
        let (feat_none, _) = detect_bin(
            &engine,
            &params,
            &ppm,
            &[plus1_peak],
            HashSet::new(),
            None,
            RejectStopCfg::DISABLED,
        );
        assert!(
            feat_none.is_empty() || (feat_none[0].monoisotopic_mass - 1200.0).abs() > 0.1,
            "without re-anchoring the minor-tooth seed should mis-place the mono, got {:?}",
            feat_none.first().map(|f| f.monoisotopic_mass)
        );
    }

    /// Eight well-separated clean z=2 elutions spread across RT and m/z — wide/long enough to yield a
    /// ≥ 2×2 margin-respecting tile grid. Each envelope is comfortably isolated so there is no
    /// co-elution or ambiguous boundary; the 2-D path should match serial exactly here.
    fn tiling_fixture() -> Vec<Scan> {
        let specs = [
            (800.0, 1.0),
            (1000.0, 4.0),
            (1300.0, 7.0),
            (1600.0, 10.0),
            (900.0, 13.0),
            (1100.0, 16.0),
            (1400.0, 19.0),
            (1700.0, 22.0),
        ];
        let mut scans: Vec<Scan> = Vec::new();
        for (i, &(mass, rt0)) in specs.iter().enumerate() {
            scans.extend(elution_scans(mass, 2, 15, 7, 0.15, rt0, (i as i32) * 100 + 1));
        }
        scans
    }

    #[test]
    fn tile2d_is_deterministic_and_matches_serial_on_clean_data() {
        let scans = tiling_fixture();
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let params = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.15,
            rt_half_window_minutes: 0.3,
            // Trace cap wide enough to claim each ~1.4-min elution whole (no fragmentation), so every
            // envelope is exactly one feature and there is no fragment to straddle a tile boundary.
            trace_max_half_width_minutes: 1.5,
            ..TraceKernelParameters::default()
        };

        let tiled = detect_features_tile2d(&engine, &params)
            .expect("fixture is wide/long enough for a ≥ 2×2 tile grid");
        // Deterministic: a second run reproduces the first exactly (order included).
        let tiled2 = detect_features_tile2d(&engine, &params).expect("second tiled run");
        assert_eq!(tiled.len(), tiled2.len(), "tiled run must be deterministic in count");
        for (a, b) in tiled.iter().zip(tiled2.iter()) {
            assert_eq!(a.charge, b.charge);
            approx(a.monoisotopic_mass, b.monoisotopic_mass, 1e-12);
            approx(a.apex_rt, b.apex_rt, 1e-12);
        }

        // Equivalent to serial on this clean, well-separated data: same feature set (charge, mono).
        let serial = detect_features_serial(&engine, &params);
        assert_eq!(
            tiled.len(),
            serial.len(),
            "tiled feature count {} must match serial {} on clean data",
            tiled.len(),
            serial.len()
        );
        let key = |f: &DetectedFeature| (f.charge, (f.monoisotopic_mass * 100.0).round() as i64);
        let mut tk: Vec<_> = tiled.iter().map(key).collect();
        let mut sk: Vec<_> = serial.iter().map(key).collect();
        tk.sort();
        sk.sort();
        assert_eq!(tk, sk, "tiled and serial must agree on the (charge, mono) multiset");
        // Sanity: at least the eight seeded envelopes are recovered.
        assert!(
            serial.len() >= 8,
            "the fixture holds at least eight detectable envelopes, got {}",
            serial.len()
        );
    }

    #[test]
    fn knee_default_params_do_not_change_detection() {
        // Default params have the knee disabled: detection on the synthetic z=2 fixture must be
        // identical whether the (present-but-disabled) knee fields hold their defaults or arbitrary
        // values — proving the added fields are byte-identical no-ops on the default path.
        let (scans, _mono_mz) = synthetic_envelope_scans();
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        let p_default = TraceKernelParameters {
            ppm_tolerance: 5.0,
            rt_sigma_minutes: 0.15,
            half_window_scans: 4,
            ..TraceKernelParameters::default()
        };
        let p_knee_off = TraceKernelParameters {
            knee_window_seeds: 1,
            knee_slope_frac: 0.9,
            knee_abs_eps: 1.0,
            ..p_default
        };
        let a = detect_features(&engine, &p_default);
        let b = detect_features(&engine, &p_knee_off);
        assert_eq!(a.len(), b.len(), "disabled knee must not change feature count");
        for (fa, fb) in a.iter().zip(b.iter()) {
            approx(fa.monoisotopic_mass, fb.monoisotopic_mass, 1e-9);
            assert_eq!(fa.charge, fb.charge);
        }
        // Sanity: the default run still detects the known charge-2, mass-1000 envelope.
        assert!(!a.is_empty());
        assert_eq!(a[0].charge, 2);
        approx(a[0].monoisotopic_mass, 1000.0, 0.01);
    }
}
