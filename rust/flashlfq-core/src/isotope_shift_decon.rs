//! Untargeted isotope-shift deconvolution — the **FlashLFQ** monoisotope-placement strategy adapted
//! to de-novo (no-identification) feature detection.
//!
//! ## What this is
//! Base FlashLFQ ([`crate::isotopic_envelope::get_isotopic_envelopes`] /
//! `CheckIsotopicEnvelopeCorrelation`) resolves the monoisotope off-by-one by comparing the
//! **theoretical isotopic distribution** to the experimental data at ¹³C shifts of `-1, 0, +1` and
//! rejecting an envelope that a shifted neighbour explains better. That path is *targeted*: it is
//! handed the identified peptide's exact isotope distribution, its monoisotopic + peakfinding masses,
//! and the m/z index, so it never has to guess where the monoisotope is.
//!
//! In the untargeted pipeline there is no identification, so this module supplies the two missing
//! pieces from the **averagine** model and a plain spectrum slice:
//! - the theoretical distribution is [`crate::deconvolution::averagine_comb_weights`], keyed by the
//!   observed **most-intense mass** (see below), and
//! - peak lookup is a nearest-peak search over an ascending `(mz, intensity)` slice rather than the
//!   indexing engine.
//!
//! ## Anchor on the most-intense peak, not the monoisotope (deliberate)
//! The averagine is looked up by the observed **most-abundant** (tallest) peak's neutral mass, not by
//! a monoisotope guess. The tallest peak is the highest-SNR, most accurately-centroided point in the
//! envelope and is the hardest to misassign; the monoisotope — especially for heavier peptides whose
//! envelope mode sits one or more ¹³C above it — is exactly the fragile quantity we are trying to
//! *recover*, so anchoring on it would be circular. This mirrors `get_most_intense_mass_index` in the
//! parity-gated classic deconvolution and FlashLFQ's use of the observed peak as the reference.
//!
//! The averagine template (indexed from the monoisotope, `k = 0`) has its **mode** — its tallest
//! tooth — at the isotope index that the most-intense mass represents. That mode index locates the
//! monoisotope relative to the anchor (`anchor − mode·spacing` under the accurate hypothesis). The
//! shift search then perturbs that placement by whole ¹³C units and keeps the shift whose envelope
//! best explains the observed intensities.
//!
//! ## Per-shift template + shared-grid scoring (removes the shift-0 bias)
//! Each shift hypothesis is scored with its **own** averagine template, rebuilt (mono-keyed) at that
//! shift's candidate monoisotope, and every shift is correlated against **one shared observed grid**
//! spanning all the candidate tooth positions. This matters: an earlier version pinned a single
//! template's mode onto the (tallest) anchor and slid only the observed sampling, which handed shift
//! 0 a guaranteed max-on-max covariance term and biased the argmax toward shift 0 regardless of the
//! rest of the envelope. Scoring on a shared grid penalises a shift for observed peaks it fails to
//! predict, so the winner is the shift that explains the *whole* envelope, not just the tallest peak.
//!
//! ## Charge-dependent shift range
//! Light peptides only need `-1 / 0 / +1`. Heavier peptides — which carry higher charge — spread
//! their envelope over more isotopes and can be seeded up to two ¹³C units off, so at
//! `|charge| >= WIDE_SHIFT_MIN_CHARGE` the search widens to `-2 .. +2` (see [`shifts_for_charge`]).
//!
//! This is **new algorithm work** (an untargeted adaptation, not a line-by-line port), so it carries
//! its own unit tests rather than a C# golden. Shifts are scored by **cosine similarity** ([`cosine`])
//! rather than FlashLFQ's Pearson: on sparse non-negative envelope vectors Pearson mean-centres, which
//! rewards co-absent teeth (shared zeros read as agreement) and blurs the miss penalty, whereas cosine
//! treats a zero as no contribution so a predicted-but-absent tooth and an unexplained observed peak
//! each lower the score — the standard MS spectral-angle behaviour.

use crate::deconvolution::{averagine_comb_weights, averagine_intensities_from_mono};
use crate::isotopic_envelope::{C13_MINUS_C12, PROTON_MASS};

/// Relative weight below which the averagine template's descending high-mass tail is dropped
/// (passed to [`averagine_comb_weights`]). Teeth up to and including the mode are always kept.
const TEMPLATE_MIN_WEIGHT: f64 = 1e-3;

/// Hard cap on averagine template length (isotope teeth) requested from the model.
const TEMPLATE_MAX_ISOTOPES: usize = 20;

/// The default ¹³C shift hypotheses (ascending): the monoisotope is one ¹³C too low (`-1`), correctly
/// placed (`0`), or one ¹³C too high (`+1`).
pub const BASE_SHIFTS: [i32; 3] = [-1, 0, 1];

/// The wider shift set used at high charge, where a heavier peptide's envelope is spread over more
/// isotopes and the seed can be misplaced by up to two ¹³C units.
pub const WIDE_SHIFTS: [i32; 5] = [-2, -1, 0, 1, 2];

/// Minimum `|charge|` at which the wider [`WIDE_SHIFTS`] (`-2..+2`) set is searched instead of the
/// default [`BASE_SHIFTS`] (`-1..+1`). Charges of this magnitude and above correspond to the heavier
/// peptides whose off-by-one can span two isotopes.
pub const WIDE_SHIFT_MIN_CHARGE: i32 = 4;

/// Number of ¹³C units to walk the monoisotope back when double-checking a high-charge placement
/// (see [`walkback_mono_high_charge`]). Two covers the common heavy-peptide off-by-{one,two}-too-high.
pub const WALKBACK_MAX_C13: i32 = 2;

/// A lower-mass monoisotope hypothesis is adopted during the walk-back check if its envelope fit is
/// within this margin of the best-fitting candidate. Small and positive: the walk-back should only
/// prefer the lower mono when it explains the data *at least as well*, not merely differently.
pub const WALKBACK_FIT_MARGIN: f64 = 0.01;

/// Monoisotopic-mass floor (Da) for the walk-back check. Mass is the physical driver of the
/// off-by-one-too-high error — only heavy peptides carry an envelope mode well above the mono — so the
/// gate keys on mass, with [`WALKBACK_MIN_CHARGE`] as a companion floor.
pub const WALKBACK_MIN_MASS_DA: f64 = 2200.0;

/// Minimum `|charge|` for the walk-back check, paired with [`WALKBACK_MIN_MASS_DA`]. A heavy peptide's
/// low-charge (z=1/z=2) observations are sparse and their envelopes are best-resolved at z>=3; walking
/// those back added no recall and risked the metric's below-mono blind spot.
pub const WALKBACK_MIN_CHARGE: i32 = 3;


/// The shift hypotheses to search for a given charge: [`WIDE_SHIFTS`] once `|charge|` reaches
/// [`WIDE_SHIFT_MIN_CHARGE`], else [`BASE_SHIFTS`].
pub fn shifts_for_charge(charge: i32) -> &'static [i32] {
    if charge.abs() >= WIDE_SHIFT_MIN_CHARGE {
        &WIDE_SHIFTS
    } else {
        &BASE_SHIFTS
    }
}

/// Result of an untargeted isotope-shift deconvolution around one anchor peak.
#[derive(Debug, Clone, PartialEq)]
pub struct ShiftDeconResult {
    /// Monoisotopic neutral mass implied by the **winning** shift (`anchor − (mode − best_shift)·
    /// spacing`, converted to neutral mass).
    pub monoisotopic_mass: f64,
    /// Charge state searched (as supplied by the caller).
    pub charge: i32,
    /// The ¹³C shift whose rebuilt template best correlated with the observed intensities (drawn from
    /// [`shifts`](Self::shifts)). Ties resolve to the lower (smaller) shift.
    pub best_shift: i32,
    /// The ¹³C shift hypotheses evaluated, ascending — [`BASE_SHIFTS`] (`-1..+1`) normally, or
    /// [`WIDE_SHIFTS`] (`-2..+2`) at `|charge| >= WIDE_SHIFT_MIN_CHARGE`. Parallel to
    /// [`shift_correlations`](Self::shift_correlations).
    pub shifts: Vec<i32>,
    /// Per-shift **cosine similarity** ([`cosine`]) of that shift's rebuilt averagine template to the
    /// observed intensities, scored on a grid shared by all shifts (so a shift is penalised for
    /// observed peaks it fails to predict). Parallel to [`shifts`](Self::shifts). Degenerate windows
    /// (zero-norm observed vector) are mapped to `-1.0` so they rank last. Use
    /// [`shift0_correlation`](Self::shift0_correlation) to read the accurate-hypothesis value.
    pub shift_correlations: Vec<f64>,
    /// Whether the accurate (shift-0) hypothesis passes the acceptance gate: its cosine similarity
    /// exceeds [`GATE_MIN_CORRELATION`] and neither ±1 neighbour beats it by [`GATE_NEIGHBOUR_MARGIN`]
    /// or more. When `false`, shift 0 is *not* confidently the monoisotope placement (a neighbour
    /// explains the data as well or better) — the signal the four-way disagreement router keys on.
    pub shift0_passes_gate: bool,
    /// Neutral mass of the anchor (observed most-abundant) peak, `to_mass(anchor_mz, charge)`.
    pub most_intense_mass: f64,
    /// Isotope index (relative to the monoisotope) of the averagine template's mode — i.e. how many
    /// ¹³C units above the monoisotope the anchor sits under the shift-0 hypothesis.
    pub mode_index: usize,
    /// Number of observed teeth that matched within tolerance for the winning shift (excluding the
    /// below-mono padding tooth).
    pub matched_isotopes: usize,
}

impl ShiftDeconResult {
    /// Correlation of the accurate (shift-0) hypothesis. Shift 0 is always evaluated, so this is
    /// always present; returns `-1.0` only in the impossible case that it is absent.
    pub fn shift0_correlation(&self) -> f64 {
        self.shifts
            .iter()
            .position(|&s| s == 0)
            .map(|i| self.shift_correlations[i])
            .unwrap_or(-1.0)
    }
}

/// Acceptance-gate thresholds, structurally from FlashLFQ's `CheckIsotopicEnvelopeCorrelation`
/// (min-correlation + neighbour-margin). NOTE: these values were calibrated for FlashLFQ's *Pearson*
/// correlation; this module now scores with **cosine**, which on non-negative vectors runs higher,
/// so the `0.7` floor is looser than under Pearson. They are reasonable starting points but should be
/// re-tuned against the `obo_envelope_probe` scorecard.
const GATE_MIN_CORRELATION: f64 = 0.7;
const GATE_NEIGHBOUR_MARGIN: f64 = 0.1;

/// `mz.ToMass(charge)` in all-`f64` (the deconvolution/`ClassExtensions` convention):
/// `|charge|·mz − charge·ProtonMass`.
#[inline]
fn mz_to_mass(mz: f64, charge: i32) -> f64 {
    charge.abs() as f64 * mz - charge as f64 * PROTON_MASS
}

/// Uncentered cosine similarity of two equal-length non-negative vectors: `⟨a,b⟩ / (‖a‖·‖b‖)`.
///
/// Preferred over [`crate::isotopic_envelope::pearson`] for scoring an averagine template against
/// observed isotope intensities. Pearson mean-centres, which on sparse non-negative envelopes turns
/// co-absent teeth into positive covariance (shared zeros read as agreement) and makes every miss a
/// negative deviation that blurs the penalty. Cosine treats a zero as no contribution, so a
/// predicted-but-absent tooth (adds to `‖a‖` only) and an unexplained observed peak (adds to `‖b‖`
/// only) each lower the score. Returns `-1.0` when either vector has zero norm (degenerate), so the
/// shift argmax and the gate rank it last.
fn cosine(a: &[f64], b: &[f64]) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    let mut dot = 0.0f64;
    let mut norm_a = 0.0f64;
    let mut norm_b = 0.0f64;
    for (&x, &y) in a.iter().zip(b.iter()) {
        dot += x * y;
        norm_a += x * x;
        norm_b += y * y;
    }
    if norm_a <= 0.0 || norm_b <= 0.0 {
        return -1.0;
    }
    dot / (norm_a.sqrt() * norm_b.sqrt())
}

/// Index of the maximum element (first on ties). Empty slice yields 0 (callers guard emptiness).
fn argmax(v: &[f64]) -> usize {
    v.iter()
        .enumerate()
        .fold(0usize, |best, (i, &x)| if x > v[best] { i } else { best })
}

/// Nearest peak to `target_mz` in the ascending `mz` slice; returns its intensity iff it is within
/// `tol_ppm`, else `None`. `mz`/`intensity` are parallel and `mz` is ascending.
fn nearest_within_ppm(mz: &[f64], intensity: &[f64], target_mz: f64, tol_ppm: f64) -> Option<f64> {
    if mz.is_empty() {
        return None;
    }
    // First element >= target, then compare it and its predecessor.
    let ip = mz.partition_point(|&m| m < target_mz);
    let mut best_idx = ip.min(mz.len() - 1);
    if ip > 0 {
        let prev = ip - 1;
        if (target_mz - mz[prev]).abs() <= (mz[best_idx] - target_mz).abs() {
            best_idx = prev;
        }
    }
    let ppm = (mz[best_idx] - target_mz).abs() / target_mz * 1e6;
    if ppm <= tol_ppm {
        Some(intensity[best_idx])
    } else {
        None
    }
}

/// Deconvolutes the isotope envelope around an anchor (observed most-abundant) peak by the FlashLFQ
/// `-1 / 0 / +1` shift comparison, using the **averagine** template keyed by the anchor's mass.
///
/// - `mz` / `intensity`: an ascending-by-m/z, parallel spectrum slice covering the envelope's
///   vicinity (e.g. an averaged composite or a single apex scan, sliced to the feature window).
/// - `anchor_mz`: the m/z of the observed most-abundant peak of the envelope (the detector's seed).
/// - `charge`: the charge state to interpret the spacing at.
/// - `tol_ppm`: ppm tolerance for matching a predicted tooth to an observed peak.
///
/// Returns `None` when the slice is empty, the averagine template is empty, or every shift
/// hypothesis is degenerate (no correlation could be computed). Otherwise returns the winning shift,
/// the implied monoisotopic mass, the per-shift correlations, and the FlashLFQ acceptance-gate
/// verdict for the accurate hypothesis.
pub fn shift_decon(
    mz: &[f64],
    intensity: &[f64],
    anchor_mz: f64,
    charge: i32,
    tol_ppm: f64,
) -> Option<ShiftDeconResult> {
    assert_eq!(mz.len(), intensity.len(), "mz and intensity must be parallel");
    if mz.is_empty() || charge == 0 {
        return None;
    }

    let most_intense_mass = mz_to_mass(anchor_mz, charge);
    // Anchor-keyed averagine (indexed from the monoisotope, k = 0) used ONLY to locate the mono
    // relative to the tallest (anchor) peak: its mode sits `mode_index` ¹³C units above the mono, so
    // the shift-0 monoisotope is `mode_index` units below the anchor.
    let base = averagine_comb_weights(most_intense_mass, TEMPLATE_MIN_WEIGHT, TEMPLATE_MAX_ISOTOPES);
    if base.is_empty() {
        return None;
    }
    let mode_index = argmax(&base) as i32;
    let spacing = C13_MINUS_C12 / charge.abs() as f64;

    let shifts = shifts_for_charge(charge);

    // For each shift hypothesis, rebuild a *mono-keyed* averagine template placed at that shift's
    // candidate monoisotope. `mono_offset` is the tooth-k=0 position as an integer ¹³C offset from
    // the anchor; the whole template spans `[mono_offset, mono_offset + len - 1]`, with a padding
    // tooth at `mono_offset - 1` (theoretical 0) to catch unexpected below-mono signal.
    let mut templates: Vec<Vec<f64>> = Vec::with_capacity(shifts.len());
    let mut mono_offsets: Vec<i32> = Vec::with_capacity(shifts.len());
    let mut j_lo = i32::MAX;
    let mut j_hi = i32::MIN;
    for &shift in shifts {
        let o_s = shift - mode_index;
        let cand_mono_mass = mz_to_mass(anchor_mz + o_s as f64 * spacing, charge);
        let t = averagine_intensities_from_mono(cand_mono_mass, TEMPLATE_MIN_WEIGHT, TEMPLATE_MAX_ISOTOPES);
        j_lo = j_lo.min(o_s - 1);
        j_hi = j_hi.max(o_s + t.len() as i32 - 1).max(o_s);
        templates.push(t);
        mono_offsets.push(o_s);
    }
    if j_hi < j_lo {
        return None;
    }

    // Observed intensities on a grid SHARED by every shift (integer ¹³C offsets j = j_lo..=j_hi from
    // the anchor). Because the observed vector is identical across shifts, each shift is scored on
    // the same evidence: a shift that leaves a strong observed peak unexplained (obs large, theor 0)
    // is penalised, and one that predicts a tooth where nothing is observed is penalised too. This
    // removes the old shift-0 bias, where pinning the template's mode onto the tallest anchor peak
    // handed shift 0 a free max-on-max covariance term regardless of the rest of the envelope.
    let grid_len = (j_hi - j_lo + 1) as usize;
    let mut obs = vec![0.0f64; grid_len];
    for (gi, j) in (j_lo..=j_hi).enumerate() {
        let target_mz = anchor_mz + j as f64 * spacing;
        obs[gi] = nearest_within_ppm(mz, intensity, target_mz, tol_ppm).unwrap_or(0.0);
    }

    let mut correlations = vec![-1.0f64; shifts.len()];
    let mut matched_counts = vec![0usize; shifts.len()];
    for si in 0..shifts.len() {
        let o_s = mono_offsets[si];
        let t = &templates[si];
        let mut theor = vec![0.0f64; grid_len];
        let mut matched = 0usize;
        for (gi, j) in (j_lo..=j_hi).enumerate() {
            let k = j - o_s;
            if k >= 0 && (k as usize) < t.len() {
                theor[gi] = t[k as usize];
                if obs[gi] > 0.0 {
                    matched += 1;
                }
            }
        }
        // Cosine (not Pearson): uncentered, so it penalises predicted-but-absent teeth and
        // unexplained observed peaks symmetrically without crediting shared zeros. Degenerate
        // (zero-norm) windows already return -1.0.
        correlations[si] = cosine(&theor, &obs);
        matched_counts[si] = matched;
    }

    // Winning shift = the best-correlating hypothesis (ties resolve to the lower index, i.e. the
    // smaller shift, which prefers the correctly-placed or lower monoisotope over a higher guess).
    let best_si = argmax(&correlations);
    let best_shift = shifts[best_si];

    // FlashLFQ acceptance gate for the accurate (shift-0) hypothesis, judged against its immediate
    // ±1 neighbours (the genuine off-by-one question, independent of the wider ±2 search).
    let corr_at = |want: i32| -> f64 {
        shifts
            .iter()
            .position(|&s| s == want)
            .map(|i| correlations[i])
            .unwrap_or(f64::NEG_INFINITY)
    };
    let corr0 = corr_at(0);
    let shift0_passes_gate = corr0 > GATE_MIN_CORRELATION
        && (corr_at(-1) - corr0) < GATE_NEIGHBOUR_MARGIN
        && (corr_at(1) - corr0) < GATE_NEIGHBOUR_MARGIN;

    // Monoisotope for the winning shift: tooth k = 0 sits at anchor + (best_shift - mode)·spacing.
    let mono_mz = anchor_mz + (best_shift - mode_index) as f64 * spacing;
    let monoisotopic_mass = mz_to_mass(mono_mz, charge);

    Some(ShiftDeconResult {
        monoisotopic_mass,
        charge,
        best_shift,
        shifts: shifts.to_vec(),
        shift_correlations: correlations,
        shift0_passes_gate,
        most_intense_mass,
        mode_index: mode_index as usize,
        matched_isotopes: matched_counts[best_si],
    })
}

/// Convenience wrapper that first locates the anchor as the most-intense peak within
/// `[window_min_mz, window_max_mz]`, then runs [`shift_decon`]. Returns `None` if the window holds no
/// peak. Useful when the caller has an m/z window (the feature's envelope span) but not a specific
/// seed peak.
pub fn shift_decon_in_window(
    mz: &[f64],
    intensity: &[f64],
    window_min_mz: f64,
    window_max_mz: f64,
    charge: i32,
    tol_ppm: f64,
) -> Option<ShiftDeconResult> {
    assert_eq!(mz.len(), intensity.len(), "mz and intensity must be parallel");
    let lo = mz.partition_point(|&m| m < window_min_mz);
    let hi = mz.partition_point(|&m| m <= window_max_mz);
    if lo >= hi {
        return None;
    }
    let mut anchor_idx = lo;
    for i in lo..hi {
        if intensity[i] > intensity[anchor_idx] {
            anchor_idx = i;
        }
    }
    shift_decon(mz, intensity, mz[anchor_idx], charge, tol_ppm)
}

/// Fraction of a feature's expected **significant** isotope teeth that are actually observed.
///
/// For a feature placed at monoisotope m/z `mono_mz` and charge `charge`, the averagine template
/// (keyed by the neutral mass) gives the expected relative intensity at each isotope index `k`. A
/// tooth counts as *expected-significant* when its weight is at least `min_rel` of the mode; it is
/// *observed* when some peak lies within `tol_ppm` of `mono_mz + k·(C13/charge)`. Returns
/// observed / expected (`1.0` when nothing is expected).
///
/// Diagnoses the **charge-harmonic** artifact: a real z=1 species mislabeled z=2 has peaks only at
/// every OTHER z=2 tooth (the 0.5-Th intervening positions are empty), so ~half its expected teeth
/// are missing and completeness falls well below 1. A genuine envelope observes nearly all of them.
pub fn isotope_completeness(
    mz: &[f64],
    intensity: &[f64],
    mono_mz: f64,
    charge: i32,
    tol_ppm: f64,
    min_rel: f64,
) -> f64 {
    if mz.is_empty() || charge == 0 {
        return 1.0;
    }
    let mono_mass = mz_to_mass(mono_mz, charge);
    let template =
        averagine_intensities_from_mono(mono_mass, TEMPLATE_MIN_WEIGHT, TEMPLATE_MAX_ISOTOPES);
    if template.is_empty() {
        return 1.0;
    }
    let mode_w = template[argmax(&template)];
    if mode_w <= 0.0 {
        return 1.0;
    }
    let spacing = C13_MINUS_C12 / charge.abs() as f64;
    let mut expected = 0usize;
    let mut observed = 0usize;
    for (k, &w) in template.iter().enumerate() {
        if w / mode_w < min_rel {
            continue;
        }
        expected += 1;
        let target = mono_mz + k as f64 * spacing;
        if nearest_within_ppm(mz, intensity, target, tol_ppm).is_some() {
            observed += 1;
        }
    }
    if expected == 0 {
        1.0
    } else {
        observed as f64 / expected as f64
    }
}

/// **Envelope-fit cosine** — the unified fit / explained / completeness metric.
///
/// Cosine similarity between the theoretical averagine envelope `t` and the observed intensities `o`,
/// evaluated over the **union grid** = {predicted isotope positions} ∪ {observed peaks in the window}:
/// `S = Σ tⱼoⱼ / (‖t‖·‖o‖)`. Each of three properties maps to one behaviour of this single formula:
/// a **missing** predicted tooth (`t` large, `o = 0`) stays in `‖t‖` but not the numerator → S drops
/// (completeness); an **unexplained** observed peak (`t = 0`, `o` large) inflates `‖o‖` with no
/// numerator → S drops (fraction-explained); matching ratios raise the numerator (goodness of fit).
/// `S²` is the fraction of observed signal energy explained by the best-scaled envelope.
///
/// The window is `[mono_mz − 0.5·spacing, mono_mz + (kmax+1)·spacing]` where `kmax` is the last
/// significant tooth (weight ≥ `min_rel` of the mode); observed peaks below `noise_floor` are ignored
/// (not counted as "unexplained"). Returns a value in `[0, 1]` (`0` on degenerate input).
///
/// A wider lower reach (2.5·spacing, to see the peak below a too-high monoisotope) was tried both
/// globally and scoped to the walk-back and regressed CA/Lumos recall in every form, so the window
/// stays tight.
pub fn envelope_fit_cosine(
    mz: &[f64],
    intensity: &[f64],
    mono_mz: f64,
    charge: i32,
    tol_ppm: f64,
    min_rel: f64,
    noise_floor: f64,
) -> f64 {
    envelope_fit_cosine_masked(
        mz, intensity, mono_mz, charge, tol_ppm, min_rel, noise_floor, &[], 0.0,
    )
}

/// ppm within which a window peak is considered to sit on an `excluded` (neighbour-owned) position.
pub const NEIGHBOR_MASK_PPM: f64 = 12.0;

/// [`envelope_fit_cosine`] with a **neighbour mask**: any observed window peak within `excl_ppm` of an
/// `excluded` m/z (ascending) is made **invisible** to this fit — it can neither be matched to one of
/// this envelope's predicted teeth nor counted as unexplained signal. `excluded` should be the isotope
/// m/z of co-eluting *neighbour* features that are **not** on this feature's own grid.
///
/// This lets a low-scoring feature in a crowded window be judged on the signal that is plausibly *its
/// own*: a competing charge that only "fits" by absorbing a neighbour's peaks (the z↔2z harmonic, whose
/// intervening teeth are neighbour-filled) loses its borrowed evidence, and the walk-back cannot anchor
/// the mono on a neighbour's peak. With `excluded` empty this is exactly [`envelope_fit_cosine`].
#[allow(clippy::too_many_arguments)]
pub fn envelope_fit_cosine_masked(
    mz: &[f64],
    intensity: &[f64],
    mono_mz: f64,
    charge: i32,
    tol_ppm: f64,
    min_rel: f64,
    noise_floor: f64,
    excluded: &[f64],
    excl_ppm: f64,
) -> f64 {
    if mz.is_empty() || charge == 0 {
        return 0.0;
    }
    let mono_mass = mz_to_mass(mono_mz, charge);
    let template =
        averagine_intensities_from_mono(mono_mass, TEMPLATE_MIN_WEIGHT, TEMPLATE_MAX_ISOTOPES);
    if template.is_empty() {
        return 0.0;
    }
    let mode_w = template[argmax(&template)];
    if mode_w <= 0.0 {
        return 0.0;
    }
    let spacing = C13_MINUS_C12 / charge.abs() as f64;
    let sig: Vec<(usize, f64)> = template
        .iter()
        .enumerate()
        .filter(|(_, &w)| w / mode_w >= min_rel)
        .map(|(k, &w)| (k, w))
        .collect();
    if sig.is_empty() {
        return 0.0;
    }
    let kmax = sig.iter().map(|(k, _)| *k).max().unwrap();
    let win_lo = mono_mz - 0.5 * spacing;
    let win_hi = mono_mz + (kmax as f64 + 1.0) * spacing;
    let lo = mz.partition_point(|&m| m < win_lo);
    let hi = mz.partition_point(|&m| m <= win_hi);

    // Neighbour-owned window peaks: invisible to this fit (masked from both matching and penalty).
    let masked: Vec<bool> = (lo..hi)
        .map(|j| !excluded.is_empty() && near_any(mz[j], excluded, excl_ppm))
        .collect();

    let mut tvec: Vec<f64> = Vec::with_capacity(sig.len() + (hi - lo));
    let mut ovec: Vec<f64> = Vec::with_capacity(sig.len() + (hi - lo));
    let mut matched = vec![false; hi.saturating_sub(lo)];

    // Predicted teeth: pull the nearest non-masked in-window peak within tolerance (else observed 0).
    for &(k, w) in &sig {
        let target = mono_mz + k as f64 * spacing;
        let mut best_j: Option<usize> = None;
        let mut best_d = f64::INFINITY;
        for j in lo..hi {
            if masked[j - lo] {
                continue;
            }
            let d = (mz[j] - target).abs();
            if d < best_d {
                best_d = d;
                best_j = Some(j);
            }
        }
        let o = match best_j {
            Some(j) if (mz[j] - target).abs() / target * 1e6 <= tol_ppm => {
                matched[j - lo] = true;
                intensity[j]
            }
            _ => 0.0,
        };
        tvec.push(w);
        ovec.push(o);
    }
    // Off-grid observed peaks in the window (above noise, not neighbour-owned) → unexplained penalty.
    for j in lo..hi {
        if !matched[j - lo] && !masked[j - lo] && intensity[j] > noise_floor {
            tvec.push(0.0);
            ovec.push(intensity[j]);
        }
    }

    let nt = tvec.iter().map(|x| x * x).sum::<f64>().sqrt();
    let no = ovec.iter().map(|x| x * x).sum::<f64>().sqrt();
    if nt <= 0.0 || no <= 0.0 {
        return 0.0;
    }
    let dot: f64 = tvec.iter().zip(ovec.iter()).map(|(t, o)| t * o).sum();
    dot / (nt * no)
}

/// Envelope-fit margin by which an alternative charge must beat the detector's own charge before
/// [`best_charge_by_fit`] overrides it. A charge harmonic (z ↔ 2z) can spuriously tie or slightly beat
/// the true charge in a crowded window — the lower charge's teeth are a subset of the higher charge's
/// grid, and interferents fill the intervening higher-charge teeth — so re-charging is gated on a
/// clear improvement. Calibrated so a genuine re-charge still passes (e.g. a real z=2 the detector
/// called z=1 fits far better at z=2) while a spurious z→2z doubling in dense signal does not.
pub const RECHARGE_PREFER_MARGIN: f64 = 0.05;

/// Minimum absolute envelope fit an alternative charge must reach to override the detector's own
/// charge in [`best_charge_by_fit`]. Blocks a spurious harmonic flip whose "better" fit is only a
/// mediocre cosine earned by absorbing interferents on a denser grid; a genuine re-charge clears this
/// comfortably. See [`RECHARGE_PREFER_MARGIN`].
pub const RECHARGE_MIN_FIT: f64 = 0.85;

/// Fraction of the anchor (most-abundant) peak intensity that an observed peak on an **intervening**
/// finer-grid isotope position may reach before [`best_charge_by_fit`] refuses to *halve* the
/// detector charge onto a coarser grid.
///
/// A candidate charge that **divides** the detector's charge (e.g. z=1 winning over the detector's
/// z=2) predicts only every `ratio`-th tooth of the finer (detector-charge) grid. If strong observed
/// signal sits on the *intervening* teeth it cannot predict, its higher cosine is an artefact of
/// ignoring real signal, not a better model — so the halving is rejected and the detector charge kept.
///
/// Traced on **DSCQGDSGGPVVCSGK** (mono 1608.6508, z=2, RT 10.776): its observed envelope is
/// **mono-tallest** (k0..k3 ≈ 3.5/2.15/1.54/0.5 e6), but the averagine for a 1608-Da peptide predicts a
/// **+1-tallest** shape, so the *true* z=2 template fits the mono-tallest data worse than the **z=1**
/// comb — whose averagine mode at 804 Da *is* the mono. The z=1 comb reaches cosine 0.87 by matching
/// only the even z=2 teeth (805.34, 806.34) and **skipping** the odd ones (805.84 ≈ 2.15e6, 806.84 ≈
/// 0.5e6), which it leaves unexplained. Those intervening peaks are 2.15e6/3.5e6 = 0.61 of the anchor —
/// far above this floor — so the halving is refused and z=2 kept. The existing `RECHARGE_PREFER_MARGIN`
/// (0.05) and `RECHARGE_MIN_FIT` (0.85) did not block it: z=1 genuinely cleared 0.85 and beat z=2 by the
/// margin. A **genuine** halving (a real z=1 the detector mislabelled z=2) has *empty* intervening
/// positions (fraction ≈ 0), so it clears this gate untouched — preserving the light-z2 recovery.
pub const RECHARGE_HALVING_MAX_INTERVENING: f64 = 0.30;

/// Largest observed peak sitting on an **intervening** finer-grid isotope position — a position on the
/// `z_hi` (detector-charge) grid that the coarser divisor charge `z_lo` does **not** predict — relative
/// to the anchor peak intensity. Requires `z_lo < z_hi` and `z_hi % z_lo == 0`.
///
/// The intervening positions are scanned from the `z_lo` monoisotope upward across the `z_lo`
/// significant envelope (`ratio = z_hi / z_lo` finer teeth per `z_lo` tooth; the `j % ratio == 0`
/// positions are the `z_lo` teeth themselves and are skipped). Peaks within `excl_ppm` of an
/// `excluded` (neighbour-owned) position, and peaks at/below `noise_floor`, are ignored — they are not
/// this envelope's own unexplained signal. Returns `0.0` when the anchor peak is absent or degenerate.
#[allow(clippy::too_many_arguments)]
fn max_intervening_fraction(
    mz: &[f64],
    intensity: &[f64],
    anchor_mz: f64,
    z_hi: i32,
    z_lo: i32,
    mono_mass: f64,
    tol_ppm: f64,
    noise_floor: f64,
    excluded: &[f64],
    excl_ppm: f64,
) -> f64 {
    if z_lo == 0 || z_hi == 0 || z_lo.abs() >= z_hi.abs() || z_hi.abs() % z_lo.abs() != 0 {
        return 0.0;
    }
    let anchor_int = match nearest_within_ppm(mz, intensity, anchor_mz, tol_ppm) {
        Some(v) if v > 0.0 => v,
        _ => return 0.0,
    };
    let template =
        averagine_intensities_from_mono(mono_mass, TEMPLATE_MIN_WEIGHT, TEMPLATE_MAX_ISOTOPES);
    if template.is_empty() {
        return 0.0;
    }
    let mode_w = template[argmax(&template)];
    if mode_w <= 0.0 {
        return 0.0;
    }
    // Last significant z_lo tooth (>= 0.2 of the mode) bounds how far the intervening scan reaches.
    let kmax_lo = template
        .iter()
        .enumerate()
        .filter(|(_, &w)| w / mode_w >= 0.2)
        .map(|(k, _)| k)
        .max()
        .unwrap_or(0) as i32;
    let ratio = z_hi.abs() / z_lo.abs();
    let mono_mz = mono_mass / z_lo.abs() as f64 + PROTON_MASS;
    let spacing_hi = C13_MINUS_C12 / z_hi.abs() as f64;
    let jmax = (kmax_lo + 1) * ratio;
    let mut max_frac = 0.0f64;
    for j in 1..=jmax {
        if j % ratio == 0 {
            continue; // a z_lo tooth — z_lo already predicts it, so it is not intervening
        }
        let target = mono_mz + j as f64 * spacing_hi;
        if !excluded.is_empty() && near_any(target, excluded, excl_ppm) {
            continue; // neighbour-owned position — not this envelope's own unexplained signal
        }
        if let Some(o) = nearest_within_ppm(mz, intensity, target, tol_ppm) {
            if o > noise_floor {
                max_frac = max_frac.max(o / anchor_int);
            }
        }
    }
    max_frac
}

/// Chooses the charge from `candidates` whose shift-placed envelope best fits the slice by
/// [`envelope_fit_cosine`]. Each candidate charge is placed with [`shift_decon`] anchored on
/// `anchor_mz`, then scored. Returns `(charge, monoisotopic_mass, cosine)` of the best, or `None` if
/// no candidate yields an envelope.
///
/// `prefer_charge` (the detector's own charge) is kept unless another candidate's fit exceeds it by
/// `prefer_margin` — loyalty to the detector's charge that blocks a spurious harmonic flip (see
/// [`RECHARGE_PREFER_MARGIN`]) while still allowing a clearly-better re-charge.
#[allow(clippy::too_many_arguments)]
pub fn best_charge_by_fit(
    mz: &[f64],
    intensity: &[f64],
    anchor_mz: f64,
    candidates: &[i32],
    tol_ppm: f64,
    noise_floor: f64,
    prefer_charge: i32,
    prefer_margin: f64,
    excluded: &[f64],
    excl_ppm: f64,
) -> Option<(i32, f64, f64)> {
    // (charge, monoisotopic_mass, cosine) for every candidate that yields an envelope.
    let mut evals: Vec<(i32, f64, f64)> = Vec::with_capacity(candidates.len());
    for &z in candidates {
        if z == 0 {
            continue;
        }
        let Some(r) = shift_decon(mz, intensity, anchor_mz, z, tol_ppm) else {
            continue;
        };
        let mono_mz = r.monoisotopic_mass / z.abs() as f64 + PROTON_MASS;
        let cos = envelope_fit_cosine_masked(
            mz, intensity, mono_mz, z, tol_ppm, 0.2, noise_floor, excluded, excl_ppm,
        );
        evals.push((z, r.monoisotopic_mass, cos));
    }
    // Best by fit (ties keep the earlier/lower candidate, since candidates are ascending).
    let best = evals
        .iter()
        .copied()
        .reduce(|a, b| if b.2 > a.2 { b } else { a })?;
    // Detector-charge loyalty: keep prefer_charge unless another candidate BOTH beats it by the
    // margin AND itself fits well in absolute terms. The absolute floor is what blocks a spurious
    // harmonic (z↔2z) win in a crowded window: there the higher charge only "fits better" by
    // absorbing interferents on its denser grid, reaching a mediocre cosine (~0.69 on HVGDLGNVTADK's
    // sodium adduct), whereas a genuine re-charge fits cleanly (~0.98 on the light z=2 class). Without
    // the floor, a large but low-quality margin flips the correct charge.
    if let Some(pref) = evals.iter().copied().find(|e| e.0 == prefer_charge) {
        let mut keep_pref = best.2 <= pref.2 + prefer_margin || best.2 < RECHARGE_MIN_FIT;
        // Anti-halving guard: if the winner HALVES the detector charge (a proper divisor of
        // prefer_charge) yet strong observed signal sits on the intervening finer-grid teeth it
        // cannot predict, its higher cosine is an artefact of ignoring real signal — keep the
        // detector charge. This is the direction opposite the light-z2 doubling recovery (best > pref),
        // which is left untouched. Traced on DSCQGDSGGPVVCSGK — see RECHARGE_HALVING_MAX_INTERVENING.
        if !keep_pref && best.0.abs() < prefer_charge.abs() && prefer_charge.abs() % best.0.abs() == 0 {
            let frac = max_intervening_fraction(
                mz, intensity, anchor_mz, prefer_charge, best.0, best.1, tol_ppm, noise_floor,
                excluded, excl_ppm,
            );
            if frac > RECHARGE_HALVING_MAX_INTERVENING {
                keep_pref = true;
            }
        }
        if keep_pref {
            return Some(pref);
        }
    }
    Some(best)
}

/// Double-check a **heavy-peptide** monoisotope by walking it back up to [`WALKBACK_MAX_C13`] ¹³C units.
///
/// The default [`envelope_fit_cosine`] window begins at `mono − 0.5·spacing`, so a strong peak one or
/// two isotopes **below** the chosen mono is invisible to the fit. That is exactly the signature of an
/// off-by-one(-or-two)-**too-high** placement, which is common for heavy peptides: their envelope mode
/// sits above the monoisotope, so the detector (and even the shift search, if its anchor/slice differ
/// from the ideal) can seed the mono on a peak one or two ¹³C too high and the fit will not notice the
/// unexplained signal beneath it.
///
/// This re-scores the mono against its 1- and 2-¹³C-lower alternatives — whose windows **do** include
/// those lower teeth — and adopts the **lowest-mass** hypothesis whose fit is within
/// [`WALKBACK_FIT_MARGIN`] of the best. It is self-protecting: a correctly-placed mono's lower
/// alternatives predict a leading tooth where nothing is observed and score far worse, so a good mono
/// is never walked back.
///
/// Gated on **mass** ([`WALKBACK_MIN_MASS_DA`]) with a companion charge floor ([`WALKBACK_MIN_CHARGE`]):
/// mass is what drives the error (only heavy peptides carry an envelope mode well above the mono).
/// A pure `argmax(template) >= 1` mode gate (~1500 Da) and a lower 2000 Da all-charge floor were both
/// tried and regressed recall — at low charge / mid mass the fit metric's window cannot see below the
/// mono, so it scores the +1 placement above the true one, and the walk-back over-corrects.
#[allow(clippy::too_many_arguments)]
pub fn walkback_mono_high_charge(
    mz: &[f64],
    intensity: &[f64],
    mono_mass: f64,
    charge: i32,
    tol_ppm: f64,
    noise_floor: f64,
    excluded: &[f64],
    excl_ppm: f64,
) -> f64 {
    if charge.abs() < WALKBACK_MIN_CHARGE || mono_mass < WALKBACK_MIN_MASS_DA || mz.is_empty() {
        return mono_mass;
    }
    let z = charge.abs() as f64;
    let n = (WALKBACK_MAX_C13 + 1) as usize;
    let mut fits = vec![f64::NEG_INFINITY; n];
    let mut best_fit = f64::NEG_INFINITY;
    for b in 0..n {
        let cand_mass = mono_mass - b as f64 * C13_MINUS_C12;
        let cand_mz = cand_mass / z + PROTON_MASS;
        let cos = envelope_fit_cosine_masked(
            mz, intensity, cand_mz, charge, tol_ppm, 0.2, noise_floor, excluded, excl_ppm,
        );
        fits[b] = cos;
        if cos > best_fit {
            best_fit = cos;
        }
    }
    // Prefer the lowest-mass (largest walk-back) candidate that still fits within the margin of the
    // best — the off-by-one error is directionally "mono too high", never too low.
    for b in (0..n).rev() {
        if fits[b] >= best_fit - WALKBACK_FIT_MARGIN {
            return mono_mass - b as f64 * C13_MINUS_C12;
        }
    }
    mono_mass
}

/// Whether `mz` is within `ppm` of any position in the ascending `sorted` list.
fn near_any(mz: f64, sorted: &[f64], ppm: f64) -> bool {
    if sorted.is_empty() {
        return false;
    }
    let ip = sorted.partition_point(|&m| m < mz);
    for cand in [ip.checked_sub(1), Some(ip)].into_iter().flatten() {
        if let Some(&f) = sorted.get(cand) {
            if (mz - f).abs() / mz * 1e6 <= ppm {
                return true;
            }
        }
    }
    false
}

/// **Neighbor-aware** variant of [`shift_decon_in_window`]: the anchor (most-abundant) peak is the
/// tallest peak in the window that is *eligible*, where a peak is eligible unless it belongs to an
/// already-detected neighbouring feature. A peak is **excluded** when it lies within `grid_ppm` of a
/// `forbidden_mz` position (a neighbour's predicted isotope m/z) AND is **not** on this feature's own
/// isotope grid (`own_mono_mz + k·spacing`, `k ∈ [0, own_kmax]`). Own-grid peaks are always eligible,
/// so a co-eluting duplicate detection cannot exclude this feature's real teeth.
///
/// This directly fixes the "distant grab": the tallest peak in a wide window is often a stronger
/// co-eluting *different* species several isotopes away; anchoring there places the monoisotope on
/// the neighbour. Excluding the neighbour's grid positions makes the shift decon anchor on this
/// feature's own strongest peak instead. `forbidden_mz` must be ascending.
#[allow(clippy::too_many_arguments)]
pub fn shift_decon_gated(
    mz: &[f64],
    intensity: &[f64],
    window_min_mz: f64,
    window_max_mz: f64,
    charge: i32,
    tol_ppm: f64,
    own_mono_mz: f64,
    own_kmax: usize,
    forbidden_mz: &[f64],
    grid_ppm: f64,
) -> Option<ShiftDeconResult> {
    assert_eq!(mz.len(), intensity.len(), "mz and intensity must be parallel");
    if charge == 0 {
        return None;
    }
    let spacing = C13_MINUS_C12 / charge.abs() as f64;
    let lo = mz.partition_point(|&m| m < window_min_mz);
    let hi = mz.partition_point(|&m| m <= window_max_mz);

    let on_own_grid = |m: f64| -> bool {
        let k = ((m - own_mono_mz) / spacing).round();
        if k < 0.0 || k > own_kmax as f64 {
            return false;
        }
        let expected = own_mono_mz + k * spacing;
        (m - expected).abs() / expected * 1e6 <= grid_ppm
    };

    let mut anchor_idx: Option<usize> = None;
    let mut anchor_int = f64::NEG_INFINITY;
    for i in lo..hi {
        let m = mz[i];
        // Excluded iff attributed to a neighbour and not one of our own teeth.
        if !on_own_grid(m) && near_any(m, forbidden_mz, grid_ppm) {
            continue;
        }
        if intensity[i] > anchor_int {
            anchor_int = intensity[i];
            anchor_idx = Some(i);
        }
    }
    let ai = anchor_idx?;
    shift_decon(mz, intensity, mz[ai], charge, tol_ppm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::isotopic_envelope::mass_to_mz_f64;

    fn ppm_diff(a: f64, b: f64) -> f64 {
        (a - b).abs() / b * 1e6
    }

    /// Real apex peaks for `23_ECCHGDLLECADDRADLAK_z4` (true mono 2246.9355, z=4). The detector
    /// seeded this feature's mono on the **+2** isotope (563.247, mass 2248.96), one of the
    /// heavy-peptide off-by-two-too-high cases. Used by the walk-back regression tests below.
    fn ecc_z4_apex() -> (Vec<f64>, Vec<f64>) {
        let apex: &[(f64, f64)] = &[
            (562.24908, 2.723067e5), (562.35175, 1.424823e6), (562.50031, 3.830424e5),
            (562.74670, 7.995974e6), (562.76428, 2.456608e5), (562.97852, 2.544323e5),
            (562.99744, 8.917290e6), (563.24756, 5.830500e6), (563.28729, 2.374914e5),
            (563.29932, 2.094416e5), (563.36053, 9.135656e5), (563.49811, 3.341500e6),
            (563.74756, 1.148486e6), (563.76013, 2.952547e5), (563.99792, 6.053904e5),
            (564.10522, 3.223775e5), (564.24969, 2.600071e5), (564.28308, 2.871794e5),
            (564.36462, 2.759059e5), (564.50531, 2.771088e5),
        ];
        (apex.iter().map(|p| p.0).collect(), apex.iter().map(|p| p.1).collect())
    }

    #[test]
    fn walkback_recovers_true_mono_for_off_by_one_high_charge() {
        // Simulate refine having landed on the +1 mono (562.997 → mass 2247.9578). The fit there
        // (0.934) looks acceptable only because the strong 562.746 peak sits just below its window.
        let (mz, inten) = ecc_z4_apex();
        let wrong_mono = 2247.9578; // +1 too high
        let corrected = walkback_mono_high_charge(&mz, &inten, wrong_mono, 4, 20.0, 0.0, &[], 0.0);
        assert!(
            ppm_diff(corrected, 2246.9355) <= 30.0,
            "walk-back should recover the true mono 2246.9355, got {corrected}"
        );
    }

    #[test]
    fn walkback_leaves_correct_mono_untouched() {
        // A correctly-placed mono must not be walked back: its 1-/2-lower alternatives predict a
        // leading tooth where nothing is observed and fit far worse.
        let (mz, inten) = ecc_z4_apex();
        let true_mono = 2246.9544; // what shift_decon already recovers on clean data
        let kept = walkback_mono_high_charge(&mz, &inten, true_mono, 4, 20.0, 0.0, &[], 0.0);
        assert!(
            (kept - true_mono).abs() < 1e-9,
            "correct mono should be left untouched, got {kept}"
        );
    }

    #[test]
    fn walkback_noop_below_mass_floor() {
        // Mass gate: a light peptide (~900 Da, below WALKBACK_MIN_MASS_DA) is a no-op at any charge —
        // its monoisotope is the mode, so an off-by-one-too-high cannot hide from the fit.
        let (mz, inten) = ecc_z4_apex();
        let m = 900.0;
        assert_eq!(walkback_mono_high_charge(&mz, &inten, m, 4, 20.0, 0.0, &[], 0.0), m);
        assert_eq!(walkback_mono_high_charge(&mz, &inten, m, 2, 20.0, 0.0, &[], 0.0), m);
    }


    /// Builds a clean synthetic spectrum from the averagine template for a peptide of monoisotopic
    /// `mono_mass` at `charge`: one peak per template tooth, placed at the tooth's m/z with the
    /// template weight as intensity. Returns ascending `(mz, intensity)` and the anchor m/z (the
    /// mode/most-abundant tooth). This is the same construction the deconvolution tests use.
    fn synthetic_envelope(mono_mass: f64, charge: i32) -> (Vec<f64>, Vec<f64>, f64) {
        // Key the template by the mono's own most-intense mass so the fixture is self-consistent
        // with what shift_decon will look up from the anchor.
        let spacing = C13_MINUS_C12 / charge as f64;
        let mono_mz = mass_to_mz_f64(mono_mass, charge);
        // Approximate the template by looking it up from the mono (good enough to synthesize peaks;
        // shift_decon re-derives its own template from the anchor mass).
        let template = crate::deconvolution::averagine_intensities_from_mono(mono_mass, 1e-3, 20);
        let mut mz: Vec<f64> = Vec::new();
        let mut inten: Vec<f64> = Vec::new();
        for (k, &w) in template.iter().enumerate() {
            mz.push(mono_mz + k as f64 * spacing);
            inten.push(w * 1.0e7);
        }
        let mode = argmax(&inten);
        let anchor_mz = mz[mode];
        (mz, inten, anchor_mz)
    }

    #[test]
    fn best_charge_keeps_detector_charge_on_clean_envelope() {
        // A clean z=2 envelope: the detector's z=2 is both the best fit and the preferred charge.
        let (mz, inten, anchor) = synthetic_envelope(1246.0, 2);
        let (z, _, _) =
            best_charge_by_fit(&mz, &inten, anchor, &[1, 2, 4], 20.0, 0.0, 2, RECHARGE_PREFER_MARGIN, &[], 0.0)
                .expect("charge");
        assert_eq!(z, 2, "clean z=2 envelope should stay z=2");
    }

    #[test]
    fn best_charge_allows_genuine_recharge_to_cleanly_fitting_alt() {
        // Detector mislabeled a clean z=2 envelope as z=1. The z=2 fit is high (clears RECHARGE_MIN_FIT)
        // and clearly beats z=1, so the genuine re-charge fires despite detector-charge loyalty — this
        // is the light-z2 charge-halving recovery that must survive the loyalty margin + floor.
        let (mz, inten, anchor) = synthetic_envelope(1246.0, 2);
        let (z, _, cos) =
            best_charge_by_fit(&mz, &inten, anchor, &[1, 2], 20.0, 0.0, 1, RECHARGE_PREFER_MARGIN, &[], 0.0)
                .expect("charge");
        assert_eq!(z, 2, "a clean z=2 fit must override a wrong z=1 detector label");
        assert!(cos >= RECHARGE_MIN_FIT, "the winning z=2 fit {cos} should clear the floor");
    }

    #[test]
    fn intervening_fraction_flags_present_and_ignores_absent() {
        // z_hi=2 (detector) vs z_lo=1 (halving candidate). A z=1 envelope (mono is the mode) plus a
        // strong peak on the intervening z=2 (odd, half-Th) position should be flagged; a clean z=1
        // envelope with NO intervening peak (a genuine halving) must read ~0.
        let z1_mono = 804.3327;
        let mono_mz = mass_to_mz_f64(z1_mono, 1);
        let sp1 = C13_MINUS_C12;
        let sp2 = C13_MINUS_C12 / 2.0;
        // z1 teeth: anchor (mono) tallest.
        let base = vec![
            (mono_mz, 3.5e6),
            (mono_mz + sp1, 1.5e6),
            (mono_mz + 2.0 * sp1, 0.5e6),
        ];

        // Present: add an intervening peak at 0.6 of the anchor.
        let mut with = base.clone();
        with.push((mono_mz + sp2, 2.1e6));
        with.sort_by(|a, b| a.0.total_cmp(&b.0));
        let wm: Vec<f64> = with.iter().map(|p| p.0).collect();
        let wi: Vec<f64> = with.iter().map(|p| p.1).collect();
        let f = max_intervening_fraction(&wm, &wi, mono_mz, 2, 1, z1_mono, 20.0, 0.0, &[], 0.0);
        assert!(
            f > RECHARGE_HALVING_MAX_INTERVENING,
            "present intervening peak should exceed the floor, got {f}"
        );

        // Absent: a genuine z=1 envelope has empty intervening positions → fraction ~0.
        let am: Vec<f64> = base.iter().map(|p| p.0).collect();
        let ai: Vec<f64> = base.iter().map(|p| p.1).collect();
        let f0 = max_intervening_fraction(&am, &ai, mono_mz, 2, 1, z1_mono, 20.0, 0.0, &[], 0.0);
        assert!(
            f0 <= RECHARGE_HALVING_MAX_INTERVENING,
            "no intervening peak → fraction ~0, got {f0}"
        );
    }

    #[test]
    fn best_charge_rejects_spurious_halving_with_intervening_signal() {
        // DSCQGDSGGPVVCSGK-shaped case: detector z=2, observed envelope MONO-TALLEST with strong
        // signal on the intervening (odd) z=2 teeth. The z=1 comb fits its own (even) teeth well and
        // skips the intervening peaks, which without the anti-halving guard let it out-score z=2 and
        // halve the mass to ~804. The guard sees the intervening signal and keeps the detector's z=2.
        let z1_mono = 804.3327;
        let mono_mz = mass_to_mz_f64(z1_mono, 1); // == z2 mono m/z (805.34)
        let sp1 = C13_MINUS_C12;
        let sp2 = C13_MINUS_C12 / 2.0;
        // z1-shaped even teeth (so z=1 fits them cleanly), plus strong intervening odd teeth.
        let tmpl1 = crate::deconvolution::averagine_intensities_from_mono(z1_mono, 1e-3, 20);
        let mut peaks: Vec<(f64, f64)> = tmpl1
            .iter()
            .enumerate()
            .take(4)
            .map(|(k, &w)| (mono_mz + k as f64 * sp1, w * 1.0e7))
            .collect();
        let anchor_int = peaks[0].1;
        peaks.push((mono_mz + sp2, 0.6 * anchor_int)); // intervening between k0 and k1
        peaks.push((mono_mz + 3.0 * sp2, 0.25 * anchor_int)); // intervening between k1 and k2
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mz: Vec<f64> = peaks.iter().map(|p| p.0).collect();
        let inten: Vec<f64> = peaks.iter().map(|p| p.1).collect();

        let (z, _mass, _cos) = best_charge_by_fit(
            &mz, &inten, mono_mz, &[1, 2, 4], 20.0, 0.0, 2, RECHARGE_PREFER_MARGIN, &[], 0.0,
        )
        .expect("charge");
        assert_eq!(z, 2, "strong intervening signal must block the spurious halving to z=1");
    }

    #[test]
    fn clean_light_peptide_places_mono_at_shift_zero() {
        // Light peptide (~1200 Da): the monoisotope is the mode, so the anchor IS the mono and the
        // accurate hypothesis is shift 0.
        let mono = 1200.0;
        let charge = 2;
        let (mz, inten, anchor) = synthetic_envelope(mono, charge);
        let r = shift_decon(&mz, &inten, anchor, charge, 15.0).expect("should deconvolute");
        assert_eq!(r.best_shift, 0, "clean light envelope: shift 0 wins");
        assert!(r.shift0_passes_gate, "clean envelope should pass the gate");
        assert!(
            ppm_diff(r.monoisotopic_mass, mono) <= 20.0,
            "recovered mono {} not within 20 ppm of {} (ppm={})",
            r.monoisotopic_mass,
            mono,
            ppm_diff(r.monoisotopic_mass, mono)
        );
    }

    #[test]
    fn gated_anchor_ignores_taller_off_grid_neighbor() {
        // Reproduce the drift0 mechanism: a target envelope plus a TALLER, unrelated species ~4
        // isotopes above, off the target's grid. Ungated shift-decon grabs the neighbour → drifts;
        // the gated variant, told the neighbour's grid m/z, anchors on the target instead.
        let charge = 1;
        let target_mono = 980.0;
        let (mut mz, mut inten, target_anchor) = synthetic_envelope(target_mono, charge);
        let target_mono_mz = mass_to_mz_f64(target_mono, charge);
        let spacing = C13_MINUS_C12 / charge as f64;

        // Foreign species ~4 isotopes up (off the integer grid so it's clearly a different peptide),
        // ~1.5× the target's tallest intensity. Add its envelope teeth.
        let target_max = inten.iter().copied().fold(0.0_f64, f64::max);
        let foreign_mono_mz = target_mono_mz + 3.94 * spacing;
        let mut foreign_grid: Vec<f64> = Vec::new();
        for (k, w) in [1.0, 0.55, 0.2].into_iter().enumerate() {
            let m = foreign_mono_mz + k as f64 * spacing;
            mz.push(m);
            inten.push(target_max * 1.5 * w);
            foreign_grid.push(m);
        }
        // Re-sort ascending (parallel).
        let mut pairs: Vec<(f64, f64)> = mz.iter().copied().zip(inten.iter().copied()).collect();
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mz: Vec<f64> = pairs.iter().map(|p| p.0).collect();
        let inten: Vec<f64> = pairs.iter().map(|p| p.1).collect();
        foreign_grid.sort_by(f64::total_cmp);

        let win_lo = target_mono_mz - 1.5;
        let win_hi = target_mono_mz + 8.0 * spacing;

        // Ungated: grabs the taller foreign peak → mono drifts up by ~4 isotopes.
        let ungated = shift_decon_in_window(&mz, &inten, win_lo, win_hi, charge, 15.0).unwrap();
        let k_ungated = ((ungated.monoisotopic_mass - target_mono) / C13_MINUS_C12).round() as i32;
        assert!(k_ungated >= 3, "ungated should drift onto the neighbour (k={k_ungated})");

        // Gated with the neighbour's grid forbidden: anchors on the target → k=0.
        let gated = shift_decon_gated(
            &mz, &inten, win_lo, win_hi, charge, 15.0, target_mono_mz, 6, &foreign_grid, 15.0,
        )
        .unwrap();
        let k_gated = ((gated.monoisotopic_mass - target_mono) / C13_MINUS_C12).round() as i32;
        assert_eq!(k_gated, 0, "gated should recover the target mono (k={k_gated})");
        // Sanity: the gated anchor is the target's own most-abundant peak.
        assert!(ppm_diff(gated.most_intense_mass, mz_to_mass(target_anchor, charge)) < 5.0);
    }

    #[test]
    fn heavy_peptide_mode_above_mono_still_recovers_mono() {
        // Heavy peptide (~3200 Da): the averagine mode sits one or more ¹³C above the monoisotope,
        // so the most-abundant peak (anchor) is NOT the mono. shift_decon must still put the mono
        // `mode` teeth below the anchor.
        let mono = 3200.0;
        let charge = 3;
        let (mz, inten, anchor) = synthetic_envelope(mono, charge);
        let r = shift_decon(&mz, &inten, anchor, charge, 15.0).expect("should deconvolute");
        assert!(r.mode_index >= 1, "heavy peptide: anchor is above the mono (mode {})", r.mode_index);
        assert_eq!(r.best_shift, 0, "clean heavy envelope still fits at shift 0");
        assert!(
            ppm_diff(r.monoisotopic_mass, mono) <= 25.0,
            "recovered mono {} not within 25 ppm of {} (ppm={})",
            r.monoisotopic_mass,
            mono,
            ppm_diff(r.monoisotopic_mass, mono)
        );
    }

    #[test]
    fn unexpected_peak_below_mono_favours_the_lower_shift() {
        // Add a strong peak one ¹³C below the true monoisotope (a co-eluting lower species, or the
        // real mono the classic anchor missed). The below-mono padding tooth should let shift -1
        // correlate at least as well, so the accurate hypothesis no longer confidently passes.
        let mono = 1500.0;
        let charge = 2;
        let (mut mz, mut inten, anchor) = synthetic_envelope(mono, charge);
        let spacing = C13_MINUS_C12 / charge as f64;
        let below_mz = mass_to_mz_f64(mono, charge) - spacing;
        // Insert keeping ascending order.
        mz.insert(0, below_mz);
        inten.insert(0, 0.9e7);
        let r = shift_decon(&mz, &inten, anchor, charge, 15.0).expect("should deconvolute");
        // With strong unexpected intensity below the mono, shift 0 should not sail through the gate.
        assert!(
            !r.shift0_passes_gate || r.best_shift == -1,
            "unexpected below-mono peak should challenge the accurate hypothesis \
             (best_shift={}, gate={}, corr={:?})",
            r.best_shift,
            r.shift0_passes_gate,
            r.shift_correlations
        );
    }

    #[test]
    fn empty_slice_returns_none() {
        assert!(shift_decon(&[], &[], 500.0, 2, 15.0).is_none());
    }

    #[test]
    fn window_wrapper_finds_the_tallest_peak() {
        let mono = 1000.0;
        let charge = 2;
        let (mz, inten, anchor) = synthetic_envelope(mono, charge);
        let lo = mz[0] - 1.0;
        let hi = *mz.last().unwrap() + 1.0;
        let r = shift_decon_in_window(&mz, &inten, lo, hi, charge, 15.0).expect("window has peaks");
        // The wrapper's anchor is the tallest peak, which is the same one the direct call uses.
        let direct = shift_decon(&mz, &inten, anchor, charge, 15.0).unwrap();
        assert_eq!(r.monoisotopic_mass, direct.monoisotopic_mass);
        assert_eq!(r.best_shift, direct.best_shift);
    }
}
