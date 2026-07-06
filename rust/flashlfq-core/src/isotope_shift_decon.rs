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
//! tooth — at the isotope index that the most-intense mass represents. Placing that mode tooth at the
//! anchor pins the whole comb; the monoisotope then falls out at `anchor − mode·spacing`. The `-1 / 0
//! / +1` search perturbs that placement by one ¹³C and keeps the shift whose comb best correlates
//! with the observed intensities, which is what corrects a mis-placed monoisotope.
//!
//! This is **new algorithm work** (an untargeted adaptation, not a line-by-line port), so it carries
//! its own unit tests rather than a C# golden. The Pearson correlation it scores with is the same
//! parity-ported [`crate::isotopic_envelope::pearson`] FlashLFQ uses.

use crate::deconvolution::averagine_comb_weights;
use crate::isotopic_envelope::{pearson, C13_MINUS_C12, PROTON_MASS};

/// Relative weight below which the averagine template's descending high-mass tail is dropped
/// (passed to [`averagine_comb_weights`]). Teeth up to and including the mode are always kept.
const TEMPLATE_MIN_WEIGHT: f64 = 1e-3;

/// Hard cap on averagine template length (isotope teeth) requested from the model.
const TEMPLATE_MAX_ISOTOPES: usize = 20;

/// The ¹³C shift hypotheses tested, in ascending order: the monoisotope is one ¹³C too low (`-1`),
/// correctly placed (`0`), or one ¹³C too high (`+1`). Index into
/// [`ShiftDeconResult::shift_correlations`] is `shift + 1`.
pub const SHIFTS: [i32; 3] = [-1, 0, 1];

/// Result of an untargeted isotope-shift deconvolution around one anchor peak.
#[derive(Debug, Clone, PartialEq)]
pub struct ShiftDeconResult {
    /// Monoisotopic neutral mass implied by the **winning** shift (`anchor − (mode − best_shift)·
    /// spacing`, converted to neutral mass).
    pub monoisotopic_mass: f64,
    /// Charge state searched (as supplied by the caller).
    pub charge: i32,
    /// The ¹³C shift (`-1`, `0`, or `+1`) whose comb best correlated with the observed intensities.
    pub best_shift: i32,
    /// Per-shift Pearson correlation of the averagine template to the observed intensities, indexed
    /// `[shift + 1]` (so `[0]` = shift −1, `[1]` = shift 0, `[2]` = shift +1). `NaN` correlations
    /// (degenerate/constant windows) are mapped to `-1.0`, matching FlashLFQ's gate handling.
    pub shift_correlations: [f64; 3],
    /// Whether the accurate (shift-0) hypothesis passes FlashLFQ's acceptance gate: its correlation
    /// exceeds `0.7` and neither shifted neighbour beats it by `0.1` or more. When `false`, shift 0
    /// is *not* confidently the monoisotope placement (a neighbour explains the data as well or
    /// better) — the signal the four-way disagreement router will key on.
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

/// FlashLFQ's acceptance gate thresholds (`CheckIsotopicEnvelopeCorrelation`).
const GATE_MIN_CORRELATION: f64 = 0.7;
const GATE_NEIGHBOUR_MARGIN: f64 = 0.1;

/// `mz.ToMass(charge)` in all-`f64` (the deconvolution/`ClassExtensions` convention):
/// `|charge|·mz − charge·ProtonMass`.
#[inline]
fn mz_to_mass(mz: f64, charge: i32) -> f64 {
    charge.abs() as f64 * mz - charge as f64 * PROTON_MASS
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
    // Averagine template keyed by the *most-intense* mass, indexed from the monoisotope (k = 0),
    // normalized so the max (mode) tooth is 1.0.
    let template = averagine_comb_weights(most_intense_mass, TEMPLATE_MIN_WEIGHT, TEMPLATE_MAX_ISOTOPES);
    if template.is_empty() {
        return None;
    }
    let mode_index = argmax(&template);
    let spacing = C13_MINUS_C12 / charge.abs() as f64;

    // Build the correlation vectors for each shift. The theoretical vector prepends a below-mono
    // padding tooth (theoretical intensity 0) so a hypothesis is penalised when there is unexpected
    // signal just below where it places the monoisotope — FlashLFQ's off-by-one discriminator.
    let mut correlations = [-1.0f64; 3];
    let mut matched_counts = [0usize; 3];
    for (si, &shift) in SHIFTS.iter().enumerate() {
        // Theoretical (T) and observed (O) over teeth k = -1 (pad), 0, 1, ..., template.len()-1.
        let mut theor: Vec<f64> = Vec::with_capacity(template.len() + 1);
        let mut obs: Vec<f64> = Vec::with_capacity(template.len() + 1);
        let mut matched = 0usize;
        for k in -1..(template.len() as i32) {
            let t = if k < 0 { 0.0 } else { template[k as usize] };
            // Predicted m/z of tooth k under this shift: the mode tooth sits at the anchor for
            // shift 0, and the whole comb slides by `shift` ¹³C units.
            let predicted_mz = anchor_mz + (k - mode_index as i32 + shift) as f64 * spacing;
            let o = nearest_within_ppm(mz, intensity, predicted_mz, tol_ppm).unwrap_or(0.0);
            if k >= 0 && o > 0.0 {
                matched += 1;
            }
            theor.push(t);
            obs.push(o);
        }
        let mut corr = pearson(&theor, &obs);
        if corr.is_nan() {
            corr = -1.0;
        }
        correlations[si] = corr;
        matched_counts[si] = matched;
    }

    // Winning shift = the best-correlating hypothesis (ties resolve to the lower index, i.e. the
    // smaller shift, which prefers the correctly-placed or lower monoisotope over a higher guess).
    let best_si = argmax(&correlations);
    let best_shift = SHIFTS[best_si];

    // FlashLFQ acceptance gate for the accurate (shift-0) hypothesis.
    let corr0 = correlations[1];
    let corr_left = correlations[0];
    let corr_right = correlations[2];
    let shift0_passes_gate = corr0 > GATE_MIN_CORRELATION
        && (corr_left - corr0) < GATE_NEIGHBOUR_MARGIN
        && (corr_right - corr0) < GATE_NEIGHBOUR_MARGIN;

    // Monoisotope for the winning shift: tooth k = 0 sits at anchor + (best_shift - mode)·spacing.
    let mono_mz = anchor_mz + (best_shift - mode_index as i32) as f64 * spacing;
    let monoisotopic_mass = mz_to_mass(mono_mz, charge);

    Some(ShiftDeconResult {
        monoisotopic_mass,
        charge,
        best_shift,
        shift_correlations: correlations,
        shift0_passes_gate,
        most_intense_mass,
        mode_index,
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
