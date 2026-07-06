//! Joint multi-envelope deconvolution — the "match multiple envelopes simultaneously" advanced decon.
//!
//! When a low-scoring feature sits in a chimeric window, the single-envelope [`crate::isotope_shift_decon`]
//! fit cannot separate its signal from co-eluting neighbours. This module instead models the window as a
//! **non-negative linear combination** of several averagine envelopes (the feature plus its overlapping
//! neighbours) and solves for their amplitudes by non-negative least squares (NNLS), optionally over a
//! small search of the target's monoisotope position. The winning placement is the one whose *joint*
//! model best explains the observed window — so a peak the target used to be scored against, but which
//! actually belongs to a neighbour, is now attributed to that neighbour rather than penalising (or
//! wrongly rewarding) the target.
//!
//! This is de-novo algorithm work (no C# golden); it carries its own unit tests.

use crate::deconvolution::averagine_intensities_from_mono;
use crate::isotopic_envelope::{C13_MINUS_C12, PROTON_MASS};

const TEMPLATE_MIN_WEIGHT: f64 = 1e-3;
const TEMPLATE_MAX_ISOTOPES: usize = 20;

/// One component of a joint fit: an averagine envelope at `(mono_mz, charge)`.
#[derive(Debug, Clone, Copy)]
pub struct Component {
    pub mono_mz: f64,
    pub charge: i32,
}

/// Result of a joint fit: per-component amplitudes (parallel to the input components) and the cosine
/// similarity of the summed model to the observed window.
#[derive(Debug, Clone)]
pub struct JointFitResult {
    pub amplitudes: Vec<f64>,
    pub fit: f64,
}

fn mz_to_mass(mz: f64, charge: i32) -> f64 {
    charge.abs() as f64 * mz - charge as f64 * PROTON_MASS
}

/// Significant averagine teeth of a component as `(mz, weight)`, weight normalised to a unit-max shape.
fn component_teeth(c: &Component, min_rel: f64) -> Vec<(f64, f64)> {
    if c.charge == 0 {
        return Vec::new();
    }
    let mass = mz_to_mass(c.mono_mz, c.charge);
    let t = averagine_intensities_from_mono(mass, TEMPLATE_MIN_WEIGHT, TEMPLATE_MAX_ISOTOPES);
    if t.is_empty() {
        return Vec::new();
    }
    let mode = t.iter().cloned().fold(0.0_f64, f64::max);
    if mode <= 0.0 {
        return Vec::new();
    }
    let spacing = C13_MINUS_C12 / c.charge.abs() as f64;
    t.iter()
        .enumerate()
        .filter(|(_, &w)| w / mode >= min_rel)
        .map(|(k, &w)| (c.mono_mz + k as f64 * spacing, w))
        .collect()
}

/// Solves `min_{x >= 0} ||A x - b||` by the Lawson–Hanson active-set NNLS. `a_cols[j]` is column `j`
/// (length = `b.len()`). Small dense systems only (a handful of components); the passive least-squares
/// subproblem is solved via the normal equations with Gaussian elimination.
pub fn nnls(a_cols: &[Vec<f64>], b: &[f64], max_iter: usize) -> Vec<f64> {
    let n = a_cols.len();
    let m = b.len();
    let mut x = vec![0.0f64; n];
    if n == 0 || m == 0 {
        return x;
    }
    let mut passive = vec![false; n];
    // w = Aᵀ (b - A x)
    let residual = |x: &[f64]| -> Vec<f64> {
        let mut r = b.to_vec();
        for (j, xj) in x.iter().enumerate() {
            if *xj != 0.0 {
                for i in 0..m {
                    r[i] -= a_cols[j][i] * xj;
                }
            }
        }
        r
    };
    let grad = |r: &[f64]| -> Vec<f64> {
        (0..n).map(|j| (0..m).map(|i| a_cols[j][i] * r[i]).sum()).collect()
    };

    let tol = 1e-9;
    for _outer in 0..max_iter.max(1) {
        let r = residual(&x);
        let w = grad(&r);
        // Pick the active variable with the largest gradient; stop if none is positive.
        let mut best_j = None;
        let mut best_w = tol;
        for j in 0..n {
            if !passive[j] && w[j] > best_w {
                best_w = w[j];
                best_j = Some(j);
            }
        }
        let Some(j_in) = best_j else { break };
        passive[j_in] = true;

        loop {
            // Least squares over the passive set: solve (A_Pᵀ A_P) s_P = A_Pᵀ b.
            let p_idx: Vec<usize> = (0..n).filter(|&j| passive[j]).collect();
            let s_p = solve_passive(a_cols, b, &p_idx, m);
            let mut s = vec![0.0f64; n];
            for (k, &j) in p_idx.iter().enumerate() {
                s[j] = s_p[k];
            }
            // If all passive components are positive, accept.
            if p_idx.iter().all(|&j| s[j] > 0.0) {
                x = s;
                break;
            }
            // Otherwise move x toward s until a passive component hits zero, then drop it.
            let mut alpha = f64::INFINITY;
            for &j in &p_idx {
                if s[j] <= 0.0 {
                    let denom = x[j] - s[j];
                    if denom > 0.0 {
                        alpha = alpha.min(x[j] / denom);
                    }
                }
            }
            if !alpha.is_finite() {
                x = s;
                break;
            }
            for j in 0..n {
                x[j] += alpha * (s[j] - x[j]);
            }
            for j in 0..n {
                if passive[j] && x[j] <= tol {
                    passive[j] = false;
                    x[j] = 0.0;
                }
            }
        }
    }
    x.iter_mut().for_each(|v| *v = v.max(0.0));
    x
}

/// Solves the normal equations `(A_Pᵀ A_P) s = A_Pᵀ b` for the passive columns `p_idx` (Gaussian
/// elimination with partial pivoting). Returns `s` parallel to `p_idx`.
fn solve_passive(a_cols: &[Vec<f64>], b: &[f64], p_idx: &[usize], m: usize) -> Vec<f64> {
    let k = p_idx.len();
    if k == 0 {
        return Vec::new();
    }
    // Gram matrix G = A_Pᵀ A_P (k×k) and rhs = A_Pᵀ b.
    let mut g = vec![vec![0.0f64; k]; k];
    let mut rhs = vec![0.0f64; k];
    for a in 0..k {
        let ca = &a_cols[p_idx[a]];
        rhs[a] = (0..m).map(|i| ca[i] * b[i]).sum();
        for bb in a..k {
            let cb = &a_cols[p_idx[bb]];
            let dot: f64 = (0..m).map(|i| ca[i] * cb[i]).sum();
            g[a][bb] = dot;
            g[bb][a] = dot;
        }
    }
    // Solve G s = rhs.
    for col in 0..k {
        // Partial pivot.
        let mut piv = col;
        for r in (col + 1)..k {
            if g[r][col].abs() > g[piv][col].abs() {
                piv = r;
            }
        }
        g.swap(col, piv);
        rhs.swap(col, piv);
        let d = g[col][col];
        if d.abs() < 1e-12 {
            continue;
        }
        for r in 0..k {
            if r == col {
                continue;
            }
            let f = g[r][col] / d;
            if f != 0.0 {
                for c in col..k {
                    g[r][c] -= f * g[col][c];
                }
                rhs[r] -= f * rhs[col];
            }
        }
    }
    (0..k)
        .map(|i| if g[i][i].abs() < 1e-12 { 0.0 } else { rhs[i] / g[i][i] })
        .collect()
}

/// Fits `components` jointly to the observed `(mz, intensity)` window by NNLS over the **union grid** of
/// all components' predicted teeth ∪ observed peaks in the span, and returns the amplitudes and the
/// cosine of the summed model to the observed. Peaks below `noise_floor` are ignored.
pub fn joint_envelope_fit(
    mz: &[f64],
    intensity: &[f64],
    components: &[Component],
    tol_ppm: f64,
    min_rel: f64,
    noise_floor: f64,
) -> JointFitResult {
    let n = components.len();
    if n == 0 || mz.is_empty() {
        return JointFitResult { amplitudes: vec![0.0; n], fit: 0.0 };
    }
    let teeth: Vec<Vec<(f64, f64)>> = components.iter().map(|c| component_teeth(c, min_rel)).collect();
    // Span of the model.
    let (mut lo_mz, mut hi_mz) = (f64::INFINITY, f64::NEG_INFINITY);
    for t in &teeth {
        for &(p, _) in t {
            lo_mz = lo_mz.min(p);
            hi_mz = hi_mz.max(p);
        }
    }
    if !lo_mz.is_finite() {
        return JointFitResult { amplitudes: vec![0.0; n], fit: 0.0 };
    }
    // Row positions: every predicted tooth, plus observed peaks in [lo,hi] not already near a tooth.
    let mut rows: Vec<f64> = teeth.iter().flat_map(|t| t.iter().map(|&(p, _)| p)).collect();
    let lo = mz.partition_point(|&m| m < lo_mz - 0.02);
    let hi = mz.partition_point(|&m| m <= hi_mz + 0.02);
    for j in lo..hi {
        if intensity[j] > noise_floor {
            rows.push(mz[j]);
        }
    }
    rows.sort_by(f64::total_cmp);
    rows.dedup_by(|a, b| (*a - *b).abs() / *b * 1e6 <= tol_ppm);

    let r = rows.len();
    let b: Vec<f64> = rows
        .iter()
        .map(|&p| nearest_within_ppm(mz, intensity, p, tol_ppm).unwrap_or(0.0))
        .collect();
    // Design matrix columns.
    let mut a_cols: Vec<Vec<f64>> = vec![vec![0.0; r]; n];
    for (jc, t) in teeth.iter().enumerate() {
        for &(pos, w) in t {
            // Row index of this tooth (nearest row position).
            if let Some(ri) = nearest_row(&rows, pos, tol_ppm) {
                a_cols[jc][ri] = w;
            }
        }
    }
    let amps = nnls(&a_cols, &b, 3 * n + 10);
    // Model and cosine.
    let mut model = vec![0.0f64; r];
    for jc in 0..n {
        if amps[jc] != 0.0 {
            for i in 0..r {
                model[i] += a_cols[jc][i] * amps[jc];
            }
        }
    }
    let dot: f64 = (0..r).map(|i| model[i] * b[i]).sum();
    let nm = model.iter().map(|x| x * x).sum::<f64>().sqrt();
    let nb = b.iter().map(|x| x * x).sum::<f64>().sqrt();
    let fit = if nm > 0.0 && nb > 0.0 { dot / (nm * nb) } else { 0.0 };
    JointFitResult { amplitudes: amps, fit }
}

/// Joint fit with a small search over the **target's** (component 0) monoisotope position: tries the
/// target mono at each `shift ∈ shifts` ¹³C units (neighbours fixed) and returns the shift whose joint
/// model fits best, as `(best_shift, JointFitResult)`.
pub fn joint_fit_target_shift(
    mz: &[f64],
    intensity: &[f64],
    components: &[Component],
    shifts: &[i32],
    tol_ppm: f64,
    min_rel: f64,
    noise_floor: f64,
) -> (i32, JointFitResult) {
    let mut best: Option<(i32, JointFitResult)> = None;
    let target = components[0];
    let spacing = C13_MINUS_C12 / target.charge.max(1).abs() as f64;
    for &s in shifts {
        let mut comps = components.to_vec();
        comps[0] = Component { mono_mz: target.mono_mz + s as f64 * spacing, charge: target.charge };
        let r = joint_envelope_fit(mz, intensity, &comps, tol_ppm, min_rel, noise_floor);
        if best.as_ref().map_or(true, |(_, br)| r.fit > br.fit) {
            best = Some((s, r));
        }
    }
    best.unwrap_or((0, JointFitResult { amplitudes: vec![0.0; components.len()], fit: 0.0 }))
}

fn nearest_within_ppm(mz: &[f64], intensity: &[f64], target: f64, tol_ppm: f64) -> Option<f64> {
    let ip = mz.partition_point(|&m| m < target);
    let mut best: Option<(f64, f64)> = None;
    for c in [ip.checked_sub(1), Some(ip)].into_iter().flatten() {
        if let Some(&m) = mz.get(c) {
            let d = (m - target).abs();
            if best.map_or(true, |(bd, _)| d < bd) {
                best = Some((d, intensity[c]));
            }
        }
    }
    match best {
        Some((d, v)) if d / target * 1e6 <= tol_ppm => Some(v),
        _ => None,
    }
}

fn nearest_row(rows: &[f64], target: f64, tol_ppm: f64) -> Option<usize> {
    let ip = rows.partition_point(|&m| m < target);
    let mut best: Option<(f64, usize)> = None;
    for c in [ip.checked_sub(1), Some(ip)].into_iter().flatten() {
        if let Some(&m) = rows.get(c) {
            let d = (m - target).abs();
            if best.map_or(true, |(bd, _)| d < bd) {
                best = Some((d, c));
            }
        }
    }
    match best {
        Some((d, i)) if d / target * 1e6 <= tol_ppm => Some(i),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::isotopic_envelope::mass_to_mz_f64;

    #[test]
    fn nnls_recovers_simple_amplitudes() {
        // Two orthogonal columns; NNLS should return the exact positive coefficients.
        let a = vec![vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0]];
        let b = vec![3.0, 5.0, 0.0];
        let x = nnls(&a, &b, 20);
        assert!((x[0] - 3.0).abs() < 1e-6 && (x[1] - 5.0).abs() < 1e-6, "got {x:?}");
    }

    #[test]
    fn nnls_clamps_negative_to_zero() {
        // b demands a negative coefficient on col 1; NNLS must clamp it to 0.
        let a = vec![vec![1.0, 0.0], vec![1.0, 0.0]];
        let b = vec![-2.0, -2.0];
        let x = nnls(&a, &b, 20);
        assert!(x[0] <= 1e-9 && x[1] <= 1e-9, "negative demand must clamp to 0, got {x:?}");
    }

    /// Builds a synthetic peak list = sum of two averagine envelopes at given (mono_mass, charge, amp).
    fn two_envelope_spectrum(
        a: (f64, i32, f64),
        b: (f64, i32, f64),
    ) -> (Vec<f64>, Vec<f64>) {
        let mut peaks: Vec<(f64, f64)> = Vec::new();
        for (mass, z, amp) in [a, b] {
            let mono_mz = mass_to_mz_f64(mass, z);
            let spacing = C13_MINUS_C12 / z as f64;
            let t = averagine_intensities_from_mono(mass, 1e-3, 20);
            for (k, &w) in t.iter().enumerate() {
                peaks.push((mono_mz + k as f64 * spacing, w * amp));
            }
        }
        peaks.sort_by(|p, q| p.0.total_cmp(&q.0));
        (peaks.iter().map(|p| p.0).collect(), peaks.iter().map(|p| p.1).collect())
    }

    #[test]
    fn joint_fit_separates_two_overlapping_envelopes() {
        // Two co-eluting z=2 species 3 Da apart — their envelope tails overlap (a chimeric window) but
        // the monos are distinct, so the joint fit is well-posed. It should explain the whole window and
        // recover the amplitude ratio.
        let (mz, inten) = two_envelope_spectrum((1500.0, 2, 1.0e7), (1503.0, 2, 6.0e6));
        let comps = [
            Component { mono_mz: mass_to_mz_f64(1500.0, 2), charge: 2 },
            Component { mono_mz: mass_to_mz_f64(1503.0, 2), charge: 2 },
        ];
        let r = joint_envelope_fit(&mz, &inten, &comps, 20.0, 0.05, 0.0);
        assert!(r.fit > 0.98, "joint model should explain the chimeric window, fit={}", r.fit);
        let ratio = r.amplitudes[0] / r.amplitudes[1];
        assert!((ratio - 10.0 / 6.0).abs() < 0.3, "amplitude ratio off: {:?}", r.amplitudes);
    }

    #[test]
    fn joint_shift_search_picks_correct_target_mono() {
        // Target truly at 1500.0 with a co-eluting neighbour 3 Da up. Seed the target one 13C too high;
        // the shift search should walk it back (shift -1) to the true mono, where the joint model — with
        // the neighbour explaining its own tail — fits best.
        let (mz, inten) = two_envelope_spectrum((1500.0, 2, 1.0e7), (1503.0, 2, 6.0e6));
        let comps = [
            Component { mono_mz: mass_to_mz_f64(1500.0 + C13_MINUS_C12, 2), charge: 2 },
            Component { mono_mz: mass_to_mz_f64(1503.0, 2), charge: 2 },
        ];
        let (shift, r) = joint_fit_target_shift(&mz, &inten, &comps, &[-2, -1, 0, 1], 20.0, 0.05, 0.0);
        assert_eq!(shift, -1, "should walk the +1-seeded target back to the true 1500.0 (got {shift})");
        assert!(r.fit > 0.98, "corrected joint fit should be high, {}", r.fit);
    }
}
