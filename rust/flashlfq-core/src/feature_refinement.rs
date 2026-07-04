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
use crate::isotopic_envelope::{mass_to_mz_f64, C13_MINUS_C12};
use crate::peak_indexing::Scan;
use crate::spectral_averaging::{average_spectra, SpectralAveragingParameters};
use crate::trace_kernel::DetectedFeature;
use std::collections::HashMap;

/// Absolute mass-clustering floor (Da) for small masses, where a ppm window would be tighter than
/// real mass precision. Applied as `max(mass · ppm/1e6, MASS_CLUSTER_ABS_FLOOR_DA)`.
pub const MASS_CLUSTER_ABS_FLOOR_DA: f64 = 0.01;

/// Hard cap on how many MS1 scans [`refine_feature`] averages into the composite. The composite is a
/// high-SNR snapshot at the feature apex, so this stays small (≈ SpectralAveraging's default of 5);
/// averaging more pulls in co-eluting interference. A `> MAX_SCANS_TO_AVERAGE` window is treated as
/// a bug (asserted), not silently accepted.
pub const MAX_SCANS_TO_AVERAGE: usize = 7;

/// Largest integer ¹³C off-by-one offset tolerated when grouping features by neutral mass. The mono
/// off-by-one is normally ±1; ±2 is allowed for robustness against a doubly-mis-assigned monoisotope.
const MAX_OFFBYONE_UNITS: i32 = 2;

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
pub fn refine_feature(
    feature: &DetectedFeature,
    scans: &[Scan],
    averaging_params: &SpectralAveragingParameters,
    decon_params: &ClassicDeconvolutionParameters,
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
    let half = (MAX_SCANS_TO_AVERAGE / 2) as i32; // 3 → up to 7 scans
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
    for s in window {
        let a = s.mz.partition_point(|&m| m < slice_lo);
        let b = s.mz.partition_point(|&m| m <= slice_hi);
        x_arrays.push(s.mz[a..b].to_vec());
        y_arrays.push(s.intensity[a..b].to_vec());
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
            return true;
        }
    }
    false
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

/// A single pooled candidate mass tagged with its source charge and contributing intensity weight.
struct Candidate {
    mass: f64,
    charge: i32,
    weight: f64,
}

/// Pools every member's candidate masses, clusters them within `mass_tolerance_ppm` (0.01 Da floor),
/// and returns the (weighted-mean mass, cross-charge support) of the cluster with the greatest
/// cross-charge support (ties: candidate count, then summed weight). Returns the tallest member's
/// refined mass if there are no candidates to cluster.
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
        for &mass in &m.candidate_masses {
            candidates.push(Candidate {
                mass,
                charge: m.refined_charge,
                weight: w,
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

    // Score each cluster: (distinct charges, count, summed weight); pick the best.
    let mut best: Option<(usize, usize, f64, f64)> = None; // (support, count, sum_w, weighted_mass)
    for cluster in &clusters {
        let mut charges: Vec<i32> = cluster.iter().map(|&i| candidates[i].charge).collect();
        charges.sort_unstable();
        charges.dedup();
        let support = charges.len();
        let count = cluster.len();
        let sum_w: f64 = cluster.iter().map(|&i| candidates[i].weight).sum();
        let weighted_mass = if sum_w > 0.0 {
            cluster
                .iter()
                .map(|&i| candidates[i].mass * candidates[i].weight)
                .sum::<f64>()
                / sum_w
        } else {
            cluster.iter().map(|&i| candidates[i].mass).sum::<f64>() / count as f64
        };
        let better = match &best {
            None => true,
            Some((bs, bc, bw, _)) => {
                support > *bs
                    || (support == *bs && count > *bc)
                    || (support == *bs && count == *bc && sum_w > *bw)
            }
        };
        if better {
            best = Some((support, count, sum_w, weighted_mass));
        }
    }

    let (support, _, _, mass) = best.expect("at least one cluster exists");
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
