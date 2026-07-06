//! Spectral averaging — faithful port of mzLib's `SpectralAveraging` project, **default-config
//! subset**.
//!
//! Purpose in the untargeted feature-detection pipeline: once the trace kernel has located a
//! candidate feature and derived an RT window, we average that window's MS1 scans into a single
//! composite spectrum (SNR ~√N) before the final [`crate::deconvolution`] pass. The **m/z binning**
//! is the load-bearing step — it is what gives the composite its mass accuracy — so it is ported
//! verbatim (see `agent_info/Feature-Detection-Design.md`, "Window + averaging").
//!
//! ## Scope (default-config subset)
//! mzLib's full project supports several outlier-rejection algorithms, several weighting schemes,
//! and file-level scan windowing. This port covers the configuration the pipeline uses faithfully
//! (m/z binning, even/TIC weighting, all normalization variants, and both `NoRejection` and
//! `SigmaClipping`), and stubs the remaining rejection/weighting algorithms until needed:
//!
//! | Knob | mzLib default | Ported here |
//! |------|---------------|-------------|
//! | `SpectralAveragingType` | `MzBinning` (only variant) | ✅ full |
//! | `NormalizationType` | `RelativeToTics` | ✅ all four variants (trivial) |
//! | `SpectraWeightingType` | `WeightEvenly` | ✅ `WeightEvenly` + `TicValue`; `MrsNoiseEstimation` panics |
//! | `OutlierRejectionType` | `NoRejection` | ✅ `NoRejection` + `SigmaClipping`; other clipping variants panic |
//! | `BinSize` | `0.01` | ✅ |
//!
//! The file-level windowing (`SpectraFileAveraging`, `AverageEverynScansWithOverlap`, scan overlap,
//! output type, thread count) is **not** ported: feature detection supplies its own scan window, so
//! the relevant surface is just `AverageSpectra(double[][] xArrays, double[][] yArrays, params)`.
//!
//! ## Parity
//! Control flow, summation order, the `floor((x - minX) / binSize)` bin index, and the "divide the
//! bin's summed intensity by the *spectrum* count (not the present-peak count)" behaviour are
//! preserved so a future C# golden matches at the standard tolerance (counts exact, floats
//! rel-1e-6). The one deliberate departure is performance-motivated and applies **only to the
//! `NoRejection` path**: mzLib pads every bin with a zero-intensity peak for each spectrum that did
//! not contribute a real peak, then averages over the padded set. When nothing is rejected that
//! padding is algebraically inert — a zero peak adds nothing to the weighted numerator and sits at
//! the running m/z mean — so [`average_bin`] elides it and folds its only real effect (dividing by
//! the full spectrum count, via the summed weight) into the averaging denominator directly. Under
//! outlier rejection the padding is *not* inert (absent-spectrum zeros participate in the clip
//! statistics and shift the surviving-weight denominator), so [`average_bin_rejected`] materializes
//! it explicitly. Results are identical up to floating-point summation order, well within the
//! rel-1e-6 tolerance. Where mzLib mutates the caller's `yArrays`
//! in place during normalization, this port normalizes an internal clone instead — the arithmetic
//! is identical; only the (undesirable) side effect on the caller is dropped.

// ---------------------------------------------------------------------------
// Configuration enums (mirroring SpectralAveraging/DataStructures/Enums)
// ---------------------------------------------------------------------------

/// `SpectralAveraging.OutlierRejectionType`. [`OutlierRejectionType::NoRejection`] and
/// [`OutlierRejectionType::SigmaClipping`] are ported; the remaining clipping variants are carried
/// for config/parity completeness but panic if dispatched (see module scope).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutlierRejectionType {
    NoRejection,
    MinMaxClipping,
    PercentileClipping,
    SigmaClipping,
    WinsorizedSigmaClipping,
    AveragedSigmaClipping,
    BelowThresholdRejection,
}

/// `SpectralAveraging.SpectraWeightingType`. [`SpectraWeightingType::WeightEvenly`] and
/// [`SpectraWeightingType::TicValue`] are ported; `MrsNoiseEstimation` panics (it needs the
/// MRS noise estimator + biweight midvariance, out of the default-config subset).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpectraWeightingType {
    WeightEvenly,
    MrsNoiseEstimation,
    TicValue,
}

/// `SpectralAveraging.NormalizationType`. All four variants are ported (they are trivial).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalizationType {
    NoNormalization,
    RelativeToTics,
    AbsoluteToTic,
    RelativeIntensity,
}

/// `SpectralAveraging.SpectralAveragingType`. mzLib defines only one variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpectralAveragingType {
    MzBinning,
}

// ---------------------------------------------------------------------------
// Parameters (mirroring SpectralAveragingParameters, averaging-only subset)
// ---------------------------------------------------------------------------

/// The subset of `SpectralAveragingParameters` that governs [`average_spectra`]. The file-level
/// windowing fields (`SpectraFileAveragingType`, `NumberOfScansToAverage`, `ScanOverlap`,
/// `OutputType`, `MaxThreadsToUsePerFile`) are intentionally omitted — feature detection supplies
/// its own scan window. `MinSigmaValue`/`MaxSigmaValue`/`Percentile` are carried so the struct can
/// carry the sigma-clipping bounds (and the not-yet-ported `Percentile`) without changing shape.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralAveragingParameters {
    pub outlier_rejection_type: OutlierRejectionType,
    pub spectral_weighting_type: SpectraWeightingType,
    pub spectral_averaging_type: SpectralAveragingType,
    pub normalization_type: NormalizationType,
    pub bin_size: f64,
    pub percentile: f64,
    pub min_sigma_value: f64,
    pub max_sigma_value: f64,
}

impl Default for SpectralAveragingParameters {
    /// FlashLFQ's configured averaging standard. This mirrors `SpectralAveragingParameters
    /// .SetDefaultValues()` for weighting (`WeightEvenly`), averaging (`MzBinning`), normalization
    /// (`RelativeToTics`) and bin size (`0.01`), but **deliberately deviates on outlier rejection**:
    /// mzLib's `SetDefaultValues()` uses `NoRejection`, whereas the pipeline runs
    /// [`OutlierRejectionType::SigmaClipping`] with min/max σ `0.5`/`3.0`. The asymmetric bounds
    /// clip the low tail of each bin (dropouts and absent-spectrum zeros) aggressively while leaving
    /// real high-intensity signal essentially untouched, which sharpens the averaged composite.
    fn default() -> Self {
        SpectralAveragingParameters {
            outlier_rejection_type: OutlierRejectionType::SigmaClipping,
            spectral_weighting_type: SpectraWeightingType::WeightEvenly,
            spectral_averaging_type: SpectralAveragingType::MzBinning,
            normalization_type: NormalizationType::RelativeToTics,
            bin_size: 0.01,
            percentile: 0.1,
            min_sigma_value: 0.5,
            max_sigma_value: 3.0,
        }
    }
}

// ---------------------------------------------------------------------------
// BinnedPeak (mirroring DataStructures/BinnedPeak.cs)
// ---------------------------------------------------------------------------

/// One peak assigned to an m/z bin, tagged with the spectrum it came from. Mirrors the internal
/// `BinnedPeak` record struct. Zero-intensity instances are synthesized to pad bins that a given
/// spectrum did not contribute a real peak to (see [`get_bins`]).
#[derive(Debug, Clone, Copy)]
struct BinnedPeak {
    mz: f64,
    intensity: f64,
    spectra_id: usize,
}

// ---------------------------------------------------------------------------
// Public entry point (mirroring SpectraAveraging.AverageSpectra)
// ---------------------------------------------------------------------------

/// Averages a group of spectra into one composite `(mz, intensity)` spectrum.
///
/// Faithful port of `SpectraAveraging.AverageSpectra` → `MzBinning`. `x_arrays[i]` / `y_arrays[i]`
/// are the m/z and intensity arrays of the *i*-th spectrum (each internally ascending in m/z, as
/// mzML centroided scans are). Returns `(mz, intensity)` for the averaged spectrum, ascending in
/// m/z, with zero-intensity bins dropped.
///
/// Panics if `x_arrays`/`y_arrays` differ in outer length, if any paired inner arrays differ in
/// length, or if the input is empty — matching the fact that the C# would throw / divide-by-zero
/// on those degenerate inputs rather than return a meaningful spectrum.
pub fn average_spectra(
    x_arrays: &[Vec<f64>],
    y_arrays: &[Vec<f64>],
    parameters: &SpectralAveragingParameters,
) -> (Vec<f64>, Vec<f64>) {
    match parameters.spectral_averaging_type {
        SpectralAveragingType::MzBinning => mz_binning(x_arrays, y_arrays, parameters),
    }
}

/// Faithful port of `SpectraAveraging.MzBinning`.
fn mz_binning(
    x_arrays: &[Vec<f64>],
    y_arrays: &[Vec<f64>],
    parameters: &SpectralAveragingParameters,
) -> (Vec<f64>, Vec<f64>) {
    assert_eq!(
        x_arrays.len(),
        y_arrays.len(),
        "x and y arrays must have the same number of spectra"
    );
    assert!(!x_arrays.is_empty(), "cannot average zero spectra");
    for (x, y) in x_arrays.iter().zip(y_arrays.iter()) {
        assert_eq!(x.len(), y.len(), "each spectrum's x and y arrays must match in length");
    }

    // normalize spectra — mzLib mutates the caller's arrays here; we normalize an internal clone
    // (and skip even that allocation when there is nothing to normalize).
    let mut owned;
    let y_norm: &[Vec<f64>] = match parameters.normalization_type {
        NormalizationType::NoNormalization => y_arrays,
        nt => {
            owned = y_arrays.to_vec();
            normalize_spectra(&mut owned, nt);
            &owned
        }
    };

    // get bins (real peaks only; the zero padding is re-derived per bin where it matters — folded
    // into average_bin for NoRejection, materialized in average_bin_rejected otherwise)
    let bins = get_bins(x_arrays, y_norm, parameters.bin_size);

    // get weights. For NoRejection the averaging denominator is the summed weight over *all* spectra
    // (the padding's only effect), hoisted out of the per-bin loop since it is bin-independent; the
    // rejection path instead sums the weights of the surviving peaks per bin.
    let weights = calculate_spectra_weights(x_arrays, y_norm, parameters.spectral_weighting_type);
    let total_weight = sum(&weights);

    // reject outliers and average bins. NoRejection keeps the padding-free fast path (average_bin);
    // every other config must materialize the zero padding and reject over it (average_bin_rejected),
    // because the rejection statistics and the surviving-weight denominator both depend on the
    // padded set — the algebraic shortcut only holds when nothing is rejected.
    let num_spectra = x_arrays.len();
    let mut averaged_peaks: Vec<(f64, f64)> = Vec::with_capacity(bins.len());
    for peaks_from_bin in &bins {
        let averaged = match parameters.outlier_rejection_type {
            OutlierRejectionType::NoRejection => average_bin(peaks_from_bin, &weights, total_weight),
            _ => average_bin_rejected(peaks_from_bin, &weights, num_spectra, parameters),
        };
        averaged_peaks.push(averaged);
    }

    // return averaged: drop zero-intensity bins, order by m/z; AbsoluteToTic re-scales by averageTic
    let mut ordered: Vec<(f64, f64)> =
        averaged_peaks.into_iter().filter(|p| p.1 != 0.0).collect();
    ordered.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let mzs: Vec<f64> = ordered.iter().map(|p| p.0).collect();
    let intensities: Vec<f64> = match parameters.normalization_type {
        // Only this branch needs the average TIC, so only this branch pays for computing it.
        NormalizationType::AbsoluteToTic => {
            let average_tic =
                y_arrays.iter().map(|y| sum(y)).sum::<f64>() / y_arrays.len() as f64;
            ordered.iter().map(|p| p.1 * average_tic).collect()
        }
        _ => ordered.iter().map(|p| p.1).collect(),
    };
    (mzs, intensities)
}

// ---------------------------------------------------------------------------
// Binning (mirroring SpectraAveraging.GetBins / AverageBin)
// ---------------------------------------------------------------------------

/// Port of `SpectraAveraging.GetBins`, minus the zero-padding step. Sorts every peak of every
/// spectrum into an m/z bin (`floor((mz - minX) / binSize)`) and returns the per-bin groups of
/// **real** peaks. mzLib additionally pads each bin with a zero-intensity peak for every absent
/// spectrum; that padding is algebraically inert (see the module "Parity" note) and is folded into
/// [`average_bin`] instead of being materialized here.
///
/// Returned as a `Vec` of bins rather than a keyed map: the bin *index* never affects the output
/// (the final spectrum is re-sorted by m/z), so only the per-bin peak groups matter. Bins are
/// produced in ascending bin-index order.
fn get_bins(x_arrays: &[Vec<f64>], y_arrays: &[Vec<f64>], bin_size: f64) -> Vec<Vec<BinnedPeak>> {
    let num_spectra = x_arrays.len();
    let min_x_value = x_arrays
        .iter()
        .flat_map(|x| x.iter().copied())
        .fold(f64::INFINITY, f64::min);

    // Sort all peaks into bins keyed by bin index. BTreeMap gives deterministic ascending-index
    // iteration; the C# uses a Dictionary whose order is irrelevant post-sort.
    use std::collections::BTreeMap;
    let mut bins: BTreeMap<i64, Vec<BinnedPeak>> = BTreeMap::new();
    for i in 0..num_spectra {
        for j in 0..x_arrays[i].len() {
            let mz = x_arrays[i][j];
            let bin_index = ((mz - min_x_value) / bin_size).floor() as i64;
            bins.entry(bin_index).or_default().push(BinnedPeak {
                mz,
                intensity: y_arrays[i][j],
                spectra_id: i,
            });
        }
    }

    bins.into_values().collect()
}

/// Port of `SpectraAveraging.AverageBin` with the zero-intensity padding folded in analytically.
///
/// mzLib computes the weighted mean over the *padded* peak set (real peaks plus one zero-intensity
/// peak per absent spectrum): numerator `Σ intensity·weight`, denominator `Σ weight`, and the plain
/// arithmetic mean of every peak's m/z. Because a padded peak has zero intensity and sits exactly at
/// the bin's running m/z mean, it changes neither the numerator nor the mean — its sole effect is to
/// add its spectrum's weight to the denominator. Summed over the whole bin that denominator is just
/// the weight over *all* spectra (`total_weight`, passed in), so this operates on the real peaks
/// alone: numerator over them, denominator `total_weight`, m/z = their arithmetic mean.
fn average_bin(real_peaks: &[BinnedPeak], weights: &[f64], total_weight: f64) -> (f64, f64) {
    let mut numerator = 0.0;
    let mut mz_sum = 0.0;
    for peak in real_peaks {
        numerator += peak.intensity * weights[peak.spectra_id];
        mz_sum += peak.mz;
    }
    let mz = mz_sum / real_peaks.len() as f64;
    let intensity = numerator / total_weight;
    (mz, intensity)
}

/// `SpectraAveraging.AverageBin` composed with `OutlierRejection.RejectOutliers`, for every config
/// other than `NoRejection`. Unlike [`average_bin`], the zero-intensity padding cannot be elided:
/// mzLib rejects outliers over the *padded* peak set (so an absent spectrum's zero participates in
/// the median/σ statistics and can itself be clipped), and averages the survivors with *their*
/// summed weight as the denominator, not the weight over all spectra. So we materialize the padding
/// (mzLib does this in `GetBins`), reject, then weight-average whatever remains.
///
/// Returns `(0.0, 0.0)` when every peak is rejected; the caller drops zero-intensity bins, matching
/// mzLib's `if (!peaksFromBin.Any()) continue;`.
fn average_bin_rejected(
    real_peaks: &[BinnedPeak],
    weights: &[f64],
    num_spectra: usize,
    parameters: &SpectralAveragingParameters,
) -> (f64, f64) {
    // Materialize mzLib's zero padding: each spectrum with no real peak in this bin contributes a
    // zero-intensity peak at the bin's mean real m/z (mzLib pads at the running m/z mean, which is
    // that same value since each padded peak sits exactly on it).
    let mz_mean = real_peaks.iter().map(|p| p.mz).sum::<f64>() / real_peaks.len() as f64;
    let mut padded: Vec<BinnedPeak> = real_peaks.to_vec();
    let mut present = vec![false; num_spectra];
    for peak in real_peaks {
        present[peak.spectra_id] = true;
    }
    for (id, &seen) in present.iter().enumerate() {
        if !seen {
            padded.push(BinnedPeak {
                mz: mz_mean,
                intensity: 0.0,
                spectra_id: id,
            });
        }
    }

    let survivors = reject_outliers(padded, parameters);
    if survivors.is_empty() {
        return (0.0, 0.0);
    }

    let mut numerator = 0.0;
    let mut denominator = 0.0;
    let mut mz_sum = 0.0;
    for peak in &survivors {
        numerator += peak.intensity * weights[peak.spectra_id];
        denominator += weights[peak.spectra_id];
        mz_sum += peak.mz;
    }
    (mz_sum / survivors.len() as f64, numerator / denominator)
}

/// Port of `OutlierRejection.RejectOutliers(List<BinnedPeak>, …)`: run the configured rejection over
/// the peaks' intensities, then keep the peaks whose intensity survived. Membership is by intensity
/// value, mirroring the C# overload; this is exact because sigma clipping never splits equal values
/// (identical intensities share the same clip predicate, so they are kept or dropped together).
fn reject_outliers(
    peaks: Vec<BinnedPeak>,
    parameters: &SpectralAveragingParameters,
) -> Vec<BinnedPeak> {
    let survivors = match parameters.outlier_rejection_type {
        OutlierRejectionType::NoRejection => return peaks,
        OutlierRejectionType::SigmaClipping => sigma_clipping(
            peaks.iter().map(|p| p.intensity).collect(),
            parameters.min_sigma_value,
            parameters.max_sigma_value,
        ),
        other => unimplemented!(
            "outlier rejection {:?} is outside the ported subset; only NoRejection and \
             SigmaClipping are ported (see module scope)",
            other
        ),
    };
    peaks
        .into_iter()
        .filter(|p| survivors.contains(&p.intensity))
        .collect()
}

/// Port of `OutlierRejection.SigmaClipping`: iteratively drop values lying more than `s_min` σ below
/// or `s_max` σ above the current median, recomputing the median and (sample) σ each pass, until a
/// pass rejects nothing. The two bounds are asymmetric by design — feature detection clips the low
/// tail (dropouts, absent-spectrum zeros) aggressively while keeping the high tail (real signal).
fn sigma_clipping(mut values: Vec<f64>, s_min: f64, s_max: f64) -> Vec<f64> {
    loop {
        let med = median(&values);
        let std_dev = sample_standard_deviation(&values);
        let before = values.len();
        values.retain(|&v| !should_clip(v, med, std_dev, s_min, s_max));
        if values.len() == before {
            break;
        }
    }
    values
}

/// mzLib's `SigmaClipping` per-value predicate: reject `value` when it is more than `s_min` σ below
/// the median or more than `s_max` σ above it. When σ is zero or NaN (constant or singleton input)
/// the comparisons evaluate false, so nothing is rejected — matching C# `double` semantics and
/// letting [`sigma_clipping`] terminate.
fn should_clip(value: f64, median: f64, std_dev: f64, s_min: f64, s_max: f64) -> bool {
    (median - value) / std_dev > s_min || (value - median) / std_dev > s_max
}

/// Sample median (`MathNet.Numerics.Statistics.Median`): the central order statistic, averaging the
/// two central values for even length. Returns NaN for empty input.
fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let mut sorted: Vec<f64> = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    }
}

/// Sample standard deviation (`MathNet.Numerics.Statistics.StandardDeviation`, Bessel-corrected,
/// dividing by `n - 1`). Returns NaN for fewer than two values, matching MathNet.
fn sample_standard_deviation(values: &[f64]) -> f64 {
    let n = values.len();
    if n < 2 {
        return f64::NAN;
    }
    let mean = values.iter().sum::<f64>() / n as f64;
    let sum_sq: f64 = values.iter().map(|&v| (v - mean) * (v - mean)).sum();
    (sum_sq / (n as f64 - 1.0)).sqrt()
}

// ---------------------------------------------------------------------------
// Normalization (mirroring SpectraNormalization)
// ---------------------------------------------------------------------------

/// Faithful port of `SpectraNormalization.NormalizeSpectra`. Mutates `y_arrays` in place.
fn normalize_spectra(y_arrays: &mut [Vec<f64>], normalization_type: NormalizationType) {
    match normalization_type {
        NormalizationType::NoNormalization => {}
        NormalizationType::AbsoluteToTic => normalize_absolute_to_tic(y_arrays),
        NormalizationType::RelativeToTics => normalize_relative_to_tics(y_arrays),
        NormalizationType::RelativeIntensity => to_relative_intensity(y_arrays),
    }
}

/// Divide each y by its own TIC (sum-to-one per spectrum). `NormalizeAbsoluteToTic`.
fn normalize_absolute_to_tic(y_arrays: &mut [Vec<f64>]) {
    for y in y_arrays.iter_mut() {
        let mut total_ion_current = sum(y);
        if total_ion_current == 0.0 {
            total_ion_current = 1.0;
        }
        for v in y.iter_mut() {
            *v /= total_ion_current;
        }
    }
}

/// Divide each y by its own TIC then multiply by the average TIC. `NormalizeRelativeToTics`
/// — the default. Puts every spectrum on the same total-intensity scale while preserving the
/// overall magnitude.
fn normalize_relative_to_tics(y_arrays: &mut [Vec<f64>]) {
    let tics: Vec<f64> = y_arrays.iter().map(|y| sum(y)).collect();
    let average_tic = sum(&tics) / tics.len() as f64;
    for i in 0..y_arrays.len() {
        let tic = if tics[i] == 0.0 { 1.0 } else { tics[i] };
        for v in y_arrays[i].iter_mut() {
            *v = *v / tic * average_tic;
        }
    }
}

/// Divide each y by the spectrum's maximum intensity. `ToRelativeIntensity`.
fn to_relative_intensity(y_arrays: &mut [Vec<f64>]) {
    for y in y_arrays.iter_mut() {
        let max_value = y.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        for v in y.iter_mut() {
            *v /= max_value;
        }
    }
}

// ---------------------------------------------------------------------------
// Weighting (mirroring SpectralWeighting)
// ---------------------------------------------------------------------------

/// Faithful port of `SpectralWeighting.CalculateSpectraWeights`. Returns weights indexed by
/// spectrum id (mzLib returns a `Dictionary<int,double>`; the keys are a dense `0..count`, so a
/// `Vec` indexed by id is equivalent).
fn calculate_spectra_weights(
    x_arrays: &[Vec<f64>],
    y_arrays: &[Vec<f64>],
    spectra_weighting_type: SpectraWeightingType,
) -> Vec<f64> {
    match spectra_weighting_type {
        SpectraWeightingType::WeightEvenly => vec![1.0; x_arrays.len()],
        SpectraWeightingType::TicValue => weight_by_tic_value(y_arrays),
        SpectraWeightingType::MrsNoiseEstimation => unimplemented!(
            "MrsNoiseEstimation weighting is outside the default-config subset (needs the MRS \
             noise estimator + biweight midvariance); port it when a config requires it"
        ),
    }
}

/// Weight each spectrum by `tic_i / max_tic`. `WeightByTicValue`.
fn weight_by_tic_value(y_arrays: &[Vec<f64>]) -> Vec<f64> {
    let tics: Vec<f64> = y_arrays.iter().map(|y| sum(y)).collect();
    let max_tic = tics.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    tics.iter().map(|&t| t / max_tic).collect()
}

// ---------------------------------------------------------------------------
// Outlier rejection: only the default `NoRejection` config is ported, and it is a no-op folded
// directly into `mz_binning` (all bins pass through). The six clipping variants are dispatched up
// front in `mz_binning` and panic — they operate on the zero-padded bin representation this fast
// path elides, so they belong with the not-yet-ported alternate-config work (see module scope).
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Small numeric helpers (kept explicit for summation-order parity)
// ---------------------------------------------------------------------------

/// Left-to-right sum, matching `IEnumerable<double>.Sum()` accumulation order.
#[inline]
fn sum(values: &[f64]) -> f64 {
    values.iter().copied().sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64) {
        assert!((a - b).abs() <= 1e-9 * b.abs().max(1.0), "expected {b}, got {a}");
    }

    /// Two spectra, one shared bin (0.01 wide), default config. Hand-computed:
    /// tics 10/20 → averageTic 15 → both normalize to 15 → bin mean intensity 15,
    /// bin m/z = mean(100.00, 100.005) = 100.0025.
    #[test]
    fn two_spectra_single_shared_bin() {
        let params = SpectralAveragingParameters {
            outlier_rejection_type: OutlierRejectionType::NoRejection,
            ..SpectralAveragingParameters::default()
        };
        let x = vec![vec![100.000], vec![100.005]];
        let y = vec![vec![10.0], vec![20.0]];
        let (mz, inten) = average_spectra(&x, &y, &params);
        assert_eq!(mz.len(), 1);
        approx(mz[0], 100.0025);
        approx(inten[0], 15.0);
    }

    /// Padding path: spectrum B lacks the 200 m/z bin, so it is padded with a zero peak, and the
    /// bin intensity divides by the spectrum count (2), not the present-peak count (1).
    /// A: x=[100,200] y=[10,30] tic 40; B: x=[100] y=[20] tic 20; averageTic 30.
    /// Normalized A=[7.5,22.5], B=[30]. Bin100: (7.5+30)/2=18.75. Bin200: (22.5+0)/2=11.25.
    #[test]
    fn padding_divides_by_spectrum_count() {
        let params = SpectralAveragingParameters {
            outlier_rejection_type: OutlierRejectionType::NoRejection,
            ..SpectralAveragingParameters::default()
        };
        let x = vec![vec![100.000, 200.000], vec![100.000]];
        let y = vec![vec![10.0, 30.0], vec![20.0]];
        let (mz, inten) = average_spectra(&x, &y, &params);
        assert_eq!(mz.len(), 2);
        approx(mz[0], 100.0);
        approx(inten[0], 18.75);
        approx(mz[1], 200.0);
        approx(inten[1], 11.25);
    }

    /// With `NoNormalization` + `WeightEvenly`, a single spectrum passes through as-is (bins of one
    /// real peak each, no padding), only re-sorted by m/z with zero-intensity peaks dropped.
    #[test]
    fn no_normalization_single_spectrum_passthrough() {
        let params = SpectralAveragingParameters {
            normalization_type: NormalizationType::NoNormalization,
            outlier_rejection_type: OutlierRejectionType::NoRejection,
            ..SpectralAveragingParameters::default()
        };
        let x = vec![vec![300.0, 100.0, 200.0]];
        let y = vec![vec![3.0, 1.0, 2.0]];
        let (mz, inten) = average_spectra(&x, &y, &params);
        assert_eq!(mz, vec![100.0, 200.0, 300.0]);
        assert_eq!(inten, vec![1.0, 2.0, 3.0]);
    }

    /// Zero-intensity bins are dropped from the output entirely.
    #[test]
    fn zero_intensity_bins_dropped() {
        let params = SpectralAveragingParameters {
            normalization_type: NormalizationType::NoNormalization,
            outlier_rejection_type: OutlierRejectionType::NoRejection,
            ..SpectralAveragingParameters::default()
        };
        let x = vec![vec![100.0, 200.0]];
        let y = vec![vec![0.0, 5.0]];
        let (mz, inten) = average_spectra(&x, &y, &params);
        assert_eq!(mz, vec![200.0]);
        assert_eq!(inten, vec![5.0]);
    }

    /// `TicValue` weighting: two spectra sharing a bin, weights tic_i / max_tic.
    /// A tic 10 → w 0.5; B tic 20 → w 1.0. NoNormalization keeps raw intensities.
    /// intensity = (10·0.5 + 20·1.0)/(0.5+1.0) = 25/1.5 = 16.666…
    #[test]
    fn tic_value_weighting() {
        let params = SpectralAveragingParameters {
            normalization_type: NormalizationType::NoNormalization,
            spectral_weighting_type: SpectraWeightingType::TicValue,
            outlier_rejection_type: OutlierRejectionType::NoRejection,
            ..SpectralAveragingParameters::default()
        };
        let x = vec![vec![100.000], vec![100.005]];
        let y = vec![vec![10.0], vec![20.0]];
        let (mz, inten) = average_spectra(&x, &y, &params);
        assert_eq!(mz.len(), 1);
        approx(inten[0], 25.0 / 1.5);
    }

    /// `sigma_clipping` iteratively drops a low outlier. [10,10,10,10,2] at min σ 0.5: the 2 lies
    /// ~2.24 σ below the median (10) so it is rejected; the remaining values are identical (σ = 0),
    /// so the next pass rejects nothing and the loop terminates.
    #[test]
    fn sigma_clipping_drops_low_outlier() {
        let survivors = sigma_clipping(vec![10.0, 10.0, 10.0, 10.0, 2.0], 0.5, 3.0);
        assert_eq!(survivors, vec![10.0, 10.0, 10.0, 10.0]);
    }

    /// The min/max bounds are asymmetric: with min σ 0.5 / max σ 3.0 on [5,10,10,10,15] (median 10,
    /// σ ≈ 3.54), the low value 5 is ~1.41 σ below the median and is rejected, while the equally
    /// distant high value 15 is only ~1.41 σ above (< 3) and is kept.
    #[test]
    fn sigma_clipping_bounds_are_asymmetric() {
        let survivors = sigma_clipping(vec![5.0, 10.0, 10.0, 10.0, 15.0], 0.5, 3.0);
        assert_eq!(survivors, vec![10.0, 10.0, 10.0, 15.0]);
    }

    /// Constant and singleton inputs have σ = 0 / NaN, so nothing is rejected and the loop halts.
    #[test]
    fn sigma_clipping_terminates_on_degenerate_input() {
        assert_eq!(sigma_clipping(vec![7.0, 7.0, 7.0], 0.5, 3.0), vec![7.0, 7.0, 7.0]);
        assert_eq!(sigma_clipping(vec![7.0], 0.5, 3.0), vec![7.0]);
    }

    /// End-to-end: five spectra sharing one bin with intensities [10,10,10,10,2]. Under the default
    /// sigma-clipping config the low 2 is clipped, so the composite intensity is the mean of the
    /// four survivors (10.0); with `NoRejection` it would be (40 + 2) / 5 = 8.4.
    #[test]
    fn sigma_clipping_end_to_end_clips_bin() {
        let x = vec![vec![100.0]; 5];
        let y = vec![vec![10.0], vec![10.0], vec![10.0], vec![10.0], vec![2.0]];

        let sigma = SpectralAveragingParameters {
            normalization_type: NormalizationType::NoNormalization,
            ..SpectralAveragingParameters::default()
        };
        let (mz, inten) = average_spectra(&x, &y, &sigma);
        assert_eq!(mz.len(), 1);
        approx(mz[0], 100.0);
        approx(inten[0], 10.0);

        let no_reject = SpectralAveragingParameters {
            normalization_type: NormalizationType::NoNormalization,
            outlier_rejection_type: OutlierRejectionType::NoRejection,
            ..SpectralAveragingParameters::default()
        };
        let (_, inten_nr) = average_spectra(&x, &y, &no_reject);
        approx(inten_nr[0], 8.4);
    }

    /// The zero-intensity padding for absent spectra participates in rejection. Four spectra have a
    /// peak at 100.0 (intensity 10); a fifth contributes only at 100.5. In the 100.0 bin the fifth
    /// spectrum's padded zero is ~2.24 σ below the median and is clipped, so the intensity is the
    /// mean of the four real peaks (10.0) rather than 40 / 5 = 8.0. In the 100.5 bin the lone real
    /// peak (10) sits among four zeros but only ~2.24 σ above the median (< 3 σ), so nothing is
    /// rejected and the intensity is 10 / 5 = 2.0.
    #[test]
    fn sigma_clipping_rejects_zero_padding() {
        let x = vec![
            vec![100.0],
            vec![100.0],
            vec![100.0],
            vec![100.0],
            vec![100.5],
        ];
        let y = vec![vec![10.0]; 5];
        let params = SpectralAveragingParameters {
            normalization_type: NormalizationType::NoNormalization,
            ..SpectralAveragingParameters::default()
        };
        let (mz, inten) = average_spectra(&x, &y, &params);
        assert_eq!(mz.len(), 2);
        approx(mz[0], 100.0);
        approx(inten[0], 10.0);
        approx(mz[1], 100.5);
        approx(inten[1], 2.0);
    }
}
