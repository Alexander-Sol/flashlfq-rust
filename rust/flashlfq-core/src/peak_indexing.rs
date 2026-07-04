//! Binned m/z peak index — port of mzLib `MassSpectrometry/PeakIndexing/IndexingEngine<T>`
//! specialised to `IndexedMassSpectralPeak` (the m/z peak flavour used by FlashLFQ's
//! `PeakIndexingEngine`).
//!
//! Ports the index build (`IndexPeaks`), the point query (`GetIndexedPeak`), and the
//! helpers it leans on (`GetBinsInRange`, `GetBestPeakFromBins`, `GetPeakFromBin`,
//! `BinarySearchForIndexedPeak`) plus the `BinsPerDalton` constant. XIC tracing
//! (`GetXic*`) is the next task (P1.11) and is **not** ported here.
//!
//! ## Fidelity notes
//! - `BinsPerDalton = 100` (the `PeakIndexingEngine` default; `IndexingEngine<T>` exposes it
//!   as a virtual `protected` member).
//! - The bin for a peak is `(int)Math.Round(mz * 100, 0)` — C# `Math.Round` is *banker's*
//!   (round-half-to-even), so this uses `f64::round_ties_even`.
//! - The jagged index is `Vec<Option<Vec<_>>>`: index `i` holds every peak whose rounded
//!   m/z bin is `i`, ordered by scan index ascending (peaks are appended in scan order).
//! - `IndexedMassSpectralPeak` stores `mz`/`intensity`/`retention_time` as `f32`, matching
//!   the C# class (the `(float)` casts in its constructor). The bin computation uses the
//!   original `f64` m/z, exactly as C# rounds `MassSpectrum.XArray[j]` (a `double`).

use std::collections::HashSet;
use std::path::PathBuf;

use mzdata::io::{DetailLevel, MZReader, ThermoRawReader};
use mzdata::prelude::*;
use mzdata::spectrum::Spectrum;

use crate::tolerance::PpmTolerance;

/// Number of m/z bins per dalton. The `PeakIndexingEngine` default in mzLib.
pub const BINS_PER_DALTON: f64 = 100.0;

/// A value identity for an indexed peak, used to mark a peak as already claimed by an XIC.
///
/// C# uses object reference equality in the `Dictionary<IIndexedPeak, ExtractedIonChromatogram>`
/// that `GetAllXics` threads through `GetXicByScanIndex`. Our peaks are `Copy` value types, so we
/// stand in a structural key: `(scan index, m/z bits, intensity bits)`. Two distinct peaks in the
/// same scan with identical m/z *and* intensity are functionally indistinguishable, so this keys
/// them the same — matching the only behaviour that matters for dedup.
pub type PeakKey = (i32, u32, u32);

/// Structural identity key for an indexed peak (see [`PeakKey`]).
#[inline]
fn peak_key(p: &IndexedMassSpectralPeak) -> PeakKey {
    (
        p.zero_based_scan_index,
        p.mz.to_bits(),
        p.intensity.to_bits(),
    )
}

/// A single indexed m/z peak. Port of mzLib `MassSpectrometry.IndexedMassSpectralPeak`.
///
/// `m()` returns the m/z (the C# `M => Mz` shortcut). All masses are `f32` to mirror the
/// C# class's `(float)` casts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IndexedMassSpectralPeak {
    /// The peak m/z (C# `Mz`, also surfaced as `M`).
    pub mz: f32,
    /// Peak intensity.
    pub intensity: f32,
    /// Zero-based index of the scan (within the indexed scan array) this peak came from.
    pub zero_based_scan_index: i32,
    /// Retention time of the scan, in minutes.
    pub retention_time: f32,
}

impl IndexedMassSpectralPeak {
    /// Builds a peak, applying the C# `(float)` narrowing casts on the mass fields.
    pub fn new(mz: f64, intensity: f64, zero_based_scan_index: i32, retention_time: f64) -> Self {
        IndexedMassSpectralPeak {
            mz: mz as f32,
            intensity: intensity as f32,
            zero_based_scan_index,
            retention_time: retention_time as f32,
        }
    }

    /// The peak's representative mass (m/z here). C# `M => Mz`.
    #[inline]
    pub fn m(&self) -> f32 {
        self.mz
    }

    /// The structural identity key used to mark this peak as claimed (see [`PeakKey`]).
    /// Exposed so callers outside this module (e.g. the trace-kernel detector) can maintain
    /// their own claim sets with the same identity semantics `get_all_xics` uses internally.
    #[inline]
    pub fn key(&self) -> PeakKey {
        peak_key(self)
    }
}

/// Per-scan metadata kept alongside the index. Port of mzLib `MassSpectrometry.ScanInfo`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScanInfo {
    /// One-based scan number as reported by the source file.
    pub one_based_scan_number: i32,
    /// Zero-based position of this scan within the indexed scan array.
    pub zero_based_scan_index: i32,
    /// Retention time, in minutes.
    pub retention_time: f64,
    /// MS level (1 for survey scans).
    pub msn_order: i32,
}

/// A raw scan to be indexed: parallel `mz`/`intensity` arrays plus scan metadata.
///
/// This is the input shape `index_peaks` consumes. It mirrors the fields of an mzLib
/// `MsDataScan` that the indexing engine reads (`MassSpectrum.XArray`/`YArray`,
/// `OneBasedScanNumber`, `RetentionTime`, `MsnOrder`).
#[derive(Debug, Clone)]
pub struct Scan {
    /// m/z values (assumed ascending, as in centroided mzML).
    pub mz: Vec<f64>,
    /// Intensity values, parallel to `mz`.
    pub intensity: Vec<f64>,
    /// One-based scan number from the source file.
    pub one_based_scan_number: i32,
    /// Retention time, in minutes.
    pub retention_time: f64,
    /// MS level.
    pub msn_order: i32,
}

/// The binned m/z index plus per-scan metadata.
///
/// Port of mzLib `PeakIndexingEngine` (the `IndexingEngine<IndexedMassSpectralPeak>`
/// specialisation), limited to indexing + point queries for P1.10.
#[derive(Debug, Clone)]
pub struct PeakIndexingEngine {
    /// Jagged index: `indexed_peaks[bin]` is the list of peaks whose rounded m/z bin is
    /// `bin`, ordered by scan index ascending. `None` for empty bins.
    indexed_peaks: Vec<Option<Vec<IndexedMassSpectralPeak>>>,
    /// Per-scan metadata, parallel to the indexed scan array (by zero-based scan index).
    scan_info: Vec<ScanInfo>,
}

impl PeakIndexingEngine {
    /// Builds the index from an array of scans (which should already be the MS1-only,
    /// scan-number-ordered array — `index_peaks` does not filter).
    ///
    /// Faithful port of `IndexingEngine<T>.IndexPeaks`: returns `None` (C# `false`/`null`)
    /// when there is nothing to index (empty, or no scan carries a peak).
    pub fn index_peaks(scans: &[Scan]) -> Option<PeakIndexingEngine> {
        if scans.is_empty() {
            return None;
        }

        // Largest m/z across all scans (C# `Max(p => p.MassSpectrum.LastX.Value)`), where
        // LastX is the last (max) m/z of a scan. Scans with no peaks are skipped, matching
        // the C# `Where(p => p.MassSpectrum.LastX != null)` filter.
        let max_last_x = scans
            .iter()
            .filter_map(|s| s.mz.last().copied())
            .fold(None::<f64>, |acc, x| Some(acc.map_or(x, |a| a.max(x))));

        let max_last_x = match max_last_x {
            Some(x) => x,
            None => return None, // no peaks anywhere
        };

        let len = (max_last_x * BINS_PER_DALTON).ceil() as usize + 1;
        let mut indexed_peaks: Vec<Option<Vec<IndexedMassSpectralPeak>>> = vec![None; len];
        let mut scan_info: Vec<ScanInfo> = Vec::with_capacity(scans.len());

        for (scan_index, scan) in scans.iter().enumerate() {
            scan_info.push(ScanInfo {
                one_based_scan_number: scan.one_based_scan_number,
                zero_based_scan_index: scan_index as i32,
                retention_time: scan.retention_time,
                msn_order: scan.msn_order,
            });

            for j in 0..scan.mz.len() {
                let mz = scan.mz[j];
                // C# (int)Math.Round(mz * BinsPerDalton, 0) — banker's rounding.
                let rounded_mz = (mz * BINS_PER_DALTON).round_ties_even() as usize;
                indexed_peaks[rounded_mz]
                    .get_or_insert_with(Vec::new)
                    .push(IndexedMassSpectralPeak::new(
                        mz,
                        scan.intensity[j],
                        scan_index as i32,
                        scan.retention_time,
                    ));
            }
        }

        Some(PeakIndexingEngine {
            indexed_peaks,
            scan_info,
        })
    }

    /// Builds the index directly from the MS1 scans of a spectra file.
    ///
    /// Mirrors `PeakIndexingEngine.InitializeIndexingEngine(MsDataFile)`: read all MS1
    /// (survey) scans, ordered by scan number, and index their peaks. Returns `None` when
    /// there is nothing to index.
    ///
    /// The format (mzML or Thermo `.raw`) is detected from the file by [`read_ms1_scans`];
    /// `.raw` reading requires a .NET 8 runtime (see [`read_ms1_scans`]).
    pub fn from_spectra_file<P: Into<PathBuf> + Clone>(
        path: P,
    ) -> std::io::Result<Option<PeakIndexingEngine>> {
        let scans = read_ms1_scans(path)?;
        Ok(Self::index_peaks(&scans))
    }

    /// Backwards-compatible alias for [`from_spectra_file`](Self::from_spectra_file).
    ///
    /// Kept for the P1.x callers and the Python `PeakIndex.from_mzml` binding; despite the name
    /// it now reads any format [`read_ms1_scans`] supports (mzML and Thermo `.raw`).
    pub fn from_mzml<P: Into<PathBuf> + Clone>(path: P) -> std::io::Result<Option<PeakIndexingEngine>> {
        Self::from_spectra_file(path)
    }

    /// Per-scan metadata array, indexed by zero-based scan index.
    pub fn scan_info(&self) -> &[ScanInfo] {
        &self.scan_info
    }

    /// Every indexed peak, flattened in bin-then-scan order.
    ///
    /// The order matches `get_all_xics`' pre-sort traversal (ascending bin, then scan order
    /// within a bin). Callers that want intensity-descending seeds sort the result themselves.
    /// Exposed for the trace-kernel detector, which seeds from the full peak set.
    pub fn all_peaks(&self) -> Vec<IndexedMassSpectralPeak> {
        self.indexed_peaks
            .iter()
            .flatten()
            .flat_map(|bin| bin.iter().copied())
            .collect()
    }

    /// Finds the peak closest to `m` in scan `zero_based_scan_index` and within `ppm`
    /// tolerance. Returns `None` if no peak qualifies.
    ///
    /// Faithful port of `IndexingEngine<T>.GetIndexedPeak` (the m/z, charge-less path).
    pub fn get_indexed_peak(
        &self,
        m: f64,
        zero_based_scan_index: i32,
        ppm: &PpmTolerance,
    ) -> Option<&IndexedMassSpectralPeak> {
        let bins = self.get_bins_in_range(m, ppm);
        if bins.is_empty() {
            return None;
        }
        let peak_indices: Vec<isize> = bins
            .iter()
            .map(|b| Self::binary_search_for_indexed_peak(b, zero_based_scan_index))
            .collect();
        Self::get_best_peak_from_bins(&bins, m, zero_based_scan_index, &peak_indices, ppm)
    }

    /// Traces a peak of m/z `m` across retention time, beginning at the scan just before
    /// `retention_time`. Faithful port of `IndexingEngine<T>.GetXic`.
    ///
    /// Finds the start scan by walking `scan_info` forward while RT is strictly below the target
    /// (so it stops on the first scan whose RT is `>= retention_time`), then delegates to
    /// [`Self::get_xic_by_scan_index`]. If every scan's RT is `>= retention_time` the start index
    /// stays `-1`, exactly as in C#; the forward walk then picks up scan 0 onward.
    #[allow(clippy::too_many_arguments)]
    pub fn get_xic(
        &self,
        m: f64,
        retention_time: f64,
        ppm: &PpmTolerance,
        missed_scans_allowed: i32,
        max_peak_half_width: f64,
        matched_peaks: Option<&HashSet<PeakKey>>,
    ) -> Vec<IndexedMassSpectralPeak> {
        let mut scan_index: i32 = -1;
        for scan in &self.scan_info {
            if scan.retention_time < retention_time {
                scan_index = scan.zero_based_scan_index;
            } else {
                break;
            }
        }
        self.get_xic_by_scan_index(
            m,
            scan_index,
            ppm,
            missed_scans_allowed,
            max_peak_half_width,
            matched_peaks,
        )
    }

    /// Traces a peak of m/z `m` across retention time, beginning at `zero_based_start_index`.
    /// Faithful port of `IndexingEngine<T>.GetXicByScanIndex` (the m/z, charge-less path).
    ///
    /// Walks backward (`direction = -1`) then forward (`direction = +1`) from the start scan,
    /// following a per-bin pointer array so each step lands on the first peak of the new scan
    /// index. The walk stops once `missed_scans_allowed` consecutive misses accumulate, the scan
    /// range is exhausted, or a scan's RT is more than `max_peak_half_width` from the initial
    /// peak's RT. A peak already present in `matched_peaks` counts as a miss (so `GetAllXics` does
    /// not re-claim it). The result is sorted by retention time ascending.
    #[allow(clippy::too_many_arguments)]
    pub fn get_xic_by_scan_index(
        &self,
        m: f64,
        zero_based_start_index: i32,
        ppm: &PpmTolerance,
        missed_scans_allowed: i32,
        max_peak_half_width: f64,
        matched_peaks: Option<&HashSet<PeakKey>>,
    ) -> Vec<IndexedMassSpectralPeak> {
        let mut xic: Vec<IndexedMassSpectralPeak> = Vec::new();
        let all_bins = self.get_bins_in_range(m, ppm);
        if all_bins.is_empty() {
            return xic;
        }

        // A pointer into each bin, seeded at the start scan; copied per direction below.
        let peak_pointer_array: Vec<isize> = all_bins
            .iter()
            .map(|b| Self::binary_search_for_indexed_peak(b, zero_based_start_index))
            .collect();

        let mut initial_peak =
            Self::get_best_peak_from_bins(&all_bins, m, zero_based_start_index, &peak_pointer_array, ppm)
                .copied();

        if let Some(p) = initial_peak {
            xic.push(p);
        }

        for &direction in &[-1i32, 1i32] {
            // The legacy code does not count the (possibly missing) initial peak as a miss.
            let mut missed_peaks = 0;
            let mut current = zero_based_start_index;
            let mut pointer_copy = peak_pointer_array.clone();

            while missed_peaks <= missed_scans_allowed {
                current += direction;
                if current < 0
                    || current > self.scan_info.len() as i32 - 1
                    || (initial_peak.is_some()
                        && (self.scan_info[current as usize].retention_time
                            - initial_peak.unwrap().retention_time as f64)
                            .abs()
                            > max_peak_half_width)
                {
                    break;
                }

                // Advance every per-bin pointer to the first peak of `current`'s scan index.
                for i in 0..pointer_copy.len() {
                    match direction {
                        -1 => {
                            loop {
                                pointer_copy[i] -= 1;
                                let keep_going = pointer_copy[i] >= 0
                                    && all_bins[i][pointer_copy[i] as usize].zero_based_scan_index
                                        > current - 1;
                                if !keep_going {
                                    break;
                                }
                            }
                            pointer_copy[i] += 1;
                        }
                        1 => {
                            while pointer_copy[i] < all_bins[i].len() as isize - 1
                                && all_bins[i][pointer_copy[i] as usize].zero_based_scan_index
                                    < current
                            {
                                pointer_copy[i] += 1;
                            }
                        }
                        _ => {}
                    }
                }

                let next_peak =
                    Self::get_best_peak_from_bins(&all_bins, m, current, &pointer_copy, ppm).copied();

                let claimed = match next_peak {
                    Some(p) => matched_peaks.map_or(false, |set| set.contains(&peak_key(&p))),
                    None => false,
                };

                match next_peak {
                    Some(p) if !claimed => {
                        if initial_peak.is_none() {
                            initial_peak = Some(p);
                        }
                        xic.push(p);
                        missed_peaks = 0;
                    }
                    _ => missed_peaks += 1,
                }
            }
        }

        xic.sort_by(|a, b| a.retention_time.total_cmp(&b.retention_time));
        xic
    }

    /// Traces every peak in the index into an XIC, greedily from most intense to least, skipping
    /// peaks already claimed by an earlier XIC. Faithful port of `IndexingEngine<T>.GetAllXics`
    /// **minus** the optional `CutPeak` step, which depends on `ExtractedIonChromatogram.CutPeak`
    /// (ported in P1.14). Only XICs of at least `num_peak_threshold` peaks are kept.
    ///
    /// Peaks are visited in descending intensity order (a stable sort over the index's
    /// bin-then-scan traversal order, matching the C# `OrderByDescending`). Each surviving XIC's
    /// peaks are marked claimed so later iterations skip them.
    pub fn get_all_xics(
        &self,
        peak_finding_tolerance: &PpmTolerance,
        max_missed_scan_allowed: i32,
        max_rt_range: f64,
        num_peak_threshold: usize,
    ) -> Vec<ExtractedIonChromatogram> {
        let mut xics: Vec<ExtractedIonChromatogram> = Vec::new();
        let mut matched_peaks: HashSet<PeakKey> = HashSet::new();

        // Flatten the jagged index in bin-then-scan order, then stable-sort by intensity desc.
        let mut sorted_peaks: Vec<IndexedMassSpectralPeak> = self
            .indexed_peaks
            .iter()
            .flatten()
            .flat_map(|bin| bin.iter().copied())
            .collect();
        sorted_peaks.sort_by(|a, b| b.intensity.total_cmp(&a.intensity));

        for peak in &sorted_peaks {
            if matched_peaks.contains(&peak_key(peak)) {
                continue;
            }
            let peak_list = self.get_xic_by_scan_index(
                peak.m() as f64,
                peak.zero_based_scan_index,
                peak_finding_tolerance,
                max_missed_scan_allowed,
                max_rt_range,
                Some(&matched_peaks),
            );
            if peak_list.len() >= num_peak_threshold {
                let xic = ExtractedIonChromatogram::new(peak_list);
                for matched_peak in &xic.peaks {
                    matched_peaks.insert(peak_key(matched_peak));
                }
                xics.push(xic);
            } else {
                matched_peaks.insert(peak_key(peak));
            }
        }

        xics
    }

    /// Collects every non-empty bin overlapping `mz ± ppm`.
    ///
    /// Faithful port of `GetBinsInRange`: `floor(min*100) ..= ceil(max*100)`, skipping
    /// out-of-range and empty bins.
    fn get_bins_in_range(&self, mz: f64, ppm: &PpmTolerance) -> Vec<&Vec<IndexedMassSpectralPeak>> {
        let ceiling_mz = (ppm.get_maximum_value(mz) * BINS_PER_DALTON).ceil() as i64;
        let floor_mz = (ppm.get_minimum_value(mz) * BINS_PER_DALTON).floor() as i64;
        let mut all_bins = Vec::new();
        for j in floor_mz..=ceiling_mz {
            if j < 0 || j as usize >= self.indexed_peaks.len() {
                continue;
            }
            if let Some(bin) = &self.indexed_peaks[j as usize] {
                all_bins.push(bin);
            }
        }
        all_bins
    }

    /// Picks the peak closest to `mz` across all candidate bins.
    ///
    /// Faithful port of `GetBestPeakFromBins` (m/z path; no charge).
    fn get_best_peak_from_bins<'a>(
        all_bins: &[&'a Vec<IndexedMassSpectralPeak>],
        mz: f64,
        zero_based_scan_index: i32,
        peak_indices_in_bins: &[isize],
        ppm: &PpmTolerance,
    ) -> Option<&'a IndexedMassSpectralPeak> {
        let mut best_peak: Option<&IndexedMassSpectralPeak> = None;
        for i in 0..all_bins.len() {
            let temp_peak = Self::get_peak_from_bin(
                all_bins[i],
                mz,
                zero_based_scan_index,
                peak_indices_in_bins[i],
                ppm,
            );
            let temp_peak = match temp_peak {
                Some(p) => p,
                None => continue,
            };
            match best_peak {
                None => best_peak = Some(temp_peak),
                Some(b) => {
                    if (temp_peak.m() as f64 - mz).abs() < (b.m() as f64 - mz).abs() {
                        best_peak = Some(temp_peak);
                    }
                }
            }
        }
        best_peak
    }

    /// Picks the peak closest to `mz` within a single bin, starting from `peak_index_in_bin`
    /// and scanning forward while the scan index matches.
    ///
    /// Faithful port of `GetPeakFromBin` (m/z path; no charge).
    fn get_peak_from_bin<'a>(
        bin: &'a [IndexedMassSpectralPeak],
        mz: f64,
        zero_based_scan_index: i32,
        peak_index_in_bin: isize,
        ppm: &PpmTolerance,
    ) -> Option<&'a IndexedMassSpectralPeak> {
        if peak_index_in_bin < 0 || peak_index_in_bin as usize >= bin.len() {
            return None;
        }
        let mut best_peak: Option<&IndexedMassSpectralPeak> = None;
        for i in (peak_index_in_bin as usize)..bin.len() {
            let peak = &bin[i];

            if peak.zero_based_scan_index > zero_based_scan_index {
                break;
            }

            if ppm.within(peak.m() as f64, mz)
                && peak.zero_based_scan_index == zero_based_scan_index
                && best_peak.map_or(true, |b| (peak.m() as f64 - mz).abs() < (b.m() as f64 - mz).abs())
            {
                best_peak = Some(peak);
            }
        }
        best_peak
    }

    /// Finds the index of the first peak in `indexed_peaks` whose scan index is
    /// `>= zero_based_scan_index` (a binary search refined by a short linear backstep).
    ///
    /// Faithful port of `BinarySearchForIndexedPeak`. Returns an `isize` so the caller can
    /// distinguish "before the list" the same way the C# integer return (clamped to 0) does.
    fn binary_search_for_indexed_peak(
        indexed_peaks: &[IndexedMassSpectralPeak],
        zero_based_scan_index: i32,
    ) -> isize {
        let mut m: isize = 0;
        let mut l: isize = 0;
        let mut r: isize = indexed_peaks.len() as isize - 1;

        while l <= r {
            m = l + ((r - l) / 2);

            if r - l < 2 {
                break;
            }
            if indexed_peaks[m as usize].zero_based_scan_index < zero_based_scan_index {
                l = m + 1;
            } else {
                r = m - 1;
            }
        }

        let mut i = m;
        while i >= 0 {
            if indexed_peaks[i as usize].zero_based_scan_index < zero_based_scan_index {
                break;
            }
            m -= 1;
            i -= 1;
        }

        if m < 0 {
            m = 0;
        }

        m
    }
}

/// An extracted-ion chromatogram: a peak trace across retention time. Partial port of mzLib
/// `MassSpectrometry.ExtractedIonChromatogram` — holds the peaks (RT-ascending) plus the summary
/// info `SetXicInfo` derives. The `CutPeak` boundary trimming is deferred to P1.14.
#[derive(Debug, Clone)]
pub struct ExtractedIonChromatogram {
    /// The traced peaks, ordered by retention time ascending.
    pub peaks: Vec<IndexedMassSpectralPeak>,
    /// The most intense peak (the apex).
    pub apex_peak: IndexedMassSpectralPeak,
    /// Minimum retention time across the peaks.
    pub start_rt: f64,
    /// Maximum retention time across the peaks.
    pub end_rt: f64,
    /// Minimum zero-based scan index across the peaks.
    pub start_scan_index: i32,
    /// Maximum zero-based scan index across the peaks.
    pub end_scan_index: i32,
    /// Intensity-weighted average m/z (C# `AverageM`/`AveragedMassOrMz`).
    pub averaged_mass_or_mz: f64,
}

impl ExtractedIonChromatogram {
    /// Builds an XIC from a non-empty peak list, computing the summary info (C# constructor +
    /// `SetXicInfo`). The apex is the max-intensity peak; RT/scan bounds are the min/max over the
    /// peaks; `averaged_mass_or_mz` is the intensity-weighted mean m/z.
    pub fn new(peaks: Vec<IndexedMassSpectralPeak>) -> ExtractedIonChromatogram {
        assert!(
            !peaks.is_empty(),
            "ExtractedIonChromatogram requires at least one peak"
        );
        let apex_peak = *peaks
            .iter()
            .max_by(|a, b| a.intensity.total_cmp(&b.intensity))
            .unwrap();
        let start_rt = peaks
            .iter()
            .map(|p| p.retention_time as f64)
            .fold(f64::INFINITY, f64::min);
        let end_rt = peaks
            .iter()
            .map(|p| p.retention_time as f64)
            .fold(f64::NEG_INFINITY, f64::max);
        let start_scan_index = peaks.iter().map(|p| p.zero_based_scan_index).min().unwrap();
        let end_scan_index = peaks.iter().map(|p| p.zero_based_scan_index).max().unwrap();
        let averaged_mass_or_mz = Self::average_m(&peaks);
        ExtractedIonChromatogram {
            peaks,
            apex_peak,
            start_rt,
            end_rt,
            start_scan_index,
            end_scan_index,
            averaged_mass_or_mz,
        }
    }

    /// Intensity-weighted average m/z over the peaks. Port of C# `AverageM`.
    fn average_m(peaks: &[IndexedMassSpectralPeak]) -> f64 {
        let sum_intensity: f64 = peaks.iter().map(|p| p.intensity as f64).sum();
        let mut averaged_m = 0.0;
        for peak in peaks {
            let weight = peak.intensity as f64 / sum_intensity;
            averaged_m += weight * peak.m() as f64;
        }
        averaged_m
    }

    /// Apex retention time (C# `ApexRT`).
    #[inline]
    pub fn apex_rt(&self) -> f64 {
        self.apex_peak.retention_time as f64
    }

    /// Apex zero-based scan index (C# `ApexScanIndex`).
    #[inline]
    pub fn apex_scan_index(&self) -> i32 {
        self.apex_peak.zero_based_scan_index
    }
}

/// Reads the MS1 scans of a spectra file into the `Scan` shape the index consumes.
///
/// Mirrors `PeakIndexingEngine.InitializeIndexingEngine(MsDataFile)`: keep only spectra whose MS
/// level is 1 (survey scans), in file order (ascending scan number). The one-based scan number is
/// derived from the file index, and the retention time from `start_time()` (minutes).
///
/// **Thermo `.raw` is centroided on read** using the vendor (Thermo) algorithm. Orbitrap MS1 survey
/// scans are acquired in *profile* mode (~9k sample points per scan), but the trace-kernel detector's
/// isotope-comb model assumes **one peak per isotope**. Feeding it raw profile points inflates the
/// peak count ~13×, samples each isotope off-apex, and double-counts profile shoulders in the summed
/// intensity. So the Thermo reader is opened with `centroiding = true`, which returns the vendor
/// centroid stream — the same peaks MetaMorpheus / FlashLFQ consume. mzML is read as stored (the
/// calibrated Lumos mzML is already centroided; a profile mzML would come through as profile).
///
/// Peaks are extracted via [`SpectrumLike::peaks`], which yields the centroid list for a centroided
/// spectrum and the raw profile arrays otherwise. Its iteration order is **not** guaranteed
/// m/z-ascending, so the pairs are sorted (the invariant [`PeakIndexingEngine::index_peaks`] and
/// [`crate::feature_refinement::refine_feature`] rely on) before the zero-intensity filter.
///
/// `.raw` support rides the `thermorawfilereader` crate, which hosts a self-contained **.NET 8
/// runtime** — reading a `.raw` therefore requires a .NET 8 runtime on the machine (mzML stays
/// pure-Rust).
pub fn read_ms1_scans<P: Into<PathBuf> + Clone>(path: P) -> std::io::Result<Vec<Scan>> {
    let path: PathBuf = path.into();
    let is_thermo_raw = path
        .extension()
        .map(|e| e.eq_ignore_ascii_case("raw"))
        .unwrap_or(false);

    if is_thermo_raw {
        // Vendor-centroid the profile Thermo scans on read (see doc above).
        let reader =
            ThermoRawReader::new_with_detail_level_and_centroiding(path, DetailLevel::Full, true)?;
        Ok(collect_ms1_scans(reader))
    } else {
        let reader = MZReader::open_path(path)?;
        Ok(collect_ms1_scans(reader))
    }
}

/// Extracts MS1 scans from any spectrum iterator (mzML or Thermo). Shared by both branches of
/// [`read_ms1_scans`]; see its docs for the centroiding and ordering rationale.
fn collect_ms1_scans<I: Iterator<Item = Spectrum>>(reader: I) -> Vec<Scan> {
    let mut scans = Vec::new();
    for spectrum in reader {
        if spectrum.ms_level() != 1 {
            continue;
        }
        // `peaks()` yields the centroid list for a centroided spectrum, raw profile points otherwise.
        let peaks = spectrum.peaks();
        let mut pairs: Vec<(f64, f64)> = Vec::with_capacity(peaks.len());
        for p in peaks.iter() {
            pairs.push((p.mz, p.intensity as f64));
        }
        // Not guaranteed m/z-ordered; sort to satisfy the ascending-m/z invariant of the index.
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut mz = Vec::with_capacity(pairs.len());
        let mut intensity = Vec::with_capacity(pairs.len());
        for (m, i) in pairs {
            mz.push(m);
            intensity.push(i);
        }
        let (mz, intensity) = remove_zero_intensity_peaks(mz, intensity);
        scans.push(Scan {
            mz,
            intensity,
            one_based_scan_number: spectrum.index() as i32 + 1,
            retention_time: spectrum.start_time(),
            msn_order: spectrum.ms_level() as i32,
        });
    }
    scans
}

/// Intensity below which mzLib's mzML reader treats a peak as "zero" and drops it.
/// Verbatim from `Readers/MzML/Mzml.cs` (`zeroEquivalentIntensity`).
const ZERO_EQUIVALENT_INTENSITY: f64 = 0.01;

/// Drops peaks whose intensity is below [`ZERO_EQUIVALENT_INTENSITY`], faithfully mirroring the
/// "Remove Zero Intensity Peaks" step of mzLib's MzML reader (`Mzml.cs`).
///
/// This matters for **profile** data: the raw mzML carries many zero-intensity sample points, and
/// `mzdata` returns them verbatim. mzLib strips them on read, which both shrinks the index and —
/// critically — makes XIC tracing stop at real gaps (a dropped zero is a genuine missed scan, not
/// a spurious in-tolerance hit that keeps the walk alive). Without this filter the Rust XIC walks
/// straight through zero-filled gaps and picks up co-eluting interference peaks the C# engine never
/// sees, diverging on the L5/L6 parity gate.
///
/// Matches mzLib's exact guard: peaks are dropped only when at least one, but not all, fall below
/// the threshold (an all-"zero" scan is kept verbatim, as in C#). Input arrays are assumed parallel
/// and mass-ascending; the relative order of surviving peaks is preserved.
fn remove_zero_intensity_peaks(mz: Vec<f64>, intensity: Vec<f64>) -> (Vec<f64>, Vec<f64>) {
    let total = intensity.len();
    let zero_count = intensity
        .iter()
        .filter(|&&v| v < ZERO_EQUIVALENT_INTENSITY)
        .count();
    if zero_count == 0 || zero_count == total {
        return (mz, intensity);
    }
    let mut kept_mz = Vec::with_capacity(total - zero_count);
    let mut kept_intensity = Vec::with_capacity(total - zero_count);
    for (m, i) in mz.into_iter().zip(intensity.into_iter()) {
        if i >= ZERO_EQUIVALENT_INTENSITY {
            kept_mz.push(m);
            kept_intensity.push(i);
        }
    }
    (kept_mz, kept_intensity)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Builds the synthetic 9-scan fixture from the C# `IndexingEngineTests` setup: five
    /// peaks at 500/500.5/501/501.5/502 in every scan, intensities scaled per scan.
    fn synthetic_scans() -> Vec<Scan> {
        let intensity = 1e6;
        let multipliers = [1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0];
        let mut scans = Vec::new();
        for (s, mult) in multipliers.iter().enumerate() {
            scans.push(Scan {
                mz: vec![500.0, 500.5, 501.0, 501.5, 502.0],
                intensity: vec![mult * intensity; 5],
                one_based_scan_number: (s + 1) as i32,
                retention_time: 1.0 + s as f64 / 10.0,
                msn_order: 1,
            });
        }
        scans
    }

    #[test]
    fn bins_per_dalton_is_100() {
        assert_eq!(BINS_PER_DALTON, 100.0);
    }

    #[test]
    fn index_peaks_returns_none_on_empty_input() {
        assert!(PeakIndexingEngine::index_peaks(&[]).is_none());
        // a scan with no peaks anywhere -> nothing to index
        let scans = vec![Scan {
            mz: vec![],
            intensity: vec![],
            one_based_scan_number: 1,
            retention_time: 1.0,
            msn_order: 1,
        }];
        assert!(PeakIndexingEngine::index_peaks(&scans).is_none());
    }

    #[test]
    fn rounds_mz_into_expected_bins() {
        // Each m/z lands in bin round(mz * 100). m/z must be ascending so the last value
        // is the max (the C# `LastX` assumption that sizes the index). The banker's-rounding
        // tie behaviour itself is covered by the round_ties_even test in
        // theoretical_isotope_distribution; here we check plain bin assignment.
        let scans = vec![Scan {
            mz: vec![500.4, 500.5, 500.51],
            intensity: vec![1.0, 2.0, 3.0],
            one_based_scan_number: 1,
            retention_time: 1.0,
            msn_order: 1,
        }];
        let engine = PeakIndexingEngine::index_peaks(&scans).expect("indexed");
        assert!(engine.indexed_peaks[50040].is_some()); // 500.4  -> 50040
        assert!(engine.indexed_peaks[50050].is_some()); // 500.5  -> 50050
        assert!(engine.indexed_peaks[50051].is_some()); // 500.51 -> 50051
    }

    #[test]
    fn get_indexed_peak_finds_exact_mz_in_scan() {
        // Mirrors C# TestGetIndexedPeakWithoutChargeState: m=500, scan 5, 20 ppm.
        let engine = PeakIndexingEngine::index_peaks(&synthetic_scans()).expect("indexed");
        let peak = engine
            .get_indexed_peak(500.0, 5, &PpmTolerance::new(20.0))
            .expect("peak should be found");
        assert!((peak.m() as f64 - 500.0).abs() < 0.001);
        assert_eq!(peak.zero_based_scan_index, 5);
    }

    #[test]
    fn get_indexed_peak_resolves_closest_of_adjacent_peaks() {
        // 500 and 500.5 are both in range of a wide tolerance around 500.2; the closer one
        // (500.0, 400 ppm away) must win over 500.5 (600 ppm away).
        let engine = PeakIndexingEngine::index_peaks(&synthetic_scans()).expect("indexed");
        let peak = engine
            .get_indexed_peak(500.2, 3, &PpmTolerance::new(1000.0))
            .expect("peak should be found");
        assert!((peak.m() as f64 - 500.0).abs() < 0.001);
        assert_eq!(peak.zero_based_scan_index, 3);
    }

    #[test]
    fn get_indexed_peak_returns_none_when_out_of_range() {
        // 400 has no peak anywhere (mirrors C# TestMissingXic's empty result).
        let engine = PeakIndexingEngine::index_peaks(&synthetic_scans()).expect("indexed");
        assert!(engine
            .get_indexed_peak(400.0, 0, &PpmTolerance::new(20.0))
            .is_none());
    }

    #[test]
    fn get_indexed_peak_returns_none_when_scan_absent() {
        // A peak at 500 exists, but not in scan 100 (beyond the 9 scans).
        let engine = PeakIndexingEngine::index_peaks(&synthetic_scans()).expect("indexed");
        assert!(engine
            .get_indexed_peak(500.0, 100, &PpmTolerance::new(20.0))
            .is_none());
    }

    #[test]
    fn get_indexed_peak_finds_peak_in_every_scan() {
        let engine = PeakIndexingEngine::index_peaks(&synthetic_scans()).expect("indexed");
        for scan in 0..9 {
            let peak = engine
                .get_indexed_peak(501.0, scan, &PpmTolerance::new(20.0))
                .unwrap_or_else(|| panic!("peak missing in scan {scan}"));
            assert_eq!(peak.zero_based_scan_index, scan);
            assert!((peak.m() as f64 - 501.0).abs() < 0.001);
        }
    }

    /// Resolves a path under the mzLib FlashLFQ golden test data directory
    /// (see [`crate::mzlib_test_data`]).
    fn test_data(relative: &str) -> PathBuf {
        crate::mzlib_test_data(relative)
    }

    #[test]
    fn indexes_real_mzml_and_round_trips_a_point_query() {
        // Build the index straight from the golden mzML (10 MS1 scans, per P0.2), then
        // pick a real peak out of the index and query it back by its own m/z + scan.
        let engine = PeakIndexingEngine::from_mzml(test_data("sliced-mzml.mzML"))
            .expect("sliced-mzml.mzML should be readable")
            .expect("the file has indexable MS1 peaks");

        // 10 MS1 scans were indexed (matches count_ms1_scans in lib.rs).
        assert_eq!(engine.scan_info().len(), 10);

        // Grab the first peak from the first non-empty bin as a known-present target.
        let target = engine
            .indexed_peaks
            .iter()
            .flatten()
            .flat_map(|bin| bin.iter())
            .next()
            .copied()
            .expect("at least one indexed peak");

        let found = engine
            .get_indexed_peak(
                target.m() as f64,
                target.zero_based_scan_index,
                &PpmTolerance::new(5.0),
            )
            .expect("the peak we pulled from the index must be queryable");

        // The query should return the same scan and an m/z within tolerance.
        assert_eq!(found.zero_based_scan_index, target.zero_based_scan_index);
        assert!((found.m() as f64 - target.m() as f64).abs() < 1e-3);
    }

    /// Builds a 10-scan fixture mirroring the C# `TestXicStops` shape: a peak at m/z 500 in every
    /// scan **except** scans 2 and 3, where it sits at m/z 250 (a charge-2 shift, out of the 20 ppm
    /// window around 500). Tracing 500 must therefore stop after the two-scan gap.
    fn xic_stops_scans() -> Vec<Scan> {
        let intensity = 1e6;
        let multipliers = [1.0, 3.0, 0.0, 0.0, 3.0, 5.0, 10.0, 5.0, 3.0, 1.0];
        let mut scans = Vec::new();
        for (s, mult) in multipliers.iter().enumerate() {
            // Scans 2 and 3 carry the peak at a shifted m/z so it is invisible at 500 ± 20 ppm.
            let mz = if s == 2 || s == 3 { 250.0 } else { 500.0 };
            scans.push(Scan {
                mz: vec![mz],
                intensity: vec![mult * intensity + 1.0], // +1 so scans 2/3 still hold a (250) peak
                one_based_scan_number: (s + 1) as i32,
                retention_time: 1.0 + s as f64 / 10.0,
                msn_order: 1,
            });
        }
        scans
    }

    #[test]
    fn get_xic_by_scan_index_traces_all_scans() {
        // Mirrors C# TestGetXicWithStartIndex: GetXicByScanIndex(500, 0, 20 ppm) -> 9 peaks.
        let engine = PeakIndexingEngine::index_peaks(&synthetic_scans()).expect("indexed");
        let xic = engine.get_xic_by_scan_index(500.0, 0, &PpmTolerance::new(20.0), 1, f64::MAX, None);
        assert_eq!(xic.len(), 9);
        // Sorted by RT ascending, covering scans 0..8.
        for (i, peak) in xic.iter().enumerate() {
            assert_eq!(peak.zero_based_scan_index, i as i32);
            assert!((peak.m() as f64 - 500.0).abs() < 0.001);
        }
    }

    #[test]
    fn get_xic_by_retention_time_traces_all_scans() {
        // Mirrors C# TestGetXicWithRetentionTime: GetXic(500, 1.0, 20 ppm) -> 9. RT 1.0 is scan 0's
        // RT, so the start index stays -1 and the forward walk picks up scans 0..8.
        let engine = PeakIndexingEngine::index_peaks(&synthetic_scans()).expect("indexed");
        let xic = engine.get_xic(500.0, 1.0, &PpmTolerance::new(20.0), 1, f64::MAX, None);
        assert_eq!(xic.len(), 9);
    }

    #[test]
    fn get_xic_returns_empty_for_absent_mz() {
        // Mirrors C# TestMissingXic: GetXicByScanIndex(400, 0, 20 ppm) -> empty.
        let engine = PeakIndexingEngine::index_peaks(&synthetic_scans()).expect("indexed");
        let xic = engine.get_xic_by_scan_index(400.0, 0, &PpmTolerance::new(20.0), 1, f64::MAX, None);
        assert!(xic.is_empty());
    }

    #[test]
    fn get_xic_stops_after_missed_scans() {
        // Mirrors C# TestXicStops: starting at scan 7, tracing 500 with 1 missed scan allowed
        // yields 6 peaks (scans 4,5,6,7,8,9) — the scan 2/3 gap halts the backward walk.
        let engine = PeakIndexingEngine::index_peaks(&xic_stops_scans()).expect("indexed");
        let xic = engine.get_xic_by_scan_index(500.0, 7, &PpmTolerance::new(20.0), 1, f64::MAX, None);
        assert_eq!(xic.len(), 6);
        let scans: Vec<i32> = xic.iter().map(|p| p.zero_based_scan_index).collect();
        assert_eq!(scans, vec![4, 5, 6, 7, 8, 9]);
    }

    #[test]
    fn get_all_xics_traces_every_mz_trace() {
        // The synthetic fixture has 5 distinct m/z (500/500.5/501/501.5/502) across 9 scans.
        // At 20 ppm each is its own bin, so GetAllXics yields 5 XICs of 9 peaks each (45 total).
        let engine = PeakIndexingEngine::index_peaks(&synthetic_scans()).expect("indexed");
        let xics = engine.get_all_xics(&PpmTolerance::new(20.0), 1, f64::MAX, 2);
        assert_eq!(xics.len(), 5);
        let total: usize = xics.iter().map(|x| x.peaks.len()).sum();
        assert_eq!(total, 45);
        for xic in &xics {
            assert_eq!(xic.peaks.len(), 9);
            // The apex is scan 4 (intensity multiplier 5, the max).
            assert_eq!(xic.apex_scan_index(), 4);
        }
    }

    #[test]
    fn scan_info_tracks_each_scan() {
        let engine = PeakIndexingEngine::index_peaks(&synthetic_scans()).expect("indexed");
        let info = engine.scan_info();
        assert_eq!(info.len(), 9);
        assert_eq!(info[0].zero_based_scan_index, 0);
        assert_eq!(info[0].one_based_scan_number, 1);
        assert_eq!(info[8].zero_based_scan_index, 8);
        assert!((info[5].retention_time - 1.5).abs() < 1e-12);
    }
}
