//! `flashlfq-core` — pure-Rust port of the FlashLFQ quantification core.
//!
//! This crate knows nothing about Python. All algorithm logic, unit tests, and the
//! golden-file parity harness live here so they run with zero Python in the loop
//! (see `agent_info/FlashLFQ-Rust-Rewrite-Feasibility.md`, "Python binding layer").
//!
//! At this stage the chemistry layers (L0–L3) are not yet ported; the crate carries the
//! Phase 0 binding/reader spike. mzML reading is provided by `mzdata` (P0.2).

pub mod chemical_formula;
pub mod chromatographic_peak;
pub mod deconvolution;
pub mod detection_type;
pub mod engine;
pub mod feature_refinement;
pub mod isotope_shift_decon;
pub mod isotopic_distribution;
pub mod joint_fit;
pub mod isotopic_envelope;
pub mod mbr;
pub mod mbr_chromatographic_peak;
pub mod mbr_scorer;
pub mod mbr_search;
pub mod parquet_output;
pub mod peak_indexing;
pub mod peptide;
pub mod periodic_table;
pub mod psm_tsv;
pub mod results;
pub mod special_functions;
pub mod spectral_averaging;
pub mod theoretical_isotope_distribution;
pub mod tolerance;
pub mod trace_kernel;

use std::path::PathBuf;
use std::sync::Arc;

use arrow::array::Float64Array;
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use mzdata::io::MzMLReader;
use mzdata::prelude::*;

/// Returns the crate's semantic version string. Placeholder used by the Phase 0
/// binding spike to confirm the core links and is callable.
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// Root of the external **mzLib** checkout this port references (but does not vendor).
///
/// The port and mzLib are separate projects; parity gates and golden regeneration read
/// mzLib's test data and C# source live from here. The location is read from the
/// `MZLIB_DIR` environment variable, falling back to `F:\mzLib` (the standard checkout on
/// the development machine). Set `MZLIB_DIR` to your own mzLib checkout on other machines.
/// See the project README for details.
pub fn mzlib_dir() -> PathBuf {
    match std::env::var_os("MZLIB_DIR") {
        Some(v) => PathBuf::from(v),
        None => PathBuf::from(r"F:\mzLib"),
    }
}

/// Resolves a path under mzLib's FlashLFQ golden test-data directory,
/// i.e. `<MZLIB_DIR>/mzLib/Test/FlashLFQ/TestData/<relative>`.
///
/// This is the single source of truth for locating the shared FlashLFQ fixtures; every
/// parity test and example resolves fixtures through it so the mzLib reference lives in
/// exactly one place (see [`mzlib_dir`]).
pub fn mzlib_test_data(relative: impl AsRef<std::path::Path>) -> PathBuf {
    mzlib_dir()
        .join("mzLib")
        .join("Test")
        .join("FlashLFQ")
        .join("TestData")
        .join(relative)
}

/// Produces a contiguous `Vec<f64>` of length `n` holding `[0.0, 1.0, ..., (n-1)]`.
///
/// Phase 0 binding spike (P0.3): a pure-Rust array producer with no Python in sight.
/// The bindings crate moves this `Vec` into a NumPy array via `into_pyarray` (the
/// no-copy path), so the core stays plain and independently testable. A deterministic
/// ramp makes the length and contents trivial to assert on either side of the FFI.
pub fn iota_f64(n: usize) -> Vec<f64> {
    (0..n).map(|i| i as f64).collect()
}

/// Builds a small two-column Arrow `RecordBatch` `(mz: f64, intensity: f64)`.
///
/// Phase 0 Arrow spike (P0.4): a pure-Rust producer of a columnar Arrow table with no
/// Python in sight. The bindings crate hands this batch to pyarrow over the Arrow C Data
/// Interface (zero-copy), so the core stays plain and independently testable. The shape
/// (a peak-list–like `mz`/`intensity` pair) previews the real quant output columns; the
/// fixed contents make the row/column counts and values trivial to assert on either side
/// of the FFI. Building the columnar arrays is the only real cost — the handoff itself is
/// free (see feasibility doc, "Zero-copy Arrow is doing some lifting").
pub fn demo_record_batch() -> RecordBatch {
    let mz = Float64Array::from(vec![100.0, 200.5, 300.25]);
    let intensity = Float64Array::from(vec![1.0e6, 2.5e6, 7.0e5]);
    let schema = Arc::new(Schema::new(vec![
        Field::new("mz", DataType::Float64, false),
        Field::new("intensity", DataType::Float64, false),
    ]));
    // try_new only fails on schema/column mismatch; our literals are statically consistent.
    RecordBatch::try_new(schema, vec![Arc::new(mz), Arc::new(intensity)])
        .expect("demo schema and columns are statically consistent")
}

/// Opens an mzML file and counts the MS1 (survey) scans.
///
/// Phase 0 reader spike (P0.2): proves `mzdata` can open a real, zlib-compressed mzML
/// from the FlashLFQ golden test data and that we can interrogate per-scan MS level.
/// Returns the number of spectra whose MS level is 1.
pub fn count_ms1_scans<P: Into<PathBuf> + Clone>(path: P) -> std::io::Result<usize> {
    let reader = MzMLReader::open_path(path)?;
    let ms1_count = reader.filter(|spectrum| spectrum.ms_level() == 1).count();
    Ok(ms1_count)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Resolves a path under the mzLib FlashLFQ golden test data directory
    /// (see [`super::mzlib_test_data`]).
    fn test_data(relative: &str) -> PathBuf {
        super::mzlib_test_data(relative)
    }

    #[test]
    fn version_is_nonempty() {
        assert!(!version().is_empty());
    }

    #[test]
    fn iota_f64_has_expected_length_and_values() {
        let v = iota_f64(5);
        assert_eq!(v.len(), 5);
        assert_eq!(v, vec![0.0, 1.0, 2.0, 3.0, 4.0]);
        assert!(iota_f64(0).is_empty());
    }

    #[test]
    fn demo_record_batch_has_expected_shape_and_values() {
        let batch = demo_record_batch();
        assert_eq!(batch.num_columns(), 2);
        assert_eq!(batch.num_rows(), 3);
        assert_eq!(batch.schema().field(0).name(), "mz");
        assert_eq!(batch.schema().field(1).name(), "intensity");

        let mz = batch
            .column(0)
            .as_any()
            .downcast_ref::<Float64Array>()
            .expect("column 0 is Float64");
        assert_eq!(mz.values(), &[100.0, 200.5, 300.25]);

        let intensity = batch
            .column(1)
            .as_any()
            .downcast_ref::<Float64Array>()
            .expect("column 1 is Float64");
        assert_eq!(intensity.values(), &[1.0e6, 2.5e6, 7.0e5]);
    }

    #[test]
    fn counts_ms1_scans_in_sliced_mzml() {
        let path = test_data("sliced-mzml.mzML");
        let ms1 = count_ms1_scans(path).expect("sliced-mzml.mzML should be readable");
        // Golden parity: the file holds 100 spectra (10 MS1 + 90 MS2), per the mzML
        // `ms level` cvParams. The C# mzLib reader sees the same survey scans.
        assert_eq!(ms1, 10, "expected 10 MS1 scans in sliced-mzml.mzML");
    }
}
