//! MS1 feature export in the TopFD / msDeconv **`.msalign`** (MS1) interchange format.
//!
//! `.msalign` is the one MS1 feature/deconvolution format that mzLib
//! (github.com/smith-chem-wisc/mzLib) can natively read — its `Readers/MsAlign/Ms1Align`
//! reader parses MS1 deconvoluted-spectrum blocks into `MsDataScan[]` — and it is a
//! widely-used top-down interchange format (TopFD writes it; TopPIC / TopMG / FLASHDeconv
//! consume it). See `agent_info/Feature-Output-Formats.md` for the decision + full spec.
//!
//! Layout (exactly what mzLib's `MsAlign.cs` parses):
//! ```text
//! ##### Parameters #####          # optional file-level header, delimited by this literal
//! #File name: <name>              # `#key : value` lines
//! ##### Parameters #####
//! BEGIN IONS                      # one block == one deconvoluted MS1 scan
//! ID=<n>                          # KEY=VALUE entry headers
//! SCANS=<one-based scan number>
//! RETENTION_TIME=<seconds>        # mzLib divides by 60 to store minutes
//! LEVEL=1
//! <monoMass>\t<intensity>\t<charge>   # tab-delimited peak lines, one per feature
//! END IONS
//! ```
//! Features are grouped by apex scan so each block is a genuine per-scan spectrum and no
//! two blocks share a `SCANS` number.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::feature_refinement::ResolvedFeature;

/// One deconvoluted MS1 feature as a single `.msalign` peak entry: a neutral monoisotopic
/// mass, its summed intensity, the charge, and the apex scan it is written under. A neutral
/// intermediate so the writer stays decoupled from the pipeline's feature types (and is
/// trivially constructible in tests).
#[derive(Debug, Clone, PartialEq)]
pub struct Ms1AlignFeature {
    /// Neutral monoisotopic mass (peak column 0).
    pub monoisotopic_mass: f64,
    /// Summed feature intensity (peak column 1).
    pub intensity: f64,
    /// Charge state (peak column 2).
    pub charge: i32,
    /// One-based scan number the feature apexes at — becomes the block's `SCANS`.
    pub apex_scan_number: i32,
    /// Apex retention time in **minutes** (converted to seconds on write for `RETENTION_TIME`).
    pub apex_retention_time_minutes: f64,
}

/// Writes `features` as an MS1 `.msalign` document to `writer`.
///
/// Features are grouped by [`Ms1AlignFeature::apex_scan_number`] into one `BEGIN IONS`
/// block each (ascending scan number); within a block, peaks are sorted by ascending
/// monoisotopic mass (TopFD convention). `source_file_name` is recorded in the optional
/// file-level header (`#File name:`); pass the spectra file's basename.
pub fn write_ms1_align<W: Write>(
    writer: &mut W,
    features: &[Ms1AlignFeature],
    source_file_name: &str,
) -> io::Result<()> {
    // Group by scan → one deconvoluted-spectrum block per apex scan (BTreeMap = ascending SCANS).
    let mut blocks: BTreeMap<i32, Vec<&Ms1AlignFeature>> = BTreeMap::new();
    for f in features {
        blocks.entry(f.apex_scan_number).or_default().push(f);
    }

    // File-level parameter header (optional per the spec, but faithful to TopFD output).
    writeln!(writer, "##### Parameters #####")?;
    writeln!(writer, "#File name: {source_file_name}")?;
    writeln!(writer, "#Number of MS1 scans: {}", blocks.len())?;
    writeln!(writer, "#Software: flashlfq-rust untargeted feature detector")?;
    writeln!(writer, "##### Parameters #####")?;

    for (id, (scan, feats)) in blocks.iter().enumerate() {
        // All features in a block share the apex scan, so any member's RT represents it.
        let rt_seconds = feats
            .first()
            .map(|f| f.apex_retention_time_minutes * 60.0)
            .unwrap_or(0.0);
        writeln!(writer, "BEGIN IONS")?;
        writeln!(writer, "ID={id}")?;
        writeln!(writer, "SCANS={scan}")?;
        writeln!(writer, "RETENTION_TIME={rt_seconds:.4}")?;
        writeln!(writer, "LEVEL=1")?;
        let mut sorted = feats.clone();
        sorted.sort_by(|a, b| a.monoisotopic_mass.total_cmp(&b.monoisotopic_mass));
        for f in sorted {
            // Tab-delimited: monoMass \t intensity \t charge (mzLib parses exactly these 3).
            writeln!(
                writer,
                "{:.5}\t{:.4}\t{}",
                f.monoisotopic_mass, f.intensity, f.charge
            )?;
        }
        writeln!(writer, "END IONS")?;
        writeln!(writer)?;
    }
    Ok(())
}

/// Convenience wrapper: writes an MS1 `.msalign` file at `path`. See [`write_ms1_align`].
pub fn write_ms1_align_file<P: AsRef<Path>>(
    path: P,
    features: &[Ms1AlignFeature],
    source_file_name: &str,
) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    write_ms1_align(&mut writer, features, source_file_name)?;
    writer.flush()
}

/// Maps the pipeline's resolved features to [`Ms1AlignFeature`] entries.
///
/// Each resolved feature becomes one entry at its neutral monoisotopic mass and summed
/// intensity. The charge and apex scan are taken from the feature's **primary** (tallest)
/// member — the same primary-member convention the resolved TSV uses. Features with no
/// members are skipped (nothing to place).
pub fn resolved_to_ms1align_features(resolved: &[ResolvedFeature]) -> Vec<Ms1AlignFeature> {
    resolved
        .iter()
        .filter_map(|r| {
            let primary = r.members.iter().max_by(|a, b| {
                a.detected
                    .summed_intensity
                    .total_cmp(&b.detected.summed_intensity)
            })?;
            Some(Ms1AlignFeature {
                monoisotopic_mass: r.monoisotopic_mass,
                intensity: r.summed_intensity,
                charge: primary.refined_charge,
                // Detector apex scan index is zero-based; `.msalign` SCANS is one-based.
                apex_scan_number: primary.detected.apex_scan_index.max(0) + 1,
                apex_retention_time_minutes: r.apex_rt,
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_to_string(features: &[Ms1AlignFeature], name: &str) -> String {
        let mut buf: Vec<u8> = Vec::new();
        write_ms1_align(&mut buf, features, name).expect("write to Vec never fails");
        String::from_utf8(buf).expect("output is valid UTF-8")
    }

    #[test]
    fn writes_spec_correct_ms1_msalign() {
        // Two features apexing at scan 5 (grouped into one block) + one at scan 9.
        let features = vec![
            Ms1AlignFeature {
                monoisotopic_mass: 1500.75000,
                intensity: 1.2e7,
                charge: 2,
                apex_scan_number: 5,
                apex_retention_time_minutes: 10.0, // -> 600 s
            },
            Ms1AlignFeature {
                monoisotopic_mass: 900.12345,
                intensity: 3.4e6,
                charge: 1,
                apex_scan_number: 5,
                apex_retention_time_minutes: 10.0,
            },
            Ms1AlignFeature {
                monoisotopic_mass: 2500.50000,
                intensity: 5.6e5,
                charge: 3,
                apex_scan_number: 9,
                apex_retention_time_minutes: 12.5, // -> 750 s
            },
        ];
        let out = write_to_string(&features, "sample.raw");
        let lines: Vec<&str> = out.lines().collect();

        // --- file-level header: the delimiter literal appears exactly twice. ---
        assert_eq!(
            out.matches("##### Parameters #####").count(),
            2,
            "header must open and close with the exact delimiter"
        );
        assert!(out.contains("#File name: sample.raw"));
        assert!(out.contains("#Number of MS1 scans: 2"));

        // --- block structure: one BEGIN/END per distinct scan (2 scans -> 2 blocks). ---
        assert_eq!(out.matches("BEGIN IONS").count(), 2);
        assert_eq!(out.matches("END IONS").count(), 2);

        // --- first block is scan 5 (ascending), carries both scan-5 features. ---
        let b1 = lines.iter().position(|l| *l == "BEGIN IONS").unwrap();
        assert_eq!(lines[b1 + 1], "ID=0");
        assert_eq!(lines[b1 + 2], "SCANS=5");
        assert_eq!(lines[b1 + 3], "RETENTION_TIME=600.0000"); // 10 min -> 600 s
        assert_eq!(lines[b1 + 4], "LEVEL=1");
        // Peaks sorted by ascending mass: 900.12345 before 1500.75.
        assert_eq!(lines[b1 + 5], "900.12345\t3400000.0000\t1");
        assert_eq!(lines[b1 + 6], "1500.75000\t12000000.0000\t2");
        assert_eq!(lines[b1 + 7], "END IONS");

        // Every peak line has exactly 3 tab-separated columns (mzLib's contract).
        for l in &lines {
            if l.contains('\t') {
                assert_eq!(l.split('\t').count(), 3, "peak line must be mass/intensity/charge");
            }
        }

        // --- second block is scan 9 with the single high-mass feature. ---
        assert!(out.contains("SCANS=9"));
        assert!(out.contains("RETENTION_TIME=750.0000"));
        assert!(out.contains("2500.50000\t560000.0000\t3"));
    }

    #[test]
    fn empty_features_emit_header_only() {
        let out = write_to_string(&[], "empty.raw");
        assert_eq!(out.matches("##### Parameters #####").count(), 2);
        assert!(out.contains("#Number of MS1 scans: 0"));
        assert_eq!(out.matches("BEGIN IONS").count(), 0);
    }
}
