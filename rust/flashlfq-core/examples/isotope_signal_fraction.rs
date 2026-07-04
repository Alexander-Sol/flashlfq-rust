//! Quick diagnostic: what fraction of the total ion current is *isotopically structured*?
//!
//! For every MS1 scan we keep a peak's intensity only if that peak belongs to a run of at least
//! `MIN_ISOTOPES` peaks (default 2, to match the detector's minimum) spaced at `C13_MINUS_C12 / z`
//! Th for some charge state `z` in `[MIN_Z, MAX_Z]`. Summing the intensity of the kept peaks and
//! dividing by ΣTIC yields
//! the fraction of signal that even *has* isotope-envelope structure — an upper bound on what any
//! deconvolution-based feature detector could ever attribute to real peptide isotope patterns.
//!
//! This is deliberately per-scan and cheap (no RT tracing, no deconvolution). It is a ceiling, not
//! a detection: a peak counts as isotopic if it lines up under *any* allowed charge, regardless of
//! envelope-shape plausibility, so the true explainable fraction is somewhat below this number.
//!
//! Usage:
//!   cargo run --release --example isotope_signal_fraction -- <spectra_file>
//!
//! Env knobs:
//!   PPM=10            m/z match tolerance for isotope spacing (default 10 ppm)
//!   MIN_Z=2 MAX_Z=6   charge-state range to test (default 2..=6)
//!   MIN_ISOTOPES=2    minimum peaks in a run for the run to count (default 2)
//!   MIN_INTENSITY=0   ignore peaks below this intensity when forming/counting chains
//!
//! Peaks arrive vendor-centroided from `read_ms1_scans` (Thermo `.raw` is centroided on read), so
//! the classification runs on one peak per isotope, not profile samples.
//!
//! Per-scan export (for manual validation / plotting — see `analysis/plot_isotope_classification.py`):
//!   EXPORT_TSV=path   write a per-peak classification TSV (scan_index, one_based_scan, rt, mz,
//!                     intensity, charge; charge 0 = excluded) for the scans named by EXPORT_SCANS
//!   EXPORT_SCANS=0    which 0-based scan indices to export: comma list and/or dash ranges,
//!                     e.g. "0", "100-105", "100,250,1000-1003" (default "0")

use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufWriter, Write};

use flashlfq_core::isotopic_envelope::C13_MINUS_C12;
use flashlfq_core::peak_indexing::read_ms1_scans;

fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key).ok().and_then(|s| s.parse().ok()).unwrap_or(default)
}
fn env_i32(key: &str, default: i32) -> i32 {
    std::env::var(key).ok().and_then(|s| s.parse().ok()).unwrap_or(default)
}

/// Parses an `EXPORT_SCANS` spec (comma-separated indices and/or `a-b` inclusive ranges) into a
/// sorted, de-duplicated set of 0-based scan indices. Malformed tokens are skipped.
fn parse_scan_selection(spec: &str) -> BTreeSet<usize> {
    let mut out = BTreeSet::new();
    for tok in spec.split(',') {
        let tok = tok.trim();
        if tok.is_empty() {
            continue;
        }
        match tok.split_once('-') {
            Some((a, b)) => {
                if let (Ok(a), Ok(b)) = (a.trim().parse::<usize>(), b.trim().parse::<usize>()) {
                    for i in a..=b.max(a) {
                        out.insert(i);
                    }
                }
            }
            None => {
                if let Ok(i) = tok.parse::<usize>() {
                    out.insert(i);
                }
            }
        }
    }
    out
}

/// Index of the peak whose m/z is closest to `target`, provided it lies within `tol_th` Th.
/// `mz` is assumed ascending (as centroided spectra are), so a binary search + two-neighbour check
/// suffices.
fn nearest_within(mz: &[f64], target: f64, tol_th: f64) -> Option<usize> {
    if mz.is_empty() {
        return None;
    }
    let idx = mz.partition_point(|&m| m < target);
    let mut best: Option<usize> = None;
    let mut best_d = tol_th;
    for c in [idx.wrapping_sub(1), idx] {
        if c < mz.len() {
            let d = (mz[c] - target).abs();
            if d <= best_d {
                best_d = d;
                best = Some(c);
            }
        }
    }
    best
}

/// Assigns each peak to a single charge: the one whose isotope run through that peak is longest
/// (ties broken toward the lower charge, which is tested first). Returns a parallel `Vec<i32>` over
/// the peaks: the assigned charge, or 0 if the peak is not part of any run of at least
/// `min_isotopes` peaks. Crediting each peak once keeps the per-charge intensity sums
/// non-overlapping, so they add up to the union total. `mz` must be ascending.
fn assign_charge(
    mz: &[f64],
    intensity: &[f64],
    min_z: i32,
    max_z: i32,
    ppm: f64,
    min_isotopes: usize,
    min_intensity: f64,
) -> Vec<i32> {
    let n = mz.len();
    let mut assigned = vec![0i32; n];
    let mut best_len = vec![0usize; n];
    if n < min_isotopes {
        return assigned;
    }
    for z in min_z..=max_z {
        let delta = C13_MINUS_C12 / z as f64;
        // next[i] = index of the peak one isotope step up from peak i (usize::MAX = none).
        let mut next = vec![usize::MAX; n];
        let mut has_prev = vec![false; n];
        for i in 0..n {
            if intensity[i] < min_intensity {
                continue;
            }
            let target = mz[i] + delta;
            let tol = target * ppm / 1e6;
            if let Some(j) = nearest_within(mz, target, tol) {
                // `next` links strictly upward in m/z (delta > tol), so chains never cycle.
                if j != i && j != usize::MAX && intensity[j] >= min_intensity {
                    next[i] = j;
                    has_prev[j] = true;
                }
            }
        }
        // Walk each chain from its head (a peak nothing points to); credit peaks in runs that reach
        // `min_isotopes`, keeping the charge with the longest supporting run for each peak.
        for i in 0..n {
            if has_prev[i] {
                continue;
            }
            let mut chain = Vec::new();
            let mut cur = i;
            while cur != usize::MAX {
                chain.push(cur);
                cur = next[cur];
            }
            let len = chain.len();
            if len >= min_isotopes {
                for &c in &chain {
                    if len > best_len[c] {
                        best_len[c] = len;
                        assigned[c] = z;
                    }
                }
            }
        }
    }
    assigned
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: isotope_signal_fraction <spectra_file>");
        std::process::exit(2);
    }
    let spectra_path = &args[1];

    let ppm = env_f64("PPM", 10.0);
    let min_z = env_i32("MIN_Z", 2);
    let max_z = env_i32("MAX_Z", 6);
    let min_isotopes = env_i32("MIN_ISOTOPES", 2).max(2) as usize;
    let min_intensity = env_f64("MIN_INTENSITY", 0.0);

    // Optional per-scan classification export for manual validation / plotting.
    let export_path = std::env::var("EXPORT_TSV").ok();
    let export_scans =
        parse_scan_selection(&std::env::var("EXPORT_SCANS").unwrap_or_else(|_| "0".into()));
    let mut export_writer: Option<BufWriter<File>> = export_path.as_ref().map(|p| {
        let mut w = BufWriter::new(File::create(p).expect("could not create EXPORT_TSV"));
        writeln!(w, "scan_index\tone_based_scan\trt\tmz\tintensity\tcharge").unwrap();
        w
    });

    eprintln!("reading MS1 scans from {spectra_path} ...");
    let scans = read_ms1_scans(spectra_path).expect("failed to read spectra file");

    let mut total_tic = 0.0f64;
    let mut isotopic_tic = 0.0f64;
    let mut total_peaks = 0usize;
    let mut isotopic_peaks = 0usize;
    // Per-charge TIC / peak tallies, indexed by charge (unused low indices stay zero).
    let mut tic_by_z = vec![0.0f64; (max_z + 1).max(0) as usize];
    let mut peaks_by_z = vec![0usize; (max_z + 1).max(0) as usize];

    // Peaks arrive vendor-centroided from `read_ms1_scans` (one peak per isotope), so classify
    // directly on the scan arrays.
    for (scan_idx, scan) in scans.iter().enumerate() {
        let charge = assign_charge(
            &scan.mz,
            &scan.intensity,
            min_z,
            max_z,
            ppm,
            min_isotopes,
            min_intensity,
        );
        let do_export = export_writer.is_some() && export_scans.contains(&scan_idx);
        for (i, &inten) in scan.intensity.iter().enumerate() {
            total_tic += inten;
            total_peaks += 1;
            let z = charge[i];
            if z != 0 {
                isotopic_tic += inten;
                isotopic_peaks += 1;
                tic_by_z[z as usize] += inten;
                peaks_by_z[z as usize] += 1;
            }
            if do_export {
                let w = export_writer.as_mut().unwrap();
                writeln!(
                    w,
                    "{}\t{}\t{:.5}\t{:.6}\t{:.4e}\t{}",
                    scan_idx, scan.one_based_scan_number, scan.retention_time, scan.mz[i], inten, z
                )
                .unwrap();
            }
        }
    }

    if let Some(mut w) = export_writer.take() {
        w.flush().unwrap();
        eprintln!(
            "exported {} scan(s) to {} (scans: {:?})",
            export_scans.len(),
            export_path.as_deref().unwrap_or(""),
            export_scans
        );
    }

    let tic_frac = if total_tic > 0.0 { 100.0 * isotopic_tic / total_tic } else { 0.0 };
    let peak_frac = if total_peaks > 0 {
        100.0 * isotopic_peaks as f64 / total_peaks as f64
    } else {
        0.0
    };

    println!("spectra file:        {spectra_path}");
    println!("MS1 scans:           {}   ({total_peaks} centroided peaks)", scans.len());
    println!(
        "charge range tested: {min_z}..={max_z}   tol: {ppm} ppm   min isotopes: {min_isotopes}   min intensity: {min_intensity}"
    );
    println!("ΣTIC (all peaks):    {total_tic:.4e}");
    println!("ΣTIC (isotopic):     {isotopic_tic:.4e}");
    println!(
        "isotope-structured signal: {tic_frac:.1}% of ΣTIC   ({isotopic_peaks} / {total_peaks} peaks = {peak_frac:.1}%)"
    );
    println!("\nper-charge breakdown (each peak credited to its strongest charge):");
    println!("  {:>5}  {:>12}  {:>9}  {:>12}", "z", "ΣTIC", "% of ΣTIC", "peaks");
    for z in min_z..=max_z {
        let zi = z as usize;
        let frac = if total_tic > 0.0 { 100.0 * tic_by_z[zi] / total_tic } else { 0.0 };
        println!("  {:>5}  {:>12.4e}  {:>8.1}%  {:>12}", z, tic_by_z[zi], frac, peaks_by_z[zi]);
    }
}
