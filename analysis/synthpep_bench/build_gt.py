#!/usr/bin/env python
"""Per-file peak-level ground truth from the FragPipe MSFragger psm.tsv for the
PXD001091 synthetic-peptide dilution series.

psm.tsv is one row per PSM (per MS2 scan). We dedupe to unique
(Spectrum-File base, Modified Peptide, Charge) peaks -- each such key is one
expected MS1 feature. rt = MEDIAN of the PSM retention times for that key
(Retention is in SECONDS in psm.tsv; we convert to MINUTES to match the
detectors' RT-apex units). mono_mass = MSFragger "Observed Mass" (neutral,
isotope-corrected monoisotopic, uncalibrated -- consistent with running the
detectors on uncalibrated data).

Output schema mirrors score_any.py's loader: mono_mass  mz  charge  rt  n_psm
One TSV per file base in gt/<base>.tsv, plus a manifest gt_manifest.tsv.

Usage: build_gt.py <psm.tsv> <out_gt_dir>
"""
import csv, sys, os, re, statistics

PROTON = 1.007276466812
BASE_RE = re.compile(r"interact-(.+?)\.pep\.xml$", re.IGNORECASE)


def base_of(spectrum_file):
    m = BASE_RE.search(spectrum_file)
    if m:
        return m.group(1)
    # fallback: strip dir + extension
    b = os.path.basename(spectrum_file)
    for ext in (".pepXML", ".pep.xml", ".mzML", ".raw"):
        if b.lower().endswith(ext.lower()):
            b = b[: -len(ext)]
    return b


def main():
    psm, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    # key -> list of (mono_mass, mz, charge, rt_min)
    per_file = {}   # base -> { (pep,z): [masses, mzs, rts] }
    n_rows = decoy = 0
    with open(psm, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            n_rows += 1
            if row.get("Is Decoy", "").strip().lower() == "true":
                decoy += 1
                continue
            try:
                z = int(row["Charge"])
                rt_min = float(row["Retention"]) / 60.0
                mass = float(row["Observed Mass"])
                mz = float(row["Observed M/Z"])
            except (ValueError, KeyError):
                continue
            base = base_of(row["Spectrum File"])
            pep = row.get("Modified Peptide") or row.get("Peptide") or ""
            d = per_file.setdefault(base, {})
            rec = d.setdefault((pep, z), [[], [], []])
            rec[0].append(mass); rec[1].append(mz); rec[2].append(rt_min)

    manifest = []
    for base, d in sorted(per_file.items()):
        path = os.path.join(outdir, base + ".tsv")
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["mono_mass", "mz", "charge", "rt", "n_psm"])
            rows = []
            for (pep, z), (masses, mzs, rts) in d.items():
                rows.append((statistics.median(masses), statistics.median(mzs),
                             z, statistics.median(rts), len(rts)))
            rows.sort()
            for (mass, mz, z, rt, n) in rows:
                w.writerow([f"{mass:.5f}", f"{mz:.5f}", z, f"{rt:.4f}", n])
        manifest.append((base, len(d), sum(len(v[2]) for v in d.values())))

    with open(os.path.join(outdir, "gt_manifest.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["base", "unique_peaks", "n_psm"])
        for row in manifest:
            w.writerow(row)

    tot_peaks = sum(m[1] for m in manifest)
    tot_psm = sum(m[2] for m in manifest)
    print(f"{n_rows} PSM rows ({decoy} decoy dropped) -> {len(manifest)} files, "
          f"{tot_peaks} unique peaks, {tot_psm} PSMs kept -> {outdir}")


if __name__ == "__main__":
    main()
