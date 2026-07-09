#!/usr/bin/env python
"""Build a per-feature table (detected + refined scores + matched flag) for the
interactive low-quality-feature notebook. One parquet per file base.

detected/refined TSVs are row-aligned (verified: refined 'Detector Mono Mass' ==
detected 'Monoisotopic Mass' per row), so we merge positionally.

matched = does the feature fall within +-20 ppm mass AND +-0.5 min RT of any GT
peak (charge-blind, same as the recall scorer). We also record the nearest GT
peak's ppm/RT deltas for context.

Usage: prep_notebook_data.py [base ...]   (default: 3 representative files)
"""
import csv, os, sys, bisect
import numpy as np
import pandas as pd

BENCH = r"F:\flashlfq-rust\analysis\synthpep_bench"
DATA = r"D:\SyntheticPeptides_PXD001091"
MASS_PPM = 20.0
RT_DELTA = 0.5
DEFAULT_BASES = ["130124_dilA_10_01", "130124_dilA_2_01", "130124_dilA_12_01"]


def load_detected(path):
    cols = ["det_mass", "charge", "mono_mz", "apex_rt", "rt_start", "rt_end",
            "intensity", "detector_score", "num_isotopes", "num_peaks"]
    d = {c: [] for c in cols}
    with open(path, encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            d["det_mass"].append(float(p[0])); d["charge"].append(int(p[1]))
            d["mono_mz"].append(float(p[2])); d["apex_rt"].append(float(p[3]))
            d["rt_start"].append(float(p[4])); d["rt_end"].append(float(p[5]))
            d["intensity"].append(float(p[6])); d["detector_score"].append(float(p[7]))
            d["num_isotopes"].append(int(p[8])); d["num_peaks"].append(int(p[9]))
    return pd.DataFrame(d)


def load_refined_cols(path):
    mass, decon, ncand = [], [], []
    with open(path, encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            mass.append(float(p[0])); decon.append(float(p[4])); ncand.append(int(p[5]))
    return np.asarray(mass), np.asarray(decon), np.asarray(ncand)


def load_gt(path):
    g = []
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g.append((float(row["mono_mass"]), float(row["rt"]), int(row["charge"])))
    g.sort()
    return g


def annotate_matches(df, gt):
    gmass = [x[0] for x in gt]
    grt = np.array([x[1] for x in gt])
    matched = np.zeros(len(df), bool)
    best_ppm = np.full(len(df), np.nan)
    best_drt = np.full(len(df), np.nan)
    masses = df["mass"].to_numpy()
    rts = df["apex_rt"].to_numpy()
    for i in range(len(df)):
        rmass, rrt = masses[i], rts[i]
        tol = rmass * MASS_PPM / 1e6
        lo = bisect.bisect_left(gmass, rmass - tol)
        hi = bisect.bisect_right(gmass, rmass + tol)
        if hi > lo:
            drt = np.abs(grt[lo:hi] - rrt)
            j = int(np.argmin(drt))
            best_ppm[i] = abs(rmass - gmass[lo + j]) / rmass * 1e6
            best_drt[i] = drt[j]
            if drt[j] <= RT_DELTA:
                matched[i] = True
    df["matched"] = matched
    df["nearest_gt_ppm"] = best_ppm
    df["nearest_gt_drt"] = best_drt
    return df


def main():
    bases = sys.argv[1:] or DEFAULT_BASES
    outdir = os.path.join(BENCH, "nb_data")
    os.makedirs(outdir, exist_ok=True)
    manifest = []
    for base in bases:
        od = os.path.join(BENCH, "ours", base)
        df = load_detected(os.path.join(od, "feat.detected.tsv"))
        rmass, decon, ncand = load_refined_cols(os.path.join(od, "feat.refined.tsv"))
        assert len(rmass) == len(df), f"row mismatch {base}"
        df["mass"] = rmass                 # refined (shifted) mass
        df["decon_score"] = decon
        df["num_candidate_masses"] = ncand
        df["log10_intensity"] = np.log10(np.clip(df["intensity"], 1.0, None))
        gt = load_gt(os.path.join(BENCH, "gt", base + ".tsv"))
        df = annotate_matches(df, gt)
        df["feat_idx"] = np.arange(len(df))
        out = os.path.join(outdir, base + ".parquet")
        df.to_parquet(out, index=False)
        n_un = int((~df["matched"]).sum())
        manifest.append((base, len(df), n_un, len(gt)))
        print(f"{base}: {len(df)} features, {n_un} unmatched, {len(gt)} gt -> {out}")
    with open(os.path.join(outdir, "manifest.tsv"), "w", encoding="utf-8") as fh:
        fh.write("base\tn_features\tn_unmatched\tn_gt\tmzml\n")
        for (base, nf, nu, ng) in manifest:
            fh.write(f"{base}\t{nf}\t{nu}\t{ng}\t{os.path.join(DATA, base + '_uncalibrated.mzML')}\n")
    print("done ->", outdir)


if __name__ == "__main__":
    main()
