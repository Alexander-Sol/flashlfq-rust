#!/usr/bin/env python
"""Batch recall scorer for the PXD001091 synthetic-peptide benchmark.

For every file base in gt/gt_manifest.tsv, score BOTH detectors against the
per-file ground truth with identical tolerances (mass +-20 ppm, RT window,
charge-among-feature-charges), exactly as analysis/dinosaur_bench/score_any.py.

Emits:
  results_perfile.tsv  -- one row per (base) with recall/charge for both tools
  prints a micro-averaged aggregate (sum matched / sum refs over all files)

Usage: score_all.py <bench_dir> [rt_delta ...]   (default deltas 0.5 1.0)
"""
import csv, sys, os, bisect

MASS_PPM = 20.0


def load_features(path):
    feats = []
    if not os.path.exists(path):
        return feats
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        cols = set(r.fieldnames or [])
        is_dino = "rtApex" in cols and "charge" in cols
        for row in r:
            try:
                if is_dino:
                    mass = float(row["mass"]); rt = float(row["rtApex"])
                    charges = {int(float(row["charge"]))}
                else:
                    mass = float(row["Monoisotopic Mass"]); rt = float(row["RT Apex"])
                    charges = {int(c) for c in row.get("Charge States", "").split(";") if c.strip()}
            except (ValueError, KeyError):
                continue
            feats.append((mass, rt, charges))
    feats.sort(key=lambda x: x[0])
    return feats


def load_refs(path):
    refs = []
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            try:
                refs.append((float(row["mono_mass"]), float(row["rt"]), int(row["charge"])))
            except (ValueError, KeyError):
                continue
    return refs


def score(feats, refs, rt_delta):
    masses = [f[0] for f in feats]
    matched = matched_charge = 0
    for (rmass, rrt, rz) in refs:
        tol = rmass * MASS_PPM / 1e6
        lo = bisect.bisect_left(masses, rmass - tol)
        hi = bisect.bisect_right(masses, rmass + tol)
        hit = None
        for i in range(lo, hi):
            if abs(feats[i][1] - rrt) <= rt_delta:
                hit = feats[i]; break
        if hit is not None:
            matched += 1
            if rz in hit[2]:
                matched_charge += 1
    return matched, matched_charge


def main():
    bench = sys.argv[1]
    deltas = [float(x) for x in sys.argv[2:]] or [0.5, 1.0]
    gtdir = os.path.join(bench, "gt")
    bases = []
    with open(os.path.join(gtdir, "gt_manifest.tsv"), newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            bases.append(row["base"])

    # aggregate accumulators: tool -> delta -> [matched, matched_charge]; refs total; features total
    agg = {t: {d: [0, 0] for d in deltas} for t in ("ours", "dino")}
    tot_refs = 0
    tot_feat = {"ours": 0, "dino": 0}

    out = os.path.join(bench, "results_perfile.tsv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        hdr = ["base", "ref_peaks", "ours_feat", "dino_feat"]
        for d in deltas:
            hdr += [f"ours_recall_{d}", f"ours_charge_{d}", f"dino_recall_{d}", f"dino_charge_{d}"]
        w.writerow(hdr)
        for base in bases:
            refs = load_refs(os.path.join(gtdir, base + ".tsv"))
            ours = load_features(os.path.join(bench, "ours", base, "feat.tsv"))
            dino = load_features(os.path.join(bench, "dino", base + ".features.tsv"))
            tot_refs += len(refs)
            tot_feat["ours"] += len(ours); tot_feat["dino"] += len(dino)
            nref = max(len(refs), 1)
            row = [base, len(refs), len(ours), len(dino)]
            for d in deltas:
                om, omc = score(ours, refs, d)
                dm, dmc = score(dino, refs, d)
                agg["ours"][d][0] += om; agg["ours"][d][1] += omc
                agg["dino"][d][0] += dm; agg["dino"][d][1] += dmc
                row += [f"{100*om/nref:.1f}", f"{100*omc/nref:.1f}",
                        f"{100*dm/nref:.1f}", f"{100*dmc/nref:.1f}"]
            w.writerow(row)

    print(f"files: {len(bases)}   total ref peaks: {tot_refs}")
    print(f"total features   ours: {tot_feat['ours']:,}   dino: {tot_feat['dino']:,}")
    nref = max(tot_refs, 1)
    print("\nMICRO-AVERAGED RECALL (sum matched / sum refs across all files):")
    for d in deltas:
        om, omc = agg["ours"][d]
        dm, dmc = agg["dino"][d]
        print(f"  RT +-{d:.2f} min:")
        print(f"    OURS     recall {om}/{tot_refs} = {100*om/nref:.2f}%   "
              f"charge {100*omc/nref:.2f}%  ({100*omc/max(om,1):.1f}% of matched)")
        print(f"    DINOSAUR recall {dm}/{tot_refs} = {100*dm/nref:.2f}%   "
              f"charge {100*dmc/nref:.2f}%  ({100*dmc/max(dm,1):.1f}% of matched)")
    print(f"\nper-file table -> {out}")


if __name__ == "__main__":
    main()
