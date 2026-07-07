#!/usr/bin/env python
"""Score recall + charge-match of a detector resolved-feature TSV against an adapted
peak-level ground-truth table (TASK 3).

Mirrors the matching logic of the example's built-in `compare_to_reference`:
  - mass match: |feat_mono - ref_mono| / ref_mono * 1e6 <= MASS_PPM  (default 20 ppm)
  - RT match:   |feat_apex_rt - ref_rt| <= RT_DELTA  (min)
  - charge match (of matched refs): ref charge in the feature's Charge States list

Resolved-feature TSV columns (from detect_features_tsv write_tsv):
    ... RT Apex ... Charge States (semicolon) ... Monoisotopic Mass ...
Ground-truth TSV columns (build_ground_truth.py):
    mono_mass  mz  charge  rt  intensity  detection

Usage:
  score_recall.py <features.tsv> <ground_truth.tsv> [rt_delta_min ...]
Prints one block per rt_delta.
"""
import csv, sys, bisect

MASS_PPM = 20.0


def load_features(path):
    feats = []  # (mono_mass, apex_rt, set(charges))
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            try:
                mass = float(row["Monoisotopic Mass"])
                rt = float(row["RT Apex"])
            except (ValueError, KeyError):
                continue
            charges = set()
            for c in row.get("Charge States", "").split(";"):
                c = c.strip()
                if c:
                    try:
                        charges.add(int(c))
                    except ValueError:
                        pass
            feats.append((mass, rt, charges))
    feats.sort(key=lambda x: x[0])
    return feats


def load_refs(path):
    refs = []  # (mono_mass, rt, charge)
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            try:
                mass = float(row["mono_mass"]); rt = float(row["rt"]); z = int(row["charge"])
            except (ValueError, KeyError):
                continue
            refs.append((mass, rt, z))
    return refs


def score(feats, refs, rt_delta):
    masses = [f[0] for f in feats]
    matched = 0
    matched_charge = 0
    for (rmass, rrt, rz) in refs:
        tol = rmass * MASS_PPM / 1e6
        lo = bisect.bisect_left(masses, rmass - tol)
        hi = bisect.bisect_right(masses, rmass + tol)
        hit = None
        for i in range(lo, hi):
            fmass, frt, fchg = feats[i]
            if abs(frt - rrt) <= rt_delta:
                hit = feats[i]
                break
        if hit is not None:
            matched += 1
            if rz in hit[2]:
                matched_charge += 1
    return matched, matched_charge


def main():
    fpath, rpath = sys.argv[1], sys.argv[2]
    deltas = [float(x) for x in sys.argv[3:]] or [0.5]
    feats = load_features(fpath)
    refs = load_refs(rpath)
    print(f"features: {len(feats)}    ref peaks: {len(refs)}")
    for d in deltas:
        m, mc = score(feats, refs, d)
        n = max(len(refs), 1)
        print(f"  RT +-{d:.2f} min:  recall {m}/{len(refs)} = {100*m/n:.1f}%   "
              f"charge-match {mc}/{len(refs)} = {100*mc/n:.1f}%  "
              f"({100*mc/max(m,1):.1f}% of matched)")


if __name__ == "__main__":
    main()
