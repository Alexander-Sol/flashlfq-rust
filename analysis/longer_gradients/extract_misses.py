#!/usr/bin/env python
"""Extract the ground-truth peaks a resolved-feature TSV MISSED (recall failures),
with diagnostics for manual investigation.

Match logic mirrors score_recall.py / compare_to_reference:
  mass: |feat_mono - ref_mono|/ref_mono*1e6 <= 20 ppm ; RT: |apex - ref_rt| <= 0.3 min.
A ref is a MISS if no feature satisfies BOTH.

For each miss we report *why*:
  - category = no_mass_candidate : no feature within 20 ppm at all (pure detection miss)
  - category = rt_miss           : a 20 ppm mass match exists but all are RT-outside 0.3 min
For rt_miss we give the nearest-RT mass-matching feature (its RT delta + charge).
We also give, for every miss, the single closest feature by mass (ppm) regardless of RT,
and its RT delta + charge, so a manual reviewer can eyeball near-misses.

Intensity decile is computed over ALL refs (D1 = faintest 10%, D10 = brightest).

Usage: extract_misses.py <features.tsv> <ground_truth.tsv> <out_misses.tsv> [rt_delta=0.3] [ppm=20]
"""
import csv, sys, bisect

def load_features(path):
    feats = []  # (mono, rt, charges_set)
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                mass = float(row["Monoisotopic Mass"]); rt = float(row["RT Apex"])
            except (ValueError, KeyError):
                continue
            chs = set()
            for c in row.get("Charge States", "").split(";"):
                c = c.strip()
                if c:
                    try: chs.add(int(c))
                    except ValueError: pass
            feats.append((mass, rt, chs))
    feats.sort(key=lambda x: x[0])
    return feats

def load_refs(path):
    refs = []
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                mass = float(row["mono_mass"]); mz = float(row["mz"])
                z = int(row["charge"]); rt = float(row["rt"]); inten = float(row["intensity"])
            except (ValueError, KeyError):
                continue
            refs.append((mass, mz, z, rt, inten))
    return refs

def main():
    fpath, rpath, opath = sys.argv[1], sys.argv[2], sys.argv[3]
    rt_delta = float(sys.argv[4]) if len(sys.argv) > 4 else 0.3
    ppm = float(sys.argv[5]) if len(sys.argv) > 5 else 20.0

    feats = load_features(fpath)
    refs  = load_refs(rpath)
    masses = [f[0] for f in feats]

    # intensity deciles over all refs
    order = sorted(range(len(refs)), key=lambda i: refs[i][4])
    decile = [0]*len(refs)
    n = len(refs)
    for rank, i in enumerate(order):
        decile[i] = min(10, 1 + rank*10//n)  # D1 faintest .. D10 brightest

    misses = []
    n_no_mass = n_rt = 0
    for i, (rmass, rmz, rz, rrt, rinten) in enumerate(refs):
        tol = rmass * ppm / 1e6
        lo = bisect.bisect_left(masses, rmass - tol)
        hi = bisect.bisect_right(masses, rmass + tol)
        # matched?
        matched = any(abs(feats[j][1] - rrt) <= rt_delta for j in range(lo, hi))
        if matched:
            continue
        # miss: diagnose
        mass_cands = list(range(lo, hi))
        if mass_cands:
            category = "rt_miss"; n_rt += 1
            # nearest-RT mass-matching feature
            j = min(mass_cands, key=lambda k: abs(feats[k][1] - rrt))
            near_rt_delta = feats[j][1] - rrt
            near_rt_charge = ";".join(str(c) for c in sorted(feats[j][2]))
        else:
            category = "no_mass_candidate"; n_no_mass += 1
            near_rt_delta = ""; near_rt_charge = ""
        # closest feature by mass overall (regardless of RT) for eyeballing
        cj = None; cbest = None
        for k in (lo-2, lo-1, lo, hi-1, hi, hi+1):
            if 0 <= k < len(feats):
                d = abs(feats[k][0] - rmass)
                if cbest is None or d < cbest:
                    cbest = d; cj = k
        if cj is not None:
            closest_ppm = (feats[cj][0] - rmass)/rmass*1e6
            closest_rt_delta = feats[cj][1] - rrt
            closest_charge = ";".join(str(c) for c in sorted(feats[cj][2]))
        else:
            closest_ppm = closest_rt_delta = ""; closest_charge = ""
        misses.append((rmass, rmz, rz, rrt, rinten, decile[i], category,
                       near_rt_delta, near_rt_charge, closest_ppm, closest_rt_delta, closest_charge))

    # sort: faintest first within category grouping? Just sort by intensity ascending (tail first).
    misses.sort(key=lambda r: r[4])
    with open(opath, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["mono_mass","mz","charge","rt","intensity","intensity_decile","category",
                    "nearest_rtmatch_rt_delta_min","nearest_rtmatch_charge",
                    "closest_by_mass_ppm","closest_by_mass_rt_delta_min","closest_by_mass_charge"])
        for m in misses:
            row = list(m)
            # format floats
            row[0] = f"{row[0]:.5f}"; row[1] = f"{row[1]:.5f}"; row[3] = f"{row[3]:.4f}"
            row[4] = f"{row[4]:.4e}"
            if row[7] != "": row[7] = f"{row[7]:.4f}"
            if row[9] != "": row[9] = f"{row[9]:.2f}"
            if row[10] != "": row[10] = f"{row[10]:.4f}"
            w.writerow(row)

    tot = len(misses)
    print(f"refs {len(refs)}  features {len(feats)}")
    print(f"misses {tot} ({100*tot/len(refs):.1f}%)  ->  {opath}")
    print(f"  no_mass_candidate: {n_no_mass} ({100*n_no_mass/max(tot,1):.1f}% of misses)")
    print(f"  rt_miss (mass ok, RT off): {n_rt} ({100*n_rt/max(tot,1):.1f}% of misses)")
    # decile breakdown of misses
    dc = [0]*11
    for m in misses: dc[m[5]] += 1
    dtot = [0]*11
    for i in range(len(refs)): dtot[decile[i]] += 1
    print("  decile  misses/refs  (miss-rate)")
    for d in range(1, 11):
        mr = 100*dc[d]/max(dtot[d],1)
        print(f"    D{d:<2} {dc[d]:5d}/{dtot[d]:<5d}  {mr:5.1f}%")

if __name__ == "__main__":
    main()
