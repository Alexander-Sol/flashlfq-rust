"""99%-coverage merge: reference Peak MZ vs feature Most-Abundant m/z, carrying RT start/apex/end
from BOTH sides. For misses, also report the nearest feature (within NEAR_PPM, any RT) so the close
non-matching peak is visible.

Usage: python merge_99_rt.py
"""
import csv, bisect

REF = r"D:\SP_Tutorial\Lumos\CA_HCD_GPTMD_Search_WideTol\Task2-SearchTask\AllQuantifiedPeaks.tsv"
FEAT = r"D:\SP_Tutorial\Lumos\untargeted_features_10min_avg_cov99.tsv"
OUT = r"D:\SP_Tutorial\Lumos\untargeted_features_10min_avg_cov99.merged_mz_rt.tsv"

MZ_PPM = 10.0    # match: within this ppm on m/z ...
RT_MIN = 0.3     # ... AND within this many minutes on apex RT
NEAR_PPM = 50.0  # for misses: report the nearest feature within this ppm (any RT)


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


# --- load reference -------------------------------------------------------------------------
refs = []
with open(REF, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        mz = fnum(row.get("Peak MZ"))
        apex = fnum(row.get("Peak RT Apex"))
        if mz is None or apex is None:
            continue
        refs.append({
            "mz": mz,
            "start": fnum(row.get("Peak RT Start")),
            "apex": apex,
            "end": fnum(row.get("Peak RT End")),
            "charge": int(fnum(row.get("Peak Charge")) or 0),
            "intensity": fnum(row.get("Peak intensity")),
            "seq": (row.get("Full Sequence") or "").strip(),
        })

# --- load features --------------------------------------------------------------------------
feats = []
with open(FEAT, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        mz = fnum(row.get("Most-Abundant m/z"))
        apex = fnum(row.get("RT Apex"))
        if mz is None or apex is None:
            continue
        feats.append({
            "mz": mz,
            "start": fnum(row.get("RT Start")),
            "apex": apex,
            "end": fnum(row.get("RT End")),
            "charge": int(fnum(row.get("Primary Charge")) or 0),
            "mono_mz": fnum(row.get("Mono m/z (primary)")),
            "intensity": fnum(row.get("Summed Intensity")),
        })
feats.sort(key=lambda x: x["mz"])
fmz = [x["mz"] for x in feats]


def window(rmz, ppm):
    lo = bisect.bisect_left(fmz, rmz * (1 - ppm * 1e-6))
    hi = bisect.bisect_right(fmz, rmz * (1 + ppm * 1e-6))
    return lo, hi


def best_match(rp):
    lo, hi = window(rp["mz"], MZ_PPM)
    best, bppm = None, None
    for fe in feats[lo:hi]:
        if abs(fe["apex"] - rp["apex"]) <= RT_MIN:
            ppm = (fe["mz"] - rp["mz"]) / rp["mz"] * 1e6
            if best is None or abs(ppm) < abs(bppm):
                best, bppm = fe, ppm
    return best, bppm


def nearest_any(rp):
    lo, hi = window(rp["mz"], NEAR_PPM)
    best, bppm = None, None
    for fe in feats[lo:hi]:
        ppm = (fe["mz"] - rp["mz"]) / rp["mz"] * 1e6
        if best is None or abs(ppm) < abs(bppm):
            best, bppm = fe, ppm
    return best, bppm


# --- merge ----------------------------------------------------------------------------------
rows = []
n_match = n_charge = n_miss = 0
used = set()
for rp in refs:
    fe, ppm = best_match(rp)
    if fe is not None:
        n_match += 1
        used.add(id(fe))
        charge_ok = rp["charge"] and fe["charge"] == rp["charge"]
        n_charge += 1 if charge_ok else 0
        rows.append({"rp": rp, "fe": fe, "ppm": ppm,
                     "status": "OK" if charge_ok else "CHARGE-DIFF", "near": False})
    else:
        n_miss += 1
        near, nppm = nearest_any(rp)
        rows.append({"rp": rp, "fe": near, "ppm": nppm,
                     "status": "MISSED", "near": True})

extra = sum(1 for fe in feats if id(fe) not in used)

# --- write TSV ------------------------------------------------------------------------------
with open(OUT, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow([
        "Full Sequence", "Status",
        "Ref z", "Ref Peak MZ", "Ref RT Start", "Ref RT Apex", "Ref RT End",
        "Feat z", "Feat MostAbund m/z", "Feat Mono m/z", "Feat RT Start", "Feat RT Apex", "Feat RT End",
        "ppm(feat-ref)", "dRT Apex(min)", "Note",
    ])

    def g(v, p):
        return f"{v:.{p}f}" if isinstance(v, float) else ""

    for r in rows:
        rp, fe = r["rp"], r["fe"]
        if fe is not None:
            drt = fe["apex"] - rp["apex"]
            note = "nearest (no apex-RT match within %.1f min)" % RT_MIN if r["near"] else ""
        else:
            drt = None
            note = "no feature within %.0f ppm" % NEAR_PPM
        w.writerow([
            rp["seq"], r["status"],
            rp["charge"], g(rp["mz"], 5), g(rp["start"], 4), g(rp["apex"], 4), g(rp["end"], 4),
            (fe["charge"] if fe else ""), g(fe["mz"], 5) if fe else "", g(fe["mono_mz"], 5) if fe else "",
            g(fe["start"], 4) if fe else "", g(fe["apex"], 4) if fe else "", g(fe["end"], 4) if fe else "",
            g(r["ppm"], 2) if r["ppm"] is not None else "", g(drt, 3) if drt is not None else "", note,
        ])

# --- console: summary + samples -------------------------------------------------------------
n = len(refs)
print(f"99% coverage merge  (Ref Peak MZ vs Feat Most-Abundant m/z; match <= {MZ_PPM:.0f} ppm & {RT_MIN} min RT)")
print(f"  reference peaks : {n}")
print(f"  features        : {len(feats)}")
print(f"  matched         : {n_match} ({100*n_match/n:.1f}%)   charge also OK: {n_charge}")
print(f"  missed          : {n_miss} ({100*n_miss/n:.1f}%)")
print(f"  extra features  : {extra}")
print(f"  merged TSV -> {OUT}\n")

hdr = f'{"Status":>11} {"z":>2} {"RefPeakMZ":>10} | ref RT s/a/e            | feat RT s/a/e          {"ppm":>7} {"dRTapx":>7}  Seq'
print("--- sample: MATCHED (top 12 by ref intensity) ---")
print(hdr)
matched = sorted([r for r in rows if r["status"] != "MISSED"], key=lambda r: -(r["rp"]["intensity"] or 0))


def rtstr(d):
    def f(v):
        return f"{v:6.3f}" if isinstance(v, float) else "   -  "
    return f'{f(d["start"])} {f(d["apex"])} {f(d["end"])}'


for r in matched[:12]:
    rp, fe = r["rp"], r["fe"]
    print(f'{r["status"]:>11} {rp["charge"]:>2} {rp["mz"]:10.4f} | {rtstr(rp)} | {rtstr(fe)} '
          f'{r["ppm"]:+7.2f} {fe["apex"]-rp["apex"]:+7.2f}  {rp["seq"][:34]}')

print("\n--- sample: MISSED with nearest feature (top 15 by ref intensity) ---")
print(hdr)
missed = sorted([r for r in rows if r["status"] == "MISSED"], key=lambda r: -(r["rp"]["intensity"] or 0))
for r in missed[:15]:
    rp, fe = r["rp"], r["fe"]
    if fe is not None:
        print(f'{"MISS(near)":>11} {rp["charge"]:>2} {rp["mz"]:10.4f} | {rtstr(rp)} | {rtstr(fe)} '
              f'{r["ppm"]:+7.2f} {fe["apex"]-rp["apex"]:+7.2f}  {rp["seq"][:34]}')
    else:
        print(f'{"MISS(none)":>11} {rp["charge"]:>2} {rp["mz"]:10.4f} | {rtstr(rp)} | {"   -      -      -   ":>22} '
              f'{"":>7} {"":>7}  {rp["seq"][:34]}')
