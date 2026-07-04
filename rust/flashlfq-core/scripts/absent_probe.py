"""For the genuinely-absent reference peaks: are they truly undetected, or detected at a different
charge/mass? Check (a) does any feature's Most-Abundant m/z match the ref's observed Peak MZ (+-10 ppm)
at +-0.3 min RT — i.e. the ISOTOPE is there but neutral mass/charge differs; (b) the nearest feature by
RT at that observed m/z, with its assigned charge."""
import csv, bisect
from collections import Counter

FEAT = r"D:\SP_Tutorial\Lumos\untargeted_features_10min_centroid_dd_cov99_minscan2.tsv"
REF = r"D:\SP_Tutorial\Lumos\CA_HCD_GPTMD_Search_WideTol\Task2-SearchTask\AllQuantifiedPeaks.tsv"
C13 = 1.00335483810
MASS_PPM, MZ_PPM, RT_MIN = 20.0, 10.0, 0.3


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


feats = []
with open(FEAT, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        feats.append({
            "m": fnum(row.get("Monoisotopic Mass")),
            "abmz": fnum(row.get("Most-Abundant m/z")),
            "rt": fnum(row.get("RT Apex")),
            "z": int(fnum(row.get("Primary Charge")) or 0),
        })
by_mass = sorted([f for f in feats if f["m"] is not None], key=lambda x: x["m"])
mkeys = [f["m"] for f in by_mass]
by_mz = sorted([f for f in feats if f["abmz"] is not None], key=lambda x: x["abmz"])
zkeys = [f["abmz"] for f in by_mz]


def match_mass(mass, rt):
    lo = bisect.bisect_left(mkeys, mass * (1 - MASS_PPM * 1e-6))
    hi = bisect.bisect_right(mkeys, mass * (1 + MASS_PPM * 1e-6))
    return any(abs(fe["rt"] - rt) <= RT_MIN for fe in by_mass[lo:hi])


def match_mz(mz, rt):
    lo = bisect.bisect_left(zkeys, mz * (1 - MZ_PPM * 1e-6))
    hi = bisect.bisect_right(zkeys, mz * (1 + MZ_PPM * 1e-6))
    hits = [fe for fe in by_mz[lo:hi] if abs(fe["rt"] - rt) <= RT_MIN]
    return hits


def nearest_rt(rt, tol=0.3):
    return [fe for fe in feats if fe["rt"] is not None and abs(fe["rt"] - rt) <= tol]


refs = []
with open(REF, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        m = fnum(row.get("Peptide Monoisotopic Mass"))
        rt = fnum(row.get("Peak RT Apex"))
        pkmz = fnum(row.get("Peak MZ"))
        z = int(fnum(row.get("Peak Charge")) or 0)
        seq = (row.get("Full Sequence") or "").strip()
        if m is None or rt is None:
            continue
        refs.append({"m": m, "rt": rt, "pkmz": pkmz, "z": z, "seq": seq})

# absent = no mass match at ref +- k*C13
absent = []
for r in refs:
    if match_mass(r["m"], r["rt"]):
        continue
    if any(match_mass(r["m"] + k * C13, r["rt"]) for k in (-1, 1, -2, 2)):
        continue
    absent.append(r)

print(f"absent reference peaks: {len(absent)}\n")
cat = Counter()
print(f"{'mass':>9} {'RT':>6} {'z':>2} {'PeakMZ':>10}  m/z-match?  nearest-feat-at-RT (z:count)   seq")
for r in absent:
    mzhits = match_mz(r["pkmz"], r["rt"]) if r["pkmz"] else []
    near = nearest_rt(r["rt"])
    zc = Counter(fe["z"] for fe in near)
    if mzhits:
        cat["isotope present, mass/charge differs"] += 1
        tag = f"YES z{Counter(fe['z'] for fe in mzhits).most_common(1)[0][0]}"
    elif near:
        cat["no m/z match, but features co-elute"] += 1
        tag = "no"
    else:
        cat["nothing co-elutes at that RT"] += 1
        tag = "no"
    zc_s = ",".join(f"z{z}:{c}" for z, c in sorted(zc.items()))
    print(f"{r['m']:9.3f} {r['rt']:6.2f} {r['z']:2} {(r['pkmz'] or 0):10.4f}  {tag:>8}   {zc_s:28} {r['seq'][:34]}")

print("\nsummary:")
for k, v in cat.most_common():
    print(f"  {v:2}  {k}")
