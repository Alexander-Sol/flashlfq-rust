"""Characterize the reference-peak misses of the default (RawSum) pipeline: classify each of the 620
reference peaks as matched / recoverable-off-by-one / genuinely-absent, then break the misses down by
intensity, charge, mass, RT, off-by-one direction, charge multiplicity, and sequence modifications.
Goal: reveal WHAT the algorithm misses, to ground a discriminator plan."""
import csv, bisect
from collections import Counter

FEAT = r"D:\SP_Tutorial\Lumos\untargeted_features_10min_centroid_dd_cov99_minscan2.tsv"  # RawSum baseline
REF = r"D:\SP_Tutorial\Lumos\CA_HCD_GPTMD_Search_WideTol\Task2-SearchTask\AllQuantifiedPeaks.tsv"
C13 = 1.00335483810
MASS_PPM, RT_MIN = 20.0, 0.3


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


# --- features (RawSum resolved) ---
feats = []
with open(FEAT, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        m, rt = fnum(row.get("Monoisotopic Mass")), fnum(row.get("RT Apex"))
        ncs = fnum(row.get("Num Charge States"))
        if m is not None and rt is not None:
            feats.append({"m": m, "rt": rt, "ncs": int(ncs or 0)})
feats.sort(key=lambda x: x["m"])
keys = [x["m"] for x in feats]


def find(mass, rt):
    lo = bisect.bisect_left(keys, mass * (1 - MASS_PPM * 1e-6))
    hi = bisect.bisect_right(keys, mass * (1 + MASS_PPM * 1e-6))
    for fe in feats[lo:hi]:
        if abs(fe["rt"] - rt) <= RT_MIN:
            return fe
    return None


# --- reference ---
refs = []
with open(REF, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        m = fnum(row.get("Peptide Monoisotopic Mass"))
        rt = fnum(row.get("Peak RT Apex"))
        z = fnum(row.get("Peak Charge"))
        inten = fnum(row.get("Peak intensity"))
        seq = (row.get("Full Sequence") or "").strip()
        if m is None or rt is None:
            continue
        refs.append({"m": m, "rt": rt, "z": int(z or 0), "int": inten, "seq": seq,
                     "mods": seq.count("["), "len": sum(1 for c in seq if c.isalpha() and c.isupper())})

# --- classify ---
for r in refs:
    if find(r["m"], r["rt"]):
        r["cls"] = "matched"
        r["k"] = 0
        r["carrier"] = None
    else:
        hit = None
        for k in (-1, 1, -2, 2):
            fe = find(r["m"] + k * C13, r["rt"])
            if fe:
                hit = (k, fe)
                break
        if hit:
            r["cls"] = "offbyone"
            r["k"] = hit[0]
            r["carrier"] = hit[1]
        else:
            r["cls"] = "absent"
            r["k"] = None
            r["carrier"] = None

n = len(refs)
by = Counter(r["cls"] for r in refs)
print(f"reference peaks: {n}")
print(f"  matched    {by['matched']:3} ({100*by['matched']/n:.1f}%)")
print(f"  offbyone   {by['offbyone']:3} ({100*by['offbyone']/n:.1f}%)")
print(f"  absent     {by['absent']:3} ({100*by['absent']/n:.1f}%)")

# --- miss rate by reference intensity decile ---
with_int = [r for r in refs if r["int"] is not None]
with_int.sort(key=lambda r: r["int"])
print(f"\nmiss rate by reference-intensity decile ({len(with_int)} refs with intensity):")
d = len(with_int) // 10 or 1
for i in range(10):
    chunk = with_int[i * d: (i + 1) * d] if i < 9 else with_int[9 * d:]
    if not chunk:
        continue
    miss = sum(1 for r in chunk if r["cls"] != "matched")
    lo, hi = chunk[0]["int"], chunk[-1]["int"]
    print(f"  D{i+1:2} int[{lo:9.2e}..{hi:9.2e}] n={len(chunk):3}  miss {miss:2} ({100*miss/len(chunk):4.0f}%)")

# --- miss breakdown by charge and by mass ---
def breakdown(key, label, bins):
    print(f"\nmiss rate by {label}:")
    for lo, hi in bins:
        sub = [r for r in refs if lo <= r[key] < hi]
        if not sub:
            continue
        miss = sum(1 for r in sub if r["cls"] != "matched")
        obo = sum(1 for r in sub if r["cls"] == "offbyone")
        print(f"  {label} [{lo},{hi}): n={len(sub):3}  miss {miss:2} ({100*miss/len(sub):4.0f}%)  of which obo {obo}")

breakdown("z", "charge", [(1, 2), (2, 3), (3, 4), (4, 5), (5, 99)])
breakdown("m", "mass", [(0, 1500), (1500, 2500), (2500, 3500), (3500, 10000)])

# --- sequence modifications: matched vs missed ---
mm = [r for r in refs if r["cls"] == "matched"]
ms = [r for r in refs if r["cls"] != "matched"]
def avg(xs, k):
    xs = [r[k] for r in xs if r[k] is not None]
    return sum(xs) / len(xs) if xs else 0
print(f"\nmodifications / length (matched vs missed):")
print(f"  matched: avg mods {avg(mm,'mods'):.2f}  avg len {avg(mm,'len'):.1f}  (n={len(mm)})")
print(f"  missed : avg mods {avg(ms,'mods'):.2f}  avg len {avg(ms,'len'):.1f}  (n={len(ms)})")
mod_miss = sum(1 for r in refs if r["mods"] > 0 and r["cls"] != "matched")
mod_tot = sum(1 for r in refs if r["mods"] > 0)
unmod_miss = sum(1 for r in refs if r["mods"] == 0 and r["cls"] != "matched")
unmod_tot = sum(1 for r in refs if r["mods"] == 0)
print(f"  modified refs   miss {mod_miss}/{mod_tot} ({100*mod_miss/max(1,mod_tot):.0f}%)")
print(f"  unmodified refs miss {unmod_miss}/{unmod_tot} ({100*unmod_miss/max(1,unmod_tot):.0f}%)")

# --- off-by-one detail ---
obo = [r for r in refs if r["cls"] == "offbyone"]
print(f"\noff-by-one ({len(obo)}): k distribution {dict(Counter(r['k'] for r in obo))}")
print(f"  carrier single-charge: {sum(1 for r in obo if r['carrier']['ncs']<=1)}  multi: {sum(1 for r in obo if r['carrier']['ncs']>1)}")

# --- absent detail (top by intensity) ---
absent = [r for r in refs if r["cls"] == "absent"]
absent.sort(key=lambda r: -(r["int"] or 0))
print(f"\nabsent ({len(absent)}) — top 12 by intensity (mass, RT, z, intensity, mods, seq):")
for r in absent[:12]:
    print(f"  {r['m']:9.3f} RT{r['rt']:6.2f} z{r['z']} int{(r['int'] or 0):9.2e} mods{r['mods']}  {r['seq'][:42]}")
