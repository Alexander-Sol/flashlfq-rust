"""For each reference peak NOT matched on neutral mass (+-20ppm, +-0.3min), check whether a resolved
feature sits at the right RT with mass = ref_mass + k*C13 for k in {-2,-1,+1,+2} -- i.e. the miss is a
recoverable OFF-BY-ONE, vs genuinely absent. Tells us if fixing the corrector will move recall."""
import csv, bisect

FEAT = r"D:\SP_Tutorial\Lumos\untargeted_features_10min_centroid_dd_cov99_minscan2.tsv"
REF = r"D:\SP_Tutorial\Lumos\CA_HCD_GPTMD_Search_WideTol\Task2-SearchTask\AllQuantifiedPeaks.tsv"
C13 = 1.00335483810
MASS_PPM, RT_MIN = 20.0, 0.3


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


feats = []
with open(FEAT, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        m, rt = fnum(row.get("Monoisotopic Mass")), fnum(row.get("RT Apex"))
        z = fnum(row.get("Primary Charge"))
        if m is not None and rt is not None:
            feats.append({"m": m, "rt": rt, "z": int(z) if z else 0})
feats.sort(key=lambda x: x["m"])
keys = [x["m"] for x in feats]


def find(mass, rt):
    lo = bisect.bisect_left(keys, mass * (1 - MASS_PPM * 1e-6))
    hi = bisect.bisect_right(keys, mass * (1 + MASS_PPM * 1e-6))
    for fe in feats[lo:hi]:
        if abs(fe["rt"] - rt) <= RT_MIN:
            return fe
    return None


refs = []
with open(REF, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        m = fnum(row.get("Peptide Monoisotopic Mass"))
        rt = fnum(row.get("Peak RT Apex"))
        z = fnum(row.get("Peak Charge"))
        if m is not None and rt is not None:
            refs.append({"m": m, "rt": rt, "z": int(z) if z else 0})

matched = obo = absent = 0
obo_by_k = {}
absent_examples = []
for r in refs:
    if find(r["m"], r["rt"]):
        matched += 1
        continue
    hit_k = None
    for k in (-1, 1, -2, 2):
        fe = find(r["m"] + k * C13, r["rt"])
        if fe:
            hit_k = k
            break
    if hit_k is not None:
        obo += 1
        obo_by_k[hit_k] = obo_by_k.get(hit_k, 0) + 1
    else:
        absent += 1
        if len(absent_examples) < 10:
            absent_examples.append((r["m"], r["rt"], r["z"]))

n = len(refs)
print(f"reference peaks: {n}")
print(f"  matched on mass         : {matched} ({100*matched/n:.1f}%)")
print(f"  RECOVERABLE off-by-one  : {obo} ({100*obo/n:.1f}%)   by k(*C13): {dict(sorted(obo_by_k.items()))}")
print(f"  genuinely absent        : {absent} ({100*absent/n:.1f}%)")
print(f"\n  ceiling if off-by-one fully fixed: {100*(matched+obo)/n:.1f}% recall")
print("\n  absent examples (mass, RT, z):")
for m, rt, z in absent_examples:
    print(f"    {m:.4f}  RT {rt:.3f}  z{z}")
