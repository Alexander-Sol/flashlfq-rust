"""Emit analysis/gains/gain_targets.tsv — a curated sample of the features CLASSIC refine drops but
SHIFT refine keeps (from <run>.refinediff.tsv, written by detect_features_tsv DIFF_REFINE=1).

Two groups, both with reported_k=0 (these are retained features, not off-by-ones):
  * ref_<seq>  — a gain whose shift-placed mono matches a FlashLFQ reference peak (20 ppm + 0.3 min):
                 a REAL peptide classic was too strict to keep (the recall win).
  * top_<n>    — the highest-intensity gains with NO reference match: strong envelopes classic still
                 dropped (worth seeing why).

Columns feed examples/obo_envelope_probe.rs. Paths per lumos-test-data-paths.
"""
import csv, bisect, os, re

DIFF = r"D:\SP_Tutorial\Lumos\untargeted_cov90_diff.refinediff.tsv"
REF = r"D:\SP_Tutorial\Lumos\CA_HCD_GPTMD_Search_WideTol\Task2-SearchTask\AllQuantifiedPeaks.tsv"
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gains")
OUT = os.path.join(OUTDIR, "gain_targets.tsv")
PROTON = 1.007276466879
MASS_PPM, RT_MIN = 20.0, 0.3
N_TOP = 12  # how many strong non-reference gains to include


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def clean_seq(seq):
    return re.sub(r"\[[^\]]*\]", "", seq or "")


# Reference peaks, indexed by mono for fast ppm lookup.
refs = []
with open(REF, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        m, rt, z = fnum(row.get("Peptide Monoisotopic Mass")), fnum(row.get("Peak RT Apex")), fnum(row.get("Peak Charge"))
        if m is None or rt is None or z is None:
            continue
        refs.append((m, rt, int(z), clean_seq(row.get("Full Sequence"))))
refs.sort()
ref_keys = [m for m, *_ in refs]


def ref_match(mass, rt, z):
    lo = bisect.bisect_left(ref_keys, mass * (1 - MASS_PPM * 1e-6))
    hi = bisect.bisect_right(ref_keys, mass * (1 + MASS_PPM * 1e-6))
    for i in range(lo, hi):
        rm, rrt, rz, seq = refs[i]
        if abs(rrt - rt) <= RT_MIN and rz == z:
            return seq
    return None


gains = []
with open(DIFF, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        m, rt, z = fnum(row.get("Mono")), fnum(row.get("RT")), fnum(row.get("Charge"))
        pkmz, inten = fnum(row.get("Pk MZ")), fnum(row.get("Summed Intensity"))
        if m is None or rt is None or z is None or pkmz is None:
            continue
        gains.append({"mono": m, "rt": rt, "z": int(z), "pkmz": pkmz, "inten": inten or 0.0})

matched = []
for g in gains:
    seq = ref_match(g["mono"], g["rt"], g["z"])
    if seq:
        matched.append((g, seq))

unmatched = [g for g in gains if ref_match(g["mono"], g["rt"], g["z"]) is None]
unmatched.sort(key=lambda g: -g["inten"])
top = unmatched[:N_TOP]

os.makedirs(OUTDIR, exist_ok=True)
rows = []
for g, seq in matched:
    rows.append((f"ref_{seq[:20]}_z{g['z']}", g))
for i, g in enumerate(top):
    rows.append((f"top{i:02d}_strong_z{g['z']}", g))

with open(OUT, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["label", "ref_mono", "rt", "z", "pk_mz", "reported_k"])
    for label, g in rows:
        w.writerow([label, f"{g['mono']:.4f}", f"{g['rt']:.3f}", g["z"], f"{g['pkmz']:.4f}", 0])

print(f"wrote {OUT}")
print(f"  total gains (classic-dropped, shift-kept): {len(gains)}")
print(f"  reference-matched (real peptides classic dropped): {len(matched)}")
print(f"  + top {len(top)} strong non-reference gains")
