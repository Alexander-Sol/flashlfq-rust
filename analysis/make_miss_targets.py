"""Emit analysis/misses/miss_targets.tsv — every reference peak the (best) shift-apex untargeted run
did NOT rediscover (no resolved feature within 20 ppm mass + 0.3 min RT of the ref monoisotope).

For each miss it records `reported_k`: the nearest integer 13C offset (in -1,+1,-2,+2) at which a
resolved feature DOES exist, or 0 if the peptide is genuinely absent from the resolved set. That lets
the explorer separate "off-by-one remainders" from "truly missing" peptides.

Columns feed examples/obo_envelope_probe.rs (OBO_TARGETS=...). Data paths per lumos-test-data-paths.
"""
import csv, bisect, os, re, sys

# Resolved-features file (best method) as argv[1]; default = cov90 shift-apex.
FEAT = sys.argv[1] if len(sys.argv) > 1 else r"D:\SP_Tutorial\Lumos\untargeted_cov90_shiftapex.tsv"
REF = r"D:\SP_Tutorial\Lumos\CA_HCD_GPTMD_Search_WideTol\Task2-SearchTask\AllQuantifiedPeaks.tsv"
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "misses")
OUT = os.path.join(OUTDIR, "miss_targets.tsv")
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
        if m is not None and rt is not None:
            feats.append((m, rt))
feats.sort()
keys = [m for m, _ in feats]


def has_feature(mass, rt):
    lo = bisect.bisect_left(keys, mass * (1 - MASS_PPM * 1e-6))
    hi = bisect.bisect_right(keys, mass * (1 + MASS_PPM * 1e-6))
    return any(abs(feats[i][1] - rt) <= RT_MIN for i in range(lo, hi))


def clean_seq(seq):
    return re.sub(r"\[[^\]]*\]", "", seq or "")


rows = []
n_absent = n_obo = 0
with open(REF, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        m = fnum(row.get("Peptide Monoisotopic Mass"))
        rt = fnum(row.get("Peak RT Apex"))
        z = fnum(row.get("Peak Charge"))
        pkmz = fnum(row.get("Peak MZ"))
        seq = clean_seq(row.get("Full Sequence"))
        if m is None or rt is None or z is None:
            continue
        if pkmz is None:
            pkmz = m / z + 1.007276466879  # fall back to mono m/z if Peak MZ absent
        if has_feature(m, rt):
            continue  # rediscovered at the mono — not a miss
        rk = 0
        for k in (-1, 1, -2, 2):
            if has_feature(m + k * C13, rt):
                rk = k
                break
        if rk == 0:
            n_absent += 1
        else:
            n_obo += 1
        rows.append((m, rt, int(z), pkmz, rk, seq))

os.makedirs(OUTDIR, exist_ok=True)
rows.sort(key=lambda r: (r[2], r[0]))
with open(OUT, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["label", "ref_mono", "rt", "z", "pk_mz", "reported_k"])
    for i, (m, rt, z, pkmz, rk, seq) in enumerate(rows):
        label = f"{i:02d}_{seq[:22] or 'peptide'}_z{z}"
        w.writerow([label, f"{m:.4f}", f"{rt:.3f}", z, f"{pkmz:.4f}", rk])

print(f"wrote {OUT}  ({len(rows)} misses: {n_absent} absent, {n_obo} off-by-one remainders)")
