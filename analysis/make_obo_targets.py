"""Emit analysis/obo_targets.tsv — every reference peak the untargeted detector got off by one 13C
unit (a feature exists at ref +- k*C13 but not at the ref mono). Columns feed
examples/obo_envelope_probe.rs (OBO_TARGETS=...), which then extracts the real envelopes.

Same join as offbyone_diag.py / char_misses.py. Data paths hard-coded (see lumos-test-data-paths).
"""
import csv, bisect, os, re

FEAT = r"D:\SP_Tutorial\Lumos\untargeted_features_10min_centroid_dd_cov99_minscan2.tsv"
REF = r"D:\SP_Tutorial\Lumos\CA_HCD_GPTMD_Search_WideTol\Task2-SearchTask\AllQuantifiedPeaks.tsv"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "obo_targets.tsv")
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
    return re.sub(r"\[[^\]]*\]", "", seq or "")  # drop [mod] annotations for a compact label


rows = []
with open(REF, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        m = fnum(row.get("Peptide Monoisotopic Mass"))
        rt = fnum(row.get("Peak RT Apex"))
        z = fnum(row.get("Peak Charge"))
        pkmz = fnum(row.get("Peak MZ"))
        seq = clean_seq(row.get("Full Sequence"))
        if m is None or rt is None or pkmz is None or z is None:
            continue
        if has_feature(m, rt):
            continue  # matched on the mono; not an off-by-one
        rk = None
        for k in (-1, 1, -2, 2):
            if has_feature(m + k * C13, rt):
                rk = k
                break
        if rk is None:
            continue  # genuinely absent, not an off-by-one
        rows.append((m, rt, int(z), pkmz, rk, seq))

# Sort by charge then mass for a sensible browsing order; assign unique labels.
rows.sort(key=lambda r: (r[2], r[0]))
with open(OUT, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["label", "ref_mono", "rt", "z", "pk_mz", "reported_k"])
    for i, (m, rt, z, pkmz, rk, seq) in enumerate(rows):
        label = f"{i:02d}_{seq[:22] or 'peptide'}_z{z}"
        w.writerow([label, f"{m:.4f}", f"{rt:.3f}", z, f"{pkmz:.4f}", rk])

print(f"wrote {OUT}  ({len(rows)} off-by-one targets)")
