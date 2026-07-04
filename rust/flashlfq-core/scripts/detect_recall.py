"""DETECTION-level recall vs FWHM: match reference peaks against each run's .detected.tsv (pre-refine,
pre-consensus), on neutral monoisotopic mass. This is the fair comparison for the narrow FWHM settings
whose refine/consensus is intractable (millions of features)."""
import csv, bisect, os

REF = r"D:\SP_Tutorial\Lumos\CA_HCD_GPTMD_Search_WideTol\Task2-SearchTask\AllQuantifiedPeaks.tsv"
TAGS = ["1p8", "3", "8", "12", "18", "36"]
LABEL = {"1p8": "1.8", "3": "3", "8": "8", "12": "12", "18": "18", "36": "36"}
DET = lambda t: rf"D:\SP_Tutorial\Lumos\untargeted_features_10min_avg_cov99_fwhm{t}.detected.tsv"
SPACING_MIN = 0.0064
MASS_PPM, RT_MIN = 20.0, 0.3


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


refs = []
with open(REF, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        m = fnum(row.get("Peptide Monoisotopic Mass"))
        rt = fnum(row.get("Peak RT Apex"))
        s, e = fnum(row.get("Peak RT Start")), fnum(row.get("Peak RT End"))
        if m is None or rt is None or s is None or e is None:
            continue
        refs.append({"mass": m, "rt": rt, "scans": (e - s) / SPACING_MIN + 1})
short = [r for r in refs if r["scans"] <= 3]
long = [r for r in refs if r["scans"] > 3]


def load(path):
    xs = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            m, rt = fnum(row.get("Monoisotopic Mass")), fnum(row.get("Apex RT"))
            if m is not None and rt is not None:
                xs.append((m, rt))
    xs.sort()
    return xs, [x[0] for x in xs]


def recall(sub, xs, keys):
    hit = 0
    for r in sub:
        lo = bisect.bisect_left(keys, r["mass"] * (1 - MASS_PPM * 1e-6))
        hi = bisect.bisect_right(keys, r["mass"] * (1 + MASS_PPM * 1e-6))
        if any(abs(rt - r["rt"]) <= RT_MIN for _, rt in xs[lo:hi]):
            hit += 1
    return hit


print(f"DETECTION-level recall (ref theoretical mass +-{MASS_PPM:.0f} ppm, +-{RT_MIN} min RT vs .detected.tsv)")
print(f"refs: {len(refs)} total, {len(short)} short (<=3 scan), {len(long)} long\n")
hdr = f'{"FWHM":>5} {"detected":>10} | {"ALL":>13} | {"SHORT":>12} | {"LONG":>12}'
print(hdr)
print("-" * len(hdr))
for t in TAGS:
    p = DET(t)
    if not os.path.exists(p):
        print(f'{LABEL[t]:>5}  (not present)')
        continue
    xs, keys = load(p)
    ra, rs, rl = recall(refs, xs, keys), recall(short, xs, keys), recall(long, xs, keys)
    print(f'{LABEL[t]:>5} {len(xs):>10} | {ra:>4}/{len(refs)} {100*ra/len(refs):>5.1f}% | '
          f'{rs:>2}/{len(short)} {100*rs/len(short):>5.1f}% | {rl:>3}/{len(long)} {100*rl/len(long):>5.1f}%')
