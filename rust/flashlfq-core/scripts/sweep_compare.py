"""Compare the FWHM sweep: for each cov99_fwhm{N}.tsv, report resolved count, RT-width percentiles,
overall recall (m/z merge + neutral-mass), and SHORT-feature (<=3 scan) recall on neutral mass."""
import csv, bisect

REF = r"D:\SP_Tutorial\Lumos\CA_HCD_GPTMD_Search_WideTol\Task2-SearchTask\AllQuantifiedPeaks.tsv"
FWHMS = ["1p8", "3", "8", "12", "18", "36"]
LABEL = {"1p8": "1.8", "3": "3", "8": "8", "12": "12", "18": "18", "36": "36"}
FEAT = lambda n: rf"D:\SP_Tutorial\Lumos\untargeted_features_10min_avg_cov99_fwhm{n}.tsv"
PROTON = 1.007276466
SPACING_MIN = 0.0064
MZ_PPM, RT_MIN, MASS_PPM = 10.0, 0.3, 20.0


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


# reference: observed Peak MZ, theoretical mass, apex RT, scan estimate
refs = []
with open(REF, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        pkmz = fnum(row.get("Peak MZ"))
        mass = fnum(row.get("Peptide Monoisotopic Mass"))
        rt = fnum(row.get("Peak RT Apex"))
        s, e = fnum(row.get("Peak RT Start")), fnum(row.get("Peak RT End"))
        if rt is None or s is None or e is None:
            continue
        refs.append({"pkmz": pkmz, "mass": mass, "rt": rt, "scans": (e - s) / SPACING_MIN + 1})

short_refs = [r for r in refs if r["scans"] <= 3]
long_refs = [r for r in refs if r["scans"] > 3]


def load(path):
    feats = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            feats.append({
                "abmz": fnum(row.get("Most-Abundant m/z")),
                "mass": fnum(row.get("Monoisotopic Mass")),
                "rt": fnum(row.get("RT Apex")),
                "w": (fnum(row.get("RT End")) or 0) - (fnum(row.get("RT Start")) or 0),
            })
    return feats


def indexed(feats, key):
    xs = sorted((f[key], f["rt"]) for f in feats if f[key] is not None and f["rt"] is not None)
    return xs, [x[0] for x in xs]


def recall(refs_sub, key_ref, xs, keys, ppm):
    hit = 0
    for r in refs_sub:
        v = r[key_ref]
        if v is None:
            continue
        lo = bisect.bisect_left(keys, v * (1 - ppm * 1e-6))
        hi = bisect.bisect_right(keys, v * (1 + ppm * 1e-6))
        if any(abs(rt - r["rt"]) <= RT_MIN for _, rt in xs[lo:hi]):
            hit += 1
    return hit


def pct(xs, ps):
    xs = sorted(xs)
    n = len(xs)
    return [round(xs[min(n - 1, int(p / 100 * (n - 1)))], 3) for p in ps]


print(f"refs: {len(refs)} total, {len(short_refs)} short (<=3 scan), {len(long_refs)} long\n")
hdr = (f'{"FWHM":>5} {"feats":>7} {"wMed":>6} {"wP99":>6} {"wMax":>7} | '
       f'{"mz%":>6} {"mzChg":>6} | {"mass%":>6} | {"SHORT mass%":>11} {"LONG mass%":>10}')
print(hdr)
print("-" * len(hdr))
for n in FWHMS:
    lbl = LABEL[n]
    try:
        feats = load(FEAT(n))
    except FileNotFoundError:
        print(f'{lbl:>5}  (not written yet)')
        continue
    widths = [f["w"] for f in feats if f["w"] is not None]
    wmed, wp99, wmax = pct(widths, [50, 99, 100])
    xs_mz, k_mz = indexed(feats, "abmz")
    xs_ma, k_ma = indexed(feats, "mass")
    mz_all = recall(refs, "pkmz", xs_mz, k_mz, MZ_PPM)
    mass_all = recall(refs, "mass", xs_ma, k_ma, MASS_PPM)
    mass_short = recall(short_refs, "mass", xs_ma, k_ma, MASS_PPM)
    mass_long = recall(long_refs, "mass", xs_ma, k_ma, MASS_PPM)
    # charge match on m/z: recount with charge check would need charge; skip, report mz recall only
    print(f'{lbl:>5} {len(feats):>7} {wmed:>6.3f} {wp99:>6.3f} {wmax:>7.3f} | '
          f'{100*mz_all/len(refs):>5.1f}% {"":>6} | {100*mass_all/len(refs):>5.1f}% | '
          f'{mass_short:>3}/{len(short_refs)} {100*mass_short/len(short_refs):>5.1f}% '
          f'{100*mass_long/len(long_refs):>9.1f}%')
