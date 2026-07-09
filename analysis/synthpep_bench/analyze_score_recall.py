#!/usr/bin/env python
"""Score -> recall tradeoff for our refined features on the PXD001091 benchmark.

Question: we emit ~7.5x more features than Dinosaur. Are the extra features a
low-score tail we could REJECT IN THE REFINER without losing recall?

We work at the *refined per-charge* level (feat.refined.tsv) -- that is where a
rejection threshold would be applied, before charge-consensus merge. Two candidate
rejection axes:
    Decon Score      -- envelope-fit cosine from refine (col 'Decon Score')
    Summed Intensity -- col 'Summed Intensity'

For every GT peak (177,759 across 57 files) we find all refined features within
+-20 ppm mass and +-0.5 min RT (charge-blind, matching score_all's recall), and
record the MAX score among them: a peak survives a rejection threshold t iff its
best matching feature scores >= t. recall(t) = frac of peaks whose best match >= t.
features_retained(t) = frac of ALL refined features with score >= t.

Also: every refined feature is labelled matched / unmatched (does it match ANY GT
peak) so we can see whether unmatched "junk" clusters at low score.

Outputs: score_recall.png  +  score_recall_curves.tsv  + printed threshold tables.
"""
import csv, os, sys, bisect
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BENCH = r"F:\flashlfq-rust\analysis\synthpep_bench"
MASS_PPM = 20.0
RT_DELTA = 0.5

# --- bins ---
DEC_BINS = np.linspace(0.0, 1.0, 501)                 # decon score in [0,1]
LOGI_LO, LOGI_HI = 2.0, 12.0
INT_BINS = np.linspace(LOGI_LO, LOGI_HI, 501)         # log10(intensity)


def read_refined(path):
    m, z, rt, inten, dec = [], [], [], [], []
    with open(path, newline="", encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            m.append(p[0]); z.append(p[1]); rt.append(p[2]); inten.append(p[3]); dec.append(p[4])
    return (np.asarray(m, float), np.asarray(z, int), np.asarray(rt, float),
            np.asarray(inten, float), np.asarray(dec, float))


def read_gt(path):
    mass, rt = [], []
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            mass.append(float(row["mono_mass"])); rt.append(float(row["rt"]))
    return np.asarray(mass, float), np.asarray(rt, float)


def main():
    bases = [r["base"] for r in csv.DictReader(
        open(os.path.join(BENCH, "gt", "gt_manifest.tsv"), encoding="utf-8"), delimiter="\t")]

    # global accumulators
    feat_dec_all = np.zeros(len(DEC_BINS) - 1)     # all refined features, decon hist
    feat_dec_match = np.zeros(len(DEC_BINS) - 1)   # matched features, decon hist
    feat_int_all = np.zeros(len(INT_BINS) - 1)
    feat_int_match = np.zeros(len(INT_BINS) - 1)
    peak_best_dec = []                              # per GT peak: best-match decon (-1 if unmatched)
    peak_best_int = []                              # per GT peak: best-match log10 intensity (-1 if unmatched)
    n_feat_total = 0
    n_gt_total = 0

    for i, base in enumerate(bases, 1):
        rp = os.path.join(BENCH, "ours", base, "feat.refined.tsv")
        gp = os.path.join(BENCH, "gt", base + ".tsv")
        if not (os.path.exists(rp) and os.path.exists(gp)):
            continue
        fm, fz, frt, fi, fd = read_refined(rp)
        order = np.argsort(fm)
        fm, fz, frt, fi, fd = fm[order], fz[order], frt[order], fi[order], fd[order]
        fli = np.log10(np.clip(fi, 10 ** LOGI_LO, None))
        gm, grt = read_gt(gp)
        n_feat_total += len(fm)
        n_gt_total += len(gm)

        matched_mask = np.zeros(len(fm), dtype=bool)
        for k in range(len(gm)):
            rmass, rrt = gm[k], grt[k]
            tol = rmass * MASS_PPM / 1e6
            lo = bisect.bisect_left(fm, rmass - tol)
            hi = bisect.bisect_right(fm, rmass + tol)
            if hi > lo:
                sl = slice(lo, hi)
                rt_ok = np.abs(frt[sl] - rrt) <= RT_DELTA
                if rt_ok.any():
                    idx = np.arange(lo, hi)[rt_ok]
                    matched_mask[idx] = True
                    peak_best_dec.append(fd[idx].max())
                    peak_best_int.append(fli[idx].max())
                    continue
            peak_best_dec.append(-1.0)
            peak_best_int.append(-1.0)

        feat_dec_all += np.histogram(fd, bins=DEC_BINS)[0]
        feat_dec_match += np.histogram(fd[matched_mask], bins=DEC_BINS)[0]
        feat_int_all += np.histogram(fli, bins=INT_BINS)[0]
        feat_int_match += np.histogram(fli[matched_mask], bins=INT_BINS)[0]
        print(f"  [{i}/{len(bases)}] {base}: {len(fm)} refined, {len(gm)} gt", flush=True)

    peak_best_dec = np.asarray(peak_best_dec)
    peak_best_int = np.asarray(peak_best_int)
    full_recall = np.mean(peak_best_dec >= 0.0)  # matched at all

    # ---- curves ----
    def retained_curve(hist, edges):
        # fraction of features with score >= edge (right-cumulative), evaluated at left edges
        rev = np.cumsum(hist[::-1])[::-1]
        total = hist.sum()
        return edges[:-1], rev / max(total, 1)

    def recall_curve(best, thresholds):
        return np.array([np.mean(best >= t) for t in thresholds])

    dec_edges = DEC_BINS[:-1]
    dec_ret = retained_curve(feat_dec_all, DEC_BINS)[1]
    dec_rec = recall_curve(peak_best_dec, dec_edges)

    int_edges = INT_BINS[:-1]
    int_ret = retained_curve(feat_int_all, INT_BINS)[1]
    int_rec = recall_curve(peak_best_int, int_edges)

    # ---- threshold tables: max threshold keeping recall >= full - loss ----
    def pick(thresholds, rec, ret, loss):
        ok = np.where(rec >= full_recall - loss)[0]
        j = ok.max()  # largest threshold index still satisfying
        return thresholds[j], rec[j], ret[j]

    print(f"\nTotal refined features: {n_feat_total:,}   GT peaks: {n_gt_total:,}")
    print(f"Full recall (any matching refined feature, +-{RT_DELTA} min): {100*full_recall:.2f}%\n")
    print("DECON SCORE rejection (reject features below threshold):")
    for loss in (0.001, 0.005, 0.010):
        t, rc, rt_ = pick(dec_edges, dec_rec, dec_ret, loss)
        print(f"  keep recall >= {100*(full_recall-loss):.2f}% : thr Decon>={t:.3f}  "
              f"recall {100*rc:.2f}%  features kept {100*rt_:.1f}%  (dropped {100*(1-rt_):.1f}%)")
    print("\nSUMMED INTENSITY rejection (reject features below threshold):")
    for loss in (0.001, 0.005, 0.010):
        t, rc, rt_ = pick(int_edges, int_rec, int_ret, loss)
        print(f"  keep recall >= {100*(full_recall-loss):.2f}% : thr Intensity>=1e{t:.2f}  "
              f"recall {100*rc:.2f}%  features kept {100*rt_:.1f}%  (dropped {100*(1-rt_):.1f}%)")

    # ---- dump curves ----
    with open(os.path.join(BENCH, "score_recall_curves.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["decon_thr", "decon_recall", "decon_feat_kept",
                    "log10int_thr", "int_recall", "int_feat_kept"])
        for a in range(len(dec_edges)):
            w.writerow([f"{dec_edges[a]:.4f}", f"{dec_rec[a]:.5f}", f"{dec_ret[a]:.5f}",
                        f"{int_edges[a]:.4f}", f"{int_rec[a]:.5f}", f"{int_ret[a]:.5f}"])

    # ---- plot ----
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25,
                         "figure.dpi": 130})
    C_REC, C_RET, C_MATCH, C_JUNK = "#1b6ca8", "#c0392b", "#2e8b57", "#d98c00"
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    def tradeoff(a, thr, rec, ret, xlabel, title, logx=False):
        a.plot(thr, 100 * rec, color=C_REC, lw=2, label="recall")
        a.set_xlabel(xlabel); a.set_ylabel("recall (%)", color=C_REC)
        a.tick_params(axis="y", labelcolor=C_REC)
        a.set_ylim(0, 100)
        a2 = a.twinx()
        a2.plot(thr, 100 * ret, color=C_RET, lw=2, ls="--", label="features kept")
        a2.set_ylabel("refined features kept (%)", color=C_RET)
        a2.tick_params(axis="y", labelcolor=C_RET); a2.set_ylim(0, 100); a2.grid(False)
        a.axhline(100 * full_recall, color="gray", ls=":", lw=1)
        a.set_title(title)
        if logx:
            a.set_xlim(thr[0], thr[-1])

    tradeoff(ax[0, 0], dec_edges, dec_rec, dec_ret, "Decon Score rejection threshold",
             "Reject on Decon Score (envelope-fit cosine)")
    tradeoff(ax[0, 1], int_edges, int_rec, int_ret, "log10(Summed Intensity) rejection threshold",
             "Reject on Summed Intensity")

    # histograms matched vs unmatched
    dc = DEC_BINS[:-1]
    ax[1, 0].plot(dc, feat_dec_match, color=C_MATCH, lw=1.6, label="matches a GT peak")
    ax[1, 0].plot(dc, feat_dec_all - feat_dec_match, color=C_JUNK, lw=1.6, label="no GT match")
    ax[1, 0].set_yscale("log"); ax[1, 0].set_xlabel("Decon Score")
    ax[1, 0].set_ylabel("refined features (log)"); ax[1, 0].legend()
    ax[1, 0].set_title("Decon Score: matched vs unmatched features")

    ic = INT_BINS[:-1]
    ax[1, 1].plot(ic, feat_int_match, color=C_MATCH, lw=1.6, label="matches a GT peak")
    ax[1, 1].plot(ic, feat_int_all - feat_int_match, color=C_JUNK, lw=1.6, label="no GT match")
    ax[1, 1].set_yscale("log"); ax[1, 1].set_xlabel("log10(Summed Intensity)")
    ax[1, 1].set_ylabel("refined features (log)"); ax[1, 1].legend()
    ax[1, 1].set_title("Summed Intensity: matched vs unmatched features")

    fig.suptitle(f"Our refined features: score -> recall  (57 files, {n_feat_total:,} refined, "
                 f"{n_gt_total:,} GT peaks, full recall {100*full_recall:.2f}%)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = os.path.join(BENCH, "score_recall.png")
    fig.savefig(out)
    print(f"\nplot -> {out}")


if __name__ == "__main__":
    main()
