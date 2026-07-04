"""Plot the real (observed) MS1 isotope envelopes for the monoisotope off-by-one cases against the
two competing interpretations:

  * FlashLFQ (ref)  — averagine anchored at the MetaMorpheus/FlashLFQ monoisotopic mass (k=0).
  * Pipeline        — averagine anchored where OUR untargeted detector placed the mono (k=reported_k).

Gray stems are the real centroids. Whichever colored template hugs the gray stems is where the
monoisotope actually is; the mismatch shows the off-by-one.

Consumes obo_series.tsv + obo_meta.tsv (written by examples/obo_envelope_probe.rs, EXPORT_DIR=...).
Produces obo_composite_envelopes.png (averaged composite) and obo_apex_envelopes.png (apex scan).

Run:  python analysis/plot_obo_envelopes.py analysis
"""
import csv, os, sys, math
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))

REF_COLOR = "#2ca02c"   # FlashLFQ reference
PIPE_COLOR = "#c2185b"  # our pipeline
OBS_COLOR = "0.55"

meta, order = {}, []
with open(os.path.join(DIR, "obo_meta.tsv"), newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        order.append(row["case"])
        meta[row["case"]] = {k: (int(row[k]) if k in ("z", "reported_k") else float(row[k]))
                             for k in ("ref_mono", "z", "rt", "pk_mz", "reported_k", "spacing", "mono_mz")}

series = defaultdict(lambda: defaultdict(list))  # case -> series -> [(mz, inten)]
with open(os.path.join(DIR, "obo_series.tsv"), newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        series[row["case"]][row["series"]].append((float(row["mz"]), float(row["intensity"])))


def template(case, k):
    """Exported averagine template for shift k in {-1, 0, +1}."""
    tag = "avg_mono" if k == 0 else ("avg_mono_m1" if k < 0 else "avg_mono_p1")
    return series[case][tag]


def draw(ax, case, obs_key):
    m = meta[case]
    z, sp, mono_mz, rk = m["z"], m["spacing"], m["mono_mz"], m["reported_k"]
    xmin, xmax = mono_mz - 1.6 * sp, mono_mz + 9.0 * sp

    obs = [(mz, i) for mz, i in series[case][obs_key] if xmin <= mz <= xmax]
    omax = max((i for _, i in obs), default=1.0) or 1.0
    for mz, i in obs:
        ax.vlines(mz, 0, i / omax, color=OBS_COLOR, linewidth=3.0, zorder=1)
    ax.vlines([], [], [], color=OBS_COLOR, linewidth=3.0, label="observed centroids")  # legend proxy

    # FlashLFQ reference template (k=0) and pipeline template (k=reported_k).
    for k, color, marker, name in [
        (0, REF_COLOR, "o", "FlashLFQ (ref) mono"),
        (rk, PIPE_COLOR, "X", f"pipeline mono (k={rk:+d})"),
    ]:
        pts = sorted((mz, w) for mz, w in template(case, k) if xmin - sp <= mz <= xmax + sp)
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "-", color=color, marker=marker,
                    ms=7, lw=1.7, alpha=0.9, label=name, zorder=3)

    # monoisotope m/z markers for each interpretation.
    ax.axvline(mono_mz, color=REF_COLOR, ls="--", lw=1.1, alpha=0.8)
    pipe_mono_mz = mono_mz + rk * sp
    ax.axvline(pipe_mono_mz, color=PIPE_COLOR, ls=":", lw=1.4, alpha=0.9)
    ax.annotate("FlashLFQ\nmono", (mono_mz, 1.11), color=REF_COLOR, fontsize=7.5,
                ha="center", va="top")
    ax.annotate("pipeline\nmono", (pipe_mono_mz, 1.11), color=PIPE_COLOR, fontsize=7.5,
                ha="center", va="top")

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(0, 1.22)
    ax.set_title(f"{case}   z{z}   ref_mono {m['ref_mono']:.3f} Da   "
                 f"(pipeline off by {rk:+d} 13C)", fontsize=10)
    ax.set_xlabel("m/z")
    ax.set_ylabel("rel. intensity")
    ax.legend(fontsize=7.5, loc="upper right", framealpha=0.92)
    ax.grid(True, alpha=0.15)


def make_figure(obs_key, title, outfile, ncols=4):
    n = len(order)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 3.0 * nrows), squeeze=False)
    flat = axes.flat
    for ax, case in zip(flat, order):
        draw(ax, case, obs_key)
    for ax in list(flat)[n:]:
        ax.axis("off")
    fig.suptitle(title, fontsize=14, y=0.997)
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    path = os.path.join(DIR, outfile)
    fig.savefig(path, dpi=110)
    plt.close(fig)
    print("wrote", path, f"({n} cases)")


make_figure("observed_composite",
            "Off-by-one cases - observed AVERAGED composite (gray) vs FlashLFQ-ref (green) and pipeline (pink) monoisotope",
            "obo_composite_envelopes.png")
make_figure("observed_apex",
            "Off-by-one cases - observed APEX scan (gray) vs FlashLFQ-ref (green) and pipeline (pink) monoisotope",
            "obo_apex_envelopes.png")
