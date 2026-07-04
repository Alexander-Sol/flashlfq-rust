"""Builds analysis/obo_envelope_explorer.ipynb from the exported off-by-one TSVs.

Dependency-free authoring: constructs the notebook JSON with the stdlib `json` module (nbformat is
not required for authoring). The notebook itself needs only numpy + matplotlib + stdlib csv to run;
ipywidgets (installed in the repo .venv) powers the interactive dropdown when present.

Run:  python analysis/build_obo_notebook.py
"""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "obo_envelope_explorer.ipynb")


_ID = [0]


def _next_id():
    _ID[0] += 1
    return f"obo-{_ID[0]:02d}"


def md(text):
    return {"cell_type": "markdown", "id": _next_id(), "metadata": {},
            "source": text.splitlines(keepends=True)}


def code(text):
    return {"cell_type": "code", "id": _next_id(), "metadata": {}, "execution_count": None,
            "outputs": [], "source": text.strip("\n").splitlines(keepends=True)}


cells = []

cells.append(md(r"""# Off-by-one envelope explorer

Explore **all ~47** monoisotope **off-by-one** misses of the untargeted detector. For each one,
compare the real observed MS1 isotope envelope against the two competing interpretations:

- **FlashLFQ (ref)** — averagine anchored at the MetaMorpheus/FlashLFQ monoisotopic mass (`k=0`).
- **Pipeline** — averagine anchored where our untargeted detector placed the mono (`k = reported_k`).

Whichever template hugs the gray observed centroids is where the monoisotope really is.

**Data** comes from `obo_meta.tsv`, `obo_series.tsv`, `obo_scans.tsv` in this folder, written by
`rust/flashlfq-core/examples/obo_envelope_probe.rs`. To (re)generate for all 47 cases:

```
python analysis/make_obo_targets.py                       # -> analysis/obo_targets.tsv (the 47)
cargo build --release --example obo_envelope_probe
# then, with env set:  EXPORT_DIR=analysis  OBO_TARGETS=analysis/obo_targets.tsv
#   ./target/release/examples/obo_envelope_probe
python analysis/build_obo_notebook.py                     # rebuild this notebook
```

Use the **dropdown** near the bottom for interactive browsing, or `plot_case(order[i])` directly.
"""))

cells.append(code(r"""
import csv, os, math
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

DATA_DIR = os.path.dirname(os.path.abspath("__file__")) if "__file__" in globals() else os.getcwd()
for cand in (DATA_DIR, os.path.join(DATA_DIR, "analysis"), "."):
    if os.path.exists(os.path.join(cand, "obo_meta.tsv")):
        DATA_DIR = cand
        break
print("data dir:", DATA_DIR)

meta, order = {}, []
with open(os.path.join(DATA_DIR, "obo_meta.tsv"), newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        order.append(row["case"])
        meta[row["case"]] = {k: (int(row[k]) if k in ("z", "reported_k") else float(row[k]))
                             for k in ("ref_mono","z","rt","pk_mz","reported_k","spacing","mono_mz")}

series = defaultdict(lambda: defaultdict(list))     # case -> series -> [(mz, inten)]
with open(os.path.join(DATA_DIR, "obo_series.tsv"), newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        series[row["case"]][row["series"]].append((float(row["mz"]), float(row["intensity"])))

scans = defaultdict(lambda: defaultdict(list))      # case -> scan_offset -> [(mz, inten)]
scan_rt = defaultdict(dict)                          # case -> scan_offset -> rt
with open(os.path.join(DATA_DIR, "obo_scans.tsv"), newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        off = int(row["scan_offset"])
        scans[row["case"]][off].append((float(row["mz"]), float(row["intensity"])))
        scan_rt[row["case"]][off] = float(row["rt"])

print(len(order), "off-by-one cases loaded")
"""))

cells.append(md("## The cases\n\nFull list with charge, reference mass, and the direction the pipeline was off (`reported_k`):"))

cells.append(code(r"""
print(f"{'idx':>3}  {'case':34} {'z':>2} {'ref_mono':>10} {'reported_k':>10}")
for i, c in enumerate(order):
    m = meta[c]
    print(f"{i:>3}  {c:34} {m['z']:>2} {m['ref_mono']:>10.3f} {m['reported_k']:>+10d}")
"""))

cells.append(md(r"""## A self-contained averagine model

So you can overlay a template at **any** shift (not just the three the Rust probe exported), here is a
small pure-Python averagine: pick the averagine composition for a target mass, then compute its
isotope envelope by convolving the elemental isotope patterns. It reproduces the *shape* (and the
mode-shift for heavy peptides) of the Rust `deconvolution::averagine_*` model; the validation cell
below overlays it on the exported Rust templates so you can confirm."""))

cells.append(code(r"""
PROTON = 1.007276466

# Monoisotopic masses and isotope-abundance vectors (probability by neutron offset from monoisotope).
_EL = {
    "C": (12.0,            [0.9893, 0.0107]),
    "H": (1.0078250319,    [0.999885, 0.000115]),
    "N": (14.0030740052,   [0.99636, 0.00364]),
    "O": (15.9949146221,   [0.99757, 0.00038, 0.00205]),
    "S": (31.97207069,     [0.9499, 0.0075, 0.0425, 0.0, 0.0001]),
}
# Averagine residue coefficients per unit multiplier (mzLib / deconvolution.rs magic numbers).
_AVG = {"C": 4.9384, "H": 7.7583, "O": 1.4773, "N": 1.3577, "S": 0.0417}
_RES_MASS = sum(_AVG[e] * _EL[e][0] for e in _AVG)   # ~111.05 Da per multiplier

def _conv_pow(base, n, maxlen=18):
    # Isotope pattern of n identical atoms: base convolved n times, truncated.
    dist = np.array([1.0])
    b = np.array(base)
    for _ in range(int(n)):
        dist = np.convolve(dist, b)[:maxlen]
    return dist

def averagine_weights(mono_mass, maxlen=16):
    # Per-isotope-index weights (index 0 = monoisotope) for a peptide of the given mono mass,
    # normalised so the max weight is 1.0.
    mult = max(mono_mass / _RES_MASS, 0.5)
    dist = np.array([1.0])
    for e, coeff in _AVG.items():
        n = round(coeff * mult)
        if n <= 0:
            continue
        dist = np.convolve(dist, _conv_pow(_EL[e][1], n, maxlen))[:maxlen]
    dist = dist[:maxlen]
    return dist / dist.max()

def to_mz(mass, z):
    return mass / z + PROTON

def template(mono_mass, z, kshift=0, maxlen=14, min_w=1e-3):
    # (mz, weight) teeth for the averagine envelope of mono_mass shifted by kshift 13C units, at
    # charge z. Weights max-normalised to 1.0; trailing teeth below min_w dropped.
    C13 = 1.0033548381
    anchor = mono_mass + kshift * C13
    w = averagine_weights(anchor, maxlen)
    mode = int(np.argmax(w))
    end = len(w)
    while end > mode + 1 and w[end-1] < min_w:
        end -= 1
    w = w[:end]
    mz = [to_mz(anchor + k * C13, z) for k in range(len(w))]
    return list(zip(mz, w))

print("averagine residue mass ~ %.3f Da/multiplier" % _RES_MASS)
"""))

cells.append(code(r"""
# Validation: Python averagine vs the exact Rust templates the probe exported (first few cases).
show = order[:4]
fig, axes = plt.subplots(1, len(show), figsize=(4*len(show), 3.2))
for ax, case in zip(np.atleast_1d(axes), show):
    m = meta[case]
    rust = sorted(series[case]["avg_mono"])
    py = template(m["ref_mono"], m["z"], 0)
    ax.plot([p[0] for p in rust], [p[1] for p in rust], "o-", label="Rust", alpha=.7)
    ax.plot([p[0] for p in py], [p[1] for p in py], "x--", label="Python", alpha=.7)
    ax.set_title(case, fontsize=7); ax.legend(fontsize=7)
fig.suptitle("Averagine template: Rust (exact) vs Python (this notebook)"); fig.tight_layout()
"""))

cells.append(md(r"""## Plotting functions

- `plot_case(case, observed='observed_composite', extra_shifts=(), legend=True, ax=None)` — observed
  envelope vs FlashLFQ-ref (green) and pipeline (pink); pass `observed='observed_apex'` for the apex
  scan, or `extra_shifts=[+2, -2]` to overlay more anchors.
- `plot_elution(case)` — every scan in the apex ±3 window as small multiples."""))

cells.append(code(r"""
REF_C, PIPE_C, OBS_C = "#2ca02c", "#c2185b", "0.55"

def plot_case(case, observed="observed_composite", extra_shifts=(), legend=True, ax=None):
    m = meta[case]; z, sp, mono_mz, rk = m["z"], m["spacing"], m["mono_mz"], m["reported_k"]
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 5))
    xmin, xmax = mono_mz - 1.6*sp, mono_mz + 9*sp
    obs = [(mz, i) for mz, i in series[case][observed] if xmin <= mz <= xmax]
    omax = max((i for _, i in obs), default=1.0) or 1.0
    for mz, i in obs:
        ax.vlines(mz, 0, i/omax, color=OBS_C, lw=3, zorder=1)
    ax.vlines([], [], [], color=OBS_C, lw=3, label="observed")

    curves = [(0, REF_C, "o", "FlashLFQ (ref)"),
              (rk, PIPE_C, "X", f"pipeline (k={rk:+d})")]
    for k in extra_shifts:
        curves.append((k, None, ".", f"k={k:+d}"))
    for k, color, marker, name in curves:
        pts = [(mz, w) for mz, w in template(m["ref_mono"], z, k) if xmin-sp <= mz <= xmax+sp]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "-", marker=marker, ms=6, lw=1.6,
                alpha=.9, color=color, label=name, zorder=3)
    ax.axvline(mono_mz, color=REF_C, ls="--", lw=1.0, alpha=.8)
    ax.axvline(mono_mz + rk*sp, color=PIPE_C, ls=":", lw=1.3, alpha=.9)
    ax.set_xlim(xmin, xmax); ax.set_ylim(0, 1.22)
    ax.set_title(f"{case}\nz{z}  {m['ref_mono']:.3f} Da  (off by {rk:+d} 13C)", fontsize=8)
    ax.set_xlabel("m/z", fontsize=8); ax.set_ylabel("rel. int", fontsize=8)
    if legend:
        ax.legend(fontsize=7, loc="upper right")
    ax.grid(True, alpha=.15)
    return ax

def plot_elution(case):
    m = meta[case]; z, sp, mono_mz, rk = m["z"], m["spacing"], m["mono_mz"], m["reported_k"]
    offs = sorted(scans[case]); n = len(offs)
    fig, axes = plt.subplots(1, n, figsize=(3.0*n, 3.2), sharey=True)
    xmin, xmax = mono_mz - 1.6*sp, mono_mz + 9*sp
    for ax, off in zip(np.atleast_1d(axes), offs):
        pts = [(mz, i) for mz, i in scans[case][off] if xmin <= mz <= xmax]
        omax = max((i for _, i in pts), default=1.0) or 1.0
        for mz, i in pts:
            ax.vlines(mz, 0, i/omax, color=OBS_C, lw=2.4)
        ax.axvline(mono_mz, color=REF_C, ls="--", lw=1.0)
        ax.axvline(mono_mz + rk*sp, color=PIPE_C, ls=":", lw=1.2)
        ax.set_xlim(xmin, xmax); ax.set_ylim(0, 1.15)
        tag = "APEX" if off == 0 else f"apex{off:+d}"
        ax.set_title(f"{tag}\nRT {scan_rt[case][off]:.3f}", fontsize=8)
        ax.set_xlabel("m/z", fontsize=8)
    np.atleast_1d(axes)[0].set_ylabel("rel. int")
    fig.suptitle(f"{case}  z{z}  — green=FlashLFQ mono, pink=pipeline mono", y=1.03, fontsize=10)
    fig.tight_layout()

# Quick look at the first case:
plot_case(order[0]);
"""))

cells.append(md(r"""## Interactive explorer

Pick any case from the dropdown. Shows the composite + apex + full elution. (Needs ipywidgets — the
repo `.venv` has it. If it is missing, use `plot_case(order[i])` / `plot_elution(order[i])` manually.)"""))

cells.append(code(r"""
try:
    from ipywidgets import interact, Dropdown
    def _explore(case):
        fig, ax = plt.subplots(1, 2, figsize=(15, 4.6))
        plot_case(case, "observed_composite", ax=ax[0]); ax[0].set_title("composite\n" + ax[0].get_title(), fontsize=8)
        plot_case(case, "observed_apex", ax=ax[1]);      ax[1].set_title("apex scan\n" + ax[1].get_title(), fontsize=8)
        fig.tight_layout(); plt.show()
        plot_elution(case); plt.show()
    interact(_explore, case=Dropdown(options=order, value=order[0], description="case"))
except Exception as e:
    print("ipywidgets not available (%s) — use plot_case(order[i]) manually." % e)
"""))

cells.append(md("## Overview grid — every case at a glance (averaged composite)"))

cells.append(code(r"""
n = len(order); ncols = 4; nrows = math.ceil(n / ncols)
fig, axes = plt.subplots(nrows, ncols, figsize=(4.4*ncols, 2.9*nrows))
for ax, case in zip(axes.flat, order):
    plot_case(case, "observed_composite", legend=False, ax=ax)
for ax in axes.flat[n:]:
    ax.axis("off")
fig.tight_layout()
"""))

cells.append(md("## Free exploration\n\nOverlay whatever shifts you like, on any case index:"))

cells.append(code(r"""
i = 0                       # <- change the index
plot_case(order[i], "observed_composite", extra_shifts=[-2, +1]);
plot_elution(order[i]);
"""))

nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "flashlfq (.venv)", "language": "python", "name": "flashlfq"},
        "language_info": {"name": "python", "version": "3"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

with open(OUT, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1)
print("wrote", OUT, "(%d cells)" % len(cells))
