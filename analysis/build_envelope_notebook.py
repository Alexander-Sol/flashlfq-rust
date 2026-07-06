"""Builds an interactive envelope-explorer notebook from probe-exported TSVs (obo_meta/series/scans).

One builder, three modes (argv[1]): `obo` | `miss` | `gain`. Each reads its own data subfolder and
writes a notebook there. The plotting/averagine machinery is identical across modes; only the intro,
title, and data location differ.

Run:  python analysis/build_envelope_notebook.py miss
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
MODE = sys.argv[1] if len(sys.argv) > 1 else "obo"

INTRO = {
    "obo": r"""# Off-by-one envelope explorer

Explore the monoisotope **off-by-one** misses of the untargeted detector. For each, compare the real
observed MS1 isotope envelope against **FlashLFQ (ref, k=0)** and the **pipeline** anchor
(`k = reported_k`). Whichever template hugs the gray observed centroids is where the mono really is.""",
    "miss": r"""# Miss envelope explorer

Every reference peak the **best (shift-apex) untargeted run still does not rediscover** (no resolved
feature within 20 ppm + 0.3 min of the ref mono). `reported_k` is the nearest 13C offset at which a
feature *does* exist (**0 = genuinely absent**; non-zero = an **off-by-one remainder**).

For each miss: gray = observed MS1 centroids at the reference RT; **green** = averagine at the
FlashLFQ mono (`k=0`); **pink** = averagine at `k=reported_k` (where a feature was actually found).
Use it to see *why* each peptide is missed — no signal at all, a co-eluting interferent, or a
surviving off-by-one.""",
    "gain": r"""# Shift-retained feature explorer

Features the detector found that **classic deconvolution drops** (finds no envelope) but the
**detector-anchored shift decon keeps** — the completeness gain behind the +2.4% recall. Each target
is a retained feature (mono = its shift-placed monoisotope). Gray = observed MS1 centroids; **green**
= the averagine template at that mono. Use it to judge whether the retained envelopes are real
peptides classic was too strict to keep, or noise.""",
}[MODE]

TITLE_TABLE = {
    "obo": "the off-by-one direction (`reported_k`)",
    "miss": "the nearest feature offset (`reported_k`; 0 = absent)",
    "gain": "the shift-placed mono (`reported_k` is 0 here)",
}[MODE]

SUBDIR = {"obo": "", "miss": "misses", "gain": "gains"}[MODE]
FNAME = {"obo": "obo_envelope_explorer.ipynb",
         "miss": "miss_envelope_explorer.ipynb",
         "gain": "gain_envelope_explorer.ipynb"}[MODE]
DATA_LOCAL = os.path.join(HERE, SUBDIR) if SUBDIR else HERE
OUT = os.path.join(DATA_LOCAL, FNAME)

_ID = [0]
def _next_id():
    _ID[0] += 1
    return f"env-{_ID[0]:02d}"

def md(text):
    return {"cell_type": "markdown", "id": _next_id(), "metadata": {},
            "source": text.splitlines(keepends=True)}

def code(text):
    return {"cell_type": "code", "id": _next_id(), "metadata": {}, "execution_count": None,
            "outputs": [], "source": text.strip("\n").splitlines(keepends=True)}

cells = []
cells.append(md(INTRO + "\n\n**Data**: `obo_meta.tsv`, `obo_series.tsv`, `obo_scans.tsv` in this "
                "folder (written by `examples/obo_envelope_probe.rs`). Use the **dropdown** below."))

cells.append(code(r"""
import csv, os, math
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

# Find the folder holding the probe TSVs (this notebook's folder, or a known subfolder).
DATA_DIR = os.getcwd()
for cand in (DATA_DIR, os.path.join(DATA_DIR, "misses"), os.path.join(DATA_DIR, "gains"),
             os.path.join(DATA_DIR, "analysis"), os.path.join(DATA_DIR, "analysis", "misses"),
             os.path.join(DATA_DIR, "analysis", "gains"), "."):
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

series = defaultdict(lambda: defaultdict(list))
with open(os.path.join(DATA_DIR, "obo_series.tsv"), newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        series[row["case"]][row["series"]].append((float(row["mz"]), float(row["intensity"])))

scans = defaultdict(lambda: defaultdict(list))
scan_rt = defaultdict(dict)
with open(os.path.join(DATA_DIR, "obo_scans.tsv"), newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        off = int(row["scan_offset"])
        scans[row["case"]][off].append((float(row["mz"]), float(row["intensity"])))
        scan_rt[row["case"]][off] = float(row["rt"])

print(len(order), "cases loaded")
"""))

cells.append(md("## The cases\n\nList with charge, reference mass, and " + TITLE_TABLE + ":"))

cells.append(code(r"""
print(f"{'idx':>3}  {'case':34} {'z':>2} {'ref_mono':>10} {'reported_k':>10}")
for i, c in enumerate(order):
    m = meta[c]
    print(f"{i:>3}  {c:34} {m['z']:>2} {m['ref_mono']:>10.3f} {m['reported_k']:>+10d}")
"""))

cells.append(md(r"""## A self-contained averagine model

Overlay a template at any shift: a small pure-Python averagine (convolves elemental isotope patterns),
reproducing the shape of the Rust `deconvolution::averagine_*` model. The validation cell overlays it
on the exact Rust templates the probe exported."""))

cells.append(code(r"""
PROTON = 1.007276466
_EL = {
    "C": (12.0,          [0.9893, 0.0107]),
    "H": (1.0078250319,  [0.999885, 0.000115]),
    "N": (14.0030740052, [0.99636, 0.00364]),
    "O": (15.9949146221, [0.99757, 0.00038, 0.00205]),
    "S": (31.97207069,   [0.9499, 0.0075, 0.0425, 0.0, 0.0001]),
}
_AVG = {"C": 4.9384, "H": 7.7583, "O": 1.4773, "N": 1.3577, "S": 0.0417}
_RES_MASS = sum(_AVG[e] * _EL[e][0] for e in _AVG)

def _conv_pow(base, n, maxlen=18):
    dist = np.array([1.0]); b = np.array(base)
    for _ in range(int(n)):
        dist = np.convolve(dist, b)[:maxlen]
    return dist

def averagine_weights(mono_mass, maxlen=16):
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
    C13 = 1.0033548381
    anchor = mono_mass + kshift * C13
    w = averagine_weights(anchor, maxlen)
    mode = int(np.argmax(w)); end = len(w)
    while end > mode + 1 and w[end-1] < min_w:
        end -= 1
    w = w[:end]
    mz = [to_mz(anchor + k * C13, z) for k in range(len(w))]
    return list(zip(mz, w))

print("averagine residue mass ~ %.3f Da/multiplier" % _RES_MASS)
"""))

cells.append(code(r"""
show = order[:4]
fig, axes = plt.subplots(1, len(show), figsize=(4*len(show), 3.2))
for ax, case in zip(np.atleast_1d(axes), show):
    m = meta[case]
    rust = sorted(series[case].get("avg_mono", []))
    py = template(m["ref_mono"], m["z"], 0)
    if rust:
        ax.plot([p[0] for p in rust], [p[1] for p in rust], "o-", label="Rust", alpha=.7)
    ax.plot([p[0] for p in py], [p[1] for p in py], "x--", label="Python", alpha=.7)
    ax.set_title(case, fontsize=7); ax.legend(fontsize=7)
fig.suptitle("Averagine template: Rust (exact) vs Python"); fig.tight_layout()
"""))

cells.append(md(r"""## Plotting functions

- `plot_case(case, observed='observed_composite', extra_shifts=(), legend=True, ax=None)`
- `plot_elution(case)` — every scan in the apex window as small multiples."""))

cells.append(code(r"""
REF_C, PIPE_C, OBS_C = "#2ca02c", "#c2185b", "0.55"

def plot_case(case, observed="observed_composite", extra_shifts=(), legend=True, ax=None):
    m = meta[case]; z, sp, mono_mz, rk = m["z"], m["spacing"], m["mono_mz"], m["reported_k"]
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 5))
    xmin, xmax = mono_mz - 1.6*sp, mono_mz + 9*sp
    obs = [(mz, i) for mz, i in series[case].get(observed, []) if xmin <= mz <= xmax]
    omax = max((i for _, i in obs), default=1.0) or 1.0
    for mz, i in obs:
        ax.vlines(mz, 0, i/omax, color=OBS_C, lw=3, zorder=1)
    ax.vlines([], [], [], color=OBS_C, lw=3, label="observed")
    curves = [(0, REF_C, "o", "FlashLFQ (ref)")]
    if rk != 0:
        curves.append((rk, PIPE_C, "X", f"feature (k={rk:+d})"))
    for k in extra_shifts:
        curves.append((k, None, ".", f"k={k:+d}"))
    for k, color, marker, name in curves:
        pts = [(mz, w) for mz, w in template(m["ref_mono"], z, k) if xmin-sp <= mz <= xmax+sp]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "-", marker=marker, ms=6, lw=1.6,
                alpha=.9, color=color, label=name, zorder=3)
    ax.axvline(mono_mz, color=REF_C, ls="--", lw=1.0, alpha=.8)
    if rk != 0:
        ax.axvline(mono_mz + rk*sp, color=PIPE_C, ls=":", lw=1.3, alpha=.9)
    ax.set_xlim(xmin, xmax); ax.set_ylim(0, 1.22)
    ax.set_title(f"{case}\nz{z}  {m['ref_mono']:.3f} Da  (k={rk:+d})", fontsize=8)
    ax.set_xlabel("m/z", fontsize=8); ax.set_ylabel("rel. int", fontsize=8)
    if legend:
        ax.legend(fontsize=7, loc="upper right")
    ax.grid(True, alpha=.15)
    return ax

def plot_elution(case):
    m = meta[case]; z, sp, mono_mz, rk = m["z"], m["spacing"], m["mono_mz"], m["reported_k"]
    offs = sorted(scans[case]); n = len(offs)
    if n == 0:
        print("no per-scan data for", case); return
    fig, axes = plt.subplots(1, n, figsize=(3.0*n, 3.2), sharey=True)
    xmin, xmax = mono_mz - 1.6*sp, mono_mz + 9*sp
    for ax, off in zip(np.atleast_1d(axes), offs):
        pts = [(mz, i) for mz, i in scans[case][off] if xmin <= mz <= xmax]
        omax = max((i for _, i in pts), default=1.0) or 1.0
        for mz, i in pts:
            ax.vlines(mz, 0, i/omax, color=OBS_C, lw=2.4)
        ax.axvline(mono_mz, color=REF_C, ls="--", lw=1.0)
        if rk != 0:
            ax.axvline(mono_mz + rk*sp, color=PIPE_C, ls=":", lw=1.2)
        ax.set_xlim(xmin, xmax); ax.set_ylim(0, 1.15)
        tag = "APEX" if off == 0 else f"apex{off:+d}"
        ax.set_title(f"{tag}\nRT {scan_rt[case][off]:.3f}", fontsize=8)
        ax.set_xlabel("m/z", fontsize=8)
    np.atleast_1d(axes)[0].set_ylabel("rel. int")
    fig.suptitle(f"{case}  z{z}  — green=FlashLFQ mono", y=1.03, fontsize=10)
    fig.tight_layout()

plot_case(order[0]);
"""))

cells.append(md(r"""## Interactive explorer

Pick any case from the dropdown — composite + apex + full elution. (Needs ipywidgets from the repo
`.venv`; else use `plot_case(order[i])` / `plot_elution(order[i])`.)"""))

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
n = len(order); ncols = 4; nrows = max(1, math.ceil(n / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(4.4*ncols, 2.9*nrows), squeeze=False)
flat = axes.flat
for ax, case in zip(flat, order):
    plot_case(case, "observed_composite", legend=False, ax=ax)
for ax in list(flat)[n:]:
    ax.axis("off")
fig.tight_layout()
"""))

cells.append(md("## Free exploration\n\nOverlay whatever shifts you like, on any case index:"))

cells.append(code(r"""
i = 0                       # <- change the index
plot_case(order[i], "observed_composite", extra_shifts=[-2, -1, +1]);
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

os.makedirs(DATA_LOCAL, exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1)
print(f"wrote {OUT}  ({len(cells)} cells, mode={MODE})")
