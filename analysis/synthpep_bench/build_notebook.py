#!/usr/bin/env python
"""Emit low_quality_features.ipynb (nbformat v4). Interactive visualizer for our
low-quality refined features: explore decon/intensity, then pull the raw MS1
envelope + XIC behind any selected feature from the mzML."""
import json, os

BENCH = r"F:\flashlfq-rust\analysis\synthpep_bench"

CELLS = []
def md(s):   CELLS.append(("markdown", s))
def code(s): CELLS.append(("code", s))

md(r"""# Low-quality feature visualizer -- PXD001091

Explore the ~7.5x feature tail our detector emits and *look at the raw MS1 signal*
behind individual low-quality features.

Data source: `nb_data/<base>.parquet` (one row per **refined per-charge feature**),
built by `prep_notebook_data.py`. Each feature carries its detector score, refine
`decon_score` (envelope-fit cosine), intensity, isotope/peak counts, and a
`matched` flag = whether it lands within +-20 ppm / +-0.5 min of a FragPipe PSM peak.

**How to use:** run all cells. Cell 5 is an interactive plotly scatter -- hover any
point to read its `feat_idx`, then call `show_feature(feat_idx)` (cell 8) to see its
raw XIC + isotope envelope. Cell 9 wires up a click/selector if `ipywidgets` is present.""")

code(r"""# --- dependencies (installs only what's missing; you run this, so you consent) ---
import sys, subprocess
def ensure(pkg, imp=None):
    try:
        __import__(imp or pkg)
    except ImportError:
        print("installing", pkg, "...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", pkg])
ensure("pyarrow")            # parquet engine for pd.read_parquet (this kernel may lack it)
ensure("pyteomics")          # pure-python mzML reader (raw signal)
ensure("psims")              # required by pyteomics to parse mzML CV params
# pandas / numpy / plotly / matplotlib / pyarrow already present
try:
    import ipywidgets; HAS_WIDGETS = True
except ImportError:
    HAS_WIDGETS = False
print("ipywidgets available:", HAS_WIDGETS)""")

code(r"""import os, numpy as np, pandas as pd
import plotly.graph_objects as go
import matplotlib.pyplot as plt
try:
    from IPython.display import display
except Exception:
    display = print

BENCH = r"{BENCH}"
NB = os.path.join(BENCH, "nb_data")
manifest = pd.read_csv(os.path.join(NB, "manifest.tsv"), sep="\t")
display(manifest)

# >>> pick any base from the manifest <<<
BASE = "130124_dilA_10_01"
MZML = manifest.set_index("base").loc[BASE, "mzml"]
df = pd.read_parquet(os.path.join(NB, BASE + ".parquet"))
print(f"{{BASE}}: {{len(df):,}} features | {{(~df.matched).sum():,}} unmatched | mzML: {{MZML}}")
df.head()""".replace("{BENCH}", BENCH.replace("\\", "\\\\")))

md(r"""## 1. Score distributions

`matched` features (correspond to a PSM peak) vs `unmatched` (the tail). Note most
features are unmatched -- a mix of genuine junk and real-but-unsequenced species.""")

code(r"""fig, ax = plt.subplots(1, 3, figsize=(14, 3.6))
for a, col, title in [(ax[0],"decon_score","Decon Score (envelope-fit cosine)"),
                      (ax[1],"log10_intensity","log10 Summed Intensity"),
                      (ax[2],"num_isotopes","Num Isotopes")]:
    a.hist(df.loc[df.matched, col], bins=60, alpha=.7, label="matched", color="#2e8b57", density=True)
    a.hist(df.loc[~df.matched, col], bins=60, alpha=.5, label="unmatched", color="#d98c00", density=True)
    a.set_title(title); a.set_yscale("log"); a.legend()
plt.tight_layout(); plt.show()
df.groupby("matched")[["decon_score","log10_intensity","num_isotopes","num_peaks","num_candidate_masses"]].describe().T""")

md(r"""## 2. Interactive scatter: decon vs intensity

Hover any point for its `feat_idx`. Green = matched, orange = unmatched. (Unmatched
are downsampled to 50k for responsiveness; all matched are shown.)""")

code(r"""n_un = (~df.matched).sum()
samp = pd.concat([df[df.matched],
                  df[~df.matched].sample(min(50000, n_un), random_state=0)])
fig = go.Figure()
for lab, sub, color in [("unmatched", samp[~samp.matched], "#d98c00"),
                        ("matched",   samp[samp.matched],  "#2e8b57")]:
    fig.add_trace(go.Scattergl(
        x=sub.decon_score, y=sub.log10_intensity, mode="markers", name=lab,
        marker=dict(size=3, color=color, opacity=0.45),
        customdata=np.stack([sub.feat_idx, sub.charge, sub.num_isotopes,
                             sub.num_peaks, sub.mass, sub.apex_rt], axis=1),
        hovertemplate=("feat_idx=%{customdata[0]}<br>decon=%{x:.3f}  logI=%{y:.2f}"
                       "<br>z=%{customdata[1]}  iso=%{customdata[2]}  peaks=%{customdata[3]}"
                       "<br>mass=%{customdata[4]:.4f}  rt=%{customdata[5]:.2f} min<extra></extra>")))
fig.update_layout(height=600, xaxis_title="Decon Score", yaxis_title="log10 Summed Intensity",
                  title=f"{BASE}: decon vs intensity  (hover -> feat_idx)", legend=dict(itemsizing="constant"))
fig.show()""")

md("## 3. Worst offenders (unmatched, lowest score)")

code(r"""cols = ["feat_idx","mass","charge","apex_rt","log10_intensity","decon_score",
        "detector_score","num_isotopes","num_peaks","num_candidate_masses","nearest_gt_ppm","nearest_gt_drt"]
print("Lowest decon_score (unmatched):");  display(df[~df.matched].nsmallest(15,"decon_score")[cols])
print("Lowest intensity (unmatched):");    display(df[~df.matched].nsmallest(15,"log10_intensity")[cols])
print("Single-isotope features (iso==1):", int((df.num_isotopes==1).sum()),
      "of which unmatched:", int(((df.num_isotopes==1)&(~df.matched)).sum()))""")

md(r"""## 4. Raw-signal viewer

`show_feature(feat_idx)` reads the mzML MS1 scans (cached on first call) and plots:
* **left** -- XIC of the mono m/z (+-ppm) across RT, with the feature's RT window shaded;
* **right** -- the MS1 spectrum nearest the apex, with expected isotope positions (green).""")

code(r"""from pyteomics import mzml
PROTON, ISO = 1.007276466812, 1.0033548378
_MS1 = {}
def load_ms1(path):
    if path in _MS1: return _MS1[path]
    scans = []
    with mzml.read(path) as rdr:
        for s in rdr:
            if s.get("ms level") != 1: continue
            rt = float(s["scanList"]["scan"][0]["scan start time"])   # minutes
            scans.append((rt, np.asarray(s["m/z array"]), np.asarray(s["intensity array"])))
    scans.sort(key=lambda x: x[0])
    _MS1[path] = (np.array([s[0] for s in scans]), scans)
    print(f"cached {len(scans)} MS1 scans from {os.path.basename(path)}")
    return _MS1[path]

def show_feature(idx, ppm=15, rt_win=1.5):
    r = df.loc[df.feat_idx == idx].iloc[0]
    mono_mz, z, niso = float(r.mono_mz), int(r.charge), int(r.num_isotopes)
    rts, scans = load_ms1(MZML)
    tol = mono_mz * ppm / 1e6
    xr, xi = [], []
    for rt, mz, it in scans:
        if abs(rt - r.apex_rt) > rt_win: continue
        lo, hi = np.searchsorted(mz, mono_mz - tol), np.searchsorted(mz, mono_mz + tol)
        xr.append(rt); xi.append(it[lo:hi].sum() if hi > lo else 0.0)
    j = int(np.argmin(np.abs(rts - r.apex_rt))); srt, mz, it = scans[j]
    fig, ax = plt.subplots(1, 2, figsize=(13, 4))
    ax[0].plot(xr, xi, marker="."); ax[0].axvline(r.apex_rt, color="k", ls=":")
    ax[0].axvspan(r.rt_start, r.rt_end, color="orange", alpha=0.15)
    ax[0].set_title(f"XIC  mono m/z {mono_mz:.4f} +-{ppm}ppm"); ax[0].set_xlabel("RT (min)"); ax[0].set_ylabel("intensity")
    m0, m1 = mono_mz - 1.0, mono_mz + (niso + 3) * ISO / z
    sel = (mz >= m0) & (mz <= m1)
    ax[1].vlines(mz[sel], 0, it[sel], color="#555")
    for k in range(niso + 3):
        ax[1].axvline(mono_mz + k * ISO / z, color="green", ls=":", alpha=0.5)
    ax[1].set_title(f"MS1 envelope @ {srt:.2f} min  (z={z})"); ax[1].set_xlabel("m/z"); ax[1].set_ylabel("intensity")
    fig.suptitle(f"feat_idx={idx}  decon={r.decon_score:.3f}  logI={r.log10_intensity:.2f}  "
                 f"iso={niso}  peaks={int(r.num_peaks)}  matched={bool(r.matched)}", fontsize=11)
    plt.tight_layout(); plt.show()""")

md("## 5. Examples: a junk feature vs a good one")

code(r"""lowq = int(df[~df.matched].nsmallest(1, "decon_score").feat_idx.iloc[0])
hiq  = int(df[df.matched].nlargest(1, "log10_intensity").feat_idx.iloc[0])
print("low-quality (unmatched, min decon):", lowq); show_feature(lowq)
print("high-quality (matched, max intensity):", hiq); show_feature(hiq)""")

md("## 6. Optional: pick a feature interactively (needs ipywidgets)")

code(r"""if HAS_WIDGETS:
    import ipywidgets as w
    sel = w.IntText(value=lowq, description="feat_idx")
    btn = w.Button(description="show feature", button_style="info")
    out = w.Output()
    def _go(_):
        out.clear_output()
        with out: show_feature(int(sel.value))
    btn.on_click(_go); display(w.HBox([sel, btn]), out); _go(None)
else:
    print("ipywidgets not installed -> just call show_feature(<feat_idx>) directly.")""")

md(r"""## Notes / caveats
* **`matched` is charge-blind** (mass+RT), same as the recall scorer. `unmatched`
  mixes genuine noise/artifacts with real-but-unsequenced species -- don't read every
  orange point as an error.
* Analysis is at the **refined per-charge** level (where a refiner reject threshold
  would apply); the headline 19.9M "features" is the charge-merged resolved view.
* From the full-dataset sweep: rejecting `decon_score < 0.66` drops ~83% of features
  for a ~0.5 pp recall cost (see `score_recall.png`). Use this notebook to eyeball
  what those low-decon features actually look like before choosing a threshold.
* To switch files, change `BASE` in cell 3 (prep more with
  `prep_notebook_data.py <base> ...`).""")

nb = {"cells": [], "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                                "language_info": {"name": "python"}},
      "nbformat": 4, "nbformat_minor": 5}
for ctype, src in CELLS:
    lines = src.split("\n")
    src_lines = [l + "\n" for l in lines[:-1]] + [lines[-1]]
    cell = {"cell_type": ctype, "metadata": {}, "source": src_lines}
    if ctype == "code":
        cell["outputs"] = []; cell["execution_count"] = None
    nb["cells"].append(cell)

out = os.path.join(BENCH, "low_quality_features.ipynb")
with open(out, "w", encoding="utf-8") as fh:
    json.dump(nb, fh, indent=1)
print("wrote", out, "with", len(CELLS), "cells")
