#!/usr/bin/env python
"""Parse a DETECT_PROGRESS log from a 95%-coverage long-gradient run and plot the
two curves the stopping-point analysis needs:

  1. wall-time (s)        vs  %TIC explained
  2. seeds considered     vs  %TIC explained

Plus a marginal-return panel (d%TIC / dseed, smoothed) that makes the knee visible,
and a printed table of candidate stopping points from several heuristics.

Usage: plot_progress_cutoff.py [path/to/progress_cov95.log]
"""
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "progress_cov95.log")

# [DETECT_PROGRESS] seeds 150000/482113 (pool_remaining 332113) | accepted 4211 | TIC 8.30e9 (88.4% of ...) | 42.3s
# NB: PowerShell's `2>` redirect wraps native stderr at ~120 cols, so a record's trailing
# "Ns" often spills onto the next physical line. We normalise all whitespace to single spaces
# and finditer, so a wrap (now just a space) is transparent.
LINE = re.compile(
    r"\[DETECT_PROGRESS\]\s+seeds\s+(\d+)/(\d+)\s+\(pool_remaining\s+(\d+)\)\s+\|\s+"
    r"accepted\s+(\d+)\s+\|\s+TIC\s+([\d.eE+-]+)\s+\(([\d.]+)%[^|]*\|\s+([\d.]+)s"
)

def read_text(path):
    """PowerShell '2>' redirect writes UTF-16 LE (BOM); older logs may be UTF-8."""
    with open(path, "rb") as fh:
        raw = fh.read()
    if raw[:2] in (b"\xff\xfe", b"\xfe\xff"):
        return raw.decode("utf-16")
    return raw.decode("utf-8", errors="replace")

blob = re.sub(r"\s+", " ", read_text(LOG))

seeds, pool_total, accepted, tic_abs, tic_pct, t_s = [], [], [], [], [], []
for m in LINE.finditer(blob):
    seeds.append(int(m.group(1)))
    pool_total.append(int(m.group(2)))
    accepted.append(int(m.group(4)))
    tic_abs.append(float(m.group(5)))
    tic_pct.append(float(m.group(6)))
    t_s.append(float(m.group(7)))

if len(seeds) < 3:
    sys.exit(f"only {len(seeds)} progress lines parsed from {LOG} - is the run done / format right?")

seeds = np.array(seeds, float)
accepted = np.array(accepted, float)
tic_pct = np.array(tic_pct, float)
t_s = np.array(t_s, float)
POOL = pool_total[-1]

print(f"parsed {len(seeds)} progress samples from {os.path.basename(LOG)}")
print(f"pool_total seeds = {POOL:,}   final: {tic_pct[-1]:.2f}% TIC, "
      f"{int(accepted[-1]):,} features, {seeds[-1]:,.0f} seeds, {t_s[-1]:.1f}s")

# de-dup any repeated x so gradients are well defined
uniq_s, idx = np.unique(seeds, return_index=True)
uniq_p = tic_pct[idx]
uniq_t = t_s[idx]

def rolling_marginal(x, p, win_frac=0.03):
    """Rolling slope dp/dx over a window that is `win_frac` of the x-range,
    evaluated at each sample. Robust to the fine-grained jitter of raw gradients."""
    span = (x[-1] - x[0]) * win_frac
    out = np.full_like(p, np.nan, dtype=float)
    for i in range(len(x)):
        lo = np.searchsorted(x, x[i] - span / 2)
        hi = np.searchsorted(x, x[i] + span / 2)
        lo = min(lo, i); hi = max(hi, i + 1)
        if hi - lo >= 2 and (x[hi - 1] - x[lo]) > 0:
            out[i] = (p[hi - 1] - p[lo]) / (x[hi - 1] - x[lo])
    return out

marg_seed = rolling_marginal(uniq_s, uniq_p)   # %TIC per seed
marg_time = rolling_marginal(uniq_t, uniq_p)   # %TIC per second

# ---- coverage-target reference points ---------------------------------------
def at_pct(target):
    j = int(np.searchsorted(tic_pct, target))
    j = min(j, len(seeds) - 1)
    return seeds[j], t_s[j], tic_pct[j], accepted[j]

print("\n--- fixed coverage targets (cost to reach) ---")
for tgt in (80, 85, 88, 90, 92, 93):
    s, t, p, a = at_pct(tgt)
    print(f"  {p:5.1f}% TIC : {s:>10,.0f} seeds ({100*s/seeds[-1]:4.1f}% of walk) | "
          f"{t:7.1f}s ({100*t/t_s[-1]:4.1f}% of time) | {int(a):>8,} feats")

# ---- Kneedle (max distance to chord) on each axis ---------------------------
def chord_knee(x, p):
    xn = (x - x.min()) / (x.max() - x.min())
    yn = (p - p.min()) / (p.max() - p.min())
    chord = yn[0] + (yn[-1] - yn[0]) * (xn - xn[0]) / (xn[-1] - xn[0])
    k = int(np.argmax(np.abs(yn - chord)))
    return k

def report(label, k):
    s = uniq_s[k]; j = min(int(np.searchsorted(seeds, s)), len(seeds) - 1)
    print(f"  {label:28s}: {uniq_p[k]:5.1f}% TIC | {s:>10,.0f} seeds "
          f"({100*s/seeds[-1]:4.1f}%) | {uniq_t[k]:7.1f}s ({100*uniq_t[k]/t_s[-1]:4.1f}%) "
          f"| {int(accepted[j]):>8,} feats")

ks_seed = chord_knee(uniq_s, uniq_p)
ks_time = chord_knee(uniq_t, uniq_p)
print("\n--- parameter-free knees (Kneedle / max-distance-to-chord) ---")
report("chord knee  (seeds axis)", ks_seed)
report("chord knee  (time axis)",  ks_time)

# ---- relative-marginal-slope rule (self-calibrating; a la DETECT_KNEE) ------
# early rate = median marginal over the first `early_frac` of the walk; stop the
# first time the rolling marginal falls below `frac` * early rate.
def slope_rule(x, marg, frac, early_frac=0.05):
    early_hi = x[0] + (x[-1] - x[0]) * early_frac
    m0 = np.nanmedian(marg[x <= early_hi])
    thresh = frac * m0
    below = np.where(marg < thresh)[0]
    return int(below[0]) if len(below) else len(x) - 1

print("\n--- relative marginal-slope rule: stop when rate < frac x early rate ---")
print("    (per-SEED marginal)")
for frac in (0.25, 0.10, 0.05, 0.02):
    report(f"  seed-slope frac={frac}", slope_rule(uniq_s, marg_seed, frac))
print("    (per-TIME marginal)")
for frac in (0.25, 0.10, 0.05, 0.02):
    report(f"  time-slope frac={frac}", slope_rule(uniq_t, marg_time, frac))

# ---- plots ------------------------------------------------------------------
knee_time_s = uniq_t[ks_time]
knee_seed_s = uniq_s[ks_seed]

fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

ax = axes[0]
ax.plot(t_s, tic_pct, "-", color="#1f77b4", lw=1.8)
ax.axvline(knee_time_s, color="#888", ls="--", lw=1,
           label=f"time-knee {uniq_p[ks_time]:.1f}% @ {knee_time_s:.0f}s")
ax.set_xlabel("Wall time (s)")
ax.set_ylabel("%TIC explained")
ax.set_title("Time vs %TIC explained")
ax.grid(True, alpha=0.3); ax.legend(fontsize=8, loc="lower right")

ax = axes[1]
ax.plot(seeds, tic_pct, "-", color="#2ca02c", lw=1.8)
ax.axvline(knee_seed_s, color="#888", ls="--", lw=1,
           label=f"seed-knee {uniq_p[ks_seed]:.1f}% @ {knee_seed_s/1e6:.2f}M")
ax.set_xlabel("Seeds considered")
ax.set_ylabel("%TIC explained")
ax.set_title("Seeds vs %TIC explained")
ax.grid(True, alpha=0.3); ax.legend(fontsize=8, loc="lower right")

ax = axes[2]
ax.plot(uniq_s, marg_seed, "-", color="#d62728", lw=1.4)
ax.axvline(knee_seed_s, color="#888", ls="--", lw=1)
ax.set_xlabel("Seeds considered")
ax.set_ylabel("Marginal %TIC per seed (rolling)")
ax.set_title("Marginal return per seed")
ax.set_yscale("log")
ax.grid(True, alpha=0.3, which="both")

fig.tight_layout()
outp = os.path.join(HERE, "progress_cutoff.png")
fig.savefig(outp, dpi=130)
print(f"\nwrote {outp}")
