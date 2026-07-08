#!/usr/bin/env python
"""Seed-rejection-rate curve: how often a *considered* (scored) seed fails to yield a
viable envelope, as detection descends into noise. Hypothesis: the marginal rejection
rate has a sharp knee that makes a good stopping point.

Reads a DETECT_PROGRESS log carrying the tally fields (scored / rej_env / rej_pers),
where for every scored seed:  scored == accepted + rej_env + rej_pers.
  rej_env  = no viable envelope (isotope-count / response gate)
  rej_pers = had an envelope but failed the chromatographic-persistence gate

Plots the MARGINAL rejection rate (per scored seed, rolling window) vs %TIC and vs
seeds considered, and locates the knee.

Usage: plot_reject_rate.py [progress_cov95_reject.log]
"""
import os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "progress_cov95_reject.log")

LINE = re.compile(
    r"\[DETECT_PROGRESS\]\s+seeds\s+(\d+)/(\d+)\s+\(pool_remaining\s+(\d+)\)\s+\|\s+"
    r"accepted\s+(\d+)\s+\|\s+seed_int\s+([\d.]+)\s+\|\s+scored\s+(\d+)\s+\|\s+"
    r"rej_env\s+(\d+)\s+\|\s+rej_pers\s+(\d+)\s+\|\s+TIC\s+([\d.eE+-]+)\s+"
    r"\(([\d.]+)%[^|]*\|\s+([\d.]+)s"
)

def read_text(path):
    with open(path, "rb") as fh:
        raw = fh.read()
    if raw[:2] in (b"\xff\xfe", b"\xfe\xff"):
        return raw.decode("utf-16")
    return raw.decode("utf-8", errors="replace")

blob = re.sub(r"\s+", " ", read_text(LOG))
cols = {k: [] for k in ("seeds", "acc", "sint", "scored", "renv", "rpers", "pct", "t")}
for m in LINE.finditer(blob):
    cols["seeds"].append(int(m.group(1)))
    cols["acc"].append(int(m.group(4)))
    cols["sint"].append(float(m.group(5)))
    cols["scored"].append(int(m.group(6)))
    cols["renv"].append(int(m.group(7)))
    cols["rpers"].append(int(m.group(8)))
    cols["pct"].append(float(m.group(10)))
    cols["t"].append(float(m.group(11)))

n = len(cols["seeds"])
if n < 5:
    sys.exit(f"only {n} samples parsed from {LOG} - has the tally landed / is the run done?")
for k in cols:
    cols[k] = np.array(cols[k], float)

seeds, scored = cols["seeds"], cols["scored"]
renv, rpers, acc = cols["renv"], cols["rpers"], cols["acc"]
pct, t = cols["pct"], cols["t"]
rej = renv + rpers

print(f"parsed {n} samples; final: scored {scored[-1]:,.0f}, accepted {acc[-1]:,.0f}, "
      f"rej_env {renv[-1]:,.0f}, rej_pers {rpers[-1]:,.0f}")
print(f"identity check scored == acc+rej ? {np.allclose(scored, acc+rej)}")
print(f"cumulative envelope-reject rate at end: {100*renv[-1]/scored[-1]:.1f}%  "
      f"total-reject rate: {100*rej[-1]/scored[-1]:.1f}%")

# ---- marginal rate per scored seed, over a rolling window in *scored* ----
def marginal(num, den, win_frac=0.03):
    span = (den[-1] - den[0]) * win_frac
    out = np.full(len(den), np.nan)
    for i in range(len(den)):
        lo = np.searchsorted(den, den[i] - span/2)
        hi = np.searchsorted(den, den[i] + span/2)
        lo = min(lo, i); hi = max(hi, i+1)
        dd = den[hi-1] - den[lo]
        if hi-lo >= 2 and dd > 0:
            out[i] = (num[hi-1] - num[lo]) / dd
    return out

marg_rej   = marginal(rej,  scored)   # fraction of considered seeds rejected (any reason)
marg_env   = marginal(renv, scored)   # fraction rejected for no viable envelope
marg_pers  = marginal(rpers, scored)

# ---- knee of the marginal-reject curve vs cost (time) -----------------------
good = ~np.isnan(marg_rej)
xt, yr = t[good], marg_rej[good]
xn = (xt - xt.min())/(xt.max()-xt.min())
yn = (yr - yr.min())/(yr.max()-yr.min())
chord = yn[0] + (yn[-1]-yn[0])*(xn-xn[0])/(xn[-1]-xn[0])
k = int(np.argmax(np.abs(yn - chord)))
knee_t = xt[k]
kj = min(int(np.searchsorted(t, knee_t)), n-1)
print(f"\nmarginal-reject knee (time axis): {marg_rej[kj]*100:.0f}% reject | "
      f"{pct[kj]:.1f}% TIC | {t[kj]:.0f}s | seed_int {cols['sint'][kj]:,.0f}")

# report marginal reject rate at each coverage crossing
print("\ncov%   marg_reject%   cum_reject%   seed_int")
for c in (85, 88, 90, 92, 95):
    j = min(int(np.searchsorted(pct, c)), n-1)
    mr = marg_rej[j]*100 if not np.isnan(marg_rej[j]) else float('nan')
    print(f"  {c}     {mr:5.1f}        {100*rej[j]/scored[j]:5.1f}       {cols['sint'][j]:,.0f}")

# ---- plots ------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

ax = axes[0]
ax.plot(pct, 100*marg_rej, "-", color="#d62728", lw=1.8, label="any reject")
ax.plot(pct, 100*marg_env, "-", color="#ff7f0e", lw=1.2, label="no-envelope")
ax.plot(pct, 100*marg_pers, "-", color="#9467bd", lw=1.0, label="persistence")
ax.axvline(pct[kj], color="#888", ls="--", lw=1, label=f"knee {pct[kj]:.0f}% TIC")
ax.set_xlabel("%TIC explained"); ax.set_ylabel("Marginal reject rate (%)")
ax.set_title("Seed-reject rate vs coverage"); ax.set_ylim(0, 100)
ax.grid(True, alpha=0.3); ax.legend(fontsize=8, loc="upper left")

ax = axes[1]
ax.plot(seeds, 100*marg_rej, "-", color="#d62728", lw=1.8)
ax.axvline(seeds[kj], color="#888", ls="--", lw=1)
ax.set_xlabel("Seeds considered"); ax.set_ylabel("Marginal reject rate (%)")
ax.set_title("Seed-reject rate vs seeds"); ax.set_ylim(0, 100)
ax.grid(True, alpha=0.3)

ax = axes[2]
ax.plot(t, 100*marg_rej, "-", color="#d62728", lw=1.8)
ax.axvline(t[kj], color="#888", ls="--", lw=1, label=f"knee {t[kj]:.0f}s")
ax.set_xlabel("Wall time (s)"); ax.set_ylabel("Marginal reject rate (%)")
ax.set_title("Seed-reject rate vs cost"); ax.set_ylim(0, 100)
ax.grid(True, alpha=0.3); ax.legend(fontsize=8, loc="lower right")

fig.tight_layout()
outp = os.path.join(HERE, "reject_rate.png")
fig.savefig(outp, dpi=130)
print(f"\nwrote {outp}")
