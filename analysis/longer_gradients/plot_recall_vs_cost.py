#!/usr/bin/env python
"""Recall-vs-cost curve for the long-gradient coverage sweep.

Reads recall_sweep_summary.txt (coverage, wall_s, feats, recall_pct) produced by
recall_sweep.ps1, splices in the 0.95 point measured in-session, and plots:

  1. recall vs coverage target        (does recall plateau?)
  2. recall vs total wall-time (cost)  (the real question: recall per second)

The %TIC knee (~87%, from the DETECT_PROGRESS analysis) is marked so we can see
whether recall is still climbing where the TIC curve says "stop".

Usage: plot_recall_vs_cost.py [summary.txt]
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
SUMM = sys.argv[1] if len(sys.argv) > 1 else \
    r"C:\Users\Alex\.claude\jobs\75877171\tmp\recall_sweep_summary.txt"

# in-session 0.95 point (long_cov95.tsv): total wall 336.8s, 416,787 resolved feats, 90.7% recall
REUSE = {0.95: (336.8, 416787, 90.7)}

rows = {}
with open(SUMM, encoding="utf-8-sig") as fh:
    next(fh, None)  # header
    for ln in fh:
        p = ln.strip().split("\t")
        if len(p) < 4:
            continue
        cov = float(p[0])
        rows[cov] = (float(p[1]), int(p[2]), float(p[3]))
rows.update(REUSE)

cov = sorted(rows)
wall = [rows[c][0] for c in cov]
feats = [rows[c][1] for c in cov]
recall = [rows[c][2] for c in cov]

print("coverage  wall_s   feats     recall%")
for c in cov:
    w, f, r = rows[c]
    print(f"  {c:.2f}   {w:7.1f}  {f:>8,}   {r:5.1f}")

# marginal recall per second between successive points (recall pp gained per extra second)
print("\nsegment    d_recall(pp)  d_wall(s)   pp/min")
for i in range(1, len(cov)):
    dr = recall[i] - recall[i-1]
    dw = wall[i] - wall[i-1]
    ppmin = dr / (dw/60) if dw > 0 else float("nan")
    print(f"  {cov[i-1]:.2f}->{cov[i]:.2f}   {dr:6.2f}      {dw:7.1f}   {ppmin:6.2f}")

TIC_KNEE = 87.4  # % TIC, time-axis Kneedle from progress analysis

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

ax = axes[0]
ax.plot([100*c for c in cov], recall, "o-", color="#1f77b4", lw=2, ms=6)
for c, r in zip(cov, recall):
    ax.annotate(f"{r:.1f}", (100*c, r), textcoords="offset points", xytext=(0, 7),
                ha="center", fontsize=8)
ax.axvline(TIC_KNEE, color="#d62728", ls="--", lw=1.2, label=f"%TIC knee ~{TIC_KNEE:.0f}%")
ax.set_xlabel("Coverage target (% of ΣTIC)")
ax.set_ylabel("Recall vs MSMS ground truth (%)")
ax.set_title("Recall vs coverage target (RT +-0.3 min)")
ax.grid(True, alpha=0.3); ax.legend(fontsize=8, loc="lower right")

ax = axes[1]
ax.plot(wall, recall, "s-", color="#2ca02c", lw=2, ms=6)
for c, w, r in zip(cov, wall, recall):
    ax.annotate(f"{100*c:.0f}%", (w, r), textcoords="offset points", xytext=(6, -10),
                ha="left", fontsize=8, color="#555")
ax.set_xlabel("Total wall time (s)")
ax.set_ylabel("Recall vs MSMS ground truth (%)")
ax.set_title("Recall vs cost (labels = coverage target)")
ax.grid(True, alpha=0.3)

fig.tight_layout()
outp = os.path.join(HERE, "recall_vs_cost.png")
fig.savefig(outp, dpi=130)
print(f"\nwrote {outp}")
