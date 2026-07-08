#!/usr/bin/env python
"""Plot recall vs coverage target and total runtime vs coverage target for the
long-gradient sweep. Data is hard-coded from the sweep results table so the plot
is reproducible without re-running the (expensive) detector.

Usage: plot_coverage_sweep.py   (writes two PNGs into this directory)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

HERE = os.path.dirname(os.path.abspath(__file__))

# coverage, achieved_tic, n_features, recall_pct, charge_of_all_pct, detect_s, total_s
ROWS = [
    (0.80, 80.0,  34259,  57.7, 55.7,   20.21,   47.82),
    (0.90, 90.0,  96058,  81.2, 76.2,   47.92,   78.23),
    (0.95, 95.0,  416787, 90.7, 79.8,  276.01,  319.37),
    (0.98, 96.5,  929797, 91.7, 77.5, 1396.41, 1464.41),  # target 0.98 unreached (seed floor)
]

cov      = [r[0] for r in ROWS]
recall   = [r[3] for r in ROWS]
total_s  = [r[6] for r in ROWS]
detect_s = [r[5] for r in ROWS]

# --- Panel 1: recall vs coverage ---
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.plot(cov, recall, "o-", color="#1f77b4", lw=2, ms=6)
for x, y in zip(cov, recall):
    ax.annotate(f"{y:.1f}", (x, y), textcoords="offset points", xytext=(0, 7),
                ha="center", fontsize=8)
ax.set_xlabel("Coverage target (fraction of ΣTIC to explain)")
ax.set_ylabel("Recall vs MSMS ground truth (%)")
ax.set_title("Long-gradient recall vs coverage target (±0.3 min)")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "long_recall_vs_coverage.png"), dpi=130)

# --- Panel 2: runtime vs coverage (dual: detect + total) ---
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.plot(cov, total_s, "s-", color="#d62728", lw=2, ms=6, label="total")
ax.plot(cov, detect_s, "^--", color="#ff7f0e", lw=1.6, ms=6, label="detect")
ax.set_xlabel("Coverage target (fraction of ΣTIC to explain)")
ax.set_ylabel("Runtime (s)")
ax.set_title("Long-gradient runtime vs coverage target")
ax.grid(True, alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(HERE, "long_runtime_vs_coverage.png"), dpi=130)
print("wrote long_recall_vs_coverage.png and long_runtime_vs_coverage.png")
