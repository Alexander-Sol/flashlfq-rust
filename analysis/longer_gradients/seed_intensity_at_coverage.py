#!/usr/bin/env python
"""Map each coverage point to the seed-intensity floor that produced it.

Reads a DETECT_PROGRESS log that includes the `seed_int` field (added to the
detector's progress line). Because seeds are visited intensity-descending, the
seed_int at the sample where %TIC first crosses a target IS the effective
seed-intensity floor at that coverage. Cross-references the recall sweep so the
final table is: coverage -> seed_int floor -> recall -> wall cost.

Usage: seed_intensity_at_coverage.py [progress_cov95_seedint.log]
"""
import os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "progress_cov95_seedint.log")

LINE = re.compile(
    r"\[DETECT_PROGRESS\]\s+seeds\s+(\d+)/(\d+)\s+\(pool_remaining\s+(\d+)\)\s+\|\s+"
    r"accepted\s+(\d+)\s+\|\s+seed_int\s+([\d.]+)\s+\|\s+TIC\s+([\d.eE+-]+)\s+"
    r"\(([\d.]+)%[^|]*\|\s+([\d.]+)s"
)

def read_text(path):
    with open(path, "rb") as fh:
        raw = fh.read()
    if raw[:2] in (b"\xff\xfe", b"\xfe\xff"):
        return raw.decode("utf-16")
    return raw.decode("utf-8", errors="replace")

blob = re.sub(r"\s+", " ", read_text(LOG))
seeds, seed_int, tic_pct, t_s, accepted = [], [], [], [], []
for m in LINE.finditer(blob):
    seeds.append(int(m.group(1)))
    accepted.append(int(m.group(4)))
    seed_int.append(float(m.group(5)))
    tic_pct.append(float(m.group(7)))
    t_s.append(float(m.group(8)))

if len(seeds) < 3:
    sys.exit(f"only {len(seeds)} samples parsed from {LOG} - has seed_int landed / is the run done?")

seeds = np.array(seeds, float); seed_int = np.array(seed_int, float)
tic_pct = np.array(tic_pct, float); t_s = np.array(t_s, float)
print(f"parsed {len(seeds)} samples; final {tic_pct[-1]:.1f}% TIC, "
      f"seed_int floor {seed_int[-1]:.0f}, {t_s[-1]:.1f}s")

# recall (from the recall sweep; 0.95 measured in-session). recall @ RT +-0.3
RECALL = {85: 69.4, 88: 76.4, 90: 81.2, 92: 85.8, 95: 90.7}
WALL   = {85: 62.6, 88: 69.2, 90: 81.0, 92: 106.4, 95: 336.8}

def floor_at(target):
    j = int(np.searchsorted(tic_pct, target))
    j = min(j, len(seeds) - 1)
    return seed_int[j], seeds[j], t_s[j]

print("\ncov%   seed_int floor   seeds@stop     recall%   wall(s)")
targets = [85, 88, 90, 92, 95]
rows = []
for c in targets:
    si, s, t = floor_at(c)
    rows.append((c, si, s))
    print(f"  {c}    {si:>12,.0f}   {s:>11,.0f}     {RECALL[c]:5.1f}    {WALL[c]:6.1f}")

# ---- plots: seed_int vs %TIC, and recall vs seed_int floor ------------------
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

ax = axes[0]
ax.plot(tic_pct, seed_int, "-", color="#1f77b4", lw=1.8)
ax.set_yscale("log")
for c, si, s in rows:
    ax.plot(c, si, "o", color="#d62728", ms=6)
    ax.annotate(f"{c}% -> {si:,.0f}", (c, si), textcoords="offset points",
                xytext=(-4, 8), ha="right", fontsize=8)
ax.set_xlabel("%TIC explained")
ax.set_ylabel("Seed-intensity floor (log)")
ax.set_title("Seed-intensity floor vs coverage")
ax.grid(True, alpha=0.3, which="both")

ax = axes[1]
si_pts = [r[1] for r in rows]
rc_pts = [RECALL[r[0]] for r in rows]
ax.plot(si_pts, rc_pts, "s-", color="#2ca02c", lw=2, ms=6)
ax.set_xscale("log")
for (c, si, s), rc in zip(rows, rc_pts):
    ax.annotate(f"{c}%", (si, rc), textcoords="offset points", xytext=(6, -4),
                fontsize=8, color="#555")
ax.set_xlabel("Seed-intensity floor (log)")
ax.set_ylabel("Recall vs MSMS (%)")
ax.set_title("Recall vs seed-intensity floor")
ax.grid(True, alpha=0.3, which="both")
ax.invert_xaxis()  # lower floor = deeper = to the right

fig.tight_layout()
outp = os.path.join(HERE, "seed_intensity_at_coverage.png")
fig.savefig(outp, dpi=130)
print(f"\nwrote {outp}")
