#!/usr/bin/env python3
"""Plot MS1 spectra with isotope-pattern classification color-coded.

Reads the per-peak TSV exported by the Rust example
`isotope_signal_fraction` (run it with `EXPORT_TSV=... EXPORT_SCANS=...`),
and draws one stem plot per scan: peaks that the detector's isotope-chain
logic marked as isotope-structured are colored by their assigned charge;
peaks with no isotopic partner ("excluded") are drawn in light gray.

Because the classification is exported straight from the Rust `assign_charge`
routine, what you see is exactly what the algorithm decided — not a Python
re-implementation. That is the point: eyeball whether the colored peaks really
sit inside isotope envelopes and whether the gray ones are genuinely noise.

Usage:
  python plot_isotope_classification.py <export.tsv> [options]

Options:
  --scan N          only plot scan_index N (default: all scans in the file)
  --mz-min X        lower m/z bound (zoom)
  --mz-max X        upper m/z bound (zoom)
  --log             log-scale the intensity axis
  --out PATH        save to PATH (.png/.pdf) instead of showing interactively;
                    for multiple scans a "-scanN" suffix is inserted
  --dpi N           output DPI when saving (default 140)
"""
import argparse
import sys
from collections import defaultdict

import matplotlib

# Charge -> color. Gray is reserved for excluded (charge 0). Colors chosen to be
# distinguishable in the default (light) theme and colorblind-friendlier than a rainbow.
CHARGE_COLORS = {
    1: "#8c564b",  # brown
    2: "#1f77b4",  # blue
    3: "#d62728",  # red
    4: "#2ca02c",  # green
    5: "#9467bd",  # purple
    6: "#ff7f0e",  # orange
}
EXCLUDED_COLOR = "#c8c8c8"  # light gray


def read_export(path):
    """Reads the TSV into {scan_index: dict(one_based, rt, peaks=[(mz, inten, charge)])}."""
    scans = {}
    with open(path, "r", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        required = ["scan_index", "one_based_scan", "rt", "mz", "intensity", "charge"]
        missing = [c for c in required if c not in idx]
        if missing:
            sys.exit(f"export file is missing columns: {missing}\nheader was: {header}")
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            si = int(f[idx["scan_index"]])
            rec = scans.get(si)
            if rec is None:
                rec = {
                    "one_based": int(f[idx["one_based_scan"]]),
                    "rt": float(f[idx["rt"]]),
                    "peaks": [],
                }
                scans[si] = rec
            rec["peaks"].append(
                (
                    float(f[idx["mz"]]),
                    float(f[idx["intensity"]]),
                    int(f[idx["charge"]]),
                )
            )
    return scans


def plot_scan(ax, scan_index, rec, mz_min=None, mz_max=None, log=False):
    peaks = rec["peaks"]
    if mz_min is not None:
        peaks = [p for p in peaks if p[0] >= mz_min]
    if mz_max is not None:
        peaks = [p for p in peaks if p[0] <= mz_max]

    # Group by charge so each class is one vlines call with a legend entry.
    by_charge = defaultdict(list)
    for mz, inten, z in peaks:
        by_charge[z].append((mz, inten))

    total_tic = sum(inten for _, inten, _ in peaks)
    iso_tic = sum(inten for _, inten, z in peaks if z != 0)
    iso_pct = 100.0 * iso_tic / total_tic if total_tic > 0 else 0.0

    # Draw excluded first (behind), then each charge on top so envelopes stand out.
    order = [0] + sorted(z for z in by_charge if z != 0)
    for z in order:
        pts = by_charge.get(z)
        if not pts:
            continue
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        if z == 0:
            color, label, lw, zorder = EXCLUDED_COLOR, f"excluded (n={len(pts)})", 0.8, 1
        else:
            color = CHARGE_COLORS.get(z, "#000000")
            label = f"z={z} (n={len(pts)})"
            lw, zorder = 1.4, 3
        ax.vlines(xs, 0, ys, color=color, linewidth=lw, label=label, zorder=zorder)

    if log:
        ax.set_yscale("log")
    ax.set_xlabel("m/z")
    ax.set_ylabel("intensity")
    ax.set_title(
        f"scan_index {scan_index}  (scan #{rec['one_based']}, RT {rec['rt']:.3f} min)   "
        f"isotope-structured: {iso_pct:.1f}% of scan TIC   |   {len(peaks)} peaks shown",
        fontsize=10,
    )
    ax.legend(fontsize=8, loc="upper right")
    ax.margins(x=0.01)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", help="per-peak export TSV from the Rust example")
    ap.add_argument("--scan", type=int, default=None, help="only plot this scan_index")
    ap.add_argument("--mz-min", type=float, default=None)
    ap.add_argument("--mz-max", type=float, default=None)
    ap.add_argument("--log", action="store_true", help="log-scale intensity")
    ap.add_argument("--out", default=None, help="save instead of show (png/pdf)")
    ap.add_argument("--dpi", type=int, default=140)
    args = ap.parse_args()

    # Use a non-interactive backend when saving so it works headless.
    if args.out:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    scans = read_export(args.tsv)
    if not scans:
        sys.exit("no scans found in export file")

    selected = [args.scan] if args.scan is not None else sorted(scans)
    for si in selected:
        if si not in scans:
            print(f"  (scan_index {si} not in export; skipping)", file=sys.stderr)
            continue
        fig, ax = plt.subplots(figsize=(13, 5))
        plot_scan(ax, si, scans[si], args.mz_min, args.mz_max, args.log)
        fig.tight_layout()
        if args.out:
            if len(selected) > 1:
                stem, dot, ext = args.out.rpartition(".")
                out = f"{stem}-scan{si}.{ext}" if dot else f"{args.out}-scan{si}"
            else:
                out = args.out
            fig.savefig(out, dpi=args.dpi)
            print(f"wrote {out}")
            plt.close(fig)
        else:
            print(f"showing scan_index {si} (close window for next)")

    if not args.out:
        plt.show()


if __name__ == "__main__":
    main()
