#!/usr/bin/env python3
"""Interactive (Plotly) MS1 spectrum viewer with isotope-pattern classification.

Reads the per-peak TSV exported by the Rust example `isotope_signal_fraction`
(run it with `EXPORT_TSV=... EXPORT_SCANS=...`) and builds a self-contained HTML
plot: peaks the detector's isotope-chain logic marked isotope-structured are
colored by assigned charge; peaks with no isotopic partner ("excluded") are gray.

Unlike the static matplotlib version, this is fully interactive — box-zoom (drag),
pan, wheel-zoom, hover readouts (m/z, intensity, charge), and legend toggles.
If the TSV holds several scans, a dropdown switches between them.

Because the classification comes straight from the Rust `assign_charge` routine,
what you see is exactly what the algorithm decided — not a Python re-derivation.

Usage:
  python plot_isotope_classification_interactive.py <export.tsv> [options]

Options:
  --out PATH     write the HTML here (default: alongside the TSV, same stem .html)
  --no-open      write the HTML but do not open a browser
  --log          start with a log-scaled intensity axis
  --mz-min X     initial lower m/z bound (you can still zoom out)
  --mz-max X     initial upper m/z bound
"""
import argparse
import os
import sys
from collections import defaultdict

import plotly.graph_objects as go

# Charge -> color. Gray is reserved for excluded (charge 0).
CHARGE_COLORS = {
    1: "#8c564b",  # brown
    2: "#1f77b4",  # blue
    3: "#d62728",  # red
    4: "#2ca02c",  # green
    5: "#9467bd",  # purple
    6: "#ff7f0e",  # orange
}
EXCLUDED_COLOR = "#b0b0b0"


def read_export(path):
    """Reads the TSV into {scan_index: {one_based, rt, peaks:[(mz, inten, charge)]}}."""
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
                rec = {"one_based": int(f[idx["one_based_scan"]]), "rt": float(f[idx["rt"]]), "peaks": []}
                scans[si] = rec
            rec["peaks"].append(
                (float(f[idx["mz"]]), float(f[idx["intensity"]]), int(f[idx["charge"]]))
            )
    return scans


def scan_title(scan_index, rec):
    peaks = rec["peaks"]
    total = sum(p[1] for p in peaks)
    iso = sum(p[1] for p in peaks if p[2] != 0)
    pct = 100.0 * iso / total if total > 0 else 0.0
    return (
        f"scan_index {scan_index}  (scan #{rec['one_based']}, RT {rec['rt']:.3f} min)   "
        f"isotope-structured: {pct:.1f}% of scan TIC   |   {len(peaks)} peaks"
    )


def add_scan_traces(fig, scan_index, rec, visible):
    """Adds one stem trace + one hover-marker trace per charge class for a scan.
    Returns the number of traces added (all sharing `visible`)."""
    by_charge = defaultdict(list)
    for mz, inten, z in rec["peaks"]:
        by_charge[z].append((mz, inten))

    order = [0] + sorted(z for z in by_charge if z != 0)  # excluded first (drawn behind)
    n_added = 0
    for z in order:
        pts = by_charge.get(z)
        if not pts:
            continue
        color = EXCLUDED_COLOR if z == 0 else CHARGE_COLORS.get(z, "#000000")
        label = f"excluded (n={len(pts)})" if z == 0 else f"z={z} (n={len(pts)})"
        group = f"z{z}"

        # Vertical stems: (mz,0)->(mz,inten)->None per peak, one trace for the whole class.
        xs, ys = [], []
        for mz, inten in pts:
            xs += [mz, mz, None]
            ys += [0, inten, None]
        fig.add_trace(
            go.Scattergl(
                x=xs, y=ys, mode="lines",
                line=dict(color=color, width=1.2),
                name=label, legendgroup=group, visible=visible,
                hoverinfo="skip",
            )
        )
        # Markers at peak tops carry the hover readout.
        mzs = [p[0] for p in pts]
        ins = [p[1] for p in pts]
        zlabel = "excluded" if z == 0 else f"z={z}"
        fig.add_trace(
            go.Scattergl(
                x=mzs, y=ins, mode="markers",
                marker=dict(color=color, size=4),
                name=label, legendgroup=group, visible=visible, showlegend=False,
                customdata=[zlabel] * len(mzs),
                hovertemplate="m/z %{x:.5f}<br>intensity %{y:.3e}<br>%{customdata}<extra></extra>",
            )
        )
        n_added += 2
    return n_added


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("tsv", help="per-peak export TSV from the Rust example")
    ap.add_argument("--out", default=None, help="output HTML path")
    ap.add_argument("--no-open", action="store_true", help="do not open a browser")
    ap.add_argument("--log", action="store_true", help="start with log-scaled intensity")
    ap.add_argument("--mz-min", type=float, default=None)
    ap.add_argument("--mz-max", type=float, default=None)
    args = ap.parse_args()

    scans = read_export(args.tsv)
    if not scans:
        sys.exit("no scans found in export file")
    scan_ids = sorted(scans)

    fig = go.Figure()
    # Track which scan each trace belongs to, so the dropdown can toggle visibility.
    trace_scan = []
    for pos, si in enumerate(scan_ids):
        visible = pos == 0  # only the first scan visible initially
        n = add_scan_traces(fig, si, scans[si], visible)
        trace_scan += [si] * n

    # Scan-selector dropdown (only if more than one scan).
    if len(scan_ids) > 1:
        buttons = []
        for si in scan_ids:
            vis = [ts == si for ts in trace_scan]
            buttons.append(
                dict(
                    label=f"scan {si}", method="update",
                    args=[{"visible": vis}, {"title.text": scan_title(si, scans[si])}],
                )
            )
        fig.update_layout(
            updatemenus=[dict(buttons=buttons, direction="down", x=1.0, xanchor="right",
                              y=1.14, yanchor="top", showactive=True)]
        )

    fig.update_layout(
        title=dict(text=scan_title(scan_ids[0], scans[scan_ids[0]])),
        xaxis_title="m/z",
        yaxis_title="intensity",
        template="plotly_white",
        hovermode="closest",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
        margin=dict(t=90),
    )
    fig.update_xaxes(showspikes=True, spikemode="across", spikesnap="cursor", spikethickness=1)
    if args.log:
        fig.update_yaxes(type="log")
    if args.mz_min is not None or args.mz_max is not None:
        lo = args.mz_min if args.mz_min is not None else min(
            p[0] for s in scans.values() for p in s["peaks"]
        )
        hi = args.mz_max if args.mz_max is not None else max(
            p[0] for s in scans.values() for p in s["peaks"]
        )
        fig.update_xaxes(range=[lo, hi])

    out = args.out or (os.path.splitext(args.tsv)[0] + ".html")
    # config: enable scroll-to-zoom; keep the full modebar (box zoom, pan, autoscale, reset).
    fig.write_html(out, auto_open=not args.no_open, config={"scrollZoom": True, "displaylogo": False})
    print(f"wrote {out}" + ("" if args.no_open else "  (opening in browser)"))


if __name__ == "__main__":
    main()
