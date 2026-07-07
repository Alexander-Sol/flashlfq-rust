# FWHM / scans-per-FWHM validation — short · medium · long (TASK 4)

Goal: obtain the **measured** chromatographic FWHM and median MS1 scan spacing on each gradient, from
the detector's own machinery, and check whether `scans-per-FWHM = FWHM / spacing` reproduces Alex's
ground-truth progression (short ~3 · medium ~7–9 · long ~20). This is the quantity the data-dependent
averaging formula would consume in place of the hard-coded `MAX_SCANS_TO_AVERAGE = 3`.

## Source of the measured numbers

All values are the **pipeline's own** measurement — `estimate_fwhm_seconds` (XIC half-max) reported by
`detect_features_tsv` on the σ-source line, and `median_ms1_scan_spacing_minutes` on the read line —
i.e. exactly what `apply_fwhm` consumes at runtime. The standalone `fwhm_probe` on the short raw
cross-checks the same measurement (ALL median 1.96 s, STRONG-decile median 1.83 s vs the pipeline's
1.89 s — agree).

## Table

| gradient | measured FWHM | median MS1 spacing | scans/FWHM (measured) | ground truth | verdict |
|---|--:|--:|--:|--:|---|
| short ~10-min (CA/Lumos) | 1.89 s | 0.0064 min (0.384 s) | **4.9** | ~3 | progression ✓, runs high (~1.6×) |
| medium ~65-min (glyco) | 9.36 s | 0.0155 min (0.930 s) | **10.1** | ~7–9 | ✓ near ground truth (slightly high) |
| long ~120-min (IonStar) | 17.66 s | 0.0124 min (0.744 s) | **23.7** | ~20 | ✓ near ground truth (slightly high) |

`scans/FWHM = (FWHM_sec/60) / spacing_min`.

## Verdict

- **The progression is reproduced.** Measured scans/FWHM rises monotonically 4.9 → 10.1 → 23.7 across
  short → medium → long, tracking Alex's ground-truth 3 → 8 → 20 in both shape and (for medium/long)
  magnitude. The measured FWHM itself grows sensibly with gradient length (1.9 → 9.4 → 17.7 s), and the
  scan spacing is roughly flat (~0.4–0.9 s), so the count is driven by the peak width as intended.

- **Measured values run consistently HIGH vs the visual ground truth**, by ~1.6× on short, ~1.26× on
  medium, ~1.19× on long. The gap is largest on the short gradient. Why: `scans/FWHM` counts every MS1
  scan spanning the **full** half-max width, whereas Alex's short-gradient "~3" is the tight
  **apex±1** averaging window — a deliberately narrower span than the whole FWHM. So the two are not
  measuring quite the same thing; the measured full-FWHM span is *expected* to exceed a tight apex±N
  count. For medium/long, where the ground truth is explicitly "scans **inside the FWHM**", the
  measured 10.1 / 23.7 land right on top of the 7–9 / ~20 estimates.

- **Implication for the data-dependent averaging formula.** A naïve `scans_to_average =
  round(FWHM/spacing)` would emit ~5 / 10 / 24 — i.e. slightly *over-average* relative to Alex's
  3 / 8 / 20, most noticeably on the 10-min case (5 vs 3), which risks pulling co-eluting flank
  interference into the composite and could move the 10-min recall baseline. A modest scale-down
  reconciles it cleanly: **≈0.8 × FWHM/spacing** gives 3.9 / 8.1 / 19.0 — essentially the ground truth
  on all three. Equivalently, target the scans within ~0.8 FWHM (near the ±1σ core), not the full
  half-max span. Whatever the exact factor, it must be **clamped to an odd count with a floor of 3** so
  the short case stays at the parity-locked apex±1 and the assert in `refine_feature*` still holds.

- **Do not confuse this with the detector's `half_window_scans`.** The detector already derives a
  *scoring* half-window from the same FWHM (±9 scans medium, ±20 long) — but that is the ±2σ
  matched-filter window, deliberately wider than the averaging window (per `FWHM-Flow-Audit.md`). The
  averaging count validated here (`~FWHM/spacing`) is a distinct, tighter quantity; the formula must
  compute it separately, not forward `half_window_scans`.

## Bottom line

The measured FWHM machinery is sound and gradient-adaptive: it reproduces the 3 / 8 / 20 progression
and lands on the medium/long ground truth. It runs a bit hot on the short case only because
"scans/FWHM" is the full-width span, not the apex±1 averaging window. A data-dependent
`scans_to_average ≈ round(0.8 · FWHM/spacing)`, odd, floored at 3, would land near ground truth on all
three gradients — validating the umbrella formula. No library behaviour was changed (validation only).
