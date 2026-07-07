# FWHM / scans-to-average flow audit

Scope: confirm the measured chromatographic FWHM (and the scans-to-average / peak-width it
implies) is fed **consistently** into every downstream step that assumes a peak width, with no
step left reading a stale hard-coded value. Feeds the umbrella goal "replace hard-coded
`MAX_SCANS_TO_AVERAGE = 3` with a data-dependent scans-to-average".

Ground-truth scans-per-FWHM by gradient (from `TODO.md`): short ~10-min ≈ 3 · medium ~65-min ≈ 7–9 ·
long ~120-min ≈ 20.

## TL;DR key finding

FWHM **is measured at runtime** and **flows correctly into the DETECTOR** (RT σ, the matched-filter
scoring window, and the legacy scan half-window are all derived from it). It **does NOT flow into the
REFINEMENT averaging window** — that is still the hard-coded `MAX_SCANS_TO_AVERAGE = 3` (apex ± 1),
independent of the measured FWHM — nor into the example pipeline's neighbour / consensus / dedup RT
tolerances (hard-coded 0.05 / 0.1 min, tuned to the 10-min case). The detector adapts to gradient
length today; the refine + resolve stages do not.

No purely behaviour-preserving wiring fix exists: threading the measured width into refinement/resolve
necessarily changes the window on every gradient (including the 10-min recall baseline), and the
data-dependent scans-to-average formula is not yet implemented. Everything actionable is therefore
**documented as a prioritized plan below, not applied** (per the "don't change 10-min numeric behavior"
constraint). No source files changed.

## Where FWHM is MEASURED at runtime

- `estimate_fwhm_seconds(engine, ppm)` — `rust/flashlfq-core/src/trace_kernel.rs:258`. Traces the
  tallest clean XICs with the production `get_all_xics`, measures each peak's FWHM by
  linear-interpolated half-max crossings (`xic_fwhm_minutes`, `trace_kernel.rs:307`), returns the
  **median seconds**. Returns `None` if fewer than `FWHM_PROBE_MIN_SAMPLE` (12) clean XICs are found.
- `median_ms1_scan_spacing_minutes(scan_info)` — `trace_kernel.rs:471`. The MS1 scan spacing, the
  other half of any "scans per peak" computation.
- The standalone `examples/fwhm_probe.rs` reports the same measurement as a distribution (p10..p90)
  for eyeballing per gradient — this is the validation tool, not part of the pipeline.

## Where the measurement is CONSUMED — `apply_fwhm`

`TraceKernelParameters::apply_fwhm(scan_info, fwhm_seconds)` — `trace_kernel.rs:223`. Single choke
point that turns a FWHM into detector windows:

- `rt_sigma_minutes = (fwhm/60)/FWHM_TO_SIGMA`  (`:226`)
- `half_window_scans = round(2σ / spacing)`     (`:227`, legacy — see table)
- `rt_half_window_minutes = 2σ`                 (`:229`, the live scoring window)

Reached two ways:
- `with_rt_from_index(engine, fallback)` — `trace_kernel.rs:211`. **Measured** path: calls
  `estimate_fwhm_seconds`, clamps to `[FWHM_FLOOR_SEC=1.0, FWHM_CEIL_SEC=60.0]`, then `apply_fwhm`.
- `with_rt_from_scans(scan_info, assumed)` — `trace_kernel.rs:198`. **Assumed-FWHM** fallback (same
  `apply_fwhm`, no measurement).

Pipeline wiring (the untargeted pipeline currently lives only in `examples/detect_features_tsv.rs`;
there is no library-level detect→refine→resolve entrypoint yet): `detect_features_tsv.rs:388` selects
`with_rt_from_index` by default (measured), or `with_rt_from_scans` when `FIXED_SIGMA` is set. So the
**detector is data-driven by default**.

## Every peak-width-dependent site

| # | Site (file:line) | Value / meaning | Source |
|---|---|---|---|
| 1 | `trace_kernel.rs:226` `rt_sigma_minutes` | RT Gaussian σ | **DERIVED** from measured FWHM |
| 2 | `trace_kernel.rs:229` `rt_half_window_minutes = 2σ` | **live** matched-filter scoring half-window (time-bounded); consumed at `seed_rt_window` `trace_kernel.rs:513` | **DERIVED** |
| 3 | `trace_kernel.rs:227` `half_window_scans = round(2σ/spacing)` | legacy scan-count window | **DERIVED**, but per its own doc (`:98`) "no longer consumed by the detector" — dead-ish |
| 4 | `feature_refinement.rs:58` `const MAX_SCANS_TO_AVERAGE = 3` | cap on scans averaged into the composite | **HARD-CODED** |
| 5 | `feature_refinement.rs:353` `half = MAX_SCANS_TO_AVERAGE/2` (`refine_feature_inner` window) | averaging window apex ± 1 | **HARD-CODED** (from #4) |
| 6 | `feature_refinement.rs:729` `half = MAX_SCANS_TO_AVERAGE/2` (`build_feature_slices` window) | averaging window apex ± 1 (shift/neighbor path) | **HARD-CODED** (from #4) |
| 7 | `trace_kernel.rs:151/156` `trace_missed_scans_allowed` / `trace_max_half_width_minutes` (defaults 1 / 0.5 min) | claim-extent XIC-walk knobs; example overrides via `TRACE_MISSED_SCANS` / `TRACE_MAX_HALF_WIDTH_SEC` (`detect_features_tsv.rs:322/326`) | **HARD-CODED** (deliberately decoupled from σ; but the 0.5 min = 30 s half-width guard is a 10-min-scale value — see gaps) |
| 8 | `detect_features_tsv.rs:514/519/806` `NeighborIndex::build*(…, 0.05)` | co-elution RT bin/tol for neighbour masking (min) | **HARD-CODED** 0.05 min (3 s) |
| 9 | `detect_features_tsv.rs:547` `LockedGrids::new(0.05)` | co-elution RT tol for grid locking (min) | **HARD-CODED** 0.05 min |
| 10 | `detect_features_tsv.rs:131` `rt_tol = 0.05` (dedup helper) | pass-1 dedup RT tol (min) | **HARD-CODED** 0.05 min |
| 11 | `detect_features_tsv.rs:651` `resolve_charge_state_consensus(&refined, 10.0, 0.1)` | cross-charge consensus RT tol (min) | **HARD-CODED** 0.1 min (6 s) |

Not in scope (classic FlashLFQ targeted quant, separate pipeline, not driven by chromatographic FWHM):
`engine.rs:58/66` `MISSED_SCANS_ALLOWED=1`, `MAX_PEAK_HALF_WIDTH=i32::MAX`; same in `mbr_search.rs`.

## Gaps (stale hard-codes that should track measured FWHM)

1. **Averaging window (#4–#6) — the primary gap.** `refine_feature*` compute their apex±1 window from
   `MAX_SCANS_TO_AVERAGE`, a const tuned to the 10-min case. The `refine_feature*` signatures take a
   `SpectralAveragingParameters` (sigma-clipping/binning config only — it has **no** scan-count field),
   so the measured FWHM is not even reachable from refinement today. On a 65-min gradient (~7–9
   scans/FWHM) and 120-min (~20) the composite still averages only 3 scans — far short of what the peak
   supports.
   - NB: the detector's derived `half_window_scans` / `rt_half_window_minutes` are the **±2σ scoring**
     window, a deliberately *wider* scale than the averaging window (the composite intentionally stays
     tight — apex±1 — to avoid pulling in co-eluting interference; see `refine_feature_inner` comment
     `feature_refinement.rs:347`). So you cannot just forward `half_window_scans` into refinement; the
     averaging window needs its own FWHM-derived count (scans-within-~1 FWHM, i.e. `FWHM/spacing`),
     not the 2σ scoring count.

2. **Example pipeline RT tolerances (#8–#11).** `0.05` min neighbour/dedup and `0.1` min consensus
   tolerances are 10-min-scale (peaks ~3 scans, FWHM ~1.9 s). On longer gradients wider peaks co-elute
   over a longer RT span, so a 3 s co-elution window and 6 s consensus window are likely too tight and
   will under-link charge states / neighbours. These live in the example driver, not the library.

3. **`trace_max_half_width_minutes` default 0.5 min (#7).** A 30 s half-width runaway guard. A 120-min
   peak (~20 scans/FWHM) can have a full width approaching/exceeding 60 s, so on the long gradient this
   guard can clip a real elution's claim extent before it ends. It is overridable
   (`TRACE_MAX_HALF_WIDTH_SEC`) and is a guard not a width, but its default is 10-min-scale.

## Already correct (no action)

- Detector RT σ and the live matched-filter scoring window (#1, #2) are fully data-driven and
  consumed correctly (`seed_rt_window` reads `rt_half_window_minutes`).
- The two averaging-window computations (#5, #6) are at least **consistent with each other** (both
  read the same `MAX_SCANS_TO_AVERAGE`), so there is no drift between the two refine paths — they are
  jointly stale, not divergently stale.

## What was fixed vs deferred

- **Fixed:** nothing. There is no behaviour-preserving wiring fix: refinement has no access to the
  measured width (a signature change is required), the correct averaging count is a *different*
  quantity from any value the detector already computes (scans-per-FWHM vs 2σ scoring window), and the
  data-dependent formula does not yet exist. Any change here moves the 10-min numbers, which the task
  reserves for the separate validation pass. Per instructions, all items are documented, not applied.

## Prioritized wiring plan (feeds the data-dependent-averaging umbrella)

1. **Compute a canonical scans-to-average once, from the measured FWHM.** Add a single helper (e.g.
   `scans_per_fwhm(fwhm_seconds, scan_info) = round((fwhm/60)/median_spacing)`, clamped to a sane odd
   band) next to `apply_fwhm` in `trace_kernel.rs`, reusing `estimate_fwhm_seconds` +
   `median_ms1_scan_spacing_minutes` (both already exist). Target: reproduce ~3 / ~7–9 / ~20 on
   short/medium/long. This is the umbrella's core formula — validate it against `fwhm_probe` output on
   all three gradients before wiring (TODO item "Test FWHM + n-scans-to-average").

2. **Thread that count into refinement.** Replace the `MAX_SCANS_TO_AVERAGE`-derived `half` at
   `feature_refinement.rs:353` and `:729` with a `scans_to_average` value carried on a params struct
   passed to `refine_feature*` / `build_feature_slices`. Keep `MAX_SCANS_TO_AVERAGE` as the default /
   floor so the 10-min case (≈3) is unchanged and the `n_avg <= cap` safety assert still holds (raise
   the assert cap to the derived count). Simplest signature: extend the existing
   `SpectralAveragingParameters` with a `scans_to_average: usize`, or add a small `RefineWindow`
   argument. This is where the 10-min recall must be re-checked.

3. **Make the example's neighbour / consensus / dedup RT tolerances (#8–#11) FWHM-relative.** Scale
   `0.05` / `0.1` min off the measured FWHM (e.g. ~1×FWHM for co-elution, ~1.5–2×FWHM for consensus)
   instead of constants, so linkage widens with the gradient. Re-tune on 10-min to hold recall.

4. **Reconsider `trace_max_half_width_minutes` default (#7)** as a multiple of measured FWHM (e.g.
   ≥3×FWHM) so the long-gradient claim extent isn't clipped, while staying a loose guard.

5. **Retire or repurpose the legacy `half_window_scans` (#3)** once the above lands — it is derived but
   unconsumed; either delete it or make it the scans-to-average carrier.
