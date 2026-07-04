# Isotope-Structured Signal: Ceiling Measurement & Coverage-Stop Plan

How much of the total ion current has genuine isotope-envelope structure, what that implies for the
untargeted detector's coverage target, and a plan to act on it. Companion to
`agent_info/Feature-Detection-Design.md` (original design) and
`agent_info/Detector-Improvement-Plan.md` (FWHM / scaling work); tracked by the auto-memory
`untargeted-feature-detection`.

## Question

The detector reports "% of ΣTIC explained" against **raw** ΣTIC, and the coverage-stop target
(`TraceKernelParameters::coverage_target`, default runner value 0.90) is a fraction of raw ΣTIC. But
raw ΣTIC includes chemical/electronic noise and lone peaks that have no isotope structure at all —
signal no isotope-comb detector can ever legitimately claim. So: **what fraction of ΣTIC is even
isotopically structured?** That fraction is the true ceiling on detector coverage, and the honest
denominator for any recall/coverage number.

## Method

`rust/flashlfq-core/examples/isotope_signal_fraction.rs` (standalone example; uses only
`read_ms1_scans` and `C13_MINUS_C12` from the crate). Per MS1 scan, in one binary-searched linear
pass per charge:

- For each charge `z ∈ [MIN_Z, MAX_Z]` (default 2..=6), spacing `Δ = C13_MINUS_C12 / z`.
- Link peak `i → j` when a peak sits within `PPM` (default 10) of `mz_i + Δ`. Links point strictly
  upward in m/z, so chains never cycle.
- A peak is **isotope-structured** if it belongs to a chain of at least `MIN_ISOTOPES` peaks
  (default **2**, to match the detector's minimum envelope length).
- `assign_charge()` credits each isotopic peak to a **single** charge — the one whose supporting run
  through that peak is longest (ties → lower charge). This keeps per-charge intensity sums
  non-overlapping, so they add to the union total.

Env knobs: `PPM`, `MIN_Z`, `MAX_Z`, `MIN_ISOTOPES` (floor 2), `MIN_INTENSITY`. This is a per-scan
*ceiling*, not a detection: no RT tracing, no deconvolution, no envelope-shape scoring — a peak
counts if it lines up under *any* allowed charge, so the true recoverable fraction is somewhat below
this number.

**Peaks are vendor-centroided.** `read_ms1_scans` opens Thermo `.raw` with
`ThermoRawReader::new_with_detail_level_and_centroiding(path, Full, true)`, so the classification runs
on the vendor centroid stream (one peak per isotope) — the same peaks MetaMorpheus/FlashLFQ consume.
This was **not** the case in the first pass of this analysis: the reader returned raw Orbitrap
*profile* data (~9k samples/scan), and manual validation with `analysis/plot_isotope_classification.py`
caught it (see "Validation finding" below). All numbers below are post-centroiding.

## Findings — CA/Lumos `04-17-23_CA_Tryp_HCD_10min.raw`

2689 MS1 scans, **2.77 M centroided peaks**, ΣTIC (centroid intensities) = 1.4316e13.

| min isotopes | isotope-structured signal | peaks counted |
|---|---|---|
| 3 (peak + 2 others) | **69.0%** of ΣTIC | 17.5% |
| 2 (peak + 1) — detector-aligned | **78.9%** of ΣTIC | 38.4% |

Per-charge breakdown at `MIN_ISOTOPES=2`, 10 ppm (each peak credited to its strongest charge):

| z | ΣTIC | % of ΣTIC | peaks |
|---|---|---|---|
| 2 | 5.09e12 | 35.6% | 397 k |
| 3 | 3.14e12 | 21.9% | 252 k |
| 4 | 1.71e12 | 11.9% | 164 k |
| 5 | 8.86e11 | 6.2% | 136 k |
| 6 | 4.70e11 | 3.3% | 114 k |
| **all** | **1.13e13** | **78.9%** | **1.06 M** |

At `MIN_ISOTOPES=3` (the trustworthy charge prior): z2 32.6% / z3 19.5% / z4 10.3% / z5 4.8% /
z6 1.9%, summing to 69.0%.

**Reading of the numbers.**

- **~21% of ΣTIC has no isotope structure** (unstructured at min 2; ~31% at min 3). This is a hard
  ceiling: no isotope-comb detector can attribute it to peptide envelopes.
- **The signal is concentrated in a minority of peaks.** At min 2, 38.4% of centroids carry 78.9% of
  the intensity; the other ~62% of peaks carry ~21%. The long peak tail is mostly noise by intensity.
- **Charge distribution is the expected tryptic profile** — monotonic decay, z2/z3 dominant (57.5%
  combined). The z5/z6 tail deflated sharply from the profile-data measurement (z6 5.3%→3.3% at min 2)
  once real centroids replaced profile shoulders, confirming the earlier high-z inflation was a
  sampling artifact.

**Validation finding (profile → centroid).** The first measurement ran on raw profile data and
reported 84.8% (min 2) over 17.6 M "peaks". The plotter showed each real peak was ~6–7 profile
samples spaced ~2 mΘ apart (median inter-peak gap 0.0019 Th, 85% of gaps <0.01 Th), and the dominant
peak of scan 1200 was mis-assigned **z=6** because its profile shoulders chained at the 0.17 Th z6
spacing. After wiring vendor centroiding into `read_ms1_scans`, the same peak is correctly **z=3**,
the scan collapses from 4959 samples to 835 centroids, and the ceiling settles at **78.9%**. Lesson:
run this class of metric on centroids, and prefer the vendor/mzdata centroiding to any hand-rolled
peak-picking.

**Implication.** The runner's default `COVERAGE_TARGET=0.90` is measured against total ΣTIC, but the
reachable ceiling is ~79% (min 2) / ~69% (min 3). Targeting 90% therefore forces the detector to keep
seeding past the isotope-structured signal into the noise floor — consistent with the observation
(memory `untargeted-feature-detection`) that high coverage over-fragments and admits spurious
features. The coverage objective is denominated wrong. (The detector now also consumes the centroided
stream, since `read_ms1_scans` changed, so its own ΣTIC/coverage numbers shifted too — re-baseline
before comparing to pre-centroiding runs.)

## Implementation plan

### Change 1 — Denominate the coverage stop against isotope-structured TIC (top priority)

**Where.** `trace_kernel.rs:404-407` computes `coverage_stop = coverage_target * total_intensity`
where `total_intensity` is the sum over **all** seed intensities; the accumulator/stop are at
`:414`, `:465`, `:468-469`. `peak_indexing.rs` is where the per-scan pass would be cheapest (peaks
are already in hand at index time).

**Problem.** `total_intensity` is raw ΣTIC. `coverage_target = 0.90` thus asks the detector to
explain 90% of a quantity that is only ~85% explainable, so the stop never triggers on structure —
it triggers (if at all) only after the detector has chased noise.

**Design.** Compute an **isotope-structured TIC** once (the `assign_charge`/mask pass, ~one linear
scan per spectrum with binary search — negligible vs. the matched filter) and use it as the
denominator:

```
coverage_stop = coverage_target * isotopic_tic     // was: * total_intensity
```

`coverage_target = 0.90` then means "explain 90% of the 84.8% that is reachable," and the detector
stops when real signal is claimed rather than when it has over-reached into noise. Keep the raw-TIC
number in the report for continuity, but drive the stop off the isotopic denominator.

**Wiring.** Add the isotopic-TIC computation behind a flag (`ISOTOPIC_COVERAGE=1` in the runner /
a `TraceKernelParameters` bool) so the change is A/B-measurable against the current behavior. Report
both `% of raw ΣTIC` and `% of isotope-structured ΣTIC` explained.

**Expected effect.** Fewer, cleaner features at a given target; the coverage sweep should plateau at
a target that corresponds to real structure instead of climbing indefinitely into the noise tail.
Measure resolved-count and PSM recall/precision vs. the `AllQuantifiedPeaks.tsv` reference at
matched targets.

### Change 2 — Isotope-mask seed prefilter

**Where.** Seed selection in `detect_features` (`trace_kernel.rs:391+`, the unclaimed-seed loop).

**Problem.** ~52% of peaks (min 2) have no isotopic partner yet are still eligible seeds. The comb
filter rejects them per-hypothesis, but only after forming and scoring the hypothesis.

**Design.** Precompute the per-scan isotope mask (same pass as Change 1) and skip seeding on peaks
the mask marks non-isotopic. This is a cheap **structural** gate that avoids forming doomed
hypotheses; it is *not* a replacement for the comb's envelope-shape scoring, which stays the real
discriminator. Signal cost is ~15% of TIC in the pruned peaks, most of it noise.

**Expected effect.** Fewer hypotheses → faster detection, fewer spurious seeds. Guard against
pruning real monoisotopic peaks whose sole partner fell below `min_intensity`; validate recall does
not regress on the reference.

### Change 3 — Empirical charge prior (lower priority)

**Where.** Charge search ordering / weighting in `detect_features`.

**Design.** Use the per-charge distribution (measured at `MIN_ISOTOPES=3`, not 2) to order or
down-weight the charge search — e.g. de-prioritize z5/z6 hypotheses that are largely coincidental.
Marginal, since the comb already scores per charge; do only if the charge search is a measured cost.

### Change 4 — Extract a shared `isotope_prefilter` module

If Change 1 or 2 lands, lift `assign_charge` / the per-scan mask out of the example into a small
`isotope_prefilter` module so the runner diagnostic and the detector share one implementation rather
than duplicating the chain-walk logic.

## Validation

- Re-run `isotope_signal_fraction` at `MIN_ISOTOPES=3` to lock the trustworthy charge prior.
- For Change 1: A/B the coverage stop (raw vs. isotopic denominator) at matched `coverage_target`,
  comparing resolved count, `% ΣTIC explained` (both denominators), and PSM recall/precision vs.
  `AllQuantifiedPeaks.tsv`.
- Sensitivity: sweep `PPM` (5/10/15) to confirm the ~85% ceiling is stable and not tolerance-driven.
- Cross-file: repeat on a **calibrated** raw to confirm the ceiling and charge profile generalize
  beyond the uncalibrated CA data.

## Open questions

- Should the mask allow z1? Excluded here (peptides are ≥2+), but singly-charged contaminants carry
  isotope structure too and currently fall into the "unstructured" 15%. Worth quantifying before
  treating the 15% as pure noise.
- Interaction with the FWHM/scaling changes in `Detector-Improvement-Plan.md`: the isotopic
  denominator should be measured *after* those land, since fragmentation changes the peak population
  the mask sees.
