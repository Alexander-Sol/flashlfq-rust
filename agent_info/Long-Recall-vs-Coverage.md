# Long-gradient recall vs coverage target (+ runtime)

Sweep of the detector's **coverage target** (`COVERAGE_TARGET`, the fraction of ΣTIC the greedy
seed loop tries to explain) on the **long 120-min IonStar** file
(`B03_19_150304_..._9B.raw`), holding every other parameter at its current default (apex-only
`shift_apex` refine, 10 ppm, data-driven σ, seed floor 1000). Recall/charge-match scored by
`analysis/longer_gradients/score_recall.py` against the 24,642-peak MSMS ground truth
(`long_ground_truth_msms.tsv`), primary criteria: 20 ppm mass, RT ≤ **0.3 min**, charge ∈ detected
charge states. Detect/total seconds are the example's own timing log. Reproduce with
`analysis/longer_gradients/sweep_coverage.ps1 <cov>`.

## Breakdown

| coverage target | achieved %TIC | #features (resolved) | recall % | charge-match % (of all ref) | charge-match (of matched) | detect s | total s |
|---|--:|--:|--:|--:|--:|--:|--:|
| 0.80 | 80.0 | 34,259 | 57.7 | 55.7 | 96.6 | 20.2 | 47.8 |
| **0.90 (default)** | **90.0** | **96,058** | **81.2** | **76.2** | 93.9 | **47.9** | **78.2** |
| 0.95 | 95.0 | 416,787 | 90.7 | 79.8 | 88.0 | 276.0 | 319.4 |
| 0.98 | **96.5** † | 929,797 | 91.7 | 77.5 | 84.5 | 1396.4 | 1464.4 |
| 0.99 | — | — | — | — | — | — | *not run — runtime prohibitive* |
| 0.995 | — | — | — | — | — | — | *not run — runtime prohibitive* |
| 0.999 | — | — | — | — | — | — | *not run — runtime prohibitive* |

† **The 0.98 target was never reached.** The run exhausted every seed above the intensity floor
(`MIN_SEED_INTENSITY = 1000`) at **96.5%** of ΣTIC and stopped — it detected **1.11 M** raw features
(929 k resolved) chasing the last 1.5% of TIC and still fell short of 98%. Detect alone took **23.3
min** (total 24.4 min). Because the same seed floor caps the achievable TIC at ~96.5% regardless of
target, **0.99 / 0.995 / 0.999 cannot reach their targets either** — they would only churn the same
~1 M seeds to seed-exhaustion for no additional recall, so they were correctly skipped.

## Marginal cost between adjacent points

| step | Δachieved-TIC | Δrecall (pp) | recall per +1% TIC | Δtotal s | recall per extra second |
|---|--:|--:|--:|--:|--:|
| 0.80 → 0.90 | +10.0 | +23.5 | 2.35 pp | +30.4 | **0.77 pp/s** |
| 0.90 → 0.95 | +5.0 | +9.5 | 1.90 pp | +241.1 | 0.039 pp/s |
| 0.95 → 0.98 | +1.5 | +1.0 | 0.67 pp | +1145.0 | **0.00087 pp/s** (≈1145 s / pp) |

Efficiency collapses by ~**900×** across the sweep: the first block buys 0.77 recall-points per
second, the last buys under 0.001 — roughly **19 minutes of extra compute per single recall point**
past 0.95.

## Plots

- `analysis/longer_gradients/long_recall_vs_coverage.png` — recall vs coverage target (clear knee at ~0.95).
- `analysis/longer_gradients/long_runtime_vs_coverage.png` — detect + total runtime vs coverage target (near-exponential blow-up).

## Interpretation

**Recall plateaus at ~0.95 (90.7%).** From the 0.90 default it climbs steeply (81.2 → 90.7%, +9.5
pp) for a **4.1×** runtime cost (78 s → 319 s). Beyond that it flatlines: 0.95 → 0.98 adds only
**+1.0 pp** (90.7 → 91.7%) while runtime jumps to **18.7× the 0.90 baseline** (78 s → 1464 s), and
0.98 does not even reach its nominal target. Charge-match-of-matched also *degrades* as coverage
rises (96.6 → 93.9 → 88.0 → 84.5%): the extra features the higher targets scrape up are marginal,
crowded, lower-fidelity placements, not clean new peptides.

**Coverage target is the wrong lever for the low-abundance tail.** The misses on this file are
concentrated in the bottom intensity deciles (~41% of missed peaks are bottom-decile; per-decile
recall runs D1 ≈ 22.6% → D10 ≈ 98.6%). Raising the global coverage target attacks that tail only
indirectly and inefficiently: it keeps lowering the *global* seed bar, so it re-examines the whole
spectrum's mid-abundance clutter (hundreds of thousands of extra features) to reach a few more faint
peptides — and then hits the hard **`MIN_SEED_INTENSITY = 1000` floor** at ~96.5% TIC, which is the
real ceiling. The faint tail lives *at or below* that seed floor, so chasing it by cranking coverage
is doubly wrong: it pays a global O(seeds) cost to reach signal the seed floor then forbids anyway.

**Recommended levers instead:**
1. **Auto-stop at the knee.** Recall gain per unit runtime falls off a cliff after ~0.95; the
   detector should detect the diminishing-returns knee (recall/TIC-increment slope, or the
   seed-floor exhaustion point) and stop, rather than grinding toward an unreachable target. The
   0.90 default is already close; **0.95 is the honest ceiling of this lever** and a reasonable
   opt-in max — anything higher is pure cost.
2. **Dedicated low-abundance seeding**, not a lower global floor. Seed the faint tail *locally* —
   e.g. residual-TIC seeding in windows the greedy pass left unclaimed, or a targeted lower
   `MIN_SEED_INTENSITY` applied only in RT/m-z neighbourhoods of unexplained signal — so the tail is
   recovered without re-deconvolving the entire mid-abundance spectrum. That directly addresses the
   bottom-decile miss population the global coverage knob cannot afford to reach.

**Bottom line:** push coverage to 0.95 if the extra ~9 pp recall is worth ~4× runtime; do **not**
push past it (18.7× runtime for +1 pp, target unreachable). The low-abundance tail needs a smarter
seeding strategy, not a higher coverage target.

## Notes / reproduction

- Sweep driver: `analysis/longer_gradients/sweep_coverage.ps1 <cov>` (writes per-point TSVs to the
  scratchpad — regenerable, **not committed**). Plot regenerator: `plot_coverage_sweep.py` (data
  hard-coded from this table).
- 0.90 total here (78 s) is slightly under the historical ~85 s in `run_timings.md` — machine
  variance; the relative blow-up factors are the point.
