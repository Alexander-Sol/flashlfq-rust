# Untargeted Detector Improvement Plan

Plan for three changes to the untargeted MS1 feature-detection pipeline
(`rust/flashlfq-core/src/{trace_kernel,feature_refinement,peak_indexing}.rs`). Companion to
`agent_info/Feature-Detection-Design.md` (the original design) and the auto-memory
`untargeted-feature-detection`.

> Note on filename: the repo root already has `PLAN.md` (the ralph-loop task list). Windows is
> case-insensitive, so a root `plan.md` would overwrite it — this doc lives in `agent_info/` instead.

## Motivation (measured, this session)

- **True chromatographic FWHM ≈ 1.8 s** on the CA/Lumos data (`examples/fwhm_probe.rs`, XIC half-max
  over 243k clean peaks; median 1.8 s, p90 3 s). The hard-coded 36 s assumption is ~20× too wide.
- **Detection is near-perfect at the true width**: at 1.8–3 s assumed FWHM the trace kernel detects
  **100 % of short (≤3 scan) reference peaks** and 99.7–99.8 % of all reference peaks (DETECT_ONLY
  runs, 99 % coverage). At 36 s, short-feature detection collapses to 69 %.
- **But narrow FWHM fragments**: 3 s → 1.12 M features, 1.8 s → 3.33 M features, because the claim
  window is welded to `2σ` and a real elution wider than that is split into many adjacent seeds.
- **The downstream can't absorb that**: `refine_feature` (~35 min at that count) and the **O(n²)**
  `resolve_charge_state_consensus` (1 M² pairwise → hours) are the wall. That wall is the *only*
  reason the resolved-level FWHM sweep peaked at 18 s — it reflects what the downstream can process,
  not where detection is best.

Conclusion: fix the downstream to scale, and fix detection so it captures each peak's *true* extent
instead of splitting-then-merging. Three changes, below.

---

## Implementation status & measured results (2026-07-04)

**Done and green:** Change 1, Change A (both trace-following claim + data-driven σ), plus two changes
that emerged during implementation — an **input-centroiding fix** and a **chromatographic-persistence
gate**. **Pending:** off-by-one corrector (next), then Change B.

- **Input centroiding (new, not in the original three).** The `.raw` was being read in **profile**
  mode (~9k samples/scan) — the detector's one-peak-per-isotope model was running on profile points the
  whole time. `read_ms1_scans` now vendor-centroids Thermo `.raw` on read
  (`ThermoRawReader::new_with_detail_level_and_centroiding(path, Full, true)`, extract via
  `spectrum.peaks()`), ~732 peaks/scan. CA/Lumos: 17.6M → 2.77M peaks. **All pre-centroiding numbers
  are superseded.**
- **Change A persistence gate (new).** `TraceKernelParameters::min_feature_scans` (default 2): reject a
  feature whose *traced extent* spans fewer than N distinct scans. At the true ~1.9 s σ, ~50% of
  features were single-scan noise doublets; the gate removes them at detection (before consensus).
- **Data-driven σ** measures **1.89 s FWHM** on CA/Lumos (matches `fwhm_probe`), via
  `with_rt_from_index`.

**End-to-end baseline (centroid raw, data-driven σ = 1.89 s, cov99, gate = 2 scans):**

| stage | value |
| --- | --- |
| input peaks | 2.77 M (centroid; was 17.6 M profile) |
| detected → refined → resolved | 142,958 → 84,302 → **58,851** |
| reference recall (±20 ppm, ±0.3 min) | **88.4 %** (548/620); charge-matched 84.2 % |
| single-scan features / width p50 | **0 %** (was ~50 %) / 5.8 scans |
| timing | detect 150 s / refine 33 s / consensus 0.14 s / **total 191 s** |

The gate cut resolved features 168k→59k **with recall slightly up**; SHORT (≤3-scan) recall dips to
76.2 % (the gate costs ~2 genuinely-short reference peaks). The remaining not-rediscovered peaks are
mostly 2.5–5.9 kDa, z3–z4 — the **off-by-one** population (next lever).

---

## Change 1 — Bucketed O(n log n) `resolve_charge_state_consensus` (unblocks everything)

**Where:** `feature_refinement.rs:303` (`resolve_charge_state_consensus`), predicate at
`features_link:368`. Called from the runner as `resolve_charge_state_consensus(&refined, 10.0, 0.1)`.

**Problem.** Current grouping is single-linkage union-find over **all pairs** (`for i in 0..n { for j
in i+1..n { if features_link(i,j) union } }`) — O(n²). At 1 M+ refined features it never finishes.

**Design — spatial pruning that yields byte-identical groups, just faster.** The predicate links two
features iff (a) `|apex_rt_i − apex_rt_j| ≤ rt_tol` **and** (b) `∃ k∈[−2,2] : ppm_diff(mass_i, mass_j
+ k·C13) ≤ mass_ppm`. Both arms are *local*: partners live in a bounded (RT, mass) neighborhood.

1. **RT bins.** `rt_bin(i) = floor(apex_rt_i / rt_tol)`. A partner within `rt_tol` sits in
   `rt_bin(i) − 1 … rt_bin(i) + 1` (a shift of `rt_tol` crosses ≤1 boundary). Bucket indices into a
   `HashMap<i64, Vec<usize>>`, each bucket's Vec **sorted by mass**.
2. **Mass window.** Within those ≤3 buckets, candidate masses satisfying arm (b) lie in
   `[mass_i − W, mass_i + W]` with `W = 2·C13_MINUS_C12 + mass_ppm·1e-6·mass_i` (the widest the ±2
   ¹³C off-by-one band can reach). Binary-search that sub-range in the mass-sorted bucket and iterate
   only it.
3. **Exact predicate + union.** For each surviving candidate `j`, apply the **unchanged**
   `features_link(i, j)` and `union(i, j)`.

**Exactness invariant.** The candidate set is a *superset* of all true partners (RT bins ±1 cover all
`|Δrt| ≤ rt_tol`; the `±W` mass window covers all off-by-one hits), and we still apply the exact
predicate — so no link is missed and none is added. **Connected components are identical to the O(n²)
version.** Downstream `resolve_group`/`resolve_mass_by_cross_charge` are untouched.

**Complexity.** O(n log n) for the per-bucket mass sorts + O(n·k), k = avg local candidates. After RT
binning (a ~0.1 min slice) and the ~4 Da mass window, k is tens, not millions.

**Risks / edge cases.** Empty/degenerate RT; features with identical apex on a bin boundary (covered
by ±1); float bin keys (use integer `floor`). Keep `rt_tol`/`mass_ppm` as the same params.

**Validation.**
- New unit test: random N=2–3k refined features, assert bucketed grouping == naive O(n²) grouping
  (same component partition, e.g. compare sorted group-signature multisets).
- Existing tests must stay green: `consensus_corrects_mono_off_by_one_across_charges`,
  `singleton_group_falls_back_to_refined_mass`.
- End-to-end: the 3 s / 99 % run must now complete **the consensus step** in seconds, not hours.

**Scope — this does not make the full run fast on its own.** The motivation names *two* walls:
`refine_feature` (~35 min at 1 M+ features) **and** the O(n²) consensus. Change 1 removes only the
second. Run before Change A (as sequenced), the 3 s / 99 % pipeline still pays the ~35-min
`refine_feature` pass because that pass is O(features) at a 1 M+ feature count — Change 1 does not
touch it. So the deliverable here is "consensus is no longer the bottleneck / no longer effectively
non-terminating," **not** a fast end-to-end run. `refine_feature` stays the gating cost until **Change
A** collapses the feature count (at which point both refine and consensus are cheap). `refine_feature`
itself is deliberately not re-architected in this plan — Change A makes it a non-issue by shrinking its
input, which is cheaper than optimizing it in place.

**Acceptance.** Identical groups to naive on the random test; the **consensus step** of the 3 s run
completes in seconds (the full run may still be refine-bound until Change A lands); wall-clock for the
consensus step scales roughly linearly with feature count.

---

## Change A — Trace-following claim extent (kill fragmentation at the source)

**Where:** `trace_kernel.rs` — `seed_rt_window` / `score_hypothesis` (claim step) and `detect_features`
(the `claimed.insert(p.key())` loop); reuse the XIC boundary walk in
`peak_indexing.rs:302` (`get_xic_by_scan_index`).

**Problem.** The claim extent is a fixed `±rt_half_window_minutes` (= 2σ) window. Narrow σ (needed for
short-feature sensitivity) → real peaks wider than 2σ get split into many features; wide σ →
over-claims and steals neighbors. Splitting-then-merging in consensus is backwards.

**Design — decouple *shape scoring* from *claim extent*.**
- **Score** each hypothesis on a modest shape window using the RT Gaussian (unchanged matched filter,
  narrow σ) — this keeps short-feature sensitivity.
- **On acceptance, claim the peak's true extent**: from the apex, trace the most-abundant isotope
  outward in RT (reuse the `get_xic_by_scan_index` walk: stop on `missed_scans_allowed` consecutive
  misses, on a valley/baseline return, or a max half-width guard), then claim the peaks along that
  traced extent for each isotope tooth.
- **Feature assembly must consume the traced extent, not the scored window.** Today `build_feature`
  (`trace_kernel.rs:480`) derives *every* field — `apex`, `start_rt`/`end_rt`, `summed_intensity`,
  and the stored `peaks` — from `hyp.peaks`, i.e. the narrow ±2σ *scoring* window. If Change A widens
  only the claim while `build_feature` still reads `hyp.peaks`, the reported RT bounds diverge from the
  intensity and coverage that back them. So the traced peaks (the union of the per-tooth walks, minus
  anything already `claimed` by an earlier feature) become the feature's peak set, and:
  - `start_rt`/`end_rt` = min/max RT over the traced peaks.
  - `summed_intensity` = Σ intensity over the traced peaks — this is the value pushed to
    `explained_intensity` for the coverage cap (`detect_features` line 465) **and** the downstream
    quant, so both now reflect the true extent and stay mutually consistent.
  - `apex` = tallest traced peak (unchanged rule, wider candidate set); `score` stays the *narrow*
    normalized `response` (Change B) — score is a shape-fit, deliberately not an extent sum.
  - The claim set written into `claimed` = exactly this traced peak set, so NMS suppression and the
    feature's own bookkeeping use one and the same peaks (no third, divergent set).
  This keeps `summed_intensity`/coverage, the RT bounds, and the claim mask all derived from a single
  peak set; the *only* thing computed on the narrow window is the acceptance score.
- Because the tallest seed of an elution now claims the **whole** elution, the smaller adjacent seeds
  are already claimed and never fire → **fragments never form**, and there is nothing to merge later.

**Why better than a fixed multiple (incl. the 1/2/3×FWHM idea).** The trace follows the actual data,
so it fits sharp *and* broad *and* tailed (asymmetric) peaks with no global constant, and valley-aware
stopping splits genuinely co-eluting same-mass peaks instead of blindly merging them.

**Supporting change — data-driven σ (a new estimator, not a fold-in).** With shape and extent
decoupled, σ should come from the run's own XIC half-max (~1.8–3 s) instead of the hard-coded constant.
Note the shape of this change: it is **not** a tweak to `with_rt_from_scans(scan_info, assumed_fwhm)`.
That method's signature only has `&[ScanInfo]` — no peaks, no bins — so it *structurally cannot* measure
FWHM from half-max; `fwhm_probe` needs the built `PeakIndexingEngine` and a peak-selection strategy.
So this is a genuinely new estimator with real design surface:
- **Input:** the built index (post-construction), not just `scan_info`. This introduces an **ordering
  dependency** — build index → estimate σ from it → set the detector params — that the current one-shot
  `with_rt_from_scans` construction doesn't have. Plan the call site accordingly (params are finalized
  *after* indexing, e.g. a `TraceKernelParameters::with_rt_from_index(&engine, …)` sibling that runs the
  probe, with `with_rt_from_scans` kept for the assumed-FWHM path).
- **Which peaks to probe:** the tallest N clean XICs (as `fwhm_probe` does), not every peak — cost and
  robustness both argue for a sample, not the full set.
- **Aggregation & robustness:** median (not mean) half-max over the sample to resist tails/co-elution;
  decide a floor/ceiling so a pathological run can't drive σ to a degenerate value.
- **Cost:** the probe is an extra pass over a peak sample at startup; bound it (sample size) so it's
  negligible vs detection.
- Keep `ASSUMED_FWHM_SEC` as an override (and the fallback when the probe can't find enough clean XICs).

**Open design questions.**
- Trace the monoisotope, the most-abundant tooth, or the intensity-summed envelope? (Most-abundant =
  the seed, simplest; envelope-sum is most robust but costlier.)
- Valley/baseline stop threshold (fraction of apex) and `missed_scans_allowed` — tune against the
  measured FWHM/tails, not guessed.
- Interaction with greedy claim ordering (tallest-first): confirm no live-lock when two real peaks
  share an isotope m/z at overlapping RT (valley split must resolve it).

**Validation.** Detection-level short/long/all recall (existing `detect_recall.py`) must hold ≈100 %
short; **feature count must drop sharply** vs the fixed-window narrow runs (3.3 M → target ≪ that);
resolved RT widths (`sweep_compare.py` / `merge_99_rt.py`) should track the reference peak widths
(~0.05–0.3 min), not the old claim-window widths.

**Acceptance.** At data-driven σ: short-feature detection ≈100 %, feature count small enough that
refine+consensus run in minutes, resolved feature RT bounds match reference peak bounds within a scan
or two.

---

## Change B — Normalized matched-filter score (fair across width & charge)

**Where:** `trace_kernel.rs:368` (the `response += wk * g * peak.intensity` accumulation and the
`HypothesisScore.response` used for cross-z NMS in `detect_features`).

**Problem.** `response = Σ (wₖ · g · intensity)` is an **unnormalized** inner product, so it grows with
how many scans/teeth a hypothesis spans. A broad or higher-charge comb can out-score the correct one
just by covering more of the window — this is how broad wrong-charge hypotheses steal narrow seeds in
cross-z NMS, and it's why any multi-width scheme can't select on raw response.

**Design.** Replace the raw sum with a **normalized correlation** — `score = Σ(wₖ·gₛ·I) / ‖w·g‖`. The
question is *which grid the norm sums over*, and this is where a naïve normalization backfires:

- **Observed-only norm is wrong** — `‖w·g‖ = sqrt(Σ (wₖ·gₛ)²)` over just the slots where a peak was
  found *rewards sparse hypotheses*: a comb that matches two teeth well but misses the rest gets a
  small norm and an inflated score. That is the opposite of the goal.
- **Full-template norm over-penalizes faint real peaks** — summing the norm over *every* comb
  `(isotope k, scan s)` slot correctly forces a hypothesis to explain the teeth it predicts, so
  genuine misses deflate the score. But it also punishes a real low-abundance peak for failing to show
  teeth the model itself predicts would sit **below the noise floor** — signal the instrument could
  never have recorded. We don't want to penalize those.

**Noise-floor-truncated template norm (the fix).** Estimate a noise floor `η` (global to start — a
run-level MS1 baseline / low percentile of peak intensity; refine to a *local* per-scan or
per-m/z-region floor if the global one proves too coarse). Estimate the hypothesis apex amplitude `A`
(matched-filter least-squares amplitude `A = Σ(wₖ·gₛ·I) / Σ(wₖ·gₛ)²` over observed slots, or simply the
seed/most-abundant-tooth intensity). Define the **expected-observable support**

`S = { (k, s) : A · wₖ · gₛ ≥ η }`

— the comb slots whose *model-predicted* intensity clears the floor. Then

`score = Σ_{(k,s)∈S} (wₖ·gₛ·I_obs) / sqrt( Σ_{(k,s)∈S} (wₖ·gₛ)² )`, with `I_obs = 0` for a slot in `S`
that has no detected peak.

This threads the needle: a slot in `S` with no peak is a tooth we *should* have seen and didn't → it
stays in the denominator and correctly penalizes (a real miss). A slot predicted below `η` is dropped
from **both** numerator and denominator → a faint real peak is not penalized for teeth below noise.
Optionally divide additionally by `‖I‖` over `S` for a bounded [0,1] cosine shape-fit.

**Risks.** Changes which charge wins NMS and the acceptance gate (currently `response > 0` &
`≥2 isotopes`) — re-validate against current recall, don't assume neutral. Two new knobs (`η` estimator
and its global-vs-local scope, plus the `A` estimator) need tuning against the measured data, not
guessed — a too-high `η` shrinks `S` toward the observed-only failure mode; a too-low `η` collapses
back to the full-template penalty. Keep the raw summed intensity separate for
`summed_intensity`/coverage bookkeeping (those stay unnormalized, and now derive from Change A's traced
extent).

**Validation.** A/B the normalized vs current score at fixed config: cross-z charge-assignment
accuracy (the `charge also matched` metric), overall + short recall, and the m/z-merge charge-OK
count. Sweep `η` (and global-vs-local scope) as part of the A/B — short-feature recall is the
canary for `η` set too high (faint real peaks lost), broad-wrong-charge seed-stealing for `η` too low.
Ship only if charge accuracy improves or holds with no recall regression.

**Acceptance.** Charge-match rate up (or flat) with no recall loss; broad wrong-charge seed-stealing
demonstrably reduced on a crafted co-elution test; short-feature recall (the low-amplitude regime the
noise-floor truncation protects) holds vs the pre-normalization baseline.

**RESULT (2026-07-04) — implemented, IN but not default; fails the acceptance bar.** `ScoreModel`
(`RawSum` | `NormalizedNoiseFloor`) + `estimate_noise_floor`, gated behind `SCORE_MODEL` / `NOISE_PCT`
/ `AMP_SEED` / `SCORE_COSINE`; RawSum is byte-identical to the pre-change score and stays the default.
A/B at cov99, gate 2, centroid, data-driven σ (RawSum baseline **88.4 % recall / 84.2 % charge**):

| config (η = p10) | recall | charge-match |
| --- | --- | --- |
| RawSum (default) | 88.4 % | **84.2 %** |
| normalized, least-squares A, template-norm | **89.8 %** | 83.4 % |
| normalized, seed A, template-norm | 89.8 % | 83.4 % |
| normalized, least-squares A, cosine | 74.7 % | 63.4 % |
| normalized, seed A, cosine | 74.7 % | 63.2 % |

η-**independent** across p1–p20 (89.7–89.8 % / 83.2–83.4 %). Findings: the **A estimator is irrelevant**
(seed ≡ LS); **cosine is harmful** (amplitude-independence lets low-intensity noise-shaped combs win
cross-z NMS); normalized template-norm is a fixed **+1.4 % recall / −0.8 % charge** tradeoff — the
*opposite* of Change B's charge-improvement intent. So **RawSum kept as default**; the normalized model
is available for a recall-over-charge regime but not promoted.

---

## Sequencing

1. **Change 1 (consensus bucketing)** first — mechanical, exactness-preserving. It makes the narrow-FWHM
   pipeline *terminate* (consensus goes from effectively-never to seconds); the run is then
   refine-bound (~35 min at 1 M+ features) but finite, which is enough to measure the baseline Change A
   needs. Lowest risk.
2. **Change A (trace-following + data-driven σ)** next — the biggest correctness/scale win; collapses
   the feature count so the rest of the pipeline is cheap and the RT bounds become real. Do the
   data-driven σ as part of this (it's what makes narrow detection usable).
3. **Change B (score normalization)** last — quality refinement on top; needs careful A/B and may
   interact with A's width handling.

After each: re-run the 3 s / 99 % case end-to-end and record detection recall, resolved recall
(m/z-merge on `Peak MZ` vs `Most-Abundant m/z`), charge-match, feature count, and RT-width distribution
using the existing scratchpad scripts.

## Validation harness (committed to the repo)

- `examples/fwhm_probe.rs` — measures true FWHM from XIC half-max.
- `examples/detect_features_tsv.rs` — env knobs `ASSUMED_FWHM_SEC`, `COVERAGE_TARGET`, `COMB_MODEL`,
  `DETECT_ONLY`; emits detected/refined/resolved TSVs incl. `Most-Abundant m/z`.
- `rust/flashlfq-core/scripts/` (versioned; see its `README.md`) — `detect_recall.py`
  (detection-level short/long/all recall), `sweep_compare.py` (resolved recall + width vs FWHM),
  `merge_99_rt.py` (per-ref merge on observed m/z with RT bounds and nearest-for-misses). Stdlib-only;
  data paths are hard-coded near the top of each (the `D:\SP_Tutorial\Lumos\...` locations).

## Next (post-baseline, in order)

1. **Off-by-one monoisotope correction — now IN scope.** The heavy-peptide (2.5–5.9 kDa, z3–z4) recall
   misses are largely off-by-one; the `correct_monoisotope_offbyone` corrector exists but regressed and
   is `#[allow(dead_code)]`. Wire a strictly non-regressive gate (only shift when the shifted envelope
   clearly fits better), validate against the centroid baseline.
2. **Change B** (normalized matched-filter score) — as specified below.

## Out of scope (tracked, not in this plan)

- Coverage-target retuning post-fix (the earlier sweet-spot was chaining/scaling-limited) — now cheap
  to sweep since the full run is ~3 min.
