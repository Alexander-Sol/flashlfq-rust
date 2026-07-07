# Detector Performance + Parallelization Strategy

Profiling of the untargeted feature-detection pipeline, the one safe optimization applied, and a
concrete parallelization plan. Measured 2026-07-06/07 on this machine (24 logical / 12 physical
cores) via `examples/detect_features_tsv.rs`.

## TL;DR

- **Detect is the whole cost.** 74% of wall-clock on the 10-min file, **94% on the 65-min file**
  (read/refine/resolve are rounding error by comparison). Refinement is *not* the bottleneck — the
  apex-only default already made it ~2-5%.
- **Within detect, `score_hypothesis` is 88-95%**, and inside it the single hottest routine is
  `PeakIndexingEngine::get_indexed_peak` (called ~10^8 times).
- **Applied opt (verified byte-identical):** removed the two per-call heap allocations inside
  `get_indexed_peak`. Detect on the 10-min file **49.8s → 26.8s (1.86×)**; total **59.4s → 36.3s
  (1.63×)**. All three output TSVs are MD5-identical before/after.
- **Recommended parallelization: per-file with rayon** (embarrassingly parallel, parity-safe,
  near-linear up to core/memory limit). For a *single* long file, an intra-file scheme is possible
  but the greedy claim-as-you-go dependency makes it a real, parity-risking refactor — designed
  below as future work, not applied.

## Measured stage breakdown

Driver: `detect_features_tsv.exe <raw> <out.tsv>`, default params (Poisson comb, shift_apex refine,
data-driven σ, 90% coverage target). Times are the example's own per-stage `Instant` timings.

### 10-min file — `04-17-23_CA_Tryp_HCD_10min.raw`
2689 MS1 scans, 2.77M peaks, measured FWHM 1.89 s.

| Stage | Before opt | After opt | % of total (after) |
|-------|-----------:|----------:|-------------------:|
| read + index | 6.52 s | 6.55 s | 18.0% |
| **detect** | **49.82 s** | **26.75 s** | **73.6%** |
| refine (apex) | 1.73 s | 1.72 s | 4.7% |
| charge-state consensus | 0.16 s | 0.17 s | 0.5% |
| write resolved | 0.40 s | 0.40 s | 1.1% |
| **TOTAL** | **59.35 s** | **36.33 s** | |

96,385 features detected, 85,050 resolved. Output byte-identical before/after (detected/refined/
resolved TSVs all MD5-match).

### 65-min file — `HFX_MB_14751_5_02062022.raw` (post-opt)
7407 MS1 scans, 16.8M peaks, measured FWHM 9.36 s.

| Stage | Time | % of total |
|-------|-----:|-----------:|
| read + index | 22.38 s | 2.9% |
| **detect** | **732.53 s** | **93.8%** |
| refine (apex) | 15.46 s | 2.0% |
| charge-state consensus | 1.71 s | 0.2% |
| **TOTAL** | **780.84 s** | |

737,957 features detected, 665,486 resolved. (Timings with `DETECT_PROFILE` on; the profiling
overhead is a fraction of a second — negligible against 732 s.)

**Why detect scales worse than linearly in file length.** The 65-min file has ~6× the peaks but
detect is ~27× slower. Two multipliers compound: (1) more seeds, and (2) a **much wider RT window**
— the data-driven FWHM is 9.36 s here vs 1.89 s on the 10-min, so `seed_rt_window` returns many more
scans, and `score_hypothesis` cost is `O(window_scans × comb_teeth × charges)` per seed. Wide
chromatography is the real cost driver, not just file size.

## Sub-stage attribution (inside `detect_features`)

Instrumented behind the new `DETECT_PROFILE=1` env flag (all timing gated on a bool read once; the
default path pays only a predictable-branch check per section — no `Instant::now()` in the hot loop
unless profiling is requested).

| Sub-stage | 10-min | 65-min |
|-----------|-------:|-------:|
| seed prep (`all_peaks` + intensity sort) | 0.25 s | 1.89 s |
| `seed_rt_window` | 0.15 s | 3.55 s |
| **`score_hypothesis` (6 charge hypotheses/seed)** | **23.84 s (88%)** | **693.48 s (94.7%)** |
| `trace_claim_extent` | 2.10 s | 25.22 s |
| other (claim bookkeeping, gates, build) | 0.84 s | 8.35 s |
| seeds considered (unclaimed at their turn) | 503,445 | 4,088,536 |
| `score_hypothesis` calls | 3,020,670 | 24,531,216 |

`score_hypothesis` dominates. Its inner loop is `window_scans × comb_teeth` calls to
`get_indexed_peak` per charge hypothesis. On the 65-min file that is ~24.5M hypotheses ×
(window × teeth) ≈ **hundreds of millions of `get_indexed_peak` calls** — that function *is* the hot
path.

## Optimization applied (safe, verified byte-identical)

**`PeakIndexingEngine::get_indexed_peak` — eliminate two heap allocations per call.**

The original allocated a `Vec<&bin>` (via `get_bins_in_range`) and a parallel `Vec<isize>` of
per-bin binary-search indices on *every* call, then reduced them in `get_best_peak_from_bins`. For a
10 ppm tolerance at m/z ~500 the candidate-bin span is only 1-2 bins, so these were tiny, extremely
short-lived Vecs — but allocated on the order of 10^8 times.

The rewrite iterates the candidate bins directly (`floor..=ceil`), doing the binary search and
best-peak selection inline, with **zero intermediate allocation**. It preserves the exact traversal
order and the strict-`<` "closest to `m`" tie-break, so the selected peak is identical in every case.
The sibling helpers (`get_bins_in_range`, `get_best_peak_from_bins`, `get_peak_from_bin`,
`binary_search_for_indexed_peak`) are untouched — still used by `get_xic` / `get_xic_by_scan_index`.

**Verification:** 10-min detected/refined/resolved TSVs are MD5-identical before vs after; the
5 `get_indexed_peak` unit tests and all 211 core tests still pass (only the pre-existing, unrelated
`raw_quant_matches_mzml_quant` fails, as before). **Result: detect 49.8s → 26.8s (1.86×).**

### Recommendations NOT applied (need more than a byte-identical guarantee, or bigger blast radius)

1. **Avoid the per-call `slots: Vec<(f64,f64)>` in `score_hypothesis` for the default `RawSum`
   model.** `RawSum` only needs `Σ tᵢ·Iᵢ`, which can be accumulated inline (bit-identical if summed
   in the same window×k order); `slots` is only genuinely needed by `NormalizedNoiseFloor`. Saves one
   ~1-3 KB allocation per hypothesis (~24.5M on the 65-min file). Left as a recommendation because it
   forks a subtle scoring function; apply behind a before/after byte-identity check.
2. **Reuse the two `HashSet`s (`used`, `observed_isotopes`) and the `peaks` Vec across the 6 charge
   hypotheses of a seed** via a caller-owned scratch buffer that is `.clear()`ed rather than
   reallocated per call. Behavior-preserving but changes `score_hypothesis`'s signature; moderate
   blast radius.
3. **Cheaper `claimed` membership.** `claimed: HashSet<PeakKey>` is probed inside the hot loop; a
   faster hasher (e.g. `FxHashMap`) or a bitset keyed on a dense peak id would cut hashing cost. Adds
   a dependency / an id scheme; measure first.
4. **Shrink the effective window.** `score_hypothesis` cost is linear in `window_scans`; the
   Gaussian weight far from the apex is ~0. A tighter `rt_half_window_minutes` (or early-terminating
   the comb once the Gaussian weight drops below a threshold) would cut work — but it can change
   which teeth are observed and therefore the output, so it is an algorithmic tuning decision, not a
   free micro-opt.

## Parallelization strategy

### Where the shared mutable state is

`detect_features` is a **greedy, claim-as-you-go** loop. The state that threads would contend on:

- `claimed: HashSet<PeakKey>` — read by every `score_hypothesis`/`trace_claim_extent` (claimed peaks
  are treated as absent) and written on every acceptance. **This is the ordering dependency**: a
  feature's score depends on what all taller, earlier-accepted features already claimed. Output is a
  deterministic function of the tallest-first seed order.
- `features: Vec<DetectedFeature>` and `explained_intensity` — append-only accumulators (trivially
  handled by per-thread partials + merge).

The `PeakIndexingEngine` itself is **immutable** during detection (pure reads), so it can be shared
across threads by `&` with no synchronization — which is what makes both strategies below feasible.

### Option A — per-file on its own thread (RECOMMENDED)

Each spectra file's read → index → detect → refine → resolve is fully independent: no cross-file
shared state. The existing MS2 driver already has the exact shape — `engine.rs::run_msms` loops
`for file_name in &file_order { index; quantify; }` collecting into per-file maps; the untargeted
path is structurally identical (only the single-file example exists today).

**Plan:** add `rayon` and turn the per-file loop into `file_order.par_iter().map(|f| detect_one(f))
.collect()`, merging the per-file result maps after the join. `detect_one` already touches nothing
global.

- **Parity:** bit-identical — each file's computation is unchanged; only the *scheduling* across
  files changes. Determinism preserved by collecting results keyed by file name.
- **Speedup:** near-linear in `min(num_files, num_cores)`. A 20-file batch on 12 physical cores ≈
  **8-12×** end-to-end.
- **Memory is the constraint.** Each in-flight file holds its scans + index. The 65-min file's index
  is on the order of ~0.5-1 GB resident; N files in flight ≈ N× that. Cap concurrency with a rayon
  thread-pool size / semaphore sized to `min(cores, RAM_budget / per_file_footprint)` rather than
  defaulting to all 24 logical threads. Drop each file's `Scan`/index as soon as its features are
  resolved (already natural if `detect_one` owns them).
- **Effort:** small. One dependency, one `par_iter`, a result merge. This is the highest
  value-per-line change and should ship first.

### Option B — intra-file parallelism (single large file)

Per-file parallelism does nothing for a *single* 65-min file that alone takes ~12 min. Options,
worst-to-best:

1. **Parallelize the 6-charge loop per seed.** Read-only on `engine`+`claimed`; but only 6-way, each
   hypothesis is microseconds, and fork/join across 4M seeds would be dominated by overhead. **Not
   worth it.**
2. **Spatial tiling (m/z bands or RT blocks) detected independently, then reconciled.** Breaks the
   global tallest-first claim order, so features straddling a tile boundary can be claimed
   differently → **output changes**. Only acceptable if a small, characterized deviation from the
   serial result is tolerable. Not parity-safe.
3. **Speculative batch scoring (parity-safe, RECOMMENDED intra-file path).** Keep acceptance serial
   and tallest-first, but score in parallel: take the top-K currently-unclaimed seeds, score all K in
   parallel against a read-only snapshot of `claimed` (rayon), then walk the K in intensity order and
   accept serially — before accepting seed *i*, check whether any peak it scored was claimed by an
   earlier acceptance *within this batch*; if so, re-score just that one seed serially. Spatial
   collisions inside a small top-K batch are rare, so re-scores are infrequent and the result is
   **bit-identical to the serial greedy**. `score_hypothesis` is 88-95% of detect and is pure/read-
   only, so this parallelizes the dominant cost. Expected ~**4-8×** on 12 physical cores after
   subtracting the serial reconciliation and re-score tail; needs benchmarking of the collision rate
   to confirm. This is a real implementation effort (batch manager + re-score path + a lock-free or
   snapshotted `claimed`), hence future work, not applied here.

### Recommended path

1. **Ship Option A (per-file rayon) now** — small, parity-safe, and the common FlashLFQ workload is
   many files. Expect 8-12× on batch runs.
2. **Then Option B.3 (speculative batch scoring)** for single-long-file latency, gated by a
   collision-rate benchmark and a byte-identity check against the serial output.
3. Fold in micro-opt recommendation #1 (`slots` for `RawSum`) and #3 (faster `claimed` hashing)
   opportunistically; both are independent of the threading work.

Combined ceiling on a batch of long files: per-file parallelism (≈ cores) × the 1.86× already banked
in `get_indexed_peak` ≈ an order of magnitude over today's serial detect, memory permitting.

## How to reproduce

```
cargo build --release --example detect_features_tsv
# sub-stage attribution:
DETECT_PROFILE=1 ./target/release/examples/detect_features_tsv <raw> <out.tsv>
# stage timings are always printed + appended to <out>.log
```
