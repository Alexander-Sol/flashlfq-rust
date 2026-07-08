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

## IMPLEMENTED: intra-file red-black RT binning (Option B, variant of B.2)

Shipped instead of B.3 — chosen for single-file latency with a low, constant memory footprint (one
file resident, no per-batch `claimed` snapshots). It is a **parity-approximate** spatial scheme, not
the bit-identical B.3. Opt-in via the **`DETECT_PARALLEL`** env var; see
`trace_kernel.rs::detect_features_parallel`. The serial path (`detect_features_serial`) is untouched
and remains the default and the reference — the 88-95% hot function `score_hypothesis` and its two
siblings (`trace_seed_extent`, `gather_extent_peaks`) keep their exact single-`HashSet` signatures, so
there is **zero hot-path cost** to the serial build.

### The scheme

- **Equal-work bins with a 4-min floor.** Seeds are partitioned into contiguous RT bins that are
  simultaneously ≈ equal seed count (load balance) and ≥ `DETECT_BIN_MIN_WIDTH_MINUTES = 4.0` wide
  (correctness). A bin closes only when both hold, so dense regions give min-width bins and sparse
  regions widen to carry a full share.
- **Red-black (even/odd) two-phase.** Bins are 2-colored by index parity. Even bins are detected in
  parallel (`rayon` `par_iter`), a barrier merges their claims, then odd bins are detected in parallel
  against those claims. Every pair of same-color bins is separated by a full ≥ 4-min bin, and a seed
  reaches at most `bin_reach_minutes` = `max(rt_half_window_minutes, trace_max_half_width_minutes)` ≈
  0.5 min from its apex — so **2 × reach (~1 min) ≪ bin width (4 min)** and concurrently-processed bins
  provably cannot claim the same peak. Each bin therefore runs the ordinary serial loop with its own
  `claimed` set; **no locking on the hot path.**
- **Odd-bin seeding is a thin strip, not the whole claim set.** An odd bin is pre-seeded only with the
  even-phase claims within `reach` of its span (found by binary search on RT-sorted even claims) —
  everything it could possibly observe from the earlier phase, and nothing more. Keeps memory low.
- **Deterministic.** Results are reassembled in bin-index order; each bin is serial internally, so the
  output is identical run-to-run (only *scheduling* varies across threads).

### Constraint: does not compose with the global-stop heuristics

Coverage target < 1, the TIC-knee auto-stop, and the reject-rate auto-stop all key on a **global**
tallest-first traversal, which binning breaks. When any is engaged, `detect_features` prints a notice
and falls back to serial. So today it is **parallel-full-detect XOR coverage/knee early-stop**, not
both. (Future work: give each bin a per-bin coverage/knee target — local fractions compose to roughly
the global one — to recover early-stop under binning.)

### Measured results (24 logical / 12 physical cores, `COVERAGE_TARGET=1.0`, `DETECT_ONLY`)

**10-min file** (`04-17-23_CA_Tryp_HCD_10min.raw`) — only 5 bins (3 even, 2 odd) because the run is
short, so parallelism is capped well below core count:

| | detect time | features |
|---|---:|---:|
| serial | 13.40 s | 142,958 |
| parallel | 7.72 s | 142,981 |

**1.74×** here is the floor, not the ceiling — 5 bins can't use 24 threads.

**65-min file** (`HFX_MB_14751_5_02062022.raw`, 16.8M peaks) — **20 bins (10 even, 10 odd)**, the real
target:

| | detect time | features |
|---|---:|---:|
| serial | 220.73 s | 972,322 |
| parallel | **46.87 s** | 972,277 |

**4.71× on detect** (20 bins across 24 threads; feature counts within 0.005%). End-to-end wall incl.
the shared 18 s read/index: 243.6 s → 70.0 s (**3.5×**). The speedup is ~4.7× rather than the ~10×
the 10-bins-per-phase geometry suggests, because (1) each phase's wall is bounded by its *slowest*
bin, and the 4-min floor forces dense-region bins to carry more than an equal share, and (2) the two
phases run in series with a barrier between (`slowest_even + slowest_odd`, not their max). Squarely in
the doc's predicted 4-8× band for intra-file parallelism.

(The serial 220.73 s here is well under the 732 s recorded higher up in this doc: that earlier
measurement predates the trace-first persistence pre-gate, which retires most seeds on a cheap
charge-independent XIC check *before* the six-charge scoring — 6.3M `score_hypothesis` calls now vs
24.5M then. The parallel/serial pair above is same-settings, so the 4.71× is apples-to-apples.)

### Parity: the boundary deviation, quantified

On the 10-min file, a fuzzy join (charge + mass ≤ 5 ppm + apex RT ≤ 0.02 min) of serial vs parallel
detected features:

- **142,283 identical** (99.04%).
- 675 serial-only, 698 parallel-only → **1,373 features differ (0.96%)**, balanced add/drop.
- **Summed intensity identical to 0.0014%.**

The balance and intensity conservation confirm the mechanism is pure boundary *reassignment* (the even
phase claims a boundary peak before the odd phase regardless of which side is taller), **not**
fragmentation — a claim-sharing bug would inflate parallel counts and intensity, which is not observed.
99% of features are untouched; only those within `reach` of a bin edge can differ.

## IMPLEMENTED: intra-file 2-D m/z×RT tiling (generalizes the red-black RT path)

Opt-in via **`DETECT_TILE2D`** (takes precedence over `DETECT_PARALLEL`); see
`trace_kernel.rs::detect_features_tile2d` and the spec `agent_info/Detector-2D-Tiling-Spec.md`.
Generalizes the 1-D red-black scheme (2 colors, RT only) to a 2-D grid, 4-colored so non-adjacent tiles
run concurrently. The serial path stays the default/reference and byte-identical; the 1-D path is
untouched. The hot scorers keep their single-`HashSet` signatures (zero serial-path cost).

### Why 2-D
The RT claim reach grows with the trace cap (forcing coarse RT slices), but the m/z claim reach is
fixed by the isotope model. Tiling the m/z axis exposes independent regions exactly when RT slices go
coarse. It also simply yields more tiles: on the 2-hr file a 24-RT-band × 12-m/z-column grid = **288
tiles** vs the 1-D path's 20 bins, so each phase's "slowest tile" tail is much smaller.

### What's new over 1-D
- **Grid + 4-coloring.** Equal-work RT bands (≥ 4-min floor, the 1-D walk) × fixed 96-Th m/z columns;
  `color = (i_rt&1, i_mz&1)`, colors run in series 0→3, tiles within a color in parallel.
- **Strip seeding from king-neighbours.** Each tile pre-seeds its `claimed` from the earlier-colour
  king-neighbours' claims in a padded box — the 2-D analogue of the 1-D odd-bin strip.
- **Seed re-anchoring (§6a) — the m/z-specific correctness piece.** RT tiling never splits an envelope
  (teeth share a scan), but an m/z border cuts *through* one: a tile's local-tallest seed can be a minor
  tooth whose true apex is across the border, which would anchor a mis-placed comb (wrong charge/mono →
  bogus feature). Before seeding a border seed the detector searches an m/z halo at that scan for the
  tallest unclaimed peak and processes it first — enforcing "the taller peak goes first" locally, as
  serial does globally. Backed by the new `PeakIndexingEngine::tallest_unclaimed_in_mz_at_scan`.
- **Collision detector.** The per-colour barrier merge into `global_claimed` *is* the detector: a
  duplicate insert (impossible in margin-respecting mode) is counted and the loser dropped, with a
  warning. It is the correctness canary — a nonzero count means the margin was mis-set.

### Reach geometry (two decoupled m/z quantities — a spec correction)
- **Halo radius / border threshold (`reach_mz`, default 4 Th).** The envelope *span* — how far a minor
  tooth sits from its apex. Measured on the IonStar 2-hr file (1.1 M detections), the per-feature m/z
  reach `= observed_isotopes × (C13−C12)/z` is bounded in **thomson and shrinks with charge** (isotope
  count grows with mass but the `1.0033/z` spacing shrinks faster): z=1 is worst (median 3 Th, max 7 Th),
  everything heavier tighter. Overall p99 = 4.0 Th, p99.9 = 5.0 Th; only ~2 % exceed 4 Th. So 4 Th covers
  ~98 % of straddles directly; the rare tail is a mis-anchor backstopped by the collision detector, not a
  lost/double-claimed peak.
- **Claim bound (`claim_reach_mz` = `max_isotopes × (C13−C12)/min_charge` ≈ 12 Th).** The hard bound on
  how far *any* comb reaches — used only for the strip pad (`reach_mz + claim_reach_mz`) and the column
  floor (`2 ×` the pad), so a neighbour always strips a long z=1 comb's tail. Kept conservative because
  it costs nothing at 96-Th columns.

The first cut conflated these (a single 6 Th for halo + pad + floor) and logged **67 collisions** on the
10-min file: a z=1 comb reaches ~12 Th and a re-anchored claim reaches `halo + claim`, both past a 6-Th
pad. Decoupling drove collisions to **0** while keeping the halo cheap.

### Halo sizing (why 4 Th, not 6)
Halo scan cost grows ~radius² (scan width × border-seed fraction). Sweeping `DETECT_TILE2D_REACH_MZ`,
with feature counts flat to < 0.1 %:

| halo | 10-min detect | 2-hr detect |
|---:|---:|---:|
| 12 Th | 12.9 s | — |
| 6 Th | 7.05 s | 67.2 s |
| **4 Th (default)** | **5.90 s** | **63.0 s** |
| 2 Th | 5.55 s | 62.2 s |

Most of the win is at 6→4; below 4 it flattens. 4 Th is the default.

### Measured results (24 logical / 12 physical cores, `COVERAGE_TARGET=1.0`, `DETECT_ONLY`, halo 4 Th)

**10-min file** (`04-17-23_CA_Tryp_HCD_10min.raw`) — 5×12 = 60 tiles:

| | detect time | features | ΣTIC | collisions |
|---|---:|---:|---:|---:|
| serial | 18.72 s | 142,958 | 90.0% | — |
| 2-D tiling | 5.90 s (**3.2×**) | 143,132 | 90.0% | 0 |

**2-hr IonStar file** (`B03_19_..._2hrs_30B_9B.raw`, 10,582 MS1 scans, 22.5M peaks) — 24×12 = 288
tiles, the real target:

| | detect time | features | ΣTIC | collisions |
|---|---:|---:|---:|---:|
| serial | 413.74 s | 1,113,403 | 96.5% | — |
| 2-D tiling | **63.03 s (6.6×)** | 1,116,158 | 96.5% | 0 |

6.6× on detect — top of the doc's predicted 4-8× band, above the 1-D path's 4.71× (measured on a
*different* 65-min file, HFX_MB), because the 288-tile grid parallelizes far better than 20 RT bins.
Baselines are not cross-file comparable: this 413.7 s is B03_19 at full coverage, not the 1-D section's
220.7 s on HFX_MB.

### Parity: boundary deviation
Fuzzy join (charge + mass ≤ 10 ppm + apex RT ≤ 0.02 min), serial vs 2-D detected features:
- **10-min:** 97.4 % match, 2.5 % boundary deviation, ΣTIC to 0.015 %.
- **2-hr:** 94.8 % match, 5.2 % boundary deviation, ΣTIC to **0.0045 %**.

Higher deviation than the 1-D scheme's ~1 % is expected — 288 tiles have far more edge length than 20
bins — but ΣTIC conservation to ≤ 0.015 % confirms it is boundary *re-fragmentation*, not lost or
double-claimed signal (a claim-sharing bug would inflate counts + intensity, and collisions are 0).

## How to reproduce

```
cargo build --release --example detect_features_tsv
# sub-stage attribution:
DETECT_PROFILE=1 ./target/release/examples/detect_features_tsv <raw> <out.tsv>
# stage timings are always printed + appended to <out>.log

# 1-D red-black RT binning (needs full coverage — global-stop heuristics disable it):
COVERAGE_TARGET=1.0 DETECT_ONLY=1 DETECT_PROFILE=1 DETECT_PARALLEL=1 \
  ./target/release/examples/detect_features_tsv <raw> <out.tsv>
# 2-D m/z×RT tiling (default halo 4 Th; override with DETECT_TILE2D_REACH_MZ / DETECT_TILE2D_MZ_WIDTH):
COVERAGE_TARGET=1.0 DETECT_ONLY=1 DETECT_PROFILE=1 DETECT_TILE2D=1 \
  ./target/release/examples/detect_features_tsv <raw> <out.tsv>
```
