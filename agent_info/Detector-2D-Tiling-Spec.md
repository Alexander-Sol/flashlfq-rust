# Detector 2D m/z×RT Tiling — Implementation Spec

Parallel intra-file detection by tessellating the m/z×RT plane into tiles, 4-coloring them so
non-adjacent tiles run concurrently, and processing tiles in intensity order. Generalizes the shipped
1-D red-black RT-binning path (`trace_kernel.rs::detect_features_parallel`, 2 colors) to 2 dimensions
and 4 colors, and adds a collision detector as a hard safety net. Opt-in, serial path unchanged.

Motivation over 1-D: the RT claim reach can grow (raising `trace_max_half_width_minutes` to capture
tailing peptides), which forces coarse RT slices. The **m/z claim reach is fixed by the isotope model
and does not grow**, so tiling the second axis keeps exposing independent regions exactly when RT
slices go coarse — restoring parallelism and load balance.

## 1. Claim geometry (what can collide)

Two features collide only if they claim the same peak, which requires their claim regions to overlap
in **both** m/z and RT. A feature seeded at `(seed_mz, seed_rt)` reaches:

- **RT:** `reach_rt = max(trace_max_half_width_minutes, rt_half_window_minutes)`. The claim-extent
  trace is capped at `trace_max_half_width_minutes` (nominal 0.5 min); the scoring window reads
  `±rt_half_window_minutes` (= 2σ, ~0.13 min nominal). Take the max so the strip (§5) covers reads too.
  **This is the axis that grows.**
- **m/z:** the theoretical bound is `max_isotopes × (C13_MINUS_C12 / min_charge)` = `12 × 1.0033 ≈
  12 Da` at z=1. **This is fixed by the isotope model — it never grows with the RT trace cap.** But it
  is ~2× conservative: measured on real data (§6a), the **actual claimed-envelope m/z span is ≤ ~6 Da
  (p99 ~4 Da, max ~6 Da) on both the 10-min and 65-min files**, because real peptides don't show all 12
  isotopes at z=1. So **derive `reach_mz` from the data** — e.g. the p99.9 observed span plus a small
  guard, clamped to the theoretical bound (~6 Da typical). `reach_mz` sets three things: the strip
  padding (§5), the re-anchoring halo radius (§6a), and the hard tile-width **floor** (`2 × reach_mz`).
  It is **not** the chosen m/z tile width — §2 deliberately runs m/z tiles *wide* (96 Da) to keep the
  re-anchoring rate low.

## 2. Tiling (equal-work grid)

A clean grid, `n_rt` RT bands × `n_mz` m/z columns:

- **RT bands:** equal seed count in RT (the shipped 1-D equal-work walk), floor width `2 × reach_rt`.
- **m/z columns:** **fixed ~96-Da-wide columns** (target width, not the `2 × reach_mz ≈ 12 Da` floor).
  Wide m/z tiles are a deliberate choice: they minimize the fraction of seeds near an m/z border and
  thus how often re-anchoring (§6a) fires — measured straddle exposure is ~1 % at 96 Da vs ~4 % at
  24 Da. The cost is fewer m/z columns (~12 over a 1200-Da range) hence less m/z parallelism; that is
  acceptable because the RT axis supplies the rest, and it still clears the occupancy target below.
  Cut points are global (same for every RT band). (Equal-work-in-m/z is a v2 option if a fixed grid
  leaves a hot column; v1 uses fixed width for simplicity.)
- **Tiles = grid cells.** Equal-work per tile is approximate (co-elution hot spots are the residual
  imbalance — §8 hot-tile subdivision is the v2 lever).

Tile counts: aim for **~2–4 tiles per thread per color** so rayon work-stealing smooths the tail. With
4 colors and ~22 threads that is ~180–350 total tiles. At 96-Da m/z columns (~12 over 1200 Da) × the RT
bands (~16 at a 4-min floor on a 65-min run) ≈ ~190 tiles ≈ 48/color ≈ 2/thread — clears the target. If
`reach_rt` grows (coarser RT bands) the m/z axis carries proportionally more of the parallelism, which
is why 96 Da (not the 12 Da floor) still leaves enough columns. If a run is too short/narrow to yield
≥ 2 tiles per axis at these widths, fall back to serial (or the 1-D path) and warn.

**Hard invariant:** every tile ≥ `2 × reach` on each axis. Validate at build time; if unmet, coarsen
that axis (fewer tiles) until it holds, else fall back. This is what makes same-color tiles provably
collision-free (§4); violating it is only allowed deliberately in the sub-margin performance mode (§8),
which *relies* on the collision detector.

## 3. Coloring & schedule

Tile `(i_rt, i_mz)` → **color = (i_rt & 1, i_mz & 1)**, four classes {(0,0),(0,1),(1,0),(1,1)}.

A claim reaches into the 8 king-neighbors (a corner seed reaches the diagonal tile), so adjacency is
the 8-neighborhood, which 4-colors as a 2×2 block. Same-color tiles differ by ≥ 2 in at least one
axis → separated by a full tile ≥ `2 × reach` on that axis → their claims cannot overlap.

```
build tiles; assign seeds; sort seeds tallest-first within each tile
global_claimed = {}                     // grows across colors
for color c in [ (0,0), (0,1), (1,0), (1,1) ]:      // 4 phases
    results[c] = par_iter over tiles of color c:
        local = strip of global_claimed near this tile (§5)   // prior colors' claims it can see
        (features, claimed_peaks) = detect_bin(seeds_of_tile, local)   // reuse serial worker
    barrier: merge results[c].claimed_peaks into global_claimed,
             detecting collisions (§6)
concat features from all tiles in (color, tile) order   // deterministic
```

Colors run in series (0→3); each color sees all earlier colors' claims via the strip. Within a color,
tiles are independent and run in parallel. Optionally wrap the 4-color loop in an outer **intensity
tranche** loop (§7) for tighter parity — omit in v1.

## 4. Why the fast path has zero collisions (nominal mode) — and the m/z exception

With margin-respecting tiles and strip seeding: same-color tiles are ≥ `2×reach` apart → never claim
the same peak. Cross-color contention is resolved by processing order (a boundary peak goes to the
earlier color that claims it; later colors exclude it via their strip). So along **RT** borders, in
nominal mode, no collision fires and no re-score runs — the detector is pure safety, and the only
deviation from serial is that boundary peaks are prioritized by color order rather than global
intensity (same class as the shipped 1-D scheme).

**The m/z axis is different, and §6a is mandatory.** Unlike RT (an envelope's teeth share one scan, so
RT tiling never splits an envelope), m/z tiling *cuts through isotope envelopes* — the whole point of
the algorithm is that teeth are spread across m/z. So an m/z border can split an envelope, and a tile's
local-tallest seed may be a *minor* tooth whose true apex sits in the neighbor tile. That breaks the
detector's core assumption (seed = envelope apex) and is not a rare safety event — it is expected at
every m/z border. See §6a.

## 5. Strip seeding (keeps the hot path single-set)

The shipped scorers (`score_hypothesis`, `trace_seed_extent`, `gather_extent_peaks`) take a single
`&HashSet<PeakKey>` and must stay that way — adding a second set to the ~90 % hot path costs ~15–20 %
(an extra `contains` per found peak). So instead of passing `global_claimed` by reference (two-set
check), pre-seed each tile's local set with only the prior-color claims it can observe:

- A tile `(i,j)`'s reads/claims stay within its region padded by `reach`. Relevant prior claims come
  only from its **already-processed (earlier-color) king-neighbors**.
- Strip = those neighbors' claimed peaks lying in the padded box
  `[rt_lo(i) − reach_rt, rt_hi(i) + reach_rt] × [mz_lo(j) − reach_mz, mz_hi(j) + reach_mz]`.
- Neighbor claim lists are small (equal-work) and already in hand; filter by the box. No global spatial
  index needed. Result: each tile runs the unmodified serial worker against one combined `HashSet`.

## 6. Collision detection & resolution (the safety net)

The per-color barrier merges each tile's `claimed_peaks` into `global_claimed` **serially** (single
thread; ~100 ns/insert, ~a few % of detect total — this is part of the serial-acceptance cost). The
merge *is* the detector: a `HashSet::insert` returning `false` means two same-color tiles claimed the
same peak — impossible in nominal mode, possible only under sub-margin tiles (§8) or a mis-set reach.

On collision, the peak belongs to the **higher-intensity seed**'s feature. The loser is handled per a
flag:

- **`drop` (default, v1):** discard the loser feature, increment a logged counter. Correct and cheap;
  fine because nominal-mode collisions are zero. A nonzero count is a red flag that the margin
  invariant was violated.
- **`rescore` (v2, for sub-margin mode):** re-run the loser seed through the serial worker against the
  updated `global_claimed` and re-decide. Restores near-serial fidelity when collisions are frequent
  by design.

## 5a. Claims cross tile borders freely (by design)

A feature's **claim is never clipped at a tile edge.** `trace_seed_extent` / `gather_extent_peaks`
query the global peak index by m/z and scan with no tile restriction, so a feature seeded in tile T
traces its full elution and claims its full isotope envelope wherever the peaks physically sit —
including into the neighbouring RT and/or m/z tile, up to `reach_rt` / `reach_mz`. This is **required**,
not merely tolerated: clipping at the edge would orphan the tail and let it re-seed as a fragment (the
exact failure whole-extent claiming exists to prevent).

Two mechanisms make crossing safe rather than a collision hazard — both already running in the shipped
1-D detector (odd bins strip-seeded from even claims; verified 99 % identical to serial, ΣTIC to
0.0014 %):

- **The `2 × reach` tile-width floor (§4):** a claim reaches ≤ `reach` past a border, never as far as
  the next same-color tile, so no two concurrently-processed tiles claim into the same region.
- **Strip seeding (§5):** the neighbour tile (different colour, different phase) is pre-seeded with the
  claims that crossed into it and excludes them; whichever tile runs first claims, the other defers.

The **only** requirement is that the width floor be enforced against the *actual* `reach_rt` — so if the
trace cap is raised to capture wide tails, the RT band floor rises with it (`2 × reach_rt`). Exceed that
(cap > ½ RT tile) and same-colour RT-neighbours can collide → the collision detector (§6) is the
backstop. (The m/z axis has the additional, seed-anchor wrinkle of §6a; RT crossing itself is clean.)

## 6a. Envelopes straddling m/z borders (seed-anchor correctness) — MANDATORY

**The invariant m/z tiling breaks.** The detector assumes the seed is the envelope's most-abundant
peak: it anchors the comb as `mono = seed − i*·spacing` (i* = most-abundant index of the weight model),
i.e. it *places the seed as the tallest tooth*. In serial this always holds, because global
tallest-first + claim-as-you-go guarantees the apex of any envelope seeds first and claims its lesser
teeth, so a minor tooth is already claimed before its turn and never seeds. RT tiling preserves this
(all teeth share one scan → one RT band). **m/z tiling does not:** if an envelope's apex A sits in tile
`j+1` and a lesser tooth B in tile `j`, and B is the tallest *unclaimed* peak in tile `j` at the time
`j` runs (e.g. `j`'s color precedes `j+1`'s), then B seeds first and anchors a comb on itself — a
**mis-placed monoisotope / wrong charge → a bogus feature**, plus A's real feature is degraded because
its lower teeth were consumed by B.

**Frequency — measured** (`MZ_STRADDLE=1` in `detect_features_tsv`, both files, at full coverage). The
failure mode is a *wrong* feature (mis-anchored comb → wrong charge/mono), not a reassignment, so even
a few percent must be handled. Measured envelope m/z spans ≤ ~6 Da (p99 ~4, median ~0.7–1.0). Fraction
of features a random tile border would split *and* orphan a tooth above the seed floor, by tile width
`W`:

| `W` | exposed | mis-seed rate w/o re-anchoring (≈ ½, color order) |
|----:|--------:|--------------------------------------------------:|
| 24 Da (= 2·reach at 12 Da) | ~4.0–4.5 % | ~2.0–2.2 % |
| 48 Da | ~2.0–2.2 % | ~1.0–1.1 % |
| 96 Da | ~1.0–1.1 % | ~0.5–0.6 % |

Nearly identical on the 10-min and 65-min files (the longer gradient did not raise it). Straddle scales
~`1/W`, so tile width is itself a knob. **v1 chooses `W = 96 Da`** (§2) to sit at the ~1 % row —
minimizing how often re-anchoring must fire — and re-anchoring eliminates that residual ~1 %.

**Handling: seed re-anchoring (local tallest-first repair) — the v1 mechanism.**

Serial never mis-anchors because global tallest-first guarantees a taller peak is resolved before any
shorter peak in its neighborhood. Re-anchoring restores that invariant *locally* across the tile
border. When a seed `B` is about to fire:

1. Search `B`'s m/z halo `[B_mz − reach_mz, B_mz + reach_mz]` at `B`'s scan (± the RT window) for the
   tallest **unclaimed** peak `A` (which may be `B` itself, and may lie across an m/z border).
2. If `A ≠ B`, process `A` as the seed first (anchor the comb's most-abundant tooth on `A`, score,
   claim its envelope) — i.e. swap the anchor from `B` to `A` before the existing scoring runs.
3. Then fall through: if `B` is still unclaimed, process `B` as normal (against the now-updated claims).

This is **unconditionally safe** — it just enforces "the taller peak goes first" across the border,
exactly as serial does globally — so it needs no isotope-consistency gate:

- If `A` is `B`'s true envelope apex → the comb anchors correctly on `A` and claims `B` as a tooth; `B`
  falls through already-claimed. No bogus feature.
- If `A` is a *distinct* co-eluting species → `A` is detected (it's real) and `B` still falls through
  and is detected on its own. No lost feature, no bogus feature. `A` is merely attributed to `B`'s turn;
  when `A`'s own tile runs later it finds `A` claimed and skips it.

**Composes with the 4-coloring for free.** Every king-neighbor of a tile is a *different* color (a 2×2
checkerboard has no same-color adjacency), so the halo peak `A` always lives in a tile that runs in a
different phase — never concurrently with `B`'s tile. So the cross-border read-and-claim cannot race,
and it is correct in both color orders (A's tile earlier → `A` already claimed, seen via the strip; A's
tile later → `B` claims `A` first, A's tile then skips it). The collision detector (§6) therefore stays
a pure backstop with `drop` resolution — re-anchoring *prevents* straddle collisions rather than
resolving them afterward, so the heavier `rescore`/eviction path is not needed in v1.

**Scope & cost.** Apply only to seeds within `reach_mz` of an m/z border — interior seeds have their
whole envelope in-tile, where within-tile tallest-first already guarantees apex-first, so they skip the
halo scan. With 96-Da m/z tiles (§2) that border zone is `2·reach_mz / 96 ≈ 12 %` of seeds, and only
~1 % actually straddle (§6a table) — so the halo scan is rare. It is one bounded m/z-range lookup at
one scan, cheap against the 6× charge scoring.

**Search both sides, but the apex is usually to the *lighter* side.** The isotope envelope is
right-skewed — there are typically more teeth *above* the most-abundant peak (mono, +1, +2 … tailing to
higher mass) than below it — so a randomly-orphaned minor tooth is more often a *heavy* tooth, whose
apex sits at *lower* m/z. The halo covers `±reach_mz` (both directions) for correctness; the asymmetry
is only an ordering/early-exit hint (scan the lighter side first). Do **not** restrict to one side —
a light minor tooth (e.g. the monoisotope when the apex is +1/+2) needs the heavier direction.

**Residual edge:** cascading (`A` may itself have a taller peak just outside `B`'s halo but inside
`A`'s) — rare at this reach; one level of re-anchoring is a strong approximation and the collision
detector backstops the rest.

**Why not intensity tranching (rejected).** Tranching does *not* help the straddle: isotope ratios are
linear factors (~1.5–2× between an apex and its adjacent tooth, ~0.3 decades), while intensities span
~4–6 decades that a handful of tranches cut into ~1-decade bands. An apex and its neighbor tooth almost
always land in the *same* tranche, so tranche cuts never fall between them and give the apex no head
start over its own teeth. See §7.

## 7. Intensity tranches (optional, OFF by default)

Wrap the 4-color loop in an outer loop over descending intensity floors (equal-seed-count tranches):
process the tallest tranche across all 4 colors, then the next, etc. This approximates global
tallest-first, marginally shrinking the boundary deviation (which distinct co-eluting feature claims a
shared peak first), at the cost of `4 × n_tranches` barriers. Tranche depth does **not** change the
serial-acceptance fraction — it is purely a parity knob.

**It does not help the §6a envelope straddle** — isotope ratios are linear O(1) factors while
intensity spans several decades, so an apex and its adjacent tooth land in the same tranche and tranche
cuts never fall between them (see §6a). Re-anchoring (§6a) is what handles straddle. Tranching is left
**off in v1**; enable only if the residual boundary deviation between distinct features proves too large
for a downstream consumer.

## 8. Sub-margin performance mode (v2)

Deliberately use tiles **smaller** than `2×reach` for finer load balance / more parallelism, accepting
frequent boundary collisions caught by the detector with `rescore`. This trades re-score cost for
granularity and is the escape hatch when a grown reach makes margin-respecting tiles too coarse. Gated
behind an explicit flag; requires `rescore` resolution and benchmarking of the collision rate. Not in
v1.

## 9. Determinism, fallback, incompatibilities

- **Deterministic:** tiles are serial internally; results reassembled in fixed (color, tile) order;
  rayon `collect` preserves order. Identical run-to-run.
- **Not bit-identical to serial:** boundary deviation as in §4 (expect somewhat more than the 1-D
  scheme's ~1 %, since more tiles = more edge length; measure it).
- **Fallback to serial** when: run too short/narrow for ≥ 2 margin-respecting tiles per axis; a
  global-stop heuristic is engaged (coverage target < 1, knee, reject-stop — same incompatibility as
  the 1-D path, since binning breaks the global tallest-first traversal those need).

## 10. Data structures & flow

1. `seeds = engine.all_peaks()`.
2. `reach_rt`, `reach_mz` from params (§1). RT band + m/z column cut points (§2); validate margins.
3. Per seed: `(i_rt, i_mz)` from cut points; `color`. Sort a seed-index permutation by
   `(color, tile, intensity desc)`; materialize one reordered `Vec` so each tile is a contiguous
   slice (as the shipped 1-D path does — one seed pool resident, low memory).
4. Per-tile bounds `[rt_lo,rt_hi] × [mz_lo,mz_hi]`; king-neighbor index (grid arithmetic).
5. Run the §3 schedule. Per-tile output `(Vec<DetectedFeature>, Vec<claimed peak+intensity>)`.
6. `global_claimed: HashSet<PeakKey>` accumulator; strip seeding from earlier-color neighbors (§5);
   barrier merge + collision check (§6).
7. Concatenate features; return.

Reuses the existing `detect_bin` worker (serial per-tile loop) with **one addition**: the re-anchoring
step (§6a) for m/z-border seeds. Gate re-anchoring so it is inert on the serial path (e.g. a
`reanchor_halo: Option<f64>` param, `None` in serial) — the serial reference must stay byte-identical.
New code is: 2-D grid build + margin validation, the re-anchoring hook, the 4-color strip-seeded
schedule, and the barrier-merge collision detector.

## 13. Implementation handoff checklist

Build order (each step independently testable against the serial reference):

1. **Reach + grid.** Compute `reach_rt` (§1), `reach_mz` data-derived (§1, one pass over detected
   envelope spans or a fixed ~6 Da to start). Build RT bands (reuse the shipped 1-D equal-work walk)
   and fixed 96-Da m/z columns (§2). Validate every tile ≥ `2·reach` per axis; else coarsen / fall
   back (§9). Assign each seed `(i_rt, i_mz, color)`; sort a permutation by `(color, tile, −intensity)`;
   materialize one reordered seed `Vec` with contiguous per-tile slices (mirror the 1-D path).
2. **Re-anchoring in the worker (§6a).** Add the halo search + anchor-swap + fall-through to
   `detect_bin`, gated off for serial. Unit-test: a synthetic envelope split across an m/z border seeds
   the same feature regardless of which side's tile runs first.
3. **4-color strip-seeded schedule (§3, §5).** Colours `(i_rt&1, i_mz&1)` run 0→3; each tile
   strip-seeds its local `claimed` from earlier-colour king-neighbours' claims within the padded box;
   `par_iter` across a colour's tiles; barrier merges into `global_claimed`.
4. **Collision detector (§6).** The serial barrier-merge doubles as the detector (`insert` == false →
   collision); `drop` resolution + a logged count. In nominal (margin-respecting, re-anchored) mode the
   count must be ~0 — that count is the correctness canary.
5. **Assemble + fallbacks (§9).** Concatenate features in (color, tile) order; serial fallback on
   short/narrow runs and when a global-stop heuristic is on.
6. **Validate** exactly as the 1-D path was: on the 65-min file, fuzzy-join vs serial (expect ~1–few %
   boundary deviation, ΣTIC preserved to ~1e-3 %, collision count ~0) and record detect speedup +
   per-phase load balance. Gate the whole path behind an env flag (e.g. `DETECT_TILE2D`), serial
   default unchanged.

**Do not touch** the hot scorers' single-`HashSet` signatures (`score_hypothesis`, `trace_seed_extent`,
`gather_extent_peaks`) — strip seeding exists specifically to avoid a two-set check on the ~90 % path.
Keep the shipped 1-D `detect_features_parallel` as-is; this is a separate opt-in path.

## 11. Open knobs (defaults chosen; confirm before build)

| Knob | Default | Notes |
|---|---|---|
| m/z tile width | **96 Da** (fixed) | wide on purpose → straddle ~1 %, minimal re-anchoring; §2 |
| RT band width | equal-work, floor `2·reach_rt` | shipped 1-D walk; §2 |
| tiles per thread per color | ~2–4 | ~190 tiles at 96-Da m/z × ~16 RT bands, 4 colors; §2 |
| seed re-anchoring | ON, m/z-border seeds, **both sides** | local tallest-first repair; apex usually to lighter side; §6a |
| tranches | off | doesn't help straddle; marginal parity only; §7 |
| collision resolution | `drop` (backstop) | re-anchoring prevents straddle collisions; `rescore` only for sub-margin mode (§8) |
| reach_rt | `max(trace cap, 2σ)` | grows with the trace cap |
| reach_mz | data-derived (p99.9 span + guard, ~6 Da), clamp to `max_isotopes·1.0033/min_charge` | sets strip pad, halo radius, tile floor — **not** the tile width |
| tile min width (floor) | `2 × reach` per axis | hard invariant; §2/§4 |

## 12. Effort

Moderate — ~1–2 focused sessions. Grid build + validation (small), re-anchoring hook in the worker
(small–moderate, §6a), 4-color strip-seeded schedule (moderate, mirrors shipped 1-D), collision
detector (small, dormant in nominal mode), fallbacks (small), plus end-to-end parity + speedup
measurement on the 65-min file (as done for the 1-D path). See §13 for the build order.
