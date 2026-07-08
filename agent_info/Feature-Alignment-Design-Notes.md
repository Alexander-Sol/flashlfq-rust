# Cross-Run Feature Alignment — Design Brainstorm

**Status:** brainstorm / not implemented. Captures the design thinking for the future
"align features across two or more runs" step. Anchored to the current data model
(`ResolvedFeature`) and to machinery that already exists in the repo.

**Date:** 2026-07-07.

## Problem statement

Given the untargeted feature lists of N runs, build **consensus groups** where each group
holds at most one feature per run and every member is the same analyte/peak across runs.
The user framed it as a bipartite (2 runs) / multi-partite (N runs) matching problem. That
framing is correct, but the structure of the data makes it far cheaper than general graph
matching — see the reframe below.

### Two distinct things, don't conflate them

- **MBR (already implemented, `mbr_search.rs`)** — identification-driven donor→acceptor
  targeted transfer. Fills gaps where a run lacks an ID.
- **Feature correspondence / alignment (this note, the future step)** — ID-agnostic grouping
  of `ResolvedFeature`s across runs. This is what OpenMS calls *feature linking*, MaxQuant the
  *match/assembly*, IonQuant/mzMine *correspondence*.

They **compose**, they don't compete: align untargeted features to get the correspondence
backbone + quant matrix, then MBR/ID-transfer to *label* groups and fill holes. Identified
features are also the best warp anchors, so the two share the RT-calibration layer.

## What we're aligning

`feature_refinement::ResolvedFeature`:

```
monoisotopic_mass   // ppm-precise; observed-vs-observed m/z is ~0 ppm median (per project memory)
charge_states       // distinct charges, ascending
apex_rt             // minutes; drifts nonlinearly run-to-run
start_rt, end_rt    // elution bounds
summed_intensity    // weak alignment axis; this IS the per-run quant value
cross_charge_support// confidence proxy — good anchor selector
members             // constituent RefinedFeatures
```

Alignment axes, by usefulness: **neutral mass** (near-exact key) > **warped RT** (the hard
axis) > charge (consistency/tiebreak) > isotope/elution shape (crowded-region tiebreak) >
intensity (weak).

## The reframe that makes it efficient

Do **not** model this as general multi-partite graph matching. It is a graph problem, but one
key almost fully determines the graph.

**Neutral mass is a blocking key, not a similarity dimension.** Same-analyte features across
runs agree in mass to a few ppm. Sort all features from all runs by mass once and sweep a
ppm-width window → **O(N log N)**, and each feature only ever sees the handful of others in its
ppm block. Instead of O(N²) all-pairs, or O(runs²·n²) pairwise-then-stitch.

This is the **same shape as the existing bucketed `feature_refinement::group_features`**
(RT-bin × binary-searched mass window). Generalize it to carry a `run_id` and enforce
one-member-per-run.

Everything hard then happens *inside* a tiny mass block, where only RT (and charge)
disambiguate.

## The three real subproblems

### 1. RT alignment (the actually-hard part)

Raw `apex_rt` isn't comparable across runs (nonlinear drift). Fit a **monotone warp per run**
into a common frame before matching.

- Anchors: unambiguous, high-`cross_charge_support`, high-intensity features that mass-match
  cleanly across runs (plus MS2-identified features when available).
- Reuse `mbr_search`/`mbr` **`get_rt_cal_spline`** and the anchor-peptide RT calibration
  already written for MBR.
- It's an estimate-then-match loop (EM-ish): better correspondences → better anchors → better
  warp. One or two iterations help.

### 2. Matching within a mass block

- **Two runs:** bipartite matching per block. Full Hungarian (O(n³)) is overkill — after
  mass+RT constraint the candidate set is tiny. Use **mutual-nearest-neighbor** (A's best is B
  AND B's best is A) on a distance combining Δmass_ppm + Δwarped_RT (+ optional isotope-shape).
  Near-optimal and near-linear here.
- **N runs:** do **NOT** do all pairwise alignments and stitch by transitive closure — that is
  exactly the megagroup bug (see below). Use **reference-based (star) alignment**: build/pick
  one consensus reference frame (pooled pseudo-run, or the richest run), align every run *to the
  reference only*, features landing in the same reference slot form a group. O(N·runs), linear,
  consistent by construction.

### 3. Enforce "≤1 feature per run per group" STRUCTURALLY

A consensus group is a **partial matching, not a connected component**. If you build edges and
union-find them, you rebuild the transitivity trap. Formulate as **assignment** (each run
contributes at most one node per group) so the constraint cannot be violated.

> **Lesson already learned in this repo:** intra-run charge consensus originally used transitive
> single-linkage on RT-range *overlap* and daisy-chained distinct peptides into megagroups (up
> to 1042 members, 11 min wide), erasing real peptides — mass recall 36.5%→82.3% once the overlap
> arm was removed and linkage required apex proximity only. Cross-run alignment has the identical
> hazard one dimension up. Carry the fix forward: constrain structurally, never union edges.

## Recommended architecture

1. **Warp** — fit a monotone RT spline per run to a common frame from cross-run mass-matched
   anchors (reuse `get_rt_cal_spline`). Align on neutral mass (charge-agnostic).
2. **Block** — sort all features (all runs) by mass; two-pointer sweep a ppm window → mass
   blocks. O(N log N).
3. **Match inside each block** in (Δmass_ppm, Δwarped_RT) space — greedy mutual-nearest with the
   one-per-run constraint, against a single reference frame for N>2. Charge agreement as a hard
   filter or strong tiebreak; isotope/apex-shape as secondary distance for crowded blocks.
4. **Emit** consensus groups carrying per-run `summed_intensity` (the quant matrix) + which runs
   are missing (MBR / ID-transfer candidates).
5. **Optionally iterate** warp↔match once or twice.

Total complexity ~**O(N log N)** dominated by the sort.

## Two things to steal from the existing MBR code

- **Decoy-based FDR on the alignment.** `mbr_search::get_random_peak` already draws mass-shifted
  decoys (5–11 H away, RT far). Do the same here: a *decoy match* is a mass-incompatible pairing
  at the same RT proximity. The rate at which decoys pass the distance threshold estimates the
  false-alignment rate and **sets** the acceptance threshold — instead of hand-tuning ppm/RT
  windows. Principled answer to "how close is close enough."
- **Shared RT-calibration layer.** Identified features are the best warp anchors; the alignment
  and MBR both consume `get_rt_cal_spline`.

## The one fork worth an explicit decision

Multi-run strategy:

| Strategy | Pros | Cons | Verdict |
|---|---|---|---|
| **Reference / star** | simplest, consistent, linear | needs a good reference frame | **recommended start** |
| Incremental consensus clustering | cheap memory, streaming | order-sensitive | fallback |
| Pairwise + graph consensus | most flexible | transitivity hazard, O(runs²) | avoid unless forced |

Start with **star alignment**; reach for the others only if a specific failure mode demands it.

## Next concrete steps (when picked up)

1. Prototype the **decoy-FDR threshold estimator** first — it de-risks all the window-tuning.
2. Sketch an `align_features` module: block-sweep + one-per-run matcher against `ResolvedFeature`,
   reusing `get_rt_cal_spline`. Generalize `group_features` to carry `run_id`.
3. Build the warp↔match iteration and validate on the two-K562-file corpus already used by Phase 3.
