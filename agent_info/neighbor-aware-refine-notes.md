# Neighbor-aware refinement — exploration notes

Branch: `neighbor-aware-refine` (off `obo-envelope-probe` @133e814).
Goal: incorporate already detected/refined neighbor features into refinement, esp. to help **low-score** features (contaminated / crowded-region envelopes). Motivating case: YAAELHLVHWNTK z3 (mono pushed +1 by refine; an M-1 interferer at 547.252 belongs to a strong co-eluting z2 neighbor, 1093.51/z2, itself likely off-by-one).

Fitness: CA/Lumos cov90 benchmark. Baseline at branch point = **602/620 recall, 584 charge-matched**. Every experiment behind an env flag; default path must stay byte-identical (excluded/neighbors empty => same result).

## Existing machinery
- `NeighborIndex` (feature_refinement.rs): RT-bucketed; `forbidden_positions(self_idx, win_lo, win_hi)` = ascending isotope m/z of co-eluting *other* features, k∈[0,kmax] (mono upward only — MISSES a neighbor's M-1 / off-by-one mono).
- `shift_decon_gated` / `four_way_decon_gated`: use forbidden positions to exclude neighbor peaks from ANCHOR selection only. Behind `NEIGHBOR_AWARE=1`, four-way comparator only — NOT in default `refine_feature_shift`.
- Default refine = `refine_feature_shift` (detector anchor, best_charge_by_fit, walkback, envelope_fit_cosine). No neighbor awareness.

## Gaps (from YAA analysis)
1. Neighbor machinery not wired into default refine.
2. Gates ANCHOR, not the FIT / walk-back / charge-selection scoring.
3. forbidden_positions is mono-upward only; misses neighbor M-1 / mis-placed-neighbor mono.
4. Chicken-and-egg: neighbors themselves off-by-one => want strongest-first, refined-neighbor context.

## Approaches to try (measure each)
- **A. Neighbor-masked fit.** `envelope_fit_cosine` drops observed window peaks that are neighbor-owned (on a neighbor grid, off own grid) from the "unexplained" penalty. Thread through best_charge_by_fit (charge selection — should subsume/strengthen anti-doubling: masked interferer teeth => 2z stops winning), walkback (won't grab interferer), and decon_score. Flag `NEIGHBOR_REFINE`.
- **B. Neighbor-gated anchor** in default (reuse shift_decon_gated). Helps "distant grab"; less for YAA (own-peak anchor).
- **C. Strongest-first + refined-neighbor context.** Sort features by intensity; refine in order; build neighbor context from already-REFINED features. Targets low-score features vs corrected neighbors.
- **D. Extend neighbor grid to ±1 mono ambiguity** (cover a neighbor's M-1 / off-by-one). Needed for YAA's 547.252.
- **E. Grid-aware down-weight of shared peaks** (competitive split) instead of hard mask — the naive delete-claimed censoring regressed; down-weight may be safer.

## Log
- (start) baseline 602/584; scaffolding masked fit first.
- **Exp A1: neighbor-masked fit, detected neighbors, ALL neighbors, hard mask.** `NEIGHBOR_REFINE=1`. Result **599/575 — REGRESSED** (-3 recall, -9 charge). Cause (hypothesised): detected neighbors are themselves off-by-one/mis-charged (chicken-egg); masking their teeth removes a feature's own real peaks and makes the correct charge fit worse. Mirrors the earlier naive-censoring regression. Plumbing kept (envelope_fit_cosine_masked + refine_feature_shift_neighbor + NEIGHBOR_REFINE), invariant when mask empty (default still 602/584).
- **Exp A2: mask only STRONGER neighbors** (`NEIGHBOR_MIN_RATIO`). Sweep (detected neighbors, hard mask):
  - all (0): 599/575 · 2×: 601/582 · **5×: 603/585** · 8×: 603/585 · 12×: 603/585.
  - **Sweet spot ratio ≥5× → 603/585 (+1 recall, +1 charge over baseline 602/584).** Small but real net win: deferring only to *dominant* (≥5×) co-eluting neighbours helps a few low-score features without the collateral of masking moderate neighbours. Confirms the "weak feature in a strong neighbour's shadow" hypothesis. Plateau (5×=8×=12×) => the useful masks come from a handful of very-dominant neighbours.
  - Note: existing `RECHARGE_MIN_FIT` floor already covers the harmonic-doubling cases, so masking's charge win here is incremental, not the doubling class.
- **Exp A3: refined-neighbor context (two-pass).** `NEIGHBOR_REFINE=refined` (refine once, rebuild index from corrected placements, refine again masking against those). At 5×: **603/585 — identical to detected single-pass**, at ~2× the refine cost. So corrected neighbour grids add nothing over raw detected here (a ≥5× dominant neighbour's placement doesn't move enough to change which peaks get masked). Detected single-pass is the keeper.

## Result so far
**Neighbour-masked fit, detected context, `NEIGHBOR_MIN_RATIO=5` (default): 603/585 vs baseline 602/584 — clean +1 recall, +1 charge, 0 lost.** Gained peptide: TLNFNAEGEPELLMLANWRPAQPLK (2905 Da, z4, RT 16.34 — a heavy peptide where a dominant co-eluting neighbour was fooling the fit). Flag-gated (`NEIGHBOR_REFINE`); default path byte-identical (mask empty). 205 lib tests pass.

Assessment: **modest but real.** The neighbour-masking mostly overlaps failure modes the shipped fixes (co-elution knitting, RECHARGE_MIN_FIT floor, fit tiebreak) already cover, so the marginal gain is small on this dataset. Not promoting to default on a +1/+1 single-dataset result; left as an opt-in experiment.

- **Exp A4: low-baseline-fit gating.** Mask only features whose unmasked fit < gate (`NEIGHBOR_FIT_GATE`). Gate 0.9 at ratio 2 → 601/582 (== ungated ratio 2); at ratio 5 → 603/585 (== ungated). **No-op** — the masking already only changes low-fit features, so gating on own-fit excludes nothing. The regression at low ratio is driven by NEIGHBOR STRENGTH (masking moderate 2–5× neighbours), not the feature's own fit, so ratio is the right knob and 5× is the ceiling. Reverted the gate (it doubled the flag-path cost for no gain).

- **Exp B (user-requested): iterative score-ordered masking.** `NEIGHBOR_REFINE=iterative`: refine in descending pass-1 score; each feature that clears `NEIGHBOR_LOCK_MIN_SCORE` locks its corrected grid (`LockedGrids`, incremental RT-bucketed) so later lower-scoring features mask its peaks.
  - lock ALL (threshold 0): 600/580 — regressed (a "higher-scoring" but mediocre feature is still uncertain; masking against it adds collateral, same as mask-all).
  - lock ≥0.9: 601/582 · lock ≥0.95: 601/583. Asymptotes to baseline as the gate tightens; **never beats baseline 602/584, and below the intensity-dominance mask (603/585).**
  - **Lesson:** for masking, "who owns this shared peak" is answered by **intensity dominance**, not refinement confidence/order. A confident-but-weak neighbour doesn't own a peak the current feature might legitimately share. Intensity-ratio single-pass (≥5×) remains the masking winner.

## Strategy 2 (user-requested): joint linear model — `src/joint_fit.rs` + `JOINT_FIT`
New module: `nnls` (Lawson–Hanson active-set NNLS), `joint_envelope_fit` (models a window as a
non-negative sum of averagine envelopes over the union grid, returns amplitudes + cosine),
`joint_fit_target_shift` (searches the target mono over ¹³C shifts with neighbours fixed). 4 unit tests
pass (recovers amplitudes of overlapping envelopes 3 Da apart; shift search recovers a mis-seeded mono).
Wired as a post-refine pass (`apply_joint_fit_pass`): for each refined feature scoring < `JOINT_FIT_MAX_SCORE`,
model its window as itself + up to 3 overlapping co-eluting neighbours, shift-search the target, accept a
move only if the target's solo re-score improves by `JOINT_FIT_MIN_GAIN`.

Results (CA/Lumos cov90):
- max_score 0.7, min_gain 0.02: moved 5473, **602/584 — net zero, 0 reference peptides changed.**
- max_score 0.9, min_gain 0.02: moved 5726, 602/584 — **+1 gained (AVVQDPALKPLALVYGEATSRR z4, heavy/chimeric) / −1 lost (YENEVALR z2, over-corrected).**
- max_score 0.9, min_gain 0.2: moved 1234, 602/584 — **0/0** (strict gate rejects both the gain and the loss; they're in the same marginal-improvement regime, so min_gain can't separate them).

**Conclusion: net-zero on recall.** The primitive is correct and genuinely separates overlapping envelopes, but (a) the low-score features it moves are overwhelmingly non-reference (junk), and (b) the residual reference misses (e.g. YAA) are NOT fixable by a mono-shift — at their apex the true-mono placement's own fit is inherently poor (contaminated envelope), so no target-fit-driven method prefers it. The joint model separates signal; it cannot manufacture fit where the evidence is weak. **Kept as a validated primitive** (the foundation for advanced multi-envelope decon / overlapping-feature quantification), opt-in, default unchanged.

## Overall assessment of neighbour-aware refinement (this session)
Neighbour information gives at best a **marginal** recall improvement on this dataset (the intensity-dominance mask, +1/+1 → 603/585), because the shipped fixes (co-elution knitting, RECHARGE_MIN_FIT, fit tiebreak, walk-back) already capture most neighbour-driven failures. The remaining ~18 misses are dominated by cases neighbour-awareness can't reach: inherently poor target evidence (contaminated apex), detection-depth, or FlashLFQ being wrong on heavy modforms. The joint-fit primitive is the most reusable artifact (advanced decon foundation).

## Approaches NOT yet tried
- **Down-weight vs hard-mask** shared peaks (competitive split) — softer than the invisible mask.
- **Joint fit at the QUANTIFICATION stage** (split shared-peak intensity between overlapping features) rather than as a mono-corrector — likely the joint model's real value.
- **Off-by-one via neighbour coherence** for the sub-walk-back-gate class (YAA), accepting that some need a second (non-apex) scan where the true mono is cleaner.
