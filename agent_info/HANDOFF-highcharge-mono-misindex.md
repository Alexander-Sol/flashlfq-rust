# Hand-off: high-charge monoisotope mis-indexing locked in by cross-charge merge (misses #7, #15)

**Status:** open. Diagnosed, not fixed. Raw data in [`misses-5-7-15-data.md`](misses-5-7-15-data.md).
**Branch context:** `neighbor-aware-refine`; benchmark baseline after the anti-halving fix (`fdac6b9`) = **604/620 recall, 585 charge-matched** on CA/Lumos cov90.

## One-paragraph summary
For some peptides the per-charge deconvolution places the monoisotope **one ¹³C off the true mono onto an isotope position that has no observed peak**, even though the detector's most-abundant peak sits *on* the true mono. `resolve_charge_state_consensus` then merges the charge states and, voting by support/fit/intensity, **locks in the mis-indexed charge and overrides the charge(s) that had the mono right**. The reported mono is off by one and the reference peak is missed. This is a genuine bug (not a de-novo limitation) for peptides whose true mono is a strong peak.

## Evidence (numbers from the data report)
All three are heavy peptides; all masses ~+10 ppm high vs theory (uncalibrated raw). "k" = isotope offset of our reported mono vs the true mono.

### #15 TLNFNAE[K]GEPELLM[ox]LANWRPAQPLK — z4, RT 16.342, ref mono 2905.4199 (the cleanest example)
- **True mono is strong**: at z4 its m/z is 727.362; in the observed envelope it is **0.62 of the mode (composite) and IS the mode (1.00) at the apex scan**. Not a weak/ambiguous mono.
- Resolved answer: **mono 2904.4358, k = −0.98** (one ¹³C **low**), merging **z3;4;5** (4 members), primary z4.
- The dominant z4 (int ~6e8) placed its mono one isotope low (727.116 — **no peak there**). The z3 that had the mono correct (2905.436, k≈0, int ~2.9e8) was **outvoted**; the z3 member actually merged is its +1. Per-charge detected m/z even include **z5:582.095 = the true mono**, yet the reported mono is stepped down to 727.116.
- Net: the reported mono sits where there is **no observed peak**, while the true (strong) mono is discarded by the merge.

### #7 [Carbamyl-on-M]TLNFNAEGEPELLMLANWRPAQPLK — z3, RT 16.406, ref mono 2894.4749
- **True mono is strong**: at z3 its m/z is 965.832; **0.65 of mode (composite), 0.62 (apex)**.
- Near the reference RT we produced only off-by-one-low fragments; nearest z3 is **mono 2893.5072, k = −0.96**, reported mono m/z **965.510 which has only 0.005 of mode = noise (no real peak)**. Its detected peak 965.844 *is* the true mono, but it was labelled +1 and the mono stepped down one.
- **Extra wrinkle (modform/RT tangle):** a correctly-placed merged feature (mono 2894.5012, k≈0, z2;3;4, huge) DOES exist — but at **RT 16.786**, matched to the *other* co-eluting isobaric modform ([Carbamyl-on-X], same 2894.47 mass, different site). So #7 is an off-by-one **plus** an isobaric-modform/RT-assignment issue; treat #15 as the clean case and use #7 to check the fix survives the modform tangle.

### #5 MVNNGHSFNVEYDDSQDKAVLK (2× ammonia-loss) — z3, RT 12.66, ref mono 2475.1013 — the MIRROR case (keep for contrast)
- Here the true mono **is genuinely weak** (0.21 of mode in composite, **absent at apex**); our reported mono (826.383) sits on a **real** peak (the +1, = the reference's own Peak MZ). Off by **+1 high**.
- So #5 is *not* the same bug — it is the weak-mono/averagine case. A fix for #7/#15 should be checked to not disturb #5 (and ideally still miss-or-fix it honestly). Do **not** lump #5 in.

## Mechanism (hypothesis to confirm)
1. **Per-charge mis-index.** At the offending charge the shift decon anchors on the detector's most-abundant (mode) peak and locates the mono as `anchor − mode_index·spacing`, where `mode_index = argmax(averagine_comb_weights(most_intense_mass))`. If the **averagine mode-index disagrees with the observed mode position by one** (observed mode at +2 but averagine predicts +3 for these ~2900 Da peptides, or vice-versa), the mono lands one isotope off — onto an empty position. #7/#15 step **down** (−1); #5 steps **up** (+1). **Verify this** by instrumenting `shift_decon`/`best_charge_by_fit` for these features (print `mode_index`, `best_shift`, the per-shift cosines, and the observed intensity at the placed-mono position vs the ±1 positions).
2. **Merge locks it in.** `resolve_mass_by_cross_charge` clusters candidate masses and ranks clusters by `(distinct-charge support, count, max fit, summed intensity)`. When two clusters are 1 ¹³C apart and the mis-indexed cluster has more charges/higher fit/intensity, it wins — overriding the correctly-placed charge. This is the same lever as the ECC case, **except the mis-indexed members are the SAME peptide at a mis-indexed charge, not a different co-eluting species** — so the co-elution knitting gate (which separates non-co-eluting species) does NOT help here; these genuinely co-elute.

## Why existing machinery doesn't catch it
- **Walk-back** (`walkback_mono_high_charge`, gated mass ≥ 2200 Da & |z| ≥ 3) only corrects **too-high** placements (walks the mono *down*). #7/#15 are too-**low**, so it can't help and would make them worse if it fired.
- **envelope_fit_cosine window** starts at `mono − 0.5·spacing`; a placement whose mono tooth is empty is only lightly penalised because the averagine mono weight is small — not enough to lose to the correct placement.
- **Anti-doubling / anti-halving guards** are about charge magnitude, not isotope index — orthogonal.

## Candidate fixes (in rough order of promise / safety)
1. **Consensus: prefer a mono placement that has a real peak at its mono.** In `resolve_mass_by_cross_charge` (`feature_refinement.rs`), when the top clusters differ by ~1 ¹³C, break toward the cluster whose monoisotope m/z actually carries observed intensity (checked at each member's charge on that member's slice), rather than an empty position. This is the user-favoured "faulty charge-state merging" angle and is the most targeted. **Risk:** needs the observed slice at resolve time (currently the consensus works on refined scalars, not spectra) — either thread a light "mono has signal?" flag from refine, or recompute against the apex scan. Measure carefully; consensus changes have regressed before (the fit-weighted-*primary* variant regressed; only the fit-*tiebreak* survived).
2. **Walk-*up* analog of the walk-back.** A sibling to `walkback_mono_high_charge` that, for the same heavy/high-charge regime, checks whether the mono was placed **too low**: if the placed-mono position has little/no observed signal but the +1 position has a strong peak (inconsistent with the averagine mono/+1 ratio — i.e. a "leading tooth absent, second tooth dominant" pattern), step the mono **up** one. Must be gated tightly and symmetric-safe (must not fire on #5, where the mono is legitimately weak — the discriminator is *no peak at all* at the placed mono vs *a small peak*). Wire it into `refine_feature_shift` next to the existing walk-back.
3. **Root-cause the mode-index.** If step (1) of the mechanism confirms an averagine mode-index mismatch for ~2900 Da peptides, consider a mass-dependent correction to `mode_index` or a shift-search that also verifies the mono tooth is populated. Upstream and cleaner but requires care with the parity-locked classic path (do NOT edit `deconvolution.rs`'s parity-gated functions; the averagine helpers `averagine_intensities_from_mono` / `averagine_comb_weights` are shared — changing them ripples everywhere and must be re-validated against the classic golden).

## Reproduce / measure
- **Benchmark** (~90 s + build; PowerShell, bash is unreliable here):
  ```
  cd F:\flashlfq-rust\rust; $env:COVERAGE_TARGET="0.90"
  rtk proxy cargo run --release --example detect_features_tsv -- `
    "D:\SP_Tutorial\Lumos\04-17-23_CA_Tryp_HCD_10min.raw" "<OUT.tsv>" `
    "D:\SP_Tutorial\Lumos\CA_HCD_GPTMD_Search_WideTol\Task2-SearchTask\AllQuantifiedPeaks.tsv"
  ```
  Prints `reference peaks rediscovered ... N / 620` and `charge state also matched: M`. Baseline **604 / 585**. Write experiments to a distinct OUT path if another process is reading `untargeted_cov90_shiftapex.*`.
- **Envelope probe** for a single case (composite + apex + per-shift decon): build an `OBO_TARGETS` TSV (cols `label ref_mono rt z pk_mz reported_k`) and run `cargo run --release --example obo_envelope_probe` with `OBO_TARGETS=<file>` (and `EXPORT_DIR=<dir>` for the long-form series). See `analysis/misses/` for the current export.
- **Data**: report at `agent_info/misses-5-7-15-data.md`; outputs at `D:\SP_Tutorial\Lumos\untargeted_cov90_shiftapex.{detected,refined,}.tsv`; reference `...\AllQuantifiedPeaks.tsv` (read the header — precursor-charge vs peak-charge columns are easy to confuse).
- Python via `& "F:\flashlfq-rust\.venv\Scripts\python.exe"`; git via `git -C F:\flashlfq-rust`.

## Constraints / pitfalls
- **Do not regress 604/585.** Every change measured on the full benchmark; add a unit test for any fix.
- **Do not touch the parity-locked classic deconvolution** (`deconvolution.rs` `classic_deconvolute` and friends). The averagine helpers are shared — changing them needs classic-golden re-validation.
- Consensus is high-leverage but has bitten before: fit-weighted-*primary* regressed; the fit-*tiebreak* and the co-elution gate worked. Prefer a narrow, tie-break-level change over re-ranking.
- Watch the **#7 modform/RT tangle** — don't overfit to it; validate on #15 (clean) first, then confirm #7 doesn't get worse.
- Relevant code: `feature_refinement.rs` (`resolve_charge_state_consensus`, `resolve_group`, `resolve_mass_by_cross_charge`, `features_link`, `refine_feature_shift`), `isotope_shift_decon.rs` (`shift_decon`, `best_charge_by_fit`, `walkback_mono_high_charge`, `envelope_fit_cosine`). Related memory: `charge-state-knitting-coelution.md`, `light-z2-charge-halving.md`, `envelope-fit-metric.md`.
