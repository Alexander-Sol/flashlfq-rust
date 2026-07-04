# Miss Characterization & Discriminator Plan

Data-grounded follow-up to `agent_info/Detector-Improvement-Plan.md` (Changes 1/A/centroiding/gate/B
done; off-by-one deferred). Characterizes what the detector misses on the CA/Lumos reference and plans
the fixes. Baseline = the default **RawSum** pipeline (centroid raw, data-driven σ ≈ 1.9 s, cov99,
`min_feature_scans=2`): **88.4 % reference recall (548/620), 84.2 % charge-matched.**

Scripts (versioned, `rust/flashlfq-core/scripts/`): `char_misses.py` (full breakdown),
`absent_probe.py` (absent → wrong-charge vs gap), `offbyone_diag.py` (off-by-one recoverability). Data
paths hard-coded near the top (the `D:\SP_Tutorial\Lumos\...` locations).

## The 72 misses are THREE distinct problems

| # | class | count | signature | fix track |
| --- | --- | --- | --- | --- |
| 1 | **off-by-one** | 47 (7.6 %) | a feature exists at ref ± k·¹³C (k=−1:24, +1:20, ±2:3). High-mass/high-charge: z3 misses 20/22 obo, z4 11/13, mass>2500 → 14/15. 26 single-charge, 21 multi (support=1). | **A** |
| 2 | **wrong-charge NMS** | 10 | observed most-abundant m/z IS detected but assigned another charge: every z1 ref stolen by **z2**; low-mass z2 refs (820/907/828 Da) stolen by **z4** (harmonic). | **B** |
| 3 | **detection gap** | 15 | observed m/z claimed by NO feature at that RT (LLVVYPWTQR, YAAELHLVHWNTK family, TEWLDGK, …). | **C** |

Ruled out: **modifications and peptide length are non-factors** (matched avg 1.00 mods / 20.5 aa vs
missed 1.04 / 21.5; modified refs miss 11 %, unmodified 13 %). Misses skew low-intensity (bottom decile
23 % vs top 3 %) **but not exclusively** — the most intense absent peak is 884 Da z1 at 7.35e8 (top
decile), a wrong-charge case.

Ceiling if all three fixed: ~100 %; off-by-one alone → 96 %.

## Track A — off-by-one discriminator (the 47)

Two independent per-feature correctors have already **regressed** (cosine envelope match; detector-comb
anchor snap 88.4 %→80.8 %), because the discriminating signal is weak on **chimeric** averaged
composites. Split the population and attack the cheaper half first:

- **A1 — consensus tiebreak for the 21 multi-charge (support=1).** These features have ≥2 grouped
  charges that *disagree* on the mono (winning cluster supported by one charge). Instead of picking the
  max-support cluster blindly, break support ties toward the candidate mass whose **averagine envelope
  best matches the pooled cross-charge evidence** (each charge's composite re-anchored at the candidate).
  Cross-charge agreement is a stronger, less chimera-sensitive signal than a single feature's envelope.
  Lives in `resolve_charge_state_consensus`; testable without touching detection.
- **A2 — per-feature discriminator for the 26 single-charge.** The hard case (no cross-charge signal).
  Design constraints learned from the two failures: (i) score **only the low isotopes** (mono…mono+3),
  where off-by-one manifests and chimeric contamination at higher teeth is excluded; (ii) use a
  **least-squares residual** of the observed-vs-averagine low-isotope pattern, not a full-envelope
  cosine; (iii) a **strict margin** — only shift when a ±1/±2 anchor fits *dramatically* better, so a
  correct mono (already a good fit) is never flipped; (iv) require a real peak to exist at the shifted
  mono position. **Validate on a held-out CALIBRATED file** (`...CLEAN-calib-calib.mzML`, already
  centroid) as well as the raw, since the discriminator must be non-regressive on both.
- **Acceptance (both):** net recall up with **charge-match not down** and the 548 currently-matched not
  regressing (per-ref diff, not just aggregate). Ship A1 and A2 independently — A1 may pass alone.

## Track B — wrong-charge NMS (the 10)

Cross-z NMS lets a higher-charge comb (a superset of the lower-charge comb's teeth) steal the seed: z1
peptides → z2, low-mass z2 → z4. Change B's normalization was aimed here but **did not help the
aggregate** (charge-match −0.8 %). Re-examine specifically on these 10:

- Does `NormalizedNoiseFloor` fix *these* while breaking others? (per-ref diff on the 10, not aggregate).
- A **charge prior / Occam penalty**: prefer the lowest charge that explains the envelope (penalize
  higher z unless it adds real, non-harmonic teeth) — directly targets harmonic stealing.
- Guard: whatever is tried must hold the 84.2 % charge baseline overall.

## Track C — detection gaps (the 15)

The observed m/z is claimed by no feature at that RT. Investigate the mechanism per case:
- claimed by a co-eluting neighbour's trace extent (Change A over-claim)? — check the claim mask.
- monoisotope too weak to seed / fails the ≥2-isotope gate? — check `.detected.tsv` at that m/z.
- RT apex just outside ±0.3 min (a matching-tolerance artifact, not a real miss)?
Likely a mix; the breakdown decides whether this is a detector gap or a scoring/tolerance issue.

## Sequencing

1. **A1** (consensus tiebreak) — cheapest, likely non-regressive, addresses 21.
2. **Track C triage** — cheap diagnostics; may reclassify some "gaps" as tolerance artifacts.
3. **Track B** (charge prior) — targets 10, follow-up to Change B.
4. **A2** (single-charge discriminator) — hardest, most regression-prone; do last with calibrated-file
   validation.

After each: re-run the default pipeline end-to-end, report the per-ref diff (fixed vs newly-broken),
recall, and charge-match. Ship only on a net gain with no charge regression.
