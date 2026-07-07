# Quantification — Kick-off (scoping + benchmark design)

> **STATUS (2026-07-07): DEFERRED — do not start quant work yet (per Alex).**
> Quantification is on hold until features have been **matched to peptides** first. Once a feature
> carries a peptide/amino-acid composition, abundance should come from an **NNLS fit against the
> precise, composition-specific isotope envelope** rather than an averagine model — composition-exact
> envelopes are more accurate than averagine, so fitting before ID would bake in avoidable error.
> Treat the plan below as background scoping only; the recommended "build intensity now on CA/Lumos"
> first step is **superseded** by this decision. Revisit after feature↔peptide matching exists.

Scopes the quantification work for the untargeted MS1 feature-detection pipeline: what
"quantification" means here, the current state of the code, proposed test cases, a head-to-head
benchmark against FlashLFQ classic, and a phased plan. Companion to
`Feature-Detection-Design.md` (source of truth) and `TODO.md` (test-data reference section).

TODO item this addresses: **[med] Quantification — kick-off.**

## 1. What "quantification" means here, and the current state of the code

For this pipeline, quantification = turning a **detected feature** (neutral mass + charge + RT
extent) into an **abundance** that is comparable across peptides within a run and across runs.
The pipeline is de novo (no PSM required to detect), so quant must also be ID-free: the number
attached to a feature has to be defensible without a sequence.

There are **two distinct quant paths in the tree today**, and only one of them is mature.

### 1a. Targeted / classic FlashLFQ path — mature, this is real quant
- `src/chromatographic_peak.rs::ChromatographicPeak::calculate_intensity_for_this_feature(integrate)`
  (port of `CalculateIntensityForThisFeature`). Sets `self.intensity` to either the **apex
  isotopic-envelope intensity** (`integrate = false`) or the **sum of all isotopic-envelope
  intensities** (`integrate = true`). This is the FlashLFQ definition of a peak's abundance:
  per-scan an `IsotopicEnvelope` already aggregates the monoisotope + isotope peaks into one
  intensity; the peak's abundance is then apex-height or integrated-area over those per-scan
  envelope intensities.
- `src/chromatographic_peak.rs::cut_peak` — valley-based RT boundary trimming, so the integral
  is over a defined elution, not the whole XIC.
- `src/results.rs` — `PeptideQuant` / `PeptideResults` / `calculate_peptide_results` /
  `engine_quantify_set` roll peak intensities up to **peptide-level** abundances (per file).
- `src/mbr.rs`, `mbr_search.rs`, `mbr_scorer.rs`, `mbr_chromatographic_peak.rs` — match-between-runs:
  transfer a donor peptide's peak into an acceptor file and quantify it, with a PEP model.
- `src/parquet_output.rs` — `write_peptide_results_parquet` (peptide abundances) and the MBR
  **feature table** (`feature_table_schema`, one row per transferred MBR candidate: intensity +
  five component scores + `mbr_pep`). Note: this "feature table" is the **MBR training table**,
  not the untargeted feature map.

### 1b. Untargeted path (trace_kernel → feature_refinement) — detection is done, quant is a stub
Every untargeted feature carries a single number, `summed_intensity`:
- `src/trace_kernel.rs::DetectedFeature.summed_intensity` — built in `build_feature` as a plain
  `Σ peak.intensity` over **all claimed centroid peaks across the RT window and all isotope teeth**.
- `src/feature_refinement.rs::RefinedFeature` / `ResolvedFeature.summed_intensity` — the same sum,
  aggregated across co-eluting charge states in `resolve_charge_state_consensus`.
- `examples/detect_features_tsv.rs::write_tsv` emits that `Summed Intensity` column.

**What is missing for untargeted quant** (the actual work of this task):
1. **A principled abundance.** `summed_intensity` is a raw sum of every claimed centroid (all
   isotopes × all scans). It is not FlashLFQ-comparable: no apex-vs-integrated choice, no per-scan
   envelope intensity, no baseline, and it double-counts the RT dimension against the isotope
   dimension. Need to define feature abundance the FlashLFQ way — per-scan envelope intensity
   (monoisotope + isotopes) forming an XIC, then apex height **or** integrated area.
2. **RT integration bounds.** `cut_peak` exists but is only wired into the targeted path. Untargeted
   features integrate over the trace-kernel window (±2σ), which over-claims co-eluting neighbours
   in dense regions. Untargeted quant needs an equivalent valley/boundary step.
3. **Cross-file feature correspondence.** Detection is single-file. Any ratio/CV/linearity metric
   needs features aligned across runs (a consensus feature map keyed on neutral mass + RT, with RT
   alignment). Does not exist for the untargeted path.
4. **Normalization across runs.** FlashLFQ does median-ratio normalization between files; there is
   no untargeted equivalent.
5. **Untargeted MBR.** No match-between-runs for *unidentified* features (rediscover a real feature
   in a run where it dipped below the detection score). Out of scope for kickoff; flag as future.
6. **Output plumbing.** The untargeted feature map is TSV-only from the example; the design doc's
   "extend the Parquet writer with a feature-map table" (Feature-Detection-Design §Proposed pipeline
   step 5) is not done.

Bottom line: **detection is solved, targeted quant is solved and golden-gated, untargeted quant is
a single unvalidated sum.** This task is to give the untargeted feature a real, benchmarked abundance.

## 2. Proposed test cases (start set)

Both references share the **identical 26-column FlashLFQ QuantifiedPeaks schema** (verified), so a
single format adapter covers both. Key quant columns: `Peptide Monoisotopic Mass` (theoretical),
`Peak MZ` (observed), `Peak Charge`, `Peak RT Apex/Start/End`, `Peak intensity`,
`Peak Detection Type` (MSMS vs MBR/Imputed), `Organism`, `Decoy Peptide`, `Full Sequence`.

### Case A — CA/Lumos 10-min (primary first case)
- raw: `D:\SP_Tutorial\Lumos\04-17-23_CA_Tryp_HCD_10min.raw` (2689 MS1 scans, ~1.9 s FWHM, ~3 scans/peak).
- ref: `D:\SP_Tutorial\Lumos\CA_HCD_GPTMD_Search_WideTol\Task2-SearchTask\AllQuantifiedPeaks.tsv`
  (~620 peaks, single file — no File-Name filtering needed).
- **Validates:** the core intensity definition and the matcher, on a small, fast, already-understood
  file (detection recall here is already characterized: ~40–41% strict, ~84% loose, gated by mass
  accuracy). Quant question: *for the features we already match, how well does our abundance track
  FlashLFQ's `Peak intensity`?* Fast iteration loop — use this to build and shake out the harness.
- **Gotcha (from memory):** the raw is **uncalibrated**; the reference reports theoretical masses.
  Expect a systematic **+10 ppm** measured-vs-theoretical offset. The matcher must tolerate it (or
  fit and subtract a global ppm offset before scoring).

### Case B — IonStar 120-min spike-in (accuracy/linearity case)
- raws: 20 files at `D:\PXD003881_IonStar_SpikeIn\*.raw` (~20 MS1 scans/peak).
- ref: `...\MM_Master_noMBR_Norm\Task1-SearchTask\Individual File Results\*_QuantifiedPeaks.tsv`.
  **Note:** each per-file-named table actually contains the full **504,784-row pooled** set — filter
  to the target raw via `File Name`, drop `Decoy Peptide = True`, and split on `Peak Detection Type`.
- **Designed spike-in.** Constant **human** background; **E. coli** peptides spiked at 5 levels across
  conditions **A–E** (filenames `human_ecoli_{A..E}`, 4 replicates each). The `Organism` column is the
  discriminator: **human peptides are the null (expected 1:1 across all conditions), E. coli peptides
  carry the known fold-changes.** This is what makes IonStar ideal for **accuracy + linearity** and
  for a **false-ratio (null) distribution** — you get both the true-positive ratio recovery and the
  false-positive spread from one dataset.
- **Validates:** does the abundance recover known ratios (E. coli measured vs designed fold change,
  linearity/R² across A–E), is the human null tight (CV within replicates, median |log-ratio| ≈ 0),
  and does it hold at 120-min scale (~20 scans/peak — stresses the integration window vs the 10-min
  case). Requires cross-file correspondence (missing piece #3), so this is the second case.

(The medium D_Morgen_Glyco case in TODO is a **recall** ground truth — PSM scan-RT table, no
intensity/apex columns — so it is *not* a quant benchmark. Keep it for the averaging-window/recall
tasks, not here.)

## 3. Benchmark methodology vs FlashLFQ classic

FlashLFQ classic IS the reference quant engine; the `Peak intensity` column is its answer for the
same raw. The benchmark is a controlled head-to-head on identical inputs.

### Matching our features to reference peaks
Reuse and extend `examples/detect_features_tsv.rs::compare_to_reference` (today it matches on
`Peptide Monoisotopic Mass` ±20 ppm + `Peak RT Apex` ±0.3 min + `Peak Charge`, and reports
**recall/charge only — it does not compare intensity**). Extend it to also carry
`Peak intensity` so a match yields an **(ours, FlashLFQ) intensity pair**.
- Match key: neutral mass (ppm, with the +10 ppm offset accommodation on the uncalibrated CA raw)
  **+** charge **+** apex RT (lenient). Prefer matching on our resolved monoisotopic mass vs the
  reference **theoretical** `Peptide Monoisotopic Mass` (both are neutral, offset-comparable),
  keeping `Peak MZ` for a sanity cross-check.
- One-to-one greedy assignment (nearest in mass then RT) to avoid a tall feature grabbing several
  reference peaks.

### Metrics
- **Intensity correlation / log-log regression.** On matched pairs, fit `log10(ours)` vs
  `log10(FlashLFQ)`: report slope (≈1 = proportional), intercept (systematic scale offset), R², and
  Pearson **and** Spearman. This is the headline quant fidelity number.
- **Median relative error** `median(|ours − ref| / ref)` and the log-ratio spread (IQR of
  `log2(ours/ref)`) — robust to the intensity-scale offset a raw sum will have.
- **Detection overlap.** Precision/recall/F1 at the match tolerance (already have recall). Report the
  intensity metrics on the matched set separately from overlap, so a quant regression is not hidden
  by a detection change.
- **CV within replicates** (IonStar A–E, 4 reps each): `CV = sd/mean` of our abundance per
  peptide-condition; compare distribution to FlashLFQ's. Lower/comparable CV = usable quant.
- **Spike-in ratio accuracy + linearity** (IonStar): for E. coli peptides, measured fold change vs
  designed across A–E — report bias (median measured/designed) and linearity (R² of measured vs
  designed on log scale). For human peptides, the **null**: median |log2 ratio| ≈ 0 and its spread.
  Compare all of these head-to-head with FlashLFQ's own numbers on the same peaks.

### MBR vs MS1-only rows
The untargeted detector is **de novo MS1 with no MBR** (missing piece #5). Handle the reference's
`Peak Detection Type` explicitly:
- **Primary head-to-head:** restrict to reference rows with `Peak Detection Type = MSMS` (FlashLFQ
  quantified them from a real MS2 ID in that file) — the fair comparison for an MS1-only engine.
- **MBR / Imputed rows:** report separately. These are peaks FlashLFQ only got by transfer; a de novo
  MS1 detector that finds them *is* rediscovering real signal, so they are a **bonus recall** bucket,
  not a miss — but never fold their intensities into the primary correlation (different provenance).
- Always drop `Decoy Peptide = True`.

## 4. Phased plan

**Phase 0 — Kickoff (this doc).** Test cases + benchmark design agreed. No code.

**Phase 1 — Integration method (the core deliverable).** Replace `summed_intensity` with a real
abundance on the untargeted path:
- Aggregate claimed peaks into a **per-scan envelope intensity** (monoisotope + isotopes at that
  scan), forming an XIC — mirror what an `IsotopicEnvelope.intensity` is on the targeted side.
- Define feature abundance = **apex envelope intensity** and **integrated area** (both, selectable),
  matching `calculate_intensity_for_this_feature(integrate)` semantics so the two paths are
  comparable.
- Apply an RT-boundary step (reuse `cut_peak`'s valley logic, or an untargeted analogue) so the
  integral stops at the elution edges instead of the ±2σ window.
- Unit-test on the existing synthetic z=2 envelope fixtures (known area).
- First concrete step below.

**Phase 2 — Benchmark harness.** Extend `compare_to_reference` to emit matched (ours, FlashLFQ)
intensity pairs to TSV; write the QuantifiedPeaks adapter (shared CA + IonStar, File-Name filter +
decoy drop + detection-type split); do the stats/plots in **Python** (per project convention —
log-log scatter, correlation, relative error). Run on **CA first**.

**Phase 3 — Accuracy tuning (single file).** Iterate the integration method against CA correlation;
then run one IonStar single file (e.g. `B03_..._B_..._9B`) to confirm it holds at 120-min scale /
~20 scans/peak.

**Phase 4 — Cross-file consensus + spike-in accuracy.** Build the multi-file feature correspondence
(align on neutral mass + aligned RT) and per-run normalization; run the full IonStar A–E ratio /
linearity / CV / null benchmark head-to-head with FlashLFQ.

### First concrete step
On the CA 10-min case: (1) add an `apex envelope intensity` + `integrated area` to
`DetectedFeature`/`ResolvedFeature` alongside the existing `summed_intensity` (keep the sum for
coverage), computed from per-scan envelope intensities; (2) extend `compare_to_reference` to also
capture `Peak intensity` and dump matched intensity pairs; (3) Python log-log scatter of our apex
intensity vs FlashLFQ `Peak intensity` on the matched set. That single scatter + its slope/R² is the
go/no-go read on whether the abundance definition is sound before touching IonStar.

## 5. Open questions / decisions for Alex

1. **Apex vs integrated area as the default untargeted abundance?** FlashLFQ exposes both
   (`integrate` flag). Which becomes the untargeted default, and do we match FlashLFQ's default for a
   clean head-to-head?
2. **Per-scan envelope intensity: monoisotope-only, top-N isotopes, or all claimed isotopes?**
   FlashLFQ sums the fitted envelope. Do we sum all claimed teeth, or restrict to the isotopes the
   averagine/Poisson model says carry the mass (to avoid pulling interference into the abundance)?
3. **RT boundary for untargeted:** reuse `cut_peak` verbatim, or a matched-filter-native boundary
   (response falls below a fraction of apex)? `cut_peak` was built for targeted envelopes.
4. **CA calibration:** benchmark on the uncalibrated raw with a fitted global ppm offset, or use a
   calibrated CA raw if one exists? Affects match rate and the intensity-pair count.
5. **IonStar designed ratios:** confirm the exact A–E E. coli spike levels (PXD003881 / Shen 2018)
   so linearity is scored against the real design, not just monotonicity. The `Organism` column gives
   us the human-null / E.coli-signal split regardless.
6. **Cross-file correspondence ownership:** does untargeted quant reuse the MBR/`ResolveIdentifications`
   neutral-mass re-keying (design doc's stated plan), or a new consensus module? Determines whether
   Phase 4 is a reuse or a build.
7. **Scope of MBR for untargeted:** explicitly out for kickoff — confirm it stays a later phase.
