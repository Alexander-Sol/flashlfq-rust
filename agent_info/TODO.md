# Untargeted feature-detection — TODO

Running task list for the untargeted MS1 feature-detection work. Check items off as they land;
add a one-line result/commit note when closing one. Priority tags: **[high] / [med] / [low]**.

See also: `Feature-Detection-Design.md` (source of truth), `Detector-Improvement-Plan.md`,
`Miss-Characterization-and-Discriminator-Plan.md`.

## Test-data reference (medium + long gradients)

Existing case — **short ~10-min** (CA/Lumos): raw `D:\SP_Tutorial\Lumos\04-17-23_CA_Tryp_HCD_10min.raw`,
ref `...\CA_HCD_GPTMD_Search_WideTol\Task2-SearchTask\AllQuantifiedPeaks.tsv` (~620 peaks). ~1.9 s FWHM.

**Medium ~65-min gradient** — `D:\D_Morgen_Glyco`.
- Refs: `Byonic and MSFragger ID\ConvertedBionicResults\*_Ngly_Generic.tsv`. Each ref's filename prefix
  maps to a `.raw` in `RawFiles\`. Many pairs available; start with one.
- Example pair:
  - raw: `D:\D_Morgen_Glyco\RawFiles\HFX_MB_14751_5_02062022.raw`
  - ref: `D:\D_Morgen_Glyco\Byonic and MSFragger ID\ConvertedBionicResults\HFX_MB_14751_5_02062022.raw_20230206_WC_Ngly_Generic.tsv`
- **~7–9 MS1 scans inside each peak's FWHM** (Alex's visual estimate = the ground-truth ballpark for
  scans-to-average here).
- **Gotcha:** ref reports the **scan RT at which a PSM was collected, NOT the peak apex** → use a
  **lenient RT acceptance delta** when scoring recall.
- Ref columns (differ from AllQuantifiedPeaks — needs a format adapter): `File Name`, `Base Sequence`,
  `Full Sequence`, `Peptide Monoisotopic Mass`, `Precursor Charge`, `Protein Accession`,
  `Scan Retention Time`. No apex/intensity columns. 13,179 rows = per-PSM (dedupe to unique
  peptide+charge for a peak-level ground truth).

**Long ~120-min gradient** — `D:\PXD003881_IonStar_SpikeIn` (IonStar spike-in).
- Raws at the top level of the folder. Refs: `MM_Master_noMBR_Norm\Task1-SearchTask\Individual File
  Results\*_QuantifiedPeaks.tsv`.
- Example pair:
  - raw: `D:\PXD003881_IonStar_SpikeIn\B03_19_150304_human_ecoli_B_3ul_3um_column_95_HCD_OT_2hrs_30B_9B.raw`
  - ref: `...\Individual File Results\B03_19_150304_...9B_QuantifiedPeaks.tsv`
- Ref is a standard FlashLFQ QuantifiedPeaks table (has `Peak RT Apex`, `Peak Charge`, `Peak intensity`,
  `Peak MZ`, `Peak Detection Type`, `Decoy Peptide`, `File Name`). **504,784 rows = many files pooled** →
  filter to this raw's `File Name`, drop `Decoy Peptide`, and consider excluding MBR rows via
  `Peak Detection Type`. The runner's existing compare logic mostly works after that filtering.
- **~20 MS1 scans inside each peak's FWHM** (Alex's estimate = the ground-truth ballpark for
  scans-to-average here).

Scans-to-average ground truth by gradient: short ~10-min ≈ 3 (apex±1) · medium ~65-min ≈ 7–9 ·
long ~120-min ≈ 20. The data-dependent formula should reproduce this progression.

## High priority

### Data-dependent averaging window (umbrella goal)
Replace hard-coded `MAX_SCANS_TO_AVERAGE = 3` with a scans-to-average count derived from the measured
chromatographic FWHM / scans-per-peak, so it adapts to gradient length instead of being tuned to
CA/Lumos 10-min. Broken into the tasks below.

- [ ] **[high] Prep the medium + long test cases.**
  Pick one medium pair (D_Morgen_Glyco) and one long pair (IonStar) to start. Build a format adapter
  for each ref (medium `Ngly_Generic` PSM table; long `QuantifiedPeaks` filtered to one file). Dedupe
  the medium ref to a peak-level ground truth. Record Alex's ground-truth scans-to-average per case
  (medium ≈ 7–9).
- [ ] **[high] Run the detector on the longer files; record run times.**
  Full-pipeline runs on the medium (~65-min) and long (~120-min) raws. Capture per-stage timings
  (read/detect/refine/consensus) — these files are much larger than the 10-min case and will stress
  the detector (see perf/parallelization task).
- [ ] **[high] Test recall on medium + long.**
  Score recall vs each ref. Medium: **lenient RT delta** (ref RT = PSM scan time, not apex). Long:
  filter ref to the single file + non-decoy. Report recall + charge-match like the 10-min case.
- [ ] **[high] Test FWHM + n-scans-to-average calculations.**
  Verify the measured FWHM and the derived scans-to-average on each gradient. Sanity-check against
  ground truth (short ~3, medium ~7–9, long ~20 scans/FWHM). This is where the data-dependent formula
  is validated.
- [ ] **[high] Ensure FWHM info flows correctly to ALL downstream steps.**
  Confirm the measured FWHM / scans-to-average is fed consistently into the averaging window, the
  detector RT σ / matched-filter window, and anything else that assumes a peak width — no step left
  reading a stale hard-coded value.

### Detector performance
- [ ] **[high] Detector performance + parallelization.**
  Profile the detector (~76–79% of total wall-clock; the real cost, not refinement). Find
  optimizations and design a parallelization strategy — evaluate per-file-on-its-own-thread (likely
  simplest/effective) vs. intra-file parallelism. The longer files above make this urgent.

## Medium priority

- [ ] **[med] Turn averaging off completely via params.**
  Today `build_feature_slices` unconditionally calls `average_spectra` (`feature_refinement.rs:738`) —
  even `shift_apex` builds a composite it never reads. Add a real parameter to disable averaging so the
  composite is not computed at all when unused.
- [ ] **[med] Change the default to NOT average (apex-only).**
  A/B on CA/Lumos 10-min: apex-only beats the averaged composite — 97.4% vs 96.0% recall, 94.4% vs
  92.3% charge. Make apex-only the default in the library/param defaults. Coordinate with the
  turn-averaging-off item. NOTE: revisit once data-dependent averaging lands — averaging may help on
  the longer gradients (more scans/peak) even though it hurts on 10-min.
- [ ] **[med] Supported output formats.**
  Determine feature-output formats to emit (e.g. ms-align files, feature files). Must be supported by
  mzLib; ideally widely used. Pick target format(s) and wire an exporter.
- [ ] **[med] Quantification — kick-off.**
  Scope quantification: choose test cases and benchmark against FlashLFQ classic. (The IonStar long
  case ships full FlashLFQ QuantifiedPeaks — a natural quant benchmark.)

## Low priority

- [ ] **[low] Clean up unused / dead code.**
  Sweep `#[allow(dead_code)]` and abandoned experiment paths.
- [ ] **[low] Algorithm explainer.**
  Write up the basic workflow / how the algorithm works (index → detect → refine → resolve).
