# Dinosaur vs. our detector — full PXD001091 synthetic-peptide benchmark (2026-07-09)

Head-to-head on **all 57 raw files** of the PXD001091 synthetic-peptide dilution
series (`130124_dilA_*`), scored against the FragPipe/MSFragger PSMs
(`Fragpip_v23p1_AllFiles\psm.tsv`, 339,398 PSMs, 0 decoys).

**Both detectors run on the same input** — the FragPipe `*_uncalibrated.mzML`
files, which already carry full MS1+MS2 (no msconvert step needed). **Concurrency
equalized at 16** (Dinosaur `--concurrency=16` vs ours `RAYON_NUM_THREADS=16`) on a
12-core / 24-thread Xeon Silver 4214R. Both feature tables scored by the same
scorer (`score_all.py`): mass ±20 ppm + RT window, against the same PSM-derived
ground truth.

## Ground truth

Per-file, deduped from the PSMs to unique **(file, modified peptide, charge)**
peaks — **177,759 expected peaks** across 57 files (mean 3,119/file, range
1,970–3,484). `mono_mass` = MSFragger "Observed Mass" (neutral, isotope-corrected,
uncalibrated). `rt` = median PSM "Retention" ÷ 60 (minutes). Because PSM scan time
is not the chromatographic apex, recall is reported at two RT tolerances.

## Results (aggregate, micro-averaged over 177,759 peaks)

| Metric                | Dinosaur   | Ours        |
|-----------------------|-----------:|------------:|
| recall (RT ±0.5 min)  |   97.30 %  |   **97.76 %** |
| recall (RT ±1.0 min)  |   97.82 %  |   **98.13 %** |
| total features        | 2,659,432  | 19,896,214  |
| runtime — total       |  2,556 s (42.6 min) | **1,768 s (29.5 min)** |
| runtime — mean/file   |   44.8 s   |   **31.0 s** |
| runtime — range/file  |  39.4–62.6 s | 25.3–55.5 s |

Per-file recall @ ±0.5 min: **ours wins 56/57 files, 1 tie, 0 losses**
(mean advantage +0.46 pp, up to +1.50 pp).

## Reading the numbers

- **Recall:** ours is higher on aggregate at both tolerances (+0.46 pp @ ±0.5,
  +0.31 pp @ ±1.0) and wins essentially every file. On this clean synthetic data
  both detectors are near the ceiling (~97–98 %) — nearly every PSM'd peptide is
  findable — so the gap is small but consistent and never negative.
- **Runtime:** ours is ~1.4× faster at matched concurrency (31.0 vs 44.8 s/file),
  reading the identical mzML. (Ours reads `.raw` directly too — see preliminary —
  but mzML is ~22 % faster with identical recall, so the full run used mzML.)
- **Feature count is NOT a like-for-like quality metric.** Ours emits ~7.5× more
  features (19.9M charge-merged vs Dinosaur's 2.66M per-charge). Ours is far more
  permissive (`MIN_SEED_INTENSITY=1000` admits a large low-intensity tail); it buys
  the recall edge at the cost of much lower specificity. Dinosaur is more
  parsimonious.

## Charge — deliberately omitted

Charge accuracy is **not reported as a headline metric** on this benchmark. The GT
`mono_mass` is the neutral monoisotopic mass (identical across charge states) and
co-eluting charge states share RT, so the mass+RT matcher cannot localize to a
charge — and a PSM's charge is just whichever precursor the instrument happened to
fragment, not ground truth about which charge states exist. The raw numbers
(ours ~91 % vs Dinosaur ~67 % "charge-among-feature-charges") mostly reflect output
format (our merged multi-charge rows vs Dinosaur's one-charge-per-row) rather than
charge-assignment quality. A fair charge comparison would require per-charge
expansion of our rows plus a charge-aware matcher; left out here as low-value.

## Preliminary: raw vs mzML for our detector (2 files, 24 threads)

Recall **identical** within 1 peak / 0.1 pp; mzML **~22 % faster** (33–36 s vs
42–45 s) — the `.raw` path pays .NET Thermo-reader overhead. Full run therefore
used mzML, which also makes ours-vs-Dinosaur a true same-input comparison.

| file            | raw recall ±0.5 | mzml recall ±0.5 | raw s | mzml s |
|-----------------|----------------:|-----------------:|------:|-------:|
| dilA_10_01      |  97.5 %         |  97.5 %          | 42.4  | 33.3   |
| dilA_1_01       |  97.6 %         |  97.6 %          | 45.4  | 36.4   |

## Provenance / how to reproduce

- GT:      `build_gt.py <psm.tsv> gt\`  → `gt\<base>.tsv` + `gt_manifest.tsv`
- Dinosaur:`dino_full.ps1`  → `dino\<base>.features.tsv`, `dino_timings.tsv`
- Ours:    `ours_full.ps1`  → `ours\<base>\feat.tsv`, `ours_timings.tsv`
  (release `detect_features_tsv.exe`, defaults: 2-D tiling, coverage 1.0)
- Score:   `score_all.py <bench_dir> 0.5 1.0` → `results_perfile.tsv` + aggregate
