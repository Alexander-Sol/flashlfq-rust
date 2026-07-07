# Detector run timings — medium + long gradients (TASK 2)

Full untargeted pipeline (read+index → detect → refine → charge-consensus → write) via
`detect_features_tsv` (release), default params (shift_apex refine, coverage 90%, 10 ppm,
data-driven σ from measured FWHM). Reproduce with `run_detector.ps1`. Per-stage timings are the
example's own timing log (`*_features.log`, verbatim below).

## MEDIUM ~65-min — HFX_MB_14751_5_02062022.raw (D_Morgen_Glyco glyco)

7,407 MS1 scans · 16,825,488 peaks · ΣTIC 5.312e12 · median scan spacing 0.0155 min · measured FWHM 9.36 s

| stage | seconds | % |
|---|--:|--:|
| read + index | 18.60 | 2.6 |
| **detect** | **679.29** | **94.1** |
| refine | 14.42 | 2.0 |
| charge-state consensus | 1.55 | 0.2 |
| write resolved | 3.08 | 0.4 |
| **TOTAL** | **722.03** | |

Features: 737,957 detected → 737,957 refined → **665,486 resolved**.

## LONG ~120-min — B03_19_150304_..._9B.raw (IonStar spike-in)

10,582 MS1 scans · 22,528,910 peaks · ΣTIC 2.300e13 · median scan spacing 0.0124 min · measured FWHM 17.66 s

| stage | seconds | % |
|---|--:|--:|
| read + index | 28.50 | 33.4 |
| **detect** | **47.95** | **56.2** |
| refine | 4.27 | 5.0 |
| charge-state consensus | 0.27 | 0.3 |
| write resolved | 0.55 | 0.6 |
| **TOTAL** | **85.32** | |

Features: 114,732 detected → 114,732 refined → **96,058 resolved**.

## Notable — detect cost tracks SEED count, not peak/scan count

The long file is *larger* (22.5M peaks, 10.6k scans vs 16.8M / 7.4k) yet detect ran **14× faster**
(48 s vs 679 s) and produced **6.4× fewer** features (115k vs 738k). Detection cost is dominated by
the number of seeds needed to reach the 90% ΣTIC coverage target, not by raw peak count: the glyco
medium has a highly fragmented TIC (many low-level glycoforms) so 90% coverage demands ~738k
features; the cleaner IonStar proteomics TIC hits 90% with ~115k. This is the relevant signal for
the detector-parallelization task — the pathological case is a fragmented spectrum (glyco), not
merely a long gradient. Detect stays the dominant stage in both (94% / 56%).

## Output files (NOT committed — regenerable, large)

`*_features.tsv` (resolved), `*.detected.tsv`, `*.refined.tsv` were left unstaged. The medium set is
~150 MB total; all are reproducible via `run_detector.ps1`. Only the tiny `*_features.log` timing
logs and this summary are committed.
