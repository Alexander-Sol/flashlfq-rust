# Longer-gradient recall — medium (65-min) + long (120-min)

Recall + charge-match of the untargeted detector's **resolved** features vs each adapted peak-level
ground truth (`analysis/longer_gradients/`). Scored by `score_recall.py`, mirroring the example's
built-in `compare_to_reference` matching logic:

- **mass match**: |feature mono − ref mono| / ref mono ≤ **20 ppm**
- **RT match**: |feature apex RT − ref RT| ≤ RT delta (below)
- **charge match**: ref charge ∈ the feature's detected `Charge States`

Reference scale (existing **short 10-min** CA/Lumos case): ~**97%** recall / ~**94%** charge-match.

## MEDIUM ~65-min — D_Morgen_Glyco glyco (10,021 unique peaks)

Ground truth RT is the **PSM scan RT, not the chromatographic apex**, so it can sit anywhere on the
elution flank. A lenient RT delta is required. The measured FWHM here is ~9.4 s (~0.16 min) and a PSM
can be triggered up to roughly a half-peak off apex, so **±0.5 min** is the justified primary delta
(≈3× FWHM — enough to bridge a PSM collected on the shoulder without inviting cross-matches). ±0.3 and
±1.0 shown for sensitivity.

| RT delta | recall | charge-match (of all ref) | charge-match (of matched) |
|---|--:|--:|--:|
| ±0.30 min | 88.5% (8,866/10,021) | 78.2% | 88.4% |
| **±0.50 min (primary)** | **92.2% (9,237/10,021)** | **78.0%** | **84.6%** |
| ±1.00 min | 95.6% (9,581/10,021) | 74.5% | 77.9% |

- Recall climbs with the RT window (88.5 → 92.2 → 95.6%), as expected when the reference RT is a
  fragmentation time rather than an apex — the wider window recovers PSMs collected far off-apex.
- **Charge-match of matched (84.6%) runs below the 10-min ~94%.** Glyco precursors carry higher and
  more varied charge; the primary-charge assignment is the weaker link here, not detection. As the RT
  window widens, charge-of-matched *falls* (88.4 → 84.6 → 77.9%) because the extra matches are looser
  RT coincidences on the wrong species — evidence the ±0.5 primary is a sound cutoff.

## LONG ~120-min — IonStar B03_19 (24,642 MSMS quantified peaks)

Ground truth RT is a **real `Peak RT Apex`**, so a tight delta is appropriate; **±0.3 min** primary
(same as the 10-min case), ±0.5 for sensitivity. This is the **noMBR** search, so there are **no MBR
rows and no decoys** in this file — the MSMS-only primary set **is** the full set (both files
identical, 24,642 peaks). The "vs full set" number therefore equals the primary number.

| RT delta | recall | charge-match (of all ref) | charge-match (of matched) |
|---|--:|--:|--:|
| **±0.30 min (primary)** | **81.2% (20,009/24,642)** | **76.2%** | **93.9%** |
| ±0.50 min | 82.7% (20,371/24,642) | 76.9% | 93.0% |

- **Charge-match of matched (93.9%) matches the 10-min ~94%** — where the detector finds a peak, it
  assigns the charge as reliably on the long gradient as on the short one.
- **Recall (81–83%) is materially below the 10-min ~97%.** ~4,600 FlashLFQ-quantified MSMS peaks are
  not rediscovered. Widening RT to ±0.5 min adds only +1.5%, so this is **not** an RT-tolerance
  artifact — it is genuine missed detection/placement on the long gradient. See flag below.

## Flags for Alex

1. **Long-gradient recall gap (~81% vs ~97% on 10-min) is the headline issue.** It is real (not RT
   slack). Leading hypotheses, in order: (a) the refinement **averaging window is still hard-coded to
   apex±1 (3 scans)** while the long peaks span ~24 scans/FWHM — the composite/placement sees a tiny
   fraction of each wide peak (this is exactly the umbrella-goal gap); (b) the example pipeline's
   co-elution / consensus RT tolerances (0.05 / 0.1 min) are 10-min-scale and under-link charge states
   over wider long-gradient elution; (c) seed floor / 90% coverage may exhaust before reaching the
   ~4,600 lower-abundance IonStar peptides. Worth an A/B once data-dependent averaging lands.
2. **Medium charge-match (~85%) is lower than both the long and short cases** — a glyco-specific
   charge-assignment weakness (high/varied precursor charge), independent of the recall question.
3. Medium recall depends heavily on RT delta because the reference RT is a PSM time, not an apex;
   ±0.5 min is the defensible primary. The long reference has true apexes and is delta-insensitive.
