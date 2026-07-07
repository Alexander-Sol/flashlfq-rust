# Longer-gradient test cases — ground-truth notes (TASK 1)

Peak-level ground truth built by `build_ground_truth.py` from the two reference tables.
Normalized schema: `mono_mass  mz  charge  rt  intensity  detection`.

## Ground-truth scans-to-average per gradient (Alex's visual estimates)

| gradient | raw | scans/FWHM (ground truth) |
|---|---|---|
| short ~10-min (CA/Lumos) | `D:\SP_Tutorial\Lumos\04-17-23_CA_Tryp_HCD_10min.raw` | ~3 (apex ±1) |
| medium ~65-min (D_Morgen_Glyco) | `HFX_MB_14751_5_02062022.raw` | ~7–9 |
| long ~120-min (IonStar) | `B03_19_150304_..._9B.raw` | ~20 |

The data-dependent averaging formula should reproduce this 3 / 7–9 / 20 progression.

## MEDIUM — D_Morgen_Glyco (Byonic N-glyco PSM table)

- Ref: `...\ConvertedBionicResults\HFX_MB_14751_5_02062022.raw_20230206_WC_Ngly_Generic.tsv`
- Columns: `File Name, Base Sequence, Full Sequence, Peptide Monoisotopic Mass, Precursor Charge,
  Protein Accession, Scan Retention Time` — a **per-PSM** table, NO apex/intensity.
- Adapter: dedupe on unique **(Full Sequence + Precursor Charge)**, keep first-seen RT; compute
  `mz = (mono_mass + z·1.007276) / z`.
- **13,179 PSM rows → 10,021 unique peaks** (`medium_ground_truth.tsv`). Many distinct glycoforms
  (Full Sequence encodes the glycan mass) keep the unique count high.
- **RT gotcha:** `rt` here is the **PSM's Scan Retention Time**, i.e. the scan where the peptide was
  fragmented — NOT the chromatographic apex. It can sit anywhere on the elution flank, so recall
  scoring against this table must use a **lenient RT delta** (see recall report; we use ±0.5 min).
- `intensity` blank; `detection = PSM`.

## LONG — IonStar spike-in (FlashLFQ QuantifiedPeaks)

- Ref: `...\Individual File Results\B03_19_150304_..._9B_QuantifiedPeaks.tsv` (504,784 pooled rows).
- Adapter: filter to `File Name == B03_19_150304_human_ecoli_B_3ul_3um_column_95_HCD_OT_2hrs_30B_9B`,
  drop `Decoy Peptide == True`, split by `Peak Detection Type`.
- Findings on this file's slice: **25,086 rows, ALL `MSMS`, ALL non-decoy (0 decoy, 0 MBR).**
  The search is the **noMBR** variant (`MM_Master_noMBR_Norm`), so MBR rows do not exist — MSMS-only
  == full set here. `long_ground_truth_msms.tsv` and `long_ground_truth_all.tsv` are therefore
  identical (both 24,642 peaks), kept separate to match the recall report's primary/secondary split.
- **444 rows dropped** because `Peak Charge == "-"` — FlashLFQ's marker for a peptide that was
  identified but whose peak could **not** be quantified (no found apex). Correctly excluded: the
  detector cannot be expected to rediscover a peak FlashLFQ itself failed to quantify.
- **→ 24,642 quantified MSMS peaks** (`long_ground_truth_msms.tsv`). `rt = Peak RT Apex` (a real
  apex, RT range 19.7–159.5 min), `mz = Peak MZ`, `charge = Peak Charge`, `intensity = Peak intensity`.
