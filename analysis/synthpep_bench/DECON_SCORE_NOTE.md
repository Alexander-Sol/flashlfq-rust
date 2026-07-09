# Note: rejecting the low-value feature tail on refine Decon Score (2026-07-09)

**Status: finding + recommendation, not yet implemented. Circle back to wire a
reject threshold into the refiner.**

## The problem

On the full 57-file PXD001091 benchmark our detector emits **~7.5x more features
than Dinosaur** (19.9M resolved / 21.2M refined-per-charge vs Dinosaur's 2.66M) for
only a small recall edge (97.76% vs 97.30% @ RT ±0.5 min). The extra features are a
permissive low-quality tail (documented `MIN_SEED_INTENSITY=1000`). Question we
asked: can we **reject in the refiner** without losing recall?

## The finding

The refine stage already computes a per-feature **`Decon Score`** (envelope-fit
cosine, column in `feat.refined.tsv`, range 0–1). It is a clean, efficient rejection
lever — better than intensity. Sweeping a reject threshold over all 21.2M refined
features vs the 177,759 PSM ground-truth peaks (±20 ppm, ±0.5 min):

| Reject rule        | Refined features dropped | Recall  |
|--------------------|-------------------------:|--------:|
| none (current)     |                       0% | 98.06%  |
| `decon < 0.198`    |                    20.9% | 97.97%  |
| **`decon < 0.664`**|                **82.9%** | **97.57%** |
| `decon < 0.796`    |                    89.3% | 97.07%  |

**Headline: `decon < 0.66` drops ~83% of features for a ~0.5 pp recall cost**
(98.06 → 97.57%, still above Dinosaur's 97.30%). That would bring us from ~19.9M
resolved features down to ≈2.7M — i.e. roughly Dinosaur's parsimony at higher recall.

At equal recall cost, Decon Score prunes 83% vs Summed Intensity's 72% — decon is the
better axis. (Recall stays essentially flat until the threshold is pushed near 1.0;
see the cliff in `score_recall.png` top-left.)

## Why it works (visual confirmation)

`low_quality_features.ipynb` shows individual features. The minimum-decon unmatched
feature (`feat_idx=37501`, z=4, decon=0.000): a pure-noise XIC and MS1 peaks that do
**not** fall on the proposed charge grid — a mis-charged/noise call the cosine
correctly scores at 0. A good feature (z=2, decon=0.989): clean Gaussian XIC +
textbook isotope envelope on the z=2 grid. The score is doing exactly what we'd
threshold on.

## Caveat

"Unmatched" ≠ "junk". Even at high decon there is a large unmatched population, some
of which are **real but unsequenced** species (contaminants/background with no PSM).
So the honest claim is "≥83% of features correspond to no PSM'd peptide and can be
cut for 0.5 pp recall," not "all noise." On a complex-matrix sample the balance would
shift — re-check the sweep there before committing a fixed threshold.

## Recommended next step (TODO)

1. Add an **env-gated reject threshold in the refiner** (e.g. `REFINE_MIN_DECON`,
   default OFF / 0.0 to preserve current behavior), dropping refined features whose
   `decon_score` < threshold before charge-state consensus.
2. A/B feature-count vs recall on this benchmark; confirm the 0.5 pp curve holds
   end-to-end (resolved level, not just refined).
3. Consider a data-driven default (the matched/unmatched crossover, not a hardcoded
   0.66) and re-validate on a complex-matrix file, not just synthetic peptides.

## Reproduce

- Sweep + plot: `python analyze_score_recall.py` → `score_recall.png`,
  `score_recall_curves.tsv` (reads `ours/<base>/feat.refined.tsv` + `gt/<base>.tsv`).
- Per-feature inspection: `python prep_notebook_data.py [base ...]` → `nb_data/`,
  then open `low_quality_features.ipynb`.
- Full benchmark context: `RESULTS.md`.
