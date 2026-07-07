# Untargeted MS1 Feature Detection — Algorithm Explainer

A developer-oriented walkthrough of how the untargeted feature-detection pipeline works end
to end: **index → detect → refine → resolve**. This distills the authoritative design in
`agent_info/Feature-Detection-Design.md`; where the two disagree, the design doc wins. All
symbol references are `file.rs:symbol` relative to `rust/flashlfq-core/src/`.

## What problem this solves

Base FlashLFQ is *targeted*: it starts from a PSM (`Identification` = sequence + mass + RT +
charge) and quantifies that known species. This pipeline is the *untargeted* superset. It takes
**raw MS1 scans** — centroided `(m/z, intensity)` peak lists per retention time — and produces a
**de novo feature map**: for every isotope-pattern feature it can find, a
`(monoisotopic neutral mass, charge state(s), RT extent, summed intensity, score)` record, with
**no** identification required. PSM identity, if ever attached, is only an annotation on top.

The hard part is inference. A targeted run is *given* the mass, the charge, and the theoretical
isotope distribution. Untargeted detection is given none of those — from a bare observed peak it
must infer charge (test z = 1..N), the monoisotopic mass (walk the isotope ladder), and whether
the pattern is real (correlate against an expected envelope). That missing "ID-free
deconvolution" layer is the ported mzLib `ClassicDeconvolutionAlgorithm` plus a purpose-built
matched-filter detector.

## Pipeline at a glance

```
 raw file (.raw/.mzML)
        │  read_ms1_scans / collect_ms1_scans           peak_indexing.rs
        ▼
 ┌──────────────────────────────────────────────────────────────────────┐
 │ 1. INDEX      PeakIndexingEngine::index_peaks                          │
 │   Scan[] ──► binned m/z index (5 bins/Da) + ScanInfo[]                 │
 │   fast (m/z, scan) lookups for isotope-ladder walking                  │
 └──────────────────────────────────────────────────────────────────────┘
        ▼
 ┌──────────────────────────────────────────────────────────────────────┐
 │ 2. DETECT     trace_kernel.rs::detect_features                         │
 │   seeds tallest-first ─► for z=1..6 lay isotope comb × RT Gaussian     │
 │   (sparse 2-D matched filter) ─► keep best-z (cross-z NMS) ─► claim    │
 │   traced elution extent ─► DetectedFeature[]  (coarse mass + charge)   │
 └──────────────────────────────────────────────────────────────────────┘
        ▼
 ┌──────────────────────────────────────────────────────────────────────┐
 │ 3. REFINE     feature_refinement.rs::refine_feature_shift[_neighbor]   │
 │   per feature: apex slice (or apex±k averaged composite) ─► shift/     │
 │   classic decon ─► precise monoisotope + recharge  ─► RefinedFeature[] │
 └──────────────────────────────────────────────────────────────────────┘
        ▼
 ┌──────────────────────────────────────────────────────────────────────┐
 │ 4. RESOLVE    feature_refinement.rs::resolve_charge_state_consensus    │
 │   union-find group co-eluting same-mass charges ─► cross-charge        │
 │   consensus corrects ±1 Da mono off-by-one ─► ResolvedFeature[]        │
 └──────────────────────────────────────────────────────────────────────┘
        ▼
 peptide-level feature map (TSV / optional .msalign)
```

Driver that wires all four stages together: `examples/detect_features_tsv.rs`.

## Stage 1 — Index

**What / why.** Deconvolution and comb-walking need to ask "is there a peak near m/z X in scan
S?" millions of times. A linear scan of every centroid would be hopeless, so the raw scans are
loaded once into a binned m/z index that answers those lookups in O(1)-ish.

**Inputs / outputs.**
- In: a spectra file. `peak_indexing.rs:read_ms1_scans` (dispatches mzML vs Thermo `.raw`) →
  `Vec<Scan>`. Each `peak_indexing.rs:Scan` is one MS1 scan: parallel ascending `mz` /
  `intensity` arrays plus `retention_time`, `one_based_scan_number`, `msn_order`.
- Out: `peak_indexing.rs:PeakIndexingEngine` via `PeakIndexingEngine::index_peaks(&scans)`.

**How.** `index_peaks` is a faithful port of mzLib `IndexingEngine<T>.IndexPeaks`. It builds
`indexed_peaks: Vec<Option<Vec<IndexedMassSpectralPeak>>>` sized to
`max_mz * BINS_PER_DALTON`, bucketing every peak into bin
`round_ties_even(mz * BINS_PER_DALTON)` (banker's rounding, matching C#). `BINS_PER_DALTON = 5`
→ 0.2 Da bins. It also records a parallel `ScanInfo[]` (RT + scan numbers). Returns `None` when
there is nothing to index (mirrors the C# null/false).

**Key API surface used downstream:**
- `PeakIndexingEngine::get_indexed_peak(mz, scan_index, ppm)` — the single most-called method;
  returns the closest peak to `mz` in that scan within tolerance (the comb-tooth lookup).
- `get_xic` / `get_xic_by_scan_index` — extracted-ion chromatogram (a peak's trace across RT),
  used by the claim-extent walk and the data-driven FWHM probe.
- `get_all_xics` — the greedy, intensity-descending, claim-as-you-go seed generator (the design's
  original seed path; the detector reuses its tallest-first + claiming discipline via `all_peaks`).
- `all_peaks` / `scan_info` — flatten the index / expose RT metadata.

**Knobs.** `BINS_PER_DALTON = 5` (bin width). `ppm_tolerance` on lookups is supplied by the
caller (detector default 10 ppm). Zero-intensity peaks are dropped (`remove_zero_intensity_peaks`,
`ZERO_EQUIVALENT_INTENSITY`).

## Stage 2 — Detect

**What / why.** Turn the raw peak cloud into coarse feature candidates: charge, an approximate
monoisotopic mass, and the set of peaks the feature owns. This is a **sparse 2-D matched filter** —
an isotope comb along m/z crossed with a Gaussian along RT — evaluated only at the comb positions
the index says could exist (never rasterized). The template *is* the shape test; there is no
separate data-vs-data correlation gate at detection time. Real LC peaks tail (EMG), but a Gaussian
template detects a tailed peak fine; tailing only matters for integration bounds, handled elsewhere.

**Entry point:** `trace_kernel.rs:detect_features(engine, params) -> Vec<DetectedFeature>`.

**Algorithm (per the loop in `detect_features`):**
1. **Seed ordering.** `engine.all_peaks()` sorted intensity-descending. Busiest/tallest signal
   first maximizes explained-intensity-per-feature. A `claimed: HashSet<PeakKey>` enforces
   claim-as-you-go: already-claimed seeds are skipped.
2. **RT window.** `trace_kernel.rs:seed_rt_window` computes the charge-independent scan window +
   Gaussian weights once per seed (shared across all charge hypotheses). Bounded **in time**
   (`rt_half_window_minutes`, ≈2σ), not in scan count — DDA interleaves a variable number of MS2
   scans, so a fixed ±N-scan window spans wildly different times and over-claims.
3. **Score every charge z = min_charge..=max_charge.** `trace_kernel.rs:score_hypothesis`:
   - Lays the isotope comb, **anchored so the most-abundant tooth `i*` sits on the seed** →
     `mono_mz = seed_mz − i*·spacing/z`, `spacing = C13_MINUS_C12 / z`
     (`C13_MINUS_C12 = 1.0033548`). Comb runs **both directions**: for mass ≳1500 Da the mono is
     at *lower* m/z than the tallest peak, so a one-directional comb would misassign the mass.
   - Envelope weights `wₖ` from `trace_kernel.rs:comb_weights` — default
     `CombWeightModel::Averagine` (reuses the precomputed averagine table), alternative
     `CombWeightModel::Poisson` (`λ = 0.00048·M`, one parameter, no table).
   - Each comb `(isotope k, scan s)` slot contributes `template = wₖ·gₛ` and the observed peak
     intensity from `get_indexed_peak`. A per-hypothesis `used` set prevents two teeth resolving
     to the same physical peak at high z; the shared `claimed` set makes a peak absent for this
     hypothesis if another feature already owns it. Response via
     `trace_kernel.rs:hypothesis_response` (default `ScoreModel::RawSum` = Σ `template·observed`).
4. **Cross-z non-max suppression.** A z=2 comb is a subset of z=4/z=6, so a real z=2 partly fires
   higher-z kernels. Keep only the highest-response z; its harmonics never win.
5. **Acceptance gates.** Need `num_isotopes_observed ≥ min_isotopes_observed` (default 2 — a lone
   peak is not a feature) and `response > 0`. Rejected seeds are retired into `claimed`.
6. **Claim the true traced extent.** `trace_kernel.rs:trace_claim_extent` walks the most-abundant
   tooth outward in RT via `get_xic_by_scan_index` (missed-scan tolerance
   `trace_missed_scans_allowed`, runaway guard `trace_max_half_width_minutes`) — the *whole*
   elution, wider than the ~2σ scoring window, so smaller adjacent seeds can't re-fire as
   fragments. A chromatographic-persistence gate (`min_feature_scans`, default 2) rejects
   single-scan noise doublets.
7. **Claim + emit.** All traced peaks go into `claimed`; `trace_kernel.rs:build_feature` derives
   apex, RT bounds, and summed intensity from that traced set. The `score` deliberately stays the
   narrow matched-filter response (a shape-fit, not an extent-sum).
8. **Stop** when explained intensity reaches `coverage_target · ΣTIC`, or when seeds fall below
   `min_seed_intensity` (seeds are sorted, so this early-breaks the loop).

**Output:** `trace_kernel.rs:DetectedFeature { monoisotopic_mass (coarse), charge, mono_mz,
apex_scan_index, apex_rt, start_rt, end_rt, summed_intensity, score, num_isotopes_observed,
peaks }`.

**Key parameters / defaults** (`trace_kernel.rs:TraceKernelParameters`, `Default`):
charge **1–6**, **10 ppm**, `weight_model = Averagine`, `min_isotopes_observed = 2`,
`min_feature_scans = 2`, `max_isotopes = 12`, `min_isotope_weight = 1e-3`,
`coverage_target = 1.0` (disabled) / `min_seed_intensity = 0.0` in the struct default — the
example driver overrides these (coverage 0.90, seed floor 1000). RT σ and window are placeholders
in `Default` and **must** be set from the data: `with_rt_from_index` (measures FWHM via XIC
half-max, `ASSUMED_FWHM_SEC` fallback = 36 s) or `with_rt_from_scans`. `FWHM_TO_SIGMA`,
`FWHM_FLOOR_SEC`/`FWHM_CEIL_SEC` bound the estimate.

## Stage 3 — Refine

**What / why.** The detector's mass is coarse (comb anchored on the tallest peak; the true
monoisotope may be one or two ¹³C away). Refinement re-deconvolves each feature against a cleaner
spectrum to place a **precise** monoisotopic mass and confirm/correct the charge.

**Entry points** (selected by `REFINE_METHOD`, default `shift_apex`):
- `feature_refinement.rs:refine_feature_shift` / `refine_feature_shift_neighbor` → the default
  detector-anchored FlashLFQ-style **shift deconvolution** (`isotope_shift_decon.rs`).
- `feature_refinement.rs:refine_feature` / `refine_feature_censored` → the parity-locked
  **classic** deconvolution path (`REFINE_METHOD=classic`).

**How (shift path, `refine_feature_shift_inner`):**
1. **Build the spectrum slice.** `feature_refinement.rs:build_feature_slices` cuts the feature's
   local m/z window. Default is **apex-only** (`shift_apex`): just the apex scan's slice — on
   CA/Lumos 10-min this beats the averaged composite (97.4% vs 96.0% recall) and skips the
   averaging cost. `shift_composite` instead averages **apex ± `MAX_SCANS_TO_AVERAGE/2`** scans
   (`MAX_SCANS_TO_AVERAGE = 3`, i.e. apex ± 1) via `spectral_averaging::average_spectra` (mzLib `SpectralAveraging`
   port: TIC-normalized m/z binning, `RelativeToTics` + `WeightEvenly` + 0.01 bin, weighted-mean
   collapse) to raise SNR ~√N.
2. **Anchor.** The detector's own most-abundant claimed peak — an anchor that *cannot* grab a
   foreign peak.
3. **Charge + mono.** With `recharge` (default on for shift methods, `RECHARGE=0` to disable):
   `isotope_shift_decon.rs:best_charge_by_fit` re-scores `charge_candidates` by the envelope-fit
   cosine (fit × explained × completeness), keeping the detector's charge unless another fits
   clearly better (`RECHARGE_PREFER_MARGIN`) — this recovers charge-halved z=2→z=1 features while
   blocking spurious z↔2z harmonic flips. Without recharge it runs a single `shift_decon` at the
   detector's charge, testing integer ¹³C shifts of the monoisotope.
4. **High-charge walk-back.** `walkback_mono_high_charge` nudges the mono down up to 2 ¹³C for
   heavy/high-z peptides whose envelope mode sits above the true mono.
5. **Optional neighbor masking** (`NEIGHBOR_REFINE`): peaks belonging to co-eluting *other*
   features (off this feature's own isotope grid, and only neighbours ≥ `NEIGHBOR_MIN_RATIO`×
   more intense) are masked from the fit, so a weak feature in a strong neighbour's shadow is
   judged on plausibly-its-own signal. `NEIGHBOR_REFINE=iterative` refines in descending-score
   order, each confident feature (fit ≥ `NEIGHBOR_LOCK_MIN_SCORE`) locking its corrected grid for
   the lower-scoring ones that follow.
6. **Optional `JOINT_FIT`** (`examples/…::apply_joint_fit_pass` using `joint_fit.rs`): for still
   low-scoring features, model the window as the target's averagine envelope **plus** overlapping
   neighbours (NNLS) and search the target mono over ¹³C shifts to maximize the joint fit —
   correcting a mis-placement the single-envelope refine can't because a neighbour's peaks
   confused it. Accepted only if it clears a margin.

**Classic path (`refine_feature`):** averages the apex±k composite, runs the parity-gated
`deconvolution.rs:classic_deconvolute` (faithful mzLib `ClassicDeconvolutionAlgorithm.Deconvolute`
port: charge/mono search against averagine, ratio bounds, adjacent-charge folding), picks the
same-charge envelope closest to the seed. `DeconEnvelope.monoisotopic_mass_predictions` surfaces
the per-peak candidate masses.

**Output:** `feature_refinement.rs:RefinedFeature { detected (cloned), refined_monoisotopic_mass,
refined_charge, candidate_masses, decon_score }`. `candidate_masses` is the raw prediction list the
next stage intersects.

**Knobs / defaults:** `REFINE_METHOD=shift_apex` (default), `RECHARGE` (on), `shift_tol_ppm = 20`,
`MAX_SCANS_TO_AVERAGE = 3`, decon `ClassicDeconvolutionParameters` = charge 1–6, 10 ppm, intensity
ratio 3, positive. `NEIGHBOR_REFINE`, `NEIGHBOR_MIN_RATIO` (5.0), `NEIGHBOR_LOCK_MIN_SCORE` (0.9),
`JOINT_FIT` / `JOINT_FIT_MAX_SCORE` (0.7) / `JOINT_FIT_MIN_GAIN` (0.02), `CENSOR_CLAIMED`.

## Stage 4 — Resolve

**What / why.** One peptide appears as several co-eluting charge states, each independently
refined. Merge them into one peptide-level record — and use the redundancy to fix the
**monoisotopic ±1 Da off-by-one**. The true mono appears in *every* charge state's candidate set;
a spurious ±1 Da prediction does not consistently recur across charges. So intersecting candidates
across charges *resolves* the off-by-one, not merely averages down noise.

**Entry point:**
`feature_refinement.rs:resolve_charge_state_consensus(refined, mass_tolerance_ppm, rt_tolerance_minutes)`.

**How:**
1. **Group** (`group_features`, single-linkage union-find via `features_link`): link two refined
   features when they co-elute (`|apex_rt_i − apex_rt_j| ≤ rt_tol`, RT-range overlap aware) **and**
   their neutral masses agree within `mass_tolerance_ppm` **after allowing an integer ¹³C
   off-by-one** (`MAX_OFFBYONE_UNITS`). Grouping on RT-range *overlap* (not just apex proximity)
   collapsed residual duplicate features to ~0%.
2. **Resolve mass** (`resolve_group` → `resolve_mass_by_cross_charge`): pool `candidate_masses`
   across members, cluster within tolerance (`MASS_CLUSTER_ABS_FLOOR_DA`), and pick the cluster
   with the greatest **cross-charge support** (softened from strict "all charges" to "≥2 charges,
   weighted by support count"). The resolved mass is that cluster's weighted mean.

**Output:** `feature_refinement.rs:ResolvedFeature { monoisotopic_mass (consensus, off-by-one
corrected), charge_states, apex_rt, start_rt, end_rt, summed_intensity, cross_charge_support,
members }` — the final peptide-level feature map.

**Knobs:** example driver calls it with **10 ppm** mass tolerance and **0.1 min** RT tolerance.

## Where to look / how to run

- **End-to-end driver:** `rust/flashlfq-core/examples/detect_features_tsv.rs` — the canonical call
  order (`read_ms1_scans` → `index_peaks` → `detect_features` → `refine_feature_shift` →
  `resolve_charge_state_consensus`), plus TSV writers and PSM-recall comparison.

  ```
  cargo run --release --example detect_features_tsv -- <spectra_file> <out.tsv> [reference.tsv]
  ```

  Writes `out.detected.tsv` (raw detections), `out.refined.tsv` (post-decon), `out.tsv` (resolved),
  and appends per-step timings to `out.log`. A `reference.tsv` (base-FlashLFQ
  `AllQuantifiedPeaks.tsv`) triggers a rediscovery report (mass ±20 ppm, RT ±0.3 min).

- **Design source of truth:** `agent_info/Feature-Detection-Design.md` (goal, validation gates,
  session-by-session rationale for every default).

- **Environment knobs** (read in `examples/detect_features_tsv.rs`, unless noted):
  - Detect: `COMB_MODEL` (averagine|poisson), `COVERAGE_TARGET` (0.90), `ASSUMED_FWHM_SEC` (36),
    `MIN_SEED_INTENSITY` (1000), `MIN_FEATURE_SCANS` (2), `SCORE_MODEL` (raw|normalized),
    `NOISE_PCT`, `AMP_SEED`, `SCORE_COSINE`, `FIXED_SIGMA`, `TRACE_MISSED_SCANS`,
    `TRACE_MAX_HALF_WIDTH_SEC`, `DETECT_ONLY`, `DETECT_PROFILE` (read in `detect_features` itself —
    per-substage timing).
  - Refine: `REFINE_METHOD` (shift_apex default; classic | shift_composite), `RECHARGE` (on),
    `NEIGHBOR_REFINE` (detected|refined|iterative), `NEIGHBOR_MIN_RATIO` (5.0),
    `NEIGHBOR_LOCK_MIN_SCORE` (0.9), `JOINT_FIT` + `JOINT_FIT_MAX_SCORE`/`MIN_GAIN`,
    `CENSOR_CLAIMED`.
  - Diagnostics / output: `FOUR_WAY_DECON` (classic×shift over composite×apex disagreement
    comparator), `DIFF_REFINE`, `DETECTOR_ANCHOR` / `NEIGHBOR_AWARE` (four-way anchor gating),
    `MSALIGN_OUT` (also emit an MS1 `.msalign` readable by mzLib's `Ms1Align`).

- **Core modules:** `peak_indexing.rs` (index + XICs), `trace_kernel.rs` (detector),
  `deconvolution.rs` (parity-gated classic decon + averagine), `isotope_shift_decon.rs` (shift
  decon, recharge, walk-back, envelope-fit cosine), `joint_fit.rs` (NNLS multi-envelope fit),
  `spectral_averaging.rs` (composite building), `feature_refinement.rs` (refine + resolve),
  `feature_export.rs` (`.msalign` writer).
