# Untargeted-detector validation harness

Python analysis scripts that score the untargeted MS1 feature detector against the CA/Lumos
reference peaks (`AllQuantifiedPeaks.tsv`). They consume the TSVs emitted by
`examples/detect_features_tsv.rs` and are the acceptance-gate tooling referenced by
`agent_info/Detector-Improvement-Plan.md`.

Data paths are hard-coded near the top of each script (the `D:\SP_Tutorial\Lumos\...` locations —
see the `lumos-test-data-paths` memory) and the `+10 ppm` calibration caveat noted there. Edit those
constants if the data moves.

| Script | Level | What it reports |
| --- | --- | --- |
| `detect_recall.py` | detection (`.detected.tsv`, pre-refine/consensus) | short/long/all recall vs reference on neutral monoisotopic mass, across the FWHM sweep. The fair comparison for narrow-FWHM settings whose refine/consensus is intractable. |
| `sweep_compare.py` | resolved (`.tsv`) | per-FWHM resolved count, RT-width percentiles, and m/z-merge + neutral-mass recall (overall + short). |
| `merge_99_rt.py` | resolved (single run) | Ref `Peak MZ` vs feature `Most-Abundant m/z` merge carrying RT start/apex/end from both sides; charge-match count; nearest-feature diagnostics for misses. Writes a merged TSV. |
| `char_misses.py` | resolved (single run) | Full miss characterization: classifies each ref as matched / off-by-one / absent; breaks misses down by intensity decile, charge, mass, RT, and sequence mods. Grounds `agent_info/Miss-Characterization-and-Discriminator-Plan.md`. |
| `absent_probe.py` | resolved (single run) | For genuinely-absent refs, distinguishes wrong-charge detection (observed m/z present, charge differs) from true detection gaps (m/z claimed by no feature). |
| `offbyone_diag.py` | resolved (single run) | Off-by-one recoverability: how many misses have a feature at ref ± k·¹³C, split by k and by carrier charge-multiplicity. |

Run with the repo's Python (no third-party deps — stdlib `csv`/`bisect` only):

```
python rust/flashlfq-core/scripts/detect_recall.py
python rust/flashlfq-core/scripts/sweep_compare.py
python rust/flashlfq-core/scripts/merge_99_rt.py
```
