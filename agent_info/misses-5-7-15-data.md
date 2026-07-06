# Data report: misses #5, #7, #15 — envelope reality vs. pipeline mono placement

Read-only data extraction for manual investigation of a suspected faulty charge-state-merging
problem. No code was run; all numbers are pulled from files on disk.

**Conventions**
- Neutral mass of a peak: `M = mz*z − z*1.007276466879`. Isotope index vs a reference mono:
  `k = (M − ref_mono) / 1.00335483810`. Uncalibrated raw data runs ~+10 ppm high, so real peaks
  land at an integer k plus a small positive residual (~+0.02–0.05 here).
- "mode" = most-abundant peak inside the analysis window. "frac_mode" = intensity / mode intensity.
- Reference `Peak MZ` in AllQuantifiedPeaks is MetaMorpheus's *most-abundant* peak (apex), **not**
  the monoisotope. For all three peptides that reference `Peak MZ` sits at the **+1** isotope.

**Column disambiguation (reference file, verified from header):** column [7] = `Precursor Charge`
(the MS2 precursor), column [14] = `Peak Charge` (the quantified MS1 peak). This report uses
`Peak Charge` / `Peak MZ` / `Peak RT Apex` throughout.

## TL;DR — the pattern

All three resolved (merged) features report a monoisotopic mass that is **off by one isotope**,
and in every case a charge that had the mono placed *correctly* fails to survive as its own
resolved row:

| case | ref mono | z | our resolved mono (k) | direction | true mono in envelope | peak where we put the mono? |
|---|---|---|---|---|---|---|
| #5  | 2475.1013 | 3 | 2476.1274 (**k=+1.02**) | one isotope **high** | weak: 0.21 of mode (composite), **absent at apex** | **YES** (that's the +1 / mode peak) |
| #7  | 2894.4749 | 3 | 2893.5072 (**k=−0.96**) | one isotope **low**  | strong: **0.65** of mode (composite), 0.62 (apex) | **NO** (0.005 of mode ≈ noise) |
| #15 | 2905.4199 | 4 | 2904.4358 (**k=−0.98**) | one isotope **low**  | strong: **0.62** of mode (composite), **is the mode at apex** | **NO** (nothing there) |

**#7 and #15 exactly match the human's observation** ("mono at ~0.6 of the most-abundant peak,
NO peak at the position we reported the mono") and **contradict** the earlier "tiny mono /
averagine-ambiguity" characterization — their monos are strong, well-resolved peaks.
**#5 is a different failure** (mono reported +1 *high*, onto the real most-abundant peak; its true
mono genuinely is weak). Do not lump #5 in with #7/#15.

---

# Peptide #5 — MVNNGHSFNVEYDDSQDKAVLK, z=3, RT 12.660

## 1. Reference truth
Matched modform for ref mono **2475.1013**: `MVN[Common Artifact:Ammonia loss on N]N[Common Artifact:Ammonia loss on N]GHSFNVEYDDSQDKAVLK`,
peptide monoisotopic mass **2475.101254**, Peak Charge **3**, Peak MZ **826.3815** (this is the +1
isotope), Peak RT Apex **12.6601** (start 12.6106, end 12.6601), Num Charge States Observed 2.

Theoretical **mono** m/z for this modform (M=2475.101254):

| z | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| mono m/z | 2476.10853 | 1238.55790 | **826.04103** | 619.78259 | 496.02753 | 413.52415 |

All 12 reference entries for base sequence `MVNNGHSFNVEYDDSQDKAVLK` (Peak Charge / Peak MZ / Peak RT Apex / pep mono / mod):

| pep mono | Peak Chg | Peak MZ | Peak RT Apex | modform |
|---|---|---|---|---|
| 2509.154352 | 4 | 628.5528 | 12.269 | (unmod) |
| 2510.138368 | 4 | 628.7986 | 12.514 | Deamidation on N |
| 2510.138368 | 4 | 628.7985 | 12.414 | Deamidation on N |
| 2525.149267 | 4 | 632.5519 | 12.103 | Oxidation on M |
| 2526.133282 | 4 | 632.7975 | 12.316 | Oxidation on M + Deamidation on N |
| 2552.160166 | 3 | 852.0694 | 12.960 | Carbamyl on K (2 loci) |
| 2511.122383 | 4 | 629.0446 | 12.627 | 2× Deamidation on N |
| 2509.154352 | 3 | 837.7349 | 13.284 | (unmod) |
| 2526.133282 | 4 | 632.7979 | 12.219 | Oxidation on M + Deamidation on N |
| **2475.101254** | **3** | **826.3815** | **12.660** | **2× Ammonia loss on N  ← THIS ROW** |
| 2536.117632 | 3 | 846.7224 | 13.871 | Ammonia loss on N + Carboxylation on E |
| 2548.094249 | 4 | 638.2857 | 12.531 | Deamidation on N + Potassium on E |

## 2. What we detected (±3 Da of 2475.1013, ±0.25 min of RT 12.660; sorted by intensity)

| charge | det mono | k vs ref | mono m/z | RT start / apex / end | summed int | nIso | nPk | det score |
|---|---|---|---|---|---|---|---|---|
| 4 | 2474.1226 | **−0.975** | 619.5379 | 12.562 / 12.660 / 12.660 | 1.157e8 | 7 | 27 | 5.36e7 |
| 3 | 2476.1274 | **+1.023** | 826.3831 | 12.611 / 12.660 / 12.660 | 6.384e7 | 4 | 12 | 4.05e7 |
| 3 | 2474.1416 | −0.956 | 825.7211 | 12.467 / 12.482 / 12.482 | 2.680e7 | 6 | 11 | 1.66e7 |
| 4 | 2475.1443 | +0.043 | 619.7934 | 12.456 / 12.467 / 12.467 | 1.100e7 | 5 | 6 | 9.03e6 |
| 4 | 2476.2959 | +1.191 | 620.0812 | 12.767 / 12.767 / 12.779 | 3.330e6 | 4 | 5 | 2.02e6 |

Note the **highest-intensity** detection (z4, 1.157e8) has a *detector* mono at k=−0.975, but it
**refines to the correct mono** (next section). The z3 detection at k=+1.02 keeps its wrong (+1)
placement.

## 3. What we refined (matched by Detector Mono Mass)

| det mono | → refined mass | k vs ref | z | apex | int | decon | 
|---|---|---|---|---|---|---|
| 2474.1226 | **2475.1259** | **+0.025 (CORRECT)** | 4 | 12.660 | 1.157e8 | 0.590 |
| 2476.1274 | 2476.1274 | +1.023 | 3 | 12.660 | 6.384e7 | 0.946 |
| 2475.1443 | 2474.1409 | −0.957 | 4 | 12.467 | 1.100e7 | 0.987 |

The dominant z4 feature refines to **2475.1259 (k≈0, correct mono)**.

## 4. What we resolved (the merge)

| mono (k) | charge states | nMembers | X-chg support | primary | mono m/z (primary) | per-charge det m/z | apex | int |
|---|---|---|---|---|---|---|---|---|
| **2476.1274 (k=+1.023)** | **3;4** | **2** | 1 | 4 | 620.039 | z3:826.7175 · z4:619.7888 | 12.660 | **1.796e8** |
| 2474.1388 (k=−0.959) | 3;4 | 5 | 2 | 3 | 825.720 | z3:826.0556 · z4:620.2892 | 12.482 | 9.326e7 |
| 2476.2959 (k=+1.191) | 4 | 1 | 1 | 4 | 620.081 | z4:620.3321 | 12.767 | 3.330e6 |

**The answer for this peptide (RT 12.660) is `2476.1274`, k=+1.02 — the mono is reported one
isotope too high.** Its intensity (1.796e8 ≈ 6.38e7 + 1.16e8) shows it **merged the correct-mono z4
(refined 2475.1259, k≈0) with the wrong +1 z3 (2476.1274)** and reported the **z3's +1 value** as the
consensus mono. The correctly-placed k=0 mono (2475.126) **does not survive as its own resolved
row at RT 12.66** (a scan of all resolved rows near 2475.13 finds only off-RT features at 13.49,
15.03, 15.27, 16.92, 18.40, 19.37, 19.77 — none at 12.66). So the merge overwrote a correct
placement.

## 5. Observed envelope (z=3). ref mono m/z = **826.04104**; our reported mono m/z (2476.1274 @ z3) = **826.38308** (k=+1.02)

**observed_composite** — mode = 826.71753 (k=+2.02), int 1.131e7:

| m/z | intensity | k vs ref | frac_mode | |
|---|---|---|---|---|
| 826.05042 | 2.383e6 | +0.028 | 0.211 | ← REF MONO (k0) |
| 826.38153 | 9.938e6 | +1.018 | 0.879 | ← OUR reported mono |
| 826.71753 | 1.131e7 | +2.023 | 1.000 | mode |
| 827.05206 | 5.561e6 | +3.023 | 0.492 | |
| 827.38794 | 2.471e5 | +4.027 | 0.022 | |
| 827.43463 | 9.518e6 | +4.167 | 0.842 | co-eluting different species |
| 828.15405 | 3.842e5 | +6.318 | 0.034 | |

**observed_apex** — mode = 826.71753 (k=+2.02), int 1.850e7:

| m/z | intensity | k vs ref | frac_mode | |
|---|---|---|---|---|
| 826.38153 | 1.626e7 | +1.018 | 0.879 | ← OUR reported mono |
| 826.71753 | 1.850e7 | +2.023 | 1.000 | mode |
| 827.05206 | 9.099e6 | +3.023 | 0.492 | |

**Peak at ref mono 826.041?** composite YES @ 0.211 of mode; **apex: NONE** (true mono absent in the
apex scan). **Peak at our reported mono 826.383?** YES — 0.879 of mode (it is the +1 isotope /
second-tallest peak, and equals the reference's own `Peak MZ`). 

**#5 verdict:** we placed the mono **+1 too high onto the real +1 peak**. The true mono is genuinely
weak (0.21 composite, absent at apex) — so for *this* peptide the "weak mono" story holds. This is
NOT the same failure as #7/#15.

---

# Peptide #7 — [Carbamyl]TLNFNAEGEPELLMLANWRPAQPLK, z=3, RT 16.406

Base sequence is **`TLNFNAEGEPELLMLANWRPAQPLK`** (25 residues; the task label was truncated).

## 1. Reference truth
Ref mono **2894.4749** corresponds to two isobaric modforms that both co-elute here:
- `[Common Artifact:Carbamyl on X]TLNFNAEGEPELLMLANWRPAQPLK` — pep mono 2894.474898, Peak Chg 3, Peak MZ 966.1754, **Peak RT Apex 16.772**
- `TLNFNAEGEPELLM[Common Artifact:Carbamyl on M]LANWRPAQPLK` — pep mono 2894.474898, Peak Chg 3, Peak MZ 966.1755, **Peak RT Apex 16.406  ← the case-07 reference (RT 16.406)**

Theoretical **mono** m/z for M=2894.474898:

| z | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| mono m/z | 2895.48217 | 1448.24473 | **965.83224** | 724.62600 | 579.90226 | 483.41976 |

All 21 reference rows for base `TLNFNAEGEPELLMLANWRPAQPLK` (pep mono / Peak Chg / Peak MZ / Peak RT Apex / mod):

| pep mono | Chg | Peak MZ | RT Apex | modform |
|---|---|---|---|---|
| 2851.469085 | 3 | 951.8400 | 16.552 | (unmod) |
| 2894.500051 | 2 | 1448.7594 | 16.786 | Citrullination on R + Trimethylation on K |
| 2867.463999 | 3 | 957.1714 | 16.358 | Oxidation on M |
| 2910.469813 | 3 | 971.5070 | 16.695 | Carbamyl on X + Oxidation on M |
| 2916.456843 | 3 | 973.5031 | 16.786 | Carbamyl on X + Sodium on E |
| 2908.490548 | 4 | 728.3877 | 15.713 | Carbamyl on M + Methylation on R |
| 2965.500779 | 3 | 989.8551 | 16.503 | Dimethyl R + Malonyl/Hydroxybutyryl K |
| **2894.474898** | **3** | **966.1754** | **16.772** | **Carbamyl on X** |
| 2873.451029 | 3 | 959.1674 | 16.536 | Sodium on E |
| 2834.442536 | 3 | 946.1654 | 16.552 | Ammonia loss on N |
| 2889.445944 | 3 | 964.4990 | 16.358 | Sodium on E + Oxidation on M |
| 2883.458914 | 3 | 962.5037 | 16.438 | Hydroxylation on N + P |
| 2937.480712 | 2 | 1470.2621 | 18.186 | Carbamyl on X + Carbamyl on K |
| 2932.430780 | 3 | 978.8228 | 16.802 | Carbamyl on X + Potassium on E |
| 2889.445944 | 3 | 964.4984 | 16.342 | Sodium on E + Oxidation on M |
| **2894.474898** | **3** | **966.1755** | **16.406** | **Carbamyl on M  ← case-07 reference** |
| 2905.419881 | 4 | 727.6176 | 16.342 | Potassium on E + Oxidation on M (= peptide #15) |
| 2851.469085 | 3 | 951.8438 | 18.640 | (unmod) |
| 2851.469085 | 3 | 951.8411 | 18.252 | (unmod) |
| 2851.469085 | 3 | 951.8414 | 19.239 | (unmod) |
| 2851.469085 | – | – | – | (unmod, not quantified) |

## 2. What we detected near RT 16.4–16.8 (±4 Da; the reference RT 16.406 and the co-eluting 16.786 apex)

Key rows (charge / det mono / k / mono m/z / apex / int / nIso / nPk / score):

| chg | det mono | k | mono m/z | apex | int | nIso | nPk | score |
|---|---|---|---|---|---|---|---|---|
| 2 | 2894.5009 | **+0.026 (CORRECT)** | 1448.2577 | 16.786 | 6.687e10 | 9 | 56 | 3.93e10 |
| 3 | 2894.5014 | **+0.026 (CORRECT)** | 965.8411 | 16.786 | 5.974e10 | 9 | 50 | 3.27e10 |
| 3 | 2893.5072 | **−0.965** | 965.5097 | 16.662 | 1.748e9 | 9 | 31 | 6.77e8 |
| 4 | 2893.4998 | −0.972 | 724.3822 | 16.786 | 8.725e8 | 8 | 19 | 4.91e8 |
| 3 | 2894.8631 | +0.387 | 965.9616 | 16.230 | 5.913e7 | 3 | 3 | 3.89e7 |
| 3 | 2891.8973 | −2.569 | 964.9730 | 16.374 | 2.345e7 | 3 | 4 | 1.42e7 |
| 3 | 2893.1219 | −1.348 | 965.3813 | 16.342 | 5.289e6 | 2 | 3 | 2.67e6 |

At the **16.786 apex** the detector places the mono **correctly** (k=+0.026) on both z2 and z3.
The **z3 at RT 16.662, det mono 2893.5072 (k=−0.965)** is the off-by-one feature: its
`Detected m/z (primary) = 965.844` — i.e. it detected the **true mono peak (965.84, k0)** but labeled
it the +1 isotope and stepped the reported mono **down** to 965.510.

## 3. What we refined (near RT 16.4–16.8)

| refined mass | k | z | det mono (k) | apex | int | decon |
|---|---|---|---|---|---|---|
| 2894.5009 | +0.026 | 2 | 2894.5009 (+0.026) | 16.786 | 6.687e10 | 0.996 |
| 2894.5014 | +0.026 | 3 | 2894.5014 (+0.026) | 16.786 | 5.974e10 | 0.998 |
| **2893.5072** | **−0.965** | 3 | 2893.5072 (−0.965) | 16.662 | 1.748e9 | 0.997 |
| 2894.5031 | +0.028 | 4 | 2893.4998 (−0.972) | 16.786 | 8.725e8 | 0.905 |
| 2892.8564 | −1.613 | 3 | 2894.8631 (+0.387) | 16.230 | 5.913e7 | 0.526 |
| 2893.1219 | −1.348 | 3 | 2893.1219 (−1.348) | 16.342 | 5.289e6 | 0.158 |

## 4. What we resolved (the merge)

| mono (k) | chg states | nMem | primary | mono m/z (prim) | det m/z (prim) | per-charge det m/z | apex | int |
|---|---|---|---|---|---|---|---|---|
| **2894.5012 (k=+0.026, CORRECT)** | **2;3;4** | **3** | 2 | 1448.258 | 1448.759 | z2:1448.7594 · z3:966.1755 · z4:724.6331 | **16.786** | **1.275e11** |
| **2893.5072 (k=−0.965)** | **3** | **1** | 3 | **965.510** | 965.844 | z3:965.8441 | **16.662** | 1.748e9 |
| 2896.4732 (k=+1.99) | 3 | 1 | 3 | 966.498 | 967.502 | z3:967.5017 | 16.695 | 6.223e7 |
| 2892.8564 (k=−1.61) | 3 | 1 | 3 | 965.293 | 966.296 | z3:966.2961 | 16.230 | 5.913e7 |
| 2893.4865 (k=−0.99) | 3 | 1 | 3 | 965.503 | 965.837 | z3:965.8372 | 16.134 | 1.070e7 |
| 2893.1219 (k=−1.35) | 3 | 1 | 3 | 965.381 | 965.716 | z3:965.7157 | 16.342 | 5.289e6 |

**Two things live in this region.** (a) A **correctly-placed merged feature** `2894.5012 (k≈0)` at
**RT 16.786** (merges z2;3;4, int 1.275e11) — this is the pipeline's good detection, matched to the
RT-16.77/16.79 modform(s). (b) At/near the **reference RT 16.406** the pipeline produced **only
off-by-one-low fragments**; the nearest z3 (the reference charge) is **`2893.5072`, k=−0.96, reported
mono m/z 965.510**. Every one of these near-RT-16.4 z3 rows sits at k ≈ −1 or lower and reports a
mono m/z (965.51, 965.29, 965.50, 965.38) where the spectrum is essentially empty, even though each
one's *detected* peak (965.84, 966.30, 965.84, 965.72) sits on a real isotope. This is a
detector/refiner off-by-one that then propagates into the resolved mono.

## 5. Observed envelope (z=3). ref mono m/z = **965.83224**; our reported mono m/z (2893.5072 @ z3) = **965.50968** (k=−0.96)

**observed_composite** — mode = 966.17554 (k=+1.03), int 5.746e9 (only peaks ≥1e7 shown; the k≈0
and k≈+1..+6 isotopes dominate, everything else is <0.005 of mode):

| m/z | intensity | k vs ref | frac_mode | |
|---|---|---|---|---|
| 965.50580 | 2.600e7 | −0.976 | **0.0045** | ← OUR reported mono (essentially noise) |
| 965.84100 | 3.709e9 | +0.026 | **0.6455** | ← REF MONO (k0) |
| 966.17554 | 5.746e9 | +1.026 | 1.0000 | mode |
| 966.50970 | 4.854e9 | +2.026 | 0.8447 | |
| 966.84357 | 3.037e9 | +3.024 | 0.5284 | |
| 967.17725 | 1.588e9 | +4.022 | 0.2764 | |
| 967.50946 | 9.295e8 | +5.015 | 0.1618 | |
| 967.84444 | 3.328e8 | +6.016 | 0.0579 | |

**observed_apex** — mode = 966.17554 (k=+1.03), int 9.540e9:

| m/z | intensity | k vs ref | frac_mode | |
|---|---|---|---|---|
| 965.50269 | 2.414e7 | −0.985 | **0.0025** | ← OUR reported mono (noise) |
| 965.84131 | 5.870e9 | +0.027 | **0.6153** | ← REF MONO (k0) |
| 966.17554 | 9.540e9 | +1.026 | 1.0000 | mode |
| 966.50970 | 8.058e9 | +2.026 | 0.8447 | |
| 966.84381 | 4.839e9 | +3.025 | 0.5072 | |
| 967.17792 | 2.585e9 | +4.024 | 0.2709 | |

**Peak at ref mono 965.832?** YES, strong — **0.646 (composite) / 0.615 (apex) of the mode**.
**Peak at our reported mono 965.510?** effectively **NO** — only 2.6e7 / 2.4e7, i.e. **0.005 / 0.003
of the mode** (M−1 baseline, not a real isotope). Composite and apex **agree**.

**#7 verdict:** matches the human exactly. True mono is a strong peak (~0.6 of mode); we reported the
mono one isotope **low**, onto an m/z with no real peak. "Tiny mono / averagine ambiguity" is wrong
for #7.

---

# Peptide #15 — TLNFNAE[Potassium on E]GEPELLM[Oxidation on M]LANWRPAQPLK, z=4, RT 16.342

Same base peptide `TLNFNAEGEPELLMLANWRPAQPLK`; modform = **Potassium on E + Oxidation on M**.

## 1. Reference truth
Ref mono **2905.4199**: `TLNFNAE[Metal:Potassium on E]GEPELLM[Common Variable:Oxidation on M]LANWRPAQPLK`,
pep mono **2905.419881**, Peak Charge **4**, Peak MZ **727.6176** (the +1 isotope), Peak RT Apex
**16.342** (start = end = 16.342), Num Charge States Observed 3. (Full 21-row reference table for the
base peptide is given under #7 above.)

Theoretical **mono** m/z for M=2905.419881:

| z | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| mono m/z | 2906.42716 | 1453.71722 | 969.48057 | **727.36225** | 582.09125 | 485.24392 |

## 2. What we detected (±3 Da of 2905.4199, ±0.25 min of RT 16.342; sorted by intensity)

| chg | det mono | k vs ref | mono m/z | RT start / apex / end | int | nIso | nPk | score |
|---|---|---|---|---|---|---|---|---|
| 4 | 2904.4360 | **−0.981** | 727.1163 | 16.325 / 16.358 / 16.374 | 6.035e8 | 7 | 19 | 3.11e8 |
| 3 | 2905.4361 | **+0.016 (CORRECT)** | 969.4860 | 16.342 / 16.358 / 16.374 | 2.947e8 | 5 | 12 | 1.89e8 |
| 5 | 2904.4342 | **−0.982** | 581.8941 | 16.342 / 16.358 / 16.358 | 7.005e7 | 5 | 6 | 6.98e7 |
| 2 | 2902.5541 | −2.856 | 1452.2843 | 16.243 / 16.253 / 16.253 | 4.018e7 | 4 | 4 | 3.10e7 |
| 4 | 2906.4295 | +1.006 | 727.6146 | 16.389 / 16.389 / 16.406 | 3.261e7 | 5 | 8 | 1.87e7 |
| 3 | 2904.2815 | −1.135 | 969.1011 | 16.230 / 16.230 / 16.243 | 2.052e7 | 2 | 2 | 1.86e7 |
| 3 | 2906.9176 | +1.493 | 969.9798 | 16.268 / 16.268 / 16.295 | 1.452e7 | 2 | 4 | 6.13e6 |
| 3 | 2906.4580 | +1.035 | 969.8266 | 16.406 / 16.438 / 16.438 | 1.117e7 | 3 | 4 | 7.48e6 |

The **dominant** detection is **z4 at k=−0.981** (int 6.0e8). The **z3 detection has the mono correct
(k=+0.016)** but is smaller (2.9e8). The z5 is also at k=−0.982.

## 3. What we refined (matched by Detector Mono Mass)

| det mono | → refined mass | k vs ref | z | apex | int | decon |
|---|---|---|---|---|---|---|
| 2904.4360 | 2904.4360 | **−0.981** | 4 | 16.358 | 6.035e8 | 0.870 |
| 2905.4361 | 2905.4361 | **+0.016 (CORRECT)** | 3 | 16.358 | 2.947e8 | 0.561 |
| 2904.4342 | 2904.4342 | −0.982 | 5 | 16.358 | 7.005e7 | 0.636 |
| 2906.4295 | 2905.4261 | +0.006 | 4 | 16.389 | 3.261e7 | 0.698 |
| 2904.2815 | 2904.2815 | −1.135 | 3 | 16.230 | 2.052e7 | 0.313 |

Refinement leaves the dominant z4 at k=−0.98 and the z3 at k≈0 (correct); they disagree by one
isotope.

## 4. What we resolved (the merge)

| mono (k) | chg states | nMem | X-chg | primary | mono m/z (prim) | per-charge det m/z | apex | int |
|---|---|---|---|---|---|---|---|---|
| **2904.4358 (k=−0.981)** | **3;4;5** | **4** | 2 | 4 | **727.116** | z3:969.8204 · z4:727.3671 · z5:582.0948 | 16.358 | **1.001e9** |
| 2908.5119 (k=+3.08) | 3 | 1 | 1 | 3 | 970.511 | z3:970.8457 | 16.487 | 1.347e8 |
| 2902.5541 (k=−2.86) | 2 | 1 | 1 | 2 | 1452.284 | z2:1452.7860 | 16.253 | 4.018e7 |
| 2904.2815 (k=−1.14) | 3 | 1 | 1 | 3 | 969.101 | z3:969.4355 | 16.230 | 2.052e7 |
| 2904.9109 (k=−0.51) | 3 | 1 | 1 | 3 | 969.311 | z3:970.3143 | 16.268 | 1.452e7 |
| 2907.4614 (k=+2.04) | 3 | 1 | 1 | 3 | 970.161 | z3:970.1611 | 16.438 | 1.117e7 |

**The answer for this peptide is `2904.4358`, k=−0.98 — mono reported one isotope too low.** It
**merges z3;4;5 (4 members)** and picks **z4 as primary**. Note the per-charge detected m/z include
**z4:727.3671 = the true mono (727.362)** and **z5:582.0948 = the true mono** — the detector's
most-abundant peaks are on the real mono, but the reported mono is stepped down one to **727.116**.
The z3 member merged in is the **+1** one (2906.458, det m/z 969.8266 ≈ 969.8204), **not** the
correctly-placed z3 (2905.4361, k≈0). The correct z3 k=0 feature **does not survive as its own
resolved row at RT 16.34** (nearest resolved rows to 2905.44 are at RT 18.45, 14.43, 7.67 — different
elutions). So the correct z3 placement was outvoted/absorbed by the dominant z4 at k=−1.

## 5. Observed envelope (z=4). ref mono m/z = **727.36225**; our reported mono m/z (2904.4358 @ z4) = **727.11623** (k=−0.98)

**observed_composite** — mode = 727.61755 (k=+1.02), int 6.116e7:

| m/z | intensity | k vs ref | frac_mode | |
|---|---|---|---|---|
| 727.31921 | 3.607e6 | −0.172 | 0.059 | |
| 727.36719 | 3.810e7 | +0.020 | **0.623** | ← REF MONO (k0) |
| 727.61755 | 6.116e7 | +1.018 | 1.000 | mode |
| 727.86823 | 5.807e7 | +2.017 | 0.950 | |
| 728.11829 | 3.412e7 | +3.014 | 0.558 | |
| 728.36877 | 2.177e7 | +4.013 | 0.356 | |
| 728.61816 | 1.607e7 | +5.007 | 0.263 | |
| 729.12109 | 3.120e6 | +7.012 | 0.051 | |

*(no peak anywhere near our reported mono 727.116 — nearest is the small 727.319 at k=−0.17)*

**observed_apex** — mode = 727.36713 (k=+0.02) **= the true mono itself**, int 1.199e8:

| m/z | intensity | k vs ref | frac_mode | |
|---|---|---|---|---|
| 727.36713 | 1.199e8 | +0.019 | 1.000 | ← REF MONO (k0) = mode |
| 727.61700 | 9.544e7 | +1.016 | 0.796 | |
| 727.86847 | 8.284e7 | +2.018 | 0.691 | |
| 728.36877 | 3.965e7 | +4.013 | 0.331 | |
| 728.62024 | 3.246e7 | +5.015 | 0.271 | |

**Peak at ref mono 727.362?** YES — **0.623 of mode (composite); it IS the mode at apex (1.000)**.
**Peak at our reported mono 727.116?** **NONE** (composite and apex agree — nothing there).

**#15 verdict:** matches the human. True mono is strong (0.62 composite, the tallest peak at apex);
we reported the mono one isotope **low**, at an m/z with no peak at all. Driven by charge-state
merging: the dominant z4 (and z5) sat one isotope low and outvoted the z3 that had the mono correct.

---

# Cross-cutting notes for the charge-merge hypothesis

1. **Two distinct off-by-one failures are being conflated.** #7 and #15 report the mono **one
   isotope LOW** onto an empty m/z while the true mono is a strong peak (~0.6 of mode) — exactly the
   human's observation. #5 reports the mono **one isotope HIGH** onto the real +1/most-abundant peak,
   and its true mono genuinely is weak (0.21 composite, absent at apex). The "tiny mono /
   averagine-ambiguity" story is only defensible for #5; it is **wrong** for #7 and #15.

2. **In all three, a charge that placed the mono correctly failed to survive.** #5: z4 refined to
   the correct 2475.126 (k≈0) but the merge reported the z3's +1 value 2476.127. #15: z3 refined to
   the correct 2905.436 (k≈0) but the merge reported the z4's −1 value 2904.436. #7: z2/z3 at the
   16.786 apex are correct (k≈0, huge), but the reference-RT (16.406) region only yields the z3
   off-by-one 2893.507 (k≈−1). The correct-mono placements do not appear as their own resolved rows
   at the reference RT.

3. **The detector's most-abundant peak is usually ON the true mono, but gets labeled +1.** #7's
   off-by-one feature has `Detected m/z (primary) = 965.844` = true mono; #15's answer carries
   `z4:727.367` and `z5:582.095` = true mono, yet the *reported* mono is stepped one isotope down.
   The mis-step is in isotope indexing / mono assignment, and the cross-charge merge then locks in
   the wrong index by choosing the highest-intensity (mis-indexed) charge as primary.

4. **Composite vs. apex agree** on the key ratios for #7 and #15 (mono ≈ 0.6 of mode; nothing at our
   reported mono). For #5 they differ only in that the weak true mono is present at 0.21 in the
   composite but drops out entirely in the single apex scan. No composite-only artifact is inflating
   or hiding the mono in the disputed cases.

*(Source files: reference `AllQuantifiedPeaks.tsv`; pipeline `untargeted_cov90_shiftapex.detected.tsv`
/ `.refined.tsv` / `.tsv`; envelopes `analysis/misses/obo_series.tsv`, `obo_meta.tsv`, `obo_scans.tsv`.)*
