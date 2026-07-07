# Reproduces the medium + long detector runs behind run_timings.md (TASK 2).
# Run from F:\flashlfq-rust\rust after: cargo build --release --example detect_features_tsv
# DETECT_PROFILE=1 enables sub-stage detail on stderr. No reference arg is passed —
# recall is scored separately in Python (score_recall.py) because neither the medium
# Ngly_Generic nor the long QuantifiedPeaks table matches the example's built-in
# compare_to_reference column schema (AllQuantifiedPeaks).

$exe = "F:\flashlfq-rust\rust\target\release\examples\detect_features_tsv.exe"
$out = "F:\flashlfq-rust\analysis\longer_gradients"
$env:DETECT_PROFILE = "1"

# MEDIUM ~65-min (D_Morgen_Glyco glyco)
& $exe "D:\D_Morgen_Glyco\RawFiles\HFX_MB_14751_5_02062022.raw" "$out\medium_features.tsv"

# LONG ~120-min (IonStar spike-in)
& $exe "D:\PXD003881_IonStar_SpikeIn\B03_19_150304_human_ecoli_B_3ul_3um_column_95_HCD_OT_2hrs_30B_9B.raw" "$out\long_features.tsv"
