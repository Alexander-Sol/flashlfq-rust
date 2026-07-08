# Coverage-target sweep on the LONG 120-min IonStar file.
# Usage: sweep_coverage.ps1 <coverage>   e.g. 0.90
param([Parameter(Mandatory=$true)][string]$cov)

$exe  = "F:\flashlfq-rust\rust\target\release\examples\detect_features_tsv.exe"
$raw  = "D:\PXD003881_IonStar_SpikeIn\B03_19_150304_human_ecoli_B_3ul_3um_column_95_HCD_OT_2hrs_30B_9B.raw"
$scr  = "C:\Users\Alex\AppData\Local\Temp\claude\F--flashlfq-rust\6ffafcbe-584e-46ee-af4e-b133d003ba90\scratchpad"
$tag  = ($cov -replace '\.','p')
$out  = Join-Path $scr "long_cov$tag.tsv"
$err  = Join-Path $scr "long_cov$tag.console.log"

$env:COVERAGE_TARGET = $cov
$env:DETECT_PROFILE  = "1"
# remove any stale per-stem log so this run's timing block is the only one
$log = $out -replace '\.tsv$','.log'
if (Test-Path $log) { Remove-Item $log }

$sw = [System.Diagnostics.Stopwatch]::StartNew()
$p = Start-Process -FilePath $exe -ArgumentList @($raw, $out) -NoNewWindow -Wait -PassThru -RedirectStandardError $err -RedirectStandardOutput ($err + '.out')
$sw.Stop()
"COVERAGE_TARGET=$cov  exit=$($p.ExitCode)  wall=$([math]::Round($sw.Elapsed.TotalSeconds,2))s  out=$out"
