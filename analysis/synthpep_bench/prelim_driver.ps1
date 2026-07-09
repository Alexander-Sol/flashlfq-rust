$bin  = "F:\flashlfq-rust\rust\target\release\examples\detect_features_tsv.exe"
$data = "D:\SyntheticPeptides_PXD001091"
$root = "F:\flashlfq-rust\analysis\synthpep_bench\prelim"
$env:RAYON_NUM_THREADS = "24"
$bases = @("130124_dilA_10_01", "130124_dilA_1_01")
$timing = @()
foreach ($base in $bases) {
    foreach ($fmt in @("raw", "mzml")) {
        if ($fmt -eq "raw") { $spec = Join-Path $data "$base.raw" }
        else                { $spec = Join-Path $data "${base}_uncalibrated.mzML" }
        $outdir = Join-Path $root "${base}_$fmt"
        New-Item -ItemType Directory -Force $outdir | Out-Null
        $out = Join-Path $outdir "feat.tsv"
        $log = Join-Path $outdir "run.log"
        $sw = [System.Diagnostics.Stopwatch]::StartNew()
        $p = Start-Process -FilePath $bin -ArgumentList @($spec, $out) -NoNewWindow -Wait -PassThru `
             -RedirectStandardOutput "$log.out" -RedirectStandardError "$log.err"
        $sw.Stop()
        $secs = [math]::Round($sw.Elapsed.TotalSeconds, 1)
        $nfeat = if (Test-Path $out) { (Get-Content $out | Measure-Object -Line).Lines - 1 } else { -1 }
        $line = "{0}`t{1}`t{2}`t{3}`texit={4}" -f $base, $fmt, $secs, $nfeat, $p.ExitCode
        $timing += $line
        Write-Output $line
    }
}
$timing | Set-Content (Join-Path $root "prelim_timings.tsv")
Write-Output "PRELIM DONE"
