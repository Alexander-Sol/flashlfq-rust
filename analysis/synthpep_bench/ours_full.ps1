# Run our detector on all 57 uncalibrated mzML (chosen format = mzML), 24 threads, sequential, timed.
$bin    = "F:\flashlfq-rust\rust\target\release\examples\detect_features_tsv.exe"
$data   = "D:\SyntheticPeptides_PXD001091"
$outroot= "F:\flashlfq-rust\analysis\synthpep_bench\ours"
$timing = "F:\flashlfq-rust\analysis\synthpep_bench\ours_timings.tsv"
$env:RAYON_NUM_THREADS = "16"
New-Item -ItemType Directory -Force $outroot | Out-Null
"base`tsecs`tfeatures`texit" | Set-Content $timing
$mzmls = Get-ChildItem $data -Filter "*_uncalibrated.mzML" | Sort-Object Name
$i = 0
foreach ($f in $mzmls) {
    $i++
    $base = $f.BaseName -replace "_uncalibrated$", ""
    $od = Join-Path $outroot $base
    New-Item -ItemType Directory -Force $od | Out-Null
    $out = Join-Path $od "feat.tsv"
    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    $p = Start-Process -FilePath $bin -ArgumentList @($f.FullName, $out) -NoNewWindow -Wait -PassThru `
        -RedirectStandardOutput "$od\run.out" -RedirectStandardError "$od\run.err"
    $sw.Stop()
    $secs = [math]::Round($sw.Elapsed.TotalSeconds, 1)
    $nfeat = if (Test-Path $out) { (Get-Content $out | Measure-Object -Line).Lines - 1 } else { -1 }
    $line = "{0}`t{1}`t{2}`t{3}" -f $base, $secs, $nfeat, $p.ExitCode
    Add-Content $timing $line
    Write-Output ("[{0}/{1}] {2}" -f $i, $mzmls.Count, $line)
}
Write-Output "ALL OURS DONE"
