# Run Dinosaur 1.2.0 on all 57 uncalibrated mzML, concurrency=16, sequential, timed.
$jar    = "C:\Users\Alex\Desktop\Dinosaur\Dinosaur-1.2.0.free.jar"
$data   = "D:\SyntheticPeptides_PXD001091"
$outdir = "F:\flashlfq-rust\analysis\synthpep_bench\dino"
$timing = "F:\flashlfq-rust\analysis\synthpep_bench\dino_timings.tsv"
New-Item -ItemType Directory -Force $outdir | Out-Null
"base`tsecs`tfeatures`texit" | Set-Content $timing
$mzmls = Get-ChildItem $data -Filter "*_uncalibrated.mzML" | Sort-Object Name
$i = 0
foreach ($f in $mzmls) {
    $i++
    $base = $f.BaseName -replace "_uncalibrated$", ""
    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    $p = Start-Process -FilePath "java" -ArgumentList @(
            "-Xmx12g", "-jar", $jar, "--verbose", "--profiling",
            "--concurrency=16", "--outDir=$outdir", "--outName=$base", "--force", $f.FullName
        ) -NoNewWindow -Wait -PassThru `
        -RedirectStandardOutput "$outdir\$base.console.out" `
        -RedirectStandardError  "$outdir\$base.console.err"
    $sw.Stop()
    $secs = [math]::Round($sw.Elapsed.TotalSeconds, 1)
    $feat = "$outdir\$base.features.tsv"
    $nfeat = if (Test-Path $feat) { (Get-Content $feat | Measure-Object -Line).Lines - 1 } else { -1 }
    $line = "{0}`t{1}`t{2}`t{3}" -f $base, $secs, $nfeat, $p.ExitCode
    Add-Content $timing $line
    Write-Output ("[{0}/{1}] {2}" -f $i, $mzmls.Count, $line)
}
Write-Output "ALL DINOSAUR DONE"
