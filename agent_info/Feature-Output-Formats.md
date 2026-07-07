# MS1 feature-output format: decision + spec

## Decision

**Primary format: TopFD / msDeconv `.msalign` (MS1 flavour).**

The untargeted detector's resolved MS1 features are exported to an MS1 `.msalign`
file (one deconvoluted-spectrum block per apex scan; one peak entry per feature).

Rationale in one line: it is the **only** MS1 feature/deconvolution interchange
format that mzLib can natively read, and it is a widely-used TopFD/TopPIC/FLASHDeconv
interchange format — so the output is consumable by the same library stack this port
targets.

## Why `.msalign` (evaluation of candidates)

mzLib is the C# library (github.com/smith-chem-wisc/mzLib) backing MetaMorpheus /
FlashLFQ, and the stated requirement is that the chosen format be **readable by mzLib**.
I checked mzLib's `master` source tree directly.

| Candidate | mzLib support | Evidence | Maps to our features? |
|-----------|---------------|----------|-----------------------|
| **TopFD/msDeconv `.msalign` (MS1)** | **YES — dedicated reader** | `mzLib/Readers/MsAlign/` → `MsAlign.cs`, `Ms1Align.cs`, `Ms2Align.cs`. `Ms1Align` overrides `DefaultMsnOrder => 1` and reads MS1 deconvoluted blocks into `MsDataScan[]`. | Directly: block per scan, peak = (neutral monoisotopic mass, intensity, charge). |
| ProMex `.ms1ft` | **NO** | No `ms1ft` / `MsFeature` reader anywhere in the mzLib tree (only `MassSpectrometry/DeconvolutionFeature*.cs`, which are in-memory decon types, not a `.ms1ft` file reader). `.ms1ft` is an Informed-Proteomics / LCMS-Spectator format, not mzLib. | Would map well (feature-level), but no mzLib reader → fails the requirement. |
| Dinosaur `.features.tsv` | **NO** | Not present in the mzLib tree. | n/a |
| FlashLFQ `QuantifiedPeaks` TSV | Produced by FlashLFQ, not a general MS1-feature interchange format | (the detector already emits a bespoke TSV) | Already emitted; not an interchange target. |

`.msalign` is the clear winner: it is the single candidate mzLib actually parses, and
it is broadly used across the top-down ecosystem (TopFD writes it; TopPIC, TopMG,
pTop, and FLASHDeconv read/write it).

Source URLs (verified against `master`):
- `https://github.com/smith-chem-wisc/mzLib/tree/master/mzLib/Readers/MsAlign`
- `https://raw.githubusercontent.com/smith-chem-wisc/mzLib/master/mzLib/Readers/MsAlign/MsAlign.cs`
- `https://raw.githubusercontent.com/smith-chem-wisc/mzLib/master/mzLib/Readers/MsAlign/Ms1Align.cs`

## Format spec (as parsed by mzLib `MsAlign.cs`, verified from source)

The file is: an **optional** file-level parameter header, followed by one or more
`BEGIN IONS` … `END IONS` spectrum blocks.

### File-level header (optional, we emit it)
- Region is delimited by the exact literal `##### Parameters #####` (5 hashes each
  side of `Parameters`), which appears once to open and once to close.
- Each header line inside starts with `#` and is parsed as `key : value` by
  `line.TrimStart('#').Split(':')` (kept only if it splits into exactly 2 parts).
- The reader treats the header as entirely optional: the first `BEGIN IONS`
  terminates header reading regardless.

### Spectrum block
- Delimiters: lines containing `BEGIN IONS` and `END IONS`.
- Entry-header lines are `KEY=VALUE` (split on `=`). Recognized keys (full switch in
  `MsAlign.cs`): `ID` / `SPECTRUM ID`, `FRACTION_ID`, `FILE_NAME`, `SCANS`,
  `RETENTION_TIME`, `LEVEL`, `ACTIVATION`, `MS_ONE_ID`, `MS_ONE_SCAN`, `PRECURSOR_MZ`,
  `PRECURSOR_CHARGE`, `PRECURSOR_MASS`, `PRECURSOR_INTENSITY`,
  `PRECURSOR_WINDOW_BEGIN`, `PRECURSOR_WINDOW_END`.
  - `RETENTION_TIME` is in **seconds** (mzLib divides by 60 to store minutes).
  - `LEVEL` is the MS order; if absent it falls back to the class `DefaultMsnOrder`
    (1 for `Ms1Align`). We write `LEVEL=1` explicitly.
  - No field is mandatory; a block with `SCANS`/`RETENTION_TIME`/`LEVEL` + peak lines
    parses into a valid `MsDataScan`.
- Peak lines are **tab-delimited**, exactly 3 columns, in order:
  `monoisotopicMass` `\t` `intensity` `\t` `charge`
  (mzLib: `monoMasses[i]=double.Parse(splits[0]); intensities[i]=double.Parse(splits[1]); charges[i]=int.Parse(splits[2]);`).

### What we write (MS1 flavour)
```
##### Parameters #####
#File name: <spectra basename>
#Number of MS1 scans: <n blocks>
#Software: flashlfq-rust untargeted feature detector
##### Parameters #####
BEGIN IONS
ID=<0-based block index>
SCANS=<1-based apex scan number>
RETENTION_TIME=<apex RT in seconds>
LEVEL=1
<monoMass>\t<summedIntensity>\t<charge>      # one line per feature at this scan
...                                           # peaks sorted by ascending mass
END IONS
```
Features are grouped by apex scan so each `BEGIN IONS` block is a genuine per-scan
deconvoluted MS1 spectrum (matching TopFD's block-per-scan convention and avoiding
duplicate `SCANS` numbers across blocks).

### Field mapping (detector → `.msalign`)
| `.msalign` field | Detector source (`ResolvedFeature`) |
|------------------|-------------------------------------|
| peak col 0 `monoMass` | `monoisotopic_mass` |
| peak col 1 `intensity` | `summed_intensity` |
| peak col 2 `charge` | primary (tallest) member's `refined_charge` |
| `SCANS` | primary member's `detected.apex_scan_index` + 1 |
| `RETENTION_TIME` | `apex_rt` × 60 (min → s) |

## Wiring

- Library module: `rust/flashlfq-core/src/feature_export.rs`
  - `Ms1AlignFeature` — neutral intermediate struct.
  - `write_ms1_align(writer, &[Ms1AlignFeature], source_file_name)` — core writer.
  - `write_ms1_align_file(path, …)` — path convenience wrapper.
  - `resolved_to_ms1align_features(&[ResolvedFeature])` — maps pipeline output.
- Opt-in in `examples/detect_features_tsv.rs`: set env var `MSALIGN_OUT=1` to also
  write `<out>.ms1.msalign` alongside the existing TSVs. The default TSV output is
  unchanged.
