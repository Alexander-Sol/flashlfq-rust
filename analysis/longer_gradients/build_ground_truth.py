#!/usr/bin/env python
"""Format adapters -> peak-level ground-truth tables for the medium (65-min D_Morgen_Glyco)
and long (120-min IonStar) gradients.

Normalized schema (one row per expected peak):
    mono_mass   monoisotopic neutral mass (Da)
    mz          monoisotopic m/z at `charge`
    charge      precursor / peak charge
    rt          retention time (min) -- SEE per-case note on what RT means
    intensity   peak intensity (long only; blank for medium)
    detection   detection type (long only: MSMS | MBR)

MEDIUM (Ngly_Generic PSM table, 13,179 rows):
  * dedupe to unique (Full Sequence + Precursor Charge) -> peak-level truth
  * mz derived from Peptide Monoisotopic Mass + charge (proton mass)
  * rt = the PSM's Scan Retention Time -- NOT the chromatographic apex, so recall
    scoring against this must use a LENIENT RT delta.

LONG (FlashLFQ QuantifiedPeaks, 504,784 pooled rows):
  * filter to this raw's File Name
  * drop Decoy Peptide == True
  * split MSMS-detected (primary truth) vs MBR (match-between-runs) via Peak Detection Type
  * rt = Peak RT Apex (a real apex), mz = Peak MZ, charge = Peak Charge
"""
import csv, sys

PROTON = 1.007276466812

MED_REF = r"D:\D_Morgen_Glyco\Byonic and MSFragger ID\ConvertedBionicResults\HFX_MB_14751_5_02062022.raw_20230206_WC_Ngly_Generic.tsv"
LONG_REF = r"D:\PXD003881_IonStar_SpikeIn\MM_Master_noMBR_Norm\Task1-SearchTask\Individual File Results\B03_19_150304_human_ecoli_B_3ul_3um_column_95_HCD_OT_2hrs_30B_9B_QuantifiedPeaks.tsv"
LONG_FILE_NAME = "B03_19_150304_human_ecoli_B_3ul_3um_column_95_HCD_OT_2hrs_30B_9B"

OUT_DIR = r"F:\flashlfq-rust\analysis\longer_gradients"


def mono_mz(mass, z):
    return (mass + z * PROTON) / z


def build_medium():
    rows = {}
    n_raw = 0
    with open(MED_REF, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            n_raw += 1
            try:
                mass = float(row["Peptide Monoisotopic Mass"])
                z = int(row["Precursor Charge"])
                rt = float(row["Scan Retention Time"])
            except (ValueError, KeyError):
                continue
            key = (row["Full Sequence"], z)
            # keep first-seen RT (they cluster; median would need buffering all)
            if key not in rows:
                rows[key] = (mass, z, rt)
    out = OUT_DIR + r"\medium_ground_truth.tsv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["mono_mass", "mz", "charge", "rt", "intensity", "detection"])
        for (mass, z, rt) in sorted(rows.values()):
            w.writerow([f"{mass:.5f}", f"{mono_mz(mass, z):.5f}", z, f"{rt:.4f}", "", "PSM"])
    print(f"MEDIUM: {n_raw} PSM rows -> {len(rows)} unique (FullSeq+charge) peaks -> {out}")
    return len(rows)


def build_long():
    msms, mbr, decoy = [], [], 0
    n_file = 0
    with open(LONG_REF, newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            if row["File Name"] != LONG_FILE_NAME:
                continue
            n_file += 1
            if row["Decoy Peptide"].strip().lower() == "true":
                decoy += 1
                continue
            try:
                mass = float(row["Peptide Monoisotopic Mass"])
                z = int(row["Peak Charge"])
                rt = float(row["Peak RT Apex"])
            except (ValueError, KeyError):
                continue
            mz = row["Peak MZ"]
            inten = row["Peak intensity"]
            det = row["Peak Detection Type"].strip()
            rec = (f"{mass:.5f}", mz, z, f"{rt:.4f}", inten, det)
            if det == "MSMS":
                msms.append(rec)
            else:
                mbr.append(rec)
    for tag, data in (("msms", msms), ("all", msms + mbr)):
        out = OUT_DIR + rf"\long_ground_truth_{tag}.tsv"
        with open(out, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["mono_mass", "mz", "charge", "rt", "intensity", "detection"])
            for rec in data:
                w.writerow(rec)
        print(f"LONG {tag}: {len(data)} peaks -> {out}")
    print(f"LONG: {n_file} rows for {LONG_FILE_NAME}; {decoy} decoy dropped; "
          f"{len(msms)} MSMS, {len(mbr)} MBR non-decoy")
    return len(msms), len(mbr)


if __name__ == "__main__":
    build_medium()
    build_long()
