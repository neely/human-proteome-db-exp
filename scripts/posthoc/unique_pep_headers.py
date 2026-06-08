"""
unique_peptide_headers.py

For each database, finds peptides unique to that search result (from
results_peptide_matrix.csv), locates them in the corresponding FASTA,
and returns the protein header for each.

Outputs:
  - unique_headers/  : one CSV per database with columns [peptide, header]

Usage (run from the directory containing the FASTA files and the matrix):
    python unique_peptide_headers.py
    python unique_peptide_headers.py --matrix results_peptide_matrix.csv --fasta_dir . --out_dir ./unique_headers

Notes:
  - Folder names in the matrix must match FASTA filenames (minus .fasta extension).
    If they don't match exactly, edit the FOLDER_TO_FASTA dict at the top of the script.
  - FASTA sequences are collapsed in memory per file; large files (~200k seqs) may
    take 30-60 seconds each.
"""

import os
import csv
import argparse
from collections import defaultdict


# If folder names don't exactly match FASTA basenames, map them here.
# e.g. "Ensembl-Human-GRCh37": "Ensembl-Human-GRCh37.pep.all.fasta"
# Leave empty to use auto-matching (strips .fasta and compares).
FOLDER_TO_FASTA = {
    "Ensembl-Human-GRCh37": "Ensembl-Human-GRCh37.pep.all.fasta",
    "Ensembl-Human-GRCh38": "Ensembl-Human-GRCh38.pep.all.fasta",
    "Ensembl-Human-T2T": "Ensembl-Human-T2T_GCA_009914755.4-2022_07-pep.fasta",
    "NCBI-Human-Celera2005": "NCBI-Human-ASM2283312v2-Celera2005.fasta",
    "NCBI-Human-Hopkins2022": "NCBI-Human-Hopkins2022-06-Ash1_v2.2.fasta",
    "NCBI-RefSeq-Human-GRCh38.p14": "NCBI-RefSeq-Human-GRCh38.p14-GCF_000001405.40-2022-02-03.fasta",
    "NCBI_RefSeq-Human-T2T": "NCBI_RefSeq-Human-T2T-CHM13v2.0-T2T-GCF_009914755.1-2022-01-24.fasta",
    "Pangenome-HPRC-Human-Barbados": "Pangenome-HPRC-Human-GCA_018466835.1-Barbados-2022_07-pep.fasta",
    "Pangenome-HPRC-Human-Caucasian": "Pangenome-HPRC-Human-GCA_018852605.1-Caucasian-2022_07-pep.fasta",
    "Pangenome-HPRC-Human-Columbia": "Pangenome-HPRC-Human-GCA_018469405.1-Columbia-2022_07-pep.fasta",
    "Pangenome-HPRC-Human-HanChineseSouth": "Pangenome-HPRC-Human-GCA_018472565.1-HanChineseSouth-2022_08-pep.fasta",
    "Pangenome-HPRC-Human-Nigeria": "Pangenome-HPRC-Human-GCA_018469415.1-Nigeria-2022_07-pep.fasta",
    "Pangenome-HPRC-Human-Pakistan": "Pangenome-HPRC-Human-GCA_018505835.1-Pakistan-2022_07-pep.fasta",
    "Pangenome-HPRC-Human-Vietnam": "Pangenome-HPRC-Human-GCA_018504055.1-Vietnam-2022_07-pep.fasta",
    "Pangenome-HPRC-Human-swUSA": "Pangenome-HPRC-Human-GCA_018504625.1-swUSA-2022_07-pep.fasta",
    "UniProt-Human-UP000005640_canonical": "UniProt-Human-UP000005640_canonical-2023_05.fasta",
    "UniProt_Human_9606_SP": "UniProt_Human_9606_SP_2023_05.fasta",
    "UniProt_Human_9606_SP-Tr": "UniProt_Human_9606_SP-Tr-2023_05.fasta",
    "UniProt_Human_9606_SP-iso": "UniProt_Human_9606_SP-iso_2023_05.fasta",
    "UniProt_Human_UP000005640_All": "UniProt_Human_UP000005640_All_2023_05.fasta",
    "nextProt": "nextprot_Human-all-peff-2023-09-11.fasta",
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--unique_dir", default="results_unique",
                   help="Directory containing per-database unique peptide CSVs")
    p.add_argument("--fasta_dir", default=".",
                   help="Directory containing FASTA files")
    p.add_argument("--out_dir", default="unique_headers",
                   help="Output directory")
    return p.parse_args()


def read_unique_dir(unique_dir):
    """
    Reads all *_unique.csv files from results_unique/.
    Returns dict {folder_name: [peptide, ...]}
    """
    unique_peps = {}
    for fname in sorted(os.listdir(unique_dir)):
        if not fname.endswith("_unique.csv"):
            continue
        folder = fname.replace("_unique.csv", "")
        peps = []
        fpath = os.path.join(unique_dir, fname)
        with open(fpath, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                val = row["peptide"].strip()
                if val:
                    peps.append(val)
        unique_peps[folder] = peps
    return unique_peps


def read_fasta_collapsed(path):
    """
    Yields (header, sequence) with multi-line sequences collapsed.
    """
    header, seq_parts = None, []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header, seq_parts = line[1:], []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def find_peptides_in_fasta(fasta_path, peptides):
    """
    Returns dict {peptide: [header, ...]} for all peptides found in fasta.
    A peptide may appear in multiple entries.
    """
    pep_set = set(peptides)
    found = defaultdict(list)
    for header, seq in read_fasta_collapsed(fasta_path):
        for pep in pep_set:
            if pep in seq:
                found[pep].append(header)
    return found


def main():
    args = parse_args()

    unique_peps = read_unique_dir(args.unique_dir)
    folder_names = list(unique_peps.keys())
    print(f"Found {len(folder_names)} databases in {args.unique_dir}:")
    for fn in folder_names:
        print(f"  {fn}: {len(unique_peps[fn])} unique peptides")

    os.makedirs(args.out_dir, exist_ok=True)

    for folder in folder_names:
        peps = unique_peps[folder]
        if not peps:
            print(f"\n{folder}: no unique peptides, skipping")
            continue

        # Resolve FASTA path
        if folder in FOLDER_TO_FASTA:
            fasta_name = FOLDER_TO_FASTA[folder]
        else:
            # Auto-match: find a .fasta file whose basename starts with folder name
            candidates = [f for f in os.listdir(args.fasta_dir)
                          if f.endswith(".fasta") and f.startswith(folder)]
            if not candidates:
                print(f"\n{folder}: WARNING — no matching FASTA found, skipping")
                continue
            fasta_name = candidates[0]

        fasta_path = os.path.join(args.fasta_dir, fasta_name)
        if not os.path.isfile(fasta_path):
            print(f"\n{folder}: WARNING — FASTA not found at {fasta_path}, skipping")
            continue

        print(f"\n{folder}: searching {fasta_name} for {len(peps)} peptides...")
        found = find_peptides_in_fasta(fasta_path, peps)

        out_path = os.path.join(args.out_dir, folder + "_unique_headers.csv")
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["peptide", "header", "note"])
            for pep in sorted(peps):
                headers = found.get(pep, [])
                if headers:
                    for h in headers:
                        w.writerow([pep, h, ""])
                else:
                    w.writerow([pep, "NOT FOUND", "peptide not located in FASTA — check sequence"])
        print(f"  Written: {out_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
