"""
fasta_tryptic_compare.py

In silico tryptic digest of all FASTA files in a directory.
Outputs:
  - peptide_matrix.csv  : peptides x FASTAs presence/absence matrix
  - summary.csv         : per-FASTA counts (total peptides, unique to that FASTA, shared across all)
  - unique_peptides/    : per-FASTA file listing peptides unique to that FASTA with source protein headers

Usage:
    python fasta_tryptic_compare.py --fasta_dir ./fastas --missed 2 --min_len 7 --max_len 50
"""

import os
import re
import csv
import argparse
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta_dir", default=".", help="Directory containing FASTA files")
    p.add_argument("--missed", type=int, default=2, help="Max missed cleavages")
    p.add_argument("--min_len", type=int, default=7, help="Min peptide length (aa)")
    p.add_argument("--max_len", type=int, default=50, help="Max peptide length (aa)")
    p.add_argument("--out_dir", default=".", help="Output directory")
    return p.parse_args()


def read_fasta(path):
    """Yield (header, sequence) tuples."""
    header, seq = None, []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:], []
            else:
                seq.append(line)
    if header is not None:
        yield header, "".join(seq)


def trypsin_digest(seq, missed, min_len, max_len):
    """
    Cleave after K or R, not before P.
    Returns list of peptide strings passing length filter.
    """
    # Split at K/R not followed by P
    sites = [0]
    for m in re.finditer(r"[KR](?!P)", seq):
        sites.append(m.end())
    sites.append(len(seq))

    base = []
    for i in range(len(sites) - 1):
        pep = seq[sites[i]:sites[i+1]]
        if pep:
            base.append(pep)

    peptides = []
    for i in range(len(base)):
        for mc in range(missed + 1):
            if i + mc + 1 > len(base):
                break
            pep = "".join(base[i:i+mc+1])
            if min_len <= len(pep) <= max_len:
                peptides.append(pep)

    return peptides


def main():
    args = parse_args()
    fasta_files = sorted(f for f in os.listdir(args.fasta_dir) if f.endswith(".fasta"))
    if not fasta_files:
        print("No .fasta files found in", args.fasta_dir)
        return

    os.makedirs(args.out_dir, exist_ok=True)
    uniq_dir = os.path.join(args.out_dir, "unique_peptides")
    os.makedirs(uniq_dir, exist_ok=True)

    # peptide -> {fasta_name -> [headers]}
    pep_to_sources = defaultdict(lambda: defaultdict(list))

    print(f"Digesting {len(fasta_files)} FASTA files...")
    for fname in fasta_files:
        path = os.path.join(args.fasta_dir, fname)
        seen = set()
        for header, seq in read_fasta(path):
            for pep in trypsin_digest(seq, args.missed, args.min_len, args.max_len):
                if pep not in seen:
                    pep_to_sources[pep][fname].append(header)
                    seen.add(pep)
        print(f"  {fname}: done")

    all_peptides = sorted(pep_to_sources.keys())
    n_fastas = len(fasta_files)
    print(f"\nTotal unique peptides across all FASTAs: {len(all_peptides)}")

    # --- peptide_matrix.csv ---
    matrix_path = os.path.join(args.out_dir, "peptide_matrix.csv")
    with open(matrix_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["peptide"] + fasta_files)
        for pep in all_peptides:
            row = [pep] + [1 if fname in pep_to_sources[pep] else 0 for fname in fasta_files]
            w.writerow(row)
    print(f"Written: {matrix_path}")

    # --- summary.csv ---
    shared_all = sum(1 for p in all_peptides if len(pep_to_sources[p]) == n_fastas)
    summary_path = os.path.join(args.out_dir, "summary.csv")
    with open(summary_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["fasta", "total_peptides", "unique_to_this_fasta", "shared_across_all"])
        for fname in fasta_files:
            total = sum(1 for p in all_peptides if fname in pep_to_sources[p])
            unique = sum(1 for p in all_peptides if list(pep_to_sources[p].keys()) == [fname])
            w.writerow([fname, total, unique, shared_all])
    print(f"Written: {summary_path}")

    # --- unique_peptides/ : one file per FASTA ---
    for fname in fasta_files:
        unique_peps = [p for p in all_peptides if list(pep_to_sources[p].keys()) == [fname]]
        out_path = os.path.join(uniq_dir, fname.replace(".fasta", "_unique.csv"))
        with open(out_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["peptide", "source_headers"])
            for pep in unique_peps:
                headers = "; ".join(pep_to_sources[pep][fname])
                w.writerow([pep, headers])
        print(f"  Unique peptides for {fname}: {len(unique_peps)}")

    print("\nDone.")


if __name__ == "__main__":
    main()
