"""
results_peptide_overlap.py

Reads results.csv from each subdirectory, extracts identified peptide IDs
from the 'col_names' column, and computes overlap across all searches.

Outputs:
  - results_summary.csv        : per-folder counts (total, unique, shared across all)
  - results_peptide_matrix.csv : peptides x folders presence/absence matrix
  - results_unique/            : per-folder CSV of peptides unique to that search

Usage (run from the directory containing the search result folders):
    python results_peptide_overlap.py
    python results_peptide_overlap.py --results_file results.csv --out_dir ./overlap_out
"""

import os
import csv
import re
import argparse
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--base_dir", default=".", help="Directory containing search result folders")
    p.add_argument("--results_file", default="results.csv", help="Name of results file in each folder")
    p.add_argument("--peptide_col", default="col_names", help="Column name containing peptide IDs")
    p.add_argument("--out_dir", default=".", help="Output directory")
    return p.parse_args()


def strip_modifications(peptide):
    """Remove SAGE mass-notation modifications e.g. ..57.0215. or ..15.994."""
    return re.sub(r'\.\.\d+\.\d+\.', '', peptide)


def read_peptides(path, col):
    peptides = set()
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if col not in reader.fieldnames:
            raise ValueError(f"Column '{col}' not found in {path}. "
                             f"Available columns: {reader.fieldnames}")
        for row in reader:
            val = strip_modifications(row[col].strip())
            if val:
                peptides.add(val)
    return peptides


def main():
    args = parse_args()

    # Find all subdirectories containing the results file
    folders = {}
    for name in sorted(os.listdir(args.base_dir)):
        fpath = os.path.join(args.base_dir, name, args.results_file)
        if os.path.isdir(os.path.join(args.base_dir, name)) and os.path.isfile(fpath):
            folders[name] = fpath

    if not folders:
        print(f"No subdirectories containing '{args.results_file}' found in {args.base_dir}")
        return

    print(f"Found {len(folders)} result folders:")
    for name in folders:
        print(f"  {name}")

    os.makedirs(args.out_dir, exist_ok=True)
    uniq_dir = os.path.join(args.out_dir, "results_unique")
    os.makedirs(uniq_dir, exist_ok=True)

    # Load peptides per folder
    pep_to_folders = defaultdict(set)
    folder_peptides = {}

    for name, fpath in folders.items():
        peps = read_peptides(fpath, args.peptide_col)
        folder_peptides[name] = peps
        for p in peps:
            pep_to_folders[p].add(name)
        print(f"  {name}: {len(peps)} peptides")

    all_peptides = sorted(pep_to_folders.keys())
    folder_names = list(folders.keys())
    n_folders = len(folder_names)
    shared_all = sum(1 for p in all_peptides if len(pep_to_folders[p]) == n_folders)

    print(f"\nTotal unique peptides across all folders: {len(all_peptides)}")
    print(f"Peptides shared across all {n_folders} folders: {shared_all}")

    # --- results_peptide_matrix.csv ---
    matrix_path = os.path.join(args.out_dir, "results_peptide_matrix.csv")
    with open(matrix_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["peptide"] + folder_names)
        for pep in all_peptides:
            row = [pep] + [1 if fn in pep_to_folders[pep] else 0 for fn in folder_names]
            w.writerow(row)
    print(f"Written: {matrix_path}")

    # --- results_summary.csv ---
    summary_path = os.path.join(args.out_dir, "results_summary.csv")
    with open(summary_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["folder", "total_peptides", "unique_to_this_folder", "shared_across_all"])
        for name in folder_names:
            total = len(folder_peptides[name])
            unique = sum(1 for p in folder_peptides[name]
                         if pep_to_folders[p] == {name})
            w.writerow([name, total, unique, shared_all])
    print(f"Written: {summary_path}")

    # --- results_unique/ : unique peptides per folder ---
    for name in folder_names:
        unique_peps = sorted(p for p in folder_peptides[name]
                             if pep_to_folders[p] == {name})
        out_path = os.path.join(uniq_dir, name + "_unique.csv")
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["peptide"])
            for pep in unique_peps:
                w.writerow([pep])
        print(f"  Unique identified peptides for {name}: {len(unique_peps)}")

    print("\nDone.")


if __name__ == "__main__":
    main()
