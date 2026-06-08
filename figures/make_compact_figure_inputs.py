#!/usr/bin/env python3
"""
Create compact figure-input files from the full in silico tryptic peptide matrix.

No external Python packages required.

The full peptide_matrix.csv is intentionally not stored in GitHub because it is
too large. Temporarily copy it into:

    derived/in_silico_tryptic_digest/peptide_matrix.csv

Then run:

    python figures/make_compact_figure_inputs.py

Expected input format:
    peptide,database1,database2,...,database21
    PEPTIDESEQ,1,0,...,1

Outputs:
    derived/in_silico_tryptic_digest/pairwise_jaccard_long.csv
    derived/in_silico_tryptic_digest/pairwise_jaccard_wide.csv
    derived/in_silico_tryptic_digest/pairwise_overlap_counts_long.csv
"""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TRYPTIC_DIR = ROOT / "derived" / "in_silico_tryptic_digest"
FULL_MATRIX = TRYPTIC_DIR / "peptide_matrix.csv"


def as_present(value: str) -> bool:
    return str(value).strip() in {"1", "1.0", "TRUE", "True", "true"}


def read_presence_sets(matrix_path: Path) -> dict[str, set[str]]:
    """
    Read peptide presence/absence matrix into:
        {database_name: set(peptide sequences)}
    """
    presence_sets: dict[str, set[str]] = {}

    with matrix_path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)

        if reader.fieldnames is None:
            raise ValueError("Input matrix has no header.")

        if "peptide" not in reader.fieldnames:
            raise ValueError("Expected a column named 'peptide'.")

        databases = [c for c in reader.fieldnames if c != "peptide"]

        if not databases:
            raise ValueError("No database columns found.")

        presence_sets = {db: set() for db in databases}

        for i, row in enumerate(reader, start=1):
            peptide = row.get("peptide", "").strip()

            if not peptide:
                continue

            for db in databases:
                if as_present(row.get(db, "")):
                    presence_sets[db].add(peptide)

            if i % 500000 == 0:
                print(f"  read {i:,} peptide rows")

    return presence_sets


def write_pairwise_outputs(presence_sets: dict[str, set[str]]) -> None:
    databases = list(presence_sets.keys())

    long_path = TRYPTIC_DIR / "pairwise_jaccard_long.csv"
    wide_path = TRYPTIC_DIR / "pairwise_jaccard_wide.csv"
    overlap_path = TRYPTIC_DIR / "pairwise_overlap_counts_long.csv"

    long_rows: list[dict[str, object]] = []

    for db1 in databases:
        set1 = presence_sets[db1]
        total1 = len(set1)

        for db2 in databases:
            set2 = presence_sets[db2]
            total2 = len(set2)

            shared = len(set1 & set2)
            union = len(set1 | set2)
            jaccard = shared / union if union else 0.0

            long_rows.append(
                {
                    "database_1": db1,
                    "database_2": db2,
                    "peptides_database_1": total1,
                    "peptides_database_2": total2,
                    "shared_peptides": shared,
                    "union_peptides": union,
                    "jaccard": f"{jaccard:.8f}",
                    "percent_database_1_in_database_2": f"{(shared / total1) if total1 else 0.0:.8f}",
                    "percent_database_2_in_database_1": f"{(shared / total2) if total2 else 0.0:.8f}",
                }
            )

    with long_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "database_1",
            "database_2",
            "peptides_database_1",
            "peptides_database_2",
            "shared_peptides",
            "union_peptides",
            "jaccard",
            "percent_database_1_in_database_2",
            "percent_database_2_in_database_1",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(long_rows)

    with overlap_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "database_1",
            "database_2",
            "peptides_database_1",
            "peptides_database_2",
            "shared_peptides",
            "union_peptides",
            "jaccard",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for row in long_rows:
            writer.writerow({key: row[key] for key in fieldnames})

    jaccard_lookup = {
        (row["database_1"], row["database_2"]): row["jaccard"]
        for row in long_rows
    }

    with wide_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["database"] + databases)

        for db1 in databases:
            writer.writerow(
                [db1] + [jaccard_lookup[(db1, db2)] for db2 in databases]
            )

    print(f"Wrote: {long_path}")
    print(f"Wrote: {wide_path}")
    print(f"Wrote: {overlap_path}")


def main() -> None:
    if not FULL_MATRIX.exists():
        raise FileNotFoundError(
            f"Missing full matrix: {FULL_MATRIX}\n"
            "Copy peptide_matrix.csv into derived/in_silico_tryptic_digest/ first."
        )

    print(f"Reading full tryptic matrix: {FULL_MATRIX}")
    presence_sets = read_presence_sets(FULL_MATRIX)

    print("Database peptide counts:")
    for db, peptides in presence_sets.items():
        print(f"  {db}: {len(peptides):,}")

    write_pairwise_outputs(presence_sets)


if __name__ == "__main__":
    main()