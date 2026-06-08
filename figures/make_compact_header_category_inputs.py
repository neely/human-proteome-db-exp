#!/usr/bin/env python3
"""
Create compact Figure 2 input files from observed peptide uniqueness outputs.

No external Python packages required.

Inputs:
    derived/peptide_identification/results_summary.csv
    derived/peptide_identification/unique_headers/*_unique_headers.csv

Outputs:
    derived/peptide_identification/observed_unique_counts.csv
    derived/peptide_identification/unique_header_category_counts.csv

The header categorization is keyword-based and should be interpreted as a
descriptive summary, not a formal ontology or enrichment analysis.
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

PEPTIDE_ID_DIR = ROOT / "derived" / "peptide_identification"
RESULTS_SUMMARY = PEPTIDE_ID_DIR / "results_summary.csv"
UNIQUE_HEADERS_DIR = PEPTIDE_ID_DIR / "unique_headers"

OBSERVED_OUT = PEPTIDE_ID_DIR / "observed_unique_counts.csv"
CATEGORY_OUT = PEPTIDE_ID_DIR / "unique_header_category_counts.csv"

CATEGORY_ORDER = [
    "immunoglobulin",
    "t_cell_receptor",
    "hla_mhc",
    "keratin",
    "olfactory_receptor",
    "zinc_finger",
    "predicted_uncharacterized",
    "not_found",
    "other",
]


def short_label(name: str) -> str:
    replacements = {
        "Pangenome-HPRC-Human-": "HPRC-",
        "UniProt_Human_": "UniProt-",
        "UniProt-Human-": "UniProt-",
        "NCBI-RefSeq-Human-": "NCBI-RefSeq-",
        "NCBI_RefSeq-Human-": "NCBI-RefSeq-",
        "NCBI-Human-": "NCBI-",
        "Ensembl-Human-": "Ensembl-",
        "_unique_headers.csv": "",
        "_results.csv": "",
    }

    out = name
    for old, new in replacements.items():
        out = out.replace(old, new)
    return out


def categorize_header(header: str) -> str:
    if not header or not header.strip():
        return "not_found"

    text = header.lower()

    if "not found" in text:
        return "not_found"

    if (
        "immunoglobulin" in text
        or any(
            term in text
            for term in [
                "ighv",
                "ighj",
                "ighd",
                "ighg",
                "igha",
                "ighm",
                "igkv",
                "igkj",
                "iglv",
                "iglj",
            ]
        )
        or re.search(r"\big[hkla][vjcagm]?\d*", text)
    ):
        return "immunoglobulin"

    if (
        "t cell receptor" in text
        or "t-cell receptor" in text
        or any(
            term in text
            for term in [
                "trav",
                "traj",
                "trbv",
                "trbj",
                "trgv",
                "trgj",
                "trdv",
                "trdj",
            ]
        )
    ):
        return "t_cell_receptor"

    if (
        "hla-" in text
        or "major histocompatibility" in text
        or "mhc" in text
        or re.search(r"\bhla[A-Za-z0-9_-]*", header)
    ):
        return "hla_mhc"

    if "keratin" in text or re.search(r"\bkrt\d", text):
        return "keratin"

    if "olfactory receptor" in text or re.search(r"\bor\d+[a-z]\d", text):
        return "olfactory_receptor"

    if "zinc finger" in text or re.search(r"\bznf\d", text):
        return "zinc_finger"

    if any(
        term in text
        for term in ["predicted", "uncharacterized", "hypothetical", "novel protein"]
    ):
        return "predicted_uncharacterized"

    return "other"


def make_observed_unique_counts() -> None:
    if not RESULTS_SUMMARY.exists():
        raise FileNotFoundError(f"Missing {RESULTS_SUMMARY}")

    with RESULTS_SUMMARY.open("r", newline="", encoding="utf-8") as f_in, OBSERVED_OUT.open(
        "w", newline="", encoding="utf-8"
    ) as f_out:
        reader = csv.DictReader(f_in)

        required = {"folder", "total_peptides", "unique_to_this_folder", "shared_across_all"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns from {RESULTS_SUMMARY}: {missing}")

        fieldnames = [
            "database",
            "label",
            "total_peptides",
            "unique_to_this_folder",
            "shared_across_all",
        ]

        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            database = row["folder"]
            writer.writerow(
                {
                    "database": database,
                    "label": short_label(database),
                    "total_peptides": row["total_peptides"],
                    "unique_to_this_folder": row["unique_to_this_folder"],
                    "shared_across_all": row["shared_across_all"],
                }
            )

    print(f"Wrote: {OBSERVED_OUT}")


def database_from_header_file(path: Path) -> str:
    return path.name.replace("_unique_headers.csv", "")


def make_unique_header_category_counts() -> None:
    files = sorted(UNIQUE_HEADERS_DIR.glob("*_unique_headers.csv"))

    if not files:
        raise FileNotFoundError(f"No *_unique_headers.csv files found in {UNIQUE_HEADERS_DIR}")

    priority = {cat: i for i, cat in enumerate(CATEGORY_ORDER)}

    # Store one final category per database + peptide.
    peptide_category: dict[tuple[str, str], str] = {}

    for path in files:
        database = database_from_header_file(path)

        with path.open("r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)

            required = {"peptide", "header"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"Missing columns from {path}: {missing}")

            for row in reader:
                peptide = row["peptide"].strip()
                header = row["header"].strip()
                category = categorize_header(header)

                key = (database, peptide)

                if key not in peptide_category:
                    peptide_category[key] = category
                else:
                    old = peptide_category[key]
                    if priority[category] < priority[old]:
                        peptide_category[key] = category

    counts: dict[tuple[str, str], int] = defaultdict(int)
    databases = sorted({database for database, _peptide in peptide_category.keys()})

    for (database, _peptide), category in peptide_category.items():
        counts[(database, category)] += 1

    with CATEGORY_OUT.open("w", newline="", encoding="utf-8") as f:
        fieldnames = ["database", "label", "category", "peptide_count"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for database in databases:
            for category in CATEGORY_ORDER:
                writer.writerow(
                    {
                        "database": database,
                        "label": short_label(database),
                        "category": category,
                        "peptide_count": counts[(database, category)],
                    }
                )

    print(f"Wrote: {CATEGORY_OUT}")

    print("Category count totals:")
    totals = defaultdict(int)
    for (_database, category), count in counts.items():
        totals[category] += count

    for category in CATEGORY_ORDER:
        print(f"  {category}: {totals[category]:,}")


def main() -> None:
    make_observed_unique_counts()
    make_unique_header_category_counts()


if __name__ == "__main__":
    main()