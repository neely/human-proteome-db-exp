#!/usr/bin/env python3
"""
Create Figure 2 draft panels.

Uses only Python standard library plus matplotlib.

Inputs:
    derived/peptide_identification/observed_unique_counts.csv
    derived/peptide_identification/unique_header_category_counts.csv

Outputs:
    figures/figure_2_unique_peptide_composition/figure_2A_observed_unique_counts.png
    figures/figure_2_unique_peptide_composition/figure_2B_unique_header_categories.png
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]

PEPTIDE_ID_DIR = ROOT / "derived" / "peptide_identification"

OBSERVED_COUNTS = PEPTIDE_ID_DIR / "observed_unique_counts.csv"
CATEGORY_COUNTS = PEPTIDE_ID_DIR / "unique_header_category_counts.csv"

OUT_DIR = ROOT / "figures" / "figure_2_unique_peptide_composition"

CATEGORY_ORDER = [
    "immunoglobulin",
    "hla_mhc",
    "zinc_finger",
    "predicted_uncharacterized",
    "t_cell_receptor",
    "keratin",
    "olfactory_receptor",
    "not_found",
    "other",
]

CATEGORY_LABELS = {
    "immunoglobulin": "Immunoglobulin",
    "t_cell_receptor": "T-cell receptor",
    "hla_mhc": "HLA / MHC",
    "keratin": "Keratin",
    "olfactory_receptor": "Olfactory receptor",
    "zinc_finger": "Zinc finger",
    "predicted_uncharacterized": "Predicted / uncharacterized",
    "not_found": "Not found",
    "other": "Other",
}


def read_csv_dicts(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return list(reader)


def make_observed_unique_counts_plot() -> None:
    rows = read_csv_dicts(OBSERVED_COUNTS)

    required = {"database", "label", "unique_to_this_folder"}
    missing = required - set(rows[0].keys()) if rows else required
    if missing:
        raise ValueError(f"Missing columns from {OBSERVED_COUNTS}: {missing}")

    rows.sort(key=lambda r: float(r["unique_to_this_folder"]))

    labels = [r["label"] for r in rows]
    counts = [float(r["unique_to_this_folder"]) for r in rows]
    y = list(range(len(rows)))

    fig, ax = plt.subplots(figsize=(9, 8))

    ax.barh(y, counts)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Observed peptides unique to one database search")
    ax.set_title("Database-specific observed peptide identifications")

    fig.tight_layout()

    out = OUT_DIR / "figure_2A_observed_unique_counts.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)

    print(f"Wrote: {out}")


def make_header_category_plot() -> None:
    rows = read_csv_dicts(CATEGORY_COUNTS)

    required = {"database", "label", "category", "peptide_count"}
    missing = required - set(rows[0].keys()) if rows else required
    if missing:
        raise ValueError(f"Missing columns from {CATEGORY_COUNTS}: {missing}")

    labels_by_db = {}
    counts_by_db = {}

    for row in rows:
        db = row["database"]
        label = row["label"]
        category = row["category"]
        count = float(row["peptide_count"])

        labels_by_db[db] = label
        counts_by_db.setdefault(db, {cat: 0.0 for cat in CATEGORY_ORDER})
        counts_by_db[db][category] = count

    databases = sorted(
        counts_by_db,
        key=lambda db: sum(counts_by_db[db].values())
    )

    labels = [labels_by_db[db] for db in databases]
    y = list(range(len(databases)))

    fig, ax = plt.subplots(figsize=(10, 8))

    left = [0.0 for _ in databases]

    for category in CATEGORY_ORDER:
        values = [counts_by_db[db].get(category, 0.0) for db in databases]

        # Skip categories that are all zero to keep the legend cleaner.
        if sum(values) == 0:
            continue

        ax.barh(
            y,
            values,
            left=left,
            label=CATEGORY_LABELS.get(category, category),
        )

        left = [old + val for old, val in zip(left, values)]

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Database-specific peptides")
    ax.set_title("Header-derived composition of database-specific peptides")
    ax.legend(frameon=False, fontsize=7, loc="lower right")

    fig.tight_layout()

    out = OUT_DIR / "figure_2B_unique_header_categories.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)

    print(f"Wrote: {out}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if not OBSERVED_COUNTS.exists():
        raise FileNotFoundError(f"Missing {OBSERVED_COUNTS}")

    if not CATEGORY_COUNTS.exists():
        raise FileNotFoundError(f"Missing {CATEGORY_COUNTS}")

    make_observed_unique_counts_plot()
    make_header_category_plot()


if __name__ == "__main__":
    main()