#!/usr/bin/env python3
"""
Create Figure 1 draft panels.

Uses only Python standard library plus matplotlib.

Inputs:
    derived/in_silico_tryptic_digest/summary.csv
    derived/in_silico_tryptic_digest/pairwise_jaccard_wide.csv

Outputs:
    figures/figure_1_tryptic_overlap/figure_1A_tryptic_counts.png
    figures/figure_1_tryptic_overlap/figure_1B_pairwise_jaccard_heatmap.png
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]

SUMMARY = ROOT / "derived" / "in_silico_tryptic_digest" / "summary.csv"
JACCARD = ROOT / "derived" / "in_silico_tryptic_digest" / "pairwise_jaccard_wide.csv"
OUT_DIR = ROOT / "figures" / "figure_1_tryptic_overlap"


def short_label(name: str) -> str:
    replacements = {
        "Pangenome-HPRC-Human-": "HPRC-",
        "UniProt_Human_": "UniProt-",
        "UniProt-Human-": "UniProt-",
        "NCBI-RefSeq-Human-": "NCBI-RefSeq-",
        "NCBI_RefSeq-Human-": "NCBI-RefSeq-",
        "NCBI-Human-": "NCBI-",
        "Ensembl-Human-": "Ensembl-",
        ".pep.all.fasta": "",
        ".fasta": "",
        "_2023_05": "",
        "-2023_05": "",
        "-2022_07-pep": "",
        "-2022_08-pep": "",
    }

    out = name
    for old, new in replacements.items():
        out = out.replace(old, new)
    return out


def read_summary() -> list[dict[str, str]]:
    with SUMMARY.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    required = {"fasta", "total_peptides", "unique_to_this_fasta"}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise ValueError(f"Missing columns from {SUMMARY}: {missing}")

    return rows


def read_jaccard() -> tuple[list[str], list[list[float]]]:
    with JACCARD.open("r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)

        labels = header[1:]
        matrix = []

        for row in reader:
            matrix.append([float(x) for x in row[1:]])

    return labels, matrix


def make_counts_plot() -> None:
    rows = read_summary()

    rows.sort(key=lambda r: float(r["total_peptides"]))

    labels = [short_label(r["fasta"]) for r in rows]
    total = [float(r["total_peptides"]) for r in rows]
    unique = [float(r["unique_to_this_fasta"]) for r in rows]

    y = list(range(len(rows)))

    fig, ax = plt.subplots(figsize=(9, 8))

    ax.barh(y, total, label="Total theoretical peptides")
    ax.barh(y, unique, label="Unique to this FASTA")

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Number of in silico tryptic peptides")
    ax.set_title("In silico tryptic peptide content across FASTA databases")
    ax.legend(frameon=False)

    fig.tight_layout()

    out = OUT_DIR / "figure_1A_tryptic_counts.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Wrote: {out}")


def make_jaccard_heatmap() -> None:
    labels, matrix = read_jaccard()
    short_labels = [short_label(x) for x in labels]

    fig, ax = plt.subplots(figsize=(9, 8))

    im = ax.imshow(matrix, aspect="auto")

    ax.set_xticks(range(len(short_labels)))
    ax.set_xticklabels(short_labels, rotation=90, fontsize=6)

    ax.set_yticks(range(len(short_labels)))
    ax.set_yticklabels(short_labels, fontsize=6)

    ax.set_title("Pairwise Jaccard similarity of in silico tryptic peptide sets")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Jaccard similarity")

    fig.tight_layout()

    out = OUT_DIR / "figure_1B_pairwise_jaccard_heatmap.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Wrote: {out}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if not SUMMARY.exists():
        raise FileNotFoundError(f"Missing {SUMMARY}")

    if not JACCARD.exists():
        raise FileNotFoundError(f"Missing {JACCARD}")

    make_counts_plot()
    make_jaccard_heatmap()


if __name__ == "__main__":
    main()