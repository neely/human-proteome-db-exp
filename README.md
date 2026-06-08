# Human Proteome Database Choice Experiment

This repository contains manuscript materials, workflow documentation, analysis scripts, and derived figure-input tables for a comparison of human protein sequence databases in bottom-up proteomics.

Bottom-up proteomics experiments depend fundamentally on the protein sequence database used during spectral searching. Every peptide-spectrum match is constrained by what sequences are present in the search space; a peptide whose sequence is absent from the database cannot be identified, regardless of spectral quality. This project asks whether major human protein sequence databases are similar enough that database choice is analytically inconsequential, or whether the choice of FASTA database measurably affects peptide identification and downstream differential abundance results.

The study compares 21 human protein sequence databases spanning UniProt, neXtProt, Ensembl, NCBI RefSeq, alternative human assemblies, T2T-CHM13, and HPRC pangenome-derived resources. The same human immune-cell mass spectrometry files were searched against each database using the same SAGE parameters, and peptide-level differential abundance was evaluated between B naive and T4 naive immune cells.

## Repository scope

This GitHub repository is intended for **figure-level reproducibility** and transparent documentation of the analysis workflow.

It includes:

- manuscript draft materials
- workflow documentation and Mermaid workflow diagram
- SAGE configuration file
- scripts used for primary and post hoc analyses
- small derived tables required to reproduce manuscript figures
- database manifest describing the 21 FASTA databases

Large files are not stored in GitHub. FASTA databases, full SAGE outputs, and bulky intermediate files are archived separately on OSF.

## Data availability

Large input and intermediate files are archived on OSF, including:

- the 21 human protein sequence FASTA files
- full SAGE search outputs, including `results.sage.tsv`, `lfq.tsv`, and `results.json`
- the full in silico tryptic peptide presence/absence matrix
- other large intermediate outputs not suitable for GitHub

GitHub contains compact derived outputs needed to regenerate the manuscript figures. Full pipeline reproduction requires downloading the OSF archive.

OSF link: **to be added**

## Repository structure

```text
.
├── config/
│   ├── params.json
│   └── database_manifest.csv
├── derived/
│   ├── differential_abundance/
│   │   └── per_database_results_csv/
│   ├── peptide_identification/
│   │   ├── results_summary.csv
│   │   ├── results_peptide_matrix.csv
│   │   ├── results_unique/
│   │   └── unique_headers/
│   ├── sequence_counts/
│   │   └── output.csv
│   └── tryptic_digest/
│       ├── summary.csv
│       └── unique_peptides/
├── figures/
│   ├── figure_1_tryptic_overlap/
│   └── figure_2_unique_peptide_composition/
├── manuscript/
├── scripts/
│   ├── posthoc/
│   ├── primary_pipeline/
│   ├── protein_level/
│   └── utilities/
└── workflow/
    ├── workflow_diagram.mermaid
    └── workflow.md
```

## Main analysis workflow

The primary analysis searched eight steady-state immune-cell mzML files against 21 human protein sequence FASTA databases using SAGE and a common `params.json` configuration.

For each database, SAGE generated peptide-spectrum match, LFQ, and metadata outputs. Peptide-level LFQ outputs were processed for differential abundance testing between B naive and T4 naive cells using an unpaired Wilcoxon rank-sum test followed by Benjamini-Hochberg correction.

The workflow is summarized in:

- `workflow/workflow_diagram.mermaid`
- `workflow/workflow.md`

## Protein sequence databases

The FASTA databases represent several major types of human protein sequence resources:

- UniProt reference proteome and taxonomy-based retrievals
- neXtProt human PEFF sequence database
- Ensembl GRCh37, GRCh38, and T2T peptide FASTAs
- NCBI RefSeq GRCh38 and T2T resources
- NCBI alternative or legacy assembly-derived resources
- HPRC pangenome-derived Ensembl rapid-release peptide FASTAs

The database manifest is provided in:

```text
config/database_manifest.csv
```

The manifest records database short name, source, category, assembly or resource, FASTA filename, sequence count, release/download date, search-output folder name, and OSF location.

## Post hoc analyses

Post hoc analyses were used to support manuscript figures and to evaluate database-specific peptide content.

### Peptide identification overlap

Script:

```text
scripts/posthoc/results_overlap.py
```

Inputs:

```text
derived/differential_abundance/per_database_results_csv/
```

Outputs:

```text
derived/peptide_identification/results_summary.csv
derived/peptide_identification/results_peptide_matrix.csv
derived/peptide_identification/results_unique/
```

This analysis compares observed peptide identifications across database searches and identifies peptides unique to each search result.

### Unique peptide header lookup

Script:

```text
scripts/posthoc/unique_pep_headers.py
```

Inputs:

```text
derived/peptide_identification/results_unique/
OSF FASTA archive
```

Outputs:

```text
derived/peptide_identification/unique_headers/
```

This analysis maps database-specific identified peptides back to protein FASTA headers. These outputs are used to examine whether database-specific identifications are associated with particular protein families, annotation classes, or highly variable regions.

### In silico tryptic FASTA comparison

Script:

```text
scripts/posthoc/tryptic_compare.py
```

Inputs:

```text
OSF FASTA archive
```

Outputs included in GitHub:

```text
derived/tryptic_digest/summary.csv
derived/tryptic_digest/unique_peptides/
```

The full tryptic peptide presence/absence matrix is archived on OSF because it exceeds GitHub-friendly file size. The compact summaries retained in GitHub support figure-level reproduction.

## Reproducibility modes

### Figure-level reproduction

Clone this repository and use the files in `derived/` with the figure-generation scripts in `figures/`.

This mode does not require downloading the full FASTA or SAGE output archive.

### Full analysis reproduction

Download the OSF archive, including FASTA files and SAGE outputs, then rerun the scripts in:

```text
scripts/primary_pipeline/
scripts/posthoc/
scripts/utilities/
```

The SAGE search parameters are provided in:

```text
config/params.json
```

## File-size policy

Large input and intermediate files are intentionally excluded from GitHub, including:

- `*.fasta`
- `*.RAW`
- `*.mzML`
- `*.mzML.gz`
- `*.sage.tsv`
- `lfq.tsv`
- `results.json`
- `MSstats-Sage.csv`
- the full tryptic `peptide_matrix.csv`

These files are ignored by `.gitignore` and should be stored in the OSF archive instead.

## Manuscript

The manuscript draft is located in:

```text
manuscript/
```

The central result is that database choice changes both peptide identification counts and peptide-level differential abundance calls, even when the same mass spectrometry files, search engine, search parameters, and statistical thresholds are used throughout.

## Citation

Manuscript citation to be added after publication.

## License

See `LICENSE`.