# human-proteome-db-exp
This repository contains the manuscript, workflow documentation, analysis scripts,
and derived tables required to reproduce the manuscript figures.

Large input and intermediate files are archived separately on OSF, including:
FASTA databases, full SAGE outputs (`results.sage.tsv`, `lfq.tsv`, `results.json`),
and optional protein-level intermediate outputs.

The repository supports figure-level reproducibility directly. Full pipeline
reproducibility requires downloading the OSF archive.

```
derived/
├── peptide_identification/
│   ├── results_summary.csv
│   ├── results_peptide_matrix.csv
│   ├── results_unique/
│   └── unique_headers/
├── differential_abundance/
│   ├── table1_peptide_dap_summary.csv
│   └── per_database_results_csv/
├── tryptic_digest/
│   ├── summary.csv
│   ├── peptide_matrix.csv
│   └── unique_peptides/
└── sequence_counts/
    └── output.csv
```
