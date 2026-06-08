\# Workflow overview



This repository documents a comparison of human protein sequence databases in bottom-up proteomics.



\## Primary analysis



Eight steady-state immune-cell mzML files were searched independently against 21 human protein sequence FASTA databases using SAGE with a common `params.json` configuration. Each search produced SAGE output files including `results.sage.tsv`, `lfq.tsv`, and `results.json`.



Peptide-level LFQ outputs were filtered and processed for differential abundance analysis between B naive and T4 naive cells. Differentially abundant peptides were identified using an unpaired Wilcoxon rank-sum test followed by Benjamini-Hochberg correction.



\## Post hoc analyses



Two post hoc analyses support the manuscript figures.



1\. Search-result peptide overlap:

&#x20;  - Script: `scripts/posthoc/results\_overlap.py`

&#x20;  - Inputs: per-database `results.csv`

&#x20;  - Outputs: `derived/peptide\_identification/results\_peptide\_matrix.csv`, `results\_summary.csv`, `results\_unique/`



2\. Unique peptide header lookup:

&#x20;  - Script: `scripts/posthoc/unique\_pep\_headers.py`

&#x20;  - Inputs: `results\_unique/` and FASTA files

&#x20;  - Outputs: `derived/peptide\_identification/unique\_headers/`



3\. In silico tryptic FASTA comparison:

&#x20;  - Script: `scripts/posthoc/tryptic\_compare.py`

&#x20;  - Inputs: 21 FASTA files archived on OSF

&#x20;  - Outputs: `derived/tryptic\_digest/summary.csv` and `unique\_peptides/`

&#x20;  - The full tryptic peptide matrix is archived on OSF because it exceeds GitHub-friendly file size.



\## Data availability



This GitHub repository contains scripts, workflow documentation, manuscript materials, and small derived files needed to reproduce figures. Large input and intermediate files, including FASTA files and full SAGE outputs, are archived on OSF.

