# Codex Project Notes

## Overview
- Project: `bch709_vibe_coding`
- Purpose (from current files): analyze yeast expression data and generate a top-10 CV gene table plus heatmaps.

## Folder Snapshot
- `README.md`: minimal title file.
- `environment.yml`: Conda env with `r-base` and `r-ggplot2`.
- `data/`: input/reference files (about 15 MB total).
- `scripts/`: analysis script(s).
- `results/`: generated tables and figures.

## Key Inputs
- `data/gasch2000.txt` (6.1 MB, 6152 lines): main expression matrix used by the R script.
- `data/chrom.sizes` (229 B)
- `data/mrna.fa.gz` (112 KB)
- `data/sacCer3.fa.gz` (3.7 MB)
- `data/saccharomyces_cerevisiae.gff.gz` (5.0 MB)

## Main Script
- `scripts/top10_cv_heatmap.R`
- Reads `data/gasch2000.txt`, takes condition columns (default all; optional first N via CLI arg), computes per-gene CV, writes top-10 table, and saves a heatmap.
- Outputs:
  - default run:
    - `results/top10_cv_table.tsv`
    - `results/gene_top10_heatmap.png`
  - first N conditions run:
    - `results/top10_cv_table_first<N>_conditions.tsv`
    - `results/gene_top10_heatmap_first<N>_conditions.png`

## Current Results Present
- `results/top10_cv_table.tsv` (11 lines including header)
- `results/top10_cv_table_first30_conditions.tsv` (11 lines including header)
- `results/gene_top10_heatmap.png` (847 KB)
- `results/gene_top10_heatmap_first30_conditions.png` (331 KB)

## How To Run
```bash
Rscript scripts/top10_cv_heatmap.R
Rscript scripts/top10_cv_heatmap.R 30
```
