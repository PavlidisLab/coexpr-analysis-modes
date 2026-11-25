# coexpr-analysis-modes

Code, figures and result tables for the manuscript
**"Evaluation of analysis modes for RNA coexpression in single-cell and bulk tissue"**
(Genome Research, GENOME/2026/282552, revised version R1).

The paper compares three ways of computing gene coexpression from the same data:

| Level | Correlation is computed across | Expression input |
|---|---|---|
| **xCell** | cells, within each subject, then averaged over subjects | per-cell CPM (no log) |
| **xSubject** | subjects (pseudobulk: mean of cell CPM per subject) | log2(x + 1), quantile-normalised across subjects |
| **xBulk** | bulk tissue samples | CPM, log2(x + 1), quantile-normalised |

At every level, coexpression values are rank-normalised (rank / max rank) and the top 1% of gene
pairs are taken as edges. Agreement is the Fisher odds ratio between edge sets: across datasets
(replicability) or across levels of the same dataset (preservation).

## Repository layout

```
figures/
  main/                 Figures 1-6 as published
  supplementary/        Supplementary Figures 1-6 (5 and 6 are new in R1, with PDF versions)
  extended_data/        Extended Data Figures 1-7
  pipeline_panels/      individual panels as written by the original pipeline
code/
  original_pipeline/    the original R pipeline: networks, results, Figures 2-5, Supplementary Figure 1, Extended Data
  snap_case_study/      SNAP case study (Ling et al. 2024): Figure 6, Supplementary Figure 4
  supp_figures/         scripts that build Supplementary Figures 5-6 and Supplementary Tables 2-5
  revision_analyses/    the R1 sensitivity analyses (A1-A10) and the baseline reproduction
results/
  A1_* ... A10_*        summary tables of each sensitivity analysis (each folder has a README)
  eric_baseline/        reproduction of the submitted results from the original coexpression matrices
  eric_saved_tables/    the original pipeline's saved Fisher tables, exported to TSV
  supp_items/           numbers.json, Supplementary Tables 2-5 (TSV), manifest.json with input sha256
```

See [`FIGURES.md`](FIGURES.md) for which file and which script belongs to each figure.

## Regenerating the new supplementary figures and tables

These need only the tables in `results/` (Python 3 with matplotlib; no other data):

```bash
bash code/supp_figures/make_figures.sh
```

This rewrites `figures/supplementary/Supplementary_Figure_{5,6}.{pdf,png}` and
`results/supp_items/` (Supplementary Tables 2-5 and every number cited in the revised text,
in `numbers.json`). Each script checks the qualitative statements made in the figure captions
and stops if the data no longer support them. `results/supp_items/manifest.json` records the
sha256 of every input table.

## Re-running the sensitivity analyses

`code/revision_analyses/run_all.sh` lists the analyses in the order they were run. They start from
the coexpression and expression matrices of the original pipeline (R objects on the lab server),
exported by the scripts in `code/revision_analyses/eric/R/`. Those inputs are not included here
(about 120 GB), and the scripts keep the server paths they were run with, so they document the
analysis rather than run out of the box. Software versions are in
`code/revision_analyses/ENVIRONMENT.txt`.

`code/revision_analyses/lib/eric_engine.py` is a Python re-implementation of the original R
definitions of xCell and xSubject coexpression. It reproduces the published networks
(`results/eric_baseline/validation_vs_saved.tsv`: all 216 xCell and 216 xSubject replicability
odds ratios and all 39 preservation odds ratios match the saved tables).

| Analysis | Reviewer point | Results |
|---|---|---|
| A1 matched sample size (S1 pooled cells, S2 per-subject cap; 50 draws) | R1 M1 | `results/A1_matched_sample` |
| A2 permutation null | R1 M2 | `results/A2_permutation_null` |
| A3 expression strata | R1 M3 | `results/A3_expression_strata` |
| A4 pseudobulk variants | R1 M4 | `results/A4_xsubject_pseudobulk` |
| A5 precision / recall / F1 / Jaccard | R2 | `results/A5_prf_jaccard` |
| A6 edge thresholds; raw and signed correlations | R2 | `results/A6_thresholds`, `results/A6_raw` |
| A7 no imputation | R2 | `results/A7_no_imputation` |
| A8 weighting subjects by cell number | R1 m3 | `results/A8_weighted` |
| A9 odds ratios at or below 1; reference/test asymmetry | R2 | `results/A9_or_le1` |
| A10 common gene set | R1 | `results/A10_common_geneset` |
| Bulk gene counts and cross-dataset correlation | R1 m2 | `results/A_bulk_check` |

## Data

The single-cell datasets (Lau 2020, Lim 2022, Nagy 2020, Pineda 2021, Ramos 2022,
ROSMAP/Mathys 2019, Velmeshev 2019) and the ROSMAP and Velmeshev bulk data are public; see
Supplementary Table 1 of the paper for sources and accessions.

## Citation

If you use this code, please cite the paper (citation to be added on publication).
