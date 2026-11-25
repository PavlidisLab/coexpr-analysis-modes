# SNAP case study (Ling et al. 2024)

Figure 6 and Supplemental Fig S13.

Inputs: the PEER latent factor 4 (SNAP) gene loadings for the Ling data, and the per-dataset
xCell, xCell-xSubject and xSubject coexpression of the seven brain datasets and the Ling data.

| Step | Script | Output |
|---|---|---|
| 1 | `PEER4-genes.py` | `<study>_PEER4_gene_correlations_with_background_<mode>.csv`: correlations among the top 1,000 SNAP genes, ribosomal protein genes and a random background set |
| 2 | `PEER4_top1.py`, `PEER4-top1_overlap.py` | top-1% edges among the same genes, and their overlap across modes |
| 3 | `vis/vis-PEER4.r` | loads step 1 output (`combined1`); Supplemental Fig S13 (`make_and_plot_top_sets`) |
| 4 | `vis/odds-ratio.r` | Figure 6A |
| 5 | `vis/compare-studies.R` | Figure 6B |
| 6 | `vis/overlap-summary.R`, `vis/compare-studies-pVals.R` | statistics in the Figure 6 caption |

Steps 4-6 run in the same R session as step 3 and use the objects it creates. Figures were
exported from the R graphics device.

The PEER analysis itself (Supplemental Fig S12) was done by Rui Xiang Yu in the repository
[xnrxng/PavlabSZProject](https://github.com/xnrxng/PavlabSZProject), which also produced the
LF4 gene-loadings file used in step 1:

| Step | File in that repository | What it does |
|---|---|---|
| 1 | `R/13-eric_peer.R` | per-dataset, per-cell-type inputs for PEER (protein-coding genes) |
| 2 | `bin/peeer_eric.sh` | runs `peertool` with 5, 10 and 15 factors per dataset |
| 3 | `R/14-eric_peer_evaluation.R` | compares the factors with Ling's published LF4; draws and combines the panels in `results/Ling/eric/plots_paper/` |

Supplemental Fig S12 panels: A = `01-ling_10factors`, B = `13-all_donors_expression`,
C = `07-lim_5factors`. Axis labels were edited for the figure.
