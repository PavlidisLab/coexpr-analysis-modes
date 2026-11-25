# Figures and the code that makes them

Published images (`figures/main`, `figures/supplementary`, `figures/extended_data`) were taken
unchanged from the submitted manuscript and supplement. Most were assembled from separate panels;
the panels as written by the original pipeline are in `figures/pipeline_panels/<analysis>/`.
Files ending in `.emf` are vector images exported from Word; open them in Word or PowerPoint
(GitHub cannot preview them).

Code paths: `original_pipeline/` = `code/original_pipeline/` (R),
`snap/` = `code/snap_case_study/`, `supp/` = `code/supp_figures/`.

| Figure | Code | Panels (`figures/pipeline_panels/`) |
|---|---|---|
| Figure 1. Analysis modes (schematic) | drawn | – |
| Figure 2. xCell replicability | `original_pipeline/041_ANALYSIS_01_xcell_coex.R` | `01_xcell_coex/figure_06-08` |
| Figure 3. xSubject replicability; xCell vs xSubject | `original_pipeline/041_ANALYSIS_02_xsbj_coex.R` | `02_xsbj_coex/figure_01-06` |
| Figure 4. xBulk replicability and preservation | `original_pipeline/041_ANALYSIS_03_xbulk_coex.R` | `03_xbulk_coex/figure_01-03` |
| Figure 5. Curated regulatory interactions (DTRI) | `original_pipeline/041_ANALYSIS_07_dtri_comparison.R` | `07_dtri_comparison/figure_01-04` |
| Figure 6A. SNAP preservation odds ratios | `snap/vis/odds-ratio.r` | – |
| Figure 6B. SNAP replication heatmap | `snap/vis/compare-studies.R`; caption statistics: `snap/vis/overlap-summary.R` (mean ± sd Spearman), `snap/vis/compare-studies-pVals.R` (Wilcoxon, Bonferroni) | – |
| Supplementary Figure 1. Dataset properties | `original_pipeline/041_ANALYSIS_01_xcell_coex.R` | `01_xcell_coex/figure_01-05`, `s_figure_01-02` |
| Supplementary Figure 2. xCell approach (schematic) | drawn | – |
| Supplementary Methods schematic, xSubject construction | drawn | – |
| Supplementary Figure 3. PEER analysis of SNAP | external: [xnrxng/PavlabSZProject](https://github.com/xnrxng/PavlabSZProject) `R/14-eric_peer_evaluation.R` (see `code/snap_case_study/README.md`) | in that repo: `results/Ling/eric/plots_paper/01-ling_10factors` (A), `13-all_donors_expression` (B), `07-lim_5factors` (C) |
| Supplementary Figure 4. SNAP genes in the Ling data | `snap/vis/vis-PEER4.r` (`make_and_plot_top_sets`) | – |
| Supplementary Figure 5 (new in R1) | `supp/make_supp_items.py` | – |
| Supplementary Figure 6 (new in R1) | `supp/make_supp_fig6.py` | – |
| Supplementary Tables 2-5 (new in R1) | `supp/make_supp_items.py` | – |
| Extended Data Figure 1. CCV correction | `original_pipeline/041_ANALYSIS_03_xbulk_coex.R` | `03_xbulk_coex/figure_04-08` |
| Extended Data Figure 2. Cell-type variation vs MGP estimates | `original_pipeline/041_ANALYSIS_03_xbulk_coex.R` (likely) | `03_xbulk_coex/` |
| Extended Data Figure 3. Clusters and CCV | `original_pipeline/041_ANALYSIS_04_xcell_clusters.R` | `04_xcell_clusters/` |
| Extended Data Figure 4. Dilution | `original_pipeline/041_ANALYSIS_06_dilution.R` | `06_dilution/figure_01-08` |
| Extended Data Figure 5. RPG synchrony | `original_pipeline/041_ANALYSIS_06_dilution.R` | `06_dilution/sfigure_03.png` (identical) |
| Extended Data Figure 6. Synchrony and preservation | `original_pipeline/041_ANALYSIS_06_dilution.R` | `06_dilution/sfigure_01.png` (identical) |
| Extended Data Figure 7. Synchrony and preservation | `original_pipeline/041_ANALYSIS_06_dilution.R` | `06_dilution/sfigure_02.png` (identical) |

"Identical": the published image matches the panel pixel for pixel after resizing. The other
assignments are by panel titles and captions.
