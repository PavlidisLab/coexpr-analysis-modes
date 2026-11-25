# Figures and the code that makes them

Published images (`figures/main`, `figures/supplemental`) were taken
unchanged from the submitted manuscript and supplement. Most were assembled from separate panels;
the panels as written by the original pipeline are in `figures/pipeline_panels/<analysis>/`.
Files ending in `.emf` are vector images exported from Word (open them in Word or PowerPoint);
each has a PNG copy next to it.

Resolution. Most images are the versions embedded in the submitted Word files (about 160-260 dpi
at print size). Higher-resolution versions replace them where a source exists: Figure 6 was
rendered at 600 dpi from the vector graphics in the submitted PDF; Supplemental Figs S13 and
S6 are the images embedded in the submitted PDF (about 380 and 640 dpi);
Supplemental Figs S8-S10 are the original pipeline PNGs; Supplemental Figs S3 and S4 are
600 dpi PNGs with vector PDFs. The panels in `figures/pipeline_panels/` (mostly EPS, vector) can
be used to rebuild the other figures at any resolution.

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
| Supplemental Fig S1. xCell approach (schematic) | drawn | – |
| Supplemental Fig S2. Dataset properties | `original_pipeline/041_ANALYSIS_01_xcell_coex.R` | `01_xcell_coex/figure_01-05`, `s_figure_01-02` |
| Supplemental Fig S3 (new in R1) | `supp/make_supp_items.py` | – |
| Supplemental Fig S4 (new in R1) | `supp/make_supp_figS4.py` | – |
| Supplemental Fig S5. CCV correction | `original_pipeline/041_ANALYSIS_03_xbulk_coex.R` | `03_xbulk_coex/figure_04-08` |
| Supplemental Fig S6. Cell-type variation vs MGP estimates | `original_pipeline/041_ANALYSIS_03_xbulk_coex.R` (likely) | `03_xbulk_coex/` |
| Supplemental Fig S7. Clusters and CCV | `original_pipeline/041_ANALYSIS_04_xcell_clusters.R` | `04_xcell_clusters/` |
| Supplemental Fig S8. RPG synchrony | `original_pipeline/041_ANALYSIS_06_dilution.R` | `06_dilution/sfigure_03.png` (identical) |
| Supplemental Fig S9. Synchrony and preservation | `original_pipeline/041_ANALYSIS_06_dilution.R` | `06_dilution/sfigure_01.png` (identical) |
| Supplemental Fig S10. Synchrony and preservation | `original_pipeline/041_ANALYSIS_06_dilution.R` | `06_dilution/sfigure_02.png` (identical) |
| Supplemental Fig S11. Dilution | `original_pipeline/041_ANALYSIS_06_dilution.R` | `06_dilution/figure_01-08` |
| Supplemental Fig S12. PEER analysis of SNAP | external: [xnrxng/PavlabSZProject](https://github.com/xnrxng/PavlabSZProject) `R/14-eric_peer_evaluation.R` (see `code/snap_case_study/README.md`) | in that repo: `results/Ling/eric/plots_paper/01-ling_10factors` (A), `13-all_donors_expression` (B), `07-lim_5factors` (C) |
| Supplemental Fig S13. SNAP genes in the Ling data | `snap/vis/vis-PEER4.r` (`make_and_plot_top_sets`) | – |
| Supplemental Methods schematic, xSubject construction | drawn | – |
| Supplemental Tables S2-S5 (new in R1) | `supp/make_supp_items.py` | – |

"Identical": the published image matches the panel pixel for pixel after resizing. The other
assignments are by panel titles and captions.
