# Bulk gene counts / cross-dataset expression correlation (Methods claim check)
Claim: ROSMAP 21,518 and Velmeshev 25,326 genes; 21,239 shared; "mean Pearson r = 0.9" (Suppl. Fig 1D).
Files: `gene_counts.tsv`, `results.tsv` (stage, reading, value, n genes), `sc_dataset_pairs.tsv` (7 sc datasets), `results_meta.json`.
Saved baseline objects `bk_exprmats_{rosmap,velmeshev}.rds` are log2(CPM+1) then quantile-normalised ("qn"); recomputing the stages from the source files reproduces them to 2e-15.
Note: the ROSMAP source file is FPKM (ROSMAP_RNAseq_FPKM_gene.tsv), renormalised to per-sample 1e6 in the pipeline ("cpm" stage); Velmeshev is counts.
