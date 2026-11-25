A4 (R1-M4) on Eric's xSubject definition. results.tsv/summary.tsv: replicability (xSubject variant vs itself, ordered dataset pairs, 216) and preservation (Eric's xCell edges as truth vs the variant, 39 groups), thresholds 0.5/1/5/10%.
Variants (lib/eric_xsbj.py; each validated: 'std' reproduces Eric's saved sbj_coexmats, top-1% overlap 1.0 in every group):
- s:std               = Eric's pipeline: mean over cells of CPM -> log2(x+1) -> quantile-normalise over subjects -> Pearson.
- s:A4_sum_counts_cpm = SUM of raw counts per subject -> CPM over the group's genes -> log2(x+1) -> quantile normalisation.   (the reviewer's "pseudobulk by summation")
- s:A4_sum_cellcpm    = SUM of per-cell CPM -> log2(x+1) -> quantile normalisation.
- s:A4_log_noqn / s:A4_nolog_noqn / s:A4_mean_nolog_qn = diagnostics for the log2 and quantile-normalisation steps that the Methods text does not mention.
Results at 1% (mean OR): see summary.tsv. Short version: summing counts leaves xSubject replicability (6.8 vs 6.5) and preservation (11.0 vs 11.1) essentially unchanged; summing cell CPM is worse (4.1; 3 groups below 1 in preservation);
dropping the log2 transform raises xSubject replicability to 14.4 (preservation 13.1), i.e. a good part of xSubject's deficit is due to log2(x+1) on sparse pseudobulk, not to mean vs sum.
