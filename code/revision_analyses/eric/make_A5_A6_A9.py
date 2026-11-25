"""A5, A6 (rank-based part), A9 on Eric Chu's pipeline; plus the baseline summary that reproduces the manuscript numbers."""
import pandas as pd, numpy as np, json, hashlib, datetime, os
B=os.path.expanduser('~/xcell-xbulk-revision/'); D=pd.read_csv(B+'eric/baseline/all_comparisons.tsv',sep='\t')
keep=['analysis','cell_type','reference_level','reference_dataset','test_level','test_dataset','frac','universe','n_genes_ref','n_genes_test','a','b','c','d','N','n_ref_pos','n_test_pos','odds_ratio','p_value','log_p','precision','recall','F1','jaccard']
src={'all_comparisons.tsv':hashlib.sha256(open(B+'eric/baseline/all_comparisons.tsv','rb').read()).hexdigest()}
def meta(d,extra): json.dump(dict(script=os.path.abspath(__file__),script_sha256=hashlib.sha256(open(__file__,'rb').read()).hexdigest(),inputs=src,timestamp=datetime.datetime.now().isoformat(),source='Eric Chu 041_CH4_FINAL saved matrices (read-only)',seed='n/a (deterministic)',**extra),open(B+d+'/results_meta.json','w'),indent=1)
main=D[D.analysis.isin(['replicability','preservation'])]
# A5
a5=main[main.frac==0.01][keep]; a5.to_csv(B+'A5_prf_jaccard/results.tsv',sep='\t',index=False)
a5.groupby(['analysis','reference_level','test_level']).agg(n=('odds_ratio','size'),mean_OR=('odds_ratio','mean'),mean_precision=('precision','mean'),mean_recall=('recall','mean'),mean_F1=('F1','mean'),mean_jaccard=('jaccard','mean'),min_jaccard=('jaccard','min'),max_jaccard=('jaccard','max')).reset_index().to_csv(B+'A5_prf_jaccard/summary.tsv',sep='\t',index=False)
meta('A5_prf_jaccard',dict(definition='reference = truth, test = prediction; precision=a/(a+b); recall=a/(a+c); F1=2a/(2a+b+c); Jaccard=a/(a+b+c); universe = all pairs of the test group genes (Eric evaluationUtils$fisher)'))
# A6
main[keep].to_csv(B+'A6_thresholds/results.tsv',sep='\t',index=False)
s6=main.groupby(['analysis','reference_level','test_level','frac']).agg(n=('odds_ratio','size'),mean_OR=('odds_ratio','mean'),median_OR=('odds_ratio','median'),min_OR=('odds_ratio','min'),n_OR_le_1=('odds_ratio',lambda x:int((x<=1).sum())),mean_F1=('F1','mean'),mean_jaccard=('jaccard','mean'),mean_n_test_pos=('n_test_pos','mean')).reset_index(); s6.to_csv(B+'A6_thresholds/summary.tsv',sep='\t',index=False)
meta('A6_thresholds',dict(definition='positive = rank-normalised coexpression >= 1-frac (frac 0.005,0.01,0.05,0.10), lower triangle; raw-correlation version in A6_raw'))
# A9
le=main[(main.odds_ratio<=1)|main.odds_ratio.isna()][keep]; le.to_csv(B+'A9_or_le1/comparisons_OR_le_1.tsv',sep='\t',index=False)
rows=[]
for (an,rl,tl,fr),g in main.groupby(['analysis','reference_level','test_level','frac']):
    i=g.odds_ratio.idxmin(); rows.append(dict(analysis=an,reference_level=rl,test_level=tl,frac=fr,n=len(g),min_OR=g.odds_ratio.min(),argmin=' | '.join(g.loc[i,['cell_type','reference_dataset','test_dataset']].astype(str)),n_OR_le_1=int((g.odds_ratio<=1).sum()),mean_OR=g.odds_ratio.mean()))
pd.DataFrame(rows).to_csv(B+'A9_or_le1/min_OR_summary.tsv',sep='\t',index=False)
r1=main[(main.frac==0.01)&(main.analysis=='replicability')]
comp=r1.groupby(['reference_level','cell_type']).size().reset_index(name='n_ordered_pairs'); comp['datasets']=comp.apply(lambda x:r1[(r1.cell_type==x.cell_type)&(r1.reference_level==x.reference_level)].reference_dataset.nunique(),axis=1); comp['expected']=comp.datasets*(comp.datasets-1); comp.to_csv(B+'A9_or_le1/completeness_replicability.tsv',sep='\t',index=False)
sym=r1.merge(r1,left_on=['reference_level','cell_type','reference_dataset','test_dataset'],right_on=['reference_level','cell_type','test_dataset','reference_dataset'],suffixes=('','_rev'))
ld=np.abs(np.log(sym.odds_ratio)-np.log(sym.odds_ratio_rev)); 
bsum=dict(max_abs_log_OR_asymmetry=float(ld.max()),median_abs_log_OR_asymmetry=float(ld.median()),frac_pairs_asym_gt_1pct=float((ld>0.01).mean()))
json.dump(bsum,open(B+'A9_or_le1/asymmetry.json','w'),indent=1)
# baseline summary reproducing manuscript numbers
sc=r1[r1.reference_level=='sc']; sb=r1[r1.reference_level=='sbj']
j=sc.merge(sb,on=['cell_type','reference_dataset','test_dataset'],suffixes=('_sc','_sbj')); gm=j.groupby(['cell_type','reference_dataset'])[['odds_ratio_sc','odds_ratio_sbj']].mean(); delta=(gm.odds_ratio_sbj-gm.odds_ratio_sc)
pr=main[(main.frac==0.01)&(main.analysis=='preservation')]
nc=pd.read_csv(B+'eric/tables/saved_sc_stats_ncells.tsv',sep='\t')
summ=dict(groups=int(len(pr)),cells_total=int(nc.n_cell.sum()),subjects_rows=int(len(nc)),sc_repro_n=int(len(sc)),sc_repro_mean_OR=float(sc.odds_ratio.mean()),sc_repro_min_OR=float(sc.odds_ratio.min()),sbj_repro_mean_OR=float(sb.odds_ratio.mean()),sbj_repro_min_OR=float(sb.odds_ratio.min()),sbj_n_OR_le_1=int((sb.odds_ratio<=1).sum()),
  mean_group_delta_OR_sbj_minus_sc=float(delta.mean()),n_groups_delta=int(len(delta)),preservation_n=int(len(pr)),preservation_n_OR_gt_1=int((pr.odds_ratio>1).sum()),preservation_mean_OR=float(pr.odds_ratio.mean()),preservation_min_OR=float(pr.odds_ratio.min()),preservation_mean_OR_of_gt1=float(pr[pr.odds_ratio>1].odds_ratio.mean()))
json.dump(summ,open(B+'eric/baseline/manuscript_numbers_check.json','w'),indent=1); print(json.dumps(summ,indent=1)); print(json.dumps(bsum)); print(s6.round(3).to_string()); print(pd.DataFrame(rows).query('frac==0.01').to_string()); print(comp.to_string())
pr.sort_values('odds_ratio')[['cell_type','reference_dataset','odds_ratio','n_genes_ref','n_ref_pos','n_test_pos']].head(4).to_csv(B+'eric/baseline/preservation_lowest.tsv',sep='\t',index=False)
print(pr.sort_values('odds_ratio')[['cell_type','reference_dataset','odds_ratio']].head(4).to_string()); print(sb[sb.odds_ratio<=1][['cell_type','reference_dataset','test_dataset','odds_ratio']].to_string())
