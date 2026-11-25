"""Replicability / preservation before and after the A2 filter. Filter F_g of a group = gene pairs whose real correlation is beyond the permutation null (positive above the null maximum or negative below the null minimum; eric/a2sig/*).
Filtered evaluation: reference positives = baseline top-1% edges that are in F_ref; test positives = baseline top-1% edges in F_test; universe = F_test (pairs significant in the test group); same 2x2 definition as Eric's fisher otherwise.
Levels: sc (xCell, full sample size), sbj (xSubject). usage: eric_a2_filtered.py"""
import sys,os,numpy as np,pandas as pd
sys.path.insert(0,os.path.expanduser('~/xcell-xbulk-revision/lib')); import eric
from multiprocessing import Pool
R=os.path.expanduser('~/xcell-xbulk-revision/'); SIG=R+'eric/a2sig/'
GR=sorted({(r.split('_')[2],r.split('_')[3]) for r in os.listdir(eric.E) if r.startswith('sc_coexmats_') and r.endswith('__genes.i32')})
def F(level,ds,ct):
    a=SIG+f'{"full" if level=="sc" else "xsubject"}_{ds}_{ct}_d{0 if level=="sc" else -1}'
    return np.union1d(np.fromfile(a+'__pos.i32',dtype=np.int32),np.fromfile(a+'__neg.i32',dtype=np.int32))
def job(t):
    kind,lr,lt,ds_r,ds_t,ct=t; fr=0.01
    ref=np.fromfile(f'{eric.E}{lr}_coexmats_{ds_r}_{ct}__{fr}.i32',dtype=np.int32); test=np.fromfile(f'{eric.E}{lt}_coexmats_{ds_t}_{ct}__{fr}.i32',dtype=np.int32)
    Fr=F(lr,ds_r,ct); Ft=F(lt,ds_t,ct); ref2=np.intersect1d(ref,Fr,assume_unique=True); test2=np.intersect1d(test,Ft,assume_unique=True)
    a=int(np.intersect1d(ref2,test2,assume_unique=True).size); S=test2.size; Tin=int(np.intersect1d(ref2,Ft,assume_unique=True).size); N=int(Ft.size); b=S-a; c=Tin-a; d=N-S-c
    OR=(a/b)/(c/d) if min(b,c,d)>0 and a>0 else np.nan
    return dict(analysis=kind,reference_level=lr,test_level=lt,reference_dataset=ds_r,test_dataset=ds_t,cell_type=ct,a=a,b=b,c=c,d=d,N=N,n_ref_pos_filtered=int(ref2.size),n_test_pos_filtered=int(S),n_F_ref=int(Fr.size),n_F_test=int(Ft.size),odds_ratio_filtered=OR)
if __name__=='__main__':
    jobs=[]
    for ct in eric.CT:
        ds=[d for d,c in GR if c==ct]
        for r in ds:
            jobs.append(('preservation','sc','sbj',r,r,ct))
            for t in ds:
                if r!=t: jobs+=[('replicability','sc','sc',r,t,ct),('replicability','sbj','sbj',r,t,ct)]
    with Pool(16) as P: res=P.map(job,jobs,chunksize=8)
    d=pd.DataFrame(res); b=pd.read_csv(R+'eric/baseline/all_comparisons.tsv',sep='\t'); b=b[(b.frac==0.01)&b.analysis.isin(['replicability','preservation'])][['analysis','reference_level','test_level','reference_dataset','test_dataset','cell_type','odds_ratio']].rename(columns={'odds_ratio':'odds_ratio_unfiltered'})
    d=d.merge(b,on=['analysis','reference_level','test_level','reference_dataset','test_dataset','cell_type']); os.makedirs(R+'A2_permutation_null',exist_ok=True); d.to_csv(R+'A2_permutation_null/replicability_with_without_filter.tsv',sep='\t',index=False)
    g=d.groupby(['analysis','reference_level']).agg(n=('odds_ratio_filtered','size'),mean_OR_unfiltered=('odds_ratio_unfiltered','mean'),mean_OR_filtered=('odds_ratio_filtered','mean'),median_unfiltered=('odds_ratio_unfiltered','median'),median_filtered=('odds_ratio_filtered','median'),n_le_1_filtered=('odds_ratio_filtered',lambda x:int((x<=1).sum())),n_nan_filtered=('odds_ratio_filtered',lambda x:int(x.isna().sum())),mean_F_fraction=('n_F_test','mean')).reset_index()
    g.to_csv(R+'A2_permutation_null/replicability_with_without_filter_summary.tsv',sep='\t',index=False); print(g.round(3).to_string())
