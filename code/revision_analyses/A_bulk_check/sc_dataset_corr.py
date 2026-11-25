"""Dataset-to-dataset correlation of per-gene mean expression, single-cell groups (Eric's pipeline objects, exported raw counts).
Per cell type, for each dataset pair: gene means over cells of log2(CPM+1) (CPM over the group's own filtered genes, as pipeline)
across genes present in both datasets' group gene sets. Also mean of CPM (linear). Pearson and Spearman."""
import sys; sys.path.insert(0,'lib')
import numpy as np, pandas as pd, itertools, eric_engine as EE
from scipy.stats import spearmanr
q=pd.read_csv(EE.X+'export_qc.tsv',sep='\t'); keys=[k for k in q.key if not k.startswith('bulk')]
D={}
for k in keys:
    ds,ct=k.split('_',1); C,gi,s=EE.load(ds,ct); P=EE.cpm(C).tocsc()
    L=P.copy(); L.data=np.log2(L.data+1)
    D[(ds,ct)]=(gi,np.asarray(L.mean(1)).ravel(),np.asarray(P.mean(1)).ravel())
rows=[]
for ct in sorted({c for _,c in D}):
    for a,b in itertools.combinations(sorted([d for d,c in D if c==ct]),2):
        ga,la,pa=D[(a,ct)]; gb,lb,pb=D[(b,ct)]; sh,ia,ib=np.intersect1d(ga,gb,return_indices=True)
        rows.append(dict(cell_type=ct,ds_a=a,ds_b=b,n_shared=len(sh),pearson_log2cpm_mean=np.corrcoef(la[ia],lb[ib])[0,1],
            spearman=spearmanr(la[ia],lb[ib])[0],pearson_cpm_mean=np.corrcoef(pa[ia],pb[ib])[0,1]))
r=pd.DataFrame(rows); r.to_csv('A_bulk_check/sc_dataset_pairs.tsv',sep='\t',index=False)
print(r.groupby('cell_type')[['n_shared','pearson_log2cpm_mean','spearman','pearson_cpm_mean']].mean().round(3)); print(r[['pearson_log2cpm_mean','spearman','pearson_cpm_mean']].agg(['mean','min','max']).round(3), len(r))
