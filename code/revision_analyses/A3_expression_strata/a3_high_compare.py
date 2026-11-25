"""A3: replicability and preservation of the HIGH-expression stratum networks recomputed from scratch (genes with mean UMI/cell >= 1 in the group; ranks within the stratum). Eric's pair-universe definition."""
import sys,os,glob,numpy as np,pandas as pd
sys.path.insert(0,os.path.expanduser('~/xcell-xbulk-revision/lib')); import eric, eric_subset as SS
from multiprocessing import Pool
R=os.path.expanduser('~/xcell-xbulk-revision/'); S=SS.S+'A3_high/'
GR=sorted({(r.split('_')[2],r.split('_')[3]) for r in os.listdir(eric.E) if r.startswith('sc_coexmats_') and r.endswith('__genes.i32')})
def L(level,ds,ct,fr): return np.fromfile(f'{S}{level}_{ds}_{ct}__{fr}.i32',dtype=np.int32)
def Gn(level,ds,ct): return np.fromfile(f'{S}{level}_{ds}_{ct}__genes.i32',dtype=np.int32)
def job(t):
    kind,lr,lt,r,s,ct,fr=t; x=eric.counts(L(lr,r,ct,fr),L(lt,s,ct,fr),Gn(lr,r,ct),Gn(lt,s,ct)); x.update(analysis=kind,reference_level=lr,test_level=lt,reference_dataset=r,test_dataset=s,cell_type=ct,frac=fr,n_genes_ref=len(Gn(lr,r,ct)),n_genes_test=len(Gn(lt,s,ct))); return x
if __name__=='__main__':
    jobs=[]
    for fr in [0.005,0.01,0.05,0.10]:
        for ct in eric.CT:
            ds=[d for d,c in GR if c==ct and os.path.exists(f'{S}sc_{d}_{ct}__genes.i32')]
            for r in ds:
                jobs.append(('preservation','sc','sbj',r,r,ct,fr))
                for t in ds:
                    if r!=t: jobs+=[('replicability','sc','sc',r,t,ct,fr),('replicability','sbj','sbj',r,t,ct,fr)]
    with Pool(12) as P: res=P.map(job,jobs,chunksize=8)
    d=pd.DataFrame(res); d.to_csv(R+'A3_expression_strata/results_high_recomputed.tsv',sep='\t',index=False)
    g=d.groupby(['analysis','reference_level','frac']).agg(n=('odds_ratio','size'),mean_genes=('n_genes_ref','mean'),mean_OR=('odds_ratio','mean'),median_OR=('odds_ratio','median'),min_OR=('odds_ratio','min'),n_le_1=('odds_ratio',lambda x:int((x<=1).sum())),n_nan=('odds_ratio',lambda x:int(x.isna().sum())),mean_F1=('F1','mean')).reset_index(); g.to_csv(R+'A3_expression_strata/summary_high_recomputed.tsv',sep='\t',index=False); print(g.round(3).to_string())
