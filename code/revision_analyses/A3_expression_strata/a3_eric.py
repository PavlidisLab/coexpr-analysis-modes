"""A3 on Eric's pipeline. (1) per-group mean raw UMI/cell per gene (from exported raw counts), genes with mean >=1 ('high') and <1 ('low'); (2) replicability/preservation within strata using the
BASELINE edges restricted to stratum pairs (not re-ranked; universe = stratum pairs shared by both groups); (3) HIGH stratum networks recomputed from scratch (ranks within the stratum) for xCell and xSubject."""
import sys,os,json,glob,numpy as np,pandas as pd
sys.path.insert(0,os.path.expanduser('~/xcell-xbulk-revision/lib')); import eric, eric_engine as EE, eric_subset as SS
from multiprocessing import Pool
R=os.path.expanduser('~/xcell-xbulk-revision/'); D=R+'A3_expression_strata/'; G=eric.G
GR=sorted({tuple(os.path.basename(f).split('__')[0].split('_',1)) for f in glob.glob(EE.X+'*__p.i32')})
def umi(g):
    ds,ct=g; C,gi,sbj=EE.load(ds,ct); m=np.asarray(C.sum(1)).ravel()/C.shape[1]; tot=float(C.sum()); subs,cnt=np.unique(sbj,return_counts=True)
    np.save(f'{D}meanumi__{ds}_{ct}.npy',m.astype(np.float32)); np.save(f'{D}genes__{ds}_{ct}.npy',gi)
    return dict(dataset=ds,cell_type=ct,cells=int(C.shape[1]),subjects_all=len(subs),subjects_ge20=int((cnt>=20).sum()),subjects_gt20=int((cnt>20).sum()),total_umi=tot,mean_umi_per_cell=tot/C.shape[1],n_genes=int(C.shape[0]),genes_mean_ge1=int((m>=1).sum()),genes_mean_lt1=int((m<1).sum()))
def stratified(t):
    kind,lr,lt,ds_r,ds_t,ct,st=t; fr=0.01
    ref=np.fromfile(f'{eric.E}{lr}_coexmats_{ds_r}_{ct}__{fr}.i32',dtype=np.int32); test=np.fromfile(f'{eric.E}{lt}_coexmats_{ds_t}_{ct}__{fr}.i32',dtype=np.int32)
    gr=eric.genes(lr,ds_r,ct); gt=eric.genes(lt,ds_t,ct); mr=np.load(f'{D}meanumi__{ds_r}_{ct}.npy'); mt=np.load(f'{D}meanumi__{ds_t}_{ct}.npy')
    hr=np.zeros(G,bool); hr[gr[mr>=1]]=True; ht=np.zeros(G,bool); ht[gt[mt>=1]]=True; ar=np.zeros(G,bool); ar[gr]=True; at=np.zeros(G,bool); at[gt]=True
    m={'all':ar&at,'high':hr&ht,'low':(ar&~hr)&(at&~ht)}[st]
    ref2=ref[m[ref//G]&m[ref%G]]; test2=test[m[test//G]&m[test%G]]; k=int(m.sum()); N=k*(k-1)//2
    a=int(np.intersect1d(ref2,test2,assume_unique=True).size); S=test2.size; T=ref2.size; b=S-a; c=T-a; d=N-S-c
    OR=(a/b)/(c/d) if min(b,c,d)>0 and a>0 else np.nan
    return dict(analysis=kind,reference_level=lr,test_level=lt,reference_dataset=ds_r,test_dataset=ds_t,cell_type=ct,stratum=st,n_genes_stratum=k,a=a,b=b,c=c,d=d,N=N,n_ref_pos=int(T),n_test_pos=int(S),odds_ratio=OR,F1=2*a/(2*a+b+c) if (2*a+b+c) else np.nan,jaccard=a/(a+b+c) if (a+b+c) else np.nan)
def sub_high(g):
    ds,ct=g; m=np.load(f'{D}meanumi__{ds}_{ct}.npy'); gi=np.load(f'{D}genes__{ds}_{ct}.npy'); return SS.run_subset(ds,ct,gi[m>=1],'A3_high')
if __name__=='__main__':
    with Pool(12) as P:
        rows=P.map(umi,GR); pd.DataFrame(rows).to_csv(D+'group_strata_counts.tsv',sep='\t',index=False)
        jobs=[]
        for ct in eric.CT:
            ds=[d for d,c in GR if c==ct]
            for r in ds:
                for t in ds:
                    for st in ['all','high','low']:
                        if r!=t: jobs+=[('replicability','sc','sc',r,t,ct,st),('replicability','sbj','sbj',r,t,ct,st)]
                        else: jobs.append(('preservation','sc','sbj',r,t,ct,st))
        res=P.map(stratified,jobs,chunksize=8); d=pd.DataFrame(res); d.to_csv(D+'results_stratified_baseline_edges.tsv',sep='\t',index=False)
        g=d.groupby(['analysis','reference_level','stratum']).agg(n=('odds_ratio','size'),mean_genes=('n_genes_stratum','mean'),mean_OR=('odds_ratio','mean'),median_OR=('odds_ratio','median'),min_OR=('odds_ratio','min'),n_le_1=('odds_ratio',lambda x:int((x<=1).sum())),n_nan=('odds_ratio',lambda x:int(x.isna().sum())),mean_F1=('F1','mean')).reset_index(); g.to_csv(D+'summary_stratified_baseline_edges.tsv',sep='\t',index=False); print(g.round(3).to_string(),flush=True)
        out=P.map(sub_high,GR); json.dump(out,open(D+'high_recompute_qc.json','w'),indent=1); print('high recompute done',flush=True)
