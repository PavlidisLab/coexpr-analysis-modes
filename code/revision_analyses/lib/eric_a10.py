"""A10 (Nairuz's suggestion; Paul's ruling): common gene set per cell type used at ALL levels (xCell, xSubject, xBulk), ranks recomputed within the set.
Set(ct,rule) = genes that passed the existing per-group filters (i.e. are in the group's gene set) in >= k of the cell type's datasets (50% / 70% / 100%; n=7 or 6 datasets), intersected with the genes present in BOTH bulk matrices.
Each group then uses Set ∩ its own genes. Bulk uses Set ∩ bulk genes (bulk expression already log2+QN over the full bulk gene set; not re-normalised: flagged).
stages: sets | groups <workers> | bulk | compare <workers>"""
import sys,os,json,glob,math,numpy as np,pandas as pd
sys.path.insert(0,os.path.expanduser('~/xcell-xbulk-revision/lib')); import eric, eric_engine as EE, eric_subset as SS
R=os.path.expanduser('~/xcell-xbulk-revision/'); D=R+'A10_common_geneset/'; S=SS.S; G=eric.G; FR=[0.005,0.01,0.05,0.10]
RULES={'p50':0.5,'p70':0.7,'p100':1.0}
GR=sorted({tuple(os.path.basename(f).split('__')[0].split('_',1)) for f in glob.glob(EE.X+'*__p.i32')})
def bulk_genes(nm): return np.fromfile(f'{EE.X}bulk_{nm}__genes.i32',dtype=np.int32)
def sets():
    both=np.intersect1d(bulk_genes('rosmap'),bulk_genes('velmeshev')); rows=[]; out={}
    for ct in eric.CT:
        ds=[d for d,c in GR if c==ct]; cnt=np.zeros(G,int)
        for d in ds: cnt[eric.genes('sc',d,ct)]+=1
        for rn,thr in RULES.items():
            k=int(math.ceil(thr*len(ds)-1e-9)); sc_rule=np.where(cnt>=k)[0]; fin=np.intersect1d(sc_rule,both); out[(ct,rn)]=fin
            rows.append(dict(cell_type=ct,rule=rn,n_datasets=len(ds),k=k,genes_sc_rule=len(sc_rule),genes_final_after_bulk_intersection=len(fin),genes_in_bulk_both=len(both)))
    pd.DataFrame(rows).to_csv(D+'gene_set_sizes.tsv',sep='\t',index=False); return out
def do_group(t):
    ct,rn,ds=t; gs=sets_cache[(ct,rn)]; r=SS.run_subset(ds,ct,gs,f'A10_{rn}'); r['rule']=rn; return r
def do_bulk(t):
    ct,rn,nm=t; keep=sets_cache[(ct,rn)]; gi=bulk_genes(nm); dim=[int(x) for x in open(f'{EE.X}bulk_{nm}__dim.txt').read().split()]
    M=np.fromfile(f'{EE.X}bulk_{nm}__expr.f64',dtype=np.float64).reshape(dim[1],dim[0]).T   # genes x samples (column-major export)
    sel=np.where(np.isin(gi,keep))[0]; gs=gi[sel]; n=len(sel); i0,i1=EE.iu(n); ids=EE.ids_from(gs,i0,i1)
    with np.errstate(invalid='ignore',divide='ignore'): Cm=np.corrcoef(M[sel,:])
    v=Cm[i0,i1]; d=np.isnan(v)
    if d.any(): v[d]=np.nanmedian(v)
    SS.write(f'A10_{rn}','bk',f'{nm}_{ct}',gs,ids,EE.rank_avg_over_max(v)); return dict(cell_type=ct,rule=rn,bulk=nm,n_genes=n,n_samples=dim[1])
sets_cache={}
def path(tag,level,key,fr): return f'{S}{tag}/{level}_{key}__{fr}.i32'
def cmp_job(t):
    rule,fr,kind,lr,lt,kr,kt,ct=t; tag=f'A10_{rule}'
    ref=np.fromfile(path(tag,lr,kr,fr),dtype=np.int32); test=np.fromfile(path(tag,lt,kt,fr),dtype=np.int32)
    gr=np.fromfile(path(tag,lr,kr,'').replace('__.i32','__genes.i32'),dtype=np.int32); gt=np.fromfile(path(tag,lt,kt,'').replace('__.i32','__genes.i32'),dtype=np.int32)
    x=eric.counts(ref,test,gr,gt); x.update(rule=rule,frac=fr,analysis=kind,cell_type=ct,ref_level=lr,test_level=lt,reference=kr,test=kt,n_genes_ref=len(gr),n_genes_test=len(gt)); return x
if __name__=='__main__':
    from multiprocessing import Pool
    stage=sys.argv[1]; sets_cache.update(sets())
    if stage=='sets': print(pd.read_csv(D+'gene_set_sizes.tsv',sep='\t').to_string())
    if stage=='groups':
        with Pool(int(sys.argv[2])) as P: res=P.map(do_group,[(ct,rn,d) for rn in RULES for ct in eric.CT for d,c in GR if c==ct],chunksize=1)
        pd.DataFrame(res).to_csv(D+'group_gene_counts.tsv',sep='\t',index=False)
    if stage=='bulk':
        with Pool(6) as P: res=P.map(do_bulk,[(ct,rn,nm) for rn in RULES for ct in eric.CT for nm in ['rosmap','velmeshev']])
        pd.DataFrame(res).to_csv(D+'bulk_gene_counts.tsv',sep='\t',index=False)
    if stage=='compare':
        jobs=[]
        for rn in RULES:
            for fr in FR:
                for ct in eric.CT:
                    ds=[d for d,c in GR if c==ct]
                    for r in ds:
                        for t in ds:
                            if r!=t: jobs+=[(rn,fr,'replicability','sc','sc',f'{r}_{ct}',f'{t}_{ct}',ct),(rn,fr,'replicability','sbj','sbj',f'{r}_{ct}',f'{t}_{ct}',ct)]
                        jobs.append((rn,fr,'preservation_xsubject','sc','sbj',f'{r}_{ct}',f'{r}_{ct}',ct))
                        for nm in ['rosmap','velmeshev']:
                            if r==nm:
                                jobs+=[(rn,fr,'preservation_xbulk_from_sc','sc','bk',f'{r}_{ct}',f'{nm}_{ct}',ct),(rn,fr,'preservation_xbulk_from_sbj','sbj','bk',f'{r}_{ct}',f'{nm}_{ct}',ct)]
                    jobs+=[(rn,fr,'replicability_xbulk','bk','bk',f'rosmap_{ct}',f'velmeshev_{ct}',ct),(rn,fr,'replicability_xbulk','bk','bk',f'velmeshev_{ct}',f'rosmap_{ct}',ct)]
        jobs=[j for j in jobs if os.path.exists(path(f'A10_{j[0]}',j[3],j[5],j[1])) and os.path.exists(path(f'A10_{j[0]}',j[4],j[6],j[1]))]
        print(len(jobs),'comparisons',flush=True)
        with Pool(int(sys.argv[2])) as P: res=P.map(cmp_job,jobs,chunksize=8)
        d=pd.DataFrame(res); d.to_csv(D+'results.tsv',sep='\t',index=False)
        g=d.groupby(['analysis','rule','frac']).agg(n=('odds_ratio','size'),mean_OR=('odds_ratio','mean'),median_OR=('odds_ratio','median'),min_OR=('odds_ratio','min'),n_le_1=('odds_ratio',lambda x:int((x<=1).sum())),mean_F1=('F1','mean')).reset_index(); g.to_csv(D+'summary.tsv',sep='\t',index=False); print(g[g.frac==0.01].round(3).to_string())
