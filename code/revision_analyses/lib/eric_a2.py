"""A2 (R1-M2): signal-to-noise vs a permutation null at matched sample size, using the same samplings as A1 (identical SeedSequence -> identical cells).
Real statistic = raw Pearson coexpression (xCell: mean over subjects of per-subject r on the sampled cells; xSubject/xBulk: Pearson across samples of Eric's transformed data).
Null = each gene's expression vector permuted independently (within subject for xCell-S2, across the pooled cells for S1, across subjects/samples for xSubject/xBulk), P permutations.
Reported: fraction of real positive correlations above the null maximum (over all pairs and permutations) and of real negative correlations below the null minimum.
usage: eric_a2.py <level S2|S1|xsubject|xbulk> <workers> <ndraws> <nperm>"""
import sys,os,json,time,numpy as np,scipy.sparse as sp
sys.path.insert(0,os.path.expanduser('~/xcell-xbulk-revision/lib')); import eric_engine as EE
O=os.path.expanduser('~/xcell-xbulk-revision/A2_permutation_null/'); SEED=20261015
def stats(real,nullmax,nullmin,extra,nullsample=None):
    r=real[~np.isnan(real)]; pos=r[r>0]; neg=r[r<0]
    if nullsample is not None and nullsample.size:
        for q in (0.999,0.9999,0.99999): extra[f'frac_pos_above_null_q{q}']=float((pos>np.quantile(nullsample,q)).mean()) if pos.size else np.nan; extra[f'frac_neg_below_null_q{1-q:.5f}']=float((neg<np.quantile(nullsample,1-q)).mean()) if neg.size else np.nan
    d=dict(n_pairs=int(r.size),n_pos=int(pos.size),n_neg=int(neg.size),frac_pos_above_null_max=float((pos>nullmax.max()).mean()) if pos.size else np.nan,frac_neg_below_null_min=float((neg<nullmin.min()).mean()) if neg.size else np.nan,
           null_max_overall=float(nullmax.max()),null_max_mean_perm=float(nullmax.mean()),null_min_overall=float(nullmin.min()),null_min_mean_perm=float(nullmin.mean()),real_max=float(r.max()),real_min=float(r.min()),
           frac_pos_above_mean_perm_max=float((pos>nullmax.mean()).mean()) if pos.size else np.nan); d.update(extra); return d
SIG=os.path.expanduser('~/xcell-xbulk-revision/eric/a2sig/')
def save_sig(level,key,draw,ids,real,nmax,nmin):
    os.makedirs(SIG,exist_ok=True); ok=~np.isnan(real)
    np.sort(ids[ok&(real>nmax.max())]).astype(np.int32).tofile(f'{SIG}{level}_{key}_d{draw}__pos.i32'); np.sort(ids[ok&(real<nmin.min())]).astype(np.int32).tofile(f'{SIG}{level}_{key}_d{draw}__neg.i32')
def xcell(level,ds,ct,draw,P):
    key=f'{ds}_{ct}'; C,gi,sbj=EE.load(ds,ct); n=C.shape[0]; cp=EE.cpm(C).tocsc(); subs,cnt=np.unique(sbj,return_counts=True); keep=[s for s,c in zip(subs,cnt) if c>=20]; N=len(keep)
    gidx=sorted({tuple(os.path.basename(f).split('__')[0].split('_',1)) for f in __import__('glob').glob(EE.X+'*__p.i32')}).index((ds,ct))
    rng=np.random.default_rng(np.random.SeedSequence([SEED,1 if level=='S1' else 2,gidx,draw])); i0,i1=EE.iu(n); blocks=[]
    if level=='full':
        for s in keep: blocks.append(cp[:,np.where(sbj==s)[0]].T.toarray())
    elif level=='S2':
        for s in keep:
            idx=np.where(sbj==s)[0]
            if idx.size>N: idx=np.sort(rng.choice(idx,size=N,replace=False))
            blocks.append(cp[:,idx].T.toarray())
    else:
        pool=np.where(np.isin(sbj,keep))[0]; idx=np.sort(rng.choice(pool,size=N,replace=False)); blocks.append(cp[:,idx].T.toarray())
    prng=np.random.default_rng(np.random.SeedSequence([SEED,99,gidx,draw]))
    def agg(bl,perm):
        a=np.zeros(len(i0)); c=np.zeros(len(i0))
        for X in bl:
            if perm: X=prng.permuted(X,axis=0)
            v=EE.corr_qlc(X)[i0,i1]; d=~np.isnan(v); a[d]+=v[d]; c+=d
        with np.errstate(invalid='ignore',divide='ignore'): return np.where(c>0,a/c,np.nan)
    real=agg(blocks,False); nmax=np.zeros(P); nmin=np.zeros(P); samp=[]
    for p in range(P): z=agg(blocks,True); nmax[p]=np.nanmax(z); nmin[p]=np.nanmin(z); samp.append(z[prng.integers(0,len(z),size=min(len(z),300000))])
    ids=EE.ids_from(gi,i0,i1); save_sig(level,key,draw,ids,real,nmax,nmin)
    return stats(real,nmax,nmin,dict(level=level,dataset=ds,cell_type=ct,draw=draw,n_subjects=N,n_genes=int(n),n_perm=P,cells_per_corr=int(blocks[0].shape[0])),np.concatenate(samp))
def samples_matrix(level,ds,ct):
    from eric_xsbj import lowtri_rank  # noqa (module import only to share conventions)
def xsubject(ds,ct,P):
    C,gi,sbj=EE.load(ds,ct); cp=EE.cpm(C).tocsc(); subs,cnt=np.unique(sbj,return_counts=True)
    S_=sp.csc_matrix((np.ones(len(sbj)),(np.arange(len(sbj)),np.searchsorted(subs,sbj))),shape=(len(sbj),len(subs))); M=EE.quantile_normalize(np.log2(np.asarray((cp@S_).todense())/cnt+1))[:,cnt>20]
    return corrnull(M,P,dict(level='xsubject',dataset=ds,cell_type=ct,draw=-1,n_subjects=int(M.shape[1]),n_genes=int(M.shape[0]),n_perm=P,cells_per_corr=int(M.shape[1])),SEED+7,gi=gi,key=f'{ds}_{ct}')
def xbulk(nm,P,ndraw_match):
    dim=[int(x) for x in open(f'{EE.X}bulk_{nm}__dim.txt').read().split()]; M=np.fromfile(f'{EE.X}bulk_{nm}__expr.f64',dtype=np.float64).reshape(dim[1],dim[0]).T
    gib=np.fromfile(f'{EE.X}bulk_{nm}__genes.i32',dtype=np.int32)
    out=[corrnull(M,P,dict(level='xbulk_full',dataset=nm,cell_type='',draw=-1,n_subjects=dim[1],n_genes=dim[0],n_perm=P,cells_per_corr=dim[1]),SEED+11,gi=gib,key=nm)]
    Nm={'rosmap':48,'velmeshev':41}[nm]   # matched to the number of subjects of the dataset's single-cell groups
    for d in range(ndraw_match):
        rng=np.random.default_rng(np.random.SeedSequence([SEED,5,d])); idx=np.sort(rng.choice(dim[1],size=Nm,replace=False))
        out.append(corrnull(M[:,idx],P,dict(level='xbulk_matched',dataset=nm,cell_type='',draw=d,n_subjects=Nm,n_genes=dim[0],n_perm=P,cells_per_corr=Nm),SEED+13+d))
    return out
def corrnull(M,P,extra,seed,gi=None,key=None):
    n=M.shape[0]; i0,i1=EE.iu(n)
    with np.errstate(invalid='ignore',divide='ignore'): real=np.corrcoef(M)[i0,i1]
    rng=np.random.default_rng(seed); nmax=np.zeros(P); nmin=np.zeros(P); samp=[]
    for p in range(P):
        with np.errstate(invalid='ignore',divide='ignore'): z=np.corrcoef(rng.permuted(M,axis=1))[i0,i1]
        nmax[p]=np.nanmax(z); nmin[p]=np.nanmin(z); samp.append(z[rng.integers(0,len(z),size=min(len(z),300000))])
    if gi is not None: save_sig(extra['level'],key,extra['draw'],EE.ids_from(gi,i0,i1),real,nmax,nmin)
    return stats(real,nmax,nmin,extra,np.concatenate(samp)[~np.isnan(np.concatenate(samp))])
def _w(t):
    kind=t[0]; t0=time.time()
    r=xcell(*t) if kind in('S1','S2','full') else (xsubject(*t[1:]) if kind=='xsubject' else xbulk(*t[1:]))
    rr=r if isinstance(r,list) else [r]
    for x in rr: x['secs']=round(time.time()-t0)
    print(json.dumps(rr),flush=True); return rr
if __name__=='__main__':
    from multiprocessing import Pool
    import glob,pandas as pd
    level=sys.argv[1]; nw=int(sys.argv[2]); nd=int(sys.argv[3]); P=int(sys.argv[4])
    groups=sorted({tuple(os.path.basename(f).split('__')[0].split('_',1)) for f in glob.glob(EE.X+'*__p.i32')}); sizes={g:os.path.getsize(f'{EE.X}{g[0]}_{g[1]}__x.f64') for g in groups}; groups.sort(key=lambda g:-sizes[g])
    if level in('S1','S2','full'): tasks=[(level,g[0],g[1],d,P) for d in range(nd) for g in groups]
    elif level=='xsubject': tasks=[('xsubject',g[0],g[1],P) for g in groups]
    else: tasks=[('xbulk','rosmap',P,nd),('xbulk','velmeshev',P,nd)]
    with Pool(nw) as Pl: res=Pl.map(_w,tasks,chunksize=1)
    pd.DataFrame([x for rr in res for x in rr]).to_csv(O+f'results_{level}.tsv',sep='\t',index=False); print('ALLDONE',flush=True)
