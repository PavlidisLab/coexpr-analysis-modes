"""A1 (R1-M1) matched-sample-size resampling on Eric's xCell definition. Seed 20261015, 50 draws (default).
 S2 (primary): for every subject with >=20 cells keep at most N_subj cells (N_subj = the group's number of subjects with >=20 cells; random without replacement), then run the UNCHANGED xCell engine
               (per-subject Pearson, NaN->median, rank(min)/max, mean over subjects, rank(average)/max). Subjects with <= N_subj cells are unchanged.
 S1 (secondary): pool the cells of those subjects, draw N_subj cells at random, one Pearson matrix over them (no subject step), NaN->median, rank(average)/max.
Only top-edge lists (0.5/1/5/10%) are kept per draw: eric/a1/<S1|S2>/<ds>_<ct>__d<NN>__<frac>.i32 (+ ../a1_meta). usage: eric_a1.py <S1|S2> <workers> [ndraws]"""
import sys,os,json,time,hashlib,numpy as np
sys.path.insert(0,os.path.expanduser('~/xcell-xbulk-revision/lib')); import eric_engine as EE
A=os.path.expanduser('~/xcell-xbulk-revision/eric/a1/'); FR=[0.005,0.01,0.05,0.10]; SEED=20261015
def one(task):
    scheme,ds,ct,d,gidx=task; key=f'{ds}_{ct}'; out=f'{A}{scheme}/{key}__d{d:02d}__0.1.i32'
    if os.path.exists(out): return None
    t0=time.time(); C,gi,sbj=EE.load(ds,ct); n=C.shape[0]; cp=EE.cpm(C).tocsc(); subs,cnt=np.unique(sbj,return_counts=True); keep=[s for s,c in zip(subs,cnt) if c>=20]; N=len(keep)
    rng=np.random.default_rng(np.random.SeedSequence([SEED,1 if scheme=='S1' else 2,gidx,d]))
    i0,i1=EE.iu(n); ids=EE.ids_from(gi,i0,i1)
    if scheme=='S2':
        acc=np.zeros(len(i0)); used=0
        for s in keep:
            idx=np.where(sbj==s)[0]
            if idx.size>N: idx=np.sort(rng.choice(idx,size=N,replace=False))
            v=EE.corr_qlc(cp[:,idx].T.toarray())[i0,i1]; d_=np.isnan(v)
            if d_.any(): v[d_]=np.nanmedian(v)
            acc+=EE.rank_min_over_max(v); used+=idx.size
        r=EE.rank_avg_over_max(acc/len(keep))
    else:
        pool=np.where(np.isin(sbj,keep))[0]; idx=np.sort(rng.choice(pool,size=N,replace=False)); used=idx.size
        v=EE.corr_qlc(cp[:,idx].T.toarray())[i0,i1]; d_=np.isnan(v)
        if d_.any(): v[d_]=np.nanmedian(v)
        r=EE.rank_avg_over_max(v)
    os.makedirs(f'{A}{scheme}',exist_ok=True)
    for fr in FR: np.sort(ids[r>=1-fr]).astype(np.int32).tofile(f'{A}{scheme}/{key}__d{d:02d}__{fr}.i32')
    return dict(scheme=scheme,dataset=ds,cell_type=ct,draw=d,n_subjects=N,cells_used=int(used),secs=round(time.time()-t0))
def _w(t):
    r=one(t)
    if r: print(json.dumps(r),flush=True)
    return r
if __name__=='__main__':
    from multiprocessing import Pool
    import glob
    scheme=sys.argv[1]; nw=int(sys.argv[2]); nd=int(sys.argv[3]) if len(sys.argv)>3 else 50
    groups=sorted({tuple(os.path.basename(f).split('__')[0].split('_',1)) for f in glob.glob(EE.X+'*__p.i32')})
    gidx={g:i for i,g in enumerate(groups)}; sizes={g:os.path.getsize(f'{EE.X}{g[0]}_{g[1]}__x.f64') for g in groups}
    tasks=[(scheme,g[0],g[1],d,gidx[g]) for d in range(nd) for g in sorted(groups,key=lambda g:-sizes[g])]
    if len(sys.argv)>4 and sys.argv[4]=='rev': tasks=tasks[::-1]
    os.makedirs(f'{A}{scheme}',exist_ok=True); print(scheme,len(tasks),'tasks',flush=True)
    with Pool(nw) as P: list(P.imap_unordered(_w,tasks,chunksize=1))
    print('ALLDONE',flush=True)
