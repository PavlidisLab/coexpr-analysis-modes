"""One pass over a group's subjects (Eric's xCell definition) computing variants from the same per-subject Pearson matrices:
 std      : Eric: per-subject lower-tri vector, NaN->median, rank(ties=min)/max, mean over subjects (>=20 cells), rank(average)/max           [validate vs saved sc_coexmats]
 A7       : NO imputation: per-subject rank over defined pairs only; mean over subjects where defined; final rank over pairs defined in >=1 subject
 A8_wrank : as std but subjects weighted by their number of cells
 rawmean  : unweighted mean raw Pearson r over subjects where defined;  A8_r: weighted by cells;  A8_z: Fisher-z mean weighted by (cells-3)
 rawsign  : pairs with mean raw r > 0
Edges (value rank >= 1-frac) -> eric/edges_var/<variant>__<dataset>_<celltype>__<frac>.i32"""
import sys,os,json,time,numpy as np
sys.path.insert(0,os.path.expanduser('~/xcell-xbulk-revision/lib')); import eric_engine as EE
V=os.path.expanduser('~/xcell-xbulk-revision/eric/edges_var/'); FR=[0.005,0.01,0.05,0.10]
def save(variant,key,ids_by_frac):
    for fr,a in ids_by_frac.items(): a.astype(np.int32).tofile(f'{V}{variant}__{key}__{fr}.i32')
def run(ds,ct):
    key=f'{ds}_{ct}'; done=f'{V}done__{key}.json'
    if os.path.exists(done): return json.load(open(done))
    t0=time.time(); C,gi,sbj=EE.load(ds,ct); n=C.shape[0]; cp=EE.cpm(C).tocsc()
    subs,counts=np.unique(sbj,return_counts=True); keep=[s for s,c in zip(subs,counts) if c>=20]
    i0,i1=EE.iu(n); L=len(i0); ids=EE.ids_from(gi,i0,i1)
    acc=np.zeros(L); accw=np.zeros(L); a7=np.zeros(L); cnt=np.zeros(L,np.int16); sr=np.zeros(L); swr=np.zeros(L); swt=np.zeros(L); sz=np.zeros(L); swz=np.zeros(L); sumw=0.0; ncell=[]
    for s in keep:
        idx=np.where(sbj==s)[0]; ns=idx.size; M=cp[:,idx].T.toarray()
        Cm=EE.corr_qlc(M); v=Cm[i0,i1]; del Cm,M; d=~np.isnan(v)
        vi=v.copy(); vi[~d]=np.nanmedian(v); r=EE.rank_min_over_max(vi); acc+=r; accw+=ns*r; sumw+=ns
        if d.all(): rk=EE.rank_min_over_max(v)
        else: rk=np.zeros(L); rk[d]=EE.rank_min_over_max(v[d])
        a7[d]+=rk[d]; cnt+=d
        vz=np.where(d,v,0.0); sr+=vz; swr+=ns*vz; swt+=ns*d
        z=np.arctanh(np.clip(vz,-1+1e-9,1-1e-9)); sz+=(ns-3)*z*d; swz+=(ns-3)*d; ncell.append(int(ns)); del v,vi,r,rk,vz,z,d
    out=dict(dataset=ds,cell_type=ct,n_genes=int(n),n_pairs=int(L),n_subjects=len(keep),cells_used=int(sum(ncell)),cells_per_subject=ncell)
    out['n_pairs_never_defined']=int((cnt==0).sum())
    save('std',key,EE.edges(ids,EE.rank_avg_over_max(acc/len(keep)),FR)); save('A8_wrank',key,EE.edges(ids,EE.rank_avg_over_max(accw/sumw),FR))
    dd=cnt>0
    with np.errstate(invalid='ignore',divide='ignore'):
        rk=np.zeros(L); rk[dd]=EE.rank_avg_over_max((a7/np.maximum(cnt,1))[dd]); save('A7',key,EE.edges(ids,rk,FR,valid=dd))
        rm=np.zeros(L); mr=(sr/np.maximum(cnt,1)); rm[dd]=EE.rank_avg_over_max(mr[dd]); save('rawmean',key,EE.edges(ids,rm,FR,valid=dd))
        np.sort(ids[dd&(mr>0)]).astype(np.int32).tofile(f'{V}rawsign__{key}__pos.i32'); out['n_rawsign_pos']=int((dd&(mr>0)).sum())
        w=swt>0; x=np.zeros(L); x[w]=EE.rank_avg_over_max((swr/np.maximum(swt,1e-12))[w]); save('A8_r',key,EE.edges(ids,x,FR,valid=w))
        w=swz>0; x=np.zeros(L); x[w]=EE.rank_avg_over_max(np.tanh(sz/np.maximum(swz,1e-12))[w]); save('A8_z',key,EE.edges(ids,x,FR,valid=w))
    # validation vs Eric's saved sc edges
    ref=np.fromfile(f'{os.path.expanduser("~/xcell-xbulk-revision/eric/edges/")}sc_coexmats_{key}__0.01.i32',dtype=np.int32); mine=np.fromfile(f'{V}std__{key}__0.01.i32',dtype=np.int32)
    out['std_vs_saved_jaccard_1pct']=float(np.intersect1d(ref,mine).size/np.union1d(ref,mine).size); out['secs']=round(time.time()-t0); json.dump(out,open(done,'w')); return out
def _w(a):
    r=run(*a); print(json.dumps({k:v for k,v in r.items() if k!='cells_per_subject'}),flush=True); return r
if __name__=='__main__':
    from multiprocessing import Pool
    import glob
    keys=sorted({os.path.basename(f).split('__')[0] for f in glob.glob(EE.X+'*__p.i32')})
    groups=[tuple(k.split('_',1)) for k in keys]
    sizes={g:os.path.getsize(f'{EE.X}{g[0]}_{g[1]}__x.f64') for g in groups}; groups.sort(key=lambda g:-sizes[g])
    only=sys.argv[2:]; 
    if only: groups=[g for g in groups if f'{g[0]}_{g[1]}' in only]
    print(len(groups),'groups',flush=True)
    with Pool(int(sys.argv[1])) as P: list(P.imap_unordered(_w,groups))
    print('ALLDONE',flush=True)
