"""Recompute Eric's xCell and xSubject networks WITHIN a gene subset (ranks recomputed over the subset's pairs; CPM library size still taken over the group's full retained gene set, as in the pipeline;
for xSubject the quantile normalisation is recomputed over the subset genes). Used by A3 (high-expression stratum) and A10 (common gene sets).
run_subset(ds,ct,keep_global_idx,tag,levels) -> writes eric/edges_sub/<tag>/<level>_<ds>_<ct>__<frac>.i32 and __genes.i32 (subset genes, global idx)."""
import sys,os,json,time,numpy as np,scipy.sparse as sp
sys.path.insert(0,os.path.expanduser('~/xcell-xbulk-revision/lib')); import eric_engine as EE
S=os.path.expanduser('~/xcell-xbulk-revision/eric/edges_sub/'); FR=[0.005,0.01,0.05,0.10]
def write(tag,level,key,gsub,ids,r):
    d=f'{S}{tag}/'; os.makedirs(d,exist_ok=True); gsub.astype(np.int32).tofile(f'{d}{level}_{key}__genes.i32')
    for fr in FR: np.sort(ids[r>=1-fr]).astype(np.int32).tofile(f'{d}{level}_{key}__{fr}.i32')
def run_subset(ds,ct,keep_global,tag,levels=('sc','sbj')):
    key=f'{ds}_{ct}'; C,gi,sbj=EE.load(ds,ct); cp=EE.cpm(C).tocsc()
    sel=np.where(np.isin(gi,keep_global))[0]; out=dict(dataset=ds,cell_type=ct,tag=tag,n_group_genes=len(gi),n_subset_genes=int(len(sel)))
    if len(sel)<30: out['skipped']='too few genes'; return out
    gs=gi[sel]; n=len(sel); i0,i1=EE.iu(n); ids=EE.ids_from(gs,i0,i1); cps=cp[sel,:].tocsc(); subs,cnt=np.unique(sbj,return_counts=True)
    if 'sc' in levels:
        keep=[s for s,c in zip(subs,cnt) if c>=20]; acc=np.zeros(len(i0))
        for s in keep:
            idx=np.where(sbj==s)[0]; v=EE.corr_qlc(cps[:,idx].T.toarray())[i0,i1]; d=np.isnan(v)
            if d.any(): v[d]=np.nanmedian(v)
            acc+=EE.rank_min_over_max(v)
        write(tag,'sc',key,gs,ids,EE.rank_avg_over_max(acc/len(keep))); out['sc_subjects']=len(keep)
    if 'sbj' in levels:
        S_=sp.csc_matrix((np.ones(len(sbj)),(np.arange(len(sbj)),np.searchsorted(subs,sbj))),shape=(len(sbj),len(subs)))
        mean=np.asarray((cps@S_).todense())/cnt; M=EE.quantile_normalize(np.log2(mean+1)); use=np.where(cnt>20)[0]
        with np.errstate(invalid='ignore',divide='ignore'): Cm=np.corrcoef(M[:,use])
        v=Cm[i0,i1]; d=np.isnan(v)
        if d.any(): v[d]=np.nanmedian(v)
        write(tag,'sbj',key,gs,ids,EE.rank_avg_over_max(v)); out['sbj_subjects']=int(len(use))
    return out
