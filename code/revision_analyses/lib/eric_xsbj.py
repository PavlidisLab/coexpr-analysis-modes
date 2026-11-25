"""xSubject variants on Eric's definitions (aggregSbj -> logTrans -> quantileNormalize -> computeCoexmat -> normalizeCoexmat), from exported raw counts.
 std              : mean over cells of CPM, log2(x+1), quantile-normalise over ALL subjects of the group, keep subjects with >20 cells, Pearson across subjects, NaN->median, rank/max   [validate vs saved sbj_coexmats]
 A4_sum_cellcpm   : SUM of cell CPM instead of mean, then identical downstream
 A4_sum_counts_cpm: SUM of raw counts per subject, CPM-normalised over the group's genes, then identical downstream
 A4_nolog_noqn / A4_log_noqn / A4_mean_nolog_qn: effect of the undocumented log2 and quantile-normalisation steps on the mean-CPM pseudobulk (not requested; diagnostic)
Edges to eric/edges_var/sbj_<variant>__<ds>_<ct>__<frac>.i32 ; sign positives sbj_<variant>__..__pos.i32"""
import sys,os,json,time,numpy as np,scipy.sparse as sp
sys.path.insert(0,os.path.expanduser('~/xcell-xbulk-revision/lib')); import eric_engine as EE
V=os.path.expanduser('~/xcell-xbulk-revision/eric/edges_var/'); FR=[0.005,0.01,0.05,0.10]
def lowtri_rank(M,gi,sel):
    n=M.shape[0]; i0,i1=EE.iu(n); ids=EE.ids_from(gi,i0,i1)
    with np.errstate(invalid='ignore',divide='ignore'): Cm=np.corrcoef(M[:,sel])
    v=Cm[i0,i1]; d=np.isnan(v); 
    if d.any(): v[d]=np.nanmedian(v)
    return ids,EE.rank_avg_over_max(v),v
def run(ds,ct):
    key=f'{ds}_{ct}'; done=f'{V}done_sbj__{key}.json'
    if os.path.exists(done): return json.load(open(done))
    t0=time.time(); C,gi,sbj=EE.load(ds,ct); cp=EE.cpm(C).tocsc(); subs,cnt=np.unique(sbj,return_counts=True)
    S=sp.csc_matrix((np.ones(len(sbj)),(np.arange(len(sbj)),np.searchsorted(subs,sbj))),shape=(len(sbj),len(subs)))
    sum_cpm=np.asarray((cp@S).todense()); mean_cpm=sum_cpm/cnt; sum_counts=np.asarray((C@S).todense()); sum_counts_cpm=sum_counts/np.maximum(sum_counts.sum(0,keepdims=True),1)*1e6
    sel=np.where(cnt>20)[0]; out=dict(dataset=ds,cell_type=ct,n_genes=int(C.shape[0]),n_subjects_all=len(subs),n_subjects_used=len(sel),cells_used=int(cnt[sel].sum()))
    variants={'std':EE.quantile_normalize(np.log2(mean_cpm+1)),'A4_sum_cellcpm':EE.quantile_normalize(np.log2(sum_cpm+1)),'A4_sum_counts_cpm':EE.quantile_normalize(np.log2(sum_counts_cpm+1)),
              'A4_log_noqn':np.log2(mean_cpm+1),'A4_nolog_noqn':mean_cpm,'A4_mean_nolog_qn':EE.quantile_normalize(mean_cpm)}
    for nm,M in variants.items():
        ids,r,v=lowtri_rank(M,gi,sel)
        for fr in FR: np.sort(ids[r>=1-fr]).astype(np.int32).tofile(f'{V}sbj_{nm}__{key}__{fr}.i32')
        np.sort(ids[v>0]).astype(np.int32).tofile(f'{V}sbj_{nm}__{key}__pos.i32'); out['n_sign_pos_'+nm]=int((v>0).sum())
    ref=np.fromfile(os.path.expanduser('~/xcell-xbulk-revision/eric/edges/')+f'sbj_coexmats_{key}__0.01.i32',dtype=np.int32); mine=np.fromfile(f'{V}sbj_std__{key}__0.01.i32',dtype=np.int32)
    out['std_vs_saved_jaccard_1pct']=float(np.intersect1d(ref,mine).size/np.union1d(ref,mine).size)
    for nm in variants:
        if nm!='std':
            o=np.fromfile(f'{V}sbj_{nm}__{key}__0.01.i32',dtype=np.int32); out['jaccard_vs_std_'+nm]=float(np.intersect1d(mine,o).size/np.union1d(mine,o).size)
    out['secs']=round(time.time()-t0); json.dump(out,open(done,'w')); return out
def _w(a):
    r=run(*a); print(json.dumps(r),flush=True); return r
if __name__=='__main__':
    from multiprocessing import Pool
    import glob
    groups=sorted({tuple(os.path.basename(f).split('__')[0].split('_',1)) for f in glob.glob(EE.X+'*__p.i32')})
    only=sys.argv[2:]
    if only: groups=[g for g in groups if f'{g[0]}_{g[1]}' in only]
    print(len(groups),'groups',flush=True)
    with Pool(int(sys.argv[1])) as P: list(P.imap_unordered(_w,groups))
    print('ALLDONE',flush=True)
