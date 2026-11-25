import sys,os,numpy as np
sys.path.insert(0,os.path.expanduser('~/xcell-xbulk-revision/lib')); import eric_engine as EE
from scipy.stats import rankdata
ds,ct='nagy','microglia'; key=f'{ds}_{ct}'
C,gi,sbj=EE.load(ds,ct); n=C.shape[0]; cp=EE.cpm(C).tocsc(); i0,i1=EE.iu(n); ids=EE.ids_from(gi,i0,i1)
subs,cnt=np.unique(sbj,return_counts=True); keep=[s for s,c in zip(subs,cnt) if c>=20]
ref=np.fromfile(os.path.expanduser('~/xcell-xbulk-revision/eric/edges/')+f'sc_coexmats_{key}__0.01.i32',dtype=np.int32)
def corsparse(Xc):   # Xc sparse cells x genes  (qlcMatrix::corSparse)
    m=Xc.shape[0]; mu=np.asarray(Xc.mean(0)).ravel(); XtX=(Xc.T@Xc).toarray(); cov=(XtX-m*np.outer(mu,mu))/(m-1); sd=np.sqrt(np.diag(cov))
    with np.errstate(invalid='ignore',divide='ignore'): return cov/np.outer(sd,sd)
def run(corr,rnd=None):
    acc=np.zeros(len(i0))
    for s in keep:
        idx=np.where(sbj==s)[0]; Cm=corr(cp[:,idx].T.tocsr()); v=Cm[i0,i1]; 
        if rnd: v=np.round(v,rnd)
        v[np.isnan(v)]=np.nanmedian(v); acc+=EE.rank_min_over_max(v)
    r=EE.rank_avg_over_max(acc/len(keep)); e=np.sort(ids[r>=0.99]); return np.intersect1d(ref,e).size/np.union1d(ref,e).size
def cordense(Xc):
    X=Xc.toarray(); m=X.shape[0]; mu=X.mean(0); cov=(X.T@X-m*np.outer(mu,mu))/(m-1); sd=np.sqrt(np.diag(cov))
    with np.errstate(invalid='ignore',divide='ignore'): return cov/np.outer(sd,sd)
print('corSparse formula, dense BLAS',run(cordense))
