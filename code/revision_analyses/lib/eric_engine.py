"""Python re-implementation of Eric Chu's xCell / xSubject coexpression definitions (041_SETUP.R: computeAggScCoex, computeCoexmat, normalizeCoexmat, aggregSbj, logTrans, quantileNormalize),
run on raw counts exported by eric/R/e3_export_counts.R. CPM = counts / colSums(group genes) * 1e6 (verified in the export). Pair id convention = lib/eric.py (max(g)*G+min(g), global C-sorted gene idx).
Everything is validated against Eric's saved matrices via their exported edge sets (see eric/validation)."""
import numpy as np, scipy.sparse as sp, os, json
from scipy.stats import rankdata
X=os.path.expanduser('~/xcell-xbulk-revision/eric/expr/'); G=29426
def load(ds,ct):
    k=f'{X}{ds}_{ct}'; gi=np.fromfile(k+'__genes.i32',dtype=np.int32)
    i=np.fromfile(k+'__i.i32',dtype=np.int32); p=np.fromfile(k+'__p.i32',dtype=np.int32); x=np.fromfile(k+'__x.f64',dtype=np.float64)
    sbj=np.array(open(k+'__cellsbj.txt').read().split('\n')[:len(p)-1])
    C=sp.csc_matrix((x,i,p),shape=(len(gi),len(p)-1)); o=np.argsort(gi,kind='stable')   # sort genes by global idx (= sorted ensembl ids)
    return C[o,:].tocsc(),gi[o],sbj
def cpm(C): s=np.asarray(C.sum(0)).ravel(); return C@sp.diags(1e6/np.maximum(s,1))
def iu(n): return np.triu_indices(n,1)   # pairs r<c in row-major == R's lower.tri column-major order
def ids_from(gi,iu0,iu1): return (gi[iu1].astype(np.int64)*G+gi[iu0]).astype(np.int32)
def rank_min_over_max(v): r=rankdata(v,method='min'); return r/r.max()
def rank_avg_over_max(v): r=rankdata(v,method='average'); return r/r.max()
def edges(ids,rank_vec,fracs,valid=None):
    out={}
    for fr in fracs:
        m=rank_vec>=1-fr
        if valid is not None: m&=valid
        out[fr]=np.sort(ids[m])
    return out
def quantile_normalize(M):   # columns = samples; rank ties.min, mean of sorted columns, as workspace$utils$quantileNormalize
    S=np.sort(M,axis=0); mean=S.mean(axis=1); R=np.apply_along_axis(lambda c:rankdata(c,method='min'),0,M).astype(int)
    return mean[R-1]

def corr_qlc(X):
    """Pearson between columns as qlcMatrix::corSparse computes it: cov=(X'X - n mu mu')/(n-1); cor=cov/outer(sd,sd) (dense BLAS here). Zero-variance genes -> NaN."""
    m=X.shape[0]; mu=X.mean(0); cov=(X.T@X-m*np.outer(mu,mu))/(m-1); sd=np.sqrt(np.diag(cov))
    with np.errstate(invalid='ignore',divide='ignore'): return cov/np.outer(sd,sd)
