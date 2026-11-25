"""Comparison on Eric Chu's pipeline edges (eric/edges/*.i32, exported by eric/R/e1_export_edges.R).
Replicates his evaluationUtils$fisher exactly: truth = reference positives; predicted = test positives; universe = ALL pairs of the TEST group's own gene set
 (so reference positives involving genes absent from the test group are dropped from truth; test positives not in truth count in b). OR=(a/b)/(c/d), p one-sided Fisher.
A 'shared' restriction (both groups' common genes only) is available for A7."""
import numpy as np, os
from scipy.stats import hypergeom
E=os.path.expanduser('~/xcell-xbulk-revision/eric/edges/'); G=29426
CT=['astrocyte','excitatory','inhibitory','microglia','oligodendrocyte','opc']
DS=['lau','lim','nagy','pineda','ramos','rosmap','velmeshev']
_cache={}
def genes(level,ds,ct):
    k=(level,ds,ct,'g')
    if k not in _cache: _cache[k]=np.fromfile(f'{E}{level}_coexmats_{ds}_{ct}__genes.i32',dtype=np.int32)
    return _cache[k]
def edges(level,ds,ct,frac,pre='coexmats'):
    return np.fromfile(f'{E}{level}_{pre}_{ds}_{ct}__{frac}.i32',dtype=np.int32)
def mask(g): m=np.zeros(G,bool); m[g]=True; return m
def counts(ref,test,gref,gtest,shared=False):
    """ref,test sorted int32 pair ids; gref,gtest gene index arrays of the two groups."""
    ms=mask(gtest); mr=mask(gref)
    if shared:
        m=ms&mr; ref=ref[m[ref//G]&m[ref%G]]; test=test[m[test//G]&m[test%G]]; k=int(m.sum())
        a=int(np.intersect1d(ref,test,assume_unique=True).size); S=test.size; T=ref.size; N=k*(k-1)//2
        b=S-a; c=T-a; d=N-S-c
    else:
        k=int(ms.sum()); N=k*(k-1)//2
        a=int(np.intersect1d(ref,test,assume_unique=True).size); S=test.size
        T_in_S=int((ms[ref//G]&ms[ref%G]).sum()); b=S-a; c=T_in_S-a; d=N-S-c
    OR=(a/b)/(c/d) if min(b,c,d)>0 and a>0 else float('nan')
    lp=float(hypergeom.logsf(a-1,N,a+c,S)) if a>0 else 0.0   # P(X>=a), draws=S test positives, K=a+c truth in universe
    return dict(a=a,b=b,c=c,d=d,N=N,n_ref_pos=int(ref.size),n_test_pos=int(S),odds_ratio=OR,log_p=lp,p_value=float(np.exp(lp)),precision=a/(a+b) if a+b else np.nan,recall=a/(a+c) if a+c else np.nan,F1=2*a/(2*a+b+c),jaccard=a/(a+b+c))
