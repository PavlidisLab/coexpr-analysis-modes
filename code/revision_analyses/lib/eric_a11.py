"""A11 (XCB request 2026-10-19): xBulk like-for-like gene universe. Separates gene universe from level.
Part 1 'restrict': baseline xCell/xSubject/xBulk networks (ranks over ALL genes, edges as originally defined at top 0.5/1/5/10%) evaluated only on pairs among
        the A10 common set (rule p50/p70/p100) of the cell type; OR pair universe = (ref group genes ∩ test group genes ∩ common set). rule 'full' = baseline definition (sanity).
Part 2 'recomp': xBulk recomputed on the common set: log2(CPM+1) (library size over all cleaned bulk genes) -> subset -> quantile normalisation over the SET -> Pearson -> rank within set.
        Compared with A10's recomputed xCell/xSubject on the same set (A10 sc/sbj edges), and with A10's xBulk (bulk NOT re-normalised) for reference.
Part 3 'random': size-matched control: for microglia & excitatory, 20 random gene sets of the p70 common-set size drawn (seeded) from genes in both bulk datasets that are present in >=1 sc dataset of the cell type;
        xBulk with (a) baseline ranks restricted to the set, (b) recomputed+renormalised; replicability ROSMAP<->Velmeshev and preservation (baseline sc/sbj edges restricted to the set, shared universe).
stages: p1 <workers> | p2 | p2cmp <workers> | p3 <workers> | summary"""
import sys,os,math,numpy as np,pandas as pd
sys.path.insert(0,os.path.expanduser('~/xcell-xbulk-revision/lib')); import eric, eric_engine as EE, eric_subset as SS
R=os.path.expanduser('~/xcell-xbulk-revision/'); D=R+'A11_xbulk_like_for_like/'; G=eric.G; FR=[0.005,0.01,0.05,0.10]; FRS=['0.005','0.01','0.05','0.1']
RULES={'p50':0.5,'p70':0.7,'p100':1.0}; SEED=20261019; NDRAW=20; P3CT=['microglia','excitatory']
DSALL=eric.DS
import glob
GR=sorted({tuple(os.path.basename(f).split('__')[0].split('_',1)) for f in glob.glob(EE.X+'*__p.i32')})
def bulk_genes(nm): return np.fromfile(f'{EE.X}bulk_{nm}__genes.i32',dtype=np.int32)
def sets():
    both=np.intersect1d(bulk_genes('rosmap'),bulk_genes('velmeshev')); out={}
    for ct in eric.CT:
        ds=[d for d,c in GR if c==ct]; cnt=np.zeros(G,int)
        for d in ds: cnt[eric.genes('sc',d,ct)]+=1
        for rn,thr in RULES.items():
            k=int(math.ceil(thr*len(ds)-1e-9)); out[(ct,rn)]=np.intersect1d(np.where(cnt>=k)[0],both)
    return out
def restr(ids,m): return ids[m[ids//G]&m[ids%G]]
def fr_str(fr): return {0.005:'0.005',0.01:'0.01',0.05:'0.05',0.1:'0.1',0.10:'0.1'}[fr]
def row(x,**kw): x.update(kw); return x
def cnt(ref,test,gr,gt,shared): return eric.counts(ref,test,gr,gt,shared=shared)
def groups_of(ct): return [d for d,c in GR if c==ct]
# ---------------- part 1
def p1_task(t):
    rn,ct,gs=t; ds=groups_of(ct); full=(rn=='full'); m=None if full else eric.mask(gs); rows=[]
    gen={}
    for lv in ('sc','sbj'):
        for d in ds: gen[(lv,d)]=eric.genes(lv,d,ct) if full else np.intersect1d(eric.genes(lv,d,ct),gs)
    for nm in ('rosmap','velmeshev'): gen[('bk',nm)]=bulk_genes(nm) if full else gs   # baseline bulk net: one per dataset (genes of that bulk)
    def E(lv,d,fr):
        e=eric.edges(lv,d,ct,fr_str(fr)) if lv!='bk' else np.fromfile(f'{eric.E}bk_coexmats_{d}__{fr_str(fr)}.i32',dtype=np.int32)
        return e if full else restr(e,m)
    for fr in FR:
        ed={(lv,d):E(lv,d,fr) for lv in ('sc','sbj') for d in ds}; ed.update({('bk',nm):E('bk',nm,fr) for nm in ('rosmap','velmeshev')})
        def C(a,b,kind): return row(cnt(ed[a],ed[b],gen[a],gen[b],not full),rule=rn,frac=fr,analysis=kind,cell_type=ct,ref_level=a[0],test_level=b[0],reference=a[1],test=b[1],part='restrict',variant='baseline_ranks',n_genes_ref=len(gen[a]),n_genes_test=len(gen[b]))
        for r in ds:
            for t_ in ds:
                if r!=t_: rows+=[C(('sc',r),('sc',t_),'replicability'),C(('sbj',r),('sbj',t_),'replicability')]
            rows.append(C(('sc',r),('sbj',r),'preservation_xsubject'))
            if r in ('rosmap','velmeshev'): rows+=[C(('sc',r),('bk',r),'preservation_xbulk_from_sc'),C(('sbj',r),('bk',r),'preservation_xbulk_from_sbj')]
        rows+=[C(('bk','rosmap'),('bk','velmeshev'),'replicability_xbulk'),C(('bk','velmeshev'),('bk','rosmap'),'replicability_xbulk')]
    return rows
# ---------------- part 2
def load_l2(nm):
    gi=np.fromfile(f'{D}inputs/bulk_{nm}__log2cpm_genes.i32',dtype=np.int32); d=[int(x) for x in open(f'{D}inputs/bulk_{nm}__log2cpm_dim.txt').read().split()]
    return gi,np.fromfile(f'{D}inputs/bulk_{nm}__log2cpm.f64').reshape(d[1],d[0]).T
def bulk_net(L,gi,keep,renorm=True):
    sel=np.where(np.isin(gi,keep))[0]; M=L[sel,:]
    if renorm: M=EE.quantile_normalize(M)
    n=len(sel); i0,i1=EE.iu(n); ids=EE.ids_from(gi[sel],i0,i1)
    with np.errstate(invalid='ignore',divide='ignore'): Cm=np.corrcoef(M)
    v=Cm[i0,i1]; dd=np.isnan(v)
    if dd.any(): v[dd]=np.nanmedian(v)
    return gi[sel],ids,EE.rank_avg_over_max(v)
def p2_task(t):
    rn,ct,nm,gs=t; gi,L=load_l2(nm); g,ids,r=bulk_net(L,gi,gs,True); SS.write(f'A11_{rn}','bk',f'{nm}_{ct}',g,ids,r); return dict(rule=rn,cell_type=ct,bulk=nm,n_genes=len(g))
def pth(tag,lv,key,fr): return f'{SS.S}{tag}/{lv}_{key}__{fr}.i32'
def p2cmp_task(t):
    rn,fr,ct=t; fs=fr_str(fr); rows=[]
    def ld(tag,lv,key): return np.fromfile(pth(tag,lv,key,fs),dtype=np.int32),np.fromfile(pth(tag,lv,key,'').replace('__.i32','__genes.i32'),dtype=np.int32)
    for variant,btag in (('recomputed_renormalised',f'A11_{rn}'),('A10_not_renormalised',f'A10_{rn}')):
        for r,te in (('rosmap','velmeshev'),('velmeshev','rosmap')):
            a,ga=ld(btag,'bk',f'{r}_{ct}'); b,gb=ld(btag,'bk',f'{te}_{ct}'); rows.append(row(eric.counts(a,b,ga,gb),rule=rn,frac=fr,analysis='replicability_xbulk',cell_type=ct,ref_level='bk',test_level='bk',reference=r,test=te,part='recomp',variant=variant,n_genes_ref=len(ga),n_genes_test=len(gb)))
        for nm in ('rosmap','velmeshev'):
            b,gb=ld(btag,'bk',f'{nm}_{ct}')
            for lv in ('sc','sbj'):
                a,ga=ld(f'A10_{rn}',lv,f'{nm}_{ct}'); rows.append(row(eric.counts(a,b,ga,gb),rule=rn,frac=fr,analysis=f'preservation_xbulk_from_{lv}',cell_type=ct,ref_level=lv,test_level='bk',reference=f'{nm}_{ct}',test=f'{nm}_{ct}',part='recomp',variant=variant,n_genes_ref=len(ga),n_genes_test=len(gb)))
    return rows
# ---------------- part 3
_bk={}
def bk_base(nm,fs):
    k=(nm,fs)
    if k not in _bk: _bk[k]=np.fromfile(f'{eric.E}bk_coexmats_{nm}__{fs}.i32',dtype=np.int32)
    return _bk[k]
def p3_task(t):
    ct,draw,size,pool=t; ss=np.random.SeedSequence([SEED,eric.CT.index(ct),draw]); rs=np.random.default_rng(ss)
    rset=np.sort(rs.choice(pool,size=size,replace=False)); m=eric.mask(rset); rows=[]
    nets={}
    L={nm:load_l2(nm) for nm in ('rosmap','velmeshev')}
    for nm in L: nets[('recomputed_renormalised',nm)]=bulk_net(L[nm][1],L[nm][0],rset,True)
    truth={}
    for fr in FR:
        fs=fr_str(fr); bl={nm:restr(bk_base(nm,fs),m) for nm in L}
        vb={'baseline_ranks_restricted':{nm:(bl[nm],rset) for nm in L}}
        vb['recomputed_renormalised']={nm:(np.sort(nets[('recomputed_renormalised',nm)][1][nets[('recomputed_renormalised',nm)][2]>=1-fr]),rset) for nm in L}
        for var,nn in vb.items():
            for r,te in (('rosmap','velmeshev'),('velmeshev','rosmap')):
                rows.append(row(eric.counts(nn[r][0],nn[te][0],rset,rset,shared=True),part='random',variant=var,analysis='replicability_xbulk',frac=fr,cell_type=ct,draw=draw,n_genes=len(rset),reference=r,test=te))
            for nm in L:
                for lv in ('sc','sbj'):
                    k=(lv,nm,fs)
                    if k not in truth: truth[k]=restr(eric.edges(lv,nm,ct,fs),m)
                    g=np.intersect1d(eric.genes(lv,nm,ct),rset)
                    rows.append(row(eric.counts(truth[k],nn[nm][0],g,rset,shared=True),part='random',variant=var,analysis=f'preservation_xbulk_from_{lv}',frac=fr,cell_type=ct,draw=draw,n_genes=len(rset),reference=f'{nm}_{ct}',test=nm))
    return rows
def summ(d,keys):
    return d.groupby(keys).agg(n=('odds_ratio','size'),mean_OR=('odds_ratio','mean'),median_OR=('odds_ratio','median'),min_OR=('odds_ratio','min'),n_le_1=('odds_ratio',lambda x:int((x<=1).sum())),mean_F1=('F1','mean')).reset_index()
if __name__=='__main__':
    from multiprocessing import Pool
    st=sys.argv[1]; S_=sets()
    if st=='p1':
        jobs=[('full',ct,None) for ct in eric.CT]+[(rn,ct,S_[(ct,rn)]) for rn in RULES for ct in eric.CT]
        with Pool(int(sys.argv[2])) as P: res=P.map(p1_task,jobs,chunksize=1)
        pd.DataFrame([r for x in res for r in x]).to_csv(D+'results_part1_restrict.tsv',sep='\t',index=False)
    if st=='p2':
        jobs=[(rn,ct,nm,S_[(ct,rn)]) for rn in RULES for ct in eric.CT for nm in ('rosmap','velmeshev')]
        with Pool(6) as P: res=P.map(p2_task,jobs); pd.DataFrame(res).to_csv(D+'part2_bulk_gene_counts.tsv',sep='\t',index=False)
    if st=='p2cmp':
        jobs=[(rn,fr,ct) for rn in RULES for fr in FR for ct in eric.CT]
        with Pool(int(sys.argv[2])) as P: res=P.map(p2cmp_task,jobs)
        pd.DataFrame([r for x in res for r in x]).to_csv(D+'results_part2_recomp.tsv',sep='\t',index=False)
    if st=='p3':
        both=np.intersect1d(bulk_genes('rosmap'),bulk_genes('velmeshev')); jobs=[]
        for ct in P3CT:
            u=np.unique(np.concatenate([eric.genes('sc',d,ct) for d in groups_of(ct)])); pool=np.intersect1d(u,both); size=len(S_[(ct,'p70')])
            print(ct,'pool',len(pool),'size',size,flush=True)
            jobs+=[(ct,dr,size,pool) for dr in range(NDRAW)]
        with Pool(int(sys.argv[2])) as P: res=P.map(p3_task,jobs,chunksize=1)
        pd.DataFrame([r for x in res for r in x]).to_csv(D+'results_part3_random.tsv',sep='\t',index=False)
