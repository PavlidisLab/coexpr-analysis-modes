"""A1 evaluation: for each scheme (S1,S2), draw d and threshold: (i) replicability of the xCell draw (ordered dataset pairs within cell type; reference draw d vs test draw d of different datasets),
(ii) preservation of the xCell draw (truth) in Eric's xSubject network of the same group (predicted), (iii) preservation in Eric's xBulk networks (rosmap, velmeshev; predicted).
Also writes the full-data baseline rows. usage: eric_a1_compare.py <S1|S2> <workers>"""
import sys,os,json,hashlib,datetime,numpy as np,pandas as pd
sys.path.insert(0,os.path.expanduser('~/xcell-xbulk-revision/lib')); import eric
from multiprocessing import Pool
R=os.path.expanduser('~/xcell-xbulk-revision/'); A=R+'eric/a1/'; FR=[0.005,0.01,0.05,0.10]
GR=sorted({(r.split('_')[2],r.split('_')[3]) for r in os.listdir(eric.E) if r.startswith('sc_coexmats_') and r.endswith('__genes.i32')})
def draw_edges(scheme,ds,ct,d,fr):
    if scheme=='base': return np.fromfile(f'{eric.E}sc_coexmats_{ds}_{ct}__{fr}.i32',dtype=np.int32)
    return np.fromfile(f'{A}{scheme}/{ds}_{ct}__d{d:02d}__{fr}.i32',dtype=np.int32)
BK={'rosmap':'rosmap','velmeshev':'velmeshev'}
def job(t):
    scheme,d,fr=t; rows=[]; E={g:draw_edges(scheme,g[0],g[1],d,fr) for g in GR}
    for ct in eric.CT:
        ds=[a for a,c in GR if c==ct]
        for r in ds:
            for s in ds:
                if r==s: continue
                x=eric.counts(E[(r,ct)],E[(s,ct)],eric.genes('sc',r,ct),eric.genes('sc',s,ct)); x.update(analysis='replicability',scheme=scheme,draw=d,frac=fr,cell_type=ct,reference_dataset=r,test_dataset=s,test_level='sc'); rows.append(x)
        for r in ds:
            sb=np.fromfile(f'{eric.E}sbj_coexmats_{r}_{ct}__{fr}.i32',dtype=np.int32)
            x=eric.counts(E[(r,ct)],sb,eric.genes('sc',r,ct),eric.genes('sbj',r,ct)); x.update(analysis='preservation_xsubject',scheme=scheme,draw=d,frac=fr,cell_type=ct,reference_dataset=r,test_dataset=r,test_level='sbj'); rows.append(x)
            if r in BK:
                bk=np.fromfile(f'{eric.E}bk_coexmats_{BK[r]}__{fr}.i32',dtype=np.int32); bgn=np.fromfile(f'{eric.E}bk_coexmats_{BK[r]}__genes.i32',dtype=np.int32)
                x=eric.counts(E[(r,ct)],bk,eric.genes('sc',r,ct),bgn); x.update(analysis='preservation_xbulk',scheme=scheme,draw=d,frac=fr,cell_type=ct,reference_dataset=r,test_dataset=BK[r],test_level='bk'); rows.append(x)
    return rows
if __name__=='__main__':
    scheme=sys.argv[1]; nd=int(sys.argv[3]) if len(sys.argv)>3 else 50
    tasks=[('base',0,fr) for fr in FR]+[(scheme,d,fr) for d in range(nd) for fr in FR
        if all(os.path.exists(f'{A}{scheme}/{a}_{c}__d{d:02d}__{fr}.i32') for a,c in GR)]
    print(scheme,len(tasks),'tasks',flush=True)
    with Pool(int(sys.argv[2])) as P: res=P.map(job,tasks,chunksize=1)
    d=pd.DataFrame([r for rr in res for r in rr]); out=R+'A1_matched_sample/'; os.makedirs(out,exist_ok=True); d.to_csv(out+f'results_{scheme}.tsv',sep='\t',index=False)
    g=d.groupby(['scheme','draw','analysis','frac']).agg(mean_OR=('odds_ratio','mean'),median_OR=('odds_ratio','median'),min_OR=('odds_ratio','min'),n=('odds_ratio','size'),n_le_1=('odds_ratio',lambda x:int((x<=1).sum())),mean_F1=('F1','mean')).reset_index()
    g.to_csv(out+f'per_draw_summary_{scheme}.tsv',sep='\t',index=False)
    s=g[g.scheme!='base'].groupby(['analysis','frac']).agg(draws=('draw','nunique'),mean_of_mean_OR=('mean_OR','mean'),sd_of_mean_OR=('mean_OR','std'),min_of_mean_OR=('mean_OR','min'),max_of_mean_OR=('mean_OR','max'),mean_n_le_1=('n_le_1','mean')).reset_index()
    b=g[g.scheme=='base'][['analysis','frac','mean_OR','n_le_1']].rename(columns={'mean_OR':'baseline_mean_OR','n_le_1':'baseline_n_le_1'}); s=s.merge(b,on=['analysis','frac']); s.to_csv(out+f'summary_{scheme}.tsv',sep='\t',index=False); print(s.round(3).to_string())
    json.dump(dict(script=os.path.abspath(__file__),script_sha256=hashlib.sha256(open(__file__,'rb').read()).hexdigest(),timestamp=datetime.datetime.now().isoformat(),seed=20261015,draws=nd,scheme=scheme,lab_code='Eric Chu R (reimplemented in lib/eric_engine.py)'),open(out+f'results_{scheme}_meta.json','w'),indent=1)
