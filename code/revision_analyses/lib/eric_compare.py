"""Replicability / preservation tables for variant edge sets on Eric's definition (lib/eric.counts). usage: eric_compare.py <A7|A8|A4|A6raw> [workers]
Edge spec -> file: base_sc / base_sbj (Eric's saved matrices), x:<variant> (eric_pass.py), s:<variant> (eric_xsbj.py)."""
import sys,os,json,hashlib,datetime,numpy as np,pandas as pd
sys.path.insert(0,os.path.expanduser('~/xcell-xbulk-revision/lib')); import eric
from multiprocessing import Pool
R=os.path.expanduser('~/xcell-xbulk-revision/'); EV=R+'eric/edges_var/'; FR=[0.005,0.01,0.05,0.10]
def path(spec,ds,ct,fr):
    k=f'{ds}_{ct}'
    if spec=='base_sc': return f'{eric.E}sc_coexmats_{k}__{fr}.i32'
    if spec=='base_sbj': return f'{eric.E}sbj_coexmats_{k}__{fr}.i32'
    if spec.startswith('x:'): return f'{EV}{spec[2:]}__{k}__{fr}.i32'
    if spec.startswith('s:'): return f'{EV}sbj_{spec[2:]}__{k}__{fr}.i32'
def load(spec,ds,ct,fr): return np.fromfile(path(spec,ds,ct,fr),dtype=np.int32)
def lvl(spec): return 'sbj' if (spec=='base_sbj' or spec.startswith('s:')) else 'sc'
GR=sorted({(r.split('_')[2],r.split('_')[3]) for r in os.listdir(eric.E) if r.startswith('sc_coexmats_') and r.endswith('__genes.i32')})
def job(t):
    kind,sr,st,ds_r,ds_t,ct,fr,sh=t
    ref=load(sr,ds_r,ct,fr); test=load(st,ds_t,ct,fr)
    r=eric.counts(ref,test,eric.genes(lvl(sr),ds_r,ct),eric.genes(lvl(st),ds_t,ct),shared=sh)
    r.update(analysis=kind,cell_type=ct,ref_spec=sr,test_spec=st,reference_dataset=ds_r,test_dataset=ds_t,frac=fr,universe='shared_genes' if sh else 'test_group_genes'); return r
def specs(which):
    S=[]   # (kind, ref spec, test spec, fracs, shared options)
    if which=='A7':
        for v in ['std','A7']:
            for sh in [False,True]:
                S+=[('replicability',f'x:{v}',f'x:{v}',FR,sh),('preservation',f'x:{v}','base_sbj',FR,sh)]
        S+=[('replicability','base_sc','base_sc',[0.01],True)]
    if which=='A8':
        for v in ['A8_wrank','A8_r','A8_z']: S+=[('replicability',f'x:{v}',f'x:{v}',FR,False),('preservation',f'x:{v}','base_sbj',FR,False)]
        S+=[('replicability','x:std','x:std',FR,False),('preservation','x:std','base_sbj',FR,False)]
    if which=='A4':
        for v in ['std','A4_sum_cellcpm','A4_sum_counts_cpm','A4_log_noqn','A4_nolog_noqn','A4_mean_nolog_qn']:
            S+=[('replicability',f's:{v}',f's:{v}',FR,False),('preservation','base_sc',f's:{v}',FR,False)]
    if which=='A6raw':
        S+=[('replicability','x:rawmean','x:rawmean',FR,False),('preservation','x:rawmean','base_sbj',FR,False),('replicability','x:std','x:std',FR,False),
            ('replicability','x:rawsign','x:rawsign',['pos'],False),('replicability','s:std','s:std',['pos'],False),('preservation','x:rawsign','s:std',['pos'],False)]
    return S
if __name__=='__main__':
    which=sys.argv[1]; jobs=[]; miss=0
    for kind,sr,st,fracs,sh in specs(which):
        for ct in eric.CT:
            ds=[d for d,c in GR if c==ct]
            for r in ds:
                for t in ds:
                    if kind=='replicability' and r==t: continue
                    if kind=='preservation' and r!=t: continue   # same-dataset group preservation (39); cross-dataset omitted here
                    for fr in fracs:
                        if os.path.exists(path(sr,r,ct,fr)) and os.path.exists(path(st,t,ct,fr)): jobs.append((kind,sr,st,r,t,ct,fr,sh))
                        else: miss+=1
    print(which,len(jobs),'jobs; missing',miss,flush=True)
    with Pool(int(sys.argv[2]) if len(sys.argv)>2 else 8) as P: res=P.map(job,jobs,chunksize=4)
    d=pd.DataFrame(res); out={'A7':'A7_no_imputation','A8':'A8_weighted','A4':'A4_xsubject_pseudobulk','A6raw':'A6_raw'}[which]; os.makedirs(R+out,exist_ok=True)
    d.to_csv(R+out+'/results.tsv',sep='\t',index=False)
    g=d.groupby(['analysis','ref_spec','test_spec','universe','frac']).agg(n=('odds_ratio','size'),mean_OR=('odds_ratio','mean'),median_OR=('odds_ratio','median'),min_OR=('odds_ratio','min'),n_OR_le_1=('odds_ratio',lambda x:int((x<=1).sum())),mean_F1=('F1','mean'),mean_jaccard=('jaccard','mean')).reset_index(); g.to_csv(R+out+'/summary.tsv',sep='\t',index=False)
    json.dump(dict(script=os.path.abspath(__file__),script_sha256=hashlib.sha256(open(__file__,'rb').read()).hexdigest(),timestamp=datetime.datetime.now().isoformat(),n_rows=len(d),n_missing_jobs=miss,seed='n/a (deterministic)',source='Eric Chu 041_CH4_FINAL (R) saved matrices + raw counts exported from sc_exprmats_*.rds',engine='lib/eric_engine.py (qlcMatrix corSparse formula)'),open(R+out+'/results_meta.json','w'),indent=1)
    print(g.round(3).to_string())
