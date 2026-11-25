"""Recompute Eric Chu's reproducibility (sc, sbj) and preservation (sc->sbj) tables from his saved coexpression matrices (via exported edge sets) and compare with his saved fisher tables.
Output: replicability.tsv, preservation.tsv (4 thresholds), validation.tsv"""
import sys,os,json,hashlib,datetime,numpy as np,pandas as pd
sys.path.insert(0,os.path.expanduser('~/xcell-xbulk-revision/lib')); import eric
from multiprocessing import Pool
B=os.path.expanduser('~/xcell-xbulk-revision/eric/'); FR=[0.005,0.01,0.05,0.10]
qc=pd.concat([pd.read_csv(B+f'edges/qc_{l}.tsv',sep='\t') for l in ['sc','sbj']]); qc['key']=qc.file.str.replace('.rds','',regex=False)
groups=sorted({(r.file.split('_')[2],r.file.split('_')[3].replace('.rds','')) for r in qc.itertuples() if r.file.startswith('sc_')})
def job(t):
    kind,lr,lt,ds_r,ds_t,ct,fr,sh=t
    ref=eric.edges(lr,ds_r,ct,fr); test=eric.edges(lt,ds_t,ct,fr)
    r=eric.counts(ref,test,eric.genes(lr,ds_r,ct),eric.genes(lt,ds_t,ct),shared=sh)
    r.update(analysis=kind,cell_type=ct,reference_level=lr,reference_dataset=ds_r,test_level=lt,test_dataset=ds_t,frac=fr,universe='shared_genes' if sh else 'test_group_genes',n_genes_ref=len(eric.genes(lr,ds_r,ct)),n_genes_test=len(eric.genes(lt,ds_t,ct))); return r
if __name__=='__main__':
    jobs=[]
    for ct in eric.CT:
        ds=[d for c,d in sorted(groups) if c==ct] if False else [d for d,c in groups if c==ct]
        for fr in FR:
            for lvl in ['sc','sbj']:
                for r in ds:
                    for t in ds:
                        if r!=t: jobs.append(('replicability',lvl,lvl,r,t,ct,fr,False))
            for r in ds: jobs.append(('preservation','sc','sbj',r,r,ct,fr,False))
            for r in ds:
                for t in ds:
                    if r!=t: jobs.append(('preservation_cross','sc','sbj',r,t,ct,fr,False))
    print(len(jobs),'comparisons',flush=True)
    with Pool(16) as P: res=P.map(job,jobs,chunksize=4)
    d=pd.DataFrame(res); d.to_csv(B+'baseline/all_comparisons.tsv',sep='\t',index=False)
    d[d.analysis=='replicability'].to_csv(B+'baseline/replicability.tsv',sep='\t',index=False); d[d.analysis!='replicability'].to_csv(B+'baseline/preservation.tsv',sep='\t',index=False)
    # validate vs saved
    val=[]
    def chk(name,mine,saved,keys_m,keys_s):
        m=mine.merge(saved,left_on=keys_m,right_on=keys_s,suffixes=('','_saved'))
        val.append(dict(table=name,n_saved=len(saved),n_matched=len(m),max_abs_log_OR_diff=float(np.abs(np.log(m.odds_ratio)-np.log(m.or_saved)).max()),max_abs_a_diff=float(np.abs(m.a-m.recovered_n).max()),n_OR_gt1_mine=int((m.odds_ratio>1).sum()),n_OR_gt1_saved=int((m.or_saved>1).sum())))
    S=lambda f:pd.read_csv(B+f'tables/{f}.tsv',sep='\t')
    one=d[d.frac==0.01]
    sc=S('saved_sc_repro_fisher').rename(columns={'or':'or_saved'}); chk('sc_repro',one[(one.analysis=='replicability')&(one.reference_level=='sc')],sc,['cell_type','reference_dataset','test_dataset'],['cell_type','dat_truth','dat_sample'])
    sb=S('saved_sbj_repro_fisher').rename(columns={'or':'or_saved'}); chk('sbj_repro',one[(one.analysis=='replicability')&(one.reference_level=='sbj')],sb,['cell_type','reference_dataset','test_dataset'],['cell_type','dat_truth','dat_sample'])
    pr=S('saved_sc_sbj_preservation_fisher').rename(columns={'or':'or_saved'}); chk('sc_sbj_preservation',one[one.analysis=='preservation'],pr,['cell_type','reference_dataset'],['cell_type','dat_truth'])
    v=pd.DataFrame(val); v.to_csv(B+'baseline/validation_vs_saved.tsv',sep='\t',index=False); print(v.to_string())
    json.dump(dict(script=os.path.abspath(__file__),script_sha256=hashlib.sha256(open(__file__,'rb').read()).hexdigest(),timestamp=datetime.datetime.now().isoformat(),source='/space/scratch/ericchu/r_cache/041_CH4_FINAL/output (read-only); R code /home/echu113/ws_rstudio/_workspaces/041_CH4_FINAL',matrices_sha256=dict(zip(qc.key,qc.sha256)),definition='positive = value>=1-frac, lower triangle; fisher as evaluationUtils$fisher'),open(B+'baseline/baseline_meta.json','w'),indent=1)
