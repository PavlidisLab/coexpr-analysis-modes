"""xBulk baseline on Eric's matrices: bulk replicability (truth = one bulk dataset, predicted = the other) and preservation of sc/sbj edges (truth) in the bulk network (predicted), recomputed and compared with bk_repro_fisher / bk_preservation_fisher."""
import sys,os,numpy as np,pandas as pd,json
sys.path.insert(0,os.path.expanduser('~/xcell-xbulk-revision/lib')); import eric
E=eric.E; B=os.path.expanduser('~/xcell-xbulk-revision/eric/')
bk=sorted({f.split('__')[0].replace('bk_coexmats_','') for f in os.listdir(E) if f.startswith('bk_coexmats_') and f.endswith('__genes.i32')})
def be(n,fr): return np.fromfile(f'{E}bk_coexmats_{n}__{fr}.i32',dtype=np.int32)
def bg(n): return np.fromfile(f'{E}bk_coexmats_{n}__genes.i32',dtype=np.int32)
FR=[0.005,0.01,0.05,0.10]; rows=[]
for fr in FR:
    for t in bk:
        for p in bk:
            if t[:3]==p[:3]: continue
            if t.endswith('pres')!=p.endswith('pres'): continue
            r=eric.counts(be(t,fr),be(p,fr),bg(t),bg(p)); r.update(analysis='bulk_replicability',frac=fr,reference=t,test=p,cell_type='',level='bk'); rows.append(r)
    for p in bk:
        ds=p.split('-')[0]
        for lvl in ['sc','sbj']:
            for ct in eric.CT:
                f=f'{E}{lvl}_coexmats_{ds}_{ct}__{fr}.i32'
                if not os.path.exists(f): continue
                r=eric.counts(np.fromfile(f,dtype=np.int32),be(p,fr),eric.genes(lvl,ds,ct),bg(p)); r.update(analysis='bulk_preservation',frac=fr,reference=f'{lvl}:{ds}',test=p,cell_type=ct,level=lvl); rows.append(r)
d=pd.DataFrame(rows); d.to_csv(B+'baseline/bulk_comparisons.tsv',sep='\t',index=False)
S=pd.read_csv(B+'tables/saved_bk_repro_fisher.tsv',sep='\t'); one=d[(d.frac==0.01)&(d.analysis=='bulk_replicability')]
m=one.merge(S,left_on=['reference','test'],right_on=['dat_truth','dataset']); print('bulk repro matched',len(m),'of',len(S),'max |log OR diff|',float(np.abs(np.log(m.odds_ratio)-np.log(m['or'])).max()))
P=pd.read_csv(B+'tables/saved_bk_preservation_fisher.tsv',sep='\t'); pr=d[(d.frac==0.01)&(d.analysis=='bulk_preservation')]
m=pr.merge(P,left_on=['test','level','cell_type'],right_on=['dataset','level','cell_type'],suffixes=('','_saved')); print('bulk preservation matched',len(m),'of',len(P),'max |log OR diff|',float(np.abs(np.log(m.odds_ratio)-np.log(m['or'])).max()))
print(one[['reference','test','odds_ratio']].to_string()); print(pr.groupby(['test','level']).odds_ratio.agg(['size','mean','min']).round(2).to_string())
