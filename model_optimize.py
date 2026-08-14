import ast, re, numpy as np, pandas as pd, cobra
from cobra.sampling import sample
from scipy.stats import spearmanr, pearsonr

MODEL="C:/Users/jccas/OneDrive/Desktop/work_sess/pa_fba/iPAE1146.xml"; MEDIA="C:/Users/jccas/OneDrive/Desktop/work_sess/pa_fba/scfm2.tsv"; PROT="prot_normalized.tsv" #temp
ALL=['AS1','ASA1','AT4','SAS1','SASA1','TS4','TT4']
CCM={'Citrate cycle (TCA cycle)','Glycolysis / Gluconeogenesis','Oxidative phosphorylation',
     'Pentose phosphate pathway','Pyruvate metabolism','Pentose and glucuronate interconversions'}
ESS=['cpd00001','cpd00009','cpd00048','cpd00067','cpd00007','cpd00011','cpd00013',
 'cpd00205','cpd00254','cpd00971','cpd00030','cpd00034','cpd00058','cpd00063','cpd00099',
 'cpd10515','cpd10516','cpd00149','cpd00244','cpd00021','cpd11574','cpd00028','cpd00531']

def find_ex(m,c):
    if c in m.reactions: return m.reactions.get_by_id(c)
    cs=[r for r in m.exchanges if re.search(re.escape(c),r.id)]
    if not cs: return None
    ext=[r for r in cs if re.search(r"_e\d?\b|\[e|_e_",r.id)]
    return (ext or cs)[0]

def media(m):
    md=pd.read_csv(MEDIA,sep='\t'); cc='compounds' if 'compounds' in md.columns else 'compound'
    for ex in m.exchanges: ex.lower_bound=0.0
    for _,r in md.iterrows():
        rx=find_ex(m,str(r[cc]))
        if rx is not None: rx.lower_bound=float(r['minflux']); rx.upper_bound=float(r['maxflux'])
    for c in ESS:
        rx=find_ex(m,c)
        if rx is not None and rx.lower_bound==0.0: rx.lower_bound=-1000.0

def or_branches(rule):
    tree=ast.parse(rule,mode='eval')
    def rec(n):
        if isinstance(n,ast.Expression): return rec(n.body)
        if isinstance(n,ast.Name): return [frozenset([n.id])]
        if isinstance(n,ast.BoolOp):
            if isinstance(n.op,ast.Or):
                o=[]; [o.extend(rec(v)) for v in n.values]; return o
            combos=[frozenset()]
            for v in n.values:
                sub=rec(v); combos=[c|s for c in combos for s in sub]
            return combos
        return [frozenset()]
    return rec(tree)

def map_flux_to_genes(model, act, restrict=None):
    """restrict: set of reaction ids to include (None=all)."""
    gf={}
    for rxn in model.reactions:
        if restrict is not None and rxn.id not in restrict: continue
        a=act.get(rxn.id,0.0)
        if a==0.0 or not rxn.gene_reaction_rule.strip(): continue
        try: br=or_branches(rxn.gene_reaction_rule)
        except Exception: br=[frozenset(g.id for g in rxn.genes)]
        if not br: continue
        sh=a/len(br)
        for b in br:
            for g in b: gf[g]=gf.get(g,0.0)+sh
    return gf

print("[load]"); model,_=cobra.io.sbml.validate_sbml_model(MODEL)
media(model); gmax=model.slim_optimize(); print(f"[media] growth {gmax:.2f}")
ccm_rxns={r.id for r in model.reactions if r.subsystem in CCM}
print(f"[ccm] {len(ccm_rxns)} central-carbon reactions")

bio=[r for r in model.reactions if r.objective_coefficient!=0]
if bio: bio[0].lower_bound=0.3*gmax   # modest growth floor
print("[sampling] n=4000 thinning=60")
S=sample(model,n=4000,thinning=60,processes=1)
act=S.abs().mean(axis=0).to_dict()

prot=pd.read_csv(PROT,sep='\t').set_index('Protein')

def report(gf_dict, cond, label):
    gf=pd.Series(gf_dict)
    d=prot[[cond]].dropna().copy(); d['flux']=d.index.map(gf).fillna(0.0)
    nz=d[d.flux>1e-9]                       # drop blocked/zero-flux genes
    x=np.log10(nz[cond]); y=np.log10(nz['flux'])
    rho,p=spearmanr(x,y) if len(nz)>3 else (np.nan,np.nan)
    pear,_=pearsonr(x,y) if len(nz)>3 else (np.nan,np.nan)
    print(f"  {label:30s} n={len(nz):4d}  Spearman={rho:+.3f} (p={p:.2g})  Pearson(log)={pear:+.3f}")
    return {'label':label,'cond':cond,'n':len(nz),'spearman':round(float(rho),4),
            'p':float(p),'pearson_log':round(float(pear),4)}, nz.reset_index()[['Protein',cond,'flux']]

gf_all=map_flux_to_genes(model,act)
gf_ccm=map_flux_to_genes(model,act,restrict=ccm_rxns)

rows=[]; scat={}
print("\n=== SAS1 (base growth) ===")
r,s=report(gf_all,'SAS1','SAS1 genome-wide'); rows.append(r); scat['SAS1_all']=s
r,s=report(gf_ccm,'SAS1','SAS1 central-carbon'); rows.append(r); scat['SAS1_ccm']=s

print("\n=== all conditions, genome-wide ===")
for c in ALL:
    r,_=report(gf_all,c,f'{c} genome-wide'); rows.append(r)
print("\n=== all conditions, central-carbon ===")
for c in ALL:
    r,_=report(gf_ccm,c,f'{c} central-carbon'); rows.append(r)

pd.DataFrame(rows).to_csv('normalized_correlation.tsv',sep='\t',index=False)
for k,v in scat.items(): v.to_csv(f'scatter_{k}.tsv',sep='\t',index=False)
print("\nsaved normalized_correlation.tsv")
