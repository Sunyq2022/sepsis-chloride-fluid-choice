"""Fixed prior-LR diagnostic: diagnose first, then explicitly estimate 8 endpoints."""
import argparse, hashlib, json
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
import n09_corrected_estimand as mod
import canonical_rerun_n02_n07_n09 as m

OUT=mod.OUT/'review_audit'
COHORT=mod.OUT/'cohort_private.parquet'

def run(mode):
    d=pd.read_parquet(COHORT).sort_values(['subject_id','stay_id']).reset_index(drop=True)
    assert len(d)==552 and d.subject_id.is_unique and d.lr_ml.notna().all()
    d['prior_LR_any']=d.lr_ml.gt(0).astype(int)
    plan={'cohort_sha256':hashlib.sha256(COHORT.read_bytes()).hexdigest(),
          'n':552,'added_feature':'prior_LR_any = I(lr_ml > 0), 24-hour pre-trigger window',
          'numeric':list(mod.NUM)+['prior_LR_any'],'categorical':mod.CAT,
          'models':['Existing L2 plus prior LR indicator','Unpenalized plus prior LR indicator'],
          'endpoints':mod.ENDPOINTS,'bootstrap':500,'seed':mod.PROTOCOL['seed'],
          'trigger':'External review identified weighted binary treatment-history imbalance',
          'interpretation':'Post-result diagnostic; original main estimates retained; both models and all eight outcomes reported regardless of direction',
          'design_gate':'Assess balance before outcomes; do not choose specifications by outcome. Report residual SMD >0.1 and ESS without tuning.'}
    path=OUT/'prior_lr_plan.json'
    if mode=='diagnose':
        if path.exists():assert json.loads(path.read_text())==plan
        else:path.write_text(json.dumps(plan,indent=2),encoding='utf-8')
    else:assert json.loads(path.read_text())==plan
    rows=[];diagnostics=[]
    for label,unpen,added in [('Original main',False,False)]+[(x,i==1,True) for i,x in enumerate(plan['models'])]:
        if mode=='estimate' and not added:continue
        mod.NUM=plan['numeric'] if added else plan['numeric'][:-1]
        dd,prep,x,ps,w=mod.fit(d,m,unpen);a=dd.exposure.to_numpy()
        dense=x.toarray() if hasattr(x,'toarray') else x
        values=[(k,dense[:,j]) for j,k in enumerate(prep.get_feature_names_out())]
        values.append(('clinical_prior_LR_any_raw',dd.prior_LR_any.to_numpy()))
        for k,v in values:
            sd=np.sqrt((v[a==0].var()+v[a==1].var())/2)
            means=[np.average(v[a==arm],weights=w[a==arm]) for arm in [0,1]]
            rows.append(dict(analysis=label,variable=k,weighted_ns=means[0],weighted_lr=means[1],smd_after=abs(means[1]-means[0])/sd if sd else 0))
        diagnostics.append(dict(analysis=label,ess_ns=w[a==0].sum()**2/(w[a==0]**2).sum(),ess_lr=w[a==1].sum()**2/(w[a==1]**2).sum(),ps_min=float(ps.min()),ps_max=float(ps.max()),weight_min=float(w.min()),weight_max=float(w.max()),max_smd=max(r['smd_after'] for r in rows if r['analysis']==label)))
    if mode=='diagnose':
        pd.DataFrame(rows).to_csv(OUT/'prior_lr_balance.csv',index=False)
        (OUT/'prior_lr_design_diagnostics.json').write_text(json.dumps(diagnostics,indent=2),encoding='utf-8')
        print(json.dumps(diagnostics,indent=2),flush=True)
        print(pd.DataFrame(rows).query("variable.str.contains('prior_LR_any')",engine='python').to_string(index=False),flush=True)
        return
    assert (OUT/'prior_lr_design_diagnostics.json').exists()
    mod.NUM=plan['numeric'];mod.PROTOCOL['bootstrap']=500
    results=[];diags=[]
    for i,label in enumerate(plan['models']):
        rr,_,_,diag,_=mod.estimate(d,m,label,unpenalized=i==1)
        results.extend(rr);diags.append(diag)
    assert len(results)==16 and all(x['bootstrap_success']==500 for x in diags)
    pd.DataFrame(results).to_csv(OUT/'prior_lr_outcomes.csv',index=False)
    (OUT/'prior_lr_estimation_checks.json').write_text(json.dumps(diags,indent=2),encoding='utf-8')
    print('Completed both fixed models and all eight endpoints.',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['diagnose','estimate']);args=p.parse_args()
    with threadpool_limits(limits=1):run(args.mode)

