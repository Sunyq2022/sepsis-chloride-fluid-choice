"""Reviewer-requested N09 revision. Original code, protocol and data are read-only.

Run with the existing icu_causal_audit_py313 environment, --mode core or source.
All exports are aggregate. Source-row caches are restricted to work/private.
The six-hour landmark is conditional on post-treatment survival, ICU residence
and freedom from the composite; it is a sensitivity analysis, not a causal fix.
"""
from pathlib import Path
import argparse, hashlib, json, sys
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

sys.dont_write_bytecode = True
REV = Path(__file__).resolve().parents[1]
ROOT = Path(__file__).resolve().parents[5]
OLD = REV.parent / 'rebuilt_analysis'
CODE = ROOT / '02代码/01分析脚本'
sys.path.insert(0, str(CODE))
import n09_corrected_estimand as n
import canonical_rerun_n02_n07_n09 as m

def save(name, value):
    (REV/'results'/name).write_text(json.dumps(value, ensure_ascii=False, indent=2, default=str), encoding='utf-8')

def hashes():
    paths = [OLD/'cohort_private.parquet', OLD/'outcomes.csv', OLD/'protocol.json',
             CODE/'n09_corrected_estimand.py', CODE/'canonical_rerun_n02_n07_n09.py',
             REV.parent/'Critical_Care_Correspondence_submission/Manuscript.docx']
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def orient(df):
    df = df.copy()
    df['rd_ns_lr_pp'] = -100*df.rd
    df['rd_ns_lr_low_pp'] = -100*df.rd_high
    df['rd_ns_lr_high_pp'] = -100*df.rd_low
    df['rr_ns_lr'] = 1/df.rr
    df['rr_ns_lr_low'] = 1/df.rr_high
    df['rr_ns_lr_high'] = 1/df.rr_low
    return df

def baseline(d, w):
    """Means/SDs and proportions among observed values; missingness is separate.
    Binary prior-LR use is an extra diagnostic, not a frozen model feature.
    Standardization denominator is the pooled unweighted population SD.
    """
    d=d.copy();d['prior_lr_any']=(d.lr_ml>0).astype(float)
    rows=[];missing=[];arm=d.exposure.to_numpy(int)
    def row(col,lev,v):
        z={'variable':col,'level':lev}
        av=[];wv=[];vs=[]
        for a,label in [(0,'NS'),(1,'LR')]:
            ix=(arm==a)&np.isfinite(v);x=v[ix];ww=w[ix]
            av.append(float(x.mean()));wv.append(float(np.average(x,weights=ww)));vs.append(float(x.var()))
            sd=float(np.sqrt(np.average((x-wv[-1])**2,weights=ww)))
            z.update({label+'_n_observed':int(ix.sum()),label+'_raw_mean':av[-1],label+'_raw_sd':float(x.std()),label+'_weighted_mean':wv[-1],label+'_weighted_sd':sd})
        denom=np.sqrt(np.mean(vs));z['smd_before']=abs(av[1]-av[0])/denom if denom else 0.;z['smd_after']=abs(wv[1]-wv[0])/denom if denom else 0.;rows.append(z)
    for col in n.NUM+['prior_lr_any']:
        row(col,'',pd.to_numeric(d[col]).to_numpy(float))
    for col in n.CAT:
        for level in sorted(d[col].dropna().unique()):
            v=np.where(d[col].isna(),np.nan,d[col].eq(level).astype(float));row(col,str(level),v)
    for col in n.NUM+n.CAT:
        z={'variable':col,'method':'median + missing indicator' if col in n.NUM else 'mode + one-hot encoding'}
        for a,label in [(0,'NS'),(1,'LR')]:
            ix=arm==a;z[label+'_missing']=int(d.loc[ix,col].isna().sum());z[label+'_n']=int(ix.sum())
        missing.append(z)
    pd.DataFrame(rows).to_csv(REV/'results/baseline_before_after.csv',index=False)
    pd.DataFrame(missing).to_csv(REV/'results/missingness.csv',index=False)

def core(d):
    plan={'date':'2026-10-24','frozen_seed':n.PROTOCOL['seed'],'bootstrap':500,
          'new_analysis':'Six-hour landmark, retain patients alive, in index ICU and composite-event-free at six hours; follow to original day 7 (162 subsequent hours). Refit original baseline PS in landmark sample and in each patient bootstrap.',
          'interpretation':'Post-treatment selection; sensitivity analysis only; neither early-event removal nor timing audit proves exchangeability.',
          'unchanged':['cohort definition','original main endpoint','original main results','original scripts','original seed'],
          'additional_components':['creatinine3','KRT'],'prior_LR_window':'24 hours before chloride trigger; not confirmed LR-naive'}
    save('revision_plan.json',plan)
    dd,prep,x,ps,w=n.fit(d,m,False);a=dd.exposure.to_numpy(int)
    main=pd.read_csv(OLD/'outcomes.csv').query("analysis == 'Main'")
    for z in main.itertuples():
        for arm,key in [(0,'risk_ns'),(1,'risk_lr')]:
            ix=a==arm
            risk=m.weighted_aj(dd.loc[ix,z.outcome+'_time'].to_numpy(float),dd.loc[ix,z.outcome+'_event'].to_numpy(int),w[ix],n.ENDPOINTS[z.outcome])
            assert abs(risk-getattr(z,key))<1e-10, (z.outcome,key,risk)
    baseline(dd,w)
    # Model-encoded balance, including imputation indicators and one-hot terms.
    dense=x.toarray() if hasattr(x,'toarray') else x;rows=[]
    for j,feature in enumerate(prep.get_feature_names_out()):
        v=dense[:,j];sd=np.sqrt((v[a==0].var()+v[a==1].var())/2)
        rows.append({'variable':feature,'smd_before':abs(v[a==1].mean()-v[a==0].mean())/sd if sd else 0,'smd_after':abs(np.average(v[a==1],weights=w[a==1])-np.average(v[a==0],weights=w[a==0]))/sd if sd else 0})
    pd.DataFrame(rows).to_csv(REV/'results/model_balance.csv',index=False)
    hist=[];support=[]
    for arm,label in [(0,'NS'),(1,'LR')]:
        ix=a==arm;pa=ps[ix];wa=w[ix]
        support.append({'arm':label,'n':int(ix.sum()),'ess':wa.sum()**2/(wa**2).sum(),'ps_min':pa.min(),'ps_p01':np.quantile(pa,.01),'ps_p50':np.median(pa),'ps_p99':np.quantile(pa,.99),'ps_max':pa.max(),'weight_min':wa.min(),'weight_max':wa.max(),'ps_under_005':int((pa<.05).sum()),'ps_over_095':int((pa>.95).sum())})
        raw,edges=np.histogram(pa,bins=np.linspace(0,1,21));weighted,_=np.histogram(pa,bins=edges,weights=wa)
        for j in range(20):hist.append({'arm':label,'left':edges[j],'right':edges[j+1],'count':int(raw[j]),'proportion':raw[j]/ix.sum(),'weighted_proportion':weighted[j]/wa.sum()})
    pd.DataFrame(support).to_csv(REV/'results/propensity_support.csv',index=False)
    pd.DataFrame(hist).to_csv(REV/'results/propensity_distribution.csv',index=False)
    orient(main).to_csv(REV/'results/table1.csv',index=False)
    frozen=[]
    for path in [OLD/'outcomes.csv',OLD/'review_audit/review_sensitivity.csv',OLD/'review_audit/prior_lr_outcomes.csv',OLD/'review_audit/round2_adjusted_outcomes.csv']:
        z=pd.read_csv(path);z['source_file']=str(path.relative_to(REV.parent));frozen.append(z)
    orient(pd.concat(frozen,ignore_index=True)).to_csv(REV/'results/frozen_estimates_reoriented.csv',index=False)
    early=d.primary_event.eq(1)&d.primary_time.le(6)
    landmark=d[d.primary_time.gt(6)].copy()
    counts=[]
    for arm,label in [(0,'NS'),(1,'LR')]:
        ix=d.exposure.eq(arm);ee=early&ix
        krt6=d.KRT_event.eq(1)&d.KRT_time.le(6)
        counts.append({'arm':label,'original_n':int(ix.sum()),'early_composite':int(ee.sum()),'early_KRT':int((krt6&ix).sum()),'KRT_first_component_early':int((ee&d.KRT_stamp.le(d.creatinine3_stamp.fillna(pd.Timestamp.max))).sum()),'early_competing':int((ix&d.primary_event.eq(2)&d.primary_time.le(6)).sum()),'landmark_n':int(landmark.exposure.eq(arm).sum())})
    pd.DataFrame(counts).to_csv(REV/'results/landmark_flow.csv',index=False)
    for ep in ['primary','creatinine3','KRT']:
        landmark[ep+'_time']-=6
        assert landmark[ep+'_time'].gt(0).all()
    rr,cc,bb,diag,hh=n.estimate(landmark,m,'Six-hour landmark',endpoints={'primary':162.,'creatinine3':162.,'KRT':162.})
    assert diag['bootstrap_success']==500 and not diag['bootstrap_errors']
    orient(pd.DataFrame(rr)).to_csv(REV/'results/landmark_estimates.csv',index=False)
    pd.DataFrame(bb).to_csv(REV/'results/landmark_balance.csv',index=False)
    save('landmark_diagnostics.json',diag)
    save('core_checks.json',{'all_main_risks_reproduced':True,'main_n':len(d),'source_hashes':hashes(),'landmark_bootstrap_success':diag['bootstrap_success'],'model_columns':n.NUM+n.CAT,'main_max_modeled_smd':max(z['smd_after'] for z in rows)})

def source(d):
    import duckdb
    private=REV/'work/private';private.mkdir(exist_ok=True)
    con=duckdb.connect();con.execute("set threads=2");con.execute("set memory_limit='4GB'")
    con.register('cohort',d)
    raw=Path('E:/ICUdatabase/MIMIC/mimic-iv')
    queries=[]
    labfile=private/'trigger_labs.parquet'
    q=f"SELECT l.labevent_id,l.itemid,l.charttime,l.storetime,l.valuenum FROM read_csv_auto('{(raw/'hosp/labevents.csv.gz').as_posix()}') l SEMI JOIN cohort c USING(labevent_id)"
    queries.append(q)
    if not labfile.exists():
        print('Reading trigger records from source labevents',flush=True);con.execute(q).df().to_parquet(labfile,index=False)
    labs=pd.read_parquet(labfile);z=d.merge(labs,on='labevent_id',validate='one_to_one',suffixes=('','_source'))
    assert len(z)==len(d) and z.itemid_source.eq(50902).all()
    available=z[['charttime','storetime']].max(axis=1)
    audit={'matched_trigger_records':len(z),'missing_storetime':int(z.storetime.isna().sum()),'trigger_mismatch':int(z.trigger_time.ne(available).sum()),'bolus_before_availability':int(z.t0.lt(available).sum()),'sepsis_derived_time_after_trigger':int(z.sepsis_known.gt(z.trigger_time).sum()),'chloride_value_mismatch':int(z.chloride.ne(z.valuenum).sum()),'sample_to_store_minutes_median':float((z.storetime-z.charttime).dt.total_seconds().median()/60),'trigger_to_bolus_minutes_median':float((z.t0-z.trigger_time).dt.total_seconds().median()/60),'status':'source-row audit; storage is a database availability proxy, not proof of clinician viewing'}
    save('availability_audit.json',audit)
    assert all(audit[k]==0 for k in ['missing_storetime','trigger_mismatch','bolus_before_availability','sepsis_derived_time_after_trigger','chloride_value_mismatch'])
    inpfile=private/'unspecified_prior_fluids.parquet'
    q=f"""SELECT c.stay_id,c.exposure,c.lr_ml,c.trigger_time,c.t0,i.itemid,i.starttime,i.endtime,i.amount,i.amountuom,i.ordercategoryname
    FROM read_csv_auto('{(raw/'icu/inputevents.csv.gz').as_posix()}') i JOIN cohort c USING(stay_id)
    WHERE i.itemid IN (226361,226364,226375) AND i.starttime<c.t0
    AND i.statusdescription IS DISTINCT FROM 'Rewritten'"""
    queries.append(q)
    if not inpfile.exists():
        print('Reading unspecified pre-ICU / OR / PACU fluid records',flush=True);con.execute(q).df().to_parquet(inpfile,index=False)
    inputs=pd.read_parquet(inpfile);rows=[]
    for label,sub in [('All pre-t0',inputs),('Recorded in 24h before trigger',inputs[(inputs.endtime<=inputs.trigger_time)&(inputs.endtime>inputs.trigger_time-pd.Timedelta(hours=24))])]:
        for arm,armname in [(0,'NS'),(1,'LR')]:
            x=sub[sub.exposure.eq(arm)]
            rows.append({'window':label,'arm':armname,'records':len(x),'patients':x.stay_id.nunique(),'patients_no_identified_prior_LR':x.loc[x.lr_ml.eq(0),'stay_id'].nunique()})
    pd.DataFrame(rows).to_csv(REV/'results/unspecified_fluids.csv',index=False)
    aug=pd.read_parquet(OLD/'review_audit/round2_private/augmented_cohort.parquet')
    rows=[]
    for arm,name in [(0,'NS'),(1,'LR')]:
        x=aug[aug.exposure.eq(arm)]
        rows.append({'arm':name,'no_LR_pretrigger':int(x.lr_ml.eq(0).sum()),'of_these_gap_LR_positive':int((x.lr_ml.eq(0)&x.gap_LR_ml.gt(0)).sum()),'pretrigger_prior_LR':int(x.lr_ml.gt(0).sum()),'any_gap_LR':int(x.gap_LR_ml.gt(0).sum())})
    pd.DataFrame(rows).to_csv(REV/'results/prior_LR_window_audit.csv',index=False)
    (REV/'code/revision_source_audit.sql').write_text(';\n\n'.join(queries)+';\n',encoding='utf-8')
    con.close()

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--mode',choices=['core','source'],required=True);args=parser.parse_args()
    before=hashes()
    d=pd.read_parquet(OLD/'cohort_private.parquet').sort_values(['subject_id','stay_id']).reset_index(drop=True)
    assert len(d)==552 and d.subject_id.is_unique
    with threadpool_limits(limits=1):
        (core if args.mode=='core' else source)(d)
    assert before==hashes(),'An immutable source changed during execution'
    save(args.mode+'_source_integrity.json',{'unchanged':True,'sha256':before})
    print(args.mode.upper(),'COMPLETE; original inputs and scripts unchanged',flush=True)
