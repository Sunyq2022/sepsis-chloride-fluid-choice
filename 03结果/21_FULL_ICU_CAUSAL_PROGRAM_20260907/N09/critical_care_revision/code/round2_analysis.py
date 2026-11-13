"""Second-review corrections; historical analyses remain read-only. Aggregate exports only."""
from pathlib import Path
import json, hashlib
import pandas as pd
import numpy as np
import duckdb
from threadpoolctl import threadpool_limits
import revision_analysis as a
n,m=a.n,a.m
REV=a.REV
OUT=REV/'results/round2'
OUT.mkdir(exist_ok=True)
PRIVATE=REV/'work/private/round2'
PRIVATE.mkdir(parents=True,exist_ok=True)
def dump(name,obj):
    (OUT/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,default=str),encoding='utf-8')
def run():
    def source_hashes():
        paths=[a.OLD/'cohort_private.parquet',a.OLD/'outcomes.csv',a.OLD/'protocol.json',
               a.CODE/'n09_corrected_estimand.py',a.CODE/'canonical_rerun_n02_n07_n09.py']
        return {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    before=source_hashes()
    dump('plan.json',{'date':'2026-11-13','seed':n.PROTOCOL['seed'],'bootstrap':n.PROTOCOL['bootstrap'],
         'changes':['Version-specific prior-admission ESRD/dialysis codes','Exclude stage 3 in pre-bolus specimens returned later'],
         'requested_strata':['Absolute only','Rise only','Both'],
         'unchanged':['exposure definition','endpoints','baseline propensity features','bootstrap seed and count'],
         'index_diagnoses':'Audit only: discharge diagnoses are not proof of baseline maintenance dialysis.'})
    d=pd.read_parquet(a.OLD/'cohort_private.parquet').sort_values(['subject_id','stay_id']).reset_index(drop=True)
    raw=Path('E:/ICUdatabase/MIMIC/mimic-iv/hosp')
    con=duckdb.connect();con.execute("SET threads=2");con.register('cohort',d)
    codes="((dx.icd_version=9 AND replace(upper(dx.icd_code),'.','') IN ('5856','V451','V4511')) OR (dx.icd_version=10 AND replace(upper(dx.icd_code),'.','') IN ('N186','Z992')))"
    sql=f"""SELECT c.stay_id, dx.hadm_id, dx.icd_code, dx.icd_version,
        dx.hadm_id=c.hadm_id AS index_admission,
        dx.hadm_id<>c.hadm_id AND ad.dischtime<=c.t0 AS completed_prior_admission
        FROM cohort c JOIN read_csv_auto('{raw.as_posix()}/admissions.csv.gz') ad ON ad.subject_id=c.subject_id
        JOIN read_csv_auto('{raw.as_posix()}/diagnoses_icd.csv.gz') dx ON dx.hadm_id=ad.hadm_id
        WHERE {codes}"""
    dx=con.execute(sql).df()
    prior=set(dx.loc[dx.completed_prior_admission,'stay_id'])
    labs=con.execute(f"""SELECT c.stay_id,l.labevent_id,l.charttime,
        greatest(l.charttime,l.storetime) AS available_time,l.valuenum
        FROM cohort c JOIN read_csv_auto('{raw.as_posix()}/labevents.csv.gz') l USING (hadm_id)
        WHERE l.itemid=50912 AND l.valuenum>0
          AND l.charttime BETWEEN c.t0-INTERVAL '7 days' AND c.t0""").df()
    con.close()
    late=set()
    for sid,g in labs.groupby('stay_id'):
        t=d.loc[d.stay_id.eq(sid),'t0'].iloc[0];g=n.creat_stages(g)
        if ((g.charttime<=t)&(g.available_time>t)&g.stage_new.eq(3)).any():late.add(sid)
    assert len(late)==2
    cleaned=d[~d.stay_id.isin(prior|late)].copy()
    cleaned.to_parquet(PRIVATE/'cohort.parquet',index=False)
    audit={'historical_n':len(d),'prior_dialysis_new_exclusions':len(prior),'late_prebolus_stage3_exclusions':len(late),
           'overlap_exclusions':len(prior&late),'corrected_n':len(cleaned),
           'prior_codes':dx[dx.completed_prior_admission].groupby(['icd_version','icd_code']).size().to_dict().__str__(),
           'index_diagnosis_patients':int(dx.loc[dx.index_admission,'stay_id'].nunique()),
           'index_diagnosis_patients_retained':int(dx.loc[dx.index_admission&dx.stay_id.isin(cleaned.stay_id),'stay_id'].nunique())}
    dump('eligibility_audit.json',audit);print(audit,flush=True)
    absolute=cleaned.chloride.ge(110);rise=(cleaned.chloride-cleaned.prior_chloride).ge(5)
    assert (absolute|rise).all()
    strata=[('Absolute only',cleaned[absolute&~rise]),('Rise only',cleaned[~absolute&rise]),('Both',cleaned[absolute&rise])]
    profiles=[]
    for name,g in strata:
        profiles.append({'trigger':name,'n':len(g),'ns_n':int(g.exposure.eq(0).sum()),'lr_n':int(g.exposure.eq(1).sum()),
          'chloride_min':float(g.chloride.min()),'chloride_median':float(g.chloride.median()),'chloride_max':float(g.chloride.max()),
          'chloride_below100':int(g.chloride.lt(100).sum()),'chloride_below105':int(g.chloride.lt(105).sum())})
    pd.DataFrame(profiles).to_csv(OUT/'trigger_profiles.csv',index=False)
    lm=cleaned[cleaned.primary_time.gt(6)].copy()
    for ep in ['primary','creatinine3','KRT']:lm[ep+'_time']-=6
    hospital=cleaned.copy();n.set_renal_event(hospital,'primary',hospital.primary_stamp,hospital[['deathtime','dischtime']].min(axis=1))
    jobs=[('Corrected overall',cleaned,n.ENDPOINTS)]+[(name,g,{'primary':168.}) for name,g in strata]+[
      ('Absolute chloride >=110',cleaned[absolute],{'primary':168.}),
      ('No identified prior LR',cleaned[cleaned.lr_ml.eq(0)],{'primary':168.}),
      ('Six-hour landmark',lm,{'primary':162.,'creatinine3':162.,'KRT':162.}),
      ('Hospital follow-up',hospital,{'primary':168.})]
    rows=[];curves=[];balances=[];diags=[]
    for name,g,eps in jobs:
        r,c,b,diag,h=n.estimate(g,m,name,endpoints=eps)
        assert diag['bootstrap_success']==n.PROTOCOL['bootstrap'] and not diag['bootstrap_errors']
        rows+=r;curves+=c;balances+=b;diags.append(diag)
        a.orient(pd.DataFrame(rows)).to_csv(OUT/'estimates.csv',index=False)
        pd.DataFrame(curves).to_csv(OUT/'curves.csv',index=False)
        pd.DataFrame(balances).to_csv(OUT/'balance.csv',index=False)
        dump('diagnostics.json',diags)
    summaries=[]
    for arm,label in [(0,'NS'),(1,'LR')]:
        g=cleaned[cleaned.exposure.eq(arm)]
        summaries.append({'arm':label,'n':len(g),'early6_composite':int((g.primary_event.eq(1)&g.primary_time.le(6)).sum()),
          'early6_krt':int((g.KRT_event.eq(1)&g.KRT_time.le(6)).sum()),'dose_median':float(g.action_amount.median()),
          'baseline_uo_complete':int(g.baseline_uo_complete.sum()),'no_post_creatinine':int(g.creatinine_tests7.eq(0).sum())})
    pd.DataFrame(summaries).to_csv(OUT/'clinical_summary.csv',index=False)
    assert before==source_hashes()
    dump('checks.json',{'historical_sources_unchanged':True,'completed_analyses':len(jobs),'all_bootstraps_successful':True,'source_hashes':before})
if __name__=='__main__':
    with threadpool_limits(limits=1):run()
