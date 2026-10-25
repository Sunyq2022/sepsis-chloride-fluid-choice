"""Authorized-data review of the retained cohort (currently 552 patients); aggregate exports only.
The entry point runs source restrictions and additional estimation, including 2000 bootstrap refits.
It does not replace the original main outputs; it is not a credential-free read-only check.
"""
from pathlib import Path
import importlib.util,json,hashlib,warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'03结果/21_FULL_ICU_CAUSAL_PROGRAM_20260907/N09'
DATA=BASE/'rebuilt_analysis';OUT=DATA/'review_audit'
def load(name):
 p=ROOT/'02代码/01分析脚本'/f'{name}.py';s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def main():
 OUT.mkdir(exist_ok=True);mod=load('n09_corrected_estimand');m=load('canonical_rerun_n02_n07_n09')
 d=pd.read_parquet(DATA/'cohort_private.parquet').reset_index(drop=True);assert d.stay_id.is_unique
 d,prep,x,ps,w=mod.fit(d,m);a=d.exposure.to_numpy(int)
 export=lambda name,rows:pd.DataFrame(rows).to_csv(OUT/name,index=False)
 rows=[]
 for arm in (0,1):
  ix=a==arm;g=d[ix];ww=w[ix]
  def add(domain,measure,value):rows.append(dict(arm='LR' if arm else 'NS',domain=domain,measure=measure,value=value))
  for k in [0,1,2]:add('Baseline','Creatinine stage '+str(k),int(g.aki_stage_creat.eq(k).sum()))
  add('Baseline','Unknown urine stage',int(g.baseline_uo_stage.isna().sum()));add('Baseline','Complete urine window',int(g.baseline_uo_complete.sum()))
  absolute=g.chloride.ge(110);delta=(g.chloride-g.prior_chloride).ge(5)
  assert (absolute|delta).all()
  for title,z in [('Absolute only',absolute&~delta),('Rise only',delta&~absolute),('Both triggers',absolute&delta)]:add('Trigger',title,int(z.sum()))
  add('Prior fluid','No LR before trigger',int(g.lr_ml.eq(0).sum()));add('Prior fluid','LR recorded before trigger',int(g.lr_ml.gt(0).sum()))
  for k in ['age','creat','chloride','sofa_score','ns_ml','lr_ml']:
   ok=g[k].notna();add('Weighted mean',k,float(np.average(g.loc[ok,k],weights=ww[ok])))
  for cat in ['decision_careunit','anchor_year_group']:
   for val in sorted(d[cat].dropna().unique()):add('Weighted percentage',cat+': '+val,float(100*ww[g[cat].eq(val)].sum()/ww.sum()))
  for q in [0,.25,.5,.75,1]:add('Weight quantile',str(q),float(np.quantile(ww,q)))
  for k in mod.NUM+mod.CAT:add('Missing count',k,int(g[k].isna().sum()))
  c=g.creatinine3_event.eq(1);r=g.KRT_event.eq(1)
  for label,z in [('Creatinine only',c&~r),('KRT only',r&~c),('Both components',c&r),('Neither',~(c|r))]:add('Components',label,int(z.sum()))
  assert int((c|r).sum())==int(g.primary_event.eq(1).sum())
  event=g.primary_event.eq(1);early=event&g.primary_time.le(6)
  for lab,z in [('Primary events',event),('Primary events within 6h',early),('KRT first',event&g.KRT_stamp.eq(g.primary_stamp)),('Creatinine first',event&g.creatinine3_stamp.eq(g.primary_stamp)),('Simultaneous components',event&g.creatinine3_stamp.eq(g.KRT_stamp))]:add('Timing',lab,int(z.sum()))
  for label,z in [('Date-only death fallback',g.deathtime.isna()&g.dod.notna()),('Fallback clamped to t0',g.deathtime.isna()&g.dod.notna()&(g.dod.dt.normalize()+pd.Timedelta(hours=12)).lt(g.t0)),('Renal event tied with competing event',g.primary_stamp.eq(g.competing_stamp)),('No creatinine after start',g.creatinine_tests7.eq(0))]:add('Timing',label,int(z.sum()))
  lag=(g.action_storetime-g.t0).dt.total_seconds()/60
  for q in [.25,.5,.75]:add('Recording delay minutes',str(q),float(lag.quantile(q)))
  at=(g[['outtime','deathtime']].min(axis=1)-g.t0).dt.total_seconds().div(3600).clip(upper=168)
  add('Observation','Patient-days alive in index ICU to day7',float(at.sum()/24))
 export('clinical_reporting.csv',rows)
 mix=[]
 for (unit,era),g in d.groupby(['decision_careunit','anchor_year_group']):
  mix.append(dict(unit=unit,anchor_era=era,NS=int(g.exposure.eq(0).sum()),LR=int(g.exposure.eq(1).sum())))
 # Protect small cells in new public support matrix, without silently replacing them by zero.
 for row in mix:
  for k in ['NS','LR']:
   if 0<row[k]<5:row[k]='<5'
 export('unit_anchor_era_support.csv',mix)
 triggers=[]
 for title,ix in [('Absolute only',d.chloride.ge(110)&~(d.chloride-d.prior_chloride).ge(5)),('Rise only',d.chloride.lt(110)),('Both',d.chloride.ge(110)&(d.chloride-d.prior_chloride).ge(5))]:
  g=d[ix];triggers.append(dict(trigger=title,n=len(g),NS=int(g.exposure.eq(0).sum()),LR=int(g.exposure.sum()),age_median=g.age.median(),chloride_median=g.chloride.median(),creatinine_median=g.creat.median(),chloride_samples48h_median=g.chloride_n48.median()))
 export('trigger_characteristics.csv',triggers)
 # Point-estimator identity checks on the frozen cohort, no changed estimates.
 result=pd.read_csv(DATA/'outcomes.csv');errs=[]
 for ep,h in mod.ENDPOINTS.items():
  target=result.query("analysis == 'Main' and outcome == @ep").iloc[0]
  for arm in (0,1):
   ix=a==arm;t=d.loc[ix,ep+'_time'].to_numpy();ev=d.loc[ix,ep+'_event'].to_numpy();ww=w[ix]
   aj=m.weighted_aj(t,ev,ww,h);prop=float(np.average(ev==1,weights=ww));early_censor=bool(((ev==0)&(t<h)).any())
   assert abs(aj-target['risk_lr' if arm else 'risk_ns'])<1e-10
   if not early_censor:assert abs(aj-prop)<1e-10
   errs.append(dict(outcome=ep,arm=arm,early_censor=early_censor,AJ=aj,weighted_event_proportion=prop))
 export('estimator_identity.csv',errs)
 # Main fit convergence is checked now; historic bootstrap convergence was not fully logged.
 with warnings.catch_warnings(record=True) as warn:
  model=LogisticRegression(C=.5,max_iter=2000);model.fit(x,a)
 metadata={'n':len(d),'NS':int((a==0).sum()),'LR':int((a==1).sum()),'encoded_dimension':x.shape[1], 'solver':model.solver,'tol':model.tol,'iterations':model.n_iter_.tolist(),'warnings':[str(v.message) for v in warn],'seed':mod.PROTOCOL['seed'],'outcomes_sha256':hashlib.sha256((DATA/'outcomes.csv').read_bytes()).hexdigest(),'historical_bootstrap_not_rerun':True,'clinical_adjudication_completed':False}
 # Source-level lab timing audit of this exact cohort. No new clinical adjudication claim.
 con=m.connection();con.execute('set threads=1');con.register('review_people',d[['stay_id','hadm_id','t0','intime','outtime']])
 sql="""select c.stay_id,l.labevent_id,l.charttime,greatest(l.charttime,l.storetime) available_time,l.valuenum from review_people c join b.mimiciv_hosp.labevents l on l.hadm_id=c.hadm_id where l.itemid=50912 and l.valuenum>0 and l.charttime<=c.t0+interval '7d'"""
 raw=con.execute(sql).df();con.close();(OUT/'lab_timing_audit.sql').write_text(sql,encoding='utf-8')
 flags=[];late_ids=[]
 for sid,g in raw.groupby('stay_id'):
  r=d[d.stay_id.eq(sid)].iloc[0];v=mod.creat_stages(g)
  gap=v.charttime.le(r.t0)&v.available_time.gt(r.t0)&v.stage_new.eq(3)
  pre=v.charttime.le(r.t0)&v.charttime.ge(r.intime)&v.stage_new.eq(3)
  if gap.any():late_ids.append(sid)
  flags.append(dict(arm=int(r.exposure),pre_sample_late_result_stage3=bool(gap.any()),retrospective_pre_stage3=bool(pre.any()),early_primary=bool(r.primary_event==1 and r.primary_time<=6)))
 flags=pd.DataFrame(flags);agg=flags.groupby('arm').sum().reset_index();agg.to_csv(OUT/'lab_timing_counts.csv',index=False)
 metadata['lab_timing_counts']=agg.to_dict('records');(OUT/'audit_checks.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8');print(json.dumps(metadata),flush=True)
 jobs=[('Main 2000 bootstrap',d,2000),('Exclude pre-sampled late-result stage3',d[~d.stay_id.isin(late_ids)].reset_index(drop=True),500),('Exclude zero recorded volume',d[d.action_amount.gt(0)].reset_index(drop=True),500)]
 results=[];diagnostics=[]
 for name,frame,reps in jobs:
  mod.PROTOCOL['bootstrap']=reps
  r,_,_,diag,_=mod.estimate(frame,m,name,False,{'primary':168.})
  results+=r;diagnostics.append(diag)
 export('review_sensitivity.csv',results);(OUT/'review_sensitivity_checks.json').write_text(json.dumps(diagnostics,indent=2),encoding='utf-8')
if __name__=='__main__':
 with threadpool_limits(limits=1):main()
