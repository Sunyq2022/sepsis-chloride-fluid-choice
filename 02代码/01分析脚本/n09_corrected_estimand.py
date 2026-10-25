"""N09 complete reconstruction. Explicit extract/analyse modes; sources opened read-only.
Historical implementation and results remain in Git and their original result directories.
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, time, warnings
warnings.filterwarnings("ignore", message="The .generic. unit for NumPy timedelta is deprecated", category=DeprecationWarning)
from pathlib import Path
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'03结果/21_FULL_ICU_CAUSAL_PROGRAM_20260907/N09'
OUT=BASE/'rebuilt_analysis'
CODE=Path(__file__).parent
PROTOCOL=json.loads((OUT/'protocol.json').read_text(encoding='utf-8'))
NUM=PROTOCOL['numeric'];CAT=PROTOCOL['categorical']
ENDPOINTS={'primary':168.,'creatinine3':168.,'KRT':168.,'urine3':168.,'broader_renal':168.,'death7':168.,'death28':672.,'death90':2160.}
LABELS={'primary':'Creatinine stage 3 or treatment-supported KRT','creatinine3':'Creatinine-defined stage 3','KRT':'Treatment-supported KRT','urine3':'Urine-output-defined stage 3','broader_renal':'Broader renal composite including urine output','death7':'7-day mortality','death28':'28-day mortality','death90':'90-day mortality'}
PROCEDURES='225441,225802,225803,225805,225809,225955'
FLOW_ITEMS='224144,224154,226499,226457,224191'

def load(name):
 s=importlib.util.spec_from_file_location(name,CODE/(name+'.py'));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def write_json(name,value):
 (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,default=str),encoding='utf-8',newline='\n')

def stage(creat,low7,low48):
 if not np.isfinite(creat) or not np.isfinite(low7):return np.nan
 if creat>=3*low7 or (creat>=4 and (creat>=1.5*low7 or (np.isfinite(low48) and creat>=low48+.3))):return 3
 if creat>=2*low7:return 2
 if creat>=1.5*low7 or (np.isfinite(low48) and creat>=low48+.3):return 1
 return 0

def creat_stages(g):
 g=g.sort_values(['charttime','available_time','labevent_id']).drop_duplicates('charttime',keep='last').copy()
 t=g.charttime.to_numpy('datetime64[ns]');v=g.valuenum.to_numpy(float)
 lo7=[];lo48=[];st=[]
 for j in range(len(g)):
  ref=v[(t<t[j])&(t>=t[j]-np.timedelta64(7,'D'))];ref48=v[(t<t[j])&(t>=t[j]-np.timedelta64(48,'h'))]
  a=float(np.min(ref)) if len(ref) else np.nan;b=float(np.min(ref48)) if len(ref48) else np.nan
  lo7.append(a);lo48.append(b);st.append(stage(v[j],a,b))
 g['low7_new']=lo7;g['low48_new']=lo48;g['stage_new']=st
 return g

def extract(audit_only=False):
 m=load('canonical_rerun_n02_n07_n09');bias=load('n09_bias_audit');conf=load('n09_pre_t0_confounder_audit');term=load('n09_terminal_bias_audit')
 con=m.connection();con.execute('set threads=1')
 queries=[]
 def q(sql):queries.append(sql);return con.execute(sql)
 risk=bias.RISKSET_SQL.replace("i.statusdescription='FinishedRunning'","i.statusdescription is distinct from 'Rewritten'")
 q('create temp table rebuild_risk as '+risk)
 pre=risk.rsplit('select * from prior_fluid where ns_ml>=500',1)[0]
 flow=q(pre+"select 'Adult first ICU source' step,count(*) n from s.main.people where age>=18 union all select 'Sepsis source',count(distinct p.stay_id) from s.main.people p join b.mimiciv_derived.sepsis3 z using(stay_id,subject_id) where p.age>=18 and z.sepsis3 union all select 'Eligible chloride measurement',count(distinct stay_id) from chloride_results union all select 'First chloride trigger',count(*) from trigger union all select 'Prior saline at least 500mL',count(*) from prior_fluid where ns_ml>=500").df().to_dict('records')
 q("""create temp table rebuild_bolus as select i.* from b.mimiciv_icu.inputevents i semi join rebuild_risk r using(stay_id)
 where itemid in(225158,225828) and ordercategoryname='03-IV Fluid Bolus' and statusdescription is distinct from 'Rewritten'
 qualify row_number() over(partition by stay_id,itemid,coalesce(linkorderid,orderid) order by starttime,orderid)=1""")
 q("""create temp table rebuild_actions as select r.*,i.itemid,i.starttime t0,i.endtime action_end,i.amount action_amount,i.originalamount,i.rate,i.storetime action_storetime,i.statusdescription,
 (select count(distinct x.itemid) from rebuild_bolus x where x.stay_id=r.stay_id and x.starttime=i.starttime) simultaneous_fluids
 from rebuild_risk r join lateral(select * from rebuild_bolus i where i.stay_id=r.stay_id and i.starttime>=r.trigger_time and i.starttime<r.trigger_time+interval '4h' order by starttime,itemid,orderid limit 1)i on true""")
 n=q('select count(*) from rebuild_actions').fetchone()[0];flow.append({'step':'First documented bolus of either fluid, any size','n':n})
 feature=bias.FEATURE_SQL.split('where coalesce(r.aki_stage_creat,0)<3')[0];parts=[]
 for offset in range(0,n,32):
  q(f'create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset {offset}')
  parts.append(q(feature.format(action_table='rebuild_batch')).df())
 d=pd.concat(parts,ignore_index=True);con.register('rebuild_people',d)
 q("create temp table rebuild_transfers as select t.* from read_csv_auto('E:/ICUdatabase/MIMIC/mimic-iv/hosp/transfers.csv.gz') t semi join rebuild_people c using(hadm_id)")
 timing=q("""select c.stay_id,a.admittime,a.dischtime,a.deathtime,p.dod,
 (select arg_max(careunit,intime) from rebuild_transfers t where t.hadm_id=c.hadm_id and t.intime<=c.t0 and (t.outtime>c.t0 or t.outtime is null)) decision_careunit,
 exists(select 1 from b.mimiciv_hosp.admissions old join read_csv_auto('E:/ICUdatabase/MIMIC/mimic-iv/hosp/diagnoses_icd.csv.gz') dx on dx.hadm_id=old.hadm_id where old.subject_id=c.subject_id and old.dischtime<=c.t0 and old.hadm_id<>c.hadm_id and replace(upper(dx.icd_code),'.','') in('5856','V451','N186','Z992')) known_prior_ESRD
 from rebuild_people c join b.mimiciv_hosp.admissions a using(hadm_id) join b.mimiciv_hosp.patients p on p.subject_id=c.subject_id""").df()
 d=d.merge(timing,on='stay_id',validate='one_to_one')
 cr=q("""select c.stay_id,l.labevent_id,l.charttime,greatest(l.charttime,l.storetime) available_time,l.valuenum from rebuild_people c join b.mimiciv_hosp.labevents l on l.hadm_id=c.hadm_id where l.itemid=50912 and l.valuenum>0 and l.charttime<=c.t0+interval '7d'""").df()
 uo=q("""select c.stay_id,k.charttime,k.aki_stage_uo,u.uo_tm_24hr,u.uo_rt_24hr from rebuild_people c join b.mimiciv_derived.kdigo_stages k on k.stay_id=c.stay_id left join b.mimiciv_derived.kdigo_uo u on u.stay_id=k.stay_id and u.charttime=k.charttime where k.charttime<=c.t0+interval '7d' and k.aki_stage_uo is not null""").df()
 kr=q(f"""select c.stay_id,p.starttime event_time,'procedure' evidence,p.itemid from rebuild_people c join b.mimiciv_icu.procedureevents p on p.stay_id=c.stay_id where p.itemid in({PROCEDURES}) and p.value is not null and p.statusdescription is distinct from 'Rewritten' and p.starttime<=c.t0+interval '7d'
 union all select c.stay_id,p.charttime,'positive flow or output',p.itemid from rebuild_people c join b.mimiciv_icu.chartevents p on p.stay_id=c.stay_id where p.itemid in({FLOW_ITEMS}) and p.valuenum>0 and p.charttime<=c.t0+interval '7d'""").df()
 groups={k:{int(i):g for i,g in f.groupby('stay_id')} for k,f in [('cr',cr),('uo',uo),('kr',kr)]}
 rows=[]
 for row in d.itertuples():
  sid=int(row.stay_id);t0=pd.Timestamp(row.t0)
  raw=groups['cr'].get(sid,cr.iloc[:0]);known=raw[(raw.available_time<=t0)&(raw.charttime<=t0)]
  precr=creat_stages(known);allcr=creat_stages(raw)
  recent=precr[precr.charttime>=t0-pd.Timedelta(hours=24)]
  latest=recent.iloc[-1] if len(recent) else None
  ug=groups['uo'].get(sid,uo.iloc[:0]);preu=ug[ug.charttime<=t0].sort_values('charttime');lastu=preu.iloc[-1] if len(preu) else None
  kg=groups['kr'].get(sid,kr.iloc[:0]);postk=kg[kg.event_time>t0];postc=allcr[(allcr.charttime>t0)&allcr.stage_new.eq(3)];postu=ug[(ug.charttime>t0)&ug.aki_stage_uo.ge(3)]
  rows.append({'stay_id':sid,'baseline_available':latest is not None and np.isfinite(latest.low7_new),'creat_rebuilt':float(latest.valuenum) if latest is not None else np.nan,'low7_rebuilt':float(latest.low7_new) if latest is not None else np.nan,'stage_rebuilt':float(latest.stage_new) if latest is not None else np.nan,
   'prior_creat3':bool(precr[(precr.charttime>=row.intime)].stage_new.eq(3).any()),'prior_uo3':bool(preu.aki_stage_uo.ge(3).any()),'prior_krt':bool((kg.event_time<=t0).any()),
   'baseline_uo_stage':float(lastu.aki_stage_uo) if lastu is not None else np.nan,'baseline_uo_complete':bool(lastu is not None and pd.notna(lastu.uo_tm_24hr) and lastu.uo_tm_24hr>=24 and pd.notna(lastu.uo_rt_24hr)),
   'creatinine3_stamp':postc.charttime.min(),'KRT_stamp':postk.event_time.min(),'urine3_stamp':postu.charttime.min(),'procedure_stamp':postk.loc[postk.evidence.eq('procedure'),'event_time'].min(),
   'creatinine_tests7':int(((raw.charttime>t0)&(raw.charttime<=min(pd.Timestamp(row.outtime),t0+pd.Timedelta(days=7)))).sum())})
 d=d.merge(pd.DataFrame(rows),on='stay_id',validate='one_to_one')
 for a,b in [('creat','creat_rebuilt'),('low7','low7_rebuilt'),('aki_stage_creat','stage_rebuilt')]:d[a]=d[b]
 alive=(d.deathtime.isna()|d.deathtime.gt(d.t0)) & d.outtime.gt(d.t0) & d.dischtime.gt(d.t0) & (d.dod.isna()|pd.to_datetime(d.dod).dt.normalize().ge(pd.to_datetime(d.t0).dt.normalize()))
 masks=[('Exclude ambiguous simultaneous fluids',d.simultaneous_fluids.ne(1)),('Exclude impossible alive or ICU chronology',~alive),('Exclude records without valid administration timing',d.action_end.lt(d.t0)),('Exclude unavailable baseline creatinine/comparator',~d.baseline_available),('Exclude documented prior creatinine stage 3',d.prior_creat3),('Exclude documented prior urine-output stage 3',d.prior_uo3),('Exclude prior treatment-supported KRT',d.prior_krt),('Exclude prior documented ESRD/dialysis dependence',d.known_prior_ESRD)]
 keep=pd.Series(True,index=d.index)
 for label,bad in masks:
  removed=int((keep&bad).sum());keep&=~bad;flow.append({'step':label,'n':removed,'remaining':int(keep.sum())})
 if audit_only:
  cohort_migration(d,keep,masks)
  (OUT/'cohort_migration.sql').write_text('\n;\n'.join(queries),encoding='utf-8')
  con.close();return
 candidates=len(d);d=d[keep].copy().reset_index(drop=True)
 assert d.stay_id.is_unique and d.subject_id.is_unique and d.baseline_available.all() and not d[['prior_creat3','prior_uo3','prior_krt']].any().any()
 d=term.add_parsimonious_features(con,d,conf)
 # Reuse strictly pre-t0 measurements; no amount of the index treatment enters the model.
 d['primary_stamp']=d[['creatinine3_stamp','KRT_stamp']].min(axis=1)
 d['broader_renal_stamp']=d[['primary_stamp','urine3_stamp']].min(axis=1)
 death=pd.to_datetime(d.deathtime);fallback=pd.to_datetime(d.dod).dt.normalize()+pd.Timedelta(hours=12)
 d['death_stamp']=death.fillna(fallback);d.loc[d.death_stamp.lt(d.t0)&d.deathtime.isna(),'death_stamp']=d.loc[d.death_stamp.lt(d.t0)&d.deathtime.isna(),'t0']
 d['competing_stamp']=d[['deathtime','outtime']].min(axis=1)
 for ep,horizon in ENDPOINTS.items():
  if ep.startswith('death'):
   t=(d.death_stamp-d.t0).dt.total_seconds()/3600;event=t.notna()&t.ge(0)&t.le(horizon);d[ep+'_time']=np.where(event,t,horizon);d[ep+'_event']=event.astype(int)
  else:set_renal_event(d,ep,d[ep+'_stamp'],d.competing_stamp,horizon)
 assert all(d[e+'_time'].ge(0).all() for e in ENDPOINTS)
 d.to_parquet(OUT/'cohort_private.parquet',index=False)
 pd.DataFrame(flow).to_csv(OUT/'flow.csv',index=False)
 items=pd.read_csv('E:/ICUdatabase/MIMIC/mimic-iv/icu/d_items.csv.gz');items[items.itemid.isin([int(i) for i in (PROCEDURES+','+FLOW_ITEMS).split(',')])].to_csv(OUT/'KRT_items.csv',index=False)
 kr_summary=kr.groupby(['evidence','itemid']).size().rename('source_records').reset_index();kr_summary.to_csv(OUT/'KRT_source_dictionary.csv',index=False)
 profile=q("select statusdescription,itemid,count(*) n,median(action_amount) median_amount,min(action_amount) min_amount,median(datediff('second',t0,action_end)) encoded_seconds from rebuild_actions group by 1,2").df();profile.to_csv(OUT/'action_profile.csv',index=False)
 write_json('extraction_checks.json',{'candidate_actions':candidates,'n':len(d),'lr_n':int(d.exposure.sum()),'ns_n':int(d.exposure.eq(0).sum()),'baseline_uo_unknown':int(d.baseline_uo_stage.isna().sum()),'baseline_uo_complete':int(d.baseline_uo_complete.sum()),'careunit_different_from_first':int(d.decision_careunit.ne(d.first_careunit).sum()),'unique_patient':True,'no_negative_followup':True,'no_prior_documented_stage3_or_KRT':True,'primary_events':int(d.primary_event.eq(1).sum()),'source_databases_read_only':True,'protocol_sha256':hashlib.sha256((OUT/'protocol.json').read_bytes()).hexdigest()})
 (OUT/'extraction.sql').write_text('\n;\n'.join(queries),encoding='utf-8')
 con.close();print('EXTRACTION COMPLETE',len(d),'LR',int(d.exposure.sum()),flush=True)

def set_renal_event(df,ep,stamp,comp,horizon=168.):
 t=(pd.to_datetime(stamp)-df.t0).dt.total_seconds()/3600;c=(pd.to_datetime(comp)-df.t0).dt.total_seconds()/3600
 df[ep+'_time']=pd.concat([t,c,pd.Series(horizon,index=df.index)],axis=1).min(axis=1)
 df[ep+'_event']=np.select([t.notna()&t.gt(0)&t.le(c.fillna(np.inf))&t.le(horizon),c.notna()&c.lt(t.fillna(np.inf))&c.le(horizon)],[1,2],default=0)

def cif_summary(t,e,w,horizon,grid):
 # Weighted Aalen-Johansen at exact event times; integral uses all event times.
 times,inv=np.unique(np.asarray(t,float),return_inverse=True);allw=np.bincount(inv,weights=w);risk=np.cumsum(allw[::-1])[::-1]
 ev1=np.bincount(inv,weights=w*(e==1),minlength=len(times));evall=np.bincount(inv,weights=w*(e>0),minlength=len(times))
 hazards=np.divide(evall,risk,out=np.zeros_like(risk),where=risk>0);surv=np.cumprod(1-hazards);before=np.r_[1.,surv[:-1]]
 cif=np.cumsum(before*np.divide(ev1,risk,out=np.zeros_like(risk),where=risk>0))
 idx=np.searchsorted(times,grid,side='right')-1;values=np.where(idx>=0,cif[np.maximum(idx,0)],0)
 sel=times<=horizon;ts=times[sel];cs=cif[sel]
 area=float(np.sum(np.diff(np.r_[0,ts,horizon])*np.r_[0,cs]))
 return float(cs[-1] if len(cs) else 0),area,values

def fit(df,m,unpenalized=False,mult=None):
 if not unpenalized:return m.ato_fit(df,NUM,CAT,mult)
 from sklearn.linear_model import LogisticRegression
 if mult is not None:df=df.loc[df.index.repeat(mult)].reset_index(drop=True)
 prep=m.ato_design(df,NUM,CAT);x=prep.fit_transform(df[NUM+CAT]);model=LogisticRegression(penalty=None,max_iter=5000);model.fit(x,df.exposure)
 ps=np.clip(model.predict_proba(x)[:,1],1e-6,1-1e-6);w=np.where(df.exposure.eq(1),1-ps,ps);return df,prep,x,ps,w

def estimate(frame,m,name,unpenalized=False,endpoints=None):
 frame=frame.reset_index(drop=True);eps=endpoints or ENDPOINTS
 d,prep,x,ps,w=fit(frame,m,unpenalized);a=d.exposure.to_numpy(int)
 grids={ep:np.linspace(0,h,57) for ep,h in eps.items()}
 def effects(dd,ww):
  aa=dd.exposure.to_numpy(int);out={}
  for ep,h in eps.items():
   val=[];ar=[];cv=[]
   for arm in (0,1):
    ix=aa==arm;r,area,curve=cif_summary(dd.loc[ix,ep+'_time'].to_numpy(float),dd.loc[ix,ep+'_event'].to_numpy(int),ww[ix],h,grids[ep]);val.append(r);ar.append(h-area if ep.startswith('death') else area);cv.append(curve)
   out[ep]={'risk_ns':val[0],'risk_lr':val[1],'rd':val[1]-val[0],'rr':val[1]/val[0] if val[0]>0 else np.nan,'time_ns':ar[0],'time_lr':ar[1],'time_difference':ar[1]-ar[0],'curves':cv}
  return out
 point=effects(d,w);boots=[];curveboots={ep:[] for ep in eps};errors=[];rng=np.random.default_rng(PROTOCOL['seed'])
 for i in range(PROTOCOL['bootstrap']):
  mult=np.bincount(rng.integers(0,len(frame),len(frame)),minlength=len(frame))
  try:
   dd,_,_,_,ww=fit(frame,m,unpenalized,mult);b=effects(dd,ww)
   boots.append({ep:{k:v for k,v in z.items() if k!='curves'} for ep,z in b.items()})
   for ep in eps:curveboots[ep].append(b[ep]['curves'])
  except Exception as exc:errors.append({'replicate':i,'error':repr(exc)})
  if (i+1)%100==0:print(name,'bootstrap',i+1,'success',len(boots),flush=True)
 results=[];curves=[]
 for ep,h in eps.items():
  p=point[ep];r={'analysis':name,'outcome':ep,'label':LABELS.get(ep,ep),'n':len(d),'ns_n':int((a==0).sum()),'lr_n':int((a==1).sum()),'ns_events':int(((a==0)&d[ep+'_event'].eq(1)).sum()),'lr_events':int(((a==1)&d[ep+'_event'].eq(1)).sum()),'bootstrap_success':len(boots)}
  for key,value in p.items():
   if key=='curves':continue
   r[key]=value;vv=np.array([b[ep][key] for b in boots],float);finite=vv[np.isfinite(vv)];r[key+'_low'],r[key+'_high']=np.quantile(finite,[.025,.975]) if len(finite) else (np.nan,np.nan)
  results.append(r)
  cb=np.array(curveboots[ep]);lo,hi=np.quantile(cb,[.025,.975],axis=0) if len(cb) else (np.full((2,57),np.nan),np.full((2,57),np.nan))
  for arm in (0,1):
   for j,t in enumerate(grids[ep]):curves.append({'analysis':name,'outcome':ep,'arm':'LR' if arm else 'NS','hour':t,'risk':p['curves'][arm][j],'lower':lo[arm,j],'upper':hi[arm,j],'at_risk':int(((a==arm)&d[ep+'_time'].ge(t)).sum())})
 dense=x.toarray() if hasattr(x,'toarray') else x;bal=[]
 for j,feature in enumerate(prep.get_feature_names_out()):
  v=dense[:,j];sd=np.sqrt((v[a==0].var()+v[a==1].var())/2)
  raw=abs(v[a==1].mean()-v[a==0].mean())/sd if sd else 0
  adj=abs(np.average(v[a==1],weights=w[a==1])-np.average(v[a==0],weights=w[a==0]))/sd if sd else 0
  bal.append({'analysis':name,'variable':feature,'smd_before':raw,'smd_after':adj})
 diag={'analysis':name,'n':len(d),'lr_n':int(a.sum()),'ns_n':int((a==0).sum()),'ess_lr':float(w[a==1].sum()**2/(w[a==1]**2).sum()),'ess_ns':float(w[a==0].sum()**2/(w[a==0]**2).sum()),'max_smd':float(max(z['smd_after'] for z in bal)),'bootstrap_success':len(boots),'bootstrap_errors':errors}
 # No individual propensity scores are exported.
 hist=[]
 for arm in (0,1):
  cnt,edges=np.histogram(ps[a==arm],bins=np.linspace(0,1,21))
  hist += [{'analysis':name,'arm':'LR' if arm else 'NS','left':edges[j],'right':edges[j+1],'count':int(v)} for j,v in enumerate(cnt)]
 return results,curves,bal,diag,hist

def analyse():
 m=load('canonical_rerun_n02_n07_n09');d=pd.read_parquet(OUT/'cohort_private.parquet').sort_values(['subject_id','stay_id']).reset_index(drop=True)
 # Runnable estimator checks independent of this study's expected effect direction.
 r,area,_=cif_summary(np.array([1.,2.,3.]),np.array([1,2,0]),np.ones(3),3.,np.array([0.,1.,3.]));assert np.isclose(r,1/3) and np.isclose(area,2/3)
 assert not set(NUM+CAT)&{'action_amount','originalamount','first_careunit'}
 jobs=[('Main',d,False,ENDPOINTS),('Bolus within 2h',d[(d.t0-d.trigger_time).dt.total_seconds()<7200],False,{'primary':168.}),('No prior LR',d[d.lr_ml.eq(0)],False,{'primary':168.}),('Absolute chloride >=110',d[d.chloride.ge(110)],False,{'primary':168.}),('Complete baseline urine window',d[d.baseline_uo_complete],False,{'primary':168.,'broader_renal':168.}),('Unpenalized propensity model',d,True,{'primary':168.})]
 hospital=d.copy();set_renal_event(hospital,'primary',hospital.primary_stamp,hospital[['deathtime','dischtime']].min(axis=1));jobs.append(('Hospital creatinine or recorded KRT',hospital,False,{'primary':168.}))
 procedure=d.copy();set_renal_event(procedure,'primary',procedure[['creatinine3_stamp','procedure_stamp']].min(axis=1),procedure.competing_stamp);jobs.append(('Procedure-only KRT component',procedure,False,{'primary':168.}))
 date=d.copy();delta=(pd.to_datetime(date.dod).dt.normalize()-date.t0.dt.normalize()).dt.total_seconds()/3600
 for ep,h in ENDPOINTS.items():
  if ep.startswith('death'):
   valid=delta.notna()&delta.ge(0)&delta.le(h);date[ep+'_time']=np.where(valid,delta,h);date[ep+'_event']=valid.astype(int)
 jobs.append(('Date-resolved mortality',date,False,{k:h for k,h in ENDPOINTS.items() if k.startswith('death')}))
 allr=[];allc=[];allb=[];diags=[];allh=[]
 for name,frame,unpen,eps in jobs:
  if min(frame.exposure.value_counts())<5:raise ValueError('Insufficient support in fixed analysis '+name)
  results,curves,bal,diag,hist=estimate(frame,m,name,unpen,eps);allr+=results;allc+=curves;allb+=bal;diags.append(diag);allh+=hist
  pd.DataFrame(allr).to_csv(OUT/'outcomes.csv',index=False);pd.DataFrame(allc).to_csv(OUT/'curves.csv',index=False);pd.DataFrame(allb).to_csv(OUT/'balance.csv',index=False);pd.DataFrame(allh).to_csv(OUT/'propensity_histogram.csv',index=False);write_json('diagnostics.json',diags)
 reporting_summaries(d)
 write_json('run_checks.json',{'status':'complete','n':len(d),'all_fixed_analyses_completed':len(jobs),'bootstrap_per_analysis':PROTOCOL['bootstrap'],'all_bootstraps_successful':all(x['bootstrap_success']==500 for x in diags),'same_main_cohort_for_all_main_endpoints':True,'endpoint_and_treatment_revisions_disclosed':True,'primary_protocol_sha256':hashlib.sha256((OUT/'protocol.json').read_bytes()).hexdigest(),'analysis_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'historical706_not_used_as_rebuilt_main':True})
 print('ALL FIXED ANALYSES COMPLETE',flush=True)

def reporting_summaries(d):
 # Complete baseline table, with missing counts separated for main/supp presentation.
 rows=[]
 for col in NUM+['action_amount','baseline_uo_stage','baseline_uo_complete']+CAT:
  for arm in (0,1):
   x=d.loc[d.exposure.eq(arm),col]
   if col in CAT or col in ('baseline_uo_complete','pressor','ventilation'):
    for level,count in x.value_counts(dropna=False).items():rows.append({'variable':col,'level':str(level),'arm':'LR' if arm else 'NS','summary':f'{count} ({100*count/len(x):.1f}%)','missing':int(x.isna().sum()),'n':len(x)})
   else:rows.append({'variable':col,'level':'','arm':'LR' if arm else 'NS','summary':f'{x.median():.2f} ({x.quantile(.25):.2f}, {x.quantile(.75):.2f})','missing':int(x.isna().sum()),'n':len(x)})
 pd.DataFrame(rows).to_csv(OUT/'baseline.csv',index=False)
 # Descriptive competing events, first-event timing and post-t0 observation coverage.
 obs=[];first=[];components=[]
 for arm in (0,1):
  a=d.exposure.eq(arm);dd=d[a];death=pd.to_datetime(dd.deathtime).notna()&pd.to_datetime(dd.deathtime).le(dd.outtime)
  event=np.select([dd.primary_event.eq(1),dd.primary_event.eq(2)&death,dd.primary_event.eq(2)&~death],[1,2,3],default=0)
  for cause,label in [(1,'Primary kidney event'),(2,'Death before kidney event'),(3,'ICU exit before kidney event')]:
   ev=np.where(event==cause,1,np.where(event>0,2,0));grid=np.linspace(0,168,57);risk,_,values=cif_summary(dd.primary_time.to_numpy(float),ev,np.ones(len(dd)),168,grid)
   for t,v in zip(grid,values):obs.append({'arm':'LR' if arm else 'NS','cause':label,'hour':t,'risk':v,'count7':int((event==cause).sum()),'n':len(dd)})
  for ep in ['primary','broader_renal','creatinine3','KRT','urine3']:
   for lo,hi,label in [(0,6,'0–6h'),(6,24,'6–24h'),(24,168,'24–168h')]:first.append({'arm':'LR' if arm else 'NS','outcome':ep,'time_bin':label,'events':int((dd[ep+'_event'].eq(1)&dd[ep+'_time'].gt(lo)&dd[ep+'_time'].le(hi)).sum())})
  for one in ['creatinine3','KRT','urine3']:
   for two in ['creatinine3','KRT','urine3']:components.append({'arm':'LR' if arm else 'NS','component1':one,'component2':two,'overlap':int((dd[one+'_event'].eq(1)&dd[two+'_event'].eq(1)).sum())})
 pd.DataFrame(obs).to_csv(OUT/'observation_curves.csv',index=False);pd.DataFrame(first).to_csv(OUT/'event_timing.csv',index=False);pd.DataFrame(components).to_csv(OUT/'component_overlap.csv',index=False)
 rows=[]
 for arm in (0,1):
  dd=d[d.exposure.eq(arm)]
  for col in ['primary_time','creatinine_tests7','action_amount']:
   x=dd[col];rows.append({'arm':'LR' if arm else 'NS','measure':col,'n':len(x),'median':float(x.median()),'q25':float(x.quantile(.25)),'q75':float(x.quantile(.75)),'zero':int(x.eq(0).sum()),'missing':int(x.isna().sum())})
 pd.DataFrame(rows).to_csv(OUT/'monitoring.csv',index=False)

def validate():
 m=load('canonical_rerun_n02_n07_n09');d=pd.read_parquet(OUT/'cohort_private.parquet').sort_values(['subject_id','stay_id']).reset_index(drop=True)
 dd,_,_,_,w=fit(d,m);a=dd.exposure.to_numpy(int);r=pd.read_csv(OUT/'outcomes.csv');main=r[r.analysis.eq('Main')]
 assert len(d)>0 and d.subject_id.is_unique and d.stay_id.is_unique
 assert d.exposure.eq(d.itemid.eq(225828).astype(int)).all()
 assert not d[['prior_creat3','prior_uo3','prior_krt']].any().any()
 assert d.baseline_available.all() and d.outtime.gt(d.t0).all() and (d.deathtime.isna()|d.deathtime.gt(d.t0)).all()
 assert d.primary_stamp.equals(d[['creatinine3_stamp','KRT_stamp']].min(axis=1))
 maxerror=0.
 for ep,h in ENDPOINTS.items():
  z=main[main.outcome.eq(ep)].iloc[0]
  for arm,key in [(0,'risk_ns'),(1,'risk_lr')]:
   ix=a==arm;t=dd.loc[ix,ep+'_time'].to_numpy(float);ev=dd.loc[ix,ep+'_event'].to_numpy(int)
   independent=m.weighted_aj(t,ev,w[ix],h);err=abs(independent-z[key]);maxerror=max(maxerror,err);assert err<1e-10
   assert (t>=0).all() and (t<=h).all()
  assert z.rd_low<=z.rd_high and z.bootstrap_success==500
 # Independent event-time integration and known competing/censoring test.
 t=np.array([1.,2.,2.,3.]);e=np.array([1,2,0,1]);w0=np.ones(4);risk,area,_=cif_summary(t,e,w0,4.,np.array([0.,1.,2.,3.,4.]))
 assert np.isclose(risk,.75) and np.isclose(area,1.25)
 checks={'cohort_and_baseline_checks':'pass','all_8_main_endpoint_risks_match_independent_AJ':True,'max_absolute_error':maxerror,'known_competing_and_censoring_example':'pass','zero_recorded_LR_dose':int((d.exposure.eq(1)&d.action_amount.eq(0)).sum()),'complete_urine_window':int(d.baseline_uo_complete.sum()),'limitations_are_not_validation_passes':True,'analysis_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'outcomes_sha256':hashlib.sha256((OUT/'outcomes.csv').read_bytes()).hexdigest()}
 write_json('independent_checks.json',checks);print('INDEPENDENT POINT-ESTIMATE AND COHORT CHECKS PASS',maxerror,flush=True)

def cohort_migration(candidates,keep,masks):
 old=pd.read_parquet(BASE/'publication_figures/cohort_private.parquet');current=pd.read_parquet(OUT/'cohort_private.parquet')
 assert len(old)==706 and old.stay_id.is_unique and len(current)>0
 new=candidates[keep];assert set(new.stay_id)==set(current.stay_id)
 oi=set(old.stay_id);ni=set(new.stay_id);ci=set(candidates.stay_id)
 reason=pd.Series('Retained',index=candidates.index);alive=pd.Series(True,index=candidates.index)
 for label,bad in masks:reason.loc[alive&bad]=label;alive&=~bad
 candidates=candidates.copy();candidates['first_exclusion']=reason
 rows=[]
 for membership,subset in [('Old 706 retained',candidates[candidates.stay_id.isin(oi&ni)]),('Old 706 excluded',candidates[candidates.stay_id.isin(oi-ni)]),('New entrants',candidates[candidates.stay_id.isin(ni-oi)])]:
  for label,g in subset.groupby('first_exclusion'):rows.append({'membership':membership,'reason':label,'n':len(g),'NS_at_new_t0':int(g.exposure.eq(0).sum()),'LR_at_new_t0':int(g.exposure.eq(1).sum())})
 if oi-ci:rows.append({'membership':'Old 706 excluded','reason':'No qualifying first-start candidate in rebuilt source','n':len(oi-ci),'NS_at_new_t0':None,'LR_at_new_t0':None})
 pd.DataFrame(rows).to_csv(OUT/'cohort_migration.csv',index=False)
 both=old[['stay_id','t0','exposure']].merge(new[['stay_id','t0','exposure']],on='stay_id',suffixes=('_old','_new'));trans=both.groupby(['exposure_old','exposure_new']).size().rename('n').reset_index();trans.to_csv(OUT/'cohort_treatment_transition.csv',index=False)
 entrants=new[new.stay_id.isin(ni-oi)]
 report={'old_n':706,'new_n':len(current),'retained_from_old':len(oi&ni),'removed_from_old':len(oi-ni),'new_entrants':len(ni-oi),'old_without_rebuilt_candidate':len(oi-ci),'retained_with_changed_t0':int(both.t0_old.ne(both.t0_new).sum()),'retained_with_changed_fluid':int(both.exposure_old.ne(both.exposure_new).sum()),'new_entrants_recorded_volume_below250':int(entrants.action_amount.lt(250).sum()),'new_entrants_NS':int(entrants.exposure.eq(0).sum()),'new_entrants_LR':int(entrants.exposure.eq(1).sum()),'equation':f'706 - {len(oi-ni)} + {len(ni-oi)} = {len(current)}','patient_identifiers_exported':False,'effect_estimates_rerun':False,'classification':'Sequential first failing rule; not an additive decomposition of treatment-effect changes'}
 write_json('cohort_migration.json',report);print(json.dumps(report),flush=True)

def reporting_audit():
 """Read frozen cohorts and old source rules; export aggregate risk tables and entrant reasons."""
 d=pd.read_parquet(OUT/'cohort_private.parquet');old=pd.read_parquet(BASE/'publication_figures/cohort_private.parquet')
 rows=[]
 for ep,h in ENDPOINTS.items():
  for day in ([0,1,2,3,4,5,6,7] if h==168 else ([0,7,14,21,28] if h==672 else [0,30,60,90])):
   for arm,x in [('NS',0),('LR',1)]:rows.append(dict(outcome=ep,arm=arm,day=day,at_risk=int((d.exposure.eq(x)&d[ep+'_time'].ge(day*24)).sum())))
 pd.DataFrame(rows).to_csv(OUT/'risk_tables.csv',index=False)
 entrants=d[~d.stay_id.isin(old.stay_id)][['stay_id','exposure']];assert entrants.stay_id.is_unique
 m=load('canonical_rerun_n02_n07_n09');con=m.connection();con.execute('set threads=1');con.register('entrants',entrants)
 sql=m.N09_CTE+" select p.stay_id,p.ns_ml,a.t0 from prior_fluid p join entrants e using(stay_id) left join action a using(stay_id)"
 con.execute('create temp table entry_old_rules as '+sql)
 renal="""select e.stay_id,e.exposure,o.ns_ml,o.t0,r.old_stage,
 exists(select 1 from b.mimiciv_derived.rrt k where k.stay_id=e.stay_id and k.charttime<=o.t0 and k.dialysis_present=1) old_rrt
 from entrants e left join entry_old_rules o using(stay_id)
 left join lateral(select case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3 else 0 end old_stage
 from b.mimiciv_derived.kdigo_creatinine k where k.stay_id=e.stay_id and k.charttime<=o.t0 order by charttime desc limit 1) r on true"""
 q=con.execute(renal).df();con.close()
 q['reason']=np.select([q.ns_ml.isna(),q.ns_ml.lt(500),q.t0.isna(),q.old_stage.ge(3).fillna(False).to_numpy(dtype=bool),q.old_rrt.to_numpy(dtype=bool)],[
 'No historical trigger','Historical completed-status prior saline below 500 mL','No historical completed bolus of at least 250 mL in window','Historical creatinine stage 3 exclusion at old t0','Historical dialysis-present exclusion'],'Unexplained')
 assert not q.reason.eq('Unexplained').any(),q.reason.value_counts().to_dict()
 agg=q.groupby('reason').agg(n=('reason','size'),LR=('exposure','sum')).reset_index();agg['NS']=agg.n-agg.LR
 agg.to_csv(OUT/'new_entrant_reasons.csv',index=False)
 (OUT/'new_entrant_audit.sql').write_text(sql+';\n'+renal,encoding='utf-8')
 print(agg.to_string(index=False),flush=True)

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['extract','analyse','all','check','bridge','reporting']);args=parser.parse_args()
 with threadpool_limits(limits=1):
  if args.mode=='bridge':extract(audit_only=True)
  if args.mode=='reporting':reporting_audit()
  if args.mode in ('extract','all'):extract()
  if args.mode in ('analyse','all'):analyse()
  if args.mode in ('check','all'):validate()
