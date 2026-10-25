"""Bounded second-review source audit. Explicit CLI, read-only databases, private rows stay local."""
from pathlib import Path
import argparse, json, hashlib, re
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
import n09_corrected_estimand as mod
import canonical_rerun_n02_n07_n09 as m

BASE=mod.BASE; DATA=mod.OUT; OUT=DATA/'review_audit'; PRIVATE=OUT/'round2_private'
def export(name,rows): pd.DataFrame(rows).to_csv(OUT/(name+'.csv'),index=False)
def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def support_audit():
    d=pd.read_parquet(DATA/'cohort_private.parquet');con=m.connection();con.execute('set threads=1');con.register('r2',d)
    source=Path('E:/ICUdatabase/MIMIC/mimic-code-main/mimic-iv/concepts/treatment/ventilation.sql').read_text(encoding='utf-8')
    prefix=source.split(', vd0 AS')[0].replace('`physionet-data.', 'b.').replace('`','')
    sql=prefix+" select c.stay_id,v.charttime,v.ventilation_status from r2 c left join lateral(select * from vs where vs.stay_id=c.stay_id and vs.charttime<=c.t0 and vs.charttime>c.t0-interval '14h' and vs.ventilation_status is not null order by charttime desc limit 1)v on true"
    v=con.execute(sql).df();v.to_parquet(PRIVATE/'pre_t0_vent_state.parquet',index=False)
    # Raw initiation evidence supports retrospective baseline state; storage delay is reported separately.
    sql2="""select c.stay_id,x.itemid,x.starttime,x.storetime from r2 c join b.mimiciv_icu.inputevents x using(stay_id)
    where x.itemid in(221906,222315) and x.starttime<=c.t0 and x.endtime>c.t0 and x.statusdescription is distinct from 'Rewritten'"""
    p=con.execute(sql2).df();p.to_parquet(PRIVATE/'pressor_start_evidence.parquet',index=False);con.close()
    v=v.merge(d[['stay_id','t0','ventilation']],on='stay_id');p=p.merge(d[['stay_id','t0']],on='stay_id')
    res={'ventilation_positive':int(d.ventilation.sum()),'positive_vent_without_pre_t0_invasive_classification_14h':int((v.ventilation.eq(1)&~v.ventilation_status.eq('InvasiveVent')).sum()),'latest_pre_t0_vent_missing':int(v.charttime.isna().sum()),'raw_pressor_active_patients':int(p.stay_id.nunique()),'derived_pressor_positive':int(d.pressor.sum()),'derived_pressor_without_raw_active_start':int((d.pressor.eq(1)&~d.stay_id.isin(p.stay_id)).sum()),'pressor_start_after_t0':int(p.starttime.gt(p.t0).sum()),'pressor_active_record_stored_after_t0_patients':int(p.loc[p.storetime.gt(p.t0),'stay_id'].nunique()),'note':'Retrospective pre-treatment state; neither intervals nor delayed storage prove clinician-visible information. Full raw SOFA component availability remains unverified.'}
    (OUT/'round2_support_checks.json').write_text(json.dumps(res,indent=2),encoding='utf-8');(OUT/'round2_support.sql').write_text(sql+';\n'+sql2+';',encoding='utf-8');print(json.dumps(res),flush=True)

def extract():
    PRIVATE.mkdir(exist_ok=True)
    d=pd.read_parquet(DATA/'cohort_private.parquet'); assert d.stay_id.is_unique
    con=m.connection(); con.execute('set threads=1'); con.register('r2',d); sqls=[]
    def q(name,sql):
        sqls.append('-- '+name+'\n'+sql+';');f=con.execute(sql).df();f.to_parquet(PRIVATE/(name+'.parquet'),index=False);return f
    q('demography',"""select c.stay_id,p.anchor_year,p.anchor_age,p.anchor_year_group,
    year(c.t0)-p.anchor_year treatment_offset,year(c.admittime)-p.anchor_year admission_offset,
    p.anchor_age+year(c.admittime)-p.anchor_year age_admission,
    cast(substr(p.anchor_year_group,1,4) as int)+1+year(c.t0)-p.anchor_year era_midpoint,
    exists(select 1 from b.mimiciv_hosp.admissions a join read_csv_auto('E:/ICUdatabase/MIMIC/mimic-iv/hosp/diagnoses_icd.csv.gz') dx using(hadm_id)
    where a.subject_id=c.subject_id and a.hadm_id<>c.hadm_id and a.dischtime<c.t0
    and (replace(dx.icd_code,'.','') like '585%' or replace(dx.icd_code,'.','') like 'N18%')) prior_CKD
    from r2 c join b.mimiciv_hosp.patients p using(subject_id)""")
    q('sofa_timing',"""select c.stay_id,z.sofa_time,z.sofa_score,s.starttime,s.endtime,s.hr,s.sofa_24hours,
    (select max(x.endtime) from b.mimiciv_derived.sofa x where x.stay_id=c.stay_id and x.hr between s.hr-23 and s.hr) last_component_hour
    from r2 c join b.mimiciv_derived.sepsis3 z using(stay_id)
    left join b.mimiciv_derived.sofa s on s.stay_id=c.stay_id and s.endtime=z.sofa_time""")
    q('labs',"""select c.stay_id,l.labevent_id,l.itemid,l.charttime,l.storetime,greatest(l.charttime,l.storetime) available_time,l.valuenum
    from r2 c join b.mimiciv_hosp.labevents l using(hadm_id)
    where l.itemid in(50912,50971,50820,50882,50813) and l.valuenum is not null
    and l.charttime between c.t0-interval '7d' and c.t0+interval '7d'""")
    q('rrt_derived',"""select c.stay_id,x.charttime,x.dialysis_present,x.dialysis_active,x.dialysis_type
    from r2 c join b.mimiciv_derived.rrt x using(stay_id) where x.charttime<=c.t0+interval '7d'""")
    concepts=Path('E:/ICUdatabase/MIMIC/mimic-code-main/mimic-iv/concepts')
    rrt=(concepts/'treatment/rrt.sql').read_text(encoding='utf-8')
    ids=sorted(set(re.findall(r'\b22\d{4}\b',rrt))); idlist=','.join(ids)
    q('rrt_chart',f"""select c.stay_id,x.itemid,x.charttime,x.storetime,x.value,x.valuenum
    from r2 c join b.mimiciv_icu.chartevents x using(stay_id)
    where x.itemid in({idlist}) and x.charttime<=c.t0+interval '7d' and x.value is not null""")
    q('rrt_procedure',f"""select c.stay_id,x.itemid,x.starttime,x.endtime,x.storetime,x.value,x.statusdescription
    from r2 c join b.mimiciv_icu.procedureevents x using(stay_id)
    where x.itemid in({mod.PROCEDURES},224270,225436) and x.starttime<=c.t0+interval '7d'
    and x.statusdescription is distinct from 'Rewritten'""")
    q('inputs',"""select c.stay_id,x.itemid,x.starttime,x.endtime,x.storetime,x.amount,x.rate,x.rateuom,x.ordercategoryname
    from r2 c join b.mimiciv_icu.inputevents x using(stay_id)
    where x.itemid in(225158,225828,221906,222315,227536,227525) and x.starttime<=c.t0+interval '7d'
    and x.statusdescription is distinct from 'Rewritten'""")
    q('urine_stages',"""select c.stay_id,x.charttime,x.aki_stage_uo from r2 c join b.mimiciv_derived.kdigo_stages x using(stay_id)
    where x.charttime<=c.t0 and x.charttime>c.t0-interval '24h'""")
    q('pressor_rates',"""select c.stay_id,n.starttime,n.endtime,n.vaso_rate from r2 c join b.mimiciv_derived.norepinephrine n using(stay_id)
    where n.starttime<c.t0 and n.endtime>c.t0-interval '12h'""")
    q('services',"""select c.stay_id,s.transfertime,s.curr_service from r2 c join read_csv_auto('E:/ICUdatabase/MIMIC/mimic-iv/hosp/services.csv.gz') s using(hadm_id) where s.transfertime<=c.t0""")
    q('ventilation',"""select c.stay_id,v.starttime,v.endtime,v.ventilation_status from r2 c join b.mimiciv_derived.ventilation v using(stay_id)
    where v.starttime<=c.t0 and v.endtime>c.t0""")
    # Preserve the original first-failing-rule definition of the ten dialysis-flag entrants.
    old=pd.read_parquet(BASE/'publication_figures/cohort_private.parquet')
    ent=d[~d.stay_id.isin(old.stay_id)][['stay_id','exposure']];con.register('entrants',ent)
    con.execute('create temp table old_entry as '+m.N09_CTE+' select p.stay_id,p.ns_ml,a.t0 from prior_fluid p join entrants e using(stay_id) left join action a using(stay_id)')
    q('entrants',"""select e.stay_id,o.t0 old_t0 from entrants e join old_entry o using(stay_id)
    left join lateral(select creat,creat_low_past_7day lo7,creat_low_past_48hr lo48 from b.mimiciv_derived.kdigo_creatinine k where k.stay_id=e.stay_id and k.charttime<=o.t0 order by charttime desc limit 1)k on true
    where o.ns_ml>=500 and o.t0 is not null and not coalesce(k.creat>=3*k.lo7 or (k.creat>=4 and k.creat>=k.lo48+.3),false)
    and exists(select 1 from b.mimiciv_derived.rrt x where x.stay_id=e.stay_id and x.charttime<=o.t0 and x.dialysis_present=1)""")
    con.close();(OUT/'round2_extraction.sql').write_text('\n\n'.join(sqls),encoding='utf-8')
    for source in ['treatment/rrt.sql','score/sofa.sql','sepsis/sepsis3.sql','treatment/ventilation.sql']:
        p=concepts/source;(OUT/('source_'+p.name)).write_bytes(p.read_bytes())
    print('ROUND2 SOURCE EXTRACTION COMPLETE',flush=True)

def analyze():
    d=pd.read_parquet(DATA/'cohort_private.parquet').reset_index(drop=True);d,_,_,ps,w=mod.fit(d,m)
    a=d.exposure.to_numpy(int); load=lambda name:pd.read_parquet(PRIVATE/(name+'.parquet'))
    demo=load('demography'); labs=load('labs'); rr=load('rrt_derived'); ce=load('rrt_chart'); pe=load('rrt_procedure');inp=load('inputs');uo=load('urine_stages');ne=load('pressor_rates');svc=load('services');sofa=load('sofa_timing');vent=load('ventilation');entrants=load('entrants')
    entrants=entrants[entrants.stay_id.isin(d.stay_id)]
    f=d.merge(demo.drop(columns='anchor_year_group'),on='stay_id',validate='one_to_one')
    f['age_changed']=f.age.ne(f.age_admission)
    offsets=f.groupby(['exposure','treatment_offset','admission_offset']).size().rename('n').reset_index();export('round2_year_offsets',offsets)
    st=sofa.merge(d[['stay_id','t0']],on='stay_id');assert len(st)==len(d)
    summary={'n':len(d),'age_changed':int(f.age_changed.sum()),'treatment_year_offset_nonzero':int(f.treatment_offset.ne(0).sum()),'sofa_missing_matched_hour':int(st.endtime.isna().sum()),'sofa_endpoint_after_t0':int(st.endtime.gt(st.t0).sum()),'sofa_component_hour_after_t0':int(st.last_component_hour.gt(st.t0).sum()),'full_raw_sofa_storage_time_adjudication':False}
    procids=set(map(int,mod.PROCEDURES.split(',')));flowids=set(map(int,mod.FLOW_ITEMS.split(',')))
    # Extended positive treatment evidence; access/pressure/goal fields alone are not proof of therapy.
    extflow={225806,225807,228004,228005,228006,224153};rows=[];features=[];monitor=[]
    dictionaries=pd.read_csv('E:/ICUdatabase/MIMIC/mimic-iv/icu/d_items.csv.gz').set_index('itemid').label.to_dict()
    for i,r in d.iterrows():
        sid=r.stay_id;t0=r.t0;lg=labs[labs.stay_id.eq(sid)];cg=ce[ce.stay_id.eq(sid)];pg=pe[pe.stay_id.eq(sid)];ig=inp[inp.stay_id.eq(sid)];ug=uo[uo.stay_id.eq(sid)]
        cr=mod.creat_stages(lg[lg.itemid.eq(50912)]); pre=cr[cr.charttime.le(t0)]; late=pre[pre.available_time.gt(t0)&pre.stage_new.eq(3)]
        cpre=cg[cg.charttime.le(t0)];ppre=pg[pg.starttime.le(t0)];ipre=ig[ig.starttime.le(t0)]
        active_ext=cpre[cpre.itemid.isin(extflow)&cpre.valuenum.gt(0)]
        active_primary=cpre[cpre.itemid.isin(flowids)&cpre.valuenum.gt(0)]
        active_proc=ppre[ppre.itemid.isin(procids)&ppre.value.notna()]
        adjunct=ipre[ipre.itemid.isin([227536,227525])&ipre.amount.gt(0)]
        ktime=r.KRT_stamp
        near_proc=pg[pg.itemid.isin(procids)&pg.value.notna()&pg.starttime.sub(ktime).abs().le(pd.Timedelta(hours=6))]
        near_flow=cg[cg.itemid.isin(flowids)&cg.valuenum.gt(0)&cg.charttime.sub(ktime).abs().le(pd.Timedelta(hours=6))]
        before_first=cg[cg.charttime.lt(ktime)] if pd.notna(ktime) else cg.iloc[:0]
        prior_filter=bool((ppre.itemid.eq(225436)).any())
        extended=bool(len(active_ext) or len(active_primary) or len(active_proc))
        classification='Prior active-treatment evidence' if extended else ('Procedure and positive-flow corroboration; initiation unresolved' if len(near_proc) and len(near_flow) else 'Single-source treatment record; initiation unresolved')
        labels=sorted(set(dictionaries.get(x,str(x)) for x in cpre.itemid)|set(dictionaries.get(x,str(x)) for x in ppre.itemid))
        rows.append(dict(stay_id=sid,arm='LR' if r.exposure else 'NS',KRT_event=r.KRT_event==1,early6=r.primary_event==1 and r.primary_time<=6,entrant=sid in set(entrants.stay_id),late_pre3=len(late)>0,pre_retrospective3=pre.stage_new.eq(3).any(),extended_prior_active=extended,prior_filter_change=prior_filter,prior_CRRT_adjunct=len(adjunct)>0,prior_documentation='; '.join(labels),classification=classification,procedure_flow_corroboration=bool(len(near_proc) and len(near_flow)),first_component='KRT' if pd.notna(r.KRT_stamp) and (pd.isna(r.creatinine3_stamp) or r.KRT_stamp<r.creatinine3_stamp) else 'Creatinine',pre_creat_last=pre.valuenum.iloc[-1] if len(pre) else np.nan,pre_creat_change=float(pre.valuenum.iloc[-1]-pre.valuenum.iloc[0]) if len(pre)>1 else np.nan,pre_uo_stage_max=ug.aki_stage_uo.max(),primary_hour=r.primary_time,indication='Unknown; no independent clinical adjudication'))
        known=lg[lg.available_time.le(t0)&lg.charttime.gt(t0-pd.Timedelta(hours=24))].sort_values(['charttime','available_time'])
        ph=known[known.itemid.eq(50820)].valuenum
        ng=ne[ne.stay_id.eq(sid)];last=ng[ng.starttime.lt(t0)&ng.endtime.gt(t0-pd.Timedelta(hours=6))].vaso_rate.max();prior=ng[ng.starttime.lt(t0-pd.Timedelta(hours=6))&ng.endtime.gt(t0-pd.Timedelta(hours=12))].vaso_rate.max()
        sg=svc[svc.stay_id.eq(sid)].sort_values('transfertime');service=sg.curr_service.iloc[-1] if len(sg) else 'Unknown'
        gap=ig[ig.starttime.ge(r.trigger_time)&ig.endtime.le(t0)&ig.amount.gt(0)&ig.itemid.isin([225158,225828])]
        features.append(dict(stay_id=sid,ph_pre=ph.iloc[-1] if len(ph) else np.nan,ne_max6=0 if pd.isna(last) else last,ne_change6=(0 if pd.isna(last) else last)-(0 if pd.isna(prior) else prior),surgical_service=service in ['CSURG','NSURG','ORTHO','PSURG','SURG','TSURG','VSURG'],gap_NS_ml=gap.loc[gap.itemid.eq(225158),'amount'].sum(),gap_LR_ml=gap.loc[gap.itemid.eq(225828),'amount'].sum()))
        end=min(t0+pd.Timedelta(days=7),r.outtime,r.deathtime if pd.notna(r.deathtime) else r.outtime)
        n=int(cr.charttime.gt(t0).mul(cr.charttime.le(end)).sum());days=(end-t0).total_seconds()/86400
        monitor.append(dict(stay_id=sid,arm='LR' if r.exposure else 'NS',tests=n,patient_days=days,weight=w[i]))
    cases=pd.DataFrame(rows);cases.to_parquet(PRIVATE/'case_review.parquet',index=False)
    feature=pd.DataFrame(features);f=f.merge(feature,on='stay_id');f.to_parquet(PRIVATE/'augmented_cohort.parquet',index=False)
    relevant=cases.KRT_event|cases.early6|cases.entrant
    summary.update({'reviewed_unique_cases':int(relevant.sum()),'KRT_cases':int(cases.KRT_event.sum()),'early_cases':int(cases.early6.sum()),'entrants':len(entrants),'early_with_late_pre3':int((cases.early6&cases.late_pre3).sum()),'prior_active_extended_all':int(cases.extended_prior_active.sum()),'prior_active_extended_KRT':int((cases.extended_prior_active&cases.KRT_event).sum())})
    report=[]
    for group,mask in [('All KRT events',cases.KRT_event),('Early primary events',cases.early6),('Historical dialysis-flag entrants',cases.entrant)]:
        for arm in ['NS','LR']:
            g=cases[mask&cases.arm.eq(arm)]
            for field in ['KRT_event','early6','late_pre3','pre_retrospective3','extended_prior_active','prior_filter_change','prior_CRRT_adjunct','procedure_flow_corroboration']:
                report.append(dict(group=group,arm=arm,measure=field,n=int(g[field].sum()),denominator=len(g)))
            report.append(dict(group=group,arm=arm,measure='KRT first component',n=int(g.first_component.eq('KRT').sum()),denominator=len(g)))
    export('round2_case_summary',report)
    labels=[]
    for arm in ['NS','LR']:
        g=cases[cases.entrant&cases.arm.eq(arm)]
        for label in sorted(set('; '.join(g.prior_documentation).split('; '))):
            if label: labels.append(dict(arm=arm,source_label=label,patients=int(g.prior_documentation.apply(lambda s:label in s.split('; ')).sum())))
    export('round2_entrant_sources',labels)
    # Deidentified source-feature summaries only; individual review packets stay private.
    er=cases[cases.early6];export('round2_early_summary',[dict(measure='Pre-start last creatinine median',value=er.pre_creat_last.median()),dict(measure='KRT first',value=er.first_component.eq('KRT').sum()),dict(measure='Pre-start late-result stage3',value=er.late_pre3.sum()),dict(measure='Pre-start urine stage unknown',value=er.pre_uo_stage_max.isna().sum()),dict(measure='Prior24h urine stage >=2',value=er.pre_uo_stage_max.ge(2).sum())])
    mr=pd.DataFrame(monitor);obs=[]
    for arm in ['NS','LR']:
        g=mr[mr.arm.eq(arm)];obs.append(dict(arm=arm,tests=int(g.tests.sum()),patient_days=g.patient_days.sum(),tests_per_day=g.tests.sum()/g.patient_days.sum(),weighted_tests_per_day=np.sum(g.tests*g.weight)/np.sum(g.patient_days*g.weight),patients_no_test=int(g.tests.eq(0).sum())))
    export('round2_monitoring',obs)
    # Same original weights, four mutually exclusive first-event states.
    typ=np.where(d.primary_event.eq(1),1,np.where(d.primary_event.eq(2),np.where(d.deathtime.notna()&d.deathtime.le(d.outtime),2,3),0));curves=[];states=[]
    for arm in [0,1]:
        ix=a==arm;ww=w[ix];ev=typ[ix];tt=d.loc[ix,'primary_time'].to_numpy()
        for code,label in [(1,'Kidney event'),(2,'Death before kidney event'),(3,'ICU exit before kidney event'),(0,'Event-free through day 7')]:
            risk=np.average(ev==code,weights=ww);states.append(dict(arm='LR' if arm else 'NS',state=label,n=int((ev==code).sum()),weighted_percent=100*risk))
            if code:
                _,_,cv=mod.cif_summary(tt,np.where(ev==code,1,np.where(ev>0,2,0)),ww,168,np.arange(0,169,3))
                curves.extend(dict(arm='LR' if arm else 'NS',state=label,hour=h,risk=v) for h,v in zip(np.arange(0,169,3),cv))
    export('round2_competing_states',states);export('round2_competing_curves',curves)
    dist=[];ecdf=[]
    for col in ['lr_ml','creat','urine6_ml','net24_ml']:
        for arm in [0,1]:
            ix=(a==arm)&d[col].notna().to_numpy();v=d.loc[ix,col].to_numpy();ww=w[ix];o=np.argsort(v);v=v[o];ww=ww[o];cum=np.cumsum(ww)/ww.sum()
            for prob in [0,.25,.5,.75,.95,1]:dist.append(dict(variable=col,arm='LR' if arm else 'NS',quantile=prob,value=v[min(np.searchsorted(cum,prob),len(v)-1)]))
            for val in np.unique(np.quantile(v,np.linspace(0,1,101))):ecdf.append(dict(variable=col,arm='LR' if arm else 'NS',value=val,cdf=np.sum(ww[v<=val])/ww.sum()))
    export('round2_weighted_quantiles',dist);export('round2_weighted_ecdf',ecdf)
    clinical=[]
    for col in ['prior_CKD','potassium','ph_pre','ne_max6','ne_change6','surgical_service','gap_NS_ml','gap_LR_ml']:
        for arm in [0,1]:
            v=f.loc[f.exposure.eq(arm),col].astype(float);clinical.append(dict(variable=col,arm='LR' if arm else 'NS',observed=int(v.notna().sum()),median=v.median(),q25=v.quantile(.25),q75=v.quantile(.75),positive=int(v.gt(0).sum())))
    export('round2_clinical_context',clinical)
    base_rd=np.average(d.loc[a==1,'primary_event'].eq(1),weights=w[a==1])-np.average(d.loc[a==0,'primary_event'].eq(1),weights=w[a==0]);influence=[]
    for i in range(len(d)):
        dd,_,_,_,ww=mod.fit(d.drop(index=i).reset_index(drop=True),m);aa=dd.exposure.to_numpy();rd=np.average(dd.loc[aa==1,'primary_event'].eq(1),weights=ww[aa==1])-np.average(dd.loc[aa==0,'primary_event'].eq(1),weights=ww[aa==0]);influence.append(dict(arm='LR' if a[i] else 'NS',early=bool(cases.iloc[i].early6),KRT=bool(cases.iloc[i].KRT_event),rd=rd,change=rd-base_rd))
    inf=pd.DataFrame(influence); export('round2_influence_summary',[dict(group=label,n=len(g),rd_min=g.rd.min(),rd_max=g.rd.max(),max_abs_change=g.change.abs().max()) for label,g in [('All patients',inf),('Early events',inf[inf.early]),('KRT events',inf[inf.KRT])]])
    # Deterministic ambiguity envelope for KRT-only cases: explicitly not another causal CI.
    ambiguity=[]
    for arm in [0,1]:
        ix=a==arm;konly=d.primary_event.eq(1)&d.KRT_event.eq(1)&~d.creatinine3_event.eq(1)
        ambiguity.append(np.sum(w[ix&konly])/np.sum(w[ix]))
    summary.update({'original_rd':base_rd,'ambiguous_KRT_only_envelope_low':base_rd-ambiguity[1],'ambiguous_KRT_only_envelope_high':base_rd+ambiguity[0],'leave_one_out_rd_min':inf.rd.min(),'leave_one_out_rd_max':inf.rd.max()})
    # Fixed clinical/time sensitivity; keep all endpoints for coherent outcome reporting.
    original_num=mod.NUM[:];original_cat=mod.CAT[:]
    mod.NUM=[x for x in original_num if x!='age']+['age_admission','era_midpoint','prior_CKD','potassium','ph_pre','ne_max6','ne_change6','surgical_service']
    mod.CAT=[x for x in original_cat if x!='anchor_year_group'];mod.PROTOCOL['bootstrap']=500
    result,cur,bal,diag,hist=mod.estimate(f,m,'Time and clinical-context adjustment')
    export('round2_adjusted_outcomes',result);export('round2_adjusted_balance',bal)
    diags=[diag];mod.NUM=original_num;mod.CAT=original_cat
    if cases.extended_prior_active.any():
        keep=~d.stay_id.isin(cases.loc[cases.extended_prior_active,'stay_id'])
        res,_,_,dg,_=mod.estimate(d[keep].reset_index(drop=True),m,'Exclude extended pre-start active KRT evidence',endpoints={'primary':168.});export('round2_classification_sensitivity',res);diags.append(dg)
    summary['diagnostics']=diags;summary['protocol_sha256']=digest(OUT/'round2_plan.json');summary['main_outcomes_sha256']=digest(DATA/'outcomes.csv');summary['clinical_adjudication']=False
    (OUT/'round2_checks.json').write_text(json.dumps(summary,indent=2,default=str),encoding='utf-8');print(json.dumps(summary,indent=2,default=str),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['extract','analyze','support','all']);args=p.parse_args()
    with threadpool_limits(limits=1):
        if args.mode in ['extract','all']:extract()
        if args.mode in ['analyze','all']:analyze()
        if args.mode in ['support','all']:support_audit()
