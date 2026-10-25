-- demography
select c.stay_id,p.anchor_year,p.anchor_age,p.anchor_year_group,
    year(c.t0)-p.anchor_year treatment_offset,year(c.admittime)-p.anchor_year admission_offset,
    p.anchor_age+year(c.admittime)-p.anchor_year age_admission,
    cast(substr(p.anchor_year_group,1,4) as int)+1+year(c.t0)-p.anchor_year era_midpoint,
    exists(select 1 from b.mimiciv_hosp.admissions a join read_csv_auto('E:/ICUdatabase/MIMIC/mimic-iv/hosp/diagnoses_icd.csv.gz') dx using(hadm_id)
    where a.subject_id=c.subject_id and a.hadm_id<>c.hadm_id and a.dischtime<c.t0
    and (replace(dx.icd_code,'.','') like '585%' or replace(dx.icd_code,'.','') like 'N18%')) prior_CKD
    from r2 c join b.mimiciv_hosp.patients p using(subject_id);

-- sofa_timing
select c.stay_id,z.sofa_time,z.sofa_score,s.starttime,s.endtime,s.hr,s.sofa_24hours,
    (select max(x.endtime) from b.mimiciv_derived.sofa x where x.stay_id=c.stay_id and x.hr between s.hr-23 and s.hr) last_component_hour
    from r2 c join b.mimiciv_derived.sepsis3 z using(stay_id)
    left join b.mimiciv_derived.sofa s on s.stay_id=c.stay_id and s.endtime=z.sofa_time;

-- labs
select c.stay_id,l.labevent_id,l.itemid,l.charttime,l.storetime,greatest(l.charttime,l.storetime) available_time,l.valuenum
    from r2 c join b.mimiciv_hosp.labevents l using(hadm_id)
    where l.itemid in(50912,50971,50820,50882,50813) and l.valuenum is not null
    and l.charttime between c.t0-interval '7d' and c.t0+interval '7d';

-- rrt_derived
select c.stay_id,x.charttime,x.dialysis_present,x.dialysis_active,x.dialysis_type
    from r2 c join b.mimiciv_derived.rrt x using(stay_id) where x.charttime<=c.t0+interval '7d';

-- rrt_chart
select c.stay_id,x.itemid,x.charttime,x.storetime,x.value,x.valuenum
    from r2 c join b.mimiciv_icu.chartevents x using(stay_id)
    where x.itemid in(224135,224139,224144,224145,224146,224149,224150,224151,224152,224153,224154,224191,224270,224404,224406,225183,225323,225436,225441,225725,225740,225776,225802,225803,225805,225806,225807,225809,225810,225951,225952,225953,225954,225955,225956,225958,225959,225961,225963,225965,225976,225977,226118,226457,226499,227124,227290,227357,227438,227525,227536,227638,227639,227640,227753,228004,228005,228006) and x.charttime<=c.t0+interval '7d' and x.value is not null;

-- rrt_procedure
select c.stay_id,x.itemid,x.starttime,x.endtime,x.storetime,x.value,x.statusdescription
    from r2 c join b.mimiciv_icu.procedureevents x using(stay_id)
    where x.itemid in(225441,225802,225803,225805,225809,225955,224270,225436) and x.starttime<=c.t0+interval '7d'
    and x.statusdescription is distinct from 'Rewritten';

-- inputs
select c.stay_id,x.itemid,x.starttime,x.endtime,x.storetime,x.amount,x.rate,x.rateuom,x.ordercategoryname
    from r2 c join b.mimiciv_icu.inputevents x using(stay_id)
    where x.itemid in(225158,225828,221906,222315,227536,227525) and x.starttime<=c.t0+interval '7d'
    and x.statusdescription is distinct from 'Rewritten';

-- urine_stages
select c.stay_id,x.charttime,x.aki_stage_uo from r2 c join b.mimiciv_derived.kdigo_stages x using(stay_id)
    where x.charttime<=c.t0 and x.charttime>c.t0-interval '24h';

-- pressor_rates
select c.stay_id,n.starttime,n.endtime,n.vaso_rate from r2 c join b.mimiciv_derived.norepinephrine n using(stay_id)
    where n.starttime<c.t0 and n.endtime>c.t0-interval '12h';

-- services
select c.stay_id,s.transfertime,s.curr_service from r2 c join read_csv_auto('E:/ICUdatabase/MIMIC/mimic-iv/hosp/services.csv.gz') s using(hadm_id) where s.transfertime<=c.t0;

-- ventilation
select c.stay_id,v.starttime,v.endtime,v.ventilation_status from r2 c join b.mimiciv_derived.ventilation v using(stay_id)
    where v.starttime<=c.t0 and v.endtime>c.t0;

-- entrants
select e.stay_id,o.t0 old_t0 from entrants e join old_entry o using(stay_id)
    left join lateral(select creat,creat_low_past_7day lo7,creat_low_past_48hr lo48 from b.mimiciv_derived.kdigo_creatinine k where k.stay_id=e.stay_id and k.charttime<=o.t0 order by charttime desc limit 1)k on true
    where o.ns_ml>=500 and o.t0 is not null and not coalesce(k.creat>=3*k.lo7 or (k.creat>=4 and k.creat>=k.lo48+.3),false)
    and exists(select 1 from b.mimiciv_derived.rrt x where x.stay_id=e.stay_id and x.charttime<=o.t0 and x.dialysis_present=1);