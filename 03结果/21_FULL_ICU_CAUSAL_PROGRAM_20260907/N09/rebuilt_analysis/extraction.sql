create temp table rebuild_risk as 
with chloride_results as (
  select p.subject_id,p.hadm_id,p.stay_id,p.age,p.gender,p.first_careunit,
         p.anchor_year_group,p.intime,p.outtime,
         l.labevent_id,greatest(l.charttime,l.storetime) trigger_time,
         l.valuenum chloride,
         (select min(x.valuenum) from b.mimiciv_hosp.labevents x
           where x.hadm_id=l.hadm_id and x.itemid=50902 and x.valuenum between 70 and 150
             and x.charttime<l.charttime and x.charttime>=l.charttime-interval '48h'
             and greatest(x.charttime,x.storetime)<=greatest(l.charttime,l.storetime)) prior_chloride,
         z.sofa_score,greatest(z.antibiotic_time,z.culture_time,z.sofa_time) sepsis_known
  from s.main.people p join b.mimiciv_derived.sepsis3 z using(stay_id,subject_id)
  join b.mimiciv_hosp.labevents l using(hadm_id)
  where p.age>=18 and z.sepsis3 and l.itemid=50902 and l.valuenum between 70 and 150
    and l.storetime is not null and greatest(l.charttime,l.storetime)>=p.intime
    and greatest(l.charttime,l.storetime)<least(p.outtime,p.intime+interval '72h')
), trigger as (
  select * from chloride_results where trigger_time>=sepsis_known
    and (chloride>=110 or chloride-prior_chloride>=5)
  qualify row_number() over(partition by stay_id order by trigger_time,labevent_id)=1
), prior_fluid as (
  select t.*,coalesce(f.ns_ml,0) ns_ml,coalesce(f.lr_ml,0) lr_ml
  from trigger t left join lateral(
    select sum(amount) filter(where itemid=225158) ns_ml,
           sum(amount) filter(where itemid=225828) lr_ml
    from b.mimiciv_icu.inputevents i where i.stay_id=t.stay_id
      and i.endtime<=t.trigger_time and i.endtime>t.trigger_time-interval '24h'
      and i.amount>0 and i.statusdescription is distinct from 'Rewritten'
      and i.ordercategoryname in('03-IV Fluid Bolus','02-Fluids (Crystalloids)')) f on true
)
select * from prior_fluid where ns_ml>=500

;

with chloride_results as (
  select p.subject_id,p.hadm_id,p.stay_id,p.age,p.gender,p.first_careunit,
         p.anchor_year_group,p.intime,p.outtime,
         l.labevent_id,greatest(l.charttime,l.storetime) trigger_time,
         l.valuenum chloride,
         (select min(x.valuenum) from b.mimiciv_hosp.labevents x
           where x.hadm_id=l.hadm_id and x.itemid=50902 and x.valuenum between 70 and 150
             and x.charttime<l.charttime and x.charttime>=l.charttime-interval '48h'
             and greatest(x.charttime,x.storetime)<=greatest(l.charttime,l.storetime)) prior_chloride,
         z.sofa_score,greatest(z.antibiotic_time,z.culture_time,z.sofa_time) sepsis_known
  from s.main.people p join b.mimiciv_derived.sepsis3 z using(stay_id,subject_id)
  join b.mimiciv_hosp.labevents l using(hadm_id)
  where p.age>=18 and z.sepsis3 and l.itemid=50902 and l.valuenum between 70 and 150
    and l.storetime is not null and greatest(l.charttime,l.storetime)>=p.intime
    and greatest(l.charttime,l.storetime)<least(p.outtime,p.intime+interval '72h')
), trigger as (
  select * from chloride_results where trigger_time>=sepsis_known
    and (chloride>=110 or chloride-prior_chloride>=5)
  qualify row_number() over(partition by stay_id order by trigger_time,labevent_id)=1
), prior_fluid as (
  select t.*,coalesce(f.ns_ml,0) ns_ml,coalesce(f.lr_ml,0) lr_ml
  from trigger t left join lateral(
    select sum(amount) filter(where itemid=225158) ns_ml,
           sum(amount) filter(where itemid=225828) lr_ml
    from b.mimiciv_icu.inputevents i where i.stay_id=t.stay_id
      and i.endtime<=t.trigger_time and i.endtime>t.trigger_time-interval '24h'
      and i.amount>0 and i.statusdescription is distinct from 'Rewritten'
      and i.ordercategoryname in('03-IV Fluid Bolus','02-Fluids (Crystalloids)')) f on true
)
select 'Adult first ICU source' step,count(*) n from s.main.people where age>=18 union all select 'Sepsis source',count(distinct p.stay_id) from s.main.people p join b.mimiciv_derived.sepsis3 z using(stay_id,subject_id) where p.age>=18 and z.sepsis3 union all select 'Eligible chloride measurement',count(distinct stay_id) from chloride_results union all select 'First chloride trigger',count(*) from trigger union all select 'Prior saline at least 500mL',count(*) from prior_fluid where ns_ml>=500
;
create temp table rebuild_bolus as select i.* from b.mimiciv_icu.inputevents i semi join rebuild_risk r using(stay_id)
 where itemid in(225158,225828) and ordercategoryname='03-IV Fluid Bolus' and statusdescription is distinct from 'Rewritten'
 qualify row_number() over(partition by stay_id,itemid,coalesce(linkorderid,orderid) order by starttime,orderid)=1
;
create temp table rebuild_actions as select r.*,i.itemid,i.starttime t0,i.endtime action_end,i.amount action_amount,i.originalamount,i.rate,i.storetime action_storetime,i.statusdescription,
 (select count(distinct x.itemid) from rebuild_bolus x where x.stay_id=r.stay_id and x.starttime=i.starttime) simultaneous_fluids
 from rebuild_risk r join lateral(select * from rebuild_bolus i where i.stay_id=r.stay_id and i.starttime>=r.trigger_time and i.starttime<r.trigger_time+interval '4h' order by starttime,itemid,orderid limit 1)i on true
;
select count(*) from rebuild_actions
;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 0
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 32
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 64
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 96
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 128
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 160
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 192
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 224
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 256
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 288
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 320
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 352
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 384
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 416
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 448
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 480
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 512
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 544
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 576
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 608
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 640
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 672
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 704
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 736
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 768
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 800
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create or replace temp table rebuild_batch as select * from rebuild_actions order by subject_id,stay_id limit 32 offset 832
;

select a.*,case when itemid=225828 then 1 else 0 end exposure,
  r.creat,r.low7,r.low48,r.aki_stage_creat,
  lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
  vit.map_mean,vit.map_min,vit.heart_rate,
  case when exists(select 1 from b.mimiciv_derived.norepinephrine n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
         or exists(select 1 from b.mimiciv_derived.vasopressin n
                    where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0)
       then 1 else 0 end pressor,
  case when exists(select 1 from b.mimiciv_derived.ventilation v
                    where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0
                      and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation
from rebuild_batch a
left join lateral(
  select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
    case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3
         when creat>=2*creat_low_past_7day then 2
         when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
  from b.mimiciv_derived.kdigo_creatinine k
  where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1) r on true
left join lateral(
  select arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50983) sodium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50882) bicarbonate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50813) lactate,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=50971) potassium,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51222) hemoglobin,
         arg_max(valuenum,(charttime,storetime,labevent_id)) filter(where itemid=51265) platelet
  from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
    and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h') lab on true
left join lateral(
  select avg(mbp) map_mean,min(mbp) map_min,avg(heart_rate) heart_rate
  from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id
    and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h') vit on true

;
create temp table rebuild_transfers as select t.* from read_csv_auto('E:/ICUdatabase/MIMIC/mimic-iv/hosp/transfers.csv.gz') t semi join rebuild_people c using(hadm_id)
;
select c.stay_id,a.admittime,a.dischtime,a.deathtime,p.dod,
 (select arg_max(careunit,intime) from rebuild_transfers t where t.hadm_id=c.hadm_id and t.intime<=c.t0 and (t.outtime>c.t0 or t.outtime is null)) decision_careunit,
 exists(select 1 from b.mimiciv_hosp.admissions old join read_csv_auto('E:/ICUdatabase/MIMIC/mimic-iv/hosp/diagnoses_icd.csv.gz') dx on dx.hadm_id=old.hadm_id where old.subject_id=c.subject_id and old.dischtime<=c.t0 and old.hadm_id<>c.hadm_id and replace(upper(dx.icd_code),'.','') in('5856','V451','N186','Z992')) known_prior_ESRD
 from rebuild_people c join b.mimiciv_hosp.admissions a using(hadm_id) join b.mimiciv_hosp.patients p on p.subject_id=c.subject_id
;
select c.stay_id,l.labevent_id,l.charttime,greatest(l.charttime,l.storetime) available_time,l.valuenum from rebuild_people c join b.mimiciv_hosp.labevents l on l.hadm_id=c.hadm_id where l.itemid=50912 and l.valuenum>0 and l.charttime<=c.t0+interval '7d'
;
select c.stay_id,k.charttime,k.aki_stage_uo,u.uo_tm_24hr,u.uo_rt_24hr from rebuild_people c join b.mimiciv_derived.kdigo_stages k on k.stay_id=c.stay_id left join b.mimiciv_derived.kdigo_uo u on u.stay_id=k.stay_id and u.charttime=k.charttime where k.charttime<=c.t0+interval '7d' and k.aki_stage_uo is not null
;
select c.stay_id,p.starttime event_time,'procedure' evidence,p.itemid from rebuild_people c join b.mimiciv_icu.procedureevents p on p.stay_id=c.stay_id where p.itemid in(225441,225802,225803,225805,225809,225955) and p.value is not null and p.statusdescription is distinct from 'Rewritten' and p.starttime<=c.t0+interval '7d'
 union all select c.stay_id,p.charttime,'positive flow or output',p.itemid from rebuild_people c join b.mimiciv_icu.chartevents p on p.stay_id=c.stay_id where p.itemid in(224144,224154,226499,226457,224191) and p.valuenum>0 and p.charttime<=c.t0+interval '7d'
;
select statusdescription,itemid,count(*) n,median(action_amount) median_amount,min(action_amount) min_amount,median(datediff('second',t0,action_end)) encoded_seconds from rebuild_actions group by 1,2