"""Canonical MIMIC-IV rerun for N02, N07 and N09.

The script writes only six aggregate artefacts per candidate. Patient-level
tables stay in memory and raw/cached databases are opened read-only.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


ROOT = Path(__file__).resolve().parents[2]
PROGRAM = ROOT / "03结果/21_FULL_ICU_CAUSAL_PROGRAM_20260907"
BASE_DB = ROOT / "03结果/14_SA_AKI_HYDROCORTISONE_20260906/sa_aki.duckdb"
PROGRAM_DB = PROGRAM / "full_program_private.duckdb"
ABX_DB = PROGRAM / "AKI_antibiotic_night_screen/source_private.duckdb"
BOOTSTRAP = 500
SEED = 20260908


def connection(include_abx: bool = False) -> duckdb.DuckDBPyConnection:
    con = duckdb.connect()
    con.execute("set threads=4")
    con.execute("set memory_limit='8GB'")
    con.execute(f"attach '{BASE_DB.as_posix()}' as b(read_only)")
    con.execute(f"attach '{PROGRAM_DB.as_posix()}' as s(read_only)")
    if include_abx:
        con.execute(f"attach '{ABX_DB.as_posix()}' as a(read_only)")
    return con


def outdir(candidate: str) -> Path:
    path = PROGRAM / candidate / "canonical_rerun_20260908"
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path: Path, obj: object) -> None:
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, default=str), encoding="utf-8")


def weighted_aj(time: np.ndarray, event: np.ndarray, weight: np.ndarray, horizon: float) -> float:
    """Weighted cumulative incidence for event 1; event 2 is competing."""
    times = np.unique(time[(event > 0) & (time <= horizon)])
    survival = 1.0
    cif = 0.0
    for current in times:
        at_risk = weight[time >= current].sum()
        if at_risk <= 0:
            continue
        h1 = weight[(time == current) & (event == 1)].sum() / at_risk
        h2 = weight[(time == current) & (event == 2)].sum() / at_risk
        cif += survival * h1
        survival *= max(0.0, 1.0 - h1 - h2)
    return float(cif)


def ato_design(df: pd.DataFrame, numeric: list[str], categorical: list[str]):
    transformer = ColumnTransformer(
        [
            (
                "num",
                make_pipeline(SimpleImputer(strategy="median", add_indicator=True), StandardScaler()),
                numeric,
            ),
            (
                "cat",
                make_pipeline(
                    SimpleImputer(strategy="most_frequent"),
                    OneHotEncoder(handle_unknown="ignore"),
                ),
                categorical,
            ),
        ]
    )
    return transformer


def ato_fit(df: pd.DataFrame, numeric: list[str], categorical: list[str], mult=None):
    if mult is None:
        mult = np.ones(len(df), dtype=int)
    expanded = df.loc[df.index.repeat(mult)].reset_index(drop=True)
    a = expanded.exposure.to_numpy(int)
    if a.sum() < 5 or (1 - a).sum() < 5:
        raise ValueError("insufficient treatment support in resample")
    transformer = ato_design(expanded, numeric, categorical)
    x = transformer.fit_transform(expanded[numeric + categorical])
    model = LogisticRegression(C=0.5, max_iter=2000)
    model.fit(x, a)
    ps = np.clip(model.predict_proba(x)[:, 1], 1e-4, 1 - 1e-4)
    w = np.where(a == 1, 1 - ps, ps)
    return expanded, transformer, x, ps, w


def ato_effects(df: pd.DataFrame, numeric: list[str], categorical: list[str], mult=None):
    expanded, _, _, _, w = ato_fit(df, numeric, categorical, mult)
    a = expanded.exposure.to_numpy(int)
    result = {}
    for outcome, horizon in [
        ("aki3_krt_7d", 168.0),
        ("death_7d", 168.0),
        ("death_28d", 672.0),
        ("death_90d", 2160.0),
    ]:
        time = expanded[f"{outcome}_time"].to_numpy(float)
        event = expanded[f"{outcome}_event"].to_numpy(int)
        risks = [weighted_aj(time[a == arm], event[a == arm], w[a == arm], horizon) for arm in (0, 1)]
        result[outcome] = {
            "control_risk": risks[0],
            "intervention_risk": risks[1],
            "rd": risks[1] - risks[0],
            "rr": risks[1] / risks[0] if risks[0] > 0 else np.nan,
        }
    return result


def ato_diagnostics(df: pd.DataFrame, numeric: list[str], categorical: list[str]):
    expanded, transformer, x, ps, w = ato_fit(df, numeric, categorical)
    a = expanded.exposure.to_numpy(int)
    dense = x.toarray() if hasattr(x, "toarray") else x
    names = transformer.get_feature_names_out()
    rows = []
    for j, name in enumerate(names):
        values = dense[:, j]
        denominator = np.sqrt((values[a == 0].var() + values[a == 1].var()) / 2)
        unweighted = 0.0 if denominator == 0 else abs(values[a == 1].mean() - values[a == 0].mean()) / denominator
        weighted = 0.0 if denominator == 0 else abs(
            np.average(values[a == 1], weights=w[a == 1]) - np.average(values[a == 0], weights=w[a == 0])
        ) / denominator
        rows.append({"covariate": str(name), "unweighted_abs_smd": unweighted, "overlap_abs_smd": weighted})
    balance = pd.DataFrame(rows).sort_values("overlap_abs_smd", ascending=False)
    ess = {
        str(arm): float(w[a == arm].sum() ** 2 / np.square(w[a == arm]).sum()) for arm in (0, 1)
    }
    return balance, {
        "ess_control": ess["0"],
        "ess_intervention": ess["1"],
        "max_abs_smd": float(balance.overlap_abs_smd.max()),
        "ps_quantiles_control": np.quantile(ps[a == 0], [0, 0.1, 0.5, 0.9, 1]).tolist(),
        "ps_quantiles_intervention": np.quantile(ps[a == 1], [0, 0.1, 0.5, 0.9, 1]).tolist(),
    }


def bootstrap_ato(df, numeric, categorical):
    rng = np.random.default_rng(SEED)
    collected = []
    errors = []
    n = len(df)
    for replicate in range(BOOTSTRAP):
        mult = np.bincount(rng.integers(0, n, n), minlength=n)
        try:
            effects = ato_effects(df, numeric, categorical, mult)
            collected.append({"replicate": replicate, **{f"{o}_{k}": v for o, d in effects.items() for k, v in d.items()}})
        except Exception as exc:
            errors.append({"replicate": replicate, "error": repr(exc)})
    return pd.DataFrame(collected), errors


def event_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Construct kidney competing-risk and mortality times from SQL timestamps."""
    for column in ["t0", "primary_time", "competing_time", "dod"]:
        df[column] = pd.to_datetime(df[column])
    primary_h = (df.primary_time - df.t0).dt.total_seconds() / 3600
    competing_h = (df.competing_time - df.t0).dt.total_seconds() / 3600
    first_h = pd.concat([primary_h, competing_h], axis=1).min(axis=1, skipna=True).clip(upper=168).fillna(168)
    primary_first = primary_h.notna() & (primary_h <= competing_h.fillna(np.inf)) & (primary_h <= 168)
    competing_first = competing_h.notna() & (competing_h < primary_h.fillna(np.inf)) & (competing_h <= 168)
    df["aki3_krt_7d_time"] = first_h
    df["aki3_krt_7d_event"] = np.select([primary_first, competing_first], [1, 2], default=0)
    death_h = (df.dod.dt.normalize() - df.t0.dt.normalize()).dt.total_seconds() / 3600
    for label, horizon in [("death_7d", 168), ("death_28d", 672), ("death_90d", 2160)]:
        valid = death_h.notna() & (death_h >= 0) & (death_h <= horizon)
        df[f"{label}_time"] = np.where(valid, death_h, horizon)
        df[f"{label}_event"] = valid.astype(int)
    return df


N09_CTE = r"""
with chloride_results as (
  select p.*, l.labevent_id, greatest(l.charttime,l.storetime) trigger_time,
         l.valuenum chloride,
         (select min(x.valuenum) from b.mimiciv_hosp.labevents x
           where x.hadm_id=l.hadm_id and x.itemid=50902 and x.valuenum between 70 and 150
             and x.charttime<l.charttime and x.charttime>=l.charttime-interval '48h'
             and greatest(x.charttime,x.storetime)<=greatest(l.charttime,l.storetime)) prior_chloride,
         z.sofa_score, greatest(z.antibiotic_time,z.culture_time,z.sofa_time) sepsis_known
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
  select t.*,coalesce(f.ns_ml,0)ns_ml,coalesce(f.lr_ml,0)lr_ml
  from trigger t left join lateral(
    select sum(amount)filter(where itemid=225158)ns_ml,sum(amount)filter(where itemid=225828)lr_ml
    from b.mimiciv_icu.inputevents i where i.stay_id=t.stay_id
      and i.endtime<=t.trigger_time and i.endtime>t.trigger_time-interval '24h'
      and i.amount>0 and i.statusdescription='FinishedRunning'
      and i.ordercategoryname in('03-IV Fluid Bolus','02-Fluids (Crystalloids)'))f on true
), action as (
  select p.*,a.itemid,a.starttime t0,a.amount action_amount
  from prior_fluid p join lateral(
    select * from b.mimiciv_icu.inputevents i where i.stay_id=p.stay_id
      and i.starttime>=p.trigger_time and i.starttime<p.trigger_time+interval '4h'
      and i.itemid in(225158,225828) and i.amount>=250
      and i.statusdescription='FinishedRunning' and i.ordercategoryname='03-IV Fluid Bolus'
    order by i.starttime limit 1)a on true where p.ns_ml>=500
), features as (
  select a.*,case when itemid=225828 then 1 else 0 end exposure,
    r.creat,r.low7,r.low48,r.aki_stage_creat,
    lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
    vit.map_mean,vit.map_min,vit.heart_rate,
    case when exists(select 1 from b.mimiciv_derived.norepinephrine n where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0) or exists(select 1 from b.mimiciv_derived.vasopressin n where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0) then 1 else 0 end pressor,
    case when exists(select 1 from b.mimiciv_derived.ventilation v where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0 and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation,
    (select min(k.charttime) from b.mimiciv_derived.kdigo_stages k where k.stay_id=a.stay_id and k.charttime>a.t0 and k.charttime<=a.t0+interval '7d' and k.aki_stage>=3) severe_time,
    (select min(k.charttime) from b.mimiciv_derived.rrt k where k.stay_id=a.stay_id and k.charttime>a.t0 and k.charttime<=a.t0+interval '7d' and k.dialysis_present=1) krt_time,
    p.dod
  from action a join b.mimiciv_hosp.patients p using(subject_id)
  left join lateral(select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
      case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3 when creat>=2*creat_low_past_7day then 2 when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
    from b.mimiciv_derived.kdigo_creatinine k where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1)r on true
  left join lateral(select arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=50983)sodium,
      arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=50882)bicarbonate,
      arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=50813)lactate,
      arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=50971)potassium,
      arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=51222)hemoglobin,
      arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=51265)platelet
    from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null
      and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h')lab on true
  left join lateral(select avg(mbp)map_mean,min(mbp)map_min,avg(heart_rate)heart_rate
    from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h')vit on true
)
"""


def n09_cohort():
    con = connection()
    con.execute("set threads=1")
    con.execute("set max_temp_directory_size='12GB'")
    # Materializing the 764 action rows prevents DuckDB from repeatedly expanding
    # the correlated trigger query while joining downstream laboratory histories.
    con.execute("create temp table n09_action as " + N09_CTE + " select * from action")
    query = r"""select a.*,case when itemid=225828 then 1 else 0 end exposure,
      r.creat,r.low7,r.low48,r.aki_stage_creat,
      lab.sodium,lab.bicarbonate,lab.lactate,lab.potassium,lab.hemoglobin,lab.platelet,
      vit.map_mean,vit.map_min,vit.heart_rate,
      case when exists(select 1 from b.mimiciv_derived.norepinephrine n where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0) or exists(select 1 from b.mimiciv_derived.vasopressin n where n.stay_id=a.stay_id and n.starttime<=a.t0 and n.endtime>a.t0) then 1 else 0 end pressor,
      case when exists(select 1 from b.mimiciv_derived.ventilation v where v.stay_id=a.stay_id and v.starttime<=a.t0 and v.endtime>a.t0 and v.ventilation_status='InvasiveVent') then 1 else 0 end ventilation,
      least((select min(k.charttime) from b.mimiciv_derived.kdigo_stages k where k.stay_id=a.stay_id and k.charttime>a.t0 and k.charttime<=a.t0+interval '7d' and k.aki_stage>=3),
            (select min(k.charttime) from b.mimiciv_derived.rrt k where k.stay_id=a.stay_id and k.charttime>a.t0 and k.charttime<=a.t0+interval '7d' and k.dialysis_present=1)) primary_time,
      case when a.deathtime is not null and a.deathtime<=a.outtime then a.deathtime else a.outtime end competing_time,
      p.dod
      from n09_action a join b.mimiciv_hosp.patients p using(subject_id)
      left join lateral(select creat,creat_low_past_7day low7,creat_low_past_48hr low48,
        case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3) then 3 when creat>=2*creat_low_past_7day then 2 when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat
        from b.mimiciv_derived.kdigo_creatinine k where k.stay_id=a.stay_id and k.charttime<=a.t0 order by charttime desc limit 1)r on true
      left join lateral(select arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=50983)sodium,
        arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=50882)bicarbonate,
        arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=50813)lactate,
        arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=50971)potassium,
        arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=51222)hemoglobin,
        arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=51265)platelet
        from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.valuenum is not null and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h')lab on true
      left join lateral(select avg(mbp)map_mean,min(mbp)map_min,avg(heart_rate)heart_rate from b.mimiciv_derived.vitalsign v where v.stay_id=a.stay_id and v.charttime<=a.t0 and v.charttime>a.t0-interval '6h')vit on true
      where coalesce(r.aki_stage_creat,0)<3 and not exists(select 1 from b.mimiciv_derived.rrt x where x.stay_id=a.stay_id and x.charttime<=a.t0 and x.dialysis_present=1)"""
    # Bound correlated-query expansion without changing any patient-level SQL rules.
    count = con.execute("select count(*) from n09_action").fetchone()[0]
    parts = []
    for offset in range(0, count, 24):
        con.execute(f"create or replace temp table n09_action_batch as select * from n09_action order by subject_id,stay_id limit 24 offset {offset}")
        parts.append(con.execute(query.replace("from n09_action a", "from n09_action_batch a")).df())
        print("N09 source batch", min(offset + 24, count), "/", count, flush=True)
    df = pd.concat(parts, ignore_index=True)
    df = event_columns(df)
    numeric = ["age","chloride","prior_chloride","ns_ml","lr_ml","action_amount","sofa_score","creat","low7","aki_stage_creat","sodium","bicarbonate","lactate","potassium","hemoglobin","platelet","map_mean","map_min","heart_rate","pressor","ventilation"]
    categorical = ["gender","first_careunit","anchor_year_group"]
    return df, numeric, categorical, con


def run_n09() -> None:
    df, numeric, categorical, con = n09_cohort()
    point = ato_effects(df, numeric, categorical)
    balance, diagnostic = ato_diagnostics(df, numeric, categorical)
    boot, errors = bootstrap_ato(df, numeric, categorical)
    rows = []
    for outcome, values in point.items():
        row = {"outcome": outcome, **values}
        for key in ["control_risk","intervention_risk","rd","rr"]:
            row[f"{key}_ci_low"], row[f"{key}_ci_high"] = np.quantile(boot[f"{outcome}_{key}"], [.025,.975])
        rows.append(row)
    output = outdir("N09")
    protocol = {"database":"MIMIC-IV","population":"Adult first ICU sepsis; first known chloride >=110 mmol/L or delta chloride >=5 in 48h within 72h ICU; >=500 mL completed NaCl 0.9% in prior 24h; next >=250 mL crystalloid bolus within 4h; no baseline AKI3/KRT.","time_zero":"Start of first qualifying post-trigger crystalloid bolus.","intervention":"Lactated Ringer's bolus","control":"NaCl 0.9% bolus","estimand":"ATO among patients with clinical equipoise for the recorded next bolus; subsequent care natural.","imputation":"Median plus missing indicator/most-frequent category, refit in each patient bootstrap.","bootstrap":BOOTSTRAP,"limits":["This is next-bolus fluid choice, not a sustained chloride-restriction policy.","KDIGO outcome timing uses derived charttime; baseline treatment covariates use data available by t0 where possible."]}
    flow = pd.DataFrame([
        {"step":"final eligible next-bolus decisions","n":len(df)},
        {"step":"NaCl 0.9%","n":int((df.exposure==0).sum())},
        {"step":"Lactated Ringer's","n":int((df.exposure==1).sum())},
        {"step":"7d AKI3/KRT events","n":int((df.aki3_krt_7d_event==1).sum())},
    ])
    outcomes = pd.DataFrame(rows)
    diagnostics = {**diagnostic,"bootstrap_requested":BOOTSTRAP,"bootstrap_success":len(boot),"bootstrap_errors":errors,"era_by_arm":df.groupby(["anchor_year_group","exposure"]).size().rename("n").reset_index().to_dict("records"),"missingness":df[numeric+categorical].isna().mean().to_dict(),"status":"FULL_MIMIC_ATO_COMPLETE"}
    write_json(output/"protocol.json",protocol);flow.to_csv(output/"flow.csv",index=False);balance.to_csv(output/"balance.csv",index=False);outcomes.to_csv(output/"outcomes.csv",index=False);write_json(output/"diagnostics.json",diagnostics)
    report = f"""# N09 MIMIC-IV canonical rerun\n\nN={len(df)}（LR {(df.exposure==1).sum()} / NaCl {(df.exposure==0).sum()}），采用重叠加权估计下一次晶体液推注选择。500次患者级bootstrap成功{len(boot)}次。\n\n最大加权后绝对SMD={diagnostic['max_abs_smd']:.3f}；ESS LR={diagnostic['ess_intervention']:.1f}、NaCl={diagnostic['ess_control']:.1f}。该结果只解释为“高氯触发后的下一次LR vs NaCl推注”，不能扩大为持续限制氯策略。主要和次要结果见outcomes.csv。\n"""
    (output/"report.md").write_text(report,encoding="utf-8")


def load_ccw_module():
    path = ROOT / "02代码/01分析脚本/48_full_ccw_initiation.py"
    spec = importlib.util.spec_from_file_location("ccw", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


N02_CTE = r"""
with base as (
 select p.*,p.intime+interval '24h' t0,z.sofa_score,greatest(z.antibiotic_time,z.culture_time,z.sofa_time)sepsis_known
 from s.main.people p join b.mimiciv_derived.sepsis3 z using(stay_id,subject_id)
 where p.age>=18 and z.sepsis3 and p.outtime>p.intime+interval '24h'
   and greatest(z.antibiotic_time,z.culture_time,z.sofa_time)<=p.intime+interval '24h'
), history as (
 select x.*,f.input_ml,o.output_ml,coalesce(f.input_ml,0)-coalesce(o.output_ml,0)balance_ml,
   v.map_n,v.map_min,v.map_mean,v.heart_rate,r.creat,r.low7,r.low48,r.aki_stage_creat,
   lab.sodium,lab.potassium,lab.bicarbonate,lab.chloride,lab.lactate,lab.hemoglobin,lab.platelet,p.dod
 from base x join b.mimiciv_hosp.patients p using(subject_id)
 left join lateral(select sum(amount)input_ml from b.mimiciv_icu.inputevents i where i.stay_id=x.stay_id and i.endtime<=x.t0 and i.endtime>x.t0-interval '24h' and i.amount>0 and i.amountuom='mL' and i.statusdescription='FinishedRunning')f on true
 left join lateral(select sum(value)output_ml from b.mimiciv_icu.outputevents o where o.stay_id=x.stay_id and o.charttime<=x.t0 and o.charttime>x.t0-interval '24h' and o.value>0 and lower(o.valueuom)='ml')o on true
 left join lateral(select count(mbp)map_n,min(mbp)map_min,avg(mbp)map_mean,avg(heart_rate)heart_rate from b.mimiciv_derived.vitalsign v where v.stay_id=x.stay_id and v.charttime<=x.t0 and v.charttime>x.t0-interval '6h' and mbp between 20 and 200)v on true
 left join lateral(select creat,creat_low_past_7day low7,creat_low_past_48hr low48,case when creat>=3*creat_low_past_7day or (creat>=4 and creat>=creat_low_past_48hr+.3)then 3 when creat>=2*creat_low_past_7day then 2 when creat>=1.5*creat_low_past_7day or creat>=creat_low_past_48hr+.3 then 1 else 0 end aki_stage_creat from b.mimiciv_derived.kdigo_creatinine k where k.stay_id=x.stay_id and k.charttime<=x.t0 order by charttime desc limit 1)r on true
 left join lateral(select arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=50983)sodium,arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=50971)potassium,arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=50882)bicarbonate,arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=50902)chloride,arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=50813)lactate,arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=51222)hemoglobin,arg_max(valuenum,(charttime,storetime,labevent_id))filter(where itemid=51265)platelet from b.mimiciv_hosp.labevents l where l.hadm_id=x.hadm_id and l.valuenum is not null and greatest(l.charttime,l.storetime)<=x.t0 and l.charttime>x.t0-interval '24h')lab on true
), eligible as (
 select * from history x where aki_stage_creat in(1,2) and balance_ml>500 and map_n>=2 and map_min>=65
 and not exists(select 1 from b.mimiciv_derived.norepinephrine n where n.stay_id=x.stay_id and n.starttime<=x.t0 and n.endtime>x.t0)
 and not exists(select 1 from b.mimiciv_derived.vasopressin n where n.stay_id=x.stay_id and n.starttime<=x.t0 and n.endtime>x.t0)
 and not exists(select 1 from b.mimiciv_derived.rrt n where n.stay_id=x.stay_id and n.charttime<=x.t0 and n.dialysis_active=1)
 and not exists(select 1 from b.mimiciv_icu.inputevents i where i.stay_id=x.stay_id and i.itemid in(221794,228340,229639) and i.starttime<x.t0 and i.endtime>x.t0-interval '6h' and i.amount>0 and i.statusdescription is distinct from 'Rewritten')
), action as (
 select e.*,a.starttime action_time,a.itemid action_itemid
 from eligible e left join lateral(select * from b.mimiciv_icu.inputevents i where i.stay_id=e.stay_id and i.itemid in(221794,228340,229639) and i.starttime>=e.t0 and i.starttime<e.t0+interval '4h' and i.amount>0 and i.statusdescription is distinct from 'Rewritten' order by i.starttime limit 1)a on true
)
"""


def n02_cohort_and_panel(con):
    con.execute("create temp table n02_action as "+N02_CTE+" select * from action")
    df = con.execute("""select a.*,
      least((select min(k.charttime)from b.mimiciv_derived.kdigo_stages k where k.stay_id=a.stay_id and k.charttime>a.t0 and k.charttime<=a.t0+interval '7d' and k.aki_stage>=3),
            (select min(k.charttime)from b.mimiciv_derived.rrt k where k.stay_id=a.stay_id and k.charttime>a.t0 and k.charttime<=a.t0+interval '7d' and k.dialysis_present=1))primary_time,
      case when deathtime is not null and deathtime<=outtime then deathtime else outtime end competing_time
      from n02_action a""").df()
    df = event_columns(df)
    df["treatment_hour"]=(pd.to_datetime(df.action_time)-pd.to_datetime(df.t0)).dt.total_seconds()/3600
    df["followup"]=df.aki3_krt_7d_time
    df["event"]=df.aki3_krt_7d_event
    df=df.reset_index(drop=True);df["pid"]=np.arange(len(df))
    basecols=["pid","gender","first_careunit","anchor_year_group","age","sofa_score","input_ml","output_ml","balance_ml","map_min","map_mean","heart_rate","creat","low7","low48","aki_stage_creat","sodium","potassium","bicarbonate","chloride","lactate","hemoglobin","platelet"]
    panel=df[basecols].loc[df.index.repeat(4)].reset_index(drop=True)
    panel["k"]=np.tile(np.arange(4),len(df))
    return df,panel


def ccw_outcome_frame(df, outcome):
    co=df[["pid","treatment_hour"]].copy()
    co["followup"]=df[f"{outcome}_time"]
    co["event"]=df[f"{outcome}_event"]
    return co


def _run_n02_full_not_current() -> None:
    con=connection();df,panel=n02_cohort_and_panel(con);ccw=load_ccw_module()
    outcomes=["aki3_krt_7d","death_7d","death_28d","death_90d"]
    point={};main_diag=None;main_balance=None
    for outcome in outcomes:
        result=ccw.run_effect(ccw_outcome_frame(df,outcome),panel,np.ones(len(df),int),outcome=="aki3_krt_7d")
        if outcome=="aki3_krt_7d":
            values,main_diag,main_balance,_,_=result
        else: values=result
        point[outcome]=values
    rng=np.random.default_rng(SEED);records=[];errors=[]
    for replicate in range(BOOTSTRAP):
        mult=np.bincount(rng.integers(0,len(df),len(df)),minlength=len(df))
        try:
            row={"replicate":replicate}
            for outcome in outcomes:
                values=ccw.run_effect(ccw_outcome_frame(df,outcome),panel,mult)
                row.update({f"{outcome}_{k}":v for k,v in values.items()})
            records.append(row)
        except Exception as exc: errors.append({"replicate":replicate,"error":repr(exc)})
    boot=pd.DataFrame(records);rows=[]
    for outcome,values in point.items():
        row={"outcome":outcome,"control_risk":values["defer_risk"],"intervention_risk":values["initiate_risk"],"rd":values["rd"],"rr":values["initiate_risk"]/values["defer_risk"] if values["defer_risk"]>0 else np.nan}
        for key,source in [("control_risk","defer_risk"),("intervention_risk","initiate_risk"),("rd","rd")]:
            row[f"{key}_ci_low"],row[f"{key}_ci_high"]=np.quantile(boot[f"{outcome}_{source}"],[.025,.975])
        ratio=boot[f"{outcome}_initiate_risk"]/boot[f"{outcome}_defer_risk"]
        row["rr_ci_low"],row["rr_ci_high"]=np.nanquantile(ratio,[.025,.975]);rows.append(row)
    output=outdir("N02")
    protocol={"database":"MIMIC-IV","population":"Adult first ICU sepsis with creatinine AKI1/2 at 24h landmark; completed-input minus recorded-output >500mL in preceding 24h; >=2 MAP values and minimum MAP>=65 in prior 6h; no active NE/VP, prior KRT, or IV loop in prior 6h.","time_zero":"ICU admission +24h common landmark.","intervention":"Initiate actual IV furosemide/bumetanide within 4h","control":"Defer IV loop for 4h; care thereafter natural","primary":"7d recorded AKI3/new KRT cumulative incidence","secondary":["7d","28d","90d recorded all-cause mortality"],"estimator":"Continuous initiation CCW with hourly treatment-intensity model and patient bootstrap.","imputation":"Median plus missing indicators inside treatment model and each bootstrap.","bootstrap":BOOTSTRAP,"limits":["Positive balance uses completed mL input and recorded mL output; unrecorded output is not zero.","Only 24h landmark estimates this precise policy; it does not represent all later de-resuscitation decisions."]}
    flow=pd.DataFrame([{"step":"final eligibility at 24h landmark","n":len(df)},{"step":"IV loop initiated within 4h","n":int(df.treatment_hour.notna().sum())},{"step":"deferred >=4h","n":int(df.treatment_hour.isna().sum())},{"step":"7d AKI3/KRT events","n":int((df.aki3_krt_7d_event==1).sum())}])
    write_json(output/"protocol.json",protocol);flow.to_csv(output/"flow.csv",index=False);main_balance.to_csv(output/"balance.csv",index=False);pd.DataFrame(rows).to_csv(output/"outcomes.csv",index=False)
    diagnostic={**main_diag,"bootstrap_requested":BOOTSTRAP,"bootstrap_success":len(boot),"bootstrap_errors":errors,"missingness":panel.drop(columns=["pid","k"]).isna().mean().to_dict(),"status":"FULL_MIMIC_CCW_COMPLETE"};write_json(output/"diagnostics.json",diagnostic)
    (output/"report.md").write_text(f"# N02 MIMIC-IV canonical rerun\n\n24h共同零点最终纳入{len(df)}人；4h内实际启动IV袢利尿{df.treatment_hour.notna().sum()}人，暂缓{df.treatment_hour.isna().sum()}人。500次患者级bootstrap成功{len(boot)}次。\n\n宽限期末ESS为{main_diag['arms'][0]['ess']:.1f}/{main_diag['arms'][1]['ess']:.1f}，最大绝对SMD={main_diag['max_abs_smd']:.3f}。主要与7/28/90天死亡结果见outcomes.csv；该估计只适用于24h正平衡且稳定的清洁启动决策。\n",encoding="utf-8")


def run_n02() -> None:
    """Write the frozen lightweight gate; do not repeat the aborted large join."""
    output=outdir("N02")
    protocol={
        "database":"MIMIC-IV",
        "population":"Adult first ICU sepsis with creatinine AKI1/2 at ICU+24h, preceding 24h completed-input minus recorded-output >500mL, >=2 MAP values and minimum MAP>=65 in prior 6h, no active NE/VP/KRT and no IV loop in prior 6h.",
        "time_zero":"ICU admission +24h common landmark.",
        "candidate_intervention":"Actual IV furosemide/bumetanide initiation during [0,4h).",
        "candidate_control":"No IV loop initiation during [0,4h); later care natural.",
        "decision":"GATE STOP before effect estimation.",
        "reason":"The clean initiation action occurred in only 84/1237 eligible patients and the attempted time-varying source join caused excessive DuckDB temporary-disk expansion. No outcome, weight, ESS, SMD or bootstrap result is claimed.",
        "future_requirement":"Rebuild a bounded patient-keyed panel from cached source tables and verify that a first IV dose represents active de-resuscitation rather than diagnostic/PRN administration before any causal fit."
    }
    flow=pd.DataFrame([
        {"step":"eligible at ICU+24h after all stated clinical restrictions and no prior-6h IV loop","n":1237},
        {"step":"actual IV loop initiation within 4h","n":84},
        {"step":"no initiation within 4h","n":1153},
    ])
    balance=pd.DataFrame([{"status":"not_fitted","ess_initiate":np.nan,"ess_defer":np.nan,"max_abs_smd":np.nan,"reason":"time-varying source panel not completed; treatment support 84/1237"}])
    outcomes=pd.DataFrame([{"outcome":"not_estimated","effect":np.nan,"ci_low":np.nan,"ci_high":np.nan,"reason":"gate stopped before outcomes and weighting"}])
    diagnostics={"status":"GATE_STOP_NO_EFFECT_ESTIMATE","eligible_n":1237,"initiate_n":84,"defer_n":1153,"initiate_fraction":84/1237,"ess":None,"max_smd":None,"bootstrap_requested":0,"temporary_database_created":False,"aborted_compute":"An unbounded correlated time-varying join expanded DuckDB temporary storage; the process was terminated and is not part of the result.","reusable_evidence":"The 1237/84 counts came from the bounded aggregate landmark query completed before the failed panel construction."}
    write_json(output/"protocol.json",protocol);flow.to_csv(output/"flow.csv",index=False);balance.to_csv(output/"balance.csv",index=False);outcomes.to_csv(output/"outcomes.csv",index=False);write_json(output/"diagnostics.json",diagnostics)
    (output/"report.md").write_text("# N02 MIMIC-IV canonical gate\n\n24小时共同复评点共有1,237名患者满足脓毒症、SCr-AKI1/2、前24小时正平衡>500mL、血流动力学稳定、无活动NE/VP/KRT且前6小时无静脉袢利尿；其中84人在随后4小时内实际新启动静脉袢利尿，1,153人暂缓。\n\n当前仅确认动作源计数。因启动比例低，且未完成有界的逐时混杂面板，不能报告ESS、SMD、主要/次要结局或治疗效应。一次无界相关连接造成DuckDB临时磁盘膨胀，已终止且未作为研究结果。后续若重启，应先按这1,237个patient key构建有界缓存，再核首剂是否真正代表主动去复苏。\n",encoding="utf-8")


def run_n07() -> None:
    con=connection(include_abx=True)
    sql=r"""with aki as(select * from a.main.expanded_baseline qualify row_number()over(partition by stay_id,drug order by available_time)=1),recovery as(select k.*,r.available_time recovery_time,r.scr recovery_scr from aki k join lateral(select * from a.main.aki_results r where r.hadm_id=k.hadm_id and r.available_time>k.available_time and r.available_time<least(k.outtime,k.available_time+interval '7d') and r.scr<=k.scr-.3 and r.scr<=1.25*k.low7 order by r.available_time limit 1)r on true where not exists(select 1 from b.mimiciv_derived.rrt x where x.stay_id=k.stay_id and x.charttime<=r.available_time and x.dialysis_active=1)),prev as(select r.*,p.charttime base_time,p.dose_g base_dose,p.doses_day base_day,p.frequency base_freq,p.pharmacy_id base_pharm,pp.charttime previous_time,pp.dose_g previous_dose,pp.doses_day previous_day from recovery r join lateral(select * from a.main.expanded_actual x where x.hadm_id=r.hadm_id and x.drug=r.drug and x.charttime<r.recovery_time and x.charttime>=r.recovery_time-interval '24h' and greatest(x.charttime,x.storetime)<=r.recovery_time and x.dose_g>0 and x.doses_day>0 and upper(coalesce(x.route,''))like'IV%' order by x.charttime desc limit 1)p on true join lateral(select * from a.main.expanded_actual x where x.hadm_id=r.hadm_id and x.drug=r.drug and x.charttime<p.charttime and x.charttime>=p.charttime-interval '24h' and greatest(x.charttime,x.storetime)<=r.recovery_time and x.dose_g>0 and x.doses_day>0 and upper(coalesce(x.route,''))like'IV%' order by x.charttime desc limit 1)pp on true where p.dose_g=pp.dose_g and p.doses_day=pp.doses_day and p.frequency<>'ONCE'),nxt as(select p.*,n.charttime next_time,n.dose_g next_dose,n.doses_day next_day,n.frequency next_frequency,n.pharmacy_id next_pharm from prev p left join lateral(select * from a.main.expanded_actual x where x.hadm_id=p.hadm_id and x.drug=p.drug and x.charttime>=p.recovery_time and x.charttime<p.recovery_time+interval '24h' and x.dose_g>0 and x.doses_day>0 and upper(coalesce(x.route,''))like'IV%' order by x.charttime limit 1)n on true)select *,case when next_dose*next_day>base_dose*base_day and next_pharm<>base_pharm then 1 when next_dose=base_dose and next_day=base_day then 0 else null end exposure from nxt"""
    df=con.execute(sql).df();output=outdir("N07")
    flow=pd.DataFrame([{"step":"AKI episodes with qualifying creatinine recovery","n":int(df.stay_id.nunique())},{"step":"stable two-dose IV beta-lactam regimen at recovery","n":len(df)},{"step":"executed up-titration/shorter interval with changed order within 24h","n":int((df.exposure==1).sum())},{"step":"same regimen execution","n":int((df.exposure==0).sum())},{"step":"unclassified/no next action","n":int(df.exposure.isna().sum())}])
    bydrug=df.groupby("drug").agg(stable_recovery=("stay_id","size"),up=("exposure",lambda x:int((x==1).sum())),maintain=("exposure",lambda x:int((x==0).sum())),unclassified=("exposure",lambda x:int(x.isna().sum()))).reset_index()
    protocol={"database":"MIMIC-IV","population":"Sepsis-AKI patients with subsequent known creatinine decline >=0.3 mg/dL and <=1.25 times pre-AKI 7d reference, still receiving a stable two-dose IV cefepime/piperacillin-tazobactam/meropenem regimen, no active prior KRT.","time_zero":"First qualifying recovery result availability.","intervention":"Executed higher daily dose intensity or shorter interval supported by a changed pharmacy order within 24h.","control":"Executed unchanged dose and frequency.","decision":"STOP before outcome estimation because only four supported up-adjustments were observed; pooling drug-specific regimens would violate consistency and cannot repair positivity."}
    write_json(output/"protocol.json",protocol);flow.to_csv(output/"flow.csv",index=False);bydrug.to_csv(output/"balance.csv",index=False);pd.DataFrame([{"outcome":"not_estimated","reason":"treatment positivity failure: four supported up-adjustments across all three drugs"}]).to_csv(output/"outcomes.csv",index=False)
    diagnostics={"status":"ACTION_IDENTIFIABLE_BUT_POSITIVITY_FAIL_NO_EFFECT_ESTIMATE","stable_recovery_n":len(df),"up_n":int((df.exposure==1).sum()),"maintain_n":int((df.exposure==0).sum()),"unclassified_n":int(df.exposure.isna().sum()),"by_drug":bydrug.to_dict("records"),"ess":None,"max_smd":None,"bootstrap_requested":0,"reason":"No credible overlap model or bootstrap can be fitted with four supported intervention actions."};write_json(output/"diagnostics.json",diagnostics)
    (output/"report.md").write_text(f"# N07 MIMIC-IV canonical rerun\n\n从真实eMAR剂量、频率和药房医嘱重建后，{len(df)}个稳定维持方案的肾恢复节点中，仅{int((df.exposure==1).sum())}例存在有新医嘱支持的实际上调，{int((df.exposure==0).sum())}例维持原方案，{int(df.exposure.isna().sum())}例不可分类。分药后每个药物上调仅0–2例。\n\n因此动作可以定义，但MIMIC-IV不存在足够的治疗支持；本轮不读取结局、不拟合权重、不bootstrap。任何合并药物或放宽为“下一剂间隔变短”都会混入自然给药节奏，不能形成可信TTE。\n",encoding="utf-8")


def main():
    parser=argparse.ArgumentParser();parser.add_argument("candidate",choices=["N02","N07","N09","all"]);args=parser.parse_args()
    selected=["N02","N07","N09"] if args.candidate=="all" else [args.candidate]
    for candidate in selected:
        globals()[f"run_{candidate.lower()}"]()


if __name__ == "__main__":
    main()
