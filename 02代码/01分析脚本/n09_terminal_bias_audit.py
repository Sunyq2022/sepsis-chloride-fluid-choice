"""Terminal bias audit for the frozen 706-patient N09 cohort.

No new subgroup is searched. The script tests whether future LR/NaCl choice is
associated with three pre-time-zero placebo outcomes and reruns the frozen
primary outcome in two prespecified support/clinical-context restrictions.
Only aggregate files are written.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT=Path(__file__).resolve().parents[2]
N09=ROOT/"03结果"/"21_FULL_ICU_CAUSAL_PROGRAM_20260907"/"N09"
COHORT=N09/"publication_figures"/"cohort_private.parquet"
OUT=N09/"final_yellow_audit_20260909"
CANON=ROOT/"02代码"/"01分析脚本"/"canonical_rerun_n02_n07_n09.py"
CONF=ROOT/"02代码"/"01分析脚本"/"n09_pre_t0_confounder_audit.py"
BALANCE=ROOT/"02代码"/"01分析脚本"/"n09_balance_optimized_sensitivity.py"
CSV_OUT=OUT/"N09_terminal_bias_audit.csv"
JSON_OUT=OUT/"N09_terminal_bias_audit.json"
MD_OUT=OUT/"N09_terminal_bias_audit.md"
SEED=20260908
BOOT=500

BASE_NUMERIC=["age","chloride","prior_chloride","ns_ml","lr_ml","sofa_score","creat","low7",
 "aki_stage_creat","sodium","bicarbonate","lactate","potassium","hemoglobin","platelet",
 "map_mean","map_min","heart_rate","pressor","ventilation"]
PARSIMONIOUS=["urine6_ml","urine24_ml","urine24_hours","net6_ml","net24_ml","known_ns24_ml",
 "known_lr24_ml","hours_since_prior_ns","prior_urine_response_delta","creat_delta24",
 "chloride_delta48","lactate_delta6"]
CATEGORICAL=["gender","first_careunit","anchor_year_group"]
PSEUDO=[
 ("pre24_aki_worsening","t0前24小时AKI分期恶化"),
 ("pre6_pressor_escalation","t0前6小时升压药升级代理"),
 ("pre6_lactate_worsening","t0前6小时乳酸升高至少0.5 mmol/L"),
]


def load(path:Path,name:str):
    spec=importlib.util.spec_from_file_location(name,path); mod=importlib.util.module_from_spec(spec)
    assert spec.loader; spec.loader.exec_module(mod); return mod


def add_parsimonious_features(con,cohort:pd.DataFrame,conf) -> pd.DataFrame:
    extras=conf.extra_features(con,cohort)
    frame=cohort.merge(extras,on="stay_id",how="left",validate="one_to_one")
    frame[["known_ns24_ml","known_lr24_ml"]]=frame[["known_ns24_ml","known_lr24_ml"]].fillna(0)
    frame["net6_ml"]=frame.input6_ml-frame.output6_ml
    frame["net24_ml"]=frame.input24_ml-frame.output24_ml
    frame["hours_since_prior_ns"]=(pd.to_datetime(frame.t0)-pd.to_datetime(frame.prior_ns_end)).dt.total_seconds()/3600
    frame["pre_urine_rate"]=frame.pre_urine_ml/frame.prior_response_pre_hours.where(frame.prior_response_pre_hours>=1)
    frame["post_urine_rate"]=frame.post_urine_ml/frame.prior_response_post_hours.where(frame.prior_response_post_hours>=1)
    frame["prior_urine_response_delta"]=frame.post_urine_rate-frame.pre_urine_rate
    frame["creat_delta24"]=frame.creat_last24-frame.creat_first24
    frame["chloride_delta48"]=frame.chloride_last48-frame.chloride_first48
    frame["lactate_delta6"]=frame.lactate_last6-frame.lactate_first6
    return frame


def placebo_features(con,frame:pd.DataFrame) -> pd.DataFrame:
    con.register("n09_fixed_keys",frame[["stay_id","hadm_id","t0"]])
    result=con.execute("""
      WITH pressor AS (
        SELECT stay_id,starttime,endtime,vaso_rate,'NE' drug FROM b.mimiciv_derived.norepinephrine
        UNION ALL SELECT stay_id,starttime,endtime,vaso_rate,'VP' FROM b.mimiciv_derived.vasopressin
        UNION ALL SELECT stay_id,starttime,endtime,vaso_rate,'EPI' FROM b.mimiciv_derived.epinephrine
        UNION ALL SELECT stay_id,starttime,endtime,vaso_rate,'DOPA' FROM b.mimiciv_derived.dopamine
      )
      SELECT c.stay_id,
        k.first_stage,k.max_stage,k.stage_n,
        l.first_lactate,l.last_lactate,l.lactate_n,
        p.agents_prior6,p.agents_last6,p.ne_prior6,p.ne_last6
      FROM n09_fixed_keys c
      LEFT JOIN LATERAL (
        SELECT arg_min(aki_stage,charttime) first_stage,max(aki_stage) max_stage,count(*) stage_n
        FROM b.mimiciv_derived.kdigo_stages x
        WHERE x.stay_id=c.stay_id AND x.charttime>=c.t0-INTERVAL '24 hours' AND x.charttime<c.t0
      ) k ON true
      LEFT JOIN LATERAL (
        SELECT arg_min(valuenum,greatest(charttime,storetime)) first_lactate,
               arg_max(valuenum,greatest(charttime,storetime)) last_lactate,count(*) lactate_n
        FROM b.mimiciv_hosp.labevents x
        WHERE x.hadm_id=c.hadm_id AND x.itemid=50813 AND x.valuenum BETWEEN 0.1 AND 40
          AND greatest(x.charttime,x.storetime)>=c.t0-INTERVAL '6 hours'
          AND greatest(x.charttime,x.storetime)<c.t0
      ) l ON true
      LEFT JOIN LATERAL (
        SELECT count(DISTINCT drug) FILTER(WHERE starttime<c.t0-INTERVAL '6 hours' AND endtime>c.t0-INTERVAL '12 hours') agents_prior6,
               count(DISTINCT drug) FILTER(WHERE starttime<c.t0 AND endtime>c.t0-INTERVAL '6 hours') agents_last6,
               max(vaso_rate) FILTER(WHERE drug='NE' AND starttime<c.t0-INTERVAL '6 hours' AND endtime>c.t0-INTERVAL '12 hours') ne_prior6,
               max(vaso_rate) FILTER(WHERE drug='NE' AND starttime<c.t0 AND endtime>c.t0-INTERVAL '6 hours') ne_last6
        FROM pressor x WHERE x.stay_id=c.stay_id AND x.starttime<c.t0 AND x.endtime>c.t0-INTERVAL '12 hours'
      ) p ON true
    """).df()
    con.unregister("n09_fixed_keys")
    result["pre24_aki_worsening"]=np.where(result.stage_n>=2,(result.max_stage>result.first_stage).astype(float),np.nan)
    prior_agents=result.agents_prior6.fillna(0); last_agents=result.agents_last6.fillna(0)
    prior_ne=result.ne_prior6.fillna(0); last_ne=result.ne_last6.fillna(0)
    result["pre6_pressor_escalation"]=((last_agents>prior_agents)|((last_ne>=prior_ne*1.25+0.01)&(last_ne>0))).astype(float)
    result["pre6_lactate_worsening"]=np.where(result.lactate_n>=2,((result.last_lactate-result.first_lactate)>=0.5).astype(float),np.nan)
    return frame.merge(result[["stay_id"]+[x[0] for x in PSEUDO]],on="stay_id",how="left",validate="one_to_one")


def effects(frame:pd.DataFrame,canonical,balance_mod,numeric:list[str],mult=None) -> dict[str,dict[str,float]]:
    expanded,_,_,_,w,converged,_,warnings=balance_mod.unpenalized_fit(frame,numeric,CATEGORICAL,mult)
    if not converged: raise RuntimeError("unpenalized PS did not converge: "+" | ".join(warnings))
    a=expanded.exposure.to_numpy(int); out={}
    time=expanded.aki3_krt_7d_time.to_numpy(float); event=expanded.aki3_krt_7d_event.to_numpy(int)
    risk=[canonical.weighted_aj(time[a==g],event[a==g],w[a==g],168.0) for g in (0,1)]
    out["aki3_krt_7d"]={"risk_ns":risk[0],"risk_lr":risk[1],"rd":risk[1]-risk[0],"rr":risk[1]/risk[0] if risk[0]>0 else np.nan}
    for key,_ in PSEUDO:
        y=pd.to_numeric(expanded[key],errors="coerce").to_numpy(float); valid=np.isfinite(y); risk=[]
        for g in (0,1):
            use=valid&(a==g); risk.append(float(np.average(y[use],weights=w[use])))
        out[key]={"risk_ns":risk[0],"risk_lr":risk[1],"rd":risk[1]-risk[0],"rr":risk[1]/risk[0] if risk[0]>0 else np.nan,
                  "evaluable_n":int(valid.sum()),"events":int(np.nansum(y))}
    return out


def analyse(name:str,frame:pd.DataFrame,canonical,balance_mod,numeric:list[str]) -> tuple[list[dict],dict]:
    point=effects(frame,canonical,balance_mod,numeric)
    diag,balance=balance_mod.diagnose(frame,numeric,CATEGORICAL,name)
    rng=np.random.default_rng(SEED); collected={k:{"rd":[],"rr":[]} for k in point}; errors=[]
    n=len(frame)
    for b in range(BOOT):
        mult=np.bincount(rng.integers(0,n,n),minlength=n)
        try:
            one=effects(frame,canonical,balance_mod,numeric,mult)
            for key in collected:
                collected[key]["rd"].append(one[key]["rd"]); collected[key]["rr"].append(one[key]["rr"])
        except Exception as exc: errors.append(type(exc).__name__)
    labels={"aki3_krt_7d":"7天AKI3或新KRT",**dict(PSEUDO)}; rows=[]
    for key,value in point.items():
        rd=np.asarray(collected[key]["rd"]); rr=np.asarray(collected[key]["rr"])
        rows.append({"analysis":name,"endpoint":key,"endpoint_label":labels[key],"endpoint_type":"主要结局" if key=="aki3_krt_7d" else "治疗前伪结局",
          "n":len(frame),"lr_n":int(frame.exposure.sum()),"ns_n":int((1-frame.exposure).sum()),
          "evaluable_n":int(value.get("evaluable_n",len(frame))),"events":int(value.get("events",(frame.aki3_krt_7d_event==1).sum())),
          "risk_ns_pct":100*value["risk_ns"],"risk_lr_pct":100*value["risk_lr"],"rd_pp":100*value["rd"],
          "rd_ci_low_pp":100*np.quantile(rd,.025),"rd_ci_high_pp":100*np.quantile(rd,.975),"rr":value["rr"],
          "rr_ci_low":np.quantile(rr,.025),"rr_ci_high":np.quantile(rr,.975),"bootstrap_success":len(rd)})
    meta={"analysis":name,"n":len(frame),"lr_n":int(frame.exposure.sum()),"ns_n":int((1-frame.exposure).sum()),
          "ess_lr":diag["ess_lr"],"ess_ns":diag["ess_ns"],"max_smd":diag["max_abs_smd"],
          "bootstrap_success":BOOT-len(errors),"bootstrap_errors":sorted(set(errors)),"worst_balance":balance.head(5).to_dict("records")}
    return rows,meta


def main() -> None:
    OUT.mkdir(parents=True,exist_ok=True)
    canonical=load(CANON,"n09_terminal_canonical"); conf=load(CONF,"n09_terminal_conf")
    balance_mod=load(BALANCE,"n09_terminal_balance")
    frame=pd.read_parquet(COHORT)
    if len(frame)!=706 or int(frame.exposure.sum())!=259: raise RuntimeError("frozen key-set mismatch")
    con=canonical.connection(); con.execute("set threads=1")
    frame=add_parsimonious_features(con,frame,conf); frame=placebo_features(con,frame); con.close()
    numeric=BASE_NUMERIC+PARSIMONIOUS
    cell=frame.groupby(["first_careunit","anchor_year_group","exposure"],dropna=False).size().unstack(fill_value=0)
    valid_cells=cell[(cell.get(0,0)>=5)&(cell.get(1,0)>=5)].index
    cell_keys=pd.MultiIndex.from_frame(frame[["first_careunit","anchor_year_group"]])
    same_cell=frame[cell_keys.isin(valid_cells)].copy().reset_index(drop=True)
    surgical=frame.first_careunit.astype(str).str.contains("Surgical|Trauma|Cardiac Vascular|CVICU|PACU|Surgery",case=False,regex=True)
    surgical |= frame.admission_type.astype(str).str.contains("SURG",case=False,regex=True)
    nonsurgical=frame[~surgical].copy().reset_index(drop=True)
    scenarios=[("fixed_706",frame),("icu_era_cells_min5_each_arm",same_cell),("exclude_surgical_trauma_cvicu",nonsurgical)]
    all_rows=[]; metas=[]
    for name,data in scenarios:
        rows,meta=analyse(name,data,canonical,balance_mod,numeric); all_rows.extend(rows); metas.append(meta)
    table=pd.DataFrame(all_rows); table.to_csv(CSV_OUT,index=False,encoding="utf-8-sig")
    pseudo_fixed=table[(table.analysis=="fixed_706")&(table.endpoint_type=="治疗前伪结局")]
    placebo_flag=bool(((pseudo_fixed.rd_ci_low_pp>0)|(pseudo_fixed.rd_ci_high_pp<0)).any())
    primary=table[table.endpoint.eq("aki3_krt_7d")]
    sign_stable=bool((np.sign(primary.rd_pp)==np.sign(primary.iloc[0].rd_pp)).all())
    decision="RED_RESIDUAL_CONFOUNDING_SIGNAL" if placebo_flag else ("YELLOW_DIRECTION_STABLE" if sign_stable else "YELLOW_DIRECTION_UNSTABLE")
    payload={"candidate":"N09","fixed_keyset":"N=706, LR=259, NaCl=447","new_subgroups_searched":False,
      "placebo_definitions":{"pre24_aki_worsening":"at least two KDIGO-stage records in [t0-24h,t0), with maximum stage greater than the earliest stage",
       "pre6_pressor_escalation":"last-6h distinct pressor count exceeds prior-6h count, or maximum norepinephrine rises by >=25% plus 0.01; proxy only",
       "pre6_lactate_worsening":"at least two result-available lactates in [t0-6h,t0), last minus first >=0.5 mmol/L"},
      "support_restriction":"ICU-type x anchor-year cells require >=5 patients in each arm","clinical_restriction":"exclude surgical/trauma/CVICU/PACU care-unit or surgical admission proxy",
      "models":"selected unpenalized pre-t0 parsimonious ATO specification; median/missing-indicator imputation and PS refit in every bootstrap",
      "scenarios":metas,"placebo_signal":placebo_flag,"primary_direction_stable":sign_stable,"decision":decision,
      "interpretation":"A pre-treatment placebo association is evidence of residual treatment-selection bias, not an effect of fluid. Absence of association cannot prove exchangeability.",
      "patient_identifiers_written":False,"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    JSON_OUT.write_text(json.dumps(payload,ensure_ascii=False,indent=2,default=str),encoding="utf-8")
    lines=["# N09黄色终审：治疗前伪结局与限定敏感性","",f"固定队列706人（LR 259、NaCl 447）；未搜索新亚组。当前判定：**{decision}**。","",
      "|分析|类型|结局|N（LR/NS）|可评价N|LR/NS风险|RD百分点（500次95%区间）|RR（500次95%区间）|ESS LR/NS|最大SMD|",
      "|---|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    meta_lookup={x["analysis"]:x for x in metas}
    for r in all_rows:
        m=meta_lookup[r["analysis"]]
        lines.append(f"|{r['analysis']}|{r['endpoint_type']}|{r['endpoint_label']}|{r['n']}（{r['lr_n']}/{r['ns_n']}）|{r['evaluable_n']}|{r['risk_lr_pct']:.2f}%/{r['risk_ns_pct']:.2f}%|{r['rd_pp']:+.2f}（{r['rd_ci_low_pp']:+.2f}至{r['rd_ci_high_pp']:+.2f}）|{r['rr']:.3f}（{r['rr_ci_low']:.3f}至{r['rr_ci_high']:.3f}）|{m['ess_lr']:.1f}/{m['ess_ns']:.1f}|{m['max_smd']:.3f}|")
    lines += ["","## 解释","","治疗前伪结局发生在液体选择之前，因此任何组间差异都不能由后续LR或NaCl引起；它提示治疗选择仍携带病情轨迹信息。伪结局无关联也不能证明不存在未测混杂。升压药升级仅为可重复的记录代理。","",
      "同ICU类型×年代限制要求每个保留单元两臂各至少5人。临床情境限制排除手术、创伤、CVICU和PACU代理。这些均为固定偏倚审计，不是寻找显著亚组。","",
      "若任一治疗前伪结局的500次区间排除0，N09不能升级为绿色；若主要肾结局在两个限制分析中改变方向，也维持黄色或转红。","" ]
    MD_OUT.write_text("\n".join(lines),encoding="utf-8")
    for path in (CSV_OUT,JSON_OUT,MD_OUT):
        if not path.exists() or path.stat().st_size==0: raise RuntimeError(f"missing {path}")
    with CSV_OUT.open(encoding="utf-8-sig",newline="") as f: list(csv.DictReader(f))
    json.loads(JSON_OUT.read_text(encoding="utf-8"))


if __name__=="__main__": main()
