"""Outcome-blind exposure and fluid-switching audit for N09.

This script reuses the frozen N09 eligibility/action logic but never queries
post-time-zero kidney outcomes or mortality. Only aggregate files are written.
"""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
CANONICAL_SCRIPT = ROOT / "02代码/01分析脚本/canonical_rerun_n02_n07_n09.py"
OUT = ROOT / "03结果/21_FULL_ICU_CAUSAL_PROGRAM_20260907/N09/bias_audit_20260908"
NS = 225158
LR = 225828
NUMERIC = [
    "age", "chloride", "prior_chloride", "ns_ml", "lr_ml", "action_amount",
    "sofa_score", "creat", "low7", "aki_stage_creat", "sodium",
    "bicarbonate", "lactate", "potassium", "hemoglobin", "platelet",
    "map_mean", "map_min", "heart_rate", "pressor", "ventilation",
]
CATEGORICAL = ["gender", "first_careunit", "anchor_year_group"]


RISKSET_SQL = r"""
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
      and i.amount>0 and i.statusdescription='FinishedRunning'
      and i.ordercategoryname in('03-IV Fluid Bolus','02-Fluids (Crystalloids)')) f on true
)
select * from prior_fluid where ns_ml>=500
"""


FEATURE_SQL = r"""
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
from {action_table} a
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
where coalesce(r.aki_stage_creat,0)<3
  and not exists(select 1 from b.mimiciv_derived.rrt x
                 where x.stay_id=a.stay_id and x.charttime<=a.t0 and x.dialysis_present=1)
"""


def load_canonical_module():
    spec = importlib.util.spec_from_file_location("n09_canonical", CANONICAL_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path.name}")
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def create_actions(con) -> None:
    con.execute("create temp table n09_riskset as " + RISKSET_SQL)
    con.execute(
        f"""create temp table action_completed as
            select r.*,i.itemid,i.starttime t0,i.endtime,i.amount action_amount,
                   i.statusdescription,i.ordercategoryname
            from n09_riskset r join lateral(
              select * from b.mimiciv_icu.inputevents i where i.stay_id=r.stay_id
                and i.starttime>=r.trigger_time and i.starttime<r.trigger_time+interval '4h'
                and i.itemid in({NS},{LR}) and i.amount>=250
                and i.statusdescription='FinishedRunning'
                and i.ordercategoryname='03-IV Fluid Bolus'
              order by i.starttime,i.itemid limit 1) i on true"""
    )
    con.execute(
        f"""create temp table action_initiated as
            select r.*,i.itemid,i.starttime t0,i.endtime,i.amount action_amount,
                   i.statusdescription,i.ordercategoryname
            from n09_riskset r join lateral(
              select * from b.mimiciv_icu.inputevents i where i.stay_id=r.stay_id
                and i.starttime>=r.trigger_time and i.starttime<r.trigger_time+interval '4h'
                and i.itemid in({NS},{LR}) and i.amount>=250
                and i.statusdescription is distinct from 'Rewritten'
                and i.ordercategoryname='03-IV Fluid Bolus'
              order by i.starttime,i.itemid limit 1) i on true"""
    )
    con.execute(
        f"""create temp table action_completed_lr_tie as
            select r.*,i.itemid,i.starttime t0,i.endtime,i.amount action_amount,
                   i.statusdescription,i.ordercategoryname
            from n09_riskset r join lateral(
              select * from b.mimiciv_icu.inputevents i where i.stay_id=r.stay_id
                and i.starttime>=r.trigger_time and i.starttime<r.trigger_time+interval '4h'
                and i.itemid in({NS},{LR}) and i.amount>=250
                and i.statusdescription='FinishedRunning'
                and i.ordercategoryname='03-IV Fluid Bolus'
              order by i.starttime,i.itemid desc limit 1) i on true"""
    )


def feature_frame(con, table: str) -> pd.DataFrame:
    count = con.execute(f"select count(*) from {table}").fetchone()[0]
    parts: list[pd.DataFrame] = []
    query = FEATURE_SQL.format(action_table="action_batch")
    for offset in range(0, count, 32):
        con.execute(
            f"create or replace temp table action_batch as "
            f"select * from {table} order by subject_id,stay_id limit 32 offset {offset}"
        )
        parts.append(con.execute(query).df())
    return pd.concat(parts, ignore_index=True)


def arm_name(value: int) -> str:
    return "LR" if value == 1 else "NaCl"


def ess(weights: np.ndarray) -> float:
    return float(weights.sum() ** 2 / np.square(weights).sum()) if len(weights) else 0.0


def summarize_point_definitions(completed: pd.DataFrame, initiated: pd.DataFrame) -> list[dict]:
    rows: list[dict] = []
    for name, frame in (("first_completed_bolus", completed), ("first_recorded_initiation", initiated)):
        rows.append(
            {
                "comparison": name,
                "n": len(frame),
                "NaCl_n": int((frame.exposure == 0).sum()),
                "LR_n": int((frame.exposure == 1).sum()),
                "metric": "eligible action definition",
                "value": "",
            }
        )
    left = completed[["stay_id", "exposure", "t0"]].rename(
        columns={"exposure": "completed_arm", "t0": "completed_t0"}
    )
    right = initiated[["stay_id", "exposure", "t0"]].rename(
        columns={"exposure": "initiated_arm", "t0": "initiated_t0"}
    )
    matched = left.merge(right, on="stay_id", how="inner")
    same_arm = matched.completed_arm.eq(matched.initiated_arm)
    delta = (pd.to_datetime(matched.completed_t0) - pd.to_datetime(matched.initiated_t0)).dt.total_seconds() / 3600
    rows.extend(
        [
            {"comparison": "definition_concordance", "n": len(matched), "NaCl_n": "", "LR_n": "", "metric": "same_arm_n", "value": int(same_arm.sum())},
            {"comparison": "definition_concordance", "n": len(matched), "NaCl_n": "", "LR_n": "", "metric": "discordant_arm_n", "value": int((~same_arm).sum())},
            {"comparison": "definition_concordance", "n": len(matched), "NaCl_n": "", "LR_n": "", "metric": "completed_t0_later_than_initiation_n", "value": int((delta > 0).sum())},
            {"comparison": "definition_concordance", "n": len(matched), "NaCl_n": "", "LR_n": "", "metric": "median_t0_delay_hours", "value": float(delta.median())},
        ]
    )
    return rows


def fluid_records(con, cohort: pd.DataFrame, scope: str) -> pd.DataFrame:
    con.register("audit_cohort", cohort[["stay_id", "t0", "exposure"]])
    restriction = (
        "i.ordercategoryname='03-IV Fluid Bolus' and i.amount>=250"
        if scope == "bolus_250ml"
        else "i.ordercategoryname in('03-IV Fluid Bolus','02-Fluids (Crystalloids)') and i.amount>0"
    )
    result = con.execute(
        f"""select c.stay_id,c.t0,c.exposure,i.itemid,i.starttime,i.endtime,i.amount,
                   i.ordercategoryname,
                   datediff('second',c.t0,i.starttime)/3600.0 hours_after_t0
              from audit_cohort c join b.mimiciv_icu.inputevents i
                on i.stay_id=c.stay_id and i.starttime>=c.t0
               and i.starttime<c.t0+interval '24h'
             where i.itemid in({NS},{LR}) and i.statusdescription='FinishedRunning'
               and {restriction}"""
    ).df()
    con.unregister("audit_cohort")
    return result


def consistency_rows(
    cohort: pd.DataFrame, records: pd.DataFrame, weights: np.ndarray, scope: str
) -> tuple[list[dict], list[dict]]:
    base = cohort[["stay_id", "exposure", "t0"]].copy()
    base["weight"] = weights
    summary_rows: list[dict] = []
    switch_rows: list[dict] = []
    for horizon in (6, 24):
        window = records[records.hours_after_t0 < horizon].copy()
        for arm in (0, 1):
            group = base[base.exposure == arm].copy()
            arm_records = window[window.stay_id.isin(group.stay_id)].copy()
            assigned_item = LR if arm == 1 else NS
            opposite_item = NS if arm == 1 else LR
            by_id = arm_records.groupby("stay_id")
            opposite = by_id.itemid.apply(lambda x: bool((x == opposite_item).any()))
            same_after = arm_records[
                (arm_records.itemid == assigned_item) & (arm_records.hours_after_t0 > 0)
            ].groupby("stay_id").size().gt(0)
            any_after = arm_records[arm_records.hours_after_t0 > 0].groupby("stay_id").size().gt(0)
            volumes = arm_records.assign(
                assigned_volume=np.where(arm_records.itemid == assigned_item, arm_records.amount, 0.0)
            ).groupby("stay_id")[["assigned_volume", "amount"]].sum()
            majority80 = (volumes.assigned_volume / volumes.amount).ge(0.80)
            group["crossover"] = group.stay_id.map(opposite).fillna(False).astype(bool)
            group["strict_consistent"] = ~group.crossover
            group["majority80"] = group.stay_id.map(majority80).fillna(False).astype(bool)
            group["same_additional"] = group.stay_id.map(same_after).eq(True)
            group["no_additional"] = ~group.stay_id.map(any_after).eq(True)
            retained = group[group.strict_consistent]
            weighted_cross = np.average(group.crossover, weights=group.weight)
            weighted_consistent = np.average(group.strict_consistent, weights=group.weight)
            summary_rows.append(
                {
                    "record_scope": scope,
                    "horizon_hours": horizon,
                    "point_arm": arm_name(arm),
                    "point_n": len(group),
                    "crossover_n": int(group.crossover.sum()),
                    "crossover_pct": round(100 * group.crossover.mean(), 2),
                    "strict_consistent_n": int(group.strict_consistent.sum()),
                    "strict_consistent_pct": round(100 * group.strict_consistent.mean(), 2),
                    "majority80_n": int(group.majority80.sum()),
                    "majority80_pct": round(100 * group.majority80.mean(), 2),
                    "same_fluid_additional_n": int(group.same_additional.sum()),
                    "no_additional_fluid_n": int(group.no_additional.sum()),
                    "ATO_weighted_crossover_pct": round(100 * weighted_cross, 2),
                    "ATO_weighted_strict_consistent_pct": round(100 * weighted_consistent, 2),
                    "strict_consistent_subset_ESS": round(ess(retained.weight.to_numpy()), 2),
                }
            )
            switches = arm_records[arm_records.itemid == opposite_item].groupby("stay_id").hours_after_t0.min()
            if len(switches):
                switch_rows.append(
                    {
                        "record_scope": scope,
                        "horizon_hours": horizon,
                        "transition": f"{arm_name(arm)} to {arm_name(1-arm)}",
                        "switch_n": len(switches),
                        "median_hours": round(float(switches.median()), 2),
                        "q25_hours": round(float(switches.quantile(0.25)), 2),
                        "q75_hours": round(float(switches.quantile(0.75)), 2),
                    }
                )
    return summary_rows, switch_rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for path in OUT.iterdir():
        if path.is_file():
            path.unlink()
    module = load_canonical_module()
    con = module.connection()
    con.execute("set threads=1")
    create_actions(con)
    completed = feature_frame(con, "action_completed")
    initiated = feature_frame(con, "action_initiated")
    completed_lr_tie = feature_frame(con, "action_completed_lr_tie")
    if len(completed) != 706 or int((completed.exposure == 1).sum()) != 259:
        raise RuntimeError(
            "Canonical point-treatment cohort was not reproduced: "
            f"n={len(completed)}, LR={int((completed.exposure == 1).sum())}"
        )

    fitted, _, _, ps, weights = module.ato_fit(completed, NUMERIC, CATEGORICAL)
    if not fitted.index.equals(completed.reset_index(drop=True).index):
        raise RuntimeError("Unexpected ATO row expansion")
    completed = completed.reset_index(drop=True)
    completed["ato_weight"] = weights
    _, diagnostic_with_amount = module.ato_diagnostics(completed, NUMERIC, CATEGORICAL)
    numeric_without_action_amount = [x for x in NUMERIC if x != "action_amount"]
    _, diagnostic_without_amount = module.ato_diagnostics(
        completed, numeric_without_action_amount, CATEGORICAL
    )

    definition_rows = summarize_point_definitions(completed, initiated)
    definition_rows.extend(
        [
            {
                "comparison": "completed_bolus_tie_break_LR_first",
                "n": len(completed_lr_tie),
                "NaCl_n": int((completed_lr_tie.exposure == 0).sum()),
                "LR_n": int((completed_lr_tie.exposure == 1).sum()),
                "metric": "eligible action definition",
                "value": "",
            },
            {
                "comparison": "propensity_model",
                "n": len(completed),
                "NaCl_n": "",
                "LR_n": "",
                "metric": "canonical_includes_post_assignment_action_amount",
                "value": 1,
            },
        ]
    )
    all_consistency: list[dict] = []
    all_switches: list[dict] = []
    for scope in ("bolus_250ml", "all_completed_lr_ns_fluid"):
        records = fluid_records(con, completed, scope)
        consistency, switches = consistency_rows(completed, records, weights, scope)
        all_consistency.extend(consistency)
        all_switches.extend(switches)
        if scope == "bolus_250ml":
            dual_at_t0 = (
                records[np.isclose(records.hours_after_t0, 0.0)]
                .groupby("stay_id").itemid.nunique().gt(1)
            )
            same_time_dual_fluid_n = int(dual_at_t0.sum())

    definition_rows.append(
        {
            "comparison": "same_time_tie_audit",
            "n": len(completed),
            "NaCl_n": "",
            "LR_n": "",
            "metric": "same_time_LR_and_NaCl_n",
            "value": same_time_dual_fluid_n,
        }
    )

    duration = (
        pd.to_datetime(completed.endtime) - pd.to_datetime(completed.t0)
    ).dt.total_seconds() / 3600
    index_summary = {}
    for arm in (0, 1):
        values = duration[completed.exposure == arm]
        index_summary[arm_name(arm)] = {
            "n": int(len(values)),
            "median_completion_hours": float(values.median()),
            "q25_completion_hours": float(values.quantile(0.25)),
            "q75_completion_hours": float(values.quantile(0.75)),
            "completion_after_6h_n": int((values > 6).sum()),
            "completion_after_24h_n": int((values > 24).sum()),
        }

    write_csv(OUT / "action_definition.csv", definition_rows)
    write_csv(OUT / "fluid_consistency.csv", all_consistency)
    write_csv(OUT / "switch_timing.csv", all_switches)
    audit = {
        "status": "COMPLETE_OUTCOME_BLIND_BIAS_AUDIT",
        "canonical_point_cohort": {
            "n": len(completed),
            "NaCl_n": int((completed.exposure == 0).sum()),
            "LR_n": int((completed.exposure == 1).sum()),
        },
        "point_propensity_support": {
            "NaCl_ESS": ess(weights[completed.exposure.to_numpy() == 0]),
            "LR_ESS": ess(weights[completed.exposure.to_numpy() == 1]),
            "NaCl_ps_quantiles": np.quantile(ps[completed.exposure.to_numpy() == 0], [0, 0.1, 0.5, 0.9, 1]).tolist(),
            "LR_ps_quantiles": np.quantile(ps[completed.exposure.to_numpy() == 1], [0, 0.1, 0.5, 0.9, 1]).tolist(),
        },
        "tie_break_audit": {
            "archived_frozen_result": {"n": 706, "NaCl_n": 447, "LR_n": 259},
            "deterministic_NaCl_first": {"n": len(completed), "NaCl_n": int((completed.exposure == 0).sum()), "LR_n": int((completed.exposure == 1).sum())},
            "deterministic_LR_first": {"n": len(completed_lr_tie), "NaCl_n": int((completed_lr_tie.exposure == 0).sum()), "LR_n": int((completed_lr_tie.exposure == 1).sum())},
            "same_time_dual_fluid_n": same_time_dual_fluid_n,
            "interpretation": "The archived 259-LR result is retained. The original SQL ordered only by start time; one dual-fluid tie can move one patient between arms. No effect estimate was rerun in this audit.",
        },
        "propensity_model_audit": {
            "canonical_with_action_amount": diagnostic_with_amount,
            "baseline_only_without_action_amount": diagnostic_without_amount,
            "decision": "action_amount is part of delivered treatment and must not be in the baseline propensity model when the treatment is fluid type; the frozen result is preserved but requires re-estimation under a separately frozen baseline-only model before causal publication claims.",
        },
        "index_completion_duration": index_summary,
        "predefined_definitions": {
            "point_primary": "First >=250 mL LR or 0.9% saline record starting within 4 h after the chloride trigger, IV Fluid Bolus category, retrospectively marked FinishedRunning.",
            "point_initiation_audit": "First >=250 mL LR or 0.9% saline bolus record within 4 h, any status other than Rewritten.",
            "subsequent_primary_scope": "All completed LR/0.9% saline records in bolus or crystalloid categories with positive amount.",
            "subsequent_narrow_scope": "Completed IV Fluid Bolus records >=250 mL.",
            "strict_consistency": "No opposite fluid record from point time through the stated horizon.",
            "majority_consistency": "At least 80% of recorded LR plus saline volume matches the point treatment.",
        },
        "interpretation": "Subsequent consistency is descriptive. Restricting to natural adherers would condition on post-treatment care and is not a valid sensitivity effect estimate without a separately specified per-protocol emulation and censoring model.",
        "outcomes_queried": False,
    }
    (OUT / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")

    def table_lines(rows: list[dict], scope: str) -> str:
        chosen = [r for r in rows if r["record_scope"] == scope]
        lines = ["|窗口|起始液体|N|交叉换液|严格一致|ATO加权交叉|一致子集ESS|", "|---:|---|---:|---:|---:|---:|---:|"]
        for row in chosen:
            lines.append(
                f"|{row['horizon_hours']}h|{row['point_arm']}|{row['point_n']}|"
                f"{row['crossover_n']}（{row['crossover_pct']:.2f}%）|"
                f"{row['strict_consistent_n']}（{row['strict_consistent_pct']:.2f}%）|"
                f"{row['ATO_weighted_crossover_pct']:.2f}%|{row['strict_consistent_subset_ESS']:.2f}|"
            )
        return "\n".join(lines)

    point = {row["comparison"]: row for row in definition_rows if row["metric"] == "eligible action definition"}
    concordance = {row["metric"]: row["value"] for row in definition_rows if row["comparison"] == "definition_concordance"}
    report = f"""# N09暴露操作化与后续液体切换审计

## 结论

本审计完整复现原点治疗队列 **706人（LR 259，NaCl 447）**，且未查询任何肾脏或死亡结局。点治疗的共同支持仍可用，但它代表的是**下一次已完成推注的选择**，不能直接解释为持续6小时或24小时的限氯策略。

## 点治疗定义

- 原定义“首个完成的推注”：N={point['first_completed_bolus']['n']}，LR={point['first_completed_bolus']['LR_n']}，NaCl={point['first_completed_bolus']['NaCl_n']}。
- 预设启动记录审计：N={point['first_recorded_initiation']['n']}，LR={point['first_recorded_initiation']['LR_n']}，NaCl={point['first_recorded_initiation']['NaCl_n']}。
- 两定义共同患者{concordance['same_arm_n'] + concordance['discordant_arm_n']}人；液体类别一致{concordance['same_arm_n']}人，不一致{concordance['discordant_arm_n']}人；完成定义的零点晚于首次启动记录者{concordance['completed_t0_later_than_initiation_n']}人。
- 同一最早时刻同时存在LR和NaCl记录者{same_time_dual_fluid_n}人。原SQL只按开始时间排序，因此冻结存档为LR 259，而无稳定tie-break的重建曾得到LR 260；显式NaCl优先为LR 259，显式LR优先为LR {int((completed_lr_tie.exposure == 1).sum())}。本审计保留原259，不静默覆盖。

`FinishedRunning`在输注结束后才能确定，因此原队列存在按完成状态选择的风险。该问题不是靠后续敏感性结果显著与否决定，而应在正文中明确限制因果解释。

原倾向模型还纳入了`action_amount`。该变量是所选液体的给药/完成剂量，位于分配之后，不应作为“LR还是NaCl”选择的基线混杂因素。纳入时ESS为NaCl {diagnostic_with_amount['ess_control']:.1f}、LR {diagnostic_with_amount['ess_intervention']:.1f}，最大SMD {diagnostic_with_amount['max_abs_smd']:.3f}；删除后ESS为NaCl {diagnostic_without_amount['ess_control']:.1f}、LR {diagnostic_without_amount['ess_intervention']:.1f}，最大SMD {diagnostic_without_amount['max_abs_smd']:.3f}。这里只报告设计诊断，没有重估结局。

## 后续所有LR/NaCl液体记录：主要一致性审计

{table_lines(all_consistency, 'all_completed_lr_ns_fluid')}

## 仅≥250 mL推注：较窄审计

{table_lines(all_consistency, 'bolus_250ml')}

“严格一致”只表示未记录相反液体；没有后续液体者也归为一致。表中另保留了同类追加液体、无追加液体和80%体积一致性，见`fluid_consistency.csv`。

## 对研究解释的影响

1. 现有N09适合表述为**高氯触发后的下一次LR与生理盐水推注选择**，不适合写成6小时或24小时持续液体策略。
2. 6小时/24小时换液比例用于判断自然治疗交叉。不能把“始终一致者”直接筛出来重估效应，因为这会按治疗后的临床过程选择患者。
3. 若以后要估计持续策略，应预先建立克隆—人工删失—加权或序贯试验，纳入每个时间段的血流动力学、补液需求和新化验；本轮没有制造这种因果估计。
4. 本次固定了两种记录范围和两个窗口，没有读取结局，也没有为显著性挑选定义。
5. 发表前需要另行冻结“显式tie-break＋仅t0前协变量”的点治疗模型，再重估原主要和次要结局；这一步不能用本审计的换液结果决定是否执行。
"""
    (OUT / "report.md").write_text(report, encoding="utf-8")
    con.close()

    expected = {"action_definition.csv", "fluid_consistency.csv", "switch_timing.csv", "audit.json", "report.md"}
    actual = {p.name for p in OUT.iterdir() if p.is_file()}
    if actual != expected:
        raise RuntimeError(f"Unexpected output files: {sorted(actual)}")
    json.loads((OUT / "audit.json").read_text(encoding="utf-8"))
    for name in ("action_definition.csv", "fluid_consistency.csv", "switch_timing.csv"):
        with (OUT / name).open("r", encoding="utf-8-sig", newline="") as handle:
            list(csv.DictReader(handle))
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
