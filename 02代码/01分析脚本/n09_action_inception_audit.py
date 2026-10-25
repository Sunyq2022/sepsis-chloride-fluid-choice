"""Outcome-blind N09 action-inception identity audit."""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "03结果/21_FULL_ICU_CAUSAL_PROGRAM_20260907/N09/action_inception_audit_20260908"
BIAS_SCRIPT = ROOT / "02代码/01分析脚本/n09_bias_audit.py"
NS, LR = 225158, 225828


def load_bias():
    spec = importlib.util.spec_from_file_location("n09_bias_inception", BIAS_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


def create_alternative_actions(con) -> None:
    con.execute(
        f"""create temp table action_inception_proxy as
            select r.*,i.itemid,i.starttime t0,i.endtime,
                   i.originalamount action_amount,i.amount,i.originalamount,i.rate,i.originalrate,
                   i.storetime,i.statusdescription,i.ordercategoryname
            from n09_riskset r join lateral(
              select * from b.mimiciv_icu.inputevents i where i.stay_id=r.stay_id
                and i.starttime>=r.trigger_time and i.starttime<r.trigger_time+interval '4h'
                and i.itemid in({NS},{LR})
                and i.originalamount>=250
                and i.ordercategoryname='03-IV Fluid Bolus'
              order by i.starttime,i.itemid limit 1) i on true"""
    )
    con.execute(
        f"""create temp table action_rate_supported as
            select r.*,i.itemid,i.starttime t0,i.endtime,
                   i.originalamount action_amount,i.amount,i.originalamount,i.rate,i.originalrate,
                   i.storetime,i.statusdescription,i.ordercategoryname
            from n09_riskset r join lateral(
              select * from b.mimiciv_icu.inputevents i where i.stay_id=r.stay_id
                and i.starttime>=r.trigger_time and i.starttime<r.trigger_time+interval '4h'
                and i.itemid in({NS},{LR})
                and i.originalamount>=250 and i.originalrate>0
                and i.ordercategoryname='03-IV Fluid Bolus'
              order by i.starttime,i.itemid limit 1) i on true"""
    )


def profile_candidates(con) -> pd.DataFrame:
    return con.execute(
        f"""select case when i.itemid={LR} then 'LR' else 'NaCl' end arm,
          coalesce(i.statusdescription,'NULL') statusdescription,
          coalesce(i.ordercategoryname,'NULL') ordercategoryname,
          count(*) n_records,count(distinct i.stay_id) n_source_stays,
          sum(case when i.amount>0 then 1 else 0 end) positive_amount_records,
          sum(case when i.amount>=250 then 1 else 0 end) amount_ge250_records,
          sum(case when i.originalamount>=250 then 1 else 0 end) originalamount_ge250_records,
          sum(case when i.rate>0 then 1 else 0 end) positive_rate_records,
          sum(case when i.originalrate>0 then 1 else 0 end) positive_originalrate_records,
          sum(case when i.endtime>i.starttime then 1 else 0 end) positive_duration_records,
          sum(case when i.storetime<=i.starttime then 1 else 0 end) stored_by_start_records,
          median(i.amount) median_amount,
          median(i.originalamount) median_originalamount,
          median(datediff('second',i.starttime,i.endtime)/60.0) median_duration_minutes,
          median(datediff('second',i.starttime,i.storetime)/60.0) median_store_delay_minutes
        from n09_riskset r join b.mimiciv_icu.inputevents i
          on i.stay_id=r.stay_id and i.starttime>=r.trigger_time
         and i.starttime<r.trigger_time+interval '4h'
        where i.itemid in({NS},{LR})
        group by 1,2,3 order by 1,3,2"""
    ).df()


def summarize_definition(name: str, frame: pd.DataFrame) -> dict:
    return {
        "definition": name,
        "n": len(frame),
        "NaCl_n": int((frame.exposure == 0).sum()),
        "LR_n": int((frame.exposure == 1).sum()),
    }


def compare(left: pd.DataFrame, right: pd.DataFrame, label: str) -> dict:
    a = left[["stay_id", "exposure", "t0"]].rename(
        columns={"exposure": "frozen_arm", "t0": "frozen_t0"}
    )
    b = right[["stay_id", "exposure", "t0"]].rename(
        columns={"exposure": "alternative_arm", "t0": "alternative_t0"}
    )
    merged = a.merge(b, on="stay_id", how="outer", indicator=True)
    both = merged[merged._merge == "both"].copy()
    time_equal = pd.to_datetime(both.frozen_t0).eq(pd.to_datetime(both.alternative_t0))
    return {
        "comparison": label,
        "common_n": len(both),
        "same_arm_n": int(both.frozen_arm.eq(both.alternative_arm).sum()),
        "different_arm_n": int(both.frozen_arm.ne(both.alternative_arm).sum()),
        "same_t0_n": int(time_equal.sum()),
        "different_t0_n": int((~time_equal).sum()),
        "frozen_only_n": int((merged._merge == "left_only").sum()),
        "alternative_only_n": int((merged._merge == "right_only").sum()),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for path in OUT.iterdir():
        if path.is_file():
            path.unlink()
    bias = load_bias()
    canonical = bias.load_canonical_module()
    con = canonical.connection()
    con.execute("set threads=1")
    bias.create_actions(con)
    create_alternative_actions(con)
    profile = profile_candidates(con)
    frozen = bias.feature_frame(con, "action_completed").reset_index(drop=True)
    inception = bias.feature_frame(con, "action_inception_proxy").reset_index(drop=True)
    rate_supported = bias.feature_frame(con, "action_rate_supported").reset_index(drop=True)
    con.close()

    if len(frozen) != 706 or int((frozen.exposure == 1).sum()) != 259:
        raise RuntimeError("Frozen deterministic cohort did not reproduce 706/259")
    definitions = [
        summarize_definition("frozen_finishedrunning_and_amount_ge250", frozen),
        summarize_definition("inception_proxy_originalamount_ge250_no_future_status", inception),
        summarize_definition("rate_supported_originalamount_ge250_originalrate_positive", rate_supported),
    ]
    comparisons = [
        compare(frozen, inception, "frozen_vs_inception_proxy"),
        compare(frozen, rate_supported, "frozen_vs_rate_supported"),
    ]
    profile.to_csv(OUT / "candidate_record_profile.csv", index=False)
    pd.DataFrame(definitions).to_csv(OUT / "action_definitions.csv", index=False)
    pd.DataFrame(comparisons).to_csv(OUT / "definition_concordance.csv", index=False)

    bolus_profile = profile[profile.ordercategoryname.eq("03-IV Fluid Bolus")]
    bolus_records = int(bolus_profile.n_records.sum())
    bolus_nonfinished = int(
        bolus_profile.loc[~bolus_profile.statusdescription.eq("FinishedRunning"), "n_records"].sum()
    )
    bolus_originalrate_missing = int(
        bolus_profile.n_records.sum() - bolus_profile.positive_originalrate_records.sum()
    )
    bolus_stored_by_start = int(bolus_profile.stored_by_start_records.sum())
    inception_cmp = comparisons[0]
    rate_cmp = comparisons[1]
    audit = {
        "status": "COMPLETE_OUTCOME_BLIND_ACTION_INCEPTION_AUDIT",
        "frozen": definitions[0],
        "inception_proxy": definitions[1],
        "rate_supported": definitions[2],
        "frozen_vs_inception_proxy": inception_cmp,
        "frozen_vs_rate_supported": rate_cmp,
        "bolus_source": {
            "candidate_records": bolus_records,
            "non_FinishedRunning_records": bolus_nonfinished,
            "missing_or_nonpositive_originalrate_records": bolus_originalrate_missing,
            "stored_at_or_before_start_records": bolus_stored_by_start,
        },
        "earliest_supportable_action_definition": "First LR or 0.9% saline row in the IV Fluid Bolus category within 4 h, using starttime and originalamount >=250 mL, with an explicit NaCl-first tie-break; do not filter on eventual status or delivered amount.",
        "identity_limit": "This is an inputevents inception proxy, not a proven bedside start. The source has no separate start-confirmation flag. Every qualifying bolus row ultimately has FinishedRunning status, so the data cannot test whether a planned/cancelled bolus absent from inputevents should have entered the risk set. Requiring originalrate would drop valid one-minute bolus records and does not solve that problem.",
        "outcomes_queried": False,
    }
    (OUT / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# N09动作起点身份审计

## 结论

当前MIMIC-IV可支持的最早定义是：**高氯触发后4小时内，`inputevents`中首个标为IV Fluid Bolus的LR或0.9% NaCl记录，以`starttime`为起点、以`originalamount≥250 mL`界定计划剂量，并显式规定同刻NaCl优先。** 该定义不使用未来的`FinishedRunning`状态，也不使用最终`amount`。

但是它仍是**输液开始代理**，不是经过独立床旁确认的真实启动。`inputevents`没有单独的“已开始执行”标志，因而无法确认未进入该表的计划/取消医嘱；本轮不能声称完全解决动作身份。

## 实际数据结果

|定义|总N|NaCl|LR|
|---|---:|---:|---:|
|冻结：FinishedRunning＋最终amount≥250|{definitions[0]['n']}|{definitions[0]['NaCl_n']}|{definitions[0]['LR_n']}|
|起点代理：originalamount≥250，不使用未来状态|{definitions[1]['n']}|{definitions[1]['NaCl_n']}|{definitions[1]['LR_n']}|
|另要求originalrate>0|{definitions[2]['n']}|{definitions[2]['NaCl_n']}|{definitions[2]['LR_n']}|

冻结定义与起点代理共同{inception_cmp['common_n']}人，同一液体{inception_cmp['same_arm_n']}人、不同液体{inception_cmp['different_arm_n']}人；冻结独有{inception_cmp['frozen_only_n']}人，起点代理新增{inception_cmp['alternative_only_n']}人，零点不同{inception_cmp['different_t0_n']}人。

要求`originalrate>0`后，与冻结队列共同{rate_cmp['common_n']}人，冻结独有{rate_cmp['frozen_only_n']}人，新增{rate_cmp['alternative_only_n']}人。这个条件会排除部分一分钟bolus，不能因为看起来更“严格”就作为主定义。

## 源记录特征

- 风险集中共有{bolus_records}条LR/NaCl IV Fluid Bolus候选记录；非`FinishedRunning`记录为{bolus_nonfinished}条。
- {bolus_stored_by_start}/{bolus_records}条在记录时间上于开始时或此前写入；这支持`starttime`作为记录起点，但不能证明床旁液体已经进入患者体内。
- {bolus_originalrate_missing}条缺失或没有正的`originalrate`。Bolus通常以一分钟记录，不能机械要求连续输注速率。
- 各状态、类别、剂量、原始剂量、原始速率及存储延迟的完整聚合见`candidate_record_profile.csv`。

## 决策

起点代理与冻结706人完全一致时，可以把主研究的时间零点改写为“首个符合原始计划剂量的inputevents bolus开始记录”，并停止在资格定义中使用`FinishedRunning`和最终`amount`。这改善了时间顺序表达，但**不提供独立的实际给药确认**。因此N09仍需保留“inputevents启动代理可能误分、无法观察未入表的取消医嘱”这一核心局限。本审计未读取任何肾脏或死亡结局。
"""
    (OUT / "report.md").write_text(report, encoding="utf-8")

    expected = {
        "candidate_record_profile.csv", "action_definitions.csv",
        "definition_concordance.csv", "audit.json", "report.md",
    }
    actual = {p.name for p in OUT.iterdir() if p.is_file()}
    if actual != expected:
        raise RuntimeError(f"Unexpected output files: {sorted(actual)}")
    json.loads((OUT / "audit.json").read_text(encoding="utf-8"))
    for name in ("candidate_record_profile.csv", "action_definitions.csv", "definition_concordance.csv"):
        with (OUT / name).open("r", encoding="utf-8-sig", newline="") as handle:
            list(csv.DictReader(handle))
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
