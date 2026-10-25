"""Outcome-blind pre-time-zero confounder augmentation for N09."""
from __future__ import annotations

import importlib.util
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


ROOT = Path(r"D:\2026project\04AKI课题基础\05因果推断")
ACTION_AUDIT = ROOT / "02代码" / "01分析脚本" / "n09_bias_audit.py"
OUT = ROOT / "03结果" / "21_FULL_ICU_CAUSAL_PROGRAM_20260907" / "N09" / "pre_t0_confounder_audit_20260908"
EXPECTED = {"diagnostics.json", "model_comparison.csv", "balance.csv", "missingness.csv", "feature_definitions.csv", "report.md"}
NS = 225158
LR = 225828
URINE_ITEMS = (226557, 226558, 226559, 226560, 226561, 226563, 226564, 226565, 226584, 226627, 226631, 226713)
SEED = 20260908


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def cv_auc(df: pd.DataFrame, numeric: list[str], categorical: list[str]) -> float:
    prep = ColumnTransformer(
        [
            ("num", make_pipeline(SimpleImputer(strategy="median", add_indicator=True), StandardScaler()), numeric),
            ("cat", make_pipeline(SimpleImputer(strategy="most_frequent"), OneHotEncoder(handle_unknown="ignore")), categorical),
        ]
    )
    model = make_pipeline(prep, LogisticRegression(C=0.5, max_iter=3000))
    folds = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    prediction = cross_val_predict(model, df[numeric + categorical], df["exposure"], cv=folds, method="predict_proba")[:, 1]
    return float(roc_auc_score(df["exposure"], prediction))


def extra_features(con, cohort: pd.DataFrame) -> pd.DataFrame:
    urine = ",".join(map(str, URINE_ITEMS))
    query = f"""select a.stay_id,
      u.urine6_ml,u.urine24_ml,u.urine6_records,u.urine24_records,u.urine6_hours,u.urine24_hours,
      io.input6_ml,io.input24_ml,io.known_ns24_ml,io.known_lr24_ml,oo.output6_ml,oo.output24_ml,
      pb.prior_ns_amount,pb.prior_ns_end,
      pur.pre_urine_ml,pur.post_urine_ml,pur.pre_urine_records,pur.post_urine_records,
      datediff('second',greatest(a.intime,pb.prior_ns_end-interval '4h'),pb.prior_ns_end)/3600.0 prior_response_pre_hours,
      datediff('second',pb.prior_ns_end,least(a.t0,pb.prior_ns_end+interval '4h'))/3600.0 prior_response_post_hours,
      cr.creat_first24,cr.creat_last24,cr.creat_n24,cr.creat_span_hours,
      cl.chloride_first48,cl.chloride_last48,cl.chloride_n48,cl.chloride_span_hours,
      la.lactate_first6,la.lactate_last6,la.lactate_n6,la.lactate_span_hours
    from conf_batch a
    left join lateral(
      select sum(value) filter(where charttime>a.t0-interval '6h') urine6_ml,
             sum(value) urine24_ml,
             count(*) filter(where charttime>a.t0-interval '6h') urine6_records,
             count(*) urine24_records,
             count(distinct date_trunc('hour',charttime)) filter(where charttime>a.t0-interval '6h') urine6_hours,
             count(distinct date_trunc('hour',charttime)) urine24_hours
      from b.mimiciv_icu.outputevents o where o.stay_id=a.stay_id and o.itemid in({urine})
        and o.value between 0 and 5000 and lower(coalesce(o.valueuom,''))='ml'
        and o.charttime<=a.t0 and coalesce(o.storetime,o.charttime)<=a.t0 and o.charttime>a.t0-interval '24h'
    ) u on true
    left join lateral(
      select sum(amount) filter(where endtime>a.t0-interval '6h') input6_ml,
             sum(amount) input24_ml,
             sum(amount) filter(where itemid={NS}) known_ns24_ml,
             sum(amount) filter(where itemid={LR}) known_lr24_ml
      from b.mimiciv_icu.inputevents i where i.stay_id=a.stay_id and i.amount>0 and lower(coalesce(i.amountuom,''))='ml'
        and i.statusdescription='FinishedRunning' and i.endtime<=a.t0 and coalesce(i.storetime,i.endtime)<=a.t0
        and i.endtime>a.t0-interval '24h'
    ) io on true
    left join lateral(
      select sum(value) filter(where charttime>a.t0-interval '6h') output6_ml,
             sum(value) output24_ml
      from b.mimiciv_icu.outputevents o where o.stay_id=a.stay_id and o.value>0
        and lower(coalesce(o.valueuom,''))='ml' and o.charttime<=a.t0 and coalesce(o.storetime,o.charttime)<=a.t0
        and o.charttime>a.t0-interval '24h'
    ) oo on true
    left join lateral(
      select amount prior_ns_amount,endtime prior_ns_end
      from b.mimiciv_icu.inputevents i where i.stay_id=a.stay_id and i.itemid={NS}
        and i.amount>=250 and i.statusdescription='FinishedRunning' and i.ordercategoryname='03-IV Fluid Bolus'
        and i.endtime<=a.trigger_time and coalesce(i.storetime,i.endtime)<=a.t0
        and i.endtime>a.trigger_time-interval '24h'
      order by i.endtime desc,i.orderid desc limit 1
    ) pb on true
    left join lateral(
      select sum(value) filter(where charttime>greatest(a.intime,pb.prior_ns_end-interval '4h')
                                and charttime<=pb.prior_ns_end) pre_urine_ml,
             sum(value) filter(where charttime>pb.prior_ns_end
                                and charttime<=least(a.t0,pb.prior_ns_end+interval '4h')) post_urine_ml,
             count(*) filter(where charttime>greatest(a.intime,pb.prior_ns_end-interval '4h')
                              and charttime<=pb.prior_ns_end) pre_urine_records,
             count(*) filter(where charttime>pb.prior_ns_end
                              and charttime<=least(a.t0,pb.prior_ns_end+interval '4h')) post_urine_records
      from b.mimiciv_icu.outputevents o where o.stay_id=a.stay_id and o.itemid in({urine})
        and o.value between 0 and 5000 and lower(coalesce(o.valueuom,''))='ml'
        and pb.prior_ns_end is not null and o.charttime<=a.t0 and coalesce(o.storetime,o.charttime)<=a.t0
    ) pur on true
    left join lateral(
      select arg_min(valuenum,greatest(charttime,storetime)) creat_first24,
             arg_max(valuenum,greatest(charttime,storetime)) creat_last24,
             count(*) creat_n24,
             datediff('second',min(greatest(charttime,storetime)),max(greatest(charttime,storetime)))/3600.0 creat_span_hours
      from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.itemid=50912 and l.valuenum between 0.1 and 30
        and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '24h'
    ) cr on true
    left join lateral(
      select arg_min(valuenum,greatest(charttime,storetime)) chloride_first48,
             arg_max(valuenum,greatest(charttime,storetime)) chloride_last48,
             count(*) chloride_n48,
             datediff('second',min(greatest(charttime,storetime)),max(greatest(charttime,storetime)))/3600.0 chloride_span_hours
      from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.itemid=50902 and l.valuenum between 70 and 150
        and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '48h'
    ) cl on true
    left join lateral(
      select arg_min(valuenum,greatest(charttime,storetime)) lactate_first6,
             arg_max(valuenum,greatest(charttime,storetime)) lactate_last6,
             count(*) lactate_n6,
             datediff('second',min(greatest(charttime,storetime)),max(greatest(charttime,storetime)))/3600.0 lactate_span_hours
      from b.mimiciv_hosp.labevents l where l.hadm_id=a.hadm_id and l.itemid=50813 and l.valuenum between 0.1 and 40
        and greatest(l.charttime,l.storetime)<=a.t0 and l.charttime>a.t0-interval '6h'
    ) la on true"""
    source = cohort.copy()
    con.register("conf_source", source)
    parts: list[pd.DataFrame] = []
    for offset in range(0, len(source), 24):
        con.execute(
            "create or replace temp table conf_batch as "
            f"select * from conf_source order by subject_id,stay_id limit 24 offset {offset}"
        )
        parts.append(con.execute(query).df())
    con.unregister("conf_source")
    return pd.concat(parts, ignore_index=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    unexpected = {path.name for path in OUT.iterdir() if path.is_file()} - EXPECTED
    if unexpected:
        raise RuntimeError(f"Unexpected files already present: {sorted(unexpected)}")

    audit = load_module(ACTION_AUDIT, "n09_action_audit")
    canonical = audit.load_canonical_module()
    con = canonical.connection()
    con.execute("set threads=1")
    con.execute("create temp table n09_riskset as " + audit.RISKSET_SQL)
    con.execute(
        f"""create temp table action_fixed as
        select r.*,i.itemid,i.starttime t0,i.endtime,i.amount action_amount,
               i.statusdescription,i.ordercategoryname
        from n09_riskset r join lateral(
          select * from b.mimiciv_icu.inputevents i where i.stay_id=r.stay_id
            and i.starttime>=r.trigger_time and i.starttime<r.trigger_time+interval '4h'
            and i.itemid in({NS},{LR}) and i.amount>=250
            and i.statusdescription='FinishedRunning' and i.ordercategoryname='03-IV Fluid Bolus'
          order by i.starttime,case when i.itemid={NS} then 0 else 1 end,i.endtime,i.orderid limit 1
        ) i on true"""
    )
    cohort = audit.feature_frame(con, "action_fixed").reset_index(drop=True)
    extras = extra_features(con, cohort)
    con.close()
    cohort = cohort.merge(extras, on="stay_id", how="left", validate="one_to_one")
    # Inputevents is the complete recorded administration table for this ICU
    # source. No matching recorded LR/NS in a fixed lookback is a recorded zero,
    # unlike absent urine or laboratory measurements, which remain missing.
    cohort[["known_ns24_ml", "known_lr24_ml"]] = cohort[["known_ns24_ml", "known_lr24_ml"]].fillna(0)
    cohort["net6_ml"] = cohort["input6_ml"] - cohort["output6_ml"]
    cohort["net24_ml"] = cohort["input24_ml"] - cohort["output24_ml"]
    cohort["hours_since_prior_ns"] = (
        pd.to_datetime(cohort["t0"]) - pd.to_datetime(cohort["prior_ns_end"])
    ).dt.total_seconds() / 3600
    cohort["pre_urine_rate"] = cohort["pre_urine_ml"] / cohort["prior_response_pre_hours"].where(cohort["prior_response_pre_hours"] >= 1)
    cohort["post_urine_rate"] = cohort["post_urine_ml"] / cohort["prior_response_post_hours"].where(cohort["prior_response_post_hours"] >= 1)
    cohort["prior_urine_response_delta"] = cohort["post_urine_rate"] - cohort["pre_urine_rate"]
    cohort["creat_delta24"] = cohort["creat_last24"] - cohort["creat_first24"]
    cohort["chloride_delta48"] = cohort["chloride_last48"] - cohort["chloride_first48"]
    cohort["lactate_delta6"] = cohort["lactate_last6"] - cohort["lactate_first6"]

    if len(cohort) != 706 or int(cohort["exposure"].sum()) != 259:
        raise RuntimeError(f"Fixed cohort mismatch: N={len(cohort)}, LR={int(cohort['exposure'].sum())}")

    base_numeric = [
        "age", "chloride", "prior_chloride", "ns_ml", "lr_ml", "sofa_score",
        "creat", "low7", "aki_stage_creat", "sodium", "bicarbonate", "lactate",
        "potassium", "hemoglobin", "platelet", "map_mean", "map_min", "heart_rate",
        "pressor", "ventilation",
    ]
    categorical = ["gender", "first_careunit", "anchor_year_group"]
    augmented = [
        "urine6_ml", "urine24_ml", "urine6_records", "urine24_records", "urine6_hours", "urine24_hours",
        "input6_ml", "input24_ml", "known_ns24_ml", "known_lr24_ml", "output6_ml", "output24_ml", "net6_ml", "net24_ml",
        "prior_ns_amount", "hours_since_prior_ns", "pre_urine_rate", "post_urine_rate", "prior_urine_response_delta",
        "creat_delta24", "creat_n24", "creat_span_hours", "chloride_delta48", "chloride_n48", "chloride_span_hours",
        "lactate_delta6", "lactate_n6", "lactate_span_hours",
    ]
    parsimonious = [
        "urine6_ml", "urine24_ml", "urine24_hours", "net6_ml", "net24_ml",
        "known_ns24_ml", "known_lr24_ml", "hours_since_prior_ns",
        "prior_urine_response_delta", "creat_delta24", "chloride_delta48", "lactate_delta6",
    ]
    specs = {
        "legal_original_without_action_amount": base_numeric,
        "pre_t0_parsimonious": base_numeric + parsimonious,
        "pre_t0_augmented": base_numeric + augmented,
    }
    comparisons: list[dict[str, object]] = []
    balances: list[pd.DataFrame] = []
    for name, numeric in specs.items():
        _, _, _, ps, weights = canonical.ato_fit(cohort, numeric, categorical)
        balance, diag = canonical.ato_diagnostics(cohort, numeric, categorical)
        y = cohort["exposure"].to_numpy(int)
        lower = max(float(np.quantile(ps[y == 0], .01)), float(np.quantile(ps[y == 1], .01)))
        upper = min(float(np.quantile(ps[y == 0], .99)), float(np.quantile(ps[y == 1], .99)))
        comparisons.append(
            {
                "model": name,
                "cv_auc": cv_auc(cohort, numeric, categorical),
                "ess_ns": diag["ess_control"],
                "ess_lr": diag["ess_intervention"],
                "max_abs_smd": diag["max_abs_smd"],
                "common_ps_q01_lower": lower,
                "common_ps_q99_upper": upper,
                "outside_common_q01_q99_n": int(((ps < lower) | (ps > upper)).sum()),
                "ps_below_0_05_n": int((ps < .05).sum()),
                "ps_above_0_95_n": int((ps > .95).sum()),
                "max_weight": float(weights.max()),
            }
        )
        balance.insert(0, "model", name)
        balances.append(balance)

    missing = pd.DataFrame(
        [
            {
                "feature": feature,
                "missing_fraction": float(cohort[feature].isna().mean()),
                "available_n": int(cohort[feature].notna().sum()),
                "unweighted_mean_ns": float(cohort.loc[cohort.exposure == 0, feature].mean()),
                "unweighted_mean_lr": float(cohort.loc[cohort.exposure == 1, feature].mean()),
            }
            for feature in augmented
        ]
    )
    definitions = pd.DataFrame(
        [
            ["urine6_ml/urine24_ml", "Sum of urinary output item values whose max(charttime, storetime) was at or before t0", "raw outputevents; 12 urine item IDs"],
            ["urine*_records/hours", "Number of records and distinct availability hours", "measurement coverage, not assumed zero"],
            ["input*_ml", "Completed mL inputs with endtime and storetime at or before t0", "all inputevents measured in mL"],
            ["known_ns24_ml/known_lr24_ml", "Completed NS/LR input in prior 24 h and stored by t0", "availability-time sensitivity to original ns_ml/lr_ml"],
            ["output*_ml", "All positive mL outputs available by t0", "all outputevents measured in mL"],
            ["net*_ml", "Known completed input minus known output", "6 h and 24 h"],
            ["prior_ns_amount/hours_since_prior_ns", "Most recent >=250 mL completed NS bolus known before trigger", "strict store/end availability"],
            ["prior_urine_response_delta", "Post-prior-NS urine rate minus pre-prior-NS urine rate; each observable before t0", "missing when either window has insufficient duration/records"],
            ["creat_delta24", "Last minus first creatinine result available in prior 24 h", "labevents max(charttime, storetime)"],
            ["chloride_delta48", "Last minus first chloride result available in prior 48 h", "labevents max(charttime, storetime)"],
            ["lactate_delta6", "Last minus first lactate result available in prior 6 h", "labevents max(charttime, storetime)"],
        ],
        columns=["feature", "definition", "source"],
    )

    comparison = pd.DataFrame(comparisons)
    balance_all = pd.concat(balances, ignore_index=True)
    diagnostics = {
        "candidate": "N09",
        "database": "MIMIC-IV",
        "n": 706,
        "lr_n": 259,
        "ns_n": 447,
        "tie_rule": "Earliest qualifying starttime; if simultaneous LR and NS, NS takes priority; then endtime and orderid.",
        "outcomes_queried": False,
        "post_t0_features_used": False,
        "models": comparisons,
        "interpretation": "Augmentation covers more prespecified clinical indications but does not prove exchangeability. CV AUC is a treatment-selection description, not a gate for whether a clinical confounder belongs in the model. Missingness indicators are part of the fitted pipeline; missing urine output is not treated as no urine output.",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "patient_level_output": False,
    }
    comparison.to_csv(OUT / "model_comparison.csv", index=False)
    balance_all.to_csv(OUT / "balance.csv", index=False)
    missing.to_csv(OUT / "missingness.csv", index=False)
    definitions.to_csv(OUT / "feature_definitions.csv", index=False)
    (OUT / "diagnostics.json").write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    base = comparison.iloc[0]
    pars = comparison.iloc[1]
    aug = comparison.iloc[2]
    worst_new = (
        balance_all[(balance_all.model == "pre_t0_augmented") & balance_all.covariate.str.contains("|".join(augmented), regex=True)]
        .sort_values("unweighted_abs_smd", ascending=False)
        .head(8)
    )
    new_lines = "\n".join(
        f"- {r.covariate}: raw |SMD| {r.unweighted_abs_smd:.3f}; augmented ATO |SMD| {r.overlap_abs_smd:.3f}"
        for r in worst_new.itertuples()
    )
    report = f"""# N09 pre-t0 confounder audit

## Result

The fixed outcome-blind cohortnant was reproduced: **N=706, LR=259, saline=447**. The action rule now explicitly assigns saline priority when LR and saline share the earliest timestamp. No post-t0 variable and no clinical outcome was queried.

|Model|CV AUC|ESS LR|ESS saline|Maximum weighted |SMD||Outside common 1st-99th percentile PS range|PS <0.05 / >0.95|
|---|---:|---:|---:|---:|---:|---:|
|Legal original, action amount removed|{base.cv_auc:.3f}|{base.ess_lr:.1f}|{base.ess_ns:.1f}|{base.max_abs_smd:.3f}|{int(base.outside_common_q01_q99_n)}|{int(base.ps_below_0_05_n)} / {int(base.ps_above_0_95_n)}|
|Pre-t0 parsimonious|{pars.cv_auc:.3f}|{pars.ess_lr:.1f}|{pars.ess_ns:.1f}|{pars.max_abs_smd:.3f}|{int(pars.outside_common_q01_q99_n)}|{int(pars.ps_below_0_05_n)} / {int(pars.ps_above_0_95_n)}|
|Pre-t0 augmented|{aug.cv_auc:.3f}|{aug.ess_lr:.1f}|{aug.ess_ns:.1f}|{aug.max_abs_smd:.3f}|{int(aug.outside_common_q01_q99_n)}|{int(aug.ps_below_0_05_n)} / {int(aug.ps_above_0_95_n)}|

Changes in cross-validated treatment-model AUC show whether the extra variables add stable out-of-sample prediction of fluid choice. A lower AUC does not prove that clinically required confounders may be omitted, but it warns that a large redundant feature set can add noise. ESS and sparse-support changes show the information cost of indication modeling.

The full augmented model had CV AUC {aug.cv_auc:.3f}, compared with {base.cv_auc:.3f} for the legal original model. Thus augmentation **did not improve cross-validated discrimination**. It nevertheless covers prespecified clinical indications that the original model omitted; AUC is not a causal-model quality gate. The accompanying loss of ESS and increase in sparse support must remain visible in interpretation.

## Strongest newly measured arm differences

{new_lines}

## Decision

The augmented model is the more defensible pre-outcome specification if it does not cause severe ESS or positivity collapse. It captures urine amount and recording coverage, known net balance, timing and response after the preceding saline bolus, and availability-time creatinine/chloride/lactate trajectories. Missingness is modeled rather than interpreted as normal physiology.

This audit still cannot certify no unmeasured confounding. The earlier audit found no credible negative-control outcome, and that limitation remains. Clinical outcomes should be run only once after reviewing the aggregate missingness and support diagnostics in this folder.
"""
    (OUT / "report.md").write_text(report, encoding="utf-8")
    actual = {path.name for path in OUT.iterdir() if path.is_file()}
    if actual != EXPECTED:
        raise RuntimeError(f"Expected six aggregate files, found {sorted(actual)}")
    print(json.dumps(diagnostics, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
