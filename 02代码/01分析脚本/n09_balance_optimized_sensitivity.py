"""Outcome-blind N09 unpenalized-overlap selection, then one sensitivity run."""
from __future__ import annotations

import importlib.util
import json
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.exceptions import ConvergenceWarning
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


ROOT = Path(r"D:\2026project\04AKI课题基础\05因果推断")
PRE_SCRIPT = ROOT / "02代码" / "01分析脚本" / "n09_pre_t0_confounder_audit.py"
OUT = ROOT / "03结果" / "21_FULL_ICU_CAUSAL_PROGRAM_20260907" / "N09" / "balance_optimized_sensitivity_20260908"
EXPECTED = {"protocol.json", "diagnostics.json", "model_diagnostics.csv", "balance.csv", "outcomes.csv", "report.md"}
NS, LR = 225158, 225828
BOOTSTRAP, SEED, MAX_ITER = 500, 20260908, 5000


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def prep(numeric: list[str], categorical: list[str]):
    return ColumnTransformer(
        [
            ("num", make_pipeline(SimpleImputer(strategy="median", add_indicator=True), StandardScaler()), numeric),
            (
                "cat",
                make_pipeline(
                    SimpleImputer(strategy="most_frequent"),
                    OneHotEncoder(handle_unknown="ignore", drop="first"),
                ),
                categorical,
            ),
        ]
    )


def unpenalized_fit(df: pd.DataFrame, numeric: list[str], categorical: list[str], mult=None):
    if mult is None:
        mult = np.ones(len(df), dtype=int)
    expanded = df.loc[df.index.repeat(mult)].reset_index(drop=True)
    y = expanded.exposure.to_numpy(int)
    transformer = prep(numeric, categorical)
    x = transformer.fit_transform(expanded[numeric + categorical])
    model = LogisticRegression(penalty=None, solver="lbfgs", max_iter=MAX_ITER, tol=1e-9)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        model.fit(x, y)
    convergence_warnings = [str(w.message) for w in caught if issubclass(w.category, ConvergenceWarning)]
    ps = np.clip(model.predict_proba(x)[:, 1], 1e-8, 1 - 1e-8)
    weight = np.where(y == 1, 1 - ps, ps)
    converged = not convergence_warnings and int(model.n_iter_[0]) < MAX_ITER
    return expanded, transformer, x, ps, weight, converged, int(model.n_iter_[0]), convergence_warnings


def diagnose(df: pd.DataFrame, numeric: list[str], categorical: list[str], name: str):
    expanded, transformer, x, ps, weight, converged, iterations, warning_text = unpenalized_fit(df, numeric, categorical)
    dense = x.toarray() if hasattr(x, "toarray") else x
    y = expanded.exposure.to_numpy(int)
    rows: list[dict[str, object]] = []
    for index, feature in enumerate(transformer.get_feature_names_out()):
        value = dense[:, index]
        pooled = np.sqrt((value[y == 0].var() + value[y == 1].var()) / 2)
        raw = 0 if pooled == 0 else abs(value[y == 1].mean() - value[y == 0].mean()) / pooled
        weighted = 0 if pooled == 0 else abs(
            np.average(value[y == 1], weights=weight[y == 1])
            - np.average(value[y == 0], weights=weight[y == 0])
        ) / pooled
        rows.append({"model": name, "covariate": str(feature), "unweighted_abs_smd": raw, "overlap_abs_smd": weighted})
    balance = pd.DataFrame(rows).sort_values("overlap_abs_smd", ascending=False)
    lower = max(float(np.quantile(ps[y == 0], .01)), float(np.quantile(ps[y == 1], .01)))
    upper = min(float(np.quantile(ps[y == 0], .99)), float(np.quantile(ps[y == 1], .99)))
    result = {
        "model": name,
        "converged": converged,
        "iterations": iterations,
        "convergence_warning_count": len(warning_text),
        "ess_ns": float(weight[y == 0].sum() ** 2 / np.square(weight[y == 0]).sum()),
        "ess_lr": float(weight[y == 1].sum() ** 2 / np.square(weight[y == 1]).sum()),
        "max_abs_smd": float(balance.overlap_abs_smd.max()),
        "ps_min": float(ps.min()),
        "ps_max": float(ps.max()),
        "ps_below_0_01_n": int((ps < .01).sum()),
        "ps_above_0_99_n": int((ps > .99).sum()),
        "outside_common_q01_q99_n": int(((ps < lower) | (ps > upper)).sum()),
        "weight_max": float(weight.max()),
        "weight_q99": float(np.quantile(weight, .99)),
        "weighted_mass_ns": float(weight[y == 0].sum()),
        "weighted_mass_lr": float(weight[y == 1].sum()),
        "warnings": warning_text,
    }
    return result, balance


def prepare_cohort(include_outcomes: bool):
    pre = load(PRE_SCRIPT, "n09_pre_balance")
    action = pre.load_module(pre.ACTION_AUDIT, "n09_action_balance")
    canonical = action.load_canonical_module()
    con = canonical.connection()
    con.execute("set threads=1")
    con.execute("create temp table n09_riskset as " + action.RISKSET_SQL)
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
    cohort = action.feature_frame(con, "action_fixed").reset_index(drop=True)
    cohort = cohort.merge(pre.extra_features(con, cohort), on="stay_id", how="left", validate="one_to_one")
    cohort[["known_ns24_ml", "known_lr24_ml"]] = cohort[["known_ns24_ml", "known_lr24_ml"]].fillna(0)
    cohort["net6_ml"] = cohort.input6_ml - cohort.output6_ml
    cohort["net24_ml"] = cohort.input24_ml - cohort.output24_ml
    cohort["hours_since_prior_ns"] = (pd.to_datetime(cohort.t0) - pd.to_datetime(cohort.prior_ns_end)).dt.total_seconds() / 3600
    cohort["pre_urine_rate"] = cohort.pre_urine_ml / cohort.prior_response_pre_hours.where(cohort.prior_response_pre_hours >= 1)
    cohort["post_urine_rate"] = cohort.post_urine_ml / cohort.prior_response_post_hours.where(cohort.prior_response_post_hours >= 1)
    cohort["prior_urine_response_delta"] = cohort.post_urine_rate - cohort.pre_urine_rate
    cohort["creat_delta24"] = cohort.creat_last24 - cohort.creat_first24
    cohort["chloride_delta48"] = cohort.chloride_last48 - cohort.chloride_first48
    cohort["lactate_delta6"] = cohort.lactate_last6 - cohort.lactate_first6
    if include_outcomes:
        # action_fixed contains 764 pre-baseline candidates.  feature_frame()
        # applies the renal baseline exclusions and returns the final 706.
        # Register those exact keys before extracting outcomes; positional LIMIT
        # on action_fixed previously reintroduced 49 excluded stays and omitted
        # 49 eligible stays.
        con.register("final_eligible_keys", cohort[["stay_id"]])
        con.execute(
            "create temp table final_action_fixed as "
            "select a.* from action_fixed a inner join final_eligible_keys k using(stay_id)"
        )
        final_n, final_unique = con.execute(
            "select count(*),count(distinct stay_id) from final_action_fixed"
        ).fetchone()
        if int(final_n) != len(cohort) or int(final_unique) != len(cohort):
            raise RuntimeError(
                f"Final outcome key mismatch: rows={final_n}, unique={final_unique}, cohort={len(cohort)}"
            )
        parts = []
        for offset in range(0, len(cohort), 32):
            con.execute(
                "create or replace temp table outcome_batch as "
                f"select * from final_action_fixed order by subject_id,stay_id limit 32 offset {offset}"
            )
            parts.append(
                con.execute(
                    """select a.stay_id,
                      least(
                        (select min(k.charttime) from b.mimiciv_derived.kdigo_stages k where k.stay_id=a.stay_id and k.charttime>a.t0 and k.charttime<=a.t0+interval '7d' and k.aki_stage>=3),
                        (select min(k.charttime) from b.mimiciv_derived.rrt k where k.stay_id=a.stay_id and k.charttime>a.t0 and k.charttime<=a.t0+interval '7d' and k.dialysis_present=1)
                      ) primary_time,
                      case when ad.deathtime is not null and ad.deathtime<=a.outtime then ad.deathtime else a.outtime end competing_time,
                      p.dod
                    from outcome_batch a join b.mimiciv_hosp.patients p using(subject_id)
                    join b.mimiciv_hosp.admissions ad using(hadm_id)"""
                ).df()
            )
        cohort = cohort.merge(pd.concat(parts, ignore_index=True), on="stay_id", how="left", validate="one_to_one")
        cohort = canonical.event_columns(cohort)
    con.close()
    if len(cohort) != 706 or int(cohort.exposure.sum()) != 259:
        raise RuntimeError(f"Frozen cohort mismatch: N={len(cohort)}, LR={int(cohort.exposure.sum())}")
    return cohort, canonical


def effects(df: pd.DataFrame, numeric: list[str], categorical: list[str], canonical, mult=None):
    expanded, _, _, _, weight, converged, _, warning_text = unpenalized_fit(df, numeric, categorical, mult)
    if not converged:
        raise RuntimeError("Unpenalized bootstrap model did not converge: " + " | ".join(warning_text))
    y = expanded.exposure.to_numpy(int)
    result = {}
    for outcome, horizon in (("aki3_krt_7d", 168.0), ("death_7d", 168.0), ("death_28d", 672.0), ("death_90d", 2160.0)):
        time = expanded[f"{outcome}_time"].to_numpy(float)
        event = expanded[f"{outcome}_event"].to_numpy(int)
        risk = [canonical.weighted_aj(time[y == arm], event[y == arm], weight[y == arm], horizon) for arm in (0, 1)]
        result[outcome] = {
            "control_risk": risk[0], "intervention_risk": risk[1],
            "rd": risk[1] - risk[0], "rr": risk[1] / risk[0] if risk[0] > 0 else np.nan,
        }
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    unexpected = {p.name for p in OUT.iterdir() if p.is_file()} - EXPECTED
    if unexpected:
        raise RuntimeError(f"Unexpected pre-existing files: {sorted(unexpected)}")
    cohort, canonical = prepare_cohort(include_outcomes=False)
    base = [
        "age", "chloride", "prior_chloride", "ns_ml", "lr_ml", "sofa_score", "creat", "low7",
        "aki_stage_creat", "sodium", "bicarbonate", "lactate", "potassium", "hemoglobin", "platelet",
        "map_mean", "map_min", "heart_rate", "pressor", "ventilation",
    ]
    pars = [
        "urine6_ml", "urine24_ml", "urine24_hours", "net6_ml", "net24_ml", "known_ns24_ml", "known_lr24_ml",
        "hours_since_prior_ns", "prior_urine_response_delta", "creat_delta24", "chloride_delta48", "lactate_delta6",
    ]
    full = [
        "urine6_ml", "urine24_ml", "urine6_records", "urine24_records", "urine6_hours", "urine24_hours",
        "input6_ml", "input24_ml", "known_ns24_ml", "known_lr24_ml", "output6_ml", "output24_ml", "net6_ml", "net24_ml",
        "prior_ns_amount", "hours_since_prior_ns", "pre_urine_rate", "post_urine_rate", "prior_urine_response_delta",
        "creat_delta24", "creat_n24", "creat_span_hours", "chloride_delta48", "chloride_n48", "chloride_span_hours",
        "lactate_delta6", "lactate_n6", "lactate_span_hours",
    ]
    categorical = ["gender", "first_careunit", "anchor_year_group"]
    specs = {
        "unpenalized_legal_base": base,
        "unpenalized_pre_t0_parsimonious": base + pars,
        "unpenalized_pre_t0_full": base + full,
    }
    diagnostics, balance_frames = [], []
    for name, numeric in specs.items():
        diag, balance = diagnose(cohort, numeric, categorical, name)
        diagnostics.append(diag)
        balance_frames.append(balance)
    diagnostic_df = pd.DataFrame([{k: v for k, v in row.items() if k != "warnings"} for row in diagnostics])
    balance_df = pd.concat(balance_frames, ignore_index=True)

    # Outcome-blind rule: prefer the most comprehensive augmented model that
    # converges, balances every fitted term, retains ESS >=100 in each arm, and
    # has fewer than 10% of patients at PS <.01 or >.99.
    selected = None
    for candidate in ("unpenalized_pre_t0_full", "unpenalized_pre_t0_parsimonious"):
        row = next(item for item in diagnostics if item["model"] == candidate)
        extreme_fraction = (row["ps_below_0_01_n"] + row["ps_above_0_99_n"]) / len(cohort)
        if row["converged"] and row["max_abs_smd"] <= .10 and min(row["ess_ns"], row["ess_lr"]) >= 100 and extreme_fraction < .10:
            selected = candidate
            break

    protocol = {
        "candidate": "N09",
        "database": "MIMIC-IV",
        "cohort": "Fixed N=706; LR=259; saline=447; explicit saline-priority tie rule",
        "outcome_key_source": "Exact stay_id keys returned by feature_frame after renal baseline exclusions; no positional truncation of the 764 pre-baseline actions.",
        "selection_was_outcome_blind": True,
        "models_checked": list(specs),
        "selection_rule": "Prefer full then parsimonious pre-t0 set if unpenalized logistic converges, all fitted SMD <=0.10, arm ESS >=100 and <10% PS below .01 or above .99.",
        "selected_model": selected,
        "outcome_analysis_label": "post-result method sensitivity",
        "bootstrap": BOOTSTRAP if selected else 0,
    }
    outcomes_rows: list[dict[str, object]] = []
    bootstrap_errors: list[dict[str, object]] = []
    if selected:
        # Only now, after outcome-blind selection, rebuild the same fixed cohort
        # with its four prespecified outcome fields.
        outcome_df, canonical = prepare_cohort(include_outcomes=True)
        chosen_numeric = specs[selected]
        point = effects(outcome_df, chosen_numeric, categorical, canonical)
        rng = np.random.default_rng(SEED)
        collected = []
        for replicate in range(BOOTSTRAP):
            mult = np.bincount(rng.integers(0, len(outcome_df), len(outcome_df)), minlength=len(outcome_df))
            try:
                estimate = effects(outcome_df, chosen_numeric, categorical, canonical, mult)
                collected.append({"replicate": replicate, **{f"{o}_{k}": v for o, values in estimate.items() for k, v in values.items()}})
            except Exception as exc:
                bootstrap_errors.append({"replicate": replicate, "error": repr(exc)})
        boot = pd.DataFrame(collected)
        for outcome, estimate in point.items():
            row = {"model": selected, "outcome": outcome, **estimate}
            for metric in ("control_risk", "intervention_risk", "rd", "rr"):
                if len(boot):
                    row[f"{metric}_ci_low"], row[f"{metric}_ci_high"] = np.quantile(boot[f"{outcome}_{metric}"], [.025, .975])
                else:
                    row[f"{metric}_ci_low"], row[f"{metric}_ci_high"] = np.nan, np.nan
            outcomes_rows.append(row)
    else:
        outcomes_rows.append({"model": "none", "outcome": "not_run", "reason": "No augmented unpenalized design passed outcome-blind diagnostics."})

    diagnostic_json = {
        "models": diagnostics,
        "selected_model": selected,
        "outcomes_read": bool(selected),
        "bootstrap_requested": BOOTSTRAP if selected else 0,
        "bootstrap_success": BOOTSTRAP - len(bootstrap_errors) if selected else 0,
        "bootstrap_errors": bootstrap_errors,
        "patient_level_output": False,
        "created_at": datetime.now().isoformat(timespec="seconds"),
    }
    (OUT / "protocol.json").write_text(json.dumps(protocol, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (OUT / "diagnostics.json").write_text(json.dumps(diagnostic_json, ensure_ascii=False, indent=2, default=str) + "\n", encoding="utf-8")
    diagnostic_df.to_csv(OUT / "model_diagnostics.csv", index=False)
    balance_df.to_csv(OUT / "balance.csv", index=False)
    pd.DataFrame(outcomes_rows).to_csv(OUT / "outcomes.csv", index=False)

    lines = []
    for row in diagnostics:
        lines.append(
            f"|{row['model']}|{row['converged']}|{row['iterations']}|{row['ess_lr']:.1f}|{row['ess_ns']:.1f}|"
            f"{row['max_abs_smd']:.6f}|{row['ps_below_0_01_n']} / {row['ps_above_0_99_n']}|{row['outside_common_q01_q99_n']}|"
        )
    result_lines = []
    for row in outcomes_rows:
        if row.get("outcome") == "not_run":
            result_lines.append("No outcome model was run because no augmented design passed.")
        else:
            result_lines.append(
                f"|{row['outcome']}|{row['control_risk']*100:.2f}%|{row['intervention_risk']*100:.2f}%|"
                f"{row['rd']*100:+.2f} ({row['rd_ci_low']*100:+.2f}, {row['rd_ci_high']*100:+.2f})|"
                f"{row['rr']:.3f} ({row['rr_ci_low']:.3f}, {row['rr_ci_high']:.3f})|"
            )
    report = """# N09 unpenalized overlap-weighting sensitivity

## Outcome-blind model selection

|Model|Converged|Iterations|ESS LR|ESS saline|Maximum weighted abs(SMD)|PS <.01 / >.99|Outside common 1st-99th range|
|---|---|---:|---:|---:|---:|---:|---:|
""" + "\n".join(lines) + f"\n\nSelected before reading outcomes: **{selected or 'none'}**. Standard unpenalized logistic overlap weights should exactly balance fitted score terms when the model converges without separation; `balance.csv` retains every fitted term rather than only the worst values.\n\n## Post-result method sensitivity\n\n|Outcome|Saline risk|LR risk|RD percentage points (95% CI)|RR (95% CI)|\n|---|---:|---:|---:|---:|\n" + "\n".join(result_lines) + f"\n\nBootstrap success: {BOOTSTRAP-len(bootstrap_errors) if selected else 0}/{BOOTSTRAP if selected else 0}. This is explicitly a result-after-method sensitivity analysis, not a replacement preregistered primary model. The earlier failed negative-control search and the retrospective treatment-initiation proxy remain limitations.\n"
    (OUT / "report.md").write_text(report, encoding="utf-8")
    actual = {p.name for p in OUT.iterdir() if p.is_file()}
    if actual != EXPECTED:
        raise RuntimeError(f"Expected six aggregate files, found {sorted(actual)}")
    print(json.dumps({"selected": selected, "diagnostics": diagnostics, "bootstrap_success": diagnostic_json["bootstrap_success"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
