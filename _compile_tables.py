#!/usr/bin/env python
"""Compile Global interpreter tables into ExplainingConclusions/tables."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/autofs/nas8/tywang/tjzou/TME_modeling_benchmark")
OUT = Path("/autofs/nas8/tywang/tjzou/ExplainingConclusions/tables")
PANEL_DIR = ROOT / "results" / "pseudo_label_explanations_panel"
HNC_DIR = ROOT / "results" / "pseudo_label_explanations"
sys.path.insert(0, str(ROOT))

from benchmark.motifs.panel import load_selected_panel  # noqa: E402
from benchmark.motifs.recovery import RULES, UTAG_RULES, rank_recovery  # noqa: E402

OUT.mkdir(parents=True, exist_ok=True)
CONTROLS = ("motif_tumor_high", "motif_cd8_high")
KIND = {
    "motif_tumor_high": "control",
    "motif_cd8_high": "control",
}


def kind_of(task: str) -> str:
    return "control" if task in CONTROLS else "spatial"


def load_panel() -> pd.DataFrame:
    panel = load_selected_panel(ROOT / "results" / "pseudo_labels_all" / "selected_catalog.csv")
    panel = panel.copy()
    panel["task"] = "motif_" + panel["motif_id"].astype(str)
    panel["kind"] = panel["task"].map(kind_of)
    return panel


def selected_pairs(panel: pd.DataFrame) -> set[tuple[str, str]]:
    return {(str(r.source_dataset), str(r.task)) for r in panel.itertuples(index=False)}


def dataset_controls_ok(rows: pd.DataFrame, pass_col: str) -> dict[str, bool]:
    ok = {}
    for dataset, sub in rows.groupby("dataset"):
        ctrl = sub[sub["task"].isin(CONTROLS)]
        ok[str(dataset)] = bool(not ctrl.empty and ctrl[pass_col].all())
    return ok


def score_method(method: str, rows: pd.DataFrame, pass_col: str, panel: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    selected = selected_pairs(panel)
    controls_ok = dataset_controls_ok(rows, pass_col)
    out = []
    n_ctrl = n_ctrl_pass = n_spatial = n_spatial_pass = n_blocked = 0
    n_spatial_raw = 0
    for dataset, task in sorted(selected):
        sub = rows[(rows["dataset"].astype(str) == dataset) & (rows["task"] == task)]
        kind = kind_of(task)
        ds_ok = controls_ok.get(dataset, False)
        if sub.empty:
            raw = False
            verdict = "missing"
        else:
            raw = bool(sub[pass_col].iloc[0])
            if kind == "control":
                verdict = "pass_control" if raw else "fail_control"
            elif not ds_ok:
                verdict = "fail_control_block"
            else:
                verdict = "pass_spatial" if raw else "fail_spatial"
        counted = raw if kind == "control" else (raw and ds_ok)
        extra = {}
        if not sub.empty:
            for col in sub.columns:
                if col not in {"dataset", "task", pass_col}:
                    extra[col] = sub.iloc[0][col]
        out.append(dict(
            method=method, dataset=dataset, task=task, kind=kind,
            raw_passed=raw, protocol_passed=counted, verdict=verdict,
            dataset_controls_ok=ds_ok, **{k: extra[k] for k in extra if k not in {
                "method", "dataset", "task", "kind", "raw_passed", "protocol_passed",
                "verdict", "dataset_controls_ok", "passed",
            }},
        ))
        if kind == "control":
            n_ctrl += 1
            n_ctrl_pass += int(raw)
        else:
            n_spatial += 1
            n_spatial_raw += int(raw)
            if not ds_ok:
                n_blocked += 1
            else:
                n_spatial_pass += int(raw)
    scored = n_spatial - n_blocked
    n_overall = n_ctrl + n_spatial
    n_overall_raw = n_ctrl_pass + n_spatial_raw
    n_overall_protocol = n_ctrl_pass + n_spatial_pass
    denom = f"{n_overall}" if n_overall else "0"
    summary = dict(
        method=method,
        n_control=n_ctrl,
        n_control_passed=n_ctrl_pass,
        n_spatial=n_spatial,
        n_spatial_blocked=n_blocked,
        n_spatial_scored=scored,
        n_spatial_passed_protocol=n_spatial_pass,
        n_spatial_passed_raw=n_spatial_raw,
        control_ok=n_ctrl_pass == n_ctrl and n_ctrl > 0,
        protocol_spatial=f"{n_spatial_pass}/{scored}" if scored else f"0/{n_spatial}",
        raw_spatial=f"{n_spatial_raw}/{n_spatial}",
        protocol_control=f"{n_ctrl_pass}/{n_ctrl}",
        # Fixed denominator: blocked spatial tasks count as failures, not dropped.
        raw_overall=f"{n_overall_raw}/{denom}",
        protocol_overall=f"{n_overall_protocol}/{denom}",
    )
    return pd.DataFrame(out), summary


def assert_overall_monotone(frame: pd.DataFrame, name: str) -> None:
    """raw_overall ≥ protocol_overall: protocol never credits a raw failure."""
    if frame is None or frame.empty:
        return
    raw_n = frame["n_control_passed"] + frame["n_spatial_passed_raw"]
    prot_n = frame["n_control_passed"] + frame["n_spatial_passed_protocol"]
    bad = frame.loc[raw_n.to_numpy() < prot_n.to_numpy()]
    if not bad.empty:
        cols = ["method", "raw_overall", "protocol_overall"]
        cols = [c for c in cols if c in bad.columns]
        raise RuntimeError(f"{name}: raw_overall < protocol_overall\n{bad[cols].to_string(index=False)}")


def top1_freq(ranks: pd.DataFrame, method: str) -> pd.DataFrame:
    sub = ranks[(ranks["method"] == method) & (ranks["rank"] == 1)].copy()
    if "selected" in sub.columns:
        sub = sub[sub["selected"].astype(str).isin(["True", "true", "1"]) | (sub["selected"] == True)]
    grp = (
        sub.groupby(["dataset", "task", "feature"], dropna=False)
        .size()
        .reset_index(name="n_folds_top1")
    )
    tot = sub.groupby(["dataset", "task"]).size().reset_index(name="n_folds")
    out = grp.merge(tot, on=["dataset", "task"])
    out["frac_top1"] = out["n_folds_top1"] / out["n_folds"]
    out = out.sort_values(["dataset", "task", "frac_top1"], ascending=[True, True, False])
    out["method"] = method
    return out


def rescore_utag_linear() -> pd.DataFrame:
    ranks = pd.read_csv(HNC_DIR / "utag_message_passing" / "feature_ranks.csv")
    folds = pd.read_csv(HNC_DIR / "utag_message_passing" / "fold_recovery.csv")
    rows = []
    for (task, seed, fold, method), sub in ranks.groupby(["task", "seed", "fold", "method"]):
        rule = UTAG_RULES.get(task)
        if rule is None:
            continue
        ranked = sub.sort_values("rank")["feature"].tolist()
        rec = rank_recovery(ranked, rule, k=5)
        rows.append(dict(
            dataset="hnc_wu2022", task=task, seed=seed, fold=fold, method=method,
            kind=rule.kind, **rec,
        ))
    rec_df = pd.DataFrame(rows)
    auc = folds[["task", "seed", "fold", "auc_lasso", "auc_ridge",
                 "faithfulness_drop_top", "faithfulness_drop_random",
                 "faithfulness_passed"]].drop_duplicates()
    rec_df = rec_df.merge(auc, on=["task", "seed", "fold"], how="left")
    return rec_df


FAMILY_DIRS = {
    "composition": PANEL_DIR / "tabular_composition",
    "density": PANEL_DIR / "tabular_density",
    "mixing": PANEL_DIR / "tabular_mixing",
    "point_pattern": PANEL_DIR / "tabular_point_pattern",
}

PAIR_DIRS = {
    "composition_density": PANEL_DIR / "tabular_composition_density",
    "composition_mixing": PANEL_DIR / "tabular_composition_mixing",
    "density_mixing": PANEL_DIR / "tabular_density_mixing",
    "composition_point_pattern": PANEL_DIR / "tabular_composition_point_pattern",
    "density_point_pattern": PANEL_DIR / "tabular_density_point_pattern",
    "mixing_point_pattern": PANEL_DIR / "tabular_mixing_point_pattern",
}


def compile_named_feature_runs(
    panel: pd.DataFrame,
    run_dirs: dict[str, Path],
    out_stem: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score Lasso/SHAP runs that each used a fixed named-feature subset."""
    fold_frames = []
    verdict_frames = []
    summaries = []
    compact = []
    auc_rows = []
    selected = selected_pairs(panel)

    joint = pd.read_csv(PANEL_DIR / "tabular" / "summary.csv")
    joint_idx = joint.set_index(["dataset", "task"])

    for family, path in run_dirs.items():
        if not (path / "summary.csv").is_file():
            print(f"WARN: missing {path / 'summary.csv'}; skip {family}")
            continue
        summary = pd.read_csv(path / "summary.csv")
        summary = summary.copy()
        summary["family"] = family
        summary["lasso_ok"] = summary["frac_lasso_passed"] >= 0.5
        summary["shap_ok"] = summary["frac_shap_passed"] >= 0.5
        fold_frames.append(summary)

        score_input = summary.drop(columns=["family"])
        lasso_tasks, lasso_sum = score_method(f"lasso_{family}", score_input, "lasso_ok", panel)
        shap_tasks, shap_sum = score_method(f"shap_{family}", score_input, "shap_ok", panel)
        lasso_tasks["family"] = family
        shap_tasks["family"] = family
        lasso_tasks["explainer"] = "lasso"
        shap_tasks["explainer"] = "shap"
        verdict_frames.extend([lasso_tasks, shap_tasks])
        lasso_sum["family"] = family
        shap_sum["family"] = family
        lasso_sum["explainer"] = "lasso"
        shap_sum["explainer"] = "shap"
        summaries.extend([lasso_sum, shap_sum])

        ranks = pd.read_csv(path / "feature_ranks.csv")
        for explainer, top in (("lasso", top1_freq(ranks, "lasso")), ("shap", top1_freq(ranks, "shap"))):
            for (dataset, task), sub in top.groupby(["dataset", "task"]):
                row = sub.iloc[0]
                compact.append(dict(
                    family=family, method=explainer, dataset=dataset, task=task,
                    top1_feature=row.feature, frac_top1=row.frac_top1, n_folds=row.n_folds,
                ))

        for dataset, task in sorted(selected):
            sub = summary[(summary["dataset"].astype(str) == dataset) & (summary["task"] == task)]
            jsub = joint_idx.loc[(dataset, task)] if (dataset, task) in joint_idx.index else None
            if sub.empty:
                continue
            row = sub.iloc[0]
            auc_rows.append(dict(
                family=family, dataset=dataset, task=task, kind=kind_of(task),
                mean_auc_lasso=row.mean_auc_lasso,
                mean_auc_ridge=row.mean_auc_ridge,
                frac_lasso_passed=row.frac_lasso_passed,
                frac_shap_passed=row.frac_shap_passed,
                frac_faithfulness_passed=row.frac_faithfulness_passed,
                joint_mean_auc_lasso=None if jsub is None else float(jsub.mean_auc_lasso),
                joint_mean_auc_ridge=None if jsub is None else float(jsub.mean_auc_ridge),
                joint_frac_lasso_passed=None if jsub is None else float(jsub.frac_lasso_passed),
                joint_frac_shap_passed=None if jsub is None else float(jsub.frac_shap_passed),
            ))

    if not fold_frames:
        return pd.DataFrame(), pd.DataFrame()

    folds = pd.concat(fold_frames, ignore_index=True)
    verdicts = pd.concat(verdict_frames, ignore_index=True)
    comparison = pd.DataFrame(summaries)
    folds.to_csv(OUT / f"{out_stem}_fold_summary.csv", index=False)
    verdicts.to_csv(OUT / f"{out_stem}_selected_task_verdicts.csv", index=False)
    comparison.to_csv(OUT / f"{out_stem}_method_comparison.csv", index=False)
    pd.DataFrame(compact).to_csv(OUT / f"{out_stem}_dominant_top1.csv", index=False)
    pd.DataFrame(auc_rows).to_csv(OUT / f"{out_stem}_auc_compare.csv", index=False)

    wide_rows = []
    for dataset, task in sorted(selected):
        row = dict(dataset=dataset, task=task, kind=kind_of(task))
        for family in run_dirs:
            for explainer in ("lasso", "shap"):
                sub = verdicts[
                    (verdicts.family == family)
                    & (verdicts.explainer == explainer)
                    & (verdicts.dataset == dataset)
                    & (verdicts.task == task)
                ]
                prefix = f"{explainer}_{family}"
                if sub.empty:
                    row[f"{prefix}_raw"] = "missing"
                    row[f"{prefix}_verdict"] = "missing"
                else:
                    row[f"{prefix}_raw"] = "pass" if bool(sub.raw_passed.iloc[0]) else "fail"
                    row[f"{prefix}_verdict"] = sub.verdict.iloc[0]
        wide_rows.append(row)
    pd.DataFrame(wide_rows).to_csv(OUT / f"{out_stem}_selected_tasks_wide.csv", index=False)
    return comparison, verdicts


def compile_tabular_families(panel: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    return compile_named_feature_runs(panel, FAMILY_DIRS, "tabular_family")


def compile_tabular_pairs(panel: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    return compile_named_feature_runs(panel, PAIR_DIRS, "tabular_pair")


def summarize_utag_rescore(rec_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (task, method), sub in rec_df.groupby(["task", "method"]):
        kind = sub["kind"].iloc[0]
        rows.append(dict(
            dataset="hnc_wu2022",
            method=f"utag_linear_{method}",
            task=task,
            kind=kind,
            n_folds=len(sub),
            frac_recovery_passed=float(sub["passed"].mean()),
            mean_auc_lasso=float(sub["auc_lasso"].mean()) if "auc_lasso" in sub else float("nan"),
            mean_auc_ridge=float(sub["auc_ridge"].mean()) if "auc_ridge" in sub else float("nan"),
            frac_faithfulness=float(sub["faithfulness_passed"].mean()) if "faithfulness_passed" in sub else float("nan"),
            most_common_top1=sub["top1"].value_counts().index[0] if len(sub) else "",
            recovery_ok=float(sub["passed"].mean()) >= 0.5,
        ))
    return pd.DataFrame(rows)


def export_recovery_regex() -> None:
    rows = []
    wide = []
    for task, rule in RULES.items():
        for i, pat in enumerate(rule.hit, 1):
            rows.append(dict(
                task=task, kind=rule.kind, role="hit", pattern_index=i, pattern=pat,
                pass_note="top-5 feature name must match (case-insensitive after normalize)",
            ))
        if not rule.miss_as_top:
            rows.append(dict(
                task=task, kind=rule.kind, role="miss_as_top", pattern_index=0, pattern="",
                pass_note="no abundance-leak Top-1 rule",
            ))
        for i, pat in enumerate(rule.miss_as_top, 1):
            rows.append(dict(
                task=task, kind=rule.kind, role="miss_as_top", pattern_index=i, pattern=pat,
                pass_note="if Top-1 matches, spatial fold fails even if hit later in top-5",
            ))
        wide.append(dict(
            task=task,
            kind=rule.kind,
            hit_patterns=" | ".join(rule.hit),
            miss_as_top_patterns=" | ".join(rule.miss_as_top) if rule.miss_as_top else "(none)",
            fold_pass_rule=(
                "top-5 hits any hit-pattern"
                if rule.kind == "control"
                else "top-5 hits any hit-pattern AND Top-1 does not match miss_as_top"
            ),
            task_pass_rule=">=50% of 15 folds pass",
        ))
    pd.DataFrame(rows).to_csv(OUT / "tabular_recovery_regex.csv", index=False)
    pd.DataFrame(wide).to_csv(OUT / "tabular_recovery_regex_by_task.csv", index=False)


MIL_NOISY_DIRS = {
    "composition": PANEL_DIR / "mil_noisy_composition",
    "mixing": PANEL_DIR / "mil_noisy_mixing",
    "celltype_density": PANEL_DIR / "mil_noisy_celltype_density",
    "composition_mixing": PANEL_DIR / "mil_noisy_composition_mixing",
    "composition_celltype_density": PANEL_DIR / "mil_noisy_composition_celltype_density",
    "mixing_celltype_density": PANEL_DIR / "mil_noisy_mixing_celltype_density",
}

MIL_IG_DIRS = {
    "composition_mixing": PANEL_DIR / "mil_ig_composition_mixing",
    "mixing": PANEL_DIR / "mil_ig_mixing",
    "mixing_celltype_density": PANEL_DIR / "mil_ig_mixing_celltype_density",
}


def compile_mil_noisy_combos(panel: pd.DataFrame) -> pd.DataFrame:
    """Score AttnMIL feature-group combos (attention + instance explainers)."""
    return compile_mil_combos(panel, MIL_NOISY_DIRS, "mil_noisy_combo", "attention")


def compile_mil_ig_combos(panel: pd.DataFrame) -> pd.DataFrame:
    """Score already-landed noisy IG feature-attribution runs with the same overall columns."""
    return compile_mil_combos(panel, MIL_IG_DIRS, "mil_ig_combo", "ig")


def compile_mil_combos(
    panel: pd.DataFrame,
    combo_dirs: dict[str, Path],
    out_stem: str,
    wide_explainer: str,
) -> pd.DataFrame:
    fold_frames = []
    verdict_frames = []
    summaries = []
    metric_rows = []
    selected = selected_pairs(panel)

    for combo, path in combo_dirs.items():
        summary_path = path / "summary.csv"
        if not summary_path.is_file():
            print(f"WARN: missing {summary_path}; skip MIL combo {combo}")
            continue
        summary = pd.read_csv(summary_path)
        summary = summary.copy()
        summary["combo"] = combo
        fold_frames.append(summary)
        for explainer, sub in summary.groupby("explainer"):
            sub = sub.copy()
            mil_ok = []
            for r in sub.itertuples(index=False):
                if r.kind == "control":
                    mil_ok.append(str(r.verdict) == "pass_control")
                else:
                    mil_ok.append(
                        float(r.frac_loc_passed) >= 0.5 and float(r.frac_faith_passed) >= 0.5
                    )
            score_in = sub.rename(columns={"verdict": "native_verdict"}).copy()
            score_in["mil_ok"] = mil_ok
            keep = [
                "dataset", "task", "mil_ok", "mean_auc", "frac_auc_passed",
                "frac_loc_passed", "frac_faith_passed", "mean_loc_enrichment",
                "mean_delta_lerf_morf", "native_verdict", "loc_verdict", "faith_verdict",
            ]
            keep = [c for c in keep if c in score_in.columns]
            tasks, summ = score_method(
                f"mil_{combo}_{explainer}", score_in[keep], "mil_ok", panel
            )
            tasks["combo"] = combo
            tasks["explainer"] = explainer
            summ["combo"] = combo
            summ["explainer"] = explainer
            verdict_frames.append(tasks)
            summaries.append(summ)
            for dataset, task in sorted(selected):
                row = sub[(sub["dataset"].astype(str) == dataset) & (sub["task"] == task)]
                if row.empty:
                    continue
                r = row.iloc[0]
                metric_rows.append(dict(
                    combo=combo, explainer=explainer, dataset=dataset, task=task,
                    kind=kind_of(task), mean_auc=r.mean_auc,
                    frac_auc_passed=r.frac_auc_passed,
                    frac_loc_passed=r.frac_loc_passed,
                    frac_faith_passed=r.frac_faith_passed,
                    mean_loc_enrichment=r.mean_loc_enrichment,
                    mean_delta_lerf_morf=r.mean_delta_lerf_morf,
                    native_verdict=r.verdict,
                ))

    if not fold_frames:
        return pd.DataFrame()

    folds = pd.concat(fold_frames, ignore_index=True)
    verdicts = pd.concat(verdict_frames, ignore_index=True)
    comparison = pd.DataFrame(summaries)
    folds.to_csv(OUT / f"{out_stem}_fold_summary.csv", index=False)
    verdicts.to_csv(OUT / f"{out_stem}_selected_task_verdicts.csv", index=False)
    comparison.to_csv(OUT / f"{out_stem}_method_comparison.csv", index=False)
    pd.DataFrame(metric_rows).to_csv(OUT / f"{out_stem}_selected_metrics.csv", index=False)

    att = verdicts[verdicts["explainer"] == wide_explainer]
    wide_rows = []
    for dataset, task in sorted(selected):
        row = dict(dataset=dataset, task=task, kind=kind_of(task))
        for combo in combo_dirs:
            sub = att[
                (att.combo == combo) & (att.dataset == dataset) & (att.task == task)
            ]
            if sub.empty:
                row[f"{combo}_raw"] = "missing"
                row[f"{combo}_verdict"] = "missing"
            else:
                row[f"{combo}_raw"] = "pass" if bool(sub.raw_passed.iloc[0]) else "fail"
                row[f"{combo}_verdict"] = sub.verdict.iloc[0]
        wide_rows.append(row)
    pd.DataFrame(wide_rows).to_csv(OUT / f"{out_stem}_{wide_explainer}_wide.csv", index=False)
    return comparison


def main() -> None:
    export_recovery_regex()
    panel = load_panel()
    panel[[
        "motif_id", "source_dataset", "source_dataset_name", "task_type", "panel",
        "definition", "n_labeled", "n_0", "n_1", "auc_composition", "auc_expression",
        "auc_density", "reason",
    ]].to_csv(OUT / "selected_motif_catalog.csv", index=False)

    tabular = pd.read_csv(PANEL_DIR / "tabular" / "summary.csv")
    tabular["lasso_ok"] = tabular["frac_lasso_passed"] >= 0.5
    tabular["shap_ok"] = tabular["frac_shap_passed"] >= 0.5
    tabular.to_csv(OUT / "lasso_shap_fold_summary_all_tasks.csv", index=False)

    lasso_tasks, lasso_sum = score_method("lasso", tabular, "lasso_ok", panel)
    shap_tasks, shap_sum = score_method("shap", tabular, "shap_ok", panel)
    lasso_tasks.to_csv(OUT / "lasso_selected_task_verdicts.csv", index=False)
    shap_tasks.to_csv(OUT / "shap_selected_task_verdicts.csv", index=False)

    ranks = pd.read_csv(PANEL_DIR / "tabular" / "feature_ranks.csv")
    top1_lasso = top1_freq(ranks, "lasso")
    top1_shap = top1_freq(ranks, "shap")
    pd.concat([top1_lasso, top1_shap], ignore_index=True).to_csv(
        OUT / "lasso_shap_top1_features.csv", index=False
    )
    # compact: top feature per selected task
    compact = []
    for method, df in (("lasso", top1_lasso), ("shap", top1_shap)):
        for (dataset, task), sub in df.groupby(["dataset", "task"]):
            row = sub.iloc[0]
            compact.append(dict(
                method=method, dataset=dataset, task=task,
                top1_feature=row.feature, frac_top1=row.frac_top1,
                n_folds=row.n_folds,
            ))
    pd.DataFrame(compact).to_csv(OUT / "lasso_shap_dominant_top1.csv", index=False)

    # embeddings
    kronos = pd.read_csv(PANEL_DIR / "kronos_embedding" / "embedding_probe_summary.csv")
    eva = pd.read_csv(PANEL_DIR / "eva_embedding" / "embedding_probe_summary.csv")
    for name, table in (("kronos", kronos), ("eva", eva)):
        table = table.rename(columns={"mean": "auc"}) if "auc" not in table.columns else table
        table["kind"] = table["task"].map(kind_of)
        table["embed_ok"] = [
            float(r.auc) >= (0.90 if r.kind == "control" else 0.60)
            for r in table.itertuples(index=False)
        ]
        table["method"] = name
        table.to_csv(OUT / f"{name}_embedding_probe.csv", index=False)
        tasks, summary = score_method(name, table, "embed_ok", panel)
        tasks.to_csv(OUT / f"{name}_selected_task_verdicts.csv", index=False)
        if name == "kronos":
            kronos_sum = summary
            kronos_tasks = tasks
        else:
            eva_sum = summary
            eva_tasks = tasks

    # UTAG portraits
    frames = []
    for path in sorted((PANEL_DIR / "utag_native_portraits").glob("*/summary.csv")):
        frame = pd.read_csv(path)
        if "dataset" not in frame.columns:
            frame["dataset"] = path.parent.name
        frames.append(frame)
    utag = pd.concat(frames, ignore_index=True)
    utag["utag_ok"] = utag["passed"].astype(bool)
    utag.to_csv(OUT / "utag_native_portraits_all_tasks.csv", index=False)
    utag_tasks, utag_sum = score_method("utag_cellular", utag, "utag_ok", panel)
    utag_tasks.to_csv(OUT / "utag_cellular_selected_task_verdicts.csv", index=False)

    # UTAG linear (HNC message-passing), rescored with UTAG_RULES
    rec = rescore_utag_linear()
    rec.to_csv(OUT / "utag_linear_hnc_rescored_folds.csv", index=False)
    utag_lin = summarize_utag_rescore(rec)
    utag_lin.to_csv(OUT / "utag_linear_hnc_rescored_summary.csv", index=False)
    orig = pd.read_csv(HNC_DIR / "utag_message_passing" / "summary.csv")
    orig.insert(0, "dataset", "hnc_wu2022")
    orig.insert(0, "note", "original run used tabular recovery rules on UTAG column names")
    orig.to_csv(OUT / "utag_linear_hnc_original_tabular_rules.csv", index=False)

    # GNN if present
    gnn_sum = None
    gnn_path = PANEL_DIR / "gnn_explainer" / "gnn_summary.csv"
    gnn_frames = list((PANEL_DIR / "gnn_explainer").glob("*/gnn_summary.csv"))
    if gnn_path.is_file() or gnn_frames:
        gnn_parts = []
        if gnn_path.is_file():
            gnn_parts.append(pd.read_csv(gnn_path))
        for p in gnn_frames:
            frame = pd.read_csv(p)
            if "dataset" not in frame.columns:
                frame["dataset"] = p.parent.name
            gnn_parts.append(frame)
        gnn = pd.concat(gnn_parts, ignore_index=True).drop_duplicates(
            subset=["dataset", "task"], keep="last"
        )
        if "frac_passed" in gnn.columns:
            gnn["gnn_ok"] = gnn["frac_passed"] >= 0.5
        else:
            gnn["gnn_ok"] = gnn["passed"].astype(bool)
        gnn.to_csv(OUT / "gnn_explainer_task_summary.csv", index=False)
        gnn_tasks, gnn_sum = score_method("gnn_explainer", gnn, "gnn_ok", panel)
        gnn_tasks.to_csv(OUT / "gnn_explainer_selected_task_verdicts.csv", index=False)

    summaries = [lasso_sum, shap_sum, kronos_sum, eva_sum, utag_sum]
    if gnn_sum:
        summaries.append(gnn_sum)
    global_cmp = pd.DataFrame(summaries)
    global_cmp.to_csv(OUT / "global_method_comparison.csv", index=False)
    assert_overall_monotone(global_cmp, "global_method_comparison")

    # selected-task wide table for markdown
    method_frames = [
        ("lasso", lasso_tasks),
        ("shap", shap_tasks),
        ("utag_cellular", utag_tasks),
        ("kronos", kronos_tasks),
        ("eva", eva_tasks),
    ]
    if gnn_sum:
        method_frames.append(("gnn_explainer", gnn_tasks))

    wide_rows = []
    for dataset, task in sorted(selected_pairs(panel)):
        row = dict(dataset=dataset, task=task, kind=kind_of(task))
        for name, tasks in method_frames:
            sub = tasks[(tasks.dataset == dataset) & (tasks.task == task)]
            if sub.empty:
                row[f"{name}_raw"] = "missing"
                row[f"{name}_verdict"] = "missing"
            else:
                row[f"{name}_raw"] = "pass" if bool(sub.raw_passed.iloc[0]) else "fail"
                row[f"{name}_verdict"] = sub.verdict.iloc[0]
        wide_rows.append(row)
    pd.DataFrame(wide_rows).to_csv(OUT / "selected_tasks_all_methods.csv", index=False)

    # Clinical Global Lasso / SHAP (real endpoints, not pseudo-labels)
    clin_dir = ROOT / "results"
    lasso_clin = pd.read_csv(clin_dir / "stability_lasso_feature_summary.csv")
    lasso_clin = lasso_clin.copy()
    lasso_clin["abs_coef"] = lasso_clin["coefficient_mean"].abs()
    lasso_top = (
        lasso_clin.sort_values(["dataset", "task", "class", "abs_coef"], ascending=[True, True, True, False])
        .groupby(["dataset", "task", "class"], as_index=False)
        .head(5)
    )
    lasso_top.to_csv(OUT / "clinical_stability_lasso_top5.csv", index=False)

    shap_cox = pd.read_csv(clin_dir / "shap_global_composition_expression_cox_feature_summary.csv")
    shap_cox = shap_cox[shap_cox["direction_consistency"] == 1].copy()
    shap_cox_top = (
        shap_cox.sort_values(["dataset", "task", "scheme", "class", "mean_abs_shap"], ascending=[True, True, True, True, False])
        .groupby(["dataset", "task", "scheme", "class"], as_index=False)
        .head(5)
    )
    shap_cox_top.to_csv(OUT / "clinical_shap_comp_expr_cox_top5.csv", index=False)

    tcn = clin_dir / "shap_cytocommunity_selected_tcn_group_features_ordered_biological.csv"
    if tcn.is_file():
        pd.read_csv(tcn).to_csv(OUT / "clinical_shap_cytocommunity_tcn_biological.csv", index=False)

    family_cmp, _ = compile_tabular_families(panel)
    pair_cmp, _ = compile_tabular_pairs(panel)
    mil_cmp = compile_mil_noisy_combos(panel)
    mil_ig_cmp = compile_mil_ig_combos(panel)
    assert_overall_monotone(family_cmp, "tabular_family_method_comparison")
    assert_overall_monotone(pair_cmp, "tabular_pair_method_comparison")
    assert_overall_monotone(mil_cmp, "mil_noisy_combo_method_comparison")
    assert_overall_monotone(mil_ig_cmp, "mil_ig_combo_method_comparison")

    overall_cols = [
        "method", "protocol_control", "raw_spatial", "protocol_spatial",
        "raw_overall", "protocol_overall",
    ]
    print("Wrote", OUT)
    print(global_cmp[overall_cols].to_string(index=False))
    print("\nSingle-family Lasso/SHAP")
    print(family_cmp[overall_cols + ["family", "explainer"]].to_string(index=False))
    print("\nPairwise Lasso/SHAP")
    print(pair_cmp[overall_cols + ["family", "explainer"]].to_string(index=False))
    if mil_cmp is not None and not mil_cmp.empty:
        print("\nMIL noisy combos (attention)")
        print(
            mil_cmp[mil_cmp["explainer"] == "attention"][
                ["combo", "protocol_control", "raw_spatial", "protocol_spatial",
                 "raw_overall", "protocol_overall"]
            ].to_string(index=False)
        )
    if mil_ig_cmp is not None and not mil_ig_cmp.empty:
        print("\nMIL IG combos")
        print(
            mil_ig_cmp[
                ["combo", "explainer", "protocol_control", "raw_spatial",
                 "protocol_spatial", "raw_overall", "protocol_overall"]
            ].to_string(index=False)
        )


if __name__ == "__main__":
    main()
