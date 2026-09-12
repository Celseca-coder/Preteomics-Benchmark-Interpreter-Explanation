#!/usr/bin/env python
"""Offline top-k precision / recall / MRR from tabular feature_ranks.csv."""
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path("/autofs/nas8/tywang/tjzou/TME_modeling_benchmark")
OUT = Path("/autofs/nas8/tywang/tjzou/ExplainingConclusions/tables")
PANEL = ROOT / "results" / "pseudo_label_explanations_panel"
REGEX = Path("/autofs/nas8/tywang/tjzou/ExplainingConclusions/tables/tabular_recovery_regex.csv")

FAMILY_DIRS = {
    "composition": PANEL / "tabular_composition",
    "density": PANEL / "tabular_density",
    "mixing": PANEL / "tabular_mixing",
    "point_pattern": PANEL / "tabular_point_pattern",
    "composition_density": PANEL / "tabular_composition_density",
    "composition_mixing": PANEL / "tabular_composition_mixing",
    "density_mixing": PANEL / "tabular_density_mixing",
    "composition_point_pattern": PANEL / "tabular_composition_point_pattern",
    "density_point_pattern": PANEL / "tabular_density_point_pattern",
    "mixing_point_pattern": PANEL / "tabular_mixing_point_pattern",
}

CONTROLS = {"motif_tumor_high", "motif_cd8_high"}


def _norm(name: str) -> str:
    return re.sub(r"\s+", " ", str(name).strip().lower())


def load_rules() -> dict[str, dict]:
    table = pd.read_csv(REGEX)
    rules: dict[str, dict] = {}
    for task, sub in table.groupby("task"):
        hits = tuple(
            str(p) for p in sub.loc[sub["role"].eq("hit"), "pattern"] if isinstance(p, str) and p
        )
        misses = tuple(
            str(p)
            for p in sub.loc[sub["role"].eq("miss_as_top"), "pattern"]
            if isinstance(p, str) and p
        )
        kind = str(sub["kind"].iloc[0])
        rules[str(task)] = dict(hit=hits, miss_as_top=misses, kind=kind)
    return rules


def hits(name: str, patterns: tuple[str, ...]) -> bool:
    if not patterns:
        return False
    text = _norm(name)
    return any(re.search(pat, text) for pat in patterns)


def vocab(family_dir: Path, dataset: str) -> list[str]:
    path = family_dir / dataset / "tabular_features.csv"
    if not path.exists():
        return []
    cols = pd.read_csv(path, nrows=0).columns.tolist()
    return [c for c in cols if c != "region_id"]


def score_ranking(ranked: list[str], rule: dict) -> dict:
    top1 = ranked[0] if ranked else ""
    hit_ranks = [i + 1 for i, name in enumerate(ranked) if hits(name, rule["hit"])]
    n1 = int(bool(ranked) and hits(ranked[0], rule["hit"]))
    n3 = sum(1 for name in ranked[:3] if hits(name, rule["hit"]))
    n5 = sum(1 for name in ranked[:5] if hits(name, rule["hit"]))
    best = hit_ranks[0] if hit_ranks else None
    return dict(
        top1=top1,
        hits_at_1=n1,
        hits_at_3=n3,
        hits_at_5=n5,
        hit_at_1=n1,
        precision_at_3=n3 / 3.0,
        precision_at_5=n5 / 5.0,
        mrr=(1.0 / best) if best else 0.0,
        best_hit_rank=best,
        miss_as_top1=bool(ranked) and hits(ranked[0], rule["miss_as_top"]),
    )


def main() -> None:
    rules = load_rules()
    fold_rows = []
    cand_rows = []
    cand_cache: dict[tuple[str, str, str], int] = {}

    for family, dest in FAMILY_DIRS.items():
        ranks_path = dest / "feature_ranks.csv"
        if not ranks_path.exists():
            print(f"skip missing {ranks_path}")
            continue
        ranks = pd.read_csv(ranks_path)
        if "selected" in ranks.columns:
            ranks = ranks[ranks["selected"].astype(str).isin(["True", "true", "1"])]
        for (dataset, task, method, seed, fold), sub in ranks.groupby(
            ["dataset", "task", "method", "seed", "fold"]
        ):
            rule = rules.get(str(task))
            if rule is None:
                continue
            ordered = (
                sub.sort_values("rank")["feature"].astype(str).tolist()
            )
            rec = score_ranking(ordered, rule)
            key = (family, str(dataset), str(task))
            if key not in cand_cache:
                names = vocab(dest, str(dataset))
                n_h = sum(1 for name in names if hits(name, rule["hit"]))
                cand_cache[key] = n_h
                cand_rows.append(
                    dict(
                        family=family,
                        dataset=dataset,
                        task=task,
                        kind=rule["kind"],
                        n_features=len(names),
                        n_hit_candidates=n_h,
                    )
                )
            n_h = cand_cache[key]
            rec.update(
                family=family,
                dataset=dataset,
                task=task,
                kind=rule["kind"],
                method=method,
                seed=seed,
                fold=fold,
                n_hit_candidates=n_h,
                recall_at_1=rec["hits_at_1"] / min(1, n_h) if n_h else 0.0,
                recall_at_3=rec["hits_at_3"] / min(3, n_h) if n_h else 0.0,
                recall_at_5=rec["hits_at_5"] / min(5, n_h) if n_h else 0.0,
            )
            fold_rows.append(rec)

    folds = pd.DataFrame(fold_rows)
    cands = pd.DataFrame(cand_rows)
    cands.to_csv(OUT / "tabular_topk_hit_candidates.csv", index=False)
    folds.to_csv(OUT / "tabular_topk_precision_folds.csv", index=False)

    keys = [
        "hit_at_1",
        "precision_at_3",
        "precision_at_5",
        "mrr",
        "recall_at_1",
        "recall_at_3",
        "recall_at_5",
        "miss_as_top1",
    ]
    task = (
        folds.groupby(["family", "dataset", "task", "kind", "method", "n_hit_candidates"], sort=False)[keys]
        .mean()
        .reset_index()
    )
    task["n_folds"] = (
        folds.groupby(["family", "dataset", "task", "kind", "method"], sort=False).size().to_numpy()
    )
    task.to_csv(OUT / "tabular_topk_precision_by_task.csv", index=False)

    method = (
        task.groupby(["family", "method", "kind"], sort=False)[
            ["hit_at_1", "precision_at_3", "precision_at_5", "mrr", "recall_at_1", "recall_at_3", "recall_at_5"]
        ]
        .mean()
        .reset_index()
    )
    n_task = task.groupby(["family", "method", "kind"], sort=False).size().rename("n_tasks")
    method = method.merge(n_task, on=["family", "method", "kind"])
    overall = (
        task.groupby(["family", "method"], sort=False)[
            ["hit_at_1", "precision_at_3", "precision_at_5", "mrr", "recall_at_1", "recall_at_3", "recall_at_5"]
        ]
        .mean()
        .reset_index()
    )
    overall["kind"] = "all"
    overall["n_tasks"] = task.groupby(["family", "method"], sort=False).size().to_numpy()
    summary = pd.concat([method, overall], ignore_index=True)
    summary.to_csv(OUT / "tabular_topk_precision_method_summary.csv", index=False)

    print("families", folds["family"].nunique(), "fold-rows", len(folds), "task-rows", len(task))
    print(summary[summary.kind.eq("all")].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
