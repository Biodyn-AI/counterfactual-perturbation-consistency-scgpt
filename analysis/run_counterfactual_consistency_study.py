#!/usr/bin/env python3
"""Run Idea 05 counterfactual perturbation consistency study end-to-end.

This script intentionally reuses existing single_cell_mechinterp scripts for:
- causal intervention generation,
- perturbation baseline metrics (AUPR/AUROC).

On top of those baselines, it computes idea-specific consistency signals:
- rank consistency (Spearman) between model effects and perturbation deltas,
- sign agreement,
- covariate-controlled consistency (baseline expression + composition shift),
- shuffled perturbation-label null controls.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]
SUBPROJECT_ROOT = REPO_ROOT / "subproject_24_counterfactual_perturbation_consistency"
SINGLE_CELL_ROOT = REPO_ROOT / "single_cell_mechinterp"

# Reuse symbol normalization logic from the existing project code.
if str(SINGLE_CELL_ROOT) not in sys.path:
    sys.path.insert(0, str(SINGLE_CELL_ROOT))

from src.eval.gene_symbols import (  # noqa: E402
    canonical_symbol,
    load_hgnc_alias_map,
    normalize_edges,
    normalize_gene_names,
)


@dataclass(frozen=True)
class PerturbEvaluation:
    label: str
    config_path: str


@dataclass(frozen=True)
class StudyRun:
    name: str
    causal_config_path: str
    perturbation_evaluations: tuple[PerturbEvaluation, ...]


DEFAULT_RUNS: tuple[StudyRun, ...] = (
    StudyRun(
        name="adamson",
        causal_config_path="configs/causal_intervention_adamson.yaml",
        perturbation_evaluations=(
            PerturbEvaluation(
                label="adamson",
                config_path="configs/perturbation_validation_adamson.yaml",
            ),
        ),
    ),
    StudyRun(
        name="dixit_tiny_neg",
        causal_config_path="configs/causal_intervention_dixit_targets_tiny_neg.yaml",
        perturbation_evaluations=(
            PerturbEvaluation(
                label="dixit_13d",
                config_path="configs/perturbation_validation_dixit_causal_dixit_targets_tiny_neg.yaml",
            ),
            PerturbEvaluation(
                label="dixit_7d",
                config_path="configs/perturbation_validation_dixit7_causal_dixit_targets_tiny_neg.yaml",
            ),
        ),
    ),
    StudyRun(
        name="shifrut",
        causal_config_path="configs/causal_intervention_shifrut_targets_tiny_neg.yaml",
        perturbation_evaluations=(
            PerturbEvaluation(
                label="shifrut",
                config_path="configs/perturbation_validation_shifrut_causal_shifrut_targets_tiny_neg.yaml",
            ),
        ),
    ),
)


def _resolve_single_cell_path(path_value: str | Path) -> Path:
    path = Path(path_value)
    if path.is_absolute():
        return path
    return SINGLE_CELL_ROOT / path


def _load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _mean_expression(adata: sc.AnnData) -> np.ndarray:
    X = adata.X
    mean = X.mean(axis=0)
    if hasattr(mean, "A1"):
        return mean.A1
    return np.asarray(mean).ravel()


def _parse_source(
    label: str,
    delimiter: str | None,
    allow_multi: bool,
    alias_map: Dict[str, str],
    control_labels: Iterable[str],
) -> str:
    raw = str(label).strip()
    parts = [raw]
    if delimiter and delimiter in raw:
        parts = [part.strip() for part in raw.split(delimiter) if part.strip()]

    control_set = {str(item).strip().lower() for item in control_labels}
    filtered = [part for part in parts if part.lower() not in control_set]
    if not filtered:
        return ""
    if len(filtered) > 1 and not allow_multi:
        return ""
    return canonical_symbol(filtered[0], alias_map)


def _choose_composition_key(
    obs_columns: Iterable[str],
    explicit_key: str | None,
) -> str | None:
    columns = set(obs_columns)
    if explicit_key and explicit_key in columns:
        return explicit_key
    for candidate in ("cell_type", "celltype", "cell_line", "sample", "patient"):
        if candidate in columns:
            return candidate
    return None


def _l1_composition_shift(
    control_distribution: pd.Series,
    group_values: pd.Series,
) -> float:
    group_distribution = group_values.value_counts(normalize=True)
    universe = control_distribution.index.union(group_distribution.index)
    control_aligned = control_distribution.reindex(universe, fill_value=0.0)
    group_aligned = group_distribution.reindex(universe, fill_value=0.0)
    return float((control_aligned - group_aligned).abs().sum())


def _spearman(values_a: np.ndarray, values_b: np.ndarray) -> float:
    if values_a.size < 3 or values_b.size < 3:
        return float("nan")
    series_a = pd.Series(values_a)
    series_b = pd.Series(values_b)
    if series_a.nunique(dropna=True) < 2 or series_b.nunique(dropna=True) < 2:
        return float("nan")
    return float(series_a.corr(series_b, method="spearman"))


def _sign_agreement(effect: np.ndarray, delta: np.ndarray) -> float:
    if effect.size == 0:
        return float("nan")
    effect_sign = np.sign(effect)
    delta_sign = np.sign(delta)
    valid = (effect_sign != 0) & (delta_sign != 0)
    if valid.sum() == 0:
        return float("nan")
    return float(np.mean(effect_sign[valid] == delta_sign[valid]))


def _residualize(y: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    if y.ndim != 1:
        raise ValueError("y must be one-dimensional")
    if covariates.ndim != 2:
        raise ValueError("covariates must be two-dimensional")

    X = np.column_stack([np.ones(len(y), dtype=float), covariates.astype(float)])
    mask = (~np.isnan(y)) & (~np.any(np.isnan(X), axis=1))
    resid = np.full(len(y), np.nan, dtype=float)
    if mask.sum() <= X.shape[1]:
        return resid

    beta, *_ = np.linalg.lstsq(X[mask], y[mask], rcond=None)
    resid[mask] = y[mask] - X[mask] @ beta
    return resid


def _compute_observed_metrics(matched_df: pd.DataFrame) -> dict[str, float]:
    effect = matched_df["effect_mean"].to_numpy(dtype=float)
    delta = matched_df["delta"].to_numpy(dtype=float)
    effect_abs = np.abs(effect)
    delta_abs = np.abs(delta)

    covariates = np.column_stack(
        [
            np.log1p(np.maximum(matched_df["baseline_expr"].to_numpy(dtype=float), 0.0)),
            matched_df["composition_shift"].to_numpy(dtype=float),
        ]
    )

    resid_effect_abs = _residualize(effect_abs, covariates)
    resid_delta_abs = _residualize(delta_abs, covariates)
    valid_resid = (~np.isnan(resid_effect_abs)) & (~np.isnan(resid_delta_abs))

    controlled_spearman = float("nan")
    if valid_resid.sum() >= 3:
        controlled_spearman = _spearman(resid_effect_abs[valid_resid], resid_delta_abs[valid_resid])

    observed = {
        "raw_spearman_abs": _spearman(effect_abs, delta_abs),
        "raw_spearman_signed": _spearman(effect, delta),
        "raw_sign_agreement": _sign_agreement(effect, delta),
        "controlled_spearman_abs": controlled_spearman,
        "n_matched_pairs": float(len(matched_df)),
    }
    observed["controlled_minus_raw_abs"] = (
        observed["controlled_spearman_abs"] - observed["raw_spearman_abs"]
        if not np.isnan(observed["controlled_spearman_abs"]) and not np.isnan(observed["raw_spearman_abs"])
        else float("nan")
    )
    return observed


def _bootstrap_metric_intervals(
    matched_df: pd.DataFrame,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> dict[str, float]:
    if matched_df.empty or n_bootstrap <= 0:
        return {
            "raw_spearman_abs_ci_low": float("nan"),
            "raw_spearman_abs_ci_high": float("nan"),
            "raw_sign_agreement_ci_low": float("nan"),
            "raw_sign_agreement_ci_high": float("nan"),
        }

    spearman_values: list[float] = []
    sign_values: list[float] = []
    n = len(matched_df)

    for _ in range(n_bootstrap):
        indices = rng.integers(0, n, size=n)
        sample = matched_df.iloc[indices]
        observed = _compute_observed_metrics(sample)
        spearman_values.append(observed["raw_spearman_abs"])
        sign_values.append(observed["raw_sign_agreement"])

    spearman_arr = np.asarray(spearman_values, dtype=float)
    sign_arr = np.asarray(sign_values, dtype=float)
    spearman_valid = spearman_arr[~np.isnan(spearman_arr)]
    sign_valid = sign_arr[~np.isnan(sign_arr)]

    return {
        "raw_spearman_abs_ci_low": float(np.percentile(spearman_valid, 2.5))
        if spearman_valid.size > 0
        else float("nan"),
        "raw_spearman_abs_ci_high": float(np.percentile(spearman_valid, 97.5))
        if spearman_valid.size > 0
        else float("nan"),
        "raw_sign_agreement_ci_low": float(np.percentile(sign_valid, 2.5))
        if sign_valid.size > 0
        else float("nan"),
        "raw_sign_agreement_ci_high": float(np.percentile(sign_valid, 97.5))
        if sign_valid.size > 0
        else float("nan"),
    }


def _compute_source_shuffle_null(
    matched_df: pd.DataFrame,
    n_shuffle: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    if matched_df.empty:
        return pd.DataFrame()

    unique_sources = np.array(sorted(matched_df["source"].unique()))
    if unique_sources.size < 2 or n_shuffle <= 0:
        return pd.DataFrame()

    source_to_delta: dict[str, dict[str, float]] = {}
    for source, source_df in matched_df.groupby("source"):
        source_to_delta[source] = dict(zip(source_df["target"], source_df["delta"]))

    effect = matched_df["effect_mean"].to_numpy(dtype=float)
    targets = matched_df["target"].to_numpy(dtype=object)
    sources = matched_df["source"].to_numpy(dtype=object)

    null_rows: list[dict[str, float]] = []
    for perm_idx in range(n_shuffle):
        shuffled_sources = rng.permutation(unique_sources)
        mapping = dict(zip(unique_sources, shuffled_sources, strict=True))

        perm_delta = np.full(len(matched_df), np.nan, dtype=float)
        for row_idx, (source, target) in enumerate(zip(sources, targets, strict=True)):
            mapped_source = mapping[source]
            perm_delta[row_idx] = source_to_delta.get(mapped_source, {}).get(target, np.nan)

        valid = ~np.isnan(perm_delta)
        if valid.sum() < 3:
            continue

        perm_effect = effect[valid]
        perm_delta_valid = perm_delta[valid]
        null_rows.append(
            {
                "shuffle_index": float(perm_idx),
                "n_pairs": float(valid.sum()),
                "null_spearman_abs": _spearman(np.abs(perm_effect), np.abs(perm_delta_valid)),
                "null_spearman_signed": _spearman(perm_effect, perm_delta_valid),
                "null_sign_agreement": _sign_agreement(perm_effect, perm_delta_valid),
            }
        )

    return pd.DataFrame(null_rows)


def _empirical_p_value(null_values: np.ndarray, observed_value: float) -> float:
    valid = null_values[~np.isnan(null_values)]
    if valid.size == 0 or np.isnan(observed_value):
        return float("nan")
    count = np.sum(valid >= observed_value)
    return float((count + 1) / (valid.size + 1))


def _null_z_score(null_values: np.ndarray, observed_value: float) -> float:
    valid = null_values[~np.isnan(null_values)]
    if valid.size < 2 or np.isnan(observed_value):
        return float("nan")
    std = float(np.std(valid, ddof=1))
    if std <= 0:
        return float("nan")
    return float((observed_value - float(np.mean(valid))) / std)


def _build_perturbation_profiles(
    perturb_cfg: dict,
    alias_map: Dict[str, str],
    allowed_sources: set[str],
    allowed_targets: set[str],
) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, float]], dict[str, float], dict[str, float]]:
    paths = perturb_cfg["paths"]
    validation_cfg = perturb_cfg.get("perturbation_validation", {})

    adata_path = _resolve_single_cell_path(paths["perturbation_h5ad"])
    adata = sc.read_h5ad(adata_path)

    obs_key = validation_cfg.get("obs_key", "perturbation")
    control_labels = validation_cfg.get("control_labels", ["control", "ctrl"])
    delimiter = validation_cfg.get("delimiter", "+")
    allow_multi = bool(validation_cfg.get("allow_multi", False))
    min_cells = int(validation_cfg.get("min_cells", 10))

    if obs_key not in adata.obs:
        raise ValueError(f"{adata_path}: obs key '{obs_key}' not found")

    obs_values = adata.obs[obs_key].astype(str)
    control_set = {str(item).strip().lower() for item in control_labels}
    control_mask = obs_values.str.lower().isin(control_set).to_numpy()
    if int(control_mask.sum()) == 0:
        raise ValueError(f"{adata_path}: no control cells found for {obs_key}")

    gene_names_norm = normalize_gene_names(adata.var_names.values, alias_map)
    gene_to_idx: dict[str, int] = {}
    for idx, gene in enumerate(gene_names_norm):
        if gene and gene not in gene_to_idx:
            gene_to_idx[gene] = idx

    target_to_idx = {target: gene_to_idx[target] for target in allowed_targets if target in gene_to_idx}

    control_mean = _mean_expression(adata[control_mask])
    baseline_by_target = {target: float(control_mean[idx]) for target, idx in target_to_idx.items()}

    composition_key = _choose_composition_key(
        adata.obs.columns,
        validation_cfg.get("composition_key"),
    )
    control_distribution = None
    if composition_key is not None:
        control_distribution = adata.obs.loc[control_mask, composition_key].astype(str).value_counts(normalize=True)

    delta_by_source: dict[str, dict[str, float]] = {}
    source_covariates: dict[str, dict[str, float]] = {}

    for label in sorted(obs_values.unique()):
        if str(label).strip().lower() in control_set:
            continue
        source = _parse_source(label, delimiter, allow_multi, alias_map, control_labels)
        if not source or source not in allowed_sources:
            continue

        group_mask = obs_values == label
        group_count = int(group_mask.sum())
        if group_count < min_cells:
            continue

        group_mean = _mean_expression(adata[group_mask.to_numpy()])
        delta = group_mean - control_mean

        delta_by_source[source] = {
            target: float(delta[idx]) for target, idx in target_to_idx.items()
        }

        composition_shift = 0.0
        if composition_key is not None and control_distribution is not None:
            composition_shift = _l1_composition_shift(
                control_distribution,
                adata.obs.loc[group_mask.to_numpy(), composition_key].astype(str),
            )

        source_covariates[source] = {
            "n_cells": float(group_count),
            "composition_shift": float(composition_shift),
            "obs_label": str(label),
        }

    adata.file.close() if getattr(adata, "isbacked", False) else None
    return delta_by_source, source_covariates, baseline_by_target, {
        "n_total_sources": float(len(delta_by_source)),
        "n_target_overlap": float(len(target_to_idx)),
    }


def _build_matched_table(
    causal_df: pd.DataFrame,
    delta_by_source: dict[str, dict[str, float]],
    source_covariates: dict[str, dict[str, float]],
    baseline_by_target: dict[str, float],
) -> pd.DataFrame:
    rows: list[dict[str, float | str]] = []
    for row in causal_df.itertuples(index=False):
        source = str(row.source)
        target = str(row.target)

        source_table = delta_by_source.get(source)
        if source_table is None or target not in source_table:
            continue

        baseline_expr = baseline_by_target.get(target)
        if baseline_expr is None:
            continue

        cov = source_covariates.get(source, {})
        rows.append(
            {
                "source": source,
                "target": target,
                "effect_mean": float(row.effect_mean),
                "delta": float(source_table[target]),
                "baseline_expr": float(baseline_expr),
                "composition_shift": float(cov.get("composition_shift", 0.0)),
                "n_cells": float(cov.get("n_cells", np.nan)),
            }
        )

    return pd.DataFrame(rows)


def _compute_timepoint_consistency(
    matched_tables: dict[tuple[str, str], pd.DataFrame],
) -> pd.DataFrame:
    key_13d = ("dixit_13d", "ablation")
    key_7d = ("dixit_7d", "ablation")
    if key_13d not in matched_tables or key_7d not in matched_tables:
        return pd.DataFrame()

    table_13d = matched_tables[key_13d][["source", "target", "delta"]].rename(
        columns={"delta": "delta_13d"}
    )
    table_7d = matched_tables[key_7d][["source", "target", "delta"]].rename(
        columns={"delta": "delta_7d"}
    )
    merged = table_13d.merge(table_7d, on=["source", "target"], how="inner")
    if merged.empty:
        return pd.DataFrame()

    delta_13d = merged["delta_13d"].to_numpy(dtype=float)
    delta_7d = merged["delta_7d"].to_numpy(dtype=float)
    sign_valid = (np.sign(delta_13d) != 0) & (np.sign(delta_7d) != 0)
    sign_agreement = (
        float(np.mean(np.sign(delta_13d[sign_valid]) == np.sign(delta_7d[sign_valid])))
        if sign_valid.sum() > 0
        else float("nan")
    )

    row = {
        "dataset_pair": "dixit_13d_vs_7d",
        "intervention": "ablation",
        "n_overlap_pairs": float(len(merged)),
        "spearman_abs_delta": _spearman(np.abs(delta_13d), np.abs(delta_7d)),
        "spearman_signed_delta": _spearman(delta_13d, delta_7d),
        "sign_agreement_delta": sign_agreement,
    }
    return pd.DataFrame([row])


def _run_command(command: list[str], cwd: Path) -> None:
    print(f"[run] cwd={cwd} :: {' '.join(command)}", flush=True)
    env = dict(os.environ)
    env["PYTHONPATH"] = "."
    subprocess.run(command, cwd=cwd, check=True, env=env)


def _ensure_causal_scores(causal_config_path: Path, skip_causal: bool) -> tuple[Path, dict]:
    causal_cfg = _load_yaml(causal_config_path)
    ci_cfg = causal_cfg.get("causal_intervention", {})
    output_dir = _resolve_single_cell_path(
        ci_cfg.get("output_dir") or causal_cfg.get("paths", {}).get("causal_output_dir", "outputs/causal")
    )
    scores_path = output_dir / "causal_scores.tsv"

    if scores_path.exists() or skip_causal:
        return scores_path, causal_cfg

    _run_command(
        [
            "python",
            "scripts/run_causal_interventions.py",
            "--config",
            str(causal_config_path.relative_to(SINGLE_CELL_ROOT)),
        ],
        cwd=SINGLE_CELL_ROOT,
    )
    return scores_path, causal_cfg


def _run_perturbation_baseline(perturb_config_path: Path, skip_eval: bool) -> tuple[Path, dict]:
    perturb_cfg = _load_yaml(perturb_config_path)
    output_dir = _resolve_single_cell_path(perturb_cfg["paths"].get("output_dir", "outputs/perturb_validation"))
    metrics_path = output_dir / "perturbation_metrics.tsv"

    if not skip_eval:
        _run_command(
            [
                "python",
                "scripts/evaluate_perturbation_validation.py",
                "--config",
                str(perturb_config_path.relative_to(SINGLE_CELL_ROOT)),
                "--permutations",
                "300",
            ],
            cwd=SINGLE_CELL_ROOT,
        )
    return metrics_path, perturb_cfg


def _select_runs(requested_runs: list[str] | None) -> list[StudyRun]:
    if not requested_runs:
        return list(DEFAULT_RUNS)

    requested = {item.strip() for item in requested_runs if item.strip()}
    selected = [run for run in DEFAULT_RUNS if run.name in requested]
    missing = requested - {run.name for run in selected}
    if missing:
        raise ValueError(f"Unknown run names: {sorted(missing)}")
    return selected


def _write_report(
    metrics_df: pd.DataFrame,
    baseline_df: pd.DataFrame,
    timepoint_df: pd.DataFrame,
    report_path: Path,
) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)

    primary = metrics_df[metrics_df["intervention"] == "ablation"].copy()
    if primary.empty:
        primary = metrics_df.copy()

    if not primary.empty:
        primary = primary.sort_values(["dataset_label", "intervention"])

    lines: list[str] = []
    lines.append("# Counterfactual Perturbation Consistency Report")
    lines.append("")
    lines.append("## Question")
    lines.append(
        "Do scGPT causal intervention effects align with real perturbation outcomes at matched gene-pair level?"
    )
    lines.append("")

    lines.append("## Evidence")
    if primary.empty:
        lines.append("No evaluable dataset/intervention rows were produced.")
    else:
        for row in primary.itertuples(index=False):
            lines.append(
                "- "
                f"{row.dataset_label} ({row.intervention}): "
                f"raw |rho|={row.raw_spearman_abs:.3f}, "
                f"95% CI [{row.raw_spearman_abs_ci_low:.3f}, {row.raw_spearman_abs_ci_high:.3f}], "
                f"controlled |rho|={row.controlled_spearman_abs:.3f}, "
                f"sign agreement={row.raw_sign_agreement:.3f}, "
                f"sign 95% CI [{row.raw_sign_agreement_ci_low:.3f}, {row.raw_sign_agreement_ci_high:.3f}], "
                f"null p(|rho|)={row.null_p_spearman_abs:.4f}, "
                f"matched pairs={int(row.n_matched_pairs)}"
            )

    lines.append("")
    lines.append("## Baseline Perturbation Benchmark")
    if baseline_df.empty:
        lines.append("No perturbation benchmark rows were available.")
    else:
        baseline_primary = baseline_df[baseline_df["intervention"] == "ablation"].copy()
        if baseline_primary.empty:
            baseline_primary = baseline_df.copy()
        for row in baseline_primary.sort_values(["dataset_label", "intervention"]).itertuples(index=False):
            lines.append(
                "- "
                f"{row.dataset_label} ({row.intervention}): "
                f"AUPR={row.aupr:.3f}, AUROC={row.auroc:.3f}, "
                f"n_pairs={int(row.n_pairs)}, n_pos={int(row.n_pos)}"
            )

    lines.append("")
    lines.append("## Timepoint Robustness")
    if timepoint_df.empty:
        lines.append("No timepoint-overlap robustness rows were available.")
    else:
        for row in timepoint_df.itertuples(index=False):
            lines.append(
                "- "
                f"{row.dataset_pair} ({row.intervention}): "
                f"n_overlap={int(row.n_overlap_pairs)}, "
                f"Spearman |delta|={row.spearman_abs_delta:.3f}, "
                f"Spearman signed delta={row.spearman_signed_delta:.3f}, "
                f"delta sign agreement={row.sign_agreement_delta:.3f}"
            )

    lines.append("")
    lines.append("## Interpretation")
    if primary.empty:
        lines.append("Insufficient matched pairs to interpret consistency.")
    else:
        mean_raw = float(primary["raw_spearman_abs"].mean(skipna=True))
        mean_ctrl = float(primary["controlled_spearman_abs"].mean(skipna=True))
        improved = int((primary["controlled_minus_raw_abs"] > 0).sum())
        total = int(primary["controlled_minus_raw_abs"].notna().sum())
        lines.append(
            f"Across evaluated rows, mean raw |rho|={mean_raw:.3f} and mean controlled |rho|={mean_ctrl:.3f}."
        )
        lines.append(
            f"Control adjustment improved |rho| in {improved}/{total} evaluable rows."
        )
        lines.append(
            "This supports partial counterfactual consistency while indicating that source/target coverage and confounding still limit alignment."
        )

    lines.append("")
    lines.append("## Limitations")
    lines.append("- Causal scores are from low-compute pair subsets, not full perturbation universes.")
    lines.append("- Composition control is coarse (group-level distribution shift), not cell-level causal adjustment.")
    lines.append("- Shuffled-label nulls are source-permutation controls, not full generative null models.")

    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _plot_consistency(metrics_df: pd.DataFrame, output_path: Path) -> None:
    if metrics_df.empty:
        return

    primary = metrics_df[metrics_df["intervention"] == "ablation"].copy()
    if primary.empty:
        primary = metrics_df.copy()

    plot_df = primary[["dataset_label", "raw_spearman_abs", "null_mean_spearman_abs", "null_std_spearman_abs"]].copy()
    plot_df = plot_df.sort_values("dataset_label").reset_index(drop=True)
    if plot_df.empty:
        return

    x = np.arange(len(plot_df))
    width = 0.36

    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.bar(x - width / 2, plot_df["raw_spearman_abs"], width=width, label="Observed |rho|", color="#3264A8")
    ax.bar(
        x + width / 2,
        plot_df["null_mean_spearman_abs"],
        width=width,
        yerr=plot_df["null_std_spearman_abs"],
        label="Null mean |rho| (source-shuffle)",
        color="#B85C38",
        capsize=4,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(plot_df["dataset_label"], rotation=20, ha="right")
    ax.set_ylabel("Spearman |rho|")
    ax.set_title("Counterfactual Consistency vs Shuffled-Label Null")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def run_study(
    selected_runs: list[StudyRun],
    shuffle_n: int,
    bootstrap_n: int,
    seed: int,
    skip_causal: bool,
    skip_eval: bool,
) -> None:
    outputs_dir = SUBPROJECT_ROOT / "outputs"
    tables_dir = outputs_dir / "tables"
    figures_dir = outputs_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    metrics_rows: list[dict[str, float | str]] = []
    null_rows: list[dict[str, float | str]] = []
    baseline_rows: list[dict[str, float | str]] = []
    matched_tables: dict[tuple[str, str], pd.DataFrame] = {}

    rng = np.random.default_rng(seed)

    for run in selected_runs:
        print(f"[study] run={run.name}", flush=True)
        causal_config_path = SINGLE_CELL_ROOT / run.causal_config_path
        causal_scores_path, causal_cfg = _ensure_causal_scores(causal_config_path, skip_causal=skip_causal)
        if not causal_scores_path.exists():
            raise FileNotFoundError(f"Missing causal scores after run: {causal_scores_path}")

        alias_path = _resolve_single_cell_path(causal_cfg["paths"].get("hgnc_alias_tsv"))
        alias_map = load_hgnc_alias_map(alias_path)

        causal_scores_df = pd.read_csv(causal_scores_path, sep="\t")
        causal_scores_df = normalize_edges(causal_scores_df, alias_map)

        allowed_sources = set(causal_scores_df["source"].dropna().astype(str))
        allowed_targets = set(causal_scores_df["target"].dropna().astype(str))

        for perturb_eval in run.perturbation_evaluations:
            print(f"[study]   dataset={perturb_eval.label}", flush=True)
            perturb_config_path = SINGLE_CELL_ROOT / perturb_eval.config_path
            baseline_metrics_path, perturb_cfg = _run_perturbation_baseline(
                perturb_config_path,
                skip_eval=skip_eval,
            )

            if baseline_metrics_path.exists():
                baseline_df = pd.read_csv(baseline_metrics_path, sep="\t")
                for row in baseline_df.itertuples(index=False):
                    baseline_rows.append(
                        {
                            "run_name": run.name,
                            "dataset_label": perturb_eval.label,
                            "intervention": str(row.intervention),
                            "aupr": float(row.aupr),
                            "auroc": float(row.auroc),
                            "perm_p_value": float(row.perm_p_value),
                            "n_pairs": float(row.n_pairs),
                            "n_pos": float(row.n_pos),
                            "n_sources": float(row.n_sources),
                        }
                    )

            (
                delta_by_source,
                source_covariates,
                baseline_by_target,
                profile_stats,
            ) = _build_perturbation_profiles(
                perturb_cfg,
                alias_map,
                allowed_sources=allowed_sources,
                allowed_targets=allowed_targets,
            )

            for intervention in sorted(causal_scores_df["intervention"].dropna().unique()):
                subset = causal_scores_df[causal_scores_df["intervention"] == intervention].copy()
                matched_df = _build_matched_table(
                    subset,
                    delta_by_source,
                    source_covariates,
                    baseline_by_target,
                )
                matched_df = matched_df.sort_values(["source", "target"], ignore_index=True)
                matched_tables[(perturb_eval.label, str(intervention))] = matched_df.copy()
                matched_export_path = (
                    tables_dir
                    / f"matched_pairs_{run.name}_{perturb_eval.label}_{intervention}.tsv"
                )
                matched_df.to_csv(matched_export_path, sep="\t", index=False)

                observed = _compute_observed_metrics(matched_df)
                bootstrap_intervals = _bootstrap_metric_intervals(
                    matched_df=matched_df,
                    n_bootstrap=bootstrap_n,
                    rng=rng,
                )
                null_df = _compute_source_shuffle_null(matched_df, n_shuffle=shuffle_n, rng=rng)

                null_mean_abs = float(null_df["null_spearman_abs"].mean()) if not null_df.empty else float("nan")
                null_std_abs = (
                    float(null_df["null_spearman_abs"].std(ddof=1)) if len(null_df) > 1 else float("nan")
                )
                null_p_abs = _empirical_p_value(
                    null_df["null_spearman_abs"].to_numpy(dtype=float)
                    if not null_df.empty
                    else np.array([]),
                    observed["raw_spearman_abs"],
                )
                null_p_sign = _empirical_p_value(
                    null_df["null_sign_agreement"].to_numpy(dtype=float)
                    if not null_df.empty
                    else np.array([]),
                    observed["raw_sign_agreement"],
                )
                null_z_abs = _null_z_score(
                    null_df["null_spearman_abs"].to_numpy(dtype=float)
                    if not null_df.empty
                    else np.array([]),
                    observed["raw_spearman_abs"],
                )

                metrics_rows.append(
                    {
                        "run_name": run.name,
                        "dataset_label": perturb_eval.label,
                        "intervention": str(intervention),
                        "n_causal_pairs": float(len(subset)),
                        "n_profile_sources": profile_stats["n_total_sources"],
                        "n_target_overlap": profile_stats["n_target_overlap"],
                        "n_matched_pairs": observed["n_matched_pairs"],
                        "raw_spearman_abs": observed["raw_spearman_abs"],
                        "raw_spearman_signed": observed["raw_spearman_signed"],
                        "raw_sign_agreement": observed["raw_sign_agreement"],
                        "raw_spearman_abs_ci_low": bootstrap_intervals["raw_spearman_abs_ci_low"],
                        "raw_spearman_abs_ci_high": bootstrap_intervals["raw_spearman_abs_ci_high"],
                        "raw_sign_agreement_ci_low": bootstrap_intervals["raw_sign_agreement_ci_low"],
                        "raw_sign_agreement_ci_high": bootstrap_intervals["raw_sign_agreement_ci_high"],
                        "controlled_spearman_abs": observed["controlled_spearman_abs"],
                        "controlled_minus_raw_abs": observed["controlled_minus_raw_abs"],
                        "null_mean_spearman_abs": null_mean_abs,
                        "null_std_spearman_abs": null_std_abs,
                        "null_p_spearman_abs": null_p_abs,
                        "null_p_sign_agreement": null_p_sign,
                        "null_z_spearman_abs": null_z_abs,
                        "null_n": float(len(null_df)),
                    }
                )

                if not null_df.empty:
                    null_export = null_df.copy()
                    null_export["run_name"] = run.name
                    null_export["dataset_label"] = perturb_eval.label
                    null_export["intervention"] = str(intervention)
                    null_rows.extend(null_export.to_dict(orient="records"))

    metrics_df = pd.DataFrame(metrics_rows)
    baseline_df = pd.DataFrame(baseline_rows)
    null_df = pd.DataFrame(null_rows)
    timepoint_df = _compute_timepoint_consistency(matched_tables)

    metrics_path = tables_dir / "counterfactual_consistency_metrics.tsv"
    baseline_path = tables_dir / "perturbation_baseline_metrics.tsv"
    null_path = tables_dir / "counterfactual_consistency_null.tsv"
    timepoint_path = tables_dir / "timepoint_consistency.tsv"

    metrics_df.to_csv(metrics_path, sep="\t", index=False)
    baseline_df.to_csv(baseline_path, sep="\t", index=False)
    null_df.to_csv(null_path, sep="\t", index=False)
    timepoint_df.to_csv(timepoint_path, sep="\t", index=False)

    _plot_consistency(metrics_df, figures_dir / "consistency_vs_null.png")
    _write_report(
        metrics_df=metrics_df,
        baseline_df=baseline_df,
        timepoint_df=timepoint_df,
        report_path=SUBPROJECT_ROOT / "reports" / "counterfactual_perturbation_consistency_report.md",
    )

    print(f"[done] metrics={metrics_path}", flush=True)
    print(f"[done] baseline={baseline_path}", flush=True)
    print(f"[done] null={null_path}", flush=True)
    print(f"[done] timepoint={timepoint_path}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run counterfactual perturbation consistency study.")
    parser.add_argument(
        "--runs",
        nargs="*",
        default=None,
        help="Optional run names subset: adamson dixit_tiny_neg shifrut",
    )
    parser.add_argument("--shuffle-n", type=int, default=200, help="Number of source-shuffle null draws")
    parser.add_argument("--bootstrap-n", type=int, default=500, help="Number of bootstrap draws for confidence intervals")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--skip-causal", action="store_true", help="Do not run causal interventions")
    parser.add_argument(
        "--skip-eval",
        action="store_true",
        help="Do not rerun perturbation baseline evaluation",
    )
    args = parser.parse_args()

    selected_runs = _select_runs(args.runs)
    run_study(
        selected_runs=selected_runs,
        shuffle_n=int(args.shuffle_n),
        bootstrap_n=int(args.bootstrap_n),
        seed=int(args.seed),
        skip_causal=bool(args.skip_causal),
        skip_eval=bool(args.skip_eval),
    )


if __name__ == "__main__":
    main()
