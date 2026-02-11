from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_counterfactual_consistency_study as study


def test_residualize_removes_linear_signal() -> None:
    x1 = np.arange(20, dtype=float)
    x2 = np.linspace(-1.0, 1.0, 20)
    cov = np.column_stack([x1, x2])
    y = 3.0 + 2.0 * x1 - 0.5 * x2

    resid = study._residualize(y, cov)
    assert np.nanmax(np.abs(resid)) < 1e-8


def test_parse_source_filters_controls_and_multigene() -> None:
    alias_map: dict[str, str] = {}

    parsed = study._parse_source(
        label="ctrl+STAT1",
        delimiter="+",
        allow_multi=False,
        alias_map=alias_map,
        control_labels=["ctrl"],
    )
    assert parsed == "STAT1"

    multi_blocked = study._parse_source(
        label="STAT1+IRF1",
        delimiter="+",
        allow_multi=False,
        alias_map=alias_map,
        control_labels=["ctrl"],
    )
    assert multi_blocked == ""


def test_source_shuffle_null_runs() -> None:
    matched = pd.DataFrame(
        {
            "source": ["A", "A", "B", "B"],
            "target": ["X", "Y", "X", "Y"],
            "effect_mean": [1.0, 0.5, -0.8, -0.4],
            "delta": [0.9, 0.4, -0.7, -0.3],
            "baseline_expr": [1.0, 1.0, 1.0, 1.0],
            "composition_shift": [0.1, 0.1, 0.2, 0.2],
            "n_cells": [30, 30, 40, 40],
        }
    )
    rng = np.random.default_rng(123)

    null = study._compute_source_shuffle_null(matched, n_shuffle=10, rng=rng)
    assert not null.empty
    assert {"null_spearman_abs", "null_sign_agreement"}.issubset(null.columns)
