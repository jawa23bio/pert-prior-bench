"""Metrics. Kept deliberately simple so every number is easy to explain.

Per held-out perturbation:
  pearson_all   correlation of predicted vs true effect over all genes
  pearson_top20 same, over the 20 genes the perturbation moves most (true |delta|)
  mse           mean squared error over all genes
  direction_top20  fraction of those top-20 genes whose sign is predicted correctly
The perturbed genes themselves are excluded (their change is trivially large).

Calibration (for models that report a std):
  coverage80    fraction of (perturbation, gene) truths inside the predicted 80%
                interval; a calibrated model gives ~0.80
  err_vs_std    Spearman correlation, across perturbations, between mean predicted
                std and actual error; positive means "says unsure when wrong"
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from .data import parse_targets

Z80 = 1.2816


def per_perturbation(true: pd.DataFrame, pred: pd.DataFrame, std: pd.DataFrame | None, top: int = 20):
    rows = []
    for p in true.index:
        cols = [g for g in true.columns if g not in set(parse_targets(p))]
        t, y = true.loc[p, cols].values, pred.loc[p, cols].values
        order = np.argsort(-np.abs(t))[:top]
        r = dict(perturbation=p,
                 pearson_all=_safe_r(t, y),
                 pearson_top20=_safe_r(t[order], y[order]),
                 mse=float(np.mean((t - y) ** 2)),
                 direction_top20=float(np.mean(np.sign(t[order]) == np.sign(y[order]))))
        if std is not None:
            s = std.loc[p, cols].values
            r["coverage80"] = float(np.mean(np.abs(t - y) <= Z80 * s))
            r["mean_std"] = float(np.mean(s))
        rows.append(r)
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame) -> dict:
    out = {c: float(df[c].mean()) for c in ["pearson_all", "pearson_top20", "mse", "direction_top20"]}
    if "coverage80" in df:
        out["coverage80"] = float(df["coverage80"].mean())
        rho = spearmanr(df["mean_std"], df["mse"]).correlation
        out["err_vs_std"] = float(rho) if np.isfinite(rho) else float("nan")
    return out


def _safe_r(a, b) -> float:
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(pearsonr(a, b)[0])
