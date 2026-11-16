"""Metrics. Kept deliberately simple so every number is easy to explain.

Per held-out perturbation:
  pearson_all   correlation of predicted vs true effect over all genes
  pearson_top20 same, over the 20 genes the perturbation moves most (true |delta|)
  mse           mean squared error over all genes
  direction_top20  fraction of those top-20 genes whose sign is predicted correctly
The perturbed genes themselves are excluded (their change is trivially large).

Centered metrics (only when `center` is given):
  pearson_ctr_all    correlation of (true - center) vs (pred - center), all genes
  pearson_ctr_top20  same, over the 20 genes with the largest |true - center|
`center` is a reference every model gets for free: the mean effect of the TRAINING
perturbations (unseen singles) or the additive expectation (doubles). Most of a
perturbation's effect vs control is shared by all perturbations, so the uncentered
metrics reward a model for reproducing that shared response. The centered metrics
ask the harder question: did the model predict what is specific to THIS
perturbation? A prediction equal to the center (the mean baseline, or a network
model for a gene the network does not cover) says nothing specific and scores 0.
The center must come from training data only; centering on a mean that includes
the held-out perturbations makes every model look anti-correlated with the truth.

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


def per_perturbation(true: pd.DataFrame, pred: pd.DataFrame, std: pd.DataFrame | None, top: int = 20,
                     center: pd.Series | pd.DataFrame | None = None):
    """`center`: optional reference effect. A Series (genes) is used for every
    perturbation; a DataFrame (perturbations x genes) gives one reference each."""
    # work on arrays: label-based lookups per perturbation dominated the runtime
    T = true.values
    Y = pred.reindex(index=true.index, columns=true.columns).values
    S = None if std is None else std.reindex(index=true.index, columns=true.columns).values
    C = None
    if center is not None:
        C = (center.reindex(index=true.index, columns=true.columns).values if isinstance(center, pd.DataFrame)
             else np.broadcast_to(center.reindex(true.columns).values, T.shape))
    col = {g: j for j, g in enumerate(true.columns)}
    rows = []
    for i, p in enumerate(true.index):
        keep = np.ones(T.shape[1], dtype=bool)
        for g in parse_targets(p):
            if g in col:
                keep[col[g]] = False
        t, y = T[i, keep], Y[i, keep]
        order = np.argsort(-np.abs(t))[:top]
        r = dict(perturbation=p,
                 pearson_all=_safe_r(t, y),
                 pearson_top20=_safe_r(t[order], y[order]),
                 mse=float(np.mean((t - y) ** 2)),
                 direction_top20=float(np.mean(np.sign(t[order]) == np.sign(y[order]))))
        if C is not None:
            tc, yc = t - C[i, keep], y - C[i, keep]
            oc = np.argsort(-np.abs(tc))[:top]
            r["pearson_ctr_all"] = _ctr_r(tc, yc)
            r["pearson_ctr_top20"] = _ctr_r(tc[oc], yc[oc])
        if S is not None:
            s = S[i, keep]
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
    # a constant vector has no correlation; the relative test also catches a
    # prediction that equals the center up to floating-point noise
    sa, sb = np.std(a), np.std(b)
    if sa <= 1e-12 * max(1.0, np.abs(a).max()) or sb <= 1e-12 * max(1.0, np.abs(b).max()):
        return float("nan")
    return float(pearsonr(a, b)[0])



def _ctr_r(tc, yc) -> float:
    """Centered correlation. No deviation predicted on these genes -> 0, not NaN."""
    r = _safe_r(tc, yc)
    if np.isnan(r) and np.std(tc) > 0:
        return 0.0
    return r
