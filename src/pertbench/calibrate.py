"""Prediction intervals that are meant to contain the truth.

The original benchmark built intervals from the bootstrap spread of the fitted
mean (prediction +/- 1.28 * bootstrap std) and found they covered ~10% of
outcomes instead of 80%. That number is real, but it is not a discovery about
the data: a bootstrap over training perturbations measures how much the FIT
moves when the training set changes. It leaves out everything the model cannot
explain (misspecification and measurement noise), which is almost all of the
error here. It is a confidence interval for the model's mean, not a prediction
interval for a new perturbation, so it under-covers by construction. The mean
baseline looked calibrated (0.83) only because its "std" is the spread of
training effects, which IS a predictive spread.

This module builds intervals from out-of-fold residuals instead:

  1. inner K-fold on the training perturbations gives a residual for every
     training perturbation from a model that never saw it;
  2. sigma_j = RMS residual of gene j (some genes are simply noisier);
  3. q = the 80th percentile of |residual| / sigma_j, pooled over perturbations
     and genes (split-conformal, CV+ style);
  4. interval for a new perturbation: prediction +/- q * sigma_j.

Coverage is then ~80% whenever held-out perturbations resemble the training
ones (exchangeability), whatever the model. With coverage fixed, the question
worth asking changes from "is it 80%?" to:

  - width: how wide do honest intervals have to be? A model that knows more
    earns narrower ones.
  - ranking: does the model's own uncertainty say which perturbations it will
    get wrong (err_vs_unc)? Per-gene sigma cannot answer this because it is the
    same for every perturbation, so the per-perturbation signal has to come from
    the model (here: bootstrap spread, used only to RANK perturbations).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .data import kfold, parse_targets

Z80 = 1.2816


def oof_residuals(make_model, train: pd.DataFrame, k: int = 5, seed: int = 0) -> pd.DataFrame:
    """Residual (truth - prediction) for every training perturbation, each from a
    model fitted without it. `make_model()` must return a fresh, unfitted model."""
    out = []
    for tr, te in kfold(list(train.index), k=k, seed=seed):
        m = make_model().fit(train.loc[tr])
        mu, _ = m.predict(te)
        out.append(train.loc[te] - mu.values)
    return pd.concat(out).loc[train.index]


def _mask_own_gene(R: pd.DataFrame) -> np.ndarray:
    """True where a residual may be used: not the perturbed gene's own column."""
    ok = np.ones(R.shape, dtype=bool)
    for i, p in enumerate(R.index):
        for t in parse_targets(p):
            if t in R.columns:
                ok[i, R.columns.get_loc(t)] = False
    return ok


def conformal_scale(resid: pd.DataFrame, level: float = 0.8) -> tuple[pd.Series, float]:
    """(sigma per gene, q). Interval half-width for gene j is q * sigma_j."""
    ok = _mask_own_gene(resid)
    r = np.where(ok, resid.values, np.nan)
    sigma = np.sqrt(np.nanmean(r**2, axis=0))
    floor = np.nanmedian(sigma) * 1e-3 + 1e-12
    sigma = np.maximum(np.nan_to_num(sigma, nan=floor), floor)
    scores = (np.abs(r) / sigma)[ok]
    n = scores.size
    q = float(np.quantile(scores, min(1.0, np.ceil((n + 1) * level) / n)))
    return pd.Series(sigma, index=resid.columns), q


def interval_report(true: pd.DataFrame, pred: pd.DataFrame, half_width: pd.DataFrame | pd.Series,
                    unc: pd.Series | None = None) -> pd.DataFrame:
    """Per perturbation: coverage, mean half-width, RMSE, and (optionally) the
    model's own per-perturbation uncertainty score."""
    ok = _mask_own_gene(true)
    err = np.abs(true.values - pred.values)
    hw = (half_width.reindex(true.columns).values[None, :] * np.ones_like(err)
          if isinstance(half_width, pd.Series) else half_width.values)
    rows = []
    for i, p in enumerate(true.index):
        m = ok[i]
        r = dict(perturbation=p,
                 coverage80=float(np.mean(err[i, m] <= hw[i, m])),
                 half_width=float(np.mean(hw[i, m])),
                 rmse=float(np.sqrt(np.mean(err[i, m] ** 2))))
        if unc is not None:
            r["unc"] = float(unc.loc[p])
        rows.append(r)
    return pd.DataFrame(rows)


def err_vs_unc(df: pd.DataFrame) -> float:
    """Spearman correlation across perturbations between the model's uncertainty
    score and its actual error. Positive = it is less sure where it is more wrong."""
    if "unc" not in df or df["unc"].nunique() < 3:
        return float("nan")
    rho = spearmanr(df["unc"], df["rmse"]).correlation
    return float(rho) if np.isfinite(rho) else float("nan")
