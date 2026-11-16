"""The v2 evaluation protocol, shared by run_v2.py and spike_in.py.

What changes relative to scripts/run_benchmark.py (which is left untouched so
every committed result still reproduces):

  splits     repeated K-fold, stratified so each fold holds the same share of
             network-covered perturbations (the original used one split; the
             covered TFs' mean score ranged from about 0.1 in one fold to 0.76
             in another, so which split you draw matters)
  metric     train-mean-centered Pearson next to the original metrics
  controls   many draws of two nulls that keep every TF a regulator, giving a
             null DISTRIBUTION and an empirical p-value instead of one shuffle
  models     adds direct[prior], which uses the network row itself
  stability  block feature scaling, no shrinkage in the output projection,
             ridge chosen by held-out MSE over a wider grid
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import partial

import numpy as np
import pandas as pd

from .data import parse_targets
from .evaluate import per_perturbation
from .models import AdditiveBaseline, DirectPriorModel, LinearModel, MeanBaseline
from .priors import permute_regulators, permute_targets, prior_svd, regulators, relabel_all

NULLS = {"regulator-permuted": permute_regulators, "target-permuted": permute_targets,
         "all-relabelled (v1 control)": relabel_all}
METRICS = ["pearson_top20", "pearson_ctr_top20", "pearson_ctr_all", "mse"]


@dataclass
class Config:
    k: int = 10
    k_prior: int = 50
    lam_grid: list[float] = field(default_factory=lambda: [0.01, 0.1, 1.0, 10.0, 100.0, 1000.0])
    folds: int = 5
    repeats: int = 5
    n_null: int = 20
    nulls: tuple[str, ...] = ("regulator-permuted", "target-permuted")
    n_boot: int = 0          # uncertainty is handled by the calibration stage
    with_linear_prior: bool = True


def stratified_kfold(items: list[str], strata: list[bool], k: int, seed: int):
    """K folds with each stratum dealt out evenly."""
    rng = np.random.default_rng(seed)
    folds = [[] for _ in range(k)]
    offset = 0
    for s in sorted(set(strata)):
        group = [x for x, t in zip(items, strata) if t == s]
        for i, x in enumerate(rng.permutation(group)):
            folds[(i + offset) % k].append(str(x))
        offset += len(group)
    return [([x for j, f in enumerate(folds) if j != i for x in f], folds[i]) for i in range(k)]


def covered(D: pd.DataFrame | None, perts: list[str]) -> list[bool]:
    """Is (any of) the perturbed gene(s) a regulator with measured targets?"""
    if D is None:
        return [False] * len(perts)
    regs = set(regulators(D))
    return [any(t in regs for t in parse_targets(p)) for p in perts]


def _linear(embedding, cfg, D=None, emb=None):
    return LinearModel(embedding, k=cfg.k, k_prior=cfg.k_prior, lam_grid=cfg.lam_grid, n_boot=cfg.n_boot,
                       prior_matrix=D, prior_embedding=emb, scaling="block", lam_out=1e-6, select="mse")


def prior_models(D: pd.DataFrame, cfg: Config, emb=None) -> list:
    """The models that read a network, built for one (real or null) network."""
    ms = [DirectPriorModel(D, n_boot=cfg.n_boot)]
    if cfg.with_linear_prior:
        ms.append(_linear("both", cfg, D=D, emb=emb if emb is not None else prior_svd(D, cfg.k_prior)))
    return ms


def make_nulls(D: pd.DataFrame, cfg: Config, pool: list[str] | None = None) -> dict[str, list]:
    """Recipes for the null networks: {null name: [zero-argument builders]}.

    Each null is a full genes x genes matrix (about 35 MB at 2,000 genes), so
    they are built one at a time when needed instead of being held in memory
    together; 99 draws of two nulls would otherwise need several GB.
    """
    out = {}
    for name in cfg.nulls:
        fn = NULLS[name]
        out[name] = [partial(fn, D, seed=1000 + i, pool=pool) if name == "regulator-permuted"
                     else partial(fn, D, seed=1000 + i) for i in range(cfg.n_null)]
    return out


def run_singles(deltas: pd.DataFrame, singles: list[str], D: pd.DataFrame | None, cfg: Config,
                nulls: dict[str, list] | None = None, log=print) -> pd.DataFrame:
    """Per-perturbation scores for every model, repeat and network.

    Columns: perturbation, the metrics, model, network ('none', 'real', or
    '<null name>'), draw (null index, -1 otherwise), rep, fold, covered.
    """
    cov = dict(zip(singles, covered(D, singles)))
    folds = [(rep, fold, tr, te) for rep in range(cfg.repeats)
             for fold, (tr, te) in enumerate(stratified_kfold(singles, [cov[s] for s in singles], cfg.folds, seed=rep))]
    rows = []

    def score(network, draw, models):
        """Fit and score one set of models on every fold of every repeat."""
        for rep, fold, tr, te in folds:
            train, test = deltas.loc[tr], deltas.loc[te]
            center = train.mean(axis=0)
            for m in models:
                mu, _ = m.fit(train).predict(te)
                df = per_perturbation(test, mu, None, center=center)
                df["model"], df["network"], df["draw"], df["rep"], df["fold"] = m.name, network, draw, rep, fold
                rows.append(df)

    score("none", -1, [MeanBaseline(), _linear("expression", cfg)])
    if D is not None:
        score("real", -1, prior_models(D, cfg))
        log("  real network done")
        for name, builders in (nulls or {}).items():
            for i, build in enumerate(builders):
                score(name, i, prior_models(build(), cfg))   # one null in memory at a time
            log(f"  {name}: {len(builders)} draws done")
    out = pd.concat(rows, ignore_index=True)
    out["covered"] = out["perturbation"].map(cov)
    return out


def run_doubles(deltas: pd.DataFrame, singles: list[str], doubles: list[str], D: pd.DataFrame | None,
                cfg: Config) -> pd.DataFrame:
    """Train on all singles, predict every double (same experiment as v1, with the
    stable model settings). The additive baseline is the bar to clear."""
    train = deltas.loc[singles]
    models = [MeanBaseline(), AdditiveBaseline(), _linear("expression", cfg)]
    if D is not None:
        models += prior_models(D, cfg)
    rows = []
    for m in models:
        mu, _ = m.fit(train).predict(doubles)
        mu.columns = deltas.columns
        df = per_perturbation(deltas.loc[doubles], mu, None)
        df["model"] = m.name
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


# ---------------------------------------------------------------- summaries
def per_pert_means(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """perturbations x (model, network): mean over repeats (and null draws)."""
    return df.pivot_table(index="perturbation", columns=["model", "network"], values=metric, aggfunc="mean")


def boot_ci(x: np.ndarray, n: int = 5000, seed: int = 0) -> tuple[float, float, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    bs = rng.choice(x, (n, len(x))).mean(axis=1)
    return float(x.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def table(df: pd.DataFrame, metric: str, subsets: dict[str, pd.Series]) -> pd.DataFrame:
    """Mean [95% CI over perturbations] per model/network and subset."""
    w = per_pert_means(df, metric)
    out = {}
    for (model, network) in w.columns:
        label = model if network in {"none", "real"} else f"{model}, {network}"
        out[label] = {}
        for name, mask in subsets.items():
            m, lo, hi = boot_ci(w.loc[mask.reindex(w.index).fillna(False).values, (model, network)].values)
            f = "{:.4f}" if metric == "mse" else "{:.3f}"
            out[label][name] = "n/a" if np.isnan(m) else f"{f.format(m)} [{f.format(lo)}, {f.format(hi)}]"
    return pd.DataFrame(out).T


def null_test(df: pd.DataFrame, metric: str, model: str, subset: pd.Series | None = None) -> pd.DataFrame:
    """Real network vs each null: is the real score outside the null distribution?

    Statistic = mean of `metric` over perturbations (in `subset`) and repeats.
    p = (1 + #null draws with statistic >= real) / (1 + #draws), one-sided.
    """
    d = df[df["model"] == model]
    if subset is not None:
        d = d[d["perturbation"].isin(subset.index[subset.values])]
    real = d[d["network"] == "real"][metric].mean()
    rows = []
    for name in sorted(set(d["network"]) - {"real", "none"}):
        stats = d[d["network"] == name].groupby("draw")[metric].mean().values
        rows.append(dict(model=model, null=name, metric=metric, real=real, null_mean=stats.mean(),
                         null_lo=np.percentile(stats, 2.5), null_hi=np.percentile(stats, 97.5),
                         n_draws=len(stats), p=(1 + int(np.sum(stats >= real))) / (1 + len(stats))))
    return pd.DataFrame(rows)
