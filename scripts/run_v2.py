"""v2 benchmark: same question, a design that can answer it.

    python scripts/run_v2.py --h5ad data/norman_pb.h5ad --already-log --min-cells 1 \
        --pert-key perturbation --control control --sep _ --n-hvg 2085 \
        --prior data/collectri.csv --out results_v2

Writes to --out:
  summary.md                   every table below, ready to paste
  singles_per_perturbation.csv one row per perturbation x model x network x draw x repeat (large)
  singles_mean.csv             the same, averaged over repeats and null draws (small; commit this one)
  null_tests.csv               real network vs each null distribution
  null_draws.csv               mean score of every null draw (covered perturbations)
  diagnostics.csv              per covered TF: measured targets, retained, transferable
  calibration.csv              per perturbation interval coverage / width (with --calibrate)
  doubles_per_perturbation.csv unseen doubles

Read src/pertbench/protocol.py for what changed relative to run_benchmark.py and why.
"""
import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pertbench import protocol as P  # noqa: E402
from pertbench.calibrate import Z80, conformal_scale, err_vs_unc, interval_report, oof_residuals  # noqa: E402
from pertbench.data import parse_targets, pseudobulk_deltas, split_singles_doubles  # noqa: E402
from pertbench.models import DirectPriorModel, MeanBaseline  # noqa: E402
from pertbench.priors import (load_collectri, prior_svd, propagation_matrix, regulators,  # noqa: E402
                              relabel_all, retained, transferability)


def load(args):
    adata = ad.read_h5ad(args.h5ad)
    if args.gene_col:
        adata.var_names = adata.var[args.gene_col].astype(str).values
        adata.var_names_make_unique()
    deltas, _ = pseudobulk_deltas(adata, args.pert_key, args.control, args.n_hvg,
                                  args.min_cells, args.sep, args.already_log)
    net = None
    if args.prior != "none":
        net = load_collectri() if args.prior == "collectri" else pd.read_csv(args.prior)
    return deltas, net


def diagnostics(D, net, deltas, singles, cfg) -> tuple[pd.DataFrame, dict]:
    """Can the network matter at all here? Computed before any model is fitted."""
    genes = set(deltas.columns)
    pert = [parse_targets(p)[0] for p in singles]
    regs = set(regulators(D))
    cov = [g for g in pert if g in regs]
    direct = net[net["source"].isin(cov) & net["target"].isin(genes) & (net["source"] != net["target"])]
    n_targets = direct.groupby("source")["target"].nunique().reindex(cov).fillna(0).astype(int)
    Rp, _ = prior_svd(D, cfg.k_prior)
    keep = retained(D, cfg.k_prior, cov)
    # transferable share, averaged over the folds of the first repeat
    tf_emb, tf_raw = [], []
    strata = [g in regs for g in pert]
    for tr, te in P.stratified_kfold(singles, strata, cfg.folds, seed=0):
        trg, teg = [parse_targets(p)[0] for p in tr], [parse_targets(p)[0] for p in te]
        tf_emb.append(transferability(Rp, trg, teg))
        tf_raw.append(transferability(D, trg, teg))
    tf_emb, tf_raw = pd.concat(tf_emb), pd.concat(tf_raw)
    per = pd.DataFrame({"measured_direct_targets": n_targets, "retained_by_rank_k_svd": keep.reindex(cov),
                        "transferable_in_svd_space": tf_emb.reindex(cov),
                        "transferable_in_network_space": tf_raw.reindex(cov)})
    per.index.name = "perturbed_tf"
    old = relabel_all(D, seed=1)
    old_regs = set(regulators(old))
    info = {
        "n_genes": int(deltas.shape[1]),
        "n_singles": len(singles),
        "regulators_among_measured_genes": len(regs),
        "perturbed_genes_covered": len(cov),
        "median_measured_direct_targets_per_covered_tf": float(n_targets.median()) if len(cov) else 0.0,
        "covered_tfs_with_5plus_targets": int((n_targets >= 5).sum()),
        "median_share_of_prior_row_kept_by_svd": round(float(keep.median()), 3) if len(keep) else None,
        "median_transferable_share_svd_space": round(float(tf_emb.median()), 3) if len(tf_emb) else None,
        "median_transferable_share_network_space": round(float(tf_raw.median()), 3) if len(tf_raw) else None,
        "v1_control_covered_tfs_given_any_targets": f"{sum(g in old_regs for g in cov)} of {len(cov)}",
    }
    return per, info


def calibration(deltas, singles, D, cfg, level=0.8) -> pd.DataFrame:
    """Bootstrap intervals (v1) vs out-of-fold conformal intervals, first repeat only."""
    cfg_b = P.Config(**{**cfg.__dict__, "n_boot": 20})
    makers = {"mean": lambda: MeanBaseline(), "linear[expression]": lambda: P._linear("expression", cfg_b)}
    if D is not None:
        makers["direct[prior]"] = lambda: DirectPriorModel(D, n_boot=20)
    cov = P.covered(D, singles)
    rows = []
    for fold, (tr, te) in enumerate(P.stratified_kfold(singles, cov, cfg.folds, seed=0)):
        train, test = deltas.loc[tr], deltas.loc[te]
        for name, make in makers.items():
            mu, sd = make().fit(train).predict(te)
            sigma, q = conformal_scale(oof_residuals(make, train, k=cfg.folds, seed=fold), level)
            unc = sd.mean(axis=1)
            for kind, hw in [("v1 (1.28 x model std)", Z80 * sd), ("out-of-fold conformal", q * sigma)]:
                df = interval_report(test, mu, hw, unc)
                df["model"], df["interval"], df["fold"] = name, kind, fold
                rows.append(df)
    return pd.concat(rows, ignore_index=True)


def md(df: pd.DataFrame) -> str:
    if isinstance(df.index, pd.MultiIndex):  # "a / b" reads better than a tuple
        df = df.copy()
        df.index = [" / ".join(map(str, t)) for t in df.index]
    return df.to_markdown()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--pert-key", default="perturbation")
    ap.add_argument("--control", default="control")
    ap.add_argument("--sep", default="+")
    ap.add_argument("--prior", default="collectri", help="'collectri', a CSV path, or 'none'")
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--min-cells", type=int, default=30)
    ap.add_argument("--already-log", action="store_true")
    ap.add_argument("--gene-col", default="")
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--k-prior", type=int, default=50)
    ap.add_argument("--lam-grid", default="0.01,0.1,1,10,100,1000")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--repeats", type=int, default=5, help="independent CV splits")
    ap.add_argument("--n-null", type=int, default=20, help="draws per null network")
    ap.add_argument("--no-linear-prior", action="store_true",
                    help="skip linear[both] (the slow model) and keep only direct[prior]")
    ap.add_argument("--calibrate", action="store_true", help="also compare interval types")
    ap.add_argument("--out", default="results_v2")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = P.Config(k=args.k, k_prior=args.k_prior, lam_grid=[float(x) for x in args.lam_grid.split(",")],
                   folds=args.folds, repeats=args.repeats, n_null=args.n_null,
                   with_linear_prior=not args.no_linear_prior)

    deltas, net = load(args)
    singles, doubles = split_singles_doubles(deltas)
    print(f"{len(singles)} singles, {len(doubles)} doubles, {deltas.shape[1]} genes")
    lines = ["# Results (v2 protocol)\n"]

    D = nulls = None
    if net is not None:
        D = propagation_matrix(net, list(deltas.columns))
        per, info = diagnostics(D, net, deltas, singles, cfg)
        per.to_csv(out / "diagnostics.csv")
        lines += ["## Before fitting anything: can the network matter here?\n", "```", json.dumps(info, indent=2), "```\n"]
        print(json.dumps(info, indent=2))
        pert_regs = [g for g in dict.fromkeys(parse_targets(p)[0] for p in singles) if g in set(regulators(D))]
        nulls = P.make_nulls(D, cfg, pool=pert_regs)

    df = P.run_singles(deltas, singles, D, cfg, nulls)
    df.to_csv(out / "singles_per_perturbation.csv", index=False)
    # compact version for the repo: mean over repeats and null draws
    (df.groupby(["perturbation", "covered", "model", "network"], sort=False)[P.METRICS + ["pearson_all"]]
       .mean().reset_index().to_csv(out / "singles_mean.csv", index=False))
    cov = df.drop_duplicates("perturbation").set_index("perturbation")["covered"]
    subsets = {f"all {len(cov)}": cov | True}
    if D is not None:
        subsets.update({f"{int(cov.sum())} covered": cov, f"{int((~cov).sum())} not covered": ~cov})

    lines.append(f"## Unseen singles ({cfg.repeats} x {cfg.folds}-fold CV; mean [95% CI over perturbations])\n")
    for metric, note in [
        ("pearson_ctr_top20", "**Primary.** Perturbation-specific effect: (truth - training mean) vs "
                              "(prediction - training mean), on the 20 genes that deviate most. 0 = no better than the mean."),
        ("pearson_ctr_all", "Same, all genes."),
        ("pearson_top20", "The v1 metric (effect vs control). Dominated by the response shared by all perturbations."),
        ("mse", "Mean squared error, all genes."),
    ]:
        lines += [f"### {metric}\n", note + "\n", md(P.table(df, metric, subsets)), ""]

    if D is not None:
        tests = pd.concat([P.null_test(df, m, model, cov) for model in sorted(set(df[df.network == "real"].model))
                           for m in ["pearson_ctr_top20", "pearson_top20"]], ignore_index=True)
        tests.to_csv(out / "null_tests.csv", index=False)
        # one row per null draw (mean over covered perturbations and repeats), for figures
        d = df[df["perturbation"].isin(cov.index[cov.values])]
        (d.groupby(["model", "network", "draw"])[["pearson_ctr_top20", "pearson_top20"]].mean()
           .reset_index().to_csv(out / "null_draws.csv", index=False))
        lines += ["## Is the real network better than a wrong one? (covered perturbations)\n",
                  "Statistic: mean score over covered perturbations. `p` = share of null draws scoring at least "
                  "as high as the real network (one-sided; the smallest possible value is 1/(1+draws)).\n",
                  md(tests.round(4).set_index(["model", "null", "metric"])), ""]
        # paired: does reading the network beat not reading it?
        w = P.per_pert_means(df, "pearson_top20")
        base = w[("linear[expression]", "none")]
        rows = []
        for model in sorted(set(df[df.network == "real"].model)):
            d = (w[(model, "real")] - base)[cov.reindex(w.index).values].dropna()
            m, lo, hi = P.boot_ci(d.values)
            nz = d[d != 0]
            rows.append(dict(comparison=f"{model} - linear[expression]", mean_diff=m, ci_lo=lo, ci_hi=hi,
                             wins=f"{int((d > 0).sum())}/{len(d)}",
                             wilcoxon_p=wilcoxon(nz).pvalue if len(nz) > 5 else np.nan))
        lines += ["## Paired differences on covered perturbations (pearson_top20)\n",
                  md(pd.DataFrame(rows).round(4).set_index("comparison")), ""]

    if args.calibrate:
        cal = calibration(deltas, singles, D, cfg)
        cal.to_csv(out / "calibration.csv", index=False)
        tab = cal.groupby(["model", "interval"], sort=False).agg(
            coverage80=("coverage80", "mean"), half_width=("half_width", "mean"), rmse=("rmse", "mean"))
        tab["err_vs_unc"] = [err_vs_unc(g) for _, g in cal.groupby(["model", "interval"], sort=False)]
        lines += ["## Uncertainty (first CV repeat)\n",
                  "`coverage80` should be 0.80. `half_width` is the price of honest coverage (narrower = the model "
                  "knows more). `err_vs_unc`: does the model's own per-perturbation uncertainty rank its errors? "
                  "(Spearman; identical for both interval types because both rank by the model's own std. For the mean "
                  "baseline that std is the spread of training effects, which is why its v1 interval was never far off.)\n",
                  md(tab.round(3)), ""]

    if doubles:
        dd = P.run_doubles(deltas, singles, doubles, D, cfg)
        dd.to_csv(out / "doubles_per_perturbation.csv", index=False)
        tab = dd.groupby("model", sort=False)[["pearson_top20", "pearson_all", "mse"]].mean()
        lines += ["## Unseen doubles (trained on all singles)\n",
                  "The additive baseline has seen both singles, so it is the bar to clear.\n", md(tab.round(3)), ""]

    (out / "summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
