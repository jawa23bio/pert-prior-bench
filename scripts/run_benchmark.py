"""Run the full comparison and write tables to --out.

Experiment 1 (unseen single perturbations): 5-fold CV over single-gene
perturbations. Train on 4 folds, predict the held-out fold.

Experiment 2 (unseen combinations): train on all singles, predict every double.

Models: mean baseline, additive (doubles only), and the linear model with three
perturbation embeddings (expression / prior / both) plus a SHUFFLED-prior control.
The shuffled control matters: if 'both' beats 'expression' only because it has
more parameters, the shuffled prior will help just as much.
"""
import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pertbench.data import kfold, parse_targets, pseudobulk_deltas, split_singles_doubles  # noqa: E402
from pertbench.evaluate import per_perturbation, summarize  # noqa: E402
from pertbench.models import AdditiveBaseline, LinearModel, MeanBaseline  # noqa: E402
from pertbench.priors import coverage, load_collectri, propagation_matrix  # noqa: E402


def shuffled(D: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Same network structure, gene labels permuted: keeps size, destroys meaning."""
    perm = np.random.default_rng(seed).permutation(D.index)
    S = D.copy()
    S.index, S.columns = perm, perm
    return S.loc[D.index, D.columns]


def build_models(D, Dshuf, args, doubles=False):
    grid = [float(x) for x in args.lam_grid.split(",")] if args.lam_grid else None
    kw = dict(k=args.k, k_prior=args.k_prior, lam=args.lam, n_boot=args.n_boot, lam_grid=grid)
    ms = [MeanBaseline()]
    if doubles:
        ms.append(AdditiveBaseline())
    ms.append(LinearModel("expression", **kw))
    if D is not None:
        ms.append(LinearModel("prior", prior_matrix=D, **kw))
        ms.append(LinearModel("both", prior_matrix=D, **kw))
        shuf = LinearModel("both", prior_matrix=Dshuf, **kw)
        shuf.name = "linear[both, shuffled prior]"
        ms.append(shuf)
    return ms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--pert-key", default="perturbation")
    ap.add_argument("--control", default="control")
    ap.add_argument("--sep", default="+")
    ap.add_argument("--prior", default="collectri", help="'collectri', a CSV path, or 'none'")
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--min-cells", type=int, default=30)
    ap.add_argument("--k", type=int, default=10, help="expression PCs")
    ap.add_argument("--k-prior", type=int, default=50, help="prior SVD components")
    ap.add_argument("--lam", type=float, default=1.0, help="ridge penalty if --lam-grid is empty")
    ap.add_argument("--lam-grid", default="0.1,1,10",
                    help="comma list; picked per model by inner CV on training perts only ('' to disable)")
    ap.add_argument("--n-boot", type=int, default=20)
    ap.add_argument("--already-log", action="store_true")
    ap.add_argument("--gene-col", default="", help="var column holding gene symbols, if var_names are IDs")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.h5ad)
    if args.gene_col:
        adata.var_names = adata.var[args.gene_col].astype(str).values
        adata.var_names_make_unique()
    deltas, n_cells = pseudobulk_deltas(adata, args.pert_key, args.control, args.n_hvg,
                                        args.min_cells, args.sep, args.already_log)
    singles, doubles = split_singles_doubles(deltas)
    print(f"{len(singles)} singles, {len(doubles)} doubles, {deltas.shape[1]} genes")

    D = Dshuf = None
    info = {"n_singles": len(singles), "n_doubles": len(doubles), "n_genes": deltas.shape[1]}
    if args.prior != "none":
        net = load_collectri() if args.prior == "collectri" else pd.read_csv(args.prior)
        D = propagation_matrix(net, list(deltas.columns))
        Dshuf = shuffled(D, seed=1)
        info["prior_coverage_singles"] = coverage(net, [parse_targets(p)[0] for p in singles])
        print(f"prior covers {info['prior_coverage_singles']:.0%} of perturbed genes")

    # Experiment 1: unseen singles
    rows = []
    for fold, (tr, te) in enumerate(kfold(singles, k=5)):
        for m in build_models(D, Dshuf, args):
            m.fit(deltas.loc[tr])
            mu, sd = m.predict(te)
            df = per_perturbation(deltas.loc[te], mu, sd)
            df["model"], df["fold"] = m.name, fold
            rows.append(df)
    single_df = pd.concat(rows)
    single_df.to_csv(out / "singles_per_perturbation.csv", index=False)

    # Experiment 2: unseen doubles
    double_df = pd.DataFrame()
    if doubles:
        rows = []
        for m in build_models(D, Dshuf, args, doubles=True):
            m.fit(deltas.loc[singles])
            mu, sd = m.predict(doubles)
            df = per_perturbation(deltas.loc[doubles], mu, sd)
            df["model"] = m.name
            rows.append(df)
        double_df = pd.concat(rows)
        double_df.to_csv(out / "doubles_per_perturbation.csv", index=False)

    # summaries
    lines = ["# Results\n", "```", json.dumps(info, indent=2), "```\n"]
    summary = {}
    for label, df in [("Unseen single perturbations (5-fold CV)", single_df),
                      ("Unseen double perturbations (trained on singles)", double_df)]:
        if df.empty:
            continue
        tab = pd.DataFrame({m: summarize(g) for m, g in df.groupby("model", sort=False)}).T
        summary[label] = tab
        lines += [f"## {label}\n", tab.round(3).to_markdown(), ""]
    if D is not None:
        lines += _prior_split(single_df, D)
    (out / "summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


def _prior_split(df, D):
    """Does the prior help more for perturbed genes the prior actually knows about?"""
    has_edges = (D.abs().sum(axis=1) > 0)
    df = df.copy()
    df["in_prior"] = [bool(has_edges.get(parse_targets(p)[0], False)) for p in df["perturbation"]]
    tab = df.groupby(["in_prior", "model"], sort=False)["pearson_top20"].mean().unstack(0)
    tab.columns = [f"in_prior={c}" for c in tab.columns]
    return ["## Singles: pearson_top20 split by whether the perturbed gene is a regulator in the prior\n",
            tab.round(3).to_markdown(), ""]


if __name__ == "__main__":
    main()
