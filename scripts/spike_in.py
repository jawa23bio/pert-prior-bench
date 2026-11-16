"""How big would a real network effect have to be for the benchmark to see it?

A null result ("the network did not help") is only informative if the analysis
would have said "it helps" had the network mattered. The synthetic positive
control shows the code can find an effect in a dense toy network. This script
asks the same thing on YOUR data, with YOUR network and its real sparsity:

  1. take the real pseudobulk effects;
  2. for every perturbed TF the network covers, ADD a planted effect that
     points along that TF's network row:
         planted(g) = beta * |specific(g)| * unit(D[g, :])
     where |specific(g)| is the size of g's real perturbation-specific effect
     (its deviation from the mean effect). beta = 0.5 plants a network-shaped
     effect half as large as everything specific the perturbation really did;
  3. hand the analyses the same network and see who notices.

For each beta it reports, on the covered perturbations:
  v1      linear[both] minus linear[expression], pearson_top20, single split,
          original model settings (what the README concluded from)
  v2      linear[both] and direct[prior]: real network vs the regulator-permuted
          null, centered Pearson, with the empirical p-value

The smallest beta an analysis detects is its detection limit. Report it next to
the null result: "the network did not help, and an effect of relative size
beta >= X would have been detected".

    python scripts/spike_in.py --h5ad data/norman_pb.h5ad --already-log --min-cells 1 \
        --sep _ --n-hvg 2085 --prior data/collectri.csv --out results_spike_in

--edge-dropout 0.3 gives the analyses a degraded copy of the network (30% of
each TF's targets removed), which is the honest setting if you believe the
curated network is incomplete.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pertbench import protocol as P  # noqa: E402
from pertbench.data import kfold, parse_targets, split_singles_doubles  # noqa: E402
from pertbench.evaluate import per_perturbation  # noqa: E402
from pertbench.models import LinearModel  # noqa: E402
from pertbench.priors import prior_svd, propagation_matrix, regulators, relabel_all  # noqa: E402
from run_v2 import load  # noqa: E402


def plant(deltas: pd.DataFrame, singles: list[str], D: pd.DataFrame, beta: float) -> pd.DataFrame:
    """Add a network-shaped effect of relative size beta to every covered single."""
    Y = deltas.values.astype(float).copy()
    row = {p: i for i, p in enumerate(deltas.index)}
    col = {g: j for j, g in enumerate(deltas.columns)}
    mean = Y[[row[p] for p in singles]].mean(axis=0)
    Dv = D.reindex(columns=deltas.columns).fillna(0.0)
    regs = set(regulators(D))
    for p in singles:
        g = parse_targets(p)[0]
        if g not in regs:
            continue
        d = Dv.loc[g].values.copy()
        spec = Y[row[p]] - mean
        if g in col:
            d[col[g]], spec[col[g]] = 0.0, 0.0
        n = np.linalg.norm(d)
        if n > 0:
            Y[row[p]] += beta * np.linalg.norm(spec) * d / n
    return pd.DataFrame(Y, index=deltas.index, columns=deltas.columns)


def degrade(D: pd.DataFrame, frac: float, seed: int = 0) -> pd.DataFrame:
    """Zero a random share of each row's non-zero entries (an incomplete network)."""
    if frac <= 0:
        return D
    rng = np.random.default_rng(seed)
    V = D.values.copy()
    V[(V != 0) & (rng.random(V.shape) < frac)] = 0.0
    return pd.DataFrame(V, index=D.index, columns=D.columns)


def v1_verdict(deltas, singles, D, cov, k=10, k_prior=50) -> dict:
    """The original comparison: one 5-fold split, original model settings."""
    emb = {"real": prior_svd(D, k_prior), "shuffled": prior_svd(relabel_all(D, seed=1), k_prior)}
    rows = []
    for tr, te in kfold(singles, k=5):
        models = {"expression": LinearModel("expression", k=k, lam_grid=[0.1, 1, 10], n_boot=0),
                  "both": LinearModel("both", k=k, k_prior=k_prior, lam_grid=[0.1, 1, 10], n_boot=0,
                                      prior_embedding=emb["real"]),
                  "shuffled": LinearModel("both", k=k, k_prior=k_prior, lam_grid=[0.1, 1, 10], n_boot=0,
                                          prior_embedding=emb["shuffled"])}
        for name, m in models.items():
            mu, _ = m.fit(deltas.loc[tr]).predict(te)
            df = per_perturbation(deltas.loc[te], mu, None)[["perturbation", "pearson_top20"]]
            df["model"] = name
            rows.append(df)
    w = pd.concat(rows).pivot(index="perturbation", columns="model", values="pearson_top20").loc[cov]
    d, lo, hi = P.boot_ci((w["both"] - w["expression"]).values)
    return {"v1 both-expr": d, "v1 ci_lo": lo, "v1 ci_hi": hi,
            "v1 both-shuffled": float((w["both"] - w["shuffled"]).mean()),
            "v1 detects": bool(lo > 0)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--pert-key", default="perturbation")
    ap.add_argument("--control", default="control")
    ap.add_argument("--sep", default="+")
    ap.add_argument("--prior", default="collectri")
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--min-cells", type=int, default=30)
    ap.add_argument("--already-log", action="store_true")
    ap.add_argument("--gene-col", default="")
    ap.add_argument("--betas", default="0,0.1,0.25,0.5,1")
    ap.add_argument("--edge-dropout", type=float, default=0.0)
    ap.add_argument("--repeats", type=int, default=2)
    ap.add_argument("--n-null", type=int, default=19)
    ap.add_argument("--no-linear-prior", action="store_true")
    ap.add_argument("--out", default="results_spike_in")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    if args.n_null < 19:
        print(f"note: with {args.n_null} null draws the smallest possible p is {1 / (1 + args.n_null):.2f}; "
              "v2 cannot reach p <= 0.05, so it will report 'not detected'. Use --n-null 19 or more.")
    deltas, net = load(args)
    singles, _ = split_singles_doubles(deltas)
    D_true = propagation_matrix(net, list(deltas.columns))
    D_seen = degrade(D_true, args.edge_dropout)
    cov = [s for s, c in zip(singles, P.covered(D_seen, singles)) if c]
    cfg = P.Config(repeats=args.repeats, n_null=args.n_null, nulls=("regulator-permuted",),
                   with_linear_prior=not args.no_linear_prior)
    pool = [parse_targets(p)[0] for p in cov]
    nulls = P.make_nulls(D_seen, cfg, pool=pool)
    cov_mask = pd.Series([s in set(cov) for s in singles], index=singles)

    rows = []
    for beta in [float(b) for b in args.betas.split(",")]:
        print(f"beta = {beta}")
        Y = plant(deltas, singles, D_true, beta)
        r = {"beta": beta, **v1_verdict(Y, singles, D_seen, cov)}
        df = P.run_singles(Y, singles, D_seen, cfg, nulls, log=lambda *_: None)
        for model in sorted(set(df[df.network == "real"].model)):
            t = P.null_test(df, "pearson_ctr_top20", model, cov_mask).iloc[0]
            r[f"v2 {model} real"], r[f"v2 {model} null"], r[f"v2 {model} p"] = t["real"], t["null_mean"], t["p"]
            r[f"v2 {model} detects"] = bool(t["p"] <= 0.05)
        rows.append(r)
    tab = pd.DataFrame(rows).set_index("beta")
    tab.to_csv(out / "spike_in.csv")

    def limit(col):
        hit = tab.index[tab[col].values.astype(bool)]
        hit = [b for b in hit if b > 0]
        return f"beta >= {min(hit)}" if hit else "not detected at any beta tried"

    lines = ["# Spike-in sensitivity\n",
             f"{len(cov)} covered perturbations; edge dropout {args.edge_dropout}; "
             f"{cfg.repeats} x {cfg.folds}-fold CV; {cfg.n_null} null draws (smallest possible p = {1 / (1 + cfg.n_null):.3f}).\n",
             "beta = size of the planted network-shaped effect relative to each perturbation's real specific effect. "
             "A `detects` at beta = 0 is a false positive.\n",
             tab.round(3).to_markdown(), "", "## Detection limits\n"]
    for col in [c for c in tab.columns if c.endswith("detects")]:
        lines.append(f"- {col.replace(' detects', '')}: {limit(col)}")
    (out / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
