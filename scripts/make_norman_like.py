"""A stand-in for data/norman_pb.h5ad that needs no download.

Same SHAPE as the real analysis, simulated responses:
  - the same perturbation names as the committed Norman results (105 singles,
    131 doubles), so CollecTRI covers the same 43 perturbed TFs;
  - about 2,085 genes, sampled so that a covered TF has a median of ~4 measured
    CollecTRI targets, as reported for the real gene set;
  - effects = a response shared by all perturbations + a few low-rank programs
    + noise. NO network effect is simulated: the network is uninformative here
    by construction.

Use it to (a) run the whole v2 pipeline end to end without the 700 MB screen,
and (b) with scripts/spike_in.py, ask how large a real network effect would
have to be before each analysis can see it, on data with CollecTRI's actual
sparsity. It is not a substitute for the real screen: the gene set is sampled,
not the real one, and the responses are made up.

    python scripts/make_norman_like.py --out data/norman_like_pb.h5ad
"""
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--prior", default=str(ROOT / "data/collectri.csv"))
    ap.add_argument("--names-from", default=str(ROOT / "results_norman"),
                    help="folder with singles_/doubles_per_perturbation.csv (for the perturbation names)")
    ap.add_argument("--n-genes", type=int, default=2085)
    ap.add_argument("--target-median", type=float, default=4, help="median measured targets per covered TF")
    ap.add_argument("--shared", type=float, default=1.0, help="size of the response shared by all perturbations")
    ap.add_argument("--specific", type=float, default=0.8, help="size of the perturbation-specific programs")
    ap.add_argument("--noise", type=float, default=0.35)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=str(ROOT / "data/norman_like_pb.h5ad"))
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    net = pd.read_csv(a.prior)
    singles = list(dict.fromkeys(pd.read_csv(Path(a.names_from) / "singles_per_perturbation.csv")["perturbation"]))
    doubles = list(dict.fromkeys(pd.read_csv(Path(a.names_from) / "doubles_per_perturbation.csv")["perturbation"]))
    pert = sorted({g for p in singles + doubles for g in p.split("+")})
    cov = [g for g in pert if g in set(net["source"])]

    # gene set: perturbed genes + a mix of network genes and anonymous genes, with the
    # mix tuned so covered TFs have about --target-median measured targets
    universe = np.array(sorted((set(net["source"]) | set(net["target"])) - set(pert)))
    order = rng.permutation(universe)
    n_free = a.n_genes - len(pert)

    def median_targets(n_real):
        genes = set(pert) | set(order[:n_real])
        sub = net[net["source"].isin(cov) & net["target"].isin(genes) & (net["source"] != net["target"])]
        return sub.groupby("source")["target"].nunique().reindex(cov).fillna(0).median()

    lo, hi = 0, min(n_free, len(order))
    while lo < hi:  # smallest n_real reaching the target median
        mid = (lo + hi) // 2
        if median_targets(mid) >= a.target_median:
            hi = mid
        else:
            lo = mid + 1
    n_real = lo
    picked = list(order[:n_real])
    # in the real gene set every one of the perturbed TFs has at least one measured target
    have = set(pert) | set(picked)
    for tf in cov:
        t = net[(net["source"] == tf) & (net["target"] != tf)]["target"]
        if not (set(t) & have) and len(t):
            g = sorted(t)[int(rng.integers(len(t)))]
            picked.append(g); have.add(g)
    n_real = len(picked)
    genes = pert + picked + [f"ANON{i:04d}" for i in range(n_free - n_real)]
    G = len(genes)

    shared = rng.normal(0, 1, G) * (rng.random(G) < 0.3)          # a third of genes respond to anything
    programs = rng.normal(0, 1, (5, G)) * (rng.random((5, G)) < 0.1)
    scale = 0.08                                                   # log-expression units, similar to real deltas

    def effect(_name):
        return scale * (a.shared * rng.uniform(0.5, 1.5) * shared
                        + a.specific * rng.normal(0, 1, 5) @ programs
                        + a.noise * rng.normal(0, 1, G))

    single_fx = {p: effect(p) for p in singles}
    rows, labels = [np.zeros(G)], ["control"]
    for p in singles:
        rows.append(single_fx[p]); labels.append(p)
    for d in doubles:
        x, y = d.split("+")
        base = single_fx.get(x, effect(x)) + single_fx.get(y, effect(y))
        rows.append(0.85 * base + scale * a.noise * rng.normal(0, 1, G)); labels.append(f"{x}_{y}")
    X = np.vstack(rows) + rng.uniform(0.2, 2.5, G)                 # add a baseline expression level
    # the perturbed gene itself goes up (CRISPRa)
    gi = {g: i for i, g in enumerate(genes)}
    for r, lab in enumerate(labels):
        if lab != "control":
            for g in lab.split("_"):
                X[r, gi[g]] += 1.5

    res = ad.AnnData(X=X.astype(np.float32),
                     obs=pd.DataFrame({"perturbation": labels, "n_cells": 300}, index=[f"pb{i}" for i in range(len(labels))]),
                     var=pd.DataFrame(index=genes))
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    res.write_h5ad(a.out)
    print(f"wrote {a.out}: {len(singles)} singles, {len(doubles)} doubles, {G} genes "
          f"({n_real} network genes, {n_free - n_real} anonymous); {len(cov)} perturbed TFs in the prior, "
          f"median measured targets "
          f"{net[net['source'].isin(cov) & net['target'].isin(set(genes)) & (net['source'] != net['target'])].groupby('source')['target'].nunique().reindex(cov).fillna(0).median():.0f}")


if __name__ == "__main__":
    main()
