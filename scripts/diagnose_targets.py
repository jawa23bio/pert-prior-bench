"""Do CollecTRI targets respond when their TF is perturbed, and in the right direction?

For each single perturbation of a TF in the network:
- location: how many of its measured targets are among the 100 genes it moves most,
  against the number expected by chance;
- direction: how often a target moves the way the edge sign says (CRISPRa raises the
  TF, so an activating edge predicts the target goes up), against a baseline of
  random genes given the same signs.

The baseline matters: in Norman, most genes go down after a perturbation and about
90% of CollecTRI edges are activating, so raw agreement is low even for random genes.

    python scripts/diagnose_targets.py --h5ad data/norman_pb.h5ad --prior data/collectri.csv
"""
import argparse

import anndata as ad
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--h5ad", default="data/norman_pb.h5ad")
    ap.add_argument("--prior", default="data/collectri.csv")
    ap.add_argument("--pert-key", default="perturbation")
    ap.add_argument("--control", default="control")
    ap.add_argument("--sep", default="_")
    ap.add_argument("--top", type=int, default=100)
    ap.add_argument("--n-random", type=int, default=200)
    a = ap.parse_args()

    A = ad.read_h5ad(a.h5ad)
    lab = list(A.obs[a.pert_key].astype(str))
    X = np.asarray(A.X)
    genes = np.asarray(A.var_names)
    gi = {g: i for i, g in enumerate(genes)}
    ctrl = X[lab.index(a.control)]
    net = pd.read_csv(a.prior)
    singles = [x for x in lab if x != a.control and a.sep not in x]
    rng = np.random.default_rng(0)

    obs = exp = 0
    agree, rand_agree, rows = [], [], []
    for tf in singles:
        sub = net[(net.source == tf) & net.target.isin(gi) & (net.target != tf)]
        if sub.empty:
            continue
        d = X[lab.index(tf)] - ctrl
        d[gi[tf]] = 0
        idx = np.array([gi[t] for t in sub.target])
        w = np.sign(sub.weight.values)
        top = set(np.argsort(-np.abs(d))[: a.top])
        hits = sum(i in top for i in idx)
        obs += hits
        exp += len(idx) * a.top / (len(genes) - 1)
        moved = d[idx] != 0
        agree.append((np.sign(d[idx]) == w)[moved])
        others = np.setdiff1d(np.arange(len(genes)), [gi[tf]])
        for _ in range(a.n_random):
            r = rng.choice(others, len(idx), replace=False)
            m = d[r] != 0
            if m.any():
                rand_agree.append((np.sign(d[r]) == w)[m].mean())
        rows.append((tf, len(idx), hits))

    per = pd.DataFrame(rows, columns=["tf", "measured_targets", f"in_top{a.top}"])
    agree = np.concatenate(agree)
    up = np.mean([(X[lab.index(t)] - ctrl > 0).mean() for t in singles])
    print(f"TFs with measured targets: {len(per)}")
    print(f"location: {obs} targets in top {a.top}, {exp:.1f} expected, {obs / exp:.2f}-fold")
    top3 = per[f"in_top{a.top}"].nlargest(3).sum() / obs
    print(f"  TFs contributing any hit: {(per[f'in_top{a.top}'] > 0).sum()}; top 3 TFs give {top3:.0%} of hits")
    print(f"direction: {agree.mean():.1%} of {len(agree)} moved targets agree with the edge sign")
    print(f"  random genes with the same signs: {np.mean(rand_agree):.1%}")
    print(f"  genes going up per perturbation: {up:.1%}; activating edges: {(net.weight > 0).mean():.1%}")
    print(per.sort_values(f"in_top{a.top}", ascending=False).head(10).to_string(index=False))


if __name__ == "__main__":
    main()
