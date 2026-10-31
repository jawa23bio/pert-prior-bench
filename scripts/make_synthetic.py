"""Simulate a small Perturb-seq-like dataset with a KNOWN regulatory network.

Used as a positive control: if the true network drives responses, a prior that
resembles it should help, and a shuffled prior should not. If the pipeline does
not show that here, something is broken before you ever touch real data.

Writes:
  synthetic.h5ad     cells x genes counts, obs['perturbation']
  synthetic_prior.csv  noisy version of the true network (source,target,weight)
"""
import argparse

import anndata as ad
import numpy as np
import pandas as pd


def main(out_dir: str, seed: int = 0):
    rng = np.random.default_rng(seed)
    n_genes, n_tf, n_nontf, n_double = 300, 60, 20, 30
    genes = [f"G{i:03d}" for i in range(n_genes)]
    tfs = genes[:n_tf]

    # true network: each TF regulates ~15 genes (including other TFs) with signs
    edges = []
    for tf in tfs:
        for t in rng.choice(genes, size=15, replace=False):
            if t != tf:
                edges.append((tf, t, rng.choice([-1.0, 1.0], p=[0.3, 0.7])))
    true = pd.DataFrame(edges, columns=["source", "target", "weight"])

    idx = {g: i for i, g in enumerate(genes)}
    A = np.zeros((n_genes, n_genes))
    for s, t, w in edges:
        A[idx[s], idx[t]] = w
    A = A / np.maximum(np.abs(A).sum(1, keepdims=True), 1)

    common = rng.normal(0, 1, n_genes) * 0.05  # generic response to any perturbation

    def effect(targets):
        x = np.zeros(n_genes)
        for t in targets:
            x[idx[t]] = 6.0  # CRISPRa-like activation
        d, cur = np.zeros(n_genes), x.copy()
        for k in range(3):
            cur = cur @ A
            d += (0.6**k) * cur
        d = 1.5 * np.tanh(d)  # saturation makes doubles non-additive
        return d + x * 0.5 + common

    base = rng.uniform(0.5, 3.0, n_genes)
    singles = tfs + genes[n_tf:n_tf + n_nontf]
    pairs = set()
    while len(pairs) < n_double:
        a, b = sorted(rng.choice(tfs[:40], 2, replace=False))
        pairs.add((a, b))
    conditions = {"control": np.zeros(n_genes)}
    conditions.update({g: effect([g]) for g in singles})
    conditions.update({f"{a}+{b}": effect([a, b]) for a, b in pairs})

    X, labels = [], []
    for name, d in conditions.items():
        n = 600 if name == "control" else 60
        lam = np.exp(base + d)[None, :] * rng.gamma(5, 1 / 5, (n, 1))
        X.append(rng.poisson(lam))
        labels += [name] * n
    adata = ad.AnnData(np.vstack(X).astype(np.float32))
    adata.var_names = genes
    adata.obs["perturbation"] = labels
    adata.write_h5ad(f"{out_dir}/synthetic.h5ad")

    # noisy prior: drop 30% of true edges, add spurious edges equal to 20% of the rest
    keep = true.sample(frac=0.7, random_state=seed)
    n_fake = int(0.2 * len(keep))
    fake = pd.DataFrame({"source": rng.choice(tfs, n_fake), "target": rng.choice(genes, n_fake),
                         "weight": rng.choice([-1.0, 1.0], n_fake)})
    pd.concat([keep, fake]).to_csv(f"{out_dir}/synthetic_prior.csv", index=False)
    print(f"wrote {adata.n_obs} cells, {len(conditions) - 1} perturbations to {out_dir}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=".")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    main(a.out, a.seed)
