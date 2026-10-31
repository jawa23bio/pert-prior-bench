"""Shrink a large screen so it fits in limited memory (e.g. 3 GB).

Reads the .h5ad in column chunks (backed mode), log-normalizes with each cell's
FULL library size, keeps the top-N variable genes plus every perturbed gene, and
writes a small .h5ad of log-normalized values. Run the benchmark on the output
with --already-log.

    python scripts/prep_large_h5ad.py --h5ad data/norman.h5ad --out data/norman_small.h5ad
"""
import argparse
import numpy as np
import scipy.sparse as sp
import anndata as ad
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from pertbench.data import parse_targets


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pert-key", default="perturbation")
    ap.add_argument("--control", default="control")
    ap.add_argument("--sep", default="_")
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--chunk", type=int, default=2000)
    a = ap.parse_args()

    A = ad.read_h5ad(a.h5ad, backed="r")
    n, g = A.shape
    cache = a.out + ".stage1.npz"
    if os.path.exists(cache):
        z = np.load(cache)
        totals, cols = z["totals"], z["cols"]
        scale = sp.diags((1e4 / totals).astype(np.float32))
        return write(A, a, scale, cols)
    totals = np.zeros(n, dtype=np.float64)
    for s in range(0, g, a.chunk):
        totals += np.asarray(A.X[:, s:s + a.chunk].sum(axis=1)).ravel()
    totals[totals == 0] = 1.0
    scale = sp.diags((1e4 / totals).astype(np.float32))
    print("library sizes done", flush=True)

    mean = np.zeros(g); sq = np.zeros(g)
    for s in range(0, g, a.chunk):
        X = scale @ sp.csr_matrix(A.X[:, s:s + a.chunk], dtype=np.float32)
        X.data = np.log1p(X.data)
        mean[s:s + a.chunk] = np.asarray(X.mean(axis=0)).ravel()
        sq[s:s + a.chunk] = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    var = sq - mean**2
    keep = set(np.argsort(var)[::-1][:a.n_hvg].tolist())
    genes = np.asarray(A.var_names)
    idx = {x: i for i, x in enumerate(genes)}
    labels = A.obs[a.pert_key].astype(str).values
    missing = set()
    for lab in np.unique(labels[labels != a.control]):
        for t in parse_targets(lab, a.sep):
            if t in idx: keep.add(idx[t])
            else: missing.add(t)
    cols = np.array(sorted(keep))
    print(f"keeping {len(cols)} genes; perturbed genes not measured: {sorted(missing)}", flush=True)
    np.savez(cache, totals=totals, cols=cols)
    print("stage 1 saved; run the same command again to write the small file", flush=True)


def write(A, a, scale, cols, min_cells=30):
    """Write one row per perturbation: mean log-normalized expression of its cells.

    Averaging here (instead of keeping every cell) keeps memory small. The
    benchmark then subtracts the control row, which gives the same pseudobulk
    effects it would compute from the cells. Run it with --already-log --min-cells 1.
    """
    import pandas as pd
    labels = A.obs[a.pert_key].astype(str).values
    groups, inv = np.unique(labels, return_inverse=True)
    counts = np.bincount(inv)
    keep_g = counts >= min_cells
    G = sp.csr_matrix((np.ones(len(inv), dtype=np.float32), (inv, np.arange(len(inv)))),
                      shape=(len(groups), len(inv)))
    G = sp.diags(1.0 / counts.astype(np.float32)) @ G
    out = np.zeros((len(groups), len(cols)), dtype=np.float32)
    for k, c in enumerate(cols):
        x = scale @ sp.csc_matrix(A.X[:, int(c):int(c) + 1], dtype=np.float32)
        x.data = np.log1p(x.data)
        out[:, k] = np.asarray((G @ x).todense()).ravel()
    obs = pd.DataFrame({a.pert_key: groups[keep_g], "n_cells": counts[keep_g]},
                       index=[f"pb{i}" for i in range(keep_g.sum())])
    res = ad.AnnData(X=out[keep_g], obs=obs, var=A.var.iloc[cols].copy())
    res.write_h5ad(a.out)
    print("wrote", a.out, res.shape, "dropped (<%d cells):" % min_cells, list(groups[~keep_g]))


if __name__ == "__main__":
    main()
