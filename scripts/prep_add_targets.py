"""Rebuild the pseudobulk with the top variable genes PLUS every measured
target of the perturbed TFs in the prior network.

Why: with 2,085 variable genes, a covered TF had a median of only 4 of its
CollecTRI targets among the measured genes, which gives the prior little to
work with. This adds every target of a perturbed TF that is detected in at
least --min-frac of cells (788 genes for Norman; median rises to 11).

Needs the stage-1 cache written by prep_large_h5ad.py (library sizes + genes).

    python scripts/prep_add_targets.py --h5ad data/norman.h5ad \
        --stage1 data/norman_pb.h5ad.stage1.npz --prior data/collectri.csv \
        --out data/norman_pb_targets.h5ad
"""
import argparse
import os
import sys
import types

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(__file__))
from prep_large_h5ad import write  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--h5ad", default="data/norman.h5ad")
    ap.add_argument("--stage1", default="data/norman_pb.h5ad.stage1.npz")
    ap.add_argument("--prior", default="data/collectri.csv")
    ap.add_argument("--out", default="data/norman_pb_targets.h5ad")
    ap.add_argument("--pert-key", default="perturbation")
    ap.add_argument("--control", default="control")
    ap.add_argument("--sep", default="_")
    ap.add_argument("--min-frac", type=float, default=0.01)
    a = ap.parse_args()

    A = ad.read_h5ad(a.h5ad, backed="r")
    z = np.load(a.stage1)
    totals, cols = z["totals"], z["cols"]
    net = pd.read_csv(a.prior)
    genes = np.asarray(A.var_names)
    idx = {g: i for i, g in enumerate(genes)}
    labels = A.obs[a.pert_key].astype(str).unique()
    pert = {g for lab in labels if lab != a.control for g in lab.split(a.sep)}
    targets = set(net[net["source"].isin(pert)]["target"]) & set(idx)
    cand = sorted({idx[g] for g in targets} - set(cols.tolist()))
    det = np.array([A.X[:, c:c + 1].nnz for c in cand]) / A.shape[0]
    add = [c for c, d in zip(cand, det) if d >= a.min_frac]
    print(f"candidate targets {len(cand)}, kept (>= {a.min_frac:.0%} of cells) {len(add)}")

    allc = np.array(sorted(set(cols.tolist()) | set(add)))
    scale = sp.diags((1e4 / totals).astype(np.float32))
    write(A, types.SimpleNamespace(pert_key=a.pert_key, out=a.out), scale, allc)


if __name__ == "__main__":
    main()
