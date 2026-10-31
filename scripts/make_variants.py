"""Inputs for two robustness checks.

unsigned prior: every CollecTRI edge gets weight +1 (tests whether wrong signs hurt).
noshared data:  subtract the mean effect of all single perturbations from every
                perturbed row, removing the response they share (about a third of
                the variance in Norman). Diagnostic only: the mean includes the
                held-out perturbations.

    python scripts/make_variants.py --prior data/collectri.csv --h5ad data/norman_pb.h5ad
"""
import argparse

import anndata as ad
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--prior", default="data/collectri.csv")
    ap.add_argument("--h5ad", default="data/norman_pb.h5ad", help="pseudobulk file (one row per perturbation)")
    ap.add_argument("--pert-key", default="perturbation")
    ap.add_argument("--control", default="control")
    ap.add_argument("--sep", default="_")
    a = ap.parse_args()

    net = pd.read_csv(a.prior)
    net["weight"] = 1
    out = a.prior.replace(".csv", "_unsigned.csv")
    net.to_csv(out, index=False)
    print("wrote", out)

    A = ad.read_h5ad(a.h5ad)
    lab = A.obs[a.pert_key].astype(str).values
    X = np.asarray(A.X, dtype=np.float32).copy()
    ctrl = X[lab == a.control][0]
    singles = (lab != a.control) & np.array([a.sep not in x for x in lab])
    X[lab != a.control] -= X[singles].mean(0) - ctrl
    A.X = X
    out = a.h5ad.replace(".h5ad", "_noshared.h5ad")
    A.write_h5ad(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
