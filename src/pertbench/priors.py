"""Regulatory prior: turn a TF -> target network into a vector per perturbed gene.

The idea: if gene g is a transcription factor, perturbing it should move its
targets, and their targets, and so on. We propagate a signed signal a few steps
through the prior network and use the result as g's embedding.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def load_collectri(organism: str = "human") -> pd.DataFrame:
    """CollecTRI TF -> target network via decoupler (handles 1.x and 2.x APIs).

    Returns columns: source, target, weight (sign of regulation).
    """
    import decoupler as dc

    if hasattr(dc, "op") and hasattr(dc.op, "collectri"):
        net = dc.op.collectri(organism=organism)
    else:
        net = dc.get_collectri(organism=organism, split_complexes=False)
    return net[["source", "target", "weight"]]


def propagation_matrix(
    net: pd.DataFrame, genes: list[str], steps: int = 3, decay: float = 0.5
) -> pd.DataFrame:
    """Rows = all genes in `genes`; row g = signed influence of perturbing g.

    D = sum_{k=1..steps} decay^(k-1) * A^k, where A is the signed adjacency
    restricted to measured genes and scaled so no row sums above 1 in magnitude.
    Genes with no outgoing edges get an all-zero row (no prior information).
    """
    idx = {g: i for i, g in enumerate(genes)}
    n = len(genes)
    A = np.zeros((n, n))
    sub = net[net["source"].isin(idx) & net["target"].isin(idx)]
    for s, t, w in sub[["source", "target", "weight"]].itertuples(index=False):
        A[idx[s], idx[t]] += np.sign(w) if w != 0 else 1.0
    scale = np.abs(A).sum(axis=1, keepdims=True)
    scale[scale == 0] = 1.0
    A = A / scale
    D, Ak = np.zeros_like(A), np.eye(n)
    for k in range(steps):
        Ak = Ak @ A
        D += (decay**k) * Ak
    return pd.DataFrame(D, index=genes, columns=genes)


def coverage(net: pd.DataFrame, perturbed_genes: list[str]) -> float:
    """Fraction of perturbed genes that are regulators in the prior."""
    sources = set(net["source"])
    return float(np.mean([g in sources for g in perturbed_genes])) if perturbed_genes else 0.0
