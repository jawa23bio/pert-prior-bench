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


# --------------------------------------------------------------------------
# Null networks
#
# A control network should be wrong about WHICH genes are linked while keeping
# everything else. Relabelling all genes at once (the original `shuffled` in
# run_benchmark.py) does not do that for the perturbed genes: only a few percent
# of measured genes are regulators, so after a full relabel almost every perturbed
# TF lands on a non-regulator and receives an all-zero row. That control compares
# "real network" with "no network", not with "wrong network". The two nulls below
# keep every regulator a regulator.
# --------------------------------------------------------------------------
def regulators(D: pd.DataFrame) -> list[str]:
    """Genes with at least one outgoing (propagated) edge among the measured genes."""
    return list(D.index[(D.abs().sum(axis=1) > 0).values])


def permute_regulators(D: pd.DataFrame, seed: int, pool: list[str] | None = None) -> pd.DataFrame:
    """Swap target sets between regulators: each regulator receives another
    regulator's row. Columns are untouched, so the same genes are targets and
    every row keeps its size. Destroys only 'this TF regulates these genes'.

    `pool` restricts the swap to a subset (e.g. the perturbed regulators), which
    keeps the null from importing the much larger target sets of unperturbed hubs.
    No regulator keeps its own row.
    """
    regs = [g for g in (pool if pool is not None else regulators(D)) if g in D.index]
    if len(regs) < 2:
        return D.copy()
    rng = np.random.default_rng(seed)
    while True:  # derangement: nobody keeps their own targets
        perm = rng.permutation(len(regs))
        if not np.any(perm == np.arange(len(regs))):
            break
    S = D.copy()
    S.loc[regs] = D.loc[[regs[i] for i in perm]].values
    return S


def permute_targets(D: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Relabel the target side only: every regulator keeps its number of targets,
    its weights and its overlap with every other regulator, but the targets are
    now random genes. Destroys 'the effect lands on these genes' while leaving
    the similarity between regulators intact."""
    rng = np.random.default_rng(seed)
    S = D.copy()
    S.columns = rng.permutation(D.columns)
    return S.loc[:, D.columns]


def relabel_all(D: pd.DataFrame, seed: int) -> pd.DataFrame:
    """The original control (kept for comparison): one permutation applied to
    rows and columns. See the note above for why it is a weak null."""
    perm = np.random.default_rng(seed).permutation(D.index)
    S = D.copy()
    S.index, S.columns = perm, perm
    return S.loc[D.index, D.columns]


def prior_svd(D: pd.DataFrame, k: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Same decomposition LinearModel uses (D = U S V^T; regulator side U S, target
    side V), computed on the non-zero rows only. All-zero rows contribute nothing
    to an SVD, so the result is identical and much faster (a few hundred
    regulators instead of every measured gene)."""
    genes = list(D.index)
    nz = (D.abs().sum(axis=1) > 0).values
    U, S, Vt = np.linalg.svd(D.values[nz], full_matrices=False)
    k = min(k, len(S))
    Rp = np.zeros((len(genes), k))
    Rp[nz] = U[:, :k] * S[:k]
    return pd.DataFrame(Rp, index=genes), pd.DataFrame(Vt[:k].T, index=D.columns)


def retained(D: pd.DataFrame, k: int, genes: list[str]) -> pd.Series:
    """Share of each regulator's prior row (squared norm) that survives the rank-k
    SVD used by LinearModel. The SVD is dominated by hub TFs with many targets; a
    perturbed TF with a handful of measured targets can be almost entirely
    truncated away, in which case the model never sees its targets."""
    Rp, _ = prior_svd(D, k)
    full = (D.values ** 2).sum(axis=1)
    kept = (Rp.values ** 2).sum(axis=1)
    out = pd.Series(np.divide(kept, full, out=np.zeros_like(kept), where=full > 0), index=D.index)
    return out.loc[[g for g in genes if g in D.index and full[D.index.get_loc(g)] > 0]]


def transferability(E: pd.DataFrame, train_regs: list[str], test_regs: list[str]) -> pd.Series:
    """For each held-out regulator: the share of its embedding (squared norm) that
    lies in the span of the TRAINING regulators' embeddings. `E` is any
    regulator x feature matrix: the propagation matrix itself, or the regulator
    side of `prior_svd` (what LinearModel actually uses).

    A model that is linear in the embedding can only predict a held-out regulator
    through directions the training regulators also use. Near 0 means the model
    has nothing to transfer, however correct the network is. Check this before
    reading a null result as evidence against the network.
    """
    tr = [g for g in train_regs if g in E.index]
    A = E.loc[tr].values
    A = A[np.abs(A).sum(axis=1) > 0]
    Q = None
    if len(A):
        _, sv, Vt = np.linalg.svd(A, full_matrices=False)
        Q = Vt[sv > 1e-10 * sv.max()]
    out = {}
    for g in test_regs:
        if g not in E.index:
            continue
        d = E.loc[g].values
        n = float(d @ d)
        if n == 0:
            continue
        out[g] = float(np.sum((Q @ d) ** 2) / n) if Q is not None else 0.0
    return pd.Series(out, dtype=float)
