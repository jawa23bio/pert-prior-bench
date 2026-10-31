"""Turn a single-cell perturbation screen into pseudobulk effect vectors.

Each perturbation p becomes one vector: mean log-expression in cells carrying p
minus mean log-expression in control cells, over a fixed gene set. Every model in
this project predicts these vectors for perturbations it has not seen.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def parse_targets(name: str, sep: str = "+") -> list[str]:
    """'KLF1+CEBPA' -> ['KLF1', 'CEBPA']. Drops placeholders like 'ctrl'."""
    parts = [p.strip() for p in str(name).replace("_", sep).split(sep)]
    return [p for p in parts if p and p.lower() not in {"ctrl", "control", "nan"}]


def pseudobulk_deltas(
    adata,
    pert_key: str = "perturbation",
    control_label: str = "control",
    n_hvg: int = 2000,
    min_cells: int = 30,
    sep: str = "+",
    already_log: bool = False,
) -> tuple[pd.DataFrame, pd.Series]:
    """Return (deltas: perturbations x genes, n_cells per perturbation).

    Genes = top `n_hvg` most variable genes, plus every perturbed gene that was
    measured (the linear model needs a row for the perturbed gene itself).
    """
    import scipy.sparse as sp

    # stay sparse until the gene set is chosen (a full screen does not fit densely)
    X = sp.csr_matrix(adata.X, dtype=np.float32)
    if not already_log:
        totals = np.asarray(X.sum(axis=1)).ravel()
        totals[totals == 0] = 1.0
        X = sp.diags(1e4 / totals) @ X
        X.data = np.log1p(X.data)
    X = X.tocsr()

    genes = np.asarray(adata.var_names)
    labels = adata.obs[pert_key].astype(str).values
    is_ctrl = labels == control_label
    if is_ctrl.sum() < min_cells:
        raise ValueError(f"Only {is_ctrl.sum()} control cells with label '{control_label}'.")

    mean = np.asarray(X.mean(axis=0)).ravel()
    var = np.asarray(X.multiply(X).mean(axis=0)).ravel() - mean**2
    keep = set(np.argsort(var)[::-1][:n_hvg])
    gene_index = {g: i for i, g in enumerate(genes)}
    for lab in np.unique(labels[~is_ctrl]):
        for t in parse_targets(lab, sep):
            if t in gene_index:
                keep.add(gene_index[t])
    cols = np.array(sorted(keep))
    X = X[:, cols].toarray()

    ctrl_mean = X[is_ctrl].mean(axis=0)
    rows, counts, names = [], [], []
    for lab in np.unique(labels[~is_ctrl]):
        mask = labels == lab
        if mask.sum() < min_cells:
            continue
        rows.append(X[mask].mean(axis=0) - ctrl_mean)
        counts.append(int(mask.sum()))
        names.append(_canonical(lab, sep))
    deltas = pd.DataFrame(np.vstack(rows), index=names, columns=genes[cols])
    deltas = deltas.groupby(level=0).mean()  # 'A+B' and 'B+A' collapse together
    n_cells = pd.Series(counts, index=names).groupby(level=0).sum()
    return deltas, n_cells.loc[deltas.index]


def _canonical(name: str, sep: str) -> str:
    return "+".join(sorted(parse_targets(name, sep)))


def split_singles_doubles(deltas: pd.DataFrame) -> tuple[list[str], list[str]]:
    singles = [p for p in deltas.index if len(parse_targets(p)) == 1]
    doubles = [p for p in deltas.index if len(parse_targets(p)) == 2]
    return singles, doubles


def kfold(items: list[str], k: int = 5, seed: int = 0) -> list[tuple[list[str], list[str]]]:
    rng = np.random.default_rng(seed)
    order = list(rng.permutation(items))
    folds = [order[i::k] for i in range(k)]
    return [([x for j, f in enumerate(folds) if j != i for x in f], folds[i]) for i in range(k)]
