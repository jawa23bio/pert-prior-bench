"""How well could ANY model do? Split-half reliability of each perturbation's effect.

A score of 0.52 means nothing without knowing the best achievable score. Each
perturbation's effect is a mean over a finite number of cells, so the "truth" a
model is scored against is itself noisy. This script splits the cells of every
perturbation (and the control cells) into two random halves, computes the
effect from each half separately, and scores one half against the other with
the benchmark's own metrics. That is what a perfect model would score against
half the data ("technical duplicate" bound).

Reported per metric:
  half_vs_half   mean agreement between the two halves
  ceiling        sqrt(2r / (1 + r)): the Spearman-Brown estimate of the best
                 score against the full-data effect (r = half_vs_half)

Read the benchmark tables against `ceiling`. If the ceiling of the centered
metric is itself low, the perturbation-specific effects are mostly noise at this
cell count and no model or prior could have been shown to predict them.

Needs the cell-level file. It is read once, in column chunks, so it fits in a few GB.

    python scripts/noise_ceiling.py --h5ad data/norman.h5ad --sep _ \
        --genes-from data/norman_pb.h5ad --stage1 data/norman_pb.h5ad.stage1.npz \
        --prior data/collectri.csv --out results_noise_ceiling
"""
import argparse
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pertbench.data import _canonical, parse_targets  # noqa: E402
from pertbench.evaluate import per_perturbation  # noqa: E402


def half_pseudobulks(A, labels, cols, totals, chunk, seeds, checkpoint=None):
    """Mean log-normalized expression per (split, label, half) over the genes in `cols`.

    The file is read ONCE, in column chunks, for all random splits together, so
    extra splits cost almost nothing. Returns (groups, M, counts) with
    M[split, group, half, gene] and counts[split, group, half].

    `checkpoint`: optional .npz path. Progress is saved there every few chunks
    and picked up on the next run, so an interrupted read does not start over.
    """
    groups, inv = np.unique(labels, return_inverse=True)
    n_rows = 2 * len(groups)
    Gs, counts = [], []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        half = np.zeros(len(labels), dtype=int)
        for k in range(len(groups)):  # alternate after shuffling: halves differ by at most one cell
            idx = rng.permutation(np.where(inv == k)[0])
            half[idx[1::2]] = 1
        gid = inv * 2 + half
        c = np.bincount(gid, minlength=n_rows).astype(np.float32)
        G = sp.csr_matrix((np.ones(len(gid), dtype=np.float32), (gid, np.arange(len(gid)))),
                          shape=(n_rows, len(gid)))
        Gs.append(sp.diags(1.0 / np.maximum(c, 1)) @ G)
        counts.append(c.reshape(len(groups), 2))
    G = sp.vstack(Gs).tocsr()
    scale = sp.diags((1e4 / totals).astype(np.float32))
    out = np.zeros((len(seeds) * n_rows, len(cols)), dtype=np.float64)
    want = {int(c): j for j, c in enumerate(cols)}
    n_chunks = -(-A.shape[1] // chunk)
    key = np.array([len(seeds), chunk, len(cols), int(np.sum(cols)), A.shape[0]])   # must match to resume
    done = 0
    if checkpoint and Path(checkpoint).exists():
        z = np.load(checkpoint)
        if np.array_equal(z["key"], key):
            out, done = z["out"], int(z["done"])
            print(f"  resuming from chunk {done + 1}/{n_chunks}", flush=True)
    for i, s in enumerate(range(0, A.shape[1], chunk)):
        if i < done:
            continue
        local = [c for c in range(s, min(s + chunk, A.shape[1])) if c in want]
        if local:
            X = sp.csc_matrix(A.X[:, s:s + chunk], dtype=np.float32)[:, [c - s for c in local]]
            X = (scale @ X).tocsc()
            X.data = np.log1p(X.data)
            out[:, [want[c] for c in local]] = (G @ X).toarray()
        print(f"  read chunk {i + 1}/{n_chunks}", flush=True)
        if checkpoint and ((i + 1) % 4 == 0 or i + 1 == n_chunks):
            tmp = str(checkpoint) + ".tmp.npz"
            np.savez(tmp, out=out, done=i + 1, key=key)
            os.replace(tmp, checkpoint)
    return groups, out.reshape(len(seeds), len(groups), 2, len(cols)), np.stack(counts)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--h5ad", required=True, help="cell-level screen")
    ap.add_argument("--pert-key", default="perturbation")
    ap.add_argument("--control", default="control")
    ap.add_argument("--sep", default="+")
    ap.add_argument("--genes-from", default="", help="pseudobulk .h5ad whose genes to use (default: all genes)")
    ap.add_argument("--stage1", default="", help="stage-1 cache from prep_large_h5ad.py (library sizes)")
    ap.add_argument("--prior", default="", help="optional source,target,weight CSV to split the summary by coverage")
    ap.add_argument("--min-cells", type=int, default=30)
    ap.add_argument("--splits", type=int, default=5, help="independent random half-splits to average")
    ap.add_argument("--chunk", type=int, default=2000)
    ap.add_argument("--checkpoint", default="", help=".npz file to save/resume the read (optional)")
    ap.add_argument("--out", default="results_noise_ceiling")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    A = ad.read_h5ad(a.h5ad, backed="r")
    genes = np.asarray(A.var_names)
    labels = A.obs[a.pert_key].astype(str).values
    if a.stage1:
        totals = np.load(a.stage1)["totals"]
    else:
        totals = np.zeros(A.shape[0])
        for s in range(0, A.shape[1], a.chunk):
            totals += np.asarray(sp.csr_matrix(A.X[:, s:s + a.chunk]).sum(axis=1)).ravel()
        totals[totals == 0] = 1.0
    if a.genes_from:
        keep = set(ad.read_h5ad(a.genes_from, backed="r").var_names)
        cols = np.array([i for i, g in enumerate(genes) if g in keep])
    else:
        cols = np.arange(len(genes))
    print(f"{A.shape[0]} cells, {len(cols)} genes", flush=True)

    rows = []
    groups, M_all, counts_all = half_pseudobulks(A, labels, cols, totals, a.chunk, seeds=range(a.splits),
                                                    checkpoint=a.checkpoint or None)
    for split in range(a.splits):
        M, counts = M_all[split], counts_all[split]
        ctrl = list(groups).index(a.control)
        singles = [k for k, g in enumerate(groups)
                   if g != a.control and len(parse_targets(g, a.sep)) == 1 and counts[k].sum() >= a.min_cells]
        names = [_canonical(groups[k], a.sep) for k in singles]
        d = [pd.DataFrame(M[singles, h] - M[ctrl, h], index=names, columns=genes[cols]) for h in (0, 1)]
        n = len(names)
        for x, y in [(0, 1), (1, 0)]:  # each half takes a turn as "truth"
            # center = mean effect of the OTHER perturbations, computed within each half
            cx = (d[x].sum(axis=0).values[None, :] - d[x].values) / (n - 1)
            cy = (d[y].sum(axis=0).values[None, :] - d[y].values) / (n - 1)
            raw = per_perturbation(d[x], d[y], None)[["perturbation", "pearson_top20", "pearson_all"]]
            ctr = per_perturbation(d[x] - cx, d[y] - cy, None)[["perturbation", "pearson_top20", "pearson_all"]]
            ctr = ctr.rename(columns={"pearson_top20": "pearson_ctr_top20", "pearson_all": "pearson_ctr_all"})
            df = raw.merge(ctr, on="perturbation")
            df["split"] = split
            df["n_cells"] = [int(counts[k].sum()) for k in singles]
            rows.append(df)
    df = pd.concat(rows)
    per = df.groupby("perturbation").mean(numeric_only=True).drop(columns="split")
    per.to_csv(out / "noise_ceiling_per_perturbation.csv")

    metrics = ["pearson_top20", "pearson_all", "pearson_ctr_top20", "pearson_ctr_all"]
    subsets = {f"all {len(per)} singles": per}
    if a.prior:
        # covered = has at least one measured target in the network (the v2 definition)
        net = pd.read_csv(a.prior)
        used = set(genes[cols])
        src = set(net[net["target"].isin(used) & (net["source"] != net["target"])]["source"])
        inp = per.index.isin(src)
        subsets[f"{int(inp.sum())} covered"] = per[inp]
        subsets[f"{int((~inp).sum())} not covered"] = per[~inp]

    def sb(r):
        return float(np.sqrt(2 * r / (1 + r))) if r > 0 else float("nan")

    tab = pd.DataFrame({name: {f"{m} {kind}": (sub[m].mean() if kind == "half_vs_half" else sb(sub[m].mean()))
                               for m in metrics for kind in ("half_vs_half", "ceiling")}
                        for name, sub in subsets.items()})
    lines = ["# Noise ceiling (split-half reliability)\n",
             f"{a.splits} random half-splits of every perturbation's cells and of the control cells; "
             f"median cells per perturbation {int(per['n_cells'].median())}.\n",
             "`half_vs_half`: agreement between the effects computed from the two halves. "
             "`ceiling`: Spearman-Brown estimate of the best achievable score against the full-data effect.\n",
             tab.round(3).to_markdown(), ""]
    if per["n_cells"].nunique() > 1:
        r = np.corrcoef(np.log(per["n_cells"]), per["pearson_ctr_top20"])[0, 1]
        lines += [f"Correlation between log(cells per perturbation) and centered reliability: {r:.2f} "
                  "(positive = the low-cell perturbations are the noisy ones).", ""]
    (out / "summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
