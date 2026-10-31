"""Rebuild the pseudobulk with the top variable genes PLUS every measured
CollecTRI target of the perturbed TFs (detected in >=1% of cells)."""
import sys, os, types, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
sys.path.insert(0, os.path.dirname(__file__))
from prep_large_h5ad import write
A = ad.read_h5ad('data/norman.h5ad', backed='r')
z = np.load('data/norman_pb.h5ad.stage1.npz'); totals, cols = z['totals'], z['cols']
net = pd.read_csv('data/collectri.csv')
genes = np.asarray(A.var_names); idx = {g: i for i, g in enumerate(genes)}
labels = A.obs['perturbation'].astype(str).unique()
pert = set(g for l in labels if l != 'control' for g in l.split('_'))
tg = set(net[net.source.isin(pert)].target) & set(idx)
cand = sorted(set(idx[g] for g in tg) - set(cols.tolist()))
det = np.array([A.X[:, c:c+1].nnz for c in cand]) / A.shape[0]
add = [c for c, d in zip(cand, det) if d >= 0.01]
print('candidate targets', len(cand), 'kept (>=1% cells)', len(add))
allc = np.array(sorted(set(cols.tolist()) | set(add)))
scale = sp.diags((1e4 / totals).astype(np.float32))
a = types.SimpleNamespace(pert_key='perturbation', out='data/norman_pb_targets.h5ad')
write(A, a, scale, allc)
