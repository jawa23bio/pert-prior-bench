"""Models that predict a perturbation's effect vector.

All models share one interface:
    fit(train: DataFrame perts x genes) -> self
    predict(perts: list[str]) -> (mean: DataFrame, std: DataFrame or None)

The two linear models follow the structure used in Ahlmann-Eltze et al.:
    Y (genes x perts) ~ G @ W @ P.T + b
G is a gene embedding from PCA of the training effects; P embeds each
perturbation through the gene it targets. The only thing that changes between
LinearModel variants is WHERE P comes from:
    "expression": the perturbed gene's row of G      (data only)
    "prior":      the perturbed gene's row of an SVD of the propagation matrix
    "both":       the two concatenated
That makes the comparison a clean ablation: does knowledge from the regulatory
network add anything beyond what the expression data already says?
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .data import parse_targets


class MeanBaseline:
    """Predict the average training effect for every perturbation."""

    name = "mean"

    def fit(self, train: pd.DataFrame):
        self.mu = train.mean(axis=0)
        self.sd = train.std(axis=0).replace(0, 1e-6)
        return self

    def predict(self, perts):
        mean = pd.DataFrame([self.mu.values] * len(perts), index=perts, columns=self.mu.index)
        std = pd.DataFrame([self.sd.values] * len(perts), index=perts, columns=self.sd.index)
        return mean, std


class AdditiveBaseline:
    """For a double A+B, predict effect(A) + effect(B) from training singles."""

    name = "additive"

    def fit(self, train: pd.DataFrame):
        self.singles = {p: train.loc[p] for p in train.index if len(parse_targets(p)) == 1}
        self.fallback = MeanBaseline().fit(train)
        return self

    def predict(self, perts):
        rows = []
        for p in perts:
            parts = parse_targets(p)
            if all(t in self.singles for t in parts):
                rows.append(sum(self.singles[t] for t in parts))
            else:
                rows.append(self.fallback.mu)
        return pd.DataFrame(rows, index=perts), None


def _ridge_solve(M: np.ndarray, Y: np.ndarray, lam: float) -> np.ndarray:
    return np.linalg.solve(M.T @ M + lam * np.eye(M.shape[1]), M.T @ Y)


class LinearModel:
    """Ahlmann-Eltze-style bilinear model with a choice of perturbation embedding.

    Uncertainty comes from a bootstrap ensemble: refit on resampled training
    perturbations and report the spread of predictions.
    """

    def __init__(self, embedding="expression", k=10, k_prior=50, lam=1.0,
                 prior_matrix: pd.DataFrame | None = None, n_boot=20, seed=0,
                 lam_grid: list[float] | None = None,
                 scaling: str = "max", lam_out: float | None = None,
                 select: str = "pearson_top20", prior_embedding=None):
        """The last four arguments default to the original behaviour.

        scaling  "max": divide each perturbation feature by its largest training
                 value (original). A feature that is tiny in training but not in a
                 held-out perturbation is then multiplied by a huge factor, which
                 produced MSEs in the tens to hundreds in two robustness runs.
                 "block": divide the expression block and the prior block each by
                 one number (their RMS), so no single feature can blow up.
        lam_out  ridge penalty for the projection onto the output basis. None
                 reuses `lam` (original), which shrinks every prediction toward
                 the mean by 1/(1+lam) when the basis is orthonormal: harmless for
                 Pearson, but it distorts MSE and interval coverage. Use ~1e-6.
        select   what the inner CV maximises: "pearson_top20" (original) or "mse".
        prior_embedding  precomputed (regulator side, target side) from
                 priors.prior_svd, to avoid redoing the SVD in every fold.
        """
        assert embedding in {"expression", "prior", "both"}
        assert scaling in {"max", "block"} and select in {"pearson_top20", "mse"}
        if embedding != "expression" and prior_matrix is None and prior_embedding is None:
            raise ValueError("prior_matrix is required for 'prior' or 'both'.")
        self.embedding, self.k, self.k_prior, self.lam = embedding, k, k_prior, lam
        self.lam_grid = lam_grid
        self.prior_matrix, self.n_boot, self.seed = prior_matrix, n_boot, seed
        self.scaling, self.lam_out, self.select = scaling, lam_out, select
        if prior_embedding is not None:
            self._prior_cache = prior_embedding
        self.name = f"linear[{embedding}]"

    # -- embeddings ---------------------------------------------------------
    def _gene_embedding(self, train: pd.DataFrame) -> pd.DataFrame:
        Yc = (train - train.mean(axis=0)).T.values  # genes x perts
        U, S, _ = np.linalg.svd(Yc, full_matrices=False)
        k = min(self.k, U.shape[1])
        return pd.DataFrame(U[:, :k] * S[:k], index=train.columns)

    def _prior_embedding(self, genes):
        """SVD of the propagation matrix D (regulator x target) = U S V^T.

        Returns (regulator side, target side):
          Rp: row g = U S of regulator g  -> describes what perturbing g does
          Rt: row j = V of target gene j  -> where in gene space that effect lands
        so that D[g, :] = Rp[g] @ Rt.T exactly when k is large enough.
        """
        D = self.prior_matrix.reindex(index=genes, columns=genes).fillna(0.0).values
        U, S, Vt = np.linalg.svd(D, full_matrices=False)
        k = min(self.k_prior, U.shape[1])
        return (pd.DataFrame(U[:, :k] * S[:k], index=genes),
                pd.DataFrame(Vt[:k].T, index=genes))

    def _pert_features(self, perts, G, R) -> np.ndarray:
        use_g = self.embedding in {"expression", "both"}
        use_p = self.embedding in {"prior", "both"}
        Gv, gi = G.values, {g: i for i, g in enumerate(G.index)}
        if use_p:
            Rp = R[0]  # regulator side of the prior
            Rv, ri = Rp.values, {g: i for i, g in enumerate(Rp.index)}
        rows = []
        for p in perts:
            vecs = []
            for t in parse_targets(p):
                parts = []
                if use_g:
                    parts.append(Gv[gi[t]] if t in gi else np.zeros(Gv.shape[1]))
                if use_p:
                    parts.append(Rv[ri[t]] if t in ri else np.zeros(Rv.shape[1]))
                vecs.append(np.concatenate(parts))
            rows.append(np.sum(vecs, axis=0))  # doubles: sum of single embeddings
        F = np.vstack(rows)
        return F

    # -- fit / predict -------------------------------------------------------
    def _output_basis(self, G, R) -> np.ndarray:
        """Gene-side basis the prediction lives in. Columns unit-normalized so the
        expression and prior parts compete on equal footing."""
        blocks = []
        if self.embedding in {"expression", "both"}:
            blocks.append(G.values)
        if self.embedding in {"prior", "both"}:
            blocks.append(R[1].values)  # target side of the prior
        B = np.hstack(blocks)
        n = np.linalg.norm(B, axis=0)
        n[n == 0] = 1.0
        return B / n

    def _fit_once(self, train: pd.DataFrame):
        G = self._gene_embedding(train)
        R = self._prior_cache if self.embedding != "expression" else None
        b = train.mean(axis=0).values
        P = self._pert_features(list(train.index), G, R)
        if self.scaling == "block":
            scale = self._block_scale(P, G.shape[1])
        else:
            scale = np.abs(P).max(axis=0)
            # A feature that is ~0 for every training perturbation carries no usable
            # information; dividing by its tiny max would blow up test predictions.
            # Switch it off (scale = inf -> feature 0) instead.
            scale[scale <= 1e-8 * max(scale.max(), 1e-12)] = np.inf
        P = P / scale
        Y = train.values - b  # perts x genes
        B = self._output_basis(G, R)
        # Y.T ~ B W P.T  ->  solve W in two ridge steps
        lam_out = self.lam if self.lam_out is None else self.lam_out
        A = _ridge_solve(B, Y.T, lam_out)              # basis x perts  (B A ~ Y.T)
        W = _ridge_solve(P, A.T, self.lam).T            # basis x feats  (P W.T ~ A.T)
        return dict(G=G, R=R, B=B, b=b, W=W, scale=scale)

    def _block_scale(self, P: np.ndarray, k_expr: int) -> np.ndarray:
        """One scale per block (RMS over training perturbations and features)."""
        scale = np.ones(P.shape[1])
        blocks = {"expression": [(0, P.shape[1])], "prior": [(0, P.shape[1])],
                  "both": [(0, k_expr), (k_expr, P.shape[1])]}[self.embedding]
        for a, z in blocks:
            rms = np.sqrt(np.mean(P[:, a:z] ** 2)) if z > a else 0.0
            scale[a:z] = rms if rms > 0 else 1.0
        return scale

    def fit(self, train: pd.DataFrame):
        self.columns = train.columns
        # the prior's decomposition depends only on the gene list, so compute it once
        if self.embedding != "expression":
            cache = getattr(self, "_prior_cache", None)
            if cache is None or not cache[0].index.equals(train.columns):
                self._prior_cache = self._prior_embedding(list(train.columns))
        if self.lam_grid:
            self.lam = self._select_lam(train)
        self.full = self._fit_once(train)
        rng = np.random.default_rng(self.seed)
        self.boots = []
        for _ in range(self.n_boot):
            pick = rng.choice(len(train), size=len(train), replace=True)
            boot = train.iloc[pick]
            boot = boot[~boot.index.duplicated()]
            self.boots.append(self._fit_once(boot))
        return self

    def _select_lam(self, train: pd.DataFrame, k: int = 3) -> float:
        """Pick the ridge penalty by inner cross-validation on the TRAINING set only,
        scoring Pearson over each held-out perturbation's top-20 moved genes."""
        from .data import kfold
        from .evaluate import per_perturbation

        scores = {}
        for lam in self.lam_grid:
            self.lam, vals = lam, []
            for tr, te in kfold(list(train.index), k=k, seed=self.seed + 1):
                m = self._fit_once(train.loc[tr])
                pred = pd.DataFrame(self._predict_with(m, te), index=te, columns=train.columns)
                if self.select == "mse":
                    vals.append(-float(np.mean((train.loc[te].values - pred.values) ** 2)))
                else:
                    vals.append(per_perturbation(train.loc[te], pred, None)["pearson_top20"].mean())
            scores[lam] = np.nanmean(vals)
        self.lam_scores = scores
        return max(scores, key=scores.get)

    def _predict_with(self, m, perts):
        P = self._pert_features(perts, m["G"], m["R"]) / m["scale"]
        return (m["B"] @ m["W"] @ P.T).T + m["b"]

    def predict(self, perts):
        mean = pd.DataFrame(self._predict_with(self.full, perts), index=perts, columns=self.columns)
        if not self.boots:
            return mean, None
        stack = np.stack([self._predict_with(m, perts) for m in self.boots])
        std = pd.DataFrame(stack.std(axis=0) + 1e-6, index=perts, columns=self.columns)
        return mean, std


class DirectPriorModel:
    """Use the network the way it is usually read: perturbing TF g moves g's targets.

        effect(g) = mean training effect + alpha * D[g, :]

    D[g, :] is g's row of the propagation matrix (its signed, propagated target
    set) and alpha is ONE number, fitted on the training perturbations the network
    covers. Perturbations the network does not cover get the mean. Doubles use the
    sum of the two rows.

    Why it exists: LinearModel sees the network only through a truncated SVD and
    learns a k x k map from a few dozen covered TFs, so it can fail to use a
    network that is exactly right (see priors.retained / priors.transferability).
    This model has nothing to learn except a scale, so if it cannot beat the mean
    and a regulator-permuted network cannot be told apart from the real one, the
    network's edges themselves are not predictive in this screen.

    alpha is clipped at 0 by default: a negative alpha would mean "targets move
    against their annotated sign", which is a different claim from the network
    being right. Set allow_negative=True to let the data choose.
    """

    def __init__(self, prior_matrix: pd.DataFrame, n_boot=20, seed=0, allow_negative=False):
        self.D, self.n_boot, self.seed, self.allow_negative = prior_matrix, n_boot, seed, allow_negative
        self.name = "direct[prior]"

    def _rows(self, perts, cols) -> np.ndarray:
        if getattr(self, "_cols", None) is None or not self._cols.equals(cols):
            self._cols = cols
            self._Dv = self.D.reindex(columns=cols).fillna(0.0).values
            self._di = {g: i for i, g in enumerate(self.D.index)}
        out = np.zeros((len(perts), len(cols)))
        for i, p in enumerate(perts):
            for t in parse_targets(p):
                if t in self._di:
                    out[i] += self._Dv[self._di[t]]
        return out

    def _alpha(self, Yc: np.ndarray, F: np.ndarray) -> float:
        den = float((F * F).sum())
        a = float((Yc * F).sum()) / den if den > 0 else 0.0
        return a if self.allow_negative else max(a, 0.0)

    def fit(self, train: pd.DataFrame):
        self.columns = train.columns
        self.mu = train.mean(axis=0).values
        F = self._rows(list(train.index), train.columns)
        # the perturbed gene's own column is excluded from scoring, so exclude it here too
        for i, p in enumerate(train.index):
            for t in parse_targets(p):
                if t in train.columns:
                    F[i, train.columns.get_loc(t)] = 0.0
        Yc = train.values - self.mu
        self.alpha = self._alpha(Yc, F)
        rng = np.random.default_rng(self.seed)
        self.alphas = []
        for _ in range(self.n_boot):
            pick = rng.choice(len(train), size=len(train), replace=True)
            self.alphas.append(self._alpha(Yc[pick], F[pick]))
        return self

    def predict(self, perts):
        F = self._rows(perts, self.columns)
        mean = pd.DataFrame(self.mu + self.alpha * F, index=perts, columns=self.columns)
        if not self.alphas:
            return mean, None
        std = pd.DataFrame(np.std(self.alphas) * np.abs(F) + 1e-6, index=perts, columns=self.columns)
        return mean, std
