"""Run with:  pytest -q     (about a minute; no data download needed)

Each test pins one property the conclusions depend on. If you change the
models or the protocol and one of these fails, the result tables no longer mean
what the README says they mean.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pertbench import protocol as P
from pertbench.calibrate import Z80, conformal_scale, interval_report, oof_residuals
from pertbench.data import kfold
from pertbench.evaluate import per_perturbation
from pertbench.models import LinearModel, MeanBaseline
from pertbench.priors import (permute_regulators, permute_targets, prior_svd, propagation_matrix,
                              regulators, relabel_all, transferability)

ROOT = Path(__file__).resolve().parents[1]


# ------------------------------------------------------------------ fixtures
def toy(n_genes=400, n_tf=40, n_other=40, targets=6, beta=1.0, seed=0):
    """Effects = shared response + beta * (network row) + noise. Sparse network:
    each TF has `targets` targets among non-TF genes, so target sets barely overlap."""
    rng = np.random.default_rng(seed)
    genes = [f"G{i:03d}" for i in range(n_genes)]
    tfs, others = genes[:n_tf], genes[n_tf:n_tf + n_other]
    net = pd.DataFrame([(tf, t, rng.choice([-1.0, 1.0])) for tf in tfs
                        for t in rng.choice(genes[n_tf + n_other:], targets, replace=False)],
                       columns=["source", "target", "weight"])
    D = propagation_matrix(net, genes)
    shared = rng.normal(0, 0.3, n_genes)
    Y = {p: shared * rng.uniform(0.7, 1.3) + beta * D.loc[p].values + rng.normal(0, 0.05, n_genes)
         for p in tfs + others}
    return pd.DataFrame(Y, index=genes).T, D, tfs, others


# ------------------------------------------------- 1. the original still runs
def test_original_pipeline_reproduces_committed_synthetic_results(tmp_path):
    run = lambda *a: subprocess.run([sys.executable, "-W", "ignore", *a], check=True, capture_output=True)
    run(str(ROOT / "scripts/make_synthetic.py"), "--out", str(tmp_path))
    run(str(ROOT / "scripts/run_benchmark.py"), "--h5ad", str(tmp_path / "synthetic.h5ad"),
        "--prior", str(tmp_path / "synthetic_prior.csv"), "--n-hvg", "300", "--out", str(tmp_path / "out"))
    for f in ["singles_per_perturbation.csv", "doubles_per_perturbation.csv"]:
        new, old = pd.read_csv(tmp_path / "out" / f), pd.read_csv(ROOT / "results_synthetic" / f)
        # float32 pseudobulk: values agree to ~1e-8 across numpy versions, not bit for bit
        pd.testing.assert_frame_equal(new, old, rtol=1e-4, atol=1e-5)


# ------------------------------------------------------------ 2. the controls
def test_v1_control_gives_most_perturbed_tfs_no_targets():
    """Documents WHY the original shuffled control is weak: after relabelling
    every gene, a perturbed TF almost always lands on a non-regulator."""
    _, D, tfs, _ = toy()
    left = [len(set(tfs) & set(regulators(relabel_all(D, s)))) for s in range(20)]
    assert np.mean(left) < 0.25 * len(tfs)


def test_regulator_permutation_keeps_every_tf_a_regulator_with_wrong_targets():
    _, D, tfs, _ = toy()
    S = permute_regulators(D, seed=3, pool=tfs)
    assert set(regulators(S)) == set(regulators(D))
    assert not any(np.allclose(S.loc[g], D.loc[g]) for g in tfs)          # nobody keeps their own row
    assert np.allclose(np.sort(np.abs(S.values).sum(1)), np.sort(np.abs(D.values).sum(1)))


def test_target_permutation_keeps_tf_similarity():
    _, D, _, _ = toy()
    T = permute_targets(D, seed=3)
    assert np.allclose(T.values @ T.values.T, D.values @ D.values.T)       # same overlaps between TFs
    assert not np.allclose(T.values, D.values)


def test_fast_prior_svd_matches_full_svd():
    _, D, _, _ = toy()
    Rp, Rt = prior_svd(D, 30)
    U, S, Vt = np.linalg.svd(D.values, full_matrices=False)
    assert np.allclose(Rp.values @ Rt.values.T, (U[:, :30] * S[:30]) @ Vt[:30], atol=1e-8)


# -------------------------------------------------------------- 3. the metric
def test_centered_metric_scores_the_mean_baseline_zero_and_a_perfect_model_one():
    Y, _, tfs, others = toy()
    tr, te = (tfs + others)[:60], (tfs + others)[60:]
    center = Y.loc[tr].mean(axis=0)
    mean_pred, _ = MeanBaseline().fit(Y.loc[tr]).predict(te)
    m = per_perturbation(Y.loc[te], mean_pred, None, center=center)
    assert (m["pearson_ctr_top20"] == 0).all() and m["pearson_top20"].mean() > 0.5   # v1 metric rewards the mean
    perfect = per_perturbation(Y.loc[te], Y.loc[te], None, center=center)
    assert np.allclose(perfect["pearson_ctr_top20"], 1)


def test_centering_on_a_mean_that_includes_the_test_set_biases_scores_negative():
    """The 'shared response removed' robustness check subtracted the mean of ALL
    singles. This shows the artefact on data with no specific signal at all."""
    rng = np.random.default_rng(0)
    Y = pd.DataFrame(rng.normal(0, 1, (50, 300)) + rng.normal(0, 1, 300), index=[f"P{i}" for i in range(50)],
                     columns=[f"G{i}" for i in range(300)])
    Yg = Y - Y.mean(axis=0)                                # global centering (leaky)
    leaky, clean = [], []
    for tr, te in kfold(list(Y.index), k=5):
        pred, _ = MeanBaseline().fit(Yg.loc[tr]).predict(te)
        leaky.append(per_perturbation(Yg.loc[te], pred, None)["pearson_all"].mean())
        expr = LinearModel("expression", lam=1.0, n_boot=0).fit(Y.loc[tr])
        mu, _ = expr.predict(te)
        clean.append(per_perturbation(Y.loc[te], mu, None, center=Y.loc[tr].mean(axis=0))["pearson_ctr_all"].mean())
    assert np.mean(leaky) < -0.1          # looks anti-predictive, purely by construction
    assert abs(np.mean(clean)) < 0.05     # train-only centering: correctly ~0


# ------------------------------------------------------------- 4. the models
def test_direct_model_recovers_a_planted_network_effect_and_nulls_do_not():
    Y, D, tfs, others = toy(beta=1.0)
    cfg = P.Config(repeats=1, n_null=5, with_linear_prior=False)
    df = P.run_singles(Y, tfs + others, D, cfg, P.make_nulls(D, cfg, pool=tfs), log=lambda *_: None)
    cov = df.drop_duplicates("perturbation").set_index("perturbation")["covered"]
    t = P.null_test(df, "pearson_ctr_top20", "direct[prior]", cov).set_index("null")
    assert t.loc["regulator-permuted", "real"] > 0.4
    assert abs(t.loc["regulator-permuted", "null_mean"]) < 0.1
    assert t["p"].max() <= 1 / 6 + 1e-9                    # smallest possible p with 5 draws


def test_no_planted_effect_means_no_detection():
    Y, D, tfs, others = toy(beta=0.0)
    cfg = P.Config(repeats=1, n_null=9, with_linear_prior=False)
    df = P.run_singles(Y, tfs + others, D, cfg, P.make_nulls(D, cfg, pool=tfs), log=lambda *_: None)
    cov = df.drop_duplicates("perturbation").set_index("perturbation")["covered"]
    t = P.null_test(df, "pearson_ctr_top20", "direct[prior]", cov)
    assert (t["real"].abs() < 0.1).all()


def test_block_scaling_cannot_blow_up_but_max_scaling_can():
    """A prior feature that is tiny for every training TF and large for a held-out
    one: 'max' scaling divides by the tiny number."""
    Y, D, tfs, others = toy()
    Rp, Rt = prior_svd(D, 30)
    Rp = Rp.copy()
    tr, te = tfs[:30] + others, tfs[30:]
    Rp.loc[tr, 0] *= 1e-6                                 # feature 0 nearly absent from training
    worst = {}
    for scaling in ["max", "block"]:
        m = LinearModel("both", k_prior=30, lam=0.1, n_boot=0, prior_embedding=(Rp, Rt), scaling=scaling)
        mu, _ = m.fit(Y.loc[tr]).predict(te)
        worst[scaling] = float(np.abs(mu.values).max())
    assert worst["block"] < 10 * np.abs(Y.values).max()
    assert worst["max"] > 100 * worst["block"]


def test_stratified_folds_are_a_partition_and_balanced():
    items = [f"P{i}" for i in range(105)]
    strata = [i < 43 for i in range(105)]
    folds = P.stratified_kfold(items, strata, 5, seed=1)
    tests = [x for _, te in folds for x in te]
    assert sorted(tests) == sorted(items)
    assert all(not (set(tr) & set(te)) for tr, te in folds)
    n_cov = [sum(int(x[1:]) < 43 for x in te) for _, te in folds]
    assert max(n_cov) - min(n_cov) <= 1


def test_transferability_is_zero_for_disjoint_targets_and_one_for_shared():
    E = pd.DataFrame(np.eye(4), index=list("ABCD"))
    E.loc["D"] = E.loc["A"] + E.loc["B"]
    assert transferability(E, ["A", "B"], ["C"]).loc["C"] == pytest.approx(0)
    assert transferability(E, ["A", "B"], ["D"]).loc["D"] == pytest.approx(1)


# -------------------------------------------------------- 5. the uncertainty
def test_conformal_intervals_cover_and_bootstrap_intervals_do_not():
    Y, _, tfs, others = toy(beta=0.0)
    perts = tfs + others
    tr, te = perts[:60], perts[60:]
    make = lambda: LinearModel("expression", lam=1.0, n_boot=20)
    mu, sd = make().fit(Y.loc[tr]).predict(te)
    sigma, q = conformal_scale(oof_residuals(make, Y.loc[tr], k=5))
    boot = interval_report(Y.loc[te], mu, Z80 * sd)["coverage80"].mean()
    conf = interval_report(Y.loc[te], mu, q * sigma)["coverage80"].mean()
    assert 0.74 < conf < 0.86
    assert boot < 0.5
