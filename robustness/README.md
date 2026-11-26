# Robustness checks on the Norman et al. 2019 screen

Metric: Pearson correlation on the top 20 changed genes (pearson_top20), 5-fold CV over 105 unseen single perturbations.
Numbers below are for the perturbed TFs with at least one CollecTRI target among the measured genes:
34 in the original 2,085-gene set and 41 in the 2,873-gene `targets` set. (43 perturbed genes are
CollecTRI regulators; 9 of them have no measured target in the original set. An earlier version of
this file labelled these numbers as the 43.)

| Run | Folder | linear[expression] | linear[both] (real CollecTRI) | linear[both, shuffled] |
|---|---|---|---|---|
| Original (2,085 genes, k_prior=50) | base | 0.598 | 0.550 | 0.592 |
| k_prior=20 | k20 | 0.598 | 0.573 | 0.587 |
| k_prior=100 | k100 | 0.598 | 0.555 | 0.592 |
| Unsigned network (all weights +1) | unsigned | 0.598 | 0.543 | 0.596 |
| + measured CollecTRI targets of perturbed TFs (2,873 genes; 41 TFs) | targets | 0.527 | 0.497 | 0.494 |
| Shared response removed | noshared | not interpretable (see note below) | | |

Conclusion: the real network never beat the shuffled network or the expression-only baseline.
Real > shuffled for 9 to 14 of the 34 TFs in the four original-gene-set runs, and 16 of 41 in `targets`.
Real vs expression-only was slightly worse in every run (paired Wilcoxon p = 0.0008 in base, 0.011 to
0.017 for k20 and k100, 0.001 for unsigned, 0.002 for targets).

The shuffled control relabels every gene, which leaves only 2 of the 34 TFs with any targets, so it
behaves like no network rather than a wrong one. v2 replaces it with permuted networks (see the main README).

The "targets" data: `scripts/prep_add_targets.py` adds every CollecTRI target of a perturbed TF that is detected in >= 1% of cells
(788 genes added; 41 of the 43 regulators then have at least one measured target, and the median among
those with any goes from 5 to 11).

Commands: see "Data and reproducing every result" in the main README. The inputs
for `unsigned` and `noshared` come from `scripts/make_variants.py`; the `targets`
data from `scripts/prep_add_targets.py`.

Note for noshared: the shared response was the mean effect over all single perturbations, including the
held-out ones. Inside each fold that makes the training mean anti-correlated with the test perturbations,
which pushes every model below zero by construction, so this run says nothing about the models. v2's
centered metric removes the shared response correctly, using the training perturbations only.
