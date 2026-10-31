# Robustness checks on the Norman et al. 2019 screen

Metric: Pearson correlation on the top 20 changed genes (pearson_top20), 5-fold CV over 105 unseen single perturbations.
Numbers below are for the 43 perturbed genes that are TFs in CollecTRI.

| Run | Folder | linear[expression] | linear[both] (real CollecTRI) | linear[both, shuffled] |
|---|---|---|---|---|
| Original (2,085 genes, k_prior=50) | base | 0.598 | 0.550 | 0.592 |
| k_prior=20 | k20 | 0.598 | 0.573 | 0.587 |
| k_prior=100 | k100 | 0.598 | 0.555 | 0.592 |
| Unsigned network (all weights +1) | unsigned | 0.598 | 0.543 | 0.596 |
| + measured CollecTRI targets of perturbed TFs (2,873 genes) | targets | 0.527 | 0.497 | 0.494 |
| Shared response removed | noshared | uninformative (all models near or below 0) | | |

Conclusion: the real network never beat the shuffled network or the expression-only baseline.
Real > shuffled for only 13 to 18 of 43 TFs. Real vs expression-only was slightly worse (Wilcoxon p <= 0.001 in base and targets).

The "targets" data: `scripts/prep_add_targets.py` adds every CollecTRI target of a perturbed TF that is detected in >= 1% of cells
(788 genes added; median measured targets per covered TF went from 4 to 11).

Commands: see "Data and reproducing every result" in the main README. The inputs
for `unsigned` and `noshared` come from `scripts/make_variants.py`; the `targets`
data from `scripts/prep_add_targets.py`.

Note for noshared: the shared response was the mean effect over all single perturbations (including held-out ones), so it is a diagnostic only.
