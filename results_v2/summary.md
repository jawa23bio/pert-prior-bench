# Results (v2 protocol)

## Before fitting anything: can the network matter here?

```
{
  "n_genes": 2085,
  "n_singles": 105,
  "regulators_among_measured_genes": 123,
  "perturbed_genes_covered": 34,
  "median_measured_direct_targets_per_covered_tf": 5.0,
  "covered_tfs_with_5plus_targets": 18,
  "median_share_of_prior_row_kept_by_svd": 0.614,
  "median_transferable_share_svd_space": 0.899,
  "median_transferable_share_network_space": 0.269,
  "v1_control_covered_tfs_given_any_targets": "2 of 34"
}
```

## Unseen singles (5 x 5-fold CV; mean [95% CI over perturbations])

### pearson_ctr_top20

**Primary.** Perturbation-specific effect: (truth - training mean) vs (prediction - training mean), on the 20 genes that deviate most. 0 = no better than the mean.

|                                   | all 105                | 34 covered             | 71 not covered        |
|:----------------------------------|:-----------------------|:-----------------------|:----------------------|
| direct[prior]                     | -0.001 [-0.007, 0.004] | -0.002 [-0.021, 0.011] | 0.000 [0.000, 0.000]  |
| direct[prior], regulator-permuted | 0.000 [-0.005, 0.005]  | 0.000 [-0.015, 0.015]  | 0.000 [0.000, 0.000]  |
| direct[prior], target-permuted    | -0.000 [-0.002, 0.001] | -0.001 [-0.006, 0.004] | 0.000 [0.000, 0.000]  |
| linear[both]                      | 0.027 [-0.073, 0.123]  | -0.030 [-0.177, 0.121] | 0.054 [-0.066, 0.174] |
| linear[both], regulator-permuted  | 0.061 [-0.028, 0.149]  | 0.082 [-0.007, 0.175]  | 0.051 [-0.068, 0.171] |
| linear[both], target-permuted     | 0.028 [-0.073, 0.126]  | -0.028 [-0.179, 0.127] | 0.055 [-0.065, 0.176] |
| linear[expression]                | 0.100 [0.002, 0.197]   | 0.209 [0.055, 0.357]   | 0.048 [-0.070, 0.168] |
| mean                              | 0.000 [0.000, 0.000]   | 0.000 [0.000, 0.000]   | 0.000 [0.000, 0.000]  |

### pearson_ctr_all

Same, all genes.

|                                   | all 105                | 34 covered             | 71 not covered        |
|:----------------------------------|:-----------------------|:-----------------------|:----------------------|
| direct[prior]                     | -0.000 [-0.002, 0.001] | -0.001 [-0.007, 0.002] | 0.000 [0.000, 0.000]  |
| direct[prior], regulator-permuted | 0.000 [-0.001, 0.001]  | 0.000 [-0.002, 0.002]  | 0.000 [0.000, 0.000]  |
| direct[prior], target-permuted    | -0.000 [-0.000, 0.000] | -0.000 [-0.001, 0.000] | 0.000 [0.000, 0.000]  |
| linear[both]                      | 0.008 [-0.056, 0.070]  | -0.011 [-0.116, 0.093] | 0.018 [-0.058, 0.096] |
| linear[both], regulator-permuted  | 0.030 [-0.026, 0.086]  | 0.060 [-0.005, 0.126]  | 0.016 [-0.060, 0.093] |
| linear[both], target-permuted     | 0.009 [-0.056, 0.071]  | -0.010 [-0.118, 0.096] | 0.018 [-0.058, 0.096] |
| linear[expression]                | 0.048 [-0.012, 0.106]  | 0.117 [0.027, 0.205]   | 0.015 [-0.060, 0.091] |
| mean                              | 0.000 [0.000, 0.000]   | 0.000 [0.000, 0.000]   | 0.000 [0.000, 0.000]  |

### pearson_top20

The v1 metric (effect vs control). Dominated by the response shared by all perturbations.

|                                   | all 105              | 34 covered           | 71 not covered       |
|:----------------------------------|:---------------------|:---------------------|:---------------------|
| direct[prior]                     | 0.514 [0.426, 0.600] | 0.566 [0.446, 0.674] | 0.489 [0.367, 0.600] |
| direct[prior], regulator-permuted | 0.514 [0.426, 0.600] | 0.566 [0.446, 0.674] | 0.489 [0.367, 0.600] |
| direct[prior], target-permuted    | 0.514 [0.426, 0.600] | 0.566 [0.446, 0.674] | 0.489 [0.367, 0.600] |
| linear[both]                      | 0.515 [0.427, 0.599] | 0.556 [0.432, 0.667] | 0.495 [0.372, 0.604] |
| linear[both], regulator-permuted  | 0.519 [0.432, 0.604] | 0.576 [0.464, 0.677] | 0.492 [0.370, 0.602] |
| linear[both], target-permuted     | 0.515 [0.427, 0.599] | 0.556 [0.432, 0.667] | 0.495 [0.373, 0.604] |
| linear[expression]                | 0.517 [0.430, 0.601] | 0.574 [0.457, 0.680] | 0.489 [0.368, 0.599] |
| mean                              | 0.514 [0.426, 0.600] | 0.566 [0.446, 0.674] | 0.489 [0.367, 0.600] |

### mse

Mean squared error, all genes.

|                                   | all 105                 | 34 covered              | 71 not covered          |
|:----------------------------------|:------------------------|:------------------------|:------------------------|
| direct[prior]                     | 0.0062 [0.0051, 0.0074] | 0.0087 [0.0061, 0.0120] | 0.0050 [0.0042, 0.0058] |
| direct[prior], regulator-permuted | 0.0062 [0.0051, 0.0074] | 0.0087 [0.0061, 0.0120] | 0.0050 [0.0042, 0.0058] |
| direct[prior], target-permuted    | 0.0062 [0.0051, 0.0074] | 0.0087 [0.0061, 0.0120] | 0.0050 [0.0042, 0.0058] |
| linear[both]                      | 0.0063 [0.0052, 0.0075] | 0.0090 [0.0063, 0.0123] | 0.0050 [0.0042, 0.0059] |
| linear[both], regulator-permuted  | 0.0062 [0.0051, 0.0074] | 0.0087 [0.0062, 0.0120] | 0.0050 [0.0042, 0.0059] |
| linear[both], target-permuted     | 0.0063 [0.0052, 0.0075] | 0.0090 [0.0063, 0.0123] | 0.0050 [0.0042, 0.0059] |
| linear[expression]                | 0.0064 [0.0053, 0.0077] | 0.0087 [0.0061, 0.0120] | 0.0053 [0.0044, 0.0064] |
| mean                              | 0.0062 [0.0051, 0.0074] | 0.0087 [0.0061, 0.0120] | 0.0050 [0.0042, 0.0058] |

## Is the real network better than a wrong one? (covered perturbations)

Statistic: mean score over covered perturbations. `p` = share of null draws scoring at least as high as the real network (one-sided; the smallest possible value is 1/(1+draws)).

|                                                        |    real |   null_mean |   null_lo |   null_hi |   n_draws |      p |
|:-------------------------------------------------------|--------:|------------:|----------:|----------:|----------:|-------:|
| direct[prior] / regulator-permuted / pearson_ctr_top20 | -0.0021 |      0.0005 |   -0.0313 |    0.0342 |        20 | 0.5238 |
| direct[prior] / target-permuted / pearson_ctr_top20    | -0.0021 |     -0.001  |   -0.0303 |    0.0222 |        20 | 0.7143 |
| direct[prior] / regulator-permuted / pearson_top20     |  0.5656 |      0.5656 |    0.5655 |    0.5659 |        20 | 0.8571 |
| direct[prior] / target-permuted / pearson_top20        |  0.5656 |      0.5656 |    0.5656 |    0.5657 |        20 | 1      |
| linear[both] / regulator-permuted / pearson_ctr_top20  | -0.0303 |      0.0822 |   -0.0463 |    0.1692 |        20 | 0.9048 |
| linear[both] / target-permuted / pearson_ctr_top20     | -0.0303 |     -0.0281 |   -0.0359 |   -0.0242 |        20 | 0.9524 |
| linear[both] / regulator-permuted / pearson_top20      |  0.556  |      0.5759 |    0.5557 |    0.6216 |        20 | 0.9524 |
| linear[both] / target-permuted / pearson_top20         |  0.556  |      0.5562 |    0.5554 |    0.5566 |        20 | 0.9048 |

## Paired differences on covered perturbations (pearson_top20)

| comparison                         |   mean_diff |   ci_lo |   ci_hi | wins   |   wilcoxon_p |
|:-----------------------------------|------------:|--------:|--------:|:-------|-------------:|
| direct[prior] - linear[expression] |     -0.0084 | -0.0222 |  0.0049 | 12/34  |       0.0907 |
| linear[both] - linear[expression]  |     -0.018  | -0.0523 |  0.0058 | 13/34  |       0.2415 |

## Uncertainty (first CV repeat)

`coverage80` should be 0.80. `half_width` is the price of honest coverage (narrower = the model knows more). `err_vs_unc`: does the model's own per-perturbation uncertainty rank its errors? (Spearman; identical for both interval types because both rank by the model's own std. For the mean baseline that std is the spread of training effects, which is why its v1 interval was never far off.)

|                                            |   coverage80 |   half_width |   rmse |   err_vs_unc |
|:-------------------------------------------|-------------:|-------------:|-------:|-------------:|
| mean / v1 (1.28 x model std)               |        0.831 |        0.086 |  0.072 |       -0.042 |
| mean / out-of-fold conformal               |        0.801 |        0.079 |  0.072 |       -0.042 |
| linear[expression] / v1 (1.28 x model std) |        0.102 |        0.007 |  0.072 |        0.03  |
| linear[expression] / out-of-fold conformal |        0.804 |        0.079 |  0.072 |        0.03  |
| direct[prior] / v1 (1.28 x model std)      |        0     |        0     |  0.072 |        0.204 |
| direct[prior] / out-of-fold conformal      |        0.801 |        0.079 |  0.072 |        0.204 |

## Unseen doubles (trained on all singles)

The additive baseline has seen both singles, so it is the bar to clear.

| model              |   pearson_top20 |   pearson_all |   mse |
|:-------------------|----------------:|--------------:|------:|
| mean               |           0.633 |         0.604 | 0.013 |
| additive           |           0.966 |         0.881 | 0.005 |
| linear[expression] |           0.686 |         0.63  | 0.012 |
| direct[prior]      |           0.633 |         0.604 | 0.013 |
| linear[both]       |           0.709 |         0.646 | 0.012 |
