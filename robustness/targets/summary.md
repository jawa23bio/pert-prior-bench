# Results

```
{
  "n_singles": 105,
  "n_doubles": 131,
  "n_genes": 2873,
  "prior_coverage_singles": 0.4095238095238095
}
```

## Unseen single perturbations (5-fold CV)

|                              |   pearson_all |   pearson_top20 |   mse |   direction_top20 |   coverage80 |   err_vs_std |
|:-----------------------------|--------------:|----------------:|------:|------------------:|-------------:|-------------:|
| mean                         |         0.528 |           0.506 | 0.005 |             0.786 |        0.835 |       -0.115 |
| linear[expression]           |         0.529 |           0.514 | 0.005 |             0.794 |        0.109 |        0.018 |
| linear[prior]                |         0.528 |           0.507 | 0.005 |             0.787 |        0.1   |        0.063 |
| linear[both]                 |         0.524 |           0.502 | 0.005 |             0.788 |        0.113 |        0.152 |
| linear[both, shuffled prior] |         0.522 |           0.498 | 0.351 |             0.782 |        0.12  |        0.004 |

## Unseen double perturbations (trained on singles)

|                              |   pearson_all |   pearson_top20 |   mse |   direction_top20 |   coverage80 |   err_vs_std |
|:-----------------------------|--------------:|----------------:|------:|------------------:|-------------:|-------------:|
| mean                         |         0.609 |           0.628 | 0.01  |             0.845 |        0.689 |       -0.07  |
| additive                     |         0.881 |           0.966 | 0.004 |             0.998 |      nan     |      nan     |
| linear[expression]           |         0.666 |           0.719 | 0.009 |             0.887 |        0.113 |        0.025 |
| linear[prior]                |         0.61  |           0.629 | 0.01  |             0.845 |        0.066 |        0.338 |
| linear[both]                 |         0.619 |           0.65  | 0.01  |             0.853 |        0.078 |        0.312 |
| linear[both, shuffled prior] |         0.677 |           0.735 | 0.009 |             0.893 |        0.161 |        0.033 |

## Singles: pearson_top20 split by whether the perturbed gene is a regulator in the prior

| model                        |   in_prior=False |   in_prior=True |
|:-----------------------------|-----------------:|----------------:|
| mean                         |            0.506 |           0.507 |
| linear[expression]           |            0.506 |           0.527 |
| linear[prior]                |            0.506 |           0.508 |
| linear[both]                 |            0.506 |           0.497 |
| linear[both, shuffled prior] |            0.501 |           0.494 |
