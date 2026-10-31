# Results

```
{
  "n_singles": 105,
  "n_doubles": 131,
  "n_genes": 2085,
  "prior_coverage_singles": 0.4095238095238095
}
```

## Unseen single perturbations (5-fold CV)

|                              |   pearson_all |   pearson_top20 |   mse |   direction_top20 |   coverage80 |   err_vs_std |
|:-----------------------------|--------------:|----------------:|------:|------------------:|-------------:|-------------:|
| mean                         |         0.519 |           0.513 | 0.006 |             0.788 |        0.83  |       -0.121 |
| linear[expression]           |         0.519 |           0.52  | 0.006 |             0.799 |        0.113 |        0.094 |
| linear[prior]                |         0.517 |           0.513 | 0.006 |             0.787 |        0.098 |       -0.009 |
| linear[both]                 |         0.519 |           0.515 | 0.006 |             0.791 |        0.098 |       -0.059 |
| linear[both, shuffled prior] |         0.517 |           0.521 | 0.006 |             0.79  |        0.134 |        0.047 |

## Unseen double perturbations (trained on singles)

|                              |   pearson_all |   pearson_top20 |   mse |   direction_top20 |   coverage80 |   err_vs_std |
|:-----------------------------|--------------:|----------------:|------:|------------------:|-------------:|-------------:|
| mean                         |         0.604 |           0.633 | 0.013 |             0.846 |        0.681 |       -0.059 |
| additive                     |         0.881 |           0.966 | 0.005 |             0.998 |      nan     |      nan     |
| linear[expression]           |         0.662 |           0.724 | 0.012 |             0.885 |        0.113 |        0.026 |
| linear[prior]                |         0.604 |           0.636 | 0.013 |             0.846 |        0.065 |        0.119 |
| linear[both]                 |         0.611 |           0.648 | 0.013 |             0.85  |        0.066 |        0.119 |
| linear[both, shuffled prior] |         0.664 |           0.726 | 0.012 |             0.886 |        0.131 |        0.007 |

## Singles: pearson_top20 split by whether the perturbed gene is a regulator in the prior

| model                        |   in_prior=False |   in_prior=True |
|:-----------------------------|-----------------:|----------------:|
| mean                         |            0.487 |           0.569 |
| linear[expression]           |            0.483 |           0.598 |
| linear[prior]                |            0.487 |           0.569 |
| linear[both]                 |            0.487 |           0.573 |
| linear[both, shuffled prior] |            0.49  |           0.587 |
