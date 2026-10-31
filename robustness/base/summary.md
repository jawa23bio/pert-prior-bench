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
| linear[prior]                |         0.516 |           0.508 | 0.006 |             0.786 |        0.097 |       -0.084 |
| linear[both]                 |         0.516 |           0.507 | 0.007 |             0.788 |        0.109 |       -0.02  |
| linear[both, shuffled prior] |         0.517 |           0.521 | 0.006 |             0.79  |        0.131 |        0.088 |

## Unseen double perturbations (trained on singles)

|                              |   pearson_all |   pearson_top20 |   mse |   direction_top20 |   coverage80 |   err_vs_std |
|:-----------------------------|--------------:|----------------:|------:|------------------:|-------------:|-------------:|
| mean                         |         0.604 |           0.633 | 0.013 |             0.846 |        0.681 |       -0.059 |
| additive                     |         0.881 |           0.966 | 0.005 |             0.998 |      nan     |      nan     |
| linear[expression]           |         0.662 |           0.724 | 0.012 |             0.885 |        0.113 |        0.026 |
| linear[prior]                |         0.605 |           0.637 | 0.013 |             0.846 |        0.064 |        0.197 |
| linear[both]                 |         0.614 |           0.656 | 0.013 |             0.851 |        0.07  |        0.133 |
| linear[both, shuffled prior] |         0.666 |           0.73  | 0.011 |             0.885 |        0.14  |       -0.011 |

## Singles: pearson_top20 split by whether the perturbed gene is a regulator in the prior

| model                        |   in_prior=False |   in_prior=True |
|:-----------------------------|-----------------:|----------------:|
| mean                         |            0.487 |           0.569 |
| linear[expression]           |            0.483 |           0.598 |
| linear[prior]                |            0.487 |           0.554 |
| linear[both]                 |            0.487 |           0.55  |
| linear[both, shuffled prior] |            0.487 |           0.592 |
