# Results

```
{
  "n_singles": 80,
  "n_doubles": 30,
  "n_genes": 300,
  "prior_coverage_singles": 0.75
}
```

## Unseen single perturbations (5-fold CV)

|                              |   pearson_all |   pearson_top20 |   mse |   direction_top20 |   coverage80 |   err_vs_std |
|:-----------------------------|--------------:|----------------:|------:|------------------:|-------------:|-------------:|
| mean                         |         0.396 |           0.331 | 0.025 |             0.628 |        0.917 |       -0.187 |
| linear[expression]           |         0.396 |           0.332 | 0.025 |             0.631 |        0.188 |       -0.099 |
| linear[prior]                |         0.345 |           0.545 | 0.05  |             0.771 |        0.23  |        0.558 |
| linear[both]                 |         0.401 |           0.537 | 0.032 |             0.731 |        0.264 |        0.627 |
| linear[both, shuffled prior] |         0.396 |           0.331 | 0.025 |             0.629 |        0.188 |       -0.114 |

## Unseen double perturbations (trained on singles)

|                              |   pearson_all |   pearson_top20 |   mse |   direction_top20 |   coverage80 |   err_vs_std |
|:-----------------------------|--------------:|----------------:|------:|------------------:|-------------:|-------------:|
| mean                         |         0.318 |           0.403 | 0.051 |             0.617 |        0.829 |        0.263 |
| additive                     |         0.783 |           0.96  | 0.031 |             0.998 |      nan     |      nan     |
| linear[expression]           |         0.441 |           0.687 | 0.048 |             0.858 |        0.221 |       -0.175 |
| linear[prior]                |         0.652 |           0.832 | 0.035 |             0.922 |        0.485 |        0.237 |
| linear[both]                 |         0.651 |           0.868 | 0.035 |             0.95  |        0.325 |       -0.187 |
| linear[both, shuffled prior] |         0.356 |           0.51  | 0.05  |             0.693 |        0.112 |       -0.493 |

## Singles: pearson_top20 split by whether the perturbed gene is a regulator in the prior

| model                        |   in_prior=True |   in_prior=False |
|:-----------------------------|----------------:|-----------------:|
| mean                         |           0.264 |            0.532 |
| linear[expression]           |           0.265 |            0.532 |
| linear[prior]                |           0.55  |            0.532 |
| linear[both]                 |           0.545 |            0.514 |
| linear[both, shuffled prior] |           0.265 |            0.531 |
