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
| linear[prior]                |         0.338 |           0.539 | 0.054 |             0.771 |        0.236 |        0.618 |
| linear[both]                 |         0.41  |           0.508 | 0.026 |             0.734 |        0.251 |        0.487 |
| linear[both, shuffled prior] |         0.396 |           0.332 | 0.025 |             0.628 |        0.188 |       -0.092 |

## Unseen double perturbations (trained on singles)

|                              |   pearson_all |   pearson_top20 |   mse |   direction_top20 |   coverage80 |   err_vs_std |
|:-----------------------------|--------------:|----------------:|------:|------------------:|-------------:|-------------:|
| mean                         |         0.314 |           0.448 | 0.052 |             0.617 |        0.832 |        0.463 |
| additive                     |         0.783 |           0.969 | 0.03  |             0.998 |      nan     |      nan     |
| linear[expression]           |         0.437 |           0.716 | 0.049 |             0.863 |        0.224 |       -0.285 |
| linear[prior]                |         0.656 |           0.852 | 0.035 |             0.923 |        0.498 |        0.172 |
| linear[both]                 |         0.652 |           0.884 | 0.035 |             0.953 |        0.336 |        0.087 |
| linear[both, shuffled prior] |         0.455 |           0.747 | 0.048 |             0.872 |        0.249 |       -0.49  |

## Singles: pearson_top20 split by whether the perturbed gene is a regulator in the prior

| model                        |   in_prior=True |   in_prior=False |
|:-----------------------------|----------------:|-----------------:|
| mean                         |           0.264 |            0.532 |
| linear[expression]           |           0.265 |            0.532 |
| linear[prior]                |           0.541 |            0.532 |
| linear[both]                 |           0.504 |            0.52  |
| linear[both, shuffled prior] |           0.266 |            0.531 |
