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
| mean                         |        -0.191 |          -0.259 | 0.006 |             0.356 |        0.83  |       -0.121 |
| linear[expression]           |        -0.073 |          -0.069 | 0.007 |             0.452 |        0.168 |        0.152 |
| linear[prior]                |        -0.15  |          -0.195 | 0.199 |             0.368 |        0.106 |        0.337 |
| linear[both]                 |        -0.09  |          -0.096 | 2.301 |             0.445 |        0.293 |        0.563 |
| linear[both, shuffled prior] |        -0.086 |          -0.082 | 0.008 |             0.451 |        0.208 |        0.218 |

## Unseen double perturbations (trained on singles)

|                              |   pearson_all |   pearson_top20 |   mse |   direction_top20 |   coverage80 |   err_vs_std |
|:-----------------------------|--------------:|----------------:|------:|------------------:|-------------:|-------------:|
| mean                         |        -0     |          -0.036 | 0.013 |             0.478 |        0.681 |       -0.059 |
| additive                     |         0.715 |           0.903 | 0.007 |             0.968 |      nan     |      nan     |
| linear[expression]           |         0.293 |           0.424 | 0.011 |             0.692 |        0.237 |       -0.007 |
| linear[prior]                |         0.098 |           0.14  | 0.013 |             0.552 |        0.075 |        0.145 |
| linear[both]                 |         0.44  |           0.545 | 0.009 |             0.793 |        0.535 |       -0.134 |
| linear[both, shuffled prior] |         0.315 |           0.443 | 0.011 |             0.711 |        0.288 |       -0.057 |

## Singles: pearson_top20 split by whether the perturbed gene is a regulator in the prior

| model                        |   in_prior=False |   in_prior=True |
|:-----------------------------|-----------------:|----------------:|
| mean                         |           -0.247 |          -0.284 |
| linear[expression]           |           -0.115 |           0.027 |
| linear[prior]                |           -0.247 |          -0.086 |
| linear[both]                 |           -0.096 |          -0.096 |
| linear[both, shuffled prior] |           -0.135 |           0.027 |
