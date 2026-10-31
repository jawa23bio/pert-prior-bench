# Results

```
{
  "n_singles": 105,
  "n_doubles": 131,
  "n_genes": 2085
}
```

## Unseen single perturbations (5-fold CV)

|                    |   pearson_all |   pearson_top20 |   mse |   direction_top20 |   coverage80 |   err_vs_std |
|:-------------------|--------------:|----------------:|------:|------------------:|-------------:|-------------:|
| mean               |         0.519 |           0.513 | 0.006 |             0.788 |        0.83  |       -0.121 |
| linear[expression] |         0.519 |           0.52  | 0.006 |             0.799 |        0.113 |        0.094 |

## Unseen double perturbations (trained on singles)

|                    |   pearson_all |   pearson_top20 |   mse |   direction_top20 |   coverage80 |   err_vs_std |
|:-------------------|--------------:|----------------:|------:|------------------:|-------------:|-------------:|
| mean               |         0.604 |           0.633 | 0.013 |             0.846 |        0.681 |       -0.059 |
| additive           |         0.881 |           0.966 | 0.005 |             0.998 |      nan     |      nan     |
| linear[expression] |         0.662 |           0.724 | 0.012 |             0.885 |        0.113 |        0.026 |
