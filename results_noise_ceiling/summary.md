# Noise ceiling (split-half reliability)

5 random half-splits of every perturbation's cells and of the control cells; median cells per perturbation 495.

`half_vs_half`: agreement between the effects computed from the two halves. `ceiling`: Spearman-Brown estimate of the best achievable score against the full-data effect.

|                                |   all 105 singles |   34 covered |   71 not covered |
|:-------------------------------|------------------:|-------------:|-----------------:|
| pearson_top20 half_vs_half     |             0.959 |        0.987 |            0.946 |
| pearson_top20 ceiling          |             0.99  |        0.997 |            0.986 |
| pearson_all half_vs_half       |             0.785 |        0.858 |            0.75  |
| pearson_all ceiling            |             0.938 |        0.961 |            0.926 |
| pearson_ctr_top20 half_vs_half |             0.969 |        0.981 |            0.963 |
| pearson_ctr_top20 ceiling      |             0.992 |        0.995 |            0.991 |
| pearson_ctr_all half_vs_half   |             0.77  |        0.794 |            0.759 |
| pearson_ctr_all ceiling        |             0.933 |        0.941 |            0.929 |

Correlation between log(cells per perturbation) and centered reliability: 0.44 (positive = the low-cell perturbations are the noisy ones).
