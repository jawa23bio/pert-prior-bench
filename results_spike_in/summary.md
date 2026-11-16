# Spike-in sensitivity

34 covered perturbations; edge dropout 0.0; 2 x 5-fold CV; 19 null draws (smallest possible p = 0.050).

beta = size of the planted network-shaped effect relative to each perturbation's real specific effect. A `detects` at beta = 0 is a false positive.

|   beta |   v1 both-expr |   v1 ci_lo |   v1 ci_hi |   v1 both-shuffled | v1 detects   |   v2 direct[prior] real |   v2 direct[prior] null |   v2 direct[prior] p | v2 direct[prior] detects   |   v2 linear[both] real |   v2 linear[both] null |   v2 linear[both] p | v2 linear[both] detects   |
|-------:|---------------:|-----------:|-----------:|-------------------:|:-------------|------------------------:|------------------------:|---------------------:|:---------------------------|-----------------------:|-----------------------:|--------------------:|:--------------------------|
|   0    |         -0.048 |     -0.089 |     -0.016 |             -0.042 | False        |                   0     |                  -0.002 |                 0.5  | False                      |                 -0.037 |                  0.068 |                0.9  | False                     |
|   0.1  |         -0.052 |     -0.092 |     -0.018 |             -0.047 | False        |                   0.095 |                  -0.003 |                 0.05 | True                       |                  0.05  |                  0.091 |                0.8  | False                     |
|   0.25 |         -0.039 |     -0.078 |     -0.003 |             -0.042 | False        |                   0.333 |                   0.008 |                 0.05 | True                       |                  0.159 |                  0.106 |                0.3  | False                     |
|   0.5  |         -0.032 |     -0.09  |      0.021 |             -0.032 | False        |                   0.603 |                   0.024 |                 0.05 | True                       |                  0.238 |                  0.082 |                0.05 | True                      |
|   1    |          0.139 |      0.061 |      0.223 |              0.127 | True         |                   0.826 |                   0.031 |                 0.05 | True                       |                  0.353 |                  0.052 |                0.05 | True                      |

## Detection limits

- v1: beta >= 1.0
- v2 direct[prior]: beta >= 0.1
- v2 linear[both]: beta >= 0.5
