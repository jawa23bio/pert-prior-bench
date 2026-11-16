# Spike-in sensitivity

34 covered perturbations; edge dropout 0.3; 2 x 5-fold CV; 19 null draws (smallest possible p = 0.050).

beta = size of the planted network-shaped effect relative to each perturbation's real specific effect. A `detects` at beta = 0 is a false positive.

|   beta |   v1 both-expr |   v1 ci_lo |   v1 ci_hi |   v1 both-shuffled | v1 detects   |   v2 direct[prior] real |   v2 direct[prior] null |   v2 direct[prior] p | v2 direct[prior] detects   |   v2 linear[both] real |   v2 linear[both] null |   v2 linear[both] p | v2 linear[both] detects   |
|-------:|---------------:|-----------:|-----------:|-------------------:|:-------------|------------------------:|------------------------:|---------------------:|:---------------------------|-----------------------:|-----------------------:|--------------------:|:--------------------------|
|   0    |         -0.031 |     -0.091 |      0.032 |             -0.022 | False        |                  -0.003 |                  -0.001 |                 0.6  | False                      |                  0.147 |                  0.096 |                0.25 | False                     |
|   0.1  |         -0.032 |     -0.088 |      0.026 |             -0.025 | False        |                   0.107 |                   0.003 |                 0.05 | True                       |                  0.203 |                  0.12  |                0.15 | False                     |
|   0.25 |         -0.016 |     -0.06  |      0.029 |             -0.017 | False        |                   0.298 |                   0.013 |                 0.05 | True                       |                  0.245 |                  0.124 |                0.05 | True                      |
|   0.5  |         -0.008 |     -0.061 |      0.039 |             -0.005 | False        |                   0.515 |                   0.023 |                 0.05 | True                       |                  0.261 |                  0.079 |                0.05 | True                      |
|   1    |         -0.031 |     -0.118 |      0.049 |             -0.02  | False        |                   0.674 |                   0.03  |                 0.05 | True                       |                  0.278 |                  0.045 |                0.05 | True                      |

## Detection limits

- v1: not detected at any beta tried
- v2 direct[prior]: beta >= 0.1
- v2 linear[both]: beta >= 0.25
