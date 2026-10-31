# pert-prior-bench

**Does knowing the regulatory network help predict what an unseen genetic perturbation will do?**

Recent benchmarks found that simple linear baselines match or beat many single-cell
foundation models at predicting perturbation effects (Ahlmann-Eltze et al., 2025).
This project asks a narrower question with a clean ablation: take the same linear
model, and change only *where the perturbation's description comes from*:

| Model | Perturbation embedding | Output space |
|---|---|---|
| `mean` | none, predicts the average training effect | |
| `additive` | doubles only: effect(A) + effect(B) | |
| `linear[expression]` | perturbed gene's row in a PCA of training effects | expression PCs |
| `linear[prior]` | perturbed gene's row in an SVD of the TF network propagation matrix | prior target space |
| `linear[both]` | both, concatenated | both |
| `linear[both, shuffled prior]` | same as `both`, with the network's gene labels shuffled | |

The shuffled control is the point: if the prior "helps" only because the model
has more parameters, the shuffled prior will help just as much.

Every linear model also reports uncertainty (bootstrap ensemble), and the
benchmark checks whether that uncertainty means anything.

## Results on Norman et al. 2019

**Short answer: no.** On the Norman et al. 2019 CRISPRa Perturb-seq screen (K562;
105 single and 131 double perturbations, 2,085 genes), adding the CollecTRI network
did not improve prediction of unseen perturbations, even for the perturbed genes it
covers. A network with shuffled gene labels did as well.

CollecTRI covers 43 of the 105 perturbed genes (41%). Unseen singles, `pearson_top20`:

| model | all 105 | 43 in prior | 62 not in prior |
|---|---|---|---|
| mean baseline | 0.51 | 0.57 | 0.49 |
| **linear[expression]** | **0.52** | **0.60** | 0.48 |
| linear[prior] | 0.51 | 0.55 | 0.49 |
| linear[both] | 0.51 | 0.55 | 0.49 |
| linear[both, shuffled prior] | 0.52 | 0.59 | 0.49 |

For the 43 covered TFs, the real network beat the shuffled one for only 16, and
`linear[both]` was slightly but consistently worse than `linear[expression]`
(paired Wilcoxon p < 0.01). For doubles, the additive baseline won by a wide margin
(0.97 vs 0.73 for the best linear model), as in the synthetic control.

### The result holds across settings

Same table for the 43 covered TFs (full outputs in [`robustness/`](robustness/README.md)):

| check | linear[expression] | linear[both] (real) | linear[both, shuffled] |
|---|---|---|---|
| original (`--k-prior 50`) | 0.60 | 0.55 | 0.59 |
| `--k-prior 20` | 0.60 | 0.57 | 0.59 |
| `--k-prior 100` | 0.60 | 0.56 | 0.59 |
| unsigned network (all edges +1) | 0.60 | 0.54 | 0.60 |
| + measured targets of perturbed TFs (2,873 genes) | 0.53 | 0.50 | 0.49 |

The last row tests the most likely objection. In the original gene set, a covered TF
had a median of only 4 of its CollecTRI targets among the measured genes, so the
network had little to work with. `scripts/prep_add_targets.py` adds every target of a
perturbed TF detected in at least 1% of cells (788 genes), raising the median to 11.
The network still did not help.

A sixth check removed the response shared by all perturbations (about a third of the
variance). After that, no model, including the expression baseline, predicted the
remaining perturbation-specific effects, so that check cannot separate the models.

### What it means

Either curated networks carry less causal information about these perturbations
than is usually assumed, or a linear model cannot use it. Some caveats point to next
steps:

- **Cell type.** CollecTRI is not specific to K562, and many of its links may not be
  active there. Context-specific networks are the obvious next test.
- **Coverage.** Only 43 perturbed genes are TFs in the network, and 18 of them have
  5 or more measured targets in the original gene set. A screen that perturbs more
  TFs, such as a genome-scale Perturb-seq screen, would give the network more chances.
- **Signal.** Known targets are enriched among each TF's most-changed genes (about
  2.2-fold), but their direction agrees with the network's sign only 43% of the time.
  The links point to the right genes more often than they predict how those genes move.
- **Uncertainty.** As on synthetic data, the bootstrap intervals are overconfident
  (`coverage80` about 0.1 instead of 0.8).

## What gets measured

For each held-out perturbation (perturbed genes excluded from scoring):

- `pearson_top20`: correlation over the 20 genes the perturbation moves most. **Primary metric.**
- `pearson_all`, `mse`, `direction_top20`: secondary.
- `coverage80`: fraction of truths inside the predicted 80% interval (calibrated ≈ 0.80).
- `err_vs_std`: does the model say "unsure" on the perturbations it gets wrong? (Spearman, positive is good.)

Two experiments: **unseen singles** (5-fold CV over single-gene perturbations) and
**unseen doubles** (train on all singles, predict every combination).
The ridge penalty is chosen per model by inner cross-validation on training
perturbations only, so the test fold is never used for tuning.

## Positive control: synthetic data with a known network

`scripts/make_synthetic.py` simulates a screen where a known network drives the
responses, then hands the model a *noisy* copy of that network (30% of true edges
dropped, 20% spurious edges added). If the pipeline is sound, the prior should
help, only for genes the prior covers, and the shuffled prior should not.

Unseen singles, `pearson_top20` (run `bash run_synthetic.sh` to reproduce):

| model | mean | in prior | not in prior |
|---|---|---|---|
| mean baseline | 0.33 | 0.26 | 0.53 |
| linear[expression] | 0.33 | 0.27 | 0.53 |
| **linear[prior]** | **0.55** | **0.56** | 0.53 |
| linear[both] | 0.53 | 0.54 | 0.52 |
| linear[both, shuffled prior] | 0.33 | 0.27 | 0.53 |

It behaves as expected. Two things it also shows, worth carrying into the real analysis:

- For doubles, the plain **additive baseline wins** (0.97). Combinations are hard to beat when effects are close to additive.
- The bootstrap intervals are **overconfident**: `coverage80` is 0.2 to 0.5 instead of 0.8. Ensembles of the same model underestimate uncertainty. Fixing that (for example with conformal calibration) is a natural extension.

## Running it on real data

```bash
pip install -r requirements.txt

# 1. sanity check on synthetic data (about 20 seconds)
bash run_synthetic.sh

# 2. get the Norman et al. 2019 CRISPRa screen (K562, singles + doubles).
#    One source is scPerturb (scperturb.org), which distributes it as an .h5ad.
#    Put it at data/norman.h5ad, then look inside to set the flags:
python scripts/inspect_data.py --h5ad data/norman.h5ad

# 3. run (set --pert-key / --control / --sep / --gene-col from step 2)
python scripts/run_benchmark.py --h5ad data/norman.h5ad \
    --pert-key perturbation --control control --sep _ \
    --prior collectri --out results_norman
python scripts/plot_results.py --results results_norman
```

**Low memory?** The full Norman file is about 700 MB. `scripts/prep_large_h5ad.py`
reads it in column chunks and writes one mean log-normalized row per perturbation
(run it twice: the first pass picks genes, the second writes the file). Then run the
benchmark on that file with `--already-log --min-cells 1`. The results above were
produced this way, with CollecTRI saved as `data/collectri.csv`.

`--prior collectri` downloads the CollecTRI TF network through `decoupler`.
You can also pass any CSV with `source,target,weight` columns.

Outputs in `--out`: `summary.md` (tables), per-perturbation CSVs, and one figure per experiment.

## Four-week plan

**Week 1: baselines on real data.** Set up, run the synthetic control, download
Norman, run with `--prior none`. Expected: `linear[expression]` is close to the
mean baseline for singles and `additive` is strong for doubles, matching the
published benchmarks. Reproducing a known result first is what makes everything
after it trustworthy.

**Week 2: add the prior.** Run with CollecTRI. Record prior coverage (what
fraction of perturbed genes are TFs in the network). Read the in-prior /
not-in-prior split: that split is where any real effect should show up.

**Week 3: robustness.** Three seeds; `--k-prior` 20/50/100; `--k` 5/10/20.
Check that conclusions hold across settings, not just one. Stretch: replace the
bootstrap intervals with split-conformal intervals and compare `coverage80`.

**Week 4: write it up.** One-page summary in this README: the question, the
table, the figure, what you'd do next. Push to GitHub.

## Reading the result honestly

Any outcome is a usable result:

- **Prior helps, shuffled doesn't, and the gain is concentrated in covered genes:**
  evidence that curated regulation carries information the expression data lacks.
- **Prior doesn't help:** also informative. The next question is why: low
  coverage, wrong network for K562, or effects dominated by a shared stress response.
- **Prior helps and shuffled helps too:** the gain is capacity, not biology. Say so.

The thing to avoid is reporting only the configuration that looks best.

## Layout

```
src/pertbench/data.py       pseudobulk effects, singles/doubles, folds
src/pertbench/priors.py     CollecTRI loader, network propagation
src/pertbench/models.py     baselines and the linear model (+ bootstrap, nested ridge)
src/pertbench/evaluate.py   metrics and calibration
scripts/                    synthetic data, inspection, benchmark, plots
```
