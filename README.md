# pert-prior-bench

**Does knowing the regulatory network help predict what an unseen genetic perturbation will do?**

![Left: scores on what is specific to each perturbation for the 34 perturbed TFs with a measured network target. No network 0.21; with the CollecTRI network -0.03; with permuted networks 0.08 on average; reading the network edges directly about 0, p = 0.55 against 99 permuted networks; two halves of the same cells agree at 0.98. Right: the direct read-out of the real network detects a planted network effect one tenth the size of each TF's own effect, while permuted networks stay near 0.](docs/figures/overview.png)

*Norman et al. 2019 CRISPRa Perturb-seq screen (K562). Left: a model with no network predicts part of each
TF's specific response; the curated network adds nothing and is no better than permuted networks. Right: the
same test detects a planted network effect one tenth the size of each TF's own effect. Details below;
`python scripts/make_figure.py` redraws it from the committed results.*

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
| `direct[prior]` | none: the training mean plus one fitted scale times the TF's own row of the network | genes |

`direct[prior]` reads the network the way it is usually read ("activating this TF
moves these genes"), with a single number to fit. If it cannot beat the mean, the
edges themselves are not predictive, whatever model is built on top of them.

Two kinds of control:

- **Wrong networks.** Every model that reads the network is rerun on networks of the
  same shape that are wrong about which genes each TF regulates: TFs swap target sets
  (*regulator-permuted*), or each TF keeps its number of targets and its overlap with
  other TFs but on random genes (*target-permuted*). Repeated draws give a null
  distribution and an empirical p-value.
- **A known answer.** A synthetic screen generated from a known network, and a
  spike-in that plants a network effect of chosen size into the real data, show how
  large an effect the benchmark can detect.

Every model also reports 80% prediction intervals, and the benchmark checks whether
they hold and whether a model knows which perturbations it will get wrong.

The first version of this benchmark (v1, kept for reproducibility) had weaker
controls and a metric that hid most of the signal. What changed and why is in
[`docs/V2.md`](docs/V2.md); the v1 results are at the end of the next section.

## Results on Norman et al. 2019

**Short answer: no.** On the Norman et al. 2019 CRISPRa Perturb-seq screen (K562;
105 single and 131 double perturbations, 2,085 genes), the CollecTRI network did not
help predict unseen perturbations, and a model that reads its edges directly could
not be told apart from a randomly rewired network.

The screen can only test this for some perturbations. 43 of the 105 perturbed genes
are regulators in CollecTRI, and 34 of those have at least one CollecTRI target among
the measured genes (median 5; 18 have five or more). The other 9 (CBFA2T3, FOXF1,
HES7, HOXB9, ISL2, MAML2, OSR2, RUNX1T1, ZBTB10) get nothing from the network, and the
remaining 62 perturbed genes are not regulators in it. Network comparisons below use
the 34.

### What is specific to each perturbation

Most of a perturbation's effect is a response shared by all perturbations, so
predicting the average training effect already scores 0.51 on a control-referenced
metric (see v1 below). The primary metric therefore scores only what is specific:
the Pearson correlation between (truth minus training mean) and (prediction minus
training mean) over the 20 genes that deviate most. The mean baseline scores 0 by
definition. 5-fold cross-validation repeated 5 times, folds stratified by network
coverage; brackets are 95% confidence intervals over perturbations.

| model | 34 TFs with a measured target | 71 other perturbations |
|---|---|---|
| mean baseline | 0 | 0 |
| **linear[expression]** | **0.21** [0.06, 0.36] | 0.05 [-0.07, 0.17] |
| linear[both], real network | -0.03 [-0.18, 0.12] | 0.05 [-0.07, 0.17] |
| linear[both], regulator-permuted (20 draws) | 0.08 [-0.01, 0.18] | 0.05 |
| linear[both], target-permuted (20 draws) | -0.03 | 0.05 |
| direct[prior], real network | -0.00 [-0.02, 0.01] | 0 (no network row) |
| direct[prior], regulator-permuted (99 draws) | -0.00 | 0 |

- **A model with no network predicts some of the specific effect for TFs.** 20 of the
  34 score above zero (Wilcoxon p = 0.019), led by the myeloid TFs (CEBPA, FOSB,
  CEBPB, CEBPE, SPI1). It is modest, it is close to nothing for the 71 other
  perturbations, and no model beats the mean on mean squared error.
- **Adding the network removes that signal.** With the real network, `linear[both]`
  scores below 18 of 20 regulator-permuted networks: the real network cannot be told
  apart from a wrong one and, if anything, sits at the low end. It scores the same as
  the target-permuted networks, so the drop comes from the network features
  themselves, not from which genes they point to.
- **Read directly, the edges predict nothing.** `direct[prior]` scores -0.002 against a
  null mean of -0.005 (p = 0.55, 99 draws). Its fitted scale is zero in 92% of fits;
  left free to go negative, it is slightly negative.

### Could the benchmark have seen an effect?

A negative result only means something if a positive one was possible. Two checks:

- **Detection limit** (`scripts/spike_in.py`). For every covered TF, a planted effect
  pointing along its network row is added to the real data, sized relative to that
  perturbation's own specific effect (beta). `direct[prior]` detects beta = 0.1, also
  when 30% of the network's edges are removed from the copy it is given.
  `linear[both]` needs 0.5 (0.25 with edges removed). The original v1 comparison
  needs beta = 1, and with edges removed it detects nothing up to 1. Nothing fires at
  beta = 0. The planted effect follows the same network rows the models see, so this
  is the favourable case.
- **Noise ceiling** (`scripts/noise_ceiling.py`). Splitting each perturbation's cells
  in half and scoring one half against the other gives 0.97 on the same centered
  metric (0.98 for the 34 TFs), an estimated ceiling of about 0.99 against the full
  data. These are cells from one experiment, so this is not a between-experiment
  replicate.

So within this screen, a network effect one tenth the size of a perturbation's own
specific effect would have been detected, the effects being predicted are almost
noise-free, and the best model reaches 0.21 of a possible 0.99.

### Uncertainty

Intervals built from out-of-fold residuals (`src/pertbench/calibrate.py`) cover 80%
of outcomes for every model, at the same half-width (0.079). No model earns narrower
intervals than the mean baseline, and a model's own uncertainty does not rank which
perturbations it gets wrong (Spearman 0.03 for `linear[expression]`). The v1
intervals, the bootstrap spread of the fitted model, covered about 10%: that spread
measures how much the fit moves when the training set changes, not how far a new
perturbation lands from it, so it under-covers by construction.

### Doubles

Trained on all singles, the additive baseline wins by a wide margin (pearson_top20
0.97, against 0.71 for the best linear model). It has already seen both halves of
every pair, so what is left to predict is the genetic interaction, which no model
here attempts.

### Original protocol (v1)

v1 used one 5-fold split and scored each prediction against control
(`pearson_top20`), so every score includes the shared response.

| model | all 105 | 34 TFs with a measured target | 71 others |
|---|---|---|---|
| mean baseline | 0.51 | 0.57 | 0.49 |
| **linear[expression]** | **0.52** | **0.60** | 0.48 |
| linear[prior] | 0.51 | 0.55 | 0.49 |
| linear[both] | 0.51 | 0.55 | 0.49 |
| linear[both], all gene labels shuffled | 0.52 | 0.59 | 0.49 |

No model beat the mean baseline by more than 0.01 overall, which is what the centered
metric above corrects. For the 34 TFs, `linear[both]` scored below
`linear[expression]` for 26 (mean difference -0.048, 95% CI -0.088 to -0.015;
Wilcoxon p = 0.0008), and 73% of that gap comes from five TFs (KMT2A, KLF1, CEBPE,
ETS2, RREB1). v1's control shuffled every gene label in the network. Because only a
few percent of measured genes are regulators, that left 2 of the 34 TFs with any
targets, so the "shuffled" model behaved like a model with no network. The permuted
networks above replace it.

The same comparison held across settings (full outputs in
[`robustness/`](robustness/README.md)):

| check | linear[expression] | linear[both] (real) | linear[both], shuffled |
|---|---|---|---|
| original (`--k-prior 50`), 34 TFs | 0.60 | 0.55 | 0.59 |
| `--k-prior 20`, 34 TFs | 0.60 | 0.57 | 0.59 |
| `--k-prior 100`, 34 TFs | 0.60 | 0.56 | 0.59 |
| unsigned network (all edges +1), 34 TFs | 0.60 | 0.54 | 0.60 |
| + measured targets of perturbed TFs (2,873 genes), 41 TFs | 0.53 | 0.50 | 0.49 |

The last row tests the most likely objection, that the network had little to work
with. `scripts/prep_add_targets.py` adds every CollecTRI target of a perturbed TF
detected in at least 1% of cells (788 genes). That gives 41 of the 43 regulators at
least one measured target and raises the median from 5 to 11 among those that have
any. The network still did not help.

A further run removed the response shared by all perturbations before fitting. It
subtracted a mean that included the held-out perturbations, which pushes every
model below zero by construction (`tests/` reproduces this on pure noise), so its
numbers say nothing about the models. The centered metric above does this
correctly, with the mean taken from training perturbations only.

### What it means

In this screen, CollecTRI's edges do not predict how genes respond when their TF is
activated: the one-parameter model that uses them directly is at chance, an effect a
tenth as large as each perturbation's own would have been seen, and the responses
are reproducible. That is a result about this network in this experiment, not about
regulatory networks in general. Caveats and next steps:

- **Cell type.** CollecTRI is not specific to K562, and many of its links may not be
  active there. A K562-specific network is the obvious next test.
- **Coverage.** Only 34 perturbed TFs have a measured target, and 18 have five or
  more. A screen that perturbs more TFs, such as a genome-scale Perturb-seq screen,
  would give the network more chances.
- **Which genes, not which way.** Known targets are enriched among each TF's 100
  most-changed genes (57 vs about 26 expected, 2.2-fold), but over half of these hits
  come from three TFs (JUN, FOSB, SPI1). The direction is not predicted: targets move
  in the direction of the edge sign 43% of the time, against 38% for random genes
  given the same signs (most genes go down after a perturbation while most edges say
  "activates"). The network carries some information about which genes respond,
  concentrated in a few TFs, and none that improves prediction.
  `python scripts/diagnose_targets.py` reproduces these numbers.
- **The gap.** The best model reaches 0.21 against a ceiling near 0.99. What closes
  that gap is open; `docs/V2.md` lists candidates.

## What gets measured

For each held-out perturbation (perturbed genes excluded from scoring):

- `pearson_ctr_top20`: correlation of (truth minus training mean) with (prediction
  minus training mean) over the 20 genes that deviate most. **Primary metric.**
- `pearson_ctr_all`: the same over all genes.
- `pearson_top20`, `pearson_all`, `mse`, `direction_top20`: scored against control
  (`pearson_top20` was the v1 primary metric).
- `coverage80`: fraction of truths inside the predicted 80% interval (calibrated ≈ 0.80).
- `err_vs_unc` (`err_vs_std` in v1): does the model say "unsure" on the
  perturbations it gets wrong? (Spearman, positive is good.)

Two experiments: **unseen singles** (5-fold cross-validation over single-gene
perturbations, repeated 5 times with folds stratified by network coverage) and
**unseen doubles** (train on all singles, predict every combination). The ridge
penalty is chosen per model by inner cross-validation on training perturbations
only, so the test fold is never used for tuning.

## Positive control: synthetic data with a known network

`scripts/make_synthetic.py` simulates a screen where a known network drives the
responses, then hands the model a *noisy* copy of that network (30% of true edges
dropped, 20% spurious edges added). If the pipeline is sound, the prior should
help, and only for genes the prior covers.

Unseen singles, `pearson_top20`, v1 protocol (run `bash run_synthetic.sh` to reproduce):

| model | all 80 | 60 in prior | 20 not in prior |
|---|---|---|---|
| mean baseline | 0.33 | 0.26 | 0.53 |
| linear[expression] | 0.33 | 0.27 | 0.53 |
| **linear[prior]** | **0.55** | **0.55** | 0.53 |
| linear[both] | 0.54 | 0.55 | 0.51 |
| linear[both], all gene labels shuffled | 0.33 | 0.27 | 0.53 |

It behaves as expected. Two limits: the shuffled row is close to "no network" here
too (only 16 of the 60 covered TFs keep any targets after the shuffle), and the
synthetic network is denser than CollecTRI and generated by the same propagation
rule the model uses, so it is an easier problem than the real one. The spike-in on
the real data is the matched positive control. For doubles, the additive baseline
wins here as well (0.96).

## Running it on real data

```bash
pip install -r requirements.txt        # or: pip install -e ".[test]"
pytest -q                              # about 15 seconds, no data needed

# 1. sanity check on synthetic data (about 20 seconds)
bash run_synthetic.sh

# 2. get the Norman et al. 2019 CRISPRa screen (K562, singles + doubles).
#    One source is scPerturb (scperturb.org), which distributes it as an .h5ad.
#    Put it at data/norman.h5ad, then look inside to set the flags:
python scripts/inspect_data.py --h5ad data/norman.h5ad

# 3. run (set --pert-key / --control / --sep / --gene-col from step 2)
python scripts/run_v2.py --h5ad data/norman.h5ad \
    --pert-key perturbation --control control --sep _ \
    --prior collectri --calibrate --out results_v2
```

**Low memory?** The full Norman file is about 700 MB. `scripts/prep_large_h5ad.py`
reads it in column chunks and writes one mean log-normalized row per perturbation
(run it twice: the first pass picks genes, the second writes the file). Then run the
benchmark on that file with `--already-log --min-cells 1`. The results above were
produced this way, with CollecTRI saved as `data/collectri.csv`.

`--prior collectri` downloads the CollecTRI TF network through `decoupler`.
You can also pass any CSV with `source,target,weight` columns.

Outputs in `--out`: `summary.md` (tables) and per-perturbation CSVs.

## Data and reproducing every result

The screen is not in the repo (about 700 MB); the CollecTRI files are. To rebuild the rest:

| file | source | how |
|---|---|---|
| `data/norman.h5ad` | Norman et al. 2019, from scPerturb ([Zenodo record 10044268](https://zenodo.org/records/10044268), v1.4, `NormanWeissman2019_filtered.h5ad`) | download and rename |
| `data/collectri.tsv` | CollecTRI from OmniPath: `https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes` | download (included) |
| `data/collectri.csv` | the TSV above | `python scripts/make_collectri_csv.py` (included) |
| `data/norman_pb.h5ad` (+ `.stage1.npz`) | the screen, one mean row per perturbation, 2,085 genes | `python scripts/prep_large_h5ad.py --h5ad data/norman.h5ad --out data/norman_pb.h5ad --n-hvg 2000`, run twice |
| `data/norman_pb_targets.h5ad` | same, plus measured CollecTRI targets of the perturbed TFs | `python scripts/prep_add_targets.py` |
| `data/collectri_unsigned.csv`, `data/norman_pb_noshared.h5ad` | inputs for two v1 robustness runs | `python scripts/make_variants.py` |

The v2 results:

```bash
A="--h5ad data/norman_pb.h5ad --already-log --min-cells 1 --sep _ --n-hvg 2085 --prior data/collectri.csv"
python scripts/run_v2.py   $A --calibrate                         --out results_v2                # ~6 min
python scripts/run_v2.py   $A --n-null 99 --no-linear-prior       --out results_v2_direct         # ~3 min
python scripts/spike_in.py $A                                     --out results_spike_in          # ~6 min
python scripts/spike_in.py $A --edge-dropout 0.3                  --out results_spike_in_dropout  # ~6 min
python scripts/noise_ceiling.py --h5ad data/norman.h5ad --sep _ --genes-from data/norman_pb.h5ad \
    --stage1 data/norman_pb.h5ad.stage1.npz --prior data/collectri.csv --out results_noise_ceiling  # <1 min
python scripts/make_figure.py                                     # docs/figures/overview.png
```

The v1 results (each run about two minutes):

```bash
B="python scripts/run_benchmark.py --already-log --min-cells 1"
$B --h5ad data/norman_pb.h5ad --n-hvg 2085 --prior none                 --out results_none
$B --h5ad data/norman_pb.h5ad --n-hvg 2085 --prior data/collectri.csv   --out results_norman
$B --h5ad data/norman_pb.h5ad --n-hvg 2085 --prior data/collectri.csv --k-prior 20  --out robustness/k20
$B --h5ad data/norman_pb.h5ad --n-hvg 2085 --prior data/collectri.csv --k-prior 100 --out robustness/k100
$B --h5ad data/norman_pb.h5ad --n-hvg 2085 --prior data/collectri_unsigned.csv      --out robustness/unsigned
$B --h5ad data/norman_pb_targets.h5ad  --n-hvg 2873 --prior data/collectri.csv      --out robustness/targets
$B --h5ad data/norman_pb_noshared.h5ad --n-hvg 2085 --prior data/collectri.csv      --out robustness/noshared
python scripts/plot_results.py --results results_norman
```

`robustness/base` is the same run as `results_norman`. `bash run_synthetic.sh`
reproduces `results_synthetic` (fixed seed; identical up to floating-point
differences of about 1e-8 across numpy versions).

## Layout

```
src/pertbench/data.py         pseudobulk effects, singles/doubles, folds
src/pertbench/priors.py       CollecTRI loader, network propagation, permuted null networks
src/pertbench/models.py       baselines, the linear model (+ bootstrap, nested ridge), direct[prior]
src/pertbench/evaluate.py     metrics, including the centered metric
src/pertbench/protocol.py     v2: repeated stratified CV, null networks, summaries
src/pertbench/calibrate.py    v2: prediction intervals from out-of-fold residuals
scripts/run_v2.py             the benchmark (v2)
scripts/spike_in.py           how large an effect could the benchmark detect?
scripts/noise_ceiling.py      split-half reliability (best achievable score)
scripts/make_figure.py        the overview figure at the top of this README
scripts/run_benchmark.py      the original benchmark (v1)
scripts/plot_results.py       figures (v1)
scripts/make_synthetic.py     synthetic positive control
scripts/make_norman_like.py   data-free stand-in with CollecTRI's real sparsity
scripts/inspect_data.py       look inside an .h5ad to set the flags
scripts/prep_large_h5ad.py    low-memory pseudobulk of a large screen
scripts/prep_add_targets.py   pseudobulk with the prior's targets added
scripts/make_collectri_csv.py OmniPath download -> source,target,weight CSV
scripts/make_variants.py      unsigned prior and shared-response-removed data (v1)
scripts/diagnose_targets.py   are network targets enriched and correctly signed?
tests/                        pytest: controls, metric, models, intervals, v1 reproduction
docs/V2.md                    what changed from v1 and why
docs/figures/                 overview figure (PNG and SVG)
results_v2/                   Norman, v2 (main result)
results_v2_direct/            Norman, direct[prior] against 99 null draws
results_spike_in/             detection limit (and _dropout: 30% of edges removed)
results_noise_ceiling/        split-half reliability
results_synthetic/            positive control (v1)
results_none/                 Norman, no prior, v1
results_norman/               Norman with CollecTRI, v1
robustness/                   four v1 variations that held, plus the flawed shared-response run
```

## References

- Norman TM, Horlbeck MA, Replogle JM, et al. Exploring genetic interaction manifolds
  constructed from rich single-cell phenotypes. *Science* 365, 786-793 (2019).
- Müller-Dott S, et al. Expanding the coverage of regulons from high-confidence prior
  knowledge for accurate estimation of transcription factor activities. *Nucleic Acids
  Research* 51, 10934-10949 (2023). https://doi.org/10.1093/nar/gkad841
- Ahlmann-Eltze C, Huber W, Anders S. Deep-learning-based gene perturbation effect
  prediction does not yet outperform simple linear baselines. *Nature Methods* 22,
  1657 (2025). https://doi.org/10.1038/s41592-025-02772-6
- Kernfeld EM, Yang Y, Weinstock JS, Battle A, Cahan P. A systematic comparison of
  computational methods for expression forecasting. *Genome Biology* (2025).
  https://doi.org/10.1186/s13059-025-03840-y
- Viñas Torné R, et al. Systema: a framework for evaluating genetic perturbation
  response prediction beyond systematic variation. *Nature Biotechnology* (2025).
  https://doi.org/10.1038/s41587-025-02777-8
