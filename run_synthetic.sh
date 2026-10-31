#!/usr/bin/env bash
# Positive control: simulate a screen with a known network and run the full benchmark.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data
python scripts/make_synthetic.py --out data
python scripts/run_benchmark.py --h5ad data/synthetic.h5ad --prior data/synthetic_prior.csv \
    --n-hvg 300 --out results_synthetic
python scripts/plot_results.py --results results_synthetic
