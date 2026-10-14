#!/usr/bin/env bash
# Run the whole pipeline from a fresh checkout: download data, compute the expected spread datasets, make the figures.
# Needs about 20 GB of disk for data/ and about 4 GB of RAM. Usage: ./run_all.sh
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=${PYTHON:-python}

$PY scripts/01_download_data.py # a) download all required data
$PY scripts/02_prepare_scatter.py        #    cache the scattered profiles as one npz
$PY scripts/03_compute_expected_spread.py  # b) compute the expected std / RMSE datasets -> outputs/expected_std/*.nc
$PY scripts/04_make_plots.py             # c) overview figures -> reports/expected_std/*.png
