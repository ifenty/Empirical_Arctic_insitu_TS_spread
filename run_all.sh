#!/usr/bin/env bash
# Run the whole pipeline from a fresh checkout: download data, compute the expected-spread datasets, make figures and reports.
# Needs about 20 GB of disk for data/ and about 4 GB of RAM. Usage: [PYTHON=.venv/bin/python] ./run_all.sh
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=${PYTHON:-python}

"$PY" scripts/01_download_data.py            # download scatter, AOTD climatology, WOA23 -> data/
"$PY" scripts/02_prepare_scatter.py          # cache the scattered profiles -> data/scatter_all.npz
"$PY" scripts/03_compute_expected_spread.py  # expected std / RMSE / bias / n -> outputs/expected_spread/*.nc
"$PY" scripts/04_make_plots.py               # figures -> reports/expected_spread/figures/
"$PY" scripts/05_report_expected_spread.py   # expected-spread report -> reports/expected_spread/
"$PY" scripts/06_report_aotd_vs_woa23.py      # AOTD vs WOA23 report -> reports/woa_compare/

# Scatter-data report: needs the NSIDC EASE2 25 km files (Earthdata login), which cannot be downloaded automatically.
if [ -f data/ease2/NSIDC0772_LatLon_EASE2_N25km_v1.1.nc ] &&
   [ -f data/ease2/EASE2_N25km.LOCImask_land50_coast0km.720x720.bin ]; then
    "$PY" scripts/07_report_scatter_data.py  # -> reports/scatter/
else
    echo "Skipping scripts/07_report_scatter_data.py: required EASE2 latitude/longitude or land-mask file is missing (see README)."
fi
