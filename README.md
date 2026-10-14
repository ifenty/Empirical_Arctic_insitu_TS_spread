# Expected spread of Arctic in situ temperature and salinity

**Goal:** estimate how much individual in situ profile values of potential temperature and salinity are expected to scatter around reference monthly climatologies (WOA23 and the published AOTD climatology) in the Arctic Ocean, as a function of location (EASE-Grid 2.0 North, 100 km) and depth (the 102 WOA standard levels). The estimates are data-honest (NaN where there is no observation), are provided as expected **RMSE** and expected **standard deviation**, and are computed against three references: WOA23, the published AOTD climatology, and the larger of the two.

The method, data, provenance and limitations are described in **[`technical_note.md`](technical_note.md)**.

## Data and credit
The scattered profiles are the AOTD *scatter* dataset provided by the authors of the Arctic Ocean Thermohaline Dataset:

> Li, J., Wu, X. & Wang, X. An Arctic Ocean Thermohaline Dataset. *Scientific Data* 12, 1607 (2025). https://doi.org/10.1038/s41597-025-05855-3
> Li, J. AOTD_scatter. figshare (2025). https://doi.org/10.6084/m9.figshare.28825490.v1 (CC BY 4.0)

Reference climatologies: WOA23 decadal (1991-2020) 0.25° monthly and seasonal T and S, NOAA NCEI; and the published AOTD monthly climatology, https://doi.org/10.6084/m9.figshare.27645738.v1.

## Reproduce
Requirements: internet access, about 20 GB free disk (WOA23 files are large), about 4 GB RAM, and **a Python (>= 3.10) environment with every package listed in `requirements.txt`** (numpy, scipy, pandas, xarray, netCDF4, pyproj, matplotlib, cartopy; gsw and pytest are only for the tests). The scripts will fail with `ImportError` if any is missing, so create the environment first.

### 1. Create the Python environment
Option A: venv + pip (works for people and automated agents; run from the repo root):
```bash
git clone https://github.com/ifenty/AOTD_v2.git && cd AOTD_v2
python3 -m venv .venv
source .venv/bin/activate          # every new shell needs this, or call .venv/bin/python directly
pip install --upgrade pip
pip install -r requirements.txt
python -c "import numpy, scipy, pandas, xarray, netCDF4, pyproj, matplotlib, cartopy; print('environment OK')"
```
Option B: conda/mamba:
```bash
conda create -n aotd_spread -c conda-forge python=3.12 --file requirements.txt
conda activate aotd_spread
```
(Cartopy and netCDF4 need compiled libraries; if `pip install` fails on Linux, use Option B or install the system packages `libgeos-dev libproj-dev libnetcdf-dev`.)
Automated agents: no shell state persists between commands, so prefix commands with `.venv/bin/python` (or set `PYTHON=.venv/bin/python` for `run_all.sh`) instead of relying on `activate`.

### 2. Run
```bash
./run_all.sh            # runs scripts 01-06 (and 07 if data/ease2 exists); the download takes the longest; uses $PYTHON, default `python`
```

Or step by step (always from the repo root):

| script | command | what it does | output |
|---|---|---|---|
| 01 | `python scripts/01_download_data.py` | downloads AOTD scatter + climatology (figshare) and WOA23 monthly/seasonal T, S (NCEI); resumable | `data/` |
| 02 | `python scripts/02_prepare_scatter.py` | concatenates yearly scatter files, normalises missing values | `data/scatter_all.npz` |
| 03 | `python scripts/03_compute_expected_spread.py` | deviations from the references in potential temperature and salinity, per-cell std/RMSE/bias/n, 3 products (about 1.5 min) | `outputs/expected_spread/*.nc` (12 files) |
| 04 | `python scripts/04_make_plots.py` | maps at 10/100/1500 m, log10 profile counts, median spread and counts vs depth, for each statistic and product | `reports/expected_spread/figures/*.png` (11 files) |
| 05 | `python scripts/05_report_expected_spread.py` | writes the expected-spread report (tables and captioned figures from 03/04) | `reports/expected_spread/expected_spread_report.md` |

Further reports (not needed for the products; 06 is run by `run_all.sh`, 07 only if `data/ease2/` exists):

| script | what it does | output |
|---|---|---|
| `python scripts/06_report_aotd_vs_woa23.py` | report comparing the published AOTD climatology with WOA23, the two references (RMSE, bias, maps) | `reports/woa_compare/woa_compare_report.md` |
| `python scripts/07_report_scatter_data.py` | report describing the scattered profiles (coverage, time, platform classes inferred heuristically). Needs the NSIDC EASE-Grid 2.0 N 25 km files in `data/ease2/` (NSIDC-0772 lat/lon, NSIDC-0609 LOCI land mask; free NASA Earthdata login, not downloaded automatically) | `reports/scatter/scatter_report.md` |

Each `reports/<name>/` folder has the report `.md`, `figures/` and `tables/`.

Tests: `pytest tests`.

## Output files
`outputs/expected_spread/expected_{std,rmse}_from_{woa23,aotd,woa23_and_aotd}_{potential_temperature,salinity}_ease2_100km.nc`; dimensions `(depth 102, y 180, x 180)`, variables `Theta_*` or `S_*` for `std`/`rmse`, plus `*_bias`, `*_n`, 2-D `lat`/`lon` and the EPSG:6931 `crs`. NaN means no observation. See the technical note for definitions.

## Layout
- `scripts/01..07_*.py` pipeline and diagnostics (above); `run_all.sh`
- `scripts/scatter_io.py` (loading the scatter files) and `scripts/potential_temperature.py` (UNESCO 1983 potential temperature): helper modules imported by the scripts
- `technical_note.md` technical note; `reports/` figures and reports; `tests/`
- `data/` (incl. `data/aotd_climatology/`, `data/woa23_raw/`, `data/scatter/`) and `outputs/` are generated and untracked

## History
This repository began as an attempt to reproduce the AOTD climatology from the authors' released MATLAB code. That code proved incomplete (only two of twelve published monthly files could be reproduced) and a replacement optimal-interpolation climatology did not improve on WOA23, so the project was refocused on the in situ spread. The reproduction and OI code was removed; it is available in the git history before the commit "Refocus repository on in situ spread estimation".
