# Expected spread of Arctic in situ temperature and salinity

**Goal:** estimate how much individual in situ profile values of potential temperature and salinity are expected to scatter around reference monthly climatologies (WOA23 and the published AOTD climatology) in the Arctic Ocean, as a function of location (EASE-Grid 2.0 North, 100 km) and depth (the 102 WOA standard levels). The estimates are data-honest (NaN where there is no observation), are provided as expected **RMSE** and expected **standard deviation**, and are computed against three references: WOA23, the published AOTD climatology, and the larger of the two.

The method, data, provenance and limitations are described in **[`docs/technical_note.md`](docs/technical_note.md)**.

## Data and credit
The scattered profiles are the AOTD *scatter* dataset provided by the authors of the Arctic Ocean Thermohaline Dataset:

> Li, J., Wu, X. & Wang, X. An Arctic Ocean Thermohaline Dataset. *Scientific Data* 12, 1607 (2025). https://doi.org/10.1038/s41597-025-05855-3
> Li, J. AOTD_scatter. figshare (2025). https://doi.org/10.6084/m9.figshare.28825490.v1 (CC BY 4.0)

Reference climatologies: WOA23 decadal (1991-2020) 0.25° monthly and seasonal T and S, NOAA NCEI; and the published AOTD monthly climatology, https://doi.org/10.6084/m9.figshare.27645738.v1.

## Reproduce
Requirements: Python >= 3.10, about 20 GB free disk (WOA23 files are large), about 4 GB RAM, internet access.

```bash
git clone https://github.com/ifenty/AOTD_v2.git && cd AOTD_v2
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
./run_all.sh            # runs scripts 01-05 below (download takes the longest)
```

Or step by step (always from the repo root):

| script | command | what it does | output |
|---|---|---|---|
| 01 | `python scripts/01_download_data.py` | downloads AOTD scatter + climatology (figshare) and WOA23 monthly/seasonal T, S (NCEI); resumable | `data/` |
| 02 | `python scripts/02_prepare_scatter.py` | concatenates yearly scatter files, normalises missing values | `data/scatter_all.npz` |
| 03 | `python scripts/03_compute_expected_spread.py` | deviations from the references in potential temperature and salinity, per-cell std/RMSE/bias/n, 3 products (about 1.5 min) | `outputs/expected_spread/*.nc` (12 files) |
| 04 | `python scripts/04_make_plots.py` | maps at 10/100/1500 m, log10 profile counts, median spread and counts vs depth, for each statistic and product | `reports/expected_spread/figures/*.png` (11 files) |
| 05 | `python scripts/05_make_report.py` | writes the report with tables and captions | `reports/expected_spread/expected_spread_report.md` |

Optional diagnostics (not needed for the products):
* `python scripts/06_compare_aotd_woa23.py` compares the published AOTD climatology with WOA23 (-> `reports/woa_compare/`).
* `python scripts/07_scatter_data_report.py` describes the scattered data, including heuristically inferred platform classes (-> `reports/scatter/`). It also needs the NSIDC EASE-Grid 2.0 N 25 km ancillary files in `data/ease2/` (NSIDC-0772 lat/lon and NSIDC-0609 LOCI land mask; these require a free NASA Earthdata login, so they are not downloaded automatically).

Tests: `pytest tests`.

## Output files
`outputs/expected_spread/expected_{std,rmse}_from_{woa23,aotd,woa23_and_aotd}_{potential_temperature,salinity}_ease2_100km.nc`; dimensions `(depth 102, y 180, x 180)`, variables `Theta_*` or `S_*` for `std`/`rmse`, plus `*_bias`, `*_n`, 2-D `lat`/`lon` and the EPSG:6931 `crs`. NaN means no observation. See the technical note for definitions.

## Layout
- `scripts/01..07_*.py` pipeline and diagnostics (above); `run_all.sh`
- `scripts/scatter_io.py` (loading the scatter files) and `scripts/potential_temperature.py` (UNESCO 1983 potential temperature): helper modules imported by the scripts
- `docs/technical_note.md` technical note; `reports/` figures and reports; `tests/`
- `data/` (incl. `data/aotd_climatology/`, `data/woa23_raw/`, `data/scatter/`) and `outputs/` are generated and untracked

## History
This repository began as an attempt to reproduce the AOTD climatology from the authors' released MATLAB code. That code proved incomplete (only two of twelve published monthly files could be reproduced) and a replacement optimal-interpolation climatology did not improve on WOA23, so the project was refocused on the in situ spread. The reproduction and OI code was removed; it is available in the git history before the commit "Refocus repository on in situ spread estimation".
