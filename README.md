# Empirical Arctic in situ temperature and salinity profile spread relative to climatologies

This repository summarizes Arctic in situ temperature and salinity differences from
two published reference climatologies: WOA23 and AOTD. It produces per-cell **empirical
RMSE, population standard deviation, mean offset and comparison count**, on a 100 km
EASE-Grid 2.0 North grid (EPSG:6931) at 102 WOA standard depths.

The analysis pools profiles from 1992 onward, comparing each value with its monthly
reference (seasonal WOA23 below 1500 m). Missing comparisons remain NaN; no spatial
gap filling is applied. A third product takes the larger finite metric from the two
references. This maximum is a sensitivity envelope, **not a calibrated confidence
bound**. The spreads include natural variability and reference mismatch as well as
observational uncertainty; they are not independently validated observation-error or
forecast-uncertainty estimates.

## Documentation and results

* [Technical note](technical_note.md): data, formulas, assumptions and limitations.
* [Empirical-spread report](reports/expected_spread/expected_spread_report.md):
  depth summaries, maps, offsets and sampling counts.
* [Climatology comparison](reports/woa_compare/woa_compare_report.md):
  published AOTD versus WOA23, including grid-alignment diagnostics.
* [Scatter-data report](reports/scatter/scatter_report.md):
  profile coverage and exploratory, unvalidated sampling classifications.

Reports distinguish measured differences from possible explanations. AOTD scatter
and climatology share input data; overlap with WOA23 is not quantified here, so the
comparisons are not independent holdout validation.

## Data and attribution

The AOTD authors provided the scattered profiles and published climatology:

* Li, J., Wu, X. & Wang, X. (2025). *An Arctic Ocean Thermohaline Dataset*.
  *Scientific Data* 12, 1607.
  [doi:10.1038/s41597-025-05855-3](https://doi.org/10.1038/s41597-025-05855-3).
* Li, J. (2025). *AOTD_scatter*, figshare, version 1 (CC BY 4.0).
  [doi:10.6084/m9.figshare.28825490.v1](https://doi.org/10.6084/m9.figshare.28825490.v1).
* Li, J. (2024). *An Arctic Ocean Thermohaline Dataset*, published climatology,
  figshare, version 1.
  [doi:10.6084/m9.figshare.27645738.v1](https://doi.org/10.6084/m9.figshare.27645738.v1).

WOA23 reference fields are NOAA NCEI's 1991–2020 (`decav91C0`) 0.25° objectively
analysed monthly and seasonal temperature and salinity:
[World Ocean Atlas](https://www.ncei.noaa.gov/products/world-ocean-atlas).
Full scientific references are in the [technical note](technical_note.md#references).

## Reproduction

Run commands from the repository root. Use Python ≥ 3.10 and all packages in
[`requirements.txt`](requirements.txt). Initial downloads require internet access.
Plan for **at least 30 GB of free disk space and several gigabytes of RAM**, with
more space if retaining additional archives or caches. Runtime, peak memory and
storage use depend on the machine and available inputs.

### 1. Install dependencies

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
```

Cartopy and netCDF4 require compiled dependencies; if suitable wheels are unavailable,
use a conda-forge environment or install the appropriate system libraries.
`gsw` provides the potential-temperature test reference; `pytest` runs the tests.

### 2. Run the pipeline

```bash
PYTHON=.venv/bin/python ./run_all.sh
```

The runner executes scripts 01–06. Script 07 also runs when the required NSIDC files
are available in `data/ease2/`; those ancillary files are not downloaded automatically.
Downloads can be rerun to reuse verified local files.

To run individual stages, use `.venv/bin/python scripts/<filename>`:

| stage | filename | purpose | output |
|---|---|---|---|
| 01 | `01_download_data.py` | Download AOTD scatter/climatology and WOA23 monthly/seasonal T and S | `data/` |
| 02 | `02_prepare_scatter.py` | Concatenate annual profiles and normalize missing values | `data/scatter_all.npz` |
| 03 | `03_compute_expected_spread.py` | Compute potential-temperature and salinity differences, then cellwise statistics | `outputs/expected_spread/` (12 NetCDF files) |
| 04 | `04_make_plots.py` | Map metrics, mean offsets and counts; summarize depth dependence | `reports/expected_spread/figures/` |
| 05 | `05_report_expected_spread.py` | Regenerate the spread report and summary CSV from stage 03 products | `reports/expected_spread/` |
| 06 | `06_report_aotd_vs_woa23.py` | Compare reference climatologies; regenerate report, tables and maps | `reports/woa_compare/`, cache in `outputs/woa/` |
| 07 | `07_report_scatter_data.py` | Describe profile sampling and exploratory classifications | `reports/scatter/` |

Stage 07 additionally needs NSIDC EASE-Grid 2.0 North 25 km latitude/longitude files
(NSIDC-0772) and a LOCI land mask (NSIDC-0609) in `data/ease2/`. Access requires a
free NASA Earthdata account. The expected filenames are
`NSIDC0772_LatLon_EASE2_N25km_v1.1.nc` and
`EASE2_N25km.LOCImask_land50_coast0km.720x720.bin`.
Stages 05 and 06 can be rerun without recomputing the
spread products, provided their inputs exist; stage 05 also links figures from 04.

```bash
.venv/bin/python -m pytest tests
```

## Product format

Files in `outputs/expected_spread/` follow this pattern:

```text
expected_{std,rmse}_from_{woa23,aotd,woa23_and_aotd}_{potential_temperature,salinity}_ease2_100km.nc
```

Dimensions are `(depth=102, y=180, x=180)`. Variables are `Theta_std`/`Theta_rmse`
or `S_std`/`S_rmse`, plus the associated `*_bias`, `*_n`, geographic coordinates and
EPSG:6931 CRS metadata. AOTD comparisons exist only at 0–1500 m; below that, the
combined product equals WOA23. Its bias and count fields follow the reference
selected for **that file's metric**; std and RMSE can select different references.
Counts are not effective independent sample sizes, and a singleton has zero
population std. See the [technical note](technical_note.md#3-calculation) before use.

## Repository layout and scope

`scripts/01..07_*.py` implement the pipeline and diagnostics. `scatter_io.py` handles
profile loading; `potential_temperature.py` implements the UNESCO surface-referenced
temperature conversion. `reports/` contains report text, tables and figures;
`tests/` contains validation checks. Downloaded `data/` and computed `outputs/` are
generated locally and are not tracked.

Earlier climatology-reproduction and optimal-interpolation experiments remain in
the git history. The current project focuses on empirical profile–reference spread,
not reconstruction of the published AOTD analysis.
