# Expected spread of Arctic in situ temperature and salinity about climatology

*Technical note. Code and instructions: see the repository `README.md`. All numbers below are produced by `scripts/03_compute_expected_spread.py` and can be regenerated with `./run_all.sh`.*

## 1. Purpose and philosophy

Ocean models and state estimates that assimilate in situ profiles need to know how far a single observation is expected to differ from the model's (or a climatology's) mean state, so that observations can be weighted sensibly. This note describes a **data-driven estimate of that spread for the Arctic Ocean**: for each location and depth, how much do individual in situ profile values scatter around a reference (monthly) climatology? Two reference climatologies are used: the World Ocean Atlas 2023 (WOA23) and the published AOTD climatology, and a third estimate combines them.

Design principles:

1. **Data-honest.** The estimate is built only from observations. Where a cell has no observation at a given depth the value is **NaN**; nothing is extrapolated, smoothed in from neighbours, or filled from a background field.
2. **Fixed, published reference climatologies.** The spread is measured about an existing climatology rather than a field we constructed from the same observations. The first estimate uses WOA23 (1991-2020, monthly to 1500 m and seasonal below), the second the published AOTD monthly climatology (0-1500 m), and the third takes the larger of the two (a conservative envelope). The method itself is reference-agnostic: any monthly climatology on the same grid can be substituted.
3. **Standard grids.** Output is on the WOA standard depth levels (102 levels, 0-5500 m) and the 100 km EASE-Grid 2.0 North grid (EPSG:6931), a regular equal-area grid suited to Arctic modelling. There is no time dimension: the estimate pools all years and seasons, but each observation is compared with the climatology of its own month (or season at depth).
4. **Potential temperature.** Model temperature state variables are potential temperature, so the deviations are formed in potential temperature, not in situ temperature.
5. **Two complementary measures** are given: the *expected RMSE* (spread about the reference, including any systematic offset) and the *expected standard deviation* (spread about the cell's own mean offset, i.e. with the bias removed).

History of the project: it started as an attempt to reproduce the AOTD climatology (Li et al., 2025) from the authors' released MATLAB code. The released code turned out to be incomplete: only two of the twelve published monthly files (March, September) could be reproduced. An optimal-interpolation alternative was prototyped but was not an improvement over WOA23, and the project was redirected to the spread estimates described here. The reproduction and OI code were removed from this repository (they remain in the git history).

## 2. Data

### 2.1 Scattered in situ profiles (AOTD scatter dataset)

* **Source:** the AOTD *scatter* dataset, 414,221 quality-controlled temperature and salinity profiles for 1983-2023 north of 60°N, interpolated to the 102 WOA standard depth levels (0-5500 m), one NetCDF file per year. Li, J. (2025), *AOTD_scatter*, figshare, https://doi.org/10.6084/m9.figshare.28825490.v1 (CC BY 4.0). We thank the AOTD authors (J. Li, X. Wu and X. Wang) for making these data available. Paper: Li, J., Wu, X. & Wang, X. (2025), An Arctic Ocean Thermohaline Dataset, *Scientific Data* 12, 1607, https://doi.org/10.1038/s41597-025-05855-3.
* **Variables used:** latitude, longitude, time, temperature `temp` (treated as in situ, ITS-90), practical salinity `salt`. The files also carry `pres` and `dens`, which are not used.
* **Cleaning on load** (`aotd/io.py`): in the files only `temp` marks missing values with NaN; salinity and pressure use exactly 0.0 and density holds arbitrary values at missing levels. All are converted to NaN, with salinity/pressure/density invalidated wherever temperature is missing or salinity is 0. The level at 0 m is empty in every profile; the 5 m level is only 67% filled.
* **Used here:** profiles from 1992 onward (400,281 of 414,221) that fall inside the 100 km grid and the WOA23 domain (60-90°N): 400,280 profiles.

**Provenance (from the paper, plus our best guesses).** The paper reports that 1,453,232 profiles were collected from several programmes and databases (WOD, CMEMS, Argo, Ice-Tethered Profilers, Bolin Centre and NABOS expedition data, and Chinese Arctic Research Expeditions on Xuelong / Xuelong-2), then de-duplicated and passed through a unified QC sequence: position/time checks, depth-reversal and repeatability checks, temperature/salinity range checks, T/S gradient checks, a density-gradient check, interpolation to the WOA standard levels with the Barker & McDougall (2020) interpolation, a second range/gradient check after interpolation, and a 3σ screening. 414,221 profiles survive. The released files carry **no source or platform identifier**, so the mix of sources is not documented. From the sampling pattern (track linking on time, position and depth reach; `scripts/06_scatter_data_report.py`, `reports/scatter/`) our *best guess* is: ~29% ITP-like drifting profilers in the ice-covered basin, ~6% Argo-like 10-day-cycle floats, ~52% ship CTD-like (stations, moorings, gliders and other), and ~14% unclassified. These classes are heuristic and unvalidated; near-duplicates and repeated profiles may remain. Because many profiles come from a few drifting platforms, **profiles are not independent samples** and the profile counts overstate the number of independent observations.

### 2.2 Reference climatologies

* **WOA23** (World Ocean Atlas 2023, decadal 1991-2020, `decav91C0`, 0.25° grid), objectively analysed fields `t_an` (in situ temperature) and `s_an` (practical salinity), from NOAA NCEI: 12 monthly files (57 levels, 0-1500 m) and 4 seasonal files (102 levels; winter = Jan-Mar, spring = Apr-Jun, summer = Jul-Sep, autumn = Oct-Dec). References: Locarnini et al. (2024) and Reagan et al. (2024), NOAA Atlas NESDIS 89 and 90. https://www.ncei.noaa.gov/products/world-ocean-atlas
* **AOTD climatology**: the published monthly climatology of the AOTD authors (57 levels, 0-1500 m, 0.25°), https://doi.org/10.6084/m9.figshare.27645738.v1. Its grid is labelled 60.0-90.0°N, -180 to 179.75°E, but its ocean mask agrees best with the WOA23 grid shifted by half a cell (cell centres 59.875-89.875°N), the labels are probably a notational slip, and we treat it as the WOA23 grid (WOA rows 599-719); this is an assumption (`reports/woa_compare/`). The published monthly files other than March and September have a deep salinity anomaly (RMSE ~0.6 psu relative to WOA23 below 800 m), so the AOTD-based estimate is the less trustworthy of the two references.

## 3. Calculation

For every profile value (profile *p*, standard level *k*) from 1992-2023:

1. **Location and cell.** The profile position is projected to EASE-Grid 2.0 North (EPSG:6931) and assigned to a 100 km × 100 km cell (180 × 180 grid). It is also assigned to the WOA23 0.25° cell containing it (nearest-cell lookup, no interpolation).
2. **Reference.** For levels 0-1500 m (indices 0-56) the reference is the WOA23 *monthly* field of the profile's calendar month; for deeper levels it is the WOA23 *seasonal* field of its season. (WOA23 monthly files stop at 1500 m; the switch is at 1500 m.) For the AOTD reference only the monthly field to 1500 m exists, so that product is NaN below 1500 m.
3. **Potential temperature.** Observed and reference in situ temperature are both converted to potential temperature referenced to the surface (UNESCO 1983 Bryden adiabatic lapse rate with the fourth-order Runge-Kutta integration used in the standard routine; pressure from depth and latitude by Saunders & Fofonoff). Salinity is the practical salinity of the observation, or, where the observed salinity is missing but temperature is valid, the reference salinity (used only inside the conversion, never stored). The code agrees with `gsw.pt_from_t` to about 1e-3 °C. The difference is
   `d_theta = theta(T_obs, S_obs, z, lat) - theta(T_ref, S_ref, z, lat)`.
   Salinity: `d_S = S_obs - S_ref`.
   The conversion is done per value, before differencing, so that T-S correlated errors are retained. (Scaling the in situ T statistic instead would give nearly the same answer: median differences are below 0.002 °C.)
4. **Accumulation per cell and level.** For each 100 km cell and level we accumulate `n = count(d)`, `Σd`, `Σd²` over all profiles, months and years.
5. **Statistics.**
   * `bias = Σd / n`
   * `RMSE = sqrt(Σd² / n)` — spread about the reference, bias included.
   * `std = sqrt(Σd²/n − bias²)` — spread about the cell's mean offset (population form; 0 if n = 1).
   * NaN wherever `n = 0` (no observation); no minimum `n` is imposed, so users should look at `n`.
6. **Three products** per variable and statistic: `from_woa23`; `from_aotd`; `from_woa23_and_aotd` = elementwise `max` of the two fields (ignoring NaN, so below 1500 m it equals the WOA23-based value). The bias and n written with the combined product belong to whichever reference gave the larger value.

### Output files (`outputs/expected_std/`, NetCDF4, CF-1.8)

`expected_{std,rmse}_from_{woa23,aotd,woa23_and_aotd}_{potential_temperature,salinity}_ease2_100km.nc`, 12 files. Dimensions `(depth=102, y=180, x=180)`; variables `Theta_std|Theta_rmse` (°C) or `S_std|S_rmse` (practical salinity), `*_bias`, `*_n` (number of values), plus 2-D `lat`, `lon` and a `crs` grid-mapping variable (EPSG:6931, EASE-Grid 2.0 North 100 km).

## 4. Results in brief

Medians over cells with data (WOA23-based, then AOTD-based, then combined):

| | level | WOA23 | AOTD | max |
|---|---|---|---|---|
| Theta RMSE (°C) | 10 m / 100 m / 1500 m | 0.84 / 0.51 / 0.076 | 0.76 / 0.43 / 0.050 | 0.92 / 0.53 / 0.077 |
| Theta std (°C) | 10 m / 100 m / 1500 m | 0.71 / 0.43 / 0.047 | 0.69 / 0.40 / 0.041 | 0.75 / 0.45 / 0.050 |
| S RMSE | 10 m / 100 m / 1500 m | 0.60 / 0.135 / 0.0078 | 0.43 / 0.103 / 0.0057 | 0.64 / 0.144 / 0.0081 |
| S std | 10 m / 100 m / 1500 m | 0.44 / 0.107 / 0.0059 | 0.36 / 0.091 / 0.0048 | 0.47 / 0.112 / 0.0060 |
| cells with data | 10 m / 100 m / 1500 m | 1666 / 1379 / 579 | 1672 / 1366 / 572 | 1684 / 1388 / 583 |

The spread is largest in the upper ocean over the shelves and the Atlantic-side seas and smallest in the ice-covered central basin and at depth. RMSE is 5-35% larger than std, i.e. the cell-mean offset from WOA23 is not negligible compared with the spread (largest for temperature at depth and salinity at 10-50 m). Data are densest (hundreds to thousands of values per cell) on the Atlantic side and in the Beaufort Sea; at 1500 m, 31% of cells with data have fewer than 5 values. Figures: `reports/expected_std/expected_{std,rmse}_overview_from_*.png` (std/RMSE maps at 10, 100 and 1500 m, log10 of the number of values, median spread and counts versus depth).

## 5. Interpretation and limitations

* The quantity is a **total** spread: it contains instrument and QC error, representation error (profiles are points, WOA23 is a 0.25° objective analysis), real synoptic and mesoscale variability, interannual variability and trends (1992-2023 pooled against a 1991-2020 climatology), and the error of the climatology itself. It is not a pure observation error.
* The reference is not unbiased with respect to these observations (WOA23 may not include them, uses different data, and is smoothed), hence the non-zero bias; this is why both RMSE and std are provided.
* Pooling over years and months: each value is compared with its own month/season, but the statistic is pooled, so it mixes seasons within a season-scale reference below 1500 m and years of different state.
* Counts are not independent: ITP/Argo profiles along a drift track are strongly correlated. Cells with small `n` (std is exactly 0 for n = 1) are noisy; no thresholding is applied.
* WOA23 is used at the nearest 0.25° cell; no spatial or temporal interpolation of the reference is performed.
* The AOTD-based product depends on an assumed grid alignment and on a published climatology whose ten non-March/September months have a deep-salinity anomaly; the combined product inherits any resulting inflation of the spread in the top 1500 m.
* Source mix (Section 2.1) is only a best guess. QC of the AOTD profiles was not re-done; this note takes the published quality control as given.

## Appendix A: code excerpts

Full code: `scripts/03_compute_expected_spread.py`, `aotd/potential_temperature.py`.

**A.1 Cell indices (EASE-Grid 2.0 North 100 km and WOA23 0.25°).**

```python
CELL, HALF, NC = 100e3, 9000e3, 180
x, y = Transformer.from_crs(4326, 6931, always_xy=True).transform(lon, lat)
ix = np.floor((x + HALF) / CELL).astype(int)
iy = np.floor((HALF - y) / CELL).astype(int)          # y runs downward
cell = iy * NC + ix
row = np.floor((lat + 90.0) / 0.25).astype(int) - R0   # R0 = 599: WOA23 row of 59.875 N
col = np.floor((lon + 180.0) / 0.25).astype(int) % 1440
```

**A.2 Reference field for each profile (WOA23).**

```python
NMONTHLY = 57                    # levels 0-1500 m have monthly fields
SEASON_FILE = np.array([13, 13, 13, 14, 14, 14, 15, 15, 15, 16, 16, 16])   # WOA23 seasonal files
Tr[:NMONTHLY] = monT[m - 1][:, r, c]                       # monthly, 0-1500 m
Tr[NMONTHLY:] = sesT[SEASON_FILE[m - 1] - 13][:, r, c]     # seasonal, deeper
```

**A.3 Potential temperature of the observation and of the reference, and the deviation.**

```python
S_for_obs = np.where(np.isnan(s_o), Sr, s_o)     # reference S only inside the conversion, only where S_obs is missing
d_theta = calculate_potential_T(t_o, S_for_obs, z, la) - calculate_potential_T(Tr, Sr, z, la)
d_S = s_o - Sr
```

```python
# aotd/potential_temperature.py (excerpt): pressure from depth and latitude, UNESCO 1983 / Saunders & Fofonoff
def depth_to_pressure(depth, lat):
    magic_number = 5.92e-3 + np.sin(np.abs(lat) * np.pi / 180) ** 2 * 5.25e-3
    return ((1 - magic_number) - np.sqrt(((1 - magic_number) ** 2) - (8.84e-6 * depth))) / 4.42e-6
# calculate_potential_T(T, S, depth, lat) then integrates the adiabatic lapse rate from P to 0 dbar with a 4th-order Runge-Kutta scheme.
```

**A.4 Accumulation of n, Σd and Σd² per 100 km cell.**

```python
good = np.isfinite(dev)
for k in range(nlev):
    g = good[k]
    n[k]  += np.bincount(ce[g], minlength=NC * NC)
    s1[k] += np.bincount(ce[g], weights=dev[k, g], minlength=NC * NC)
    s2[k] += np.bincount(ce[g], weights=dev[k, g] ** 2, minlength=NC * NC)
```

**A.5 Statistics and the combined product.**

```python
mean = s1 / n
rmse = np.where(n > 0, np.sqrt(s2 / n), np.nan)
std  = np.where(n > 0, np.sqrt(np.maximum(s2 / n - mean ** 2, 0)), np.nan)
bias = np.where(n > 0, mean, np.nan)
combined = np.fmax(from_woa23, from_aotd)          # NaN-ignoring elementwise max
```

**A.6 Reading the output.**

```python
import xarray as xr
ds = xr.open_dataset("outputs/expected_std/expected_rmse_from_woa23_potential_temperature_ease2_100km.nc")
ds.Theta_rmse.sel(depth=100, method="nearest").plot()   # 100 m map; NaN where there is no observation
```

## References

* Li, J., Wu, X. & Wang, X. An Arctic Ocean Thermohaline Dataset. *Sci. Data* 12, 1607 (2025). https://doi.org/10.1038/s41597-025-05855-3
* Li, J. AOTD_scatter. figshare (2025). https://doi.org/10.6084/m9.figshare.28825490.v1 ; Li, J. An Arctic Ocean Thermohaline Dataset (climatology). figshare (2024). https://doi.org/10.6084/m9.figshare.27645738.v1
* Locarnini, R. A. et al. World Ocean Atlas 2023, Volume 1: Temperature. NOAA Atlas NESDIS 89 (2024); Reagan, J. R. et al. World Ocean Atlas 2023, Volume 2: Salinity. NOAA Atlas NESDIS 90 (2024).
* UNESCO. Algorithms for computation of fundamental properties of seawater. UNESCO Technical Papers in Marine Science 44 (1983).
* Barker, P. M. & McDougall, T. J. Two interpolation methods using multiply-rotated piecewise cubic Hermite interpolating polynomials. *J. Atmos. Oceanic Technol.* 37 (2020).
* Brodzik, M. J. et al. EASE-Grid 2.0: incremental but significant improvements for Earth-gridded data sets. *ISPRS Int. J. Geo-Inf.* 1, 32-45 (2012).
