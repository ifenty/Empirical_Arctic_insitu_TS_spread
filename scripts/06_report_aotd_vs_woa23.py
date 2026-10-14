"""Compare AOTD and WOA23 monthly climatologies under an adopted grid alignment.

AOTD array indices are assigned to WOA23 centres (rows 599:720), without regridding.
Mask agreement supports this assumption but does not prove the intended geolocation.
Run from the repository root to regenerate the report, statistics and figures.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cf



def published():
    """(12, 2, 57, 121, 1440) float32 published AOTD monthly T,S."""
    arr = np.full((12, 2, 57, 121, 1440), np.nan, dtype="float32")
    for m in range(12):
        with xr.open_dataset(f"data/aotd_climatology/climatology_{m + 1:02d}.nc") as ds:
            arr[m, 0] = ds["temp"].values
            arr[m, 1] = ds["salt"].values
    return arr


RAW, OUT, REP = "data/woa23_raw", "outputs/woa", "reports/woa_compare"
NL = 57
R0 = 599  # WOA23 lat row of 59.875 N; rows 599..719 cover 59.875..89.875


def build():
    f = f"{OUT}/woa_on_aotd_grid.npy"
    if Path(f).exists():
        return np.load(f, mmap_mode="r")
    arr = np.lib.format.open_memmap(f, mode="w+", dtype="float32", shape=(12, 2, NL, 121, 1440))
    for v, (letter, var) in enumerate((("t", "t_an"), ("s", "s_an"))):
        for m in range(12):
            with xr.open_dataset(f"{RAW}/woa23_decav91C0_{letter}{m + 1:02d}_04.nc", decode_times=False) as ds:
                arr[m, v] = ds[var].isel(time=0, lat=slice(R0, 720)).values
    arr.flush()
    return np.load(f, mmap_mode="r")


def stats(woa, pub, depth, lat=None):
    if lat is None:
        lat = -89.875 + np.arange(R0, 720) * 0.25
    w_lat = np.cos(np.deg2rad(lat)).astype("float64")[:, None]
    w_2d = np.broadcast_to(w_lat, (len(lat), woa.shape[-1]))
    rows = []
    for k in range(len(depth)):
        for v, name in enumerate(("T", "S")):
            for m in range(12):
                a, b = np.asarray(pub[m, v, k]), np.asarray(woa[m, v, k])
                both = np.isfinite(a) & np.isfinite(b)
                row = dict(level=k, depth=float(depth[k]), var=name, month=m + 1, n_aotd=int(np.isfinite(a).sum()),
                           n_woa=int(np.isfinite(b).sum()), n_both=int(both.sum()),
                           bias=np.nan, rmse=np.nan, bias_unweighted=np.nan, rmse_unweighted=np.nan,
                           maxabs=np.nan, corr=np.nan)
                if both.sum() > 10:
                    av, bv = a[both].astype("float64"), b[both].astype("float64")
                    wv = w_2d[both]
                    wsum = wv.sum()
                    d = av - bv
                    w_bias = float(np.sum(wv * d) / wsum)
                    w_rmse = float(np.sqrt(np.sum(wv * (d ** 2)) / wsum))
                    da = av - np.sum(wv * av) / wsum
                    db = bv - np.sum(wv * bv) / wsum
                    denom = np.sqrt(np.sum(wv * da ** 2) * np.sum(wv * db ** 2))
                    w_corr = float(np.sum(wv * da * db) / denom) if denom > 0 and np.ptp(av) > 0 and np.ptp(bv) > 0 else np.nan
                    row.update(
                        bias=w_bias,
                        rmse=w_rmse,
                        bias_unweighted=float(d.mean()),
                        rmse_unweighted=float(np.sqrt((d ** 2).mean())),
                        maxabs=float(np.abs(d).max()),
                        corr=w_corr,
                    )
                rows.append(row)
    return pd.DataFrame(rows)


def mask_diagnostics(woa, pub, depth, levels=(0, 20), max_shift=3):
    """Compare January masks, using interior rows and cyclic longitude for shift tests."""
    rows = []
    nrow = woa.shape[-2]
    if nrow <= 2 * max_shift:
        raise ValueError("Latitude grid is too small for the requested mask-shift test")
    interior = slice(max_shift, nrow - max_shift)
    for k in levels:
        for v, name in enumerate(("T", "S")):
            a, b = np.isfinite(pub[0, v, k]), np.isfinite(woa[0, v, k])
            trials = []
            for dr in range(-max_shift, max_shift + 1):
                for dc in range(-max_shift, max_shift + 1):
                    shifted = np.roll(b[max_shift + dr:nrow - max_shift + dr], dc, axis=1)
                    trials.append((float(np.mean(a[interior] != shifted)), dr, dc))
            best, dr, dc = min(trials)
            rows.append(dict(depth=float(depth[k]), var=name,
                             mismatch_pct=100 * float(np.mean(a != b)),
                             interior_zero_pct=100 * float(np.mean(a[interior] != b[interior])),
                             interior_best_pct=100 * best, row_offset=dr, column_roll=dc))
    return pd.DataFrame(rows)


def comparison_table(df):
    rows = [
        "| depth (m) | T RMSE Mar/Sep | T RMSE other 10 | T offset other 10 | S RMSE Mar/Sep | S RMSE other 10 |",
        "|---|---|---|---|---|---|",
    ]
    for d in (0, 10, 50, 100, 300, 500, 750, 800, 1000, 1500):
        q = df[df.depth == d]
        if q.empty:
            continue
        med = q.groupby(["var", q.month.isin((3, 9))])[["rmse", "bias"]].median()
        values = [med.loc[("T", True), "rmse"], med.loc[("T", False), "rmse"],
                  med.loc[("T", False), "bias"], med.loc[("S", True), "rmse"],
                  med.loc[("S", False), "rmse"]]
        cells = ["–" if not np.isfinite(v) else f"{v:+.2f}" if i == 2 else f"{v:.2f}"
                 for i, v in enumerate(values)]
        rows.append(f"| {d} | " + " | ".join(cells) + " |")
    return "\n".join(rows)


def report(df, alignment):
    mask_rows = [
        "| variable / depth | mask mismatch, full grid (%) | interior mismatch, zero shift (%) | best tested interior mismatch (%) | best row offset / column roll |",
        "|---|---|---|---|---|",
    ]
    for row in alignment.itertuples():
        mask_rows.append(f"| {row.var} / {row.depth:g} m | {row.mismatch_pct:.2f} | "
                         f"{row.interior_zero_pct:.2f} | {row.interior_best_pct:.2f} | "
                         f"{row.row_offset:+d} / {row.column_roll:+d} |")
    return f"""# Published AOTD and WOA23: monthly climatology comparison

Generated by `scripts/06_report_aotd_vs_woa23.py` from the local published climatology
files. See the [README](../../README.md) for reproduction, the
[technical note](../../technical_note.md) for assumptions, and the
[empirical-spread report](../expected_spread/expected_spread_report.md) for the
observation-minus-reference products.

## 1. Purpose and method

This is a **comparison of two analyses**, not an independent accuracy assessment.
Differences are AOTD minus WOA23. Temperature (T) here is treated as **in situ**
temperature, unlike the potential-temperature differences in the spread products.
Practical salinity (S) is conventionally labelled psu.

* WOA23: `decav91C0`, objectively analysed `t_an` and `s_an`, monthly, 0.25°,
  57 levels from 0 to 1500 m, reference period 1991–2020.
* AOTD: published monthly fields for the dataset spanning 1983–2023.
* Each monthly statistic uses cells finite in **both** fields at the selected depth.
  Because 0.25° latitude–longitude cells converge toward the pole, primary summary
  statistics (`rmse`, `bias`, `corr`) weight each cell by **cos(latitude)**; unweighted
  grid-index means (`rmse_unweighted`, `bias_unweighted`) are also recorded in the CSV.
  Masks may differ between months, variables and depths.
* Area-weighted RMSE = √(Σ w(AOTD − WOA23)² / Σ w); area-weighted offset (`bias` in the CSV) = Σ w(AOTD − WOA23) / Σ w, with w = cos(lat).
  Correlation is the cos(lat)-weighted spatial Pearson correlation on the same finite-cell intersection.
  Statistics require more than ten shared cells. These summaries are not sampling
  confidence intervals.

### Adopted grid alignment

The local AOTD files expose integer index coordinates rather than geographic cell
centres. For this analysis their 121 × 1440 arrays are assigned to WOA23 rows 599–719:
latitude centres 59.875–89.875°N and longitude centres −179.875–179.875°E.
The fields are then compared cell for cell, without regridding. Maps use these same
adopted coordinates.

The following **January mask diagnostic** compares finite-value masks at the adopted
alignment and tests integer row/column offsets from −3 to +3. Shift tests omit the
outer three latitude rows on each side and wrap longitude cyclically, so their
percentages are not directly comparable to the full-grid percentage.
Offsets describe sampling a displaced WOA row and rolling its columns, not a
geolocation correction.

{chr(10).join(mask_rows)}

Mask agreement supports the adopted correspondence, but cannot resolve a fractional
cell shift or establish the intended AOTD geolocation. No independent confirmation
of that alignment is provided here.

## 2. Monthly statistics by depth

Each entry below is the **median of monthly statistics** in its group, not an RMSE
pooled across months. March and September are shown separately to describe the
observed monthly contrast; this is a descriptive, post hoc grouping, not a test of
its cause. Masks are determined separately for every monthly comparison.
T quantities are in °C; S quantities are in psu.

{comparison_table(df)}

The unrounded values, per-field counts (`n_aotd`, `n_woa`) and intersection counts
(`n_both`) are in [the level-statistics CSV](tables/level_stats.csv).

## 3. Figures

**Figure 1. Monthly comparison statistics at 0–1500 m.** Rows: T and S.
Columns: difference RMSE, mean offset and valid-cell counts. Lines show January,
March, June and September; the dashed count curve is WOA23 in January.
Counts in the last column describe each field separately, whereas RMSE and offsets
use their intersection.

![Monthly climatology statistics versus depth](figures/aotd_vs_woa23_vs_depth.png)

**Figure 2. AOTD, WOA23 and their difference at selected depths.** Each figure
has rows for January, March, June and September. The first two columns share
a colour range within each row (pooled 1st–99th percentiles); difference ranges
are symmetric about zero (98th percentile of absolute differences). Scales can
differ between rows. All fields are plotted at the adopted WOA23 centres.

| depth | temperature | salinity |
|---|---|---|
| 10 m | ![T at 10 m](figures/maps/aotd_vs_woa23_T_L2_10m.png) | ![S at 10 m](figures/maps/aotd_vs_woa23_S_L2_10m.png) |
| 100 m | ![T at 100 m](figures/maps/aotd_vs_woa23_T_L20_100m.png) | ![S at 100 m](figures/maps/aotd_vs_woa23_S_L20_100m.png) |
| 1500 m | ![T at 1500 m](figures/maps/aotd_vs_woa23_T_L56_1500m.png) | ![S at 1500 m](figures/maps/aotd_vs_woa23_S_L56_1500m.png) |

## 4. Findings and interpretation

1. **Differences depend strongly on depth and month.** The table shows larger
   near-surface RMSEs than at depth for temperature. Surface temperature offsets
   also vary by monthly group. A pooled group median should not be interpreted
   as an offset present in every month or location.
2. **Deep salinity merits targeted investigation.** In the ten-month group,
   salinity RMSE increases sharply between 750 and 800 m and remains much larger
   at 1000 and 1500 m than in March/September. This is a documented disagreement
   between products. Its cause and the relative accuracy of the two fields have
   not been established by this comparison.
3. **Consequences for the spread products are conditional.** Reference differences
   can change the empirical observation-minus-reference metrics where profiles
   sample the affected cells. A large climatology-to-climatology RMSE does not
   directly determine the magnitude of the profile-based spread or its combined
   maximum.

Potential explanations include different input sampling, reference periods, analysis
methods and grid interpretation. These are hypotheses, not attributions established
here. WOA23 is an analysed reference rather than truth. AOTD scatter and climatology
share source observations, and the extent of overlap with WOA23 is not quantified;
neither comparison constitutes an independent holdout validation.

## 5. Reproduction and attribution

Run `python scripts/06_report_aotd_vs_woa23.py` from the repository root after
downloading the reference files. It regenerates this report, the
[level statistics](tables/level_stats.csv),
[mask diagnostics](tables/mask_alignment.csv), and linked figures. The extracted
WOA23 cache is `outputs/woa/woa_on_aotd_grid.npy`; remove it if changing the source
files or extraction convention.

AOTD: Li, J., Wu, X. & Wang, X. (2025), *An Arctic Ocean Thermohaline Dataset*,
[doi:10.1038/s41597-025-05855-3](https://doi.org/10.1038/s41597-025-05855-3);
published climatology,
[doi:10.6084/m9.figshare.27645738.v1](https://doi.org/10.6084/m9.figshare.27645738.v1).
WOA23: Locarnini et al. (2024), *Temperature*, NOAA Atlas NESDIS 89; Reagan et al.
(2024), *Salinity*, NOAA Atlas NESDIS 90,
[NOAA NCEI](https://www.ncei.noaa.gov/products/world-ocean-atlas).
"""


def plots(df):
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    for v, name in enumerate(("T", "S")):
        s = df[df["var"] == name]
        for j, (col, lab) in enumerate((("rmse", "RMSE AOTD-WOA23"), ("bias", "bias AOTD-WOA23"), ("n_aotd", "valid cells"))):
            a = ax[v, j]
            for m in (1, 3, 6, 9):
                q = s[s.month == m]
                a.plot(q[col], q.depth, label=f"month {m}")
            if col == "n_aotd":
                a.plot(s[s.month == 1].n_woa, s[s.month == 1].depth, "k--", label="WOA23 (Jan)")
            a.set_ylim(1500, 0); a.set_title(f"{name} {lab}"); a.set_ylabel("depth (m)"); a.grid(alpha=.3)
        ax[v, 0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{REP}/figures/aotd_vs_woa23_vs_depth.png", dpi=100); plt.close(fig)


def maps(woa, pub, depth, lat, lon, levels=(2, 20, 56), months=((0, "Jan"), (2, "Mar"), (5, "Jun"), (8, "Sep"))):
    LON, LAT = np.meshgrid(lon, lat)
    tr = ccrs.PlateCarree()
    for k in levels:
        for v, name in enumerate(("T", "S")):
            fig = plt.figure(figsize=(12, 15))
            for r, (m, mn) in enumerate(months):
                a, b = np.asarray(pub[m, v, k]), np.asarray(woa[m, v, k])
                allv = np.concatenate([a[np.isfinite(a)], b[np.isfinite(b)]])
                vl = np.percentile(allv, [1, 99]) if len(allv) else (0, 1)
                dd = a - b
                dl = np.nanpercentile(np.abs(dd), 98) if np.isfinite(dd).any() else 1
                for c, (f, t, cm, lim) in enumerate(((a, "AOTD published", "viridis", vl), (b, "WOA23", "viridis", vl),
                                                     (dd, "AOTD - WOA23", "RdBu_r", (-dl, dl)))):
                    ax = fig.add_subplot(4, 3, 3 * r + c + 1, projection=ccrs.NorthPolarStereo()); ax.set_extent([-180, 180, 62, 90], tr)
                    h = ax.pcolormesh(LON, LAT, f, transform=tr, cmap=cm, vmin=lim[0], vmax=lim[1], shading="auto", rasterized=True)
                    ax.add_feature(cf.LAND, facecolor=".85"); ax.coastlines(lw=.3); ax.set_title(f"{name} {mn} {depth[k]:.0f} m: {t}", fontsize=8)
                    cb = fig.colorbar(h, ax=ax, shrink=.7, pad=.02); cb.ax.tick_params(labelsize=7)
            fig.tight_layout(); fig.savefig(f"{REP}/figures/maps/aotd_vs_woa23_{name}_L{k}_{depth[k]:.0f}m.png", dpi=80); plt.close(fig)


def main():
    for d in (OUT, REP, f"{REP}/figures/maps", f"{REP}/tables"):
        Path(d).mkdir(parents=True, exist_ok=True)
    pub = published()
    woa = build()
    with xr.open_dataset(f"{RAW}/woa23_decav91C0_t01_04.nc", decode_times=False) as ds:
        depth = ds.depth.values
        lat, lon = ds.lat.values[R0:720], ds.lon.values
    df = stats(woa, pub, depth, lat=lat)
    alignment = mask_diagnostics(woa, pub, depth)
    df.to_csv(f"{REP}/tables/level_stats.csv", index=False)
    alignment.to_csv(f"{REP}/tables/mask_alignment.csv", index=False)
    plots(df)
    maps(woa, pub, depth, lat, lon)
    path = Path(REP) / "woa_compare_report.md"
    path.write_text(report(df, alignment), encoding="utf-8")
    print("wrote", path)
    print(df.groupby("var")[["rmse", "bias", "corr"]].median())


if __name__ == "__main__":
    main()
