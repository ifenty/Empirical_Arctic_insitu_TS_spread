"""Compare the published AOTD monthly climatology with WOA23 (decav91C0, 0.25 deg, 57 levels, objectively analysed t_an/s_an).

The published grid is labelled 60.0..90.0 N / -180..179.75 E, but it is taken to BE the WOA23 grid (cell centres 59.875..89.875 N,
-179.875..179.875 E, WOA rows 599..719): the ocean mask agrees best at zero shift (level 0: 2.4% mismatch, worse at every other
integer row/column shift). So the comparison is cell for cell, no regridding. WOA23 is 1991-2020, AOTD is 1983-2023.

Outputs: outputs/woa/woa_on_aotd_grid.npy (untracked), reports/woa_compare/{tables/level_stats.csv, figures/*.png, figures/maps/*.png}
Usage (from the repo root): python scripts/06_report_aotd_vs_woa23.py
"""
import os, sys
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
for d in (OUT, REP, f"{REP}/figures/maps", f"{REP}/tables"):
    os.makedirs(d, exist_ok=True)
NL = 57
R0 = 599  # WOA23 lat row of 59.875 N; rows 599..719 cover 59.875..89.875


def build():
    f = f"{OUT}/woa_on_aotd_grid.npy"
    if os.path.exists(f):
        return np.load(f, mmap_mode="r")
    arr = np.lib.format.open_memmap(f, mode="w+", dtype="float32", shape=(12, 2, NL, 121, 1440))
    for v, (letter, var) in enumerate((("t", "t_an"), ("s", "s_an"))):
        for m in range(12):
            ds = xr.open_dataset(f"{RAW}/woa23_decav91C0_{letter}{m + 1:02d}_04.nc", decode_times=False)
            a = ds[var].isel(time=0, lat=slice(R0, 720)).values
            arr[m, v] = a
            ds.close()
    arr.flush()
    return np.load(f, mmap_mode="r")


def stats(woa, pub):
    rows = []
    depth = xr.open_dataset(f"{RAW}/woa23_decav91C0_t01_04.nc", decode_times=False).depth.values
    for k in range(NL):
        for v, name in enumerate(("T", "S")):
            for m in range(12):
                a, b = np.asarray(pub[m, v, k]), np.asarray(woa[m, v, k])
                both = np.isfinite(a) & np.isfinite(b)
                row = dict(level=k, depth=float(depth[k]), var=name, month=m + 1, n_aotd=int(np.isfinite(a).sum()),
                           n_woa=int(np.isfinite(b).sum()), n_both=int(both.sum()))
                if both.sum() > 10:
                    d = a[both] - b[both]
                    row.update(bias=float(d.mean()), rmse=float(np.sqrt((d**2).mean())), maxabs=float(np.abs(d).max()),
                               corr=float(np.corrcoef(a[both], b[both])[0, 1]))
                rows.append(row)
    return pd.DataFrame(rows), depth


def plots(df):
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    for v, name in enumerate(("T", "S")):
        s = df[df["var"] == name]
        g = s.groupby("depth")
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


def maps(woa, pub, depth, levels=(2, 20, 56), months=((0, "Jan"), (2, "Mar"), (5, "Jun"), (8, "Sep"))):
    lat = np.linspace(60, 90, 121); lon = np.arange(-180, 180, .25); LON, LAT = np.meshgrid(lon, lat)
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


if __name__ == "__main__":
    pub = published()
    woa = build()
    df, depth = stats(woa, pub)
    df.to_csv(f"{REP}/tables/level_stats.csv", index=False)
    plots(df)
    maps(woa, pub, depth)
    print(df.groupby("var")[["rmse", "bias", "corr"]].median())
