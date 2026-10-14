"""One overview figure per product (from_woa23, from_aotd, from_woa23_and_aotd): std maps, log10 profile-count maps, and depth profiles.
Usage: PYTHONPATH=. python scripts/expected_std_overview.py
"""
import numpy as np
import xarray as xr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cf

D = "outputs/expected_std"
LEVELS = (2, 20, 56)
tr = ccrs.PlateCarree()
def make(prod):
    global dep
    ds = {v: xr.open_dataset(f"{D}/expected_std_from_{prod}_{nm}_ease2_100km.nc") for v, nm in (("T", "temperature"), ("S", "salinity"))}
    dep = ds["T"].depth.values
    lat, lon = ds["T"].lat.values, ds["T"].lon.values
    same_n = bool(np.array_equal(ds["T"]["T_n"].values, ds["S"]["S_n"].values))

    fig = plt.figure(figsize=(14, 15))
    rows = [("T_std", "T", "magma_r", None), ("S_std", "S", "magma_r", None), ("T_n", "T", "viridis", "log10 N")]
    if not same_n:
        rows.append(("S_n", "S", "viridis", "log10 N"))
    nr = len(rows) + 1
    for r, (var, v, cmap, lab) in enumerate(rows):
        a = ds[v][var].values.astype("float64")
        for c, k in enumerate(LEVELS):
            f = a[k]
            if lab:
                f = np.where(f > 0, np.log10(np.where(f > 0, f, 1)), np.nan)
                vmin, vmax = 0, np.nanmax(np.log10(np.maximum(a[a > 0], 1))) if (a > 0).any() else 1
                vmax = np.nanmax(f) if np.isfinite(f).any() else 1
            else:
                vmin, vmax = 0, np.nanpercentile(f, 98) if np.isfinite(f).any() else 1
            ax = fig.add_subplot(nr, 3, 3 * r + c + 1, projection=ccrs.NorthPolarStereo()); ax.set_extent([-180, 180, 62, 90], tr)
            m = ax.pcolormesh(lon, lat, f, transform=tr, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto", rasterized=True)
            ax.add_feature(cf.LAND, facecolor=".85"); ax.coastlines(lw=.3)
            ax.set_title(f"[{prod}] {var if not lab else 'log10 N (' + v + ')'} at {dep[k]:.0f} m", fontsize=9)
            cb = fig.colorbar(m, ax=ax, shrink=.7); cb.ax.tick_params(labelsize=7)
    ax1 = fig.add_subplot(nr, 2, 2 * nr - 1)
    import warnings; warnings.simplefilter("ignore")
    for v in "TS":
        s = ds[v][f"{v}_std"].values.reshape(len(dep), -1)
        ax1.plot(dep, np.nanmedian(s, axis=1), label=f"{v}_std, median over cells")
    ax1.set_xscale("symlog", linthresh=50); ax1.set_xlabel("depth (m)"); ax1.legend(fontsize=8); ax1.grid(alpha=.3)
    ax2 = fig.add_subplot(nr, 2, 2 * nr)
    n = ds["T"]["T_n"].values.reshape(len(dep), -1)
    ax2.plot(dep, np.maximum(n.sum(1), 1), label="profiles (all cells)")
    ax2.plot(dep, np.maximum((n > 0).sum(1), 1), label="cells with data")
    ax2.set_yscale("log"); ax2.set_xscale("symlog", linthresh=50); ax2.set_xlim(0, 5500); ax1.set_xlim(0, 5500); ax2.set_xlabel("depth (m)"); ax2.legend(fontsize=8); ax2.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(f"reports/expected_std/expected_std_overview_from_{prod}.png", dpi=90); plt.close(fig)
    print("T and S counts identical:", same_n)


for prod in ("woa23", "aotd", "woa23_and_aotd"):
    make(prod)
