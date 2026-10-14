"""Overview figures of the expected-std and expected-RMSE datasets (potential temperature and salinity): metric maps, log10 profile-count
maps, and median metric / counts vs depth. One figure per metric (std, rmse) and product (from_woa23, from_aotd, from_woa23_and_aotd).
Usage (from the repo root): python scripts/04_make_plots.py
"""
import warnings
import numpy as np
import xarray as xr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cf

warnings.simplefilter("ignore")
D = "outputs/expected_std"
LEVELS = (2, 20, 56)
tr = ccrs.PlateCarree()
VARS = (("Theta", "potential_temperature"), ("S", "salinity"))


def make(kind, prod):
    ds = {v: xr.open_dataset(f"{D}/expected_{kind}_from_{prod}_{nm}_ease2_100km.nc") for v, nm in VARS}
    dep = ds["Theta"].depth.values
    lat, lon = ds["Theta"].lat.values, ds["Theta"].lon.values
    same_n = bool(np.array_equal(ds["Theta"]["Theta_n"].values, ds["S"]["S_n"].values))
    rows = [(f"Theta_{kind}", "Theta", "magma_r", None), (f"S_{kind}", "S", "magma_r", None), ("Theta_n", "Theta", "viridis", "log10 N")]
    if not same_n:
        rows.append(("S_n", "S", "viridis", "log10 N"))
    nr = len(rows) + 1
    fig = plt.figure(figsize=(14, 15))
    for r, (var, v, cmap, lab) in enumerate(rows):
        a = ds[v][var].values.astype("float64")
        for c, k in enumerate(LEVELS):
            f = a[k]
            if lab:
                f = np.where(f > 0, np.log10(np.where(f > 0, f, 1)), np.nan)
                vmin, vmax = 0, (np.nanmax(f) if np.isfinite(f).any() else 1)
            else:
                vmin, vmax = 0, (np.nanpercentile(f, 98) if np.isfinite(f).any() else 1)
            ax = fig.add_subplot(nr, 3, 3 * r + c + 1, projection=ccrs.NorthPolarStereo()); ax.set_extent([-180, 180, 62, 90], tr)
            m = ax.pcolormesh(lon, lat, f, transform=tr, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto", rasterized=True)
            ax.add_feature(cf.LAND, facecolor=".85"); ax.coastlines(lw=.3)
            ax.set_title(f"[{kind} from {prod}] {var if not lab else 'log10 N (' + v + ')'} at {dep[k]:.0f} m", fontsize=8)
            cb = fig.colorbar(m, ax=ax, shrink=.7); cb.ax.tick_params(labelsize=7)
    ax1 = fig.add_subplot(nr, 2, 2 * nr - 1)
    for v, _ in VARS:
        s = ds[v][f"{v}_{kind}"].values.reshape(len(dep), -1)
        ax1.plot(dep, np.nanmedian(s, axis=1), label=f"{v}_{kind}, median over cells")
    ax1.set_xscale("symlog", linthresh=50); ax1.set_xlim(0, 5500); ax1.set_xlabel("depth (m)"); ax1.legend(fontsize=8); ax1.grid(alpha=.3)
    ax2 = fig.add_subplot(nr, 2, 2 * nr)
    n = ds["Theta"]["Theta_n"].values.reshape(len(dep), -1)
    ax2.plot(dep, np.maximum(n.sum(1), 1), label="profiles (all cells)")
    ax2.plot(dep, np.maximum((n > 0).sum(1), 1), label="cells with data")
    ax2.set_yscale("log"); ax2.set_xscale("symlog", linthresh=50); ax2.set_xlim(0, 5500); ax2.set_xlabel("depth (m)"); ax2.legend(fontsize=8); ax2.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(f"reports/expected_std/expected_{kind}_overview_from_{prod}.png", dpi=90); plt.close(fig)


for kind in ("std", "rmse"):
    for prod in ("woa23", "aotd", "woa23_and_aotd"):
        make(kind, prod)
