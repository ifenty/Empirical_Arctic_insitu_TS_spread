"""Figures of the expected-std and expected-RMSE datasets (potential temperature and salinity), written to reports/expected_spread/figures/:
  expected_{std,rmse}_overview_from_{product}.png   metric maps, log10 profile-count maps, median metric / counts vs depth
  spread_vs_depth.png, spread_vs_depth_full.png      median and 10-90% range of RMSE, std and RMSE/std versus depth, all products
  bias_maps_from_woa23.png                            mean offset (obs - WOA23) maps
  product_difference_maps.png                         AOTD-based minus WOA23-based RMSE maps
  counts_vs_depth.png                                 number of values, cells with data, fraction of sparse cells
Usage (from the repo root): python scripts/04_make_plots.py
"""
import os
import numpy as np
import xarray as xr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cf

D = "outputs/expected_spread"
LEVELS = (2, 20, 56)
tr = ccrs.PlateCarree()
VARS = (("Theta", "potential_temperature"), ("S", "salinity"))


def cell_percentiles(values, percentiles):
    values = values.reshape(len(values), -1)
    supported = np.isfinite(values).any(axis=1)
    result = np.full((len(percentiles), len(values)), np.nan)
    result[:, supported] = np.nanpercentile(values[supported], percentiles, axis=1)
    return result


def rmse_std_ratio(rmse, std, rmse_n, std_n):
    ratio = np.full_like(rmse, np.nan, dtype="float64")
    np.divide(rmse, std, out=ratio, where=(rmse_n >= 5) & (std_n >= 5) & (std > 0))
    return ratio


def make(kind, prod):
    ds = {v: load(kind, prod, v) for v, _ in VARS}
    dep = ds["Theta"].depth.values
    lat, lon = ds["Theta"].lat.values, ds["Theta"].lon.values
    same_n = bool(np.array_equal(ds["Theta"]["Theta_n"].values, ds["S"]["S_n"].values))
    rows = [(f"Theta_{kind}", "Theta", "magma_r", None), (f"S_{kind}", "S", "magma_r", None), ("Theta_n", "Theta", "viridis", "log10 N")]
    if not same_n:
        rows.append(("S_n", "S", "viridis", "log10 N"))
    nr = len(rows) + 1
    fig = plt.figure(figsize=(14, 15))
    grid = fig.add_gridspec(nr, 6)
    for r, (var, v, cmap, lab) in enumerate(rows):
        a = ds[v][var].values.astype("float64")
        for c, k in enumerate(LEVELS):
            f = a[k]
            if lab:
                f = np.where(f > 0, np.log10(np.where(f > 0, f, 1)), np.nan)
                vmin, vmax = 0, (np.nanmax(f) if np.isfinite(f).any() else 1)
            else:
                vmin, vmax = 0, (np.nanpercentile(f, 98) if np.isfinite(f).any() else 1)
            ax = fig.add_subplot(grid[r, 2 * c:2 * c + 2], projection=ccrs.NorthPolarStereo()); ax.set_extent([-180, 180, 62, 90], tr)
            m = ax.pcolormesh(lon, lat, f, transform=tr, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto", rasterized=True)
            ax.add_feature(cf.LAND, facecolor=".85"); ax.coastlines(lw=.3)
            ax.set_title(f"[{kind} from {prod}] {var if not lab else 'log10 N (' + v + ')'} at {dep[k]:.0f} m", fontsize=8)
            cb = fig.colorbar(m, ax=ax, shrink=.7); cb.ax.tick_params(labelsize=7)
    ax1 = fig.add_subplot(grid[-1, :3])
    for v, _ in VARS:
        s = ds[v][f"{v}_{kind}"].values.reshape(len(dep), -1)
        ax1.plot(dep, cell_percentiles(s, (50,))[0], label=f"{v}_{kind}, median over cells")
    ax1.set_xscale("symlog", linthresh=50); ax1.set_xlim(0, 5500); ax1.set_xlabel("depth (m)"); ax1.legend(fontsize=8); ax1.grid(alpha=.3)
    ax2 = fig.add_subplot(grid[-1, 3:])
    n = ds["Theta"]["Theta_n"].values.reshape(len(dep), -1)
    ax2.plot(dep, np.maximum(n.sum(1), 1), label="profiles (all cells)")
    ax2.plot(dep, np.maximum((n > 0).sum(1), 1), label="cells with data")
    ax2.set_yscale("log"); ax2.set_xscale("symlog", linthresh=50); ax2.set_xlim(0, 5500); ax2.set_xlabel("depth (m)"); ax2.legend(fontsize=8); ax2.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(f"reports/expected_spread/figures/expected_{kind}_overview_from_{prod}.png", dpi=90); plt.close(fig)


PRODS = ("woa23", "aotd", "woa23_and_aotd")
COL = {"woa23": "C0", "aotd": "C3", "woa23_and_aotd": "k"}
NAMES = {"Theta": "potential_temperature", "S": "salinity"}
UNIT = {"Theta": "degC", "S": "psu"}
OUTD = "reports/expected_spread/figures"


def load(kind, prod, v):
    with xr.open_dataset(f"{D}/expected_{kind}_from_{prod}_{NAMES[v]}_ease2_100km.nc") as ds:
        return ds.load()


def polar_axes(fig, nr, nc, i):
    ax = fig.add_subplot(nr, nc, i, projection=ccrs.NorthPolarStereo()); ax.set_extent([-180, 180, 62, 90], tr)
    ax.add_feature(cf.LAND, facecolor=".85"); ax.coastlines(lw=.3)
    return ax


def vs_depth(full):
    fig, ax = plt.subplots(2, 3, figsize=(15, 9), sharey=True)
    for r, v in enumerate(("Theta", "S")):
        data = {}
        for k in ("rmse", "std"):
            for p in (("woa23",) if full else PRODS):
                data[(k, p)] = load(k, p, v)[f"{v}_{k}"].values
        dep = load("rmse", "woa23", v).depth.values
        for c, k in enumerate(("rmse", "std", "ratio")):
            a = ax[r, c]
            for p in (("woa23",) if full else PRODS):
                if k == "ratio":
                    rmse_n = load("rmse", p, v)[f"{v}_n"].values
                    std_n = load("std", p, v)[f"{v}_n"].values
                    f = rmse_std_ratio(data[("rmse", p)], data[("std", p)], rmse_n, std_n)
                else:
                    f = data[(k, p)]
                lo, med, hi = cell_percentiles(f, (10, 50, 90))
                a.plot(med, dep, color=COL[p], label=f"from {p}" + (" (median)" if p == "woa23" else ""))
                if p == "woa23":
                    a.fill_betweenx(dep, lo, hi, color=COL[p], alpha=.15, label="10-90% of cells")
            a.set_ylim((5500 if full else 1500), 0)
            if full:
                a.set_yscale("symlog", linthresh=100); a.set_ylim(5500, 0)
            a.set_xlabel({"rmse": f"expected RMSE ({UNIT[v]})", "std": f"expected std ({UNIT[v]})", "ratio": "RMSE / std"}[k])
            a.set_title(f"{v}: " + {"rmse": "RMSE", "std": "std", "ratio": "RMSE / std (both counts >= 5)"}[k], fontsize=10); a.grid(alpha=.3)
            if k == "ratio":
                a.set_xlim(1, 3)
            else:
                a.set_xlim(left=0)
            a.set_ylabel("depth (m)")
    ax[0, 0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f"{OUTD}/spread_vs_depth{'_full' if full else ''}.png", dpi=90); plt.close(fig)


def bias_maps():
    fig = plt.figure(figsize=(14, 9))
    for r, v in enumerate(("Theta", "S")):
        ds = load("rmse", "woa23", v); b = ds[f"{v}_bias"].values; dep = ds.depth.values
        for c, k in enumerate(LEVELS):
            f = b[k]; lim = np.nanpercentile(np.abs(f), 98)
            ax = polar_axes(fig, 2, 3, 3 * r + c + 1)
            m = ax.pcolormesh(ds.lon.values, ds.lat.values, f, transform=tr, cmap="RdBu_r", vmin=-lim, vmax=lim, shading="auto", rasterized=True)
            ax.set_title(f"{v}_bias (obs - WOA23) at {dep[k]:.0f} m ({UNIT[v]})", fontsize=9)
            cb = fig.colorbar(m, ax=ax, shrink=.7); cb.ax.tick_params(labelsize=7)
    fig.tight_layout(); fig.savefig(f"{OUTD}/bias_maps_from_woa23.png", dpi=90); plt.close(fig)


def product_difference_maps():
    fig = plt.figure(figsize=(14, 9))
    for r, v in enumerate(("Theta", "S")):
        a, w = load("rmse", "aotd", v), load("rmse", "woa23", v); dep = a.depth.values
        for c, k in enumerate(LEVELS):
            f = a[f"{v}_rmse"].values[k] - w[f"{v}_rmse"].values[k]; lim = np.nanpercentile(np.abs(f), 98)
            ax = polar_axes(fig, 2, 3, 3 * r + c + 1)
            m = ax.pcolormesh(a.lon.values, a.lat.values, f, transform=tr, cmap="RdBu_r", vmin=-lim, vmax=lim, shading="auto", rasterized=True)
            ax.set_title(f"{v} RMSE: from AOTD - from WOA23, {dep[k]:.0f} m ({UNIT[v]})", fontsize=9)
            cb = fig.colorbar(m, ax=ax, shrink=.7); cb.ax.tick_params(labelsize=7)
    fig.tight_layout(); fig.savefig(f"{OUTD}/product_difference_maps.png", dpi=90); plt.close(fig)


def counts_vs_depth():
    ds = load("rmse", "woa23", "Theta"); dep = ds.depth.values
    n = ds["Theta_n"].values.reshape(len(dep), -1); has = n > 0
    neff = ds["Theta_n_eff"].values.reshape(len(dep), -1) if "Theta_n_eff" in ds else n
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    ax[0].plot(dep, np.maximum(n.sum(1), 1), label="all values (n)")
    ax[0].plot(dep, np.maximum(neff.sum(1), 1), "--", label="sampling days (n_eff)")
    ax[0].set_yscale("log"); ax[0].set_title("comparisons across all cells"); ax[0].legend(fontsize=8)
    ax[1].plot(dep, has.sum(1), label="n >= 1 (RMSE)")
    ax[1].plot(dep, (n >= 2).sum(1), "--", label="n >= 2 (sample std)")
    ax[1].set_title("cells with data (of 180x180)"); ax[1].legend(fontsize=8)
    for th in (5, 30):
        ax[2].plot(dep, 100 * ((n < th) & has).sum(1) / np.maximum(has.sum(1), 1), label=f"n < {th}")
    ax[2].plot(dep, 100 * ((neff < 5) & has).sum(1) / np.maximum(has.sum(1), 1), "--", label="n_eff < 5")
    ax[2].set_title("% of cells with data that are sparse"); ax[2].legend(fontsize=8)
    for a in ax:
        a.set_xscale("symlog", linthresh=100); a.set_xlim(0, 5500); a.set_xlabel("depth (m)"); a.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(f"{OUTD}/counts_vs_depth.png", dpi=90); plt.close(fig)


def main():
    os.makedirs(OUTD, exist_ok=True)
    for kind in ("std", "rmse"):
        for prod in PRODS:
            make(kind, prod)
    vs_depth(False); vs_depth(True); bias_maps(); product_difference_maps(); counts_vs_depth()
    print("figures written to", OUTD)


if __name__ == "__main__":
    main()
