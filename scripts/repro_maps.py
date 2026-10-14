"""Maps comparing the faithful-port reproduction with the published monthly climatology.

Outputs (reports/reproduction/maps/):
  cmp_{T|S}_L{level}_{depth}m.png   rows = months (Jan, Mar, Jun, Sep); cols = published, port, port - published
  diff_all_months_{T|S}_L{level}_{depth}m.png   port - published for all 12 months
  diff_depth_{T|S}_{month}.png      port - published at several depths for one month
Usage: PYTHONPATH=. .venv/bin/python scripts/repro_maps.py   (single process, ~1 GB)
"""
import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

OUT = "reports/reproduction/maps"
os.makedirs(OUT, exist_ok=True)
pub = np.load("outputs/repro/pub_month.npy", mmap_mode="r")
depth = np.load("data/scatter_all.npz")["depth"]
lat = np.linspace(60, 90, 121)
lon = np.arange(-180, 180, 0.25)
LON, LAT = np.meshgrid(lon, lat)
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
VARS = {0: ("T", "degC", "RdYlBu_r", 1.0), 1: ("S", "psu", "viridis", 0.5)}  # name, unit, cmap, diff color limit
proj = ccrs.NorthPolarStereo()
tr = ccrs.PlateCarree()
_cache = {}


def port(k):
    if k not in _cache:
        _cache[k] = np.load(f"outputs/repro/month_{k}.npy")
    return _cache[k]


def get(k, m, v):
    return port(k)[m, v], np.array(pub[m, v, k])


def panel(ax, field, cmap, vmin, vmax, title):
    ax.set_extent([-180, 180, 62, 90], tr)
    ax.set_boundary(plt.matplotlib.path.Path(np.c_[np.cos(np.linspace(0, 2 * np.pi, 100)) * .5 + .5,
                                                   np.sin(np.linspace(0, 2 * np.pi, 100)) * .5 + .5]), transform=ax.transAxes)
    im = ax.pcolormesh(LON, LAT, field, transform=tr, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto", rasterized=True)
    ax.add_feature(cfeature.LAND, facecolor="0.85", zorder=2)
    ax.coastlines(linewidth=0.3, zorder=3)
    ax.set_title(title, fontsize=8)
    return im


def lim(a):
    a = a[~np.isnan(a)]
    return np.percentile(a, 2), np.percentile(a, 98)


def cmp_figure(k, v, months=(0, 2, 5, 8)):
    name, unit, cmap, dl = VARS[v]
    fig = plt.figure(figsize=(10, 3.2 * len(months)))
    for r, m in enumerate(months):
        a, b = get(k, m, v)
        vmin, vmax = lim(b)
        d = a - b
        dd = max(dl, 0.0) if False else np.nanpercentile(np.abs(d), 98) or dl
        for c, (f, t, cm, lo, hi) in enumerate(((b, "published", cmap, vmin, vmax), (a, "port", cmap, vmin, vmax),
                                                 (d, "port - published", "RdBu_r", -dd, dd))):
            ax = fig.add_subplot(len(months), 3, r * 3 + c + 1, projection=proj)
            im = panel(ax, f, cm, lo, hi, f"{name} {MONTHS[m]} {depth[k]:.0f} m: {t}")
            fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02, label=unit if c < 2 else f"{unit} (98th pct of |diff|)")
    fig.tight_layout()
    fig.savefig(f"{OUT}/cmp_{name}_L{k}_{depth[k]:.0f}m.png", dpi=110)
    plt.close(fig)


def diff_all_months(k, v):
    name, unit, _, dl = VARS[v]
    ds = [get(k, m, v) for m in range(12)]
    dd = np.nanpercentile(np.abs(np.concatenate([(a - b).ravel() for a, b in ds])), 98) or dl
    fig = plt.figure(figsize=(14, 10.5))
    for m in range(12):
        ax = fig.add_subplot(3, 4, m + 1, projection=proj)
        im = panel(ax, ds[m][0] - ds[m][1], "RdBu_r", -dd, dd, f"{MONTHS[m]}")
    fig.suptitle(f"{name} port - published, {depth[k]:.0f} m (level {k})", fontsize=11)
    fig.colorbar(im, ax=fig.axes, shrink=0.5, pad=0.02, label=unit)
    fig.savefig(f"{OUT}/diff_all_months_{name}_L{k}_{depth[k]:.0f}m.png", dpi=100)
    plt.close(fig)


def diff_depth(v, m, levels=(2, 10, 20, 30, 40, 50, 56)):
    name, unit, _, dl = VARS[v]
    fig = plt.figure(figsize=(16, 7))
    for i, k in enumerate(levels):
        a, b = get(k, m, v)
        d = a - b
        dd = np.nanpercentile(np.abs(d), 98) or dl
        ax = fig.add_subplot(2, 4, i + 1, projection=proj)
        im = panel(ax, d, "RdBu_r", -dd, dd, f"{name} {MONTHS[m]} {depth[k]:.0f} m (limit {dd:.2f} {unit})")
        fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    fig.savefig(f"{OUT}/diff_depth_{name}_{MONTHS[m]}.png", dpi=100)
    plt.close(fig)


if __name__ == "__main__":
    for v in (0, 1):
        for k in (2, 20, 40, 56):
            cmp_figure(k, v)
        for k in (2, 20, 56):
            diff_all_months(k, v)
        for m in (0, 2, 5, 8):
            diff_depth(v, m)
    print("done")
