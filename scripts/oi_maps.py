"""Maps of the OI climatology (masked field + calibrated std) next to the published field.
Usage: PYTHONPATH=. python scripts/oi_maps.py 2 20 56
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cf

sys.path.insert(0, os.path.dirname(__file__))
from reproduce import published

OUT = "reports/oi/maps"; os.makedirs(OUT, exist_ok=True)
lat = np.linspace(60, 90, 121); lon = np.arange(-180, 180, .25); LON, LAT = np.meshgrid(lon, lat)
tr = ccrs.PlateCarree()
MONTHS = ((0, "Jan"), (2, "Mar"), (5, "Jun"), (8, "Sep"))


def panel(fig, pos, f, title, cmap, vlim, label):
    ax = fig.add_subplot(*pos, projection=ccrs.NorthPolarStereo()); ax.set_extent([-180, 180, 62, 90], tr)
    m = ax.pcolormesh(LON, LAT, f, transform=tr, cmap=cmap, vmin=vlim[0], vmax=vlim[1], shading="auto", rasterized=True)
    ax.add_feature(cf.LAND, facecolor=".85"); ax.coastlines(lw=.3); ax.set_title(title, fontsize=8)
    c = fig.colorbar(m, ax=ax, shrink=.7, pad=.02); c.set_label(label, fontsize=7); c.ax.tick_params(labelsize=7)


def main(levels):
    pub = published()
    for k in levels:
        d = np.load(f"outputs/oi/level_{k}.npz")
        for v, name in enumerate(("T", "S")):
            fig = plt.figure(figsize=(12, 15))
            for r, (mo, mn) in enumerate(MONTHS):
                p, f, s = pub[mo, v, k], d[f"{name}_field"][mo], d[f"{name}_std"][mo]
                vl = np.nanpercentile(np.concatenate([p[np.isfinite(p)], f[np.isfinite(f)]]), [1, 99])
                panel(fig, (4, 3, 3 * r + 1), p, f"{name} {mn} L{k}: published", "viridis", vl, name)
                panel(fig, (4, 3, 3 * r + 2), f, f"{name} {mn} L{k}: OI (NaN = no data)", "viridis", vl, name)
                panel(fig, (4, 3, 3 * r + 3), s, f"{name} {mn} L{k}: OI calibrated std", "magma", (0, np.nanpercentile(s, 99)), name)
            fig.tight_layout(); fig.savefig(f"{OUT}/oi_{name}_L{k}.png", dpi=80); plt.close(fig)


if __name__ == "__main__":
    main([int(x) for x in sys.argv[1:]])
