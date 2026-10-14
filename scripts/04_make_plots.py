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
    nep = ds["Theta_n_epochs"].values.reshape(len(dep), -1) if "Theta_n_epochs" in ds else neff
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    ax[0].plot(dep, np.maximum(n.sum(1), 1), label="all values (n)")
    ax[0].plot(dep, np.maximum(neff.sum(1), 1), "--", label="sampling days (n_eff)")
    ax[0].plot(dep, np.maximum(nep.sum(1), 1), ":", label="1-mo epochs (n_epochs)")
    ax[0].set_yscale("log"); ax[0].set_title("comparisons across all cells"); ax[0].legend(fontsize=8)
    ax[1].plot(dep, has.sum(1), label="n >= 1 (RMSE & SE)")
    ax[1].plot(dep, (n >= 2).sum(1), "--", label="n >= 2 (raw std)")
    ax[1].plot(dep, (nep >= 2).sum(1), ":", label="n_epochs >= 2 (declust std)")
    ax[1].set_title("cells with data (of 180x180)"); ax[1].legend(fontsize=8)
    for th in (5, 30):
        ax[2].plot(dep, 100 * ((n < th) & has).sum(1) / np.maximum(has.sum(1), 1), label=f"n < {th}")
    ax[2].plot(dep, 100 * ((nep < 5) & has).sum(1) / np.maximum(has.sum(1), 1), "--", label="n_epochs < 5")
    ax[2].set_title("% of cells with data that are sparse"); ax[2].legend(fontsize=8)
    for a in ax:
        a.set_xscale("symlog", linthresh=100); a.set_xlim(0, 5500); a.set_xlabel("depth (m)"); a.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(f"{OUTD}/counts_vs_depth.png", dpi=90); plt.close(fig)


def declustered_vs_depth():
    fig, ax = plt.subplots(2, 3, figsize=(15, 9), sharey=True)
    for r, v in enumerate(("Theta", "S")):
        ds_s = load("std", "woa23", v)
        ds_r = load("rmse", "woa23", v)
        dep = ds_s.depth.values
        # Col 0: std comparison (profile-weighted vs epoch-declustered vs EB 500km vs EB pan-Arctic)
        for name, col, ls, lab in ((f"{v}_std", "C0", "-", "profile-weighted std (n >= 2)"),
                                   (f"{v}_std_declustered", "C2", "-", "epoch-declustered std (M >= 2)"),
                                   (f"{v}_std_se", "C1", "-", "EB regularized (500 km prior, M >= 1)"),
                                   (f"{v}_std_se_pan_arctic", "C5", "--", "EB regularized (pan-Arctic prior)")):
            if name in ds_s:
                lo, med, hi = cell_percentiles(ds_s[name].values, (10, 50, 90))
                ax[r, 0].plot(med, dep, color=col, ls=ls, label=lab)
                if name == f"{v}_std_se":
                    ax[r, 0].fill_betweenx(dep, lo, hi, color=col, alpha=.15, label="EB (500 km) 10-90% of cells")
        ax[r, 0].set_title(f"{v}: profile-weighted, declustered & EB regularized std", fontsize=9.5)
        ax[r, 0].set_xlabel(f"std ({UNIT[v]})")
        # Col 1: variance components (within vs between vs total declustered)
        for name, col, lab in ((f"{v}_std_within", "C4", "within-epoch std (sub-monthly)"),
                               (f"{v}_std_between", "C3", "between-epoch std (synoptic/interannual)"),
                               (f"{v}_std_declustered", "C2", "epoch-declustered total std")):
            if name in ds_s:
                _, med, _ = cell_percentiles(ds_s[name].values, (10, 50, 90))
                ax[r, 1].plot(med, dep, color=col, label=lab)
        ax[r, 1].set_title(f"{v}: within- vs. between-epoch std (WOA23)", fontsize=9.5)
        ax[r, 1].set_xlabel(f"std component ({UNIT[v]})")
        # Col 2: RMSE comparison (profile-weighted vs epoch-declustered vs EB 500km vs EB pan-Arctic)
        for name, col, ls, lab in ((f"{v}_rmse", "C0", "-", "profile-weighted RMSE (n >= 1)"),
                                   (f"{v}_rmse_declustered", "C2", "-", "epoch-declustered RMSE"),
                                   (f"{v}_rmse_se", "C1", "-", "EB regularized RMSE (500 km prior)"),
                                   (f"{v}_rmse_se_pan_arctic", "C5", "--", "EB regularized RMSE (pan-Arctic prior)")):
            if name in ds_r:
                _, med, _ = cell_percentiles(ds_r[name].values, (10, 50, 90))
                ax[r, 2].plot(med, dep, color=col, ls=ls, label=lab)
        ax[r, 2].set_title(f"{v}: profile-weighted, declustered & EB regularized RMSE", fontsize=9.5)
        ax[r, 2].set_xlabel(f"RMSE ({UNIT[v]})")
        for c in range(3):
            ax[r, c].set_ylim(1500, 0); ax[r, c].set_xlim(left=0); ax[r, c].grid(alpha=.3); ax[r, c].set_ylabel("depth (m)")
    for c in range(3):
        ax[0, c].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{OUTD}/declustered_spread_vs_depth.png", dpi=90); plt.close(fig)


def declustered_maps():
    fig = plt.figure(figsize=(14, 15))
    for r_base, v in enumerate(("Theta", "S")):
        ds = load("std", "woa23", v)
        dep = ds.depth.values
        if f"{v}_std_se" not in ds:
            plt.close(fig)
            return
        for sub, (field_expr, cmap, sym, title_prefix) in enumerate((
            (ds[f"{v}_std_se"].values, "magma_r", False, f"{v}_std_se (EB regularized std)"),
            (ds[f"{v}_std_se"].values - ds[f"{v}_std"].values, "RdBu_r", True, f"{v}: EB regularized - profile-weighted std"),
        )):
            row_idx = 2 * r_base + sub
            for c, k in enumerate(LEVELS):
                f = field_expr[k]
                lim = np.nanpercentile(np.abs(f) if sym else f, 98) if np.isfinite(f).any() else 1.0
                vmin, vmax = (-lim, lim) if sym else (0.0, lim)
                ax = polar_axes(fig, 4, 3, 3 * row_idx + c + 1)
                m = ax.pcolormesh(ds.lon.values, ds.lat.values, f, transform=tr, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto", rasterized=True)
                ax.set_title(f"{title_prefix} at {dep[k]:.0f} m ({UNIT[v]})", fontsize=8.5)
                cb = fig.colorbar(m, ax=ax, shrink=.7); cb.ax.tick_params(labelsize=7)
    fig.tight_layout(); fig.savefig(f"{OUTD}/declustered_std_maps_from_woa23.png", dpi=90); plt.close(fig)


def method_comparison_maps(kind, v):
    ds = load(kind, "woa23", v)
    if f"{v}_{kind}_declustered" not in ds:
        return
    dep = ds.depth.values
    depths = (10, 50, 100, 300, 1000)
    k_list = [int(np.argmin(np.abs(dep - d))) for d in depths]
    raw_all = ds[f"{v}_{kind}"].values
    dec_all = ds[f"{v}_{kind}_declustered"].values
    se_all = ds[f"{v}_{kind}_se"].values
    fig = plt.figure(figsize=(15, 16))
    for r, k in enumerate(k_list):
        f_raw, f_dec, f_se = raw_all[k], dec_all[k], se_all[k]
        f_diff = f_se - f_raw
        pooled = np.concatenate([x[np.isfinite(x)] for x in (f_raw, f_dec, f_se) if np.isfinite(x).any()])
        vmax = float(np.percentile(pooled, 98)) if pooled.size else 1.0
        dlim = float(np.nanpercentile(np.abs(f_diff), 98)) if np.isfinite(f_diff).any() else 1.0
        panels = (
            (f_raw, "magma_r", (0.0, vmax), f"Profile-weighted {v}_{kind}"),
            (f_dec, "magma_r", (0.0, vmax), f"Epoch-declustered {v}_{kind}"),
            (f_se, "magma_r", (0.0, vmax), f"EB regularized {v}_{kind}_se"),
            (f_diff, "RdBu_r", (-dlim, dlim), f"EB regularized - profile-weighted"),
        )
        for c, (field, cmap, (vmin, vm), title) in enumerate(panels):
            ax = polar_axes(fig, len(k_list), 4, 4 * r + c + 1)
            m = ax.pcolormesh(ds.lon.values, ds.lat.values, field, transform=tr, cmap=cmap, vmin=vmin, vmax=vm, shading="auto", rasterized=True)
            ax.set_title(f"{dep[k]:.0f} m: {title} ({UNIT[v]})", fontsize=8)
            cb = fig.colorbar(m, ax=ax, shrink=.68, pad=.02); cb.ax.tick_params(labelsize=6.5)
    fig.tight_layout(); fig.savefig(f"{OUTD}/method_comparison_{kind}_{v}_from_woa23.png", dpi=90); plt.close(fig)


def variance_components_maps(v):
    ds = load("std", "woa23", v)
    if f"{v}_std_within" not in ds:
        return
    dep = ds.depth.values
    depths = (10, 50, 100, 300, 1000)
    k_list = [int(np.argmin(np.abs(dep - d))) for d in depths]
    w_all = ds[f"{v}_std_within"].values
    b_all = ds[f"{v}_std_between"].values
    t_all = ds[f"{v}_std_declustered"].values
    fig = plt.figure(figsize=(15, 16))
    for r, k in enumerate(k_list):
        fw, fb, ft = w_all[k], b_all[k], t_all[k]
        with np.errstate(invalid="ignore", divide="ignore"):
            frac_b = np.where(ft > 0, np.clip((fb ** 2) / (ft ** 2), 0.0, 1.0), np.nan)
        pooled = np.concatenate([x[np.isfinite(x)] for x in (fw, fb, ft) if np.isfinite(x).any()])
        vmax = float(np.percentile(pooled, 98)) if pooled.size else 1.0
        panels = (
            (fw, "magma_r", (0.0, vmax), f"{v}_std_within (<1 mo)", UNIT[v]),
            (fb, "magma_r", (0.0, vmax), f"{v}_std_between (>1 mo)", UNIT[v]),
            (ft, "magma_r", (0.0, vmax), f"Epoch-declustered {v}_std", UNIT[v]),
            (frac_b, "viridis", (0.0, 1.0), "between / total variance", "fraction"),
        )
        for c, (field, cmap, (vmin, vm), title, u_str) in enumerate(panels):
            ax = polar_axes(fig, len(k_list), 4, 4 * r + c + 1)
            m = ax.pcolormesh(ds.lon.values, ds.lat.values, field, transform=tr, cmap=cmap, vmin=vmin, vmax=vm, shading="auto", rasterized=True)
            ax.set_title(f"{dep[k]:.0f} m: {title} ({u_str})", fontsize=8)
            cb = fig.colorbar(m, ax=ax, shrink=.68, pad=.02); cb.ax.tick_params(labelsize=6.5)
    fig.tight_layout(); fig.savefig(f"{OUTD}/variance_components_{v}_from_woa23.png", dpi=90); plt.close(fig)


def method_distributions_by_depth():
    depths = (10, 50, 100, 300, 1000)
    fig, ax = plt.subplots(4, len(depths), figsize=(16, 12), sharey=False)
    rows_spec = (
        ("std", "Theta", "Theta std (°C)"),
        ("rmse", "Theta", "Theta RMSE (°C)"),
        ("std", "S", "S std (psu)"),
        ("rmse", "S", "S RMSE (psu)"),
    )
    for r, (kind, v, rlab) in enumerate(rows_spec):
        ds = load(kind, "woa23", v)
        if f"{v}_{kind}_declustered" not in ds:
            plt.close(fig)
            return
        dep = ds.depth.values
        raw_all = ds[f"{v}_{kind}"].values
        dec_all = ds[f"{v}_{kind}_declustered"].values
        se_all = ds[f"{v}_{kind}_se"].values
        pan_all = ds[f"{v}_{kind}_se_pan_arctic"].values if f"{v}_{kind}_se_pan_arctic" in ds else se_all
        for c, d in enumerate(depths):
            k = int(np.argmin(np.abs(dep - d)))
            a = ax[r, c]
            vals = [x[np.isfinite(x)] for x in (raw_all[k], dec_all[k], se_all[k], pan_all[k])]
            pooled = np.concatenate([x for x in vals if x.size])
            xmax = float(np.percentile(pooled, 98)) if pooled.size else 1.0
            grid_x = np.linspace(0.0, xmax, 200)
            for arr, col, ls, lab in zip(vals, ("C0", "C2", "C1", "C5"), ("-", "-", "-", "--"),
                                         ("Profile-wtd", "Epoch-declust", "EB (500km)", "EB (pan-Arc)")):
                if arr.size:
                    s_arr = np.sort(arr)
                    cdf = np.searchsorted(s_arr, grid_x, side="right") / float(arr.size)
                    a.plot(grid_x, cdf, color=col, ls=ls, lw=1.5, label=f"{lab} ({np.median(arr):.2g})")
            a.set_xlim(0.0, xmax); a.set_ylim(0.0, 1.02); a.grid(alpha=.3)
            if r == 0:
                a.set_title(f"{dep[k]:.0f} m", fontsize=10, fontweight="bold")
            if c == 0:
                a.set_ylabel(f"{rlab}\nEmpirical CDF", fontsize=8.5)
            a.set_xlabel(UNIT[v], fontsize=8)
            a.legend(fontsize=6, loc="lower right")
    fig.tight_layout(); fig.savefig(f"{OUTD}/method_distributions_by_depth.png", dpi=90); plt.close(fig)


def prior_comparison_maps():
    depths = (10, 50, 100, 300, 1000)
    for v in ("Theta", "S"):
        ds = load("std", "woa23", v)
        if f"{v}_std_se_pan_arctic" not in ds:
            return
        dep = ds.depth.values
        k_list = [int(np.argmin(np.abs(dep - d))) for d in depths]
        dec_all = ds[f"{v}_std_declustered"].values
        se500_all = ds[f"{v}_std_se"].values
        sepan_all = ds[f"{v}_std_se_pan_arctic"].values
        fig = plt.figure(figsize=(15, 16))
        for r, k in enumerate(k_list):
            f_dec, f_500, f_pan = dec_all[k], se500_all[k], sepan_all[k]
            f_diff = f_500 - f_pan
            pooled = np.concatenate([x[np.isfinite(x)] for x in (f_dec, f_500, f_pan) if np.isfinite(x).any()])
            vmax = float(np.percentile(pooled, 98)) if pooled.size else 1.0
            dlim = float(np.nanpercentile(np.abs(f_diff), 98)) if np.isfinite(f_diff).any() else 1.0
            panels = (
                (f_dec, "magma_r", (0.0, vmax), f"Epoch-declustered {v}_std"),
                (f_500, "magma_r", (0.0, vmax), f"EB (500 km prior) {v}_std_se"),
                (f_pan, "magma_r", (0.0, vmax), f"EB (pan-Arctic prior)"),
                (f_diff, "RdBu_r", (-dlim, dlim), f"500 km prior - pan-Arctic prior"),
            )
            for c, (field, cmap, (vmin, vm), title) in enumerate(panels):
                ax = polar_axes(fig, len(k_list), 4, 4 * r + c + 1)
                m = ax.pcolormesh(ds.lon.values, ds.lat.values, field, transform=tr, cmap=cmap, vmin=vmin, vmax=vm, shading="auto", rasterized=True)
                ax.set_title(f"{dep[k]:.0f} m: {title} ({UNIT[v]})", fontsize=8)
                cb = fig.colorbar(m, ax=ax, shrink=.68, pad=.02); cb.ax.tick_params(labelsize=6.5)
        fig.tight_layout(); fig.savefig(f"{OUTD}/prior_comparison_std_{v}_from_woa23.png", dpi=90); plt.close(fig)


def main():
    os.makedirs(OUTD, exist_ok=True)
    for kind in ("std", "rmse"):
        for prod in PRODS:
            make(kind, prod)
    vs_depth(False); vs_depth(True); bias_maps(); product_difference_maps(); counts_vs_depth()
    declustered_vs_depth(); declustered_maps()
    for kind in ("std", "rmse"):
        for v in ("Theta", "S"):
            method_comparison_maps(kind, v)
    for v in ("Theta", "S"):
        variance_components_maps(v)
    method_distributions_by_depth()
    prior_comparison_maps()
    print("figures written to", OUTD)


if __name__ == "__main__":
    main()
