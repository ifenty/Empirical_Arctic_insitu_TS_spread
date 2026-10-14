#!/usr/bin/env python
"""Descriptive report of the AOTD scattered (profile) data.

Run from the repo root:  python scripts/07_report_scatter_data.py   (optional diagnostic; needs the NSIDC EASE2 25 km ancillary files in data/ease2, see README)
Writes ONLY to reports/scatter/ (md report, figures/*.png, tables/*.csv).
Diagnostic counts are computed here; source attribution is external to the arrays.
Sampling-pattern classes and coverage ratings are unvalidated heuristics.
"""
import datetime as dtm
import hashlib
import importlib.metadata
import os
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import numpy as np
import pandas as pd
import pyproj
import gsw
from matplotlib.colors import LogNorm, BoundaryNorm
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)
OUT = ROOT / "reports" / "scatter"
FIG = OUT / "figures"
TAB = OUT / "tables"
for p in (FIG, TAB):
    p.mkdir(parents=True, exist_ok=True)

T_START = time.time()
plt.rcParams.update({"figure.dpi": 100, "savefig.dpi": 110, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.3})


def log(msg):
    print(f"[{time.time() - T_START:6.0f}s] {msg}", flush=True)


# ----------------------------------------------------------------------------- helpers
SEC = {i: [] for i in range(1, 9)}  # markdown per section
K = {}  # key numbers for the executive summary / fitness table
TABLES_WRITTEN = []


def savefig(fig, name, caption=None):
    fig.savefig(FIG / f"{name}.png", bbox_inches="tight")
    plt.close(fig)
    return f"![{caption or name}](figures/{name}.png)\n\n*Figure: {caption or name}.*\n"


def savetab(df, name, index=False):
    df.to_csv(TAB / f"{name}.csv", index=index)
    TABLES_WRITTEN.append(name + ".csv")


def fmtv(v, dec):
    if v is None or pd.isna(v):
        return "–"
    if isinstance(v, (int, np.integer)):
        return f"{v:,}"
    if isinstance(v, (float, np.floating)):
        return f"{v:,.{dec}f}"
    return str(v)


def mdtable(df, dec=2, index=False):
    if index:
        df = df.reset_index()
    cols = list(df.columns)
    decs = dec if isinstance(dec, dict) else {}
    def escape(value):
        return str(value).replace("|", "\\|").replace("\n", "<br>")

    lines = ["| " + " | ".join(escape(c) for c in cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for r in df.itertuples(index=False, name=None):
        cells = []
        for c, value in zip(cols, r):
            d = decs.get(c, dec if not isinstance(dec, dict) else 2)
            formatted = (str(int(value)) if c in {"year", "first year", "year_first", "year_last"} and pd.notna(value)
                         else fmtv(value, d))
            cells.append(escape(formatted))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def pct(a, b=1.0, dec=1):
    return f"{100.0 * a / b:.{dec}f}%"


def add(sec, *txt):
    for t in txt:
        SEC[sec].append(t if t.endswith("\n") else t + "\n")


# ----------------------------------------------------------------------------- data
log("loading data")
NPZ = ROOT / "data" / "scatter_all.npz"
if NPZ.exists():
    INPUT_FILES = [NPZ]
    z = np.load(NPZ)
    lat = z["lat"].astype("float64")
    lon = z["lon"].astype("float64")
    tsec = z["time"].astype("int64")
    T = z["temp"].astype("float32")
    S = z["salt"].astype("float32")
    P = z["pres"].astype("float32")
    DN = z["dens"].astype("float32")
    depth = z["depth"].astype("float64")
    year = z["year"].astype(int)
    month = z["month"].astype(int)
    z.close()
else:
    from scatter_io import load_scatter

    d = load_scatter()
    lat, lon = d["lat"], d["lon"]
    tsec = d["time"].astype("datetime64[s]").astype("int64")
    T, S, P, DN, depth = d["temp"], d["salt"], d["pres"], d["dens"], d["depth"]
    year, month = d["year"], d["month"]
    INPUT_FILES = sorted((ROOT / "data" / "scatter").glob("AOTD*.nc"))
N, NZ = T.shape
tdays = tsec / 86400.0
vT = ~np.isnan(T)
vS = ~np.isnan(S)
nlevT = vT.sum(1)
nlevS = vS.sum(1)
hasT = nlevT > 0
hasS = nlevS > 0
lastT = NZ - 1 - np.argmax(vT[:, ::-1], axis=1)
firstT = np.argmax(vT, axis=1)
maxdepth = np.where(hasT, depth[lastT], np.nan)
firstdepth = np.where(hasT, depth[firstT], np.nan)
lastS = NZ - 1 - np.argmax(vS[:, ::-1], axis=1)
maxdepthS = np.where(hasS, depth[lastS], np.nan)
IDX1500 = int(np.where(depth == 1500)[0][0])  # last of the 57 WOA monthly levels
NL57 = IDX1500 + 1
GEN_DATE = dtm.date.today().isoformat()
DECADES = [(1983, 1992), (1993, 2002), (2003, 2012), (2013, 2023)]
DEC_LAB = [f"{a}-{b}" for a, b in DECADES]
dec_idx = np.full(N, -1)
for i, (a, b) in enumerate(DECADES):
    dec_idx[(year >= a) & (year <= b)] = i
SEASONS_DJF = {"DJF": (12, 1, 2), "MAM": (3, 4, 5), "JJA": (6, 7, 8), "SON": (9, 10, 11)}
SEASONS_WOA = {"JFM": (1, 2, 3), "AMJ": (4, 5, 6), "JAS": (7, 8, 9), "OND": (10, 11, 12)}
MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
lon180 = ((lon + 180.0) % 360.0) - 180.0  # lon == 180 maps to -180 (same meridian)
K.update(N=N, year_min=int(year.min()), year_max=int(year.max()))
log(f"loaded {N} profiles")

# ----------------------------------------------------------------------------- projection / boxes
# EASE-Grid 2.0 North (EPSG:6931, Lambert azimuthal equal area): 100 km boxes are exact 4x4 blocks of 25 km cells
import xarray as xr

PROJ_CRS = pyproj.CRS.from_epsg(6931)
TR = pyproj.Transformer.from_crs("EPSG:4326", PROJ_CRS, always_xy=True)
EDIR = ROOT / "data" / "ease2"
# land mask: uint8 720x720, rows top (y decreasing) to bottom. Verified against the shipped PNG and the profile positions:
# 255 = ocean/water, 0 = land, 101 = Greenland ice sheet / land-ice, 254 = off-earth.
LMASK = np.fromfile(EDIR / "EASE2_N25km.LOCImask_land50_coast0km.720x720.bin", dtype="uint8").reshape(720, 720)
_l = xr.open_dataset(EDIR / "NSIDC0772_LatLon_EASE2_N25km_v1.1.nc")
ELAT = _l["latitude"].values
ELON = _l["longitude"].values
_l.close()
CELL = 25e3
E0 = -9.0e6
NE = 720
BOX = 100e3
X0 = E0
NB = 180
px, py = TR.transform(lon, lat)
px = np.asarray(px)
py = np.asarray(py)
assert px.min() > X0 and px.max() < X0 + NB * BOX and py.min() > X0 and py.max() < X0 + NB * BOX
bi = np.floor((px - X0) / BOX).astype(int)
bj = np.floor((py - X0) / BOX).astype(int)  # ascending y
cell = bj * NB + bi  # flat id, row=y ascending, col=x
NC = NB * NB
xe = X0 + BOX * np.arange(NB + 1)
xc = 0.5 * (xe[1:] + xe[:-1])
INV = pyproj.Transformer.from_crs(PROJ_CRS, "EPSG:4326", always_xy=True)
CXX, CYY = np.meshgrid(xc, xc)
CLON, CLAT = INV.transform(CXX, CYY)
CLON = np.asarray(CLON)
CLAT = np.asarray(CLAT)
ece = E0 + CELL * np.arange(NE + 1)  # 25 km cell edges, ascending
ecc = 0.5 * (ece[1:] + ece[:-1])
EOCEAN = (LMASK == 255) & (np.nan_to_num(ELAT, nan=-90) >= 60)  # ocean cells north of 60N (rows y descending)
ELAND = np.isin(LMASK, [0, 101])


def bcount(mask=None, idx=None):
    c = cell if idx is None else cell[idx]
    if mask is not None:
        c = c[mask]
    return np.bincount(c, minlength=NC)


def ndistinct(values, mask=None):
    c = cell if mask is None else cell[mask]
    v = values if mask is None else values[mask]
    pairs = np.unique(np.column_stack((c, v)), axis=0)
    return np.bincount(pairs[:, 0].astype(np.int64), minlength=NC)


def ocean_grid_coverage(counts, ocean, threshold=1):
    """Fraction of ocean-mask cells meeting the count threshold (not area-weighted)."""
    return float(((counts >= threshold) & ocean).sum() / ocean.sum())


def level_mean_std(arr):
    """Population moments per level, retaining NaN for unobserved levels."""
    count = (~np.isnan(arr)).sum(axis=0)
    dtype = np.result_type(arr.dtype, np.float32)
    mean = np.full(arr.shape[1], np.nan, dtype=dtype)
    np.divide(np.nansum(arr, axis=0), count, out=mean, where=count > 0)
    variance = np.full(arr.shape[1], np.nan, dtype=dtype)
    np.divide(np.nansum((arr - mean) ** 2, axis=0), count, out=variance, where=count > 0)
    return mean, np.sqrt(variance)


# WOA23 0.25 deg mask: only used for the lat/lon-grid tables in section 3e
WOA_MASK_FILE = ROOT / "data/woa23_raw/woa23_decav91C0_t00_04.nc"
_w = xr.open_dataset(WOA_MASK_FILE, decode_times=False)
_wl = _w["lat"].values
_sel = _wl >= 60
WLAT = _wl[_sel]
WLON = _w["lon"].values
WOCEAN = ~np.isnan(_w["t_an"].isel(time=0, depth=0).values[_sel, :])
_w.close()

# ocean area per 100 km box = number of EASE2 ocean cells (>=60N) x 625 km2 (equal-area grid)
_cnt_cells = EOCEAN.reshape(NB, 4, NB, 4).sum(axis=(1, 3))[::-1]  # flip rows -> y ascending
OCEAN_AREA = _cnt_cells.ravel().astype(float) * (CELL * CELL)  # m2
OCEAN_BOX = OCEAN_AREA > 0
TOT_OCEAN_AREA = OCEAN_AREA.sum()
K["ocean_area_mkm2"] = TOT_OCEAN_AREA / 1e12
K["n_ocean_boxes"] = int(OCEAN_BOX.sum())


def area_frac(boolbox):
    """Fraction of the 60-90N EASE2 ocean area lying in boxes where boolbox is True."""
    return float(OCEAN_AREA[boolbox].sum() / TOT_OCEAN_AREA)


# ----------------------------------------------------------------------------- map helpers
import cartopy.crs as ccrs

MAPCRS = ccrs.LambertAzimuthalEqualArea(central_longitude=0, central_latitude=90)
PC = ccrs.PlateCarree()
_xt = MAPCRS.transform_points(PC, np.array([0.0, 100.0]), np.array([75.0, 62.0]))
_xr, _yr = TR.transform(np.array([0.0, 100.0]), np.array([75.0, 62.0]))
assert np.allclose(_xt[:, 0], _xr, atol=50) and np.allclose(_xt[:, 1], _yr, atol=50), "cartopy/pyproj mismatch"
_LAND_PLOT = np.ma.masked_where(~ELAND, np.ones(ELAND.shape))[::-1]


def polar_ax(fig, nr, nc, k, title=None, lat_min=60):
    ax = fig.add_subplot(nr, nc, k, projection=MAPCRS)
    R_ = float(np.hypot(*MAPCRS.transform_point(0.0, lat_min, PC)))
    ax.set_extent([-R_, R_, -R_, R_], MAPCRS)
    ax.set_autoscale_on(False)
    th = np.linspace(0, 2 * np.pi, 120)
    circ = mpath.Path(np.vstack([np.sin(th), np.cos(th)]).T * 0.5 + 0.5)
    ax.set_boundary(circ, transform=ax.transAxes)
    # Land and land ice are codes 0 and 101; code 255 is water.
    ax.pcolormesh(ece, ece, _LAND_PLOT, cmap=matplotlib.colors.ListedColormap(["0.8"]), transform=MAPCRS, zorder=1,
                  rasterized=True)
    ax.gridlines(draw_labels=False, linewidth=0.3, color="k", alpha=0.4, linestyle=":",
                 ylocs=[70, 80], xlocs=np.arange(-180, 181, 45))
    ax.grid(False)
    if title:
        ax.set_title(title, fontsize=9)
    return ax


def box_pcolor(ax, arr, norm=None, cmap="viridis", vmin=None, vmax=None):
    a = np.ma.masked_where(~np.isfinite(arr) | (arr <= 0) if norm is not None else ~np.isfinite(arr), arr)
    kw = dict(cmap=cmap, transform=MAPCRS, zorder=2, rasterized=True)
    if norm is not None:
        kw["norm"] = norm
    else:
        kw.update(vmin=vmin, vmax=vmax)
    return ax.pcolormesh(xe, xe, a, **kw)


def cbar(fig, m, ax, label, **kw):
    cb = fig.colorbar(m, ax=ax, orientation="horizontal", pad=0.03, fraction=0.045, shrink=0.8, **kw)
    cb.set_label(label)
    return cb


def grid2(a):
    return a.reshape(NB, NB)


# ============================================================================= SECTION 1: overview
log("section 1")
nlev_edges = [(0, 0), (1, 1), (2, 5), (6, 10), (11, 20), (21, 40), (41, 57), (58, 80), (81, 102)]
rows = []
for a, b in nlev_edges:
    m = (nlevT >= a) & (nlevT <= b)
    rows.append({"levels with valid T": f"{a}" if a == b else f"{a}-{b}", "profiles": int(m.sum()),
                 "fraction (%)": 100 * m.mean(), "cumulative (%)": 100 * (nlevT <= b).mean()})
t_nlev = pd.DataFrame(rows)
savetab(t_nlev, "s1_levels_per_profile_distribution")

TH = [0, 10, 50, 100, 300, 700, 1000, 1500, 2000, 3000, 4000, 5500]
rows = []
for th in TH:
    rows.append({"max valid depth >= (m)": th,
                 "profiles with T (N)": int((maxdepth >= th).sum()),
                 "fraction of profiles with T (%)": 100 * (maxdepth >= th).sum() / hasT.sum(),
                 "profiles with S (N)": int((maxdepthS >= th).sum()),
                 "fraction of profiles with S (%)": 100 * (maxdepthS >= th).sum() / max(hasS.sum(), 1)})
t_md = pd.DataFrame(rows)
savetab(t_md, "s1_max_depth_exceedance")
md_frac = {th: float((maxdepth >= th).sum() / hasT.sum()) for th in TH}

n_noS = int((hasT & ~hasS).sum())
n_noT = int((~hasT & hasS).sum())
n_both = int((hasT & hasS).sum())
n_single = int((nlevT == 1).sum())
n_leq3 = int(((nlevT >= 1) & (nlevT <= 3)).sum())
gaps = hasT & ((lastT - firstT + 1) > nlevT)
cellsT = int(vT.sum())
cellsS = int(vS.sum())
cells_TnoS = int((vT & ~vS).sum())
cells_SnoT = int((vS & ~vT).sum())
K.update(hasT=int(hasT.sum()), hasS=int(hasS.sum()), both=n_both, noS=n_noS, single=n_single,
         frac_noS=n_noS / max(hasT.sum(), 1), med_nlev=float(np.median(nlevT[hasT])),
         med_maxdepth=float(np.nanmedian(maxdepth)), md_frac=md_frac,
         nan_frac_T=1 - cellsT / vT.size)


def desc(name, arr, unit):
    v = arr[np.isfinite(arr)]
    q = np.percentile(v, [0.1, 1, 50, 99, 99.9])
    return {"variable": name, "unit": unit, "N values": int(v.size), "min": float(v.min()), "p0.1": q[0], "p1": q[1],
            "median": q[2], "mean": float(v.mean()), "std": float(v.std()), "p99": q[3], "p99.9": q[4],
            "max": float(v.max())}


t_rng = pd.DataFrame([desc("temp (in-situ T)", T, "degC"), desc("salt (practical S)", S, "-"),
                      desc("dens (in-situ density)", DN, "kg m-3"), desc("pres", P, "dbar")])
savetab(t_rng, "s1_value_ranges")
K["T_min"], K["T_max"] = float(np.nanmin(T)), float(np.nanmax(T))
K["S_min"], K["S_max"] = float(np.nanmin(S)), float(np.nanmax(S))

fig, axs = plt.subplots(1, 3, figsize=(15, 4))
axs[0].bar(np.arange(0, NZ + 1), np.bincount(nlevT, minlength=NZ + 1), width=0.9, color="C0")
axs[0].set_yscale("log")
axs[0].set_xlabel("number of levels with valid T")
axs[0].set_ylabel("profiles")
axs[0].set_title("(a) levels per profile")
mdc = np.bincount(lastT[hasT], minlength=NZ)
axs[1].bar(np.arange(NZ), mdc, width=0.9, color="C1")
axs[1].set_yscale("log")
axs[1].set_xticks(np.arange(0, NZ, 10))
axs[1].set_xticklabels([f"{int(depth[i])}" for i in range(0, NZ, 10)], rotation=45)
axs[1].set_xlabel("deepest valid T level (m; index axis, WOA levels)")
axs[1].set_ylabel("profiles")
axs[1].set_title("(b) maximum valid depth per profile")
dd = np.sort(maxdepth[hasT])
zz = np.linspace(0, 5500, 600)
axs[2].plot(zz, 100 * (1 - np.searchsorted(dd, zz, side="left") / dd.size), lw=2)
axs[2].set_xscale("symlog", linthresh=50)
axs[2].set_xlim(0, 5600)
for th in (100, 300, 700, 1000, 1500, 2000, 3000):
    axs[2].axvline(th, color="grey", lw=0.5)
axs[2].set_xlabel("depth z (m)")
axs[2].set_ylabel("% of profiles with valid T at z or deeper")
axs[2].set_title("(c) fraction of profiles reaching depth z")
f1 = savefig(fig, "s1_levels_and_maxdepth",
             "Profile completeness and deepest non-missing temperature level; the exceedance curve describes depth reach, not continuous coverage")

add(1, "## 1. Overview of the profile data\n",
    f"The scatter files hold **{N:,}** profiles from {K['year_min']} to {K['year_max']}, each on the {NZ} WOA standard "
    f"depths (0-{int(depth[-1])} m). Counts below refer to a profile as having a *valid T* if at least one level is "
    f"non-NaN. Here *valid* means present after missing-value normalisation, not verified by an independent quality assessment. "
    f"A profile is one published record; a level-value is one temperature or salinity entry at a standard depth.\n",
    mdtable(pd.DataFrame([
        ["Profiles in total", f"{N:,}", "100%"],
        ["Profiles with valid T", f"{hasT.sum():,}", pct(hasT.sum(), N)],
        ["Profiles with valid S", f"{hasS.sum():,}", pct(hasS.sum(), N)],
        ["Profiles with both T and S", f"{n_both:,}", pct(n_both, N)],
        ["T valid but S entirely missing", f"{n_noS:,}", pct(n_noS, hasT.sum()) + " of profiles with T"],
        ["S valid but T entirely missing", f"{n_noT:,}", pct(n_noT, N)],
        ["Profiles with a single valid T level", f"{n_single:,}", pct(n_single, hasT.sum()) + " of profiles with T"],
        ["Profiles with 1-3 valid T levels", f"{n_leq3:,}", pct(n_leq3, hasT.sum()) + " of profiles with T"],
        ["Profiles with internal gaps (missing levels between first and last valid)", f"{int(gaps.sum()):,}",
         pct(gaps.sum(), hasT.sum())],
        ["T cells valid (of N x 102)", f"{cellsT:,}", pct(cellsT, vT.size)],
        ["S cells valid", f"{cellsS:,}", pct(cellsS, vT.size)],
        ["T valid where S missing (cells)", f"{cells_TnoS:,}", pct(cells_TnoS, cellsT) + " of T cells"],
        ["S valid where T missing (cells)", f"{cells_SnoT:,}", pct(cells_SnoT, cellsS) + " of S cells"],
        ["Median number of valid T levels per profile", f"{K['med_nlev']:.0f}", ""],
        ["Median maximum valid depth (m)", f"{K['med_maxdepth']:.0f}", ""],
    ], columns=["quantity", "value", "share"])),
    "**Levels with valid T per profile**\n", mdtable(t_nlev, dec=2),
    "**Maximum valid depth** (deepest non-NaN T level; S analogous)\n", mdtable(t_md, dec=2),
    "**Value ranges** (non-missing level-values; percentiles and population standard deviation)\n", mdtable(t_rng, dec=3),
    f1,
    f"**Interpretation.** About {pct(K['nan_frac_T'], 1, 0)} of all temperature entries on the standard-depth grid are missing. "
    f"This reflects limited depth reach and missing levels, not an estimated rejection rate. Only "
    f"{pct(md_frac[1000], 1)} of profiles with T reach 1000 m and {pct(md_frac[2000], 1)} reach 2000 m. "
    f"Profiles reaching 5500 m: {int((maxdepth >= 5500).sum()):,}.\n")
log("section 1 done")

# ============================================================================= SECTION 2: time
log("section 2")
ycount = pd.Series(year).value_counts().sort_index()
years_all = np.arange(K["year_min"], K["year_max"] + 1)
ycount = ycount.reindex(years_all, fill_value=0)
mcount = np.bincount(month, minlength=13)[1:]
ym = np.zeros((len(years_all), 12), int)
np.add.at(ym, (year - years_all[0], month - 1), 1)
lat_band = np.digitize(lat, [70, 80])  # 0: 60-70, 1: 70-80, 2: 80-90
yb = np.zeros((len(years_all), 3), int)
np.add.at(yb, (year - years_all[0], lat_band), 1)
t_year = pd.DataFrame({"year": years_all, "profiles": ycount.values,
                       "lat 60-70N": yb[:, 0], "lat 70-80N": yb[:, 1], "lat 80-90N": yb[:, 2],
                       "profiles with T>=700 m": np.bincount(year[maxdepth >= 700] - years_all[0],
                                                             minlength=len(years_all)),
                       "% with S": [100 * hasS[year == y].mean() if (year == y).any() else np.nan for y in years_all]})
savetab(t_year, "s2_profiles_per_year")
t_month = pd.DataFrame({"month": MONTH_NAMES, "profiles": mcount, "fraction (%)": 100 * mcount / N,
                        "profiles lat>80N": np.bincount(month[lat > 80], minlength=13)[1:],
                        "profiles lat>70N": np.bincount(month[lat > 70], minlength=13)[1:]})
savetab(t_month, "s2_profiles_per_month_of_year")
savetab(pd.DataFrame(ym, index=years_all, columns=MONTH_NAMES), "s2_profiles_year_x_month", index=True)
decc = np.bincount(dec_idx[dec_idx >= 0], minlength=4)
t_dec = pd.DataFrame({"decade": DEC_LAB, "profiles": decc, "fraction (%)": 100 * decc / N,
                      "profiles per year": decc / np.array([b - a + 1 for a, b in DECADES]),
                      "median max depth (m)": [np.nanmedian(maxdepth[dec_idx == i]) for i in range(4)],
                      "% reaching 1000 m": [100 * np.mean(maxdepth[dec_idx == i] >= 1000) for i in range(4)],
                      "% with S": [100 * hasS[dec_idx == i].mean() for i in range(4)]})
savetab(t_dec, "s2_profiles_per_decade")
rows = []
for nm, grp in list(SEASONS_DJF.items()) + list(SEASONS_WOA.items()):
    m = np.isin(month, grp)
    rows.append({"season": nm, "scheme": "calendar (DJF..)" if nm in SEASONS_DJF else "WOA-style (JFM..)",
                 "profiles": int(m.sum()), "fraction (%)": 100 * m.mean(),
                 "profiles lat>80N": int((m & (lat > 80)).sum()),
                 "% of lat>80N total": 100 * (m & (lat > 80)).sum() / (lat > 80).sum(),
                 "profiles lat>70N": int((m & (lat > 70)).sum())})
t_seas = pd.DataFrame(rows)
savetab(t_seas, "s2_profiles_per_season")

pk = int(ycount.idxmax())
share_post05 = float(ycount.loc[2005:].sum() / N)
share_post13 = float(ycount.loc[2013:].sum() / N)
jas = np.isin(month, (7, 8, 9)).mean()
jfm = np.isin(month, (1, 2, 3)).mean()
jas80 = np.isin(month[lat > 80], (7, 8, 9)).mean()
jfm80 = np.isin(month[lat > 80], (1, 2, 3)).mean()
pre05 = ycount.loc[1983:2004].mean()
post05 = ycount.loc[2005:].mean()
K.update(peak_year=pk, peak_year_n=int(ycount.loc[pk]), share_post05=share_post05, jas=jas, jfm=jfm, jas80=jas80,
         jfm80=jfm80, pre05_mean=float(pre05), post05_mean=float(post05), ymin_year=int(ycount.idxmin()),
         ymin_n=int(ycount.min()), mmax=MONTH_NAMES[int(np.argmax(mcount))], mmin=MONTH_NAMES[int(np.argmin(mcount))],
         mratio=float(mcount.max() / mcount.min()), n80=int((lat > 80).sum()))

fig, axs = plt.subplots(2, 2, figsize=(14, 8))
ax = axs[0, 0]
bot = np.zeros(len(years_all))
for k_, (lab, c) in enumerate(zip(["60-70N", "70-80N", "80-90N"], ["C0", "C1", "C3"])):
    ax.bar(years_all, yb[:, k_], bottom=bot, color=c, label=lab)
    bot += yb[:, k_]
ax.set_xlabel("year")
ax.set_ylabel("profiles")
ax.set_title("(a) profiles per year by latitude band")
ax.legend()
ax = axs[0, 1]
ax.bar(range(1, 13), mcount, color="C2")
ax.set_xticks(range(1, 13))
ax.set_xticklabels(MONTH_NAMES)
ax.set_ylabel("profiles")
ax.set_title("(b) profiles per month of year (all years)")
ax2 = ax.twinx()
ax2.plot(range(1, 13), np.bincount(month[lat > 80], minlength=13)[1:], "k.-", label="lat>80N (right)")
ax2.set_ylabel("profiles N of 80N")
ax2.grid(False)
ax2.legend(loc="upper left")
ax = axs[1, 0]
im = ax.pcolormesh(np.arange(1, 14) - 0.5, np.append(years_all, years_all[-1] + 1) - 0.5, np.ma.masked_equal(ym, 0),
                   norm=LogNorm(), cmap="magma_r")
ax.set_xticks(range(1, 13))
ax.set_xticklabels(MONTH_NAMES)
ax.set_ylabel("year")
ax.set_title("(c) profiles per year x month")
ax.grid(False)
fig.colorbar(im, ax=ax, label="profiles")
ax = axs[1, 1]
ymm = np.zeros((len(years_all), 12), int)
m80 = lat > 80
np.add.at(ymm, (year[m80] - years_all[0], month[m80] - 1), 1)
im = ax.pcolormesh(np.arange(1, 14) - 0.5, np.append(years_all, years_all[-1] + 1) - 0.5, np.ma.masked_equal(ymm, 0),
                   norm=LogNorm(), cmap="magma_r")
ax.set_xticks(range(1, 13))
ax.set_xticklabels(MONTH_NAMES)
ax.set_ylabel("year")
ax.set_title("(d) profiles north of 80N per year x month")
ax.grid(False)
fig.colorbar(im, ax=ax, label="profiles")
f2 = savefig(fig, "s2_time_distribution",
             "Profile counts by year, latitude band and calendar month; year-by-month panels use logarithmic colour scales")

add(2, "## 2. Temporal distribution\n",
    f2,
    "**Profiles per decade**\n", mdtable(t_dec, dec={"fraction (%)": 1, "profiles per year": 0, "median max depth (m)": 0,
                                                    "% reaching 1000 m": 1, "% with S": 1}),
    "**Profiles per season** (month-of-year bins, all years)\n", mdtable(t_seas, dec=1),
    "**Profiles per month of year**\n", mdtable(t_month, dec=1),
    f"The full year table is `tables/s2_profiles_per_year.csv` and the year x month matrix is "
    f"`tables/s2_profiles_year_x_month.csv`.\n",
    f"**Interpretation.** The record peaks in {pk} ({int(ycount.loc[pk]):,} profiles) and has its minimum in "
    f"{K['ymin_year']} ({K['ymin_n']:,}). Mean annual sampling is {pre05:,.0f} profiles/yr for 1983-2004 and "
    f"{post05:,.0f} profiles/yr for 2005-2023 (x{post05 / pre05:.1f}); {pct(share_post05, 1)} of all profiles are from "
    f"2005 onward. Sampling is seasonally uneven: {pct(jas, 1)} of profiles fall in Jul-Sep against {pct(jfm, 1)} in "
    f"Jan-Mar; north of 80N the corresponding shares are {pct(jas80, 1)} and {pct(jfm80, 1)}. "
    f"Across all profiles, the most-sampled month is {K['mmax']} and the least-sampled is {K['mmin']} "
    f"(count ratio {K['mratio']:.1f}). These are sampling distributions, not evidence of a particular platform, "
    f"sea-ice mechanism or environmental trend. The four reporting periods are not equal in length: "
    f"2013-2023 contains eleven years, the others ten. Quarters pool calendar months across years; DJF is not grouped into winter-years.\n")
log("section 2 done")

# ============================================================================= SECTION 3: space
log("section 3")
THR = [1, 5, 10, 50, 100, 1000]


def box_stats(c, label):
    occ = c > 0
    d = {"slice": label, "boxes >=1": int(occ.sum())}
    for t in THR[1:]:
        d[f"boxes >={t}"] = int((c >= t).sum())
    d["median per occupied"] = float(np.median(c[occ])) if occ.any() else np.nan
    d["mean per occupied"] = float(c[occ].mean()) if occ.any() else np.nan
    d["max"] = int(c.max())
    d["% ocean area with 0 profiles"] = 100 * (1 - area_frac(occ & OCEAN_BOX))
    d["% ocean area with >=10"] = 100 * area_frac((c >= 10) & OCEAN_BOX)
    return d


cnt = bcount()
nyr = ndistinct(year)
nmo = ndistinct(month)
nym = ndistinct(year * 12 + month)
occ_all = cnt > 0
occ_nonocean = int((occ_all & ~OCEAN_BOX).sum())
t_boxdist = pd.DataFrame([box_stats(cnt, "all years, all months")])
savetab(t_boxdist, "s3_box_count_distribution")
rows = [box_stats(cnt, "all")]
dec_cnt = []
for i, lab in enumerate(DEC_LAB):
    c = bcount(dec_idx == i)
    dec_cnt.append(c)
    rows.append(box_stats(c, "decade " + lab))
seas_cnt = {}
for nm, grp in SEASONS_WOA.items():
    c = bcount(np.isin(month, grp))
    seas_cnt[nm] = c
    rows.append(box_stats(c, "season " + nm))
for nm, grp in SEASONS_DJF.items():
    rows.append(box_stats(bcount(np.isin(month, grp)), "season " + nm))
t_boxslices = pd.DataFrame(rows)
savetab(t_boxslices, "s3_box_count_distribution_by_slice")

# CSV of 100 km boxes
oi = np.nonzero(occ_all)[0]
t_boxes = pd.DataFrame({"x": xc[oi % NB], "y": xc[oi // NB], "lon": CLON.ravel()[oi], "lat": CLAT.ravel()[oi],
                        "count": cnt[oi], "n_years": nyr[oi], "n_months": nmo[oi], "n_year_months": nym[oi],
                        "ocean_area_km2": OCEAN_AREA[oi] / 1e6})
savetab(t_boxes, "s3_box100km_counts")

zero_area = 1 - area_frac(occ_all & OCEAN_BOX)
K.update(nbox1=int((cnt >= 1).sum()), nbox10=int((cnt >= 10).sum()), nbox100=int((cnt >= 100).sum()),
         zero_area=zero_area, cov10_all=area_frac((cnt >= 10) & OCEAN_BOX), cnt_med=float(np.median(cnt[occ_all])),
         cnt_mean=float(cnt[occ_all].mean()), cnt_max=int(cnt.max()), occ_nonocean=occ_nonocean)
top = np.argsort(cnt)[::-1][:5]
K["topbox"] = [(float(CLAT.ravel()[i]), float(CLON.ravel()[i]), int(cnt[i])) for i in top]

# maps
fig = plt.figure(figsize=(7.5, 8))
ax = polar_ax(fig, 1, 1, 1, "Profiles per 100 km x 100 km box (1983-2023)")
m = box_pcolor(ax, grid2(cnt), norm=LogNorm(1, cnt.max()), cmap="viridis")
cbar(fig, m, ax, "profiles per box (log)")
f3a = savefig(fig, "s3a_box_counts_all", "Profiles per equal-area 100 km box over the full record; logarithmic colour scale")

fig = plt.figure(figsize=(13, 13))
axl = []
for i in range(4):
    ax = polar_ax(fig, 2, 2, i + 1, f"{DEC_LAB[i]}  ({decc[i]:,} profiles)")
    m = box_pcolor(ax, grid2(dec_cnt[i]), norm=LogNorm(1, cnt.max()), cmap="viridis")
    axl.append(ax)
cbar(fig, m, axl, "profiles per 100 km box (log)")
f3b = savefig(fig, "s3b_box_counts_by_decade", "Profile counts per box by decade")

fig = plt.figure(figsize=(13, 13))
axl = []
for i, (nm, grp) in enumerate(SEASONS_WOA.items()):
    ax = polar_ax(fig, 2, 2, i + 1, f"{nm}  ({int(seas_cnt[nm].sum()):,} profiles)")
    m = box_pcolor(ax, grid2(seas_cnt[nm]), norm=LogNorm(1, cnt.max()), cmap="viridis")
    axl.append(ax)
cbar(fig, m, axl, "profiles per 100 km box (log)")
f3c = savefig(fig, "s3c_box_counts_by_season", "Profile counts per box by season")

fig = plt.figure(figsize=(14, 7.5))
ax1 = polar_ax(fig, 1, 2, 1, "Number of distinct years sampled per box")
m1 = box_pcolor(ax1, np.where(grid2(cnt) > 0, grid2(nyr), np.nan), norm=None, cmap="YlGnBu", vmin=1, vmax=41)
cbar(fig, m1, ax1, "distinct years (1-41)")
ax2 = polar_ax(fig, 1, 2, 2, "Number of distinct calendar months sampled per box")
m2 = box_pcolor(ax2, np.where(grid2(cnt) > 0, grid2(nmo), np.nan), norm=None, cmap="YlOrRd", vmin=1, vmax=12)
cbar(fig, m2, ax2, "distinct months of year (1-12)")
f3d = savefig(fig, "s3d_years_months_per_box",
             "Distinct years and calendar months represented in occupied boxes; blank boxes have no profiles, and twelve months need not occur in the same year")

# lat/lon grids
GRID = {}
for res in (1.0, 0.25):
    nlat_ = int(30 / res)
    nlon_ = int(360 / res)
    ilat = np.clip(np.floor((lat - 60) / res).astype(int), 0, nlat_ - 1)
    ilon = np.floor((lon180 + 180) / res).astype(int) % nlon_
    g = np.bincount(ilat * nlon_ + ilon, minlength=nlat_ * nlon_).reshape(nlat_, nlon_)
    GRID[res] = g
    nz = np.nonzero(g.ravel())[0]
    cl = 60 + res * (nz // nlon_ + 0.5)
    co = -180 + res * (nz % nlon_ + 0.5)
    savetab(pd.DataFrame({"lat": cl, "lon": co, "count": g.ravel()[nz]}), f"s3_grid_{res:g}deg_counts")
rows = []
for res, g in GRID.items():
    occ = g > 0
    # ocean cells on this grid from the WOA mask (0.25 deg blocks)
    d = {"grid": f"{res:g} deg lat/lon", "cells with data": int(occ.sum())}
    for t in THR[1:]:
        d[f"cells >={t}"] = int((g >= t).sum())
    d["median per occupied"] = float(np.median(g[occ]))
    d["mean per occupied"] = float(g[occ].mean())
    d["max"] = int(g.max())
    rows.append(d)
wo = WOCEAN
n025_ocean = int(wo.sum())
wo1 = wo.reshape(30, 4, 360, 4).any(axis=(1, 3))
n1_ocean = int(wo1.sum())
t_grids = pd.DataFrame(rows)
t_grids["ocean cells in grid (WOA23 0.25 deg mask)"] = [n1_ocean, n025_ocean]
t_grids["occupied cells outside ocean mask"] = [
    int(((GRID[1.0] > 0) & ~wo1).sum()), int(((GRID[0.25] > 0) & ~wo).sum())]
t_grids["% ocean cells with data"] = [
    100 * ocean_grid_coverage(GRID[1.0], wo1), 100 * ocean_grid_coverage(GRID[0.25], wo)]
savetab(t_grids, "s3_latlon_grid_distribution")
K["g1_frac"] = ocean_grid_coverage(GRID[1.0], wo1)
K["g025_frac"] = ocean_grid_coverage(GRID[0.25], wo)
K["g025_ge10"] = ocean_grid_coverage(GRID[0.25], wo, threshold=10)
K["g1_ge10"] = ocean_grid_coverage(GRID[1.0], wo1, threshold=10)

fig = plt.figure(figsize=(14, 7.5))
for k_, res in enumerate((1.0, 0.25)):
    ax = polar_ax(fig, 1, 2, k_ + 1, f"Profiles per {res:g} deg lat/lon cell")
    g = GRID[res]
    nlat_, nlon_ = g.shape
    m = ax.pcolormesh(np.linspace(-180, 180, nlon_ + 1), 60 + res * np.arange(nlat_ + 1), np.ma.masked_equal(g, 0),
                      norm=LogNorm(1, g.max()), cmap="viridis", transform=PC, zorder=2, rasterized=True)
    cbar(fig, m, ax, "profiles per cell (log)")
f3e = savefig(fig, "s3e_latlon_grid_counts", "Counts on 1 and 0.25 degree grids")

tb = t_boxdist.iloc[0]
add(3, "## 3. Spatial distribution\n",
    f"Profiles were binned on a 100 km x 100 km grid that is an exact 4x4 block aggregation of the NSIDC EASE-Grid 2.0 "
    f"North 25 km grid (EPSG:6931, Lambert azimuthal equal area, so every box has the same area; {NB}x{NB} boxes cover the "
    f"grid domain, all profiles are north of 60N). Ocean area is approximated from the EASE2 25 km land mask "
    f"(`EASE2_N25km.LOCImask_land50_coast0km`, code 255 = water; 625 km2 per cell whose centre is north of 60N): {K['ocean_area_mkm2']:.2f} million km2 in "
    f"{K['n_ocean_boxes']:,} boxes containing any ocean. This includes marginal seas such as Hudson Bay, the Baltic, "
    f"Gulf of Alaska and the parts of the Sea of Okhotsk north of 60N, so the 'ocean without data' fraction counts every one of "
    f"them. {occ_nonocean} occupied boxes contain no ocean cell of that mask. This coastal/mask mismatch is "
    f"retained in profile counts but contributes no ocean area. A box is counted as covered when its profile "
    f"count meets the threshold; this does not establish coverage throughout its ocean area, at every depth or in every year.\n",
    "### 3a. Counts per 100 km box\n", f3a,
    "**Distribution of box counts, and by decade and season** (`tables/s3_box_count_distribution*.csv`)\n",
    mdtable(t_boxslices[["slice", "boxes >=1", "boxes >=5", "boxes >=10", "boxes >=50", "boxes >=100", "boxes >=1000",
                         "median per occupied", "mean per occupied", "max", "% ocean area with 0 profiles",
                         "% ocean area with >=10"]], dec=1),
    f"**Interpretation.** {int(tb['boxes >=1']):,} boxes hold at least one profile; the median occupied box has "
    f"{tb['median per occupied']:.0f} profiles, the mean {tb['mean per occupied']:.0f}, the maximum "
    f"{int(tb['max']):,}. The difference between median and mean indicates strongly concentrated sampling; "
    f"counts alone cannot identify its platform sources or number of independent observations. "
    f"{pct(zero_area, 1)} of the 60-90N ocean area lies in boxes with no profile at all, and only "
    f"{pct(K['cov10_all'], 1)} lies in boxes with >=10 profiles over the whole 41 years. The densest boxes (centre "
    f"lat, lon, count): " + "; ".join(f"({a:.1f}N, {b:.1f}E, {c:,})" for a, b, c in K['topbox']) + ".\n",
    "### 3b. By decade\n", f3b,
    "### 3c. By season (WOA-style quarters)\n", f3c,
    "### 3d. Temporal diversity of each box\n", f3d,
    f"Boxes sampled in >=10 distinct years: {int((nyr >= 10).sum()):,}; in >=20 years: {int((nyr >= 20).sum()):,}; "
    f"in all 12 calendar months: {int((nmo == 12).sum()):,}; in >=6 calendar months: {int((nmo >= 6).sum()):,}. "
    f"Ocean-area fraction with data in >=10 distinct years: {pct(area_frac((nyr >= 10) & OCEAN_BOX), 1)}; "
    f"in >=6 calendar months: {pct(area_frac((nmo >= 6) & OCEAN_BOX), 1)}.\n",
    "### 3e. Lat/lon grids\n", f3e,
    mdtable(t_grids, dec=1),
    f"On the 1 deg grid {pct(K['g1_frac'], 1)} of WOA23 ocean cells north of 60N contain a profile "
    f"({pct(K['g1_ge10'], 1)} contain >=10); on the 0.25 deg grid {pct(K['g025_frac'], 1)} "
    f"({pct(K['g025_ge10'], 1)} with >=10). Numerators and denominators both use the WOA23 surface ocean mask "
    f"(a 1 deg cell is ocean if any of its sixteen 0.25 deg cells is ocean). Occupied cells outside that mask "
    f"are listed separately. These are fractions of cells, not ocean area; latitude/longitude cells shrink towards the pole "
    f"and should not be compared directly with the EASE2 area fractions. "
    f"Cell lists: `tables/s3_grid_1deg_counts.csv`, `tables/s3_grid_0.25deg_counts.csv`; "
    f"100 km box list: `tables/s3_box100km_counts.csv` (columns x, y [EPSG:6931 metres], lon, lat, count, n_years, "
    f"n_months, n_year_months, ocean_area_km2).\n")
K["nyr10_frac"] = area_frac((nyr >= 10) & OCEAN_BOX)
log("section 3 done")

# ============================================================================= SECTION 4: depth
log("section 4")
nT_lev = vT.sum(0)
nS_lev = vS.sum(0)
t_lev = pd.DataFrame({"level": np.arange(NZ), "depth_m": depth, "T_obs": nT_lev, "S_obs": nS_lev,
                      "T_fraction_of_profiles_%": 100 * nT_lev / N, "S_fraction_of_profiles_%": 100 * nS_lev / N})
for i, lab in enumerate(DEC_LAB):
    t_lev[f"T_obs_{lab}"] = vT[dec_idx == i].sum(0)
# box coverage per level (T)
cov1 = np.zeros(NZ, int)
cov10 = np.zeros(NZ, int)
cov1a = np.zeros(NZ)
cov10a = np.zeros(NZ)
for lv in range(NZ):
    c = np.bincount(cell[vT[:, lv]], minlength=NC)
    cov1[lv] = (c >= 1).sum()
    cov10[lv] = (c >= 10).sum()
    cov1a[lv] = area_frac((c >= 1) & OCEAN_BOX)
    cov10a[lv] = area_frac((c >= 10) & OCEAN_BOX)
t_lev["boxes_with_>=1_T"] = cov1
t_lev["boxes_with_>=10_T"] = cov10
t_lev["ocean_area_%_with_>=1"] = 100 * cov1a
t_lev["ocean_area_%_with_>=10"] = 100 * cov10a
savetab(t_lev, "s4_observations_per_level")
SEL = [0, 100, 500, 1000, 2000, 3000, 5000]
SELI = [int(np.where(depth == z)[0][0]) for z in SEL]
t_selcov = pd.DataFrame({"depth (m)": SEL, "T obs": nT_lev[SELI], "% of profiles": 100 * nT_lev[SELI] / N,
                         "boxes >=1": cov1[SELI], "boxes >=10": cov10[SELI],
                         "ocean area % >=1": 100 * cov1a[SELI], "ocean area % >=10": 100 * cov10a[SELI]})
savetab(t_selcov, "s4_box_coverage_selected_depths")
K["cov_sel"] = {z: (float(cov1a[i]), float(cov10a[i])) for z, i in zip(SEL, SELI)}
rng = [("0-1500 m (WOA monthly & seasonal, 57 levels)", slice(0, NL57)),
       (">1500-5500 m (45 levels; WOA annual only)", slice(NL57, NZ)),
       ("0-5500 m (WOA annual, 102 levels)", slice(0, NZ))]
rows = []
for nm, sl in rng:
    rows.append({"range": nm, "levels": len(range(*sl.indices(NZ))), "T obs": int(vT[:, sl].sum()),
                 "S obs": int(vS[:, sl].sum()),
                 "profiles with any T": int(vT[:, sl].any(1).sum()),
                 "% of all profiles": 100 * vT[:, sl].any(1).mean()})
t_rng_obs = pd.DataFrame(rows)
savetab(t_rng_obs, "s4_observations_by_depth_range")
K["obs_0_1500"] = int(vT[:, :NL57].sum())
K["obs_deep"] = int(vT[:, NL57:].sum())
K["prof_deep"] = int(vT[:, NL57:].any(1).sum())

fig, axs = plt.subplots(1, 3, figsize=(15, 6))
for ax, zmax in zip(axs[:2], (1500, 5500)):
    mk = depth <= zmax
    ax.plot(nT_lev[mk], depth[mk], "C3-", lw=2, label="T")
    ax.plot(nS_lev[mk], depth[mk], "C0--", lw=2, label="S")
    ax.set_xscale("log")
    ax.set_ylim(zmax, 0)
    ax.set_xlabel("number of valid observations")
    ax.set_ylabel("depth (m)")
    ax.set_title(f"observations per standard level, 0-{zmax} m")
    ax.legend()
ax = axs[2]
for i, lab in enumerate(DEC_LAB):
    ax.plot(100 * t_lev[f"T_obs_{lab}"] / decc[i], depth, label=lab)
ax.set_ylim(5500, 0)
ax.set_xlabel("% of the decade's profiles with valid T at level")
ax.set_ylabel("depth (m)")
ax.set_title("coverage with depth, by decade")
ax.legend()
f4a = savefig(fig, "s4a_observations_per_level", "Observations per level")

fig, axs = plt.subplots(1, 2, figsize=(12, 5.5))
ax = axs[0]
ax.plot(100 * cov1a, depth, label="boxes with >=1 profile")
ax.plot(100 * cov10a, depth, label="boxes with >=10 profiles")
ax.set_ylim(5500, 0)
ax.set_xlabel("% of 60-90N ocean area (100 km boxes)")
ax.set_ylabel("depth (m)")
ax.set_title("horizontal coverage vs depth (T)")
ax.legend()
ax = axs[1]
ax.plot(cov1, depth, label=">=1 profile")
ax.plot(cov10, depth, label=">=10 profiles")
ax.set_xscale("log")
ax.set_ylim(5500, 0)
ax.set_xlabel("number of 100 km boxes")
ax.set_ylabel("depth (m)")
ax.set_title("occupied boxes vs depth")
ax.legend()
f4b = savefig(fig, "s4b_box_coverage_vs_depth", "Box coverage vs depth")

fig = plt.figure(figsize=(18, 9))
axl = []
for k_, (z, li) in enumerate(zip(SEL, SELI)):
    c = np.bincount(cell[vT[:, li]], minlength=NC)
    ax = polar_ax(fig, 2, 4, k_ + 1, f"T observations at {z} m")
    m = box_pcolor(ax, grid2(c), norm=LogNorm(1, cnt.max()), cmap="viridis")
    axl.append(ax)
cbar(fig, m, axl, "profiles per 100 km box (log)")
f4c = savefig(fig, "s4c_box_counts_by_depth", "Counts per box at selected depths")

dd_tab = t_lev[["depth_m"] + [f"T_obs_{lab}" for lab in DEC_LAB]].iloc[[0, 4, 10, 20, 24, 28, 36, 40, 46, 56, 66, 76, 86, 101]]
add(4, "## 4. Depth distribution\n",
    "Counts of non-missing level-values at each of the 102 standard depths (full table `tables/s4_observations_per_level.csv`, "
    "including counts by reporting period and box coverage at every level). Reaching a depth and having a value at that "
    "depth are different criteria; Section 1 reports depth reach, whereas this section reports values at individual levels.\n", f4a,
    "**Observations by depth range**\n", mdtable(t_rng_obs, dec=1),
    f"There are {K['obs_0_1500']:,} T observations in 0-1500 m but only {K['obs_deep']:,} below 1500 m "
    f"({pct(K['obs_deep'], K['obs_0_1500'] + K['obs_deep'], 2)} of all T observations) from {K['prof_deep']:,} profiles.\n",
    "**Horizontal coverage at selected depths** (T; boxes of 100 km; ocean area from the EASE2 mask)\n",
    mdtable(t_selcov, dec=1), f4b,
    "**Profiles with valid T per level and decade** (selected levels; all levels in the CSV)\n", mdtable(dd_tab, dec=0),
    "**Maps of counts at selected depths**\n", f4c)
fr0, fr5, fr10 = (float(vT[:, i].mean()) for i in (0, 1, 2))
n_mismatch = int((vT != vS).sum())
K.update(fill0=fr0, fill5=fr5, fill10=fr10, TS_mismatch=n_mismatch, n_valid_TS=int((vT & vS).sum()))
add(4, "### Surface-level availability and T/S consistency\n",
    f"The 0 m level contains **{int(vT[:, 0].sum()):,}** valid T values and {int(vS[:, 0].sum()):,} valid S values "
    f"({pct(1 - fr0, 1, 1)} empty). The 5 m level is {pct(fr5, 1)} filled and the 10 m level {pct(fr10, 1)} "
    f"(15 m: {pct(float(vT[:, 3].mean()), 1)}). These are observed availability patterns in the released arrays. "
    f"They do not by themselves identify the upstream interpolation behaviour that produced them. "
    f"Values at 5 or 10 m must not be presented as observations at 0 m or as full mixed-layer coverage. "
    f"T and S are present together in {int((vT & vS).sum()):,} level-values; entries where exactly one of "
    f"them is valid: {n_mismatch:,}. Missing values are NaN after `scatter_io.load_scatter` / `data/scatter_all.npz` "
    f"normalisation. If the cache is absent, the yearly files are loaded through the same normalising loader.\n")
log("section 4 done")

# ============================================================================= SECTION 5: regions
log("section 5")
REG_NAMES = ["Canadian Arctic Archipelago", "Hudson Bay / Foxe Basin", "Central Arctic (>80N: Nansen/Amundsen/Makarov)",
             "Baffin Bay / Davis Strait", "Fram Strait / Greenland Sea", "Iceland Sea", "Norwegian Sea", "Barents Sea",
             "Kara Sea", "Laptev Sea", "East Siberian Sea", "Chukchi Sea", "Beaufort Sea / Canada Basin",
             "Bering Sea / Bering Strait", "Labrador Sea / subpolar N Atlantic", "Other"]
REG_DEF = ["polygon (-128,70) (-120,68.5) (-100,67) (-82,67.5) (-80,73.5) (-72,76) (-66,78) (-60,81) (-62,83) (-85,83) "
           "(-105,80) (-120,77) (-126,74) [lon E, lat N]",
           "lon -95..-76, lat 51..66", "lat >= 80 (anywhere not already in the CAA polygon)",
           "lon -80..-50, lat 64..80", "lon -35..15, lat 72..80", "lon -40..-8, lat 64.5..72",
           "lon -12..18, lat 62..72", "lon 18..55, lat 66..80 plus lon 15..18, lat 72..80", "lon 55..100, lat 66..80",
           "lon 100..140, lat 66..80", "lon 140..175, lat 66..80",
           "lon 175..180 and -180..-157, lat 66..73", "lon -157..-100, lat 68..80; lon 175..180 / -180..-157, lat 73..80",
           "lat < 66, lon >= 160 or <= -150", "lat < 64.5, lon -75..-12; lat < 62, lon -12..12",
           "everything else (Baltic/Scandinavian coast, Okhotsk side, unassigned strips)"]
CAA_POLY = mpath.Path(np.array([(-128, 70), (-120, 68.5), (-100, 67), (-82, 67.5), (-80, 73.5), (-72, 76), (-66, 78),
                                (-60, 81), (-62, 83), (-85, 83), (-105, 80), (-120, 77), (-126, 74), (-128, 70)], float))


def assign_regions(lo, la):
    lo = ((np.asarray(lo) + 180.0) % 360.0) - 180.0
    la = np.asarray(la)
    reg = np.full(lo.shape, len(REG_NAMES) - 1)
    done = np.zeros(lo.shape, bool)

    def put(i, m):
        nonlocal done
        mm = m & ~done
        reg[mm] = i
        done |= mm

    put(0, CAA_POLY.contains_points(np.c_[lo.ravel(), la.ravel()]).reshape(lo.shape))
    put(1, (lo >= -95) & (lo <= -76) & (la >= 51) & (la <= 66))
    put(2, la >= 80)
    put(3, (lo >= -80) & (lo < -50) & (la >= 64) & (la < 80))
    put(4, (lo >= -35) & (lo < 15) & (la >= 72) & (la < 80))
    put(5, (lo >= -40) & (lo < -8) & (la >= 64.5) & (la < 72))
    put(6, (lo >= -12) & (lo < 18) & (la >= 62) & (la < 72))
    put(7, ((lo >= 18) & (lo < 55) & (la >= 66) & (la < 80)) | ((lo >= 15) & (lo < 18) & (la >= 72) & (la < 80)))
    put(8, (lo >= 55) & (lo < 100) & (la >= 66) & (la < 80))
    put(9, (lo >= 100) & (lo < 140) & (la >= 66) & (la < 80))
    put(10, (lo >= 140) & (lo < 175) & (la >= 66) & (la < 80))
    put(11, ((lo >= 175) | (lo <= -157)) & (la >= 66) & (la < 73))
    put(12, ((lo >= -157) & (lo < -100) & (la >= 68) & (la < 80)) | (((lo >= 175) | (lo <= -157)) & (la >= 73) & (la < 80)))
    put(13, (la < 66) & ((lo >= 160) | (lo <= -150)))
    put(14, ((la < 64.5) & (lo >= -75) & (lo < -12)) | ((la < 62) & (lo >= -12) & (lo < 12)))
    return reg


REG = assign_regions(lon, lat)
NREG = len(REG_NAMES)
EREG = assign_regions(np.nan_to_num(ELON, nan=0.0), np.nan_to_num(ELAT, nan=0.0))  # region of every EASE2 cell
_rr, _cc = np.meshgrid(np.arange(NE), np.arange(NE), indexing="ij")
ECELLBOX = ((NE - 1 - _rr) // 4) * NB + _cc // 4  # 100 km box id of each 25 km cell (y ascending)
NORDIC = np.isin(REG, [4, 5, 6])
CENTRAL = REG == 2
K["nordic_n"] = int(NORDIC.sum())
rows = []
for i, nm in enumerate(REG_NAMES):
    m = REG == i
    wm = (EREG == i) & EOCEAN
    ar = float(wm.sum() * CELL * CELL)
    d = {"region": nm, "profiles": int(m.sum()), "% of profiles": 100 * m.mean()}
    for k_, lab in enumerate(DEC_LAB):
        d[lab] = int((m & (dec_idx == k_)).sum())
    d["median max depth (m)"] = float(np.nanmedian(maxdepth[m])) if m.any() else np.nan
    d["% reaching 1000 m"] = 100 * float(np.mean(maxdepth[m] >= 1000)) if m.any() else np.nan
    d["% with S"] = 100 * float(hasS[m].mean()) if m.any() else np.nan
    d["ocean area (1e6 km2)"] = ar / 1e12
    d["% area in boxes >=1"] = 100 * float((wm & (cnt[ECELLBOX] >= 1)).sum() * CELL * CELL / ar) if ar > 0 else np.nan
    d["% area in boxes >=10"] = 100 * float((wm & (cnt[ECELLBOX] >= 10)).sum() * CELL * CELL / ar) if ar > 0 else np.nan
    rows.append(d)
t_reg = pd.DataFrame(rows)
savetab(t_reg, "s5_regions")
savetab(pd.DataFrame({"region": REG_NAMES, "definition": REG_DEF}), "s5_region_definitions")
K["reg"] = t_reg.set_index("region")

fig = plt.figure(figsize=(10, 9))
ax = polar_ax(fig, 1, 1, 1, "Regions used in this report (EASE2 ocean cells)")
cm = plt.get_cmap("tab20", NREG)
wr = np.ma.masked_where(~EOCEAN, EREG)[::-1]
ax.pcolormesh(ece, ece, wr, cmap=cm, vmin=-0.5, vmax=NREG - 0.5, transform=MAPCRS, zorder=2, rasterized=True)
from matplotlib.patches import Patch

ax.legend(handles=[Patch(color=cm(i), label=REG_NAMES[i]) for i in range(NREG)], loc="center left",
          bbox_to_anchor=(1.02, 0.5), fontsize=8)
f5a = savefig(fig, "s5a_regions_map", "Region definitions")
rm = pd.DataFrame(np.array([[np.nanmedian(maxdepth[(REG == i) & (dec_idx == k_)]) if ((REG == i) & (dec_idx == k_)).any()
                              else np.nan for k_ in range(4)] for i in range(NREG)]), index=REG_NAMES, columns=DEC_LAB)
savetab(rm, "s5_median_max_depth_region_decade", index=True)
add(5, "## 5. Regional distribution\n",
    "Regions are simple lon/lat boxes/polygons (applied in priority order, first match wins) defined in this script; "
    "they are approximate reporting regions, not validated oceanographic boundaries. Every profile is assigned once, "
    "and the central-Arctic region uses lat >=80N after excluding the Canadian Arctic Archipelago polygon; "
    "the historical '>80N' label in tables is shorthand for that rule. All regions are clipped by the available "
    "north-of-60N dataset and ocean mask. Regional area coverage uses the count in each whole 100 km box, "
    "including profiles in a neighbouring region if the box crosses a boundary.\n",
    mdtable(pd.DataFrame({"region": REG_NAMES, "definition (lon E, lat N)": REG_DEF})), f5a,
    "**Profiles per region and decade; median maximum depth; coverage of the region's ocean area by 100 km boxes**\n",
    mdtable(t_reg, dec={"% of profiles": 1, "median max depth (m)": 0, "% reaching 1000 m": 1, "% with S": 1,
                        "ocean area (1e6 km2)": 2, "% area in boxes >=1": 1, "% area in boxes >=10": 1}),
    "**Median maximum depth (m) per region and decade**\n", mdtable(rm.reset_index().rename(columns={"index": "region"}), dec=0),
    f"The Nordic-seas subset (Fram Strait/Greenland Sea + Iceland Sea + Norwegian Sea) used in section 7 has "
    f"{K['nordic_n']:,} profiles; the central Arctic (>80N) has {int(CENTRAL.sum()):,}.\n")
log("section 5 done")

# ============================================================================= SECTION 6: inferred provenance
log("section 6: track linking")
order = np.argsort(tdays, kind="stable")
xs = px[order] / 1e3
ys = py[order] / 1e3
ts = tdays[order]
TMAX, VMAX, DCAP, D0 = 13.0, 45.0, 100.0, 2.0  # d, km/d (0.52 m/s), km, km


def dmax(dt):
    return D0 + np.minimum(VMAX * dt, DCAP)


CH = 2.0
edges = np.arange(ts[0], ts[-1] + CH, CH)
li = np.searchsorted(ts, edges)
KEEP = 6
PASS_A_DT = 3.0  # pass A: short steps (drifting platforms, cost = dt + travel time); pass B: 3-13 d (profiling floats)


def topk(i, j, c, key, keep=KEEP):
    k1, k2 = (i, j) if key == "i" else (j, i)
    o = np.lexsort((c, k1))
    i, j, c, k1 = i[o], j[o], c[o], k1[o]
    first = np.r_[0, np.nonzero(np.diff(k1))[0] + 1]
    rank = np.arange(i.size) - np.repeat(first, np.diff(np.r_[first, i.size]))
    kk = rank < keep
    return i[kk], j[kk], c[kk]


LA, LB = [], []
for k_ in range(len(edges) - 1):
    a, b = li[k_], li[k_ + 1]
    if a == b:
        continue
    c_end = np.searchsorted(ts, edges[k_ + 1] + TMAX)
    pr = np.c_[xs[a:c_end], ys[a:c_end]]
    tl = cKDTree(pr[: b - a])
    trr = cKDTree(pr)
    sd = tl.sparse_distance_matrix(trr, DCAP + D0, output_type="ndarray")
    i = sd["i"].astype(np.int64) + a
    j = sd["j"].astype(np.int64) + a
    d = sd["v"]
    dt = ts[j] - ts[i]
    ok = (dt > 1.0 / 1440) & (dt <= TMAX) & (d <= dmax(dt))
    i, j, d, dt = i[ok], j[ok], d[ok], dt[ok]
    pa = dt <= PASS_A_DT
    if pa.any():
        LA.append(topk(i[pa], j[pa], dt[pa] + d[pa] / VMAX, "i"))
    if (~pa).any():
        LB.append(topk(i[~pa], j[~pa], d[~pa] / dmax(dt[~pa]) + 0.001 * dt[~pa], "i"))
cands = []
for L in (LA, LB):
    I_ = np.concatenate([x[0] for x in L]); J_ = np.concatenate([x[1] for x in L]); C_ = np.concatenate([x[2] for x in L])
    I_, J_, C_ = topk(I_, J_, C_, "j")
    o = np.argsort(C_, kind="stable")
    cands.append((I_[o], J_[o]))
CI = np.concatenate([cands[0][0], cands[1][0]])
CJ = np.concatenate([cands[0][1], cands[1][1]])
log(f"  candidate links: {CI.size:,}")
succ = np.full(N, -1)
has_pred = np.zeros(N, bool)
for i, j in zip(CI.tolist(), CJ.tolist()):
    if succ[i] < 0 and not has_pred[j]:
        succ[i] = j
        has_pred[j] = True
sidx = np.nonzero(succ >= 0)[0]
nlinks = sidx.size
g = coo_matrix((np.ones(nlinks), (sidx, succ[sidx])), shape=(N, N))
ncomp, tlab_s = connected_components(g, directed=False)
trk = np.empty(N, int)
trk[order] = tlab_s  # track label per original profile
log(f"  links: {nlinks:,}, tracks (incl. singletons): {ncomp:,}")

# per-profile step to successor (sorted space)
step_dt_s = np.full(N, np.nan)
step_d_s = np.full(N, np.nan)
step_dt_s[sidx] = ts[succ[sidx]] - ts[sidx]
step_d_s[sidx] = np.hypot(xs[succ[sidx]] - xs[sidx], ys[succ[sidx]] - ys[sidx])
PS = pd.DataFrame({"trk": tlab_s, "t": ts, "x": xs, "y": ys, "dt": step_dt_s, "d": step_d_s,
                   "md": maxdepth[order], "lat": lat[order], "S": hasS[order].astype(float),
                   "fd": firstdepth[order]})
PS["reg"] = (dt_ok := (PS["dt"] >= 8) & (PS["dt"] <= 12.5)).astype(float)
PS["reg"] = np.where(PS["dt"].notna(), PS["reg"], np.nan)
PS["short"] = np.where(PS["dt"].notna(), ((PS["dt"] > 0.02) & (PS["dt"] <= 3)).astype(float), np.nan)
PS["v"] = PS["d"] / PS["dt"]  # km/d
gb = PS.groupby("trk")
TS_ = pd.DataFrame({"n": gb.size(), "t0": gb["t"].min(), "t1": gb["t"].max(), "lat": gb["lat"].median(),
                    "med_maxdepth": gb["md"].median(), "frac_S": gb["S"].mean(), "med_dt": gb["dt"].median(),
                    "frac_regular_8_12d": gb["reg"].mean(), "frac_short_le3d": gb["short"].mean(),
                    "med_speed_km_d": gb["v"].median(), "path_km": gb["d"].sum(),
                    "xr": gb["x"].agg(lambda s: s.max() - s.min()), "yr": gb["y"].agg(lambda s: s.max() - s.min())})
TS_["duration_d"] = TS_["t1"] - TS_["t0"]
TS_["extent_km"] = np.hypot(TS_["xr"], TS_["yr"])
TS_["med_speed_cm_s"] = TS_["med_speed_km_d"] * 1e5 / 86400.0
TS_["median_lon"] = PS.assign(lon=lon180[order]).groupby("trk")["lon"].median()
TS_ = TS_.drop(columns=["xr", "yr"])
n_ = TS_["n"]
a_ = TS_["frac_short_le3d"].fillna(0)
b_ = np.clip((TS_["lat"] - 66) / 6, 0, 1)
md_ = TS_["med_maxdepth"]
c_ = np.where((md_ >= 100) & (md_ <= 1000), 1.0, np.where(((md_ >= 50) & (md_ < 100)) | ((md_ > 1000) & (md_ <= 1500)), 0.5, 0.0))
d_ = (TS_["extent_km"] >= 20).astype(float)
TS_["score_ITP"] = 0.35 * a_ + 0.25 * b_ + 0.25 * c_ + 0.15 * d_
a2 = TS_["frac_regular_8_12d"].fillna(0)
c2 = np.where((md_ >= 900) & (md_ <= 2300), 1.0, np.where(((md_ >= 600) & (md_ < 900)) | ((md_ > 2300) & (md_ <= 4000)), 0.5, 0.0))
e2 = ((TS_["med_dt"] >= 8) & (TS_["med_dt"] <= 12.5)).astype(float)
TS_["score_Argo"] = 0.5 * a2 + 0.3 * c2 + 0.2 * e2
is_fixed = (n_ >= 20) & (TS_["extent_km"] <= 10) & (TS_["duration_d"] >= 5)
is_itp = (n_ >= 10) & (TS_["score_ITP"] >= 0.80) & ~is_fixed
is_argo = (n_ >= 5) & (a2 >= 0.5) & (TS_["score_Argo"] >= 0.70) & ~is_fixed
both = is_itp & is_argo
is_itp &= ~(both & (TS_["score_Argo"] > TS_["score_ITP"]))
is_argo &= ~(both & (TS_["score_ITP"] >= TS_["score_Argo"]))
is_dshal = (n_ >= 10) & (a_ >= 0.7) & (TS_["lat"] >= 70) & (md_ < 100) & ~is_fixed & ~is_itp & ~is_argo
TS_["track_class"] = np.select([is_fixed, is_itp, is_argo, is_dshal],
                               ["fixed-point series", "ITP-like", "Argo-like", "drifting shallow platform"], "none")
TS_["track_score"] = np.select([is_fixed, is_itp, is_argo, is_dshal],
                               [np.minimum(1, 0.5 + 0.5 * np.minimum(1, n_ / 200)), TS_["score_ITP"], TS_["score_Argo"],
                                0.6 * a_], 0.0)
savetab(TS_[TS_["n"] >= 5].reset_index().rename(columns={"trk": "track_id"}), "s6_tracks_n_ge_5")
tc_s = TS_["track_class"].reindex(tlab_s).to_numpy()
tsc_s = TS_["track_score"].reindex(tlab_s).to_numpy()
track_class = np.empty(N, object)
track_score = np.empty(N)
track_class[order] = tc_s
track_score[order] = tsc_s
tn = TS_["n"].reindex(trk).to_numpy()
log(f"  tracks n>=5: {(TS_['n'] >= 5).sum():,}; ITP-like {int(is_itp.sum()):,}, Argo-like {int(is_argo.sum()):,}, fixed {int(is_fixed.sum()):,}")

# neighbourhood counts for the remaining profiles (cruise-like clusters): 50 km and 2 d
log("section 6: cluster counts")
SC = 25.0  # km per day in the scaled space (2 d = 50 km)
tree3 = cKDTree(np.c_[px / 1e3, py / 1e3, tdays * SC])
nn_cnt = tree3.query_ball_point(np.c_[px / 1e3, py / 1e3, tdays * SC], r=50.0, return_length=True) - 1

# Profile-level labels and heuristic scores (not calibrated probabilities).
CLASSES = ["ITP-like", "Argo-like", "ship-CTD-like", "other/unclassified"]
pclass = np.full(N, "other/unclassified", dtype=object)
psub = np.full(N, "unclassified", dtype=object)
pconf = np.zeros(N)
m = track_class == "ITP-like"
ok_md = (maxdepth >= 100) & (maxdepth <= 1000)
pclass[m] = "ITP-like"
psub[m] = "high-latitude candidate track, short steps"
pconf[m] = track_score[m] * np.where(ok_md[m], 1.0, 0.7)
m = track_class == "Argo-like"
ok_md2 = (maxdepth >= 900) & (maxdepth <= 2300)
pclass[m] = "Argo-like"
psub[m] = "regular ~10 d cycle track"
pconf[m] = track_score[m] * np.where(ok_md2[m], 1.0, 0.7)
m = track_class == "fixed-point series"
psub[m] = "low-extent candidate series"
pconf[m] = track_score[m]
m = track_class == "drifting shallow platform"
psub[m] = "drifting high-latitude track, shallow (<100 m median max depth)"
pconf[m] = track_score[m]
rest = (track_class == "none") | (track_class == None)
sc_ship = 0.25 * (nn_cnt >= 2) + 0.25 * (maxdepth >= 300) + 0.25 * np.isin(month, [6, 7, 8, 9, 10]) + 0.25 * hasS
m = rest & hasS & ((nn_cnt >= 2) | (maxdepth >= 300))
pclass[m] = "ship-CTD-like"
psub[m] = np.where(nn_cnt[m] >= 2, "station cluster (>=3 profiles within 50 km and 2 d)", "isolated deep cast")
pconf[m] = sc_ship[m]
m = rest & ~((pclass == "ship-CTD-like"))
psub[m & ~hasS] = "T-only (no salinity)"
psub[m & hasS] = "isolated shallow profile"
pconf[m] = 0.2 + 0.0 * pconf[m]
K["class_counts"] = {c: int((pclass == c).sum()) for c in CLASSES}
log(f"  classes: {K['class_counts']}")


# ---- class tables
log("section 6: class tables")
conf_hi = 0.6
yrs = years_all
rows = []
for c in CLASSES:
    m = pclass == c
    rows.append({"class (INFERRED)": c, "profiles": int(m.sum()), "fraction (%)": 100 * m.mean(),
                 "mean heuristic score": float(pconf[m].mean()), "profiles with score >= 0.6": int((m & (pconf >= conf_hi)).sum()),
                 "% with score >= 0.6": 100 * float((pconf[m] >= conf_hi).mean()),
                 "median max depth (m)": float(np.nanmedian(maxdepth[m])), "% with S": 100 * float(hasS[m].mean()),
                 "first year": int(year[m].min()) if m.any() else np.nan,
                 "median lat (N)": float(np.median(lat[m]))})
t_cls = pd.DataFrame(rows)
savetab(t_cls, "s6_class_summary")
t_sub = (pd.DataFrame({"class": pclass, "subclass": psub, "conf": pconf}).groupby(["class", "subclass"])
         .agg(profiles=("conf", "size"), mean_score=("conf", "mean")).reset_index())
t_sub["fraction (%)"] = 100 * t_sub["profiles"] / N
savetab(t_sub, "s6_subclass_summary")
ct_year = pd.crosstab(year, pd.Categorical(pclass, categories=CLASSES)).reindex(yrs, fill_value=0)
savetab(ct_year.reset_index().rename(columns={"row_0": "year", "index": "year"}), "s6_class_per_year")
ct_reg = pd.crosstab(pd.Categorical([REG_NAMES[i] for i in REG], categories=REG_NAMES), pd.Categorical(pclass, categories=CLASSES))
ct_reg.index.name = "region"
savetab(ct_reg.reset_index(), "s6_class_per_region")
rows = []
for c in CLASSES:
    m = pclass == c
    d = {"class (INFERRED)": c}
    for th in (100, 300, 700, 1000, 1500, 2000, 3000):
        d[f"% >= {th} m"] = 100 * float((maxdepth[m] >= th).mean())
    rows.append(d)
t_clsdepth = pd.DataFrame(rows)
savetab(t_clsdepth, "s6_class_maxdepth_exceedance")

# track statistics per class (tracks with n>=5, plus the fixed-point series)
TSn = TS_[(TS_["n"] >= 5)]
rows = []
for c in ["ITP-like", "Argo-like", "fixed-point series", "drifting shallow platform", "none"]:
    q = TSn[TSn["track_class"] == c]
    rows.append({"track class (INFERRED)": c if c != "none" else "other tracks (n>=5, unclassified)",
                 "tracks": len(q), "profiles in tracks": int(q["n"].sum()), "median profiles/track": float(q["n"].median()),
                 "median duration (d)": float(q["duration_d"].median()), "median step dt (d)": float(q["med_dt"].median()),
                 "median drift speed (cm/s)": float(q["med_speed_cm_s"].median()),
                 "median extent (km)": float(q["extent_km"].median()),
                 "median max depth (m)": float(q["med_maxdepth"].median()), "median lat": float(q["lat"].median())})
t_trk = pd.DataFrame(rows)
savetab(t_trk, "s6_track_statistics_by_class")
n_trk5 = int((TS_["n"] >= 5).sum())
n_in_trk5 = int(TSn["n"].sum())
K.update(n_trk5=n_trk5, frac_in_trk5=n_in_trk5 / N, nlinks=int(nlinks))

# validation
first_itp = int(year[pclass == "ITP-like"].min())
first_argo = int(year[pclass == "Argo-like"].min())
itp_pre04 = int(((pclass == "ITP-like") & (year < 2004)).sum())
argo_pre00 = int(((pclass == "Argo-like") & (year < 2000)).sum())
itp_tr = TS_[TS_["track_class"] == "ITP-like"]
arg_tr = TS_[TS_["track_class"] == "Argo-like"]
frac10 = float(((arg_tr["med_dt"] >= 9) & (arg_tr["med_dt"] <= 11)).mean()) if len(arg_tr) else np.nan
argo_lat = float(np.median(lat[pclass == "Argo-like"]))
argo_nordic_share = float(np.isin(REG[pclass == "Argo-like"], [4, 5, 6, 7, 14]).mean())
itp_lat72 = float((lat[pclass == "ITP-like"] >= 72).mean())
K.update(first_itp=first_itp, first_argo=first_argo, itp_pre04=itp_pre04, argo_pre00=argo_pre00,
         argo_lat=argo_lat, argo_nordic_share=argo_nordic_share, itp_lat72=itp_lat72)
cc = K["class_counts"]
K["frac_class"] = {c: cc[c] / N for c in CLASSES}

# ---- figures
cols = {"ITP-like": "C3", "Argo-like": "C0", "ship-CTD-like": "C2", "other/unclassified": "0.5"}
fig, axs = plt.subplots(1, 2, figsize=(14, 4.5))
bot = np.zeros(len(yrs))
for c in CLASSES:
    axs[0].bar(yrs, ct_year[c].values, bottom=bot, color=cols[c], label=c)
    bot += ct_year[c].values
axs[0].set_xlabel("year"); axs[0].set_ylabel("profiles"); axs[0].legend(); axs[0].set_title("INFERRED class per year (counts)")
bot = np.zeros(len(yrs))
tot = ct_year.sum(axis=1).values.clip(1)
for c in CLASSES:
    axs[1].bar(yrs, 100 * ct_year[c].values / tot, bottom=bot, color=cols[c], label=c)
    bot += 100 * ct_year[c].values / tot
axs[1].set_xlabel("year"); axs[1].set_ylabel("% of the year's profiles"); axs[1].set_title("INFERRED class per year (share)")
axs[1].axvline(2003.5, color="k", ls=":")
f6a = savefig(fig, "s6a_class_per_year",
             "Unvalidated sampling-pattern labels by year; neither counts nor shares identify platform contributions")

fig = plt.figure(figsize=(13, 13))
axl = []
for k_, c in enumerate(CLASSES):
    m = pclass == c
    ax = polar_ax(fig, 2, 2, k_ + 1, f"{c} (INFERRED): {int(m.sum()):,} profiles")
    mm = box_pcolor(ax, grid2(bcount(m)), norm=LogNorm(1, cnt.max()), cmap="viridis")
    axl.append(ax)
cbar(fig, mm, axl, "profiles per 100 km box (log)")
f6b = savefig(fig, "s6b_class_maps", "Spatial distribution of heuristic labels, with a common logarithmic count scale")

fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
for c in CLASSES:
    v = np.sort(maxdepth[(pclass == c) & hasT])
    axs[0].plot(v, 100 * (1 - np.arange(v.size) / v.size), color=cols[c], label=c)
    h = np.bincount(lastT[(pclass == c) & hasT], minlength=NZ)
    axs[1].plot(depth, 100 * h / h.sum(), color=cols[c], label=c)
axs[0].set_xscale("symlog", linthresh=50); axs[0].set_xlabel("max valid depth z (m)"); axs[0].set_ylabel("% of class with max depth >= z")
axs[0].legend(); axs[0].set_title("maximum depth per INFERRED class (exceedance)")
axs[1].set_xscale("symlog", linthresh=50); axs[1].set_yscale("log"); axs[1].set_xlabel("deepest valid level (m)")
axs[1].set_ylabel("% of class profiles per level"); axs[1].set_title("distribution of deepest level")
f6c = savefig(fig, "s6c_class_maxdepth", "Depth reach by heuristic label; depth itself is an input to these rules")

sdt = step_dt_s[sidx]
fig, axs = plt.subplots(1, 3, figsize=(16, 4.5))
axs[0].hist(np.log10(sdt), bins=120, color="C1")
axs[0].set_xlabel("log10 time step between linked profiles (days)"); axs[0].set_ylabel("links")
axs[0].set_title("all linked pairs (pass A + pass B)")
for v_ in (1.0 / 24, 0.25, 1, 10):
    axs[0].axvline(np.log10(v_), color="k", lw=0.5, ls=":")
if len(arg_tr):
    axs[1].hist(arg_tr["med_dt"].clip(0, 14), bins=np.arange(0, 14.01, 0.25), color="C0")
axs[1].set_xlabel("median time step of Argo-like tracks (d)"); axs[1].set_title("regular-cycle signal (expect ~10 d)")
ci_ = ct_year["ITP-like"].values; ca_ = ct_year["Argo-like"].values
axs[2].semilogy(yrs, ci_.clip(0.5), "C3o-", label="ITP-like"); axs[2].semilogy(yrs, ca_.clip(0.5), "C0s-", label="Argo-like")
axs[2].axvline(2004, color="C3", ls=":"); axs[2].axvline(2000, color="C0", ls=":")
axs[2].set_xlabel("year"); axs[2].set_ylabel("profiles"); axs[2].legend(); axs[2].set_title("first appearance (year not used by classifier)")
f6d = savefig(fig, "s6d_validation",
             "Linked-profile timing and historical benchmark years; these are consistency checks, not independent platform validation")

# ---- near duplicates and identical copies
log("section 6: duplicates")
o_t = np.argsort(tsec, kind="stable")
tt = tsec[o_t].astype(float)
prs = cKDTree(tt[:, None]).query_pairs(21600.0, output_type="ndarray")
pa_ = o_t[prs[:, 0]]
pb_ = o_t[prs[:, 1]]
del prs
dts = np.abs(tsec[pa_] - tsec[pb_]).astype(float)
dla = np.abs(lat[pa_] - lat[pb_])
dlo = np.abs(lon[pa_] - lon[pb_])
dlo_w = np.minimum(dlo, 360 - dlo)
la1, la2 = np.deg2rad(lat[pa_]), np.deg2rad(lat[pb_])
hav = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin(np.deg2rad(dlo_w) / 2) ** 2
dkm = 2 * 6371.0 * np.arcsin(np.sqrt(np.clip(hav, 0, 1)))
up6 = (dts < 21600) & (dla < 0.1) & (dlo < 0.1)
up5 = (dts <= 18000) & (dla < 0.1) & (dlo < 0.1)
up6w = (dts < 21600) & (dla < 0.1) & (dlo_w < 0.1)
p1km1h = (dkm <= 1.0) & (dts <= 3600)
p1km1h_dt0 = p1km1h & (dts == 0)
same_xyt = (dts == 0) & (dla == 0) & (dlo == 0)


def nprof(mask):
    return int(np.unique(np.concatenate([pa_[mask], pb_[mask]])).size) if mask.any() else 0


# exact copies
hashT = np.empty(N, dtype="S8")
Tn = np.nan_to_num(T, nan=-999.0)
for i in range(N):
    hashT[i] = hashlib.blake2b(Tn[i].tobytes(), digest_size=8).digest()
del Tn
hashT = hashT.view("<u8")
dfh = pd.DataFrame({"h": hashT, "n": nlevT, "lat": lat, "lon": lon, "t": tsec, "y": year, "cls": pclass})
rows = []
for minlev in (1, 5):
    q = dfh[(dfh["n"] >= minlev)]
    gsz = q.groupby("h")["n"].transform("size")
    qq = q[gsz >= 2]
    ng = qq["h"].nunique()
    # copies sharing exact position and time vs other
    gg = qq.groupby("h").agg(size=("n", "size"), nlat=("lat", "nunique"), nlon=("lon", "nunique"), nt=("t", "nunique"))
    same_all = int(((gg["nlat"] == 1) & (gg["nlon"] == 1) & (gg["nt"] == 1)).sum())
    same_pos_difft = int(((gg["nlat"] == 1) & (gg["nlon"] == 1) & (gg["nt"] > 1)).sum())
    diff_pos = int(((gg["nlat"] > 1) | (gg["nlon"] > 1)).sum())
    rows.append({"min valid levels": minlev, "profiles considered": len(q), "identical-T groups (size>=2)": ng,
                 "profiles in such groups": len(qq), "profiles beyond one per T group": len(qq) - ng,
                 "largest group": int(gg["size"].max()) if ng else 0,
                 "groups: same lat,lon,time": same_all, "groups: same lat,lon, different time": same_pos_difft,
                 "groups: different position": diff_pos})
t_ident = pd.DataFrame(rows)
savetab(t_ident, "s6_identical_T_summary")
q = dfh[dfh["n"] >= 1]
gsz = q.groupby("h")["n"].transform("size")
qq = q[gsz >= 2].copy()
for members in qq.groupby("h").groups.values():
    reference = T[members[0]]
    if not all(np.array_equal(reference, T[i], equal_nan=True) for i in members[1:]):
        raise RuntimeError("Temperature hash collision: identical-T counts cannot be reported")
gg = qq.groupby("h").agg(size=("n", "size"), n_levels=("n", "max"), lat_min=("lat", "min"), lat_max=("lat", "max"),
                         lon_min=("lon", "min"), lon_max=("lon", "max"), t_first=("t", "min"), t_last=("t", "max"),
                         year_first=("y", "min")).sort_values("size", ascending=False)
gg["t_first"] = pd.to_datetime(gg["t_first"], unit="s")
gg["t_last"] = pd.to_datetime(gg["t_last"], unit="s")
t_identtop = gg.head(25).reset_index(drop=True)
savetab(t_identtop, "s6_identical_T_top_groups")
ident_pair = np.zeros(pa_.size, bool)
ident_pair = (hashT[pa_] == hashT[pb_]) & (nlevT[pa_] >= 1)
# S identical as well?
rows = [
    ["pairs within 5 h (all, input to the tests)", int((dts <= 18000).sum()), nprof(dts <= 18000)],
    ["pairs |dlat|<0.1, |dlon|<0.1, dt<6 h  (criterion as coded in upstream duplicate.m)", int(up6.sum()), nprof(up6)],
    ["pairs |dlat|<0.1, |dlon|<0.1, dt<=5 h  (5 h variant motivated by the manuscript)", int(up5.sum()), nprof(up5)],
    ["  same but with dateline-aware dlon", int(up6w.sum()), nprof(up6w)],
    ["pairs within 1 km (haversine) and 1 h", int(p1km1h.sum()), nprof(p1km1h)],
    ["pairs within 1 km and identical time stamp", int(p1km1h_dt0.sum()), nprof(p1km1h_dt0)],
    ["pairs with identical lat, lon and time", int(same_xyt.sum()), nprof(same_xyt)],
    ["pairs (6 h/0.1 deg criterion) with identical T arrays", int((up6 & ident_pair).sum()), nprof(up6 & ident_pair)],
    ["pairs (1 km/1 h) with identical T arrays", int((p1km1h & ident_pair).sum()), nprof(p1km1h & ident_pair)],
]
t_dup = pd.DataFrame(rows, columns=["criterion", "pairs", "profiles involved"])
t_dup["% of all profiles involved"] = 100 * t_dup["profiles involved"] / N
savetab(t_dup, "s6_near_duplicate_pairs")
# which inferred classes make up the near-duplicate pairs (6 h / 0.1 deg)
cls_a = pclass[pa_[up6]]
cls_b = pclass[pb_[up6]]
t_dupcls = pd.crosstab(cls_a, cls_b)
t_dupcls = t_dupcls.reindex(index=CLASSES, columns=CLASSES, fill_value=0)
t_dupcls.index.name = "class A / class B"
savetab(t_dupcls.reset_index(), "s6_near_duplicate_pairs_by_class")
yr_dup = pd.Series(year[pa_[up6]]).value_counts().sort_index()
savetab(yr_dup.rename("pairs").reset_index().rename(columns={"index": "year", "row_0": "year"}), "s6_near_duplicate_pairs_per_year")
# clusters
gup = coo_matrix((np.ones(int(up6.sum())), (pa_[up6], pb_[up6])), shape=(N, N))
ncl, lab = connected_components(gup, directed=False)
csz = np.bincount(lab)
big = np.nonzero(csz >= 2)[0]
sz = csz[big]
edges_ = [(2, 2), (3, 3), (4, 4), (5, 9), (10, 49), (50, 10 ** 7)]
t_clsz = pd.DataFrame([{"cluster size": f"{a}" if a == b else (f">={a}" if b > 10 ** 6 else f"{a}-{b}"),
                        "clusters": int(((sz >= a) & (sz <= b)).sum()),
                        "profiles": int(sz[(sz >= a) & (sz <= b)].sum())} for a, b in edges_])
savetab(t_clsz, "s6_near_duplicate_cluster_sizes")
top_c = big[np.argsort(sz)[::-1][:15]]
rows = []
for cidx in top_c:
    mem = np.nonzero(lab == cidx)[0]
    hh = pd.Series(hashT[mem]).value_counts()
    rows.append({"cluster size": len(mem), "centre lat": float(lat[mem].mean()), "centre lon": float(lon180[mem].mean()),
                 "start": str(np.datetime64(int(tsec[mem].min()), "s")), "end": str(np.datetime64(int(tsec[mem].max()), "s")),
                 "lat span (deg)": float(np.ptp(lat[mem])), "profiles with an identical-T match": int((hh[hh >= 2]).sum()),
                 "dominant inferred class": pd.Series(pclass[mem]).value_counts().index[0],
                 "dominant region": REG_NAMES[int(pd.Series(REG[mem]).value_counts().index[0])]})
t_cltop = pd.DataFrame(rows)
savetab(t_cltop, "s6_near_duplicate_top_clusters")
K.update(dup6_pairs=int(up6.sum()), dup6_prof=nprof(up6), dup5_pairs=int(up5.sum()), dup1km_pairs=int(p1km1h.sum()),
         dup1km_prof=nprof(p1km1h), dup_ident_pairs=int((up6 & ident_pair).sum()),
         ident_redundant1=int(t_ident.iloc[0]["profiles beyond one per T group"]),
         ident_redundant5=int(t_ident.iloc[1]["profiles beyond one per T group"]),
         ident_groups5=int(t_ident.iloc[1]["identical-T groups (size>=2)"]),
         ncl2=int(len(big)), maxcl=int(sz.max()) if len(sz) else 0)
fig, axs = plt.subplots(1, 3, figsize=(16, 4.3))
axs[0].hist(dts[up6] / 3600, bins=48, color="C4"); axs[0].set_yscale("log")
axs[0].set_xlabel("time difference of near-duplicate pairs (h)"); axs[0].set_ylabel("pairs"); axs[0].set_title("pairs meeting the 0.1 deg / 6 h test")
axs[1].hist(np.log10(dkm[up6] + 1e-3), bins=60, color="C5"); axs[1].set_xlabel("log10 distance (km)"); axs[1].set_title("separation of those pairs")
axs[2].bar(range(len(t_clsz)), t_clsz["clusters"], color="C6"); axs[2].set_yscale("log")
axs[2].set_xticks(range(len(t_clsz))); axs[2].set_xticklabels(t_clsz["cluster size"]); axs[2].set_xlabel("cluster size (profiles chained by the 0.1 deg/6 h test)")
axs[2].set_ylabel("clusters")
f6e = savefig(fig, "s6e_near_duplicates",
             "Space-time proximity counts and connected-component sizes; repeated sampling and duplicate records are not distinguished")
del pa_, pb_, dts, dla, dlo, dlo_w, dkm, hav, la1, la2

# ---- markdown
LIM = (
    "**Limits of interpretation.** No profile-level platform or source identifiers are available, and no labelled "
    "validation sample was used. The rules can connect unrelated records or split a real platform's series, especially "
    "where several sampling systems operate nearby. Irregular cycles, gaps and isolated records can prevent an Argo-like "
    "assignment. *Ship-CTD-like* is a residual label for salinity-bearing clustered or sufficiently deep profiles: it may "
    "include other observations such as moorings, gliders, bottle data or fragments of profiler series. "
    "Neither this label nor ITP-like and Argo-like establishes instrument type. Scores in [0,1] quantify the stated "
    "rules only; they are not probabilities, confidence intervals or estimates of classification accuracy. "
    "No precision, recall or source-specific quality ranking can be established from these arrays."
)
add(6, "## 6. Heuristic sampling-pattern labels and proximity screening\n",
    "> **Labels are inferred sampling patterns, not identified platforms.** All profile membership and class counts "
    "below come from this script's unvalidated rules, not the source identifiers of the original observing systems.\n",
    "**Externally reported provenance.** The dataset and its compilation are described by "
    "[Li, Wu and Wang (2025)](https://doi.org/10.1038/s41597-025-05855-3). An earlier repository review of the released "
    "`duplicate.m` reports the source numbering 1 ITP, 2 Argo, 3 BOL (Bolin Centre), 4 NABOS, "
    "5 CNarc (Xuelong/Xuelong-2), 6 WOD18 and 7 CMEMS. Those identifiers are not fields in the scatter arrays used here; "
    "this source inventory must not be confused with the heuristic membership below. "
    "This diagnostic run does not re-execute or independently validate the upstream compilation.\n",
    "### 6.1 Method\n",
    f"1. **Track linking.** Profiles are sorted by time; for every profile, candidate successors within {TMAX:.0f} d are searched "
    f"within {DCAP + D0:.0f} km (KD-tree in EPSG:6931 km). A link i->j requires dt in (1 min, {TMAX:.0f} d] and distance <= "
    f"{D0:.0f} km + min({VMAX:.0f} km/d x dt, {DCAP:.0f} km). These are projected distances, not geodesic distances; "
    f"the additive allowance means this is not a strict speed bound. Links are assigned greedily, one successor and one "
    f"predecessor per profile, in two passes: pass A for steps <= {PASS_A_DT:.0f} d (cost = dt + distance/{VMAX:.0f} km d-1), then pass B for "
    f"{PASS_A_DT:.0f}-{TMAX:.0f} d (cost = distance / allowed distance + 0.001 x dt in days). "
    f"Each pass retains at most {KEEP} lowest-cost outgoing and incoming candidates per profile. "
    f"Linked chains are candidate tracks, not recovered platform trajectories: {int(nlinks):,} links, "
    f"{n_trk5:,} tracks with >= 5 profiles containing {pct(n_in_trk5, N, 1)} of all profiles.\n"
    "2. **Track statistics** (number of profiles, duration, median time step, fraction of steps in 8-12.5 d and <= 3 d, drift speed, extent, "
    "median maximum depth, median latitude) are in `tables/s6_tracks_n_ge_5.csv`. "
    "Extent is the diagonal of the track's projected bounding box; step statistics exclude terminal profiles with no successor.\n"
    "3. **Track classes.** *ITP-like*: >= 10 profiles, score >= 0.80 where score = 0.35 x fraction of steps <= 3 d + 0.25 x latitude factor "
    "(linear from 0 at 66N to 1 at >=72N) + 0.25 x depth factor + 0.15 x [extent >=20 km]. "
    "The depth factor is 1 for median maximum depth 100-1000 m, 0.5 for 50-<100 m or >1000-1500 m, and 0 otherwise. "
    "The short-step fraction includes steps >0.02 d and <=3 d. *Argo-like*: >=5 profiles, >=50% of steps "
    "8-12.5 d, score >=0.70 with score = 0.5 x fraction of regular steps + 0.3 x depth factor + 0.2 x [median step 8-12.5 d]. "
    "Its depth factor is 1 for 900-2300 m, 0.5 for 600-<900 m or >2300-4000 m, and 0 otherwise. "
    "*Fixed-point series*: >= 20 profiles, extent <= 10 km, duration >= 5 d (reported under 'other'). *Drifting shallow platform*: "
    ">=10 profiles, >=70% short steps, median lat >=70N, median max depth <100 m (under 'other'). "
    "Fixed-point rules take priority. If both ITP-like and Argo-like pass, the larger score wins, with ties assigned ITP-like; "
    "shallow-drift rules apply only to tracks not assigned those classes.\n"
    "4. **Profile labels and scores.** ITP-like and Argo-like track members inherit the label and score, multiplied by 0.7 "
    "if the individual depth reach lies outside 100-1000 m or 900-2300 m respectively. Fixed-point and shallow-drift "
    "members remain *other/unclassified*. Fixed-point score = 0.5 + 0.5 x min(1, track length/200); "
    "shallow-drift score = 0.6 x short-step fraction. Among remaining profiles, those with salinity and either "
    ">=2 neighbours in the search neighbourhood or depth reach >=300 m are *ship-CTD-like*: "
    "score = 0.25 x [>=2 neighbours] + 0.25 x [>=300 m] + 0.25 x [Jun-Oct] + 0.25 x [has S]. "
    "The neighbourhood is a radius-50 sphere in (projected x km, projected y km, time days x25): "
    "spatial distance squared + (25 x time separation) squared <=50 squared. Thus 50 km and 2 d are "
    "separate maxima, not a rectangular space-time search. Other residual profiles receive score 0.2. "
    "Square brackets denote indicator values 0 or 1; no year or sea-ice field is used.\n",
    LIM + "\n",
    "### 6.2 Results: classes\n",
    mdtable(t_cls, dec={"fraction (%)": 1, "mean heuristic score": 2, "% with score >= 0.6": 1, "median max depth (m)": 0, "% with S": 1,
                        "median lat (N)": 1, "first year": 0}),
    "**Rule-based subclasses** (the cluster wording uses the ellipsoidal neighbourhood defined above):\n",
    mdtable(t_sub, dec={"mean_score": 2, "fraction (%)": 1}),
    f6a,
    "**Track statistics per inferred class** (tracks with >= 5 profiles)\n", mdtable(t_trk, dec=1),
    f6b, "**Maximum depth by inferred class**\n", mdtable(t_clsdepth, dec=1), f6c,
    "**Inferred class per region** (profiles; `tables/s6_class_per_region.csv`)\n",
    mdtable(ct_reg.reset_index().rename(columns={"index": "region"})),
    "### 6.3 Internal consistency and historical benchmarks\n", f6d,
    f"* The classifier does not use the year. ITP-like profiles first appear in **{first_itp}**, with "
    f"{itp_pre04:,} before the historical 2004 ITP deployment benchmark. Argo-like profiles first appear in **{first_argo}**, "
    f"with {argo_pre00:,} before the 2000 Argo benchmark. These early labels show why a sampling-pattern name "
    f"cannot be interpreted as an identified platform.\n"
    f"* {pct(frac10, 1)} of Argo-like tracks have a median step of 9-11 d (regularity is built into the rule, so this is not independent); "
    f"the all-link histogram also depends on the candidate search and greedy matching. Neither timing panel independently "
    f"validates the labels. Argo-like profiles have median latitude {argo_lat:.1f}N and {pct(argo_nordic_share, 1)} lie in the "
    f"Nordic seas / Barents / subpolar North Atlantic reporting regions; {pct(itp_lat72, 1)} of ITP-like profiles are north of 72N. "
    f"Latitude is itself an input to the ITP-like rule.\n",
    "### 6.4 Space-time proximity and identical temperature arrays\n",
    "**Reported upstream rule.** The earlier review of `duplicate.m` describes a comparison of sources 2-7 with previously "
    "accepted profiles using |dlat|<0.1 deg, |dlon|<0.1 deg and dt<6 h, without dateline wrapping; "
    "the publication describes a 5 h window. That review reports that source 1 (ITP) and pairs within one input file "
    "were not compared by this routine. These statements concern the released workflow, not the identity of individual "
    "records in this report.\n",
    "**Observed proximity counts.** The tests below operate on all published profiles, including repeated observations "
    "from the same sampling series. Angular windows are not constant physical distances, particularly at high latitudes. "
    "Each unordered pair is counted once; a profile can occur in many pairs. "
    "Passing a proximity test is a review flag, not proof of erroneous duplication or a reason to delete a record.\n",
    mdtable(t_dup, dec=1),
    f"{K['dup6_pairs']:,} pairs involve {K['dup6_prof']:,} distinct profiles ({pct(K['dup6_prof'], N, 1)}) under the 0.1 deg / 6 h test; "
    f"{K['dup1km_pairs']:,} pairs ({K['dup1km_prof']:,} profiles) lie within 1 km and 1 h. Of the 6 h pairs, {K['dup_ident_pairs']:,} have "
    f"matching temperature values and missing-level patterns. Salinity, pressure, density and original observation identifiers "
    f"are not part of that equality test.\n",
    "**Proximity pairs by inferred label of the two members** (A is the earlier member in stable time order; "
    "the matrix is not symmetrised)\n", mdtable(t_dupcls.reset_index()),
    "**Clusters** (chains of profiles linked by the 0.1 deg / 6 h test)\n", mdtable(t_clsz),
    f"There are {K['ncl2']:,} connected components of >=2 profiles; the largest has {K['maxcl']:,} profiles. "
    f"Connections are transitive: two members of a component may be much farther apart than the pairwise thresholds. "
    f"These components are not duplicate-event counts or validated platform tracks. Largest components:\n",
    mdtable(t_cltop, dec=2), f6e,
    "**Identical temperature arrays** (stored values and missing-level pattern; hash matches are checked against the arrays)\n",
    mdtable(t_ident),
    f"The number of profiles beyond one per temperature group is {K['ident_redundant1']:,} for all lengths "
    f"and {K['ident_redundant5']:,} for profiles with >=5 valid levels (in {K['ident_groups5']:,} groups). "
    f"These are not counts of confirmed duplicate observations. Identical short or rounded arrays can occur in "
    f"distinct observations; equality of temperature alone does not establish equality of the complete record. "
    f"Source identifiers and original records would be needed to resolve their meaning. "
    f"Largest groups: `tables/s6_identical_T_top_groups.csv`.\n")
log("section 6 done")

# ============================================================================= SECTION 7: data quality
log("section 7: gsw quality checks (all profiles, chunked)")
BANDS = [("0-50 m", 0, 50), ("50-200 m", 50, 200), ("200-700 m", 200, 700), ("700-1500 m", 700, 1500),
         ("1500-3000 m", 1500, 3000), ("3000-5500 m", 3000, 5501)]
fig, axs = plt.subplots(2, 6, figsize=(17, 6))
rows = []
for k, (lab, a, b) in enumerate(BANDS):
    lv = (depth >= a) & (depth < b)
    tv = T[:, lv][vT[:, lv]]
    sv = S[:, lv][vS[:, lv]]
    axs[0, k].hist(tv, bins=100, color="C0")
    axs[0, k].set_yscale("log")
    axs[0, k].set_title(lab)
    axs[0, k].set_xlabel("temperature (degC)")
    axs[1, k].hist(sv, bins=100, color="C1")
    axs[1, k].set_yscale("log")
    axs[1, k].set_xlabel("salinity (PSS-78)")
    rows.append({"depth band": lab, "T obs": tv.size, "T mean": tv.mean(), "T std": tv.std(), "T min": tv.min(),
                 "T p1": np.percentile(tv, 1), "T p99": np.percentile(tv, 99), "T max": tv.max(),
                 "S obs": sv.size, "S mean": sv.mean(), "S std": sv.std(), "S min": sv.min(),
                 "S p1": np.percentile(sv, 1), "S p99": np.percentile(sv, 99), "S max": sv.max()})
axs[0, 0].set_ylabel("count (log)")
axs[1, 0].set_ylabel("count (log)")
fig.suptitle("Histograms of T (top) and S (bottom) per depth band, all profiles")
f7a = savefig(fig, "s7a_hist_TS_per_band", "T and S histograms per depth band")
t_band = pd.DataFrame(rows)
savetab(t_band, "s7_TS_stats_per_depth_band")

THR_INV = [0.001, 0.01, 0.03, 0.1]
inv_cnt = np.zeros((len(THR_INV), NZ - 1), np.int64)
pair_cnt = np.zeros(NZ - 1, np.int64)
prof_inv = np.zeros((len(THR_INV), N), bool)
prof_has_pair = np.zeros(N, bool)
frz_within = np.zeros(NZ, np.int64)
frz_below = np.zeros(NZ, np.int64)
frz_valid = np.zeros(NZ, np.int64)
SAE = np.linspace(30.0, 38.0, 321)
CTE = np.linspace(-2.2, 10.0, 489)
H_TS = np.zeros((len(SAE) - 1, len(CTE) - 1), np.int64)
CH = 20000
for i0 in range(0, N, CH):
    sl = slice(i0, min(N, i0 + CH))
    Tc, Sc, Pc = T[sl].astype("float64"), S[sl].astype("float64"), P[sl].astype("float64")
    SA = gsw.SA_from_SP(Sc, Pc, lon180[sl, None], lat[sl, None])
    CT = gsw.CT_from_t(SA, Tc, Pc)
    pm = 0.5 * (Pc[:, :-1] + Pc[:, 1:])
    d_ = gsw.rho(SA[:, :-1], CT[:, :-1], pm) - gsw.rho(SA[:, 1:], CT[:, 1:], pm)  # >0: denser water above
    ok = np.isfinite(d_)
    pair_cnt += ok.sum(0)
    prof_has_pair[sl] = ok.any(1)
    for j, th in enumerate(THR_INV):
        m = ok & (d_ > th)
        inv_cnt[j] += m.sum(0)
        prof_inv[j, sl] = m.any(1)
    tf = gsw.t_freezing(SA, Pc, 0.0)
    dd = Tc - tf
    okf = np.isfinite(dd)
    frz_valid += okf.sum(0)
    frz_within += (okf & (np.abs(dd) <= 0.01)).sum(0)
    frz_below += (okf & (dd < -0.01)).sum(0)
    g = np.isfinite(SA) & np.isfinite(CT)
    H_TS += np.histogram2d(SA[g], CT[g], bins=[SAE, CTE])[0].astype(np.int64)
    del Tc, Sc, Pc, SA, CT, pm, d_, tf, dd
log("  gsw done")
n_pairs = int(pair_cnt.sum())
t_inv = pd.DataFrame({"threshold rho(k)-rho(k+1) (kg/m3)": THR_INV,
                      "adjacent-level pairs flagged": inv_cnt.sum(1),
                      "% of pairs": 100 * inv_cnt.sum(1) / n_pairs,
                      "profiles with >=1 flagged pair": prof_inv.sum(1),
                      "% of profiles with S": 100 * prof_inv.sum(1) / hasS.sum(),
                      "profiles with assessable pairs": int(prof_has_pair.sum()),
                      "% of assessable profiles": 100 * prof_inv.sum(1) / prof_has_pair.sum()})
savetab(t_inv, "s7_static_instability")
bandsI = [("0-50 m", 0, 50), ("50-200 m", 50, 200), ("200-700 m", 200, 700), ("700-1500 m", 700, 1500), ("1500-5500 m", 1500, 5500)]
rows = []
for lab, a, b in bandsI:
    lv = np.where((depth[:-1] >= a) & (depth[:-1] < b))[0]
    rows.append({"upper level in band": lab, "pairs": int(pair_cnt[lv].sum()),
                 **{f"% pairs > {th}": 100 * inv_cnt[j, lv].sum() / max(1, pair_cnt[lv].sum()) for j, th in enumerate(THR_INV)}})
t_invb = pd.DataFrame(rows)
savetab(t_invb, "s7_static_instability_per_band")
t_invcls = pd.DataFrame({"class (INFERRED)": CLASSES,
                         "profiles with S": [int((hasS & (pclass == c)).sum()) for c in CLASSES],
                         "profiles with assessable pairs": [int((prof_has_pair & (pclass == c)).sum()) for c in CLASSES],
                         "profiles flagged (>0.03 kg/m3)": [int((prof_inv[2] & (pclass == c)).sum()) for c in CLASSES],
                         "% of assessable profiles flagged": [
                             100 * (prof_inv[2] & (pclass == c)).sum() / max(1, int(((pclass == c) & prof_has_pair).sum()))
                             for c in CLASSES]})
savetab(t_invcls, "s7_static_instability_per_class")
fig, axs = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
for j, th in enumerate(THR_INV):
    axs[0].plot(100 * inv_cnt[j] / np.maximum(pair_cnt, 1) + 1e-4, depth[:-1], label=f"> {th} kg/m3")
axs[0].set_xscale("log")
axs[0].invert_yaxis()
axs[0].set_xlabel("% of assessable adjacent-level pairs flagged (+1e-4)")
axs[0].set_ylabel("depth of upper level (m)")
axs[0].set_title("(a) common-pressure density inversion screening")
axs[0].legend()
axs[1].plot(100 * frz_within / np.maximum(frz_valid, 1), depth, label="|T - T_freeze| <= 0.01")
axs[1].plot(100 * frz_below / np.maximum(frz_valid, 1), depth, label="T < T_freeze - 0.01")
axs[1].set_xscale("symlog", linthresh=0.01)
axs[1].set_xlabel("% of levels with assessable T, S and pressure")
axs[1].set_title("(b) values at / below the freezing point")
axs[1].legend()
f7b = savefig(fig, "s7b_instability_freezing",
             "Density-inversion and freezing-point screening at standard depths; flags are not validated error diagnoses")
fz = pd.DataFrame({"level": np.arange(NZ), "depth_m": depth, "valid_T_with_S": frz_valid, "within_0.01_of_Tfreeze": frz_within,
                   "below_Tfreeze_minus_0.01": frz_below})
savetab(fz, "s7_freezing_point_per_level")
n_within, n_below, n_fv = int(frz_within.sum()), int(frz_below.sum()), int(frz_valid.sum())
n_Tm2 = int((T == -2.0).sum())

# T-S diagram
fig, ax = plt.subplots(figsize=(8, 7))
Hm = np.ma.masked_equal(H_TS, 0)
mm = ax.pcolormesh(SAE, CTE, Hm.T, norm=LogNorm(vmin=1, vmax=H_TS.max()), cmap="viridis", rasterized=True)
cs_s, cs_t = np.meshgrid(np.linspace(30, 38, 120), np.linspace(-2.2, 10, 120))
sg = gsw.sigma0(cs_s, cs_t)
cc_ = ax.contour(cs_s, cs_t, sg, levels=np.arange(22, 29, 0.5), colors="k", linewidths=0.5)
ax.clabel(cc_, fmt="%.1f", fontsize=7)
ax.set_xlim(SAE[0], SAE[-1])
ax.set_ylim(-2.2, 10)
ax.set_xlabel("Absolute salinity SA (g/kg)")
ax.set_ylabel("Conservative temperature CT (degC)")
ax.set_title("SA-CT histogram in the displayed window; contours: sigma0 (kg/m3)")
cb = fig.colorbar(mm, ax=ax)
cb.set_label("number of observations")
f7c = savefig(fig, "s7c_TS_diagram",
             "Absolute salinity–Conservative Temperature histogram (log counts) within SA 30–38 g/kg and CT -2.2–10 degC; contours show sigma0")
K["TS_n"] = int(H_TS.sum())

# spikes
dz = np.diff(depth)
gT = np.abs(np.diff(T, axis=1)) / dz
gS = np.abs(np.diff(S, axis=1)) / dz
nP_T = int((~np.isnan(gT)).sum())
nP_S = int((~np.isnan(gS)).sum())
spT, spS = gT > 0.3, gS > 0.2
K["spike_T_prof"], K["spike_S_prof"] = int(spT.any(1).sum()), int(spS.any(1).sum())
t_spk = pd.DataFrame([
    {"test": "|dT/dz| > 0.3 degC/m between adjacent valid levels", "pairs flagged": int(spT.sum()), "% of pairs": 100 * spT.sum() / nP_T,
     "profiles flagged": K["spike_T_prof"], "% of profiles (with T)": 100 * K["spike_T_prof"] / hasT.sum()},
    {"test": "|dS/dz| > 0.2 per m between adjacent valid levels", "pairs flagged": int(spS.sum()), "% of pairs": 100 * spS.sum() / nP_S,
     "profiles flagged": K["spike_S_prof"], "% of profiles (with T)": 100 * K["spike_S_prof"] / hasS.sum()},
])
t_spk = t_spk.rename(columns={"% of profiles (with T)": "% of profiles (with T resp. S)"})
savetab(t_spk, "s7_spikes")
spk_depth = pd.DataFrame({"upper_depth_m": depth[:-1], "T_pairs_flagged": spT.sum(0), "S_pairs_flagged": spS.sum(0)})
savetab(spk_depth, "s7_spikes_per_level")
del gT, gS, spT, spS

# mean/std per level and water masses
fig, axs = plt.subplots(1, 4, figsize=(15, 6), sharey=True)
lvl_tab = {"depth_m": depth[:NL57]}
for nm, msk, col, key in (("Central Arctic (>80N)", CENTRAL, "C3", "central"), ("Nordic seas", NORDIC, "C0", "nordic")):
    for v, arr0, ai in (("T", T, 0), ("S", S, 2)):
        arr = arr0[msk][:, :NL57]
        mu, sd = level_mean_std(arr)
        lvl_tab[f"{key}_{v}_n"] = np.sum(~np.isnan(arr), 0)
        lvl_tab[f"{key}_{v}_mean"] = mu
        lvl_tab[f"{key}_{v}_std"] = sd
        axs[ai].plot(mu, depth[:NL57], color=col, label=nm)
        axs[ai].fill_betweenx(depth[:NL57], mu - sd, mu + sd, color=col, alpha=0.2)
        axs[ai + 1].plot(sd, depth[:NL57], color=col, label=nm)
axs[0].set_xlabel("mean T (degC)")
axs[1].set_xlabel("std of T (degC)")
axs[2].set_xlabel("mean S (PSS-78)")
axs[3].set_xlabel("std of S")
axs[0].invert_yaxis()
axs[0].set_ylabel("depth (m)")
axs[0].legend(fontsize=7)
fig.suptitle("Mean +/- std (shaded) of T and S at the 57 levels 0-1500 m")
f7d = savefig(fig, "s7d_mean_std_profiles",
             "Observation-weighted means and population standard deviations by depth; shading is mean ±1 standard deviation, not a confidence interval")
savetab(pd.DataFrame(lvl_tab), "s7_mean_std_57_levels_central_nordic")
WM = [("Polar Mixed Layer (0-50 m)", 0, 50), ("Halocline (50-200 m)", 50, 200), ("Atlantic layer (200-900 m)", 200, 900),
      ("Intermediate (900-2000 m)", 900, 2000), ("Deep / bottom water (>=2000 m)", 2000, 5501)]
rows = []
for nm, msk in (("Central Arctic >80N", CENTRAL), ("Nordic seas", NORDIC)):
    for wl, a, b in WM:
        lv = (depth >= a) & (depth < b)
        tv = T[msk][:, lv]
        sv = S[msk][:, lv]
        tv, sv = tv[~np.isnan(tv)], sv[~np.isnan(sv)]
        rows.append({"region": nm, "layer (nominal depth range)": wl, "T obs": tv.size,
                     "T mean": tv.mean() if tv.size else np.nan, "T std": tv.std() if tv.size else np.nan,
                     "S mean": sv.mean() if sv.size else np.nan, "S std": sv.std() if sv.size else np.nan})
t_wm = pd.DataFrame(rows)
savetab(t_wm, "s7_water_mass_layers")

add(7, "## 7. Thermodynamic screening and descriptive statistics\n",
    "All statistics use the normalised arrays (NaN = missing). Temperature is treated as in-situ temperature and salinity as practical "
    "salinity (an assumption: the files do not state it). Pressure comes from the per-profile `pres` array. Seawater properties were "
    f"computed with `gsw` across the full dataset (including **{int(hasS.sum()):,} profiles with S**), "
    f"in chunks of 20,000 profiles, without random subsampling. A level or pair is assessable only where the "
    f"required T, S and pressure values produce finite seawater properties. These diagnostics flag unusual "
    f"patterns under the stated assumptions; they do not validate the upstream QC, identify erroneous values "
    f"or justify automatic rejection. No instrument uncertainty or independent truth dataset is supplied.\n",
    "### 7.1 T and S per depth band\n", f7a, mdtable(t_band, dec=2),
    "### 7.2 Common-pressure density inversions\n",
    "For each pair of adjacent valid levels the in-situ densities of both parcels were compared at the mid-pressure "
    "(`gsw.rho(SA, CT, p_mid)`); a pair is flagged if the upper parcel is denser than the lower by more than the stated threshold. "
    "Absolute Salinity and Conservative Temperature are calculated using the recorded position and pressure. "
    "This is an adjacent-level screen, not a complete stability analysis or an error diagnosis. "
    "Physical variability, measurement uncertainty and interpolation can all affect the result. "
    "The existing four thresholds are retained; no uncertainty-calibrated threshold is inferred here.\n",
    mdtable(t_inv, dec=3), mdtable(t_invb, dec=3),
    f"Out of {n_pairs:,} adjacent-level pairs, {t_inv.iloc[1, 2]:.3f}% exceed 0.01 kg/m3 and "
    f"{t_inv.iloc[2, 2]:.3f}% exceed 0.03 kg/m3. {t_inv.iloc[1, 4]:.2f}% of profiles with S have at least one pair above 0.01 kg/m3 "
    f"and {t_inv.iloc[2, 4]:.2f}% above 0.03 kg/m3.\n",
    "**Denominators.** Overall profile percentages use all profiles with salinity, including those without an "
    "assessable adjacent-level pair. The additional assessable-profile column and the class table below use "
    "only profiles with at least one such pair; an unassessable profile is not a clean-QC result. "
    "Depth-band percentages use assessable pairs whose upper depth lies in that band.\n",
    "**By inferred sampling-pattern label** (not a source-specific quality comparison)\n", mdtable(t_invcls, dec=2),
    "### 7.3 Freezing-point and vertical-gradient screens\n", f7b,
    f"* Values within +/-0.01 degC of the in-situ freezing temperature (`gsw.t_freezing(SA, p, 0)`): **{n_within:,}** "
    f"({pct(n_within, n_fv, 2)} of {n_fv:,} assessable T/S/pressure level-values); values more than 0.01 degC below freezing: {n_below:,} "
    f"({pct(n_below, n_fv, 3)}). Exactly T = -2.0 degC: {n_Tm2:,} values (the data minimum is {K['T_min']:.2f} degC).\n",
    mdtable(t_spk, dec=3),
    "The freezing calculation uses the air-free (`saturation_fraction=0`) convention. Near-freezing values can be "
    "physically plausible; below-freezing differences require context and the stated temperature/salinity assumptions. "
    "The gradient tests use absolute differences between adjacent non-missing standard levels divided by their "
    "depth separation. They detect steep gradients, not uniquely isolated spikes: real fronts, stratification, "
    "interpolation and measurement problems can all trigger them. Thresholds are fixed screening choices, "
    "not demonstrated limits on acceptable measurements.\n",
    "Per-level counts: `tables/s7_spikes_per_level.csv`, `tables/s7_freezing_point_per_level.csv`.\n",
    "### 7.4 T-S diagram\n", f7c,
    f"The displayed histogram contains {K['TS_n']:,} level-values within its SA/CT window. "
    f"Values outside that window are not shown, so it is not a histogram of every non-missing value. "
    f"No additional quality filter is applied before binning. All depths are pooled and profiles with more levels "
    f"contribute more entries; the diagram is not weighted to equalise profiles, area or time.\n",
    "### 7.5 Mean and spread of T and S, 0-1500 m, in two basins\n", f7d,
    "Central Arctic = the 'Central Arctic (>80N ...)' region of Section 5; Nordic seas = the Fram Strait/Greenland Sea, Iceland Sea and "
    "Norwegian Sea regions. Means and standard deviations pool non-missing values at each depth without "
    "area, time or platform reweighting. The contributing sample changes with depth. Spread combines "
    "spatial, seasonal and interannual variability, sampling differences and any measurement or processing effects; "
    "it is neither a measurement-error estimate nor a confidence interval for a climatological mean.\n",
    "**Nominal water-mass layers by depth range** (fixed depth boundaries, not property-defined)\n", mdtable(t_wm, dec=2))
log("section 7 done")

# ============================================================================= SECTION 8: fitness for purpose
log("section 8")
LV = {z: int(np.where(depth == z)[0][0]) for z in (10, 100, 300, 500, 700, 2000)}
OCN = OCEAN_BOX
HIGH = (CLAT.ravel() >= 70) & OCN
QMASK = {q: np.isin(month, m) for q, m in SEASONS_WOA.items()}


def covA(pm, lv, thr=10, region=None):
    c = np.bincount(cell[pm & vT[:, lv]], minlength=NC)
    ok = (c >= thr) & OCN
    if region is None:
        return 100 * area_frac(ok)
    return 100 * OCEAN_AREA[ok & region].sum() / OCEAN_AREA[region].sum()


allp = np.ones(N, bool)
q500 = {q: covA(QMASK[q], LV[500]) for q in QMASK}
q10 = {q: covA(QMASK[q], LV[10]) for q in QMASK}
m1 = float(np.mean(list(q500.values())))
m4 = float(np.mean(list(q10.values())))
m2 = covA(allp, LV[300], region=HIGH)
m3 = covA(allp, LV[2000])
c100S = np.bincount(cell[vS[:, LV[100]]], minlength=NC)
m5 = 100 * area_frac((c100S >= 10) & OCN)
m6 = 100 * area_frac((nyr >= 10) & OCN)
cen_m = []
for mth in range(1, 13):
    pm = CENTRAL & (month == mth) & (maxdepth >= 200)
    cen_m.append((MONTH_NAMES[mth - 1], int(pm.sum()), int(len(np.unique(year[pm]))), int(len(np.unique(year[CENTRAL & (month == mth)])))))
t_cen = pd.DataFrame(cen_m, columns=["month", "central-Arctic profiles reaching >=200 m", "distinct years (those profiles)",
                                     "distinct years (any central profile)"])
savetab(t_cen, "s8_central_arctic_monthly")
cen_ok = (t_cen.iloc[:, 1] >= 100) & (t_cen.iloc[:, 2] >= 10)
m7 = 100 * float(cen_ok.sum()) / 12
EUL = (nyr >= 10) & (nmo >= 6) & (cnt >= 100)
m8 = 100 * area_frac(EUL & OCN)
n_fixed = int(is_fixed.sum())
m9 = 100 * K["g1_ge10"]
c700d = [100 * area_frac((np.bincount(cell[(maxdepth >= 700) & (dec_idx == i)], minlength=NC) >= 10) & OCN) for i in range(4)]
m10 = 100 * area_frac((np.bincount(cell[maxdepth >= 700], minlength=NC) >= 10) & OCN)


def rate(v):
    return "Good" if v >= 60 else ("Marginal" if v >= 20 else "Poor")


FIT = [
    ("Gridded climatology 0-1500 m (monthly/seasonal)", m1,
     "mean over the 4 WOA quarters of the % of 60-90N ocean area with >=10 T profiles at 500 m; by quarter: "
     + ", ".join(f"{q} {v:.0f}%" for q, v in q500.items())),
    ("Regional climatology of the Atlantic layer (north of 70N)", m2,
     f"% of mask ocean area in boxes with centres >=70N and >=10 T profiles at 300 m (all years); "
     f"{int(t_cen.iloc[:, 1].sum()):,} central-Arctic profiles reach >=200 m"),
    ("Basin-scale deep (>=2000 m) analysis", m3,
     f"% of 60-90N mask ocean area with >=10 T profiles at 2000 m (shelves are included; not a bathymetry-limited deep-water fraction); "
     f"{K['prof_deep']:,} profiles ({pct(K['prof_deep'], N, 1)}) have data below 1500 m"),
    ("Mixed layer / near-surface (0-50 m)", m4,
     "mean over quarters of % ocean area with >=10 profiles at 10 m (0 m is empty, 5 m partly filled); by quarter: "
     + ", ".join(f"{q} {v:.0f}%" for q, v in q10.items())),
    ("Halocline structure (50-200 m; needs S)", m5, "% ocean area with >=10 S profiles at 100 m (all years)"),
    ("Interannual variability / trends", m6,
     "% ocean area sampled in >=10 distinct years (not accounting for profiles per year or season)"),
    ("Seasonal cycle in the central Arctic", m7,
     f"% of calendar months with >=100 central-Arctic profiles reaching >=200 m from >=10 distinct years; {int(cen_ok.sum())} of 12 months pass"),
    ("Eulerian time series at a point", m8,
     f"% mask ocean area in 100 km boxes with >=100 profiles, >=10 years, >=6 calendar months; candidate low-extent series: {n_fixed}"),
    ("Validation of models/reanalyses (point-wise, 1 deg)", m9,
     "% of WOA-mask 1 deg ocean cells north of 60N with >=10 profiles (all years, any depth; cell-count rather than area fraction)"),
    ("Ocean heat content 0-700 m", m10,
     "% ocean area with >=10 profiles reaching >=700 m (all years); by decade: " + ", ".join(f"{l} {v:.0f}%" for l, v in zip(DEC_LAB, c700d))),
]
winter_years = t_cen.loc[t_cen["month"].isin(["Dec", "Jan", "Feb", "Mar"]), "distinct years (those profiles)"]
NOTES = {
    FIT[0][0]: "A 500 m quarterly proxy does not establish monthly support through 1500 m; interpolation choices require separate assessment.",
    FIT[1][0]: "Counts at 300 m are a coverage proxy, not property-based identification of Atlantic-layer water or regional representativeness.",
    FIT[2][0]: "Few values reach these depths; the metric does not identify their source platforms or their seasonal representativeness.",
    FIT[3][0]: f"No 0 m values; 5 m fill is {pct(fr5, 1)}. Values at 10 m do not resolve mixed-layer structure or surface conditions.",
    FIT[4][0]: "A 100 m salinity count does not establish resolved halocline gradients, measurement quality or independent sampling.",
    FIT[5][0]: "Distinct-year counts do not ensure contiguous records, consistent locations/seasons or a valid trend uncertainty model.",
    FIT[6][0]: f"Only {pct((CENTRAL & (year < 2003)).sum(), CENTRAL.sum(), 1)} of central-Arctic profiles predate 2003; Dec-Mar depth-reaching profiles cover {winter_years.min()}-{winter_years.max()} distinct years per month. Unequal year/location support can confound a seasonal comparison.",
    FIT[7][0]: "A sampled 100 km box is not a fixed-point time series; candidate low-extent chains also require original-record verification.",
    FIT[8][0]: "Evaluate at observed locations, times and depths; account for representativeness and possible overlap with model assimilation data.",
    FIT[9][0]: "Depth reach alone does not ensure complete 0-700 m profiles, valid heat integration or representative basin/time coverage.",
}
t_fit = pd.DataFrame([{"purpose": a, "coverage metric (%)": b, "coverage category": rate(b), "metric definition / evidence": c,
                       "comment": NOTES[a]} for a, b, c in FIT])
savetab(t_fit, "s8_fitness_for_purpose")
add(8, "## 8. Application-specific coverage screens\n",
    "The table retains the original heuristic categories: **Good** >=60%, **Marginal** >=20% and <60%, "
    "**Poor** <20%. They describe a single sampling-count proxy for each application, **not validated fitness "
    "for that purpose**. In particular, a 'Good' category is not an endorsement for assimilation, trend estimation, "
    "independent model validation or a resolved seasonal climatology.\n",
    f"Most metrics weight the {K['ocean_area_mkm2']:.1f} million km2 of EASE2 mask ocean area north of 60N by "
    "the fraction in 100 km boxes meeting a count rule. The Atlantic-layer proxy instead uses the area in boxes "
    "whose centres are >=70N. The central-Arctic seasonal proxy is a fraction of twelve calendar months, "
    "and the model-validation proxy is a fraction of WOA-mask 1 deg ocean cells, not area. "
    "Each row states its denominator and selection. Shelves and marginal seas are retained in the broad area denominator; "
    "no bathymetric deep-ocean denominator is calculated. The categories are not directly comparable across these different metrics.\n",
    "No proxy measures independent sample size, instrument uncertainty, full depth coverage, temporal homogeneity "
    "or causal effects of observing systems. Counts pool years unless stated otherwise, and a depth-range application "
    "is generally represented by only one selected depth or a depth-reach criterion.\n",
    mdtable(t_fit, dec=1),
    "**Central Arctic by calendar month** (profiles reaching >=200 m)\n", mdtable(t_cen),
    "### Implications for analysis\n",
    f"* **Observed strengths.** {n_both:,} profiles contain both T and S; the upper ocean is much more densely sampled than "
    "deep levels. Regional, depth-specific and time-specific counts provide a starting point for selecting supported analyses.\n"
    f"* **Observed constraints.** {pct(zero_area, 1)} of mask ocean area lies in boxes without any profile, "
    f"and {pct(K['cov10_all'], 1)} lies in boxes with >=10 over the whole record. Those counts do not quantify "
    "the effective number of independent observations. Spatial and calendar support change through time; deep and surface-level data are limited.\n"
    "* **Use with explicit safeguards.** For regional or profile-location comparisons, verify the relevant year, season, "
    "depth and QC support rather than relying on total counts. Gridding and heat-content integration require a documented "
    "missing-data and weighting method. Trends require evaluation of changing sampling and serial dependence. "
    "Model validation requires checking whether these observations contributed to the evaluated analysis or assimilation. "
    "This report does not provide an observation-error model, establish validation independence or certify reliability "
    "for assimilation.\n")
log("section 8 done")

# ============================================================================= assemble the markdown
log("assembling markdown")
yr_cnt = np.bincount(year, minlength=year.max() + 1)
post05 = int((year >= 2005).sum())
summer = int(np.isin(month, [6, 7, 8, 9]).sum())
winter = int(np.isin(month, [12, 1, 2, 3]).sum())
y_peak = int(np.argmax(yr_cnt))
cls_txt = "; ".join(f"{c} {K['class_counts'][c]:,} ({pct(K['frac_class'][c], 1, 1)})" for c in CLASSES)
EXEC = [
    f"* **Scope.** {N:,} profile records cover {K['year_min']}-{K['year_max']}; {int(hasT.sum()):,} contain temperature "
    f"and {n_both:,} contain both temperature and salinity. Profiles with temperature but no salinity: {n_noS:,}. "
    f"The median profile has {K['med_nlev']:.0f} non-missing temperature levels and reaches {K['med_maxdepth']:.0f} m.",
    f"* **Depth.** Fraction reaching >=700 m: {pct((maxdepth >= 700).sum(), N)}; >=1500 m: {pct((maxdepth >= 1500).sum(), N)}; "
    f">=3000 m: {pct((maxdepth >= 3000).sum(), N)}. The 0 m level is empty; 5 m is {pct(nT_lev[1], N, 0)} filled, 10 m {pct(nT_lev[2], N, 0)}.",
    f"* **Time.** {pct(post05, N)} of the profiles are from 2005 onwards; peak year {y_peak} ({yr_cnt[y_peak]:,}). "
    f"Jun-Sep holds {pct(summer, N)} of profiles, Dec-Mar {pct(winter, N)}.",
    f"* **Spatial support.** Of {K['ocean_area_mkm2']:.1f} million km2 of mask ocean area north of 60N, "
    f"{pct(zero_area, 1)} lies in 100 km boxes without any profile; {pct(K['cov10_all'], 1)} lies in boxes with >=10, "
    f"and {pct(K['nyr10_frac'], 1)} in boxes sampled in >=10 distinct years. These are box-based coverage proxies, not resolved ocean coverage.",
    f"* **Unvalidated sampling-pattern labels.** {cls_txt}. No platform identifiers are available; "
    f"ship-CTD-like includes residual patterns that can arise from other observing systems.",
    f"* **Proximity and temperature equality.** {K['dup6_pairs']:,} unordered pairs meet the 0.1 deg / 6 h screen "
    f"and {K['dup1km_pairs']:,} meet the 1 km / 1 h screen. For temperature arrays with >=5 valid levels, "
    f"{K['ident_redundant5']:,} records lie beyond one per identical-T group. None of these counts establishes erroneous duplication.",
    f"* **Thermodynamic screening.** {t_inv.iloc[1, 4]:.1f}% of profiles with salinity have a common-pressure "
    f"density inversion >0.01 kg/m3; {K['spike_T_prof']:,} profiles have an adjacent-level temperature gradient "
    f">0.3 degC/m; {pct(n_within, n_fv, 2)} of assessable level-values are within 0.01 degC of freezing. "
    f"These are screening flags, not validated error diagnoses.",
    "* **Analysis implications.** Section 8 reports application-specific count proxies using the existing "
    "Good/Marginal/Poor categories. Even a Good category does not establish independence, trend reliability, "
    "assimilation suitability or a resolved seasonal climatology.",
]
HEAD = [
    "# AOTD profile data: coverage, sampling patterns and screening diagnostics\n",
    f"*Generated {GEN_DATE} by [the report script](../../scripts/07_report_scatter_data.py) from "
    f"{'`data/scatter_all.npz` (normalised cache)' if NPZ.exists() else 'the yearly files loaded through `scatter_io.load_scatter`'} "
    "and the local ocean-mask inputs listed below.*\n",
    "This report describes the released profile records, their spatial, temporal and depth support, "
    "and transparent screening diagnostics. It does not reconstruct platform identities, "
    "validate upstream quality control or demonstrate fitness for a scientific application.\n",
    "## Executive summary\n", "\n".join(EXEC) + "\n",
    "## Data notes and conventions\n",
    "* **Missing versus valid.** Raw-file salinity and pressure sentinels and density entries at missing temperature "
    "levels are handled by `scatter_io.normalize_missing_values`. This report uses the normalised arrays with NaN "
    "as the missing-value marker; *valid* means present, not independently QC-verified.\n"
    f"* **T/S availability.** {K['n_valid_TS']:,} level-values contain both temperature and salinity; "
    f"entries with only one variable: {K['TS_mismatch']:,}. The 0 m level has neither variable in any profile. "
    "Depth reach is based on the deepest non-missing temperature level unless otherwise specified.\n"
    "* **Thermodynamic convention.** T is treated as in-situ temperature and S as practical salinity; "
    "this convention is an assumption because the files do not fully document it. Pressure is the recorded per-profile array. "
    "Derived TEOS-10 quantities are diagnostics under these assumptions, not a substitute for source metadata.\n"
    "* **Counts and denominators.** A profile is a published record; a level-value is one variable at one depth. "
    "Percentages use the denominator stated in each table and may round a nonzero count to 0.0%. "
    "Profile counts are not independent sample sizes. Standard deviations are population descriptors (ddof=0).\n"
    "* **Spatial conventions.** EASE2 boxes are 100 km x100 km blocks in EPSG:6931. Mask ocean area is "
    "the number of 25 km water cells with centres >=60N times 625 km2. Longitudes in tables are degrees east, "
    "so negative values are west. Regions use ordered approximate lon/lat rules (Section 5).\n"
    "* **Evidence levels.** Source inventories are externally reported provenance; counts are observed diagnostics; "
    "ITP-like, Argo-like and ship-CTD-like are unvalidated sampling-pattern inferences. "
    "Proximity, temperature equality and thermodynamic flags do not establish erroneous observations.\n",
    "## Report guide\n",
    "[Overview](#1-overview-of-the-profile-data) · [Time](#2-temporal-distribution) · "
    "[Space](#3-spatial-distribution) · [Depth](#4-depth-distribution) · "
    "[Regions](#5-regional-distribution) · "
    "[Heuristics and proximity](#6-heuristic-sampling-pattern-labels-and-proximity-screening) · "
    "[Thermodynamic diagnostics](#7-thermodynamic-screening-and-descriptive-statistics) · "
    "[Application screens](#8-application-specific-coverage-screens) · "
    "[Reproducibility](#reproducibility-and-source-attribution)\n",
]
parts = HEAD
for i in range(1, 9):
    parts += SEC[i]
log("recording input checksums")
input_rows = []
for path in INPUT_FILES + [EDIR / "EASE2_N25km.LOCImask_land50_coast0km.720x720.bin",
                           EDIR / "NSIDC0772_LatLon_EASE2_N25km_v1.1.nc", WOA_MASK_FILE]:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    input_rows.append({"input": str(path.relative_to(ROOT)), "bytes": path.stat().st_size,
                       "sha256": digest.hexdigest()})
input_manifest = pd.DataFrame(input_rows)
savetab(input_manifest, "s0_input_manifest")
versions = "; ".join(f"{name} {importlib.metadata.version(name)}" for name in
                     ("numpy", "scipy", "pandas", "xarray", "netCDF4", "pyproj", "matplotlib", "cartopy", "gsw"))
parts += [
    "## Reproducibility and source attribution\n",
    "Run from the repository root using the environment described in [README](../../README.md):\n"
    "```bash\n.venv/bin/python scripts/07_report_scatter_data.py\n```\n"
    "The script writes this report, CSV tables and PNG figures only under `reports/scatter/`. "
    "It uses an existing normalised cache when present; otherwise it loads and normalises the local yearly "
    "scatter files. It does not download inputs, modify observations, remove flagged records or rerun the main pipeline. "
    "Required ancillary files are the NSIDC 25 km land mask and coordinate file in `data/ease2/`, plus "
    "`data/woa23_raw/woa23_decav91C0_t00_04.nc` for the latitude/longitude ocean mask. "
    "The WOA23 mask comes from non-missing annual surface temperature, not from a bathymetric depth mask. "
    "The EASE2 Sea Ice Regions file is not used by this script.\n",
    f"**Runtime.** Python {sys.version.split()[0]}; {versions}. "
    "There is no random subsampling. The generation date is the local run date. Greedy track linking uses stable "
    "time sorting; tied-cost assignments can depend on input order and library behaviour, so retain the input "
    "manifest and runtime versions alongside any scientific use of these diagnostics.\n",
    "**Input identity for this run.** The SHA-256 manifest records local bytes, not independent verification "
    "of the original download or authors' upstream workflow. A pre-existing cache is not rebuilt by this report.\n",
    mdtable(input_manifest),
    "**Sources and credit.**\n"
    "* Li, J., Wu, X. & Wang, X. (2025). *An Arctic Ocean Thermohaline Dataset*. "
    "*Scientific Data* 12, 1607. [Article and compilation methods](https://doi.org/10.1038/s41597-025-05855-3).\n"
    "* Li, J. (2025). *AOTD_scatter*, figshare, version 1, CC BY 4.0. "
    "[Published scatter data](https://doi.org/10.6084/m9.figshare.28825490.v1). "
    "These observations are credited to the dataset authors and contributing observing programmes; "
    "the heuristic labels in this report are not metadata supplied by those programmes.\n"
    "* NOAA NCEI, [World Ocean Atlas 2023](https://www.ncei.noaa.gov/products/world-ocean-atlas), "
    "1991-2020, 0.25 deg annual temperature file: used only to define the latitude/longitude ocean mask here.\n"
    "* NSIDC, [NSIDC-0772](https://nsidc.org/data/nsidc-0772) coordinate ancillary data and "
    "[NSIDC-0609](https://nsidc.org/data/nsidc-0609) EASE-Grid 2.0 land/ocean/coast/ice mask: "
    "used for the equal-area grid and ocean-area accounting.\n"
    "* WHOI, [Ice-Tethered Profiler programme](https://www.whoi.edu/itp/), and "
    "[Argo programme](https://argo.ucsd.edu/): historical context for the benchmark years, "
    "not profile-identity validation.\n",
    "For the distinct expected-spread products and their assumptions, see the "
    "[technical note](../../technical_note.md). None of the screening thresholds, heuristic labels or coverage "
    "categories here should be interpreted as calibrated observation errors for those products.\n",
    "## Appendix: complete CSV table index\n",
    "\n".join(f"* [{n}](tables/{n})" for n in sorted(set(TABLES_WRITTEN))) + "\n",
]
(OUT / "scatter_report.md").write_text("\n".join(parts), encoding="utf-8")
log(f"done: wrote {OUT / 'scatter_report.md'}")
