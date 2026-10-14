"""Expected standard deviation and expected RMSE of in situ potential temperature (theta) and salinity (S), on the EASE2 100 km grid.

For every in situ profile value (scatter data, 1992 onwards) the deviation d = obs - reference is formed at the profile's location
(nearest WOA23 0.25 deg cell). Both the observation and the reference are converted to potential temperature (UNESCO 1983, surface
reference, scripts/potential_temperature.py) BEFORE differencing: d_theta = theta(T_obs, S_obs) - theta(T_ref, S_ref). WOA23 t_an, the scatter
T and the published AOTD T are all treated as in situ. Where S_obs is missing but T_obs is valid, S_ref is used transiently for theta(T_obs).
References (one product each):
    woa23           WOA23 decav91C0 monthly (0-1500 m) or seasonal (deeper) field of the profile month/season
    aotd            published AOTD monthly climatology (0-1500 m, treated as the WOA23 0.25 deg grid, rows 599-719); NaN deeper
    woa23_and_aotd  elementwise fmax of the two metric fields, ignoring NaN (WOA23 alone below 1500 m)
Per 100 km cell and standard depth level (pooled over 1992-2023, no time dimension, NaN where no observation):
    rmse = sqrt(mean(d^2))                 spread about the reference, bias included
    std  = sqrt(mean(d^2) - mean(d)^2)     spread about the cell's mean offset (population std; 0 when n=1)
    bias = mean(d), n = number of values

Usage (from the repo root): python scripts/03_compute_expected_spread.py   (about 1.5 min, about 4 GB RAM)
Outputs (untracked): outputs/expected_std/expected_{std,rmse}_from_{woa23,aotd,woa23_and_aotd}_{potential_temperature,salinity}_ease2_100km.nc
"""
import os
import numpy as np
import sys
import xarray as xr
from pyproj import Transformer

from potential_temperature import calculate_potential_T

RAW = "data/woa23_raw"
OUT = "outputs/expected_std"
os.makedirs(OUT, exist_ok=True)
NMONTHLY = 57  # levels with a monthly WOA23 field (0-1500 m); deeper levels use seasonal
R0 = 599  # WOA23 lat row of the cell centred on 59.875 N
NROW = 720 - R0
CELL, HALF = 100e3, 9000e3
NC = 180
YEAR0 = 1992  # first year of in situ data used
SEASON_FILE = np.array([13, 13, 13, 14, 14, 14, 15, 15, 15, 16, 16, 16])  # WOA23 winter, spring, summer, autumn


def load_ref(letter, var):
    """monthly (12, 57, NROW, 1440) and seasonal (4, 102-57, NROW, 1440) float32 reference fields, plus depth (102,)."""
    mon = np.empty((12, NMONTHLY, NROW, 1440), "float32")
    for m in range(12):
        ds = xr.open_dataset(f"{RAW}/woa23_decav91C0_{letter}{m + 1:02d}_04.nc", decode_times=False)
        mon[m] = ds[var].isel(time=0, lat=slice(R0, 720)).values
        ds.close()
    ses = None
    for i, f in enumerate((13, 14, 15, 16)):
        ds = xr.open_dataset(f"{RAW}/woa23_decav91C0_{letter}{f:02d}_04.nc", decode_times=False)
        a = ds[var].isel(time=0, lat=slice(R0, 720), depth=slice(NMONTHLY, None)).values
        depth = ds.depth.values
        ds.close()
        if ses is None:
            ses = np.empty((4,) + a.shape, "float32")
        ses[i] = a
    return mon, ses, depth


def load_aotd(key):
    """published AOTD monthly field (12, 57, NROW, 1440) float32."""
    out = np.empty((12, NMONTHLY, NROW, 1440), "float32")
    for m in range(12):
        ds = xr.open_dataset(f"data/aotd_climatology/climatology_{m + 1:02d}.nc")
        out[m] = ds[key].values
        ds.close()
    return out


def main():
    d = np.load("data/scatter_all.npz")
    lat, lon, month = d["lat"].astype("float64"), d["lon"].astype("float64"), d["month"].astype(int)
    depth_all = d["depth"].astype("float64")
    nlev = depth_all.shape[0]
    x, y = Transformer.from_crs(4326, 6931, always_xy=True).transform(lon, lat)
    ix = np.floor((x + HALF) / CELL).astype(int)
    iy = np.floor((HALF - y) / CELL).astype(int)
    year = d["year"].astype(int)
    inside = (year >= YEAR0) & (ix >= 0) & (ix < NC) & (iy >= 0) & (iy < NC)
    cell = iy * NC + ix
    row = np.floor((lat + 90.0) / 0.25).astype(int) - R0
    col = np.floor((lon + 180.0) / 0.25).astype(int) % 1440
    inside &= (row >= 0) & (row < NROW)
    print("profiles", len(lat), "inside grid/domain", int(inside.sum()), flush=True)

    monT, sesT, depth = load_ref("t", "t_an")
    monS, sesS, _ = load_ref("s", "s_an")
    assert nlev == len(depth)
    aT, aS = load_aotd("temp"), load_aotd("salt")
    To, So = d["temp"], d["salt"]
    z = depth_all[:, None]
    acc = {(rn, v): [np.zeros((nlev, NC * NC)) for _ in range(3)] for rn in ("woa23", "aotd") for v in ("potential_temperature", "salinity")}
    for m in range(1, 13):
        sel = np.nonzero((month == m) & inside)[0]
        r, c, ce, la = row[sel], col[sel], cell[sel], lat[sel]
        t_o, s_o = To[sel].T.astype("float64"), So[sel].T.astype("float64")
        for rn in ("woa23", "aotd"):
            Tr = np.full((nlev, len(sel)), np.nan); Sr = Tr.copy()
            if rn == "woa23":
                Tr[:NMONTHLY], Sr[:NMONTHLY] = monT[m - 1][:, r, c], monS[m - 1][:, r, c]
                k = SEASON_FILE[m - 1] - 13
                Tr[NMONTHLY:], Sr[NMONTHLY:] = sesT[k][:, r, c], sesS[k][:, r, c]
            else:
                Tr[:NMONTHLY], Sr[:NMONTHLY] = aT[m - 1][:, r, c], aS[m - 1][:, r, c]
            S_for_obs = np.where(np.isnan(s_o), Sr, s_o)
            devs = {"potential_temperature": calculate_potential_T(t_o, S_for_obs, z, la) - calculate_potential_T(Tr, Sr, z, la),
                    "salinity": s_o - Sr}
            for v, dev in devs.items():
                n, s1, s2 = acc[(rn, v)]
                good = np.isfinite(dev)
                for k in range(nlev):
                    g = good[k]
                    if g.any():
                        n[k] += np.bincount(ce[g], minlength=NC * NC)
                        s1[k] += np.bincount(ce[g], weights=dev[k, g], minlength=NC * NC)
                        s2[k] += np.bincount(ce[g], weights=dev[k, g] ** 2, minlength=NC * NC)
        print("month", m, "done", flush=True)

    shp = (nlev, NC, NC)
    metrics = {}  # (product, var) -> dict(std, rmse, bias, n)
    for (rn, v), (n, s1, s2) in acc.items():
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = s1 / n
            metrics[(rn, v)] = dict(rmse=np.where(n > 0, np.sqrt(s2 / n), np.nan),
                                    std=np.where(n > 0, np.sqrt(np.maximum(s2 / n - mean ** 2, 0)), np.nan),
                                    bias=np.where(n > 0, mean, np.nan), n=n)
    for v in ("potential_temperature", "salinity"):
        w, a = metrics[("woa23", v)], metrics[("aotd", v)]
        out = {}
        for k in ("std", "rmse"):
            out[k] = np.fmax(w[k], a[k])
            pick_w = np.isnan(a[k]) | (w[k] >= a[k])  # bias and n come from the reference giving the larger value
            out[f"bias_{k}"] = np.where(pick_w, w["bias"], a["bias"])
            out[f"n_{k}"] = np.where(pick_w, w["n"], a["n"])
        metrics[("woa23_and_aotd", v)] = out

    xc = -HALF + CELL * (np.arange(NC) + 0.5)
    yc = HALF - CELL * (np.arange(NC) + 0.5)
    X, Y = np.meshgrid(xc, yc)
    glon, glat = Transformer.from_crs(6931, 4326, always_xy=True).transform(X, Y)
    titles = {"woa23": "WOA23", "aotd": "AOTD", "woa23_and_aotd": "the larger of the WOA23-based and AOTD-based estimates"}
    for (rn, v), mt in metrics.items():
        short = "Theta" if v == "potential_temperature" else "S"
        unit = "degree_Celsius" if v == "potential_temperature" else "1e-3"
        pname = v.replace("_", " ")
        for kind, desc in (("std", "standard deviation about the cell mean offset"), ("rmse", "RMSE (about the reference, bias included)")):
            bias = mt[f"bias_{kind}"] if rn == "woa23_and_aotd" else mt["bias"]
            n = mt[f"n_{kind}"] if rn == "woa23_and_aotd" else mt["n"]
            ds = xr.Dataset(
                {
                    f"{short}_{kind}": (("depth", "y", "x"), mt[kind].reshape(shp).astype("float32"),
                                        dict(long_name=f"expected {desc} of in situ {pname} from {titles[rn]}", units=unit,
                                             grid_mapping="crs", coordinates="lat lon",
                                             comment="Deviation = obs - reference, both as potential temperature (T only), pooled per cell 1992-2023, no time dimension. "
                                                     f"Product: {rn}. NaN where no observation; AOTD reference exists only to 1500 m.")),
                    f"{short}_bias": (("depth", "y", "x"), bias.reshape(shp).astype("float32"),
                                      dict(long_name=f"mean in situ {pname} minus reference", units=unit, grid_mapping="crs", coordinates="lat lon")),
                    f"{short}_n": (("depth", "y", "x"), n.reshape(shp).astype("int32"),
                                   dict(long_name="number of in situ values used", units="1", grid_mapping="crs", coordinates="lat lon")),
                    "crs": ((), np.int32(0), dict(grid_mapping_name="lambert_azimuthal_equal_area", latitude_of_projection_origin=90.0,
                                                  longitude_of_projection_origin=0.0, false_easting=0.0, false_northing=0.0,
                                                  semi_major_axis=6378137.0, inverse_flattening=298.257223563,
                                                  long_name="EASE-Grid 2.0 North 100 km (EPSG:6931)")),
                    "lat": (("y", "x"), glat.astype("float32"), dict(standard_name="latitude", units="degrees_north")),
                    "lon": (("y", "x"), glon.astype("float32"), dict(standard_name="longitude", units="degrees_east")),
                },
                coords=dict(depth=("depth", depth.astype("float32"), dict(standard_name="depth", units="m", positive="down", axis="Z")),
                            y=("y", yc, dict(standard_name="projection_y_coordinate", units="m", axis="Y")),
                            x=("x", xc, dict(standard_name="projection_x_coordinate", units="m", axis="X"))),
                attrs=dict(title=f"Expected {kind} of in situ {pname} from {titles[rn]} on the EASE2 100 km Arctic grid",
                           Conventions="CF-1.8", history="scripts/03_compute_expected_spread.py"))
            f = f"{OUT}/expected_{kind}_from_{rn}_{v}_ease2_100km.nc"
            ds.to_netcdf(f, encoding={k: dict(zlib=True, complevel=4) for k in (f"{short}_{kind}", f"{short}_bias", f"{short}_n")})
            print("wrote", f, flush=True)


if __name__ == "__main__":
    main()
