"""Expected sample standard deviation and expected RMSE of in situ potential temperature (theta) and salinity (S), on the EASE2 100 km grid.

For every in situ profile value (scatter data, 1992 onwards) the deviation d = obs - reference is formed at the profile's location
and timestamp using bilinear spatial and linear temporal interpolation of the 0.25 deg reference climatology (with coast/mask-aware
weight renormalization over finite corner/month values). Both the observation and the reference are converted to potential temperature
(UNESCO 1983, surface reference, scripts/potential_temperature.py) BEFORE differencing: d_theta = theta(T_obs, S_obs) - theta(T_ref, S_ref).
WOA23 t_an, the scatter T and the published AOTD T are all treated as in situ. Where S_obs is missing but T_obs is valid, S_ref is used
transiently for theta(T_obs).
References (one product each):
    woa23           WOA23 decav91C0 monthly (0-1500 m) or seasonal (deeper) field interpolated to profile location and time
    aotd            published AOTD monthly climatology (0-1500 m, treated as the WOA23 0.25 deg grid, rows 599-719); NaN deeper
    woa23_and_aotd  elementwise fmax of the two metric fields, ignoring NaN (WOA23 alone below 1500 m)
Per 100 km cell and standard depth level (pooled over 1992-2023, no time dimension, NaN where no observation):
    rmse    = sqrt(mean(d^2))                                spread about the reference, bias included (NaN when n < 1)
    std     = sqrt(n / (n - 1) * (mean(d^2) - mean(d)^2))    sample standard deviation about the cell's mean offset (ddof=1; NaN when n < 2)
    bias    = mean(d), n = number of finite comparisons
    n_eff   = number of distinct (100 km cell, UTC calendar day) sampling days with valid comparisons (de-clustered count proxy)
    n_years = number of distinct calendar years (1992-2023) with valid comparisons

Usage (from the repo root): python scripts/03_compute_expected_spread.py   (about 1.5 min, about 4 GB RAM)
Outputs (untracked): outputs/expected_spread/expected_{std,rmse}_from_{woa23,aotd,woa23_and_aotd}_{potential_temperature,salinity}_ease2_100km.nc
"""
import os
import numpy as np
import xarray as xr
from pyproj import Transformer

from potential_temperature import calculate_potential_T

RAW = "data/woa23_raw"
OUT = "outputs/expected_spread"
NMONTHLY = 57  # levels with a monthly WOA23 field (0-1500 m); deeper levels use seasonal
R0 = 599  # WOA23 lat row of the cell centred on 59.875 N
NROW = 720 - R0
CELL, HALF = 100e3, 9000e3
NC = 180
YEAR0 = 1992  # first year of in situ data used
SEASON_FILE = np.array([13, 13, 13, 14, 14, 14, 15, 15, 15, 16, 16, 16])  # WOA23 winter, spring, summer, autumn


def check_depth_levels(actual, expected, source):
    if not np.array_equal(actual, expected):
        raise ValueError(f"Depth levels in {source} do not match the standard depth grid")


def profile_grid_indices(lat, lon):
    lat, lon = np.asarray(lat), np.asarray(lon)
    if not (np.isfinite(lat).all() and np.isfinite(lon).all()):
        raise ValueError("Profile latitude and longitude must be finite")
    if np.any((lat < -90) | (lat > 90) | (lon < -180) | (lon > 180)):
        raise ValueError("Profile latitude or longitude is outside geographic bounds")
    x, y = Transformer.from_crs(4326, 6931, always_xy=True).transform(lon, lat)
    ix = np.floor((x + HALF) / CELL).astype(int)
    iy = np.floor((HALF - y) / CELL).astype(int)
    # The northern edge belongs to the last WOA23 latitude cell.
    row = np.minimum(np.floor((lat + 90.0) / 0.25).astype(int), 719) - R0
    col = np.floor((lon + 180.0) / 0.25).astype(int) % 1440
    inside = (lat >= 60) & (ix >= 0) & (ix < NC) & (iy >= 0) & (iy < NC) & (row >= 0) & (row < NROW)
    return iy * NC + ix, row, col, inside


def spatial_interp_weights(lat, lon):
    """Bilinear interpolation indices and weights on the 0.25 deg reference subgrid (rows R0..719, cols 0..1439)."""
    lat, lon = np.asarray(lat, dtype="float64"), np.asarray(lon, dtype="float64")
    rf = (lat - (-89.875 + R0 * 0.25)) / 0.25
    r0 = np.clip(np.floor(rf).astype(int), 0, NROW - 1)
    r1 = np.minimum(r0 + 1, NROW - 1)
    wr = np.clip(rf - r0, 0.0, 1.0)

    cf = (lon + 179.875) / 0.25
    c0_raw = np.floor(cf).astype(int)
    wc = cf - c0_raw
    c0 = c0_raw % 1440
    c1 = (c0 + 1) % 1440
    return r0, r1, wr, c0, c1, wc


def temporal_interp_weights(time_sec):
    """Linear temporal interpolation indices and weights for 12 monthly and 4 seasonal reference fields."""
    ts = np.asarray(time_sec, dtype="int64").astype("datetime64[s]")
    y0 = ts.astype("datetime64[Y]").astype("datetime64[s]")
    y1 = (ts.astype("datetime64[Y]") + np.timedelta64(1, "Y")).astype("datetime64[s]")
    frac = (ts - y0) / (y1 - y0)

    pm = frac * 12.0 - 0.5
    m0 = np.floor(pm).astype(int) % 12
    m1 = (m0 + 1) % 12
    wm = pm - np.floor(pm)

    ps = frac * 4.0 - 0.5
    s0 = np.floor(ps).astype(int) % 4
    s1 = (s0 + 1) % 4
    ws = ps - np.floor(ps)
    return m0, m1, wm, s0, s1, ws


def interp_ref_2d_time(field, t0, t1, wt, r0, r1, wr, c0, c1, wc):
    """Bilinear spatial + linear temporal interpolation with mask-aware weight renormalization over finite corners."""
    w00 = (1.0 - wr) * (1.0 - wc)
    w01 = (1.0 - wr) * wc
    w10 = wr * (1.0 - wc)
    w11 = wr * wc
    num = np.zeros((field.shape[1], len(r0)), dtype="float64")
    den = np.zeros_like(num)
    for ti, tw in ((t0, 1.0 - wt), (t1, wt)):
        for ri, ci, sw in ((r0, c0, w00), (r0, c1, w01), (r1, c0, w10), (r1, c1, w11)):
            v = field[ti, :, ri, ci].T.astype("float64") if field.ndim == 4 else field[ti][:, ri, ci]
            w = (tw * sw)[None, :]
            fin = np.isfinite(v)
            num += np.where(fin, v * w, 0.0)
            den += np.where(fin, w, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)


def load_ref(letter, var):
    """monthly (12, 57, NROW, 1440) and seasonal (4, 102-57, NROW, 1440) float32 reference fields, plus depth (102,)."""
    mon = np.empty((12, NMONTHLY, NROW, 1440), "float32")
    monthly_depth = None
    for m in range(12):
        path = f"{RAW}/woa23_decav91C0_{letter}{m + 1:02d}_04.nc"
        with xr.open_dataset(path, decode_times=False) as ds:
            if monthly_depth is None:
                monthly_depth = ds.depth.values
            check_depth_levels(ds.depth.values, monthly_depth, path)
            mon[m] = ds[var].isel(time=0, lat=slice(R0, 720)).values
    ses = None
    depth = None
    for i, f in enumerate((13, 14, 15, 16)):
        path = f"{RAW}/woa23_decav91C0_{letter}{f:02d}_04.nc"
        with xr.open_dataset(path, decode_times=False) as ds:
            if depth is None:
                depth = ds.depth.values
            check_depth_levels(ds.depth.values, depth, path)
            a = ds[var].isel(time=0, lat=slice(R0, 720), depth=slice(NMONTHLY, None)).values
        if ses is None:
            ses = np.empty((4,) + a.shape, "float32")
        ses[i] = a
    check_depth_levels(monthly_depth, depth[:NMONTHLY], f"{letter} monthly WOA23 files")
    return mon, ses, depth


def load_aotd(key, depth):
    """published AOTD monthly field (12, 57, NROW, 1440) float32."""
    out = np.empty((12, NMONTHLY, NROW, 1440), "float32")
    for m in range(12):
        path = f"data/aotd_climatology/climatology_{m + 1:02d}.nc"
        with xr.open_dataset(path) as ds:
            check_depth_levels(ds.depth.values, depth[:NMONTHLY], path)
            out[m] = ds[key].values
    return out


def deviation_statistics(n, s1, s2, n_eff=None, n_years=None):
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = s1 / n
        pop_var = np.maximum(s2 / n - mean ** 2, 0.0)
        sample_var = np.where(n >= 2, pop_var * (n / (n - 1.0)), np.nan)
        out = dict(
            rmse=np.where(n > 0, np.sqrt(s2 / n), np.nan),
            std=np.where(n >= 2, np.sqrt(sample_var), np.nan),
            bias=np.where(n > 0, mean, np.nan),
            n=n,
        )
        if n_eff is not None:
            out["n_eff"] = n_eff
        if n_years is not None:
            out["n_years"] = n_years
        return out


def combine_statistics(w, a):
    out = {}
    for kind in ("std", "rmse"):
        out[kind] = np.fmax(w[kind], a[kind])
        # Resolve ties at the precision stored in NetCDF, so metadata is reproducible from the products.
        w_stored, a_stored = w[kind].astype("float32"), a[kind].astype("float32")
        pick_w = np.isnan(a_stored) | (w_stored >= a_stored)
        out[f"bias_{kind}"] = np.where(pick_w, w["bias"], a["bias"])
        out[f"n_{kind}"] = np.where(pick_w, w["n"], a["n"])
        if "n_eff" in w and "n_eff" in a:
            out[f"n_eff_{kind}"] = np.where(pick_w, w["n_eff"], a["n_eff"])
        if "n_years" in w and "n_years" in a:
            out[f"n_years_{kind}"] = np.where(pick_w, w["n_years"], a["n_years"])
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    d = np.load("data/scatter_all.npz")
    lat, lon, month = d["lat"].astype("float64"), d["lon"].astype("float64"), d["month"].astype(int)
    depth_all = d["depth"].astype("float64")
    nlev = depth_all.shape[0]
    cell, row, col, inside = profile_grid_indices(lat, lon)
    year = d["year"].astype(int)
    if "time" in d:
        time_sec = d["time"].astype("int64")
    else:
        time_sec = ((year - 1970) * 365.2425 * 86400 + (month - 0.5) * 30.436875 * 86400).astype("int64")
    day_id = time_sec // 86400
    r0_all, r1_all, wr_all, c0_all, c1_all, wc_all = spatial_interp_weights(lat, lon)
    m0_all, m1_all, wm_all, s0_all, s1_all, ws_all = temporal_interp_weights(time_sec)
    if np.any((month < 1) | (month > 12)):
        raise ValueError("Profile month must be in 1..12")
    inside &= year >= YEAR0
    print("profiles", len(lat), "inside grid/domain", int(inside.sum()), flush=True)

    monT, sesT, depth = load_ref("t", "t_an")
    monS, sesS, salinity_depth = load_ref("s", "s_an")
    check_depth_levels(depth_all, depth, "scatter cache")
    check_depth_levels(salinity_depth, depth, "WOA23 salinity files")
    aT, aS = load_aotd("temp", depth), load_aotd("salt", depth)
    To, So = d["temp"], d["salt"]
    z = depth_all[:, None]
    nm_actual = min(NMONTHLY, nlev)
    acc = {
        (rn, v): [
            np.zeros((nlev, NC * NC)),
            np.zeros((nlev, NC * NC)),
            np.zeros((nlev, NC * NC)),
            np.zeros((nlev, NC * NC), dtype="int32"),
            np.zeros((nlev, NC * NC), dtype="uint64"),
        ]
        for rn in ("woa23", "aotd")
        for v in ("potential_temperature", "salinity")
    }
    for m in range(1, 13):
        sel = np.nonzero((month == m) & inside)[0]
        if len(sel) == 0:
            continue
        ce, la = cell[sel], lat[sel]
        yr_bit = np.uint64(1) << np.clip(year[sel] - YEAR0, 0, 62).astype("uint64")
        cd_key = ce.astype("int64") * 100000 + day_id[sel]
        r0, r1, wr, c0, c1, wc = r0_all[sel], r1_all[sel], wr_all[sel], c0_all[sel], c1_all[sel], wc_all[sel]
        m0, m1, wm, s0, s1_i, ws = m0_all[sel], m1_all[sel], wm_all[sel], s0_all[sel], s1_all[sel], ws_all[sel]
        t_o, s_o = To[sel].T.astype("float64"), So[sel].T.astype("float64")
        for rn in ("woa23", "aotd"):
            Tr = np.full((nlev, len(sel)), np.nan)
            Sr = Tr.copy()
            if rn == "woa23":
                Tr[:nm_actual] = interp_ref_2d_time(monT[:, :nm_actual], m0, m1, wm, r0, r1, wr, c0, c1, wc)
                Sr[:nm_actual] = interp_ref_2d_time(monS[:, :nm_actual], m0, m1, wm, r0, r1, wr, c0, c1, wc)
                if nlev > nm_actual:
                    Tr[nm_actual:] = interp_ref_2d_time(sesT, s0, s1_i, ws, r0, r1, wr, c0, c1, wc)
                    Sr[nm_actual:] = interp_ref_2d_time(sesS, s0, s1_i, ws, r0, r1, wr, c0, c1, wc)
            else:
                Tr[:nm_actual] = interp_ref_2d_time(aT[:, :nm_actual], m0, m1, wm, r0, r1, wr, c0, c1, wc)
                Sr[:nm_actual] = interp_ref_2d_time(aS[:, :nm_actual], m0, m1, wm, r0, r1, wr, c0, c1, wc)
            S_for_obs = np.where(np.isnan(s_o), Sr, s_o)
            devs = {"potential_temperature": calculate_potential_T(t_o, S_for_obs, z, la) - calculate_potential_T(Tr, Sr, z, la),
                    "salinity": s_o - Sr}
            for v, dev in devs.items():
                n, s1, s2, neff, ymask = acc[(rn, v)]
                good = np.isfinite(dev)
                for k in range(nlev):
                    g = good[k]
                    if g.any():
                        cg = ce[g]
                        n[k] += np.bincount(cg, minlength=NC * NC)
                        s1[k] += np.bincount(cg, weights=dev[k, g], minlength=NC * NC)
                        s2[k] += np.bincount(cg, weights=dev[k, g] ** 2, minlength=NC * NC)
                        uniq_cd = np.unique(cd_key[g])
                        neff[k] += np.bincount(uniq_cd // 100000, minlength=NC * NC).astype("int32")
                        np.bitwise_or.at(ymask[k], cg, yr_bit[g])
        print("month", m, "done", flush=True)

    shp = (nlev, NC, NC)
    metrics = {}  # (product, var) -> dict(std, rmse, bias, n, n_eff, n_years)
    for (rn, v), (n, s1, s2, neff, ymask) in acc.items():
        nyears = np.bitwise_count(ymask).astype("int32")
        metrics[(rn, v)] = deviation_statistics(n, s1, s2, n_eff=neff, n_years=nyears)
    for v in ("potential_temperature", "salinity"):
        w, a = metrics[("woa23", v)], metrics[("aotd", v)]
        metrics[("woa23_and_aotd", v)] = combine_statistics(w, a)

    xc = -HALF + CELL * (np.arange(NC) + 0.5)
    yc = HALF - CELL * (np.arange(NC) + 0.5)
    X, Y = np.meshgrid(xc, yc)
    glon, glat = Transformer.from_crs(6931, 4326, always_xy=True).transform(X, Y)
    titles = {"woa23": "WOA23", "aotd": "AOTD", "woa23_and_aotd": "the larger of the WOA23-based and AOTD-based estimates"}
    for (rn, v), mt in metrics.items():
        short = "Theta" if v == "potential_temperature" else "S"
        unit = "degree_Celsius" if v == "potential_temperature" else "1e-3"
        pname = v.replace("_", " ")
        for kind, desc in (("std", "sample standard deviation (ddof=1) about the cell mean offset"), ("rmse", "RMSE (about the reference, bias included)")):
            bias = mt[f"bias_{kind}"] if rn == "woa23_and_aotd" else mt["bias"]
            n = mt[f"n_{kind}"] if rn == "woa23_and_aotd" else mt["n"]
            n_eff = mt[f"n_eff_{kind}"] if rn == "woa23_and_aotd" else mt["n_eff"]
            n_years = mt[f"n_years_{kind}"] if rn == "woa23_and_aotd" else mt["n_years"]
            ds = xr.Dataset(
                {
                    f"{short}_{kind}": (("depth", "y", "x"), mt[kind].reshape(shp).astype("float32"),
                                        dict(long_name=f"expected {desc} of in situ {pname} from {titles[rn]}", units=unit,
                                             grid_mapping="crs", coordinates="lat lon",
                                             comment="Deviation = obs - reference (bilinearly interpolated in space and linearly in time), both as potential temperature (T only), pooled per cell 1992-2023, no time dimension. "
                                                     f"Product: {rn}. Sample std uses ddof=1 (NaN when n < 2); RMSE is NaN when n < 1; AOTD reference exists only to 1500 m.")),
                    f"{short}_bias": (("depth", "y", "x"), bias.reshape(shp).astype("float32"),
                                      dict(long_name=f"mean in situ {pname} minus reference", units=unit, grid_mapping="crs", coordinates="lat lon")),
                    f"{short}_n": (("depth", "y", "x"), n.reshape(shp).astype("int32"),
                                   dict(long_name="number of finite in situ comparison values used", units="1", grid_mapping="crs", coordinates="lat lon")),
                    f"{short}_n_eff": (("depth", "y", "x"), n_eff.reshape(shp).astype("int32"),
                                       dict(long_name="number of distinct (cell, UTC day) sampling days with valid comparisons (de-clustered count proxy)", units="1", grid_mapping="crs", coordinates="lat lon")),
                    f"{short}_n_years": (("depth", "y", "x"), n_years.reshape(shp).astype("int32"),
                                         dict(long_name="number of distinct calendar years (1992-2023) with valid comparisons", units="1", grid_mapping="crs", coordinates="lat lon")),
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
            ds.to_netcdf(f, encoding={k: dict(zlib=True, complevel=4) for k in (f"{short}_{kind}", f"{short}_bias", f"{short}_n", f"{short}_n_eff", f"{short}_n_years")})
            print("wrote", f, flush=True)


if __name__ == "__main__":
    main()
