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
    rmse    = sqrt(mean(d^2))                                profile-weighted spread about the reference, bias included (NaN when n < 1)
    std     = sqrt(n / (n - 1) * (mean(d^2) - mean(d)^2))    profile-weighted sample standard deviation about the cell's mean offset (ddof=1; NaN when n < 2)
    bias    = mean(d), n = number of finite comparisons
    n_eff   = number of distinct (100 km cell, UTC calendar day) sampling days with valid comparisons (de-clustered count proxy)
    n_years = number of distinct calendar years (1992-2023) with valid comparisons
    n_epochs = M, number of distinct 1-month synoptic epochs (cell, year, month) with valid comparisons
    std_declustered, rmse_declustered, bias_declustered, std_within, std_between = epoch-declustered (1-month synoptic epoch) statistics
    std_se, rmse_se = Empirical Bayes regularized (nu_0 = 5 epochs, 500 km local moving-window prior) epoch-declustered spread (defined for all M >= 1)
    std_se_pan_arctic, rmse_se_pan_arctic = Empirical Bayes regularized variant using a single pan-Arctic depth-level prior (for comparison with the 500 km local prior)

Usage (from the repo root): python scripts/03_compute_expected_spread.py   (about 1.5 min, about 4 GB RAM)
Outputs (untracked): outputs/expected_spread/expected_{std,rmse}_from_{woa23,aotd,woa23_and_aotd}_{potential_temperature,salinity}_ease2_100km.nc
"""
import math
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


def local_500km_prior(var_tot_dec, M, radius_cells=5):
    """Compute 500 km circular moving-window median prior on the 2D grid at each depth, plus pan-Arctic scalar prior."""
    v2d = var_tot_dec.reshape(var_tot_dec.shape[0], -1) if var_tot_dec.ndim > 1 else var_tot_dec[None, :]
    m2d = M.reshape(v2d.shape)
    nlev, ncells = v2d.shape
    pan_prior = np.zeros((nlev, 1), dtype="float64")
    for k in range(nlev):
        well = v2d[k][(m2d[k] >= 5) & np.isfinite(v2d[k])]
        if well.size == 0:
            well = v2d[k][np.isfinite(v2d[k])]
        pan_prior[k, 0] = float(np.median(well)) if well.size > 0 else np.nan
    side = int(math.isqrt(ncells))
    if side * side != ncells or side <= 2 * radius_cells:
        out_pan = np.broadcast_to(pan_prior, v2d.shape).copy()
        return (out_pan[0] if var_tot_dec.ndim == 1 else out_pan), (out_pan[0] if var_tot_dec.ndim == 1 else out_pan)

    offsets = [(dy, dx) for dy in range(-radius_cells, radius_cells + 1)
               for dx in range(-radius_cells, radius_cells + 1) if dy * dy + dx * dx <= radius_cells * radius_cells]
    R = radius_cells
    v_grid = v2d.reshape(nlev, side, side)
    m_grid = m2d.reshape(nlev, side, side)
    v_pad = np.pad(v_grid, ((0, 0), (R, R), (R, R)), constant_values=np.nan)
    m_pad = np.pad(m_grid, ((0, 0), (R, R), (R, R)), constant_values=0)

    stack_v = np.empty((len(offsets), nlev, side, side), dtype="float64")
    stack_m = np.empty((len(offsets), nlev, side, side), dtype="float64")
    for idx, (dy, dx) in enumerate(offsets):
        stack_v[idx] = v_pad[:, R + dy:R + dy + side, R + dx:R + dx + side]
        stack_m[idx] = m_pad[:, R + dy:R + dy + side, R + dx:R + dx + side]

    has_target = (m_grid >= 1).reshape(nlev, -1)
    stack_v = stack_v.reshape(len(offsets), nlev, -1)
    stack_m = stack_m.reshape(len(offsets), nlev, -1)
    loc_prior = np.broadcast_to(pan_prior, (nlev, ncells)).copy()
    with np.errstate(invalid="ignore"):
        for k in range(nlev):
            t_idx = np.nonzero(has_target[k])[0]
            if t_idx.size == 0:
                continue
            sv = stack_v[:, k, t_idx]
            sm = stack_m[:, k, t_idx]
            v5 = np.where((sm >= 5) & np.isfinite(sv), sv, np.nan)
            c5 = np.isfinite(v5).sum(axis=0)
            med5 = np.full(t_idx.size, np.nan)
            if (c5 >= 3).any():
                med5[c5 >= 3] = np.nanmedian(v5[:, c5 >= 3], axis=0)
            v2 = np.where((sm >= 2) & np.isfinite(sv), sv, np.nan)
            c2 = np.isfinite(v2).sum(axis=0)
            med2 = np.full(t_idx.size, np.nan)
            if (c2 >= 1).any():
                med2[c2 >= 1] = np.nanmedian(v2[:, c2 >= 1], axis=0)
            loc_prior[k, t_idx] = np.where(c5 >= 3, med5, np.where(c2 >= 1, med2, pan_prior[k, 0]))
    out_pan = np.broadcast_to(pan_prior, (nlev, ncells)).copy()
    if var_tot_dec.ndim == 1:
        return loc_prior[0], out_pan[0]
    return loc_prior, out_pan


def deviation_statistics(n, s1, s2, n_eff=None, n_years=None, m_ep=None, s1_ep=None, s2_ep=None, s_inv_ne=None, s_ssw=None, nu0=5.0):
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
        if m_ep is not None:
            M = m_ep.astype("float64")
            mu_dec = np.where(M >= 1, s1_ep / M, np.nan)
            var_within = np.where(n > M, np.maximum(s_ssw / (n - M), 0.0), np.nan)
            s2_dbar = np.where(M >= 2, np.maximum((s2_ep - M * (mu_dec ** 2)) / (M - 1.0), 0.0), np.nan)
            var_between = np.where(
                M >= 2,
                np.maximum(s2_dbar - np.where(np.isfinite(var_within), var_within * (s_inv_ne / M), 0.0), 0.0),
                np.nan,
            )
            var_tot_dec = np.where(
                M >= 2,
                np.maximum(s2_dbar + np.where(np.isfinite(var_within), var_within * (1.0 - s_inv_ne / M), 0.0), 0.0),
                np.nan,
            )
            rmse_dec = np.where(
                M >= 2,
                np.sqrt(np.maximum(((M - 1.0) / M) * var_tot_dec + mu_dec ** 2, 0.0)),
                np.where(M == 1, np.sqrt(np.maximum(mu_dec ** 2 + np.where(n > 1, s_ssw / n, 0.0), 0.0)), np.nan),
            )
            # Empirical Bayes shrinkage (nu_0 = 5 epochs) with 500 km local moving-window prior and pan-Arctic scalar prior
            loc_prior_var, pan_prior_var = local_500km_prior(var_tot_dec, M, radius_cells=5)
            local_var_m1 = np.where(np.isfinite(var_within), np.maximum(var_within, loc_prior_var), loc_prior_var)
            var_se = np.where(
                M >= 2,
                (nu0 * loc_prior_var + (M - 1.0) * var_tot_dec) / (nu0 + (M - 1.0)),
                np.where(M == 1, local_var_m1, np.nan),
            )
            rmse_se = np.where(M >= 1, np.sqrt(np.maximum(var_se + mu_dec ** 2, 0.0)), np.nan)

            pan_var_m1 = np.where(np.isfinite(var_within), np.maximum(var_within, pan_prior_var), pan_prior_var)
            var_se_pan = np.where(
                M >= 2,
                (nu0 * pan_prior_var + (M - 1.0) * var_tot_dec) / (nu0 + (M - 1.0)),
                np.where(M == 1, pan_var_m1, np.nan),
            )
            rmse_se_pan = np.where(M >= 1, np.sqrt(np.maximum(var_se_pan + mu_dec ** 2, 0.0)), np.nan)
            out.update(
                n_epochs=m_ep.astype("int32"),
                bias_declustered=mu_dec,
                std_within=np.sqrt(var_within),
                std_between=np.sqrt(var_between),
                std_declustered=np.sqrt(var_tot_dec),
                rmse_declustered=rmse_dec,
                std_se=np.sqrt(np.maximum(var_se, 0.0)),
                rmse_se=rmse_se,
                std_se_pan_arctic=np.sqrt(np.maximum(var_se_pan, 0.0)),
                rmse_se_pan_arctic=rmse_se_pan,
            )
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
        for extra in ("n_eff", "n_years", "n_epochs", "bias_declustered", "std_within", "std_between"):
            if extra in w and extra in a:
                out[f"{extra}_{kind}"] = np.where(pick_w, w[extra], a[extra])
        for prefix in ("declustered", "se", "se_pan_arctic"):
            mk = f"{kind}_{prefix}"
            if mk in w and mk in a:
                out[mk] = np.fmax(w[mk], a[mk])
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
            np.zeros((nlev, NC * NC), dtype="int32"),
            np.zeros((nlev, NC * NC)),
            np.zeros((nlev, NC * NC)),
            np.zeros((nlev, NC * NC)),
            np.zeros((nlev, NC * NC)),
        ]
        for rn in ("woa23", "aotd")
        for v in ("potential_temperature", "salinity")
    }
    for m in range(1, 13):
        sel = np.nonzero((month == m) & inside)[0]
        if len(sel) == 0:
            continue
        ce, la, yr_sel = cell[sel], lat[sel], year[sel] - YEAR0
        yr_bit = np.uint64(1) << np.clip(yr_sel, 0, 62).astype("uint64")
        cd_key = ce.astype("int64") * 100000 + day_id[sel]
        ep_key = ce.astype("int64") * 64 + yr_sel.astype("int64")
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
                n, s1, s2, neff, ymask, m_ep, s1_ep, s2_ep, s_inv_ne, s_ssw = acc[(rn, v)]
                good = np.isfinite(dev)
                for k in range(nlev):
                    g = good[k]
                    if g.any():
                        cg = ce[g]
                        dg = dev[k, g]
                        n[k] += np.bincount(cg, minlength=NC * NC)
                        s1[k] += np.bincount(cg, weights=dg, minlength=NC * NC)
                        s2[k] += np.bincount(cg, weights=dg ** 2, minlength=NC * NC)
                        uniq_cd = np.unique(cd_key[g])
                        neff[k] += np.bincount(uniq_cd // 100000, minlength=NC * NC).astype("int32")
                        np.bitwise_or.at(ymask[k], cg, yr_bit[g])
                        uniq_ep, inv_ep, n_e = np.unique(ep_key[g], return_inverse=True, return_counts=True)
                        ep_cells = uniq_ep // 64
                        s1_e = np.bincount(inv_ep, weights=dg)
                        s2_e = np.bincount(inv_ep, weights=dg ** 2)
                        dbar_e = s1_e / n_e
                        ssw_e = np.maximum(s2_e - n_e * (dbar_e ** 2), 0.0)
                        m_ep[k] += np.bincount(ep_cells, minlength=NC * NC).astype("int32")
                        s1_ep[k] += np.bincount(ep_cells, weights=dbar_e, minlength=NC * NC)
                        s2_ep[k] += np.bincount(ep_cells, weights=dbar_e ** 2, minlength=NC * NC)
                        s_inv_ne[k] += np.bincount(ep_cells, weights=1.0 / n_e, minlength=NC * NC)
                        s_ssw[k] += np.bincount(ep_cells, weights=ssw_e, minlength=NC * NC)
        print("month", m, "done", flush=True)

    shp = (nlev, NC, NC)
    metrics = {}
    for (rn, v), (n, s1, s2, neff, ymask, m_ep, s1_ep, s2_ep, s_inv_ne, s_ssw) in acc.items():
        nyears = np.bitwise_count(ymask).astype("int32")
        metrics[(rn, v)] = deviation_statistics(n, s1, s2, n_eff=neff, n_years=nyears, m_ep=m_ep, s1_ep=s1_ep, s2_ep=s2_ep, s_inv_ne=s_inv_ne, s_ssw=s_ssw)
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
            sfx = f"_{kind}" if rn == "woa23_and_aotd" else ""
            bias = mt[f"bias{sfx}"]
            n = mt[f"n{sfx}"]
            n_eff = mt[f"n_eff{sfx}"]
            n_years = mt[f"n_years{sfx}"]
            n_epochs = mt[f"n_epochs{sfx}"]
            bias_dec = mt[f"bias_declustered{sfx}"]
            std_w = mt[f"std_within{sfx}"]
            std_b = mt[f"std_between{sfx}"]
            val_dec = mt[f"{kind}_declustered"]
            val_se = mt[f"{kind}_se"]
            val_se_pan = mt[f"{kind}_se_pan_arctic"]
            ds = xr.Dataset(
                {
                    f"{short}_{kind}": (("depth", "y", "x"), mt[kind].reshape(shp).astype("float32"),
                                        dict(long_name=f"profile-weighted expected {desc} of in situ {pname} from {titles[rn]}", units=unit,
                                             grid_mapping="crs", coordinates="lat lon",
                                             comment="Deviation = obs - reference (bilinearly interpolated in space and linearly in time), both as potential temperature (T only), pooled per cell 1992-2023, no time dimension. "
                                                     f"Product: {rn}. Sample std uses ddof=1 (NaN when n < 2); RMSE is NaN when n < 1; AOTD reference exists only to 1500 m.")),
                    f"{short}_{kind}_declustered": (("depth", "y", "x"), val_dec.reshape(shp).astype("float32"),
                                                    dict(long_name=f"epoch-declustered (1-month synoptic epoch) {desc} of in situ {pname} from {titles[rn]}", units=unit,
                                                         grid_mapping="crs", coordinates="lat lon",
                                                         comment="Equal weight per (cell, year, month) synoptic epoch plus random-effects within-epoch variance restoration. NaN when n_epochs < 2 for std (n_epochs < 1 for RMSE).")),
                    f"{short}_{kind}_se": (("depth", "y", "x"), val_se.reshape(shp).astype("float32"),
                                           dict(long_name=f"Empirical Bayes regularized (500 km local moving-window prior, nu_0=5) {kind} of in situ {pname} from {titles[rn]}", units=unit,
                                                grid_mapping="crs", coordinates="lat lon",
                                                comment="Shrinks epoch-declustered variance toward the 500 km circular moving-window median prior across neighboring cells with n_epochs >= 5 (prior strength nu_0 = 5 epochs). Defined for all sampled cells with n_epochs >= 1.")),
                    f"{short}_{kind}_se_pan_arctic": (("depth", "y", "x"), val_se_pan.reshape(shp).astype("float32"),
                                                      dict(long_name=f"Empirical Bayes regularized (single pan-Arctic depth-level prior, nu_0=5) {kind} of in situ {pname} from {titles[rn]}", units=unit,
                                                           grid_mapping="crs", coordinates="lat lon",
                                                           comment="Shrinks epoch-declustered variance toward the single pan-Arctic depth-level median prior across all cells with n_epochs >= 5 (nu_0 = 5 epochs). Provided for comparison with the 500 km local prior.")),
                    f"{short}_std_within": (("depth", "y", "x"), std_w.reshape(shp).astype("float32"),
                                            dict(long_name=f"within-epoch (sub-monthly / sub-grid / instrument) standard deviation of in situ {pname}", units=unit,
                                                 grid_mapping="crs", coordinates="lat lon")),
                    f"{short}_std_between": (("depth", "y", "x"), std_b.reshape(shp).astype("float32"),
                                             dict(long_name=f"between-epoch (synoptic / seasonal-anomaly / interannual) standard deviation of in situ {pname}", units=unit,
                                                  grid_mapping="crs", coordinates="lat lon")),
                    f"{short}_bias": (("depth", "y", "x"), bias.reshape(shp).astype("float32"),
                                      dict(long_name=f"profile-weighted mean in situ {pname} minus reference", units=unit, grid_mapping="crs", coordinates="lat lon")),
                    f"{short}_bias_declustered": (("depth", "y", "x"), bias_dec.reshape(shp).astype("float32"),
                                                  dict(long_name=f"epoch-weighted (1-month de-clustered) mean in situ {pname} minus reference", units=unit, grid_mapping="crs", coordinates="lat lon")),
                    f"{short}_n": (("depth", "y", "x"), n.reshape(shp).astype("int32"),
                                   dict(long_name="number of finite in situ comparison values used", units="1", grid_mapping="crs", coordinates="lat lon")),
                    f"{short}_n_eff": (("depth", "y", "x"), n_eff.reshape(shp).astype("int32"),
                                       dict(long_name="number of distinct (cell, UTC day) sampling days with valid comparisons", units="1", grid_mapping="crs", coordinates="lat lon")),
                    f"{short}_n_epochs": (("depth", "y", "x"), n_epochs.reshape(shp).astype("int32"),
                                          dict(long_name="number of distinct 1-month synoptic epochs (cell, year, month) with valid comparisons", units="1", grid_mapping="crs", coordinates="lat lon")),
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
            comp_vars = (f"{short}_{kind}", f"{short}_{kind}_declustered", f"{short}_{kind}_se", f"{short}_{kind}_se_pan_arctic",
                         f"{short}_std_within", f"{short}_std_between", f"{short}_bias", f"{short}_bias_declustered",
                         f"{short}_n", f"{short}_n_eff", f"{short}_n_epochs", f"{short}_n_years")
            ds.to_netcdf(f, encoding={k: dict(zlib=True, complevel=4) for k in comp_vars})
            print("wrote", f, flush=True)


if __name__ == "__main__":
    main()
