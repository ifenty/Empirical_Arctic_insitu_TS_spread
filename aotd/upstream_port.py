"""Faithful Python port of the upstream MATLAB climatology code (cressman_bkg.m, barnes_season.m, barnes_month.m).

This module deliberately KEEPS the upstream behaviour, including the issues listed in docs/upstream_bug_log.md
(BUG-18 ... BUG-25), so that the published gridded product can be reproduced. Use `aotd.improved` for the corrected
method.

Array convention: (lat, lon) = (121, 1440), i.e. MATLAB's (lon, lat) arrays transposed (the published NetCDF files
store temp as (depth, lat, lon) when opened in Python).
"""
from __future__ import annotations

import numpy as np

from .geo import upstream_grid
from .neighbors import SphereNeighbors

NY, NX = 121, 1440
LAT_1D, LON_1D = upstream_grid()
LON2D, LAT2D = np.meshgrid(LON_1D, LAT_1D)  # (lat, lon)

CRESSMAN_RADII_KM = (321.0, 267.0, 214.0)
BARNES_RADIUS_KM = 555.0  # barnes_month.m (barnes_season.m has an unused 999)
BARNES_N = (8.0e4, 1.6e4)  # km^2 ("alpha" in the paper), weight exp(-r^2 / (0.2 n))
SEASON_OF_MONTH = np.array([1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4])  # JFM, AMJ, JAS, OND


# ----------------------------------------------------------------------------- binning (obs_merge)
def obs_merge(lat, lon, values, nanmean: bool):
    """Average observations onto grid nodes: every obs is assigned to ALL nodes within 0.25 deg Euclidean distance
    in (lat, lon) degrees (BUG-21); only the lon == -180 column handles the dateline.

    values: list of 1-D arrays aligned with lat/lon. Returns (list of (NY,NX) means, (NY,NX) member count).
    nanmean=False reproduces MATLAB `mean` (a single NaN member makes the cell NaN, BUG-22).
    """
    lat = np.asarray(lat, dtype="float64")
    lon = np.asarray(lon, dtype="float64")
    i0 = np.clip(np.floor((lat - 60.0) / 0.25).astype(int), 0, NY - 2)
    j0 = np.clip(np.floor((lon + 180.0) / 0.25).astype(int), 0, NX - 1)
    sums = [np.zeros(NY * NX) for _ in values]
    cnts = [np.zeros(NY * NX) for _ in values]
    num = np.zeros(NY * NX)
    for di in (0, 1):
        for dj in (0, 1):
            ii = i0 + di
            jj = j0 + dj
            wrap = jj >= NX
            jj = np.where(wrap, 0, jj)
            node_lat = 60.0 + 0.25 * ii
            node_lon = -180.0 + 0.25 * jj
            dlon = np.abs(lon - node_lon)
            quirk = jj == 0  # lon == -180 node: compares with -180 for lon<0 and +180 for lon>=0
            dlon = np.where(quirk, np.abs(lon - np.where(lon >= 0, 180.0, -180.0)), dlon)
            inside = np.hypot(lat - node_lat, dlon) <= 0.25
            flat = (ii * NX + jj)[inside]
            num += np.bincount(flat, minlength=NY * NX)
            for k, v in enumerate(values):
                vv = np.asarray(v, dtype="float64")[inside]
                good = ~np.isnan(vv)
                if nanmean:
                    sums[k] += np.bincount(flat[good], weights=vv[good], minlength=NY * NX)
                    cnts[k] += np.bincount(flat[good], minlength=NY * NX)
                else:
                    sums[k] += np.bincount(flat, weights=np.where(good, vv, np.nan), minlength=NY * NX)
                    cnts[k] += np.bincount(flat, minlength=NY * NX)
    out = []
    for k in range(len(values)):
        with np.errstate(invalid="ignore", divide="ignore"):
            m = np.where(cnts[k] > 0, sums[k] / cnts[k], np.nan)
        out.append(m.reshape(NY, NX))
    return out, num.reshape(NY, NX)


# ----------------------------------------------------------------------------- weighted sums
def _kernel_cressman(r, R):
    return (R * R - r * r) / (R * R + r * r)


def _kernel_barnes(r, n):
    return np.exp(-(r * r) / (0.2 * n))


def _accumulate(obs_grid, R, kernel, kparam):
    """For every grid node: sum(w), sum(w * obs), count over obs-grid nodes (non-NaN) closer than R km."""
    ok = ~np.isnan(obs_grid)
    nb = SphereNeighbors(LON2D[ok], LAT2D[ok], LON2D, LAT2D)
    vals = obs_grid[ok]
    n = NY * NX
    sw = np.zeros(n)
    swv = np.zeros(n)
    cnt = np.zeros(n)
    for t0, t1, i, j, r in nb.pairs(R):
        w = kernel(r, kparam)
        m = t1 - t0
        sw[t0:t1] = np.bincount(j, weights=w, minlength=m)
        swv[t0:t1] = np.bincount(j, weights=w * vals[i], minlength=m)
        cnt[t0:t1] = np.bincount(j, minlength=m)
    return sw.reshape(NY, NX), swv.reshape(NY, NX), cnt.reshape(NY, NX)


def _damped_denominator(sw, cnt):
    return sw + (cnt == 1)  # BUG-19: ad hoc "+1" when exactly one obs is in range


def cressman_background(obs, ok_level, land, R=CRESSMAN_RADII_KM[0]):
    """First pass of cressman_bkg.m: distance-weighted mean at ALL valid cells (not a correction)."""
    sw, swv, cnt = _accumulate(obs, R, _kernel_cressman, R)
    z = np.full((NY, NX), np.nan)
    sel = ok_level & (cnt > 0)
    z[sel] = swv[sel] / _damped_denominator(sw, cnt)[sel]
    z[land] = np.nan
    return z


def correction_pass(z, obs, ok_level, land, R, kernel, kparam, only_obs_cells: bool):
    """z <- z + sum(w (obs_b - z_target)) / sum(w)  (BUG-18: residual taken at the TARGET node, not at obs b)."""
    sw, swv, cnt = _accumulate(obs, R, kernel, kparam)
    sel = ok_level & (cnt > 0)
    if only_obs_cells:
        sel &= ~np.isnan(obs)  # BUG-20
    znew = z.copy()
    num = swv - z * sw
    znew[sel] = z[sel] + num[sel] / _damped_denominator(sw, cnt)[sel]
    znew[land] = np.nan
    return znew


# ----------------------------------------------------------------------------- one level, all periods
def annual_field(obs, ok_level, land):
    z = cressman_background(obs, ok_level, land, CRESSMAN_RADII_KM[0])
    for R in CRESSMAN_RADII_KM[1:]:
        z = correction_pass(z, obs, ok_level, land, R, _kernel_cressman, R, only_obs_cells=True)
    z[land] = np.nan
    return z


def seasonal_field(annual, obs, ok_level, land):
    z = annual.copy()
    z[land] = np.nan
    for R in CRESSMAN_RADII_KM:
        z = correction_pass(z, obs, ok_level, land, R, _kernel_cressman, R, only_obs_cells=True)
    z[np.isnan(annual)] = np.nan
    z[land] = np.nan
    return z


def monthly_field(seasonal, obs, ok_level):
    z = seasonal.copy()
    for n in BARNES_N:
        sw, swv, cnt = _accumulate(obs, BARNES_RADIUS_KM, _kernel_barnes, n)
        znew = z.copy()
        sel = cnt > 0
        znew[sel] = z[sel] + (swv - z * sw)[sel] / _damped_denominator(sw, cnt)[sel]
        znew[~ok_level] = np.nan
        z = znew
    return z


def process_level(k, data, ok_level, land):
    """Run annual -> seasonal -> monthly for temperature and salinity at level index k.

    data: dict with lat, lon, month (1-D), temp, salt (n x >=k+1). Returns dict of arrays
    {'annual': (2,NY,NX), 'season': (4,2,NY,NX), 'month': (12,2,NY,NX)} with var order (temp, salt).
    """
    lat, lon, month = data["lat"], data["lon"], data["month"]
    t_all = data["temp"][:, k].astype("float64")
    s_all = data["salt"][:, k].astype("float64")

    sel = ~np.isnan(t_all) | ~np.isnan(s_all)  # cressman_bkg.m
    (t_m, s_m), _ = obs_merge(lat[sel], lon[sel], [t_all[sel], s_all[sel]], nanmean=True)
    annual = []
    for obs in (t_m, s_m):
        obs = obs.copy()
        obs[land] = np.nan
        annual.append(annual_field(obs, ok_level, land))
    annual = np.stack(annual)

    def merged_period(idx):
        sel = idx & ~np.isnan(t_all)  # barnes_season/month.m: only T-valid obs; S taken at the same obs
        (t_m, s_m), _ = obs_merge(lat[sel], lon[sel], [t_all[sel], s_all[sel]], nanmean=False)
        out = []
        for obs in (t_m, s_m):
            obs = obs.copy()
            obs[land] = np.nan
            out.append(obs)
        return out

    season = np.full((4, 2, NY, NX), np.nan)
    for sidx in range(4):
        obs_t, obs_s = merged_period(SEASON_OF_MONTH[month - 1] == sidx + 1)
        for v, obs in enumerate((obs_t, obs_s)):
            season[sidx, v] = seasonal_field(annual[v], obs, ok_level, land)

    mon = np.full((12, 2, NY, NX), np.nan)
    for m in range(12):
        obs_t, obs_s = merged_period(month == m + 1)
        for v, obs in enumerate((obs_t, obs_s)):
            mon[m, v] = monthly_field(season[SEASON_OF_MONTH[m] - 1, v], obs, ok_level)
    return {"annual": annual, "season": season, "month": mon}
