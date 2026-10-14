"""Data-honest optimal-interpolation (OI) climatology with formal error estimates.

Staged like the upstream method (annual -> seasonal -> monthly), but every stage is a local OI of the residuals from
the previous stage with a Gaussian covariance, an explicit observation-noise model and a posterior error variance.
Far from data the estimate relaxes to the previous stage (no nearest-observation polygons), and cells without
enough nearby observations are set to NaN. Time is pooled (all years equal weight): see HANDOFF.md, temporal item.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from .upstream_port import NX, NY, LAT2D, LON2D, SEASON_OF_MONTH

R_EARTH = 6371.0
PARAMS = dict(
    L=(300.0, 200.0, 150.0),   # km e-folding scale of exp(-(r/L)^2): annual, seasonal, monthly
    K=24,                      # nearest observation nodes used per target cell
    rmax_fac=2.0,              # only nodes within rmax_fac * L are used
    nsr_a=0.1, nsr_b=0.3,      # obs noise variance = sigma^2 * (a + b / n_profiles_in_node)
    dmax=150.0,                # km: NaN if the nearest observation node (any month) is further away
    relerr_max=0.25,           # NaN if calibrated error variance exceeds this fraction of the level's total variance
)


def _xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


TARGET_XYZ = _xyz(LAT2D.ravel(), LON2D.ravel())


def node_table(flat, val, sel):
    """Mean and profile count per grid node (nearest node) for profiles in `sel`. Returns (node idx, mean, n)."""
    f, v = flat[sel], val[sel]
    good = ~np.isnan(v)
    f, v = f[good], v[good]
    cnt = np.bincount(f, minlength=NY * NX)
    sm = np.bincount(f, weights=v, minlength=NY * NX)
    idx = np.nonzero(cnt)[0]
    return idx, sm[idx] / cnt[idx], cnt[idx].astype(float)


def _chord(km):
    return 2.0 * np.sin(km / (2.0 * R_EARTH))


def oi_stage(idx, resid, n, targets, L, p=PARAMS, chunk=4000):
    """Local OI of node residuals. Returns increment, posterior variance, sigma^2 (prior), nearest-node distance km.

    targets: indices into the flattened grid. Prior mean of the residual is 0.
    """
    nt = len(targets)
    if len(idx) < 3:
        return np.zeros(nt), np.full(nt, np.nan), np.nan, np.full(nt, np.inf)
    nsr = p["nsr_a"] + p["nsr_b"] / n
    sig2 = max(np.var(resid) / (1.0 + nsr.mean()), 1e-12)
    pos = TARGET_XYZ[idx]
    tree = cKDTree(pos)
    K = min(p["K"], len(idx))
    rmax = p["rmax_fac"] * L
    d, nb = tree.query(TARGET_XYZ[targets], k=K, distance_upper_bound=_chord(rmax))
    miss = ~np.isfinite(d)
    nb = np.where(miss, 0, nb)
    dmin = 2.0 * R_EARTH * np.arcsin(np.minimum(d[:, 0], 2.0) / 2.0)
    dmin = np.where(miss[:, 0], np.inf, dmin)
    inc = np.zeros(nt)
    post = np.full(nt, sig2)
    for a in range(0, nt, chunk):
        b = min(a + chunk, nt)
        nbc, ms = nb[a:b], miss[a:b]
        P = pos[nbc]                                            # (m,K,3)
        ch = np.linalg.norm(P[:, :, None, :] - P[:, None, :, :], axis=-1)
        rr = 2.0 * R_EARTH * np.arcsin(np.minimum(ch, 2.0) / 2.0)
        C = sig2 * np.exp(-(rr / L) ** 2)
        c = sig2 * np.exp(-((2.0 * R_EARTH * np.arcsin(np.minimum(d[a:b], 2.0) / 2.0)) / L) ** 2)
        c = np.where(ms, 0.0, c)
        C = np.where(ms[:, :, None] | ms[:, None, :], 0.0, C)
        C = C + (sig2 * nsr[nbc] + np.where(ms, 1.0, 0.0) * sig2)[:, :, None] * np.eye(K)[None]
        w = np.linalg.solve(C, c[:, :, None])[:, :, 0]
        inc[a:b] = (w * np.where(ms, 0.0, resid[nbc])).sum(1)
        post[a:b] = sig2 - (w * c).sum(1)
    return inc, np.maximum(post, 0.0), sig2, dmin


def run_variable(flat, val, month, ok, exclude=None, p=PARAMS):
    """Staged OI for one variable at one level.

    flat: node index of each profile; val: values; month: 1..12; ok: (NY*NX,) bool valid ocean cells at this level;
    exclude: optional bool per profile to withhold (cross-validation).
    Returns dict(field (12,NY,NX), std (12,NY,NX), dmin (NY,NX), tot_var float).
    """
    use = np.ones(len(val), bool) if exclude is None else ~exclude
    tgt = np.nonzero(ok)[0]
    # annual
    idx, m, n = node_table(flat, val, use)
    mu = np.average(m, weights=None)
    tot_var = float(np.var(m))
    inc, P_a, s2_a, dmin = oi_stage(idx, m - mu, n, tgt, p["L"][0], p)
    A = np.full(NY * NX, np.nan); A[tgt] = mu + inc
    # nodes may lie outside `ok` (e.g. in cells masked at this level); drop them
    A_full = A.copy()
    PA = np.full(NY * NX, np.nan); PA[tgt] = P_a
    field = np.full((12, NY * NX), np.nan)
    err2 = np.full((12, NY * NX), np.nan)
    seas = {}
    for s in range(1, 5):
        sel = np.isin(month, np.nonzero(SEASON_OF_MONTH == s)[0] + 1) & use
        idx, m, n = node_table(flat, val, sel)
        keep = ok[idx]
        idx, m, n = idx[keep], m[keep], n[keep]
        inc, P_s, s2_s, _ = oi_stage(idx, m - A_full[idx], n, tgt, p["L"][1], p)
        S = np.full(NY * NX, np.nan); S[tgt] = A_full[tgt] + (inc if len(idx) >= 3 else 0)
        Pst = np.full(NY * NX, np.nan)
        if np.isfinite(s2_s):  # seasonal error carried in proportion to how little the data reduced the prior variance
            Pst[tgt] = P_s + (P_s / s2_s) * PA[tgt]
        else:
            Pst[tgt] = PA[tgt]
        seas[s] = (S, Pst, s2_s)
    s2_months = []
    nodata = []
    for mo in range(1, 13):
        S, Pst, s2_s = seas[SEASON_OF_MONTH[mo - 1]]
        sel = (month == mo) & use
        idx, m, n = node_table(flat, val, sel)
        keep = ok[idx]
        idx, m, n = idx[keep], m[keep], n[keep]
        inc, P_m, s2_m, _ = oi_stage(idx, m - S[idx], n, tgt, p["L"][2], p)
        if len(idx) >= 3:
            field[mo - 1, tgt] = S[tgt] + inc
            err2[mo - 1, tgt] = P_m + (P_m / s2_m) * Pst[tgt]
            s2_months.append(s2_m)
        else:  # no monthly information: seasonal field with the typical monthly anomaly variance added
            field[mo - 1, tgt] = S[tgt]
            err2[mo - 1, tgt] = Pst[tgt]
            nodata.append(mo - 1)
    for mi in nodata:
        err2[mi, tgt] += np.median(s2_months) if s2_months else 0.0
    return dict(field=field.reshape(12, NY, NX), std=np.sqrt(err2).reshape(12, NY, NX),
                dmin=_scatter(dmin, tgt), tot_var=tot_var)


def _scatter(v, tgt):
    out = np.full(NY * NX, np.nan)
    out[tgt] = v
    return out.reshape(NY, NX)


def apply_mask(res, inflate=1.0, p=PARAMS):
    """NaN where the nearest observation node is beyond dmax or the calibrated error variance exceeds
    relerr_max * total variance. `inflate` is the empirical cross-validation calibration of the formal std."""
    bad = ~(res["dmin"] <= p["dmax"])
    f = np.where(bad[None], np.nan, res["field"])
    s = np.where(bad[None], np.nan, res["std"] * inflate)
    bad2 = s**2 > p["relerr_max"] * res["tot_var"]
    return np.where(bad2, np.nan, f), np.where(bad2, np.nan, s)


def block_holdout(flat, frac=0.1, di=8, dj=32, seed=0):
    """Withhold random spatial blocks (2 deg lat x 8 deg lon) of profiles for cross-validation."""
    i, j = flat // NX, flat % NX
    bid = (i // di) * 1000 + (j // dj)
    rng = np.random.default_rng(seed)
    ub = np.unique(bid)
    chosen = rng.choice(ub, size=max(1, int(frac * len(ub))), replace=False)
    return np.isin(bid, chosen)


def profile_nodes(lat, lon):
    i = np.clip(np.rint((lat - 60.0) / 0.25).astype(int), 0, NY - 1)
    j = np.rint((lon + 180.0) / 0.25).astype(int) % NX
    return i * NX + j
