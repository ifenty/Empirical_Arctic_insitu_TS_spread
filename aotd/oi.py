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
    relerr_max=5.0,            # NaN if calibrated std exceeds this fraction of the local natural variability (loc_sd)
    Kw=50, Rw=500.0, n0=10, floor=0.02,  # moving window for local variances: nearest Kw nodes within Rw km; shrinkage
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


def local_var(pos, resid, nsr, query, glob, p=PARAMS):
    """Moving-window prior variance: mean squared residual of the Kw nearest nodes within Rw km of each query point,
    deflated for observation noise and shrunk toward the level-wide value `glob` when few nodes are in the window."""
    tree = cKDTree(pos)
    Kw = min(p["Kw"], len(pos))
    d, nb = tree.query(query, k=Kw, distance_upper_bound=_chord(p["Rw"]))
    v = np.isfinite(d)
    nb = np.where(v, nb, 0)
    cnt = v.sum(1)
    loc = (np.where(v, resid[nb] ** 2, 0.0).sum(1) / np.maximum(cnt, 1)) / (1.0 + (np.where(v, nsr[nb], 0.0).sum(1) / np.maximum(cnt, 1)))
    out = (cnt * loc + p["n0"] * glob) / (cnt + p["n0"])
    return np.maximum(out, p["floor"] * glob)


def oi_stage(idx, resid, n, targets, L, p=PARAMS, chunk=4000):
    """Local OI of node residuals with a moving-window prior variance. Prior mean of the residual is 0.

    Returns increment, posterior variance, prior variance at the targets, nearest-node distance km.
    """
    nt = len(targets)
    if len(idx) < 3:
        return np.zeros(nt), np.full(nt, np.nan), np.full(nt, np.nan), np.full(nt, np.inf)
    nsr = p["nsr_a"] + p["nsr_b"] / n
    glob = max(np.var(resid) / (1.0 + nsr.mean()), 1e-12)
    pos = TARGET_XYZ[idx]
    s2_node = local_var(pos, resid, nsr, pos, glob, p)
    s2_t = local_var(pos, resid, nsr, TARGET_XYZ[targets], glob, p)
    sn, st = np.sqrt(s2_node), np.sqrt(s2_t)
    tree = cKDTree(pos)
    K = min(p["K"], len(idx))
    d, nb = tree.query(TARGET_XYZ[targets], k=K, distance_upper_bound=_chord(p["rmax_fac"] * L))
    miss = ~np.isfinite(d)
    nb = np.where(miss, 0, nb)
    dkm = 2.0 * R_EARTH * np.arcsin(np.minimum(np.where(miss, 2.0, d), 2.0) / 2.0)
    dmin = np.where(miss[:, 0], np.inf, dkm[:, 0])
    inc = np.zeros(nt)
    post = s2_t.copy()
    eye = np.eye(K)[None]
    for a in range(0, nt, chunk):
        b = min(a + chunk, nt)
        nbc, ms = nb[a:b], miss[a:b]
        P = pos[nbc]
        ch = np.linalg.norm(P[:, :, None, :] - P[:, None, :, :], axis=-1)
        rr = 2.0 * R_EARTH * np.arcsin(np.minimum(ch, 2.0) / 2.0)
        s = sn[nbc]
        C = s[:, :, None] * s[:, None, :] * np.exp(-(rr / L) ** 2)
        c = st[a:b, None] * s * np.exp(-(dkm[a:b] / L) ** 2)
        c = np.where(ms, 0.0, c)
        C = np.where(ms[:, :, None] | ms[:, None, :], 0.0, C)
        C = C + (s**2 * nsr[nbc] + np.where(ms, 1.0, 0.0) * s2_t[a:b, None])[:, :, None] * eye
        w = np.linalg.solve(C, c[:, :, None])[:, :, 0]
        inc[a:b] = (w * np.where(ms, 0.0, resid[nbc])).sum(1)
        post[a:b] = s2_t[a:b] - (w * c).sum(1)
    return inc, np.maximum(post, 0.0), s2_t, dmin


def local_sd(idx, m, targets, p=PARAMS):
    """Local natural variability: std of node means in the moving window (all months/years pooled)."""
    pos = TARGET_XYZ[idx]
    tree = cKDTree(pos)
    Kw = min(p["Kw"], len(idx))
    d, nb = tree.query(TARGET_XYZ[targets], k=Kw, distance_upper_bound=_chord(p["Rw"]))
    v = np.isfinite(d)
    nb = np.where(v, nb, 0)
    cnt = v.sum(1)
    vals = np.where(v, m[nb], np.nan)
    sd = np.sqrt(np.nanvar(vals, axis=1))
    return np.where(cnt >= 5, sd, np.nan)


def run_variable(flat, val, month, ok, exclude=None, p=PARAMS):
    """Staged OI for one variable at one level.

    flat: node index of each profile; val: values; month: 1..12; ok: (NY*NX,) bool valid ocean cells at this level;
    exclude: optional bool per profile to withhold (cross-validation).
    Returns dict(field (12,NY,NX), std (12,NY,NX), dmin (NY,NX), loc_sd (NY,NX)).
    """
    use = np.ones(len(val), bool) if exclude is None else ~exclude
    tgt = np.nonzero(ok)[0]
    idx, m, n = node_table(flat, val, use)
    mu = m.mean()
    inc, P_a, s2_a, dmin = oi_stage(idx, m - mu, n, tgt, p["L"][0], p)
    lsd = local_sd(idx, m, tgt, p)
    A = np.full(NY * NX, np.nan); A[tgt] = mu + inc
    PA = np.full(NY * NX, np.nan); PA[tgt] = P_a
    field = np.full((12, NY * NX), np.nan)
    err2 = np.full((12, NY * NX), np.nan)
    seas = {}
    for s in range(1, 5):
        sel = np.isin(month, np.nonzero(SEASON_OF_MONTH == s)[0] + 1) & use
        idx, m, n = node_table(flat, val, sel)
        keep = ok[idx]
        idx, m, n = idx[keep], m[keep], n[keep]
        inc, P_s, s2_s, _ = oi_stage(idx, m - A[idx], n, tgt, p["L"][1], p)
        S = np.full(NY * NX, np.nan); S[tgt] = A[tgt] + inc
        Pst = np.full(NY * NX, np.nan)
        if len(idx) >= 3:  # error of the previous stage is carried in proportion to how little the data reduced the prior
            Pst[tgt] = P_s + (P_s / s2_s) * PA[tgt]
        else:
            Pst[tgt] = PA[tgt]
        seas[s] = (S, Pst)
    s2_months, nodata = [], []
    for mo in range(1, 13):
        S, Pst = seas[SEASON_OF_MONTH[mo - 1]]
        sel = (month == mo) & use
        idx, m, n = node_table(flat, val, sel)
        keep = ok[idx]
        idx, m, n = idx[keep], m[keep], n[keep]
        inc, P_m, s2_m, _ = oi_stage(idx, m - S[idx], n, tgt, p["L"][2], p)
        if len(idx) >= 3:
            field[mo - 1, tgt] = S[tgt] + inc
            err2[mo - 1, tgt] = P_m + (P_m / s2_m) * Pst[tgt]
            s2_months.append(s2_m)
        else:  # no monthly information: seasonal field with the typical local monthly anomaly variance added
            field[mo - 1, tgt] = S[tgt]
            err2[mo - 1, tgt] = Pst[tgt]
            nodata.append(mo - 1)
    for mi in nodata:
        err2[mi, tgt] += np.nanmedian(s2_months, axis=0) if s2_months else 0.0
    return dict(field=field.reshape(12, NY, NX), std=np.sqrt(err2).reshape(12, NY, NX),
                dmin=_scatter(dmin, tgt), loc_sd=_scatter(lsd, tgt))


def _scatter(v, tgt):
    out = np.full(NY * NX, np.nan)
    out[tgt] = v
    return out.reshape(NY, NX)


def apply_mask(res, inflate=1.0, p=PARAMS):
    """NaN where the nearest observation node is beyond dmax, or the calibrated std exceeds relerr_max times the local
    natural variability. `inflate` is the empirical cross-validation calibration of the formal std."""
    bad = ~(res["dmin"] <= p["dmax"])
    s = res["std"] * inflate
    bad2 = bad[None] | ~(s <= p["relerr_max"] * res["loc_sd"][None])
    return np.where(bad2, np.nan, res["field"]), np.where(bad2, np.nan, s)


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
