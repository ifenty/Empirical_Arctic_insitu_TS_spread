"""Loading of the published AOTD scatter (profile) files."""
from __future__ import annotations

import glob
import os

import numpy as np
import xarray as xr

T0 = np.datetime64("1980-01-01T00:00:00")  # the files carry the malformed unit "since 1980-01-01"


def days_to_datetime64(days: np.ndarray) -> np.ndarray:
    days = np.asarray(days, dtype="float64")
    if not np.isfinite(days).all():
        raise ValueError("Scatter profile times must be finite days since 1980-01-01")
    return T0 + np.round(days * 86400.0).astype("timedelta64[s]")


def normalize_missing_values(d: dict) -> dict:
    """In the published files only `temp` marks missing data with NaN; `salt` and `pres` use exactly 0.0 and `dens`
    holds arbitrary numbers (1003-1025) at missing levels. Make NaN the universal marker, in place."""
    valid = np.isfinite(d["temp"])
    d["temp"] = np.where(valid, d["temp"], np.nan).astype(d["temp"].dtype)
    d["salt"] = np.where(valid & np.isfinite(d["salt"]) & (d["salt"] != 0), d["salt"], np.nan).astype(d["salt"].dtype)
    d["pres"] = np.where(valid & np.isfinite(d["pres"]) & (d["pres"] != 0), d["pres"], np.nan).astype(d["pres"].dtype)
    d["dens"] = np.where(valid & np.isfinite(d["salt"]) & np.isfinite(d["dens"]), d["dens"], np.nan).astype(d["dens"].dtype)
    return d


def load_scatter(scatter_dir: str = "data/scatter", years=None, dtype="float32", normalize_missing=True) -> dict:
    """Concatenate yearly AOTD files into one dict of arrays.

    Returns dict with 1-D ``lat, lon, time`` (datetime64[s]) and ``year``, ``month``
    and 2-D ``temp, salt, pres, dens`` of shape (n_profiles, n_depth), plus ``depth``.
    The yearly files store profile index under the (misleading) dimension names lat/lon/time.
    """
    files = sorted(glob.glob(os.path.join(scatter_dir, "AOTD*.nc")))
    if years is not None:
        years = set(int(y) for y in years)
        files = [f for f in files if int(os.path.basename(f)[4:8]) in years]
    if not files:
        raise FileNotFoundError(f"No matching AOTD scatter files found in {scatter_dir}")
    out = {k: [] for k in ("lat", "lon", "time", "temp", "salt", "pres", "dens")}
    depth = None
    for f in files:
        with xr.open_dataset(f, decode_times=False) as d:
            file_depth = d["depth"].values.astype("float64")
            if depth is not None and not np.array_equal(file_depth, depth):
                raise ValueError(f"Inconsistent scatter depth levels in {f}")
            depth = file_depth
            for k in ("lat", "lon"):
                out[k].append(d[k].values.astype("float64"))
            out["time"].append(days_to_datetime64(d["time"].values))
            for k in ("temp", "salt", "pres", "dens"):
                out[k].append(d[k].values.astype(dtype))
    res = {k: np.concatenate(v) for k, v in out.items()}
    res["depth"] = depth
    if normalize_missing:
        normalize_missing_values(res)
    t = res["time"]
    res["year"] = t.astype("datetime64[Y]").astype(int) + 1970
    res["month"] = t.astype("datetime64[M]").astype(int) % 12 + 1
    return res
