"""Concatenate the yearly AOTD scatter files into one cache, data/scatter_all.npz.

Missing values are normalised to NaN (the raw files mark missing salinity/pressure with 0.0 and leave arbitrary densities at missing
levels; see aotd/io.py). Keys: lat, lon, time (int64 seconds since 1970-01-01, i.e. datetime64[s] viewed as int), year, month, depth (102 WOA levels, m), and
temp, salt, pres, dens of shape (n_profiles, 102).

Usage (from the repo root): python scripts/02_prepare_scatter.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from aotd.io import load_scatter

d = load_scatter("data/scatter")
d["time"] = d["time"].astype("datetime64[s]").astype("int64")
np.savez("data/scatter_all.npz", **d)
print("profiles:", len(d["lat"]), " years:", d["year"].min(), "-", d["year"].max(), " levels:", len(d["depth"]))
