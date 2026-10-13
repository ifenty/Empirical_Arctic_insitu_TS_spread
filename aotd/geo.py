"""Grids and spherical geometry helpers."""
from __future__ import annotations

import numpy as np
from pyproj import Transformer

EARTH_R_KM = 6371.0  # radius used by the MATLAB `distance(...,6371)` calls and by this package


def lonlat_to_xyz(lon, lat):
    lon = np.deg2rad(np.asarray(lon, dtype="float64"))
    lat = np.deg2rad(np.asarray(lat, dtype="float64"))
    c = np.cos(lat)
    return np.stack([c * np.cos(lon), c * np.sin(lon), np.sin(lat)], axis=-1)


def arc_km_to_chord(r_km):
    return 2.0 * np.sin(np.asarray(r_km) / (2.0 * EARTH_R_KM))


def chord_to_arc_km(chord):
    return 2.0 * EARTH_R_KM * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0))


# --------------------------------------------------------------------------- grids
def upstream_grid():
    """The grid of the published AOTD climatology: nodes on whole multiples of 0.25 deg.

    Returns lat (121,), lon (1440,) and 2-D (lat, lon) arrays.
    """
    lat = np.arange(121) * 0.25 + 60.0
    lon = np.arange(1440) * 0.25 - 180.0
    return lat, lon


def woa_grid(lat_min=60.0):
    """WOA23 0.25 deg grid (cell centres at x.125) restricted to lat >= lat_min."""
    lat = np.arange(-89.875, 90.0, 0.25)
    lat = lat[lat > lat_min]
    lon = np.arange(-179.875, 180.0, 0.25)
    return lat, lon


EASE2_N_CELL_M = 25000.0
EASE2_N_FULL = 720
EASE2_N_ORIGIN_M = 9_000_000.0  # upper-left corner of the full 720 x 720 grid (EASE-Grid 2.0 North, 25 km)


def ease2_north_25km(half_width_m: float = 3_500_000.0):
    """EASE-Grid 2.0 North (EPSG:6931) 25 km cell centres in a window of +-half_width_m.

    Returns x (n,), y (n,) [m; y decreasing like the NSIDC grid], lon (n,n), lat (n,n) of shape (y, x).
    The window +-3500 km contains the whole region north of 60N (max radius 3298 km).
    """
    n_half = int(round(half_width_m / EASE2_N_CELL_M))
    mid = EASE2_N_FULL // 2
    i = np.arange(mid - n_half, mid + n_half)
    x = -EASE2_N_ORIGIN_M + (i + 0.5) * EASE2_N_CELL_M
    y = EASE2_N_ORIGIN_M - (i + 0.5) * EASE2_N_CELL_M
    xx, yy = np.meshgrid(x, y)
    tr = Transformer.from_crs("EPSG:6931", "EPSG:4326", always_xy=True)
    lon, lat = tr.transform(xx, yy)
    return x, y, lon, lat


def lonlat_to_ease2(lon, lat):
    tr = Transformer.from_crs("EPSG:4326", "EPSG:6931", always_xy=True)
    return tr.transform(lon, lat)
