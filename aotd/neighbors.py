"""Chunked spherical fixed-radius neighbour search (replaces MATLAB's per-cell `distance()` loops)."""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from .geo import arc_km_to_chord, chord_to_arc_km, lonlat_to_xyz


class SphereNeighbors:
    """Pairs (observation, target) closer than R km on a sphere of radius 6371 km.

    Targets are processed in chunks so memory stays bounded for 555 km radii at the surface.
    """

    def __init__(self, obs_lon, obs_lat, tgt_lon, tgt_lat, chunk: int = 1000):
        self.obs_xyz = lonlat_to_xyz(obs_lon, obs_lat)
        self.tgt_xyz = lonlat_to_xyz(np.ravel(tgt_lon), np.ravel(tgt_lat))
        self.tree = cKDTree(self.obs_xyz) if len(self.obs_xyz) else None
        self.chunk = chunk
        self.n_tgt = len(self.tgt_xyz)

    def pairs(self, R_km: float):
        """Yield (t0, t1, i_obs, j_tgt_local, r_km) per target chunk with r < R_km."""
        if self.tree is None:
            return
        maxd = float(arc_km_to_chord(R_km))
        for t0 in range(0, self.n_tgt, self.chunk):
            t1 = min(t0 + self.chunk, self.n_tgt)
            sub = cKDTree(self.tgt_xyz[t0:t1])
            m = self.tree.sparse_distance_matrix(sub, maxd, output_type="ndarray")
            if len(m) == 0:
                continue
            r = chord_to_arc_km(m["v"])
            keep = r < R_km
            yield t0, t1, m["i"][keep], m["j"][keep], r[keep]
