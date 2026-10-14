"""Expected standard deviation of in situ T and S about WOA23, on the EASE2-Grid 2.0 North 100 km grid.

For every in situ profile value (scatter data, 1992 onwards) the deviation d = obs - WOA23 reference is formed at the profile's
location (nearest WOA23 0.25 deg cell). The reference is the monthly WOA23 field of the profile's month for the 57 monthly levels
(0-1500 m) and the seasonal field of the profile's season for deeper levels. Deviations are pooled per 100 km cell and standard
depth level:

    T_std = sqrt(mean(d^2))        (spread of the observations about the reference, bias included)
    T_bias = mean(d), T_n = number of values

Cells/levels with no observation are NaN (no extrapolation). No time dimension.

Usage: PYTHONPATH=. python scripts/expected_std.py
Outputs (untracked): outputs/expected_std/expected_std_{temperature,salinity}_ease2_100km.nc
"""
import os
import numpy as np
import xarray as xr
from pyproj import Transformer

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


def main():
    d = np.load("data/scatter_all.npz")
    lat, lon, month = d["lat"].astype("float64"), d["lon"].astype("float64"), d["month"].astype(int)
    nlev = d["depth"].shape[0]
    t = Transformer.from_crs(4326, 6931, always_xy=True)
    x, y = t.transform(lon, lat)
    ix = np.floor((x + HALF) / CELL).astype(int)
    iy = np.floor((HALF - y) / CELL).astype(int)
    year = d["year"].astype(int)
    inside = (year >= YEAR0) & (ix >= 0) & (ix < NC) & (iy >= 0) & (iy < NC)
    cell = iy * NC + ix
    row = np.floor((lat + 90.0) / 0.25).astype(int) - R0
    col = np.floor((lon + 180.0) / 0.25).astype(int) % 1440
    inside &= (row >= 0) & (row < NROW)
    print("profiles", len(lat), "inside grid/domain", int(inside.sum()), flush=True)

    result = {}
    for letter, var, key, name in (("t", "t_an", "temp", "temperature"), ("s", "s_an", "salt", "salinity")):
        mon, ses, depth = load_ref(letter, var)
        assert nlev == len(depth)
        obs = d[key]
        n = np.zeros((nlev, NC * NC)); s1 = np.zeros_like(n); s2 = np.zeros_like(n)
        for m in range(1, 13):
            sel = np.nonzero((month == m) & inside)[0]
            r, c, ce = row[sel], col[sel], cell[sel]
            ref = np.empty((nlev, len(sel)), "float32")
            ref[:NMONTHLY] = mon[m - 1][:, r, c]
            ref[NMONTHLY:] = ses[SEASON_FILE[m - 1] - 13][:, r, c]
            dev = obs[sel].T.astype("float64") - ref  # (nlev, nsel)
            good = np.isfinite(dev)
            for k in range(nlev):
                g = good[k]
                if not g.any():
                    continue
                n[k] += np.bincount(ce[g], minlength=NC * NC)
                s1[k] += np.bincount(ce[g], weights=dev[k, g], minlength=NC * NC)
                s2[k] += np.bincount(ce[g], weights=dev[k, g] ** 2, minlength=NC * NC)
            print(letter, "month", m, "done", flush=True)
        with np.errstate(invalid="ignore", divide="ignore"):
            std = np.where(n > 0, np.sqrt(s2 / n), np.nan)
            bias = np.where(n > 0, s1 / n, np.nan)
        result[name] = (std.reshape(nlev, NC, NC), bias.reshape(nlev, NC, NC), n.reshape(nlev, NC, NC), depth)

    xc = -HALF + CELL * (np.arange(NC) + 0.5)
    yc = HALF - CELL * (np.arange(NC) + 0.5)
    X, Y = np.meshgrid(xc, yc)
    glon, glat = Transformer.from_crs(6931, 4326, always_xy=True).transform(X, Y)
    for name, (std, bias, n, depth) in result.items():
        short = "T" if name == "temperature" else "S"
        unit = "degree_Celsius" if name == "temperature" else "1e-3"  # practical salinity (WOA23 s_an)
        sname = "sea_water_temperature" if name == "temperature" else "sea_water_salinity"
        ds = xr.Dataset(
            {
                f"{short}_std": (("depth", "y", "x"), std.astype("float32"),
                                 dict(long_name=f"expected standard deviation of in situ {name} about WOA23",
                                      units=unit, grid_mapping="crs", coordinates="lat lon",
                                      comment="sqrt(mean((obs - WOA23)^2)) over all profile values in the cell, 1992-2023, no time dimension. "
                                              "Reference: WOA23 decav91C0 monthly (0-1500 m) or seasonal (deeper) field of the profile month. NaN where no observation.")),
                f"{short}_bias": (("depth", "y", "x"), bias.astype("float32"),
                                  dict(long_name=f"mean in situ {name} minus WOA23", units=unit, grid_mapping="crs", coordinates="lat lon")),
                f"{short}_n": (("depth", "y", "x"), n.astype("int32"),
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
            attrs=dict(title=f"Expected standard deviation of in situ {name} relative to WOA23 on the EASE2 100 km Arctic grid",
                       Conventions="CF-1.8", source="AOTD scatter profiles 1992-2023 interpolated to WOA standard depths; WOA23 decav91C0",
                       history="scripts/expected_std.py"))
        enc = {v: dict(zlib=True, complevel=4) for v in (f"{short}_std", f"{short}_bias", f"{short}_n")}
        f = f"{OUT}/expected_std_{name}_ease2_100km.nc"
        ds.to_netcdf(f, encoding=enc)
        print("wrote", f, flush=True)


if __name__ == "__main__":
    main()
