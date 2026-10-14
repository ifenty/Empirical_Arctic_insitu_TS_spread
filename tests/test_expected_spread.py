import importlib

import numpy as np
import pytest
import xarray as xr

spread = importlib.import_module("03_compute_expected_spread")


def test_pole_is_included_in_last_reference_cell():
    cell, row, col, inside = spread.profile_grid_indices(
        np.array([60, 89.99, 90, 59.9]), np.array([-180, 180, 0, 0]),
    )
    np.testing.assert_array_equal(inside, [True, True, True, False])
    np.testing.assert_array_equal(row[:3], [1, 120, 120])
    np.testing.assert_array_equal(col[:3], [0, 0, 720])
    assert np.all((cell[inside] >= 0) & (cell[inside] < spread.NC ** 2))


@pytest.mark.parametrize("lat,lon", [(np.nan, 0), (80, np.inf), (91, 0), (80, 181)])
def test_invalid_geographic_coordinates_fail_explicitly(lat, lon):
    with pytest.raises(ValueError, match="Profile latitude"):
        spread.profile_grid_indices(np.array([lat]), np.array([lon]))


def test_depth_validation_checks_values_not_only_length():
    spread.check_depth_levels(np.array([0, 5, 10]), np.array([0, 5, 10]), "test")
    with pytest.raises(ValueError, match="Depth levels in test"):
        spread.check_depth_levels(np.array([0, 5, 15]), np.array([0, 5, 10]), "test")


def test_sample_statistics_and_empty_or_singleton_cells():
    result = spread.deviation_statistics(
        np.array([0, 1, 2]), np.array([0, 3, 4]), np.array([0, 9, 10]),
        n_eff=np.array([0, 1, 2]), n_years=np.array([0, 1, 2]),
    )
    np.testing.assert_allclose(result["bias"], [np.nan, 3, 2], equal_nan=True)
    # Sample std (ddof=1): NaN for n=0 and n=1; for n=2 with values [1, 3], mean=2, sample_std=sqrt(2)
    np.testing.assert_allclose(result["std"], [np.nan, np.nan, np.sqrt(2)], equal_nan=True)
    np.testing.assert_allclose(result["rmse"], [np.nan, 3, np.sqrt(5)], equal_nan=True)
    np.testing.assert_array_equal(result["n_eff"], [0, 1, 2])
    np.testing.assert_array_equal(result["n_years"], [0, 1, 2])


def test_envelope_selects_each_statistics_own_bias_and_count():
    woa = dict(std=np.array([2, np.nan, 1, np.nan, 1]),
               rmse=np.array([3, np.nan, 2, np.nan, 2]),
               bias=np.array([10, np.nan, 30, np.nan, 50]),
               n=np.array([10, 0, 30, 0, 50]))
    aotd = dict(std=np.array([1, 3, np.nan, np.nan, 1]),
                rmse=np.array([4, 4, np.nan, np.nan, 2]),
                bias=np.array([20, 40, np.nan, np.nan, 60]),
                n=np.array([20, 40, 0, 0, 60]))
    result = spread.combine_statistics(woa, aotd)
    np.testing.assert_allclose(result["std"], [2, 3, 1, np.nan, 1], equal_nan=True)
    np.testing.assert_allclose(result["rmse"], [4, 4, 2, np.nan, 2], equal_nan=True)
    np.testing.assert_allclose(result["bias_std"], [10, 40, 30, np.nan, 50], equal_nan=True)
    np.testing.assert_allclose(result["bias_rmse"], [20, 40, 30, np.nan, 50], equal_nan=True)
    np.testing.assert_array_equal(result["n_std"], [10, 40, 30, 0, 50])
    np.testing.assert_array_equal(result["n_rmse"], [20, 40, 30, 0, 50])


def test_envelope_ties_are_resolved_at_stored_precision():
    woa = dict(std=np.array([1.0]), rmse=np.array([2.0]), bias=np.array([0.0]), n=np.array([10]))
    aotd = dict(std=np.array([1.0 + 1e-9]), rmse=np.array([2.0 + 1e-9]), bias=np.array([1.0]), n=np.array([20]))
    result = spread.combine_statistics(woa, aotd)
    np.testing.assert_array_equal(result["n_std"], woa["n"])
    np.testing.assert_array_equal(result["n_rmse"], woa["n"])


def test_main_writes_complete_products_from_synthetic_profiles(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "data").mkdir()
    depth = np.array([0, 100, 2000])
    temperature = np.array([[2, 2, 2], [4, 4, 4]], dtype="float32")
    salinity = np.full_like(temperature, 35)
    np.savez(
        "data/scatter_all.npz", lat=[80, 80], lon=[0, 0], month=[1, 1],
        year=[1992, 2023], depth=depth, temp=temperature, salt=salinity,
    )
    monkeypatch.setattr(spread, "NC", 2)
    monkeypatch.setattr(spread, "NMONTHLY", 2)
    monkeypatch.setattr(
        spread, "profile_grid_indices",
        lambda lat, lon: (np.array([0, 0]), np.array([0, 0]), np.array([0, 0]), np.array([True, True])),
    )

    def reference(letter, var):
        value = 1 if letter == "t" else 35
        return np.full((12, 2, 1, 1), value), np.full((4, 1, 1, 1), value), depth

    monkeypatch.setattr(spread, "load_ref", reference)
    monkeypatch.setattr(spread, "load_aotd", lambda key, depth: np.full((12, 2, 1, 1), 2 if key == "temp" else 34))
    monkeypatch.setattr(
        spread, "spatial_interp_weights",
        lambda lat, lon: (np.zeros(len(lat), int), np.zeros(len(lat), int), np.zeros(len(lat)),
                          np.zeros(len(lat), int), np.zeros(len(lat), int), np.zeros(len(lat))),
    )
    spread.main()
    files = list((tmp_path / "outputs" / "expected_spread").glob("*.nc"))
    assert len(files) == 12
    for path in files:
        with xr.open_dataset(path) as ds:
            assert dict(ds.sizes) == {"depth": 3, "y": 2, "x": 2}
            short = "Theta" if "potential_temperature" in path.name else "S"
            metric = "std" if "expected_std_" in path.name else "rmse"
            count, values = ds[f"{short}_n"].values, ds[f"{short}_{metric}"].values
            assert f"{short}_n_eff" in ds and f"{short}_n_years" in ds
            np.testing.assert_array_equal(np.isnan(values), count < (2 if metric == "std" else 1))
            assert (count[:, 1:, :] == 0).all()
            if "from_aotd_" in path.name:
                assert (count[2] == 0).all()
            else:
                np.testing.assert_array_equal(count[:, 0, 0], [2, 2, 2])
                np.testing.assert_array_equal(ds[f"{short}_n_eff"].values[:, 0, 0], [2, 2, 2])
                np.testing.assert_array_equal(ds[f"{short}_n_years"].values[:, 0, 0], [2, 2, 2])
            if short == "Theta" and "from_woa23_" in path.name:
                assert values[0, 0, 0] == pytest.approx(np.sqrt(2) if metric == "std" else np.sqrt(5))
