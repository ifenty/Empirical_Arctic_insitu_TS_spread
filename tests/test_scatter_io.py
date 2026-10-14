import numpy as np
import pytest
import xarray as xr

from scatter_io import days_to_datetime64, load_scatter, normalize_missing_values


def test_time_conversion_rounds_to_seconds():
    np.testing.assert_array_equal(
        days_to_datetime64(np.array([0, 1, 1.5])),
        np.array(["1980-01-01", "1980-01-02", "1980-01-02T12:00:00"], dtype="datetime64[s]"),
    )


def test_invalid_times_raise():
    with pytest.raises(ValueError, match="times must be finite"):
        days_to_datetime64(np.array([np.nan]))


def test_normalization_preserves_zero_temperature_and_masks_missing_values():
    data = {
        "temp": np.array([0, np.nan, 1, np.inf, 2], dtype="float32"),
        "salt": np.array([35, 35, 0, 35, np.inf], dtype="float32"),
        "pres": np.array([10, 10, 0, 10, 20], dtype="float32"),
        "dens": np.array([1025, 1025, 1025, 1025, 1025], dtype="float32"),
    }
    normalize_missing_values(data)
    np.testing.assert_allclose(data["temp"], [0, np.nan, 1, np.nan, 2], equal_nan=True)
    np.testing.assert_allclose(data["salt"], [35, np.nan, np.nan, np.nan, np.nan], equal_nan=True)
    np.testing.assert_allclose(data["pres"], [10, np.nan, np.nan, np.nan, 20], equal_nan=True)
    np.testing.assert_allclose(data["dens"], [1025, np.nan, np.nan, np.nan, np.nan], equal_nan=True)
    assert all(a.dtype == np.dtype("float32") for a in data.values())


def write_scatter(path, depth, day):
    values = np.ones((1, len(depth)), dtype="float32")
    xr.Dataset(
        {key: (("profile", "depth"), values) for key in ("temp", "salt", "pres", "dens")},
        coords={"lat": ("profile", [80]), "lon": ("profile", [10]), "time": ("profile", [day]), "depth": depth},
    ).to_netcdf(path)


def test_loading_valid_yearly_files(tmp_path):
    write_scatter(tmp_path / "AOTD1983.nc", [0, 5], 1096)
    write_scatter(tmp_path / "AOTD1984.nc", [0, 5], 1461)
    data = load_scatter(str(tmp_path))
    assert data["temp"].shape == (2, 2)
    np.testing.assert_array_equal(data["year"], [1983, 1984])
    np.testing.assert_array_equal(data["month"], [1, 1])
    assert load_scatter(str(tmp_path), years=[1983])["temp"].shape == (1, 2)


def test_loading_empty_directory_has_clear_error(tmp_path):
    with pytest.raises(FileNotFoundError, match="No matching AOTD scatter files"):
        load_scatter(str(tmp_path))


def test_mismatched_depths_fail_before_concatenation(tmp_path):
    write_scatter(tmp_path / "AOTD1983.nc", [0, 5], 1096)
    write_scatter(tmp_path / "AOTD1984.nc", [0, 10], 1461)
    with pytest.raises(ValueError, match="Inconsistent scatter depth levels"):
        load_scatter(str(tmp_path))
