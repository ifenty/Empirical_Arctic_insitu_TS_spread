"""Focused checks for report statistics and regeneration-safe summaries."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr


def load_script(filename):
    path = Path(__file__).resolve().parents[1] / "scripts" / filename
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


spread = load_script("05_report_expected_spread.py")
compare = load_script("06_report_aotd_vs_woa23.py")


def test_spread_table_discloses_each_finite_mask(monkeypatch):
    fields = {"woa23": [1, 3, np.nan], "aotd": [np.nan, 2, np.nan],
              "woa23_and_aotd": [1, 3, 4]}

    def fake_load(kind, product, variable):
        return xr.Dataset({f"{variable}_{kind}": (("depth", "y", "x"), [[fields[product]]])},
                          coords={"depth": [10]})

    monkeypatch.setattr(spread, "load", fake_load)
    monkeypatch.setattr(spread, "DEPTHS", (10,))
    text = spread.table("rmse", "Theta")
    assert "valid cells: WOA23 / AOTD / combined" in text
    assert "| 10 | 2 [1.2–2.8] | 2 [2–2] | 3 [1.4–3.8] | 2 / 1 / 3 |" in text


def test_sparse_percentages_condition_on_positive_counts(monkeypatch):
    ds = xr.Dataset({"Theta_n": (("depth", "y", "x"), [[[0, 1, 4, 30]], [[0, 0, 0, 0]]])},
                    coords={"depth": [10, 50]})
    monkeypatch.setattr(spread, "load", lambda *args: ds.copy(deep=True))
    monkeypatch.setattr(spread, "DEPTHS", (10, 50))
    text = spread.counts_table()
    assert "| 10 | 35 | 3 | 66.7 | 66.7 |" in text
    assert "| 50 | 0 | 0 | – | – |" in text


def test_reference_statistics_use_finite_intersection():
    woa = np.zeros((12, 2, 1, 1, 12))
    woa[:, 0, 0, 0] = np.arange(12)
    woa[:, 1] = 29
    pub = woa.copy()
    pub[:, 0] += 2
    pub[:, 1] += 1
    woa[:, 0, 0, 0, -1] = np.nan
    df = compare.stats(woa, pub, np.array([10]))
    t = df[(df["var"] == "T") & (df.month == 1)].iloc[0]
    assert (t.n_aotd, t.n_woa, t.n_both) == (12, 11, 11)
    assert t.bias == pytest.approx(2)
    assert t.rmse == pytest.approx(2)
    assert t["corr"] == pytest.approx(1)
    assert df.loc[df["var"] == "S", "corr"].isna().all()


def test_insufficient_shared_cells_have_explicit_missing_statistics():
    fields = np.zeros((12, 2, 1, 1, 10))
    df = compare.stats(fields, fields, np.array([10]))
    assert (df.n_both == 10).all()
    assert df[["bias", "rmse", "corr", "maxabs"]].isna().all().all()
    assert "| 10 | – | – | – | – | – |" in compare.comparison_table(df)


def test_group_summary_is_median_of_monthly_metrics():
    rows = []
    for v in ("T", "S"):
        for month in range(1, 13):
            rows.append(dict(depth=10, var=v, month=month, rmse=month, bias=-month))
    text = compare.comparison_table(pd.DataFrame(rows))
    assert "| 10 | 6.00 | 6.50 | -6.50 | 6.00 | 6.50 |" in text


def test_mask_shift_diagnostic_keeps_longitude_periodic():
    rng = np.random.default_rng(42)
    mask = rng.random((9, 17)) > .4
    woa = np.where(mask, 1., np.nan)[None, None, None]
    woa = np.repeat(woa, 2, axis=1)
    pub = np.roll(woa, 2, axis=-1)
    alignment = compare.mask_diagnostics(woa, pub, np.array([10]), levels=(0,), max_shift=3)
    assert (alignment.interior_best_pct == 0).all()
    assert (alignment.row_offset == 0).all()
    assert (alignment.column_roll == 2).all()
    assert (alignment.interior_zero_pct > 0).all()
