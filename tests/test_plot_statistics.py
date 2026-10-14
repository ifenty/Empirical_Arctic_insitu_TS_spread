import importlib

import numpy as np

plots = importlib.import_module("04_make_plots")


def test_percentiles_preserve_empty_levels_without_warnings():
    values = np.array([[np.nan, np.nan], [1, 3]])
    with np.errstate(all="raise"):
        result = plots.cell_percentiles(values, (10, 50, 90))
    np.testing.assert_allclose(result, [[np.nan, 1.2], [np.nan, 2], [np.nan, 2.8]], equal_nan=True)


def test_all_levels_can_be_empty():
    result = plots.cell_percentiles(np.full((3, 2, 2), np.nan), (10, 50, 90))
    assert result.shape == (3, 3)
    assert np.isnan(result).all()


def test_ratio_requires_both_counts_and_nonzero_standard_deviation():
    with np.errstate(all="raise"):
        result = plots.rmse_std_ratio(
            np.array([2, 2, 2, 2]), np.array([1, 0, 1, 1]),
            np.array([5, 5, 5, 4]), np.array([5, 5, 4, 5]),
        )
    np.testing.assert_allclose(result, [2, np.nan, np.nan, np.nan], equal_nan=True)


def test_main_creates_figure_directory_before_rendering(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    seen = []

    def render(*args):
        assert (tmp_path / plots.OUTD).is_dir()
        seen.append(args)

    monkeypatch.setattr(plots, "make", render)
    for name in ("vs_depth", "bias_maps", "product_difference_maps", "counts_vs_depth",
                 "declustered_vs_depth", "declustered_maps", "method_comparison_maps",
                 "variance_components_maps", "method_distributions_by_depth", "prior_comparison_maps"):
        monkeypatch.setattr(plots, name, render)
    plots.main()
    assert len(seen) == 21
