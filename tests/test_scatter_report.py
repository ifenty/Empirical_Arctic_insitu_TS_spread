"""Regression checks for the standalone report without loading the full dataset."""
import ast
from pathlib import Path
import re
import warnings

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "reports" / "scatter"


@pytest.fixture
def helpers():
    # Executing only helper definitions avoids this diagnostic's data/plotting work.
    source = ast.parse((ROOT / "scripts" / "07_report_scatter_data.py").read_text())
    names = {"fmtv", "mdtable", "ndistinct", "ocean_grid_coverage", "level_mean_std"}
    module = ast.Module(
        body=[node for node in source.body if isinstance(node, ast.FunctionDef) and node.name in names],
        type_ignores=[],
    )
    namespace = {"np": np, "pd": pd}
    exec(compile(module, "scatter_report_helpers", "exec"), namespace)
    return namespace


def test_distinct_year_months_stay_in_their_original_boxes(helpers):
    helpers.update(cell=np.array([0, 0, 1, 1, 2]), NC=4)
    values = np.array([24001, 24001, 24001, 24002, 24001])
    np.testing.assert_array_equal(helpers["ndistinct"](values), [1, 2, 1, 0])
    np.testing.assert_array_equal(
        helpers["ndistinct"](values, np.array([False, False, True, True, False])),
        [0, 2, 0, 0],
    )
    np.testing.assert_array_equal(helpers["ndistinct"](values, np.zeros(5, bool)), [0, 0, 0, 0])


def test_ocean_coverage_excludes_land_and_unobserved_cells(helpers):
    counts = np.array([[10, 30], [0, 100]])
    ocean = np.array([[True, False], [True, False]])
    assert helpers["ocean_grid_coverage"](counts, ocean) == 0.5
    assert helpers["ocean_grid_coverage"](counts, ocean, threshold=10) == 0.5
    assert helpers["ocean_grid_coverage"](counts, ocean, threshold=11) == 0.0


def test_markdown_escapes_threshold_pipes_and_retains_integer_counts(helpers):
    table = pd.DataFrame({"|dT/dz|": ["> 0.3\nscreen"], "count": [2], "fraction": [0.5]})
    rendered = helpers["mdtable"](table, dec=2)
    assert "\\|dT/dz\\|" in rendered
    assert "> 0.3<br>screen" in rendered
    assert "| 2 | 0.50 |" in rendered
    assert helpers["fmtv"](np.float32(np.nan), 2) == "–"
    assert helpers["fmtv"](pd.NA, 2) == "–"
    assert "| 1984 |" in helpers["mdtable"](pd.DataFrame({"first year": [1984]}))


def test_empty_standard_levels_have_nan_moments_without_warnings(helpers):
    values = np.array([[np.nan, 1.0, 5.0], [np.nan, 3.0, np.nan]])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        mean, std = helpers["level_mean_std"](values)
    np.testing.assert_allclose(mean, [np.nan, 2.0, 5.0], equal_nan=True)
    np.testing.assert_allclose(std, [np.nan, 1.0, 0.0], equal_nan=True)


def test_level_moments_preserve_float32_statistics(helpers):
    values = np.arange(24, dtype="float32").reshape(6, 4)
    values[::2, 1] = np.nan
    mean, std = helpers["level_mean_std"](values)
    np.testing.assert_array_equal(mean, np.nanmean(values, axis=0))
    np.testing.assert_array_equal(std, np.nanstd(values, axis=0))


def test_rendered_report_links_and_table_structure():
    text = (REPORT / "scatter_report.md").read_text()
    headings = {
        re.sub(r"[^\w -]", "", line.lstrip("# ").lower()).replace(" ", "-")
        for line in text.splitlines() if line.startswith("#")
    }
    for target in re.findall(r"\]\(([^)]+)\)", text):
        if target.startswith("#"):
            assert target[1:] in headings, target
        elif not target.startswith(("http://", "https://")):
            assert (REPORT / target).is_file(), target
    rows_in_table = []
    for line in text.splitlines() + [""]:
        if line.startswith("|"):
            rows_in_table.append(line)
        elif rows_in_table:
            delimiters = [len(re.findall(r"(?<!\\)\|", row)) for row in rows_in_table]
            assert len(set(delimiters)) == 1, rows_in_table
            rows_in_table = []
    table_links = set(re.findall(r"\]\((tables/[^)]+)\)", text))
    assert table_links == {f"tables/{path.name}" for path in (REPORT / "tables").glob("*.csv")}


def test_static_screen_table_reports_its_actual_denominator():
    table = pd.read_csv(REPORT / "tables" / "s7_static_instability_per_class.csv")
    expected = 100 * table["profiles flagged (>0.03 kg/m3)"] / table["profiles with assessable pairs"]
    np.testing.assert_allclose(table["% of assessable profiles flagged"], expected)
    assert (table["profiles with assessable pairs"] <= table["profiles with S"]).all()
