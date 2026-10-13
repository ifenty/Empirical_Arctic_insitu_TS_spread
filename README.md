# AOTD_v2

Python port, reproduction and improved re-build of the Arctic Ocean Thermohaline Dataset
(AOTD; Sci. Data 2025, https://www.nature.com/articles/s41597-025-05855-3), compared against WOA23 (1991-2020).

Large data files (scatter, climatologies, WOA23, new fields) are **not** tracked; only code, reports and figures.

## Layout
- `aotd/` – Python package (io, geo, neighbors, faithful port, improved method)
- `upstream_matlab/` – original MATLAB for reference
- `docs/upstream_bug_log.md` – bugs found in the MATLAB (source for the note to the author)
- `reports/` – reports and figures
- `STATUS.md` – running progress log

Run with `PYTHONPATH=. python ...` from the repo root.
