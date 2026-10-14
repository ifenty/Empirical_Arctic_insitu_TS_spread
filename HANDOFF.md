# Handoff (2026-10-13)

## Goal
Plan in `STATUS.md`: 1 faithful port + reproduction, 2 QC port, 3 improved method, 4 1992-2023 dataset (WOA grid + EASE2 25 km), 5 WOA23 comparison, 6 scatter report, 7 note to author.

## Done
- Data downloaded, grids, neighbour search, faithful port (`aotd/`).
- Scatter report: `reports/scatter/`. Bug log: `docs/upstream_bug_log.md` (30 issues).
- Port matches the published valid-cell counts exactly. Monthly T RMSE vs published: 0.6 C (50 m), 0.4 C (100 m), 0.13 C (700 m). Residual is being investigated (data/QC version).

## Incident
- At about 09:31 the WSL VM restarted (journald "uncleanly shut down") and killed the full-reproduction run (10 workers, `scripts/reproduce.py`).
- Nothing was written: no `outputs/repro/month_*.npy` and no `reports/reproduction/level_stats.csv`. Only the `outputs/repro/pub_month.npy` cache (1.9 GB) exists.
- Root cause: memory in `SphereNeighbors` (`aotd/neighbors.py`). With `chunk=8000` targets, shallow levels with dense observations built huge pair arrays (more than 3.5 GB per process). 10 workers exhausted the shared VM.
- Fix: `chunk=1000`. Results are unchanged. Peak is about 1.1 GB per process.

## Current run
- Launched as systemd unit `aotd-repro5`: 2 workers, `MemoryMax=5600M`, about 1 GB anon total. Check with `systemctl --user status aotd-repro5`; log at `outputs/repro/log.txt`.
- It resumes safely and skips levels that already have `month_{k}.npy`. The script only writes `reports/reproduction/level_stats.csv` when the run finishes. Levels take about 2 min each per worker, deeper levels faster.
- If the unit died, relaunch it with the command in `.github/copilot-instructions.md`.

## Next steps
1. Wait for the run to finish. Verify `outputs/repro/month_k.npy` for k=1..56 and `reports/reproduction/level_stats.csv`.
2. Investigate the residual vs published (data/QC version), then the QC port. Update `STATUS.md`.

## Gotchas
- Use `.venv/bin/python`. The conda python lacks numpy.
- Never exceed 2 workers or 20% RAM. Another session shares this VM.
- The bathymetry file is missing upstream. `reproduce.py` builds the level mask from the published fields as a workaround.

## New requirement (2026-10-13): uncertainty + data-honest product
The climatology will constrain an ocean model, so it must not contain extrapolated "made-up" values.
- Provide a data-honest version: NaN wherever there is no real observation within a justified distance/support threshold. Do not fill by smoothing or background fields.
- Provide formal uncertainty estimates for every valid cell (T and S, per level and period). Candidates: observation-count and nearest-obs-distance support, sampling variability of the contributing obs, and an objective-analysis error estimate (Barnes/Cressman weights, or optimal-interpolation error with a noise-to-signal ratio).
- Choose the cutoff criteria with the user's constraints: no values where support is too distant or too sparse. Document them in `docs/` and the product metadata.
- Applies to the improved method (step 3) and the 1992-2023 dataset (step 4). Keep the faithful port unchanged.

## New requirement (2026-10-13): vertical grid
New climatologies must use the original WOA standard vertical grid (not a custom or reduced one). The scatter data already carries all 102 WOA levels (0-5500 m); the published AOTD uses only the top 57 (0-1500 m). Assumption: use the full 102-level WOA grid wherever data support it, with NaN below the data-supported depth. Confirm with the user whether the 0-1500 m subset is acceptable for the first release.

## Output format
New climatologies will be written as CF-compliant NetCDF4 (matching the published AOTD and WOA23), one file per month like the published product unless the user says otherwise. NaN / `_FillValue` where unsupported; uncertainty, observation count and nearest-observation distance as separate variables. No new climatology files exist yet; the faithful-port reproduction output is `.npy` in `outputs/repro/` (untracked).

## Reproduction diagnosis (2026-10-13)
- Port reproduces published March and September exactly (level 20: T 0.01/0.08 degC, S ~0.00-0.01 psu; monthly correction fields r=1.00).
- The other 10 published months lack the monthly obs signal (correction RMS ~0.15 degC vs 0.65 in port, uncorrelated). Not a month-label shift. Logged as BUG-29; confirm with the author. See `reports/reproduction/reproduction_report.md` (Diagnosis).
- Only level 20 was analysed in depth (`outputs/repro/lvl20_full.npz` holds its annual/seasonal fields). Next: repeat the Mar/Sep vs other-month split over all levels, and add per-month columns to the report. Validate the new product against Mar/Sep only.

### Hypotheses tested for the 10 non-Mar/Sep published months (level 20 T, January) - none reproduce them
- Month-label shift (obs from each of the 12 months): max correlation of correction with published = 0.37.
- Year subsets of January obs (eras, single years, <=2012/2020/2022, >=2005): max 0.21.
- Short-radius, min-count Barnes correction on the port seasonal first guess (dead-code `Barnes()` in `cressman_bkg.m` suggests ~80 km radius, >=10 obs): correlation rises to about 0.55 at best (R=50-65 km, min count 10-20) but never approaches 1.0. Grid: kernel exp(-c r^2/(0.2*8e4)), c in {1,4}, R in {50..120} km, min count in {5,10,20}.
- RMSE vs published is not a useful test here: "no correction" already gives about 0.18 degC.
Conclusion: unexplained from released code and data. Ask the author for the code/data used for those months. Scripts used were throwaway (in /tmp).

## Temporal treatment (open design item)
Upstream mapping has no temporal decorrelation: weights are purely spatial; `Lt = 1000` is defined in all three scripts but never used (BUG-26). All years 1983-2023 are pooled equally. The new product should handle time explicitly (time-dependent weights and/or trend removal before mapping), and the uncertainty should include interannual variability. Decision pending with the user (trend removal vs weighting).
