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
