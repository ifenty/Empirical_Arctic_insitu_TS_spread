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

## Wedge/block patterns in the port (2026-10-13)
- No cap on profile count in the upstream code. Jan 1500 m S: max 2065 obs nodes within the 555 km Barnes radius; 57.5% of valid cells have 0 nodes in range (keep the seasonal first guess), 7% have 1-3.
- The wedges/blocks are the circular footprints of the hard radius cutoff (see `reports/reproduction/maps/count_vs_wedges_S_Jan_1500m.png`). Published Jan field is smooth there.
- Removing the BUG-19 `+1` damping at 1500 m raises Jan S RMSE vs published from 0.0057 to 0.64 psu, so released damping is what matches the published field (`maps/damping_test_S_Jan_1500m.png`).
- New product design: smooth taper/weights, minimum observation count and maximum distance, NaN beyond them (data-honest), no hard-radius jumps.

## 10 m salinity, Jan/Jun: smooth published vs polygon port (2026-10-13)
- Port polygons come from the normalised Barnes correction: weights exp(-r^2/(0.2n)) are ~1e-9 beyond a few hundred km, but the correction is (sum w*obs - z*sum w)/sum w, so the nearest observation node dominates (Voronoi-like cells). Published Mar/Sep show the same polygons (and match the port), so that is the released method.
- Published Jan/Jun are smooth with fine-scale structure: not produced by this method. Level 2 (10 m) Jan S vs published: seasonal field alone r=0.992, RMSE 0.35 psu; port r=0.902, RMSE 1.24; port with denominator sw+1 everywhere r=0.969, RMSE 0.69; sw+0.1 r=0.955; sw+0.01 r=0.943. Any variant of the released monthly correction is worse than no correction, consistent with BUG-29. Monthly correction fields uncorrelated with published (|r|<0.05).
- Implication for the new product: use a smooth, un-normalised-at-distance estimator (optimal interpolation / kriging-style weights with an explicit prior and error variance) so the field relaxes to the background, not to the nearest obs.

## OI climatology prototype (2026-10-13), levels 2/20/56 (10, 100, 1500 m) only
- Code: `aotd/oi.py` (staged annual -> seasonal -> monthly local OI, Gaussian covariance L=300/200/150 km, K=24 nearest nodes, node-mean obs with noise var sigma^2*(0.1+0.3/n), posterior error variance, mask), `scripts/oi_levels.py` (run + block-holdout CV, ~2 min per variable-level, 1.4 GB peak), `scripts/oi_maps.py`. Outputs `outputs/oi/level_k.npz` (untracked), `reports/oi/cv_stats.csv`, `reports/oi/maps/oi_{T,S}_L{k}.png`.
- No polygons: field relaxes to the previous stage away from data. Calibrated std maps are in the figures.
- CV (withheld 10% of 2x8 degree blocks, all years pooled): formal std is too small (z-std 2-6.5, robust 1.2-2.4); std is inflated by the robust z (`apply_mask(inflate=...)`). CV RMSE: T 1.34/1.01/0.26 degC, S 0.59/0.17/0.059 psu at 10/100/1500 m.
- Mask: NaN if nearest obs node (all months/years) > 150 km, or calibrated error variance > 0.25 * total node variance. Valid fraction of published-valid cells: T 0.74/0.86/0.90, S 0.99/0.98/0.36 (10/100/1500 m). Thresholds are tunable policy, to confirm with the user.
- Caveats: ocean mask still taken from the published fields; time pooled equally over 1983-2023 (no trend/decorrelation); nearest-node binning, no QC port yet; error propagation between stages is approximate (seasonal error carried in proportion to (P_m/sigma_m^2)); only 3 levels run. Not yet validated against Mar/Sep beyond the RMSE/corr columns in cv_stats.csv.

## Moving-window OI variances (2026-10-13)
- `aotd/oi.py`: prior variance is now local (mean squared residual of the Kw=50 nearest nodes within Rw=500 km, noise-deflated, shrunk toward the level value with n0=10, floor 2%); covariance sigma_i*sigma_j*rho; posterior variance uses the local prior. `local_sd` = std of node means in the same window = local natural variability.
- Mask: NaN if nearest node > 150 km or calibrated std > relerr_max * local_sd (relerr_max now 1.0). `scripts/oi_remask.py` re-applies a mask to saved levels and prints Jan coverage per threshold.
- CV z-std improved slightly (robust 1.0-2.2). Jan valid fraction at relerr_max=1.0: T 0.28/0.30/0.14, S 0.73/0.57/0.10 (10/100/1500 m).
- Finding: a relative rule NaNs the central Arctic at 10 m T, where T is near freezing and local_sd is tiny even though absolute errors are small (see `reports/oi/maps/oi_T_L2.png`). Probably needs an absolute floor (max of relative and absolute tolerance); thresholds are a user policy decision.
