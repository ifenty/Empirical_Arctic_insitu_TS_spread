# Draft note to the AOTD author: reproduction of the monthly climatology

Subject: AOTD monthly climatology: 10 of 12 monthly files cannot be reproduced from the released code and scatter data

## Summary
We ported `cressman_bkg.m`, `barnes_season.m` and `barnes_month.m` to Python and ran them on the released scatter data (AOTD1983-2023, 414,221 profiles), comparing with the 12 released monthly climatology files (57 levels).

- **March and September are reproduced essentially exactly.** Median RMSE over levels: T 0.006 / 0.056 degC, S 0.001 / 0.005 psu (March / September); the monthly correction fields (month minus seasonal) correlate 1.00 with the released ones.
- **The other ten months are not reproduced.** Median RMSE over levels: T 0.31 degC, S 0.23 psu (worst levels near 10 m: about 0.9 degC and 1.2 psu). At level 20 the released correction relative to the seasonal field is about 0.15 degC (0.05 psu) RMS, against 0.3-0.7 degC (0.17-0.38 psu) from the released code, and the two are uncorrelated.

Valid-cell counts match the released files exactly for temperature in all months.

## What we ruled out (level 20, temperature, January)
- Sparse data: every month has 18k-32k valid observations at this level.
- Month-label shift: corrections built from each of the 12 months' observations correlate at most 0.37 with the released January correction.
- Year subsets of the observations (eras, single years, up to 2012/2020/2022, from 2005): at most 0.21.
- A short-radius, minimum-count variant suggested by the unused `Barnes()` function in `cressman_bkg.m`: correlation rises to about 0.55 at best, never near 1.

## Questions
1. Were the released monthly files for months other than March and September produced by `barnes_month.m` as released, or by an earlier/different version (different radius, kernel width, minimum observation count, first guess)?
2. Were they computed from the same scatter data (same QC version, same `AOTD.mo` month assignment)?
3. Could you share the code or parameters, or the intermediate seasonal files (`climatology_seasonal_XX.nc`) used?
4. The bathymetry file `topographic_land.mat` and `std_data.mat` referenced by the scripts are not in the release. Can they be provided?

## Other issues found in the released code
See `docs/upstream_bug_log.md` (BUG-1 to BUG-29), e.g. the `+1` damping when a single observation is in range (BUG-19), corrections applied only at cells containing observations (BUG-20), observations assigned to several nodes by a Euclidean degree radius (BUG-21), and `mean` versus `nanmean` for salinity (BUG-22).

Details and per-level statistics: `reports/reproduction/reproduction_report.md`.
