# Upstream MATLAB bug / issue log

Running log of problems found in https://github.com/niuda-niuer/An-Arctic-Ocean-Thermohaline-Dataset
(copy in `upstream_matlab/`). Each item: file, what is wrong, evidence, how the new code handles it.
IDs are referenced from code comments (`# BUG-nn`). This is the source for the note to the paper's authors.

Status legend: CONFIRMED (verified numerically or by inspection), SUSPECTED (needs test), DESIGN (not a bug, a limitation).

## Reproducibility / packaging
- **BUG-01 (CONFIRMED)** Hard-coded absolute paths (`/data2/ljl/...`), `parpool(32/52)`, and a missing input
  `topographic_land.mat` (bathymetry `MASK`, ETOPO5-derived). Code cannot run as published. README gives no
  data-flow description; script order (duplicate -> QC -> three_std -> cressman_bkg -> barnes_season -> barnes_month)
  is implicit. `clearvars -except ...` lists variables that do not exist in the repo.
- **BUG-02 (CONFIRMED)** Published scatter NetCDF files have an invalid time unit (`since 1980-01-01 00:00:00`, no unit)
  so xarray/NCO/CF tools cannot decode time (values are days). Dimensions of the per-profile arrays are named
  `lat`, `lon`, `time` (all length n_profiles) which is not CF-compliant and breaks tools expecting `lat/lon` axes.
- **BUG-03 (CONFIRMED)** The `source` ID assigned in `duplicate.m` (1=ITP, 2=Argo, 3=BOL, 4=NABOS, 5=CNarc, 6=WOD,
  7=CMEMS) is not in the published scatter files, so users cannot filter by platform/provenance.
- **BUG-04 (CONFIRMED)** Grid is on lat = 60,60.25,...,90 and lon = -180,...,179.75 (values on whole 0.25 multiples)
  whereas WOA23 0.25 deg uses cell centres (60.125...). Published `geolon/geolat` are 2-D copies; no `lat`/`lon` 1-D
  coordinate variables, no units/attributes, 2-D coordinates only. Includes both lon=-180 and lat=90 (a single point).

## Quality control (ocean_*_quality_check.m)
- **BUG-05 (CONFIRMED)** `ocean_vars_quality_check.m`: `Tmin` literal contains `-2-2` (evaluates to -4, merging two
  elements): the vector has 32 elements while `depth_std` has 33, so `[Tmax,Tmin,depth_std]` concatenation would
  error / is misaligned in MATLAB (horzcat dimension mismatch). `Smax`/`Smin` are 33 long; `Tmax`,`Tmin` lengths differ.
  Depth-dependent thresholds therefore cannot be what the paper's Table 2 states.
- **BUG-06 (CONFIRMED)** Time test uses `time>=start & time<end`, lat/lon tests use strict `>`/`<` bounds
  (`lat>60`, `lon>-180`) so profiles exactly at 60N, 90N, +-180 are rejected although paper says 60-90N.
- **BUG-07 (CONFIRMED)** Temperature and salinity gradient test compares each level only with the *first* later level
  at least 3 m deeper (`break` after first hit). In the interpolated data (5 m spacing near surface) this is
  a 2-level difference; thresholds apply to inconsistent depth separations; no handling of unequal spacing.
- **BUG-08 (CONFIRMED)** Gradient test, last-point rule: `ind0 = max(find(local_qc~=0,1,'last')-1,1); if local_qc(ind0)==2 ... local_qc(ind0+1)+1`
  gives the last valid point a pass only based on its neighbour; and the case of exactly one valid point is
  accepted when only the range test passed.
- **BUG-09 (CONFIRMED)** `ocean_dens_quality_check.m`: the gradient `Gd` has one fewer element than `d`, but
  `local_qc(ind)` is indexed with a logical vector of length n-1 into a vector of length `wid`, and uses
  pressure thresholds on `pp(2:end)` while skipping NaNs; the resulting flags are shifted by one level relative to
  the density array (and misaligned when NaNs are interleaved). Inversion threshold is a *gradient* (kg m-4) but
  quoted in the paper as a density-step; units/sign conventions undocumented.
- **BUG-10 (CONFIRMED)** `ocean_pres_quality_check.m`: loop logic flags repeats incorrectly (only the second
  element is rescued when `j==2`; later repeats are silently `continue`d without being flagged), and the
  `length(pp)==1` branch sits inside the loop and is never reached for a one-point profile.
- **BUG-11 (CONFIRMED)** `ocean_constant_pres_vars.m` removes points only where *pressure* is NaN, then applies
  `p1(p1>=6000)=nan` *after* the subset, leaving s1,t1 with values where p1 is NaN (later counts `~isnan(s1)` are
  computed on s1, not on the valid-pair mask).
- **BUG-12 (CONFIRMED)** `ocean_vars_interp.m`: `salt1(1:idx1)=nan; salt1(idx2:end)=nan` blanks *the nearest
  standard level to the first/last observation as well as everything beyond it* - i.e. extrapolation guard
  removes valid data at the profile ends (up to half a level thickness of real data). With 100-500 m level spacing
  this is up to 250 m of valid data lost per profile at depth. In the <=3 point branch the guard uses `p1` while
  the interpolation used the unfiltered `p`, mixing arrays.
- **BUG-13 (CONFIRMED)** `ocean_vars_interp.m`: temperature interpolated as `gsw_t_from_CT(salt1,temp1,target)`
  where `salt1` is *practical salinity* (SP) but the function needs *absolute salinity* SA (`gsw_t_from_CT(SA,CT,p)`):
  in-situ temperature is converted with the wrong salinity variable, a systematic (small, ~1e-3 degC) error and
  incorrect where `salt1` is NaN (then the temperature is NaN even though CT exists).
- **BUG-14 (CONFIRMED)** `ocean_vars_interp.m`: `pres = repmat(target,1,len)` then `Sa = gsw_SA_from_SP(salt,pres,lon,lat)`
  with `pres` set to the *nominal depth in m treated as dbar* (target is depth levels from WOA in metres) - no
  depth->pressure conversion (gsw_p_from_z); error up to ~1% at 5500 m.
- **BUG-15 (CONFIRMED)** `duplicate.m`: deduplication is order dependent (first source wins; source order ITP, Argo, BOL,
  NABOS, CNarc, WOD, CMEMS is arbitrary) and not symmetric; the test is `dt<6h & |dlat|<0.1 & |dlon|<0.1` where
  `|dlon|` is *degrees of longitude* (not km: 0.1 deg of lon is 11 km at the equator but only 1 km at 85N); no
  dateline wrap (lon 179.99 vs -179.99 treated as 360 deg apart). Paper says 5 hours, code uses 6 hours.
  Comparison is O(N^2) (loops over every new profile).
- **BUG-16 (CONFIRMED)** `three_std.m`: boxes and 'seasons' for the 3-sigma test are keyed
  T=1 -> AMJ, 2 -> JAS, 3 -> OND, 4 -> JFM, while `barnes_season.m`/`barnes_month.m` use 1 -> JFM, 2 -> AMJ, 3 -> JAS, 4 -> OND
  (inconsistent season indexing). Only `k = 2` (the 5 m level!) is processed in `for k = 2%1:102` - the loop
  limit was left in a debug state so the published 3-sigma test could only have applied to a single level (or the
  published code differs from what was run). `std(~isnan(tk))` tests the std of a *logical mask* (almost always 0 or ~0.5)
  instead of the data. `rmoutliers` (MAD-based, 3 scaled MAD) is applied before the std so the '3 sigma' is
  computed on already trimmed data.
- **BUG-17 (CONFIRMED)** `three_std_check.m`: `load('ts_range_v5.mat')` replaces the arrays just computed (undocumented
  input); grid index search uses `distance` along lat/lon separately instead of true nearest neighbour; the
  `for k=1:102` loop is O(Nprof*102) in MATLAB; `t1`/`s1` are NaN at levels where T missing but S valid.

## Objective analysis (cressman_bkg.m, barnes_season.m, barnes_month.m)
- **BUG-18 (CONFIRMED)** Cressman/Barnes corrections use the residual `obs - field(i,j)` where `(i,j)` is the *target*
  grid point, not the field value interpolated to the observation location `b` as in the paper's Eq. (1)
  (`f_b^n`). Equivalent only when observation and target coincide; for every other observation in the
  influence radius the residual is wrong: a gradient across the radius becomes a spurious residual that is
  added back, producing smoothing bias and unphysical amplification where gradients are strong (shelf breaks,
  fronts).
- **BUG-19 (CONFIRMED)** First pass: weights sum `+1` in the denominator if exactly one observation is within the radius
  (`numel(tt)==1`) - an ad-hoc damping that is not in the paper, creates a discontinuity (one obs -> value
  shrunk toward 0 degC / 0 psu in pass 1 background because the background is *not* a first-guess anomaly but the
  absolute value). A single observation at 0 distance (weight 1) gives T/2; for salinity of 34 gives 17 psu in the
  Cressman background (later corrected in later passes only where an observation exists on that grid cell).
- **BUG-20 (CONFIRMED)** Cressman passes 2-3 and the seasonal/monthly corrections are applied *only at grid points
  that contain an observation* (`if ~isnan(zo(i,j))`). Cells without observation keep the background value from
  pass 1 (pure distance-weighted mean), so the final field is discontinuous between observed and unobserved cells
  (visible as speckle/ 'bullseye' structures around isolated observations).
- **BUG-21 (CONFIRMED)** `obs_merge`: observations are averaged into grid cells with a *Euclidean radius of 0.25 deg in
  (lat,lon) degrees* (`abs(del_x + del_y*1i) <= 0.25`). In km this footprint shrinks with cos(lat): at 85N a 0.25 deg
  longitude is 2.4 km; at 60N it is 14 km. Each observation is also assigned to up to 5 neighbouring grid points
  (centre and 4 neighbours at exactly 0.25) so the same data are counted several times (non-independent obs, reduces
  effective spatial resolution) and the weights are not volume-normalised. Only the `lon==-180` column is treated
  for the dateline (obs at 179.9 never reach lon=-179.75).
- **BUG-22 (CONFIRMED)** `barnes_season.m`/`barnes_month.m` `obs_merge`: `mean()` (not `nanmean`) is used for salinity with
  `s(~isnan(t))` only: any NaN salinity among temperature-valid profiles in a cell turns the whole cell's
  salinity into NaN; temperature and salinity use different valid sets. `cressman_bkg.m` uses `nanmean` -> annual
  and seasonal are not computed consistently.
- **BUG-23 (CONFIRMED)** Every grid cell/level uses the same horizontal influence radii (321/267/214 km seasonal;
  555 km Barnes for months) regardless of data density or local deformation radius; the paper itself states
  seasonal variability is weak because of this. In sparse regions (East Siberian/Kara) the first guess dominates.
- **BUG-24 (CONFIRMED)** `Barnes_sc` radius R: `barnes_season.m` has R=999 km (unused debug path), `barnes_month.m` has R=555 km,
  `Barnes()` stand-alone has `max(...,1000)`; the paper only gives alpha=8e4/1.6e4 km2 and gamma=0.2, never R=555.
  Weight in the code is `exp(-r^2/(0.2*n))` with n = 8e4 then 1.6e4 km2 (so e-folding scale 126 km then 57 km).
  The cut-off at 555 km truncates a Gaussian whose weight at the cut-off is exp(-555^2/16000)=4e-9 (pass 1) -
  harmless, but the second-pass kernel (57 km) with `n=16e3` uses the same R so it is wasted. No convergence check.
- **BUG-25 (CONFIRMED)** First guess for each month is the *seasonal* field, which is itself the annual field corrected only
  at cells with seasonal observations (BUG-20): months in a season share an identical first guess (no true
  seasonal cycle in unobserved regions - "seasonal variability is relatively weak" in the paper).
  Season mapping `JFM/AMJ/JAS/OND` is not centred on the solstices/equinoxes (WOA also uses this convention, OK) but
  the month->season mapping in barnes_month uses 1..4 indexes inconsistently with three_std.
- **BUG-26 (CONFIRMED)** Time is not weighted: all years 1983-2023 are pooled into a single climatology
  (strongly biased to 2005+ and to summer; no attempt to remove the Arctic warming trend / sea-ice change
  signal; the 1983-2023 mean is not a stationary climate state). Per-box averaging is also not weighted by profile
  count or error variance.
- **BUG-27 (SUSPECTED)** Land mask is `MASK==-1`; `mask(i,j)>=depth` uses ETOPO5 (5 arcmin, ~9 km): shallow shelves,
  straits (CAA, Fram) misrepresented; the comparison to WOA23 shows ~2% disagreement in valid points.
- **BUG-28 (CONFIRMED)** Density: `Lt=1000`, `h=fspecial` and `obs_modas`, `Barnes()` are dead code but left in scripts, with different
  radii constants; `fi`, `Lt` unused.

- **BUG-29 (CONFIRMED, reproduction finding)** The published monthly files for months other than March and September
  do not contain the monthly observation signal that `barnes_month.m` produces from the released scatter data.
  At level 20 (T and S), the port reproduces published March and September to RMSE 0.01 degC / 0.00 psu (monthly
  correction fields correlate 1.0 with the published ones). For the other ten months the published correction
  relative to the seasonal field has an RMS of about 0.15 degC (0.05 psu), against 0.3-0.7 degC (0.17-0.38 psu) in the
  port, and the two are uncorrelated (|r| < 0.4). All months have 18k-32k valid observations at this level.
  Observations from no single other month reproduce the published January or June correction (max r = 0.37), so this
  is not a month-label shift. Likely cause: the released monthly files for those months came from a different or
  incomplete run than the released code and scatter data. Needs confirmation from the author.
