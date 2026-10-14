"""Summarize the expected-spread NetCDF products and caption the figures from script 04."""
from contextlib import ExitStack
from pathlib import Path

import numpy as np
import xarray as xr

D = Path("outputs/expected_spread")
OUT = Path("reports/expected_spread/expected_spread_report.md")
PRODS = ("woa23", "aotd", "woa23_and_aotd")
NAMES = {"Theta": "potential_temperature", "S": "salinity"}
DEPTHS = (10, 50, 100, 200, 500, 1000, 1500, 2500, 4000)


def load(kind, prod, v):
    return xr.open_dataset(D / f"expected_{kind}_from_{prod}_{NAMES[v]}_ease2_100km.nc")


def fmt(x):
    return "–" if not np.isfinite(x) else f"{x:.3g}"


def table(kind, v):
    rows = [
        "| depth (m) | " + " | ".join(f"{p} median [P10–P90]" for p in PRODS)
        + " | valid cells: WOA23 / AOTD / combined |",
        "|---|---|---|---|---|",
    ]
    with ExitStack() as stack:
        ds = {p: stack.enter_context(load(kind, p, v)) for p in PRODS}
        dep = ds["woa23"].depth.values
        fields = {p: ds[p][f"{v}_{kind}"].values for p in PRODS}
        for d in DEPTHS:
            k = int(np.argmin(np.abs(dep - d)))
            cells, counts = [], []
            for p in PRODS:
                f = fields[p][k]
                f = f[np.isfinite(f)]
                counts.append(str(f.size))
                cells.append(
                    "–" if not f.size else
                    f"{np.median(f):.3g} [{np.percentile(f, 10):.2g}–{np.percentile(f, 90):.2g}]"
                )
            rows.append(f"| {dep[k]:.0f} | " + " | ".join(cells) + " | " + " / ".join(counts) + " |")
    return "\n".join(rows)


def counts_table():
    rows = [
        "| depth (m) | values (`n`) | days (`n_eff`) | 1-mo epochs (`n_epochs`) | median years (`n_years`) | cells `n >= 1` (`n >= 2` / `M >= 2`) | % with `n < 5` | % with `n_epochs < 5` |",
        "|---|---|---|---|---|---|---|---|",
    ]
    with load("rmse", "woa23", "Theta") as ds:
        dep = ds.depth.values
        counts = ds["Theta_n"].values
        neff = ds["Theta_n_eff"].values if "Theta_n_eff" in ds else counts
        nep = ds["Theta_n_epochs"].values if "Theta_n_epochs" in ds else neff
        nyrs = ds["Theta_n_years"].values if "Theta_n_years" in ds else counts
        for d in DEPTHS:
            k = int(np.argmin(np.abs(dep - d)))
            n, ne, np_ep, ny = counts[k], neff[k], nep[k], nyrs[k]
            has = n > 0
            has2 = int((n >= 2).sum())
            has_ep2 = int((np_ep >= 2).sum())
            med_yr = fmt(np.median(ny[has])) if has.any() else "–"
            frac_n5 = fmt(100 * ((n < 5) & has).sum() / has.sum()) if has.any() else "–"
            frac_ep5 = fmt(100 * ((np_ep < 5) & has).sum() / has.sum()) if has.any() else "–"
            rows.append(
                f"| {dep[k]:.0f} | {int(n.sum()):,} | {int(ne.sum()):,} | {int(np_ep.sum()):,} | {med_yr} | "
                f"{int(has.sum())} ({has2} / {has_ep2}) | {frac_n5} | {frac_ep5} |"
            )
    return "\n".join(rows)


def declustered_table(v):
    rows = [
        f"| depth (m) | Profile-weighted `{v}_std` (`n >= 2`) | Epoch-declustered `{v}_std_declustered` (`M >= 2`) | EB regularized `{v}_std_se` (500 km prior) | EB `{v}_std_se_pan_arctic` (pan-Arctic prior) | `{v}_std_within` / `between` | `{v}_rmse`: profile-wtd / declust / EB (500km) |",
        "|---|---|---|---|---|---|---|",
    ]
    with load("std", "woa23", v) as ds_s, load("rmse", "woa23", v) as ds_r:
        dep = ds_s.depth.values
        if f"{v}_std_declustered" not in ds_s:
            return ""
        for d in DEPTHS:
            k = int(np.argmin(np.abs(dep - d)))
            def m_p(arr):
                g = arr[np.isfinite(arr)]
                return "–" if not g.size else f"{np.median(g):.3g} [{np.percentile(g, 10):.2g}–{np.percentile(g, 90):.2g}]"
            def m_only(arr):
                g = arr[np.isfinite(arr)]
                return "–" if not g.size else f"{np.median(g):.3g}"
            s_raw = m_p(ds_s[f"{v}_std"].values[k])
            s_dec = m_p(ds_s[f"{v}_std_declustered"].values[k])
            s_se  = m_p(ds_s[f"{v}_std_se"].values[k])
            s_pan = m_p(ds_s[f"{v}_std_se_pan_arctic"].values[k]) if f"{v}_std_se_pan_arctic" in ds_s else "–"
            s_wb  = f"{m_only(ds_s[f'{v}_std_within'].values[k])} / {m_only(ds_s[f'{v}_std_between'].values[k])}"
            r_all = " / ".join(m_only(ds_r[nm].values[k]) for nm in (f"{v}_rmse", f"{v}_rmse_declustered", f"{v}_rmse_se"))
            rows.append(f"| {dep[k]:.0f} | {s_raw} | {s_dec} | {s_se} | {s_pan} | {s_wb} | {r_all} |")
    return "\n".join(rows)


def bias_table():
    rows = [
        "| depth (m) | Theta bias median (°C) | Theta bias P10–P90 | S bias median (psu) | S bias P10–P90 |",
        "|---|---|---|---|---|",
    ]
    with load("rmse", "woa23", "Theta") as t, load("rmse", "woa23", "S") as s:
        dep = t.depth.values
        fields = {"Theta": t["Theta_bias"].values, "S": s["S_bias"].values}
        for d in DEPTHS:
            k = int(np.argmin(np.abs(dep - d)))
            r = [f"{dep[k]:.0f}"]
            for v in ("Theta", "S"):
                f = fields[v][k]
                f = f[np.isfinite(f)]
                r += [fmt(np.median(f)) if f.size else "–",
                      f"{np.percentile(f, 10):.2g} … {np.percentile(f, 90):.2g}" if f.size else "–"]
            rows.append("| " + " | ".join(r) + " |")
    return "\n".join(rows)


def fig(name, cap):
    return f"{cap}\n\n![{name}](figures/{name})\n"


def report():
    return f"""# Empirical spread of Arctic profiles relative to climatology

Generated by `scripts/05_report_expected_spread.py` from the local NetCDF products in
`outputs/expected_spread/`. See the [technical note](../../technical_note.md) for the
method, the [README](../../README.md) for reproduction, and the
[reference comparison](../woa_compare/woa_compare_report.md) for climatology diagnostics.

## 1. Scope and definitions

The products summarize finite observation-minus-reference differences within each
100 km EASE-Grid 2.0 North cell (EPSG:6931; 180 × 180 cells) and each of the 102 WOA
standard depth levels. The computation selects AOTD scatter profiles from 1992 onward
and pools years and months; there is no time dimension in the output.
Temperature differences use surface-referenced potential temperature (Theta);
salinity differences use practical salinity (S; conventionally reported here as psu).

* **Profile-weighted fields (`*_rmse`, `*_std`, `*_bias`)**:
  assign equal weight (`1 / n`) to every finite comparison value:
  * **RMSE** = √mean(d²): profile-weighted mismatch relative to the reference, including the mean offset (defined when `n >= 1`).
  * **std** = √(n / (n − 1) · mean((d − mean(d))²)): profile-weighted **sample** standard deviation (`ddof=1`) after removing the cell's mean offset. It is `NaN` when `n < 2`, preventing zero-variance artifacts in singleton cells.
  * **bias** = mean(d), where d = observation − reference. This is a sample mean offset, not an independently established climatology error.
* **Epoch-declustered fields (`*_std_declustered`, `*_rmse_declustered`, `*_bias_declustered`, `*_std_within`, `*_std_between`)**:
  group observations in each 100 km cell and depth into **1-month synoptic epochs** `(cell, year, month)` (count `M = n_epochs`),
  assigning equal weight `1 / M` to each epoch and restoring pooled within-epoch variance via the Law of Total Variance.
* **Empirical Bayes (EB) state-estimation regularized fields (`*_std_se`, `*_rmse_se`) and pan-Arctic comparison (`*_std_se_pan_arctic`, `*_rmse_se_pan_arctic`)**:
  apply Empirical Bayes shrinkage (`nu_0 = 5` epochs) of the epoch-declustered variance toward a **500 km circular moving-window median prior** across neighboring well-sampled cells (`M >= 5`, with fallback to `M >= 2` within 500 km and then the depth-level median), providing a regime-preserving, non-degenerate error estimate for every sampled cell (`M >= 1`). For comparison, `*_se_pan_arctic` shrinks toward a single horizontally uniform pan-Arctic depth-level prior.
* **n** counts finite differences at a cell and depth, **n_eff** counts distinct `(100 km cell, UTC calendar day)` sampling days, **n_epochs** (`M`) counts distinct `(100 km cell, year, month)` 1-month synoptic epochs, and **n_years** counts distinct calendar years (1992–2023) represented in the cell.

The `woa23` product uses the WOA23 0.25° monthly climatology at 0–1500 m and seasonal
climatology below 1500 m, **bilinearly interpolated in space and linearly interpolated in time**
(between mid-months or mid-seasons, with coast/mask-aware weight renormalization over finite corners).
The `aotd` product uses the published monthly AOTD climatology at 0–1500 m only, interpolated identically.
`woa23_and_aotd` takes the larger finite metric at each cell and level;
when only one reference is available it uses that reference. Below 1500 m it therefore
equals the WOA23 product. Its `bias`, `n`, `n_eff` and `n_years` follow the selected reference **separately for
std and RMSE**, so the two combined files may have different bias and count fields.

A metric is NaN when insufficient finite comparisons are available (`n < 1` for RMSE, `n < 2` for sample std), including when the reference
is missing. No spatial filling of the output spread grid is applied. These empirical summaries are neither pure
observation-error estimates nor validated forecast uncertainties. The combined maximum
is a reference-sensitivity envelope, not a calibrated confidence bound.

## 2. Spread versus depth

Tables and depth curves summarize **each product's own finite cells**, with equal weight
per cell, not per profile. The masks can differ between references and depths; lower
medians do not by themselves establish a better reference. P10–P90 is the spatial
percentile range, not a confidence interval. Counts for all three masks are shown.

{fig("spread_vs_depth.png", "**Figure 1. Spatial median RMSE, std and cellwise RMSE/std at 0–1500 m.** Top: Theta; bottom: S. WOA23 shading shows P10–P90 across finite cells. Ratios require n ≥ 5 for both metrics and positive std; the combined ratio can divide metrics selected from different references.")}
{fig("spread_vs_depth_full.png", "**Figure 2. WOA23-based summaries over 0–5500 m.** As in Figure 1, with a symlog depth axis. Sparse deep sampling should be considered alongside the curves.")}
### Theta RMSE (°C)

{table("rmse", "Theta")}

### Theta std (°C)

{table("std", "Theta")}

### S RMSE (psu)

{table("rmse", "S")}

### S std (psu)

{table("std", "S")}

## 3. Spatial distributions

The overview figures show Theta and S metrics at 10, 100 and 1500 m, log10 counts,
and depth summaries. Count maps describe the available comparisons, not an effective
sample size. Colour scales are panel-specific; consult the colour bars before comparing
magnitudes. The maps describe spatial structure but do not separate the processes
responsible for it.

| reference | std | RMSE |
|---|---|---|
| WOA23 | ![WOA23 std overview](figures/expected_std_overview_from_woa23.png) | ![WOA23 RMSE overview](figures/expected_rmse_overview_from_woa23.png) |
| AOTD | ![AOTD std overview](figures/expected_std_overview_from_aotd.png) | ![AOTD RMSE overview](figures/expected_rmse_overview_from_aotd.png) |
| Combined maximum | ![Combined std overview](figures/expected_std_overview_from_woa23_and_aotd.png) | ![Combined RMSE overview](figures/expected_rmse_overview_from_woa23_and_aotd.png) |

## 4. Mean offsets and reference sensitivity

{fig("bias_maps_from_woa23.png", "**Figure 3. Mean observation-minus-WOA23 offsets at 10, 100 and 1500 m.** Each symmetric colour range is set by the 98th percentile of the absolute cell offset.")}
{bias_table()}

{fig("product_difference_maps.png", "**Figure 4. AOTD-based RMSE minus WOA23-based RMSE where both metrics exist.** Blue denotes a smaller empirical AOTD-based RMSE, not independent validation of AOTD. Even within a shared output cell, the contributing observations can differ because the reference masks differ.")}
## 5. Available comparisons, 1-month epoch de-clustering, and state-estimation fields

{fig("counts_vs_depth.png", "**Figure 5. WOA23-based potential-temperature comparison counts versus depth.** Panels show total finite differences (`n`), distinct sampling days (`n_eff`) and 1-month epochs (`n_epochs`), cells with data (`n >= 1`, `n >= 2` and `n_epochs >= 2`), and sparse-cell percentages conditional on `n > 0`.")}
{counts_table()}

Because drifting Ice-Tethered Profilers and dense shipboard surveys can collect hundreds of profiles within a single month, **profile-weighted** pooling (`1/n`) pins the cell mean and standard deviation to individual synoptic deployments. Grouping observations into **1-month synoptic epochs** `(cell, year, month)` (count `M = n_epochs`, assuming a ~30-day Arctic synoptic/eddy decorrelation timescale) and applying the Law of Total Variance produces the **epoch-declustered** estimator, which separates **within-epoch** sub-monthly/sub-grid spread (`*_std_within`) from **between-epoch** synoptic/interannual spread (`*_std_between`). The **Empirical Bayes (EB) regularized** estimator (`*_std_se`, `*_rmse_se`) then shrinks low-`M` cells (`nu_0 = 5` epochs) toward a **500 km circular moving-window median prior** across neighboring well-sampled cells (`M >= 5`) so cells visited in only `M = 1` or a few epochs are smoothly regularized for ocean state estimation.

{fig("declustered_spread_vs_depth.png", "**Figure 6. Profile-weighted, epoch-declustered (1-month epochs), and Empirical Bayes (EB) regularized spread versus depth (from WOA23).** Left: sample std comparison; middle: decomposition into within-epoch (sub-monthly) and between-epoch (synoptic/interannual) std; right: RMSE comparison.")}
{fig("method_distributions_by_depth.png", "**Figure 7. Empirical cumulative distribution functions (CDFs) of profile-weighted, epoch-declustered, and Empirical Bayes regularized (500 km vs. pan-Arctic prior) std and RMSE across cells at 10, 50, 100, 300 and 1000 m (from WOA23).** Epoch de-clustering shifts the distribution toward larger spreads by removing single-month ITP dominance; Empirical Bayes shrinkage regularizes sparse `M = 1–4` cells toward the local 500 km prior.")}
{fig("declustered_std_maps_from_woa23.png", "**Figure 8. Empirical Bayes regularized std (`*_std_se`, rows 1 & 3) and difference between EB regularized and profile-weighted std (`*_std_se - *_std`, rows 2 & 4) at 10, 100 and 1500 m (from WOA23).** Red indicates regions where 1-month epoch de-clustering and local 500 km prior regularization increase the spread over profile-weighted pooling (both in ITP-dominated basins and in sparsely sampled cells); blue indicates noisy small-sample cells pulled back toward the local 500 km prior.")}
### Multi-depth map comparisons: Profile-weighted vs. Epoch-declustered vs. Empirical Bayes regularized

Each five-depth figure (10, 50, 100, 300 and 1000 m) shows the profile-weighted metric (column 1), 1-month epoch-declustered metric (column 2), Empirical Bayes regularized metric (`*_se`, column 3, all three sharing a common colour scale per depth row), and the difference `EB regularized − profile-weighted` (column 4).

**Why Empirical Bayes regularized (`*_se`, 500 km local prior) differs from profile-weighted spread, and how the 500 km local prior compares with a single pan-Arctic prior (`*_se_pan_arctic`):**

1. **Why EB regularized (`*_se`) is larger than profile-weighted spread in ITP-dominated basins (Beaufort Gyre and Central Arctic) while preserving local hydrographic regimes:**
   * *Epoch de-clustering (removal of single-epoch ITP domination):* In the ice-covered Central Arctic (`> 80°N`) and Canada Basin, a single drifting Ice-Tethered Profiler often records dozens to hundreds of profiles within one 100 km cell during a single month (at 10 m north of 80°N, the median cell has `n = 51` profiles but only `M = 5` distinct 1-month epochs). Profile-weighted (`1/n`) pooling gives overwhelming weight to that single-month burst, measuring primarily sub-monthly within-epoch jitter (`*_std_within`). Weighting each 1-month epoch equally (`1/M`) restores the between-epoch variance (`*_std_between`), increasing median Central Arctic 10 m `Theta_std` from `0.044 °C` (profile-weighted) to `0.063 °C` (epoch-declustered) and `S_std` from `0.436 psu` to `0.580 psu`.
   * *500 km local prior vs. single pan-Arctic prior in the Central Arctic:* Under perennial sea ice north of 80°N, near-surface water is thermally buffered near the freezing point. With a **single pan-Arctic prior** (`*_std_se_pan_arctic`, which pools all cells north of 60°N into a single 10 m prior of `0.82 °C`), low-`M` Central Arctic cells were pulled upward to a physically unrealistic median of `0.690 °C` (`+0.65 °C` above profile-weighted). With the **500 km moving-window prior** (`*_std_se`), each Central Arctic cell draws its prior from well-sampled cells (`M >= 5`) within 500 km inside the same ice-covered regime, keeping the median 10 m Central Arctic `Theta_std_se` at a physically realistic **`0.071 °C`** while regularizing `S_std_se` to **`0.655 psu`**.
2. **Why EB regularized (`*_se`) damps isolated small-sample spikes along the Russian coast while retaining regional river-plume variability (unlike the single pan-Arctic prior):**
   * Along the Kara, Laptev and East Siberian shelves (`68–78°N, 60–175°E`), river-runoff plumes and summer heating create high variance (median 10 m profile-weighted `S_std = 1.50 psu`, epoch-declustered `1.73 psu`; profile-weighted `Theta_std = 0.727 °C`, epoch-declustered `0.984 °C`), but sampling is temporally sparse (median `n = 6` profiles across `M = 3` 1-month epochs).
   * A **single pan-Arctic prior** (`~0.55 psu` for 10 m `S`) over-damped the entire Russian shelf downward to a median of `1.03 psu` (`-0.47 psu` below profile-weighted). By contrast, the **500 km moving-window prior** draws its prior from neighboring shelf cells within 500 km, **preserving the high shelf regime** (median Russian coast 10 m `S_std_se = 1.73 psu` and `Theta_std_se = 1.30 °C`) while still smoothing out isolated single-cruise extreme outliers (`> 3–4 psu`, blue coastal cells in Column 4) toward the 500 km shelf median.

| variable | standard deviation (`std` vs. `std_declustered` vs. `std_se`) | root-mean-square error (`rmse` vs. `rmse_declustered` vs. `rmse_se`) |
|---|---|---|
| Potential temperature (`Theta`) | ![Theta std multi-depth comparison](figures/method_comparison_std_Theta_from_woa23.png) | ![Theta RMSE multi-depth comparison](figures/method_comparison_rmse_Theta_from_woa23.png) |
| Salinity (`S`) | ![S std multi-depth comparison](figures/method_comparison_std_S_from_woa23.png) | ![S RMSE multi-depth comparison](figures/method_comparison_rmse_S_from_woa23.png) |

### Within-epoch vs. between-epoch variance decomposition across depths

Each five-depth figure (10, 50, 100, 300 and 1000 m) decomposes the 1-month epoch-declustered spread into sub-monthly within-epoch spread (`*_std_within`, column 1), inter-epoch synoptic/interannual spread (`*_std_between`, column 2), total epoch-declustered spread (`*_std_declustered`, column 3), and the fraction of total variance explained by between-epoch variations (`std_between² / std_declustered²`, column 4).

| Potential temperature (`Theta`) variance components | Salinity (`S`) variance components |
|---|---|
| ![Theta variance components](figures/variance_components_Theta_from_woa23.png) | ![S variance components](figures/variance_components_S_from_woa23.png) |

### Direct comparison: EB regularized with 500 km local moving-window prior (`*_std_se`) vs. single pan-Arctic prior (`*_std_se_pan_arctic`)

The figures below compare epoch-declustered spread (`*_std_declustered`, column 1), EB regularized spread with the **500 km moving-window prior** (`*_std_se`, column 2), EB regularized spread with the **single pan-Arctic prior** (`*_std_se_pan_arctic`, column 3), and their difference `500 km prior − pan-Arctic prior` (column 4) across 10, 50, 100, 300 and 1000 m. In Column 4, blue in the Central Arctic shows where the 500 km prior avoids importing warm Atlantic-inflow temperature variance into the ice-covered basin, while red in the Barents/Norwegian Seas and Russian shelves shows where the 500 km prior preserves high boundary-current and river-plume variability.

### De-clustered and state-estimation Theta spread from WOA23 (°C)

{declustered_table("Theta")}

### De-clustered and state-estimation S spread from WOA23 (psu)

{declustered_table("S")}

## 6. Interpretation and limitations

* Mismatch combines observational uncertainty, unresolved spatial and temporal
  variability, interannual change, trends, and reference limitations. This analysis
  does not identify their individual contributions.
* For a single reference and the same samples with `n >= 2`, `RMSE² = ((n − 1) / n) · std² + bias²`.
  This identity need not hold across the separately selected combined `std` and `RMSE` products.
* AOTD scatter and climatology share 1983–2023 source data, making `from_aotd` an **in-sample**
  residual calculation; the extent of overlap with WOA23 is not established here. Smaller mismatch
  relative to AOTD reflects in-sample fitting and shared processing rather than out-of-sample skill.
* Reference values are bilinearly interpolated in space and linearly interpolated between mid-months
  (or mid-seasons below 1500 m) with mask-aware renormalization near coasts. AOTD is mapped to WOA23
  centres under an adopted grid-alignment assumption.
* The [climatology comparison](../woa_compare/woa_compare_report.md) documents a
  marked deep-salinity anomaly (800–1500 m) in ten AOTD months, which inflates `from_aotd` and
  `woa23_and_aotd` salinity spread where deep profiles sample those months.
* Applying these pooled spreads to assimilation or prediction requires minimum-count filtering
  (e.g., on `n_eff` and `n_years`), spatial regularization where cells are `NaN`, and a treatment
  of correlated errors.
* **Guidance for ocean state estimation (`*_std_se` vs. `*_rmse_se`):** When weighting profile
  observations against an ocean model state (`J_obs = Σ ((y_i − H(x)) / σ_i)²`),
  use **`*_std_se`** (`Theta_std_se`, `S_std_se` with the 500 km local prior) rather than `*_rmse_se`.
  `*_std_se` estimates the random sub-grid, synoptic and interannual spread around the local mean,
  whereas `*_rmse_se` includes the squared mean offset (`*_bias_declustered`) relative to the
  reference climatology—which is often an artifact of climatological smoothing across fronts or
  multi-decadal warming/freshening trends rather than observation error. Use `*_rmse_se` only when
  bounding total departure from a static WOA23/AOTD background state (e.g. background error around
  a climatological initial condition), and apply a modest sensor-accuracy floor at abyssal depths
  (e.g. `σ_Theta >= 0.01–0.02 °C`, `σ_S >= 0.002–0.005 psu`).

## 7. Data and provenance

[Level statistics](tables/level_statistics.csv) contains finite-cell counts, medians
and P10/P90 for each level and product. Its `bias` rows use the **RMSE file's** bias
field, including for the combined product; the combined std file may contain another
reference's bias. The figures are generated by `scripts/04_make_plots.py`.

NetCDF products follow
`expected_{{std,rmse}}_from_{{woa23,aotd,woa23_and_aotd}}_{{potential_temperature,salinity}}_ease2_100km.nc`
in `outputs/expected_spread/` (12 files). Each contains the metric, `*_bias`, `*_n`, `*_n_eff`,
`*_n_years`, coordinates and CRS metadata; regenerate them with `scripts/03_compute_expected_spread.py`.

Source credit: Li, J., Wu, X. & Wang, X. (2025), *An Arctic Ocean Thermohaline Dataset*,
[doi:10.1038/s41597-025-05855-3](https://doi.org/10.1038/s41597-025-05855-3);
AOTD scatter (CC BY 4.0),
[doi:10.6084/m9.figshare.28825490.v1](https://doi.org/10.6084/m9.figshare.28825490.v1);
AOTD climatology,
[doi:10.6084/m9.figshare.27645738.v1](https://doi.org/10.6084/m9.figshare.27645738.v1).
WOA23 temperature and salinity are NOAA NCEI products; full references are in the
[technical note](../../technical_note.md#references).
"""


def write_statistics():
    with (OUT.parent / "tables/level_statistics.csv").open("w", encoding="utf-8") as fh:
        fh.write("variable,metric,product,depth_m,cells_with_data,median,p10,p90\n")
        for v in NAMES:
            for kind in ("rmse", "std", "rmse_declustered", "std_declustered", "rmse_se", "std_se",
                         "rmse_se_pan_arctic", "std_se_pan_arctic", "std_within", "std_between", "bias", "bias_declustered"):
                base_file = "rmse" if kind.startswith(("rmse", "bias")) else "std"
                for p in PRODS:
                    with load(base_file, p, v) as ds:
                        var_name = f"{v}_{kind}"
                        if var_name not in ds:
                            continue
                        f = ds[var_name].values
                        for k, d in enumerate(ds.depth.values):
                            g = f[k][np.isfinite(f[k])]
                            if g.size:
                                fh.write(f"{v},{kind},{p},{d:g},{g.size},{np.median(g):.6g},"
                                         f"{np.percentile(g, 10):.6g},{np.percentile(g, 90):.6g}\n")


def main():
    (OUT.parent / "tables").mkdir(parents=True, exist_ok=True)
    text = report()
    write_statistics()
    OUT.write_text(text, encoding="utf-8")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
