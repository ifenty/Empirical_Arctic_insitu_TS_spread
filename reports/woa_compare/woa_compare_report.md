# AOTD (published) vs WOA23 decav91C0

Produced by `scripts/woa_compare.py`. Per-level, per-month numbers: `level_stats.csv`. Depth profiles: `aotd_vs_woa23_vs_depth.png`. Maps (Jan/Mar/Jun/Sep, published | WOA23 | difference) at 10, 100 and 1500 m: `maps/`.

## Method
- WOA23: objectively analysed `t_an` / `s_an`, monthly, 0.25 deg, 57 levels (0-1500 m), 1991-2020. AOTD: 1983-2023.
- Grid: the published grid is labelled 60.0..90.0 N, -180..179.75 E, but it is treated as the WOA23 grid (centres 59.875..89.875 N, -179.875..179.875 E, WOA rows 599-719) and compared cell for cell. The ocean mask agrees best at zero shift (level 0: 2.4% of cells differ; level 20: 2.0%) and worse at every other integer row/column shift tested (+-3). Field RMSE cannot decide this because the products differ for real reasons.
- Valid-cell counts are almost identical (AOTD has 0.4% more cells overall), so the two products have essentially the same footprint.
- Statistics are over cells valid in both products. Differences are AOTD minus WOA23.
- Mar/Sep are split out because they are the two published months the port reproduces; the other ten months are the ones with the BUG-29 defect (see `reports/reproduction/`).

## Results (median of the months in each group)

| depth (m) | T RMSE Mar/Sep | T RMSE other 10 | T bias other 10 | S RMSE Mar/Sep | S RMSE other 10 |
|---|---|---|---|---|---|
| 0 | 1.00 | 1.59 | +0.59 | 1.76 | 1.70 |
| 10 | 0.87 | 1.34 | +0.52 | 1.33 | 1.26 |
| 50 | 0.59 | 0.64 | +0.10 | 0.52 | 0.52 |
| 100 | 0.68 | 0.56 | +0.03 | 0.29 | 0.26 |
| 300 | 0.41 | 0.39 | +0.06 | 0.06 | 0.13 |
| 500 | 0.25 | 0.25 | +0.03 | 0.03 | 0.03 |
| 1000 | 0.11 | 0.11 | +0.01 | 0.02 | **0.61** |
| 1500 | 0.05 | 0.07 | +0.02 | 0.01 | **0.63** |

Units: T degC, S psu.

## Findings
1. **Upper ocean (0-100 m): the products differ by 1-2 degC / 1-2 psu RMSE, with correlation 0.87-0.97.** The difference maps show large-scale, coherent patches over the Siberian shelves, Beaufort/Canada Basin and the Atlantic-side seas, not noise. Part of this is real (1983-2023 vs 1991-2020, mapping methods, ice-covered sampling bias), and the published product has hard polygon edges in March (see 10 m S maps).
2. **Warm bias in the surface layer for ten months.** AOTD minus WOA23 in T is about +0.5-0.6 degC at 0-10 m (Mar/Sep +0.03). This is mainly a seasonal sampling and period effect, not yet attributed.
3. **Below 300 m, temperature agrees well in all months** (RMSE 0.25 degC at 500 m, 0.05-0.07 at 1500 m, r 0.995+).
4. **Deep salinity in ten months is wrong relative to WOA23.** From 800 m down, RMSE jumps from about 0.02 psu to about 0.49 (800 m), 0.61 (1000 m) and 0.63 psu (1500 m), with correlation falling from 0.97 to 0.07-0.12. The break is between 750 and 800 m. March and September stay at 0.01-0.02 psu and r 0.95-0.98. This is independent of our port and confirms BUG-29: the ten non-Mar/Sep monthly files are not consistent with a physically smooth deep salinity field. A smaller anomaly is also present at 300 m (RMSE 0.13 vs 0.06, r 0.79 vs 0.94).

## Caveats
- WOA23 is itself an objective analysis with its own smoothing and gaps, not truth; large differences are not all AOTD errors.
- Different reference periods (1991-2020 vs 1983-2023) and different input data (AOTD adds ITP, Argo, ice-tethered and other Arctic sources).
- The grid-alignment conclusion rests on the ocean mask only. The author should confirm that the published grid is the WOA23 grid.
