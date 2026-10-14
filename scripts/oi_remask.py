"""Re-apply the NaN mask to saved OI levels with a given relerr_max (no recomputation). Prints coverage per threshold.
Usage: PYTHONPATH=. python scripts/oi_remask.py 2 20 56 [relerr_max]
"""
import sys
import numpy as np
import pandas as pd
from aotd import oi

levels = [int(x) for x in sys.argv[1:] if "." not in x]
cut = [float(x) for x in sys.argv[1:] if "." in x]
cv = pd.read_csv("reports/oi/cv_stats.csv")
for k in levels:
    d = dict(np.load(f"outputs/oi/level_{k}.npz"))
    for name in ("T", "S"):
        infl = max(float(cv[(cv.level == k) & (cv["var"] == name)].z_robust.iloc[0]), 1.0)
        res = dict(field=d[f"{name}_raw_field"], std=d[f"{name}_raw_std"], dmin=d[f"{name}_dmin"], loc_sd=d[f"{name}_loc_sd"])
        okc = np.isfinite(res["field"][0])
        row = []
        for rm in (0.5, 0.75, 1.0, 1.5):
            p = dict(oi.PARAMS, relerr_max=rm)
            f, s = oi.apply_mask(res, infl, p)
            row.append(f"{rm}: {np.isfinite(f[0])[okc].mean():.2f}")
        print(k, name, "infl %.2f" % infl, "Jan valid frac by relerr_max ->", "  ".join(row))
        if cut:
            p = dict(oi.PARAMS, relerr_max=cut[0])
            d[f"{name}_field"], d[f"{name}_std"] = oi.apply_mask(res, infl, p)
    if cut:
        np.savez(f"outputs/oi/level_{k}.npz", **d)
