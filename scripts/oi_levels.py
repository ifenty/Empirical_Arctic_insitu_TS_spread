"""Run the OI climatology on a few WOA levels, cross-validate with spatial block holdout, compare with published.

Usage: PYTHONPATH=. python scripts/oi_levels.py 2 20 56
Outputs: outputs/oi/level_k.npz (masked field/std, raw field/std), reports/oi/cv_stats.csv
"""
import os, sys
import numpy as np
import pandas as pd

from aotd import oi
sys.path.insert(0, os.path.dirname(__file__))
from reproduce import published

OUT, REP = "outputs/oi", "reports/oi"
os.makedirs(OUT, exist_ok=True); os.makedirs(REP, exist_ok=True)


def main(levels):
    d = np.load("data/scatter_all.npz")
    lat, lon, month = d["lat"], d["lon"], d["month"]
    flat = oi.profile_nodes(lat, lon)
    pub = published()
    rows = []
    for k in levels:
        ok = ~np.isnan(pub[:, :, k]).all(axis=(0, 1))
        okf = ok.ravel()
        hold = oi.block_holdout(flat)
        res_all = {}
        for v, name in enumerate(("T", "S")):
            val = d["temp" if v == 0 else "salt"][:, k].astype("float64")
            full = oi.run_variable(flat, val, month, okf)
            cv = oi.run_variable(flat, val, month, okf, exclude=hold)
            # cross-validation at withheld nodes (month-specific node means vs month field)
            zs, errs, base = [], [], []
            for mo in range(1, 13):
                idx, m, n = oi.node_table(flat, val, (month == mo) & hold)
                keep = okf[idx] & np.isfinite(cv["field"][mo - 1].ravel()[idx])
                idx, m, n = idx[keep], m[keep], n[keep]
                f = cv["field"][mo - 1].ravel()[idx]; s = cv["std"][mo - 1].ravel()[idx]
                errs.append(m - f); zs.append((m - f) / s)
            e = np.concatenate(errs); z = np.concatenate(zs)
            z = z[np.isfinite(z)]
            z_rob = 1.4826 * np.median(np.abs(z))
            infl = max(z_rob, 1.0)
            fm, sm = oi.apply_mask(full, infl)
            res_all[name] = (full, fm, sm)
            cov = [np.isfinite(fm[mo]).mean() / okf.mean() * 1.0 for mo in range(12)]
            pm = pub[:, v, k]
            for mo in (2, 8):  # March, September: published months that match the released code
                a, b = fm[mo], pm[mo]
                bt = np.isfinite(a) & np.isfinite(b)
                dd = a[bt] - b[bt]
                rows.append(dict(level=k, var=name, month=mo + 1, n_cv=len(e), cv_rmse=float(np.sqrt(np.mean(e**2))),
                                 cv_bias=float(e.mean()), z_std=float(np.std(z)), z_robust=float(z_rob), frac_valid=float(np.isfinite(a).sum() / okf.sum()),
                                 pub_valid=float(np.isfinite(b).sum() / okf.sum()), rmse_vs_pub=float(np.sqrt(np.mean(dd**2))),
                                 corr_vs_pub=float(np.corrcoef(a[bt], b[bt])[0, 1])))
            print(k, name, "CV rmse %.3f bias %.3f z-std %.2f robust %.2f" % (np.sqrt(np.mean(e**2)), e.mean(), np.std(z), z_rob),
                  "mean valid frac %.2f" % np.mean(cov), flush=True)
        pd.DataFrame([r for r in rows if r["level"] == k]).to_csv(f"{REP}/cv_stats_{k}.csv", index=False)
        np.savez(f"{OUT}/level_{k}.npz",
                 **{f"{n}_{q}": a for n, (full, fm, sm) in res_all.items()
                    for q, a in (("field", fm), ("std", sm), ("raw_field", full["field"]), ("raw_std", full["std"]), ("dmin", full["dmin"]), ("loc_sd", full["loc_sd"]))})


if __name__ == "__main__":
    main([int(x) for x in sys.argv[1:]])
