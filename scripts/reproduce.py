"""Run the faithful MATLAB port on all 57 levels and compare with the published monthly climatology.

Outputs (untracked): outputs/repro/month_k.npy ; (tracked) reports/reproduction/level_stats.csv
Usage: PYTHONPATH=. python scripts/reproduce.py [nproc]
"""
import os, sys, time
import numpy as np
import pandas as pd
import xarray as xr
from multiprocessing import Pool

from aotd import upstream_port as up

OUT = "outputs/repro"
REP = "reports/reproduction"
NL = 57


def published():
    """(12, 2, 57, 121, 1440) float32 published monthly T,S (cached)."""
    f = f"{OUT}/pub_month.npy"
    if not os.path.exists(f):
        arr = np.full((12, 2, NL, 121, 1440), np.nan, dtype="float32")
        for m in range(12):
            ds = xr.open_dataset(f"data/climatology/climatology_{m+1:02d}.nc")
            arr[m, 0] = ds["temp"].values
            arr[m, 1] = ds["salt"].values
        np.save(f, arr)
    return np.load(f, mmap_mode="r")


DATA = None


def init():
    global DATA
    d = np.load("data/scatter_all.npz")
    DATA = dict(lat=d["lat"], lon=d["lon"], month=d["month"],
                temp=np.ascontiguousarray(d["temp"][:, :NL]), salt=np.ascontiguousarray(d["salt"][:, :NL]))


def work(k):
    data = DATA
    pub = published()
    ok = ~np.isnan(pub[:, :, k]).all(axis=(0, 1))  # mask workaround: bathymetry file is missing upstream
    ok0 = ~np.isnan(pub[:, :, 0]).all(axis=(0, 1))
    r = up.process_level(k, data, ok, ~ok0)
    np.save(f"{OUT}/month_{k}.npy", r["month"].astype("float32"))
    rows = []
    for m in range(12):
        for v, name in enumerate(("T", "S")):
            a = r["month"][m, v]
            b = pub[m, v, k]
            both = ~np.isnan(a) & ~np.isnan(b)
            dd = a[both] - b[both]
            rows.append(dict(level=k, month=m + 1, var=name, n_pub=int((~np.isnan(b)).sum()), n_port=int((~np.isnan(a)).sum()),
                             rmse=float(np.sqrt(np.mean(dd**2))), bias=float(dd.mean()), maxabs=float(np.abs(dd).max()),
                             corr=float(np.corrcoef(a[both], b[both])[0, 1])))
    return rows


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(REP, exist_ok=True)
    published()
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    levels = [k for k in range(1, NL) if not os.path.exists(f"{OUT}/month_{k}.npy")]
    rows = []
    t = time.time()
    with Pool(min(nproc, 2), initializer=init) as p:
        for r in p.imap_unordered(work, levels):
            rows += r
            print(r[0]["level"], "done", int(time.time() - t), flush=True)
    if rows:
        old = f"{REP}/level_stats.csv"
        df = pd.DataFrame(rows)
        if os.path.exists(old):
            df = pd.concat([pd.read_csv(old), df])
        df.sort_values(["level", "var", "month"]).to_csv(old, index=False)
