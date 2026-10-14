"""Build reports/reproduction/reproduction_report.md (+ figure) from level_stats.csv.

Usage: PYTHONPATH=. .venv/bin/python scripts/reproduction_report.py
"""
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REP = "reports/reproduction"
d = np.load("data/scatter_all.npz")
depth = d["depth"]
df = pd.read_csv(f"{REP}/level_stats.csv")
df["depth_m"] = depth[df["level"].values]
df["count_diff"] = df["n_port"] - df["n_pub"]

g = df.groupby(["var", "level", "depth_m"]).agg(
    rmse=("rmse", "mean"), rmse_max_month=("rmse", "max"), bias=("bias", "mean"),
    maxabs=("maxabs", "max"), corr=("corr", "min"), count_diff=("count_diff", lambda x: int(np.abs(x).max())),
).reset_index()

fig, ax = plt.subplots(1, 3, figsize=(12, 6), sharey=True)
for v, c in (("T", "C0"), ("S", "C1")):
    x = g[g["var"] == v]
    ax[0].plot(x.rmse, x.depth_m, c, label=v)
    ax[1].plot(x.bias, x.depth_m, c, label=v)
    ax[2].plot(x["corr"], x.depth_m, c, label=v)
ax[0].set_xlabel("monthly-mean RMSE (degC for T, psu for S)")
ax[1].set_xlabel("monthly-mean bias (port - published)")
ax[2].set_xlabel("minimum monthly correlation")
ax[0].invert_yaxis()
ax[0].set_ylabel("depth (m)")
for a in ax:
    a.grid(alpha=.3)
    a.legend()
fig.tight_layout()
fig.savefig(f"{REP}/reproduction_vs_depth.png", dpi=120)


def table(v):
    x = g[g["var"] == v][["level", "depth_m", "rmse", "rmse_max_month", "bias", "maxabs", "corr", "count_diff"]]
    return x.round(3).to_markdown(index=False)


lines = [
    "# Reproduction of the published AOTD monthly climatology (faithful MATLAB port)",
    "",
    "Compares `aotd/upstream_port.py` output with the published monthly climatology, per level, over all valid cells",
    "in both products. Source: `reports/reproduction/level_stats.csv` (12 months x T/S x levels).",
    f"Levels covered: {df.level.min()}-{df.level.max()} ({df.level.nunique()} levels; level 0 and any level missing from the CSV were not scored).",
    "",
    "## Summary",
]
for v, name, unit in (("T", "Temperature", "degC"), ("S", "Salinity", "psu")):
    x = g[g["var"] == v]
    w = x.loc[x.rmse.idxmax()]
    lines.append(
        f"- **{name}**: RMSE median {x.rmse.median():.3f} {unit}, worst {w.rmse:.3f} {unit} at level {int(w.level)} "
        f"({w.depth_m:.0f} m); minimum monthly correlation {x['corr'].min():.3f}; "
        f"levels with valid-cell count mismatch: {int((x.count_diff != 0).sum())} of {len(x)} "
        f"(largest mismatch {int(x.count_diff.max())} cells)."
    )
lines += [
    "",
    "Conclusion: the port reproduces the published fields closely (valid-cell counts match exactly for T) but not",
    "bit-for-bit. The residual is largest in the upper ocean. Known candidate causes are a different data/QC",
    "version and the missing bathymetry file upstream (see `docs/upstream_bug_log.md`, `STATUS.md`).",
    "",
    "![RMSE, bias and correlation vs depth](reproduction_vs_depth.png)",
    "",
    "Columns: rmse = mean over months; rmse_max_month = worst month; maxabs = largest absolute difference;",
    "corr = minimum over months; count_diff = largest |n_port - n_pub| valid cells.",
    "",
    "## Temperature, per level",
    "",
    table("T"),
    "",
    "## Salinity, per level",
    "",
    table("S"),
    "",
]
open(f"{REP}/reproduction_report.md", "w").write("\n".join(lines))
print("ok")
