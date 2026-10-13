# Copilot instructions: AOTD_v2

Python port, reproduction and improved rebuild of the Arctic Ocean Thermohaline Dataset (AOTD, Sci. Data 2025), compared against WOA23. See `README.md` (layout, data provenance), `STATUS.md` (newest-first progress log) and `HANDOFF.md` (current state, next steps).

## Environment
- Use the project venv: `.venv/bin/python`. The conda `python` on PATH has no numpy.
- Run from the repo root with `PYTHONPATH=.`, e.g. `PYTHONPATH=. .venv/bin/python scripts/reproduce.py 2`.
- Tracked: code, reports, figures, docs. Untracked (see `.gitignore`): `data/`, `outputs/`, `*.npy`, `*.npz`, `*.nc`.

## Resource limits (hard rules)
The host is a shared WSL VM (28 GB RAM) that also runs other sessions. A too-heavy run has crashed the VM before.
- Use at most 2 worker processes.
- Use at most 20% of total RAM (about 5.6 GB) for any job.
- Set `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`.
- Enforce the memory cap with a cgroup, not just good intentions:
  ```
  systemd-run --user --unit=<name> --collect -p MemoryMax=5600M -p MemorySwapMax=0 \
    --working-directory=$PWD --setenv=PYTHONPATH=. --setenv=OMP_NUM_THREADS=1 \
    -p StandardOutput=file:$PWD/outputs/repro/log.txt -p StandardError=file:$PWD/outputs/repro/log.txt \
    $PWD/.venv/bin/python scripts/reproduce.py 2
  ```
  Check status with `systemctl --user status <name>` and `journalctl --user -u <name>`.
- Measure before scaling up: one level of `process_level` took about 130 s and about 1.9 GB peak RSS.

## Conventions
- Faithful port code (`aotd/upstream_port.py`) must mirror the MATLAB in `upstream_matlab/`, including its quirks. Fixes go in the improved method, not the port.
- Log every MATLAB bug found in `docs/upstream_bug_log.md`.
- Update `STATUS.md` (newest entry on top) after each milestone.
- Land mask note: in the LOCI mask, 255 is ocean and 0 is land.
- Do not commit data or large outputs.
