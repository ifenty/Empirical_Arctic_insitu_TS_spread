import os
from pathlib import Path
import subprocess

import pytest


@pytest.mark.parametrize("optional_inputs", [False, True])
def test_pipeline_quotes_interpreter_and_checks_optional_files(tmp_path, optional_inputs):
    root = Path(__file__).resolve().parents[1]
    runner = tmp_path / "mock python"
    runner.write_text('#!/usr/bin/env bash\nprintf "%s\\n" "$1" >> "$CALL_LOG"\n')
    runner.chmod(0o755)
    script = tmp_path / "run_all.sh"
    script.write_text((root / "run_all.sh").read_text())
    ease = tmp_path / "data" / "ease2"
    ease.mkdir(parents=True)
    if optional_inputs:
        for name in ("NSIDC0772_LatLon_EASE2_N25km_v1.1.nc", "EASE2_N25km.LOCImask_land50_coast0km.720x720.bin"):
            (ease / name).touch()
    log = tmp_path / "calls"
    result = subprocess.run(
        ["bash", str(script)],
        env=dict(os.environ, PYTHON=str(runner), CALL_LOG=str(log)),
        capture_output=True, text=True, check=True,
    )
    calls = log.read_text().splitlines()
    assert len(calls) == (7 if optional_inputs else 6)
    assert calls[0] == "scripts/01_download_data.py"
    assert calls[-1] == ("scripts/07_report_scatter_data.py" if optional_inputs else "scripts/06_report_aotd_vs_woa23.py")
    assert ("Skipping" in result.stdout) is not optional_inputs
