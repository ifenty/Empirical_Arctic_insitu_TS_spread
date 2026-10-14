"""Download all input data needed by the pipeline into data/ (untracked).

  data/scatter/AOTD{1983..2023}.nc      AOTD scattered profiles, figshare doi:10.6084/m9.figshare.28825490 (Li et al., CC BY 4.0), ~90 MB
  data/aotd_climatology/climatology_{01..12}.nc  AOTD gridded monthly climatology, figshare doi:10.6084/m9.figshare.27645738 (CC BY 4.0), ~2 GB
  data/woa23_raw/woa23_decav91C0_{t,s}{01..16}_04.nc  WOA23 1991-2020 0.25 deg T and S: 12 monthly + 4 seasonal files, ~18 GB, NOAA NCEI

Downloads are checked before replacing local files. Existing files of the right size
(and, for figshare, matching MD5) are skipped; interrupted files are downloaded again.

Usage (from the repo root):
    python scripts/01_download_data.py                # everything
    python scripts/01_download_data.py --no-climatology   # skip the AOTD climatology (NOT possible if you want script 03, which needs it for the from_aotd and combined products)
    python scripts/01_download_data.py --only scatter     # or: climatology, woa23
"""
import argparse
import hashlib
import json
import os
import urllib.request

FIGSHARE_API = "https://api.figshare.com/v2/articles/{}/versions/1"
SCATTER_ID, CLIM_ID = 28825490, 27645738
WOA_URL = "https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/{var}/netcdf/decav91C0/0.25/woa23_decav91C0_{v}{n:02d}_04.nc"


def file_md5(path):
    checksum = hashlib.md5()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 22):
            checksum.update(chunk)
    return checksum.hexdigest()


def fetch(url, dest, size=None, md5=None):
    dest = os.fspath(dest)
    if os.path.exists(dest) and (size is None or os.path.getsize(dest) == size):
        if (size is not None or md5 is not None) and (md5 is None or file_md5(dest) == md5):
            print("have", dest)
            return
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    tmp = dest + ".part"
    with urllib.request.urlopen(url, timeout=120) as r:
        content_length = r.headers.get("Content-Length")
        response_size = int(content_length) if content_length is not None else None
        if size is None:
            if response_size is None:
                raise ValueError(f"Cannot verify download size for {url}: no Content-Length or expected size")
            size = response_size
        elif response_size is not None and response_size != size:
            raise ValueError(f"Download size metadata mismatch for {url}: expected {size}, server reports {response_size}")
        if os.path.exists(dest) and os.path.getsize(dest) == size and (md5 is None or file_md5(dest) == md5):
            print("have", dest)
            return
        print("get ", url, "->", dest, flush=True)
        with open(tmp, "wb") as f:
            while chunk := r.read(1 << 22):
                f.write(chunk)
    actual_size = os.path.getsize(tmp)
    if actual_size != size:
        raise ValueError(f"Incomplete download for {dest}: expected {size} bytes, received {actual_size}")
    if md5 is not None and file_md5(tmp) != md5:
        raise ValueError(f"MD5 checksum mismatch for {dest}")
    os.replace(tmp, dest)


def figshare(article_id, outdir):
    with urllib.request.urlopen(FIGSHARE_API.format(article_id), timeout=120) as r:
        meta = json.load(r)
    for f in meta["files"]:
        fetch(f["download_url"], os.path.join(outdir, f["name"]), f["size"], f.get("computed_md5"))


def woa23(outdir="data/woa23_raw"):
    for v, var in (("t", "temperature"), ("s", "salinity")):
        for n in range(1, 17):  # 01-12 monthly (57 levels, 0-1500 m), 13-16 seasonal (102 levels); 00 (annual) is not used
            fetch(WOA_URL.format(var=var, v=v, n=n), f"{outdir}/woa23_decav91C0_{v}{n:02d}_04.nc")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=("scatter", "climatology", "woa23"))
    ap.add_argument("--no-climatology", action="store_true")
    a = ap.parse_args()
    if a.only in (None, "scatter"):
        figshare(SCATTER_ID, "data/scatter")
    if a.only == "climatology" or (a.only is None and not a.no_climatology):
        figshare(CLIM_ID, "data/aotd_climatology")
    if a.only in (None, "woa23"):
        woa23()
    print("done")
