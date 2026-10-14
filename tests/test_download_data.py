import hashlib
import importlib
import io

import pytest

download = importlib.import_module("01_download_data")


def response(body, size=None):
    stream = io.BytesIO(body)
    stream.headers = {} if size is None else {"Content-Length": str(size)}
    return stream


def test_verified_download_replaces_file(tmp_path, monkeypatch):
    target = tmp_path / "input.nc"
    target.write_bytes(b"old")
    body = b"complete data"
    monkeypatch.setattr(download.urllib.request, "urlopen", lambda *a, **k: response(body, len(body)))
    download.fetch("https://example.test/input", target, len(body), hashlib.md5(body).hexdigest())
    assert target.read_bytes() == body
    assert not target.with_suffix(".nc.part").exists()


@pytest.mark.parametrize("body,size,md5,error", [
    (b"short", 10, None, "Incomplete download"),
    (b"wrong", 5, hashlib.md5(b"right").hexdigest(), "MD5 checksum mismatch"),
])
def test_failed_verification_preserves_existing_file(tmp_path, monkeypatch, body, size, md5, error):
    target = tmp_path / "input.nc"
    target.write_bytes(b"previous valid version")
    monkeypatch.setattr(download.urllib.request, "urlopen", lambda *a, **k: response(body, size))
    with pytest.raises(ValueError, match=error):
        download.fetch("https://example.test/input", target, size, md5)
    assert target.read_bytes() == b"previous valid version"


def test_verified_figshare_cache_does_not_download(tmp_path, monkeypatch):
    target = tmp_path / "input.nc"
    target.write_bytes(b"verified")

    def unexpected_request(*args, **kwargs):
        pytest.fail("A verified cache should not make a download request")

    monkeypatch.setattr(download.urllib.request, "urlopen", unexpected_request)
    download.fetch("https://example.test/input", target, 8, hashlib.md5(b"verified").hexdigest())


def test_woa_cache_is_checked_against_server_size(tmp_path, monkeypatch):
    target = tmp_path / "input.nc"
    target.write_bytes(b"truncated")
    body = b"complete WOA file"
    monkeypatch.setattr(download.urllib.request, "urlopen", lambda *a, **k: response(body, len(body)))
    download.fetch("https://example.test/input", target)
    assert target.read_bytes() == body


def test_missing_size_metadata_fails_explicitly(tmp_path, monkeypatch):
    monkeypatch.setattr(download.urllib.request, "urlopen", lambda *a, **k: response(b"data"))
    with pytest.raises(ValueError, match="Cannot verify download size"):
        download.fetch("https://example.test/input", tmp_path / "input.nc")


def test_server_size_must_match_figshare_metadata(tmp_path, monkeypatch):
    monkeypatch.setattr(download.urllib.request, "urlopen", lambda *a, **k: response(b"data", 4))
    with pytest.raises(ValueError, match="metadata mismatch"):
        download.fetch("https://example.test/input", tmp_path / "input.nc", size=5)


def test_figshare_uses_cited_version():
    assert download.FIGSHARE_API.format(download.SCATTER_ID).endswith("/28825490/versions/1")
