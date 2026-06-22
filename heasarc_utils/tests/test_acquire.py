import os

import pytest

from heasarc_utils.lib import heasarc
from heasarc_utils.lib.heasarc_help import acquire
from obsdb_catalog import neutral

# Relative level-2 file each mission's profile expects on disk. Fresh archive
# downloads are gzipped (the live run confirmed this); reduction later gunzips.
LEVEL2_REL = {
    "ixpe": "{obsid}/event_l2/ixpe{obsid}_det1_evt2_v01.fits.gz",
    "nicer": "{obsid}/xti/event_cl/ni{obsid}_0mpu7_cl.evt.gz",
    "nustar": "{obsid}/event_cl/nu{obsid}A01_cl.evt.gz",
}
OBSID = {"ixpe": "01002901", "nicer": "4690010101", "nustar": "30702017002"}


def _write(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(b"x")


def _fake(mission, obsid, *, found=True, create=True, nest=False):
    """Build a fake Heasarc whose download_data lands the level-2 file(s)."""
    class _Result:
        def to_table(self):
            return ["row"] if found else []

    class FakeHeasarc:
        @classmethod
        def query_tap(cls, adql):
            assert obsid in adql  # the ObsID is filtered in ADQL
            return _Result()

        @classmethod
        def locate_data(cls, tab, catalog_name=None):
            return {"access_url": ["u"]}

        @classmethod
        def download_data(cls, links, location=None, host=None):
            if not create:
                return
            root = os.path.join(location, "ftp_prefix") if nest else location
            _write(os.path.join(root, LEVEL2_REL[mission].format(obsid=obsid)))

    return FakeHeasarc


@pytest.mark.parametrize("mission", ["ixpe", "nicer", "nustar"])
def test_acquire_records_level2(tmp_path, monkeypatch, mission):
    obsid = OBSID[mission]
    monkeypatch.setattr(heasarc, "_heasarc", lambda: _fake(mission, obsid))
    summary = acquire(obsid, mission, str(tmp_path))

    assert summary["obsid"] == obsid and summary["mission"] == mission
    entry = neutral.read(tmp_path)["entries"][obsid]
    assert entry["max_level"] == 2 and entry["levels_present"] == [1, 2]
    (prod,) = entry["products"]
    assert prod["kind"] == "level2_evt" and prod["produced_by"] == "heasarc"
    assert prod["path"] == LEVEL2_REL[mission].format(obsid=obsid)

    expected_params = {
        "ixpe": {"detectors": [1, 2, 3], "calibration_version": "v01",
                 "note": heasarc._IXPE_NOTE},
        "nicer": {"mpu": 7},
        "nustar": {"fpm": "A", "note": heasarc._NUSTAR_NOTE},
    }[mission]
    assert prod["params"] == expected_params


def test_reacquire_regresses_to_level2(tmp_path, monkeypatch):
    obsid = OBSID["ixpe"]
    # pre-existing level-3 reduction product
    pcube_rel = f"{obsid}/event_l2/ixpe{obsid}_det1_evt2_v01_rej_src_pcube.fits"
    _write(os.path.join(tmp_path, pcube_rel))
    neutral.advance(tmp_path, "ixpe", obsid, level=3, products=[
        {"kind": "pcube", "path": pcube_rel, "produced_by": "x", "params": {}}])

    monkeypatch.setattr(heasarc, "_heasarc", lambda: _fake("ixpe", obsid))
    acquire(obsid, "ixpe", str(tmp_path))

    entry = neutral.read(tmp_path)["entries"][obsid]
    assert entry["max_level"] == 2                       # regressed
    assert {p["kind"] for p in entry["products"]} == {"pcube", "level2_evt"}


def test_obsid_not_found_raises_and_no_record(tmp_path, monkeypatch):
    monkeypatch.setattr(heasarc, "_heasarc", lambda: _fake("ixpe", "01002901", found=False))
    with pytest.raises(heasarc.HeasarcError):
        acquire("01002901", "ixpe", str(tmp_path))
    assert neutral.read(tmp_path)["entries"] == {}


def test_download_without_level2_raises_and_no_record(tmp_path, monkeypatch):
    monkeypatch.setattr(heasarc, "_heasarc", lambda: _fake("nicer", "4690010101", create=False))
    with pytest.raises(heasarc.HeasarcError):
        acquire("4690010101", "nicer", str(tmp_path))
    assert neutral.read(tmp_path)["entries"] == {}


def test_normalize_layout_relocates_nested(tmp_path, monkeypatch):
    obsid = OBSID["nustar"]
    monkeypatch.setattr(heasarc, "_heasarc", lambda: _fake("nustar", obsid, nest=True))
    acquire(obsid, "nustar", str(tmp_path))
    # relocated to the top level
    assert os.path.isdir(os.path.join(tmp_path, obsid, "event_cl"))
    assert obsid in neutral.read(tmp_path)["entries"]


def test_unregistered_mission_raises(tmp_path):
    with pytest.raises(heasarc.HeasarcError):
        acquire("12345", "chandra", str(tmp_path))


# ── timeout_s tests ───────────────────────────────────────────────────────


def test_timeout_threaded_to_download(tmp_path, monkeypatch):
    """Differential: timeout_s value threads through all three layers to the
    astroquery object at download time.

    A capturing fake reads ``cls.timeout`` inside ``download_data``; the value
    must equal exactly what was passed to ``acquire``.  The default (omitted
    timeout_s) must land as 1800.0, confirming the old astroquery factory
    default (30 s) is never used when this tool is the caller.
    """
    obsid = OBSID["ixpe"]
    mission = "ixpe"
    captured: list[float] = []

    class CapturingFakeHeasarc:
        @classmethod
        def query_tap(cls, adql):
            class _R:
                def to_table(self):
                    return ["row"]
            return _R()

        @classmethod
        def locate_data(cls, tab, catalog_name=None):
            return {"access_url": ["u"]}

        @classmethod
        def download_data(cls, links, location=None, host=None):
            captured.append(cls.timeout)
            _write(os.path.join(location, LEVEL2_REL[mission].format(obsid=obsid)))

    monkeypatch.setattr(heasarc, "_heasarc", lambda: CapturingFakeHeasarc)

    # custom value — must arrive unchanged
    acquire(obsid, mission, str(tmp_path), timeout_s=900.0)
    assert captured[-1] == 900.0, (
        f"expected timeout 900.0 at download_data, got {captured[-1]}"
    )

    # default (omitted) — must be 1800.0, NOT the astroquery factory 30 s
    acquire(obsid, mission, str(tmp_path))
    assert captured[-1] == 1800.0, (
        f"default timeout must be 1800.0, got {captured[-1]}"
    )


def test_timeout_does_not_suppress_obsid_error(tmp_path, monkeypatch):
    """Regression: passing timeout_s must not alter the loud-failure contract
    for a missing ObsID — HeasarcError still raised, nothing recorded."""
    monkeypatch.setattr(heasarc, "_heasarc",
                        lambda: _fake("ixpe", "01002901", found=False))
    with pytest.raises(heasarc.HeasarcError):
        acquire("01002901", "ixpe", str(tmp_path), timeout_s=1800.0)
    assert neutral.read(tmp_path)["entries"] == {}


def test_timeout_default_in_all_three_layers():
    """Regression: default timeout_s is 1800 s at every call layer, confirming
    the old astroquery 30 s bound cannot silently reappear when the argument is
    omitted anywhere in the chain."""
    import inspect
    from heasarc_utils.lib import heasarc as hm

    for fn, name in [
        (acquire, "heasarc_help.acquire"),
        (hm.acquire_obsid, "heasarc.acquire_obsid"),
        (hm.download, "heasarc.download"),
    ]:
        params = inspect.signature(fn).parameters
        assert "timeout_s" in params, f"{name} missing timeout_s parameter"
        assert params["timeout_s"].default == 1800.0, (
            f"{name}: timeout_s default is {params['timeout_s'].default}, "
            f"expected 1800.0"
        )


def test_tool_manifest_exposes_timeout_s():
    """Pin the AI-facing half: TOOL_MANIFEST must expose timeout_s in args_schema.

    A lib/ change that adds timeout_s to heasarc_help.acquire but forgets to
    surface it in skills/__init__.py would leave the LLM unable to set the
    parameter.  This test makes that gap a Layer-B failure.

    Asserts:
    - timeout_s is present in args_schema.properties
    - its type is "number" (JSON Schema — not "float" or "integer")
    - timeout_s is NOT listed in required (optional argument with a default)
    - the description records the 1800 s default so the LLM knows what to expect
    """
    from heasarc_utils.skills import TOOL_MANIFEST

    entry = TOOL_MANIFEST["research_toolkit.heasarc_utils.acquire"]
    schema = entry["args_schema"]
    props = schema["properties"]

    assert "timeout_s" in props, (
        "timeout_s missing from args_schema.properties — "
        "LLM callers cannot pass it"
    )
    assert props["timeout_s"]["type"] == "number", (
        f"timeout_s type should be 'number', got {props['timeout_s'].get('type')!r}"
    )
    assert "timeout_s" not in schema.get("required", []), (
        "timeout_s must not be required — it has a default of 1800 s"
    )
    desc = props["timeout_s"].get("description", "")
    assert "1800" in desc, (
        "timeout_s description must mention the 1800 s default so the LLM "
        f"knows the baseline; got: {desc!r}"
    )


def test_register_new_mission_profile(tmp_path, monkeypatch):
    obsid = "999"

    def _level2(obs_db_path, oid):
        full = os.path.join(obs_db_path, oid, "data", "evt.fits")
        if not os.path.exists(full):
            raise heasarc.HeasarcError("missing")
        return {"kind": "level2_evt", "path": os.path.relpath(full, obs_db_path),
                "produced_by": "heasarc", "params": {}}

    monkeypatch.setitem(heasarc.MISSIONS, "newmis",
                        heasarc.MissionProfile(catalog="newmaster", level2=_level2))

    class _Result:
        def to_table(self):
            return ["row"]

    class FakeHeasarc:
        @classmethod
        def query_tap(cls, adql):
            assert "newmaster" in adql
            return _Result()

        @classmethod
        def locate_data(cls, tab, catalog_name=None):
            return {"x": 1}

        @classmethod
        def download_data(cls, links, location=None, host=None):
            _write(os.path.join(location, obsid, "data", "evt.fits"))

    monkeypatch.setattr(heasarc, "_heasarc", lambda: FakeHeasarc)
    acquire(obsid, "newmis", str(tmp_path))
    assert neutral.read(tmp_path)["entries"][obsid]["max_level"] == 2
