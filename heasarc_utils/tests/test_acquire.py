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
