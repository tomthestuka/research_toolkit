"""acquire against a fake archive (the ``_heasarc`` seam) and a throwaway store (root conftest
``store``): what a download stores, when it is answered from the store, and every refusal."""

import gzip
import os
from pathlib import PurePosixPath

import pytest

import storelib
from heasarc_utils.lib import heasarc
from heasarc_utils.lib.heasarc_help import ALGORITHM_VERSION, KIND, TOOL, acquire

# The level-2 file each mission's download holds, below its ObsID folder. Fresh archive downloads are
# gzipped (the live run of 2026-06-17 confirmed this); acquire gunzips them, so the stored file is the
# plain twin (``_plain``).
LEVEL2_REL = {
    "ixpe": "event_l2/ixpe{obsid}_det1_evt2_v01.fits.gz",
    "nicer": "xti/event_cl/ni{obsid}_0mpu7_cl.evt.gz",
    "nustar": "event_cl/nu{obsid}A01_cl.evt.gz",
}
OBSID = {"ixpe": "01002901", "nicer": "4690010101", "nustar": "30702017002"}
# Where HEASARC's archive keeps each ObsID below the mission (the access_url locate_data returned for
# these ObsIDs, 2026-11-04).
TREE = {"ixpe": "data/obs/01", "nicer": "data/obs/2022_02", "nustar": "data/obs/07/3"}
IXPE_LOCKED_REL = "event_l2/ixpe{obsid}_det1_evt2_v03.fits.gz" + heasarc.LOCK_SUFFIX


def _url(mission, obsid, tree=None):
    """An access_url shaped as HEASARC's (its double slash included)."""
    return f"https://heasarc.gsfc.nasa.gov/FTP/{mission}/{tree or TREE[mission]}//{obsid}/"


def _plain(rel: str) -> str:
    return rel[:-3] if rel.endswith(".gz") else rel


#: The processing the fake master-table row names, and its level-2 files' headers state by default.
SOFTVER, CALDBVER = "Hea_test_V6.36", "20251101"


def _fits(softver=SOFTVER, caldbver=CALDBVER) -> bytes:
    """A minimal FITS file: a primary header carrying the processing keywords a level-2 event states."""
    cards = ["SIMPLE  =                    T", "BITPIX  =                    8", "NAXIS   =                    0",
             f"SOFTVER = '{softver}'", f"CALDBVER= '{caldbver}'", "END"]
    return "".join(c.ljust(80) for c in cards).ljust(2880).encode("ascii")


def _write(path, data=None, *, softver=SOFTVER):
    """A fixture file: an event file (``.fits``/``.evt``) is a FITS file stating ``softver``, anything
    else one byte; a ``.gz`` path gets valid gzip content (as an archive download has)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    plain = path[:-3] if path.endswith(".gz") else path
    if data is None:
        data = _fits(softver) if plain.endswith((".fits", ".evt")) else b"x"
    with (gzip.open(path, "wb") if path.endswith(".gz") else open(path, "wb")) as fh:
        fh.write(data)


def _row(processing_date, public_date=None):
    """A master-table row as astroquery's table answers ``table[column][0]``."""
    return {"processing_date": [processing_date], "software_version": [SOFTVER],
            "caldb_version": [CALDBVER], "public_date": [public_date if public_date is not None else 60000.0]}


def _fake(mission, obsid, *, processing_date, found=True, files=None, folder=None, public_date=None,
          downloads=None, tree=None, softver=SOFTVER):
    """A fake Heasarc: its row carries ``processing_date`` (unique per test, so each test's recipe is
    its own), its download lands ``files`` (default: the mission's level-2 event) below the ObsID."""
    class _Result:
        def to_table(self):
            return _row(processing_date, public_date) if found else []

    class FakeHeasarc:
        @classmethod
        def query_tap(cls, adql):
            assert obsid in adql  # the ObsID is filtered in ADQL
            return _Result()

        @classmethod
        def locate_data(cls, tab, catalog_name=None):
            return {"access_url": [_url(folder or mission, obsid, tree=tree or TREE.get(mission, "data/obs"))]}

        @classmethod
        def download_data(cls, links, location=None, host=None):
            if downloads is not None:
                downloads.append(cls.timeout)
            for rel in (files if files is not None else [LEVEL2_REL[mission]]):
                # astroquery unpacks the observation as <location>/<obsid>/ (every host)
                _write(os.path.join(location, obsid, rel.format(obsid=obsid)), softver=softver)

    return FakeHeasarc


@pytest.fixture
def archive(store, tmp_path, monkeypatch):
    """Install a fake archive: ``archive(mission, obsid, **fake)``; the working folder is the test's."""
    monkeypatch.chdir(tmp_path)
    date = 61000 + int(store["tag"][:6], 16) / 2 ** 24  # a processing date of this test alone

    def install(mission, obsid, **fake):
        fake.setdefault("processing_date", date)
        monkeypatch.setattr(heasarc, "_heasarc", lambda: _fake(mission, obsid, **fake))
    return install


def _products():
    """The catalog's product rows of kind raw, oldest first."""
    with storelib.connect() as conn:
        return conn.execute("SELECT state, reason, recipe FROM catalog.product WHERE kind = 'raw' "
                            "ORDER BY claimed_at").fetchall()


@pytest.mark.parametrize("mission", ["ixpe", "nicer", "nustar"])
def test_acquire_stores_the_observation_as_a_raw_product(archive, tmp_path, mission):
    obsid = OBSID[mission]
    archive(mission, obsid, files=[LEVEL2_REL[mission], "auxil/extra.txt"])
    summary = acquire(obsid, mission)

    record = summary["record"]
    assert summary["downloaded"] is True and summary["locked"] is False and summary["files"] == 2
    assert record["store"] == storelib.store_name() and record["bucket"] == mission
    assert PurePosixPath(record["key"]).parts[:2] == (KIND, obsid) and record["key"].endswith("/recipe.json")
    level2 = _plain(LEVEL2_REL[mission].format(obsid=obsid))
    assert summary["level2_evt"]["key"] == f"{storelib.product_folder(record)}/{level2}"
    [folder] = storelib.stage_in(record, tmp_path / "in")
    assert (folder / level2).read_bytes() == _fits(), "the stored event is the gunzipped file"
    assert sorted(str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file()) == sorted(
        [level2, "auxil/extra.txt", "recipe.json"])
    assert summary["archive"]["software_version"] == "Hea_test_V6.36" and summary["archive"]["public"] is True


def test_a_repeat_acquire_is_answered_from_the_store_and_a_reprocessing_downloads_afresh(archive):
    obsid, downloads = OBSID["ixpe"], []
    archive("ixpe", obsid, downloads=downloads)
    first = acquire(obsid, "ixpe")
    again = acquire(obsid, "ixpe")
    assert len(downloads) == 1 and again["downloaded"] is False and again["record"] == first["record"]

    archive("ixpe", obsid, downloads=downloads, processing_date=first["archive"]["processing_date"] + 1.5)
    later = acquire(obsid, "ixpe")
    assert len(downloads) == 2 and later["downloaded"] is True and later["record"] != first["record"]


def test_an_obsid_the_archive_lacks_raises_before_any_claim(archive):
    archive("ixpe", "01002901", found=False)
    before = len(_products())
    with pytest.raises(heasarc.HeasarcError, match="not found"):
        acquire("01002901", "ixpe")
    assert len(_products()) == before


def test_a_download_without_its_level2_event_stores_nothing_and_says_why(archive):
    archive("nicer", OBSID["nicer"], files=["auxil/other.txt"])
    with pytest.raises(heasarc.HeasarcError, match="NICER level-2 cleaned event"):
        acquire(OBSID["nicer"], "nicer")
    state, reason, _ = _products()[-1]
    assert state == "failed" and "NICER level-2 cleaned event" in reason


def test_a_download_older_than_the_processing_its_row_names_stores_nothing(archive):
    """The archive's catalog row can name a newer processing than the files it still serves (seen for
    IXPE 01250101: row 1.11, files 1.5): the download fails, naming both, and nothing is stored."""
    obsid = OBSID["ixpe"]
    archive("ixpe", obsid, softver="Hea_10Mar2022_V6.30_IXPE_21JUNE2022_1.5")
    with pytest.raises(heasarc.HeasarcError) as raised:
        acquire(obsid, "ixpe")
    assert "Hea_10Mar2022_V6.30_IXPE_21JUNE2022_1.5" in str(raised.value) and SOFTVER in str(raised.value)
    assert "lags its catalog" in str(raised.value)
    state, reason, recipe = _products()[-1]
    assert state == "failed" and "lags its catalog" in reason and recipe["settings"]["software_version"] == SOFTVER
    with storelib.connect() as conn:
        assert conn.execute("SELECT count(*) FROM catalog.object o JOIN catalog.product p USING (recipe_hash) "
                            "WHERE p.state = 'failed' AND p.subject = %s", (obsid,)).fetchone()[0] == 0


def test_the_file_state_reads_softver_and_caldbver_from_the_header(tmp_path):
    path = tmp_path / "evt.fits"
    path.write_bytes(_fits("Hea_X_1.11", "20260101"))
    assert heasarc.file_state(str(path)) == {"SOFTVER": "Hea_X_1.11", "CALDBVER": "20260101"}
    state = {"software_version": "Hea_X_1.11", "caldb_version": "20260101"}
    assert heasarc.state_mismatch(state, str(path)) is None
    assert "CALDBVER '20260101'" in heasarc.state_mismatch({**state, "caldb_version": "20270101"}, str(path))
    assert heasarc.state_mismatch({"software_version": None, "caldb_version": None}, str(path)) is None


def test_an_archive_folder_other_than_the_mission_is_refused(archive):
    archive("ixpe", "01002901", folder="ixpe_other")
    with pytest.raises(heasarc.HeasarcError, match="register the mission under its archive name"):
        acquire("01002901", "ixpe")


def test_unregistered_mission_raises():
    with pytest.raises(heasarc.HeasarcError, match="no MissionProfile"):
        acquire("123", "nosuchmission")


def test_timeout_threads_to_the_download(archive):
    """The value acquire is given reaches astroquery's ``timeout`` at download time."""
    captured = []
    archive("ixpe", OBSID["ixpe"], downloads=captured)
    acquire(OBSID["ixpe"], "ixpe", timeout_s=900.0)
    assert captured == [900.0]


def test_timeout_default_in_every_layer():
    """Omitted, the timeout is 1800 s at every layer, never astroquery's factory 30 s."""
    import inspect

    for fn in (acquire, heasarc.download_obsid, heasarc.download):
        assert inspect.signature(fn).parameters["timeout_s"].default == 1800.0, fn.__name__


def test_the_manifest_names_the_tool_its_version_and_its_store():
    from heasarc_utils.skills import TOOL_MANIFEST

    entry = TOOL_MANIFEST[TOOL]
    assert entry["writes"] == ["db"] and entry["algorithm_version"] == ALGORITHM_VERSION
    schema = entry["args_schema"]
    assert schema["properties"]["timeout_s"]["type"] == "number"
    assert "1800" in schema["properties"]["timeout_s"]["description"]
    assert schema["required"] == ["obsid", "mission"]


# ── encrypted (proprietary-period) observations ───────────────────────────
# HEASARC serves an observation still inside its proprietary period GPG-encrypted. That download is
# complete and correct; it simply cannot be read. It is stored, and asking again while it is still
# proprietary is the same recipe, so nothing is transferred again.


def test_an_encrypted_download_is_stored_locked_and_not_fetched_again_until_public(archive):
    obsid, downloads = "04004001", []
    archive("ixpe", obsid, files=[IXPE_LOCKED_REL], tree="data/obs/04",
            public_date=heasarc.mjd_now() + 72, downloads=downloads)
    first = acquire(obsid, "ixpe")
    assert first["locked"] is True and first["downloaded"] is True and "proprietary" in first["note"]
    assert first["level2_evt"]["key"].endswith(".fits.gz.gpg") and 71 <= first["days_until_public"] <= 73

    again = acquire(obsid, "ixpe")
    assert len(downloads) == 1 and again["downloaded"] is False and again["locked"] is True

    # the proprietary period has lapsed: the archive serves it in the clear
    archive("ixpe", obsid, tree="data/obs/04", public_date=heasarc.mjd_now() - 1, downloads=downloads,
            processing_date=first["archive"]["processing_date"])
    public = acquire(obsid, "ixpe")
    assert len(downloads) == 2 and public["downloaded"] is True and public["locked"] is False


def test_a_public_observation_still_served_encrypted_stores_nothing(archive):
    obsid = "04004002"
    archive("ixpe", obsid, files=[IXPE_LOCKED_REL], tree="data/obs/04", public_date=heasarc.mjd_now() - 3)
    with pytest.raises(heasarc.HeasarcError, match="still serves it GPG-encrypted"):
        acquire(obsid, "ixpe")
    state, reason, recipe = _products()[-1]
    assert state == "failed" and recipe["subject"] == obsid and "GPG-encrypted" in reason


def test_an_archive_row_without_a_processing_date_is_refused_before_any_claim(archive):
    archive("ixpe", "01002901", processing_date=None)
    before = len(_products())
    with pytest.raises(heasarc.HeasarcError, match="no processing_date"):
        acquire("01002901", "ixpe")
    assert len(_products()) == before


def test_a_decrypted_copy_wins_over_its_encrypted_twin(archive):
    obsid = "04004001"
    archive("ixpe", obsid, tree="data/obs/04",
            files=[IXPE_LOCKED_REL, LEVEL2_REL["ixpe"].replace("v01", "v03")])
    summary = acquire(obsid, "ixpe")
    assert summary["locked"] is False and not heasarc.is_locked(summary["level2_evt"]["key"])


def test_a_mission_profile_is_all_a_mission_needs(archive, monkeypatch):
    obsid = "999"
    monkeypatch.setitem(heasarc.MISSIONS, "ixpe", heasarc.MissionProfile(
        catalog="newmaster", level2=lambda paths, oid: heasarc._one(paths, oid, "data", "evt.fits", "event")))
    archive("ixpe", obsid, files=["data/evt.fits"], tree="data/obs")
    summary = acquire(obsid, "ixpe")
    assert summary["level2_evt"]["key"].endswith("/data/evt.fits")


# ── the archive's answers (no store) ─────────────────────────────────────


def test_the_archives_own_message_names_an_obsid_without_data():
    links = {"access_url": [""], "error_message": ["NotFoundFault: No products found for id 64"]}
    with pytest.raises(heasarc.HeasarcError, match="No products found"):
        heasarc.archive_path(links, "60464001")
    two = {"access_url": [_url("ixpe", "01002901"), _url("ixpe", "01002901", tree="data/obs/02")]}
    with pytest.raises(heasarc.HeasarcError, match="2 places"):
        heasarc.archive_path(two, "01002901")


def test_the_archive_state_reads_the_row_and_leaves_blanks_empty():
    state = heasarc.archive_state({"processing_date": [61000.25], "software_version": [b"Hea_V6"],
                                   "caldb_version": ["  "], "public_date": [float("nan")]})
    assert state == {"processing_date": 61000.25, "software_version": "Hea_V6", "caldb_version": None,
                     "public_date": None}
    assert heasarc.archive_state([]) == dict.fromkeys(state)


@pytest.mark.parametrize("mission, paths, found", [
    ("ixpe", ["event_l2/ixpe01002901_det1_evt2_v01_rej.fits", "event_l2/ixpe01002901_det1_evt2_v02.fits"],
     "event_l2/ixpe01002901_det1_evt2_v02.fits"),
    ("nicer", ["xti/event_cl/ni4690010101_0mpu7_cl.evt.gpg", "xti/event_cl/ni4690010101_0mpu7_cl.evt"],
     "xti/event_cl/ni4690010101_0mpu7_cl.evt"),
    ("nustar", ["event_cl/sub/nu30702017002A01_cl.evt", "event_cl/nu30702017002A01_cl.evt.gz"],
     "event_cl/nu30702017002A01_cl.evt.gz"),
])
def test_each_mission_names_its_level2_event_among_the_files(mission, paths, found):
    assert heasarc.MISSIONS[mission].level2(paths, OBSID[mission]) == found
