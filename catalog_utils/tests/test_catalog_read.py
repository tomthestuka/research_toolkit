import json

import pytest

from catalog_utils.lib.catalog_read import data_catalog
from obsdb_catalog import neutral


def _write_catalog(mission_dir, entries):
    mission_dir.mkdir(parents=True, exist_ok=True)
    (mission_dir / "catalog.json").write_text(json.dumps({"entries": entries}), encoding="utf-8")


def test_data_catalog_lists_obsids_levels_products(tmp_path):
    d = tmp_path / "IXPE"  # basename → mission 'ixpe'
    _write_catalog(d, {
        "01002901": {"max_level": 3, "products": [{"kind": "level2_evt"}, {"kind": "pcube"}]},
        "05250399": {"max_level": 2, "products": [{"kind": "level2_evt"}]},
    })
    out = data_catalog(str(d))
    assert out.startswith("Local ixpe data (2 ObsIDs):")
    assert "- 01002901: max_level=3, products=[level2_evt, pcube]" in out
    assert "- 05250399: max_level=2, products=[level2_evt]" in out


def test_data_catalog_min_level_filters(tmp_path):
    d = tmp_path / "IXPE"
    _write_catalog(d, {
        "01002901": {"max_level": 3, "products": [{"kind": "pcube"}]},
        "05250399": {"max_level": 2, "products": [{"kind": "level2_evt"}]},
    })
    out = data_catalog(str(d), min_level=3)
    assert out.startswith("Local ixpe data (1 ObsIDs):")
    assert "01002901" in out and "05250399" not in out


def test_data_catalog_absent_is_empty(tmp_path):
    d = tmp_path / "NICER"
    d.mkdir()
    out = data_catalog(str(d))
    assert out.startswith("Local nicer data (0 ObsIDs):")
    assert "(no local data catalogued)" in out


def test_data_catalog_malformed_raises(tmp_path):
    d = tmp_path / "IXPE"
    d.mkdir()
    (d / "catalog.json").write_text("{not json", encoding="utf-8")
    with pytest.raises(neutral.CatalogError):
        data_catalog(str(d))


def test_data_catalog_flags_locked_encrypted_after_recording(tmp_path):
    """A recorded product whose file is gone but whose .gpg twin is present → locked.

    This is local encryption of an already-catalogued file: the readable path the
    catalog knows has been replaced by its encrypted twin.
    """
    d = tmp_path / "IXPE"
    d.mkdir()
    _write_catalog(d, {
        "05250399": {"max_level": 2, "products": [
            {"kind": "level2_evt",
             "path": "05250399/event_l2/ixpe05250399_det1_evt2_v01.fits.gz"}]},
    })
    evt = d / "05250399" / "event_l2"
    evt.mkdir(parents=True)
    (evt / "ixpe05250399_det1_evt2_v01.fits.gz.gpg").write_bytes(b"x")  # only the encrypted copy
    out = data_catalog(str(d))
    assert "05250399" in out
    assert "LOCKED, NOT USABLE" in out and "GPG-encrypted" in out


def test_data_catalog_flags_locked_when_archive_served_encrypted(tmp_path):
    """A product recorded *as* its .gpg path → locked, and still levelled honestly.

    This is the proprietary-period case the archive produces: the level-2 event is
    genuinely there and genuinely level 2, so the level stays 2 and the mark — not a
    demoted level — is what says it cannot be read. It therefore survives a
    ``min_level=2`` filter, which is exactly why the mark has to be on the line.
    """
    d = tmp_path / "IXPE"
    d.mkdir()
    rel = "04004001/event_l2/ixpe04004001_det1_evt2_v03.fits.gz.gpg"
    _write_catalog(d, {
        "04004001": {"max_level": 2, "products": [{"kind": "level2_evt", "path": rel}]},
    })
    evt = d / "04004001" / "event_l2"
    evt.mkdir(parents=True)
    (evt / "ixpe04004001_det1_evt2_v03.fits.gz.gpg").write_bytes(b"x")

    out = data_catalog(str(d), min_level=2)
    assert "- 04004001: max_level=2" in out
    assert "LOCKED, NOT USABLE" in out


def test_data_catalog_not_locked_when_usable_product_present(tmp_path):
    """A present, usable product means NOT locked, even if a .gpg twin also exists."""
    d = tmp_path / "IXPE"
    d.mkdir()
    _write_catalog(d, {
        "01002901": {"max_level": 2, "products": [
            {"kind": "level2_evt", "path": "01002901/event_l2/evt2.fits"}]},
    })
    evt = d / "01002901" / "event_l2"
    evt.mkdir(parents=True)
    (evt / "evt2.fits").write_bytes(b"x")       # decrypted, usable
    (evt / "evt2.fits.gpg").write_bytes(b"x")   # an encrypted twin is also around
    out = data_catalog(str(d))
    assert "01002901" in out and "LOCKED" not in out


def test_data_catalog_marks_partially_locked_entry(tmp_path):
    """A readable product beside an encrypted one is reported as both, not collapsed."""
    d = tmp_path / "IXPE"
    d.mkdir()
    _write_catalog(d, {
        "01002901": {"max_level": 3, "products": [
            {"kind": "pcube", "path": "01002901/event_l2/pcube.fits"},
            {"kind": "level2_evt", "path": "01002901/event_l2/evt2.fits.gpg"},
        ]},
    })
    evt = d / "01002901" / "event_l2"
    evt.mkdir(parents=True)
    (evt / "pcube.fits").write_bytes(b"x")
    (evt / "evt2.fits.gpg").write_bytes(b"x")

    out = data_catalog(str(d))
    assert "LOCKED, NOT USABLE (some products)" in out


def _header(cards):
    return "".join(c.ljust(80) for c in [*cards, "END"]).ljust(2880).encode("ascii", "replace")


def _fits_with_object(path, name=None, *, extension=True, extra=(), opener=open):
    """A small IXPE-shaped event file: an empty primary, then an EVENTS header carrying OBJECT."""
    obj = [f"OBJECT  = '{name}'"] if name is not None else []
    primary = ["SIMPLE  =                    T", "BITPIX  =                    8", "NAXIS   =                    0"]
    events = ["XTENSION= 'BINTABLE'", "BITPIX  =                    8", "NAXIS   =                    2",
              "NAXIS1  =                    0", "NAXIS2  =                    0", "EXTNAME = 'EVENTS'", *obj, *extra]
    path.parent.mkdir(parents=True, exist_ok=True)
    with opener(path, "wb") as fh:
        fh.write(_header(primary) + _header(events) if extension else _header([*primary, *obj, *extra]))


def _listing(tmp_path, rel, **fits):
    d = tmp_path / "IXPE"
    _fits_with_object(d / rel, **fits)
    _write_catalog(d, {"04250201": {"max_level": 2, "products": [{"kind": "level2_evt", "path": rel}]}})
    return next(l for l in data_catalog(str(d)).splitlines() if "04250201" in l)


def test_data_catalog_names_the_source_its_event_file_records(tmp_path):
    assert _listing(tmp_path, "04250201/e.fits", name="IGR J17091-3624").startswith(
        "- 04250201 (IGR J17091-3624): max_level=2")
    assert _listing(tmp_path, "04250201/p.fits", name="GX 339-4", extension=False).startswith("- 04250201 (GX 339-4):")
    assert _listing(tmp_path, "04250201/q.fits", name="O''Brien X-1").startswith("- 04250201 (O'Brien X-1):")


def test_a_header_without_a_name_shows_its_coordinates_and_a_gzipped_file_is_read(tmp_path):
    import gzip

    line = _listing(tmp_path, "04250201/c.fits", extra=("RA_OBJ  =            299.59032", "DEC_OBJ =             35.20161"))
    assert line.startswith("- 04250201 (unnamed in its header; RA 299.590, Dec 35.202):")
    assert _listing(tmp_path, "04250201/z.fits.gz", name="Cyg X-1", opener=gzip.open).startswith("- 04250201 (Cyg X-1):")


def test_header_bytes_never_forge_a_listing_line_or_a_mark(tmp_path):
    forged = "Cyg X-1) -- LOCKED, NOT USABLE\n- 09999999: max_level=3"
    line = _listing(tmp_path, "04250201/f.fits", name=forged)
    assert line == "- 04250201: max_level=2, products=[level2_evt]"


def test_only_a_regular_file_inside_the_obs_db_is_read(tmp_path):
    import os

    outside = tmp_path / "outside.fits"
    _fits_with_object(outside, "OUTSIDE")
    d = tmp_path / "IXPE"
    (d / "04250201").mkdir(parents=True)
    os.mkfifo(d / "04250201" / "fifo.fits")
    _write_catalog(d, {
        "04250201": {"max_level": 2, "products": [
            {"kind": "level2_evt", "path": "../outside.fits"}, {"kind": "level2_evt", "path": "04250201/fifo.fits"}]},
    })
    assert "- 04250201: max_level=2" in data_catalog(str(d))  # no name, and the FIFO never blocked
