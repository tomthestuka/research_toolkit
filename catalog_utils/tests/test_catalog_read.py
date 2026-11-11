"""data_catalog against a throwaway store (root conftest ``store``), its products seeded through the
host's seed channel (``seed``): the listing of what the catalog holds, each line with its record, the
source an event names, and the lock mark."""

import base64
import json

import pytest

from catalog_utils.lib.catalog_read import data_catalog


def _header(cards):
    return "".join(c.ljust(80) for c in [*cards, "END"]).ljust(2880).encode("ascii", "replace")


def _event(name=None, *, extension=True, extra=()):
    """A small IXPE-shaped event file: an empty primary, then an EVENTS header carrying OBJECT."""
    obj = [f"OBJECT  = '{name}'"] if name is not None else []
    primary = ["SIMPLE  =                    T", "BITPIX  =                    8", "NAXIS   =                    0"]
    events = ["XTENSION= 'BINTABLE'", "BITPIX  =                    8", "NAXIS   =                    2",
              "NAXIS1  =                    0", "NAXIS2  =                    0", "EXTNAME = 'EVENTS'", *obj, *extra]
    return _header(primary) + _header(events) if extension else _header([*primary, *obj, *extra])


def _product(seed, kind, obsid, files: dict[str, bytes], mission="ixpe"):
    """A synthetic product of ``kind`` for ``obsid`` holding ``files`` (relative path: bytes), seeded; its
    record."""
    encoded = {rel: base64.b64encode(data).decode("ascii") for rel, data in files.items()}
    return seed(kind=kind, subject=obsid, bucket=mission, files=encoded)["record"]


@pytest.fixture
def obsid(store):
    """An ObsID of this test alone (a store may outlive a test)."""
    return "9" + str(int(store["tag"], 16))[:7]


def _line(obsid):
    return next(line for line in data_catalog(subject=obsid).splitlines() if line.startswith(f"- {obsid}"))


def test_each_product_is_listed_with_its_record_and_the_source_its_event_names(seed, obsid):
    record = _product(seed, "raw", obsid, {f"event_l2/ixpe{obsid}_det1_evt2_v01.fits": _event("IGR J17091-3624"),
                                 "auxil/x.txt": b"x"})
    listing = data_catalog(mission="ixpe", subject=obsid)
    assert listing.splitlines()[:2] == ["Stored products (1):", "ixpe:"]
    line = _line(obsid)
    assert line.startswith(f"- {obsid} (IGR J17091-3624) raw [verified, 2 files]: record ")
    assert json.loads(line.split("record ", 1)[1]) == record


def test_a_header_without_a_name_shows_its_coordinates(seed, obsid):
    _product(seed, "raw", obsid, {f"event_l2/ixpe{obsid}_det1_evt2_v01.fits": _event(
        extra=("RA_OBJ  =            299.59032", "DEC_OBJ =             35.20161"))})
    assert _line(obsid).startswith(f"- {obsid} (unnamed in its header; RA 299.590, Dec 35.202) raw")


def test_a_name_in_the_primary_header_is_read(seed, obsid):
    _product(seed, "raw", obsid, {f"event_l2/ixpe{obsid}_det1_evt2_v01.fits": _event("GX 339-4", extension=False)})
    assert _line(obsid).startswith(f"- {obsid} (GX 339-4) raw")


def test_header_bytes_never_forge_a_listing_line_or_a_mark(seed, obsid):
    forged = "Cyg X-1) -- LOCKED, NOT USABLE\n- 09999999: max_level=3"
    _product(seed, "raw", obsid, {f"event_l2/ixpe{obsid}_det1_evt2_v01.fits": _event(forged)})
    listing = data_catalog(subject=obsid)
    assert "LOCKED" not in listing and "09999999" not in listing
    assert _line(obsid).startswith(f"- {obsid} raw [verified, 1 files]")


def test_an_encrypted_event_marks_its_line_locked(seed, obsid):
    _product(seed, "raw", obsid, {f"event_l2/ixpe{obsid}_det1_evt2_v03.fits.gz.gpg": b"encrypted"})
    listing = data_catalog(subject=obsid)
    assert "LOCKED, NOT USABLE" in listing and "proprietary period" in listing


def test_a_filter_that_matches_nothing_lists_nothing(store):
    assert data_catalog(subject="no-such-obsid") == "Stored products (0):\n(no product stored)"


def test_a_product_of_another_kind_is_listed_beside_the_download(seed, obsid):
    _product(seed, "raw", obsid, {f"event_l2/ixpe{obsid}_det1_evt2_v01.fits": _event("Cyg X-1")})
    _product(seed, "pha", obsid, {"spec.pha": b"spectrum"})
    lines = [line for line in data_catalog(subject=obsid).splitlines() if line.startswith("- ")]
    assert [line.split(" [")[0] for line in lines] == [f"- {obsid} pha", f"- {obsid} (Cyg X-1) raw"]
