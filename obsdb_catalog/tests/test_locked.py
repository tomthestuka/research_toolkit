"""The encryption-wrapper rule: a .gpg product is seen, not read as empty ground.

HEASARC serves an observation still inside its proprietary period GPG-encrypted.
Before this rule existed, every product pattern anchored on ``.fits``/``.evt`` and
so rejected the encrypted name, which made a complete download read as an empty
directory at every layer: the scan dropped the ObsID, healthcheck compared against
that same scan and had nothing to flag, and the listing never mentioned it — so
re-acquiring was the only inference available and it could never succeed.
"""

import json
import os

from obsdb_catalog import _domain, healthcheck, neutral, repair

IXPE_EVT = "ixpe{obsid}_det1_evt2_{ver}.fits.gz"


def _write(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(b"x")


def _ixpe_obsid(root, obsid, *, ver="v03", locked=False, derived=None):
    """Lay down an IXPE ObsID directory; ``locked`` wraps its event in .gpg."""
    name = IXPE_EVT.format(obsid=obsid, ver=ver) + (neutral.LOCK_SUFFIX if locked else "")
    _write(os.path.join(root, obsid, "event_l2", name))
    for extra in derived or []:
        _write(os.path.join(root, obsid, "event_l2", extra))
    return name


# ── the rule itself ────────────────────────────────────────────────────────


def test_split_lock_strips_only_the_known_wrapper():
    assert neutral.split_lock("evt2_v03.fits.gz.gpg") == ("evt2_v03.fits.gz", True)
    assert neutral.split_lock("evt2_v03.fits.gz") == ("evt2_v03.fits.gz", False)


def test_split_lock_leaves_an_unknown_suffix_unknown():
    """An unfamiliar suffix is not guessed at.

    Deciding today that some other wrapper means "encrypted" would mislabel a file's
    usability in a way that is very hard to debug later, so an unknown suffix stays
    unrecognised and simply fails to match a product pattern, exactly as before.
    """
    name = "evt2_v03.fits.gz.xyz"
    assert neutral.split_lock(name) == (name, False)
    assert neutral.is_locked(name) is False


def test_is_locked_reads_the_path():
    assert neutral.is_locked("a/b/evt2.fits.gz.gpg") is True
    assert neutral.is_locked("a/b/evt2.fits.gz") is False


# ── the scanner ────────────────────────────────────────────────────────────


def test_scan_sees_an_all_encrypted_obsid(tmp_path):
    """The regression: a complete but encrypted download is scanned, not dropped."""
    root = tmp_path / "IXPE"
    name = _ixpe_obsid(root, "04004001", locked=True)

    scanned = _domain.scan("ixpe", str(root))

    assert "04004001" in scanned, "an encrypted ObsID must not vanish from the scan"
    (product,) = scanned["04004001"]
    assert product["kind"] == "level2_evt"
    assert product["path"].endswith(name)
    assert neutral.is_locked(product["path"])
    assert product["params"]["calibration_version"] == "v03", (
        "the version is read through the wrapper, not lost with it"
    )


def test_scan_prefers_a_decrypted_copy_over_its_encrypted_twin(tmp_path):
    root = tmp_path / "IXPE"
    plain = _ixpe_obsid(root, "01002901", ver="v01")
    _write(os.path.join(root, "01002901", "event_l2", plain + neutral.LOCK_SUFFIX))

    (product,) = _domain.scan("ixpe", str(root))["01002901"]

    assert not neutral.is_locked(product["path"]), (
        "with both copies present the readable one represents the ObsID"
    )


def test_scan_ignores_an_unknown_suffix(tmp_path):
    """An unfamiliar wrapper is not silently promoted to a usable product."""
    root = tmp_path / "IXPE"
    _write(os.path.join(root, "05250399", "event_l2",
                        IXPE_EVT.format(obsid="05250399", ver="v01") + ".xyz"))

    assert _domain.scan("ixpe", str(root)) == {}


def test_scan_skips_a_human_working_copy_of_an_obsid(tmp_path):
    """A directory holding another ObsID's files is not that ObsID.

    Quarantine and working copies sit beside the real data under a suffixed name
    (``<obsid>_bad_gpg``, ``<obsid>_picorr``) while still holding files stamped with
    the original ObsID. The scan is anchored on the directory's own ObsID — as the
    acquire-side deriver already was — so such a copy cannot enter the catalog as an
    invented observation. Live disk holds two of these for the encrypted ObsIDs, and
    they became visible for the first time when the scan learned to read `.gpg`.
    """
    root = tmp_path / "IXPE"
    _ixpe_obsid(root, "04004001", locked=True)
    # a quarantined copy: directory renamed, files still stamped 04004001
    _write(os.path.join(root, "04004001_bad_gpg", "event_l2",
                        IXPE_EVT.format(obsid="04004001", ver="v03") + neutral.LOCK_SUFFIX))

    scanned = _domain.scan("ixpe", str(root))

    assert set(scanned) == {"04004001"}


def test_scan_sees_an_encrypted_nicer_event(tmp_path):
    """The rule is mission-wide, not an IXPE special case."""
    root = tmp_path / "NICER"
    _write(os.path.join(root, "4690010101", "xti", "event_cl",
                        "ni4690010101_0mpu7_cl.evt.gz" + neutral.LOCK_SUFFIX))

    (product,) = _domain.scan("nicer", str(root))["4690010101"]

    assert product["kind"] == "level2_evt" and neutral.is_locked(product["path"])


# ── the layers that read the scan ──────────────────────────────────────────


def test_healthcheck_flags_an_uncatalogued_encrypted_obsid(tmp_path):
    """Sight restored: the check that should have caught this now can.

    ``uncatalogued`` is computed from the scan, so while the scan was blind there
    was nothing to compare and the ObsID was invisible rather than reported.
    """
    root = tmp_path / "IXPE"
    _ixpe_obsid(root, "04004001", locked=True)
    (root / "catalog.json").write_text(
        json.dumps({"mission": "ixpe", "entries": {}}), encoding="utf-8")

    report = healthcheck.healthcheck(str(root))

    assert report is not True
    assert any("04004001" in path for path in report["uncatalogued"])


def test_repair_catalogues_an_encrypted_obsid_at_its_honest_level(tmp_path):
    """Repair records it, levelled for what it is: a level-2 event that is present.

    The level stays honest — the file genuinely is the level-2 event — and the
    ``.gpg`` path is what tells every reader it cannot be opened. Nothing has to be
    demoted or hidden for that to be unambiguous.
    """
    root = tmp_path / "IXPE"
    _ixpe_obsid(root, "04004001", locked=True)

    repair.repair(str(root), log=None)

    entry = neutral.read(str(root))["entries"]["04004001"]
    assert entry["max_level"] == 2 and entry["levels_present"] == [1, 2]
    (product,) = entry["products"]
    assert neutral.is_locked(product["path"])
    assert healthcheck.healthcheck(str(root)) is True, (
        "once catalogued, the encrypted ObsID is no longer drift"
    )
