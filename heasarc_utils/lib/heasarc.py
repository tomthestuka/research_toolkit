"""HEASARC acquire engine — authoritative ObsID resolution + download.

The download mechanics are **mission-agnostic**: given a HEASARC master-table
name, ``astroquery.heasarc`` resolves an ObsID to its archive location and
downloads it, for any mission. The only per-mission knowledge — the master-table
name and which file is the representative level-2 event — lives in the
``MissionProfile`` registry (``MISSIONS``), which is the single extension point.
Adding a mission is one ``MISSIONS`` entry against an unchanged engine.

Grounded in the current ``astroquery.heasarc`` API (v0.4.x):
``query_region(catalog=..., column_filters={"obsid": ...}) -> locate_data ->
download_data(location=...)``. See
https://astroquery.readthedocs.io/en/stable/heasarc/heasarc.html

``astroquery`` is imported lazily (via :func:`_heasarc`) so this module imports
cleanly without it and tests can inject a fake.

The level-2 derivers here mirror ``obsdb_catalog._domain`` exactly — the two must
stay in sync so a freshly recorded download and a ``repair`` re-derivation agree.
"""

from __future__ import annotations

import glob
import gzip
import os
import re
import shutil
from dataclasses import dataclass
from typing import Callable

__all__ = ["HeasarcError", "MissionProfile", "MISSIONS", "acquire_obsid",
           "resolve", "download", "normalize_layout", "gunzip_tree",
           "discover_master"]

# Notes mirrored from obsdb_catalog._domain (representative-product convention).
_IXPE_NOTE = "det1 listed; det2/det3 follow same naming pattern"
_NUSTAR_NOTE = "fpmA listed; fpmB follows same naming pattern"


class HeasarcError(Exception):
    """Raised when resolution or download cannot produce the expected data."""


def _heasarc():
    """Lazy accessor for the astroquery Heasarc class (the test seam)."""
    from astroquery.heasarc import Heasarc  # imported only when actually acquiring
    return Heasarc


# --------------------------------------------------------------------------- #
# Mission-agnostic engine
# --------------------------------------------------------------------------- #
def resolve(catalog: str, obsid: str):
    """Resolve an ObsID to downloadable links via a HEASARC master table.

    Looks the ObsID up in ``catalog`` with a TAP/ADQL query and returns the
    ``locate_data`` links. ``query_tap`` is used (not ``query_region`` with
    ``column_filters``, which drops the filter for a position-less query — verified
    against the installed astroquery v0.4.x); ``SELECT *`` carries the ``__row``
    column that ``locate_data`` requires. Raises :class:`HeasarcError` if the
    ObsID is not found (no fabrication — a missing ObsID fails loudly).
    """
    obsid = str(obsid)
    if not obsid.isalnum():  # ADQL is string-formatted; keep the literal safe
        raise HeasarcError(f"invalid obsid {obsid!r} (expected alphanumeric)")
    heasarc = _heasarc()
    tab = heasarc.query_tap(
        f"SELECT * FROM {catalog} WHERE obsid='{obsid}'"
    ).to_table()
    if tab is None or len(tab) == 0:
        raise HeasarcError(
            f"ObsID {obsid!r} not found in HEASARC catalog {catalog!r}"
        )
    return heasarc.locate_data(tab, catalog_name=catalog)


def download(links, obs_db_path: str, host: str | None = None,
             timeout_s: float = 1800.0) -> None:
    """Download resolved links into ``obs_db_path`` via astroquery.

    The astroquery ``Heasarc`` singleton carries a ``timeout`` attribute
    (default 30 s) that its internal ``_download_heasarc`` passes directly to
    the underlying HTTP transfer.  Large IXPE observations can exceed 1 GB and
    need several minutes — the factory default is far too short.

    This function temporarily raises ``h.timeout`` to ``timeout_s`` for the
    duration of the call and restores the previous value on exit, so concurrent
    or subsequent callers are unaffected.
    """
    h = _heasarc()
    _prev_timeout = getattr(h, "timeout", None)
    h.timeout = timeout_s
    try:
        kwargs = {"location": obs_db_path}
        if host is not None:
            kwargs["host"] = host
        h.download_data(links, **kwargs)
    finally:
        h.timeout = _prev_timeout


def normalize_layout(obs_db_path: str, obsid: str) -> str:
    """Ensure the ObsID dir sits at ``<obs_db_path>/<obsid>``; return its path.

    ``download_data`` untars ``heasarc-data.tar``; if it nests the ObsID under an
    extra prefix, relocate it to the expected top level. Raises
    :class:`HeasarcError` if the ObsID directory cannot be located.
    """
    target = os.path.join(obs_db_path, str(obsid))
    if os.path.isdir(target):
        return target
    candidates = [
        d for d in glob.glob(os.path.join(obs_db_path, "**", str(obsid)), recursive=True)
        if os.path.isdir(d)
    ]
    if not candidates:
        raise HeasarcError(
            f"download produced no directory for ObsID {obsid!r} under {obs_db_path}"
        )
    shutil.move(candidates[0], target)
    return target


def gunzip_tree(obs_db_path: str, obsid: str) -> int:
    """Decompress every ``*.gz`` under the ObsID directory in place; return the count.

    Fresh HEASARC downloads arrive gzipped (``*.evt.gz`` / ``*.fits.gz``); acquire
    decompresses them here so downstream reduction and the level-2 derivers operate
    on plain files. Each ``<name>.gz`` becomes ``<name>`` and the compressed copy is
    removed. Decompression writes a ``.part`` sidecar first and renames on success,
    so an interrupted gunzip leaves the original ``.gz`` intact rather than a
    half-written file masquerading as data.
    """
    root = os.path.join(obs_db_path, str(obsid))
    count = 0
    for gz in sorted(glob.glob(os.path.join(root, "**", "*.gz"), recursive=True)):
        dest = gz[:-3]                      # strip the trailing ".gz"
        part = dest + ".part"
        with gzip.open(gz, "rb") as fin, open(part, "wb") as fout:
            shutil.copyfileobj(fin, fout, length=16 * 1024 * 1024)
        os.replace(part, dest)
        os.remove(gz)
        count += 1
    return count


def discover_master(mission: str):
    """List candidate HEASARC master tables for a mission name (authoring aid).

    Convenience for registering a new :class:`MissionProfile`; not used on the
    download path. Returns whatever ``Heasarc.list_catalogs`` yields.
    """
    return _heasarc().list_catalogs(master=True, keywords=mission)


# --------------------------------------------------------------------------- #
# The extension point: per-mission profiles
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class MissionProfile:
    """A mission's acquire knowledge: master table + level-2 representative."""
    catalog: str                                  # HEASARC master-table name
    level2: Callable[[str, str], dict]            # (obs_db_path, obsid) -> product dict


def _pick(matches: list[str]) -> str | None:
    """Pick a single representative file, preferring an uncompressed file over its
    ``.gz`` twin. Fresh archive downloads are gzipped (``*.evt.gz`` / ``*.fits.gz``);
    :func:`gunzip_tree` decompresses them during acquire, so the plain file is normally
    present — glob both and prefer the plain one (gz-tolerant for robustness).
    """
    if not matches:
        return None
    plain = [m for m in matches if not m.endswith(".gz")]
    return sorted(plain or matches)[0]


def _one(obs_db_path: str, obsid: str, pattern: str, what: str) -> str:
    """Return the representative file matching ``pattern`` under the ObsID; else raise."""
    full = _pick(glob.glob(os.path.join(obs_db_path, str(obsid), pattern)))
    if full is None:
        raise HeasarcError(f"no {what} for ObsID {obsid!r} (looked for {pattern!r})")
    return full


# base det1 level-2 event (compressed or not); excludes _rej/_src/_bkg/_pcube/_pha.
_IXPE_BASE_EVT = re.compile(r"_det1_evt2_(?P<ver>v\d+)\.fits(?:\.gz)?$")


def _ixpe_level2(obs_db_path: str, obsid: str) -> dict:
    cands = [
        f for f in glob.glob(os.path.join(
            obs_db_path, str(obsid), "event_l2", f"ixpe{obsid}_det1_evt2_*.fits*"))
        if _IXPE_BASE_EVT.search(os.path.basename(f))
    ]
    base = _pick(cands)
    if base is None:
        raise HeasarcError(f"no IXPE level-2 event for ObsID {obsid!r}")
    ver = _IXPE_BASE_EVT.search(os.path.basename(base)).group("ver")
    return {
        "kind": "level2_evt",
        "path": os.path.relpath(base, obs_db_path),
        "produced_by": "heasarc",
        "params": {"detectors": [1, 2, 3], "calibration_version": ver, "note": _IXPE_NOTE},
    }


def _nicer_level2(obs_db_path: str, obsid: str) -> dict:
    full = _one(obs_db_path, obsid, os.path.join("xti", "event_cl", "*mpu7_cl.evt*"),
                "NICER level-2 cleaned event")
    return {
        "kind": "level2_evt",
        "path": os.path.relpath(full, obs_db_path),
        "produced_by": "heasarc",
        "params": {"mpu": 7},
    }


def _nustar_level2(obs_db_path: str, obsid: str) -> dict:
    full = _one(obs_db_path, obsid, os.path.join("event_cl", f"nu{obsid}A01_cl.evt*"),
                "NuSTAR level-2 cleaned event (FPMA)")
    return {
        "kind": "level2_evt",
        "path": os.path.relpath(full, obs_db_path),
        "produced_by": "heasarc",
        "params": {"fpm": "A", "note": _NUSTAR_NOTE},
    }


MISSIONS: dict[str, MissionProfile] = {
    "ixpe": MissionProfile(catalog="ixmaster", level2=_ixpe_level2),
    "nicer": MissionProfile(catalog="nicermastr", level2=_nicer_level2),
    "nustar": MissionProfile(catalog="numaster", level2=_nustar_level2),
}


def acquire_obsid(mission: str, obsid: str, obs_db_path: str,
                  *, host: str | None = None,
                  timeout_s: float = 1800.0) -> dict:
    """Resolve → download → normalize → gunzip → derive the level-2 product. No recording.

    Returns the level-2 product dict (relative ``path`` confirmed on disk). The
    high-level :func:`heasarc_utils.lib.heasarc_help.acquire` adds catalog
    recording. Raises :class:`HeasarcError` for an unknown mission, a missing
    ObsID, or a download that does not yield the expected level-2 file.

    Parameters
    ----------
    timeout_s : float
        HTTP transfer timeout in seconds, passed through to :func:`download`.
        Default 1800 s (30 min). The astroquery factory default (30 s) is too
        short for large IXPE transfers.
    """
    profile = MISSIONS.get(mission.lower())
    if profile is None:
        raise HeasarcError(
            f"no MissionProfile for {mission!r}; register one in MISSIONS "
            f"(see heasarc_utils/README.md; use discover_master() to find its "
            f"master table). Known missions: {sorted(MISSIONS)}"
        )
    links = resolve(profile.catalog, obsid)
    download(links, obs_db_path, host=host, timeout_s=timeout_s)
    normalize_layout(obs_db_path, obsid)
    gunzip_tree(obs_db_path, obsid)
    return profile.level2(obs_db_path, obsid)
