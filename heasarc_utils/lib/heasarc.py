"""HEASARC acquire engine — authoritative ObsID resolution + download.

The download mechanics are **mission-agnostic**: given a HEASARC master-table name,
``astroquery.heasarc`` resolves an ObsID to its archive location and downloads it, for any mission.
The archive's own link names the mission folder and tree (``.../FTP/ixpe/data/obs/01/01002901/``);
astroquery unpacks an observation as ``<location>/<obsid>/``. The only per-mission knowledge — the
master-table name and which file is the representative level-2 event — lives in the
``MissionProfile`` registry (``MISSIONS``, keyed by HEASARC's name for the mission), the single
extension point here.

The master-table row also states the observation's processing (``processing_date``,
``software_version``, ``caldb_version``) and when it becomes public (``public_date``), columns all
three registered master tables carry (HEASARC TAP schema, checked 2026-11-10): :func:`archive_state`
reads them, and ``heasarc_help.acquire`` makes them part of what identifies a download, so a
reprocessed or newly public observation is a new product.

Grounded in the current ``astroquery.heasarc`` API (v0.4.x): ``query_tap -> locate_data ->
download_data(location=...)``. See https://astroquery.readthedocs.io/en/stable/heasarc/heasarc.html

``astroquery`` is imported lazily (via :func:`_heasarc`) so this module imports cleanly without it and
tests can inject a fake.
"""

from __future__ import annotations

import fnmatch
import glob
import gzip
import math
import os
import posixpath
import re
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Callable

#: Encryption wrapper a file may arrive in: HEASARC serves an observation still inside its
#: proprietary period GPG-encrypted (``....fits.gz.gpg``).
LOCK_SUFFIX = ".gpg"

__all__ = ["HeasarcError", "MissionProfile", "MISSIONS", "LOCK_SUFFIX", "split_lock", "is_locked",
           "profile", "resolve", "archive_state", "file_state", "state_mismatch", "archive_path", "download", "download_obsid",
           "gunzip_tree", "discover_master", "mjd_now"]


def split_lock(name: str) -> tuple[str, bool]:
    """``(name without the encryption wrapper, whether it had one)``: ``"x.fits.gz.gpg"`` ->
    ``("x.fits.gz", True)``. Exactly :data:`LOCK_SUFFIX` is recognised."""
    if name.endswith(LOCK_SUFFIX):
        return name[: -len(LOCK_SUFFIX)], True
    return name, False


def is_locked(path: str) -> bool:
    """Whether ``path`` names an encrypted (unreadable as is) file."""
    return split_lock(os.path.basename(path))[1]


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
    """Resolve an ObsID via a HEASARC master table: ``(its row, as a table; its download links)``.

    Looks the ObsID up in ``catalog`` with a TAP/ADQL query and returns the row with the
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
    return tab, heasarc.locate_data(tab, catalog_name=catalog)


def mjd_now() -> float:
    """The current time as a Modified Julian Date."""
    return (datetime.now(timezone.utc) - datetime(1858, 11, 17, tzinfo=timezone.utc)).total_seconds() / 86400.0


def _cell(tab, column):
    """The first row's ``column`` as a plain Python value; ``None`` when the table lacks it or the
    cell is masked or blank."""
    try:
        value = tab[column][0]
    except (KeyError, IndexError, TypeError, ValueError):
        return None
    if value is None or getattr(value, "mask", False) is True or type(value).__name__ == "MaskedConstant":
        return None
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, bytes):
        value = value.decode("ascii", "replace")
    if isinstance(value, str):
        return value.strip() or None
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def archive_state(tab) -> dict:
    """What the master-table row says about the archive's copy of the observation:
    ``{"processing_date", "software_version", "caldb_version", "public_date"}`` (MJDs as floats, the
    versions as the archive spells them; ``None`` where the row holds nothing). The archive replaces
    files in place when it reprocesses an observation, and these columns are what changes then."""
    out = {name: _cell(tab, name)
           for name in ("processing_date", "software_version", "caldb_version", "public_date")}
    for name in ("processing_date", "public_date"):
        if out[name] is not None:
            try:
                out[name] = float(out[name])
            except (TypeError, ValueError):
                out[name] = None
    return out


def archive_path(links, obsid: str) -> tuple[str, list[str]]:
    """``(mission folder, tree)`` of an ObsID as HEASARC's archive holds it, read from the
    ``access_url`` of ``links`` (``.../FTP/<mission>/<tree...>/<obsid>/``): the folder names the
    local copy clones. Raises :class:`HeasarcError` with the archive's own message for a row that
    carries one (an ObsID with no products yet), and for links that disagree or do not end in the
    ObsID."""
    urls = [str(u) for u in links["access_url"]]
    errors = [str(e) for e in (links["error_message"] if "error_message" in getattr(links, "colnames", links) else [])]
    found = set()
    for i, url in enumerate(urls):
        if i < len(errors) and errors[i].strip():
            raise HeasarcError(f"HEASARC has no data for ObsID {obsid!r}: {errors[i].strip()}")
        parts = [p for p in url.split("FTP/", 1)[-1].split("/") if p]
        if "FTP/" not in url or len(parts) < 2 or parts[-1] != str(obsid):
            raise HeasarcError(f"archive link {url!r} is not .../FTP/<mission>/.../{obsid}/")
        found.add((parts[0], tuple(parts[1:-1])))
    if len(found) != 1:
        raise HeasarcError(f"archive links for ObsID {obsid!r} name {len(found)} places: {sorted(found)}")
    folder, tree = found.pop()
    return folder, list(tree)


def download(links, location: str, host: str | None = None,
             timeout_s: float = 1800.0) -> None:
    """Download resolved links into ``location`` via astroquery.

    Every astroquery host unpacks an observation as ``<location>/<obsid>/`` (the archive's tar
    is ObsID-rooted; the SciServer and S3 copies take the link's last folder), so the caller passes
    the folder the ObsID's own folder belongs in. The folder is created first.

    The astroquery ``Heasarc`` singleton carries a ``timeout`` attribute
    (default 30 s) that its internal ``_download_heasarc`` passes directly to
    the underlying HTTP transfer.  Large IXPE observations can exceed 1 GB and
    need several minutes — the factory default is far too short.

    This function temporarily raises ``h.timeout`` to ``timeout_s`` for the
    duration of the call and restores the previous value on exit, so concurrent
    or subsequent callers are unaffected.
    """
    os.makedirs(location, exist_ok=True)
    h = _heasarc()
    _prev_timeout = getattr(h, "timeout", None)
    h.timeout = timeout_s
    try:
        kwargs = {"location": location}
        if host is not None:
            kwargs["host"] = host
        h.download_data(links, **kwargs)
    finally:
        h.timeout = _prev_timeout


def gunzip_tree(obs_dir: str) -> int:
    """Decompress every ``*.gz`` under the ObsID directory in place; return the count.

    Fresh HEASARC downloads arrive gzipped (``*.evt.gz`` / ``*.fits.gz``); acquire
    decompresses them here so downstream reduction and the level-2 derivers operate
    on plain files. Each ``<name>.gz`` becomes ``<name>`` and the compressed copy is
    removed. Decompression writes a ``.part`` sidecar first and renames on success,
    so an interrupted gunzip leaves the original ``.gz`` intact rather than a
    half-written file masquerading as data.
    """
    root = str(obs_dir)
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
    level2: Callable[[list[str], str], str]       # (paths below the ObsID folder, obsid) -> one of them


def _pick(matches: list[str]) -> str | None:
    """Pick a single representative file: usable over encrypted, then uncompressed over its ``.gz``
    twin (:func:`gunzip_tree` normally leaves only the plain file; the rule stays gz-tolerant). A file
    still wrapped in ``.gpg`` ranks last, so a decrypted copy wins whenever one exists."""
    if not matches:
        return None
    usable = [m for m in matches if not is_locked(m)]
    for tier in (usable, matches):
        if tier:
            plain = [m for m in tier if not m.endswith(".gz")]
            return sorted(plain or tier)[0]
    return None


def _in(paths: list[str], folder: str, pattern: str) -> list[str]:
    """The ``paths`` lying directly in ``folder`` whose name matches the glob ``pattern``."""
    return [p for p in paths if posixpath.dirname(p) == folder and fnmatch.fnmatchcase(posixpath.basename(p), pattern)]


def _one(paths: list[str], obsid: str, folder: str, pattern: str, what: str) -> str:
    """The representative file matching ``folder/pattern`` among ``paths``; else raise. The pattern is
    wrapper-agnostic (it ends in ``*``), so an encrypted copy is found when it is all that was served."""
    found = _pick(_in(paths, folder, pattern))
    if found is None:
        raise HeasarcError(f"no {what} for ObsID {obsid!r} (looked for {folder}/{pattern})")
    return found


# base det1 level-2 event (compressed or not); excludes _rej/_src/_bkg/_pcube/_pha.
# Matched against the name with any encryption wrapper split off, so a proprietary
# observation served as ``....fits.gz.gpg`` is recognised as the level-2 event it is.
_IXPE_BASE_EVT = re.compile(r"_det1_evt2_(?P<ver>v\d+)\.fits(?:\.gz)?$")


def _ixpe_level2(paths: list[str], obsid: str) -> str:
    cands = [p for p in _in(paths, "event_l2", f"ixpe{obsid}_det1_evt2_*.fits*")
             if _IXPE_BASE_EVT.search(split_lock(posixpath.basename(p))[0])]
    base = _pick(cands)
    if base is None:
        raise HeasarcError(f"no IXPE level-2 event for ObsID {obsid!r}")
    return base


def _nicer_level2(paths: list[str], obsid: str) -> str:
    return _one(paths, obsid, "xti/event_cl", "*mpu7_cl.evt*", "NICER level-2 cleaned event")


def _nustar_level2(paths: list[str], obsid: str) -> str:
    return _one(paths, obsid, "event_cl", f"nu{obsid}A01_cl.evt*", "NuSTAR level-2 cleaned event (FPMA)")


#: Keyed by HEASARC's name for the mission, which is also its archive folder
#: (``https://heasarc.gsfc.nasa.gov/FTP/<name>/``) and the store's bucket; ``acquire`` refuses a link
#: naming another. A new mission: its master table (``discover_master``) and its level-2
#: representative here. Chandra (``chanmaster``) and MAXI (``maximaster``) are on HEASARC; MAXI's
#: master table returns no ``locate_data`` link for its rows (checked 2026-11-04), so it needs a locator
#: of its own before this engine can fetch it.
MISSIONS: dict[str, MissionProfile] = {
    "ixpe": MissionProfile(catalog="ixmaster", level2=_ixpe_level2),
    "nicer": MissionProfile(catalog="nicermastr", level2=_nicer_level2),
    "nustar": MissionProfile(catalog="numaster", level2=_nustar_level2),
}


def profile(mission: str) -> MissionProfile:
    """The registered profile of ``mission`` (HEASARC's name, any case); raises for another."""
    found = MISSIONS.get(mission.lower())
    if found is None:
        raise HeasarcError(
            f"no MissionProfile for {mission!r}; register one in MISSIONS "
            f"(see heasarc_utils/README.md; use discover_master() to find its "
            f"master table). Known missions: {sorted(MISSIONS)}"
        )
    return found


def file_state(path: str) -> dict:
    """The processing keywords a level-2 file's own header states: ``{"SOFTVER", "CALDBVER"}``, each
    from the primary header or else the first extension (``None`` where neither has it). HEASARC's
    master tables copy the processing script version into ``software_version`` and the CALDB version
    into ``caldb_version`` (checked 2026-11-10/11: IXPE 01250101, NuSTAR 30702017002 and NICER 4690010101, the
    row's values against the level-2 event headers' ``SOFTVER``/``CALDBVER``)."""
    from astropy.io import fits  # the archive environment holds it; imported where it is used

    out = {"SOFTVER": None, "CALDBVER": None}
    with fits.open(path, memmap=False) as hdus:
        for hdu in hdus[:2]:
            for key in out:
                value = hdu.header.get(key)
                if out[key] is None and isinstance(value, str) and value.strip():
                    out[key] = value.strip()
    return out


def state_mismatch(state: dict, path: str) -> str | None:
    """Why the level-2 file at ``path`` is not the processing the master-table ``state`` names, or
    ``None`` when its header agrees: its ``SOFTVER`` must equal the row's ``software_version`` (and be
    present when the row names one), its ``CALDBVER`` the row's ``caldb_version`` when both name one.
    The archive replaces files after its catalog row changes, so for a while the row can name a newer
    processing than the files it serves. The processing date has no header keyword to compare."""
    found = file_state(path)
    problems = []
    if state.get("software_version") and found["SOFTVER"] != state["software_version"]:
        problems.append(f"SOFTVER {found['SOFTVER']!r} where the row names software_version "
                        f"{state['software_version']!r}")
    if state.get("caldb_version") and found["CALDBVER"] and found["CALDBVER"] != state["caldb_version"]:
        problems.append(f"CALDBVER {found['CALDBVER']!r} where the row names caldb_version "
                        f"{state['caldb_version']!r}")
    return "; ".join(problems) or None


def download_obsid(links, obsid: str, folder: str, *, host: str | None = None,
                   timeout_s: float = 1800.0) -> list[str]:
    """Download the ObsID ``links`` name into ``folder`` and return the paths of its files, relative to
    ``folder`` and in posix form. astroquery unpacks to ``<location>/<obsid>/``; the ObsID folder's
    contents are moved up into ``folder`` (everything else the download left is removed), and
    gunzipped (:func:`gunzip_tree`). Raises :class:`HeasarcError` when the download yields no ObsID
    folder."""
    location = os.path.join(folder, ".download")
    download(links, location, host=host, timeout_s=timeout_s)
    obs_dir = os.path.join(location, str(obsid))
    if not os.path.isdir(obs_dir):
        raise HeasarcError(f"the download produced no folder {obsid}/ for ObsID {obsid!r}")
    for name in os.listdir(obs_dir):
        os.replace(os.path.join(obs_dir, name), os.path.join(folder, name))
    shutil.rmtree(location)
    gunzip_tree(folder)
    out = []
    for here, _dirs, files in os.walk(folder):
        out += [os.path.relpath(os.path.join(here, f), folder).replace(os.sep, "/") for f in files
                if not os.path.islink(os.path.join(here, f))]
    return sorted(out)
