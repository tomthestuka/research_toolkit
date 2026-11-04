"""Mission layout knowledge — used only by :mod:`repair` and :mod:`healthcheck`.

:mod:`obsdb_catalog.neutral` is kept domain-free; all per-mission knowledge
(directory layout, filename conventions, which ``kind`` is which level, the
constant reduction params) lives here. This module is plain standard library
(glob/re/os); it does **not** import the science stack or ``auto_research``.

The conventions mirror the live bootstrapped catalogs:

Each ObsID's folder is found by :data:`TREES` (HEASARC's archive tree, or flat ``<obsid>/``):

* **IXPE** ``<obsid>/event_l2/`` — one representative ``det1`` entry per product
  kind. ``level2_evt`` (``produced_by:"heasarc"``), ``pha``
  (``produced_by:"filter_inst_bkg_IXPE.sh"``). A PCUBE depends on a chosen
  energy binning and regions, so it is a science artifact, not a mission-team
  product: a ``*_pcube.fits`` on disk is never derived into the catalog. ``filter`` is read from the
  filename infix (``_rej_src_`` -> ``"rej_src"``; ``_src_`` -> ``"src"``, the
  ``03002599`` special case); ``calibration_version`` from the ``v0N`` token.
* **NICER** ``<obsid>/xti/event_cl/`` — ``pha`` / ``pha_bkg`` / ``arf`` / ``rmf``
  (``produced_by:"nicer_red.sh"``), constant params ``background_model:"3C50"``,
  ``mpu:7``.
* **NuSTAR** ``<obsid>/products/`` — per FPM A/B: ``pha_grp`` (``nustar_grp.sh``,
  ``grouping:"optimal"``), ``pha_src`` / ``arf`` / ``rmf`` (``nustar_red.sh``).

Every product pattern here matches a filename with its encryption wrapper split off
(:func:`neutral.split_lock`), so an observation the archive served GPG-encrypted —
still inside its proprietary period — is scanned as the products it holds rather than
read as empty ground. The scan records the ``.gpg`` path as-is; that path is what tells
every reader the product is present but unusable (:func:`neutral.is_locked`).
"""

from __future__ import annotations

import glob
import os
import re

from .neutral import is_locked, split_lock

#: Where an ObsID's folder sits below a mission folder (``<obs_db_root>/<mission>/``), as glob
#: patterns on ``{obsid}``: HEASARC's archive tree as served — the path after ``/FTP/<mission>/`` in
#: the ``access_url`` its master table's ``locate_data`` returns (checked 2026-11-04:
#: ``ixpe/data/obs/01/01002901/``, ``nicer/data/obs/2022_02/4690010101/``,
#: ``nustar/data/obs/07/3/30702017002/``, ``chandra/data/byobsid/4/3814/``) — then the flat
#: ``<obsid>/`` of obs-dbs built before. A mission absent here is read flat. The mission folder
#: is HEASARC's own name for the mission, matched case-insensitively (``mission_from_path``).
TREES: dict[str, tuple[str, ...]] = {
    "ixpe": ("data/obs/*/{obsid}", "{obsid}"),
    "nicer": ("data/obs/*/{obsid}", "{obsid}"),
    "nustar": ("data/obs/*/*/{obsid}", "{obsid}"),
    "chandra": ("data/byobsid/*/{obsid}",),
}


def obsid_dirs(mission: str, obs_db_path: str) -> list[str]:
    """Every folder under a mission folder that sits where an ObsID's folder would (:data:`TREES`);
    a scanner keeps those holding its mission's subfolders. The folder name is the ObsID."""
    found: set[str] = set()
    for pattern in TREES.get(mission.lower(), ("{obsid}",)):
        found.update(d for d in glob.glob(os.path.join(obs_db_path, pattern.format(obsid="*")))
                     if os.path.isdir(d))
    return sorted(found)


#: Marker stamped into a param whose value cannot be re-derived during repair.
UNRECOVERABLE = "repaired session, unrecoverable"

_LEVEL2_KINDS = {"level2_evt", "raw_evt"}
_LEVEL1_KINDS = {"raw", "level1"}

#: Representative-product note, mirrored by heasarc_utils when it records level-2 events.
_NUSTAR_NOTE = "fpmA listed; fpmB follows same naming pattern"


def mission_from_path(obs_db_path: str) -> str:
    return os.path.basename(os.path.normpath(obs_db_path)).lower()


def kind_level(kind: str) -> int:
    """Map a product ``kind`` to its processing level.

    Archive download supplies up to level 2 (cleaned events); everything else
    (pcube, pha, arf, rmf, grouped spectra, ...) is a level-3 reduction product.
    """
    if kind in _LEVEL1_KINDS:
        return 1
    if kind in _LEVEL2_KINDS:
        return 2
    return 3


def _rel(obs_db_path: str, full: str) -> str:
    return os.path.relpath(full, obs_db_path)


def _prefer(matches: list[str]) -> str | None:
    """Pick one representative: usable over encrypted, then uncompressed over ``.gz``.

    Fresh archive downloads are gzipped (``*.evt.gz`` / ``*.fits.gz``) and reduction
    later gunzips them; glob both and prefer the plain file when present. A file
    still wrapped in :data:`neutral.LOCK_SUFFIX` is ranked last, so an ObsID that
    holds both an encrypted and a decrypted copy is represented by the copy that can
    actually be read. Mirrors ``heasarc_utils.lib.heasarc._pick`` — keep the two in sync.
    """
    if not matches:
        return None
    usable = [m for m in matches if not is_locked(m)]
    for tier in (usable, matches):
        if tier:
            plain = [m for m in tier if not m.endswith(".gz")]
            return sorted(plain or tier)[0]
    return None


# --------------------------------------------------------------------------- #
# Per-mission disk scanners: obs_db_path -> {obsid: [product, ...]}
# Products are fully formed (kind, path, produced_by, params); files exist
# (they came from glob). created_at is left to the writer.
# --------------------------------------------------------------------------- #
def scan(mission: str, obs_db_path: str) -> dict[str, list[dict]]:
    fn = _SCANNERS.get(mission.lower())
    return fn(obs_db_path) if fn else {}


def _scan_nicer(obs_db_path: str) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for obs_dir in obsid_dirs("nicer", obs_db_path):
        evt_cl = os.path.join(obs_dir, "xti", "event_cl")
        if not os.path.isdir(evt_cl):
            continue
        obsid = os.path.basename(obs_dir)
        products = []
        # level-2 cleaned event (archive; gzipped on download, gunzipped by reduction).
        l2 = _prefer(glob.glob(os.path.join(evt_cl, "*mpu7_cl.evt*")))
        if l2:
            products.append({
                "kind": "level2_evt",
                "path": _rel(obs_db_path, l2),
                "produced_by": "heasarc",
                "params": {"mpu": 7},
            })
        # level-3 nicerl3-spect reduction products.
        specs = [
            ("pha", "*sr.pha", {"background_model": "3C50", "mpu": 7}),
            ("pha_bkg", "*bg.pha", {"background_model": "3C50", "mpu": 7}),
            ("arf", "*mpu7.arf", {"mpu": 7}),
            ("rmf", "*mpu7.rmf", {"mpu": 7}),
        ]
        for kind, pattern, params in specs:
            for full in sorted(glob.glob(os.path.join(evt_cl, pattern))):
                products.append({
                    "kind": kind,
                    "path": _rel(obs_db_path, full),
                    "produced_by": "nicer_red.sh",
                    "params": dict(params),
                })
        if products:
            out[obsid] = products
    return out


def _scan_nustar(obs_db_path: str) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    # A downloaded ObsID has event_cl (level-2); a reduced one also has products/ (level-3).
    folders = {os.path.basename(d): d for d in obsid_dirs("nustar", obs_db_path)
               if os.path.isdir(os.path.join(d, "event_cl")) or os.path.isdir(os.path.join(d, "products"))}

    # kind suffix, glob, produced_by, extra params (besides fpm) for level-3 products
    specs = [
        ("pha_grp", "*{fpm}01_sr_grp.pha", "nustar_grp.sh", {"grouping": "optimal"}),
        ("pha_src", "*{fpm}01_sr.pha", "nustar_red.sh", {}),
        ("arf", "*{fpm}01_sr.arf", "nustar_red.sh", {}),
        ("rmf", "*{fpm}01_sr.rmf", "nustar_red.sh", {}),
    ]
    for obsid in sorted(folders):
        products = []
        # level-2 cleaned events: representative FPMA, from the archive download.
        l2 = _prefer(glob.glob(
            os.path.join(folders[obsid], "event_cl", f"nu{obsid}A01_cl.evt*")
        ))
        if l2:
            products.append({
                "kind": "level2_evt",
                "path": _rel(obs_db_path, l2),
                "produced_by": "heasarc",
                "params": {"fpm": "A", "note": _NUSTAR_NOTE},
            })
        # level-3 reduction products, per FPM.
        prod_dir = os.path.join(folders[obsid], "products")
        for fpm in ("A", "B"):
            for kind, pat, produced_by, extra in specs:
                for full in sorted(glob.glob(os.path.join(prod_dir, pat.format(fpm=fpm)))):
                    params = {"fpm": fpm}
                    params.update(extra)
                    products.append({
                        "kind": f"{kind}_fpm{fpm}",
                        "path": _rel(obs_db_path, full),
                        "produced_by": produced_by,
                        "params": params,
                    })
        if products:
            out[obsid] = products
    return out


# IXPE level-2 event file for the det1 representative: ixpe<obsid>_det1_evt2_v<NN>.fits
_IXPE_EVT = re.compile(r"_det1_evt2_(?P<ver>v\d+)\.fits(?:\.gz)?$")
# pha derived names carry the filter infix and calibration version. PCUBEs are artifacts (user
# ruling 2026-11-04) and stay out of the catalog.
_IXPE_DERIVED = re.compile(
    r"_det1_evt2_(?P<ver>v\d+)_(?P<filter>rej_src|src)_pha1\.fits$"
)
_IXPE_NOTE = "det1 listed; det2/det3 follow same naming pattern"


def _scan_ixpe(obs_db_path: str) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for obs_dir in obsid_dirs("ixpe", obs_db_path):
        evt_l2 = os.path.join(obs_dir, "event_l2")
        if not os.path.isdir(evt_l2):
            continue
        obsid = os.path.basename(obs_dir)
        # Anchored on this directory's ObsID, mirroring heasarc._ixpe_level2: the
        # directory name *is* the ObsID, so a file carrying a different one is a
        # human working copy (``<obsid>_bad_gpg``, ``<obsid>_picorr``) rather than
        # this observation's product, and cataloguing it would invent an ObsID.
        names = sorted(
            os.path.basename(f)
            for f in glob.glob(os.path.join(evt_l2, f"ixpe{obsid}_det1_evt2_*.fits*"))
        )
        # Patterns match the name with any encryption wrapper split off, so an
        # encrypted product is recognised as the product it is. Recording it keeps
        # the ObsID visible; readers see the ``.gpg`` path and know it is unusable.
        derived = [(n, m) for n in names if (m := _IXPE_DERIVED.search(split_lock(n)[0]))]
        # Per the bootstrap convention, the catalogued products use the
        # instrumental-background-rejected filter when it exists, falling back to
        # the non-rejected "src" filter (the 03002599 special case). The other
        # variant's files are intentionally left out of the catalog.
        chosen = "rej_src" if any(m.group("filter") == "rej_src" for _, m in derived) \
            else "src"

        products = []
        # base det1 level-2 event (single representative; prefer uncompressed over .gz).
        base = _prefer([os.path.join(evt_l2, n) for n in names
                        if _IXPE_EVT.search(split_lock(n)[0])])
        if base:
            ver = _IXPE_EVT.search(split_lock(os.path.basename(base))[0]).group("ver")
            products.append({
                "kind": "level2_evt",
                "path": _rel(obs_db_path, base),
                "produced_by": "heasarc",
                "params": {
                    "detectors": [1, 2, 3],
                    "calibration_version": ver,
                    "note": _IXPE_NOTE,
                },
            })
        for name, m in derived:
            if m.group("filter") != chosen:
                continue
            ver = m.group("ver")
            products.append({
                "kind": "pha",
                "path": _rel(obs_db_path, os.path.join(evt_l2, name)),
                "produced_by": "filter_inst_bkg_IXPE.sh",
                "params": {
                    "filter": chosen,
                    "stokes": "I",
                    "detectors": [1, 2, 3],
                    "calibration_version": ver,
                },
            })
        if products:
            out[obsid] = products
    return out


_SCANNERS = {
    "ixpe": _scan_ixpe,
    "nicer": _scan_nicer,
    "nustar": _scan_nustar,
}
