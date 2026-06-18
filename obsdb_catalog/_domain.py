"""Mission layout knowledge — used only by :mod:`repair` and :mod:`healthcheck`.

:mod:`obsdb_catalog.neutral` is kept domain-free; all per-mission knowledge
(directory layout, filename conventions, which ``kind`` is which level, the
constant reduction params) lives here. This module is plain standard library
(glob/re/os); it does **not** import the science stack or ``auto_research``.

The conventions mirror the live bootstrapped catalogs:

* **IXPE** ``<obsid>/event_l2/`` — one representative ``det1`` entry per product
  kind. ``level2_evt`` (``produced_by:"heasarc"``), ``pcube`` / ``pha``
  (``produced_by:"filter_inst_bkg_IXPE.sh"``). ``filter`` is read from the
  filename infix (``_rej_src_`` -> ``"rej_src"``; ``_src_`` -> ``"src"``, the
  ``03002599`` special case); ``calibration_version`` from the ``v0N`` token.
* **NICER** ``<obsid>/xti/event_cl/`` — ``pha`` / ``pha_bkg`` / ``arf`` / ``rmf``
  (``produced_by:"nicer_red.sh"``), constant params ``background_model:"3C50"``,
  ``mpu:7``.
* **NuSTAR** ``<obsid>/products/`` — per FPM A/B: ``pha_grp`` (``nustar_grp.sh``,
  ``grouping:"optimal"``), ``pha_src`` / ``arf`` / ``rmf`` (``nustar_red.sh``).
"""

from __future__ import annotations

import glob
import os
import re

#: Marker stamped into a param whose value cannot be re-derived during repair.
UNRECOVERABLE = "repaired session, unrecoverable"

_LEVEL2_KINDS = {"level2_evt", "raw_evt"}
_LEVEL1_KINDS = {"raw", "level1"}


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
    for evt_cl in sorted(glob.glob(os.path.join(obs_db_path, "*", "xti", "event_cl"))):
        obsid = evt_cl.split(os.sep)[-3]
        products = []
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
    for prod_dir in sorted(glob.glob(os.path.join(obs_db_path, "*", "products"))):
        obsid = prod_dir.split(os.sep)[-2]
        products = []
        # kind suffix, glob, produced_by, extra params (besides fpm)
        specs = [
            ("pha_grp", "*{fpm}01_sr_grp.pha", "nustar_grp.sh", {"grouping": "optimal"}),
            ("pha_src", "*{fpm}01_sr.pha", "nustar_red.sh", {}),
            ("arf", "*{fpm}01_sr.arf", "nustar_red.sh", {}),
            ("rmf", "*{fpm}01_sr.rmf", "nustar_red.sh", {}),
        ]
        for fpm in ("A", "B"):
            for kind, pat, produced_by, extra in specs:
                matches = sorted(glob.glob(os.path.join(prod_dir, pat.format(fpm=fpm))))
                for full in matches:
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
_IXPE_EVT = re.compile(r"_det1_evt2_(?P<ver>v\d+)\.fits$")
# pcube/pha derived names carry the filter infix and calibration version.
_IXPE_DERIVED = re.compile(
    r"_det1_evt2_(?P<ver>v\d+)_(?P<filter>rej_src|src)_(?P<prod>pcube|pha1)\.fits$"
)
_IXPE_NOTE = "det1 listed; det2/det3 follow same naming pattern"


def _scan_ixpe(obs_db_path: str) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for evt_l2 in sorted(glob.glob(os.path.join(obs_db_path, "*", "event_l2"))):
        obsid = evt_l2.split(os.sep)[-2]
        names = sorted(
            os.path.basename(f)
            for f in glob.glob(os.path.join(evt_l2, "ixpe*_det1_evt2_*.fits"))
        )
        derived = [(n, m) for n in names if (m := _IXPE_DERIVED.search(n))]
        # Per the bootstrap convention, the catalogued products use the
        # instrumental-background-rejected filter when it exists, falling back to
        # the non-rejected "src" filter (the 03002599 special case). The other
        # variant's files are intentionally left out of the catalog.
        chosen = "rej_src" if any(m.group("filter") == "rej_src" for _, m in derived) \
            else "src"

        products = []
        for name in names:
            m_evt = _IXPE_EVT.search(name)
            if m_evt:
                products.append({
                    "kind": "level2_evt",
                    "path": _rel(obs_db_path, os.path.join(evt_l2, name)),
                    "produced_by": "heasarc",
                    "params": {
                        "detectors": [1, 2, 3],
                        "calibration_version": m_evt.group("ver"),
                        "note": _IXPE_NOTE,
                    },
                })
        for name, m in derived:
            if m.group("filter") != chosen:
                continue
            ver, prod = m.group("ver"), m.group("prod")
            rel = _rel(obs_db_path, os.path.join(evt_l2, name))
            if prod == "pcube":
                products.append({
                    "kind": "pcube",
                    "path": rel,
                    "produced_by": "filter_inst_bkg_IXPE.sh",
                    "params": {
                        "filter": chosen,
                        "detectors": [1, 2, 3],
                        "calibration_version": ver,
                        "note": _IXPE_NOTE,
                    },
                })
            else:  # pha1
                products.append({
                    "kind": "pha",
                    "path": rel,
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
