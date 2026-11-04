"""Where an observation's folder is — the one resolver every reader of an obs-db goes through.

An obs-db root holds one folder per mission, named as HEASARC names the mission
(``<obs_db_root>/ixpe/``), and below it each ObsID's folder where HEASARC's archive keeps it
(``ixpe/data/obs/01/01002901/``; :data:`obsdb_catalog._domain.TREES`). Catalog product paths are
relative to the mission folder. Pure standard library.
"""

from __future__ import annotations

import glob
import os
from pathlib import Path

from . import _domain, neutral

__all__ = ["mission_dir", "locate"]


def mission_dir(obs_db_root, mission: str) -> Path:
    """The mission's folder under ``obs_db_root``: an existing child whose name matches ``mission``
    case-insensitively, else ``<obs_db_root>/<mission, lowercase>`` (HEASARC's spelling)."""
    root, want = Path(obs_db_root), mission.lower()
    if root.is_dir():
        for child in sorted(root.iterdir()):
            if child.is_dir() and child.name.lower() == want:
                return child
    return root / want


def locate(obs_db_root, mission: str, obsid: str) -> Path:
    """The folder holding ``obsid``'s data under ``<obs_db_root>/<mission>/``.

    The catalog's record comes first: the folder a recorded product of the entry lies in, up to
    the path part that is the ObsID (its newest level-2 event before other kinds; a path leaving
    the mission folder is passed over). Without one, the
    mission's tree rule (:data:`_domain.TREES`) is searched. Raises :class:`FileNotFoundError`
    naming where it looked, and :class:`LookupError` when the tree rule finds the ObsID twice.
    """
    folder, obsid = mission_dir(obs_db_root, mission), str(obsid)
    try:
        entry = neutral.read(folder)["entries"].get(obsid) or {}
    except neutral.CatalogError:
        entry = {}
    # The newest level-2 event first (a re-acquire into another layout wins), then the rest.
    recorded = [p or {} for p in entry.get("products") or []]
    level2 = sorted((p for p in recorded if p.get("kind") == "level2_evt"),
                    key=lambda p: str(p.get("created_at") or ""), reverse=True)
    products = level2 + [p for p in recorded if p.get("kind") != "level2_evt"]
    for product in products:
        parts = Path(product.get("path") or "").parts
        if obsid in parts and ".." not in parts and not os.path.isabs(product.get("path")):
            found = folder.joinpath(*parts[: parts.index(obsid) + 1])
            if found.is_dir():
                return found
    patterns = _domain.TREES.get(mission.lower(), ("{obsid}",))
    found = sorted({d for pattern in patterns
                    for d in glob.glob(str(folder / pattern.format(obsid=glob.escape(obsid))))
                    if os.path.isdir(d)})
    if len(found) > 1:
        raise LookupError(f"{mission} ObsID {obsid} is in more than one folder under {folder}: {found}")
    if not found:
        raise FileNotFoundError(f"no folder for {mission} ObsID {obsid} under {folder}: its catalog "
                                f"records none, and none matches {', '.join(patterns)}")
    return Path(found[0])
