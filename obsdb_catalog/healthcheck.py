"""Detect-only catalog health check. Never writes.

Returns ``True`` when the catalog is healthy (the caller does nothing).
Otherwise returns a structured report (a non-empty dict, hence truthy) listing
what is broken. Use :mod:`obsdb_catalog.repair` to fix what this finds.

Checks:

* **malformed** — the catalog file is present but does not parse / has the wrong
  shape (delegated to :func:`neutral.read`, which raises loudly).
* **missing** — a catalogued ``product.path`` is not on disk.
* **level_inconsistent** — ``levels_present`` is not ``[1..max_level]``.
* **regressed** — a simple flag: ``max_level < 3`` while the entry still lists
  level-3 products (a re-reduction is pending or was forgotten).
* **uncatalogued** — files present on disk (per the mission layout) that are not
  in the catalog.
"""

from __future__ import annotations

import os

from . import _domain, neutral

__all__ = ["healthcheck"]


def healthcheck(obs_db_path):
    """Return ``True`` if healthy, else a report dict of what is broken."""
    obs_db_path = os.fspath(obs_db_path)
    mission = _domain.mission_from_path(obs_db_path)

    try:
        catalog = neutral.read(obs_db_path)
    except neutral.CatalogError as exc:
        return {"malformed": str(exc)}

    entries = catalog.get("entries", {})
    missing: list[dict] = []
    level_inconsistent: list[dict] = []
    regressed: list[str] = []

    catalogued_paths: set[str] = set()
    for obsid, entry in entries.items():
        products = entry.get("products", [])
        max_level = entry.get("max_level")
        for product in products:
            path = product.get("path", "")
            catalogued_paths.add(path)
            if not os.path.exists(os.path.join(obs_db_path, path)):
                missing.append({"obsid": obsid, "kind": product.get("kind"), "path": path})

        if isinstance(max_level, int):
            if entry.get("levels_present") != list(range(1, max_level + 1)):
                level_inconsistent.append({
                    "obsid": obsid,
                    "max_level": max_level,
                    "levels_present": entry.get("levels_present"),
                })
            # Simple regressed flag: still carrying level-3 products below 3.
            if max_level < 3 and any(
                _domain.kind_level(p.get("kind", "")) == 3 for p in products
            ):
                regressed.append(obsid)

    uncatalogued: list[str] = []
    for obsid, products in _domain.scan(mission, obs_db_path).items():
        for product in products:
            if product["path"] not in catalogued_paths:
                uncatalogued.append(product["path"])

    report = {}
    if missing:
        report["missing"] = missing
    if level_inconsistent:
        report["level_inconsistent"] = level_inconsistent
    if regressed:
        report["regressed"] = regressed
    if uncatalogued:
        report["uncatalogued"] = sorted(uncatalogued)

    return report or True
