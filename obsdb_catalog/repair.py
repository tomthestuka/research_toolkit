"""Deliberate catalog recovery — the only write-side recovery path.

Follows the prevent-then-restore consensus: durable atomic writes
(:mod:`neutral`) make corruption rare; recovery restores the last good state
and reconciles it against disk. It never hand-fixes corrupt JSON bytes and never
fabricates data values — anything it cannot re-derive is flagged with
``_domain.UNRECOVERABLE``.

Algorithm:

1. ``read`` the catalog. If it does not parse, restore the newest *parseable*
   snapshot from ``.catalog_history/`` (walking newest -> oldest).
2. Reconcile each entry against disk: drop products whose files are gone (and
   drop entries left empty), then recompute ``max_level`` / ``levels_present``
   from the surviving products' kinds.
3. Re-derive products that are present on disk (per the mission layout) but
   absent from the catalog, additively — existing entries' richer params are
   preserved, only genuinely-missing products are added.
4. Write the result atomically through :mod:`neutral`'s shared primitive.

Invoke deliberately (an operator/agent, or a maintenance run) — never wire it to
run automatically.
"""

from __future__ import annotations

import glob
import json
import os

from . import _domain, neutral

__all__ = ["repair", "repair_all"]


def _newest_parseable_snapshot(obs_db_path: str):
    hist = os.path.join(obs_db_path, neutral.HISTORY_DIR)
    snaps = sorted(
        glob.glob(os.path.join(hist, "catalog-*.json")), reverse=True
    )
    for snap in snaps:
        try:
            with open(snap, encoding="utf-8") as fh:
                data = json.load(fh)
        except (OSError, json.JSONDecodeError, UnicodeDecodeError):
            continue
        if isinstance(data, dict) and isinstance(data.get("entries"), dict):
            return data
    return None


def _recompute_levels(entry: dict) -> None:
    products = entry.get("products", [])
    if not products:
        entry["max_level"] = 0
        entry["levels_present"] = []
        return
    max_level = max(_domain.kind_level(p.get("kind", "")) for p in products)
    entry["max_level"] = max_level
    entry["levels_present"] = list(range(1, max_level + 1))


def repair(obs_db_path, mission: str | None = None) -> dict:
    """Recover and reconcile the catalog at ``obs_db_path``; returns the result."""
    obs_db_path = os.fspath(obs_db_path)
    mission = (mission or _domain.mission_from_path(obs_db_path)).lower()

    with neutral._lock(obs_db_path):
        # 1. read, or restore the newest parseable snapshot.
        try:
            catalog = neutral.read(obs_db_path)
        except neutral.CatalogError:
            catalog = _newest_parseable_snapshot(obs_db_path) or {
                "mission": mission, "entries": {}
            }
        catalog["mission"] = mission
        entries = catalog.setdefault("entries", {})

        # 2. reconcile existing entries against disk.
        for obsid in list(entries):
            entry = entries[obsid]
            entry["obsid"] = obsid
            kept = [
                p for p in entry.get("products", [])
                if os.path.exists(os.path.join(obs_db_path, p.get("path", "")))
            ]
            if not kept:
                del entries[obsid]
                continue
            entry["products"] = kept
            _recompute_levels(entry)

        # 3. additively re-derive products present on disk but uncatalogued.
        for obsid, derived in _domain.scan(mission, obs_db_path).items():
            entry = entries.get(obsid)
            if entry is None:
                entry = {"obsid": obsid, "max_level": 0,
                         "levels_present": [], "products": []}
                entries[obsid] = entry
            present = {(p.get("kind"), p.get("path")) for p in entry["products"]}
            added = False
            for product in derived:
                if (product["kind"], product["path"]) not in present:
                    entry["products"].append(product)
                    added = True
            if added:
                _recompute_levels(entry)

        neutral._commit(obs_db_path, catalog)
        return catalog


def repair_all(root) -> dict:
    """Repair every mission catalog under a single ``root`` directory.

    Discovers missions the way the orchestration does: each immediate subdirectory of ``root``
    that holds a ``catalog.json`` is a mission (keyed by its lowercased basename); unrelated
    working directories under ``root`` (no ``catalog.json``) are ignored. Each discovered
    mission is repaired via :func:`repair`. Returns ``{mission: repaired_catalog}``.
    """
    root = os.fspath(root)
    results: dict[str, dict] = {}
    for name in sorted(os.listdir(root)):
        sub = os.path.join(root, name)
        if os.path.isdir(sub) and os.path.exists(os.path.join(sub, neutral.CATALOG_NAME)):
            results[name.lower()] = repair(sub)
    return results
