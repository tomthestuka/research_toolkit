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

Every change is reported on stdout as it happens (``-`` dropped, ``+`` added,
``~`` changed, ``=`` tally, ``!`` trouble), each line prefixed with its mission.
Because this is the one path that deletes catalog rows, the operator needs to
see *what* went rather than only how many entries survived. Pass ``log=None`` to
silence, or a list's ``append`` to capture.

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


def repair(obs_db_path, mission: str | None = None, *, log=print) -> dict:
    """Recover and reconcile the catalog at ``obs_db_path``; returns the result.

    Reports **every** change on stdout as it is made — one line per dropped product, dropped
    entry, re-derived product and level change, plus a closing per-mission tally. Repair is a
    deliberate write against data a person cares about, so a bare ``entries: 47`` leaves the
    operator unable to tell a 97-row cleanup from a no-op; these lines are the audit trail.

    ``log`` takes any ``(str) -> None`` callable and defaults to :func:`print`; pass
    ``log=None`` to silence it, or a list's ``append`` to capture. Lines carry the mission
    prefix so :func:`repair_all`'s interleaved output stays attributable, and are shaped to be
    greppable: ``-`` dropped, ``+`` added, ``~`` changed, ``=`` tally, ``!`` trouble.
    """
    obs_db_path = os.fspath(obs_db_path)
    mission = (mission or _domain.mission_from_path(obs_db_path)).lower()
    emit = log if callable(log) else (lambda _line: None)
    tag = f"[{mission}]"
    dropped = added = 0

    with neutral._lock(obs_db_path):
        # 1. read, or restore the newest parseable snapshot.
        try:
            catalog = neutral.read(obs_db_path)
        except neutral.CatalogError:
            snapshot = _newest_parseable_snapshot(obs_db_path)
            if snapshot is None:
                emit(f"{tag} ! catalog unreadable and no parseable snapshot — starting empty")
            else:
                emit(f"{tag} ~ catalog unreadable — restored newest parseable snapshot")
            catalog = snapshot or {"mission": mission, "entries": {}}
        catalog["mission"] = mission
        entries = catalog.setdefault("entries", {})
        before_products = sum(len(e.get("products", [])) for e in entries.values())

        # 2. reconcile existing entries against disk.
        for obsid in list(entries):
            entry = entries[obsid]
            entry["obsid"] = obsid
            was_level = entry.get("max_level")
            kept, gone = [], []
            for product in entry.get("products", []):
                target = os.path.join(obs_db_path, product.get("path", ""))
                (kept if os.path.exists(target) else gone).append(product)
            for product in gone:
                dropped += 1
                emit(f"{tag} - drop {obsid} {product.get('kind')} "
                     f"{product.get('path')} (file missing)")
            if not kept:
                del entries[obsid]
                emit(f"{tag} - drop entry {obsid} (no surviving products)")
                continue
            entry["products"] = kept
            _recompute_levels(entry)
            if gone and entry["max_level"] != was_level:
                emit(f"{tag} ~ level {obsid}: {was_level} -> {entry['max_level']}")

        # 3. additively re-derive products present on disk but uncatalogued.
        for obsid, derived in _domain.scan(mission, obs_db_path).items():
            entry = entries.get(obsid)
            if entry is None:
                entry = {"obsid": obsid, "max_level": 0,
                         "levels_present": [], "products": []}
                entries[obsid] = entry
                emit(f"{tag} + entry {obsid} (on disk, absent from catalog)")
            was_level = entry.get("max_level")
            present = {(p.get("kind"), p.get("path")) for p in entry["products"]}
            new = [p for p in derived if (p["kind"], p["path"]) not in present]
            for product in new:
                entry["products"].append(product)
                added += 1
                emit(f"{tag} + add {obsid} {product['kind']} "
                     f"{product['path']} (on disk, uncatalogued)")
            if new:
                _recompute_levels(entry)
                if entry["max_level"] != was_level:
                    emit(f"{tag} ~ level {obsid}: {was_level} -> {entry['max_level']}")

        after_products = sum(len(e.get("products", [])) for e in entries.values())
        if dropped or added:
            emit(f"{tag} = {len(entries)} entries, {before_products} -> {after_products} "
                 f"products (-{dropped} +{added})")
        else:
            emit(f"{tag} = no changes ({len(entries)} entries, {after_products} products)")

        neutral._commit(obs_db_path, catalog)
        return catalog


def repair_all(root, *, log=print) -> dict:
    """Repair every mission catalog under a single ``root`` directory.

    Discovers missions the way the orchestration does: each immediate subdirectory of ``root``
    that holds a ``catalog.json`` is a mission (keyed by its lowercased basename); unrelated
    working directories under ``root`` (no ``catalog.json``) are ignored. Each discovered
    mission is repaired via :func:`repair`. Returns ``{mission: repaired_catalog}``.

    ``log`` is forwarded to each :func:`repair`, so every mission's changes are reported under
    its own prefix as they happen (defaults to :func:`print`; ``None`` silences).
    """
    root = os.fspath(root)
    results: dict[str, dict] = {}
    for name in sorted(os.listdir(root)):
        sub = os.path.join(root, name)
        if os.path.isdir(sub) and os.path.exists(os.path.join(sub, neutral.CATALOG_NAME)):
            results[name.lower()] = repair(sub, log=log)
    return results
