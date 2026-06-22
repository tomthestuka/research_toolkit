"""The prisoned catalog writer — pure standard library, zero domain knowledge.

This is the *only* code path that producers use to mutate a mission's
``<obs_db_path>/catalog.json``. It knows nothing about IXPE/NICER/NuSTAR file
layouts, the science stack, or the orchestration repo; the catalog is plain
JSON and that JSON shape is the entire contract.

Reliability properties (this module is hand-built and tested for them):

* **Durable atomic write** — temp file in the same dir, ``flush`` + ``fsync``,
  ``os.replace`` (atomic rename), then ``fsync`` the parent dir. A reader never
  observes a half-written file.
* **Versioned history** — each write snapshots the superseded catalog into
  ``<obs_db_path>/.catalog_history/`` and prunes to the last ``HISTORY_KEEP``.
* **Concurrency-safe** — an advisory ``fcntl.flock`` on a sidecar lock file
  serialises read-modify-write so concurrent producers never clobber each other.
* **Idempotent** — products are upserted on ``(kind, path)``; re-recording the
  same product updates it in place, never appends a duplicate.
* **Truthful** — a product is recorded only after its file is confirmed on disk;
  a producer that records an absent file fails loudly.

``max_level`` / ``levels_present`` are set authoritatively from the caller's
``level`` (callers are hardcoded, trusted libraries). ``level`` may *regress*
(e.g. 3 -> 2 after an upstream re-run): the stale higher-level products are left
in place — to be overwritten in place by the next re-reduction — but the dropped
``max_level`` is what flags them as stale to readers. Authoritative
recomputation-from-disk is :mod:`obsdb_catalog.repair`'s job, not this module's.
"""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import shutil
import tempfile
from contextlib import contextmanager
from datetime import datetime, timezone

__all__ = ["CatalogError", "read", "new", "advance"]

CATALOG_NAME = "catalog.json"
LOCK_NAME = "catalog.json.lock"
HISTORY_DIR = ".catalog_history"
HISTORY_KEEP = 20


class CatalogError(Exception):
    """Raised on a malformed catalog, a malformed product, or an absent file.

    A reader hard-errors on a malformed catalog, so we surface these loudly
    rather than writing or returning junk.
    """


# --------------------------------------------------------------------------- #
# Paths / small helpers
# --------------------------------------------------------------------------- #
def _catalog_path(obs_db_path: str) -> str:
    return os.path.join(obs_db_path, CATALOG_NAME)


def mission_from_path(obs_db_path: str) -> str:
    """Default mission name from the obs-db directory (e.g. ``.../IXPE`` -> ``ixpe``)."""
    return os.path.basename(os.path.normpath(obs_db_path)).lower()


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------- #
# Locking
# --------------------------------------------------------------------------- #
@contextmanager
def _lock(obs_db_path: str):
    """Serialise read-modify-write with an advisory exclusive flock on a sidecar.

    We lock ``catalog.json.lock`` rather than the catalog itself, because
    ``os.replace`` swaps the catalog's inode underneath us during a write.
    """
    os.makedirs(obs_db_path, exist_ok=True)
    lock_path = os.path.join(obs_db_path, LOCK_NAME)
    fd = os.open(lock_path, os.O_CREAT | os.O_RDWR, 0o644)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        yield
    finally:
        try:
            fcntl.flock(fd, fcntl.LOCK_UN)
        finally:
            os.close(fd)


# --------------------------------------------------------------------------- #
# Durable atomic write + versioned history
# --------------------------------------------------------------------------- #
def _fsync_dir(path: str) -> None:
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _atomic_dump(path: str, data: dict) -> None:
    directory = os.path.dirname(path) or "."
    fd, tmp = tempfile.mkstemp(dir=directory, prefix=".catalog-", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(data, fh, indent=2)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
        tmp = None
    finally:
        if tmp is not None and os.path.exists(tmp):
            os.unlink(tmp)
    _fsync_dir(directory)


def _snapshot_history(obs_db_path: str, catalog_path: str) -> None:
    hist = os.path.join(obs_db_path, HISTORY_DIR)
    os.makedirs(hist, exist_ok=True)
    base = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%f")
    dst = os.path.join(hist, f"catalog-{base}.json")
    n = 1
    while os.path.exists(dst):  # avoid same-microsecond collision
        dst = os.path.join(hist, f"catalog-{base}-{n}.json")
        n += 1
    shutil.copyfile(catalog_path, dst)
    _prune_history(hist)


def _prune_history(hist: str) -> None:
    snaps = sorted(
        f for f in os.listdir(hist) if f.startswith("catalog-") and f.endswith(".json")
    )
    for stale in snaps[:-HISTORY_KEEP]:
        try:
            os.unlink(os.path.join(hist, stale))
        except FileNotFoundError:
            pass


def _commit(obs_db_path: str, catalog: dict) -> None:
    """Snapshot the existing catalog to history, then atomically write the new one.

    Caller must hold :func:`_lock`. Reused by :mod:`obsdb_catalog.repair`.
    """
    path = _catalog_path(obs_db_path)
    if os.path.exists(path):
        _snapshot_history(obs_db_path, path)
    _atomic_dump(path, catalog)


# --------------------------------------------------------------------------- #
# Read
# --------------------------------------------------------------------------- #
def read(obs_db_path: str) -> dict:
    """Return the full validated catalog for ``obs_db_path``.

    Returns ``{"mission": <derived>, "entries": {}}`` if the file is absent.
    Raises :class:`CatalogError` if the file is present but malformed — a torn
    or corrupt catalog must fail loudly, never parse to junk. Does **not** fetch
    sub-items; callers index ``catalog["entries"][obsid]``.
    """
    obs_db_path = os.fspath(obs_db_path)
    path = _catalog_path(obs_db_path)
    mission = mission_from_path(obs_db_path)
    if not os.path.exists(path):
        return {"mission": mission, "entries": {}}
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise CatalogError(f"malformed catalog at {path}: {exc}") from exc
    if not isinstance(data, dict) or not isinstance(data.get("entries"), dict):
        raise CatalogError(f"malformed catalog at {path}: missing 'entries' object")
    data.setdefault("mission", mission)
    return data


# --------------------------------------------------------------------------- #
# Product preparation + upsert
# --------------------------------------------------------------------------- #
def _prepare_product(obs_db_path: str, product: dict, compute_sha256: bool) -> dict:
    if not isinstance(product, dict):
        raise CatalogError(f"product must be a dict, got {type(product).__name__}")
    for key in ("kind", "path", "produced_by"):
        if not product.get(key):
            raise CatalogError(f"product missing required field {key!r}: {product!r}")
    path = product["path"]
    if os.path.isabs(path):
        raise CatalogError(f"product path must be relative to obs_db_path: {path!r}")
    full = os.path.join(obs_db_path, path)
    if not os.path.exists(full):
        # Truthful: never record a product whose file is not on disk.
        raise CatalogError(f"refusing to record absent product file: {full}")

    out = {
        "kind": product["kind"],
        "path": path,
        "produced_by": product["produced_by"],
        "params": dict(product.get("params") or {}),
        "created_at": product.get("created_at") or _now_iso(),
    }
    sha = product.get("sha256")
    if compute_sha256 and not sha:
        sha = _sha256(full)
    if sha:
        out["sha256"] = sha
    return out


def _upsert_product(products: list, product: dict) -> None:
    """Insert or replace ``product`` in ``products``, matched on ``(kind, path)``."""
    key = (product["kind"], product["path"])
    for i, existing in enumerate(products):
        if (existing.get("kind"), existing.get("path")) == key:
            products[i] = product
            return
    products.append(product)


# --------------------------------------------------------------------------- #
# Writers
# --------------------------------------------------------------------------- #
def _record(obs_db_path, mission, obsid, *, level, products, compute_sha256):
    obs_db_path = os.fspath(obs_db_path)
    if not isinstance(level, int) or level < 1:
        raise CatalogError(f"level must be a positive int, got {level!r}")
    if not obsid:
        raise CatalogError("obsid is required")
    prepared = [_prepare_product(obs_db_path, p, compute_sha256) for p in products]

    with _lock(obs_db_path):
        catalog = read(obs_db_path)
        catalog["mission"] = (mission or catalog.get("mission") or
                              mission_from_path(obs_db_path))
        entries = catalog.setdefault("entries", {})
        entry = entries.get(obsid)
        if entry is None:
            entry = {"obsid": obsid, "max_level": level,
                     "levels_present": [], "products": []}
            entries[obsid] = entry
        entry["obsid"] = obsid
        prods = entry.setdefault("products", [])
        for prepared_product in prepared:
            _upsert_product(prods, prepared_product)
        # `level` is authoritative and may regress; stale higher-level products
        # are left in place (overwritten by the next re-reduction).
        entry["max_level"] = level
        entry["levels_present"] = list(range(1, level + 1))
        _commit(obs_db_path, catalog)
        return entry


def advance(obs_db_path, mission, obsid, *, level, products, compute_sha256=False):
    """Record a reduction's output for ``obsid`` at target ``level``.

    Idempotently upserts each product on ``(kind, path)`` and sets the entry's
    ``max_level`` to ``level`` (which may be lower than before — that regression
    is the signal that higher-level products are now stale). Each product file
    must already exist on disk. Returns the updated entry.

    ``products`` is an iterable of dicts: required ``kind``, ``path`` (relative
    to ``obs_db_path``), ``produced_by``; optional ``params`` (dict),
    ``created_at`` (else stamped now), ``sha256`` (else computed iff
    ``compute_sha256``).
    """
    return _record(obs_db_path, mission, obsid, level=level,
                   products=products, compute_sha256=compute_sha256)


def new(obs_db_path, mission, obsid, *, products, compute_sha256=False):
    """Record a freshly acquired (downloaded) ObsID — level-2 events.

    Thin lifecycle wrapper over :func:`advance` with ``level=2``. Re-acquiring an
    ObsID that was at level 3 therefore regresses ``max_level`` to 2, flagging
    its stale level-3 products for re-reduction.
    """
    return _record(obs_db_path, mission, obsid, level=2,
                   products=products, compute_sha256=compute_sha256)
