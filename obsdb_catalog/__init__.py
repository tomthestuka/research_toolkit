"""obsdb_catalog — the single trusted writer for per-mission data catalogs.

A standalone, dependency-free package: the *only* code path permitted to mutate
a mission's ``<obs_db_path>/catalog.json`` (the inventory of what data is local
and how it was produced). Producers that land data on disk record it by calling
:func:`obsdb_catalog.neutral.advance` / :func:`obsdb_catalog.neutral.new`;
nothing else opens the catalog for writing.

Layers:

* :mod:`obsdb_catalog.neutral` — the prisoned writer. Pure standard library,
  zero domain knowledge. Every producer calls this.
* :mod:`obsdb_catalog.repair` — deliberate recovery (restore + reconcile). May
  use mission layout knowledge; persists through ``neutral``'s primitive.
* :mod:`obsdb_catalog.healthcheck` — detect-only drift report; never writes.

This package is intentionally *not* an LLM tool: no ``skills/``, no
``TOOL_MANIFEST``, no CLI, and not ``*_utils``-named (so it is invisible to the
orchestration repo's tool discovery). Catalog writes are a side effect of
data-producing libraries, never a model action.
"""

from . import neutral
from .neutral import CatalogError, advance, new, read

__all__ = ["neutral", "read", "new", "advance", "CatalogError"]
