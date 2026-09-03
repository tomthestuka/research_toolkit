"""catalog_utils — the AI-facing READ of the local data catalog.

The domain read-half of the per-mission data inventory. ``obsdb_catalog`` (its
sibling package) is the single trusted *writer* and low-level reader of a mission's
``<obs_db_path>/catalog.json``; this package renders that catalog into the listing
analysts consult BEFORE acquiring, so they acquire only genuinely-missing data.

Split (the two-halves toolkit layout):
  skills/__init__.py  — TOOL_MANIFEST (pure data; no heavy imports)
  lib/catalog_read.py — the read + listing logic (pure standard library, reusing
                        ``obsdb_catalog.neutral.read``)

Pure standard library (no science stack), so the orchestration layer can both
discover the manifest and invoke the read in-process. The orchestration layer owns
*where* data lives (it resolves each mission's ``obs_db_path`` from its own config
and passes it in); this package owns *how* the catalog is read and presented.

This package re-exports TOOL_MANIFEST for ToolkitRegistry discovery.
"""

from __future__ import annotations

from .skills import TOOL_MANIFEST

__all__ = ["TOOL_MANIFEST"]
