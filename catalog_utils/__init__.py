"""catalog_utils — the AI-facing read of the data catalog.

Renders the store's catalog (Postgres, read with the unit's own catalog login) into the listing analysts
consult BEFORE acquiring, so they acquire only genuinely missing data, and from which they take the
records analysis tools read.

Split (CLAUDE.md §5):
  skills/__init__.py  — TOOL_MANIFEST (pure data; no heavy imports)
  lib/catalog_read.py — the query and the listing (``storelib`` imported when it runs)

This package re-exports TOOL_MANIFEST for ToolkitRegistry discovery.
"""

from __future__ import annotations

from .skills import TOOL_MANIFEST

__all__ = ["TOOL_MANIFEST"]
