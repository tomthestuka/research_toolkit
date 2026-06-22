"""heasarc_utils — HEASARC archive *acquire* tools for the auto_research sandbox.

The acquire engine and per-mission profiles live in lib/ (human-authored,
importable domain code). The AI-dedicated tool interface lives in skills/:
  skills/__init__.py  — TOOL_MANIFEST (pure data; no astroquery import)

This package re-exports TOOL_MANIFEST for ToolkitRegistry discovery. astroquery
is imported lazily inside lib/heasarc.py, only when a download actually runs.

It records downloaded data through obsdb_catalog (the single trusted catalog
writer) — the first producer wired into that catalog.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Mirror xspec_utils: put lib/ on sys.path so the REPL subprocess can import the
# domain modules, and keep heavy imports out of this file.
_lib = Path(__file__).parent / "lib"
if str(_lib) not in sys.path:
    sys.path.insert(0, str(_lib))

from .skills import TOOL_MANIFEST  # noqa: E402

__all__ = ["TOOL_MANIFEST"]
