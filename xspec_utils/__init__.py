"""xspec_utils — XSPEC spectral analysis tools for the auto_research sandbox.

Domain code lives in lib/ (xspec_help, NICER, NuSTAR, red_ixpe, find_src,
src_pos) — the editable working copy a tool author may extend or refactor.

The AI-dedicated tool interface lives in skills/:
  skills/__init__.py  — TOOL_MANIFEST (pure data; no xspec import)
  skills/session.py   — xspec-dependent wrappers (REPL subprocess only)

This package re-exports TOOL_MANIFEST for ToolkitRegistry discovery. The lib/
sys.path insertion lets the REPL subprocess import xspec_help etc. as bare names.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add lib/ to sys.path so the REPL subprocess can import xspec_help,
# NICER, NuSTAR, red_ixpe, etc. as bare module names.
_lib = Path(__file__).parent / "lib"
if str(_lib) not in sys.path:
    sys.path.insert(0, str(_lib))

from .skills import TOOL_MANIFEST  # noqa: E402

__all__ = ["TOOL_MANIFEST"]
