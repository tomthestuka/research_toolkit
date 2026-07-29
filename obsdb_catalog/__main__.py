"""Command-line entry for obsdb_catalog — deliberate human/operator maintenance.

This is **not** an agent-callable tool: the package is not ``*_utils``-named and exposes no
``skills``/``TOOL_MANIFEST``, so it stays invisible to LLM tool discovery. This CLI is a
convenience for a person (or a maintenance script) to run the catalog operations by hand::

    python -m obsdb_catalog healthcheck <obs_db_path>   # detect-only; exits 1 on drift
    python -m obsdb_catalog repair      <obs_db_path>    # recover+reconcile ONE mission (writes)
    python -m obsdb_catalog repair_all  <root>           # repair EVERY mission under a root (writes)

``obs_db_path`` is one mission's directory (``<root>/<MISSION>``); ``root`` is the parent
that holds the mission subdirectories. ``repair``/``repair_all`` write the catalog and are the
only write-side recovery — invoke them deliberately, never automatically.

The repair commands narrate **every** change to stdout as it happens (``-`` dropped, ``+``
added, ``~`` changed, ``=`` tally, ``!`` trouble, each prefixed with its mission), then print
the closing JSON summary. The summary alone cannot distinguish a large cleanup from a no-op,
which matters because repair is the only path that deletes catalog rows.
"""

from __future__ import annotations

import argparse
import json
import sys

from .healthcheck import healthcheck
from .repair import repair, repair_all


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="obsdb_catalog",
        description="Inspect and recover per-mission data catalogs (deliberate human use).",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_hc = sub.add_parser("healthcheck", help="detect-only drift report for one mission (no writes)")
    p_hc.add_argument("obs_db_path", help="one mission's obs-db directory (<root>/<MISSION>)")

    p_rep = sub.add_parser("repair", help="recover+reconcile ONE mission catalog (writes)")
    p_rep.add_argument("obs_db_path", help="one mission's obs-db directory (<root>/<MISSION>)")

    p_all = sub.add_parser("repair_all", help="repair EVERY mission catalog under a root (writes)")
    p_all.add_argument("root", help="parent directory holding the mission subdirectories")

    args = parser.parse_args(argv)

    if args.command == "healthcheck":
        result = healthcheck(args.obs_db_path)
        if result is True:
            print("healthy")
            return 0
        print(json.dumps(result, indent=2, sort_keys=True))  # the drift report
        return 1

    if args.command == "repair":
        catalog = repair(args.obs_db_path)
        print(json.dumps(
            {"mission": catalog.get("mission"), "entries": len(catalog.get("entries", {}))},
            indent=2, sort_keys=True,
        ))
        return 0

    if args.command == "repair_all":
        results = repair_all(args.root)
        print(json.dumps(
            {m: {"entries": len(c.get("entries", {}))} for m, c in results.items()},
            indent=2, sort_keys=True,
        ))
        return 0

    parser.error(f"unknown command: {args.command}")  # unreachable (subparsers required)
    return 2


if __name__ == "__main__":
    sys.exit(main())
