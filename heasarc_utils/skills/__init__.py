"""heasarc_utils.skills — AI-dedicated tool interface for the auto_research sandbox.

Exposes TOOL_MANIFEST for ToolkitRegistry discovery. Pure-data module: no
astroquery/science imports here, so the registry can load it in the orchestration
environment. The astroquery-dependent code lives in lib/ (imported lazily, only
when a download actually runs in the REPL subprocess).
"""

from __future__ import annotations

TOOL_MANIFEST: dict = {
    "research_toolkit.heasarc_utils.acquire": {
        "module": "heasarc_utils.lib.heasarc_help",
        "fn": "acquire",
        "description": (
            "Acquire (download) HEASARC archive data for one observation into a "
            "mission's obs-db directory and record it in the data catalog. "
            "Resolves the ObsID's archive location authoritatively via the mission "
            "master table (no hand-built paths), downloads it, and records the "
            "level-2 event at processing level 2. Data lands at "
            "<obs_db_path>/<obsid>/. ALWAYS (re)downloads — call it only when data "
            "is genuinely missing; a re-acquire regresses the catalog to level 2, "
            "flagging any stale level-3 products for re-reduction. Fails loudly if "
            "the ObsID is not found or the download does not yield the expected "
            "level-2 file (nothing is recorded on failure). Supported missions: "
            "ixpe, nicer, nustar."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "obsid": {
                    "type": "string",
                    "description": (
                        "Observation ID in the mission's native format, e.g. "
                        "'01002901' (IXPE), '4690010101' (NICER), '30702017002' "
                        "(NuSTAR)."
                    ),
                },
                "mission": {
                    "type": "string",
                    "enum": ["ixpe", "nicer", "nustar"],
                    "description": "Mission name (lowercase).",
                },
                "obs_db_path": {
                    "type": "string",
                    "description": (
                        "The mission's obs-db directory, e.g. "
                        "'/data/fzhou/XSPEC/Cyg_X-1/IXPE'. Data is downloaded to "
                        "<obs_db_path>/<obsid>/ and recorded in "
                        "<obs_db_path>/catalog.json."
                    ),
                },
                "host": {
                    "type": "string",
                    "enum": ["heasarc", "sciserver", "aws"],
                    "description": (
                        "Optional download host. Omit to let astroquery choose "
                        "based on the environment."
                    ),
                },
            },
            "required": ["obsid", "mission", "obs_db_path"],
        },
    },
}
