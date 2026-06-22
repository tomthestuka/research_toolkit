"""heasarc_utils.skills — AI-dedicated tool interface for the auto_research sandbox.

Exposes TOOL_MANIFEST for ToolkitRegistry discovery. Pure-data module: no
astroquery/science imports here, so the registry can load it in the orchestration
environment. The astroquery-dependent code lives in lib/ (imported lazily, only
when a download actually runs in the REPL subprocess).
"""

from __future__ import annotations

TOOL_MANIFEST = {
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
            "ixpe, nicer, nustar. "
            "For large IXPE observations (> 1 GB) that time out, pass a longer "
            "timeout_s (e.g. 1800 or 3600); the default is already 1800 s. "
            "The astroquery factory default (30 s) was too short and is always "
            "overridden by this tool."
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
                "timeout_s": {
                    "type": "number",
                    "description": (
                        "HTTP transfer timeout in seconds for the archive download. "
                        "Default 1800 s (30 min). The astroquery factory default "
                        "(30 s) is too short for large IXPE observations (> 1 GB) "
                        "and is always replaced by this value. Increase to 3600 s "
                        "or more if transfers are still cut off on slow connections. "
                        "Omitting this argument uses the 1800 s default."
                    ),
                },
            },
            "required": ["obsid", "mission", "obs_db_path"],
        },
    },
}
