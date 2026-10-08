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
            "Acquire (download) HEASARC archive data for one IXPE, NICER or NuSTAR "
            "observation into a mission's obs-db directory and record it in the data "
            "catalog. Call this when an ObsID is genuinely missing locally."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "METHOD: resolves the ObsID's archive location authoritatively via the "
                "mission master table (rather than any hand-built path), downloads it, and "
                "records the level-2 event at processing level 2. Data lands at "
                "<obs_db_path>/<obsid>/.\n"
                "RE-DOWNLOADS ON EVERY CALL: reserve it for genuinely missing observations -- "
                "confirm locality with a data-catalog tool first. A re-acquire regresses the "
                "catalog to level 2, which flags any stale level-3 products for re-reduction.\n"
                "ENCRYPTED (PROPRIETARY) DATA: an observation still inside its proprietary "
                "period is served GPG-encrypted. That download is complete but unreadable, and "
                "it is recorded as present-but-locked rather than reported as a failure, so the "
                "data-catalog listing shows it as LOCKED, NOT USABLE. Asking for it again is "
                "answered from the archive's public_date without re-transferring: while it is "
                "still proprietary the download is skipped and the reply says when it becomes "
                "public ('downloaded': false, 'locked': true). Treat a locked observation as "
                "unavailable and report that honestly.\n"
                "SUPPORTED MISSIONS: ixpe, nicer, nustar.\n"
                "TIMEOUTS: this tool overrides the astroquery factory default of 30 s, which "
                "is too short for archive downloads, with a 1800 s default; for large IXPE "
                "observations (> 1 GB) that still time out, pass a longer 'timeout_s' (e.g. "
                "1800 or 3600).\n"
                "FAILS loudly, recording nothing, when the ObsID is absent from the archive "
                "or the download yields no level-2 file."
            ),
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
