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
        "writes": ["db"],
        "algorithm_version": "1",
        "description": (
            "Acquire (download) HEASARC archive data for one IXPE, NICER or NuSTAR observation into "
            "the data store, and return its record. Call this when the data catalog lists no raw "
            "product for the ObsID, or to fetch the archive's newest processing of it."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "METHOD: resolves the ObsID in the mission's master table (rather than any hand-built "
                "path), downloads the observation, and stores it as one product of kind raw in the "
                "mission's bucket, its files laid out as the ObsID folder of HEASARC's archive "
                "(event_l2/..., xti/event_cl/..., event_cl/...). The result names the product's "
                "'record' (pass it to tools that read the observation), the 'level2_evt' record of its "
                "level-2 event, and 'downloaded'.\n"
                "ONE DOWNLOAD PER ARCHIVE STATE: the archive's processing date, software and CALDB "
                "versions and public status identify the download. When the store already holds that "
                "download, the call returns its record at once ('downloaded': false); a reprocessed "
                "observation, or one that has become public, is downloaded afresh.\n"
                "ENCRYPTED (PROPRIETARY) DATA: an observation still inside its proprietary period is "
                "served GPG-encrypted. That download is complete but unreadable; it is stored and the "
                "result says 'locked': true with the days until it becomes public. Treat a locked "
                "observation as unavailable and report that honestly.\n"
                "BUSY: when another unit holds the same download, the call says so and names when its "
                "claim ends; work on other data meanwhile and call again after that time.\n"
                "SUPPORTED MISSIONS: ixpe, nicer, nustar.\n"
                "TIMEOUTS: this tool overrides the astroquery factory default of 30 s with a 1800 s "
                "default; for large IXPE observations (> 1 GB) that still time out, pass a longer "
                "'timeout_s' (e.g. 3600).\n"
                "FAILS loudly, storing nothing, when the ObsID is absent from the archive, the archive "
                "states no processing date for it, the download yields no level-2 file, the archive still "
                "serves an older processing than its catalog names (the level-2 header's SOFTVER or CALDBVER "
                "differs from the row; call again later), or a public observation is still served "
                "encrypted (call again later)."
            ),
            "properties": {
                "obsid": {
                    "kind": "setting",
                    "type": "string",
                    "description": (
                        "Observation ID in the mission's native format, e.g. "
                        "'01002901' (IXPE), '4690010101' (NICER), '30702017002' "
                        "(NuSTAR)."
                    ),
                },
                "mission": {
                    "kind": "setting",
                    "type": "string",
                    "enum": ["ixpe", "nicer", "nustar"],
                    "description": "Mission name (lowercase); also the store bucket the download goes to.",
                },
                "host": {
                    "kind": "setting",
                    "type": "string",
                    "enum": ["heasarc", "sciserver", "aws"],
                    "description": (
                        "Optional download host. Omit to let astroquery choose "
                        "based on the environment."
                    ),
                },
                "timeout_s": {
                    "kind": "setting",
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
            "required": ["obsid", "mission"],
        },
    },
}
