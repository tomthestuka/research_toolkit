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
            "observation into the obs-db root, laid out as HEASARC's archive is, and record "
            "it in that mission's data catalog. Call this when an ObsID is genuinely missing "
            "locally."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "METHOD: resolves the ObsID's archive location authoritatively via the "
                "mission master table (rather than any hand-built path), downloads it, and "
                "records the level-2 event at processing level 2. Data lands where HEASARC's "
                "archive keeps it, below the mission's folder: <obs_db_root>/<mission>/<archive "
                "tree>/<obsid>/ (e.g. ixpe/data/obs/01/01002901/), recorded in "
                "<obs_db_root>/<mission>/catalog.json; the first acquire of a mission creates both. "
                "The result names mission_dir and obsid_dir.\n"
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
                    "description": "Mission name (lowercase).",
                },
                "obs_db_root": {
                    "kind": "path",
                    "type": "string",
                    "description": (
                        "The obs-db root: one folder per mission, named as HEASARC names it "
                        "(ixpe, nicer, nustar), each laid out as HEASARC's archive tree. A new, "
                        "empty folder becomes one on the first acquire; a mission folder laid out "
                        "flat (<mission>/<obsid>/) is refused."
                    ),
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
            "required": ["obsid", "mission", "obs_db_root"],
        },
    },
}
