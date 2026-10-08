"""catalog_utils.skills — AI-facing tool interface for the local data-catalog read.

Exposes TOOL_MANIFEST for ToolkitRegistry discovery. Pure-data module: no heavy
imports here (the read logic in ``lib/catalog_read.py`` is pure standard library and
imported only when the tool runs).

Like every toolkit tool, it takes the mission's ``obs_db_path`` directly (the analyst
context lists each mission's ``obs_db_path``, exactly as ``acquire``/``pcube`` receive it),
so the tool runs uniformly in the worker sandbox and — since the read is pure standard
library — in-process for the analyst lead. ``obs_db_path`` is required; ``min_level`` is
optional.
"""

from __future__ import annotations

TOOL_MANIFEST = {
    "research_toolkit.catalog_utils.data_catalog": {
        "module": "catalog_utils.lib.catalog_read",
        "fn": "data_catalog",
        "description": (
            "List observation data already on disk — which ObsIDs are local, their "
            "processing level, the products derived from them, and whether any are "
            "encrypted and unreadable. Read it BEFORE requesting data acquisition so you "
            "only acquire what is genuinely missing."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "TAKES ONE PATH: obs_db_path selects which data to list, because each "
                "observatory's data lives in its own directory (pick the one you want from "
                "your local observation databases list). The catalog read is "
                "<obs_db_path>/catalog.json; the listing names the observatory it belongs "
                "to. min_level is the only other argument.\n"
                "LEVELS: an archive download supplies up to level 2 (cleaned events); "
                "reduction products (pcube, pha, arf, rmf, ...) are level 3.\n"
                "LOCKED ENTRIES: a line marked 'LOCKED, NOT USABLE' holds data that is "
                "present but GPG-encrypted, which happens while an observation is still "
                "inside its proprietary period. Its level is reported honestly for what the "
                "products are, so read the mark rather than the level: no analysis tool can "
                "open those bytes, and acquiring it again returns the same encrypted file "
                "until the observation becomes public. Treat such an ObsID as unavailable "
                "for analysis and say so, rather than re-acquiring it."
            ),
            "properties": {
                "obs_db_path": {
                    "type": "string",
                    "description": (
                        "Directory holding one observatory's data, e.g. "
                        "'/data/fzhou/XSPEC/Cyg_X-1/IXPE' (from your local observation "
                        "databases list). The catalog read is <obs_db_path>/catalog.json."
                    ),
                },
                "min_level": {
                    "type": "integer",
                    "description": (
                        "Only list ObsIDs at or above this processing level. A locked "
                        "(encrypted) ObsID is levelled honestly, so it can pass this "
                        "filter while still being unreadable — check each line's mark."
                    ),
                },
            },
            "required": ["obs_db_path"],
        },
    },
}
