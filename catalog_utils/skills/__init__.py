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
            "List locally-available observation data for a mission — which ObsIDs are "
            "already on disk, their processing level (HEASARC supplies up to level 2), and "
            "the data products produced from them — read from the mission data catalog at "
            "<obs_db_path>/catalog.json. Use this BEFORE requesting data acquisition so you "
            "only acquire what is genuinely missing. Pass the mission's obs_db_path (listed "
            "in your local observation databases)."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "obs_db_path": {
                    "type": "string",
                    "description": (
                        "The mission's obs-db directory, e.g. "
                        "'/data/fzhou/XSPEC/Cyg_X-1/IXPE' (from your local observation "
                        "databases list). The catalog read is <obs_db_path>/catalog.json."
                    ),
                },
                "min_level": {
                    "type": "integer",
                    "description": "Only list ObsIDs available at or above this processing level.",
                },
            },
            "required": ["obs_db_path"],
        },
    },
}
