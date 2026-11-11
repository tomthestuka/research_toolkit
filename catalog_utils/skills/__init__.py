"""catalog_utils.skills — AI-facing tool interface for the data catalog read.

Exposes TOOL_MANIFEST for ToolkitRegistry discovery. Pure-data module: no heavy imports here (the read
logic in ``lib/catalog_read.py`` imports ``storelib`` only when the tool runs).
"""

from __future__ import annotations

TOOL_MANIFEST = {
    "research_toolkit.catalog_utils.data_catalog": {
        "module": "catalog_utils.lib.catalog_read",
        "fn": "data_catalog",
        "reads": ["db"],
        "description": (
            "List the observation data and products the data store holds: one line per product with "
            "its ObsID, the source its event file names, its kind, its state and its record. Read it "
            "BEFORE acquiring, so you acquire only what is missing, and to find the records the "
            "analysis tools take."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "PRODUCTS: each line is one product: an ObsID's archive download (kind raw) or a product "
                "a tool made from it. Its 'record' (store, bucket, key, version_id) is the value a tool "
                "takes for a record argument; pass it as listed.\n"
                "STATE: 'current' marks the product the catalog names for its mission, kind and ObsID; "
                "'verified' means the runtime has checked its bytes; 'committed' means it was stored and "
                "awaits that check. The listing is grouped by mission; within a mission, for each kind and "
                "ObsID, the current product is listed with any made after it, newest first.\n"
                "LOCKED ENTRIES: a line marked 'LOCKED, NOT USABLE' holds data that is stored but "
                "GPG-encrypted, which happens while an observation is still inside its proprietary period. "
                "No analysis tool can open those bytes; acquiring it again returns the same copy until the "
                "observation becomes public. Treat such an ObsID as unavailable for analysis and say so.\n"
                "FILTERS: mission, kind and subject each keep the products that match; give none to list "
                "everything."
            ),
            "properties": {
                "mission": {
                    "kind": "setting",
                    "type": "string",
                    "description": "Keep the products of this mission's bucket, e.g. 'ixpe', 'nicer', 'nustar'.",
                },
                "kind": {
                    "kind": "setting",
                    "type": "string",
                    "description": "Keep the products of this kind, e.g. 'raw' for archive downloads.",
                },
                "subject": {
                    "kind": "setting",
                    "type": "string",
                    "description": "Keep the products of this ObsID.",
                },
            },
            "required": [],
        },
    },
}
