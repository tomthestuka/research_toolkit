"""storelib — the toolkit's door to a data store (S3 bytes + a Postgres catalog).

A plain support package (CLAUDE.md §2): no ``skills/`` and no ``TOOL_MANIFEST``, and its name lies
outside the ``*_utils`` glob, so tool discovery never lists it. Tools import it to read records and to
store products; see ``README.md``. ``boto3`` and ``psycopg`` are imported only when a function runs.
"""

from __future__ import annotations

from .store import (STORE_NAME_VARIABLE, ProductBusy, ProductRefused, StoreError, check_record, client, connect, is_product,
                    normalise, open_remote, produce, product_folder, remote, resolve, stage_in, store_name)

__all__ = ["STORE_NAME_VARIABLE", "ProductBusy", "ProductRefused", "StoreError", "check_record", "client", "connect", "is_product",
           "normalise", "open_remote", "produce", "product_folder", "remote", "resolve", "stage_in", "store_name"]
