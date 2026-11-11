# catalog_utils

The **data-catalog read**: what the data store holds, each product with the record a tool takes.
*(LLM toolkit)*

## Purpose

`data_catalog` lists the products the store's catalog holds, so an analyst acquires only what is
missing and passes the right records to the analysis tools. One line per product: its ObsID (the
subject), the observed source as a raw download's level-2 event header names it, its kind, its state
(`current`, `verified` or `committed`), its file count and its record; the mark **LOCKED, NOT
USABLE** for a download whose level-2 event is GPG-encrypted (an observation still in its proprietary
period).

## Architecture

- **`lib/catalog_read.py`** — one catalog query (the unit's own Postgres login, libpq `PG*`) and the
  listing. Per (bucket, kind, subject), the catalog's key of the current pointer, it lists the current
  product and every committed or verified product claimed after it, newest first, grouped by bucket
  (the mission). A raw download's level-2 event is
  found among its files by `heasarc_utils`' mission profiles; its source name is read from the first
  two FITS headers by one version-pinned ranged GET (`storelib.client`), and only a name of printable
  characters outside the listing's own punctuation is shown, so header bytes cannot forge a line.
- **`skills/__init__.py`** — `TOOL_MANIFEST` (pure data) wrapping `data_catalog`.
- **`tool_config.yaml`** — `env: archive` (`envs/archive/`: psycopg, boto3).

## Interfaces / Usage

- `data_catalog(mission=None, kind=None, subject=None) -> str` — the listing; each argument keeps the
  products of that bucket, kind or subject. An empty store or filter lists `(no product stored)`.
- Manifest key: `research_toolkit.catalog_utils.data_catalog`.

## Dependencies

`storelib` (the catalog connection and the S3 client), `heasarc_utils.lib.heasarc` (mission profiles,
`is_locked`). The analysts call it as a toolkit tool, in `envs/archive/`.

## Changelog

- 2026-11-09 — README added (the package briefing the runtime roles read).
- 2026-11-10 — **The catalog is the store's** (auto_research Phase 3 step 0b, S8): `data_catalog(mission,
  kind, subject)` queries the Postgres catalog (products, the current pointer) instead of a mission's
  `catalog.json` (retired with `obsdb_catalog`); each line carries its record; levels are replaced by
  kinds and states; the source name is read from the stored event by a pinned ranged GET. Tests
  against a throwaway store. `tool_config.yaml` added (`env: archive`).
- 2026-11-10 — Stateless review fixes: a product's files and versions come from one aggregate of pairs; a
  product is matched by its own `<kind>/<subject>/<hash>/recipe.json`.
- 2026-11-10 — The current pointer is keyed (bucket, kind, subject) in the catalog: the listing joins it on
  the bucket too, so two missions sharing an ObsID string each show their own current product once.
  Test: auto_research `tests/test_toolkit_store.py` (pointers need the runtime's reconciler, which a
  toolkit test's producer keys cannot act as).
- 2026-11-10 — The manifest declares `reads: ["db"]` (it reads the store without a record argument);
  its tests seed their products through the host's seed channel.
