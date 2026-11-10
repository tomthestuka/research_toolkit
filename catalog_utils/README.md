# catalog_utils

The **data-catalog read**: what each mission's catalog says is already on disk. *(LLM toolkit)*

## Purpose

`data_catalog` lists the observations a mission's observation database holds, so an analyst
acquires only what is missing. One line per catalogued ObsID: its highest processing level, the
product kinds derived from it, the observed source as its event file's header names it, and the
mark **LOCKED, NOT USABLE** for data that is present but GPG-encrypted (an observation still in
its proprietary period).

## Architecture

- **`lib/catalog_read.py`** — the read and the listing, pure standard library. It reads
  `<obs_db_path>/catalog.json` through `obsdb_catalog.neutral.read` (the single trusted reader and
  writer, [`obsdb_catalog`](../obsdb_catalog/README.md)) and derives the lock mark live from the
  recorded product paths, since a human may decrypt data out of band.
- **`skills/__init__.py`** — `TOOL_MANIFEST` (pure data) wrapping `data_catalog`.

The caller passes the mission's `obs_db_path`; this package reads no configuration.

## Interfaces / Usage

- `data_catalog(obs_db_path, min_level=None) -> str` — the listing for the mission folder
  `obs_db_path`. `min_level` keeps the ObsIDs available at or above that processing level (a
  locked ObsID keeps its level, so the mark is what a caller reads). An absent catalog lists zero
  ObsIDs.
- Manifest key: `research_toolkit.catalog_utils.data_catalog`.

## Dependencies

Standard library and `obsdb_catalog`. The analysts call it as a toolkit tool; `auto_research`'s
analyst lead runs the same read as a one-shot call in the analysis jail.

## Changelog

- 2026-11-09 — README added (the package briefing the runtime roles read).
