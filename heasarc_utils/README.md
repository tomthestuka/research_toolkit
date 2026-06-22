# heasarc_utils

The **acquire** toolkit: fetch not-yet-local archive data for an ObsID into a mission's
obs-db directory and record it in the data catalog. A `*_utils` toolkit (CLAUDE.md §5) —
LLM-discoverable via `skills/TOOL_MANIFEST`, human-half in `lib/`.

## Purpose

Given `(obsid, mission, obs_db_path)`, download the observation from the HEASARC archive into
`<obs_db_path>/<obsid>/` and record its level-2 event in `<obs_db_path>/catalog.json`. This is
the first producer wired into [`obsdb_catalog`](../obsdb_catalog/README.md) — it records
through `obsdb_catalog.neutral.new`, the single trusted catalog writer.

## Architecture

The download path is resolved **authoritatively** via `astroquery.heasarc` — never by
hand-built FTP paths. This matters because, e.g., the NICER archive path component is the
*observation date* (`/FTP/nicer/data/obs/2022_03/<obsid>/`), which is not derivable from the
ObsID (CLAUDE.md §4, "ground every line — no guessing").

- **`lib/heasarc.py`** — the **mission-agnostic engine** (`resolve` → `download` →
  `normalize_layout`) plus the **`MissionProfile` registry** (`MISSIONS`). astroquery is
  imported lazily (`_heasarc()`) so the module loads without it and tests inject a fake.
- **`lib/heasarc_help.py`** — `acquire(...)`, the abstracted human/LLM entry; wires the engine
  to `obsdb_catalog.neutral.new`.
- **`skills/__init__.py`** — `TOOL_MANIFEST` (pure data; no astroquery import) wrapping
  `acquire`, so the orchestration environment can discover the tool.

Grounded in — and **live-verified against** — the installed `astroquery.heasarc` (v0.4.x):
`query_tap("SELECT * FROM <master> WHERE obsid='…'") → locate_data(tab, catalog_name=…) →
download_data(location=…)`
([docs](https://astroquery.readthedocs.io/en/stable/heasarc/heasarc.html)). **astroquery
only — no wget fallback; failures raise and record nothing.**

Two facts the live run pinned down (don't "fix" them back):
- **Resolve uses `query_tap`, not `query_region(column_filters=…)`** — the latter silently
  drops the obsid filter for a position-less query and returns the whole catalog. `SELECT *`
  carries the `__row` column `locate_data` needs.
- **Fresh archive event files are gzipped** (`*.evt.gz` / `*.fits.gz`); reduction later
  gunzips them. The level-2 globs are gz-tolerant and prefer the uncompressed twin when both
  exist. `download_data` untars an obsid-rooted tar, so data lands cleanly at
  `<obs_db_path>/<obsid>/`.

## Interfaces / Usage

```python
from heasarc_utils.lib.heasarc_help import acquire

acquire("01002901", "ixpe", "/data/fzhou/XSPEC/Cyg_X-1/IXPE")
# -> downloads to <obs_db_path>/01002901/, records level2_evt at level 2,
#    returns {"obsid", "mission", "level2_evt": <rel path>, "catalog_entry"}

# Large/slow IXPE transfers that would otherwise time out:
acquire("03005801", "ixpe", "/data/fzhou/XSPEC/Cyg_X-1/IXPE", timeout_s=3600)
```

Signature: `acquire(obsid, mission, obs_db_path, *, host=None, chatter=False, timeout_s=1800.0)`

- **Always (re)downloads** — the tool executes; the caller decides *when* to acquire. A
  deliberate re-acquire records at level 2, regressing `max_level` to 2 and thereby flagging
  any stale level-3 products for re-reduction.
- **Fails loudly** — unknown mission, ObsID not in the master table, or a download that does
  not yield the expected level-2 file all raise `heasarc.HeasarcError`; nothing is recorded.
- `host` (`heasarc`/`sciserver`/`aws`) is optional; omit to let astroquery choose.
- `timeout_s` — HTTP transfer timeout in seconds. Defaults to 1800 s (30 min). The astroquery
  factory default is **30 s**, which is far too short for large IXPE observations (> 1 GB).
  This tool always overrides it. The three-layer threading path is:
  `acquire(timeout_s)` → `acquire_obsid(timeout_s)` → `download(timeout_s)` →
  `h.timeout = timeout_s` on the astroquery singleton before `download_data`.

Shipped missions and the level-2 representative each records (these mirror
`obsdb_catalog._domain` exactly — keep the two in sync):

| mission | master table | level-2 file (under `<obsid>/`) | `params` |
|---|---|---|---|
| ixpe | `ixmaster` | `event_l2/ixpe<obsid>_det1_evt2_v0N.fits` | `detectors:[1,2,3]`, `calibration_version`, `note` |
| nicer | `nicermastr` | `xti/event_cl/ni<obsid>_0mpu7_cl.evt` | `mpu:7` |
| nustar | `numaster` | `event_cl/nu<obsid>A01_cl.evt` | `fpm:"A"`, `note` |

## Adding a mission

The engine is mission-agnostic; adding a mission (Chandra, XMM, Swift, Suzaku, …) is **one
registry entry**:

1. Find its master table: `heasarc.discover_master("<mission>")`
   (`list_catalogs(master=True, keywords=…)`).
2. Add a `MissionProfile` to `MISSIONS` in `lib/heasarc.py`: the master `catalog` name and a
   `level2(obs_db_path, obsid) -> product` deriver that globs the representative cleaned-event
   file (grounded in that mission's archive layout docs) and **raises if absent**.
3. For `repair`/`healthcheck` coverage, add the matching scanner to
   `obsdb_catalog/_domain.py` (the catalog-side analogue) with the identical `level2_evt`
   `kind`/`params`/`produced_by`.

`acquire` on an unregistered mission raises a clear, actionable error.

## Dependencies

- `astroquery` (≥ 0.4.x; pulls `pyvo`) + `astropy` — imported lazily, only when downloading.
- `obsdb_catalog.neutral` (same repo) — the catalog writer.
- Does **not** import `auto_research`.

## Tests

```bash
python3 -m pytest heasarc_utils/tests/
```
Offline: a fake `Heasarc` (injected via the `_heasarc` seam) creates the expected (gzipped)
files, so the resolve/normalize/record/regression/extensibility paths run without network or
astroquery. The download layout, the `query_tap` resolve path, and gzip handling were pinned
by a live NICER `acquire` run on 2026-06-17 (see `logs/dev/live_acquire_nicer.md`).

## Changelog

- 2026-06-17 — Initial build: mission-agnostic astroquery acquire engine + `MissionProfile`
  registry (ixpe/nicer/nustar), `acquire` recording level-2 events via `obsdb_catalog.neutral`,
  pure-data `TOOL_MANIFEST`, offline tests.
- 2026-06-17 — Live-verified against installed astroquery; two corrections: resolve via
  `query_tap` (not `query_region(column_filters=…)`, which drops the obsid filter), and
  gz-tolerant level-2 globs (fresh downloads are `.gz`). Layout confirmed obsid-rooted.
- 2026-06-22 — `timeout_s` parameter added (default 1800 s). Threaded from
  `heasarc_help.acquire` → `heasarc.acquire_obsid` → `heasarc.download`, where it is
  applied as `h.timeout = timeout_s` on the astroquery singleton before `download_data`.
  Fixes IXPE acquisition timeouts on large (> 1 GB) observations (03005801, 02001901,
  01250401) that aborted under the astroquery factory default of 30 s. Surfaced in
  `TOOL_MANIFEST` args_schema. Three new tests: differential wiring, error-contract
  regression, default-in-all-layers guard.
