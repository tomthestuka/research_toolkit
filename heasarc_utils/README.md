# heasarc_utils

The **acquire** toolkit: fetch archive data for an ObsID and store it in the data store as a product
of kind `raw`. A `*_utils` toolkit (CLAUDE.md §5) — LLM-discoverable via `skills/TOOL_MANIFEST`,
human-half in `lib/`.

## Purpose

Given `(obsid, mission)`, download the observation from the HEASARC archive and store it, through
[`storelib`](../storelib/README.md), as the product `raw/<obsid>/<short hash>/` in the mission's bucket
of the store the process's keys reach (`storelib.store_name()`, `TOOLKIT_STORE`): every file of the ObsID folder at its path in that folder (`event_l2/...`,
`xti/event_cl/...`), beside the product's `recipe.json`. The result names the product's record, which
the analysis tools take, and the record of its level-2 event.

## Architecture

The download path is resolved **authoritatively** via `astroquery.heasarc` — never by hand-built FTP
paths. This matters because, e.g., the NICER archive path component is the *observation date*
(`/FTP/nicer/data/obs/2022_03/<obsid>/`), which is not derivable from the ObsID (CLAUDE.md §4, "ground
every line — no guessing").

- **`lib/heasarc.py`** — the **mission-agnostic engine** (`resolve` → `archive_path` → `archive_state`
  → `download_obsid` = `download` + `gunzip_tree`) plus the **`MissionProfile` registry** (`MISSIONS`).
  astroquery is imported lazily (`_heasarc()`) so the module loads without it and tests inject a fake.
- **`lib/heasarc_help.py`** — `acquire(...)`, the abstracted human/LLM entry; wires the engine to
  `storelib.produce`.
- **`skills/__init__.py`** — `TOOL_MANIFEST` (pure data; no astroquery import) wrapping `acquire`,
  declaring `writes: ["db"]` and `algorithm_version`.
- **`tool_config.yaml`** — the call timeout ceiling and `env: archive` (`envs/archive/`: astroquery,
  astropy, boto3 + awscrt, s3fs, psycopg).

Grounded in — and **live-verified against** — the installed `astroquery.heasarc` (v0.4.x):
`query_tap("SELECT * FROM <master> WHERE obsid='…'") → locate_data(tab, catalog_name=…) →
download_data(location=…)` ([docs](https://astroquery.readthedocs.io/en/stable/heasarc/heasarc.html)).
**astroquery only — no wget fallback; failures raise and store nothing.**

Facts the live runs pinned down (don't "fix" them back):
- **Resolve uses `query_tap`, not `query_region(column_filters=…)`** — the latter silently drops the
  obsid filter for a position-less query and returns the whole catalog. `SELECT *` carries the `__row`
  column `locate_data` needs, and the row's archive state (below).
- **Fresh archive event files are gzipped** (`*.evt.gz` / `*.fits.gz`); **acquire gunzips the whole
  ObsID tree after download** (`heasarc.gunzip_tree`), so downstream reduction operates on plain files
  and the stored level-2 event is the decompressed file. The level-2 matchers stay gz-tolerant.
- **The archive's link names the mission.** `archive_path` reads the mission folder and tree from the
  first `access_url` `locate_data` returns (`.../FTP/<mission>/<tree>/<obsid>/`; checked 2026-11-04:
  `ixpe/data/obs/01/01002901/`, `nicer/data/obs/2022_02/4690010101/`,
  `nustar/data/obs/07/3/30702017002/`); a link naming an archive folder other than the mission is
  refused. Every astroquery host unpacks an observation as `<location>/<obsid>/`; `download_obsid`
  moves that folder's contents up into the production's temporary folder, so the stored paths are the
  ObsID folder's own.

## Interfaces / Usage

```python
from heasarc_utils.lib.heasarc_help import acquire

acquire("01002901", "ixpe")
# -> {"obsid", "mission", "record": {"store": "db", "bucket": "ixpe",
#     "key": "raw/01002901/<short hash>/recipe.json", "version_id": ...},
#     "level2_evt": {... "key": "raw/01002901/<short hash>/event_l2/ixpe01002901_det1_evt2_v01.fits"},
#     "downloaded": True, "locked": False, "archive": {...}, "files": <n>}

acquire("03005801", "ixpe", timeout_s=3600)   # large/slow IXPE transfers
```

Signature: `acquire(obsid, mission, *, host=None, chatter=False, timeout_s=1800.0)`. The process holds
the unit's store keys (the runtime sets them); `TOOL`, `ALGORITHM_VERSION` (`"1"`), `KIND` (`"raw"`)
are the recipe's names.

- **One download per archive state.** The recipe's settings are the mission and the master-table row's
  `processing_date`, `software_version`, `caldb_version` (the archive replaces files in place when it
  reprocesses an observation, and these columns change then; all three registered master tables carry
  them, HEASARC TAP schema checked 2026-11-10) and whether the observation is public yet
  (`public_date`). A repeat acquire of the same state is answered from the store
  (`downloaded: False`, no transfer); a reprocessed observation, or one that has become public, is a
  new recipe and a fresh download. `host` and `timeout_s` change no bytes and are not in the recipe.
- **Encrypted (proprietary) observations** — HEASARC serves an observation still inside its
  proprietary period GPG-encrypted (`…fits.gz.gpg`). That download is complete and correct; it is
  stored, and the summary carries `locked: True`, `days_until_public` and a `note`. Asking again while
  it is still proprietary is the same recipe, so the stored copy answers without a transfer; once it
  is public the recipe changes and the readable copy is downloaded.
- **Fails loudly** — unknown mission, ObsID not in the master table, a row stating no
  `processing_date` (the archive state identifies a download), an archive folder other than the
  mission's, a download that does not yield the expected level-2 file, a download whose level-2
  header states another processing than the row (its `SOFTVER` differs from `software_version`, or its
  `CALDBVER` from a named `caldb_version`: the archive lags its catalog, seen for IXPE 01250101, row 1.11
  and files 1.5; `heasarc.state_mismatch`), or a public observation the
  archive still serves encrypted (a lag; storing it would answer every later acquire of that state
  with the locked copy) all raise
  `heasarc.HeasarcError`; the claim is given up with the reason (`catalog.fail`) and nothing is
  committed. `storelib.ProductBusy` while another unit downloads the same recipe.
- `host` (`heasarc`/`sciserver`/`aws`) is optional; omit to let astroquery choose.
- `timeout_s` — HTTP transfer timeout in seconds. Defaults to 1800 s (30 min). The astroquery factory
  default is **30 s**, far too short for large IXPE observations (> 1 GB); this tool always overrides
  it: `acquire(timeout_s)` → `download_obsid(timeout_s)` → `download(timeout_s)` → `h.timeout =
  timeout_s` on the astroquery singleton before `download_data`. The claim's lease (7 h) outlives the
  tool's call timeout ceiling (6 h).

Shipped missions and the level-2 representative each checks (a path below the ObsID folder):

| mission (HEASARC's folder = bucket) | master table | level-2 file |
|---|---|---|
| ixpe | `ixmaster` | `event_l2/ixpe<obsid>_det1_evt2_v0N.fits` |
| nicer | `nicermastr` | `xti/event_cl/*mpu7_cl.evt` |
| nustar | `numaster` | `event_cl/nu<obsid>A01_cl.evt` |

## Adding a mission

The engine is mission-agnostic: the mission folder comes from HEASARC's own link. A mission HEASARC's
`locate_data` serves needs:

1. Its master table: `heasarc.discover_master("<mission>")` (`list_catalogs(master=True, keywords=…)`),
   with `processing_date`, `software_version`, `caldb_version` and `public_date` columns (check its TAP
   schema; a column it lacks reads as `None`).
2. A `MissionProfile` in `MISSIONS` in `lib/heasarc.py`, keyed by HEASARC's name for the mission (its
   `FTP/<name>/` folder, and the store bucket): the master `catalog` name and a
   `level2(paths, obsid) -> path` deriver that picks the representative cleaned-event file among the
   ObsID folder's relative paths (grounded in that mission's archive layout docs) and **raises if
   absent**.
3. The mission in the manifest (`heasarc_utils/skills/__init__.py`: the `mission` enum and the
   description's mission list), a bucket of that name in the store, and readers for its data where a
   toolkit analyses it.

Chandra is on HEASARC (`chanmaster`; `chandra/data/byobsid/<last digit>/<obsid>/`): it needs step 2.
MAXI is on HEASARC (`maximaster`), but `locate_data` returns no link for its rows (checked
2026-11-04), so it needs a locator of its own before this engine can fetch it.

`acquire` on an unregistered mission raises a clear, actionable error.

## Dependencies

- `astroquery` (≥ 0.4.x; pulls `pyvo`) + `astropy` — imported lazily, only when downloading.
- `storelib` (same repo) — the store; `envs/archive/` holds its clients.
- Used by `catalog_utils` (the mission profiles find a raw product's level-2 event).
- Does **not** import `auto_research`.

## Tests

```bash
python3 -m pytest heasarc_utils/tests/
```
Offline: a fake `Heasarc` (injected via the `_heasarc` seam) answers the master-table row and
creates the expected gzipped files, so the resolve/gunzip/store/reuse/reprocessing/encryption paths run
without network or astroquery, against a throwaway store (the root conftest's `store` fixture; skipped
with its reason without one). The download layout, the `query_tap` resolve path, and
gzip handling were pinned by a live NICER `acquire` run on 2026-06-17 (see
`logs/dev/live_acquire_nicer.md`).

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
- 2026-06-25 — `acquire` now **gunzips the downloaded ObsID tree** before deriving the
  level-2 product (`heasarc.gunzip_tree`, called in `acquire_obsid` after `normalize_layout`).
  Fresh archive downloads are gzipped; decompressing at acquire time means reduction and the
  level-2 derivers see plain `.fits`/`.evt` files and the recorded level-2 path is the
  decompressed file. Decompression writes a `.part` sidecar then renames, so an interrupted
  gunzip leaves the `.gz` intact. Tests: `test_acquire_gunzips_download` (3 missions); the
  offline fakes now write valid gzip content and the level-2 path assertion expects the plain
  file. The level-2 globs remain gz-tolerant.
- 2026-08-03 — **An encrypted download is data, not an error.** HEASARC serves a still-
  proprietary observation GPG-encrypted; `_ixpe_level2` matched the name through
  `_IXPE_BASE_EVT`, which the `.gpg` suffix failed, so it **raised before `acquire` ever
  reached `neutral.new`** — a completed multi-GB download recorded nothing, the catalog could
  never learn the ObsID existed, and each later request re-downloaded it (observed live: SS 433
  `04004001`, 46 encrypted files, re-fetched for ~25 min to no effect). Now the level-2
  matchers run against the name with its wrapper split off (`neutral.split_lock`, the one
  shared rule), `_pick` prefers a readable copy over an encrypted twin, and the encrypted
  product is recorded at its honest level with `locked: True`. New `heasarc.availability(catalog,
  obsid)` reads `public_date` from the mission master table (present in all three), so a repeat
  request for a locked ObsID is answered from metadata instead of a re-transfer: skip while
  proprietary, download once public. `acquire`'s summary gains `downloaded`/`locked` (+
  `public_date_mjd`/`days_until_public`/`note` on a skip); `TOOL_MANIFEST` states the behavior.
  Five offline tests (`test_acquire.py`), no network.
- 2026-11-04 — **The local obs-db clones HEASARC's archive.** `acquire(obsid, mission, <obs-db root>)`
  (was the mission folder) downloads into `<obs-db root>/<mission>/<tree>/<obsid>/`, mission folder and
  tree read from the archive link (`heasarc.archive_path`), and records in
  `<obs-db root>/<mission>/catalog.json` (created on the first acquire); catalog paths stay
  relative to the mission folder. `normalize_layout` is gone (nothing is moved); `gunzip_tree` and
  the level-2 derivers take the ObsID folder. The result names `mission_dir` and `obsid_dir`.
  Manifest: the obs-db root (kind `path`) replaces the mission folder. Tests: an acquire into an empty
  root creates HEASARC's tree and a healthy catalog, per mission; a link naming another archive
  folder is refused.
- 2026-11-04 — A flat mission folder (`<mission>/<obsid>/`) is refused as an acquire target, the
  refusal naming both layouts (ruling: the obs-db is HEASARC's tree only).
- 2026-11-10 — **Downloads are store products** (auto_research Phase 3 step 0b, S8).
  `acquire(obsid, mission)` stores the observation through `storelib.produce` (bucket = mission,
  kind `raw`, subject = ObsID; files at their ObsID-folder paths) and returns its record and its
  level-2 event's record; the recipe holds the archive state (processing date, software and CALDB
  versions, public or not), so a repeat is answered from the store and a reprocessing or a lapsed
  proprietary period downloads afresh (`availability` and the locked-copy check retired). Engine:
  `resolve` returns the row with the links, `archive_state`, `download_obsid`, `profile`, `mjd_now`;
  level-2 derivers pick among relative paths; `split_lock`/`is_locked`/`LOCK_SUFFIX` moved here from
  `obsdb_catalog`. Manifest: no store location argument; `writes: ["db"]`, `algorithm_version: "1"`.
  `tool_config.yaml`: `env: archive`. Tests against a throwaway store.
- 2026-11-10 — Stateless review fixes: a row with no `processing_date` is refused before any claim; a
  public observation still served encrypted fails its production (nothing stored; a later acquire
  retries); manifest text says when to acquire again and what a busy answer names.
- 2026-11-10 — `acquire` stores into the store `TOOLKIT_STORE` names (`storelib.store_name()`).
- 2026-11-11 — Before anything is stored, the level-2 event's own `SOFTVER` (and `CALDBVER` when the
  row names one) must equal the row's `software_version`/`caldb_version` (`heasarc.file_state`,
  `state_mismatch`); a lagging archive fails the production naming both values. The keyword match was
  checked on IXPE 01250101, NuSTAR 30702017002 and NICER 4690010101 (its archive header, 2026-11-11); the processing date has no header keyword and is
  not compared. Fixture events are FITS headers stating the row's versions.
