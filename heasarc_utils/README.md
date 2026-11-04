# heasarc_utils

The **acquire** toolkit: fetch not-yet-local archive data for an ObsID into a mission's
obs-db directory and record it in the data catalog. A `*_utils` toolkit (CLAUDE.md §5) —
LLM-discoverable via `skills/TOOL_MANIFEST`, human-half in `lib/`.

## Purpose

Given `(obsid, mission, obs_db_root)`, download the observation from the HEASARC archive into
a local clone of the archive's own layout, `<obs_db_root>/<mission>/<HEASARC's tree>/<obsid>/`
(e.g. `ixpe/data/obs/01/01002901/`), and record its level-2 event in
`<obs_db_root>/<mission>/catalog.json`; the first acquire of a mission creates both. The mission
folder and the tree are read from the archive's own link, so they are HEASARC's names. This is
the first producer wired into [`obsdb_catalog`](../obsdb_catalog/README.md) — it records
through `obsdb_catalog.neutral.new`, the single trusted catalog writer.

## Architecture

The download path is resolved **authoritatively** via `astroquery.heasarc` — never by
hand-built FTP paths. This matters because, e.g., the NICER archive path component is the
*observation date* (`/FTP/nicer/data/obs/2022_03/<obsid>/`), which is not derivable from the
ObsID (CLAUDE.md §4, "ground every line — no guessing").

- **`lib/heasarc.py`** — the **mission-agnostic engine** (`resolve` → `archive_path` →
  `download` → `gunzip_tree`) plus the **`MissionProfile` registry** (`MISSIONS`). astroquery is
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
- **Fresh archive event files are gzipped** (`*.evt.gz` / `*.fits.gz`); **acquire gunzips the
  whole ObsID tree after download** (`heasarc.gunzip_tree`), so downstream reduction and the
  level-2 derivers operate on plain files and the recorded level-2 path is the decompressed
  `.fits`. The level-2 globs stay gz-tolerant and prefer the uncompressed twin for robustness.
- **The local tree is HEASARC's.** `archive_path` reads the mission folder and tree from the
  first `access_url` `locate_data` returns (`.../FTP/<mission>/<tree>/<obsid>/`; checked
  2026-11-04: `ixpe/data/obs/01/01002901/`, `nicer/data/obs/2022_02/4690010101/`,
  `nustar/data/obs/07/3/30702017002/`). Every astroquery host unpacks an observation as
  `<location>/<obsid>/` (the archive's tar is ObsID-rooted; the SciServer and S3 copies take the
  link's last folder), so `acquire_obsid` downloads with `location = <obs_db_root>/<mission>/<tree>`
  and nothing is moved. A link naming an archive folder other than the mission is refused.

## Interfaces / Usage

```python
from heasarc_utils.lib.heasarc_help import acquire

acquire("01002901", "ixpe", "/xray02/fzhou_obsdb")
# -> downloads to /xray02/fzhou_obsdb/ixpe/data/obs/01/01002901/, records level2_evt at level 2
#    in /xray02/fzhou_obsdb/ixpe/catalog.json, returns
#    {"obsid", "mission", "mission_dir", "obsid_dir", "level2_evt": <path relative to
#     mission_dir>, "catalog_entry", "downloaded": True, "locked": False}

# Large/slow IXPE transfers that would otherwise time out:
acquire("03005801", "ixpe", "/xray02/fzhou_obsdb", timeout_s=3600)
```

Signature: `acquire(obsid, mission, obs_db_root, *, host=None, chatter=False, timeout_s=1800.0)`.
A reader finds the ObsID's folder with `obsdb_catalog.locate(obs_db_root, mission, obsid)`.

- **Re-downloads on every call** — the tool executes; the caller decides *when* to acquire. A
  deliberate re-acquire records at level 2, regressing `max_level` to 2 and thereby flagging
  any stale level-3 products for re-reduction. The single exception is an encrypted
  observation, below.
- **Encrypted (proprietary) observations** — HEASARC serves an observation still inside its
  proprietary period GPG-encrypted (`…fits.gz.gpg`). That download is complete and correct;
  it just cannot be read without the key. So it is **recorded as present** (level 2 — the file
  genuinely is the level-2 event, and the `.gpg` path is what tells readers it is unusable)
  rather than raising, and the summary carries `locked: True`. Asking for the same ObsID again
  is then answered from the archive's `public_date` (`heasarc.availability`, one metadata row,
  no transfer): while it is still proprietary the download is **skipped** and the reply reports
  when it unlocks (`downloaded: False`); once public, acquire downloads normally and the
  readable copy replaces the encrypted record. Before this, the deriver rejected the `.gpg`
  name and raised *before* recording, so a finished multi-GB download left no trace and every
  later request re-downloaded it — a loop no number of runs could break.
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

| mission (HEASARC's folder) | master table | level-2 file (under the ObsID folder) | `params` |
|---|---|---|---|
| ixpe | `ixmaster` | `event_l2/ixpe<obsid>_det1_evt2_v0N.fits` | `detectors:[1,2,3]`, `calibration_version`, `note` |
| nicer | `nicermastr` | `xti/event_cl/ni<obsid>_0mpu7_cl.evt` | `mpu:7` |
| nustar | `numaster` | `event_cl/nu<obsid>A01_cl.evt` | `fpm:"A"`, `note` |

## Adding a mission

The engine is mission-agnostic: the local folder names and tree come from HEASARC's own link.
A mission HEASARC's `locate_data` serves needs:

1. Its master table: `heasarc.discover_master("<mission>")`
   (`list_catalogs(master=True, keywords=…)`).
2. A `MissionProfile` in `MISSIONS` in `lib/heasarc.py`, keyed by HEASARC's name for the mission
   (its `FTP/<name>/` folder): the master `catalog` name and a `level2(obsid_dir, obsid) -> product`
   deriver that globs the representative cleaned-event file (grounded in that mission's archive
   layout docs), returns its path and **raises if absent**.
3. In `obsdb_catalog/_domain.py`: its tree rule in `TREES` (the archive's path below
   `FTP/<name>/`, with `*` for the folders the ObsID alone does not give) and, for
   `repair`/`healthcheck` coverage, a scanner with the identical `level2_evt`
   `kind`/`params`/`produced_by`.

4. The mission in the manifest (`heasarc_utils/skills/__init__.py`: the `mission` enum and the
   description's mission list), and readers for its data where a toolkit analyses it.

A missing or wrong `TREES` rule leaves the mission's ObsID folders invisible to `repair` and
`healthcheck` (which then report nothing) and to `locate` without a catalog record.

Chandra is on HEASARC (`chanmaster`; `chandra/data/byobsid/<last digit>/<obsid>/`, its tree rule
already in `TREES`): it needs steps 2 and the scanner. MAXI is on HEASARC (`maximaster`; tree
`maxi/data/obs/MJD<thousands>/MJD<day>/`), but `locate_data` returns no link for its rows
(checked 2026-11-04), so it needs a locator of its own before this engine can fetch it.

`acquire` on an unregistered mission raises a clear, actionable error.

## Dependencies

- `astroquery` (≥ 0.4.x; pulls `pyvo`) + `astropy` — imported lazily, only when downloading.
- `obsdb_catalog.neutral` (same repo) — the catalog writer.
- Does **not** import `auto_research`.

## Tests

```bash
python3 -m pytest heasarc_utils/tests/
```
Offline: a fake `Heasarc` (injected via the `_heasarc` seam) creates the expected gzipped
files, so the resolve/normalize/gunzip/record/regression/extensibility paths run without
network or astroquery; `test_acquire_gunzips_download` asserts the `.gz` is decompressed and
removed and the plain file is recorded. The download layout, the `query_tap` resolve path, and
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
- 2026-11-04 — **The local obs-db clones HEASARC's archive.** `acquire(obsid, mission, obs_db_root)`
  (was `obs_db_path`) downloads into `<obs_db_root>/<mission>/<tree>/<obsid>/`, mission folder and
  tree read from the archive link (`heasarc.archive_path`), and records in
  `<obs_db_root>/<mission>/catalog.json` (created on the first acquire); catalog paths stay
  relative to the mission folder. `normalize_layout` is gone (nothing is moved); `gunzip_tree` and
  the level-2 derivers take the ObsID folder. The result names `mission_dir` and `obsid_dir`.
  Manifest: `obs_db_root` (kind `path`) replaces `obs_db_path`. Tests: an acquire into an empty
  root creates HEASARC's tree and a healthy catalog, per mission; a link naming another archive
  folder is refused.
