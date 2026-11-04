# research_toolkit -- an ai-controlled toolkit repo.

Orientation and conventions for coding agents: see [CLAUDE.md](CLAUDE.md). Each package's
`README.md` is its briefing.

## Packages

- [`xspec_utils/`](xspec_utils/README.md) — joint spectro-polarimetric analysis tools
  (IXPE/NICER/NuSTAR) on HEASoft/XSPEC. *(LLM toolkit)*
- [`heasarc_utils/`](heasarc_utils/README.md) — acquire: download archive data for an ObsID
  and record it in the catalog. *(LLM toolkit)*
- [`catalog_utils/`](catalog_utils/) — `data_catalog`: what each mission's catalog says
  is local, read-only. *(LLM toolkit)*
- [`obsdb_catalog/`](obsdb_catalog/README.md) — the single trusted writer for per-mission
  `catalog.json` (provenance/inventory). *(support package; not an LLM tool)*

## Argument kinds

Every `args_schema` property in a `TOOL_MANIFEST` declares a `kind`, which `auto_research`'s input
validator enforces before a call runs:

- `data`: a scientific artifact another tool returned, passed by reference to that result;
- `parameter`: a number that changes the science, passed by reference or typed with its reason;
  an optional one declares its `default`, which the runtime passes and records;
- `path`: a file or folder path;
- `setting`: text, a flag, an option, an identifier or an operational number (a timeout, an
  iteration count);
- `output`: where the tool writes what it produces. It names a `shape` (`folder`, `file` or
  `prefix`) and, for a file or prefix, may name a `default` file name; the runtime places it in
  the call's own folder.

A manifest with a property lacking a kind is refused at launch and at the authoring gate.

## What a tool writes

Every file a tool writes is one of three kinds:

- **a science artifact** (anything that depends on a `parameter` or on a setting outside the
  mission team's defaults): written under an `output` argument and returned under `products`
  (`{"path", "caption"}`) or, for a PNG, `figures`; the runtime receipts exactly what is returned;
- **an obs-db item** (raw archive data and the mission team's standard reductions with their
  recommended settings): written into the observation database and recorded through
  `obsdb_catalog.neutral`, never returned under `products`;
- **scratch** (a program's own bookkeeping): anywhere else, recorded and checked by nothing.

The authoring gate checks artifacts and obs-db items on what the tests call.

## Observation databases

An obs-db root holds one folder per mission, named as HEASARC names it (`ixpe`, `nicer`,
`nustar`), each with its `catalog.json`; below it every ObsID sits where HEASARC's archive keeps it
(`ixpe/data/obs/01/01002901/`). `heasarc_utils.acquire(obsid, mission, obs_db_root)` writes that
tree; the Python readers find an ObsID with `obsdb_catalog.locate(obs_db_root, mission, obsid)`
(listing ObsIDs through `obsdb_catalog._domain.obsid_dirs`). The `xspec_utils` batch shell scripts
still read a flat `<obsid>/` layout.

## Software environments

A package names the environment its tools run in, in its `tool_config.yaml` (`env: <id>`); the
definition lives at `envs/<id>/` (`env.yaml`, `env.lock`, optional `env.pip`, `activate.sh`,
`boot.py`), and `auto_research` builds it as `ar-<id>`. `envs/xspec/` serves `xspec_utils`.

## Tests

`conftest.py` defines the test classes: unmarked tests are fast and deterministic; `slow` ones run
on real data (`pytest -m "not slow"` skips them); `live` ones reach the network and run only with
`--live` or `TOOLKIT_LIVE=1`.

## Changelog

- 2026-11-04 — Argument kind `output` (with `shape` and `default`) and the three kinds of file a
  tool writes; `xspec_utils` migrated (its README).
- 2026-11-04 — Observation databases clone HEASARC's archive layout (`acquire(obs_db_root)`,
  `obsdb_catalog.locate`).
