# research_toolkit -- an ai-controlled toolkit repo.

Orientation and conventions for coding agents: see [CLAUDE.md](CLAUDE.md). Each package's
`README.md` is its briefing.

## Packages

- [`xspec_utils/`](xspec_utils/README.md) — joint spectro-polarimetric analysis tools
  (IXPE/NICER/NuSTAR) on HEASoft/XSPEC. *(LLM toolkit)*
- [`heasarc_utils/`](heasarc_utils/README.md) — acquire: download archive data for an ObsID
  and store it as a product in the data store. It reaches the HEASARC archive
  (`astroquery.heasarc`) for the missions its profile registry holds: IXPE, NICER and NuSTAR.
  *(LLM toolkit)*
- [`catalog_utils/`](catalog_utils/README.md) — `data_catalog`: the products the data store holds,
  each with its record, read-only. *(LLM toolkit)*
- [`storelib/`](storelib/README.md) — the tools' one door to the data store: read a record's version,
  stage records locally, store a product once per recipe. *(support package; not an LLM tool)*

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
  the call's own folder;
- `record`: an object in the data store, `{"store", "bucket", "key", "version_id"}` (or a list of
  them), as the data catalog lists it; the runtime pins a record given without its `version_id` to
  the current version and records the version read.

A tool entry that stores products also declares `writes` (the stores it writes, e.g. `["db"]`) and
`algorithm_version` (text bumped whenever what it stores changes). A manifest with a property lacking
a kind, or a writer without its version, is refused at launch and at the authoring gate.

## What a tool writes

Every file a tool writes is one of three kinds:

- **a science artifact** (anything that depends on a `parameter` or on a setting outside the
  mission team's defaults): written under an `output` argument and returned under `products`
  (`{"path", "caption"}`) or, for a PNG, `figures`; the runtime receipts exactly what is returned;
- **a store product** (raw archive data and the mission team's standard reductions with their
  recommended settings): stored in the data store through `storelib.produce`, once per recipe, by a
  tool whose manifest declares `writes`; its record is returned;
- **scratch** (a program's own bookkeeping, temporary files of a reduction): in the call's working
  folder, recorded and checked by nothing.

The authoring gate checks artifacts and store products on what the tests call.

## The data store

Observation data and the products made from it live in a data store (named by `TOOLKIT_STORE`, `db` in the default configuration): bytes behind
the S3 interface, and a catalog of what was made from what. A tool reaches it only through
`storelib`, with the keys its process is given. Every object is named by a record
`{"store", "bucket", "key", "version_id"}`; a product's files sit at
`<bucket>/<kind>/<subject>/<short hash>/`, beside a `recipe.json` listing them, and its record is the
one of that `recipe.json`. In this toolkit:

- **bucket** = the mission, as HEASARC names it: `ixpe`, `nicer`, `nustar`;
- **subject** = the ObsID, in the mission's native format (`01002901`, `4690010101`, `30702017002`);
- **kind** = what the product is: `raw` (an archive download, laid out as HEASARC's ObsID folder:
  `event_l2/...`, `xti/event_cl/...`, `event_cl/...`), `pha` (IXPE Stokes I/Q/U spectra);
- **short hash** = the first digits of the product's recipe hash (tool, `algorithm_version`,
  settings and the records it read): the same recipe is made once, and a new one has its own folder.

`catalog_utils.data_catalog` lists the products with their records; `heasarc_utils.acquire` stores a
download; `xspec_utils.load_data` takes the records of the products it loads.

## Software environments

A package names the environment its tools run in, in its `tool_config.yaml` (`env: <id>`); the
definition lives at `envs/<id>/` (`env.yaml`, `env.lock`, optional `env.pip`, `activate.sh`,
`boot.py`), and `auto_research` builds it as `ar-<id>`. `envs/xspec/` serves `xspec_utils`.
`envs/analysis/` is the toolkit's default environment (`auto_research`'s `default_env`): a derived
definition on the human `astropy` prefix with HEASoft and CALDB, where a unit's own REPL and shell
run. Its `env.yaml` `probe` is the launch smoke test. `envs/archive/` serves `heasarc_utils` and
`catalog_utils`: astroquery, astropy and the store's clients (boto3 with awscrt, s3fs, psycopg).

## Tests

`conftest.py` defines the test classes: unmarked tests are fast and deterministic; `slow` ones take
minutes (`pytest -m "not slow"` skips them); `live` ones reach the network and run only with
`--live` or `TOOLKIT_LIVE=1`; `data("<name>")` ones read a named data set from the folder
`TOOLKIT_DATA_<NAME>` names (e.g. `TOOLKIT_DATA_CYGX1_IXPE=/path/to/IXPE pytest`) and skip, with the
reason, without it. A test taking the `store` fixture reads and writes a throwaway store the host
gives (`TOOLKIT_STORE` and the store's keys) and skips, with the reason, without one; its inputs reach
the store only through the host's seed channel (the `seed` fixture, `TOOLKIT_STORE_SEED`: a synthetic
product of given files, or a data set's file), and only a manifest tool's call writes the store. A tool
reading the store without a `record` argument declares `reads: ["db"]` (`data_catalog`).

## Changelog

- 2026-11-04 — Argument kind `output` (with `shape` and `default`) and the three kinds of file a
  tool writes; `xspec_utils` migrated (its README).
- 2026-11-04 — Observation databases clone HEASARC's archive layout (`acquire` into the obs-db root,
  `obsdb_catalog.locate`).
- 2026-11-09 — `envs/analysis/`: the default environment, moved here from `auto_research`'s config
  (`analyze_env`, `analyst_env.sh`, `analyze_probe`); `env.yaml` takes a one-line `probe`.
- 2026-11-09 — Packages: `heasarc_utils` states the archive it reaches (HEASARC) and its missions
  (IXPE, NICER, NuSTAR); the roles read this file as the toolkit's domain briefing.
- 2026-11-09 — `catalog_utils/README.md` added; every package carries its briefing (checked at launch).
- 2026-11-10 — The data store (auto_research Phase 3 step 0b, S8): argument kind `record`, manifest
  `writes` + `algorithm_version`; store products replace obs-db items; `storelib` replaces
  `obsdb_catalog` (deleted; `catalog.json` is gone); the store briefing (bucket = mission, subject =
  ObsID, kind = `raw`/`pha`); `envs/archive/` for `heasarc_utils` and `catalog_utils`; the `store` test
  fixture.
- 2026-11-10 — The store's name is read from `TOOLKIT_STORE` (the runtime and the CI pipeline set it); a
  key already holding other bytes fails a production as an ordinary failure (`storelib`).
- 2026-11-10 — Tests stage store inputs through the seed channel (`seed` fixture); `data_catalog`
  declares `reads: ["db"]`; only a manifest tool writes the store (`red_ixpe.pha` writes locally).
