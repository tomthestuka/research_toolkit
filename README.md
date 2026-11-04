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
  iteration count).

A manifest with a property lacking a kind is refused at launch and at the authoring gate.

## Software environments

A package names the environment its tools run in, in its `tool_config.yaml` (`env: <id>`); the
definition lives at `envs/<id>/` (`env.yaml`, `env.lock`, optional `env.pip`, `activate.sh`,
`boot.py`), and `auto_research` builds it as `ar-<id>`. `envs/xspec/` serves `xspec_utils`.

## Tests

`conftest.py` defines the test classes: unmarked tests are fast and deterministic; `slow` ones run
on real data (`pytest -m "not slow"` skips them); `live` ones reach the network and run only with
`--live` or `TOOLKIT_LIVE=1`.
