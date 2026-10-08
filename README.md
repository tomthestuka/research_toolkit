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

## Software environments

A package names the environment its tools run in, in its `tool_config.yaml` (`env: <id>`); the
definition lives at `envs/<id>/` (`env.yaml`, `env.lock`, optional `env.pip`, `activate.sh`,
`boot.py`), and `auto_research` builds it as `ar-<id>`. `envs/xspec/` serves `xspec_utils`.

## Tests

`conftest.py` defines the test classes: unmarked tests are fast and deterministic; `slow` ones run
on real data (`pytest -m "not slow"` skips them); `live` ones reach the network and run only with
`--live` or `TOOLKIT_LIVE=1`.
