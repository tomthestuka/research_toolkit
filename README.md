# research_toolkit -- an ai-controlled toolkit repo.

Each package's `README.md` is its briefing: read the one for the package you are working in.
Build-time conventions live in this repository's agent-instruction file, which runtime agents
do not read.

## Packages

- [`xspec_utils/`](xspec_utils/README.md) — joint spectro-polarimetric analysis tools
  (IXPE/NICER/NuSTAR) on HEASoft/XSPEC. *(LLM toolkit)*
- [`heasarc_utils/`](heasarc_utils/README.md) — acquire: download archive data for an ObsID
  and record it in the catalog. *(LLM toolkit)*
- [`obsdb_catalog/`](obsdb_catalog/README.md) — the single trusted writer for per-mission
  `catalog.json` (provenance/inventory). *(support package; not an LLM tool)*
