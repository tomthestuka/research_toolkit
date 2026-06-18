# obsdb_catalog

The **single trusted writer** for a mission's data catalog
(`<obs_db_path>/catalog.json` — the inventory of what observation data is local
and how it was produced). Standalone, standard-library only: no science stack,
no `auto_research`, no sibling-toolkit import. The JSON shape is the entire
contract; a separate read side (`auto_research`) consumes it read-only.

**Prisoned write path.** `catalog.json` is mutated *only* through this package —
`neutral` for producers, `repair` for deliberate recovery. No producer, script,
or runtime-authored tool ever opens the catalog for writing. Not `*_utils`-named
and exposes no `skills/` / `TOOL_MANIFEST` / CLI, so it is invisible to LLM tool
discovery; catalog writes are a side effect of data-producing libraries, never a
model action.

## Layers

| Module | Role | Domain knowledge |
|---|---|---|
| `neutral` | The writer every producer calls. Durable atomic write, versioned history, `flock`-serialised, idempotent, truthful. | none (pure stdlib) |
| `repair` | Deliberate recovery: restore newest parseable snapshot, reconcile against disk, additively re-derive. Never fabricates. | yes (`_domain`) |
| `healthcheck` | Detect-only drift report; never writes. | yes (`_domain`) |
| `_domain` | Per-mission layout/filename conventions, `kind`→level, disk scanners. Imported only by `repair`/`healthcheck`. | — |

## Producer usage

```python
from obsdb_catalog import neutral

# fresh download -> level-2 events
neutral.new(obs_db_path, "ixpe", obsid, products=[
    {"kind": "level2_evt", "path": "<rel>", "produced_by": "heasarc",
     "params": {"detectors": [1, 2, 3], "calibration_version": "v01"}},
])

# reduction output at a target level (idempotent upsert on (kind, path))
neutral.advance(obs_db_path, "nicer", obsid, level=3, products=[
    {"kind": "pha", "path": "<rel>", "produced_by": "nicer_red.sh",
     "params": {"background_model": "3C50", "mpu": 7}},
])
```

`level` is authoritative and may regress (e.g. 3→2 after an upstream pipeline
re-run): stale higher-level products are left in place (overwritten in place by
the next re-reduction), and the dropped `max_level` is what flags them stale to
readers. `levels_present = [1..max_level]`. A product is recorded only after its
file is confirmed on disk; product `path` is relative to `obs_db_path`.

## Maintenance

```python
from obsdb_catalog.healthcheck import healthcheck
from obsdb_catalog.repair import repair

healthcheck(obs_db_path)   # True if healthy, else a report of what's broken
repair(obs_db_path)        # restore + reconcile + re-derive; writes atomically
```

`repair` is the only write-side recovery — invoke it deliberately, never
automatically.

The per-mission disk knowledge used by `repair`/`healthcheck` lives in `_domain.py`
(`scan`, `kind_level`). Its `level2_evt` derivation **mirrors `heasarc_utils`** (the
download producer) so a freshly recorded download and a `repair` re-derivation agree —
keep the two in sync when adding or changing a mission. Running `repair` on the live
catalogs backfills `level2_evt` for NICER/NuSTAR (the original bootstrap recorded it for
IXPE only); until then `healthcheck` reports the cleaned-event files as `uncatalogued`.

## Tests

```bash
python3 -m pytest obsdb_catalog/tests/
```

## Changelog

- 2026-06-17 — `_domain` scanners derive a representative `level2_evt` for NICER
  (`*mpu7_cl.evt*`) and NuSTAR (`event_cl/nu<obsid>A01_cl.evt*`) in addition to IXPE, so
  `repair` backfills downloaded level-2 events and stays consistent with `heasarc_utils`.
  Globs are gz-tolerant (fresh archive downloads are gzipped, reduction gunzips), preferring
  the uncompressed twin when both exist.
