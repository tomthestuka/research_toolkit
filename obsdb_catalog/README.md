# obsdb_catalog

The **single trusted writer** for a mission's data catalog
(`<obs_db_path>/catalog.json` — the inventory of what observation data is local
and how it was produced). Standalone, standard-library only: no science stack,
no `auto_research`, no sibling-toolkit import. The JSON shape is the entire
contract; a separate read side (`auto_research`) consumes it read-only.

**Prisoned write path.** `catalog.json` is mutated *only* through this package —
`neutral` for producers, `repair` for deliberate recovery. No producer, script,
or runtime-authored tool ever opens the catalog for writing. Not `*_utils`-named
and exposes no `skills/` / `TOOL_MANIFEST`, so it is invisible to LLM tool
discovery; catalog writes are a side effect of data-producing libraries, never a
model action. Its only command line is a deliberate **human/operator** maintenance
entry (`python -m obsdb_catalog …`, see Maintenance) — not an agent-callable tool.

## Layers

| Module | Role | Domain knowledge |
|---|---|---|
| `neutral` | The writer every producer calls. Durable atomic write, versioned history, `flock`-serialised, idempotent, truthful. | none (pure stdlib) |
| `repair` | Deliberate recovery: restore newest parseable snapshot, reconcile against disk, additively re-derive. Never fabricates. | yes (`_domain`) |
| `healthcheck` | Detect-only drift report; never writes. | yes (`_domain`) |
| `_domain` | Per-mission layout/filename conventions, `kind`→level, disk scanners. Imported only by `repair`/`healthcheck`. | — |

**Encrypted products.** HEASARC serves an observation still inside its proprietary
period GPG-encrypted, so a product can arrive as `<product>.gpg` — present, complete,
and unreadable. `neutral.split_lock` / `neutral.is_locked` are the one shared rule for
that wrapper: every product pattern matches the name with the wrapper split off, so an
encrypted observation is scanned as the products it holds instead of reading as empty
ground, and it is catalogued at the level its products genuinely are. **Encryption state
is derived from the recorded path, never stored** — it changes out of band, and a stored
boolean would drift out of agreement with the path beside it. A decryption therefore
surfaces like any other on-disk change: the `.gpg` path goes missing, the readable file
shows up as uncatalogued, and `repair` reconciles both. Only `.gpg` is recognised; an
unfamiliar suffix stays unknown rather than being guessed at.

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
from obsdb_catalog.repair import repair, repair_all

healthcheck(obs_db_path)   # True if healthy, else a report of what's broken
repair(obs_db_path)        # restore + reconcile + re-derive ONE mission; writes atomically
repair_all(root)           # repair EVERY mission under a root -> {mission: catalog}
```

`repair_all(root)` discovers missions the way the orchestration does — each immediate
subdirectory of `root` holding a `catalog.json` is a mission (lowercased basename);
unrelated working dirs are skipped — then runs `repair` on each. `obs_db_path` is one
mission's directory (`<root>/<MISSION>`); `root` is the parent that holds them.

`repair`/`repair_all` are the only write-side recovery — invoke them deliberately, never
automatically.

Both **narrate every change on stdout** as it is made, one line per dropped product, dropped
entry, re-derived product and level change, closed by a per-mission tally. Because repair is
the one path that deletes catalog rows, a bare `entries: 47` cannot tell an operator whether
97 rows were cleaned or nothing happened — the lines are the audit trail. They are prefixed
with the mission (so `repair_all`'s interleaved output stays attributable) and shaped to be
greppable:

```
[ixpe] - drop 02008201 qu_slices qu_slices_02008201a.csv (file missing)
[ixpe] + add  02008301 pcube     pcube_02008301.csv (on disk, uncatalogued)
[ixpe] ~ level 02008201: 3 -> 2
[ixpe] = 47 entries, 519 -> 422 products (-97 +0)
```

`-` dropped · `+` added · `~` changed · `=` tally · `!` trouble. Pass `log=None` to silence,
or any `(str) -> None` callable (e.g. a list's `append`) to capture instead of print.

Or from the command line (the same three operations, for a human/operator):

```bash
python -m obsdb_catalog healthcheck <obs_db_path>   # detect-only; exits 1 on drift
python -m obsdb_catalog repair      <obs_db_path>   # recover+reconcile ONE mission (writes)
python -m obsdb_catalog repair_all  <root>          # repair EVERY mission under a root (writes)
```

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
- 2026-06-18 — Added `repair_all(root)`: repair every mission catalog under one root, with
  missions discovered by the `catalog.json` marker (lowercased subdir basename), unrelated
  working dirs skipped. Thin loop over the per-mission `repair` primitive (unchanged); returns
  `{mission: catalog}`. Covered by `tests/test_repair.py`.
- 2026-06-18 — Added a `__main__` CLI (`python -m obsdb_catalog {healthcheck,repair,repair_all}`)
  for deliberate human/operator maintenance; `healthcheck` exits 1 on drift (prints the report).
  Pure-stdlib (argparse/json), not an agent-callable tool — the package stays invisible to LLM
  tool discovery.
- 2026-08-03 — **The scanners see GPG-encrypted products.** A proprietary-period observation
  arrives as `…fits.gz.gpg`, which every anchored product pattern rejected, so `_scan_ixpe`
  produced `products == []` and the trailing `if products:` **dropped the ObsID entirely** —
  and because `healthcheck.uncatalogued` is computed from that same scan, the check that should
  have caught it defined the problem in a way that excluded it. Verified live: `04004001` (46
  encrypted files) and `05250399` (87), both complete downloads, neither in a 47-entry catalog
  and neither seen by the scan. New `neutral.split_lock`/`is_locked`/`LOCK_SUFFIX` hold the rule
  once, in the module every layer already imports (the blindness was four independent copies of
  "what a present file looks like"); `_domain`'s patterns match through the wrapper and
  `_prefer` ranks a readable copy above an encrypted twin. NICER/NuSTAR gain a correctness fix
  in passing: their globs already *matched* a `.gpg` file and would have recorded it as if
  usable. `repair` and `healthcheck` needed no change — they read the scan. Tests:
  `tests/test_locked.py` (9).
