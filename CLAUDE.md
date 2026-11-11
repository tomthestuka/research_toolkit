# CLAUDE.md

**Single source of truth for build-time coding agents** (Claude Code, Codex, any harness)
working on `research_toolkit`. `AGENTS.md` is a symlink to this file ([§8](#8-file-linking)).

> This governs the agents that **construct and maintain this repo's code**. It is *not* for
> the runtime LLM agents that *use* the toolkit — they see only each tool's
> `skills/TOOL_MANIFEST`, never this file or any `README.md`
> ([§4](#4-coding-agent-contract)).

---

## 1. What this repo is

`research_toolkit` is the **domain tool library** for X-ray spectro-polarimetric analysis
(IXPE / NICER / NuSTAR, built on HEASoft/XSPEC/ixpeobssim). It is consumed **read-only** by
the sibling orchestration repo `auto_research`, whose runtime analysts discover and invoke
these tools in a sandbox.

**Trust boundary (one-directional, hard rule):** nothing here may `import auto_research`.
This repo is self-contained; the cross-repo contract is the toolkit contract (manifests, envs, tests)
and the data store's product contract ([§6](#6-the-data-store)), which `storelib/contract.py` mirrors. Keeping the dependency one-way is what lets the orchestration
side treat this repo as an untrusted, swappable tool source.

---

## 2. Package taxonomy

Two kinds of top-level package live here, and the distinction matters:

| Kind | Examples | Discoverable as an LLM tool? | Built by |
|---|---|---|---|
| **`<domain>_utils/` toolkits** | `xspec_utils`, `heasarc_utils`, `catalog_utils` (exist); `lib_utils` (planned) | **Yes** — via `skills/TOOL_MANIFEST` | human + (Phase-2) runtime tool-author, on a branch |
| **Plain support packages** | `storelib` | **No** — deliberately invisible | human / build-time agent, by hand |

- **`*_utils` toolkits** are the AI-facing surface. Their name and `TOOL_MANIFEST` are how
  the orchestration repo finds tools. Follow the uniform layout in
  [§5](#5-the-_utils-toolkit-layout).
- **`storelib`** is a reliability-bearing support package — the tools' *one door* to the data
  store: version-pinned reads and create-once products recorded in the catalog. It is
  intentionally **not** `*_utils`-named and exposes **no `skills/` / `TOOL_MANIFEST`**, so tool
  discovery never surfaces it. Store writes are a side effect of data-producing tools (a manifest
  `writes`), never a model action. See its `README.md`.

---

## 3. `README.md` is the module briefing (the house style)

Every package keeps a `README.md` that is its **briefing and current ground truth** — the
analogue of `DOCS.md` in `auto_research`. Read the `README.md` of any package before you
touch it; trust it over assumptions. Sections, in order:

1. **Purpose** — what the package is for, in plain language.
2. **Architecture** — how it is designed and why; how it fits the repo.
3. **Interfaces / Usage** — the public surface (functions, signatures, invariants) and how a
   caller uses it.
4. **Dependencies** — what it imports and what depends on it.
5. **Reference** — every public symbol with meaning, params, return (may be folded into §3
   for small packages).
6. **Changelog** — dated entries at the **end**, always **appended at the bottom**.

The **root `README.md`** is the repo-level briefing. Keep all of them **synced**: any code
change updates the relevant `README.md` in the **same session**, including a changelog line.
An out-of-sync README is a bug.

---

## 4. Coding-agent contract

**Start of session:** read this `CLAUDE.md`, then the `README.md` of every package you will
touch. State briefly what you intend to change and which READMEs you will update.

**Hard rules:**

- **READMEs and `CLAUDE.md` are build-time only.** Never let them enter a runtime LLM's
  context. The runtime's view of a tool is its `skills/TOOL_MANIFEST` entry (and the `lib/`
  code it wraps) — nothing else. If you want the runtime to "know" something, it belongs in
  the `TOOL_MANIFEST` description or the tool's `args_schema`, not in a README.
- **Keep READMEs synced** ([§3](#3-readmemd-is-the-module-briefing-the-house-style)).
- **Never `import auto_research`** ([§1](#1-what-this-repo-is)).
- **Tests ship with the code.** An untested tool/package is not mergeable; put tests under
  `<package>/tests/` and run `python3 -m pytest <package>/tests/`.
- **Minimal code.** Prefer the fewest modules/functions/lines that satisfy the interface;
  for well-trodden primitives prefer a stable community wheel over hand-rolling.
- **Ground every line in an authoritative source — no guessing.** Every line you write and
  every piece of logic must be backed by authoritative documentation, and by a known-good
  example where one exists. Before you call any external API (HEASoft/XSPEC/ixpeobssim/
  astropy/…) or rely on a standard-library or domain behavior, consult its official docs —
  and an existing example (the current toolkit code, an upstream cookbook/tutorial) — and
  write **only** what the docs actually support. Do **not** write from a remembered
  signature, a plausible-looking flag, or "it should work like this": this repo wraps
  fragile science-stack APIs that fail silently or diverge when misused. If you cannot ground
  a function, argument, or behavior in an authoritative source, stop — find the source or ask
  the human — rather than guess. When the grounding is non-obvious, cite it (doc URL / section
  / example path) in the `README.md` or a code comment.
- **No hardcoded secrets**; the science stack (HEASoft/XSPEC/ixpeobssim) is imported **only**
  inside `skills/session.py` / the REPL subprocess, never at package import time.
- **Git:** you run git only when the human explicitly asks. The human merges to the trusted
  branch out of band; the Phase-2 runtime tool-author lands **one atomic commit per tool on a
  working branch** and never touches the trusted branch.

---

## 5. The `*_utils` toolkit layout

`xspec_utils` is the template every new toolkit mirrors exactly:

```
<domain>_utils/
  __init__.py          # re-exports TOOL_MANIFEST for discovery; puts lib/ on sys.path.
                       #   NO heavy/science imports here — must load in the orchestration env.
  lib/                 # ── HUMAN-HALF: human-authored, importable, runnable domain code ──
    <mission>.py       #   thin per-mission strategy modules (NICER.py, NuSTAR.py, red_ixpe.py)
    <helpers>.py       #   small shared helpers (find_src.py, src_pos.py, …)
    <domain>_help.py   #   the high-level, ABSTRACTED entry API a person calls (load_data, …)
  skills/              # ── AI-HALF: the runtime tool interface ──
    __init__.py        #   TOOL_MANIFEST — PURE DATA (no science imports); each entry wraps a
                       #   lib/ fn: {module, fn, description, args_schema}
    session.py         #   science-stack-dependent wrappers; imported ONLY in the REPL subprocess
  tool_config.yaml     #   per-package runtime config: per-tool call timeouts, and `env: <id>` —
                       #   the software environment the package's tools run in (§6b)
  README.md            #   the module briefing (§3)
  tests/               #   tests that ship with the toolkit
```

Invariants (they are *why* the split exists):
- **`skills/__init__.py` (`TOOL_MANIFEST`) is pure data** — no XSPEC/HEASoft import — so the
  orchestration env (no science stack) can discover it. Entries point at `lib/` functions;
  domain code is never duplicated into `skills/`.
- **`lib/` is the human-half** — importable/runnable by a person; its `*_help.py` entry points
  are abstracted on purpose so error-prone choices (data-group order, gain linking, response
  paths) live inside tested code, not in per-call LLM arguments.
- **`session.py` holds science-stack code**, imported only inside the REPL subprocess.
- A toolkit's `lib/` may be **synced/mirrored from a live working directory** (e.g.
  `xspec_utils/lib/` via a sync script — see the package docstring). Where so, **edit the
  upstream source, not the mirror**, and re-sync.

---

## 6. The data store

Observation data and the products made from it live in a data store (`auto_research`'s store `db`):
bytes behind the S3 interface plus a Postgres catalog of products, their objects, a "current" pointer
per (kind, subject), credentials and events. The orchestration runs it, issues each unit its own keys
(`AWS_*`, `PG*`) for the run, receipts every call from the store's own records, verifies committed
products and moves "current". The toolkit's side:

- **Tools are the only interface.** A tool reaches the store only through `storelib`; no tool takes a
  path into the store. Inputs are **records** (`{"store", "bucket", "key", "version_id"}`, argument
  kind `record`); every read names a version.
- **Layout** (`storelib/README.md`): `<bucket>/<kind>/<subject>/<short hash>[.<generation>]/<file>` plus
  a `recipe.json` written last. Here bucket = the mission as HEASARC names it (`ixpe`, `nicer`,
  `nustar`), subject = the ObsID, kind = `raw` (an archive download, laid out as HEASARC's ObsID
  folder) or `pha` (IXPE Stokes spectra); the root `README.md` is the briefing.
- **A product is made once per recipe** — tool (manifest key), `algorithm_version` (manifest; bump it
  when what the tool stores changes), settings, input records — through `storelib.produce`: claim,
  create-only uploads with five metadata entries, `recipe.json`, commit (or fail with the reason). A
  tool that stores products declares `writes: ["db"]` and `algorithm_version` in its manifest entry.
- **What a tool writes is one of three classes:**
  - a **science artifact** — anything that depends on a `parameter` or a setting outside the
    mission team's defaults — goes under an argument of kind `output` (its `shape` folder, file or
    prefix) and is returned under `products` (`{"path", "caption"}`) or, for a PNG, `figures`, so the
    runtime receipts it;
  - a **store product** — raw archive data and the mission team's standard reductions — goes into
    the store through `storelib.produce`, its record returned;
  - **scratch** — a program's own bookkeeping, a reduction's temporary files — stays in the call's
    working folder (a production's temporary folder is removed with it).

---

## 6b. Software environments (`envs/`)

A package names the environment its tools run in, in its `tool_config.yaml` (`env: <id>`). The
definition is a folder of plain files at `envs/<id>/`:

```
envs/<id>/
  env.yaml      # purpose (one line), base (a human environment to layer on, or none), externals,
                #   probe (one line of Python the launch runs in the default environment; optional)
  env.lock      # conda explicit list: @EXPLICIT, then a conda-forge URL#md5 per package — exact
  env.pip       # wheels conda-forge lacks: name==version --hash=sha256:<hex> per line (optional)
  activate.sh   # variables the tools need; sourcing a declared external's init (optional)
  boot.py       # analyst REPL session setup, run once at REPL start (optional)
```

- **Built** (a lock, no base) or **derived** (a human `base`, no lock). A base is named, never a
  path: the host declares each human conda environment by name (`auto_research` config
  `human_envs`). Externals (HEASoft, CALDB) are the human's declared software outside conda
  (`auto_research` config `env_externals`), stamped on every receipt, never rebuilt; an
  `activate.sh` sources an external's init by its absolute path.
- **The default environment** (`auto_research`'s toolkit contract `default_env`): `envs/analysis/`,
  derived on the human `astropy` environment with HEASoft and CALDB. A unit's own REPL and shell run
  in it, and so do packages that name no environment. It has no
  `boot.py`: none of its tools needs an XSPEC session (`envs/xspec/boot.py` holds that), and the shell
  never runs a boot file. Its `probe` (`import xspec, heasarc_utils, catalog_utils`) runs at launch,
  in the jail when it is on.
- `auto_research` builds a definition, jailed, into `<runtime_dir>/conda/envs/ar-<id>` — never
  into this repo or a human's conda environments — runs the package's tests in it at the gate, and
  runs its tools in it for the analysts. Several packages may share one environment; a package
  names exactly one.
- The runtime coder writes definitions through its environment tools (`env_solve`,
  `env_pin_wheels`, `cand_env_activate`; it never runs conda). A human writing one by hand uses the
  same functions (`auto_research.sandbox.env_build.solve` / `pin_wheels`) so the lock is exact.
- Give paths in `activate.sh` and `boot.py` absolutely: the runtime sources a copy.
- A tool routed to a pipeline environment runs only in the jail (`agent_os_jail: true`); every
  path `activate.sh` or `boot.py` reads must be on the jail's allow-list or a declared external.
- `envs/xspec/` serves `xspec_utils`: Python 3.13 and numpy 2.x match HEASoft 6.34's PyXspec
  build (compiled against Python 3.13 / numpy 2.2), so its lock pins them; it also holds the store's
  clients (boto3 with awscrt, s3fs, psycopg).
- `envs/archive/` serves `heasarc_utils` and `catalog_utils` (built: astroquery, astropy, boto3 with
  awscrt, s3fs, psycopg, pytest). The derived default environment cannot gain packages and the human
  `astropy` environment lacks the store's clients, so these packages name a built one, named after
  what it reaches (the archive and the store).

---

## 7. Tech & conventions

- **Python ≥ 3.10**; standard scientific stack for `*_utils` toolkits; a support package imports
  its third-party clients only inside the function that uses them (`storelib`: boto3, psycopg,
  s3fs), so it imports in any Python.
- **Tests:** `pytest`, deterministic only; live science-stack runs are exercised by hand, not
  asserted. `python3 -m pytest <package>/tests/`. The root `conftest.py` defines four classes:
  unmarked (fast, deterministic, seconds), `slow` (real data or a real model, minutes; runs by
  default, `-m "not slow"` gives the quick suite), `live` (network or an external service;
  skipped unless `--live` or `TOOLKIT_LIVE=1`) and `data("<name>")` (reads a named data set from
  the folder `TOOLKIT_DATA_<NAME>` names, never a hardcoded host path; skipped with the reason when
  that is unset or absent). The authoring gate deselects `live` (it has no network), sets
  `TOOLKIT_DATA_<NAME>` for each data set its host's `test_data` names, and fails a change whose
  data test it skipped, and the launch refuses a data set no host declares. Live data sets:
  `cygx1_ixpe`, `cygx1_nustar`. A test of the change that the gate skips or marks xfail fails it too.
- **Mission data** is in the data store (§6), named by records; a tool writes only what its class
  allows (§6): artifacts under its output argument, store products through `storelib.produce`,
  scratch in its working folder. A test taking the root conftest's `store` fixture reaches a throwaway
  store the host gives (`TOOLKIT_STORE` and its keys); without one it skips with the reason, and the
  authoring gate (which gives one) fails a store test it skipped. A test stages its store inputs only
  through the host's seed channel (the `seed` fixture); only a manifest tool's call writes the store,
  and library code never calls `storelib.produce`. A tool reading the store without a `record`
  argument declares `reads: ["db"]`. No tool calls another manifest tool: shared logic lives in `lib/`.

---

## 8. File linking

- This `CLAUDE.md` at the repo root is the single source of truth.
- Other agent-instruction files are **symlinks** to it (`ln -s CLAUDE.md AGENTS.md`). Never
  fork the content — one file, many links.

---

## Changelog

- 2026-06-17 — **Created `CLAUDE.md`.** Recorded the repo orientation, the two-repo trust
  boundary, the package taxonomy (`*_utils` toolkits vs. plain support packages), the
  README-as-module-briefing + sync convention (mirroring `auto_research`'s `DOCS.md`), the
  `*_utils` layout template, and the data-catalog write boundary. Added alongside the new
  `obsdb_catalog` package (the single trusted catalog writer + `repair`/`healthcheck` + tests).
- 2026-06-17 — **§4: added the "ground every line in an authoritative source — no guessing"
  rule.** Code must be backed by official docs (and a known-good example where one exists);
  no remembered signatures or plausible-looking flags, because this repo wraps fragile
  science-stack APIs that fail silently when misused.
- 2026-06-17 — **Built `heasarc_utils`** (the acquire toolkit) and marked it in §2. A
  mission-agnostic `astroquery.heasarc` engine + `MissionProfile` registry (ixpe/nicer/nustar,
  extensible) downloads an ObsID and records its level-2 event via `obsdb_catalog.neutral.new`
  — the first producer wired into the catalog. Extended `obsdb_catalog._domain` to derive a
  representative `level2_evt` for NICER/NuSTAR (mirrors `heasarc_utils`) so `repair` backfills
  downloaded level-2 events.
- 2026-06-17 — `heasarc_utils` live-verified against installed astroquery (NICER acquire,
  `logs/dev/live_acquire_nicer.md`): resolve uses `query_tap` (not `query_region(column_filters)`,
  which drops the obsid filter); fresh archive event files are gzipped → level-2 globs made
  gz-tolerant in `heasarc_utils` + `obsdb_catalog._domain`; `download_data` tar is obsid-rooted.
- 2026-10-07 — **Environments and test classes (auto_research env plan T1).** §5 layout gains
  `tool_config.yaml`'s `env:`; new §6b, the `envs/<id>/` contract; `envs/xspec/` defined for
  `xspec_utils` (`env: xspec`): Python 3.13 / numpy 2.3 / astropy 7.1 matching HEASoft 6.34's
  PyXspec build, `ixpeobssim` 33.0.0 as a hash-pinned wheel, HEASoft/CALDB activation and the
  XSPEC session boot moved in from `auto_research`'s `analyst_env.sh` / `repl_init.py`. Root
  `conftest.py`: `slow` and `live` markers, `live` opt-in (§7).
- 2026-11-04 — **Obs-db as a clone of HEASARC's archive; the three output classes.** §2 lists
  `catalog_utils` and `obsdb_catalog`'s resolver and human CLI. §6: an obs-db root of
  HEASARC-named mission folders laid out as HEASARC's tree (`_domain.TREES`), `locate`, product
  paths relative to the mission folder, producers `heasarc_utils.acquire` (into the obs-db root) and
  `xspec_utils.red_ixpe.pha`, and what a tool writes (artifacts, obs-db items, scratch). §7: the
  mission data paths and what a tool writes there.
- 2026-11-09 — §6b: `envs/analysis/` is the default environment (derived on the human `astropy`
  environment, HEASoft and CALDB; no `boot.py`), moved here from `auto_research`'s config and its
  `analyst_env.sh`/`repl_init.py`; a base is a name the host declares (`human_envs`); `env.yaml`
  takes an optional one-line `probe`.
- 2026-11-09 — **The `data` test class (auto_research Phase 3 step 0a W3).** §7: `data("<name>")` in the root `conftest.py` (root from `TOOLKIT_DATA_<NAME>`, skipped with the reason when absent); `xspec_utils/tests/test_outputs.py` reads Cyg X-1's IXPE and NuSTAR folders through `data("cygx1_ixpe")` / `data("cygx1_nustar")` instead of hardcoded paths. The gate deselects `live` and fails a change whose data test it skipped.
- 2026-11-10 — **The data store** (auto_research Phase 3 step 0b, S8). §1: the cross-repo contract
  names the store's product contract. §2: `storelib` replaces `obsdb_catalog` (deleted; `catalog.json`
  retired, not migrated in code). §6 rewritten: the data store (tools the only interface, records,
  the layout and its meaning here, one product per recipe, the three output classes with store
  products). §6b: `envs/archive/` for `heasarc_utils` and `catalog_utils`; `envs/xspec/` gains the store
  clients. §7: mission data is in the store; the `store` test fixture.
- 2026-11-10 — §7 (auto_research rulings): test inputs only through the seed channel; only a manifest
  tool writes the store, library code never calls `storelib.produce`; manifest `reads`; no tool calls
  another manifest tool.
