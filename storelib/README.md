# storelib

The toolkit's one door to a data store: S3 bytes plus a Postgres catalog. *(support package; not an
LLM tool)*

## Purpose

Tools read observation data and store the products they make through this package, never through a
path into the store. A read names an object version (a **record**,
`{"store", "bucket", "key", "version_id"}`); a write stores a **product** once per recipe, create-only,
and records it in the catalog. The orchestration (`auto_research`) issues each unit its own keys, reads
the store's own records to receipt every call, verifies what is committed (its reconciler) and moves
each (kind, subject)'s "current" pointer; this package is the tool side of that contract.

## Architecture

```
tool process (jailed; AWS_* + PG* + TOOLKIT_STORE from the runtime, AR_TOOLKIT_COMMIT + AR_ENV_FINGERPRINT beside them)
  reads   resolve(record)            HEAD by version id: size + recorded sha256
          stage_in(records, folder)  GET by version id, checked; a recipe.json record = its whole product
          remote / open_remote       s3fs version_aware, the version id in every request
  writes  produce(...)               catalog.claim(recipe)
                                       claimed  → make(work) in a temporary folder of the working folder
                                                  → PUT If-None-Match:* each declared file, 5 metadata entries
                                                  → PUT recipe.json (LAST) → catalog.commit
                                                  (any failure → catalog.fail with the reason)
                                       reuse    → the stored product, nothing made
                                       in_progress / failed / nondeterministic → raised
```

- **Layout** (the store's, domain-neutral): `<bucket>/<kind>/<subject>/<short hash>[.<generation>]/<file>`
  plus `recipe.json`. The short hash is the first 8 hex digits of the catalog's recipe hash, lengthened
  while another product of the same (kind, subject) shares it; the generation is the catalog claim's
  (omitted for 1), so a retried or taken-over production never meets an abandoned one's keys.
- **Recipe** = exactly `{kind, subject, tool, algorithm_version, settings, inputs}`: `tool` the manifest
  key, `settings` normalised (a whole float is an integer, tuples are lists: `2.0` and `2` name one
  recipe; a non-finite number is refused), `inputs` the pinned records read, in order. The toolkit
  commit and the environment fingerprint are provenance on every object and in the commit, never part
  of identity.
- **Metadata** on every object (`x-amz-meta-*`): `recipe-hash`, `sha256` (of its bytes), `producer` (the
  claim's owner), `toolkit-commit`, `env-fingerprint` (`contract.META_NAMES`).
- **`recipe.json`** (`contract.recipe_file`), written last: `{recipe_hash, generation, recipe, objects:
  [{key, version_id, sha256, size}]}` listing every other object of the product: the completeness
  marker the orchestration's reconciler requires to adopt an uncommitted product.
- **Create-only.** Every PUT carries `If-None-Match: *`. A key that already exists (an earlier attempt
  of the same generation, after a renewed live lease) is compared by sha256: equal bytes are reused;
  other bytes are a concurrency fault, failed as an ordinary failure (`catalog.fail`), so a retry claims
  a new generation and writes its own folder. Nondeterminism (one recipe, other bytes across
  productions) is the reconciler's verdict, never this helper's.
- **The store's name** comes from `TOOLKIT_STORE`, which the runtime sets beside each unit's keys and
  the CI pipeline sets for its throwaway store (`store_name()`); a product's records carry it, and a
  process with no name claims nothing. A manifest's `writes` names the store declaratively (`["db"]`);
  the launch checks it against the configured stores.
- **Checksums.** Up to 5 GiB (S3's single-PUT limit) one PUT with the object's SHA-256, which the store
  checks; above, a multipart upload in 256 MiB parts with a CRC64NVME **full-object** checksum named at
  CreateMultipartUpload (S3 takes only CRC-family checksums as full-object checksums of a multipart
  upload; computed with `awscrt`). Reads check the stored checksum (botocore's response validation)
  and the producer's sha256 and size; a staged file whose bytes already match is kept.
- **Temporary files** of a production (selected events, links) stay in its temporary folder and are
  removed with it; only the files `make` declares are stored, none of them named `recipe.json`, each
  inside the temporary folder.
- **Before the claim** every input record is resolved (HEAD by version id), so a recipe names only
  versions the store holds; a bucket outside S3's naming rules and a kind or subject holding `/`, named
  `.`/`..` or padded with spaces are refused.
- **Failure paths.** Everything after a successful claim runs under one guard: any failure gives the
  claim up with its reason. A commit whose answer is lost, or which the reconciler's adoption of this
  very production overtook, is read back from the catalog and reported as stored. A reuse of a recipe
  stored in another bucket is refused (the bucket is not part of the recipe: a tool names it in its
  settings when it differs).
- **Reads.** A single object is read again at every `stage_in`, so the store's own record of the call
  names it; a product's `recipe.json` is read again every time and its listed files are kept when the
  staged copy holds the listed bytes (their mode follows `read_only`).
- Grounding: AWS S3 User Guide "Checking object integrity" (full-object vs composite checksums);
  botocore's S3 service model (`IfNoneMatch` on PutObject and CompleteMultipartUpload; `ChecksumType`);
  s3fs `version_aware`; auto_research step-0b probes P8, P10, Q2 (versitygw 1.11.0: SHA-256 single PUT
  checked, `If-None-Match:*` → 412, 5.51 GB CRC64NVME multipart; astropy over s3fs pinned reads); the
  catalog functions in `auto_research/src/auto_research/store/catalog.sql`.

## Interfaces / Usage

```python
import storelib

out = storelib.produce(bucket="ixpe", kind="pha", subject="02008201",
                       tool="research_toolkit.<package>.<tool>", algorithm_version="1",
                       settings={"src_rad_arcmin": 1.0}, inputs=[raw_record], make=lambda work: [...])
# -> {"status": "stored" | "reused", "recipe_hash", "record": <its recipe.json record>, "objects": [...]}
[folder] = storelib.stage_in(out["record"], "staged")   # the product's files, checked, read-only
```

- `produce(*, bucket, kind, subject, tool, algorithm_version, settings, inputs, make, lease_s=3600)`
  — see Architecture. Only a manifest tool (one declaring `writes`) calls it, inside its own call;
  library code and tests never do (a test stages its inputs through the host's seed channel); this
  package's own tests are the exception, being the helper's. `make(work)` returns the files to store (paths in `work`), or `None` for every
  regular file there (symbolic links left out). Raises `ProductBusy` (another unit holds a live claim),
  `ProductRefused` (the catalog holds the recipe as failed or nondeterministic), `StoreError` (the
  store refused, no keys, `TOOLKIT_STORE`, `AR_TOOLKIT_COMMIT` or `AR_ENV_FINGERPRINT` unset: then nothing is claimed),
  or what `make` raised (the claim given up with that reason).
- `stage_in(records, folder, *, read_only=True)` → local paths: a product record gives
  `<folder>/<bucket>/<kind>/<subject>/<short hash>/`, an object record `<folder>/<bucket>/<key>`.
- `resolve(record)` → `{store, bucket, key, version_id, size, sha256}`.
- `remote(record)` → `(url, storage_options)` for astropy's `fits.open(url, use_fsspec=True,
  fsspec_kwargs=...)`; `open_remote(record)` → a seekable binary file.
- `check_record(record, *, pinned=True)`, `is_product(record)`, `product_folder(record)`, `normalise(value)`,
  `client()` (boto3 S3, path-style, from the environment), `connect()` (psycopg, autocommit, libpq `PG*`).
- `contract` — the names the orchestration's reconciler reads, mirrored here (the toolkit imports
  nothing from `auto_research`; `auto_research`'s `tests/test_toolkit_store.py` checks they are equal).

## Dependencies

`boto3`/`botocore` (+ `awscrt` for objects above 5 GiB), `psycopg` 3, `s3fs` for `open_remote`; each
imported only when a function runs, so the package imports in any Python. Environments: `envs/archive/`,
`envs/xspec/`. Used by `heasarc_utils`, `catalog_utils`, `xspec_utils`. Never imports `auto_research`.

## Reference

| Symbol | Meaning |
|---|---|
| `STORE_NAME_VARIABLE`, `store_name()` | `TOOLKIT_STORE`: the store this process's keys reach, and its name. |
| `SINGLE_PUT_MAX`, `PART_SIZE` | 5 GiB; 256 MiB. |
| `SHORT_HASH` | 8: the hex digits of a product folder's short hash. |
| `DEFAULT_LEASE_S` | 3600: a claim's lease (the catalog accepts up to a day). |
| `TOOLKIT_COMMIT_VARIABLE`, `ENV_FINGERPRINT_VARIABLE` | `AR_TOOLKIT_COMMIT`, `AR_ENV_FINGERPRINT`. |
| `StoreError`, `ProductBusy`, `ProductRefused` | See Interfaces. |
| `contract.META_*`, `META_NAMES`, `RECIPE_FILE`, `RECIPE_KEYS`, `RECORD_KEYS`, `RECIPE_FILE_KEYS`, `RECIPE_OBJECT_KEYS` | The product contract. |
| `contract.recipe_file`, `contract.check_recipe_file`, `ContractError` | Build and check a `recipe.json` body. |

Tests: `python3 -m pytest storelib/tests/` — the value checks run anywhere; the rest take the root
conftest's `store` fixture (a throwaway store the host gives through `TOOLKIT_STORE` and its keys;
auto_research's CI pipeline, or its `tests/test_toolkit_store.py`).

## Changelog

- 2026-11-10 — Created (auto_research Phase 3 step 0b, S8): `produce`, `stage_in`, `resolve`, `remote`,
  `open_remote`, the mirrored product contract; replaces `obsdb_catalog`.
- 2026-11-10 — Stateless review fixes: inputs resolved before the claim; bucket, kind and subject
  checked; a nested `recipe.json` refused; the claim given up for any failure after it (folder and
  metadata included); a lost or overtaken commit read back; a reuse in another bucket refused; the
  product's own `recipe.json` picked by its place (`<kind>/<subject>/<hash>/recipe.json`); a single
  object read again at every staging.
- 2026-11-10 — Rulings R3/R5: other bytes at a key of the same generation fail the production as an
  ordinary failure (a retry claims a new generation); the store's name is read from `TOOLKIT_STORE`
  (`store_name()`, `STORE_NAME_VARIABLE`); `produce` takes no `store` and `DEFAULT_STORE` is gone.
- 2026-11-10 — Only a manifest tool calls `produce` (auto_research ruling); a test's store inputs come
  through the seed channel (the root conftest's `seed`).
