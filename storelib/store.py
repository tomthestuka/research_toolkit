"""The toolkit's one door to a data store: S3 bytes plus a Postgres catalog.

A tool reaches a store only through these functions, with the keys the runtime hands its process as
environment variables (``AWS_ACCESS_KEY_ID``, ``AWS_SECRET_ACCESS_KEY``, ``AWS_ENDPOINT_URL_S3``,
``AWS_REGION`` for S3; ``PGHOST``, ``PGPORT``, ``PGUSER``, ``PGPASSWORD``, ``PGDATABASE`` for the
catalog, read by libpq itself), and ``AR_TOOLKIT_COMMIT`` / ``AR_ENV_FINGERPRINT`` for the provenance a
product carries. It reads no credentials file and no configuration.

Reading (every read names a version):

* :func:`resolve` — a record's object version: its size and recorded sha256 (HEAD by version id);
* :func:`stage_in` — local, checked copies of records; a record of a product's :data:`RECIPE_FILE`
  stages every object it lists, each at its listed version and checked against its listed sha256;
* :func:`open_remote` — a version-pinned remote file (s3fs, ``version_aware``) for header reads.

Writing (:func:`produce`): claim the recipe in the catalog; on ``claimed`` run the caller's ``make``
into a temporary folder in the working folder; upload each output create-only (``If-None-Match: *``)
under ``<kind>/<subject>/<short hash>[.<generation>]/`` in the bucket, with the five metadata entries;
write :data:`RECIPE_FILE` last; commit. A failure on the way gives the claim up (``catalog.fail``).

Checksums: a single PUT (up to :data:`SINGLE_PUT_MAX`, S3's single-PUT limit of 5 GiB) carries the
object's SHA-256, which the store checks; a larger object goes up in parts with a CRC64NVME full-object
checksum named at CreateMultipartUpload (S3 allows only CRC-family checksums as full-object checksums
of a multipart upload), computed with ``awscrt``. Grounding: AWS S3 User Guide, "Checking object
integrity" (full-object and composite checksums); botocore's S3 service model (``ChecksumType``,
``IfNoneMatch`` on ``PutObject`` and ``CompleteMultipartUpload``); probe P8/Q2 in
``auto_research``'s step-0b probes (versitygw 1.11.0: SHA-256 single PUT checked, ``If-None-Match:*``
→ 412, 5.51 GB CRC64NVME multipart).
"""

from __future__ import annotations

import base64
import hashlib
import json
import math
import os
import re
import tempfile
from pathlib import Path, PurePosixPath

from . import contract as c

#: The variable naming the store this process's keys reach (the runtime and the CI pipeline set it):
#: the name a product's records carry.
STORE_NAME_VARIABLE = "TOOLKIT_STORE"
#: S3's largest single PUT (5 GiB); a larger object is uploaded in parts.
SINGLE_PUT_MAX = 5 * 1024 ** 3
#: The part size of a multipart upload (S3: 5 MiB to 5 GiB per part, at most 10,000 parts).
PART_SIZE = 256 * 1024 ** 2
#: The hex digits of a product folder's short hash (lengthened on a clash with another recipe).
SHORT_HASH = 8
#: The default claim lease, in seconds (catalog.claim accepts up to one day).
DEFAULT_LEASE_S = 3600
_CHUNK = 1 << 20
_S3_VARIABLES = ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_ENDPOINT_URL_S3")
_PG_VARIABLES = ("PGHOST", "PGUSER", "PGDATABASE")
#: An S3 bucket name (lower-case letters, digits, dots and hyphens; 3 to 63 characters).
_BUCKET = re.compile(r"^[a-z0-9][a-z0-9.-]{1,61}[a-z0-9]$")
TOOLKIT_COMMIT_VARIABLE = "AR_TOOLKIT_COMMIT"
ENV_FINGERPRINT_VARIABLE = "AR_ENV_FINGERPRINT"


class StoreError(RuntimeError):
    """The store or its catalog refused, or the environment names no store."""


class ProductBusy(StoreError):
    """Another unit holds a live claim on the same recipe."""


class ProductRefused(StoreError):
    """The catalog holds the recipe as failed or nondeterministic: a human decides."""


class _KeyTaken(StoreError):
    """The store already holds other bytes at a key this production writes: an earlier attempt of the
    same generation (a renewed live lease after a killed call). A concurrency fault, failed as an
    ordinary failure so a retry claims a new generation and its own folder."""


# ------------------------------------------------------------------ connections


def _require(names) -> None:
    missing = [n for n in names if not os.environ.get(n)]
    if missing:
        raise StoreError(f"no store is reachable from this process: {', '.join(missing)} unset (the runtime "
                         "sets them for a unit holding store keys)")


def client():
    """A boto3 S3 client from the environment: path-style addressing, the request checksum sent only
    where this module names one, every response checksum the client knows checked."""
    _require(_S3_VARIABLES)
    import boto3
    from botocore.config import Config

    return boto3.client("s3", region_name=os.environ.get("AWS_REGION") or "us-east-1",
                        config=Config(s3={"addressing_style": "path"}, retries={"max_attempts": 3},
                                      request_checksum_calculation="when_required",
                                      response_checksum_validation="when_supported"))


def connect():
    """An autocommit psycopg connection to the store's catalog, as libpq's ``PG*`` variables name it."""
    _require(_PG_VARIABLES)
    import psycopg

    return psycopg.connect("", autocommit=True)


def store_name() -> str:
    """The name of the store this process's keys reach (:data:`STORE_NAME_VARIABLE`); raises
    :class:`StoreError` when the runtime set none."""
    name = os.environ.get(STORE_NAME_VARIABLE, "").strip()
    if not name:
        raise StoreError(f"no store is named for this process: {STORE_NAME_VARIABLE} unset (the runtime sets it "
                         "beside the store keys)")
    return name


def _provenance() -> tuple[str, str]:
    values = tuple(os.environ.get(n, "").strip() for n in (TOOLKIT_COMMIT_VARIABLE, ENV_FINGERPRINT_VARIABLE))
    missing = [n for n, v in zip((TOOLKIT_COMMIT_VARIABLE, ENV_FINGERPRINT_VARIABLE), values) if not v]
    if missing:
        raise StoreError(f"a product records the toolkit commit and the environment it was made in: "
                         f"{', '.join(missing)} unset (the runtime sets them beside the store keys)")
    return values


# ------------------------------------------------------------------ values


def normalise(value):
    """``value`` as the JSON a recipe hashes: a whole float becomes an integer (``2.0`` and ``2`` name
    one setting), tuples become lists, numpy scalars their Python values. Raises ``ValueError`` for a
    non-finite number, a non-text mapping key, or a value JSON cannot hold."""
    if hasattr(value, "item") and callable(value.item) and not isinstance(value, (list, tuple, dict)):
        value = value.item()  # a numpy scalar
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"a recipe holds finite numbers only, not {value!r}")
        return int(value) if value.is_integer() else value
    if isinstance(value, dict):
        bad = [k for k in value if not isinstance(k, str)]
        if bad:
            raise ValueError(f"a recipe's mapping keys are text, not {bad[0]!r}")
        return {k: normalise(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [normalise(v) for v in value]
    raise ValueError(f"a recipe holds JSON values, not {type(value).__name__}")


def check_record(record, *, pinned: bool = True) -> dict:
    """``record`` as ``{store, bucket, key, version_id}``, each non-empty text. With ``pinned`` the
    version id is required (every read names a version; the runtime pins a record before the call)."""
    keys = set(c.RECORD_KEYS) if pinned else set(c.RECORD_KEYS) - {"version_id"}
    if not isinstance(record, dict) or not keys <= set(record) or not set(record) <= set(c.RECORD_KEYS) or any(
            not isinstance(record[k], str) or not record[k].strip() for k in record):
        raise ValueError(f"a record is {{\"store\", \"bucket\", \"key\", \"version_id\"}}, each non-empty text; "
                         f"got {record!r}")
    if not _BUCKET.match(record["bucket"]) or ".." in record["bucket"]:
        raise ValueError(f"a record's bucket is an S3 bucket name, not {record['bucket']!r}")
    return {k: record[k] for k in c.RECORD_KEYS if k in record}


def is_product(record) -> bool:
    """Whether ``record`` names a product's :data:`RECIPE_FILE` (a whole product), not one object."""
    return PurePosixPath(record["key"]).name == c.RECIPE_FILE


def product_folder(record) -> str:
    """The folder key of the product whose :data:`RECIPE_FILE` ``record`` names."""
    return str(PurePosixPath(record["key"]).parent)


def _b64(digest: bytes) -> str:
    return base64.b64encode(digest).decode("ascii")


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(_CHUNK):
            h.update(chunk)
    return h.hexdigest()


def _code(exc) -> str:
    return str(getattr(exc, "response", {}).get("Error", {}).get("Code", ""))


# ------------------------------------------------------------------ reads


def resolve(record, *, s3=None) -> dict:
    """The object version ``record`` names: ``{store, bucket, key, version_id, size, sha256}``, the
    sha256 as the producer recorded it in the object's metadata (``None`` when it carries none).
    Raises :class:`StoreError` when the store holds no such version."""
    record = check_record(record)
    s3 = s3 or client()
    from botocore.exceptions import ClientError

    try:
        head = s3.head_object(Bucket=record["bucket"], Key=record["key"], VersionId=record["version_id"])
    except ClientError as exc:
        raise StoreError(f"{record['bucket']}/{record['key']} version {record['version_id']}: "
                         f"{_code(exc) or exc}") from exc
    return {**record, "size": int(head["ContentLength"]),
            "sha256": (head.get("Metadata") or {}).get(c.META_SHA256)}


def _local_name(base: Path, bucket: str, key: str) -> Path:
    parts = PurePosixPath(key).parts
    if not parts or any(p in ("", ".", "..") for p in parts) or PurePosixPath(key).is_absolute():
        raise StoreError(f"key {key!r} does not name a file below its bucket")
    return base.joinpath(bucket, *parts)


def _download(s3, bucket: str, key: str, version_id: str, target: Path, *, sha256: str | None,
              size: int | None, read_only: bool, keep: bool = False) -> str:
    """GET one version into ``target`` (via a ``.part`` file), hashing it on the way; the hash and the
    size must equal ``sha256`` and ``size`` when given. With ``keep``, a file already there with
    exactly those bytes is kept (its mode set) instead. Returns the sha256."""
    if keep and target.is_file() and sha256 is not None and size is not None and target.stat().st_size == size:
        if _sha256_file(target) == sha256:
            os.chmod(target, 0o444 if read_only else 0o644)
            return sha256  # staged earlier with these bytes
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_name(target.name + ".part")
    body = s3.get_object(Bucket=bucket, Key=key, VersionId=version_id)["Body"]
    h, n = hashlib.sha256(), 0
    try:
        with open(part, "wb") as fh:
            for chunk in body.iter_chunks(_CHUNK):
                h.update(chunk)
                n += len(chunk)
                fh.write(chunk)
    finally:
        body.close()
    digest = h.hexdigest()
    if (sha256 is not None and digest != sha256) or (size is not None and n != size):
        part.unlink(missing_ok=True)
        raise StoreError(f"{bucket}/{key} version {version_id}: the bytes read ({n} B, sha256 {digest}) "
                         f"differ from the record ({size} B, sha256 {sha256})")
    if read_only:
        os.chmod(part, 0o444)
    os.replace(part, target)
    return digest


def stage_in(records, folder, *, read_only: bool = True, s3=None) -> list[Path]:
    """Local copies of ``records`` (one record or a list) below ``folder``, in order: a product's
    :data:`RECIPE_FILE` record gives the product's folder (``<folder>/<bucket>/<kind>/<subject>/<short
    hash>/``) holding every object it lists, each fetched at its listed version and checked against its
    listed size and sha256; any other record gives its file (``<folder>/<bucket>/<key>``), checked
    against the sha256 its producer recorded. Files are read-only unless ``read_only`` is false. A
    product's file already staged with its listed bytes is kept: the recipe file is read again every
    time, so the store's record of the call names the product and, through it, every listed version. A
    single object is read again every time."""
    records = [records] if isinstance(records, dict) else list(records)
    base = Path(folder)
    s3 = s3 or client()
    out = []
    for record in records:
        record = check_record(record)
        bucket = record["bucket"]
        if is_product(record):
            target = _local_name(base, bucket, record["key"])
            _download(s3, bucket, record["key"], record["version_id"], target, sha256=None, size=None,
                      read_only=read_only)
            try:
                listing = c.check_recipe_file(json.loads(target.read_text(encoding="utf-8")))
            except (ValueError, c.ContractError) as exc:
                raise StoreError(f"{bucket}/{record['key']}: not a {c.RECIPE_FILE} ({exc})") from exc
            prefix = product_folder(record) + "/"
            for item in listing["objects"]:
                if not item["key"].startswith(prefix):
                    raise StoreError(f"{bucket}/{record['key']} lists {item['key']}, outside its folder")
                _download(s3, bucket, item["key"], item["version_id"], _local_name(base, bucket, item["key"]),
                          sha256=item["sha256"], size=item["size"], read_only=read_only, keep=True)
            out.append(target.parent)
        else:
            known = resolve(record, s3=s3)
            target = _local_name(base, bucket, record["key"])
            _download(s3, bucket, record["key"], record["version_id"], target, sha256=known["sha256"],
                      size=known["size"], read_only=read_only)
            out.append(target)
    return out


def remote(record) -> tuple[str, dict]:
    """``(url, storage_options)`` of ``record``'s version for fsspec readers, e.g. astropy's
    ``fits.open(url, use_fsspec=True, fsspec_kwargs=storage_options)`` (s3fs ``version_aware``; the
    version id rides in the URL, so every request names it)."""
    record = check_record(record)
    _require(_S3_VARIABLES)
    options = {"version_aware": True, "endpoint_url": os.environ["AWS_ENDPOINT_URL_S3"],
               "config_kwargs": {"s3": {"addressing_style": "path"}}}
    return f"s3://{record['bucket']}/{record['key']}?versionId={record['version_id']}", options


def open_remote(record):
    """A binary, seekable file of ``record``'s version, read lazily over S3 (s3fs, ``version_aware``)."""
    import s3fs

    url, options = remote(record)
    record = check_record(record)
    fs = s3fs.S3FileSystem(**options)
    return fs.open(f"{record['bucket']}/{record['key']}", "rb", version_id=record["version_id"])


# ------------------------------------------------------------------ writes


def _short(conn, recipe_hash: str, kind: str, subject: str) -> str:
    """The shortest prefix of ``recipe_hash``, from :data:`SHORT_HASH` digits, that no other product of
    ``(kind, subject)`` shares."""
    others = [r[0] for r in conn.execute(
        "SELECT recipe_hash FROM catalog.product WHERE kind = %s AND subject = %s AND recipe_hash <> %s",
        (kind, subject, recipe_hash)).fetchall()]
    n = SHORT_HASH
    while any(o.startswith(recipe_hash[:n]) for o in others):
        n += 1
    return recipe_hash[:n]


def _existing(s3, bucket: str, key: str) -> tuple[str, str | None, int]:
    """``(version id, recorded sha256, size)`` of the latest version at ``key``."""
    listing = s3.list_object_versions(Bucket=bucket, Prefix=key)
    latest = [v for v in listing.get("Versions", []) if v["Key"] == key and v.get("IsLatest")]
    if not latest:
        raise StoreError(f"{bucket}/{key}: the store refused a create-only write and lists no version there")
    version = latest[0]["VersionId"]
    head = s3.head_object(Bucket=bucket, Key=key, VersionId=version)
    return version, (head.get("Metadata") or {}).get(c.META_SHA256), int(head["ContentLength"])


def _put(s3, bucket: str, key: str, path: Path, meta: dict) -> dict:
    """Store ``path`` at ``key`` create-only with ``meta`` (its sha256 added); ``{key, version_id,
    sha256, size}``. A key that already holds the same bytes (this production's earlier attempt) is
    reused; other bytes there raise :class:`_KeyTaken`."""
    from botocore.exceptions import ClientError

    size = path.stat().st_size
    if size <= SINGLE_PUT_MAX:
        digest = hashlib.sha256()
        with open(path, "rb") as fh:
            while chunk := fh.read(_CHUNK):
                digest.update(chunk)
        sha = digest.hexdigest()
        metadata = {**meta, c.META_SHA256: sha}
        try:
            with open(path, "rb") as fh:
                answer = s3.put_object(Bucket=bucket, Key=key, Body=fh, IfNoneMatch="*", Metadata=metadata,
                                       ChecksumAlgorithm="SHA256", ChecksumSHA256=_b64(digest.digest()))
        except ClientError as exc:
            if _code(exc) not in ("PreconditionFailed", "412"):
                raise
            return _reuse_or_refuse(s3, bucket, key, sha, size)
    else:
        sha, answer = _put_parts(s3, bucket, key, path, meta)
        if answer is None:
            return _reuse_or_refuse(s3, bucket, key, sha, size)
    return {"key": key, "version_id": answer["VersionId"], "sha256": sha, "size": size}


def _reuse_or_refuse(s3, bucket, key, sha, size) -> dict:
    version, stored, stored_size = _existing(s3, bucket, key)
    if stored != sha or stored_size != size:
        raise _KeyTaken(f"{bucket}/{key} already holds other bytes (sha256 {stored}, {stored_size} B) "
                                f"than this production made (sha256 {sha}, {size} B)")
    return {"key": key, "version_id": version, "sha256": sha, "size": size}


def _put_parts(s3, bucket: str, key: str, path: Path, meta: dict):
    """A multipart upload with a CRC64NVME full-object checksum: ``(sha256, answer)``, the answer
    ``None`` when the key already exists (412)."""
    try:
        from awscrt import checksums
    except ImportError as exc:
        raise StoreError(f"{path.name} is larger than {SINGLE_PUT_MAX} B, and an upload in parts carries a "
                         "CRC64NVME checksum, which needs the awscrt package in this environment") from exc
    from botocore.exceptions import ClientError

    def crc(value: int) -> str:
        return _b64(value.to_bytes(8, "big"))

    sha, whole = hashlib.sha256(), 0
    with open(path, "rb") as fh:
        while chunk := fh.read(_CHUNK):
            sha.update(chunk)
            whole = checksums.crc64nvme(chunk, whole)
    metadata = {**meta, c.META_SHA256: sha.hexdigest()}
    upload = s3.create_multipart_upload(Bucket=bucket, Key=key, Metadata=metadata, ChecksumAlgorithm="CRC64NVME",
                                        ChecksumType="FULL_OBJECT")["UploadId"]
    try:
        parts = []
        with open(path, "rb") as fh:
            number = 1
            while chunk := fh.read(PART_SIZE):
                value = crc(checksums.crc64nvme(chunk))
                answer = s3.upload_part(Bucket=bucket, Key=key, UploadId=upload, PartNumber=number, Body=chunk,
                                        ChecksumAlgorithm="CRC64NVME", ChecksumCRC64NVME=value)
                parts.append({"ETag": answer["ETag"], "PartNumber": number, "ChecksumCRC64NVME": value})
                number += 1
        try:
            done = s3.complete_multipart_upload(Bucket=bucket, Key=key, UploadId=upload, IfNoneMatch="*",
                                                MultipartUpload={"Parts": parts}, ChecksumCRC64NVME=crc(whole),
                                                ChecksumType="FULL_OBJECT")
        except ClientError as exc:
            if _code(exc) not in ("PreconditionFailed", "412"):
                raise
            s3.abort_multipart_upload(Bucket=bucket, Key=key, UploadId=upload)
            return sha.hexdigest(), None
        return sha.hexdigest(), done
    except BaseException:
        try:
            s3.abort_multipart_upload(Bucket=bucket, Key=key, UploadId=upload)
        except Exception:  # noqa: BLE001 - the upload's own failure is the one raised
            pass
        raise


def _outputs(work: Path, made) -> list[tuple[str, Path]]:
    """``(relative posix path, file)`` of what ``make`` declared: the listed files, or every regular
    file below ``work`` (symbolic links left out) when it returned ``None``."""
    root = work.resolve()
    if made is None:
        files = sorted(p for p in work.rglob("*") if p.is_file() and not p.is_symlink())
    else:
        files = [Path(p) if Path(p).is_absolute() else work / p for p in made]
    out = []
    for f in files:
        rel = f.resolve().relative_to(root) if f.resolve().is_relative_to(root) else None
        if rel is None or not f.is_file():
            raise StoreError(f"{f} is not a file in the production's folder {work}")
        if rel.name == c.RECIPE_FILE:
            raise StoreError(f"{c.RECIPE_FILE} names the product's own record; name the output {rel} otherwise")
        out.append((rel.as_posix(), f))
    if not out:
        raise StoreError("the production made no file to store")
    if len({r for r, _ in out}) != len(out):
        raise StoreError("the production names one file twice")
    return out


def _is_marker(key: str) -> bool:
    """Whether ``key`` is a product folder's own :data:`RECIPE_FILE` (``<kind>/<subject>/<hash>/recipe.json``)."""
    parts = PurePosixPath(key).parts
    return len(parts) == 4 and parts[-1] == c.RECIPE_FILE


def _reused(claim: dict, bucket: str) -> dict:
    """The answer for a recipe the catalog already holds; refused when it lies in another bucket."""
    objects = claim["objects"]
    markers = [o for o in objects if _is_marker(o["key"])]
    if len(markers) != 1:
        raise StoreError(f"recipe {claim['recipe_hash']} is recorded with {len(markers)} {c.RECIPE_FILE} files")
    if markers[0]["bucket"] != bucket:
        raise StoreError(f"recipe {claim['recipe_hash']} is stored in bucket {markers[0]['bucket']!r}, not {bucket!r}: "
                         "name the bucket in the recipe's settings to make it another product")
    return {"status": "reused", "recipe_hash": claim["recipe_hash"],
            "record": {k: markers[0][k] for k in c.RECORD_KEYS}, "objects": objects}


def _committed_as(recipe_hash: str, marker: dict) -> bool:
    """Whether the catalog holds ``recipe_hash`` committed (or verified) with ``marker`` among its objects,
    read on a connection of its own."""
    try:
        with connect() as conn:
            row = conn.execute(
                "SELECT p.state FROM catalog.product p JOIN catalog.object o USING (recipe_hash) "
                "WHERE p.recipe_hash = %s AND o.key = %s AND o.version_id = %s",
                (recipe_hash, marker["key"], marker["version_id"])).fetchone()
    except Exception:  # noqa: BLE001 - unknown is not committed; the production's failure is raised
        return False
    return row is not None and row[0] in ("committed", "verified")


def _answer(row) -> dict:
    value = row[0]
    return value if isinstance(value, dict) else json.loads(value)


def produce(*, bucket: str, kind: str, subject: str, tool: str, algorithm_version: str,
            settings: dict, inputs, make, lease_s: int = DEFAULT_LEASE_S) -> dict:
    """Make a product once per recipe and store it.

    The recipe is ``{kind, subject, tool, algorithm_version, settings, inputs}``: ``tool`` the manifest
    key of the tool making it, ``settings`` what changes its bytes (normalised, :func:`normalise`),
    ``inputs`` the pinned records it reads, in the order it reads them. The catalog's claim decides:

    * ``claimed`` — ``make(work)`` runs with an empty temporary folder ``work`` in the working folder and
      returns the files to store (paths in ``work``), or ``None`` for every regular file there. Each is
      stored create-only at ``<kind>/<subject>/<short hash>[.<generation>]/<its path in work>`` in
      ``bucket``, then :data:`RECIPE_FILE`, then the catalog commit. Everything else ``make`` leaves in
      ``work`` is removed with it, in the store :func:`store_name` names. Any failure gives the claim up
      with its reason (``catalog.fail``, an ordinary failure, so a retry claims a new generation: other
      bytes already at one of these keys included, an earlier attempt of the same generation) and is
      raised. Nondeterminism (one recipe, other bytes across productions) is the reconciler's verdict;
    * ``reuse`` — the stored product, nothing made;
    * ``in_progress`` — :class:`ProductBusy` (another unit is making it);
    * ``failed`` / ``nondeterministic`` — :class:`ProductRefused` with the catalog's reason.

    Returns ``{"status": "stored" | "reused", "recipe_hash", "record", "objects"}``: ``record`` the
    product's :data:`RECIPE_FILE` record (what a reader passes on), ``objects`` every stored object
    ``{store, bucket, key, version_id, sha256, size}``.
    """
    commit, fingerprint = _provenance()
    store = store_name()
    recipe = {"kind": kind, "subject": subject, "tool": tool, "algorithm_version": algorithm_version,
              "settings": normalise(dict(settings)), "inputs": [check_record(r) for r in inputs]}
    for name in ("kind", "subject", "tool", "algorithm_version"):
        value = recipe[name]
        if not isinstance(value, str) or not value.strip() or value != value.strip() or (
                name in ("kind", "subject") and ("/" in value or value in (".", ".."))):
            raise ValueError(f"recipe {name} is non-empty text without surrounding spaces"
                             f"{', /, or a name of . or ..' if name in ('kind', 'subject') else ''}; got {value!r}")
    if not _BUCKET.match(bucket) or ".." in bucket:
        raise ValueError(f"bucket is an S3 bucket name, not {bucket!r}")
    s3 = client()
    for record in recipe["inputs"]:
        resolve(record, s3=s3)  # each input names a version the store holds
    with connect() as conn:
        claim = _answer(conn.execute("SELECT catalog.claim(%s::jsonb, make_interval(secs => %s))",
                                     (json.dumps(recipe), int(lease_s))).fetchone())
        status, recipe_hash = claim["status"], claim["recipe_hash"]
        if status == "reuse":
            return _reused(claim, bucket)
        if status == "in_progress":
            raise ProductBusy(f"another unit ({claim['owner']}) is making this product, its claim held until "
                              f"{claim['lease_until']}; ask again after that (recipe {recipe_hash})")
        if status != "claimed":
            raise ProductRefused(f"the catalog holds this recipe as {status}: {claim.get('reason')} "
                                 f"(recipe {recipe_hash}); a human decides")
        objects = None
        try:
            generation = int(claim["generation"])
            folder = f"{kind}/{subject}/{_short(conn, recipe_hash, kind, subject)}"
            folder += f".{generation}" if generation > 1 else ""
            meta = {c.META_RECIPE_HASH: recipe_hash, c.META_PRODUCER: claim["owner"],
                    c.META_TOOLKIT_COMMIT: commit, c.META_ENV_FINGERPRINT: fingerprint}
            with tempfile.TemporaryDirectory(dir=os.getcwd(), prefix="produce-") as work:
                files = _outputs(Path(work), make(Path(work)))
                stored = [_put(s3, bucket, f"{folder}/{rel}", path, meta) for rel, path in files]
                body = json.dumps(c.recipe_file(recipe_hash, generation, recipe, stored), indent=1,
                                  sort_keys=True).encode("utf-8")
                marker = Path(work) / c.RECIPE_FILE
                marker.write_bytes(body)
                stored.append(_put(s3, bucket, f"{folder}/{c.RECIPE_FILE}", marker, meta))
            objects = [{"store": store, "bucket": bucket, **o} for o in stored]
            conn.execute("SELECT catalog.commit(%s, %s::jsonb, %s, %s)",
                         (recipe_hash, json.dumps(objects), commit, fingerprint))
        except BaseException as exc:
            if objects is not None and _committed_as(recipe_hash, objects[-1]):
                # the commit's answer was lost, or the reconciler adopted this very production first
                return {"status": "stored", "recipe_hash": recipe_hash,
                        "record": {k: objects[-1][k] for k in c.RECORD_KEYS}, "objects": objects}
            try:
                conn.execute("SELECT catalog.fail(%s, %s)",
                             (recipe_hash, f"{type(exc).__name__}: {exc}"[:2000] or type(exc).__name__))
            except Exception:  # noqa: BLE001 - the production's own failure is the one raised
                pass
            raise
    record = {k: objects[-1][k] for k in c.RECORD_KEYS}
    return {"status": "stored", "recipe_hash": recipe_hash, "record": record, "objects": objects}
