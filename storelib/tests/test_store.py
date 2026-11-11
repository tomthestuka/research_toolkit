"""storelib against a throwaway store (the root conftest's ``store`` fixture): what a production stores,
what a reader stages, and every refusal on the way."""

from __future__ import annotations

import json
import os
import stat

import pytest

import storelib
from storelib import contract, store as st

BUCKET = "ixpe"


# ------------------------------------------------------------------ values (no store)


def test_a_whole_float_and_an_integer_name_one_setting():
    assert storelib.normalise({"a": 2.0, "b": (1, 2.5), "c": [True, None, "x"]}) == {
        "a": 2, "b": [1, 2.5], "c": [True, None, "x"]}
    assert json.dumps(storelib.normalise({"x": 2.0})) == json.dumps({"x": 2})


@pytest.mark.parametrize("value", [float("nan"), float("inf"), {1: "a"}, object()])
def test_a_value_a_recipe_cannot_hash_is_refused(value):
    with pytest.raises(ValueError):
        storelib.normalise(value)


@pytest.mark.parametrize("record", [
    {"store": "db", "bucket": "ixpe", "key": "k"},
    {"store": "db", "bucket": "ixpe", "key": "k", "version_id": " "},
    {"store": "db", "bucket": "ixpe", "key": "k", "version_id": "v", "extra": "x"},
    "db/ixpe/k",
])
def test_a_read_names_a_whole_pinned_record(record):
    with pytest.raises(ValueError, match="a record is"):
        storelib.check_record(record)


def test_the_recipe_file_shape_is_checked():
    good = {"key": "k/s/h/out", "version_id": "v", "sha256": "0" * 64, "size": 1}
    recipe = dict.fromkeys(contract.RECIPE_KEYS, "x") | {"settings": {}, "inputs": []}
    assert contract.recipe_file("a" * 64, 1, recipe, [good])["generation"] == 1
    for bad in ([], [{**good, "size": -1}], [{**good, "sha256": "X" * 64}], [{**good, "more": 1}]):
        with pytest.raises(contract.ContractError):
            contract.recipe_file("a" * 64, 1, recipe, bad)


def test_a_process_without_store_keys_is_told_which_are_unset(monkeypatch):
    for name in ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_ENDPOINT_URL_S3"):
        monkeypatch.delenv(name, raising=False)
    with pytest.raises(storelib.StoreError, match="AWS_ACCESS_KEY_ID"):
        storelib.client()


# ------------------------------------------------------------------ production


def _writer(files: dict[str, bytes], log: list | None = None):
    def make(work):
        if log is not None:
            log.append(work)
        for rel, data in files.items():
            (work / rel).parent.mkdir(parents=True, exist_ok=True)
            (work / rel).write_bytes(data)
        (work / "scratch.tmp.link").symlink_to(work / next(iter(files)))  # a link is never stored
        return None
    return make


def _produce(store, make, *, subject=None, settings=None, inputs=()):
    return storelib.produce(bucket=BUCKET, kind="test", subject=subject or f"s{store['tag']}",
                            tool="research_toolkit.storelib_test.make", algorithm_version="1",
                            settings=settings if settings is not None else {"n": 2.0}, inputs=list(inputs),
                            make=make)


def _product_row(recipe_hash):
    with storelib.connect() as conn:
        return conn.execute("SELECT state, generation, reason FROM catalog.product WHERE recipe_hash = %s",
                            (recipe_hash,)).fetchone()


def test_a_production_stores_each_file_then_its_recipe_and_commits(store, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    made = []
    out = _produce(store, _writer({"spec.pha": b"spectrum", "sub/resp.rmf": b"response"}, made))
    assert out["status"] == "stored"
    keys = [o["key"] for o in out["objects"]]
    folder = f"test/s{store['tag']}/{out['recipe_hash'][:8]}"
    assert keys == [f"{folder}/spec.pha", f"{folder}/sub/resp.rmf", f"{folder}/recipe.json"]
    assert out["record"] == {k: out["objects"][-1][k] for k in contract.RECORD_KEYS}
    assert not made[0].exists(), "the production's folder outlives it"
    assert _product_row(out["recipe_hash"])[0] == "committed"

    s3 = storelib.client()
    for obj in out["objects"]:
        head = s3.head_object(Bucket=BUCKET, Key=obj["key"], VersionId=obj["version_id"])
        meta = head["Metadata"]
        assert sorted(meta) == sorted(contract.META_NAMES)
        assert meta[contract.META_SHA256] == obj["sha256"] and meta[contract.META_RECIPE_HASH] == out["recipe_hash"]
        assert meta[contract.META_PRODUCER] == os.environ["PGUSER"]
    marker = json.loads(s3.get_object(Bucket=BUCKET, Key=keys[-1], VersionId=out["objects"][-1]["version_id"])
                        ["Body"].read())
    contract.check_recipe_file(marker)
    assert marker["objects"] == [{k: o[k] for k in contract.RECIPE_OBJECT_KEYS} for o in out["objects"][:-1]]
    assert marker["recipe"]["settings"] == {"n": 2} and marker["generation"] == 1


def test_the_same_recipe_is_reused_and_made_once(store, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    first = _produce(store, _writer({"a.txt": b"one"}), settings={"n": 2})
    made = []
    again = _produce(store, _writer({"a.txt": b"other"}, made), settings={"n": 2.0})
    assert again["status"] == "reused" and not made
    assert again["record"] == first["record"] and again["recipe_hash"] == first["recipe_hash"]


def test_a_reader_stages_a_whole_product_and_one_object_read_only(store, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    out = _produce(store, _writer({"spec.pha": b"spectrum", "sub/resp.rmf": b"response"}))
    [folder, single] = storelib.stage_in([out["record"], {k: out["objects"][0][k] for k in contract.RECORD_KEYS}],
                                         tmp_path / "in")
    assert (folder / "spec.pha").read_bytes() == b"spectrum" and (folder / "sub" / "resp.rmf").read_bytes() == b"response"
    assert single.read_bytes() == b"spectrum" and single == folder / "spec.pha"
    assert not os.stat(folder / "spec.pha").st_mode & stat.S_IWUSR
    assert storelib.resolve({k: out["objects"][0][k] for k in contract.RECORD_KEYS})["sha256"] == out["objects"][0]["sha256"]
    # A second staging keeps the files that hold the listed bytes.
    again = storelib.stage_in(out["record"], tmp_path / "in")
    assert again == [folder]


def test_a_read_whose_bytes_differ_from_the_record_is_refused(store, tmp_path):
    s3 = storelib.client()
    key = f"test/s{store['tag']}/forged/x.dat"
    version = s3.put_object(Bucket=BUCKET, Key=key, Body=b"real", IfNoneMatch="*",
                            Metadata={contract.META_SHA256: "0" * 64})["VersionId"]
    with pytest.raises(storelib.StoreError, match="differ from the record"):
        storelib.stage_in({"store": store["store"], "bucket": BUCKET, "key": key, "version_id": version}, tmp_path)
    assert not (tmp_path / BUCKET / key).exists()


def test_a_failed_production_gives_its_claim_up_and_a_retry_writes_a_new_generation(store, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    def broken(work):
        (work / "half.dat").write_bytes(b"x")
        raise RuntimeError("the reduction stopped")

    with pytest.raises(RuntimeError, match="the reduction stopped"):
        _produce(store, broken)
    out = _produce(store, _writer({"a.txt": b"made"}))
    state, generation, _ = _product_row(out["recipe_hash"])
    assert (state, generation) == ("committed", 2)
    assert out["objects"][0]["key"].split("/")[2] == f"{out['recipe_hash'][:8]}.2"


def test_a_production_that_names_no_file_is_refused_and_failed(store, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(storelib.StoreError, match="made no file"):
        _produce(store, lambda work: [])
    with storelib.connect() as conn:
        row = conn.execute("SELECT state, reason FROM catalog.product WHERE subject = %s",
                           (f"s{store['tag']}",)).fetchone()
    assert row[0] == "failed" and "made no file" in row[1]


def test_other_bytes_at_a_key_of_this_production_fail_it_and_a_retry_writes_a_new_generation(
        store, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    subject = f"s{store['tag']}"
    recipe = {"kind": "test", "subject": subject, "tool": "research_toolkit.storelib_test.make",
              "algorithm_version": "1", "settings": {"n": 2}, "inputs": []}
    with storelib.connect() as conn:  # an earlier attempt of this production (its lease still live)
        claim = conn.execute("SELECT catalog.claim(%s::jsonb)", (json.dumps(recipe),)).fetchone()[0]
        folder = f"test/{subject}/{st._short(conn, claim['recipe_hash'], 'test', subject)}"
    storelib.client().put_object(Bucket=BUCKET, Key=f"{folder}/a.txt", Body=b"earlier", IfNoneMatch="*",
                                 Metadata={contract.META_SHA256: "1" * 64})
    with pytest.raises(storelib.StoreError, match="already holds other bytes"):
        _produce(store, _writer({"a.txt": b"now"}))
    state, _, reason = _product_row(claim["recipe_hash"])
    assert state == "failed" and "already holds other bytes" in reason  # a concurrency fault, not a verdict
    out = _produce(store, _writer({"a.txt": b"now"}))
    assert out["status"] == "stored" and out["objects"][0]["key"] == f"{folder}.2/a.txt"
    assert _product_row(claim["recipe_hash"])[:2] == ("committed", 2)


def test_a_process_whose_store_is_unnamed_claims_nothing(store, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("TOOLKIT_STORE")
    with pytest.raises(storelib.StoreError, match="TOOLKIT_STORE unset"):
        _produce(store, _writer({"a.txt": b"x"}))
    with storelib.connect() as conn:
        assert conn.execute("SELECT count(*) FROM catalog.product WHERE subject = %s",
                            (f"s{store['tag']}",)).fetchone()[0] == 0


def test_the_same_bytes_at_a_key_of_this_production_are_kept(store, tmp_path, monkeypatch):
    import hashlib

    monkeypatch.chdir(tmp_path)
    subject = f"s{store['tag']}"
    recipe = {"kind": "test", "subject": subject, "tool": "research_toolkit.storelib_test.make",
              "algorithm_version": "1", "settings": {"n": 2}, "inputs": []}
    with storelib.connect() as conn:
        claim = conn.execute("SELECT catalog.claim(%s::jsonb)", (json.dumps(recipe),)).fetchone()[0]
        folder = f"test/{subject}/{st._short(conn, claim['recipe_hash'], 'test', subject)}"
    earlier = storelib.client().put_object(
        Bucket=BUCKET, Key=f"{folder}/a.txt", Body=b"same", IfNoneMatch="*",
        Metadata={contract.META_SHA256: hashlib.sha256(b"same").hexdigest()})["VersionId"]
    out = _produce(store, _writer({"a.txt": b"same"}))
    assert out["status"] == "stored" and out["objects"][0]["version_id"] == earlier


def test_a_product_records_its_provenance_or_is_never_claimed(store, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("AR_ENV_FINGERPRINT")
    with pytest.raises(storelib.StoreError, match="AR_ENV_FINGERPRINT unset"):
        _produce(store, _writer({"a.txt": b"x"}))
    with storelib.connect() as conn:
        assert conn.execute("SELECT count(*) FROM catalog.product WHERE subject = %s",
                            (f"s{store['tag']}",)).fetchone()[0] == 0


def test_a_product_records_the_inputs_it_read(store, seed, tmp_path, monkeypatch):
    import base64

    monkeypatch.chdir(tmp_path)
    source = seed(kind="raw", subject=f"r{store['tag']}", bucket=BUCKET,
                  files={"raw.evt": base64.b64encode(b"events").decode("ascii")})
    derived = _produce(store, _writer({"spec.pha": b"spectrum"}), inputs=[source["record"]])
    with storelib.connect() as conn:
        recipe = conn.execute("SELECT recipe FROM catalog.product WHERE recipe_hash = %s",
                              (derived["recipe_hash"],)).fetchone()[0]
    assert recipe["inputs"] == [source["record"]]
    with pytest.raises(ValueError, match="a record is"):
        _produce(store, _writer({"x": b"x"}), inputs=[{"store": "db", "bucket": BUCKET, "key": "k"}])


def test_an_object_above_the_single_put_limit_goes_up_in_parts(store, tmp_path, monkeypatch):
    pytest.importorskip("awscrt", reason="the multipart path needs awscrt")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(st, "SINGLE_PUT_MAX", 1024)
    monkeypatch.setattr(st, "PART_SIZE", 5 * 1024 ** 2)
    data = os.urandom(5 * 1024 ** 2 + 4096)
    out = _produce(store, _writer({"big.dat": data}))
    head = storelib.client().head_object(Bucket=BUCKET, Key=out["objects"][0]["key"],
                                         VersionId=out["objects"][0]["version_id"], ChecksumMode="ENABLED")
    assert head["ContentLength"] == len(data) and head.get("ChecksumCRC64NVME")
    [folder] = storelib.stage_in(out["record"], tmp_path / "in")
    assert (folder / "big.dat").read_bytes() == data


def test_a_remote_read_is_pinned_to_its_version(store, tmp_path, monkeypatch):
    pytest.importorskip("s3fs", reason="remote reads need s3fs")
    monkeypatch.chdir(tmp_path)
    out = _produce(store, _writer({"a.txt": b"remote bytes"}))
    record = {k: out["objects"][0][k] for k in contract.RECORD_KEYS}
    url, options = storelib.remote(record)
    assert url.endswith(f"?versionId={record['version_id']}") and options["version_aware"] is True
    with storelib.open_remote(record) as fh:
        assert fh.read() == b"remote bytes"


def test_a_bucket_outside_s3s_naming_is_refused_before_any_read():
    for bucket in ("..", "/etc", "Ixpe", "a"):
        with pytest.raises(ValueError, match="bucket"):
            storelib.check_record({"store": "db", "bucket": bucket, "key": "k", "version_id": "v"})


def test_a_production_names_no_file_recipe_json_and_no_file_outside_its_folder(store, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(storelib.StoreError, match="names the product's own record"):
        _produce(store, _writer({"sub/recipe.json": b"{}"}))
    outside = tmp_path / "outside.txt"
    outside.write_bytes(b"x")
    with pytest.raises(storelib.StoreError, match="is not a file in the production's folder"):
        _produce(store, lambda work: [outside], subject=f"o{store['tag']}")


def test_an_input_naming_a_version_the_store_lacks_is_refused_before_the_claim(store, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    missing = {"store": store["store"], "bucket": BUCKET, "key": "test/none/x", "version_id": "no-such-version"}
    with pytest.raises(storelib.StoreError, match="no-such-version"):
        _produce(store, _writer({"a.txt": b"x"}), inputs=[missing])
    with storelib.connect() as conn:
        assert conn.execute("SELECT count(*) FROM catalog.product WHERE subject = %s",
                            (f"s{store['tag']}",)).fetchone()[0] == 0


def test_a_recipe_stored_in_another_bucket_is_refused_as_this_buckets(store, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _produce(store, _writer({"a.txt": b"x"}))
    with pytest.raises(storelib.StoreError, match="is stored in bucket 'ixpe', not 'nicer'"):
        storelib.produce(bucket="nicer", kind="test", subject=f"s{store['tag']}",
                         tool="research_toolkit.storelib_test.make", algorithm_version="1", settings={"n": 2},
                         inputs=[], make=_writer({"a.txt": b"x"}))
