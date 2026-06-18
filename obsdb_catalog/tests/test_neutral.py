import json
import os
import threading

import pytest

from obsdb_catalog import neutral


def touch(obs_db, rel, content=b"data"):
    full = os.path.join(obs_db, rel)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, "wb") as fh:
        fh.write(content)
    return rel


def prod(kind, rel, produced_by="tool", params=None):
    return {"kind": kind, "path": rel, "produced_by": produced_by, "params": params or {}}


def test_read_absent_returns_empty(tmp_path):
    cat = neutral.read(tmp_path)
    assert cat == {"mission": os.path.basename(tmp_path).lower(), "entries": {}}


def test_round_trip(tmp_path):
    p = touch(tmp_path, "obs1/products/a.pha")
    neutral.advance(tmp_path, "nustar", "obs1", level=3,
                    products=[prod("pha_src_fpmA", p, params={"fpm": "A"})])
    cat = neutral.read(tmp_path)
    entry = cat["entries"]["obs1"]
    assert entry["max_level"] == 3
    assert entry["levels_present"] == [1, 2, 3]
    assert entry["products"][0]["kind"] == "pha_src_fpmA"
    assert entry["products"][0]["params"] == {"fpm": "A"}
    assert "created_at" in entry["products"][0]


def test_idempotent_upsert(tmp_path):
    p = touch(tmp_path, "obs1/x.pha")
    neutral.advance(tmp_path, "m", "obs1", level=3,
                    products=[prod("pha", p, params={"v": 1})])
    neutral.advance(tmp_path, "m", "obs1", level=3,
                    products=[prod("pha", p, params={"v": 2})])
    products = neutral.read(tmp_path)["entries"]["obs1"]["products"]
    assert len(products) == 1               # not duplicated
    assert products[0]["params"]["v"] == 2  # updated in place


def test_truthful_absent_file_raises_and_no_write(tmp_path):
    with pytest.raises(neutral.CatalogError):
        neutral.advance(tmp_path, "m", "obs1", level=3,
                        products=[prod("pha", "obs1/missing.pha")])
    # nothing recorded
    assert neutral.read(tmp_path)["entries"] == {}
    assert not os.path.exists(os.path.join(tmp_path, neutral.CATALOG_NAME))


def test_relative_path_required(tmp_path):
    absolute = touch(tmp_path, "obs1/a.pha")
    with pytest.raises(neutral.CatalogError):
        neutral.advance(tmp_path, "m", "obs1", level=3,
                        products=[prod("pha", os.path.join(tmp_path, absolute))])


def test_regression_keeps_products_lowers_max_level(tmp_path):
    evt = touch(tmp_path, "obs1/event_l2/evt.fits")
    pcube = touch(tmp_path, "obs1/event_l2/pcube.fits")
    neutral.advance(tmp_path, "ixpe", "obs1", level=2,
                    products=[prod("level2_evt", evt)])
    neutral.advance(tmp_path, "ixpe", "obs1", level=3,
                    products=[prod("pcube", pcube)])
    # upstream re-run regresses to level 2
    entry = neutral.advance(tmp_path, "ixpe", "obs1", level=2,
                            products=[prod("level2_evt", evt)])
    assert entry["max_level"] == 2
    assert entry["levels_present"] == [1, 2]
    kinds = {p["kind"] for p in entry["products"]}
    assert kinds == {"level2_evt", "pcube"}  # stale pcube left in place
    # re-reduction restores level 3
    entry = neutral.advance(tmp_path, "ixpe", "obs1", level=3,
                            products=[prod("pcube", pcube)])
    assert entry["max_level"] == 3


def test_new_is_level_two(tmp_path):
    evt = touch(tmp_path, "obs1/event_l2/evt.fits")
    entry = neutral.new(tmp_path, "ixpe", "obs1", products=[prod("level2_evt", evt)])
    assert entry["max_level"] == 2
    assert entry["levels_present"] == [1, 2]


def test_history_snapshots_and_retention(tmp_path):
    p = touch(tmp_path, "obs1/a.pha")
    for _ in range(neutral.HISTORY_KEEP + 5):
        neutral.advance(tmp_path, "m", "obs1", level=3, products=[prod("pha", p)])
    hist = os.path.join(tmp_path, neutral.HISTORY_DIR)
    snaps = [f for f in os.listdir(hist) if f.startswith("catalog-")]
    assert len(snaps) == neutral.HISTORY_KEEP  # pruned to last N
    for snap in snaps:  # every snapshot parses
        with open(os.path.join(hist, snap)) as fh:
            json.load(fh)


def test_atomicity_failed_write_leaves_original(tmp_path, monkeypatch):
    p = touch(tmp_path, "obs1/a.pha")
    neutral.advance(tmp_path, "m", "obs1", level=3, products=[prod("pha", p)])
    before = neutral.read(tmp_path)

    def boom(*a, **k):
        raise RuntimeError("disk full")

    monkeypatch.setattr(neutral.json, "dump", boom)
    with pytest.raises(RuntimeError):
        neutral.advance(tmp_path, "m", "obs2", level=3,
                        products=[prod("pha", touch(tmp_path, "obs2/a.pha"))])
    # original intact, parseable, no temp files leaked
    assert neutral.read(tmp_path) == before
    leaked = [f for f in os.listdir(tmp_path) if f.startswith(".catalog-")]
    assert leaked == []


def test_read_malformed_raises(tmp_path):
    with open(os.path.join(tmp_path, neutral.CATALOG_NAME), "w") as fh:
        fh.write("{ this is not json")
    with pytest.raises(neutral.CatalogError):
        neutral.read(tmp_path)


def test_concurrency_no_lost_entries(tmp_path):
    n = 40
    for i in range(n):
        touch(tmp_path, f"obs{i}/a.pha")

    def worker(lo, hi):
        for i in range(lo, hi):
            neutral.advance(tmp_path, "m", f"obs{i}", level=3,
                            products=[prod("pha", f"obs{i}/a.pha")])

    threads = [
        threading.Thread(target=worker, args=(0, n // 2)),
        threading.Thread(target=worker, args=(n // 2, n)),
    ]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    entries = neutral.read(tmp_path)["entries"]
    assert set(entries) == {f"obs{i}" for i in range(n)}


def test_sha256_optional(tmp_path):
    p = touch(tmp_path, "obs1/a.pha")
    neutral.advance(tmp_path, "m", "obs1", level=3, products=[prod("pha", p)])
    assert "sha256" not in neutral.read(tmp_path)["entries"]["obs1"]["products"][0]
    neutral.advance(tmp_path, "m", "obs1", level=3,
                    products=[prod("pha", p)], compute_sha256=True)
    assert "sha256" in neutral.read(tmp_path)["entries"]["obs1"]["products"][0]
