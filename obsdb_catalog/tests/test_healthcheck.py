import json
import os

from obsdb_catalog import neutral
from obsdb_catalog.healthcheck import healthcheck
from obsdb_catalog.tests.test_neutral import prod, touch


def test_healthy_returns_true(tmp_path):
    p = touch(tmp_path, "obs1/a.pha")
    neutral.advance(tmp_path, "m", "obs1", level=3, products=[prod("pha", p)])
    assert healthcheck(tmp_path) is True


def test_missing_file_reported(tmp_path):
    p = touch(tmp_path, "obs1/a.pha")
    neutral.advance(tmp_path, "m", "obs1", level=3, products=[prod("pha", p)])
    os.unlink(os.path.join(tmp_path, p))
    report = healthcheck(tmp_path)
    assert report is not True
    assert report["missing"][0]["obsid"] == "obs1"


def test_malformed_reported(tmp_path):
    with open(os.path.join(tmp_path, neutral.CATALOG_NAME), "w") as fh:
        fh.write("not json")
    report = healthcheck(tmp_path)
    assert "malformed" in report


def test_regressed_reported(tmp_path):
    evt = touch(tmp_path, "obs1/event_l2/evt.fits")
    pcube = touch(tmp_path, "obs1/event_l2/pcube.fits")
    neutral.advance(tmp_path, "ixpe", "obs1", level=3, products=[prod("pcube", pcube)])
    neutral.advance(tmp_path, "ixpe", "obs1", level=2,
                    products=[prod("level2_evt", evt)])
    report = healthcheck(tmp_path)
    assert report["regressed"] == ["obs1"]


def test_level_inconsistent_reported(tmp_path):
    catalog = {
        "mission": "m",
        "entries": {
            "obs1": {
                "obsid": "obs1",
                "max_level": 3,
                "levels_present": [1, 3],  # not [1,2,3]
                "products": [],
            }
        },
    }
    with open(os.path.join(tmp_path, neutral.CATALOG_NAME), "w") as fh:
        json.dump(catalog, fh)
    report = healthcheck(tmp_path)
    assert report["level_inconsistent"][0]["obsid"] == "obs1"


def test_uncatalogued_reported(tmp_path):
    nicer = tmp_path / "NICER"
    nicer.mkdir()
    touch(str(nicer), "5000/xti/event_cl/ni5000mpu7_sr.pha")
    report = healthcheck(str(nicer))
    assert any("5000" in p for p in report["uncatalogued"])
