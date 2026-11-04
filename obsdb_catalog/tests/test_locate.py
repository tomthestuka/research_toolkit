"""Where an ObsID's folder is: HEASARC's archive tree, the catalog's record, the flat layout."""

import pytest

from obsdb_catalog import locate, mission_dir, neutral
from obsdb_catalog.__main__ import main
from obsdb_catalog.healthcheck import healthcheck
from obsdb_catalog.repair import repair
from obsdb_catalog.tests.test_neutral import prod, touch

# HEASARC's tree below each mission folder, as the archive's access_url gives it (2026-11-04).
CASES = {
    "ixpe": ("01002901", "data/obs/01/01002901", "event_l2/ixpe01002901_det1_evt2_v01.fits"),
    "nicer": ("4690010101", "data/obs/2022_02/4690010101", "xti/event_cl/ni4690010101_0mpu7_cl.evt"),
    "nustar": ("30702017002", "data/obs/07/3/30702017002", "event_cl/nu30702017002A01_cl.evt"),
}


@pytest.mark.parametrize("mission", sorted(CASES))
def test_a_heasarc_tree_is_found_scanned_and_healthy(tmp_path, mission):
    obsid, tree, level2 = CASES[mission]
    folder = tmp_path / mission
    touch(str(folder), f"{tree}/{level2}")

    assert locate(tmp_path, mission, obsid) == folder / tree  # by the tree rule, no catalog yet
    catalog = repair(str(folder), log=None)
    paths = [p["path"] for p in catalog["entries"][obsid]["products"]]
    assert paths == [f"{tree}/{level2}"]  # relative to the mission folder, through the tree
    assert healthcheck(str(folder)) is True
    assert locate(tmp_path, mission, obsid) == folder / tree  # now by the catalog's record


def test_the_catalog_record_comes_before_the_tree_rule(tmp_path):
    folder = tmp_path / "ixpe"
    touch(str(folder), "data/obs/01/01002901/event_l2/x.fits")
    flat = touch(str(folder), "01002901/event_l2/ixpe01002901_det1_evt2_v01.fits")
    neutral.new(str(folder), "ixpe", "01002901", products=[prod("level2_evt", flat)])
    assert locate(tmp_path, "ixpe", "01002901") == folder / "01002901"


def test_the_mission_folder_matches_case_insensitively(tmp_path):
    touch(str(tmp_path / "IXPE"), "data/obs/02/02008201/event_l2/x.fits")
    assert mission_dir(tmp_path, "ixpe") == tmp_path / "IXPE"
    assert locate(tmp_path, "ixpe", "02008201") == tmp_path / "IXPE" / "data/obs/02/02008201"
    assert mission_dir(tmp_path, "NuSTAR") == tmp_path / "nustar"  # absent: HEASARC's spelling


def test_an_absent_or_doubled_obsid_fails_naming_where_it_looked(tmp_path):
    touch(str(tmp_path / "nicer"), "data/obs/2022_02/4690010101/xti/a")
    touch(str(tmp_path / "nicer"), "data/obs/2022_03/4690010101/xti/a")
    with pytest.raises(LookupError, match="more than one"):
        locate(tmp_path, "nicer", "4690010101")
    with pytest.raises(FileNotFoundError, match="data/obs/\\*/\\{obsid\\}"):
        locate(tmp_path, "nicer", "4690010199")


def test_the_cli_prints_the_folder(tmp_path, capsys):
    touch(str(tmp_path / "nustar"), "data/obs/07/3/30702017002/event_cl/a")
    assert main(["locate", str(tmp_path), "nustar", "30702017002"]) == 0
    assert capsys.readouterr().out.strip() == str(tmp_path / "nustar" / "data/obs/07/3/30702017002")


def test_the_newest_level2_record_wins_and_a_path_leaving_the_folder_is_passed_over(tmp_path):
    folder = tmp_path / "ixpe"
    flat = touch(str(folder), "01002901/event_l2/ixpe01002901_det1_evt2_v01.fits")
    tree = touch(str(folder), "data/obs/01/01002901/event_l2/ixpe01002901_det1_evt2_v01.fits")
    neutral.new(str(folder), "ixpe", "01002901",
                products=[dict(prod("level2_evt", flat), created_at="2026-01-01T00:00:00+00:00")])
    neutral.new(str(folder), "ixpe", "01002901",
                products=[dict(prod("level2_evt", tree), created_at="2026-11-04T00:00:00+00:00")])
    assert locate(tmp_path, "ixpe", "01002901") == folder / "data/obs/01/01002901"

    elsewhere = tmp_path / "other" / "02008201"
    elsewhere.mkdir(parents=True)
    touch(str(folder), "data/obs/02/02008201/event_l2/a.fits")
    entry = {"obsid": "02008201", "max_level": 3, "levels_present": [1, 2, 3], "products": [
        {"kind": "pha", "path": "../other/02008201/a.csv", "produced_by": "x", "params": {}}]}
    catalog = neutral.read(str(folder))
    catalog["entries"]["02008201"] = entry
    with neutral._lock(str(folder)):
        neutral._commit(str(folder), catalog)
    assert locate(tmp_path, "ixpe", "02008201") == folder / "data/obs/02/02008201"


def test_the_cli_exits_1_naming_where_it_looked(tmp_path, capsys):
    assert main(["locate", str(tmp_path), "ixpe", "01002901"]) == 1
    assert "no folder for ixpe ObsID 01002901" in capsys.readouterr().err
