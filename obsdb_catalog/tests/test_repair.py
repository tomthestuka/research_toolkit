import os

from obsdb_catalog import neutral
from obsdb_catalog.healthcheck import healthcheck
from obsdb_catalog.repair import repair
from obsdb_catalog.tests.test_neutral import prod, touch


def test_restore_from_snapshot_when_unparseable(tmp_path):
    p = touch(tmp_path, "obs1/a.pha")
    neutral.advance(tmp_path, "m", "obs1", level=3, products=[prod("pha", p)])
    # a second write snapshots the state that already contains obs1
    p2 = touch(tmp_path, "obs2/a.pha")
    neutral.advance(tmp_path, "m", "obs2", level=3, products=[prod("pha", p2)])

    # corrupt the live catalog
    with open(os.path.join(tmp_path, neutral.CATALOG_NAME), "w") as fh:
        fh.write("{ broken")

    catalog = repair(tmp_path, mission="m")
    assert "obs1" in catalog["entries"]          # recovered from snapshot
    assert healthcheck(tmp_path) is True          # result is consistent on disk


def test_reconcile_drops_entries_with_missing_files(tmp_path):
    p = touch(tmp_path, "obs1/a.pha")
    neutral.advance(tmp_path, "m", "obs1", level=3, products=[prod("pha", p)])
    os.unlink(os.path.join(tmp_path, p))
    catalog = repair(tmp_path, mission="m")
    assert "obs1" not in catalog["entries"]


def test_rederive_nicer_from_disk(tmp_path):
    nicer = tmp_path / "NICER"
    nicer.mkdir()
    obs = "4690010101"
    for suffix in ("mpu7_sr.pha", "mpu7_bg.pha", "mpu7.arf", "mpu7.rmf"):
        touch(str(nicer), f"{obs}/xti/event_cl/ni{obs}{suffix}")

    catalog = repair(str(nicer))  # mission inferred from path -> "nicer"
    entry = catalog["entries"][obs]
    assert entry["max_level"] == 3
    assert entry["levels_present"] == [1, 2, 3]
    by_kind = {p["kind"]: p for p in entry["products"]}
    assert set(by_kind) == {"pha", "pha_bkg", "arf", "rmf"}
    assert by_kind["pha"]["produced_by"] == "nicer_red.sh"
    assert by_kind["pha"]["params"] == {"background_model": "3C50", "mpu": 7}
    assert by_kind["arf"]["params"] == {"mpu": 7}
    assert healthcheck(str(nicer)) is True


def test_rederive_nustar_per_fpm(tmp_path):
    nustar = tmp_path / "NuSTAR"
    nustar.mkdir()
    obs = "30702017002"
    for fpm in ("A", "B"):
        for suffix in ("_sr_grp.pha", "_sr.pha", "_sr.arf", "_sr.rmf"):
            touch(str(nustar), f"{obs}/products/nu{obs}{fpm}01{suffix}")

    catalog = repair(str(nustar))
    by_kind = {p["kind"]: p for p in catalog["entries"][obs]["products"]}
    assert "pha_grp_fpmA" in by_kind and "rmf_fpmB" in by_kind
    assert by_kind["pha_grp_fpmA"]["produced_by"] == "nustar_grp.sh"
    assert by_kind["pha_grp_fpmA"]["params"] == {"fpm": "A", "grouping": "optimal"}
    assert by_kind["pha_src_fpmB"]["params"] == {"fpm": "B"}


def test_rederive_ixpe_filter_and_special_case(tmp_path):
    ixpe = tmp_path / "IXPE"
    ixpe.mkdir()
    # normal obs uses the _rej_src_ infix
    touch(str(ixpe), "01002901/event_l2/ixpe01002901_det1_evt2_v01.fits")
    touch(str(ixpe), "01002901/event_l2/ixpe01002901_det1_evt2_v01_rej_src_pcube.fits")
    touch(str(ixpe), "01002901/event_l2/ixpe01002901_det1_evt2_v01_rej_src_pha1.fits")
    # the 03002599 special case uses the _src_ infix
    touch(str(ixpe), "03002599/event_l2/ixpe03002599_det1_evt2_v01.fits")
    touch(str(ixpe), "03002599/event_l2/ixpe03002599_det1_evt2_v01_src_pcube.fits")

    catalog = repair(str(ixpe))
    normal = {p["kind"]: p for p in catalog["entries"]["01002901"]["products"]}
    assert normal["pcube"]["params"]["filter"] == "rej_src"
    assert normal["pcube"]["params"]["calibration_version"] == "v01"
    assert normal["pha"]["params"]["stokes"] == "I"

    special = {p["kind"]: p for p in catalog["entries"]["03002599"]["products"]}
    assert special["pcube"]["params"]["filter"] == "src"


def test_repair_additive_preserves_existing_params(tmp_path):
    nicer = tmp_path / "NICER"
    nicer.mkdir()
    obs = "4690010101"
    p = touch(str(nicer), f"{obs}/xti/event_cl/ni{obs}mpu7_sr.pha")
    # pre-existing entry with a richer param the scanner would not produce
    neutral.advance(str(nicer), "nicer", obs, level=3,
                    products=[prod("pha", p, produced_by="nicer_red.sh",
                                   params={"background_model": "3C50", "mpu": 7,
                                           "extra": "keep-me"})])
    repair(str(nicer))
    pha = next(p for p in neutral.read(str(nicer))["entries"][obs]["products"]
               if p["kind"] == "pha")
    assert pha["params"].get("extra") == "keep-me"  # not clobbered by re-derive
