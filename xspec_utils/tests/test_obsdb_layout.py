"""The xspec_utils readers find an observation through ``obsdb_catalog.locate``, in an obs-db laid
out as HEASARC's archive is (``<root>/<mission>/data/obs/.../<obsid>/``)."""

from __future__ import annotations

from pathlib import Path

from xspec_utils.lib import NICER, NuSTAR, red_ixpe


def _touch(path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"x")
    return str(path)


def test_nicer_files_are_found_in_the_month_folder(tmp_path):
    obs = "4690010101"
    cl = tmp_path / "nicer" / "data/obs/2022_02" / obs / "xti/event_cl"
    src, bkg = _touch(cl / f"ni{obs}mpu7_sr.pha"), _touch(cl / f"ni{obs}mpu7_bg.pha")
    rmf, arf = _touch(cl / f"ni{obs}mpu7.rmf"), _touch(cl / f"ni{obs}mpu7.arf")
    path = str(tmp_path / "nicer")
    assert NICER.pha(obs, path) == (src, bkg)
    assert (NICER.rmf(obs, path), NICER.arf(obs, path)) == (rmf, arf)


def test_nustar_files_are_found_in_the_archive_tree(tmp_path):
    obs = "30702017002"
    products = tmp_path / "NuSTAR" / "data/obs/07/3" / obs / "products"
    src = [_touch(products / f"nu{obs}{m}01_sr_grp.pha") for m in "AB"]
    bkg = [_touch(products / f"nu{obs}{m}01_bk.pha") for m in "AB"]
    rmf = [_touch(products / f"nu{obs}{m}01_sr.rmf") for m in "AB"]
    path = str(tmp_path / "NuSTAR")
    assert NuSTAR.pha(obs, path) == (src, bkg)
    assert NuSTAR.rmf(obs, path) == rmf


def test_ixpe_level2_files_are_listed_from_the_archive_tree(tmp_path):
    files = [_touch(tmp_path / "ixpe" / f"data/obs/{obs[:2]}" / obs / "event_l2"
                    / f"ixpe{obs}_det{du}_evt2_v01_rej.fits")
             for obs in ("01002901", "02008201") for du in (1, 2, 3)]
    assert red_ixpe.get_all_file_list(str(tmp_path / "ixpe")) == files


def test_pha_refuses_inputs_outside_an_ixpe_mission_folder_before_writing(tmp_path):
    import pytest

    obs = "02008201"
    event_l2 = tmp_path / "elsewhere" / "data/obs/02" / obs / "event_l2"  # no folder named ixpe above
    files = [_touch(event_l2 / f"ixpe{obs}_det{du}_evt2_v01_rej.fits") for du in (1, 2, 3)]
    before = sorted(p.name for p in event_l2.iterdir())
    with pytest.raises(ValueError, match="obsdb_catalog.locate"):
        red_ixpe.pha(obs, files, src_rad=1.0, bkg_inner_rad=2.5, bkg_outer_rad=5.0,
                     ra=299.59, dec=35.2, irf_date="20230101")
    assert sorted(p.name for p in event_l2.iterdir()) == before
