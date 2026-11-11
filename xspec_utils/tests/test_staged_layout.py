"""The xspec_utils readers find an observation's files in its folder staged from the store (laid out as
HEASARC's ObsID folder), and load_data refuses records it cannot match to the epoch before staging."""

from __future__ import annotations

from pathlib import Path

import pytest

from xspec_utils.lib import NICER, NuSTAR, red_ixpe


def _touch(path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"x")
    return str(path)


def test_nicer_files_are_found_in_the_observation_folder(tmp_path):
    obs = "4690010101"
    cl = tmp_path / "xti/event_cl"
    src, bkg = _touch(cl / f"ni{obs}mpu7_sr.pha"), _touch(cl / f"ni{obs}mpu7_bg.pha")
    rmf, arf = _touch(cl / f"ni{obs}mpu7.rmf"), _touch(cl / f"ni{obs}mpu7.arf")
    assert NICER.pha(obs, str(tmp_path)) == (src, bkg)
    assert (NICER.rmf(obs, str(tmp_path)), NICER.arf(obs, str(tmp_path))) == (rmf, arf)


def test_nustar_files_are_found_in_the_observation_folder(tmp_path):
    obs = "30702017002"
    products = tmp_path / "products"
    src = [_touch(products / f"nu{obs}{m}01_sr_grp.pha") for m in "AB"]
    bkg = [_touch(products / f"nu{obs}{m}01_bk.pha") for m in "AB"]
    rmf = [_touch(products / f"nu{obs}{m}01_sr.rmf") for m in "AB"]
    assert NuSTAR.pha(obs, str(tmp_path)) == (src, bkg)
    assert NuSTAR.rmf(obs, str(tmp_path)) == rmf


def test_ixpe_level2_files_and_spectra_are_listed_from_the_observation_folder(tmp_path):
    obs = "02008201"
    files = [_touch(tmp_path / "event_l2" / f"ixpe{obs}_det{du}_evt2_v01_rej.fits") for du in (1, 2, 3)]
    assert red_ixpe.get_all_file_list(str(tmp_path), obs) == files
    src = [_touch(tmp_path / "event_l2" / f"ixpe{obs}_det{du}_evt2_v01_rej_src_pha1.fits") for du in (1, 2, 3)]
    bkg = [_touch(tmp_path / "event_l2" / f"ixpe{obs}_det{du}_evt2_v01_rej_bkg_pha1.fits") for du in (1, 2, 3)]
    _touch(tmp_path / "event_l2" / f"ixpe{obs}_det1_evt2_v01_rej_src_pha1q.fits")  # Q/U follow each I file
    assert red_ixpe.spectra(str(tmp_path), obs) == (src, bkg)
    Path(bkg[2]).unlink()
    with pytest.raises(ValueError, match="one bkg PHA1 spectrum of 02008201 DU3"):
        red_ixpe.spectra(str(tmp_path), obs)


def test_ixpe_spectra_of_the_rejected_events_are_read_where_plain_ones_lie_beside_them(tmp_path):
    obs = "03002599"
    rej = [_touch(tmp_path / f"ixpe{obs}_det{du}_evt2_v01_rej_{role}_pha1.fits") for role in ("src", "bkg")
           for du in (1, 2, 3)]
    for role in ("src", "bkg"):
        for du in (1, 2, 3):
            _touch(tmp_path / f"ixpe{obs}_det{du}_evt2_v01_{role}_pha1.fits")
    assert red_ixpe.spectra(str(tmp_path), obs) == (rej[:3], rej[3:])


def test_ixpe_without_level2_events_says_where_it_looked(tmp_path):
    with pytest.raises(ValueError, match="No level-2 event files of 01002901"):
        red_ixpe.get_all_file_list(str(tmp_path), "01002901")


@pytest.mark.parametrize("records, said", [
    ([], "pass one record each"),
    ({"store": "db", "bucket": "ixpe", "key": "pha/01002901/abcd1234/x.fits", "version_id": "v"},
     "pass the record of a product"),
    ({"store": "db", "bucket": "ixpe", "key": "pha/02008201/abcd1234/recipe.json", "version_id": "v"},
     "is a product of"),
    ({"store": "db", "bucket": "ixpe", "key": "pha/01002901/abcd1234/recipe.json"}, "a record is"),
])
def test_load_data_refuses_records_it_cannot_match_to_the_epoch(tmp_path, records, said):
    from xspec_utils.lib import xspec_help

    with pytest.raises(ValueError, match=said):
        xspec_help.load_data({"IXPE": ["01002901"], "STATE": "test"}, ixpe=records, staged=str(tmp_path))
