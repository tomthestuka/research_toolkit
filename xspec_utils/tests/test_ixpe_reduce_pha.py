"""Tests for xspec_utils' ixpe_reduce_pha (lib/red_ixpe.reduce_pha + skills/session.py).

Uses the real, already-local Level-2 event files for IXPE ObsID 01002901 --
one of the twelve Cyg X-1 ObsIDs that already carries real, already-reduced
per-DU PHA1 products on disk -- copied into an isolated tmp directory per
test so this tool's reduction (or reuse check) runs on a scratch copy rather
than mutating the shared analyst-env data directory.
"""

from __future__ import annotations

import os
import re
import shutil

import numpy as np
import pytest
from astropy.io import fits

from obsdb_catalog import neutral
from xspec_utils.lib import find_src, red_ixpe
from xspec_utils.skills import TOOL_MANIFEST
from xspec_utils.skills import session as ixpe_session

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
OBSID = "01002901"
VERSION = "v01"
IRF_DATE = "20230702"


def _copy_event_files(dest_dir: str) -> list[str]:
    """Copy the real, already-local Level-2 (``_rej``) DU1/2/3 event files
    for OBSID into ``<dest_dir>/<OBSID>/event_l2/`` and return the sorted
    list."""
    src_evt_dir = os.path.join(REAL_IXPE_DIR, OBSID, "event_l2")
    dst_evt_dir = os.path.join(dest_dir, OBSID, "event_l2")
    os.makedirs(dst_evt_dir, exist_ok=True)
    file_list = []
    for det in (1, 2, 3):
        src = os.path.join(src_evt_dir, f"ixpe{OBSID}_det{det}_evt2_{VERSION}_rej.fits")
        if not os.path.exists(src):
            pytest.skip(f"real fixture ObsID {OBSID} not local at {src}")
        dst = os.path.join(dst_evt_dir, os.path.basename(src))
        shutil.copyfile(src, dst)
        file_list.append(dst)
    return sorted(file_list)


def _copy_existing_pha1(dest_dir: str) -> None:
    """Copy OBSID's real, already-reduced per-DU Stokes-I PHA1 src+bkg
    products alongside the copied event files, so the reuse path has
    genuine on-disk products to find (mirrors the real state of 01002901)."""
    src_evt_dir = os.path.join(REAL_IXPE_DIR, OBSID, "event_l2")
    dst_evt_dir = os.path.join(dest_dir, OBSID, "event_l2")
    for det in (1, 2, 3):
        for suffix in ("src_pha1", "bkg_pha1"):
            name = f"ixpe{OBSID}_det{det}_evt2_{VERSION}_rej_{suffix}.fits"
            src = os.path.join(src_evt_dir, name)
            if not os.path.exists(src):
                pytest.skip(f"real fixture PHA1 product not local at {src}")
            shutil.copyfile(src, os.path.join(dst_evt_dir, name))


def _rate_sum(pha_path: str) -> float:
    with fits.open(pha_path) as hdul:
        return float(np.nansum(hdul["SPECTRUM"].data["RATE"]))


# --------------------------------------------------------------------------- #
# Silent-misuse guard -- must fail before any pha()/find_src call
# --------------------------------------------------------------------------- #

def test_missing_du_raises(tmp_path):
    file_list = _copy_event_files(str(tmp_path))
    two_du_list = [f for f in file_list if "_det3_" not in f]
    assert len(two_du_list) == 2

    with pytest.raises(ValueError, match="DU3"):
        red_ixpe.reduce_pha(OBSID, two_du_list, irf_date=IRF_DATE)


# --------------------------------------------------------------------------- #
# Reuse path -- 01002901 already has real per-DU PHA1 products on disk
# --------------------------------------------------------------------------- #

def test_capability_reuse_path_resolves_existing_pha(tmp_path):
    """The reuse path must actually resolve the real, already-reduced PHA1
    products -- not just return well-shaped metadata -- and must not derive
    a source position or write region files, since no reduction runs."""
    dest = str(tmp_path)
    file_list = _copy_event_files(dest)
    _copy_existing_pha1(dest)

    result = red_ixpe.reduce_pha(OBSID, file_list, irf_date=IRF_DATE)

    assert result["reused"] is True
    assert result["ra"] is None and result["dec"] is None
    assert result["src_reg"] is None and result["bkg_reg"] is None
    assert len(result["src_pha"]) == 3
    assert len(result["bkg_pha"]) == 3

    for f in result["src_pha"] + result["bkg_pha"]:
        assert f.endswith("_pha1.fits")
        assert os.path.exists(f)
        assert _rate_sum(f) > 0, f"{f}: real PHA1 product has zero net rate"


# --------------------------------------------------------------------------- #
# Fresh-reduction path -- force the reuse check to miss by omitting the
# existing PHA1 products, so the tool must centroid via find_src and run a
# real xpselect + xpbin PHA1/PHA1Q/PHA1U reduction.
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def fresh_reduction(tmp_path_factory):
    """One real reduce_pha call on a scratch copy of 01002901 carrying only
    the Level-2 event files (no existing PHA1 products). Shared read-only
    across the tests below to avoid re-running the real ixpeobssim pipeline
    more than necessary."""
    tmp_dir = tmp_path_factory.mktemp("reduce_pha_fresh")
    file_list = _copy_event_files(str(tmp_dir))
    result = red_ixpe.reduce_pha(OBSID, file_list, irf_date=IRF_DATE)
    return result, file_list


def test_capability_fresh_reduction_produces_real_pha_and_regions(fresh_reduction):
    result, _file_list = fresh_reduction
    assert result["reused"] is False
    assert isinstance(result["ra"], float) and isinstance(result["dec"], float)

    assert os.path.exists(result["src_reg"])
    assert os.path.exists(result["bkg_reg"])
    with open(result["src_reg"]) as f:
        src_lines = f.readlines()
    assert src_lines[0].strip() == "fk5"
    assert src_lines[1].startswith("circle(")
    with open(result["bkg_reg"]) as f:
        bkg_lines = f.readlines()
    assert bkg_lines[1].startswith("annulus(")

    assert len(result["src_pha"]) == 3 and len(result["bkg_pha"]) == 3
    for f in result["src_pha"] + result["bkg_pha"]:
        assert os.path.exists(f)
        assert _rate_sum(f) > 0, f"{f}: freshly reduced PHA1 product has zero net rate"


def test_differential_centroid_matches_independent_find_src_call(fresh_reduction, tmp_path_factory):
    """The RA/Dec fed to the reduction must match an independent direct call
    to find_src.run_extraction on the same Level-2 DU1 event file -- not the
    tool's own wrapper -- to a tight tolerance (the centroid is deterministic
    and independent of src_rad/bkg_*_rad, which only set the region sizes)."""
    result, file_list = fresh_reduction
    du1_evt = next(f for f in file_list if "_det1_" in f)

    ref_dir = str(tmp_path_factory.mktemp("find_src_ref"))
    ref_src_reg = os.path.join(ref_dir, "src.reg")
    ref_bkg_reg = os.path.join(ref_dir, "bkg.reg")
    find_src.run_extraction(du1_evt, ref_src_reg, ref_bkg_reg, None, 1.0, 3.0, 5.0)

    with open(ref_src_reg) as f:
        region_line = f.readlines()[1]
    m = re.match(r'circle\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,', region_line)
    ref_ra, ref_dec = float(m.group(1)), float(m.group(2))

    assert result["ra"] == pytest.approx(ref_ra, abs=1e-6)
    assert result["dec"] == pytest.approx(ref_dec, abs=1e-6)


def test_repro_reduce_pha_stable_across_calls(fresh_reduction, tmp_path_factory):
    """Reproducibility: two independent fresh reduce_pha calls on independent
    scratch copies with identical arguments give the same centroided RA/Dec
    and the same net PHA1 rates (the underlying xpselect/xpbin pipeline is
    deterministic on identical, unchanged input data)."""
    a, _ = fresh_reduction
    tmp_dir = tmp_path_factory.mktemp("reduce_pha_repro")
    file_list = _copy_event_files(str(tmp_dir))
    b = red_ixpe.reduce_pha(OBSID, file_list, irf_date=IRF_DATE)

    assert a["reused"] is False and b["reused"] is False
    assert a["ra"] == pytest.approx(b["ra"], abs=1e-6)
    assert a["dec"] == pytest.approx(b["dec"], abs=1e-6)

    for key in ("src_pha", "bkg_pha"):
        for fa, fb in zip(a[key], b[key]):
            assert _rate_sum(fa) == pytest.approx(_rate_sum(fb), rel=1e-6)


# --------------------------------------------------------------------------- #
# Manifest wiring (AI-facing half)
# --------------------------------------------------------------------------- #

def test_tool_manifest_exposes_class_parameters():
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.ixpe_reduce_pha"]
    assert entry["module"] == "xspec_utils.skills.session"
    assert entry["fn"] == "ixpe_reduce_pha"
    props = entry["args_schema"]["properties"]
    for name in ("obs_db_path", "obsid", "src_rad", "bkg_inner_rad", "bkg_outer_rad", "irf_date"):
        assert name in props, f"{name} missing from args_schema.properties"
    assert entry["args_schema"]["required"] == ["obs_db_path", "obsid"]
    assert "1.0" in props["src_rad"]["description"]
    assert "3.0" in props["bkg_inner_rad"]["description"]
    assert "5.0" in props["bkg_outer_rad"]["description"]
    assert "reused" in entry["description"]
    assert "find_src" in entry["description"]


# --------------------------------------------------------------------------- #
# Skill wrapper: file discovery, catalog recording, failure contract
# --------------------------------------------------------------------------- #

def test_skill_wrapper_reuse_records_catalog(tmp_path):
    obs_db_path = str(tmp_path)
    _copy_event_files(obs_db_path)
    _copy_existing_pha1(obs_db_path)

    result = ixpe_session.ixpe_reduce_pha(obs_db_path, OBSID, irf_date=IRF_DATE)

    assert result["exit_code"] == 0, result.get("stderr")
    assert result["reused"] is True
    assert result["ra"] is None and result["dec"] is None
    for f in result["src_pha"] + result["bkg_pha"]:
        assert os.path.exists(f)

    entry = neutral.read(obs_db_path)["entries"][OBSID]
    assert entry["max_level"] == 3
    pha_products = [p for p in entry["products"] if p["kind"] == "pha"]
    assert len(pha_products) == 6
    for p in pha_products:
        assert p["produced_by"] == "xspec_utils.red_ixpe.pha"
        assert p["params"]["reused"] is True
        assert not os.path.isabs(p["path"])
        assert os.path.exists(os.path.join(obs_db_path, p["path"]))


def test_skill_wrapper_fresh_reduction_records_catalog_with_derived_position(tmp_path):
    obs_db_path = str(tmp_path)
    _copy_event_files(obs_db_path)

    result = ixpe_session.ixpe_reduce_pha(obs_db_path, OBSID, irf_date=IRF_DATE)

    assert result["exit_code"] == 0, result.get("stderr")
    assert result["reused"] is False
    assert isinstance(result["ra"], float) and isinstance(result["dec"], float)
    assert os.path.exists(result["src_reg"]) and os.path.exists(result["bkg_reg"])

    entry = neutral.read(obs_db_path)["entries"][OBSID]
    assert entry["max_level"] == 3
    pha_products = [p for p in entry["products"] if p["kind"] == "pha"]
    assert len(pha_products) == 6
    dus = {p["params"]["du"] for p in pha_products}
    assert dus == {"DU1", "DU2", "DU3"}
    for p in pha_products:
        assert p["params"]["reused"] is False
        assert p["params"]["irf_date"] == IRF_DATE


def test_skill_wrapper_missing_event_files_fails_loudly_no_record(tmp_path):
    result = ixpe_session.ixpe_reduce_pha(str(tmp_path), "99999999")
    assert result["exit_code"] == 1
    assert "stderr" in result
    assert neutral.read(str(tmp_path))["entries"] == {}
