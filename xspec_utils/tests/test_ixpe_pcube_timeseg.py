"""Tests for xspec_utils' ixpe_pcube_timeseg (lib/red_ixpe.py + skills/session.py).

Uses the real, already-local Level-2 event files for IXPE ObsID 01250401
(the real validation ObsID) as fixture data, copied into an
isolated tmp directory per test so the reduction runs on a scratch copy
rather than mutating the shared analyst-env data directory. segment_duration_s
is set to split the observation into exactly 2 segments (rather than the
tool's ~170-segment production default) to keep the real ixpeobssim pipeline
fast while still exercising genuine multi-segment behaviour.
"""

from __future__ import annotations

import csv
import glob
import math
import os
import shutil

import numpy as np
import pytest
from astropy.io import fits

from obsdb_catalog import neutral
from xspec_utils.lib import red_ixpe
from xspec_utils.skills import session as ixpe_session
from xspec_utils.skills import TOOL_MANIFEST

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
OBSID = "01250401"
IRF_DATE = "20230702"
# Half the observation's elapsed span -> exactly 2 segments (3 real
# xpselect+xpbin reductions total per call, instead of ~170 at the tool's
# production default segment_duration_s=5000s).
SEGMENT_DURATION_S = 424200.0

# Published 3-DU-merged accuracy anchor for this ObsID.
PUBLISHED_PD_PCT, PUBLISHED_PD_ERR = 8.33, 0.17
PUBLISHED_PA_DEG, PUBLISHED_PA_ERR = 17.78, 0.60


def _copy_real_event_files(dest_dir: str) -> list[str]:
    """Copy the real, already-local 01250401 Level-2 DU1/2/3 event files into
    ``<dest_dir>/<OBSID>/event_l2/`` and return the sorted file list."""
    src_evt_dir = os.path.join(REAL_IXPE_DIR, OBSID, "event_l2")
    dst_evt_dir = os.path.join(dest_dir, OBSID, "event_l2")
    os.makedirs(dst_evt_dir, exist_ok=True)
    file_list = []
    for det in (1, 2, 3):
        src = os.path.join(src_evt_dir, f"ixpe{OBSID}_det{det}_evt2_v03.fits")
        if not os.path.exists(src):
            pytest.skip(f"real fixture ObsID {OBSID} not local at {src}")
        dst = os.path.join(dst_evt_dir, os.path.basename(src))
        shutil.copyfile(src, dst)
        file_list.append(dst)
    return sorted(file_list)


# --------------------------------------------------------------------------- #
# Fast, pure-logic unit tests for the segmentation / GTI-overlap helpers
# --------------------------------------------------------------------------- #

def test_segment_edges_exact_division():
    assert red_ixpe._segment_edges(0.0, 100.0, 25.0) == [0.0, 25.0, 50.0, 75.0, 100.0]


def test_segment_edges_partial_remainder():
    # 90 / 25 -> 3 full segments (0-25-50-75) + a partial final edge at 90.
    assert red_ixpe._segment_edges(0.0, 90.0, 25.0) == [0.0, 25.0, 50.0, 75.0, 90.0]


def test_segment_edges_duration_exceeds_span_gives_single_segment():
    assert red_ixpe._segment_edges(10.0, 60.0, 1000.0) == [10.0, 60.0]


def test_gti_overlap_partial_window_sums_disjoint_intervals():
    # Two GTIs [0,20] and [50,80]; window [10,60) clips each to [10,20) and
    # [50,60) -- 10s + 10s of real overlap, tightest bound [10,60).
    clip = red_ixpe._gti_overlap([0.0, 50.0], [20.0, 80.0], 10.0, 60.0)
    assert clip == (10.0, 60.0, 20.0)


def test_gti_overlap_no_overlap_returns_none():
    assert red_ixpe._gti_overlap([100.0], [110.0], 0.0, 50.0) is None


def test_gti_overlap_window_fully_inside_one_gti():
    clip = red_ixpe._gti_overlap([0.0], [100.0], 20.0, 30.0)
    assert clip == (20.0, 30.0, 10.0)


# --------------------------------------------------------------------------- #
# Manifest wiring (AI-facing half)
# --------------------------------------------------------------------------- #

def test_tool_manifest_exposes_class_parameters():
    """A caller must be able to discover and override the axes this tool's
    science varies over (segment_duration_s, energy band, region, background)
    from the manifest alone -- not just from the lib docstring."""
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.ixpe_pcube_timeseg"]
    assert entry["module"] == "xspec_utils.skills.session"
    assert entry["fn"] == "ixpe_pcube_timeseg"
    props = entry["args_schema"]["properties"]
    for name in ("segment_duration_s", "emin_kev", "emax_kev", "src_rad", "bkg",
                 "bkg_inner_rad", "bkg_outer_rad", "ra", "dec", "irf_date"):
        assert name in props, f"{name} missing from args_schema.properties"
    assert entry["args_schema"]["required"] == ["obs_db_path", "obsid"]
    assert props["bkg"]["type"] == "boolean"
    assert props["segment_duration_s"]["type"] == "number"
    assert "5000" in props["segment_duration_s"]["description"]
    assert "false" in props["bkg"]["description"].lower()


# --------------------------------------------------------------------------- #
# Silent-misuse guards -- must fail loudly before any reduction is attempted,
# so these run in milliseconds (file-list validation + the ra/dec XOR check
# both happen before find_src / fits.open / xpselect / xpbin).
# --------------------------------------------------------------------------- #

def test_missing_du_raises(tmp_path):
    """Dropping DU3 from the file list must raise rather than silently
    reduce on 2 DUs."""
    file_list = _copy_real_event_files(str(tmp_path))
    two_du_list = [f for f in file_list if "_det3_" not in f]
    assert len(two_du_list) == 2

    with pytest.raises(ValueError) as excinfo:
        red_ixpe.pcube_timeseg(OBSID, two_du_list, IRF_DATE)
    assert "DU3" in str(excinfo.value)


def test_duplicate_du_version_raises(tmp_path):
    """Two Level-2 files resolving to the same DU (e.g. two calibration
    versions) must raise rather than silently picking one."""
    file_list = _copy_real_event_files(str(tmp_path))
    evt_dir = os.path.dirname(file_list[0])
    det1_dup = os.path.join(evt_dir, f"ixpe{OBSID}_det1_evt2_v99.fits")
    shutil.copyfile(os.path.join(evt_dir, f"ixpe{OBSID}_det1_evt2_v03.fits"), det1_dup)

    with pytest.raises(ValueError) as excinfo:
        red_ixpe.pcube_timeseg(OBSID, file_list + [det1_dup], IRF_DATE)
    assert "multiple" in str(excinfo.value)
    assert "DU1" in str(excinfo.value)


def test_partial_radec_raises(tmp_path):
    """Exactly one of ra/dec given (without the other) must raise -- fast,
    because the check runs before any auto-centroid attempt or file open."""
    file_list = _copy_real_event_files(str(tmp_path))

    with pytest.raises(ValueError, match="both ra and dec"):
        red_ixpe.pcube_timeseg(OBSID, file_list, IRF_DATE, ra=299.59, dec=None)
    with pytest.raises(ValueError, match="both ra and dec"):
        red_ixpe.pcube_timeseg(OBSID, file_list, IRF_DATE, ra=None, dec=35.2)


# --------------------------------------------------------------------------- #
# Real-data reduction: lib.red_ixpe.pcube_timeseg
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def timeseg_result(tmp_path_factory):
    """One real pcube_timeseg call on a scratch copy of 01250401, split into
    exactly 2 elapsed-time segments. Shared read-only across the tests below
    to avoid re-running the real ixpeobssim pipeline more than necessary."""
    tmp_dir = tmp_path_factory.mktemp("pcube_timeseg")
    file_list = _copy_real_event_files(str(tmp_dir))
    result = red_ixpe.pcube_timeseg(
        OBSID, file_list, IRF_DATE, bkg=False, segment_duration_s=SEGMENT_DURATION_S)
    return result, file_list


def test_two_segments_produced(timeseg_result):
    result, _ = timeseg_result
    assert len(result["timeseries"]) == 2
    assert {row["segment_id"] for row in result["timeseries"]} == {0, 1}
    assert result["integrated"]["segment_id"] == -1


def test_capability_segments_are_gti_clipped_not_full_obs(timeseg_result):
    """Pins the 03002599-style defect: each segment's actual elapsed/GTI
    window must be strictly narrower than the integrated whole-observation
    window, and segments must not silently inherit the identical full-obs
    exposure_s (the copied-header defect this tool exists to fix)."""
    result, _ = timeseg_result
    integrated = result["integrated"]
    full_span = integrated["t_stop_mjd"] - integrated["t_start_mjd"]

    exposures = []
    for row in result["timeseries"]:
        seg_span = row["t_stop_mjd"] - row["t_start_mjd"]
        assert seg_span < full_span
        exposures.append(row["exposure_s"])

    assert len(set(round(e, 3) for e in exposures)) > 1, (
        "segments must not all report the identical (full-obs) exposure_s"
    )
    # Segment exposures partition the integrated exposure -- no double count,
    # no silent inflation (small tolerance for floating point).
    assert sum(exposures) <= integrated["exposure_s"] * 1.001


def test_accuracy_anchor_matches_published_value(timeseg_result):
    """Integrated (whole-obs, bkg=False, 2-8 keV, DU-merged) PD/PA agree with
    the published 3-DU-merged anchor within a 3-sigma combined tolerance --
    loose enough to be stable, tight enough to catch a wrong DU-merge,
    energy-band, or region bug."""
    result, _ = timeseg_result
    integrated = result["integrated"]

    pd_sigma = math.hypot(PUBLISHED_PD_ERR, integrated["PD_err"])
    assert abs(integrated["PD_pct"] - PUBLISHED_PD_PCT) < 3 * pd_sigma

    pa_sigma = math.hypot(PUBLISHED_PA_ERR, integrated["PA_err"])
    assert abs(integrated["PA_deg"] - PUBLISHED_PA_DEG) < 3 * pa_sigma


def test_differential_readback_from_segment_fits(timeseg_result):
    """Independently re-derive PD/PA for segment 0 straight from the per-DU
    PCUBE FITS files with astropy + ixpeobssim's own Kislat (2015) formula --
    NOT via the tool's xBinnedPolarizationCube reader -- and check agreement
    with the tool-returned segment values."""
    from ixpeobssim.evt.kislat2015 import xStokesAnalysis

    result, file_list = timeseg_result
    seg0 = next(row for row in result["timeseries"] if row["segment_id"] == 0)

    evt_dir = os.path.dirname(file_list[0])
    pcube_files = sorted(glob.glob(os.path.join(evt_dir, "*_seg000src_pcube.fits")))
    assert len(pcube_files) == 3, "expected one segment-0 PCUBE per DU"

    counts = i_tot = q_tot = u_tot = w2_tot = mu_num = 0.0
    for path in pcube_files:
        with fits.open(path) as hdul:
            data = hdul["POLARIZATION"].data
            counts += float(data["COUNTS"][0])
            i_tot += float(data["I"][0])
            q_tot += float(data["Q"][0])
            u_tot += float(data["U"][0])
            w2_tot += float(data["W2"][0])
            mu_num += float(data["MU"][0]) * float(data["I"][0])
    mu = mu_num / i_tot

    pd, _pd_err, pa, _pa_err = xStokesAnalysis.calculate_polarization(
        np.array([i_tot]), np.array([q_tot]), np.array([u_tot]),
        np.array([mu]), np.array([w2_tot]), degrees=True)

    assert pd[0] * 100.0 == pytest.approx(seg0["PD_pct"], abs=0.01)
    assert pa[0] == pytest.approx(seg0["PA_deg"], abs=0.01)


def test_csv_matches_returned_rows(timeseg_result):
    result, _ = timeseg_result
    with open(result["csv_path"], newline="") as f:
        rows = list(csv.DictReader(f))
    assert list(rows[0].keys()) == red_ixpe._TIMESEG_FIELDS
    assert len(rows) == len(result["timeseries"]) + 1
    assert int(rows[-1]["segment_id"]) == -1
    assert float(rows[-1]["PD_pct"]) == pytest.approx(result["integrated"]["PD_pct"])


def test_repro_pcube_timeseg_stable_across_calls(tmp_path_factory):
    """Reproducibility: two independent pcube_timeseg calls on independent
    scratch copies with identical arguments give the same PD/PA/MDP99 within
    tolerance, for both the integrated cube and every segment."""
    results = []
    for _ in range(2):
        tmp_dir = tmp_path_factory.mktemp("pcube_timeseg_repro")
        file_list = _copy_real_event_files(str(tmp_dir))
        results.append(red_ixpe.pcube_timeseg(
            OBSID, file_list, IRF_DATE, bkg=False, segment_duration_s=SEGMENT_DURATION_S))

    a, b = results
    assert len(a["timeseries"]) == len(b["timeseries"]) == 2
    for field in ("PD_pct", "PA_deg", "MDP99"):
        assert a["integrated"][field] == pytest.approx(b["integrated"][field], rel=1e-2)
        for row_a, row_b in zip(a["timeseries"], b["timeseries"]):
            assert row_a["segment_id"] == row_b["segment_id"]
            assert row_a[field] == pytest.approx(row_b[field], rel=1e-2)


# --------------------------------------------------------------------------- #
# Skill wrapper: file discovery, CSV placement, catalog recording
# --------------------------------------------------------------------------- #

def test_skill_wrapper_records_catalog_and_falls_back_to_plain_evt2(tmp_path):
    """01250401's real files have no '_rej' variant (mirrors 03002599) -- the
    skill wrapper's glob must fall back to the plain evt2_v*.fits pattern.
    Confirms end-to-end: file discovery, reduction, CSV on disk, and
    level-3 catalog recording."""
    obs_db_path = str(tmp_path)
    _copy_real_event_files(obs_db_path)

    result = ixpe_session.ixpe_pcube_timeseg(
        obs_db_path, OBSID, irf_date=IRF_DATE, bkg=False,
        segment_duration_s=SEGMENT_DURATION_S)

    assert result["exit_code"] == 0, result.get("stderr")
    assert os.path.exists(result["csv_path"])
    assert len(result["timeseries"]) == 2

    entry = neutral.read(obs_db_path)["entries"][OBSID]
    assert entry["max_level"] == 3
    (prod,) = [p for p in entry["products"] if p["kind"] == "pcube_timeseg"]
    assert prod["path"] == f"{OBSID}/pcube_timeseg_{OBSID}.csv"
    assert prod["produced_by"] == "xspec_utils.red_ixpe.pcube_timeseg"
    assert prod["params"]["bkg"] is False
    assert prod["params"]["segment_duration_s"] == SEGMENT_DURATION_S


def test_wrapper_reinvocation_ignores_generated_products(tmp_path):
    """Regression: the wrapper's fallback glob must match only original
    Level-2 event files, never this tool's own generated products
    (..._seg000src.fits, ..._seg000src_pcube.fits, ..._intsrc.fits, ...).
    Without the v[0-9][0-9] anchor, a second run on a no-'_rej' ObsID like
    01250401 rediscovers those generated files as inputs and
    _validate_du_files then sees multiple files per DU and fails."""
    obs_db_path = str(tmp_path)
    _copy_real_event_files(obs_db_path)

    first = ixpe_session.ixpe_pcube_timeseg(
        obs_db_path, OBSID, irf_date=IRF_DATE, bkg=False,
        segment_duration_s=SEGMENT_DURATION_S)
    assert first["exit_code"] == 0, first.get("stderr")
    assert len(first["timeseries"]) == 2

    # Direct discovery check: the generated _seg*/_pcube files now sitting in
    # event_l2/ (written by the first run) must NOT be matched by the same
    # strict pattern the wrapper uses.
    evt_dir = os.path.join(obs_db_path, OBSID, "event_l2")
    strict_glob = sorted(glob.glob(os.path.join(evt_dir, f"ixpe{OBSID}_det[123]_evt2_v[0-9][0-9].fits")))
    assert len(strict_glob) == 3, (
        f"expected exactly 3 original Level-2 files (one per DU), got {strict_glob}"
    )

    second = ixpe_session.ixpe_pcube_timeseg(
        obs_db_path, OBSID, irf_date=IRF_DATE, bkg=False,
        segment_duration_s=SEGMENT_DURATION_S)
    assert second["exit_code"] == 0, second.get("stderr")
    assert len(second["timeseries"]) == 2


def test_skill_wrapper_missing_event_files_fails_loudly_no_record(tmp_path):
    result = ixpe_session.ixpe_pcube_timeseg(str(tmp_path), "99999999")
    assert result["exit_code"] == 1
    assert "stderr" in result
    assert neutral.read(str(tmp_path))["entries"] == {}
