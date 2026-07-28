"""Tests for xspec_utils' ixpe_pcube_timeseg (lib/red_ixpe.py + skills/session.py).

Uses the real, already-local Level-2 event files for IXPE ObsIDs 01250401 and
02001901 (both real validation ObsIDs) as fixture data, copied into an
isolated tmp directory per test so the reduction runs on a scratch copy
rather than mutating the shared analyst-env data directory. Several
segment_duration_s values are used across tests -- not the tool's ~150-170
segment production default -- to keep the real ixpeobssim pipeline fast
while still exercising genuine multi-segment, gate, resume and override
behaviour on real data.
"""

from __future__ import annotations

import csv
import glob
import inspect
import json
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
OBSID2 = "02001901"
IRF_DATE = "20230702"
# Half the observation's elapsed span -> exactly 2 segments (3 real
# xpselect+xpbin reductions total per call, instead of ~170 at the tool's
# production default segment_duration_s=5000s).
SEGMENT_DURATION_S = 424200.0
# Span-minus-500s: the observation's own TSTOP is 500s short of an exact
# multiple, so this produces exactly 2 windows -- one ~848000s window
# carrying essentially the whole observation's exposure, and one genuine
# 500s-wide tail window whose real GTI overlap (494.9s, measured directly
# from the event file's GTI extension) sits just under the default
# min_seg_exposure_s=1000s. Used both to exercise the GTI gate cheaply (one
# real reduction skipped without a xpselect/xpbin call) and, in the repro
# test, to keep the real pipeline call count minimal.
GATE_SEGMENT_DURATION_S = 847900.0

# Published 3-DU-merged accuracy anchor for this ObsID.
PUBLISHED_PD_PCT, PUBLISHED_PD_ERR = 8.33, 0.17
PUBLISHED_PA_DEG, PUBLISHED_PA_ERR = 17.78, 0.60


def _copy_real_event_files(dest_dir: str, obsid: str = OBSID, version: str = "v03") -> list[str]:
    """Copy the real, already-local Level-2 DU1/2/3 event files for ``obsid``
    into ``<dest_dir>/<obsid>/event_l2/`` and return the sorted file list."""
    src_evt_dir = os.path.join(REAL_IXPE_DIR, obsid, "event_l2")
    dst_evt_dir = os.path.join(dest_dir, obsid, "event_l2")
    os.makedirs(dst_evt_dir, exist_ok=True)
    file_list = []
    for det in (1, 2, 3):
        src = os.path.join(src_evt_dir, f"ixpe{obsid}_det{det}_evt2_{version}.fits")
        if not os.path.exists(src):
            pytest.skip(f"real fixture ObsID {obsid} not local at {src}")
        dst = os.path.join(dst_evt_dir, os.path.basename(src))
        shutil.copyfile(src, dst)
        file_list.append(dst)
    return sorted(file_list)


def _real_segment_exposures(evt_path: str, duration: float) -> list[float]:
    """GTI-clipped, DEADC-scaled good exposure per elapsed-time window, using
    the production ``_segment_edges``/``_gti_overlap`` helpers directly
    against one real event file's own EVENTS/GTI/PRIMARY headers -- exactly
    the computation ``pcube_timeseg``'s gate performs before ever calling
    xpselect/xpbin."""
    with fits.open(evt_path) as hdul:
        hdr = hdul["EVENTS"].header
        tstart, tstop = float(hdr["TSTART"]), float(hdr["TSTOP"])
        gti = hdul["GTI"].data
        gti_start = [float(x) for x in gti["START"]]
        gti_stop = [float(x) for x in gti["STOP"]]
        deadc = float(hdul["PRIMARY"].header["DEADC"])
    edges = red_ixpe._segment_edges(tstart, tstop, duration)
    exposures = []
    for i in range(len(edges) - 1):
        overlap = red_ixpe._gti_overlap(gti_start, gti_stop, edges[i], edges[i + 1])
        exposures.append(overlap[2] * deadc if overlap else 0.0)
    return exposures


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
# Fast, pure-logic unit test for the incremental-CSV / provenance helpers
# --------------------------------------------------------------------------- #

def test_provenance_csv_round_trip(tmp_path):
    """Round-trips the checkpoint CSV's provenance header and rows through
    the exact helpers pcube_timeseg uses -- pure I/O plumbing, no reduction."""
    csv_path = str(tmp_path / "pcube_timeseg_roundtrip.csv")
    provenance = {"segment_duration_s": 5000.0, "emin_kev": 2.0, "emax_kev": 8.0,
                  "src_rad": 1.0, "bkg": False, "bkg_inner_rad": 3.0, "bkg_outer_rad": 5.0,
                  "ra": 299.59, "dec": 35.2, "irf_date": "20230702",
                  "min_seg_exposure_s": 1000.0, "tstart": 100.0, "tstop": 6100.0}
    red_ixpe._write_timeseg_csv_header(csv_path, provenance)

    row0 = {f: 0.0 for f in red_ixpe._TIMESEG_FIELDS}
    row0.update(segment_id=0, PD_pct=8.04, PA_deg=17.7)
    red_ixpe._append_timeseg_csv_row(csv_path, row0)

    read_back = red_ixpe._read_existing_timeseg_csv(csv_path)
    assert read_back is not None
    read_provenance, rows_by_segid, integrated = read_back
    assert read_provenance == provenance
    assert integrated is None
    assert set(rows_by_segid) == {0}
    assert rows_by_segid[0]["PD_pct"] == pytest.approx(8.04)
    assert rows_by_segid[0]["PA_deg"] == pytest.approx(17.7)


def test_read_existing_csv_missing_provenance_header_treated_as_absent(tmp_path):
    """A CSV without the provenance comment line (e.g. hand-edited, or from
    a version of this tool predating the header) must not be silently
    resumed from -- callers see it as absent and rebuild from scratch."""
    csv_path = str(tmp_path / "pcube_timeseg_no_header.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=red_ixpe._TIMESEG_FIELDS)
        writer.writeheader()
    assert red_ixpe._read_existing_timeseg_csv(csv_path) is None


# --------------------------------------------------------------------------- #
# Manifest wiring (AI-facing half)
# --------------------------------------------------------------------------- #

def test_tool_manifest_exposes_class_parameters():
    """A caller must be able to discover and override the axes this tool's
    science varies over (segment_duration_s, energy band, region, background,
    the GTI gate, and resume) from the manifest alone -- not just from the
    lib docstring."""
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.ixpe_pcube_timeseg"]
    assert entry["module"] == "xspec_utils.skills.session"
    assert entry["fn"] == "ixpe_pcube_timeseg"
    props = entry["args_schema"]["properties"]
    for name in ("segment_duration_s", "emin_kev", "emax_kev", "src_rad", "bkg",
                 "bkg_inner_rad", "bkg_outer_rad", "ra", "dec",
                 "min_seg_exposure_s", "resume"):
        assert name in props, f"{name} missing from args_schema.properties"
    # The IRF epoch is auto-selected inside the tool, not an AI-facing knob (the
    # human-half red_ixpe.pcube_timeseg keeps an override).
    assert "irf_date" not in props
    assert entry["args_schema"]["required"] == ["obs_db_path", "obsid"]
    assert props["bkg"]["type"] == "boolean"
    assert props["segment_duration_s"]["type"] == "number"
    assert props["min_seg_exposure_s"]["type"] == "number"
    assert props["resume"]["type"] == "boolean"
    assert "5000" in props["segment_duration_s"]["description"]
    assert "1000" in props["min_seg_exposure_s"]["description"]
    assert "false" in props["bkg"]["description"].lower()
    # The return contract is bounded (no full timeseries key) -- must be discoverable
    # from the manifest text. Under the brief/args-schema split that prose lives in
    # args_schema["description"], so assert over all text the model sees for this tool.
    manifest_text = " ".join([
        entry["description"],
        entry["args_schema"].get("description", ""),
        *(p.get("description", "") for p in props.values()),
    ])
    assert "n_segments" in manifest_text
    assert "csv_path" in manifest_text
    assert "timeseries" not in manifest_text


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
    assert result["n_expected"] == 2
    assert result["n_skipped_empty"] == 0
    assert result["complete"] is True
    assert result["resumed"] is False


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


def test_capability_segment_counts_partition_integrated_from_fits(timeseg_result):
    """Closes a coverage gap in test_capability_segments_are_gti_clipped_not_full_obs:
    that test only checks t_start_mjd/t_stop_mjd/exposure_s, which are computed
    by _gti_overlap directly from the window edges (t0, t1) -- independently of
    whatever xpselect actually selected. A bug that silently dropped the
    tmin/tmax kwargs from _reduce's sel_kwargs (so every "segment" reduced the
    whole observation) would leave that arithmetic untouched and pass anyway.

    This test instead reads COUNTS straight from each segment's raw per-DU
    PCUBE FITS files -- NOT via xBinnedPolarizationCube, and NOT via any
    window-edge arithmetic -- and proves the segment reductions actually used
    a strictly narrower, time-filtered event selection than the integrated
    reduction: each segment's summed COUNTS is strictly less than the
    integrated cube's, and the two segments' COUNTS partition (sum to) the
    integrated COUNTS. A dropped-tmin/tmax bug would make every segment's
    COUNTS equal the integrated COUNTS and their sum roughly double it --
    this test fails hard in that case.
    """
    result, file_list = timeseg_result
    evt_dir = os.path.dirname(file_list[0])

    def _summed_counts(glob_suffix: str) -> float:
        pcube_files = sorted(glob.glob(os.path.join(evt_dir, f"*_{glob_suffix}_pcube.fits")))
        assert len(pcube_files) == 3, f"expected one {glob_suffix} PCUBE per DU, got {pcube_files}"
        total = 0.0
        for path in pcube_files:
            with fits.open(path) as hdul:
                total += float(hdul["POLARIZATION"].data["COUNTS"][0])
        return total

    integrated_counts = _summed_counts("intsrc")

    seg_counts = {}
    for row in result["timeseries"]:
        i = row["segment_id"]
        seg_counts[i] = _summed_counts(f"seg{i:03d}src")

    for i, counts in seg_counts.items():
        assert counts > 0, f"segment {i} has zero raw PCUBE counts"
        assert counts < integrated_counts, (
            f"segment {i} counts ({counts}) are not less than integrated counts "
            f"({integrated_counts}) -- looks like it reduced the whole observation"
        )
        assert counts < integrated_counts * 0.98, (
            f"segment {i} counts ({counts}) are suspiciously close to integrated "
            f"counts ({integrated_counts}); a narrower time-filtered selection "
            f"should carry a materially smaller share of the total events"
        )

    assert sum(seg_counts.values()) == pytest.approx(integrated_counts, rel=0.02)


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
        provenance_line = f.readline()
        rows = list(csv.DictReader(f))
    assert provenance_line.startswith("# pcube_timeseg provenance: ")
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
# GTI pre-reduction gate: validated against real GTI structure
# --------------------------------------------------------------------------- #

def test_empty_slice_gate_default_separates_known_pathological_segments(tmp_path):
    """At the production default segment_duration_s=5000s, real GTI structure
    puts exactly one near-zero-exposure segment in each fixture -- a real GTI
    overlap of 1.7s (DEADC-scaled: 1.68s) at 01250401 segment 7, and 241.0s
    (DEADC-scaled: 240.18s) at 02001901 segment 148 (both hand-flagged by
    prior worker runs) -- squarely below the default min_seg_exposure_s,
    while every other real segment in both observations clears it. This
    grounds the 1000s default in the actual data rather than merely
    asserting it (computed directly from GTI/DEADC headers -- no
    xpselect/xpbin call, so it stays fast even at the ~170-segment
    production default)."""
    default = inspect.signature(red_ixpe.pcube_timeseg).parameters["min_seg_exposure_s"].default

    # known_exposure is the DEADC-scaled good exposure pcube_timeseg's gate
    # actually compares to min_seg_exposure_s (raw GTI overlap * DEADC).
    cases = [(OBSID, "v03", 7, 1.6826600471973419), (OBSID2, "v01", 148, 240.1806)]
    for obsid, version, pathological_idx, known_exposure in cases:
        file_list = _copy_real_event_files(str(tmp_path), obsid, version)
        exposures = _real_segment_exposures(file_list[0], 5000.0)

        assert exposures[pathological_idx] == pytest.approx(known_exposure, abs=0.05)
        assert exposures[pathological_idx] < default

        other_nonzero = [e for i, e in enumerate(exposures) if i != pathological_idx and e > 0]
        assert other_nonzero, f"{obsid}: expected other real segments at this duration"
        assert min(other_nonzero) >= default, (
            f"{obsid}: a genuine segment's exposure fell below min_seg_exposure_s "
            f"({min(other_nonzero)} < {default}); the default no longer separates "
            f"pathological slices from real ones for this observation"
        )


def test_gate_skips_before_reduction_and_override_replaces_stale_csv(tmp_path, monkeypatch):
    """Combines the safety (mismatched-provenance override) and capability
    (GTI gate fires before any reduction) checks on one real fixture to
    minimise real xpselect/xpbin calls: a first call at SEGMENT_DURATION_S
    (2 segments) establishes a CSV, then a second call at a different
    duration must both (a) delete and rebuild the stale CSV rather than glue
    old and new rows together, and (b) skip its own tail window -- a real
    500s-wide window with measured 494.9s of GTI overlap, just under the
    default 1000s gate -- without ever invoking xpbin for it."""
    file_list = _copy_real_event_files(str(tmp_path))
    first = red_ixpe.pcube_timeseg(
        OBSID, file_list, IRF_DATE, bkg=False, segment_duration_s=SEGMENT_DURATION_S)
    csv_path = first["csv_path"]
    with open(csv_path) as f:
        original_csv_text = f.read()

    xpbin_calls = []
    orig_xpbin = red_ixpe.pipeline.xpbin

    def _spy_xpbin(*args, **kwargs):
        xpbin_calls.append((args, kwargs))
        return orig_xpbin(*args, **kwargs)

    monkeypatch.setattr(red_ixpe.pipeline, "xpbin", _spy_xpbin)

    second = red_ixpe.pcube_timeseg(
        OBSID, file_list, IRF_DATE, bkg=False,
        segment_duration_s=GATE_SEGMENT_DURATION_S, resume=True)

    # Safety: mismatched provenance (different segment_duration_s) -> override,
    # never a silent glue of the old 2-segment binning onto the new one.
    assert second["resumed"] is False
    with open(csv_path) as f:
        new_csv_text = f.read()
    assert new_csv_text != original_csv_text
    assert {row["segment_id"] for row in second["timeseries"]} != {0, 1}

    # Capability: the gate skipped the 494.9s tail window before any
    # xpselect/xpbin call -- only 2 real reductions ran (the surviving
    # ~848000s window's source cube, plus the integrated cube), never a
    # third for the gated window.
    assert second["n_expected"] == 1
    assert second["n_skipped_empty"] == 1
    assert len(second["timeseries"]) == 1
    assert len(xpbin_calls) == 2


def test_resume_skips_completed_segments(tmp_path, monkeypatch):
    """Deleting one segment's row from an otherwise-complete, matching-
    provenance CSV and rerunning with resume=True must recompute only that
    segment -- proving resume actually skips expensive work, not just
    produces a well-shaped result."""
    file_list = _copy_real_event_files(str(tmp_path))
    first = red_ixpe.pcube_timeseg(
        OBSID, file_list, IRF_DATE, bkg=False, segment_duration_s=SEGMENT_DURATION_S)
    assert first["complete"] is True
    csv_path = first["csv_path"]

    with open(csv_path, newline="") as f:
        provenance_line = f.readline()
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        kept_rows = [r for r in reader if int(r["segment_id"]) != 1]
    with open(csv_path, "w", newline="") as f:
        f.write(provenance_line)
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(kept_rows)

    xpbin_calls = []
    orig_xpbin = red_ixpe.pipeline.xpbin

    def _spy_xpbin(*args, **kwargs):
        xpbin_calls.append((args, kwargs))
        return orig_xpbin(*args, **kwargs)

    monkeypatch.setattr(red_ixpe.pipeline, "xpbin", _spy_xpbin)

    second = red_ixpe.pcube_timeseg(
        OBSID, file_list, IRF_DATE, bkg=False,
        segment_duration_s=SEGMENT_DURATION_S, resume=True)

    assert second["resumed"] is True
    assert second["complete"] is True
    assert len(second["timeseries"]) == 2
    assert {row["segment_id"] for row in second["timeseries"]} == {0, 1}
    # Exactly one reduction (segment 1's source cube) ran; segment 0 and the
    # integrated row were re-emitted from disk, and the integrated row was
    # already present too, so xpbin fires only once.
    assert len(xpbin_calls) == 1


# --------------------------------------------------------------------------- #
# Skill wrapper: file discovery, CSV placement, catalog recording, bounded
# return contract
# --------------------------------------------------------------------------- #

def test_skill_wrapper_records_catalog_bounded_return_and_falls_back_to_plain_evt2(tmp_path):
    """01250401's real files have no '_rej' variant (mirrors 03002599) -- the
    skill wrapper's glob must fall back to the plain evt2_v*.fits pattern.
    Confirms end-to-end: file discovery, reduction, CSV on disk, level-3
    catalog recording with the completeness flag, and that the wrapper's
    return dict is bounded (no full-length 'timeseries' key) regardless --
    the differential between the bounded transport dict and the complete
    on-disk CSV pins the transport-limit defect."""
    obs_db_path = str(tmp_path)
    _copy_real_event_files(obs_db_path)

    result = ixpe_session.ixpe_pcube_timeseg(
        obs_db_path, OBSID, bkg=False,
        segment_duration_s=SEGMENT_DURATION_S)

    assert result["exit_code"] == 0, result.get("stderr")
    assert os.path.exists(result["csv_path"])
    assert result["n_segments"] == 2
    assert result["n_expected"] == 2
    assert result["n_skipped_empty"] == 0
    assert result["complete"] is True
    assert result["resumed"] is False
    assert "timeseries" not in result
    assert len(result["head"]) == 2
    assert len(result["tail"]) == 2
    assert result["integrated"]["segment_id"] == -1

    with open(result["csv_path"], newline="") as f:
        f.readline()  # provenance header
        csv_rows = list(csv.DictReader(f))
    assert len(csv_rows) == result["n_segments"] + 1  # + integrated

    entry = neutral.read(obs_db_path)["entries"][OBSID]
    assert entry["max_level"] == 3
    (prod,) = [p for p in entry["products"] if p["kind"] == "pcube_timeseg"]
    assert prod["path"] == f"{OBSID}/pcube_timeseg_{OBSID}.csv"
    assert prod["produced_by"] == "xspec_utils.red_ixpe.pcube_timeseg"
    assert prod["params"]["bkg"] is False
    assert prod["params"]["segment_duration_s"] == SEGMENT_DURATION_S
    assert prod["params"]["complete"] is True


def test_wrapper_return_is_bounded_independent_of_segment_count(tmp_path):
    """The wrapper's return dict must stay well under the sandbox transport
    chunk limit regardless of how many segments the reduction produced --
    this is what breaks on 01250401/02001901 at the ~150-170 segment
    production default. Uses GATE_SEGMENT_DURATION_S (1 real segment) to
    keep the call itself cheap while asserting the bound is structural
    (bounded head/tail, no full array), not merely small because few
    segments happened to run."""
    obs_db_path = str(tmp_path)
    _copy_real_event_files(obs_db_path)

    result = ixpe_session.ixpe_pcube_timeseg(
        obs_db_path, OBSID, bkg=False,
        segment_duration_s=GATE_SEGMENT_DURATION_S)

    assert result["exit_code"] == 0, result.get("stderr")
    assert "timeseries" not in result
    assert all(len(v) <= 3 for v in result.values() if isinstance(v, list))

    serialized = json.dumps(result)
    assert len(serialized.encode("utf-8")) < 16000

    with open(result["csv_path"], newline="") as f:
        f.readline()
        csv_rows = list(csv.DictReader(f))
    assert len(csv_rows) == result["n_segments"] + 1


def test_wrapper_reinvocation_ignores_generated_products(tmp_path):
    """Regression: the wrapper's fallback glob must match only original
    Level-2 event files, never this tool's own generated products
    (..._seg000src.fits, ..._seg000src_pcube.fits, ..._intsrc.fits, ...).
    Without the v[0-9][0-9] anchor, a second run on a no-'_rej' ObsID like
    01250401 rediscovers those generated files as inputs and
    _validate_du_files then sees multiple files per DU and fails. The
    second call also exercises resume end-to-end through the wrapper: same
    arguments -> same provenance -> resumed=True, no new reduction work."""
    obs_db_path = str(tmp_path)
    _copy_real_event_files(obs_db_path)

    first = ixpe_session.ixpe_pcube_timeseg(
        obs_db_path, OBSID, bkg=False,
        segment_duration_s=SEGMENT_DURATION_S)
    assert first["exit_code"] == 0, first.get("stderr")
    assert first["n_segments"] == 2

    # Direct discovery check: the generated _seg*/_pcube files now sitting in
    # event_l2/ (written by the first run) must NOT be matched by the same
    # strict pattern the wrapper uses.
    evt_dir = os.path.join(obs_db_path, OBSID, "event_l2")
    strict_glob = sorted(glob.glob(os.path.join(evt_dir, f"ixpe{OBSID}_det[123]_evt2_v[0-9][0-9].fits")))
    assert len(strict_glob) == 3, (
        f"expected exactly 3 original Level-2 files (one per DU), got {strict_glob}"
    )

    second = ixpe_session.ixpe_pcube_timeseg(
        obs_db_path, OBSID, bkg=False,
        segment_duration_s=SEGMENT_DURATION_S)
    assert second["exit_code"] == 0, second.get("stderr")
    assert second["n_segments"] == 2
    assert second["resumed"] is True


def test_skill_wrapper_missing_event_files_fails_loudly_no_record(tmp_path):
    result = ixpe_session.ixpe_pcube_timeseg(str(tmp_path), "99999999")
    assert result["exit_code"] == 1
    assert "stderr" in result
    assert neutral.read(str(tmp_path))["entries"] == {}


def test_repro_integrated_within_tol(tmp_path_factory):
    """Reproducibility of the wrapper's bounded return: two independent
    wrapper calls on independent scratch copies with identical arguments
    reproduce the integrated PD/PA and the surviving segment's QI/UI within
    1e-6 relative tolerance (the underlying reduction is a deterministic
    read of real, identical input data -- verified empirically bit-exact
    across independent copies). Uses GATE_SEGMENT_DURATION_S to keep this to
    one real segment reduction (+ integrated) per call."""
    results = []
    for _ in range(2):
        obs_db_path = str(tmp_path_factory.mktemp("pcube_timeseg_wrapper_repro"))
        _copy_real_event_files(obs_db_path)
        results.append(ixpe_session.ixpe_pcube_timeseg(
            obs_db_path, OBSID, bkg=False,
            segment_duration_s=GATE_SEGMENT_DURATION_S))

    a, b = results
    assert a["exit_code"] == 0 and b["exit_code"] == 0
    assert a["n_segments"] == b["n_segments"] == 1

    for field in ("PD_pct", "PA_deg"):
        assert a["integrated"][field] == pytest.approx(b["integrated"][field], rel=1e-6)
    for field in ("QI", "UI"):
        assert a["integrated"][field] == pytest.approx(b["integrated"][field], rel=1e-6)
        assert a["head"][0][field] == pytest.approx(b["head"][0][field], rel=1e-6)
