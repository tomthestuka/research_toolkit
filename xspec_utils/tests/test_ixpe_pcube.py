"""Tests for xspec_utils' ixpe_pcube (lib/red_ixpe.py's extended pcube + skills/session.py).

Uses the real, already-local Level-2 event files for IXPE ObsIDs 01002901 (Cyg
X-1, May 2022) and 01250101 (Cyg X-1, June 2022) -- the Krawczynski+2022
energy-resolved-polarimetry / native multi-ObsID-combination brief's own
example pair -- copied into isolated tmp obs_db_path directories per test so
reductions run on scratch copies and never touch the shared analyst-env data
directory or its catalog.json.

Also covers a latent _get_obs_file_list defect this brief's own procedure
(red_ixpe.get_all_file_list -> filter to the requested ObsID(s)) exposed: a
raw substring match on the mission-wide file list pulls in
'01002901_picorr' -- a differently-processed variant whose filenames still
embed the same 9-digit ObsID -- whenever a plain '01002901' is requested,
making _validate_du_files see two files per DU and raise. Fixed by matching
each file's own observation directory instead of a substring test.
"""

from __future__ import annotations

import glob
import math
import os
import shutil

import numpy as np
import pytest

from obsdb_catalog import neutral
from xspec_utils.lib import red_ixpe
from xspec_utils.skills import session as ixpe_session
from xspec_utils.skills import TOOL_MANIFEST

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
OBSID_A = "01002901"   # Cyg X-1, May 2022
OBSID_B = "01250101"   # Cyg X-1, June 2022
VERSION = "v01"

# Cyg X-1 centroid (lib/src_pos.py's 01002901 entry) -- one shared region
# centre for both ObsIDs, as a caller supplies for a combined reduction.
RA, DEC = 299.5971468, 35.1922507

# Both real ObsIDs' own midpoint MJD resolves to this installed epoch
# (verified directly against irf_epoch.resolve_irf_date); pinned explicitly
# here so every reduction in this file -- solo or combined -- uses the
# identical epoch, independent of the auto-resolution logic under test
# elsewhere (test_irf_epoch.py).
IRF_DATE = "20220702"

DEFAULT_EBINNING = [2.0, 3.0, 4.0, 6.0, 8.0]


def _copy_real_event_files(dest_dir: str, obsid: str, dus=(1, 2, 3)) -> list[str]:
    """Copy the real, already-local Level-2 _rej event files for `obsid`
    (DUs in `dus`) into `<dest_dir>/<obsid>/event_l2/` and return the sorted
    list of copied paths."""
    src_evt_dir = os.path.join(REAL_IXPE_DIR, obsid, "event_l2")
    dst_evt_dir = os.path.join(dest_dir, obsid, "event_l2")
    os.makedirs(dst_evt_dir, exist_ok=True)
    file_list = []
    for det in dus:
        src = os.path.join(src_evt_dir, f"ixpe{obsid}_det{det}_evt2_{VERSION}_rej.fits")
        if not os.path.exists(src):
            pytest.skip(f"real fixture ObsID {obsid} not local at {src}")
        dst = os.path.join(dst_evt_dir, os.path.basename(src))
        shutil.copyfile(src, dst)
        file_list.append(dst)
    return sorted(file_list)


def _cube_total_counts(cube) -> float:
    """Sum total COUNTS across every energy band straight off the returned,
    background-subtracted xBinnedPolarizationCube -- the analyzed product the
    tool actually returns, not on-disk intermediate files."""
    return float(np.sum(cube.COUNTS))


# --------------------------------------------------------------------------- #
# Regression: get_all_file_list's mission-wide file_list must not collide
# plain '01002901' with the '01002901_picorr' variant (read-only; the real
# mission directory is never written to).
# --------------------------------------------------------------------------- #

def test_get_obs_file_list_excludes_picorr_variant_on_full_mission_tree():
    if not os.path.isdir(os.path.join(REAL_IXPE_DIR, f"{OBSID_A}_picorr")):
        pytest.skip(f"{OBSID_A}_picorr fixture not local")
    file_list = red_ixpe.get_all_file_list(REAL_IXPE_DIR)
    obs_files = red_ixpe._get_obs_file_list(OBSID_A, file_list)
    assert len(obs_files) == 3
    expected_dir = os.path.join(REAL_IXPE_DIR, OBSID_A, "event_l2")
    assert all(os.path.dirname(f) == expected_dir for f in obs_files)
    # Must not raise: exactly one file per DU once the picorr variant is excluded.
    red_ixpe._validate_du_files(OBSID_A, obs_files, caller="test")


# --------------------------------------------------------------------------- #
# Fast, pure-logic / silent-misuse guards -- must fail before any real
# xpselect/xpbin call.
# --------------------------------------------------------------------------- #

def test_provenance_csv_round_trip(tmp_path):
    csv_path = str(tmp_path / "ixpe_pcube_roundtrip.csv")
    provenance = {"obsids": [OBSID_A], "ebinning": DEFAULT_EBINNING, "src_rad": 1.0,
                  "bkg_inner_rad": 3.0, "bkg_outer_rad": 5.0, "ra": RA, "dec": DEC,
                  "irf_date": IRF_DATE, "combined": False}
    bands = [
        {"E_LO": 2.0, "E_HI": 3.0, "PD": 0.05, "PD_err": 0.01, "PA": 10.0, "PA_err": 5.0,
         "QN": 0.03, "QN_err": 0.01, "UN": 0.04, "UN_err": 0.01, "MDP99": 0.06, "SIGNIF": 4.2},
        {"E_LO": 3.0, "E_HI": 4.0, "PD": 0.07, "PD_err": 0.01, "PA": 12.0, "PA_err": 4.0,
         "QN": 0.05, "QN_err": 0.01, "UN": 0.06, "UN_err": 0.01, "MDP99": 0.05, "SIGNIF": 5.1},
    ]
    red_ixpe.write_pcube_csv(csv_path, provenance, bands)

    read_provenance, read_bands = red_ixpe.read_pcube_csv(csv_path)
    assert read_provenance == provenance
    assert len(read_bands) == 2
    for expected, got in zip(bands, read_bands):
        for key in expected:
            assert got[key] == pytest.approx(expected[key])


def test_read_pcube_csv_missing_header_raises(tmp_path):
    csv_path = str(tmp_path / "no_header.csv")
    with open(csv_path, "w") as f:
        f.write("band_id,E_LO,E_HI\n0,2.0,3.0\n")
    with pytest.raises(ValueError, match="provenance"):
        red_ixpe.read_pcube_csv(csv_path)


def test_pcube_raises_on_empty_obs_list(tmp_path):
    file_list = _copy_real_event_files(str(tmp_path), OBSID_A)
    with pytest.raises(ValueError, match="empty list"):
        red_ixpe.pcube([], file_list, ebinning=DEFAULT_EBINNING, src_rad=1.0,
                        bkg_inner_rad=3.0, bkg_outer_rad=5.0, ra=RA, dec=DEC,
                        irf_date=IRF_DATE)


def test_pcube_multi_obs_changes_false_raises(tmp_path):
    file_list = (_copy_real_event_files(str(tmp_path), OBSID_A)
                 + _copy_real_event_files(str(tmp_path), OBSID_B))
    with pytest.raises(ValueError, match="multi-ObsID"):
        red_ixpe.pcube([OBSID_A, OBSID_B], file_list, changes=False)


def test_wrapper_missing_obsid_fails_loudly_no_record(tmp_path):
    obs_db_path = str(tmp_path)
    _copy_real_event_files(obs_db_path, OBSID_A)
    result = ixpe_session.ixpe_pcube(obs_db_path, "99999999", ra=RA, dec=DEC,
                                      irf_date=IRF_DATE)
    assert result["exit_code"] == 1
    assert "stderr" in result
    assert neutral.read(obs_db_path)["entries"] == {}


def test_wrapper_missing_du_fails_loudly_no_record(tmp_path):
    obs_db_path = str(tmp_path)
    _copy_real_event_files(obs_db_path, OBSID_A, dus=(1, 2))
    result = ixpe_session.ixpe_pcube(obs_db_path, OBSID_A, ra=RA, dec=DEC,
                                      irf_date=IRF_DATE)
    assert result["exit_code"] == 1
    assert "DU3" in result["stderr"]
    assert neutral.read(obs_db_path)["entries"] == {}


def test_wrapper_ebinning_too_short_fails_loudly_no_record(tmp_path):
    obs_db_path = str(tmp_path)
    _copy_real_event_files(obs_db_path, OBSID_A)
    result = ixpe_session.ixpe_pcube(obs_db_path, OBSID_A, ra=RA, dec=DEC,
                                      ebinning=[5.0], irf_date=IRF_DATE)
    assert result["exit_code"] == 1
    assert "ebinning" in result["stderr"]
    assert neutral.read(obs_db_path)["entries"] == {}
    # No partial CSV left behind either.
    assert not glob.glob(os.path.join(obs_db_path, OBSID_A, "ixpe_pcube_*.csv"))


def test_wrapper_ebinning_empty_list_fails_loudly_no_record(tmp_path):
    obs_db_path = str(tmp_path)
    _copy_real_event_files(obs_db_path, OBSID_A)
    result = ixpe_session.ixpe_pcube(obs_db_path, OBSID_A, ra=RA, dec=DEC,
                                      ebinning=[], irf_date=IRF_DATE)
    assert result["exit_code"] == 1
    assert "ebinning" in result["stderr"]
    assert neutral.read(obs_db_path)["entries"] == {}
    assert not glob.glob(os.path.join(obs_db_path, OBSID_A, "ixpe_pcube_*.csv"))


# --------------------------------------------------------------------------- #
# Manifest wiring (AI-facing half)
# --------------------------------------------------------------------------- #

def test_tool_manifest_wiring():
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.ixpe_pcube"]
    assert entry["module"] == "xspec_utils.skills.session"
    assert entry["fn"] == "ixpe_pcube"
    props = entry["args_schema"]["properties"]
    for name in ("obs_db_path", "obs", "ra", "dec", "ebinning", "src_rad",
                 "bkg_inner_rad", "bkg_outer_rad"):
        assert name in props, f"{name} missing from args_schema.properties"
    assert entry["args_schema"]["required"] == ["obs_db_path", "obs", "ra", "dec"]
    # The IRF epoch is auto-selected inside the tool -- not an AI-facing knob.
    assert "irf_date" not in props
    # 'obs' must accept either a bare ObsID string or a list of ObsID strings.
    obs_types = {opt.get("type") for opt in props["obs"]["oneOf"]}
    assert obs_types == {"string", "array"}
    assert props["ebinning"]["minItems"] == 2
    assert "2,3,4,6,8" in props["ebinning"]["description"].replace(" ", "")
    # The manifest documents the units contract and the combined-ObsID mode. The one-line
    # brief carries discovery keywords while the prose lives in args_schema["description"]
    # (deferred-loading split), so assert over the whole manifest text, not just the brief.
    manifest_text = f"{entry['description']} {entry['args_schema'].get('description', '')}"
    assert "fraction" in manifest_text.lower()
    assert "combined" in manifest_text


# --------------------------------------------------------------------------- #
# Real-data reduction: combined multi-ObsID cube (default ebinning)
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def combined_result(tmp_path_factory):
    """One real combined [OBSID_A, OBSID_B] ixpe_pcube call, default ebinning
    (4 bands). Shared read-only across the tests below to avoid re-running
    the real ixpeobssim pipeline more than once."""
    tmp_dir = tmp_path_factory.mktemp("ixpe_pcube_combined")
    obs_db_path = str(tmp_dir)
    _copy_real_event_files(obs_db_path, OBSID_A)
    _copy_real_event_files(obs_db_path, OBSID_B)
    result = ixpe_session.ixpe_pcube(obs_db_path, [OBSID_A, OBSID_B], ra=RA, dec=DEC,
                                      irf_date=IRF_DATE)
    assert result["exit_code"] == 0, result.get("stderr")
    return result, obs_db_path


def test_combined_default_edges_multi_band_increasing(combined_result):
    result, _ = combined_result
    bands = result["bands"]
    assert len(bands) == 4
    e_lo = [b["E_LO"] for b in bands]
    e_hi = [b["E_HI"] for b in bands]
    assert e_lo == sorted(e_lo) and len(set(e_lo)) == len(e_lo)
    assert e_hi == sorted(e_hi) and len(set(e_hi)) == len(e_hi)
    assert e_lo == DEFAULT_EBINNING[:-1]
    assert e_hi == DEFAULT_EBINNING[1:]
    assert result["provenance"]["combined"] is True
    assert result["provenance"]["obsids"] == sorted([OBSID_A, OBSID_B])


def test_combined_cube_counts_sum_and_exceed_solo(tmp_path_factory):
    """Capability: the natively-combined cube red_ixpe.pcube *returns* really
    carries both ObsIDs' events -- its total COUNTS must (a) strictly exceed
    either single-ObsID cube's and (b) equal their sum within tolerance.
    Compares the RETURNED, analyzed cube objects (not intermediate FITS on
    disk), so a bug that reduced only one ObsID into the returned cube while
    still writing both ObsIDs' src_pcube.fits would fail here."""
    # Three independent scratch reductions on identical region/ebinning/epoch.
    dir_a = str(tmp_path_factory.mktemp("solo_a"))
    fa = _copy_real_event_files(dir_a, OBSID_A)
    src_a, _, _ = red_ixpe.pcube(OBSID_A, fa, ebinning=DEFAULT_EBINNING,
                                 src_rad=1.0, bkg_inner_rad=3.0, bkg_outer_rad=5.0,
                                 ra=RA, dec=DEC, irf_date=IRF_DATE)

    dir_b = str(tmp_path_factory.mktemp("solo_b"))
    fb = _copy_real_event_files(dir_b, OBSID_B)
    src_b, _, _ = red_ixpe.pcube(OBSID_B, fb, ebinning=DEFAULT_EBINNING,
                                 src_rad=1.0, bkg_inner_rad=3.0, bkg_outer_rad=5.0,
                                 ra=RA, dec=DEC, irf_date=IRF_DATE)

    dir_c = str(tmp_path_factory.mktemp("combined"))
    fc = (_copy_real_event_files(dir_c, OBSID_A)
          + _copy_real_event_files(dir_c, OBSID_B))
    src_c, _, _ = red_ixpe.pcube([OBSID_A, OBSID_B], fc, ebinning=DEFAULT_EBINNING,
                                 src_rad=1.0, bkg_inner_rad=3.0, bkg_outer_rad=5.0,
                                 ra=RA, dec=DEC, irf_date=IRF_DATE)

    counts_a = _cube_total_counts(src_a)
    counts_b = _cube_total_counts(src_b)
    counts_c = _cube_total_counts(src_c)
    print(f"counts_a={counts_a} counts_b={counts_b} counts_c={counts_c}")

    assert counts_a > 0 and counts_b > 0
    # The combined cube used strictly more events than either alone, by a real
    # margin (not a rounding-level nudge): it must exceed the larger solo cube.
    assert counts_c > max(counts_a, counts_b) * 1.3
    # And it is the native sum of both, within tolerance (linear src-bkg
    # subtraction is ~additive for a bright source; small bkg-scaling slack).
    assert counts_c == pytest.approx(counts_a + counts_b, rel=0.05)


def test_combined_catalog_records_both_obsids(combined_result):
    result, obs_db_path = combined_result
    entries = neutral.read(obs_db_path)["entries"]
    for obsid in (OBSID_A, OBSID_B):
        assert entries[obsid]["max_level"] == 3
        (prod,) = [p for p in entries[obsid]["products"] if p["kind"] == "pcube"]
        assert prod["path"] == os.path.relpath(result["csv_path"], obs_db_path)
        assert prod["produced_by"] == "xspec_utils.red_ixpe.pcube"
        assert prod["params"]["combined"] is True
        assert prod["params"]["obsids"] == sorted([OBSID_A, OBSID_B])
    # Combined product lands at the mission root, not nested under either ObsID.
    assert os.path.dirname(result["csv_path"]) == os.path.normpath(obs_db_path)


def test_differential_pd_pa_from_qn_un(combined_result):
    """Independent recomputation: PD=sqrt(QN^2+UN^2), PA=0.5*atan2(UN,QN) from
    the tool's own returned QN/UN -- never referencing the tool's own PD/PA --
    must match the tool's returned PD/PA for every band."""
    result, _ = combined_result
    for band in result["bands"]:
        pd_check = math.hypot(band["QN"], band["UN"])
        pa_check_deg = math.degrees(0.5 * math.atan2(band["UN"], band["QN"]))
        assert pd_check == pytest.approx(band["PD"], rel=0.01, abs=1e-4)
        assert pa_check_deg == pytest.approx(band["PA"], rel=0.01, abs=0.5)


# --------------------------------------------------------------------------- #
# Real-data reduction: single ObsID, single-edge-pair -> single broadband row
# --------------------------------------------------------------------------- #

def test_single_obs_single_edge_pair_yields_one_broadband_row(tmp_path):
    obs_db_path = str(tmp_path)
    _copy_real_event_files(obs_db_path, OBSID_A)
    result = ixpe_session.ixpe_pcube(obs_db_path, OBSID_A, ra=RA, dec=DEC,
                                      ebinning=[2.0, 8.0], irf_date=IRF_DATE)
    assert result["exit_code"] == 0, result.get("stderr")
    assert len(result["bands"]) == 1
    band = result["bands"][0]
    assert band["E_LO"] == pytest.approx(2.0)
    assert band["E_HI"] == pytest.approx(8.0)
    assert result["provenance"]["combined"] is False

    entries = neutral.read(obs_db_path)["entries"]
    assert set(entries) == {OBSID_A}
    assert entries[OBSID_A]["max_level"] == 3
    assert os.path.dirname(result["csv_path"]) == os.path.join(obs_db_path, OBSID_A)

    read_provenance, read_bands = red_ixpe.read_pcube_csv(result["csv_path"])
    assert read_provenance == result["provenance"]
    assert len(read_bands) == 1
    assert read_bands[0]["PD"] == pytest.approx(band["PD"])
    assert read_bands[0]["PA"] == pytest.approx(band["PA"])


# --------------------------------------------------------------------------- #
# Reproducibility
# --------------------------------------------------------------------------- #

def test_repro_ixpe_pcube_stable_across_calls(tmp_path_factory):
    """Two independent ixpe_pcube calls on independent scratch copies of
    01002901, identical arguments (default ebinning [2,3,4,6,8]), reproduce
    per-band PD and PA within tolerance -- the underlying xpselect/xpbin
    pipeline is deterministic on unchanged input."""
    results = []
    for _ in range(2):
        tmp_dir = tmp_path_factory.mktemp("ixpe_pcube_repro")
        obs_db_path = str(tmp_dir)
        _copy_real_event_files(obs_db_path, OBSID_A)
        results.append(ixpe_session.ixpe_pcube(obs_db_path, OBSID_A, ra=RA, dec=DEC,
                                                irf_date=IRF_DATE))

    a, b = results
    assert a["exit_code"] == 0 and b["exit_code"] == 0
    assert len(a["bands"]) == len(b["bands"]) == 4
    for band_a, band_b in zip(a["bands"], b["bands"]):
        assert band_a["E_LO"] == pytest.approx(band_b["E_LO"])
        assert band_a["PD"] == pytest.approx(band_b["PD"], rel=1e-3)
        assert band_a["PA"] == pytest.approx(band_b["PA"], rel=1e-3)
