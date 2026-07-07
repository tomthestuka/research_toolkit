"""Tests for xspec_utils' ixpe_HR (lib/red_ixpe.HR + skills/session.py).

Uses real, already-local per-DU IXPE PHA1 (Stokes I) src+bkg products for two
ObsIDs already on disk -- 01250401 (4U 1630-47, plain evt2_v03 files, all 3 DUs)
and 01002901 (Cyg X-1, _rej variant, all 3 DUs) -- copied into isolated tmp
directories per test so this tool's reads run on scratch copies rather than the
shared analyst-env data directory.

Covers the three defects reported against the prior red_ixpe.HR: DU3 silently
dropped by an off-by-one loop bound, HR computed as H/S instead of the
community-normalized (H-S)/(H+S), and an error propagated from misreading
XSPEC's Spectrum.rate tuple -- replaced here by reading PHA RATE/STAT_ERR
directly (no XSPEC AllData loading at all) and a closed-form Gaussian
propagation.
"""

from __future__ import annotations

import math
import os
import shutil

import numpy as np
import pytest
from astropy.io import fits

from obsdb_catalog import neutral
from xspec_utils.lib import red_ixpe
from xspec_utils.skills import TOOL_MANIFEST
from xspec_utils.skills import session as ixpe_session

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
OBSID_3DU = "01250401"   # 4U 1630-47; plain evt2_v03.fits, no _rej variant
VERSION_3DU = "v03"
OBSID_REJ = "01002901"   # Cyg X-1; _rej instrumental-background-rejected variant
VERSION_REJ = "v01"

SOFT = (2.0, 4.0)
HARD = (4.0, 8.0)


def _copy_event_and_pha(dest_dir: str, obsid: str, version: str, rej: bool,
                         dus=(1, 2, 3)) -> list[str]:
    """Copy real Level-2 event files + real per-DU src/bkg Stokes-I PHA1 products
    for `obsid` into `<dest_dir>/<obsid>/event_l2/`. Returns the sorted list of
    copied Level-2 event file paths."""
    tag = f"{version}_rej" if rej else version
    src_evt_dir = os.path.join(REAL_IXPE_DIR, obsid, "event_l2")
    dst_evt_dir = os.path.join(dest_dir, obsid, "event_l2")
    os.makedirs(dst_evt_dir, exist_ok=True)
    file_list = []
    for det in dus:
        base = f"ixpe{obsid}_det{det}_evt2_{tag}"
        evt_name = f"{base}.fits"
        src_evt = os.path.join(src_evt_dir, evt_name)
        if not os.path.exists(src_evt):
            pytest.skip(f"real fixture ObsID {obsid} not local at {src_evt}")
        dst_evt = os.path.join(dst_evt_dir, evt_name)
        shutil.copyfile(src_evt, dst_evt)
        file_list.append(dst_evt)
        for suffix in ("src_pha1", "bkg_pha1"):
            name = f"{base}_{suffix}.fits"
            src = os.path.join(src_evt_dir, name)
            if not os.path.exists(src):
                pytest.skip(f"real fixture PHA1 product not local at {src}")
            shutil.copyfile(src, os.path.join(dst_evt_dir, name))
    return sorted(file_list)


def _pha_path(dest_dir, obsid, version, rej, det, region):
    tag = f"{version}_rej" if rej else version
    return os.path.join(dest_dir, obsid, "event_l2",
                         f"ixpe{obsid}_det{det}_evt2_{tag}_{region}_pha1.fits")


# --------------------------------------------------------------------------- #
# hr_from_bands -- pure closed-form math, no I/O
# --------------------------------------------------------------------------- #

def test_hr_from_bands_matches_manual_formula():
    hr, err = red_ixpe.hr_from_bands(hard_rate=3.0, hard_err=0.2, soft_rate=1.0, soft_err=0.1)
    denom = 4.0
    expected_hr = (3.0 - 1.0) / denom
    expected_err = (2.0 / denom ** 2) * math.sqrt(1.0 ** 2 * 0.2 ** 2 + 3.0 ** 2 * 0.1 ** 2)
    assert hr == pytest.approx(expected_hr)
    assert err == pytest.approx(expected_err)


def test_hr_from_bands_nonpositive_denominator_returns_nan():
    hr, err = red_ixpe.hr_from_bands(hard_rate=-1.0, hard_err=0.1, soft_rate=1.0, soft_err=0.1)
    assert math.isnan(hr) and math.isnan(err)
    hr0, err0 = red_ixpe.hr_from_bands(0.0, 0.1, 0.0, 0.1)
    assert math.isnan(hr0) and math.isnan(err0)


# --------------------------------------------------------------------------- #
# capability-du: DU3 must contribute, and IXPE flying with a DU off is graceful
# --------------------------------------------------------------------------- #

def test_capability_all_three_dus_contribute_and_det3_matters(tmp_path):
    dest = str(tmp_path)
    file_list = _copy_event_and_pha(dest, OBSID_3DU, VERSION_3DU, rej=False)

    result_3du = red_ixpe.HR(OBSID_3DU, file_list, soft_band=SOFT, hard_band=HARD, bkg=False)
    assert result_3du["dus_used"] == [1, 2, 3]

    two_du_list = [f for f in file_list if "_det3_" not in f]
    result_2du = red_ixpe.HR(OBSID_3DU, two_du_list, soft_band=SOFT, hard_band=HARD, bkg=False)
    assert result_2du["dus_used"] == [1, 2]

    # Dropping a genuinely-present DU3 must change the summed net rate -- proves
    # the old off-by-one bug (which silently omitted the last-loaded spectrum,
    # i.e. always DU3) is fixed.
    assert result_3du["soft_rate"] != pytest.approx(result_2du["soft_rate"])
    assert result_3du["hard_rate"] != pytest.approx(result_2du["hard_rate"])
    assert result_3du["soft_rate"] > result_2du["soft_rate"]
    assert result_3du["hard_rate"] > result_2du["hard_rate"]


def test_capability_two_dus_present_returns_finite_hr_no_crash(tmp_path):
    """IXPE can genuinely fly with a DU off -- a file_list carrying only 2 DUs'
    event files must yield a finite HR using exactly those 2 DUs, not a crash."""
    dest = str(tmp_path)
    file_list = _copy_event_and_pha(dest, OBSID_3DU, VERSION_3DU, rej=False)
    two_du_list = [f for f in file_list if "_det1_" not in f]

    result = red_ixpe.HR(OBSID_3DU, two_du_list, soft_band=SOFT, hard_band=HARD, bkg=False)
    assert result["dus_used"] == [2, 3]
    assert math.isfinite(result["HR"])
    assert math.isfinite(result["HR_err"])


def test_present_du_missing_pha_product_raises_not_silently_dropped(tmp_path):
    """A DU present in file_list (its Level-2 event file exists) but whose PHA1
    product was never reduced must raise -- never be silently dropped like the
    old DU3 bug, and never treated as a graceful 'DU off'."""
    dest = str(tmp_path)
    file_list = _copy_event_and_pha(dest, OBSID_3DU, VERSION_3DU, rej=False)
    os.remove(_pha_path(dest, OBSID_3DU, VERSION_3DU, False, 3, "src"))

    with pytest.raises(ValueError, match="DU3"):
        red_ixpe.HR(OBSID_3DU, file_list, soft_band=SOFT, hard_band=HARD, bkg=False)


def test_hr_raises_on_duplicate_du_versions(tmp_path):
    """More than one Level-2 event file for the same DU (e.g. two processing
    versions, or a mixed plain+_rej set) must raise -- never silently
    double-count that DU's rate into an inflated soft_rate/hard_rate."""
    dest = str(tmp_path)
    file_list = _copy_event_and_pha(dest, OBSID_3DU, VERSION_3DU, rej=False)

    # A second, distinct event-file path for DU1 -- same _det1_evt2_ tag, a
    # different processing version -- alongside the real one in file_list.
    evt_dir = os.path.join(dest, OBSID_3DU, "event_l2")
    det1_v03 = os.path.join(evt_dir, f"ixpe{OBSID_3DU}_det1_evt2_{VERSION_3DU}.fits")
    det1_v02 = os.path.join(evt_dir, f"ixpe{OBSID_3DU}_det1_evt2_v02.fits")
    shutil.copyfile(det1_v03, det1_v02)
    file_list_with_dup = file_list + [det1_v02]

    with pytest.raises(ValueError, match="more than one Level-2 event file for DU"):
        red_ixpe.HR(OBSID_3DU, file_list_with_dup, soft_band=SOFT, hard_band=HARD, bkg=False)


# --------------------------------------------------------------------------- #
# capability-ebounds: band -> channel resolved live from the RMF, never hardcoded
# --------------------------------------------------------------------------- #

def test_capability_ebounds_channels_resolved_from_rmf_and_bands_differ(tmp_path):
    dest = str(tmp_path)
    _copy_event_and_pha(dest, OBSID_3DU, VERSION_3DU, rej=False, dus=(1,))
    pha_path = _pha_path(dest, OBSID_3DU, VERSION_3DU, False, 1, "src")

    # Independent astropy read -- not the tool's own helper.
    with fits.open(pha_path) as hdul:
        spec = hdul["SPECTRUM"]
        rmf_path = spec.header["RESPFILE"]
        rate = spec.data["RATE"].astype(float)
    with fits.open(rmf_path) as rmf_hdul:
        e_lo = rmf_hdul["EBOUNDS"].data["E_MIN"].astype(float)
        e_hi = rmf_hdul["EBOUNDS"].data["E_MAX"].astype(float)
    e_mid = 0.5 * (e_lo + e_hi)
    soft_mask_ref = (e_mid >= SOFT[0]) & (e_mid < SOFT[1])
    hard_mask_ref = (e_mid >= HARD[0]) & (e_mid < HARD[1])

    assert soft_mask_ref.sum() > 0 and hard_mask_ref.sum() > 0
    assert not np.array_equal(soft_mask_ref, hard_mask_ref)

    soft_rate_ref = float(np.nansum(rate[soft_mask_ref]))
    hard_rate_ref = float(np.nansum(rate[hard_mask_ref]))

    soft_rate_tool, _ = red_ixpe._pha_band_rate(pha_path, SOFT, 1)
    hard_rate_tool, _ = red_ixpe._pha_band_rate(pha_path, HARD, 1)

    assert soft_rate_tool == pytest.approx(soft_rate_ref, rel=1e-9)
    assert hard_rate_tool == pytest.approx(hard_rate_ref, rel=1e-9)
    assert soft_rate_tool != pytest.approx(hard_rate_tool)


def test_capability_ebounds_fallback_to_installed_du_rmf_when_respfile_stale(tmp_path):
    """A PHA1 reduced under an environment that no longer exists (RESPFILE points
    at a now-missing absolute path -- the real situation for most already-reduced
    Cyg X-1 IXPE PHA1 products on this machine) must still resolve a channel grid,
    via a DU-matched installed CALDB RMF -- the identical EBOUNDS grid, since it is
    a fixed instrumental characteristic invariant across IRF date/DU -- and must
    give the same band rates as the working, non-stale RESPFILE would."""
    dest = str(tmp_path)
    _copy_event_and_pha(dest, OBSID_3DU, VERSION_3DU, rej=False, dus=(1,))
    pha_path = _pha_path(dest, OBSID_3DU, VERSION_3DU, False, 1, "src")

    ref_soft, _ = red_ixpe._pha_band_rate(pha_path, SOFT, 1)
    ref_hard, _ = red_ixpe._pha_band_rate(pha_path, HARD, 1)

    # Simulate a migrated reduction environment -- rewrite RESPFILE to a path
    # that cannot exist, in this test's scratch copy only.
    with fits.open(pha_path, mode="update") as hdul:
        hdul["SPECTRUM"].header["RESPFILE"] = (
            "/grad/fzhou/.local/lib/python3.13/site-packages/ixpeobssim/"
            "caldb/ixpe/gpd/cpf/rmf/ixpe_d1_obssim_alpha075_v012.rmf"
        )

    fallback_rmf = red_ixpe._installed_du_rmf(1)
    with fits.open(fallback_rmf) as rmf_hdul:
        e_lo_ref = rmf_hdul["EBOUNDS"].data["E_MIN"].astype(float)
        e_hi_ref = rmf_hdul["EBOUNDS"].data["E_MAX"].astype(float)

    e_lo, e_hi = red_ixpe._resolve_rmf_ebounds(pha_path, 1)
    assert np.array_equal(e_lo, e_lo_ref)
    assert np.array_equal(e_hi, e_hi_ref)

    soft_rate, _ = red_ixpe._pha_band_rate(pha_path, SOFT, 1)
    hard_rate, _ = red_ixpe._pha_band_rate(pha_path, HARD, 1)
    assert soft_rate == pytest.approx(ref_soft, rel=1e-9)
    assert hard_rate == pytest.approx(ref_hard, rel=1e-9)


def test_ebounds_fallback_channel_count_mismatch_raises(tmp_path):
    """A fallback RMF whose EBOUNDS channel count disagrees with the PHA's own
    channel count must raise -- never silently accept a mismatched grid."""
    dest = str(tmp_path)
    _copy_event_and_pha(dest, OBSID_3DU, VERSION_3DU, rej=False, dus=(1,))
    pha_path = _pha_path(dest, OBSID_3DU, VERSION_3DU, False, 1, "src")

    with fits.open(pha_path) as hdul:
        spec_hdu = hdul["SPECTRUM"]
        truncated = fits.BinTableHDU(data=spec_hdu.data[:10], header=spec_hdu.header.copy())
        truncated.header["RESPFILE"] = "/nonexistent/env/ixpe_d1_obssim_alpha075_v012.rmf"
        new_hdul = fits.HDUList([hdul["PRIMARY"].copy(), truncated])
        new_hdul.writeto(pha_path, overwrite=True)

    with pytest.raises(ValueError, match="EBOUNDS channels"):
        red_ixpe._resolve_rmf_ebounds(pha_path, 1)


# --------------------------------------------------------------------------- #
# convention: (H-S)/(H+S), bounded [-1, +1]
# --------------------------------------------------------------------------- #

def test_convention_hr_matches_manual_formula_and_bounded(tmp_path):
    dest = str(tmp_path)
    file_list = _copy_event_and_pha(dest, OBSID_3DU, VERSION_3DU, rej=False)
    result = red_ixpe.HR(OBSID_3DU, file_list, soft_band=SOFT, hard_band=HARD, bkg=False)

    h, s = result["hard_rate"], result["soft_rate"]
    manual_hr = (h - s) / (h + s)
    assert result["HR"] == pytest.approx(manual_hr, rel=1e-9)
    assert -1.0 <= result["HR"] <= 1.0


# --------------------------------------------------------------------------- #
# differential-error: closed-form Gaussian propagation, recomputed independently
# --------------------------------------------------------------------------- #

def test_differential_hr_err_matches_closed_form_propagation(tmp_path):
    dest = str(tmp_path)
    file_list = _copy_event_and_pha(dest, OBSID_3DU, VERSION_3DU, rej=False)
    result = red_ixpe.HR(OBSID_3DU, file_list, soft_band=SOFT, hard_band=HARD, bkg=False)

    h, s = result["hard_rate"], result["soft_rate"]
    h_err, s_err = result["hard_err"], result["soft_err"]
    denom = h + s
    expected_hr_err = (2.0 / denom ** 2) * math.sqrt(s ** 2 * h_err ** 2 + h ** 2 * s_err ** 2)
    assert result["HR_err"] == pytest.approx(expected_hr_err, rel=1e-9)

    # Independently recompute soft_err/hard_err: quadrature sum of STAT_ERR over
    # the RMF-resolved in-band channels and over the 3 DUs -- never via the
    # tool's own error helper.
    total_var = {SOFT: 0.0, HARD: 0.0}
    for det in (1, 2, 3):
        pha_path = _pha_path(dest, OBSID_3DU, VERSION_3DU, False, det, "src")
        with fits.open(pha_path) as hdul:
            spec = hdul["SPECTRUM"]
            rmf_path = spec.header["RESPFILE"]
            stat_err = spec.data["STAT_ERR"].astype(float)
        with fits.open(rmf_path) as rmf_hdul:
            e_lo = rmf_hdul["EBOUNDS"].data["E_MIN"].astype(float)
            e_hi = rmf_hdul["EBOUNDS"].data["E_MAX"].astype(float)
        e_mid = 0.5 * (e_lo + e_hi)
        for band in (SOFT, HARD):
            mask = (e_mid >= band[0]) & (e_mid < band[1])
            total_var[band] += float(np.nansum(stat_err[mask] ** 2))

    assert result["soft_err"] == pytest.approx(total_var[SOFT] ** 0.5, rel=1e-9)
    assert result["hard_err"] == pytest.approx(total_var[HARD] ** 0.5, rel=1e-9)


# --------------------------------------------------------------------------- #
# capability-bkg-toggle
# --------------------------------------------------------------------------- #

def test_capability_bkg_subtraction_reduces_net_rate(tmp_path):
    """bkg=True must be a genuine subtraction, not a no-op -- proven by an
    independent (not via the tool's own _du_net_band_rate/_pha_backscal helpers)
    recomputation of the total BACKSCAL-scaled background expected to be
    subtracted per band, across all 3 DUs."""
    dest = str(tmp_path)
    file_list = _copy_event_and_pha(dest, OBSID_3DU, VERSION_3DU, rej=False)

    no_bkg = red_ixpe.HR(OBSID_3DU, file_list, soft_band=SOFT, hard_band=HARD, bkg=False)
    with_bkg = red_ixpe.HR(OBSID_3DU, file_list, soft_band=SOFT, hard_band=HARD, bkg=True)

    expected_subtracted = {SOFT: 0.0, HARD: 0.0}
    for det in (1, 2, 3):
        src_path = _pha_path(dest, OBSID_3DU, VERSION_3DU, False, det, "src")
        bkg_path = _pha_path(dest, OBSID_3DU, VERSION_3DU, False, det, "bkg")

        with fits.open(src_path) as hdul:
            spec = hdul["SPECTRUM"]
            rmf_path = spec.header["RESPFILE"]
            src_backscal = float(spec.header["BACKSCAL"])
        with fits.open(bkg_path) as hdul:
            bkg_spec = hdul["SPECTRUM"]
            bkg_rate = bkg_spec.data["RATE"].astype(float)
            bkg_backscal = float(bkg_spec.header["BACKSCAL"])
        with fits.open(rmf_path) as rmf_hdul:
            e_lo = rmf_hdul["EBOUNDS"].data["E_MIN"].astype(float)
            e_hi = rmf_hdul["EBOUNDS"].data["E_MAX"].astype(float)
        e_mid = 0.5 * (e_lo + e_hi)

        scale = src_backscal / bkg_backscal
        for band in (SOFT, HARD):
            mask = (e_mid >= band[0]) & (e_mid < band[1])
            expected_subtracted[band] += scale * float(np.nansum(bkg_rate[mask]))

    if expected_subtracted[SOFT] <= 0 and expected_subtracted[HARD] <= 0:
        pytest.skip("fixture has no in-band background counts in either band -- "
                     "cannot prove subtraction is non-trivial")

    band_rate_key = {SOFT: "soft_rate", HARD: "hard_rate"}
    for band, key in band_rate_key.items():
        if expected_subtracted[band] <= 0:
            continue
        # Strict inequality -- proves subtraction actually reduced the net rate,
        # never a no-op (equality would pass even if bkg were subtracting zero).
        assert with_bkg[key] < no_bkg[key]
        # Differential proof: the observed reduction matches the independently
        # computed BACKSCAL-scaled background exactly.
        assert no_bkg[key] - with_bkg[key] == pytest.approx(expected_subtracted[band], rel=1e-9)


# --------------------------------------------------------------------------- #
# reproducibility
# --------------------------------------------------------------------------- #

def test_repro_hr_stable_across_calls(tmp_path):
    dest = str(tmp_path)
    file_list = _copy_event_and_pha(dest, OBSID_3DU, VERSION_3DU, rej=False)

    a = red_ixpe.HR(OBSID_3DU, file_list, soft_band=SOFT, hard_band=HARD, bkg=False)
    b = red_ixpe.HR(OBSID_3DU, file_list, soft_band=SOFT, hard_band=HARD, bkg=False)

    assert a["HR"] == pytest.approx(b["HR"], abs=1e-9)
    assert a["HR_err"] == pytest.approx(b["HR_err"], abs=1e-9)
    assert a["dus_used"] == b["dus_used"]


# --------------------------------------------------------------------------- #
# Manifest wiring (AI-facing half)
# --------------------------------------------------------------------------- #

def test_tool_manifest_exposes_class_parameters():
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.ixpe_HR"]
    assert entry["module"] == "xspec_utils.skills.session"
    assert entry["fn"] == "ixpe_HR"
    props = entry["args_schema"]["properties"]
    for name in ("obs_db_path", "obsids", "soft_band", "hard_band", "bkg"):
        assert name in props, f"{name} missing from args_schema.properties"
    assert entry["args_schema"]["required"] == ["obs_db_path", "obsids"]
    assert "(H-S)/(H+S)" in entry["description"]
    assert "DU3" in entry["description"]


# --------------------------------------------------------------------------- #
# Skill wrapper: survey CSV, catalog recording, graceful partial/total failure
# --------------------------------------------------------------------------- #

def test_skill_wrapper_computes_and_records_catalog(tmp_path):
    obs_db_path = str(tmp_path)
    _copy_event_and_pha(obs_db_path, OBSID_3DU, VERSION_3DU, rej=False)

    result = ixpe_session.ixpe_HR(obs_db_path, [OBSID_3DU])

    assert result["exit_code"] == 0
    assert result["errors"] == {}
    assert os.path.exists(result["csv_path"])
    assert len(result["rows"]) == 1
    row = result["rows"][0]
    assert row["obsid"] == OBSID_3DU
    assert -1.0 <= row["HR"] <= 1.0

    with open(result["csv_path"]) as f:
        header = f.readline().strip().split(",")
    assert header == ["obsid", "soft_rate", "soft_err", "hard_rate", "hard_err", "HR", "HR_err"]

    entry = neutral.read(obs_db_path)["entries"][OBSID_3DU]
    assert entry["max_level"] == 3
    hr_products = [p for p in entry["products"] if p["kind"] == "hardness_ratio"]
    assert len(hr_products) == 1
    assert hr_products[0]["produced_by"] == "xspec_utils.red_ixpe.HR"
    assert hr_products[0]["params"]["dus_used"] == [1, 2, 3]
    assert not os.path.isabs(hr_products[0]["path"])
    assert os.path.exists(os.path.join(obs_db_path, hr_products[0]["path"]))


def test_skill_wrapper_rej_variant_obsid_computed(tmp_path):
    obs_db_path = str(tmp_path)
    _copy_event_and_pha(obs_db_path, OBSID_REJ, VERSION_REJ, rej=True)

    result = ixpe_session.ixpe_HR(obs_db_path, [OBSID_REJ], bkg=True)

    assert result["exit_code"] == 0, result.get("stderr")
    assert result["errors"] == {}
    assert len(result["rows"]) == 1


def test_skill_wrapper_glob_duplicate_du_versions_recorded_as_error(tmp_path):
    """A directory carrying two processing versions of the same DU's Level-2
    event file (OBSID_3DU has no _rej variant, so the wrapper's glob falls
    back to plain evt2_v*.fits and matches both versions for DU1) must be
    reported as this ObsID's error via HR's duplicate-DU guard -- never a
    silently inflated row, and never a crash of the rest of the survey."""
    obs_db_path = str(tmp_path)
    _copy_event_and_pha(obs_db_path, OBSID_3DU, VERSION_3DU, rej=False)
    _copy_event_and_pha(obs_db_path, OBSID_REJ, VERSION_REJ, rej=True)

    evt_dir = os.path.join(obs_db_path, OBSID_3DU, "event_l2")
    det1_v03 = os.path.join(evt_dir, f"ixpe{OBSID_3DU}_det1_evt2_{VERSION_3DU}.fits")
    det1_v02 = os.path.join(evt_dir, f"ixpe{OBSID_3DU}_det1_evt2_v02.fits")
    shutil.copyfile(det1_v03, det1_v02)

    result = ixpe_session.ixpe_HR(obs_db_path, [OBSID_3DU, OBSID_REJ])

    assert result["exit_code"] == 0, result.get("stderr")
    assert OBSID_3DU in result["errors"]
    assert "DU" in result["errors"][OBSID_3DU]
    assert len(result["rows"]) == 1
    assert result["rows"][0]["obsid"] == OBSID_REJ
    assert all(row["obsid"] != OBSID_3DU for row in result["rows"])


def test_skill_wrapper_partial_failure_records_error_no_dummy_row(tmp_path):
    obs_db_path = str(tmp_path)
    _copy_event_and_pha(obs_db_path, OBSID_3DU, VERSION_3DU, rej=False)

    result = ixpe_session.ixpe_HR(obs_db_path, [OBSID_3DU, "99999999"])

    assert result["exit_code"] == 0
    assert "99999999" in result["errors"]
    assert len(result["rows"]) == 1
    assert all(row["obsid"] != "99999999" for row in result["rows"])

    with open(result["csv_path"]) as f:
        lines = f.readlines()
    assert len(lines) == 2  # header + exactly one real row, no dummy row

    assert neutral.read(obs_db_path)["entries"].get("99999999") is None


def test_skill_wrapper_missing_pha_reports_error_no_pha_produced(tmp_path):
    """A present-but-unreduced DU inside a batch call is reported as this
    ObsID's error (via red_ixpe.HR's ValueError) rather than crashing the
    whole survey or inventing a row for it."""
    obs_db_path = str(tmp_path)
    _copy_event_and_pha(obs_db_path, OBSID_3DU, VERSION_3DU, rej=False)
    os.remove(_pha_path(obs_db_path, OBSID_3DU, VERSION_3DU, False, 2, "src"))

    result = ixpe_session.ixpe_HR(obs_db_path, [OBSID_3DU])

    assert result["exit_code"] == 1
    assert "stderr" in result
    assert not os.path.exists(os.path.join(obs_db_path, "bhxrb_hardness_ratios.csv"))


def test_skill_wrapper_all_obsids_missing_fails_loudly_no_csv(tmp_path):
    obs_db_path = str(tmp_path)
    result = ixpe_session.ixpe_HR(obs_db_path, ["99999999", "88888888"])

    assert result["exit_code"] == 1
    assert "stderr" in result
    assert not os.path.exists(os.path.join(obs_db_path, "bhxrb_hardness_ratios.csv"))
    assert neutral.read(obs_db_path)["entries"] == {}
