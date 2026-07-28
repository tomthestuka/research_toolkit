"""Tests for xspec_utils' fit_tbabs_pollin_powerlaw (lib/xspec_help.py) --
the paper's actual linear energy-dependence route for Cyg X-1 IXPE
polarization (tbabs*pollin*powerlaw), as opposed to worker 0's inline
runfit.py scipy fit on ixpe_pcube's PD(E) points (a data-analysis boundary
violation the brief flags as not toolkit-grounded).

Runs a real PyXSPEC fit against the already-local IXPE Level-2 event files
for Cyg X-1 ObsID 01002901, via a fresh ixpeobssim reduction
(ixpe_changes=True) at the paper's region radii -- copied into an isolated
tmp directory per test-module run so this never mutates the shared
analyst-env data directory. Mirrors test_fit_tbabs_polconst_powerlaw.py's
fixture/test shape exactly (same load path, same ObsID/region).
"""

from __future__ import annotations

import os
import re
import shutil

import numpy as np
import pytest
from xspec import AllData, AllModels, Fit, Model

from xspec_utils.lib import red_ixpe, xspec_help
from xspec_utils.skills import TOOL_MANIFEST

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
OBSID = "01002901"
VERSION = "v01"
SRC_RAD, BKG_INNER_RAD, BKG_OUTER_RAD = 1.333, 2.5, 5.167


def _copy_event_files(dest_dir: str, obsid: str = OBSID) -> None:
    src_evt_dir = os.path.join(REAL_IXPE_DIR, obsid, "event_l2")
    dst_evt_dir = os.path.join(dest_dir, obsid, "event_l2")
    os.makedirs(dst_evt_dir, exist_ok=True)
    for det in (1, 2, 3):
        src = os.path.join(src_evt_dir, f"ixpe{obsid}_det{det}_evt2_{VERSION}_rej.fits")
        if not os.path.exists(src):
            pytest.skip(f"real fixture ObsID {obsid} not local at {src}")
        shutil.copyfile(src, os.path.join(dst_evt_dir, os.path.basename(src)))


@pytest.fixture(scope="module")
def fresh_ixpe_load(tmp_path_factory):
    """One real fresh ixpeobssim reduction (ixpe_changes=True, paper region
    radii) + load_data call for 01002901, shared read-only across the tests
    below to avoid repeating the expensive xpselect/xpbin pipeline."""
    ixpe_dir = str(tmp_path_factory.mktemp("pollin_ixpe"))
    _copy_event_files(ixpe_dir)
    AllData.clear()
    AllModels.clear()
    indices = xspec_help.load_data(
        epoch={"IXPE": [OBSID]},
        ixpe_path=ixpe_dir,
        ixpe_changes=True,
        ixpe_src_rad=SRC_RAD,
        ixpe_bkg_inner_rad=BKG_INNER_RAD,
        ixpe_bkg_outer_rad=BKG_OUTER_RAD,
    )
    return indices, ixpe_dir


@pytest.fixture(scope="module")
def first_fit(fresh_ixpe_load):
    """The first real fit call on the fresh load -- shared across the
    output-contract, capability and differential tests below; the
    reproducibility test below performs its own, second, independent call."""
    indices, _ = fresh_ixpe_load
    return xspec_help.fit_tbabs_pollin_powerlaw(indices)


# --------------------------------------------------------------------------- #
# Edge cases -- same 1-vs->1-ObsID handling as fit_tbabs_polconst_powerlaw
# --------------------------------------------------------------------------- #

def test_guardrail_missing_ixpe_raises_cleanly():
    with pytest.raises(ValueError, match="requires indices\\['IXPE'\\]"):
        xspec_help.fit_tbabs_pollin_powerlaw({})


def test_guardrail_non_triplet_du_raises_cleanly():
    with pytest.raises(ValueError, match="not an .I, Q, U. triplet"):
        xspec_help.fit_tbabs_pollin_powerlaw({"IXPE": {OBSID: {"DU1": [1, 2]}}})


# --------------------------------------------------------------------------- #
# fit_tbabs_pollin_powerlaw -- capability, output contract, differential
# --------------------------------------------------------------------------- #

def test_capability_background_applied_and_pollin_tied(fresh_ixpe_load, first_fit):
    """(a) Background genuinely applied: net Stokes-I counts < gross I counts
    for every DU (never asserted for Q/U -- those can legitimately exceed
    gross). (b) pollin's A1/Aslope/psi1/psislope are genuinely tied: exactly
    one free parameter of each across the whole session, never one per DU --
    and in particular Aslope (the PD energy-slope term) is free, not frozen
    at 0, so the linear fit genuinely exercises the energy dependence."""
    indices, _ = fresh_ixpe_load

    for du, triplet in indices["IXPE"][OBSID].items():
        i_idx = triplet[0]
        net, _var, total, _model = AllData(i_idx).rate
        assert net < total, f"{du}: net Stokes-I rate not below gross ({net} !< {total})"

    n_free = {"A1": 0, "Aslope": 0, "psi1": 0, "psislope": 0}
    for dg in range(1, AllData.nGroups + 1):
        m = AllModels(dg)
        for pname in n_free:
            p = getattr(m.pollin, pname)
            if not p.link and not p.frozen:
                n_free[pname] += 1
    for pname, n in n_free.items():
        assert n == 1, f"expected exactly one free pollin.{pname}, got {n}"


def test_capability_slope_free_improves_over_slope_frozen(fresh_ixpe_load, first_fit):
    """The energy-dependence term was actually exercised: refitting with
    pollin's Aslope frozen at 0 (a constant-PD model) strictly improves the
    fit (chi2_frozen > chi2_free, i.e. positive Delta-chi2) relative to the
    slope-free fit reported in first_fit -- proving the linear
    energy-dependence term genuinely reduces chi2, not merely perturbs it.
    Rebuilds its own model rather than mutating the shared fixture's live
    model, so it runs independently of test execution order."""
    _indices, _ixpe_dir = fresh_ixpe_load
    chi2_free = first_fit["chi2"]

    m = Model("tbabs*pollin*powerlaw")
    m.TBabs.nH = first_fit["nH"]
    m.TBabs.nH.frozen = True
    m.pollin.A1 = first_fit["PD_A"]
    m.pollin.A1.frozen = False
    m.pollin.Aslope = 0.0
    m.pollin.Aslope.frozen = True
    m.pollin.psi1 = first_fit["PA_A"]
    m.pollin.psi1.frozen = False
    m.pollin.psislope = first_fit["PA_slope"]
    m.pollin.psislope.frozen = False
    m.powerlaw.PhoIndex = first_fit["Gamma"]
    m.powerlaw.PhoIndex.frozen = False
    Fit.perform()
    chi2_frozen = Fit.statistic

    delta_chi2 = chi2_frozen - chi2_free
    assert delta_chi2 > 1e-6, (
        f"freezing pollin.Aslope at 0 did not improve-away chi2 "
        f"(free={chi2_free:.4f}, frozen={chi2_frozen:.4f}, "
        f"delta={delta_chi2:.4f}); the energy-dependence term is a no-op")


def test_capability_output_contract(first_fit):
    """Every fitted number the brief requires must be in the returned dict --
    analysts never see XSPEC's own (muted) terminal output."""
    for key in ("PD_A", "PD_A_err", "PD_slope", "PD_slope_err",
                "PA_A", "PA_A_err", "PA_slope", "PA_slope_err",
                "Gamma", "Gamma_err", "nH", "chi2", "dof", "stat_per_dof",
                "conf_level", "indices", "stdout", "exit_code"):
        assert key in first_fit, f"{key} missing from result"

    assert first_fit["exit_code"] == 0
    assert 0.0 <= first_fit["PD_A"] <= 1.0
    assert first_fit["PD_A_err"] > 0
    assert first_fit["PD_slope_err"] > 0
    assert first_fit["PA_A_err"] > 0
    assert first_fit["PA_slope_err"] > 0
    assert first_fit["Gamma_err"] > 0
    assert first_fit["dof"] > 0
    assert first_fit["stat_per_dof"] == pytest.approx(first_fit["chi2"] / first_fit["dof"], rel=1e-9)
    assert first_fit["indices"] is not None and "IXPE" in first_fit["indices"]


def test_differential_pollin_slope_matches_independent_pcube_trend(fresh_ixpe_load, first_fit):
    """The pollin PD_slope's sign and rough magnitude must agree with an
    independent inverse-variance-weighted linear regression of ixpe_pcube's
    own PD(E) points (not with this helper's own fit code) -- the trusted
    energy-resolved PCUBE route for Cyg X-1 shows PD rising roughly
    2.6%->5.0% across 2-8 keV."""
    _indices, ixpe_dir = fresh_ixpe_load

    src_reg = os.path.join(ixpe_dir, OBSID, "src.reg")
    with open(src_reg) as f:
        region_line = f.readlines()[1]
    m = re.match(r'circle\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,', region_line)
    ra, dec = float(m.group(1)), float(m.group(2))

    file_list = red_ixpe.get_all_file_list(ixpe_dir)
    src_pcube, _bkg_pcube, _irf_date = red_ixpe.pcube(
        OBSID, file_list, ebinning=[2.0, 4.0, 6.0, 8.0], src_rad=SRC_RAD,
        bkg_inner_rad=BKG_INNER_RAD, bkg_outer_rad=BKG_OUTER_RAD,
        ra=ra, dec=dec,
    )
    bands = red_ixpe.pcube_bands(src_pcube)

    e_ctr = np.array([0.5 * (b["E_LO"] + b["E_HI"]) for b in bands])
    pd = np.array([b["PD"] for b in bands])
    pd_err = np.array([b["PD_err"] for b in bands])

    # Independent inverse-variance-weighted linear regression PD = a + b*E,
    # via the closed-form weighted least squares normal equations -- never
    # this helper's own pollin fit code.
    w = 1.0 / pd_err ** 2
    Sw, Swx, Swy = w.sum(), (w * e_ctr).sum(), (w * pd).sum()
    Swxx, Swxy = (w * e_ctr ** 2).sum(), (w * e_ctr * pd).sum()
    denom = Sw * Swxx - Swx ** 2
    pcube_slope = (Sw * Swxy - Swx * Swy) / denom

    assert np.sign(first_fit["PD_slope"]) == np.sign(pcube_slope), (
        f"pollin PD_slope={first_fit['PD_slope']:.5f} sign disagrees with "
        f"independent PCUBE weighted-linear trend slope={pcube_slope:.5f}")

    # Rough magnitude agreement: within a factor of 4 -- pollin ties the
    # continuum shape via a single powerlaw across the whole band while
    # PCUBE is a model-independent per-band measurement, so exact agreement
    # is not expected, only the same order of magnitude.
    ratio = first_fit["PD_slope"] / pcube_slope
    assert 0.25 <= ratio <= 4.0, (
        f"pollin PD_slope={first_fit['PD_slope']:.5f} vs PCUBE trend "
        f"slope={pcube_slope:.5f} (ratio={ratio:.2f}) not within rough "
        f"factor-of-4 agreement")


def test_repro_stable_across_calls(fresh_ixpe_load, first_fit):
    """Reproducibility: a second, independent fit_tbabs_pollin_powerlaw call
    on the identical already-loaded fresh-reduced spectra gives the same
    PD_A, PD_slope, PA_A, PA_slope, chi2 within tight tolerance -- the
    load_data->fit sequence is deterministic on unchanged data."""
    indices, _ = fresh_ixpe_load
    second_fit = xspec_help.fit_tbabs_pollin_powerlaw(indices)

    assert second_fit["PD_A"] == pytest.approx(first_fit["PD_A"], rel=1e-3)
    assert second_fit["PD_slope"] == pytest.approx(first_fit["PD_slope"], rel=1e-3, abs=1e-6)
    assert second_fit["PA_A"] == pytest.approx(first_fit["PA_A"], rel=1e-3)
    assert second_fit["PA_slope"] == pytest.approx(first_fit["PA_slope"], rel=1e-3, abs=1e-6)
    assert second_fit["chi2"] == pytest.approx(first_fit["chi2"], rel=1e-3)


# --------------------------------------------------------------------------- #
# Manifest wiring (AI-facing half)
# --------------------------------------------------------------------------- #

def test_tool_manifest_exposes_pollin_parameters():
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.fit_tbabs_pollin_powerlaw"]
    assert entry["module"] == "xspec_utils.lib.xspec_help"
    assert entry["fn"] == "fit_tbabs_pollin_powerlaw"
    props = entry["args_schema"]["properties"]
    for name in ("indices", "gamma_init", "pd_a_init", "pd_slope_init",
                 "pa_a_init", "pa_slope_init", "nh", "freeze_nh", "erange",
                 "conf_level", "n_iterations", "log_path"):
        assert name in props, f"{name} missing from args_schema.properties"
    assert entry["args_schema"]["required"] == ["indices"]
    # Brief/args-schema split: assert over all manifest text the model sees.
    manifest_text = " ".join([
        entry["description"],
        entry["args_schema"].get("description", ""),
        *(p.get("description", "") for p in props.values()),
    ])
    assert "0.683" in manifest_text
    assert "one free" in manifest_text.lower()
    assert "ixpe_changes" in manifest_text
    assert "pollin" in manifest_text
