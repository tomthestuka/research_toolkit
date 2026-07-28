"""Tests for xspec_utils' fit_tbabs_polconst_powerlaw (lib/xspec_help.py) and
its load_data companion fix (lib/xspec_help.py's IXPE ixpe_changes=True branch).

Runs a real PyXSPEC polconst spectropolarimetric cross-check fit against the
already-local IXPE Level-2 event files for Cyg X-1 ObsID 01002901, via a
fresh ixpeobssim reduction (ixpe_changes=True) at the paper's region radii --
copied into an isolated tmp directory per test-module run so this never
mutates the shared analyst-env data directory. Before this fix, load_data's
IXPE ixpe_changes=True branch never centroided a source position and always
raised ValueError (it called red_ixpe.pha with ra=None/dec=None, which
red_ixpe.pha requires whenever changes=True); the fix routes it through the
same find_src.run_extraction centroiding red_ixpe.reduce_pha already uses.
"""

from __future__ import annotations

import os
import re
import shutil

import pytest
from xspec import AllData, AllModels, Fit

from xspec_utils.lib import find_src, red_ixpe, xspec_help
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
    ixpe_dir = str(tmp_path_factory.mktemp("polconst_ixpe"))
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
    return xspec_help.fit_tbabs_polconst_powerlaw(indices)


# --------------------------------------------------------------------------- #
# load_data companion fix -- fresh IXPE reduction now actually centroids
# --------------------------------------------------------------------------- #

def test_differential_load_data_centroid_matches_independent_find_src_call(fresh_ixpe_load):
    """The RA/Dec load_data's ixpe_changes=True branch centroids at (used to
    resolve the response/BACKSCAL for this region) must match an independent
    direct find_src.run_extraction call on the same DU1 Level-2 event file --
    not load_data's own code path -- confirming the companion fix genuinely
    centroids rather than producing well-shaped-but-arbitrary metadata."""
    _indices, ixpe_dir = fresh_ixpe_load
    du1_evt = os.path.join(ixpe_dir, OBSID, "event_l2",
                           f"ixpe{OBSID}_det1_evt2_{VERSION}_rej.fits")
    assert os.path.exists(du1_evt)

    ref_dir = os.path.join(ixpe_dir, "find_src_ref")
    os.makedirs(ref_dir, exist_ok=True)
    ref_src_reg = os.path.join(ref_dir, "src.reg")
    ref_bkg_reg = os.path.join(ref_dir, "bkg.reg")
    find_src.run_extraction(du1_evt, ref_src_reg, ref_bkg_reg, None,
                             SRC_RAD, BKG_INNER_RAD, BKG_OUTER_RAD)
    with open(ref_src_reg) as f:
        region_line = f.readlines()[1]
    m = re.match(r'circle\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,', region_line)
    ref_ra, ref_dec = float(m.group(1)), float(m.group(2))

    # load_data's own centroid is what its region file (written next to the
    # event files, one level above event_l2/) actually recorded.
    src_reg = os.path.join(ixpe_dir, OBSID, "src.reg")
    assert os.path.exists(src_reg)
    with open(src_reg) as f:
        region_line = f.readlines()[1]
    m = re.match(r'circle\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,', region_line)
    ra, dec = float(m.group(1)), float(m.group(2))

    assert ra == pytest.approx(ref_ra, abs=1e-6)
    assert dec == pytest.approx(ref_dec, abs=1e-6)


def test_guardrail_multi_obsid_no_longer_raises_not_implemented(tmp_path):
    """Regression guard for the removed single-IXPE-ObsID guard-rail: a
    2-ObsID IXPE list must no longer short-circuit on NotImplementedError.
    ``ixpe_path`` deliberately resolves to zero files, so this reaches (and
    stops at) ``red_ixpe.get_all_file_list``'s own "no files found" check --
    a ValueError -- which only happens once load_data's IXPE branch has
    actually started processing a multi-ObsID list, rather than bailing out
    at the top on ``len(epoch['IXPE']) > 1``. This never calls AllData (the
    ValueError fires before any spectrum is loaded), so it cannot disturb
    the shared XSPEC session other tests in this module rely on -- the full
    multi-ObsID load/capability/differential/gain-linking coverage on real
    data lives in test_load_data_multi_ixpe.py."""
    with pytest.raises(ValueError, match="No files found"):
        xspec_help.load_data(
            epoch={"IXPE": [OBSID, "01250101"]},
            ixpe_path=str(tmp_path),
            clear_existing=False,
        )


# --------------------------------------------------------------------------- #
# fit_tbabs_polconst_powerlaw -- capability, output contract, differential
# --------------------------------------------------------------------------- #

def test_capability_background_applied_and_pd_pa_tied(fresh_ixpe_load, first_fit):
    """(a) Background was genuinely applied: net Stokes-I counts < gross I
    counts for every DU (never asserted for Q/U -- those can legitimately
    exceed gross). (b) PD/PA are genuinely tied: exactly one free PD and one
    free PA parameter across the whole session (inspect AllModels parameter
    links), never one per DU."""
    indices, _ = fresh_ixpe_load

    for du, triplet in indices["IXPE"][OBSID].items():
        i_idx = triplet[0]
        # Spectrum.rate = (net rate w/ background subtracted, net rate variance,
        # total rate without background, predicted model rate) -- PyXspec's own
        # order (xspec_help.audit_paths' inline comment mislabels this tuple).
        net, _var, total, _model = AllData(i_idx).rate
        assert net < total, f"{du}: net Stokes-I rate not below gross ({net} !< {total})"

    n_free_pd = n_free_pa = 0
    for dg in range(1, AllData.nGroups + 1):
        m = AllModels(dg)
        if not m.polconst.A.link and not m.polconst.A.frozen:
            n_free_pd += 1
        if not m.polconst.psi.link and not m.polconst.psi.frozen:
            n_free_pa += 1
    assert n_free_pd == 1, f"expected exactly one free PD parameter, got {n_free_pd}"
    assert n_free_pa == 1, f"expected exactly one free PA parameter, got {n_free_pa}"


def test_capability_output_contract(first_fit):
    """Every fitted number the brief requires must be in the returned dict --
    analysts never see XSPEC's own (muted) terminal output. Pins the trimmed
    contract exactly: the removed error-command bookkeeping (error_method)
    must not reappear."""
    for key in ("PD", "PA", "PD_err", "PA_err", "Gamma", "Gamma_err", "nH",
                "chi2", "dof", "stat_per_dof", "conf_level", "indices", "stdout", "exit_code"):
        assert key in first_fit, f"{key} missing from result"
    assert "error_method" not in first_fit

    assert first_fit["exit_code"] == 0
    assert 0.0 <= first_fit["PD"] <= 1.0
    assert first_fit["PD_err"] > 0
    assert first_fit["PA_err"] > 0
    assert first_fit["Gamma_err"] > 0
    assert first_fit["dof"] > 0
    assert first_fit["stat_per_dof"] == pytest.approx(first_fit["chi2"] / first_fit["dof"], rel=1e-9)
    assert first_fit["indices"] is not None and "IXPE" in first_fit["indices"]


def test_differential_polconst_matches_pcube_broadband(fresh_ixpe_load, first_fit):
    """Compare the polconst broadband PD/PA against the independently-trusted
    ixpe_pcube broadband PD/PA for the same ObsID/region over 2-8 keV --
    agreement within the combined error. ixpe_pcube is a route this toolkit
    already validated at <=2.3 sigma vs the published PCUBE result, so it is
    an independent reference, not this helper's own newly-added code."""
    _indices, ixpe_dir = fresh_ixpe_load

    src_reg = os.path.join(ixpe_dir, OBSID, "src.reg")
    with open(src_reg) as f:
        region_line = f.readlines()[1]
    m = re.match(r'circle\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,', region_line)
    ra, dec = float(m.group(1)), float(m.group(2))

    file_list = red_ixpe.get_all_file_list(ixpe_dir)
    src_pcube, _bkg_pcube, _irf_date = red_ixpe.pcube(
        OBSID, file_list, ebinning=[2.0, 8.0], src_rad=SRC_RAD,
        bkg_inner_rad=BKG_INNER_RAD, bkg_outer_rad=BKG_OUTER_RAD,
        ra=ra, dec=dec,
    )
    band = red_ixpe.pcube_bands(src_pcube)[0]

    combined_pd_err = (first_fit["PD_err"] ** 2 + band["PD_err"] ** 2) ** 0.5
    assert abs(first_fit["PD"] - band["PD"]) <= 2 * combined_pd_err, (
        f"polconst PD={first_fit['PD']:.4f}+/-{first_fit['PD_err']:.4f} vs "
        f"pcube PD={band['PD']:.4f}+/-{band['PD_err']:.4f}")

    combined_pa_err = (first_fit["PA_err"] ** 2 + band["PA_err"] ** 2) ** 0.5
    assert abs(first_fit["PA"] - band["PA"]) <= 2 * combined_pa_err, (
        f"polconst PA={first_fit['PA']:.3f}+/-{first_fit['PA_err']:.3f} vs "
        f"pcube PA={band['PA']:.3f}+/-{band['PA_err']:.3f}")


def test_repro_stable_across_calls(fresh_ixpe_load, first_fit):
    """Reproducibility: a second, independent fit_tbabs_polconst_powerlaw
    call on the identical already-loaded fresh-reduced spectra gives the
    same PD, PA, Gamma, chi2 within tight tolerance -- the load_data->fit
    sequence is deterministic on unchanged data."""
    indices, _ = fresh_ixpe_load
    second_fit = xspec_help.fit_tbabs_polconst_powerlaw(indices)

    assert second_fit["PD"] == pytest.approx(first_fit["PD"], rel=1e-3)
    assert second_fit["PA"] == pytest.approx(first_fit["PA"], rel=1e-3)
    assert second_fit["Gamma"] == pytest.approx(first_fit["Gamma"], rel=1e-3)
    assert second_fit["chi2"] == pytest.approx(first_fit["chi2"], rel=1e-3)


# --------------------------------------------------------------------------- #
# Manifest wiring (AI-facing half)
# --------------------------------------------------------------------------- #

def test_tool_manifest_exposes_class_parameters():
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.fit_tbabs_polconst_powerlaw"]
    assert entry["module"] == "xspec_utils.lib.xspec_help"
    assert entry["fn"] == "fit_tbabs_polconst_powerlaw"
    props = entry["args_schema"]["properties"]
    for name in ("indices", "gamma_init", "pd_init", "pa_init", "nh", "freeze_nh",
                 "erange", "conf_level", "n_iterations", "log_path"):
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


def test_load_data_manifest_exposes_ixpe_fresh_reduction_knobs():
    """Regression guard for the load_data companion fix's AI-facing half:
    the manifest args_schema must expose ixpe_changes and the three region
    radii, or a tool-calling analyst can never reach the fixed fresh-
    reduction path at all."""
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.load_data"]
    props = entry["args_schema"]["properties"]
    for name in ("ixpe_changes", "ixpe_src_rad", "ixpe_bkg_inner_rad", "ixpe_bkg_outer_rad"):
        assert name in props, f"{name} missing from load_data's args_schema.properties"
    assert "ixpe_changes" not in entry["args_schema"]["required"]
