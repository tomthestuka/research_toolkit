"""Tests for xspec_utils' bhxrb_pd_hr_survey_plot (`bhxrb_pd_hr_survey.pd_hr_survey_plot`).

Reuses test_bhxrb_pd_hr_survey.py's fixture pattern: a real bhxrb_hardness_ratios.csv
is computed once per session by running the real `ixpe_HR` skill against real,
already-local per-DU Stokes-I PHA1 products for four Cyg X-1 IXPE ObsIDs
(01002901, 01250401, 02250901, 03005101), copied into an isolated scratch mission
directory so the real /data/fzhou/XSPEC/Cyg_X-1/IXPE catalog.json is never touched.
That HR CSV plus the real, already-local pcube_timeseg_<obsid>.csv integrated rows
feed `pd_hr_survey` to build a genuine `bhxrb_pd_hardness_survey.csv`, which is the
plotter's actual input -- exactly mirroring the real chain
(ixpe_pcube_timeseg -> ixpe_HR -> bhxrb_pd_hr_survey -> bhxrb_pd_hr_survey_plot).
"""

from __future__ import annotations

import csv
import hashlib
import math
import os
import shutil

import pytest
from scipy import stats as scipy_stats

from xspec_utils.lib import bhxrb_pd_hr_survey as survey
from xspec_utils.skills import TOOL_MANIFEST
from xspec_utils.skills import session as ixpe_session

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
OBSID_INFO = {  # obsid -> (version, rej)
    "01002901": ("v01", True),
    "01250401": ("v03", False),
    "02250901": ("v01", False),
    "03005101": ("v01", False),
}
ALL_OBSIDS = list(OBSID_INFO)
DETECTED_OBSIDS = ["01002901", "01250401", "02250901"]
NON_DETECTED_OBSIDS = ["03005101"]

SUFFIXES = ("pd_vs_hr_combined", "pd_vs_hr_per_source", "pd_distribution_comparison",
            "hr_detection_bias", "feasibility_mdp99_exposure")


def _copy_event_and_src_pha(dest_dir: str, obsid: str, version: str, rej: bool) -> None:
    tag = f"{version}_rej" if rej else version
    src_evt_dir = os.path.join(REAL_IXPE_DIR, obsid, "event_l2")
    dst_evt_dir = os.path.join(dest_dir, obsid, "event_l2")
    os.makedirs(dst_evt_dir, exist_ok=True)
    for det in (1, 2, 3):
        base = f"ixpe{obsid}_det{det}_evt2_{tag}"
        evt_name = f"{base}.fits"
        src_evt = os.path.join(src_evt_dir, evt_name)
        if not os.path.exists(src_evt):
            pytest.skip(f"real fixture ObsID {obsid} not local at {src_evt}")
        shutil.copyfile(src_evt, os.path.join(dst_evt_dir, evt_name))
        pha_name = f"{base}_src_pha1.fits"
        src_pha = os.path.join(src_evt_dir, pha_name)
        if not os.path.exists(src_pha):
            pytest.skip(f"real fixture PHA1 product not local at {src_pha}")
        shutil.copyfile(src_pha, os.path.join(dst_evt_dir, pha_name))


def _require_pcube_fixtures():
    for obsid in ALL_OBSIDS:
        path = os.path.join(REAL_IXPE_DIR, obsid, f"pcube_timeseg_{obsid}.csv")
        if not os.path.exists(path):
            pytest.skip(f"real fixture pcube_timeseg CSV not local at {path}")


@pytest.fixture(scope="module")
def real_hr_csv(tmp_path_factory):
    scratch = str(tmp_path_factory.mktemp("pd_hr_plot_hr_scratch"))
    for obsid, (version, rej) in OBSID_INFO.items():
        _copy_event_and_src_pha(scratch, obsid, version, rej)
    result = ixpe_session.ixpe_HR(scratch, ALL_OBSIDS)
    assert result["exit_code"] == 0, result.get("stderr")
    assert result["errors"] == {}
    return result["csv_path"]


@pytest.fixture(scope="module")
def real_survey_csv(tmp_path_factory, real_hr_csv):
    _require_pcube_fixtures()
    out_dir = str(tmp_path_factory.mktemp("pd_hr_plot_survey"))
    obsids_with_sources = {obsid: "Cyg X-1" for obsid in ALL_OBSIDS}
    result = survey.pd_hr_survey(out_dir, obsids_with_sources, real_hr_csv, REAL_IXPE_DIR)
    assert result["exit_code"] == 0
    return result["survey_csv"]


# --------------------------------------------------------------------------- #
# Capability: exactly five named, non-empty files, sha256 matches file bytes
# --------------------------------------------------------------------------- #

def test_produces_five_named_files(tmp_path, real_survey_csv):
    prefix = str(tmp_path / "pdhr")
    result = survey.pd_hr_survey_plot(REAL_IXPE_DIR, real_survey_csv, save_path_prefix=prefix)
    assert result["exit_code"] == 0
    assert len(result["figures"]) == 5

    seen_suffixes = set()
    for fig in result["figures"]:
        path = fig["path"]
        assert os.path.exists(path)
        data = open(path, "rb").read()
        assert len(data) > 0
        assert fig["sha256"] == hashlib.sha256(data).hexdigest()
        assert fig["mime_type"] == "image/png"
        assert fig["produced_by"] == "xspec_utils.lib.bhxrb_pd_hr_survey.pd_hr_survey_plot"
        assert fig["caption"]
        matched = [s for s in SUFFIXES if path == f"{prefix}_{s}.png"]
        assert len(matched) == 1, f"unexpected figure path {path}"
        seen_suffixes.add(matched[0])
    assert seen_suffixes == set(SUFFIXES)

    written = sorted(p for p in os.listdir(tmp_path) if p.startswith("pdhr_"))
    assert len(written) == 5


def test_pdf_format_written_with_correct_mime(tmp_path, real_survey_csv):
    prefix = str(tmp_path / "pdhr_pdf")
    result = survey.pd_hr_survey_plot(REAL_IXPE_DIR, real_survey_csv, save_path_prefix=prefix, fmt="pdf")
    assert result["exit_code"] == 0
    assert len(result["figures"]) == 5
    for fig in result["figures"]:
        assert fig["path"].endswith(".pdf")
        assert fig["mime_type"] == "application/pdf"
        assert os.path.getsize(fig["path"]) > 0


# --------------------------------------------------------------------------- #
# Reproducibility: identical calls -> byte-identical figures
# --------------------------------------------------------------------------- #

def test_repro_pd_hr_survey_plot(tmp_path_factory, real_survey_csv):
    out_a = str(tmp_path_factory.mktemp("pdhr_plot_a") / "run")
    out_b = str(tmp_path_factory.mktemp("pdhr_plot_b") / "run")
    a = survey.pd_hr_survey_plot(REAL_IXPE_DIR, real_survey_csv, save_path_prefix=out_a)
    b = survey.pd_hr_survey_plot(REAL_IXPE_DIR, real_survey_csv, save_path_prefix=out_b)

    def _by_suffix(figs):
        return {s: f["sha256"] for f in figs for s in SUFFIXES if f["path"].endswith(f"_{s}.png")}

    sha_a, sha_b = _by_suffix(a["figures"]), _by_suffix(b["figures"])
    assert set(sha_a) == set(SUFFIXES) == set(sha_b)
    for s in SUFFIXES:
        assert sha_a[s] == sha_b[s], f"figure {s} not byte-identical across repeated calls"


# --------------------------------------------------------------------------- #
# Differential: plotted quantities vs an independent recompute from the raw
# survey CSV values (never against the plotter's own helpers).
# --------------------------------------------------------------------------- #

def test_combined_points_match_independent_csv_read(tmp_path, real_survey_csv):
    with open(real_survey_csv, newline="") as f:
        independent_detected = {
            r["obsid"]: (float(r["PD_pct"]), float(r["HR"]))
            for r in csv.DictReader(f) if not r["note"] and r["detected"] == "True"
        }
    assert sorted(independent_detected) == sorted(DETECTED_OBSIDS)

    result = survey.pd_hr_survey_plot(REAL_IXPE_DIR, real_survey_csv, save_path_prefix=str(tmp_path / "diff"))
    plotted = {p["obsid"]: (p["PD_pct"], p["HR"]) for p in result["plotted"]["combined_points"]}
    assert set(plotted) == set(independent_detected)
    for obsid, (pd_pct, hr) in independent_detected.items():
        assert plotted[obsid][0] == pytest.approx(pd_pct, rel=1e-9)
        assert plotted[obsid][1] == pytest.approx(hr, rel=1e-9)


def test_detection_bias_hr_matches_independent_csv_read_and_scipy(tmp_path, real_survey_csv):
    with open(real_survey_csv, newline="") as f:
        rows = [r for r in csv.DictReader(f) if not r["note"]]
    expected_det = sorted(float(r["HR"]) for r in rows if r["detected"] == "True")
    expected_nondet = sorted(float(r["HR"]) for r in rows if r["detected"] == "False")

    result = survey.pd_hr_survey_plot(REAL_IXPE_DIR, real_survey_csv, save_path_prefix=str(tmp_path / "bias"))
    bias = result["plotted"]["detection_bias_hr"]
    assert sorted(bias["detected"]) == pytest.approx(expected_det, rel=1e-9)
    assert sorted(bias["non_detected"]) == pytest.approx(expected_nondet, rel=1e-9)

    ref = scipy_stats.mannwhitneyu(expected_det, expected_nondet, alternative="two-sided")
    mine = survey._mann_whitney(expected_det, expected_nondet)
    assert mine["statistic"] == pytest.approx(float(ref.statistic), rel=1e-9)
    assert mine["p"] == pytest.approx(float(ref.pvalue), rel=1e-9)


def test_mdp99_exposure_matches_independent_pcube_read(tmp_path, real_survey_csv):
    result = survey.pd_hr_survey_plot(REAL_IXPE_DIR, real_survey_csv, save_path_prefix=str(tmp_path / "feas"))
    plotted = {p["obsid"]: p for p in result["plotted"]["mdp99_exposure"]}
    assert set(plotted) == set(ALL_OBSIDS)
    for obsid, p in plotted.items():
        csv_path = os.path.join(REAL_IXPE_DIR, obsid, f"pcube_timeseg_{obsid}.csv")
        with open(csv_path, newline="") as f:
            f.readline()  # provenance header
            for r in csv.DictReader(f):
                if int(r["segment_id"]) == -1:
                    assert p["exposure_s"] == pytest.approx(float(r["exposure_s"]), rel=1e-9)
                    assert p["MDP99"] == pytest.approx(float(r["MDP99"]), rel=1e-9)
                    break
            else:
                pytest.fail(f"no integrated row found for {obsid}")


# --------------------------------------------------------------------------- #
# Edge cases
# --------------------------------------------------------------------------- #

def test_missing_survey_csv_fails_loudly(tmp_path):
    missing = str(tmp_path / "does_not_exist_survey.csv")
    result = survey.pd_hr_survey_plot(REAL_IXPE_DIR, missing, save_path_prefix=str(tmp_path / "nosurvey"))
    assert result["exit_code"] == 1
    assert "stderr" in result and result["stderr"]
    assert "figures" not in result
    assert not any(p.startswith("nosurvey") for p in os.listdir(tmp_path))


def test_wrong_schema_survey_csv_fails_loudly(tmp_path):
    bad = tmp_path / "bad_schema_survey.csv"
    with open(bad, "w", newline="") as f:
        f.write("not_obsid,foo\n1,2\n")
    prefix = str(tmp_path / "badsurvey")
    result = survey.pd_hr_survey_plot(REAL_IXPE_DIR, str(bad), save_path_prefix=prefix)
    assert result["exit_code"] == 1
    assert "stderr" in result and result["stderr"]
    assert "figures" not in result
    assert not any(p.startswith("badsurvey") for p in os.listdir(tmp_path))


def test_all_obsids_skipped_upstream_exit_code_1_no_figures(tmp_path):
    """Every row carries a non-empty note (e.g. every obsid failed upstream) ->
    nothing usable to plot -> exit_code=1, no files written."""
    all_missing = tmp_path / "all_missing_survey.csv"
    with open(all_missing, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=survey.SURVEY_FIELDS)
        writer.writeheader()
        writer.writerow({"obsid": "99999999", "source": "Cyg X-1", "PD_pct": "", "PD_err": "",
                          "MDP99": "", "detected": "", "mdp_flagged": "", "HR": "", "HR_err": "",
                          "note": "no readable pcube_timeseg CSV"})
    prefix = str(tmp_path / "allmissing")
    result = survey.pd_hr_survey_plot(REAL_IXPE_DIR, str(all_missing), save_path_prefix=prefix)
    assert result["exit_code"] == 1
    assert "stderr" in result and result["stderr"]
    assert "figures" not in result
    assert not any(p.startswith("allmissing") for p in os.listdir(tmp_path))


def test_empty_survey_csv_zero_rows_fails_loudly(tmp_path):
    empty = tmp_path / "empty_survey.csv"
    with open(empty, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=survey.SURVEY_FIELDS)
        writer.writeheader()
    prefix = str(tmp_path / "emptysurvey")
    result = survey.pd_hr_survey_plot(REAL_IXPE_DIR, str(empty), save_path_prefix=prefix)
    assert result["exit_code"] == 1
    assert "stderr" in result and result["stderr"]
    assert "figures" not in result
    assert not any(p.startswith("emptysurvey") for p in os.listdir(tmp_path))


def test_no_soft_hard_groups_annotates_not_crashes(tmp_path, real_survey_csv):
    """Default soft_state_sources/hard_state_sources are empty -- the
    distribution-comparison panel must still render (annotated), not crash,
    and the returned group data must both be empty lists."""
    result = survey.pd_hr_survey_plot(REAL_IXPE_DIR, real_survey_csv, save_path_prefix=str(tmp_path / "nogroups"))
    assert result["exit_code"] == 0
    assert result["plotted"]["distribution_groups"] == {"soft": [], "hard": []}
    fig = [f for f in result["figures"] if f["path"].endswith("_pd_distribution_comparison.png")][0]
    assert os.path.getsize(fig["path"]) > 0


def test_soft_hard_groups_populate_distribution_panel(tmp_path, real_survey_csv):
    result = survey.pd_hr_survey_plot(
        REAL_IXPE_DIR, real_survey_csv, save_path_prefix=str(tmp_path / "groups"),
        soft_state_sources=["Cyg X-1"], hard_state_sources=[])
    groups = result["plotted"]["distribution_groups"]
    assert len(groups["soft"]) == len(DETECTED_OBSIDS)
    assert groups["hard"] == []


def test_invalid_fmt_rejected(tmp_path, real_survey_csv):
    with pytest.raises(ValueError):
        survey.pd_hr_survey_plot(REAL_IXPE_DIR, real_survey_csv,
                                  save_path_prefix=str(tmp_path / "bad"), fmt="svg")


# --------------------------------------------------------------------------- #
# Manifest wiring (AI-facing half) + no science-stack import
# --------------------------------------------------------------------------- #

def test_tool_manifest_wired_to_lib_directly_no_science_stack():
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.bhxrb_pd_hr_survey_plot"]
    assert entry["module"] == "xspec_utils.lib.bhxrb_pd_hr_survey"
    assert entry["fn"] == "pd_hr_survey_plot"
    props = entry["args_schema"]["properties"]
    for name in ("obs_db_path", "survey_csv", "save_path_prefix", "fmt",
                 "soft_state_sources", "hard_state_sources", "mdp_flag_pct"):
        assert name in props, f"{name} missing from args_schema.properties"
    assert entry["args_schema"]["required"] == ["obs_db_path", "survey_csv"]
    assert props["fmt"]["enum"] == ["png", "pdf"]


def test_lib_module_imports_no_xspec_or_ixpeobssim():
    import ast
    import inspect
    tree = ast.parse(inspect.getsource(survey))
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    assert "xspec" not in imported
    assert "ixpeobssim" not in imported
    assert "red_ixpe" not in imported
