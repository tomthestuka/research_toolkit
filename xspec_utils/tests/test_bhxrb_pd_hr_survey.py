"""Tests for xspec_utils' bhxrb_pd_hr_survey (lib/bhxrb_pd_hr_survey.py).

PD comes from the real, already-local `pcube_timeseg_<obsid>.csv` integrated
(`segment_id=-1`) rows for four Cyg X-1 IXPE ObsIDs (01002901, 01250401,
02250901, 03005101) under `/data/fzhou/XSPEC/Cyg_X-1/IXPE` -- read-only; this
tool never writes there (`pcube_csv_dir` is a read-only input). Independently
checked before writing this test: the real integrated rows give
`PD_pct > MDP99` (detected) for 01002901/01250401/02250901 and `PD_pct < MDP99`
(non-detected) for 03005101 -- a genuine 3-detected/1-non-detected split that
also naturally exercises the `min_n_intra=4` "omit the source" edge case
without any hand-built fixture.

HR comes from running the real `ixpe_HR` skill against each ObsID's real,
already-local per-DU Stokes-I PHA1 products, copied into an isolated scratch
mission directory (mirroring `test_ixpe_HR.py`'s pattern) so the HR CSV this
tool consumes is genuine, computed once per test session, while the real
mission's `bhxrb_hardness_ratios.csv`/`catalog.json` are never touched.

`obsids_with_sources` below groups the four ObsIDs under test-chosen labels
purely to exercise the per-source / soft-vs-hard grouping mechanism the tool
is built around (the brief requires source-group membership to be a caller
argument, never hardcoded) -- all four are the same real Cyg X-1 source, so
this validates the statistics on genuine real PD/HR values, not a real
astrophysical state split.
"""

from __future__ import annotations

import csv
import math
import os
import shutil

import numpy as np
import pytest
from scipy import stats as scipy_stats

from obsdb_catalog import neutral
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
DETECTED_OBSIDS = ["01002901", "01250401", "02250901"]  # PD_pct > MDP99 in the real CSVs
NON_DETECTED_OBSIDS = ["03005101"]                       # PD_pct < MDP99 in the real CSV
NO_HR_OBSID = "02001901"  # real pcube_timeseg CSV exists, but never in real_hr_csv


def _copy_event_and_src_pha(dest_dir: str, obsid: str, version: str, rej: bool) -> list[str]:
    """Copy real Level-2 event files + real per-DU src Stokes-I PHA1 products for
    `obsid` into `<dest_dir>/<obsid>/event_l2/` (bkg PHA1 not needed: ixpe_HR's
    default bkg=False never reads it). Returns the copied event file paths."""
    tag = f"{version}_rej" if rej else version
    src_evt_dir = os.path.join(REAL_IXPE_DIR, obsid, "event_l2")
    dst_evt_dir = os.path.join(dest_dir, obsid, "event_l2")
    os.makedirs(dst_evt_dir, exist_ok=True)
    file_list = []
    for det in (1, 2, 3):
        base = f"ixpe{obsid}_det{det}_evt2_{tag}"
        evt_name = f"{base}.fits"
        src_evt = os.path.join(src_evt_dir, evt_name)
        if not os.path.exists(src_evt):
            pytest.skip(f"real fixture ObsID {obsid} not local at {src_evt}")
        dst_evt = os.path.join(dst_evt_dir, evt_name)
        shutil.copyfile(src_evt, dst_evt)
        file_list.append(dst_evt)
        pha_name = f"{base}_src_pha1.fits"
        src_pha = os.path.join(src_evt_dir, pha_name)
        if not os.path.exists(src_pha):
            pytest.skip(f"real fixture PHA1 product not local at {src_pha}")
        shutil.copyfile(src_pha, os.path.join(dst_evt_dir, pha_name))
    return sorted(file_list)


@pytest.fixture(scope="module")
def real_hr_csv(tmp_path_factory):
    """A real bhxrb_hardness_ratios.csv for the four real ObsIDs, computed by
    ixpe_HR from genuine PHA1 products copied into an isolated scratch mission
    directory -- never the real /data/fzhou/XSPEC/Cyg_X-1/IXPE catalog.json."""
    scratch = str(tmp_path_factory.mktemp("bhxrb_hr_scratch"))
    for obsid, (version, rej) in OBSID_INFO.items():
        _copy_event_and_src_pha(scratch, obsid, version, rej)
    result = ixpe_session.ixpe_HR(scratch, ALL_OBSIDS)
    assert result["exit_code"] == 0, result.get("stderr")
    assert result["errors"] == {}
    return result["csv_path"]


def _same_source(label: str = "Cyg X-1") -> dict:
    return {obsid: label for obsid in ALL_OBSIDS}


def _two_groups() -> dict:
    return {"01002901": "GroupA", "01250401": "GroupA",
            "02250901": "GroupB", "03005101": "GroupB"}


# --------------------------------------------------------------------------- #
# Pure-logic unit tests: readers + statistics helpers
# --------------------------------------------------------------------------- #

def test_read_pcube_integrated_round_trip(tmp_path):
    csv_path = str(tmp_path / "pcube_timeseg_roundtrip.csv")
    with open(csv_path, "w", newline="") as f:
        f.write(survey._PCUBE_PROVENANCE_PREFIX + '{"segment_duration_s": 5000.0}\n')
        fields = ["segment_id", "t_start_mjd", "t_stop_mjd", "exposure_s",
                  "PD_pct", "PD_err", "PA_deg", "PA_err", "QI", "QI_err",
                  "UI", "UI_err", "MDP99", "SIGNIF"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        row0 = {k: 0.0 for k in fields}
        row0.update(segment_id=0, PD_pct=8.0, MDP99=1.0)
        writer.writerow(row0)
        row_int = {k: 0.0 for k in fields}
        row_int.update(segment_id=-1, PD_pct=3.5, PD_err=0.2, MDP99=0.6)
        writer.writerow(row_int)

    result = survey._read_pcube_integrated(csv_path)
    assert result == {"PD_pct": 3.5, "PD_err": 0.2, "MDP99": 0.6}


def test_read_pcube_integrated_missing_header_treated_as_absent(tmp_path):
    csv_path = str(tmp_path / "pcube_timeseg_no_header.csv")
    with open(csv_path, "w", newline="") as f:
        f.write("segment_id,PD_pct,PD_err,MDP99\n-1,3.5,0.2,0.6\n")
    assert survey._read_pcube_integrated(csv_path) is None


def test_read_pcube_integrated_missing_file_returns_none(tmp_path):
    assert survey._read_pcube_integrated(str(tmp_path / "absent.csv")) is None


def test_read_pcube_integrated_no_integrated_row_returns_none(tmp_path):
    csv_path = str(tmp_path / "pcube_timeseg_no_int.csv")
    with open(csv_path, "w", newline="") as f:
        f.write(survey._PCUBE_PROVENANCE_PREFIX + '{}\n')
        fields = ["segment_id", "PD_pct", "PD_err", "MDP99"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerow({"segment_id": 0, "PD_pct": 1.0, "PD_err": 0.1, "MDP99": 0.5})
    assert survey._read_pcube_integrated(csv_path) is None


def test_read_hr_rows_missing_file_raises(tmp_path):
    with pytest.raises(ValueError, match="hr_csv not found"):
        survey._read_hr_rows(str(tmp_path / "absent.csv"))


def test_read_hr_rows_parses_columns(tmp_path):
    csv_path = str(tmp_path / "bhxrb_hardness_ratios.csv")
    with open(csv_path, "w", newline="") as f:
        f.write("obsid,soft_rate,soft_err,hard_rate,hard_err,HR,HR_err\n")
        f.write("01002901,1.0,0.1,2.0,0.1,0.3333333333333333,0.05\n")
    rows = survey._read_hr_rows(csv_path)
    assert rows == {"01002901": {"HR": pytest.approx(1.0 / 3.0), "HR_err": pytest.approx(0.05)}}


def test_spearman_matches_scipy_directly():
    x, y = [1.0, 2.0, 3.0, 4.0], [4.0, 1.0, 3.0, 2.0]
    ref = scipy_stats.spearmanr(x, y)
    r = survey._spearman(x, y)
    assert r["n"] == 4
    assert r["rho"] == pytest.approx(float(ref.statistic), rel=1e-9)
    assert r["p"] == pytest.approx(float(ref.pvalue), rel=1e-9)


def test_spearman_nan_below_three_pairs():
    r = survey._spearman([1.0, 2.0], [2.0, 1.0])
    assert r["n"] == 2
    assert math.isnan(r["rho"]) and math.isnan(r["p"])


def test_spearman_ignores_nonfinite_pairs():
    r = survey._spearman([1.0, 2.0, 3.0, float("nan")], [4.0, 1.0, 3.0, 2.0])
    assert r["n"] == 3


def test_ks_two_sample_matches_scipy_directly():
    a, b = [0.1, 0.2, 0.3, 0.9], [0.5, 0.6, 0.7, 0.8]
    ref = scipy_stats.ks_2samp(a, b)
    r = survey._ks_two_sample(a, b)
    assert r["n_soft"] == 4 and r["n_hard"] == 4
    assert r["statistic"] == pytest.approx(float(ref.statistic), rel=1e-9)
    assert r["p"] == pytest.approx(float(ref.pvalue), rel=1e-9)


def test_ks_two_sample_empty_group_is_nan_not_exception():
    r = survey._ks_two_sample([], [1.0, 2.0])
    assert r["n_soft"] == 0 and r["n_hard"] == 2
    assert math.isnan(r["statistic"]) and math.isnan(r["p"])
    r2 = survey._ks_two_sample([1.0], [])
    assert math.isnan(r2["statistic"]) and math.isnan(r2["p"])


def test_mann_whitney_matches_scipy_directly():
    a, b = [1.0, 2.0, 3.0], [4.0, 5.0, 6.0, 7.0]
    ref = scipy_stats.mannwhitneyu(a, b, alternative="two-sided")
    r = survey._mann_whitney(a, b)
    assert r["n_detected"] == 3 and r["n_non_detected"] == 4
    assert r["statistic"] == pytest.approx(float(ref.statistic), rel=1e-9)
    assert r["p"] == pytest.approx(float(ref.pvalue), rel=1e-9)


def test_mann_whitney_empty_group_is_nan_not_exception():
    r = survey._mann_whitney([], [1.0, 2.0])
    assert r["n_detected"] == 0 and r["n_non_detected"] == 2
    assert math.isnan(r["statistic"]) and math.isnan(r["p"])


# --------------------------------------------------------------------------- #
# Capability: real detection classification matches an independent recompute
# --------------------------------------------------------------------------- #

def test_capability_detection_matches_independent_recompute_from_real_csvs(tmp_path, real_hr_csv):
    """Proves the tool actually consumed the real PCUBE integrated rows (not
    well-shaped dummy metadata): recompute detected=PD_pct>MDP99 directly
    from each real pcube_timeseg CSV with an independent reader, and check it
    matches the tool's own classification for every obsid."""
    independent = {}
    for obsid in ALL_OBSIDS:
        csv_path = os.path.join(REAL_IXPE_DIR, obsid, f"pcube_timeseg_{obsid}.csv")
        with open(csv_path, newline="") as f:
            f.readline()  # provenance header
            for r in csv.DictReader(f):
                if int(r["segment_id"]) == -1:
                    independent[obsid] = float(r["PD_pct"]) > float(r["MDP99"])
                    break

    result = survey.pd_hr_survey(str(tmp_path), _same_source(), real_hr_csv, REAL_IXPE_DIR)
    assert result["exit_code"] == 0
    by_obsid = {r["obsid"]: r for r in result["rows"]}
    for obsid in ALL_OBSIDS:
        assert by_obsid[obsid]["note"] == ""
        assert by_obsid[obsid]["detected"] == independent[obsid]

    assert sorted(o for o, r in by_obsid.items() if r["detected"]) == sorted(DETECTED_OBSIDS)
    assert sorted(o for o, r in by_obsid.items() if not r["detected"]) == sorted(NON_DETECTED_OBSIDS)

    with open(result["survey_csv"], newline="") as f:
        header = f.readline().strip().split(",")
    assert header == survey.SURVEY_FIELDS


# --------------------------------------------------------------------------- #
# Differential: combined Spearman / KS / Mann-Whitney vs direct scipy recompute
# --------------------------------------------------------------------------- #

def test_differential_combined_spearman_ks_mwu_match_scipy_directly(tmp_path, real_hr_csv):
    result = survey.pd_hr_survey(str(tmp_path), _two_groups(), real_hr_csv, REAL_IXPE_DIR,
                                  soft_state_sources=["GroupA"], hard_state_sources=["GroupB"])
    assert result["exit_code"] == 0
    by_obsid = {r["obsid"]: r for r in result["rows"]}

    detected = [by_obsid[o] for o in DETECTED_OBSIDS]
    non_detected = [by_obsid[o] for o in NON_DETECTED_OBSIDS]

    pd_vals = np.array([r["PD_pct"] for r in detected])
    hr_vals = np.array([r["HR"] for r in detected])
    ref_spear = scipy_stats.spearmanr(pd_vals, hr_vals)
    assert result["spearman_combined"]["rho"] == pytest.approx(float(ref_spear.statistic), rel=1e-9)
    assert result["spearman_combined"]["p"] == pytest.approx(float(ref_spear.pvalue), rel=1e-9)

    soft_pd = np.array([r["PD_pct"] for r in detected if r["source"] == "GroupA"])
    hard_pd = np.array([r["PD_pct"] for r in detected if r["source"] == "GroupB"])
    ref_ks = scipy_stats.ks_2samp(soft_pd, hard_pd)
    assert result["ks_result"]["statistic"] == pytest.approx(float(ref_ks.statistic), rel=1e-9)
    assert result["ks_result"]["p"] == pytest.approx(float(ref_ks.pvalue), rel=1e-9)

    det_hr = np.array([r["HR"] for r in detected])
    nondet_hr = np.array([r["HR"] for r in non_detected])
    ref_mwu = scipy_stats.mannwhitneyu(det_hr, nondet_hr, alternative="two-sided")
    assert result["mwu_result"]["statistic"] == pytest.approx(float(ref_mwu.statistic), rel=1e-9)
    assert result["mwu_result"]["p"] == pytest.approx(float(ref_mwu.pvalue), rel=1e-9)


def test_pd_hr_survey_ks_grouping_uses_caller_supplied_groups(tmp_path, real_hr_csv):
    """The KS split is computed only when soft_state_sources/hard_state_sources
    are actually supplied by the caller -- proves the groups are consumed as an
    argument, not any hardcoded source list."""
    result = survey.pd_hr_survey(
        str(tmp_path), _two_groups(), real_hr_csv, REAL_IXPE_DIR,
        soft_state_sources=["GroupA"], hard_state_sources=["GroupB"])
    assert result["ks_result"]["n_soft"] == 2
    assert result["ks_result"]["n_hard"] == 1
    assert not math.isnan(result["ks_result"]["statistic"])


# --------------------------------------------------------------------------- #
# Edge cases: min_n_intra omission, missing CSV/HR row, empty-group NaN
# --------------------------------------------------------------------------- #

def test_source_below_min_n_intra_omitted_not_crashed(tmp_path, real_hr_csv):
    """Only 3 detected ObsIDs share the 'Cyg X-1' label -- below the default
    min_n_intra=4 -- so the source must be omitted from
    intra_source_correlations.csv (and per_source_correlations), never crash."""
    result = survey.pd_hr_survey(str(tmp_path), _same_source(), real_hr_csv, REAL_IXPE_DIR)
    assert result["exit_code"] == 0
    assert "Cyg X-1" not in result["per_source_correlations"]
    with open(result["intra_csv"], newline="") as f:
        rows = list(csv.DictReader(f))
    assert rows == []
    # But the combined correlation (n=3, no per-source floor) is still computed.
    assert result["spearman_combined"]["n"] == 3
    assert not math.isnan(result["spearman_combined"]["rho"])


def test_source_at_or_above_min_n_intra_computed_and_matches_scipy(tmp_path, real_hr_csv):
    result = survey.pd_hr_survey(
        str(tmp_path), _same_source(), real_hr_csv, REAL_IXPE_DIR, min_n_intra=3)
    assert "Cyg X-1" in result["per_source_correlations"]
    res = result["per_source_correlations"]["Cyg X-1"]
    by_obsid = {r["obsid"]: r for r in result["rows"]}
    detected = [by_obsid[o] for o in DETECTED_OBSIDS]
    ref = scipy_stats.spearmanr([r["PD_pct"] for r in detected], [r["HR"] for r in detected])
    assert res["rho"] == pytest.approx(float(ref.statistic), rel=1e-9)
    with open(result["intra_csv"], newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1 and rows[0]["source"] == "Cyg X-1"


def test_missing_pcube_csv_recorded_as_note_and_skipped_not_raised(tmp_path, real_hr_csv):
    obsids = {**_same_source(), "99999999": "Cyg X-1"}
    result = survey.pd_hr_survey(str(tmp_path), obsids, real_hr_csv, REAL_IXPE_DIR)
    assert result["exit_code"] == 0
    assert "99999999" in result["skipped"]
    by_obsid = {r["obsid"]: r for r in result["rows"]}
    assert "no readable pcube_timeseg CSV" in by_obsid["99999999"]["note"]
    assert by_obsid["99999999"]["detected"] == ""
    # And it must never reach the catalog.
    assert neutral.read(str(tmp_path))["entries"].get("99999999") is None


def test_missing_hr_row_recorded_as_note_and_skipped_not_raised(tmp_path, real_hr_csv):
    """02001901 carries a real pcube_timeseg CSV but was never included when
    building real_hr_csv -- exercises the 'HR row missing for this obsid' path
    on real, already-local data rather than a hand-built CSV."""
    csv_path = os.path.join(REAL_IXPE_DIR, NO_HR_OBSID, f"pcube_timeseg_{NO_HR_OBSID}.csv")
    if not os.path.exists(csv_path):
        pytest.skip(f"real fixture pcube_timeseg CSV not local at {csv_path}")
    obsids = {**_same_source(), NO_HR_OBSID: "Cyg X-1"}
    result = survey.pd_hr_survey(str(tmp_path), obsids, real_hr_csv, REAL_IXPE_DIR)
    assert result["exit_code"] == 0
    assert NO_HR_OBSID in result["skipped"]
    by_obsid = {r["obsid"]: r for r in result["rows"]}
    assert "no HR row" in by_obsid[NO_HR_OBSID]["note"]
    # PD/MDP99 were still read (the failure is HR-side only).
    assert not math.isnan(by_obsid[NO_HR_OBSID]["PD_pct"])
    assert neutral.read(str(tmp_path))["entries"].get(NO_HR_OBSID) is None


def test_empty_non_detected_group_yields_nan_mwu_not_exception(tmp_path, real_hr_csv):
    """Restrict the survey to only the three detected ObsIDs -- the
    non-detected group is empty -- Mann-Whitney U must come back NaN, never
    raise scipy's 'empty array' error."""
    only_detected = {o: "Cyg X-1" for o in DETECTED_OBSIDS}
    result = survey.pd_hr_survey(str(tmp_path), only_detected, real_hr_csv, REAL_IXPE_DIR)
    assert result["exit_code"] == 0
    assert result["mwu_result"]["n_detected"] == 3
    assert result["mwu_result"]["n_non_detected"] == 0
    assert math.isnan(result["mwu_result"]["statistic"])
    assert math.isnan(result["mwu_result"]["p"])


def test_empty_ks_groups_by_default_yield_nan_not_exception(tmp_path, real_hr_csv):
    """soft_state_sources/hard_state_sources default to empty -- the KS test
    must come back NaN, never raise, until the caller supplies both groups."""
    result = survey.pd_hr_survey(str(tmp_path), _same_source(), real_hr_csv, REAL_IXPE_DIR)
    assert result["ks_result"]["n_soft"] == 0 and result["ks_result"]["n_hard"] == 0
    assert math.isnan(result["ks_result"]["statistic"])
    assert math.isnan(result["ks_result"]["p"])


def test_mdp_flag_pct_flags_high_mdp99_obsids(tmp_path, real_hr_csv):
    """All four real integrated MDP99 values are < 15% (the default). Between
    them (real values ~0.50/0.59/0.66/0.82%) sits a threshold that flags some
    but not all -- proving mdp_flagged is genuinely data-driven, not a
    constant column."""
    result = survey.pd_hr_survey(
        str(tmp_path), _same_source(), real_hr_csv, REAL_IXPE_DIR, mdp_flag_pct=0.6)
    flagged = {r["obsid"]: r["mdp_flagged"] for r in result["rows"]}
    assert any(flagged.values())
    assert not all(flagged.values())


# --------------------------------------------------------------------------- #
# Reproducibility
# --------------------------------------------------------------------------- #

def test_repro_pd_hr_survey_stable_across_calls(tmp_path_factory, real_hr_csv):
    out_a = str(tmp_path_factory.mktemp("pd_hr_repro_a"))
    out_b = str(tmp_path_factory.mktemp("pd_hr_repro_b"))
    a = survey.pd_hr_survey(out_a, _two_groups(), real_hr_csv, REAL_IXPE_DIR,
                             min_n_intra=2, soft_state_sources=["GroupA"],
                             hard_state_sources=["GroupB"])
    b = survey.pd_hr_survey(out_b, _two_groups(), real_hr_csv, REAL_IXPE_DIR,
                             min_n_intra=2, soft_state_sources=["GroupA"],
                             hard_state_sources=["GroupB"])

    for key in ("rho", "p", "n"):
        assert a["spearman_combined"][key] == pytest.approx(b["spearman_combined"][key], rel=1e-10) \
            if not (isinstance(a["spearman_combined"][key], float) and math.isnan(a["spearman_combined"][key])) \
            else math.isnan(b["spearman_combined"][key])
    for key in ("statistic", "p"):
        assert a["ks_result"][key] == pytest.approx(b["ks_result"][key], rel=1e-10)
        assert a["mwu_result"][key] == pytest.approx(b["mwu_result"][key], rel=1e-10)
    assert set(a["per_source_correlations"]) == set(b["per_source_correlations"])
    for source in a["per_source_correlations"]:
        ra = a["per_source_correlations"][source]["rho"]
        rb = b["per_source_correlations"][source]["rho"]
        if isinstance(ra, float) and math.isnan(ra):
            assert math.isnan(rb)
        else:
            assert ra == pytest.approx(rb, rel=1e-10)


# --------------------------------------------------------------------------- #
# Catalog recording
# --------------------------------------------------------------------------- #

def test_catalog_records_one_product_per_obsid_with_both_pd_and_hr(tmp_path, real_hr_csv):
    obsids = {**_same_source(), "99999999": "Cyg X-1"}
    result = survey.pd_hr_survey(str(tmp_path), obsids, real_hr_csv, REAL_IXPE_DIR)

    entries = neutral.read(str(tmp_path))["entries"]
    for obsid in ALL_OBSIDS:
        entry = entries[obsid]
        assert entry["max_level"] == 3
        products = [p for p in entry["products"] if p["kind"] == "pd_hardness_survey"]
        assert len(products) == 1
        p = products[0]
        assert p["produced_by"] == "xspec_utils.lib.bhxrb_pd_hr_survey.pd_hr_survey"
        assert not os.path.isabs(p["path"])
        assert os.path.exists(os.path.join(str(tmp_path), p["path"]))
        assert p["params"]["source"] == "Cyg X-1"
    assert "99999999" not in entries


# --------------------------------------------------------------------------- #
# Manifest wiring + no-science-stack import guard
# --------------------------------------------------------------------------- #

def test_tool_manifest_wired_to_lib_directly_no_science_stack():
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.bhxrb_pd_hr_survey"]
    assert entry["module"] == "xspec_utils.lib.bhxrb_pd_hr_survey"
    assert entry["fn"] == "pd_hr_survey"
    props = entry["args_schema"]["properties"]
    for name in ("obs_db_path", "obsids_with_sources", "hr_csv", "pcube_csv_dir",
                 "mdp_flag_pct", "min_n_intra", "soft_state_sources", "hard_state_sources"):
        assert name in props, f"{name} missing from args_schema.properties"
    assert entry["args_schema"]["required"] == [
        "obs_db_path", "obsids_with_sources", "hr_csv", "pcube_csv_dir"]
    assert "15.0" in props["mdp_flag_pct"]["description"]


def test_lib_module_imports_no_xspec_or_ixpeobssim():
    """Static guard: the module must not *import* the science stack (mentions
    in comments/docstrings explaining the design are fine)."""
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
    assert "pa_pd_stability" not in imported
