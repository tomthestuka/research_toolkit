"""Tests for xspec_utils' ixpe_pa_pd_stability (lib/pa_pd_stability.py).

Uses the six real, already-local pcube_timeseg_<obsid>.csv products (already
produced by ixpe_pcube_timeseg -- see red_ixpe.pcube_timeseg /
tests/test_ixpe_pcube_timeseg.py) directly from
/data/fzhou/XSPEC/Cyg_X-1/IXPE. pa_pd_stability only *reads* those CSVs and
writes its own catalog/summary to a separate output_dir, so unlike the
pcube_timeseg tests there is no need to copy event files into scratch -- the
real obs_db_path is used as-is (read-only) and tmp_path is used only for
output_dir.

The "variable vs constant" capability requirement is covered two ways.
First, on real data: 03006001 lands p_2D just under the default
variability_p_threshold (genuinely variable at that threshold) and 01002901
lands p_2D very far above it (genuinely constant) -- both verified directly
against the real CSVs before writing this test, so the discriminator is
proven on real, local measurements. Second, on a controlled fixture: two
pcube_timeseg CSVs (`synth_pa_pd_variable`, `synth_pa_pd_constant`) are
written directly into tmp_path as test-local scaffolding -- the same pattern
`test_read_timeseg_csv_round_trip` already uses -- with QI/UI values chosen
by hand so the discriminator's verdict is pinned by construction rather than
by a real reduction's numbers, which could drift if the reduction is ever
regenerated.
"""

from __future__ import annotations

import csv
import math
import os

import numpy as np
import pytest
from scipy import stats as scipy_stats

from xspec_utils.lib import pa_pd_stability as pps
from xspec_utils.skills import TOOL_MANIFEST

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
ALL_OBSIDS = ["01002901", "01250401", "02001901", "03002599", "03005801", "03006001"]
CONSTANT_OBSID = "01002901"   # p_2D ~ 0.96 at the real, local data -- far above threshold
VARIABLE_OBSID = "03006001"   # p_2D ~ 0.0079 at the real, local data -- below default 0.01


def _require_real_csv(obsid: str) -> str:
    path = os.path.join(REAL_IXPE_DIR, obsid, f"pcube_timeseg_{obsid}.csv")
    if not os.path.exists(path):
        pytest.skip(f"real fixture pcube_timeseg CSV not local at {path}")
    return path


# --------------------------------------------------------------------------- #
# Pure-logic unit tests: wrap/weighted-stat/moment-scatter/structure-function
# --------------------------------------------------------------------------- #

def test_wrap180_no_op_inside_range():
    assert pps._wrap180(30.0) == pytest.approx(30.0)
    assert pps._wrap180(-45.0) == pytest.approx(-45.0)


def test_wrap180_handles_boundary_flip():
    # +89 and -89 are 2 degrees apart across the +/-90 branch cut.
    d = 89.0 - (-89.0)  # naive difference: 178
    assert pps._wrap180(d) == pytest.approx(-2.0)


def test_weighted_stat_hand_computed():
    # Two points, hand-computed weighted mean/chi2: x=[1,3], err=[1,1] ->
    # weights [1,1], mean=2, chi2 = 1*(1-2)^2 + 1*(3-2)^2 = 2, dof=1.
    r = pps._weighted_stat(np.array([1.0, 3.0]), np.array([1.0, 1.0]))
    assert r["mean"] == pytest.approx(2.0)
    assert r["chi2"] == pytest.approx(2.0)
    assert r["dof"] == 1
    assert r["p"] == pytest.approx(scipy_stats.chi2.sf(2.0, 1))


def test_weighted_stat_excludes_nonfinite_and_nonpositive_errors():
    x = np.array([1.0, 3.0, float("nan"), 5.0])
    err = np.array([1.0, 1.0, 1.0, -1.0])
    r = pps._weighted_stat(x, err)
    assert r["n"] == 2
    assert r["mean"] == pytest.approx(2.0)


def test_weighted_stat_single_point_chi2_nan():
    r = pps._weighted_stat(np.array([7.0]), np.array([1.0]))
    assert r["n"] == 1
    assert r["mean"] == pytest.approx(7.0)
    assert math.isnan(r["chi2"])
    assert r["dof"] == 0


def test_weighted_stat_pa_wraparound_gives_small_chi2_not_large():
    # Two segments that are physically ~2 deg apart (89, -89) must NOT be
    # scored as a ~178 deg discrepancy: the doubled-angle circular mean and
    # wrapped residuals must keep chi2 small.
    r = pps._weighted_stat_pa(np.array([89.0, -89.0]), np.array([5.0, 5.0]))
    assert r["chi2"] < 1.0, "PA wraparound not handled -- chi2 should be tiny for a ~2 deg true difference"


def test_weighted_stat_pa_matches_naive_when_no_wraparound():
    # Away from the +/-90 branch cut, the circular estimator must reduce to
    # the ordinary weighted mean/chi2 (independent hand computation).
    pa = np.array([10.0, 20.0])
    err = np.array([2.0, 2.0])
    r = pps._weighted_stat_pa(pa, err)
    w = 1.0 / err ** 2
    mean_naive = float(np.sum(w * pa) / np.sum(w))
    chi2_naive = float(np.sum(w * (pa - mean_naive) ** 2))
    assert r["mean"] == pytest.approx(mean_naive, abs=1e-6)
    assert r["chi2"] == pytest.approx(chi2_naive, abs=1e-6)


def test_moment_intrinsic_scatter_upper_limit_when_noise_dominates():
    # Values with zero real scatter (all identical) but nonzero measurement
    # error -> sigma_obs^2 = 0 < mean_var -> excess <= 0 -> upper limit.
    values = np.array([5.0, 5.0, 5.0])
    errs = np.array([1.0, 1.0, 1.0])
    mask = np.array([True, True, True])
    sigma_int, upper = pps._moment_intrinsic_scatter(values, errs, mask)
    assert sigma_int == pytest.approx(0.0)
    assert upper is True


def test_moment_intrinsic_scatter_detects_real_scatter():
    # Sample variance well above the (tiny) measurement variance -> a
    # positive, non-upper-limit sigma_int recovered via hand-computed check.
    values = np.array([0.0, 10.0, 20.0])
    errs = np.array([0.1, 0.1, 0.1])
    mask = np.array([True, True, True])
    sigma_int, upper = pps._moment_intrinsic_scatter(values, errs, mask)
    sigma_obs2 = float(np.var(values, ddof=1))
    mean_var = float(np.mean(errs ** 2))
    assert sigma_int == pytest.approx(math.sqrt(sigma_obs2 - mean_var), rel=1e-9)
    assert upper is False


def test_structure_function_known_pairs():
    # t=[0,5,10] (seconds), x=[1,2,4]; pairs: (0,5)->dt=5,dx=1 ; (5,10)->dt=5,dx=2 ;
    # (0,10)->dt=10,dx=3. Lag 5s (half-width 1s) averages the two dt=5 pairs:
    # mean(1^2, 2^2) = 2.5. Lag 10s picks up only the dt=10 pair: 3^2=9.
    # Lag 100s has no qualifying pair -> NaN.
    sf = pps.structure_function([0.0, 5.0, 10.0], [1.0, 2.0, 4.0], [5.0, 10.0, 100.0], half_width_s=1.0)
    assert sf[5.0] == pytest.approx(2.5)
    assert sf[10.0] == pytest.approx(9.0)
    assert math.isnan(sf[100.0])


def test_structure_function_pa_wraps_before_squaring():
    sf = pps.structure_function([0.0, 5.0], [89.0, -89.0], [5.0], half_width_s=1.0, wrap=True)
    assert sf[5.0] == pytest.approx(4.0)  # (-2 deg)^2, not (178 deg)^2


# --------------------------------------------------------------------------- #
# Self-contained CSV reader: round-trip + missing-header safety
# --------------------------------------------------------------------------- #

def test_read_timeseg_csv_round_trip(tmp_path):
    csv_path = str(tmp_path / "pcube_timeseg_roundtrip.csv")
    provenance = {"segment_duration_s": 5000.0, "tstart": 0.0, "tstop": 5000.0}
    with open(csv_path, "w", newline="") as f:
        f.write(pps._PROVENANCE_HEADER_PREFIX + '{"segment_duration_s": 5000.0, "tstart": 0.0, "tstop": 5000.0}\n')
        writer = csv.DictWriter(f, fieldnames=pps._TIMESEG_FIELDS)
        writer.writeheader()
        row0 = {k: 0.0 for k in pps._TIMESEG_FIELDS}
        row0.update(segment_id=0, PD_pct=8.0, PA_deg=17.0)
        writer.writerow(row0)
        row_int = {k: 0.0 for k in pps._TIMESEG_FIELDS}
        row_int.update(segment_id=-1, PD_pct=8.1, PA_deg=17.1)
        writer.writerow(row_int)

    result = pps._read_timeseg_csv(csv_path)
    assert result is not None
    read_provenance, rows_by_segid, integrated = result
    assert read_provenance == provenance
    assert set(rows_by_segid) == {0}
    assert rows_by_segid[0]["PD_pct"] == pytest.approx(8.0)
    assert integrated["PD_pct"] == pytest.approx(8.1)


def test_read_timeseg_csv_missing_header_treated_as_absent(tmp_path):
    csv_path = str(tmp_path / "pcube_timeseg_no_header.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=pps._TIMESEG_FIELDS)
        writer.writeheader()
    assert pps._read_timeseg_csv(csv_path) is None


def test_read_timeseg_csv_missing_file_returns_none(tmp_path):
    assert pps._read_timeseg_csv(str(tmp_path / "does_not_exist.csv")) is None


def test_read_timeseg_csv_differential_agreement_with_red_ixpe_reader():
    """Independent-implementation differential: pa_pd_stability's own CSV
    reader must agree with red_ixpe._read_existing_timeseg_csv (a genuinely
    separate implementation of the same file format -- pa_pd_stability never
    imports red_ixpe) on a real pcube_timeseg CSV."""
    from xspec_utils.lib import red_ixpe  # science-stack import, test-only

    csv_path = _require_real_csv(CONSTANT_OBSID)
    ours = pps._read_timeseg_csv(csv_path)
    theirs = red_ixpe._read_existing_timeseg_csv(csv_path)
    assert ours is not None and theirs is not None
    our_prov, our_rows, our_int = ours
    their_prov, their_rows, their_int = theirs
    assert our_prov == their_prov
    assert set(our_rows) == set(their_rows)
    for segid in our_rows:
        assert our_rows[segid] == their_rows[segid]
    assert our_int == their_int


# --------------------------------------------------------------------------- #
# Manifest wiring (AI-facing half)
# --------------------------------------------------------------------------- #

def test_tool_manifest_wired_to_lib_directly_no_science_stack():
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.ixpe_pa_pd_stability"]
    assert entry["module"] == "xspec_utils.lib.pa_pd_stability"
    assert entry["fn"] == "pa_pd_stability"
    props = entry["args_schema"]["properties"]
    for name in ("obs_db_path", "obsids", "output_dir", "sources",
                 "variability_p_threshold", "sf_lags_ks", "min_seg_for_sf"):
        assert name in props, f"{name} missing from args_schema.properties"
    assert entry["args_schema"]["required"] == ["obs_db_path", "obsids", "output_dir"]
    assert "0.01" in props["variability_p_threshold"]["description"]


def test_lib_module_imports_no_xspec_or_ixpeobssim():
    """Static guard: the module must not *import* the science stack (mentions
    in comments/docstrings explaining the design are fine) -- the whole point
    of keeping this analysis tool importable without HEASoft."""
    import ast
    import inspect
    tree = ast.parse(inspect.getsource(pps))
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    assert "xspec" not in imported
    assert "ixpeobssim" not in imported
    assert "red_ixpe" not in imported


# --------------------------------------------------------------------------- #
# Real-data capability + differential + reproducibility
# --------------------------------------------------------------------------- #

def test_capability_discriminates_variable_from_constant_on_real_obsids(tmp_path):
    """Proves the 2D Stokes test actually consumes the segment values (not
    just returns well-shaped metadata): a genuinely constant real ObsID must
    land p_2D above variability_p_threshold, and a genuinely variable one
    below it -- both at the tool's default threshold=0.01."""
    _require_real_csv(CONSTANT_OBSID)
    _require_real_csv(VARIABLE_OBSID)

    result = pps.pa_pd_stability(
        REAL_IXPE_DIR, [CONSTANT_OBSID, VARIABLE_OBSID], str(tmp_path))
    assert result["exit_code"] == 0

    const_stats = result["per_obs"][CONSTANT_OBSID]
    var_stats = result["per_obs"][VARIABLE_OBSID]
    assert const_stats["error"] is None and var_stats["error"] is None

    assert const_stats["p_2D"] > 0.5, "known-constant real ObsID must NOT be flagged variable"
    assert const_stats["variable"] is False

    assert var_stats["p_2D"] < 0.01
    assert var_stats["variable"] is True


def test_capability_discriminates_on_controlled_fixture_csvs(tmp_path):
    """Controlled, deterministic capability grader: build two pcube_timeseg
    CSVs directly in tmp_path -- a genuinely VARIABLE one (QI/UI values that
    scatter far beyond their errors) and a genuinely CONSTANT one (QI/UI
    identical within their errors) -- and assert the authoritative 2D Stokes
    test flips `variable` True for the first and False for the second at the
    default threshold. Proves the discriminator consumes the segment values,
    with no dependence on any real reduction's drifting p-value."""

    def _write_csv(db_dir, obsid, qi_vals, ui_vals, err):
        obsdir = db_dir / obsid
        obsdir.mkdir(parents=True, exist_ok=True)
        path = obsdir / f"pcube_timeseg_{obsid}.csv"
        with open(path, "w", newline="") as f:
            f.write(pps._PROVENANCE_HEADER_PREFIX +
                    '{"segment_duration_s": 5000.0}\n')
            writer = csv.DictWriter(f, fieldnames=pps._TIMESEG_FIELDS)
            writer.writeheader()
            for i, (qi, ui) in enumerate(zip(qi_vals, ui_vals)):
                row = {k: 0.0 for k in pps._TIMESEG_FIELDS}
                row.update(segment_id=i,
                           t_start_mjd=60000.0 + i * 5000.0 / 86400.0,
                           t_stop_mjd=60000.0 + (i + 1) * 5000.0 / 86400.0,
                           exposure_s=5000.0,
                           PD_pct=1.0, PD_err=err, PA_deg=10.0, PA_err=err,
                           QI=qi, QI_err=err, UI=ui, UI_err=err,
                           MDP99=1.0, SIGNIF=5.0)
                writer.writerow(row)
            # integrated row (segment_id == -1) so the reader/labels are exercised
            row = {k: 0.0 for k in pps._TIMESEG_FIELDS}
            row.update(segment_id=-1, PD_pct=1.0, PD_err=err, PA_deg=10.0,
                       PA_err=err, QI=float(np.mean(qi_vals)), QI_err=err,
                       UI=float(np.mean(ui_vals)), UI_err=err, MDP99=1.0, SIGNIF=5.0)
            writer.writerow(row)
        return path

    db = tmp_path / "db"
    err = 0.01
    # VARIABLE: QI/UI swing by ~10 sigma between segments
    _write_csv(db, "synth_pa_pd_variable",
               qi_vals=[0.10, 0.20, 0.10, 0.20, 0.10, 0.20],
               ui_vals=[0.05, 0.15, 0.05, 0.15, 0.05, 0.15], err=err)
    # CONSTANT: QI/UI identical to within << 1 sigma
    _write_csv(db, "synth_pa_pd_constant",
               qi_vals=[0.10, 0.10, 0.10, 0.10, 0.10, 0.10],
               ui_vals=[0.05, 0.05, 0.05, 0.05, 0.05, 0.05], err=err)

    out = pps.pa_pd_stability(str(db),
                              ["synth_pa_pd_variable", "synth_pa_pd_constant"],
                              str(tmp_path / "out"))
    var = out["per_obs"]["synth_pa_pd_variable"]
    con = out["per_obs"]["synth_pa_pd_constant"]
    assert var["error"] is None and con["error"] is None

    assert var["p_2D"] < 0.01 and var["variable"] is True
    assert con["p_2D"] > 0.5 and con["variable"] is False
    # And the discriminator must genuinely differ, not both land the same:
    assert var["variable"] != con["variable"]


def test_differential_chi2_pd_independent_numpy_recompute(tmp_path):
    """Independent differential (NOT calling the tool's own _weighted_stat
    helper): recompute the weighted mean and chi-squared for PD directly from
    the real CSV with a standalone numpy expression, and check agreement with
    the tool's chi2_PD for the same ObsID."""
    csv_path = _require_real_csv(CONSTANT_OBSID)
    with open(csv_path, newline="") as f:
        f.readline()
        rows = [r for r in csv.DictReader(f) if int(r["segment_id"]) >= 0]
    pd = np.array([float(r["PD_pct"]) for r in rows])
    pd_err = np.array([float(r["PD_err"]) for r in rows])
    mask = np.isfinite(pd) & np.isfinite(pd_err) & (pd_err > 0)
    w = np.reciprocal(np.square(pd_err[mask]))
    xbar = np.dot(w, pd[mask]) * np.reciprocal(np.sum(w))
    chi2_independent = np.dot(w, np.square(pd[mask] - xbar))

    result = pps.pa_pd_stability(REAL_IXPE_DIR, [CONSTANT_OBSID], str(tmp_path))
    chi2_tool = result["per_obs"][CONSTANT_OBSID]["PD"]["chi2"]
    assert float(chi2_independent) == pytest.approx(chi2_tool, rel=1e-8)


def test_repro_pa_pd_stability_stable_across_calls(tmp_path_factory):
    """Reproducibility: the tool is pure arithmetic over an unchanged, real,
    already-local input CSV -- two independent calls must agree near-exactly."""
    _require_real_csv(CONSTANT_OBSID)
    out_a = str(tmp_path_factory.mktemp("pa_pd_repro_a"))
    out_b = str(tmp_path_factory.mktemp("pa_pd_repro_b"))
    a = pps.pa_pd_stability(REAL_IXPE_DIR, [CONSTANT_OBSID], out_a)
    b = pps.pa_pd_stability(REAL_IXPE_DIR, [CONSTANT_OBSID], out_b)

    sa, sb = a["per_obs"][CONSTANT_OBSID], b["per_obs"][CONSTANT_OBSID]
    for key in ("chi2_2D", "dof_2D", "p_2D", "sigma_PD_int", "sigma_PA_int",
                "PD_ptp", "PA_ptp", "rchi2_PD", "rchi2_PA", "mean_MDP99"):
        va, vb = sa[key], sb[key]
        if isinstance(va, float) and math.isnan(va):
            assert math.isnan(vb)
        else:
            assert va == pytest.approx(vb, rel=1e-10, abs=1e-12)
    for quantity in ("PD", "PA", "QI", "UI"):
        for field in ("mean", "chi2", "dof", "p", "n"):
            assert sa[quantity][field] == pytest.approx(sb[quantity][field], rel=1e-10, abs=1e-12) \
                if not (isinstance(sa[quantity][field], float) and math.isnan(sa[quantity][field])) \
                else math.isnan(sb[quantity][field])
    assert sa["variable"] == sb["variable"]


# --------------------------------------------------------------------------- #
# Missing-CSV edge case -- must fail loudly? No: recorded, never raised.
# --------------------------------------------------------------------------- #

def test_missing_csv_recorded_as_error_and_skipped_not_raised(tmp_path):
    result = pps.pa_pd_stability(REAL_IXPE_DIR, ["99999999"], str(tmp_path))
    assert result["exit_code"] == 0
    assert result["per_obs"]["99999999"]["error"] is not None
    with open(result["catalog_csv"], newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1
    assert rows[0]["obsid"] == "99999999"
    assert rows[0]["error"] != ""
    assert rows[0]["variable"] == "False"


def test_mixed_missing_and_real_obsid_only_real_one_contributes(tmp_path):
    _require_real_csv(CONSTANT_OBSID)
    result = pps.pa_pd_stability(REAL_IXPE_DIR, [CONSTANT_OBSID, "99999999"], str(tmp_path))
    assert result["per_obs"][CONSTANT_OBSID]["error"] is None
    assert result["per_obs"]["99999999"]["error"] is not None
    with open(result["catalog_csv"], newline="") as f:
        rows = {r["obsid"]: r for r in csv.DictReader(f)}
    assert len(rows) == 2
    assert rows[CONSTANT_OBSID]["N_seg"] != "0"
    assert rows["99999999"]["N_seg"] == "0"


# --------------------------------------------------------------------------- #
# sources labeling + cross-source Spearman across all six real ObsIDs
# --------------------------------------------------------------------------- #

def test_sources_label_flows_into_catalog(tmp_path):
    _require_real_csv(CONSTANT_OBSID)
    result = pps.pa_pd_stability(
        REAL_IXPE_DIR, [CONSTANT_OBSID], str(tmp_path), sources={CONSTANT_OBSID: "Cyg X-1"})
    assert result["per_obs"][CONSTANT_OBSID]["source"] == "Cyg X-1"
    with open(result["catalog_csv"], newline="") as f:
        rows = list(csv.DictReader(f))
    assert rows[0]["source"] == "Cyg X-1"


def test_cross_source_spearman_across_six_real_obsids(tmp_path):
    for obsid in ALL_OBSIDS:
        _require_real_csv(obsid)
    result = pps.pa_pd_stability(REAL_IXPE_DIR, ALL_OBSIDS, str(tmp_path))
    assert all(result["per_obs"][o]["error"] is None for o in ALL_OBSIDS)

    for key in ("mdp99_vs_rchi2_pd", "mdp99_vs_rchi2_pa"):
        res = result["spearman"][key]
        assert res["n"] == 6
        assert res["note"] is None
        assert -1.0 <= res["rho"] <= 1.0
        assert 0.0 <= res["p"] <= 1.0

    assert os.path.exists(result["summary_txt"])
    with open(result["summary_txt"]) as f:
        summary = f.read()
    for obsid in ALL_OBSIDS:
        assert obsid in summary
    assert "Spearman" in summary


def test_structure_function_computed_only_above_min_seg_for_sf(tmp_path):
    """03002599 has only 44 segments (n2d above); with min_seg_for_sf set
    above that count, every SF lag must come back NaN -- proving the gate,
    not just checking the default behaviour."""
    _require_real_csv("03002599")
    result = pps.pa_pd_stability(
        REAL_IXPE_DIR, ["03002599"], str(tmp_path), min_seg_for_sf=1000)
    s = result["per_obs"]["03002599"]
    assert all(math.isnan(v) for v in s["sf_pd"].values())
    assert all(math.isnan(v) for v in s["sf_pa"].values())

    result2 = pps.pa_pd_stability(
        REAL_IXPE_DIR, ["03002599"], str(tmp_path), min_seg_for_sf=10, sf_lags_ks=[5, 10, 20])
    s2 = result2["per_obs"]["03002599"]
    assert set(s2["sf_pd"].keys()) == {5000.0, 10000.0, 20000.0}
    # At least one lag should find a qualifying pair with 44 real segments.
    assert any(not math.isnan(v) for v in s2["sf_pd"].values())
