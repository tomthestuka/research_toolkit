"""Tests for the abs(tau_peak) population correlation pairs added to
`qu_lag_dcf` (`lib/qu_lag_dcf.py`).

The census's population step correlates the peak q(t)-u(t) cross-correlation
lag against a source-brightness proxy (`r_proxy`) and a combined PD
significance (`pd_significance`). The physical prescription names the
MAGNITUDE of the lag, `abs(tau_peak_s)` -- how long the lag is -- while the
pre-existing pair correlated the SIGNED `tau_peak_s`, which conflates lag
length with lag SENSE (q-leads-u vs u-leads-q). This adds two new pairs,
`abs_tau_peak_vs_r_proxy` and `abs_tau_peak_vs_pd_significance`, alongside the
two pre-existing signed pairs -- purely additive: `POPULATION_FIELDS`, row
count of qu_dcf/ccf_master outputs, and the two signed pairs' names/values are
all unchanged.

Uses four already-local qu_slices_<subobs_key>.csv products, copied into an
isolated tmp obs_db_path per test so the tool's own writes (qu_dcf CSVs,
population_correlations.csv, catalog.json) never touch the shared
analyst-env mission directory:

- 01002901, 02250301: real, already-local Cyg X-1 IXPE ObsIDs.
- synth_qu_lag_evensamp01, synth_qu_lag_null01: marked-synthetic ObsIDs
  (ground truth read from their own `qu_slice_pcubes` provenance header,
  never hardcoded).

Combined, these four give MIXED-SIGN `tau_peak_s` across `status=='ok'`
subobs (confirmed empirically: 01002901/02250301/synth_qu_lag_null01 recover
a negative peak lag, synth_qu_lag_evensamp01 a positive one) -- exactly the
condition the brief requires to demonstrate the abs and signed statistics
are genuinely different computations, not well-shaped duplicates.
"""

from __future__ import annotations

import csv
import os
import shutil

import numpy as np
import pytest
from scipy import stats

from xspec_utils.lib import qu_lag_dcf as m

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
OBSID_A = "01002901"
OBSID_B = "02250301"
OBSID_EVENSAMP = "synth_qu_lag_evensamp01"
OBSID_NULL = "synth_qu_lag_null01"
MIXED_SIGN_OBSIDS = [OBSID_A, OBSID_B, OBSID_EVENSAMP, OBSID_NULL]

ABS_PAIRS = ("abs_tau_peak_vs_r_proxy", "abs_tau_peak_vs_pd_significance")
SIGNED_PAIRS = ("tau_peak_vs_r_proxy", "tau_peak_vs_pd_significance")
ALL_PAIRS = set(SIGNED_PAIRS) | set(ABS_PAIRS)


def _copy_qu_slices(dest_dir: str, obsid: str) -> str:
    fname = "qu_slices_" + obsid + ".csv"
    src_csv = os.path.join(REAL_IXPE_DIR, obsid, fname)
    if not os.path.exists(src_csv):
        pytest.skip(f"fixture qu_slices CSV not local at {src_csv}")
    dst_dir = os.path.join(dest_dir, obsid)
    os.makedirs(dst_dir, exist_ok=True)
    dst_csv = os.path.join(dst_dir, fname)
    shutil.copyfile(src_csv, dst_csv)
    return dst_csv


def _read_ccf_master_raw(path: str) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _read_population_csv(path: str) -> dict:
    with open(path, newline="") as f:
        return {r["pair"]: r for r in csv.DictReader(f)}


def _zero_good_exposure(csv_path: str) -> None:
    """Rewrite every data row's `good_exposure_s` cell to `0.0` in place,
    leaving the provenance header line, the column-header line, and every
    other column (`counts`, `q`, `q_err`, `u`, `u_err`, the times,
    `under_target`) byte-for-byte untouched. `r_proxy_from_rows` returns NaN
    for non-positive total exposure, so this drives exactly that subobs's
    `r_proxy` to NaN while leaving its `pd_significance` (which never reads
    `good_exposure_s`) finite."""
    with open(csv_path, newline="") as f:
        header_line = f.readline()
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)
    for row in rows:
        row["good_exposure_s"] = "0.0"
    with open(csv_path, "w", newline="") as f:
        f.write(header_line)
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


# --------------------------------------------------------------------------- #
# Capability: mixed-sign tau_peak_s -> abs pairs present and genuinely
# different from the signed pairs; signed pairs numerically unchanged.
# --------------------------------------------------------------------------- #

def test_capability_abs_pairs_present_and_differ_from_signed(tmp_path):
    for obsid in MIXED_SIGN_OBSIDS:
        _copy_qu_slices(str(tmp_path), obsid)
    out_dir = os.path.join(str(tmp_path), "out")

    res = m.qu_lag_dcf(str(tmp_path), MIXED_SIGN_OBSIDS, out_dir,
                        n_permutations=50, n_bootstrap=50, seed=0)
    assert res["exit_code"] == 0
    pop = res["population_correlations"]

    # All four rows present, in both the returned dict and the CSV.
    assert set(pop) == ALL_PAIRS
    csv_pairs = _read_population_csv(res["population_csv"])
    assert set(csv_pairs) == ALL_PAIRS

    # Confirm the fixture set genuinely carries mixed-sign tau_peak_s, so this
    # is not a vacuous check.
    master = _read_ccf_master_raw(res["ccf_master_csv"])
    ok_taus = [float(r["tau_peak_s"]) for r in master if r["status"] == "ok"]
    assert any(t > 0 for t in ok_taus) and any(t < 0 for t in ok_taus)

    # The abs pair must be a genuinely different computation from the signed
    # pair sharing the same covariate -- not a well-shaped duplicate.
    assert pop["abs_tau_peak_vs_pd_significance"]["rho"] != pytest.approx(
        pop["tau_peak_vs_pd_significance"]["rho"], abs=1e-9)
    assert pop["abs_tau_peak_vs_r_proxy"]["rho"] != pytest.approx(
        pop["tau_peak_vs_r_proxy"]["rho"], abs=1e-9)

    # Every pair reports the same n (all four ok rows survive both filters).
    for pair in ALL_PAIRS:
        assert pop[pair]["n"] == 4

    # The two pre-existing signed rows are numerically unchanged from what an
    # independent scipy.stats.spearmanr computes directly on the signed
    # tau_peak_s -- i.e. the pre-change behaviour for these two pairs.
    ok_rows = [r for r in master if r["status"] == "ok"]
    tau_signed = np.array([float(r["tau_peak_s"]) for r in ok_rows])
    r_proxy = np.array([float(r["r_proxy"]) for r in ok_rows])
    pd_sig = np.array([float(r["pd_significance"]) for r in ok_rows])
    expected_r = stats.spearmanr(tau_signed, r_proxy)
    expected_pd = stats.spearmanr(tau_signed, pd_sig)
    assert pop["tau_peak_vs_r_proxy"]["rho"] == pytest.approx(expected_r.statistic, abs=1e-12)
    assert pop["tau_peak_vs_r_proxy"]["p"] == pytest.approx(expected_r.pvalue, abs=1e-12)
    assert pop["tau_peak_vs_pd_significance"]["rho"] == pytest.approx(expected_pd.statistic, abs=1e-12)
    assert pop["tau_peak_vs_pd_significance"]["p"] == pytest.approx(expected_pd.pvalue, abs=1e-12)


# --------------------------------------------------------------------------- #
# Differential: independent scipy.stats.spearmanr recompute from ccf_master
# .csv on disk (never the module's own _spearman), for both abs and signed
# pairs.
# --------------------------------------------------------------------------- #

def test_differential_population_matches_independent_spearmanr(tmp_path):
    for obsid in MIXED_SIGN_OBSIDS:
        _copy_qu_slices(str(tmp_path), obsid)
    out_dir = os.path.join(str(tmp_path), "out")

    res = m.qu_lag_dcf(str(tmp_path), MIXED_SIGN_OBSIDS, out_dir,
                        n_permutations=50, n_bootstrap=50, seed=1)

    master = _read_ccf_master_raw(res["ccf_master_csv"])
    ok_rows = [r for r in master if r["status"] == "ok"]
    tau_signed = np.array([float(r["tau_peak_s"]) for r in ok_rows])
    tau_abs = np.abs(tau_signed)
    r_proxy = np.array([float(r["r_proxy"]) for r in ok_rows])
    pd_sig = np.array([float(r["pd_significance"]) for r in ok_rows])

    pop = _read_population_csv(res["population_csv"])

    for series, pair_signed, pair_abs, covariate in (
        (tau_signed, "tau_peak_vs_r_proxy", "abs_tau_peak_vs_r_proxy", r_proxy),
        (tau_signed, "tau_peak_vs_pd_significance", "abs_tau_peak_vs_pd_significance", pd_sig),
    ):
        expected_signed = stats.spearmanr(series, covariate)
        assert float(pop[pair_signed]["rho"]) == pytest.approx(expected_signed.statistic, abs=1e-12)
        assert float(pop[pair_signed]["p"]) == pytest.approx(expected_signed.pvalue, abs=1e-12)

        expected_abs = stats.spearmanr(tau_abs, covariate)
        assert float(pop[pair_abs]["rho"]) == pytest.approx(expected_abs.statistic, abs=1e-12)
        assert float(pop[pair_abs]["p"]) == pytest.approx(expected_abs.pvalue, abs=1e-12)


# --------------------------------------------------------------------------- #
# Differential (exclude_noise_limited path): 02250301 truncated to 6 data
# rows sits between a lowered n_min_slices=5 (so status=='ok') and the
# provenance's fixed n_min=8 (so noise_limited=True) -- the one condition
# under which exclude_noise_limited actually removes a row from the ok set.
# Confirms the abs pairs are filtered identically to the signed pairs on
# that same reduced row set, and that the filter is not silently a no-op.
# --------------------------------------------------------------------------- #

def test_noise_limited_subobs_filtered_identically_from_all_four_pairs(tmp_path):
    for obsid in MIXED_SIGN_OBSIDS:
        _copy_qu_slices(str(tmp_path), obsid)

    # Truncate 02250301 (25 data rows) to 6 data rows: keep the provenance
    # header line and the column-header line byte-for-byte, drop the rest.
    truncated_path = os.path.join(str(tmp_path), OBSID_B, f"qu_slices_{OBSID_B}.csv")
    with open(truncated_path) as f:
        lines = f.readlines()
    with open(truncated_path, "w") as f:
        f.writelines(lines[:2 + 6])  # provenance header + column header + 6 data rows

    res_true = m.qu_lag_dcf(
        str(tmp_path), MIXED_SIGN_OBSIDS, os.path.join(str(tmp_path), "out_true"),
        n_min_slices=5, n_permutations=50, n_bootstrap=50, seed=3,
        exclude_noise_limited=True)
    res_false = m.qu_lag_dcf(
        str(tmp_path), MIXED_SIGN_OBSIDS, os.path.join(str(tmp_path), "out_false"),
        n_min_slices=5, n_permutations=50, n_bootstrap=50, seed=3,
        exclude_noise_limited=False)

    # Non-vacuity: the truncated subobs must genuinely be 'ok' (6 >= 5) AND
    # noise_limited (6 < the provenance's own n_min=8) in BOTH runs, and every
    # other subobs must be 'ok' and NOT noise_limited -- fail loudly, never skip.
    for res in (res_true, res_false):
        truncated_info = res["per_subobs"][OBSID_B]
        assert truncated_info["status"] == "ok", (
            f"{OBSID_B} must be status=='ok' at n_min_slices=5 with 6 rows; got {truncated_info}")
        assert truncated_info["noise_limited"] is True, (
            f"{OBSID_B} must be noise_limited (6 rows < provenance n_min=8); got {truncated_info}")
        for obsid in (OBSID_A, OBSID_EVENSAMP, OBSID_NULL):
            info = res["per_subobs"][obsid]
            assert info["status"] == "ok", f"{obsid} must be status=='ok'; got {info}"
            assert info["noise_limited"] is False, f"{obsid} must not be noise_limited; got {info}"

    pop_true = res_true["population_correlations"]
    pop_false = res_false["population_correlations"]
    for pair in ALL_PAIRS:
        assert pop_true[pair]["n"] == 3, f"{pair} (exclude_noise_limited=True) expected n=3, got {pop_true[pair]}"
        assert pop_false[pair]["n"] == 4, f"{pair} (exclude_noise_limited=False) expected n=4, got {pop_false[pair]}"

    # Independent differential over the FILTERED set (exclude_noise_limited=True):
    # read ccf_master.csv back from disk, keep status=='ok' rows, drop the
    # truncated obsid's row by obsid, and call scipy.stats.spearmanr directly
    # (never the module's own _spearman).
    master_true = _read_ccf_master_raw(res_true["ccf_master_csv"])
    filtered_rows = [r for r in master_true if r["status"] == "ok" and r["obsid"] != OBSID_B]
    assert len(filtered_rows) == 3
    tau_signed = np.array([float(r["tau_peak_s"]) for r in filtered_rows])
    tau_abs = np.abs(tau_signed)
    r_proxy = np.array([float(r["r_proxy"]) for r in filtered_rows])
    pd_sig = np.array([float(r["pd_significance"]) for r in filtered_rows])

    expected = {
        "tau_peak_vs_r_proxy": stats.spearmanr(tau_signed, r_proxy),
        "abs_tau_peak_vs_r_proxy": stats.spearmanr(tau_abs, r_proxy),
        "tau_peak_vs_pd_significance": stats.spearmanr(tau_signed, pd_sig),
        "abs_tau_peak_vs_pd_significance": stats.spearmanr(tau_abs, pd_sig),
    }
    for pair, exp in expected.items():
        assert pop_true[pair]["rho"] == pytest.approx(exp.statistic, abs=1e-12)
        assert pop_true[pair]["p"] == pytest.approx(exp.pvalue, abs=1e-12)

    # The filter must not be a silent no-op: at least one pair's rho must
    # genuinely differ between the True (02250301 dropped) and False
    # (02250301 kept) runs by more than the numerical-noise floor.
    n_differ = sum(
        1 for pair in ALL_PAIRS
        if abs(pop_true[pair]["rho"] - pop_false[pair]["rho"]) > 1e-9
    )
    assert n_differ >= 1, (
        "exclude_noise_limited=True/False produced numerically identical rho "
        f"for all four pairs on this fixture set: true={pop_true}, false={pop_false}")


# --------------------------------------------------------------------------- #
# Edge cases: empty ok set, and exactly 2 ok rows -> all four rows emitted
# with the documented degenerate values, never an exception.
# --------------------------------------------------------------------------- #

def test_edge_case_empty_ok_set_emits_four_degenerate_rows(tmp_path):
    """Truncate a real qu_slices CSV below n_min_slices so every requested
    subobs is 'insufficient_slices' -- the ok set is empty."""
    src_path = _copy_qu_slices(str(tmp_path), OBSID_A)
    with open(src_path) as f:
        lines = f.readlines()
    with open(src_path, "w") as f:
        f.writelines(lines[:2 + 3])  # provenance + column header + 3 data rows

    out_dir = os.path.join(str(tmp_path), "out")
    res = m.qu_lag_dcf(str(tmp_path), [OBSID_A], out_dir, n_permutations=5, n_bootstrap=5)
    assert res["exit_code"] == 0

    pop = res["population_correlations"]
    assert set(pop) == ALL_PAIRS
    for pair in ALL_PAIRS:
        assert pop[pair]["n"] == 0
        assert np.isnan(pop[pair]["rho"])
        assert np.isnan(pop[pair]["p"])

    csv_pairs = _read_population_csv(res["population_csv"])
    assert set(csv_pairs) == ALL_PAIRS
    for pair in ALL_PAIRS:
        assert csv_pairs[pair]["n"] == "0"
        assert csv_pairs[pair]["rho"] == "nan"
        assert csv_pairs[pair]["p"] == "nan"


def test_edge_case_two_ok_rows_emits_four_degenerate_rows(tmp_path):
    """Exactly 2 status=='ok' subobs (both above n_min_slices): Spearman
    needs n>=3, so all four pairs must report n=2, rho/p NaN -- never raise."""
    _copy_qu_slices(str(tmp_path), OBSID_A)
    _copy_qu_slices(str(tmp_path), OBSID_B)

    out_dir = os.path.join(str(tmp_path), "out")
    res = m.qu_lag_dcf(str(tmp_path), [OBSID_A, OBSID_B], out_dir,
                        n_permutations=5, n_bootstrap=5)
    assert res["exit_code"] == 0
    assert res["per_subobs"][OBSID_A]["status"] == "ok"
    assert res["per_subobs"][OBSID_B]["status"] == "ok"

    pop = res["population_correlations"]
    assert set(pop) == ALL_PAIRS
    for pair in ALL_PAIRS:
        assert pop[pair]["n"] == 2
        assert np.isnan(pop[pair]["rho"])
        assert np.isnan(pop[pair]["p"])


# --------------------------------------------------------------------------- #
# Edge case: a NaN covariate (non-positive total exposure -> r_proxy=nan)
# drops that subobs from the r_proxy pairs (signed and abs alike) while its
# finite pd_significance keeps it in the pd_significance pairs.
# --------------------------------------------------------------------------- #

def test_nan_covariate_gives_equal_n_for_signed_and_abs_pairs(tmp_path):
    for obsid in MIXED_SIGN_OBSIDS:
        _copy_qu_slices(str(tmp_path), obsid)
    zeroed_path = os.path.join(str(tmp_path), OBSID_NULL, f"qu_slices_{OBSID_NULL}.csv")
    _zero_good_exposure(zeroed_path)

    out_dir = os.path.join(str(tmp_path), "out")
    res = m.qu_lag_dcf(str(tmp_path), MIXED_SIGN_OBSIDS, out_dir,
                        n_permutations=25, n_bootstrap=25, seed=7)
    assert res["exit_code"] == 0

    master = _read_ccf_master_raw(res["ccf_master_csv"])
    zeroed_row = next(r for r in master if r["obsid"] == OBSID_NULL)
    assert zeroed_row["status"] == "ok"
    # Not vacuous: r_proxy is genuinely non-finite here, pd_significance is not.
    assert not np.isfinite(float(zeroed_row["r_proxy"]))
    assert np.isfinite(float(zeroed_row["pd_significance"]))

    pop = res["population_correlations"]
    assert pop["tau_peak_vs_r_proxy"]["n"] == 3
    assert pop["abs_tau_peak_vs_r_proxy"]["n"] == 3
    assert pop["tau_peak_vs_pd_significance"]["n"] == 4
    assert pop["abs_tau_peak_vs_pd_significance"]["n"] == 4

    ok_rows = [r for r in master if r["status"] == "ok"]
    tau_signed = np.array([float(r["tau_peak_s"]) for r in ok_rows])
    tau_abs = np.abs(tau_signed)
    r_proxy = np.array([float(r["r_proxy"]) for r in ok_rows])
    pd_sig = np.array([float(r["pd_significance"]) for r in ok_rows])

    for series, pair, covariate in (
        (tau_signed, "tau_peak_vs_r_proxy", r_proxy),
        (tau_abs, "abs_tau_peak_vs_r_proxy", r_proxy),
        (tau_signed, "tau_peak_vs_pd_significance", pd_sig),
        (tau_abs, "abs_tau_peak_vs_pd_significance", pd_sig),
    ):
        mask = np.isfinite(series) & np.isfinite(covariate)
        expected = stats.spearmanr(series[mask], covariate[mask])
        assert pop[pair]["n"] == int(mask.sum())
        assert pop[pair]["rho"] == pytest.approx(expected.statistic, abs=1e-12)
        assert pop[pair]["p"] == pytest.approx(expected.pvalue, abs=1e-12)


# --------------------------------------------------------------------------- #
# Reproducibility: identical seed over identical inputs -> all four
# population rho/p equal within tolerance, and ccf_master.csv +
# population_correlations.csv both byte-identical.
# --------------------------------------------------------------------------- #

def test_repro_qu_lag_dcf_population(tmp_path_factory):
    results = []
    ccf_bytes = []
    population_bytes = []
    for i in range(2):
        tmp_dir = str(tmp_path_factory.mktemp(f"qu_lag_pop_repro_{i}"))
        for obsid in MIXED_SIGN_OBSIDS:
            _copy_qu_slices(tmp_dir, obsid)
        res = m.qu_lag_dcf(tmp_dir, MIXED_SIGN_OBSIDS, os.path.join(tmp_dir, "out"),
                            n_permutations=50, n_bootstrap=50, seed=42)
        results.append(res["population_correlations"])
        with open(res["ccf_master_csv"], "rb") as f:
            ccf_bytes.append(f.read())
        with open(res["population_csv"], "rb") as f:
            population_bytes.append(f.read())

    a, b = results
    assert set(a) == ALL_PAIRS
    for pair in ALL_PAIRS:
        assert a[pair]["n"] == b[pair]["n"]
        if np.isnan(a[pair]["rho"]):
            assert np.isnan(b[pair]["rho"])
        else:
            assert a[pair]["rho"] == pytest.approx(b[pair]["rho"], abs=1e-12)
            assert a[pair]["p"] == pytest.approx(b[pair]["p"], abs=1e-12)

    assert ccf_bytes[0] == ccf_bytes[1]
    assert population_bytes[0] == population_bytes[1]


# --------------------------------------------------------------------------- #
# Output contract: population_correlations.csv row order (signed rows first,
# unchanged, then the two abs rows) and column shape (POPULATION_FIELDS
# unchanged).
# --------------------------------------------------------------------------- #

def test_population_csv_row_order_and_field_shape(tmp_path):
    for obsid in MIXED_SIGN_OBSIDS:
        _copy_qu_slices(str(tmp_path), obsid)
    out_dir = os.path.join(str(tmp_path), "out")
    res = m.qu_lag_dcf(str(tmp_path), MIXED_SIGN_OBSIDS, out_dir,
                        n_permutations=5, n_bootstrap=5, seed=0)

    with open(res["population_csv"], newline="") as f:
        reader = csv.DictReader(f)
        assert reader.fieldnames == m.POPULATION_FIELDS == ["pair", "n", "rho", "p"]
        rows = list(reader)
    assert [r["pair"] for r in rows] == [
        "tau_peak_vs_r_proxy", "tau_peak_vs_pd_significance",
        "abs_tau_peak_vs_r_proxy", "abs_tau_peak_vs_pd_significance",
    ]


# --------------------------------------------------------------------------- #
# Manifest wiring: output/returns documentation mentions the four pairs.
# --------------------------------------------------------------------------- #

def test_tool_manifest_documents_four_population_pairs():
    from xspec_utils.skills import TOOL_MANIFEST
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.ixpe_qu_lag_dcf"]
    desc = entry["args_schema"]["description"]
    for pair in ALL_PAIRS:
        assert pair in desc
    assert "sense" in desc.lower()
    assert "magnitude" in desc.lower()
