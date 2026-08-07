"""Tests for xspec_utils' ixpe_qu_lag_dcf (lib/qu_lag_dcf.py).

Uses the three qu_slices_<subobs_key>.csv products already on disk under
/data/fzhou/XSPEC/Cyg_X-1/IXPE (the data catalog confirms all three local),
copied into an isolated tmp obs_db_path per test so the tool's own writes
(qu_dcf_<subobs_key>.csv, catalog.json) never touch the shared analyst-env
mission directory:

- 01002901: the real, already-local Cyg X-1 IXPE ObsID (read-only Layer-A
  differential/reproducibility grounding).
- synth_qu_lag_evensamp01: a marked-synthetic ObsID (an isolated Gaussian
  q(t) bump with an identical bump in u(t) after a 300s injected lag,
  evenly sampled at 100s cadence) -- injection-recovery + a genuine
  differential capability check.
- synth_qu_lag_null01: a marked-synthetic ObsID (pure per-slice Gaussian
  noise around a constant q0/u0, with a stated per-slice error far larger
  than the actual sample scatter) -- the guard-null check, and the
  variance-floor uncorrected-denominator fallback.

Both synthetic ObsIDs' qu_slices CSVs carry their own ground truth (injected
lag, seed, ...) in the `qu_slice_pcubes provenance` JSON header -- read from
there, never hardcoded, so the tests stay grounded in the data's own marked
provenance rather than assumptions about it.
"""

from __future__ import annotations

import csv
import json
import math
import os
import shutil

import numpy as np
import pytest

from obsdb_catalog import neutral
from xspec_utils.lib import qu_lag_dcf as m
from xspec_utils.skills import TOOL_MANIFEST

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
OBSID_REAL = "01002901"
OBSID_EVENSAMP = "synth_qu_lag_evensamp01"
OBSID_NULL = "synth_qu_lag_null01"
ALL_OBSIDS = [OBSID_REAL, OBSID_EVENSAMP, OBSID_NULL]


def _independent_denominator(q, q_err, u, u_err, variance_floor=1e-12):
    """Faithful plain-numpy reimplementation of the documented D:
    D = (var(q) - mean(q_err**2)) * (var(u) - mean(u_err**2)) with
    np.var (population, ddof=0); falls back to var(q)*var(u) when either
    corrected variance <= variance_floor; floored at variance_floor if
    still non-positive. Does NOT call m.denominator / m._corrected_var."""
    vq = float(np.var(q) - np.mean(np.asarray(q_err) ** 2))
    vu = float(np.var(u) - np.mean(np.asarray(u_err) ** 2))
    if vq <= variance_floor or vu <= variance_floor:
        D = float(np.var(q) * np.var(u))
    else:
        D = float(vq * vu)
    return D if D > 0 else variance_floor


def _copy_qu_slices(dest_dir: str, obsid: str) -> str:
    """Copy `obsid`'s real, already-local qu_slices_<obsid>.csv (only that
    file -- the tool never reads an ObsID's Level-2 event data, so copying
    the multi-GB event_l2/ tree alongside it would only slow the suite down)
    into `<dest_dir>/<obsid>/`. Returns the copied CSV path."""
    fname = "qu_slices_" + obsid + ".csv"
    src_csv = os.path.join(REAL_IXPE_DIR, obsid, fname)
    if not os.path.exists(src_csv):
        pytest.skip(f"fixture qu_slices CSV not local at {src_csv}")
    dst_dir = os.path.join(dest_dir, obsid)
    os.makedirs(dst_dir, exist_ok=True)
    dst_csv = os.path.join(dst_dir, fname)
    shutil.copyfile(src_csv, dst_csv)
    return dst_csv


@pytest.fixture()
def obs_db(tmp_path):
    """A scratch obs_db_path carrying all three fixture ObsIDs' qu_slices CSVs."""
    for obsid in ALL_OBSIDS:
        _copy_qu_slices(str(tmp_path), obsid)
    return str(tmp_path)


# --------------------------------------------------------------------------- #
# Pure-logic unit tests: denominator fallback, symmetric zero-centered bins
# --------------------------------------------------------------------------- #

def test_denominator_uses_corrected_variance_when_above_floor():
    rng = np.random.default_rng(1)
    q = rng.normal(0.0, 1.0, 200)
    u = rng.normal(0.0, 1.0, 200)
    q_err = np.full(200, 0.01)
    u_err = np.full(200, 0.01)
    D, uncorrected = m.denominator(q, q_err, u, u_err, variance_floor=1e-12)
    assert uncorrected is False
    expected = (np.var(q) - np.mean(q_err ** 2)) * (np.var(u) - np.mean(u_err ** 2))
    assert D == pytest.approx(expected)


def test_denominator_falls_back_when_corrected_variance_below_floor():
    # Constant series (zero true variance) with a stated per-point error that
    # would drive the corrected variance negative -- must trigger the fallback.
    q = np.full(50, 0.02)
    u = np.full(50, -0.01)
    q_err = np.full(50, 0.05)
    u_err = np.full(50, 0.05)
    D, uncorrected = m.denominator(q, q_err, u, u_err, variance_floor=1e-12)
    assert uncorrected is True
    assert D == pytest.approx(np.var(q) * np.var(u)) or D == pytest.approx(1e-12)


def test_dcf_bins_always_center_one_bin_on_zero_lag():
    edges, centers = m.dcf_bins(width=10.0, half_range=55.0)
    assert 0.0 in centers
    assert edges.size == centers.size + 1
    assert np.allclose(np.diff(centers), 10.0)


def test_bin_udcf_known_input():
    dt = np.array([-5.0, 0.0, 0.0, 5.0])
    udcf = np.array([1.0, 2.0, 4.0, 8.0])
    edges = np.array([-7.5, -2.5, 2.5, 7.5])
    dcf, n_pairs = m.bin_udcf(dt, udcf, edges)
    assert n_pairs.tolist() == [1, 2, 1]
    assert dcf[0] == pytest.approx(1.0)
    assert dcf[1] == pytest.approx(3.0)  # mean(2, 4)
    assert dcf[2] == pytest.approx(8.0)


def test_r_proxy_and_pd_significance_known_input():
    rows = [
        {'counts': 100.0, 'good_exposure_s': 50.0, 'q': 0.02, 'q_err': 0.01, 'u': 0.01, 'u_err': 0.01},
        {'counts': 200.0, 'good_exposure_s': 50.0, 'q': 0.03, 'q_err': 0.01, 'u': 0.02, 'u_err': 0.01},
    ]
    r_proxy = m.r_proxy_from_rows(rows, src_rad=1.0)
    assert r_proxy == pytest.approx(300.0 / (100.0 * math.pi))

    pd_sig = m.pd_significance_from_rows(rows)
    q_w, sigma_q = m._weighted_mean_and_err([0.02, 0.03], [0.01, 0.01])
    u_w, sigma_u = m._weighted_mean_and_err([0.01, 0.02], [0.01, 0.01])
    pd = math.sqrt(q_w ** 2 + u_w ** 2)
    sigma_pd = math.sqrt((q_w * sigma_q) ** 2 + (u_w * sigma_u) ** 2) / pd
    assert pd_sig == pytest.approx(pd / sigma_pd)


# --------------------------------------------------------------------------- #
# Manifest wiring (AI-facing half)
# --------------------------------------------------------------------------- #

def test_tool_manifest_exposes_class_parameters():
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.ixpe_qu_lag_dcf"]
    assert entry["module"] == "xspec_utils.lib.qu_lag_dcf"
    assert entry["fn"] == "qu_lag_dcf"
    props = entry["args_schema"]["properties"]
    for name in ("obs_db_path", "obsids", "output_dir", "sources", "lag_bin_s",
                 "n_min_slices", "n_permutations", "n_bootstrap", "variance_floor",
                 "coherent_p_threshold", "max_lag_frac", "exclude_noise_limited", "seed"):
        assert name in props, f"{name} missing from args_schema.properties"
    assert entry["args_schema"]["required"] == ["obs_db_path", "obsids", "output_dir"]


# --------------------------------------------------------------------------- #
# Edge case: fewer than n_min_slices -> insufficient_slices, excluded, no raise
# --------------------------------------------------------------------------- #

def test_insufficient_slices_excluded_not_raised(tmp_path):
    """Truncate a real qu_slices CSV to 3 rows (< default n_min_slices=8):
    must be recorded as 'insufficient_slices' in ccf_master.csv, excluded from
    the population correlations, no qu_dcf CSV written, no catalog entry --
    never raised."""
    src_path = _copy_qu_slices(str(tmp_path), OBSID_REAL)
    with open(src_path) as f:
        lines = f.readlines()
    # header + column header + 3 data rows
    with open(src_path, "w") as f:
        f.writelines(lines[:2 + 3])

    out_dir = os.path.join(str(tmp_path), "out")
    res = m.qu_lag_dcf(str(tmp_path), [OBSID_REAL], out_dir,
                        n_permutations=5, n_bootstrap=5)
    assert res["exit_code"] == 0
    info = res["per_subobs"][OBSID_REAL]
    assert info["status"] == "insufficient_slices"
    assert info["n_slices"] == 3
    assert info["csv_path"] is None
    assert not any(f.startswith("qu_dcf_") for f in os.listdir(os.path.join(str(tmp_path), OBSID_REAL)))
    assert neutral.read(str(tmp_path))["entries"] == {}

    with open(os.path.join(out_dir, "ccf_master.csv"), newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1
    assert rows[0]["status"] == "insufficient_slices"
    assert rows[0]["tau_peak_s"] == "nan"


# --------------------------------------------------------------------------- #
# Trust gap: a missing/unreadable upstream qu_slices file must never vanish
# --------------------------------------------------------------------------- #

def test_missing_qu_slices_recorded_not_dropped(tmp_path):
    """An obsid with no qu_slices file on disk at all (a failed or
    not-yet-run ixpe_qu_slice_pcubes step) must appear in ccf_master.csv/
    per_subobs as status='no_qu_slices' -- never silently omitted. Combined
    with a real ObsID in the same call to confirm both are surfaced and only
    the real, 'ok' one feeds the population correlations."""
    _copy_qu_slices(str(tmp_path), OBSID_REAL)
    missing_obsid = "09999999"  # no qu_slices_*.csv under this obsid's directory

    out_dir = os.path.join(str(tmp_path), "out")
    res = m.qu_lag_dcf(str(tmp_path), [OBSID_REAL, missing_obsid], out_dir,
                        n_permutations=5, n_bootstrap=5)
    assert res["exit_code"] == 0

    missing_info = res["per_subobs"][missing_obsid]
    assert missing_info["status"] == "no_qu_slices"
    assert missing_info["n_slices"] == 0
    assert missing_info["csv_path"] is None
    real_info = res["per_subobs"][OBSID_REAL]
    assert real_info["status"] == "ok"

    with open(os.path.join(out_dir, "ccf_master.csv"), newline="") as f:
        rows = {r["subobs_key"]: r for r in csv.DictReader(f)}
    assert len(rows) == 2
    assert rows[missing_obsid]["status"] == "no_qu_slices"
    assert rows[missing_obsid]["tau_peak_s"] == "nan"
    assert rows[OBSID_REAL]["status"] == "ok"

    with open(os.path.join(out_dir, "population_correlations.csv"), newline="") as f:
        pop_rows = list(csv.DictReader(f))
    for row in pop_rows:
        assert int(row["n"]) == 1  # only the one 'ok' row (01002901) counted

    entries = neutral.read(str(tmp_path))["entries"]
    assert missing_obsid not in entries


def test_unreadable_qu_slices_recorded_not_dropped(tmp_path):
    """A qu_slices CSV whose first line is not the expected provenance
    header must be recorded status='unreadable_qu_slices' -- never silently
    skipped -- with no qu_dcf CSV/catalog entry written for it."""
    csv_path = _copy_qu_slices(str(tmp_path), OBSID_REAL)
    with open(csv_path) as f:
        lines = f.readlines()
    lines[0] = "garbage,header\n"
    with open(csv_path, "w") as f:
        f.writelines(lines)

    out_dir = os.path.join(str(tmp_path), "out")
    res = m.qu_lag_dcf(str(tmp_path), [OBSID_REAL], out_dir, n_permutations=5, n_bootstrap=5)
    assert res["exit_code"] == 0

    info = res["per_subobs"][OBSID_REAL]
    assert info["status"] == "unreadable_qu_slices"
    assert info["n_slices"] == 0
    assert info["csv_path"] is None

    with open(os.path.join(out_dir, "ccf_master.csv"), newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1
    assert rows[0]["status"] == "unreadable_qu_slices"

    assert not any(f.startswith("qu_dcf_") for f in os.listdir(os.path.join(str(tmp_path), OBSID_REAL)))
    assert neutral.read(str(tmp_path))["entries"] == {}


# --------------------------------------------------------------------------- #
# Reproducibility: identical seed -> identical tau_peak/dcf_peak/perm_p
# --------------------------------------------------------------------------- #

def test_repro_qu_lag_dcf_stable_across_calls(tmp_path_factory):
    results = []
    for i in range(2):
        tmp_dir = str(tmp_path_factory.mktemp(f"qu_lag_repro_{i}"))
        _copy_qu_slices(tmp_dir, OBSID_REAL)
        res = m.qu_lag_dcf(tmp_dir, [OBSID_REAL], os.path.join(tmp_dir, "out"),
                            n_permutations=100, n_bootstrap=100, seed=42)
        results.append(res["per_subobs"][OBSID_REAL])
    a, b = results
    for key in ("tau_peak_s", "dcf_peak", "perm_p"):
        assert a[key] == pytest.approx(b[key], abs=1e-9)


# --------------------------------------------------------------------------- #
# Capability: injection-recovery on an evenly-sampled, injected-lag fixture
# --------------------------------------------------------------------------- #

def test_capability_injection_recovery_evensamp(tmp_path):
    csv_path = _copy_qu_slices(str(tmp_path), OBSID_EVENSAMP)
    with open(csv_path) as f:
        provenance = json.loads(f.readline()[len("# qu_slice_pcubes provenance: "):])
    injected_lag_s = provenance["injected_lag_s"]

    res = m.qu_lag_dcf(str(tmp_path), [OBSID_EVENSAMP], os.path.join(str(tmp_path), "out"),
                        n_permutations=200, n_bootstrap=200, seed=0)
    info = res["per_subobs"][OBSID_EVENSAMP]
    assert info["status"] == "ok"
    # Recovered lag matches the injected lag within one lag bin.
    with open(info["csv_path"]) as f:
        dcf_provenance = json.loads(f.readline()[len("# qu_lag_dcf provenance: "):])
    lag_bin_s_used = dcf_provenance["lag_bin_s_used"]
    assert abs(info["tau_peak_s"] - injected_lag_s) <= lag_bin_s_used
    # Fails if the tool returned well-shaped zeros/metadata without actually
    # locating the correlation: the permutation null must call it significant.
    assert info["perm_p"] < 0.01
    assert info["coherent_rotation"] is True


# --------------------------------------------------------------------------- #
# Guard: uncorrelated series must not be flagged coherent
# --------------------------------------------------------------------------- #

def test_guard_null_no_significant_lag(tmp_path):
    _copy_qu_slices(str(tmp_path), OBSID_NULL)
    res = m.qu_lag_dcf(str(tmp_path), [OBSID_NULL], os.path.join(str(tmp_path), "out"),
                        n_permutations=200, n_bootstrap=200, seed=0)
    info = res["per_subobs"][OBSID_NULL]
    assert info["status"] == "ok"
    assert info["perm_p"] > 0.1
    assert info["coherent_rotation"] is False


# --------------------------------------------------------------------------- #
# Differential: dcf_zero vs an independent zero-lag noise-corrected Pearson
# correlation, computed via a plain numpy path (not the tool's binning helper)
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("obsid", ALL_OBSIDS)
def test_differential_dcf_zero_matches_independent_pearson(tmp_path, obsid):
    csv_path = _copy_qu_slices(str(tmp_path), obsid)
    _, rows = m._read_qu_slices_csv(csv_path)
    q = np.array([r["q"] for r in rows])
    q_err = np.array([r["q_err"] for r in rows])
    u = np.array([r["u"] for r in rows])
    u_err = np.array([r["u_err"] for r in rows])
    a, b = q - np.mean(q), u - np.mean(u)
    D = _independent_denominator(q, q_err, u, u_err, variance_floor=1e-12)
    independent_zero_lag = np.mean(a * b) / math.sqrt(D)

    res = m.qu_lag_dcf(str(tmp_path), [obsid], os.path.join(str(tmp_path), "out"),
                        n_permutations=5, n_bootstrap=5)
    with open(os.path.join(str(tmp_path), "out", "ccf_master.csv"), newline="") as f:
        row = next(csv.DictReader(f))
    assert float(row["dcf_zero"]) == pytest.approx(independent_zero_lag, abs=0.01)


# --------------------------------------------------------------------------- #
# Differential: r_proxy vs an independent recompute from the same CSV rows
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("obsid", ALL_OBSIDS)
def test_differential_r_proxy_matches_independent_computation(tmp_path, obsid):
    csv_path = _copy_qu_slices(str(tmp_path), obsid)
    provenance, rows = m._read_qu_slices_csv(csv_path)
    total_counts = sum(r["counts"] for r in rows)
    total_exposure = sum(r["good_exposure_s"] for r in rows)
    src_rad = provenance["src_rad"]
    independent = total_counts / (total_exposure * math.pi * src_rad ** 2)

    res = m.qu_lag_dcf(str(tmp_path), [obsid], os.path.join(str(tmp_path), "out"),
                        n_permutations=5, n_bootstrap=5)
    with open(os.path.join(str(tmp_path), "out", "ccf_master.csv"), newline="") as f:
        row = next(csv.DictReader(f))
    assert float(row["r_proxy"]) == pytest.approx(independent, rel=1e-9)


# --------------------------------------------------------------------------- #
# Output contract: CSVs on disk, catalog recording, population correlations
# --------------------------------------------------------------------------- #

def test_output_contract_csvs_and_catalog_recording(obs_db):
    out_dir = os.path.join(obs_db, "out")
    res = m.qu_lag_dcf(obs_db, ALL_OBSIDS, out_dir, n_permutations=20, n_bootstrap=20, seed=0)
    assert res["exit_code"] == 0

    assert os.path.exists(res["ccf_master_csv"])
    assert os.path.exists(res["population_csv"])
    with open(res["ccf_master_csv"], newline="") as f:
        master_rows = list(csv.DictReader(f))
    assert len(master_rows) == 3
    assert {r["status"] for r in master_rows} == {"ok"}

    for obsid in ALL_OBSIDS:
        info = res["per_subobs"][obsid]
        assert os.path.exists(info["csv_path"])
        entry = neutral.read(obs_db)["entries"][obsid]
        assert entry["max_level"] == 3
        products = [p for p in entry["products"] if p["kind"] == "qu_dcf"]
        assert len(products) == 1
        assert products[0]["produced_by"] == "xspec_utils.lib.qu_lag_dcf.qu_lag_dcf"
        assert products[0]["params"]["subobs_key"] == obsid

    with open(res["population_csv"], newline="") as f:
        pop_rows = {r["pair"]: r for r in csv.DictReader(f)}
    assert set(pop_rows) == {
        "tau_peak_vs_r_proxy", "tau_peak_vs_pd_significance",
        "abs_tau_peak_vs_r_proxy", "abs_tau_peak_vs_pd_significance",
    }
    for row in pop_rows.values():
        assert int(row["n"]) == 3


def test_qu_dcf_csv_has_provenance_header_and_expected_columns(obs_db):
    out_dir = os.path.join(obs_db, "out")
    res = m.qu_lag_dcf(obs_db, [OBSID_EVENSAMP], out_dir, n_permutations=10, n_bootstrap=10, seed=0)
    csv_path = res["per_subobs"][OBSID_EVENSAMP]["csv_path"]
    with open(csv_path, newline="") as f:
        header_line = f.readline()
        assert header_line.startswith("# qu_lag_dcf provenance: ")
        reader = csv.DictReader(f)
        assert reader.fieldnames == m._QU_DCF_FIELDS
        rows = list(reader)
    assert len(rows) > 0
    assert any(float(r["lag_s"]) == 0.0 for r in rows)
