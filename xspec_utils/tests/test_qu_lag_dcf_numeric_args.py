"""Tests for the ixpe_qu_lag_dcf sandbox numeric/boolean argument marshalling
fix (fix-ixpe-qu-lag-dcf-numeric-args) -- second and final link of the same
chain as ``test_ixpe_qu_slice_pcubes_numeric_args.py``.

The sandbox harness invokes a tool as a bare ``result = fn(**args)`` with no
schema coercion, so a numeric/boolean TOOL_MANIFEST argument declared
integer/number/boolean can arrive at ``qu_lag_dcf``/``qu_lag_dcf_plot`` as a
Python ``str``. Before this fix an explicit ``n_min_slices`` (even at its own
documented default, 8) crashed with ``TypeError: '<' not supported between
instances of 'int' and 'str'`` at line 532's ``if n_slices < n_min_slices``,
so a caller could only ever run the tool on its built-in defaults -- exactly
by luck, since those defaults happen to match the qu-lag census prescription.

These tests pin the fix on the real, already-local qu_slices series for
Cyg X-1 IXPE ObsID 01002901 (8 slices) and the two marked-synthetic
validation ObsIDs ``test_qu_lag_dcf.py`` already uses, copied into an
isolated tmp obs_db_path per fixture so the tool's own writes never touch
the shared analyst-env mission directory. ``n_permutations``/``n_bootstrap``
are kept small (20-200), not the tool's production default (1000), except
where reproducing the already-committed ``qu_dcf_01002901.csv`` numbers
requires the real defaults (that reproduction is itself fast: 01002901 has
only 8 slices, so the O(N^2)-pair DCF binning is tiny even at 1000
permutations/bootstrap draws).
"""

from __future__ import annotations

import ast
import csv
import inspect
import json
import math
import os
import shutil

import pytest

from obsdb_catalog import neutral
from xspec_utils.lib import qu_lag_dcf as m
# Bare import, not the dotted "xspec_utils.lib.toolargs" -- xspec_utils'
# package __init__ already put lib/ on sys.path, and qu_lag_dcf.py itself
# imports toolargs the same bare way, so this is the identical module object.
import toolargs

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
OBSID_REAL = "01002901"
OBSID_EVENSAMP = "synth_qu_lag_evensamp01"
_QU_DCF_HEADER_PREFIX = "# qu_lag_dcf provenance: "
_QU_SLICE_HEADER_PREFIX = "# qu_slice_pcubes provenance: "


def _copy_qu_slices(dest_dir: str, obsid: str) -> str:
    """Copy `obsid`'s real, already-local `qu_slices_<obsid>.csv` into
    `<dest_dir>/<obsid>/`. Returns the copied CSV path."""
    fname = "qu_slices_" + obsid + ".csv"
    src_csv = os.path.join(REAL_IXPE_DIR, obsid, fname)
    if not os.path.exists(src_csv):
        pytest.skip(f"fixture qu_slices CSV not local at {src_csv}")
    dst_dir = os.path.join(dest_dir, obsid)
    os.makedirs(dst_dir, exist_ok=True)
    dst_csv = os.path.join(dst_dir, fname)
    shutil.copyfile(src_csv, dst_csv)
    return dst_csv


def _read_committed_qu_dcf_rows(obsid: str) -> list[dict]:
    """Read the already-committed, pre-fix `qu_dcf_<obsid>.csv` shipped
    alongside the real qu_slices fixture -- the science-path reference this
    fix must not move (its own provenance header carries the production
    defaults: n_permutations=1000, n_bootstrap=1000, seed=0,
    variance_floor=1e-12, coherent_p_threshold=0.01, max_lag_frac=0.5)."""
    path = os.path.join(REAL_IXPE_DIR, obsid, f"qu_dcf_{obsid}.csv")
    if not os.path.exists(path):
        pytest.skip(f"committed qu_dcf CSV not local at {path}")
    with open(path, newline="") as f:
        header = f.readline()
        assert header.startswith(_QU_DCF_HEADER_PREFIX)
        rows = list(csv.DictReader(f))
    return rows


def _rewrite_qu_slice_n_min(csv_path: str, new_n_min: int) -> None:
    """Bump only the provenance-header `n_min` field of an already-copied
    `qu_slices_<subobs_key>.csv` fixture, leaving every data row byte-
    identical. `noise_limited` is, by definition,
    `n_rows < provenance['n_min']` (the *upstream* `ixpe_qu_slice_pcubes`
    call's own threshold metadata) -- rewriting that one metadata field on a
    scratch copy is the only way to exercise the noise_limited branch from
    real data, in the same spirit as the corrupt-header/truncation tests
    already in test_qu_lag_dcf.py. No science value is invented; only the
    upstream tool's own recorded threshold is edited."""
    with open(csv_path, newline="") as f:
        lines = f.readlines()
    header = json.loads(lines[0][len(_QU_SLICE_HEADER_PREFIX):])
    header["n_min"] = new_n_min
    lines[0] = _QU_SLICE_HEADER_PREFIX + json.dumps(header, sort_keys=True) + "\n"
    with open(csv_path, "w", newline="") as f:
        f.writelines(lines)


def _nan_aware_equal(a, b) -> bool:
    """Recursive dict/float equality that treats NaN==NaN as equal --
    plain `==` on two independently-computed `population_correlations`
    dicts would otherwise fail on a NaN Spearman rho/p even when both sides
    are the identical NaN."""
    if isinstance(a, dict) and isinstance(b, dict):
        return set(a) == set(b) and all(_nan_aware_equal(a[k], b[k]) for k in a)
    if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
        return True
    return a == b


# --------------------------------------------------------------------------- #
# Differential: coercion is obtained from lib/toolargs.py, never a private
# re-implementation -- paired with a science-path differential against the
# already-committed qu_dcf_01002901.csv (computed before this fix, at the
# real production defaults).
# --------------------------------------------------------------------------- #

def test_toolargs_reused_not_reimplemented():
    source = inspect.getsource(m)
    tree = ast.parse(source)

    imports_from_toolargs = any(
        isinstance(node, ast.ImportFrom)
        and node.module == "toolargs"
        and any(alias.name == "coerce_tool_args" for alias in node.names)
        for node in ast.walk(tree)
    )
    assert imports_from_toolargs, "qu_lag_dcf.py must import coerce_tool_args from toolargs"
    assert m.coerce_tool_args is toolargs.coerce_tool_args

    called_names = {
        node.func.id for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    assert "coerce_tool_args" in called_names

    # No second, private coercion implementation living in this module --
    # that would be exactly the shared-bug-on-both-sides failure mode the
    # brief calls out.
    private_coercers = [
        node.name for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef) and "coerce" in node.name.lower()
    ]
    assert private_coercers == []


# --------------------------------------------------------------------------- #
# Contract: the brief's step 1 requires the public signature and every
# default preserved EXACTLY (value and native type) -- a future edit that
# swaps a default for a string, or reorders/renames a parameter, is caught
# here before it ever reaches a real call.
# --------------------------------------------------------------------------- #

def test_public_signature_and_defaults_unchanged():
    sig = inspect.signature(m.qu_lag_dcf)
    assert list(sig.parameters) == [
        "obs_db_path", "obsids", "output_dir", "sources", "lag_bin_s",
        "n_min_slices", "n_permutations", "n_bootstrap", "variance_floor",
        "coherent_p_threshold", "max_lag_frac", "exclude_noise_limited", "seed",
    ]
    defaults = {
        name: p.default for name, p in sig.parameters.items()
        if p.default is not inspect.Parameter.empty
    }
    assert defaults["sources"] is None
    assert defaults["lag_bin_s"] is None
    assert defaults["n_min_slices"] == 8 and type(defaults["n_min_slices"]) is int
    assert defaults["n_permutations"] == 1000 and type(defaults["n_permutations"]) is int
    assert defaults["n_bootstrap"] == 1000 and type(defaults["n_bootstrap"]) is int
    assert defaults["variance_floor"] == 1e-12 and type(defaults["variance_floor"]) is float
    assert defaults["coherent_p_threshold"] == 0.01 and type(defaults["coherent_p_threshold"]) is float
    assert defaults["max_lag_frac"] == 0.5 and type(defaults["max_lag_frac"]) is float
    assert defaults["exclude_noise_limited"] is True and type(defaults["exclude_noise_limited"]) is bool
    assert defaults["seed"] == 0 and type(defaults["seed"]) is int

    sig_plot = inspect.signature(m.qu_lag_dcf_plot)
    assert list(sig_plot.parameters) == [
        "obs_db_path", "ccf_master_csv", "output_dir", "save_path_prefix", "max_panels",
    ]
    defaults_plot = {
        name: p.default for name, p in sig_plot.parameters.items()
        if p.default is not inspect.Parameter.empty
    }
    assert defaults_plot["save_path_prefix"] is None
    assert defaults_plot["max_panels"] == 12 and type(defaults_plot["max_panels"]) is int


def test_science_path_reproduces_committed_qu_dcf_csv(tmp_path):
    """Running qu_lag_dcf with every numeric kwarg at its documented default
    (the coerced path, native-typed call) must reproduce the lag_s/dcf
    columns already committed in qu_dcf_01002901.csv -- the same computation
    this fix must never perturb for an already-native call."""
    committed_rows = _read_committed_qu_dcf_rows(OBSID_REAL)
    assert len(committed_rows) > 1  # else the loop below would pass vacuously

    _copy_qu_slices(str(tmp_path), OBSID_REAL)
    out_dir = os.path.join(str(tmp_path), "out")
    res = m.qu_lag_dcf(str(tmp_path), [OBSID_REAL], out_dir)
    assert res["exit_code"] == 0
    info = res["per_subobs"][OBSID_REAL]
    assert info["status"] == "ok"

    with open(info["csv_path"], newline="") as f:
        f.readline()  # provenance header
        fresh_rows = list(csv.DictReader(f))

    assert len(fresh_rows) == len(committed_rows)
    committed_by_lag = {
        round(float(r["lag_s"]), 3): (float(r["dcf"]), int(r["n_pairs"])) for r in committed_rows
    }
    for row in fresh_rows:
        lag = round(float(row["lag_s"]), 3)
        assert lag in committed_by_lag
        committed_dcf, committed_n_pairs = committed_by_lag[lag]
        assert float(row["dcf"]) == pytest.approx(committed_dcf, abs=1e-9)
        # n_pairs is RNG-independent (fixed by the slice timestamps alone,
        # like lag_s/dcf), so it must match exactly, not just approximately.
        assert int(row["n_pairs"]) == committed_n_pairs


# --------------------------------------------------------------------------- #
# Regression: the reported defect -- string-typed numerics, exactly as the
# sandbox harness delivers them, must succeed and agree with native-typed
# equivalents.
# --------------------------------------------------------------------------- #

def test_numeric_kwarg_string_vs_native_dcf(tmp_path_factory):
    string_kwargs = {
        "n_min_slices": "8", "n_permutations": "200", "n_bootstrap": "200",
        "seed": "0", "variance_floor": "1e-12",
    }
    native_kwargs = {
        "n_min_slices": 8, "n_permutations": 200, "n_bootstrap": 200,
        "seed": 0, "variance_floor": 1e-12,
    }

    tmp_string = str(tmp_path_factory.mktemp("qu_lag_numeric_string"))
    _copy_qu_slices(tmp_string, OBSID_REAL)
    res_string = m.qu_lag_dcf(tmp_string, [OBSID_REAL], os.path.join(tmp_string, "out"),
                               **string_kwargs)

    tmp_native = str(tmp_path_factory.mktemp("qu_lag_numeric_native"))
    _copy_qu_slices(tmp_native, OBSID_REAL)
    res_native = m.qu_lag_dcf(tmp_native, [OBSID_REAL], os.path.join(tmp_native, "out"),
                               **native_kwargs)

    assert res_string["exit_code"] == 0, res_string
    assert res_native["exit_code"] == 0, res_native

    info_s = res_string["per_subobs"][OBSID_REAL]
    info_n = res_native["per_subobs"][OBSID_REAL]
    assert info_s["status"] == info_n["status"] == "ok"
    for field in ("tau_peak_s", "dcf_peak", "perm_p"):
        assert info_s[field] == info_n[field]


def test_all_declared_args_string_vs_native_identical(tmp_path_factory):
    """All nine coerced arguments, string- vs native-typed, compared on the
    FULL products (every ccf_master.csv column, every qu_dcf data row --
    including the RNG-dependent perm_lo99/perm_hi99 columns) rather than
    three scalars, so a coercion bug that only shows up in the permutation/
    bootstrap streams (not in tau_peak_s/dcf_peak/perm_p) cannot slip past."""
    native_kwargs = {
        "lag_bin_s": None, "n_min_slices": 8, "n_permutations": 200, "n_bootstrap": 200,
        "variance_floor": 1e-12, "coherent_p_threshold": 0.01, "max_lag_frac": 0.5,
        "exclude_noise_limited": True, "seed": 0,
    }
    string_kwargs = {
        "lag_bin_s": "none", "n_min_slices": "8", "n_permutations": "200", "n_bootstrap": "200",
        "variance_floor": "1e-12", "coherent_p_threshold": "0.01", "max_lag_frac": "0.5",
        "exclude_noise_limited": "true", "seed": "0",
    }

    tmp_native = str(tmp_path_factory.mktemp("qu_lag_all_args_native"))
    _copy_qu_slices(tmp_native, OBSID_REAL)
    res_native = m.qu_lag_dcf(tmp_native, [OBSID_REAL], os.path.join(tmp_native, "out"), **native_kwargs)

    tmp_string = str(tmp_path_factory.mktemp("qu_lag_all_args_string"))
    _copy_qu_slices(tmp_string, OBSID_REAL)
    res_string = m.qu_lag_dcf(tmp_string, [OBSID_REAL], os.path.join(tmp_string, "out"), **string_kwargs)

    assert res_native["exit_code"] == 0, res_native
    assert res_string["exit_code"] == 0, res_string
    info_n = res_native["per_subobs"][OBSID_REAL]
    info_s = res_string["per_subobs"][OBSID_REAL]
    assert info_n["status"] == info_s["status"] == "ok"

    with open(res_native["ccf_master_csv"], newline="") as f:
        master_native = list(csv.DictReader(f))
    with open(res_string["ccf_master_csv"], newline="") as f:
        master_string = list(csv.DictReader(f))
    assert master_native == master_string  # every CCF_MASTER_FIELDS column, verbatim

    def _dcf_rows(info):
        with open(info["csv_path"], newline="") as f:
            f.readline()  # provenance header
            return list(csv.DictReader(f))

    dcf_rows_native = _dcf_rows(info_n)
    dcf_rows_string = _dcf_rows(info_s)
    assert len(dcf_rows_native) > 1  # else the equality below could pass vacuously
    # Verbatim row equality, including perm_lo99/perm_hi99 -- this pins the
    # coerced path to a bit-identical permutation/bootstrap draw stream,
    # which the tau_peak_s/dcf_peak/perm_p-only comparison above does not.
    assert dcf_rows_native == dcf_rows_string

    with open(info_s["csv_path"], newline="") as f:
        header = json.loads(f.readline()[len(_QU_DCF_HEADER_PREFIX):])
    # Every coerced value must have landed as its real, native JSON type --
    # a str would not survive json.dumps/json.loads as an int/float/None.
    assert header["lag_bin_s"] is None
    assert type(header["n_permutations"]) is int and header["n_permutations"] == 200
    assert type(header["n_bootstrap"]) is int and header["n_bootstrap"] == 200
    assert type(header["seed"]) is int and header["seed"] == 0
    assert type(header["variance_floor"]) is float
    assert type(header["coherent_p_threshold"]) is float
    assert type(header["max_lag_frac"]) is float

    assert _nan_aware_equal(res_native["population_correlations"], res_string["population_correlations"])


# --------------------------------------------------------------------------- #
# Capability: n_min_slices must genuinely gate, not be silently swallowed
# while the tool falls back to its own default.
# --------------------------------------------------------------------------- #

def test_n_min_slices_actually_gates(tmp_path_factory):
    tmp_a = str(tmp_path_factory.mktemp("qu_lag_gate_9"))
    _copy_qu_slices(tmp_a, OBSID_REAL)
    res_9 = m.qu_lag_dcf(tmp_a, [OBSID_REAL], os.path.join(tmp_a, "out"),
                          n_min_slices=9, n_permutations=20, n_bootstrap=20)
    assert res_9["exit_code"] == 0
    info_9 = res_9["per_subobs"][OBSID_REAL]
    assert info_9["status"] == "insufficient_slices"
    assert "tau_peak_s" not in info_9

    tmp_b = str(tmp_path_factory.mktemp("qu_lag_gate_8"))
    _copy_qu_slices(tmp_b, OBSID_REAL)
    res_8 = m.qu_lag_dcf(tmp_b, [OBSID_REAL], os.path.join(tmp_b, "out"),
                          n_min_slices=8, n_permutations=20, n_bootstrap=20)
    assert res_8["exit_code"] == 0
    info_8 = res_8["per_subobs"][OBSID_REAL]
    assert info_8["status"] == "ok"
    assert math.isfinite(info_8["tau_peak_s"])

    # Same check with n_min_slices delivered string-typed, exactly as the
    # sandbox harness would deliver it -- must gate identically.
    tmp_c = str(tmp_path_factory.mktemp("qu_lag_gate_9_str"))
    _copy_qu_slices(tmp_c, OBSID_REAL)
    res_9_str = m.qu_lag_dcf(tmp_c, [OBSID_REAL], os.path.join(tmp_c, "out"),
                              n_min_slices="9", n_permutations=20, n_bootstrap=20)
    assert res_9_str["per_subobs"][OBSID_REAL]["status"] == "insufficient_slices"


# --------------------------------------------------------------------------- #
# Capability: n_permutations, string-typed, must genuinely set the
# permutation null's own trial count (not merely be accepted and then
# ignored while the tool falls back to its default 1000).
# --------------------------------------------------------------------------- #

def test_n_permutations_actually_sets_permutation_count(tmp_path_factory):
    # perm_p is add-one smoothed: p = (count_ge + 1) / (n_permutations + 1),
    # so it must land exactly on the k/(n_permutations+1) grid. If the
    # string were silently ignored and the default 1000 used instead, p
    # would be k/1001 -- and k*8/1001 (or k*20/1001) is integral only at
    # k=1001 (p=1.0), which the p<1.0 guard below excludes.
    tmp_a = str(tmp_path_factory.mktemp("qu_lag_permcount_7"))
    _copy_qu_slices(tmp_a, OBSID_EVENSAMP)
    res_7 = m.qu_lag_dcf(tmp_a, [OBSID_EVENSAMP], os.path.join(tmp_a, "out"),
                          n_permutations="7", n_bootstrap="20", seed="0")
    assert res_7["exit_code"] == 0
    perm_p_7 = res_7["per_subobs"][OBSID_EVENSAMP]["perm_p"]
    assert 0.0 < perm_p_7 < 1.0
    scaled_7 = perm_p_7 * 8
    assert abs(scaled_7 - round(scaled_7)) < 1e-9

    tmp_b = str(tmp_path_factory.mktemp("qu_lag_permcount_19"))
    _copy_qu_slices(tmp_b, OBSID_EVENSAMP)
    res_19 = m.qu_lag_dcf(tmp_b, [OBSID_EVENSAMP], os.path.join(tmp_b, "out"),
                           n_permutations="19", n_bootstrap="20", seed="0")
    assert res_19["exit_code"] == 0
    perm_p_19 = res_19["per_subobs"][OBSID_EVENSAMP]["perm_p"]
    assert 0.0 < perm_p_19 < 1.0
    scaled_19 = perm_p_19 * 20
    assert abs(scaled_19 - round(scaled_19)) < 1e-9


# --------------------------------------------------------------------------- #
# Reproducibility: same seed (string vs native) -> identical tau_peak_s/
# dcf_peak/perm_p; different seeds must be free to disagree (else a seed
# silently ignored would pass the same-seed check by accident).
# --------------------------------------------------------------------------- #

def test_repro_qu_lag_dcf_seeded(tmp_path_factory):
    tmp_a = str(tmp_path_factory.mktemp("qu_lag_repro_string"))
    _copy_qu_slices(tmp_a, OBSID_EVENSAMP)
    res_string = m.qu_lag_dcf(tmp_a, [OBSID_EVENSAMP], os.path.join(tmp_a, "out"),
                               n_permutations="200", n_bootstrap="200", seed="0")

    tmp_b = str(tmp_path_factory.mktemp("qu_lag_repro_native"))
    _copy_qu_slices(tmp_b, OBSID_EVENSAMP)
    res_native = m.qu_lag_dcf(tmp_b, [OBSID_EVENSAMP], os.path.join(tmp_b, "out"),
                               n_permutations=200, n_bootstrap=200, seed=0)

    assert res_string["exit_code"] == 0 and res_native["exit_code"] == 0
    info_s = res_string["per_subobs"][OBSID_EVENSAMP]
    info_n = res_native["per_subobs"][OBSID_EVENSAMP]
    for field in ("tau_peak_s", "dcf_peak", "perm_p"):
        assert info_s[field] == pytest.approx(info_n[field], abs=1e-12)

    tmp_c = str(tmp_path_factory.mktemp("qu_lag_repro_seed_diff"))
    _copy_qu_slices(tmp_c, OBSID_EVENSAMP)
    res_diff_seed = m.qu_lag_dcf(tmp_c, [OBSID_EVENSAMP], os.path.join(tmp_c, "out"),
                                  n_permutations=200, n_bootstrap=200, seed=1)
    assert res_diff_seed["exit_code"] == 0
    info_diff = res_diff_seed["per_subobs"][OBSID_EVENSAMP]

    # A different seed is free to give a different empirical permutation
    # envelope -- proving the seed genuinely reaches the RNG rather than
    # being ignored, which would otherwise let the same-seed check above
    # pass by accident. perm_p itself is coarsely quantized (multiples of
    # 1/(n_permutations+1)) and can tie by chance, so compare the continuous
    # per-bin perm_lo99 envelope written to each subobs's own CSV instead.
    def _perm_lo99(info):
        with open(info["csv_path"], newline="") as f:
            f.readline()
            return [float(r["perm_lo99"]) for r in csv.DictReader(f)]

    assert _perm_lo99(info_diff) != _perm_lo99(info_n)


# --------------------------------------------------------------------------- #
# Edge case (brief edge case (d)): lag_bin_s=None, and every None-like str
# ('', 'none', 'null'), and the kwarg simply omitted, must all map to the
# SAME "use the median sampling interval" behaviour -- never to 0.0.
# --------------------------------------------------------------------------- #

def test_lag_bin_s_none_like_maps_to_median_not_zero(tmp_path_factory):
    variants = [
        ("omitted", {}),
        ("none_native", {"lag_bin_s": None}),
        ("empty_string", {"lag_bin_s": ""}),
        ("none_string", {"lag_bin_s": "none"}),
        ("null_string", {"lag_bin_s": "null"}),
    ]
    dcf_rows_by_label = {}
    for label, kwargs in variants:
        tmp_dir = str(tmp_path_factory.mktemp(f"qu_lag_none_{label}"))
        _copy_qu_slices(tmp_dir, OBSID_REAL)
        res = m.qu_lag_dcf(tmp_dir, [OBSID_REAL], os.path.join(tmp_dir, "out"),
                            n_permutations=20, n_bootstrap=20, seed=0, **kwargs)
        assert res["exit_code"] == 0, (label, res)
        info = res["per_subobs"][OBSID_REAL]
        assert info["status"] == "ok", label

        with open(info["csv_path"], newline="") as f:
            header_line = f.readline()
            data_rows = list(csv.DictReader(f))
        header = json.loads(header_line[len(_QU_DCF_HEADER_PREFIX):])
        assert header["lag_bin_s"] is None, f"{label}: lag_bin_s not coerced to None"

        with open(res["ccf_master_csv"], newline="") as f:
            master_row = next(r for r in csv.DictReader(f) if r["subobs_key"] == OBSID_REAL)
        lag_bin_used = float(master_row["lag_bin_s_used"])
        median_dt = float(master_row["median_dt_s"])
        assert lag_bin_used == pytest.approx(median_dt, rel=1e-12), (
            f"{label}: lag_bin_s_used ({lag_bin_used}) did not fall back to "
            f"median_dt_s ({median_dt}) -- never 0.0"
        )
        dcf_rows_by_label[label] = data_rows

    baseline = dcf_rows_by_label["omitted"]
    for label, rows in dcf_rows_by_label.items():
        assert rows == baseline, f"{label} produced qu_dcf rows different from the omitted-kwarg default"


# --------------------------------------------------------------------------- #
# Capability: lag_bin_s and max_lag_frac (both float-typed, string-typed
# here) must genuinely reach the binning grid, not be silently swallowed.
# --------------------------------------------------------------------------- #

def test_float_args_actually_change_binning(tmp_path_factory):
    tmp_base = str(tmp_path_factory.mktemp("qu_lag_float_base"))
    _copy_qu_slices(tmp_base, OBSID_REAL)
    res_base = m.qu_lag_dcf(tmp_base, [OBSID_REAL], os.path.join(tmp_base, "out"),
                             n_permutations=20, n_bootstrap=20, seed=0)
    assert res_base["exit_code"] == 0
    info_base = res_base["per_subobs"][OBSID_REAL]
    with open(res_base["ccf_master_csv"], newline="") as f:
        master_base = next(r for r in csv.DictReader(f) if r["subobs_key"] == OBSID_REAL)
    median_dt = float(master_base["median_dt_s"])
    assert float(master_base["lag_bin_s_used"]) == pytest.approx(median_dt, rel=1e-9)

    def _dcf_rows(info):
        with open(info["csv_path"], newline="") as f:
            f.readline()
            return list(csv.DictReader(f))

    rows_default = _dcf_rows(info_base)
    assert len(rows_default) > 1  # else the strict-inequality checks below would be vacuous

    # Coarse bin width, string-typed: 3x the median spacing must produce
    # strictly fewer bins than the default -- a silently-ignored float would
    # reproduce the default bin count exactly.
    tmp_coarse = str(tmp_path_factory.mktemp("qu_lag_float_coarse"))
    _copy_qu_slices(tmp_coarse, OBSID_REAL)
    coarse_lag_bin = 3.0 * median_dt
    res_coarse = m.qu_lag_dcf(tmp_coarse, [OBSID_REAL], os.path.join(tmp_coarse, "out"),
                               n_permutations=20, n_bootstrap=20, seed=0,
                               lag_bin_s=str(coarse_lag_bin))
    assert res_coarse["exit_code"] == 0
    info_coarse = res_coarse["per_subobs"][OBSID_REAL]
    with open(res_coarse["ccf_master_csv"], newline="") as f:
        master_coarse = next(r for r in csv.DictReader(f) if r["subobs_key"] == OBSID_REAL)
    assert float(master_coarse["lag_bin_s_used"]) == pytest.approx(coarse_lag_bin, rel=1e-9)
    rows_coarse = _dcf_rows(info_coarse)
    assert len(rows_coarse) < len(rows_default)

    # Narrower max_lag_frac, string-typed: the covered |lag| range must
    # shrink relative to the default 0.5 -- a silently-ignored float would
    # reproduce the default's max |lag| exactly.
    tmp_narrow = str(tmp_path_factory.mktemp("qu_lag_float_narrow"))
    _copy_qu_slices(tmp_narrow, OBSID_REAL)
    res_narrow = m.qu_lag_dcf(tmp_narrow, [OBSID_REAL], os.path.join(tmp_narrow, "out"),
                               n_permutations=20, n_bootstrap=20, seed=0,
                               max_lag_frac="0.25")
    assert res_narrow["exit_code"] == 0
    info_narrow = res_narrow["per_subobs"][OBSID_REAL]
    rows_narrow = _dcf_rows(info_narrow)
    assert len(rows_narrow) > 1 and len(rows_default) > 1
    max_abs_lag_default = max(abs(float(r["lag_s"])) for r in rows_default)
    max_abs_lag_narrow = max(abs(float(r["lag_s"])) for r in rows_narrow)
    assert max_abs_lag_narrow < max_abs_lag_default


# --------------------------------------------------------------------------- #
# Capability: exclude_noise_limited (boolean, string-typed) must genuinely
# reach the population-correlation filter, not be silently swallowed.
# --------------------------------------------------------------------------- #

def test_exclude_noise_limited_bool_actually_filters(tmp_path_factory):
    def _build_obs_db(label):
        tmp_dir = str(tmp_path_factory.mktemp(f"qu_lag_excl_{label}"))
        _copy_qu_slices(tmp_dir, OBSID_REAL)
        evensamp_csv = _copy_qu_slices(tmp_dir, OBSID_EVENSAMP)
        with open(evensamp_csv, newline="") as f:
            lines = f.readlines()
        n_rows = len(lines) - 2  # minus provenance header + column header
        assert n_rows >= 8  # stays 'ok' under this tool's own n_min_slices default
        # Bump the UPSTREAM qu_slice_pcubes provenance's own n_min above the
        # real row count, so THIS tool reconstructs noise_limited=True for
        # this subobs while its own n_min_slices gate (still satisfied) still
        # lets it through as status='ok' -- exactly the case exclude_noise_
        # limited is meant to filter out of the population correlation.
        _rewrite_qu_slice_n_min(evensamp_csv, n_rows + 1)
        return tmp_dir

    results = {}
    for label, kwargs in (
        ("true_native", {"exclude_noise_limited": True}),
        ("true_string", {"exclude_noise_limited": "true"}),
        ("false_native", {"exclude_noise_limited": False}),
        ("false_string", {"exclude_noise_limited": "false"}),
    ):
        tmp_dir = _build_obs_db(label)
        res = m.qu_lag_dcf(tmp_dir, [OBSID_REAL, OBSID_EVENSAMP], os.path.join(tmp_dir, "out"),
                            n_permutations=20, n_bootstrap=20, seed=0, **kwargs)
        assert res["exit_code"] == 0, (label, res)
        assert res["per_subobs"][OBSID_REAL]["status"] == "ok"
        assert res["per_subobs"][OBSID_EVENSAMP]["status"] == "ok"
        results[label] = res

    n_true_native = results["true_native"]["population_correlations"]["tau_peak_vs_r_proxy"]["n"]
    n_true_string = results["true_string"]["population_correlations"]["tau_peak_vs_r_proxy"]["n"]
    n_false_native = results["false_native"]["population_correlations"]["tau_peak_vs_r_proxy"]["n"]
    n_false_string = results["false_string"]["population_correlations"]["tau_peak_vs_r_proxy"]["n"]

    assert n_true_native == 1  # the noise_limited evensamp subobs dropped
    assert n_true_string == 1
    assert n_false_native == 2  # both subobs kept
    assert n_false_string == 2
    assert n_true_native != n_false_native  # must not be softened to a tie

    assert _nan_aware_equal(results["true_native"]["population_correlations"],
                             results["true_string"]["population_correlations"])
    assert _nan_aware_equal(results["false_native"]["population_correlations"],
                             results["false_string"]["population_correlations"])

    # The filter touches only the population correlation, never ccf_master's
    # per-subobs status.
    for res in results.values():
        assert res["per_subobs"][OBSID_REAL]["status"] == "ok"
        assert res["per_subobs"][OBSID_EVENSAMP]["status"] == "ok"


# --------------------------------------------------------------------------- #
# No-op guarantee: qu_lag_dcf_plot's max_panels coerces identically for
# string/native input and writes PNG only; qu_lag_dcf/qu_lag_dcf_plot with
# every numeric kwarg omitted still reproduce their committed baselines.
# --------------------------------------------------------------------------- #

def test_plotter_max_panels_and_defaults(tmp_path_factory):
    tmp_dir = str(tmp_path_factory.mktemp("qu_lag_plot_defaults"))
    _copy_qu_slices(tmp_dir, OBSID_REAL)
    out_dir = os.path.join(tmp_dir, "out")
    dcf_res = m.qu_lag_dcf(tmp_dir, [OBSID_REAL], out_dir, n_permutations=50, n_bootstrap=50)
    assert dcf_res["exit_code"] == 0
    assert dcf_res["per_subobs"][OBSID_REAL]["status"] == "ok"

    figs_dir = os.path.join(tmp_dir, "figs")
    res_native = m.qu_lag_dcf_plot(tmp_dir, dcf_res["ccf_master_csv"], figs_dir,
                                    save_path_prefix=os.path.join(figs_dir, "native"),
                                    max_panels=4)
    res_string = m.qu_lag_dcf_plot(tmp_dir, dcf_res["ccf_master_csv"], figs_dir,
                                    save_path_prefix=os.path.join(figs_dir, "string"),
                                    max_panels="4")
    assert res_native["exit_code"] == 0 and res_string["exit_code"] == 0
    assert len(res_native["plotted"]["dcf_panels"]) == len(res_string["plotted"]["dcf_panels"])
    assert set(res_native["plotted"]["dcf_panels"]) == set(res_string["plotted"]["dcf_panels"])

    for fig in res_native["figures"] + res_string["figures"]:
        assert fig["path"].endswith(".png")
        assert os.path.exists(fig["path"])

    # Defaults (every numeric kwarg omitted) still reproduce the committed
    # baseline: qu_lag_dcf's own per_subobs values against the pre-fix
    # committed qu_dcf CSV, and qu_lag_dcf_plot's default max_panels=12
    # covers every ok subobs (there is only one here).
    committed_rows = _read_committed_qu_dcf_rows(OBSID_REAL)
    tmp_full = str(tmp_path_factory.mktemp("qu_lag_plot_full_defaults"))
    _copy_qu_slices(tmp_full, OBSID_REAL)
    full_res = m.qu_lag_dcf(tmp_full, [OBSID_REAL], os.path.join(tmp_full, "out"))
    assert full_res["exit_code"] == 0
    with open(full_res["per_subobs"][OBSID_REAL]["csv_path"], newline="") as f:
        f.readline()
        full_rows = list(csv.DictReader(f))
    committed_by_lag = {round(float(r["lag_s"]), 3): float(r["dcf"]) for r in committed_rows}
    for row in full_rows:
        assert float(row["dcf"]) == pytest.approx(committed_by_lag[round(float(row["lag_s"]), 3)], abs=1e-9)

    default_plot = m.qu_lag_dcf_plot(tmp_full, full_res["ccf_master_csv"],
                                      os.path.join(tmp_full, "figs"))
    assert default_plot["exit_code"] == 0
    assert set(default_plot["plotted"]["dcf_panels"]) == {OBSID_REAL}


# --------------------------------------------------------------------------- #
# Argument-error-not-data-status: an un-coercible numeric argument is a
# caller error, never a per-subobs status row.
# --------------------------------------------------------------------------- #

def test_bad_argument_raises_not_statused(tmp_path):
    _copy_qu_slices(str(tmp_path), OBSID_REAL)
    out_dir = os.path.join(str(tmp_path), "out")

    with pytest.raises(ValueError, match="n_permutations"):
        m.qu_lag_dcf(str(tmp_path), [OBSID_REAL], out_dir, n_permutations="many")

    assert not os.path.exists(os.path.join(out_dir, "ccf_master.csv"))
    assert not os.path.exists(out_dir)
    assert neutral.read(str(tmp_path))["entries"] == {}


def test_bad_max_panels_raises_before_reading_ccf_master(tmp_path_factory):
    tmp_dir = str(tmp_path_factory.mktemp("qu_lag_plot_bad_max_panels"))
    _copy_qu_slices(tmp_dir, OBSID_REAL)
    out_dir = os.path.join(tmp_dir, "out")
    dcf_res = m.qu_lag_dcf(tmp_dir, [OBSID_REAL], out_dir, n_permutations=20, n_bootstrap=20)
    assert dcf_res["exit_code"] == 0

    with pytest.raises(ValueError, match="max_panels"):
        m.qu_lag_dcf_plot(tmp_dir, dcf_res["ccf_master_csv"], os.path.join(tmp_dir, "figs"),
                           max_panels="many")
    assert not os.path.exists(os.path.join(tmp_dir, "figs"))
