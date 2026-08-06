"""Tests for the ixpe_qu_slice_pcubes sandbox numeric/boolean argument
marshalling fix (fix-ixpe-qu-slice-pcubes-numeric-args).

The sandbox harness invokes a tool as a bare ``result = fn(**args)`` with no
schema coercion, so a numeric/boolean TOOL_MANIFEST argument declared
integer/number/boolean can arrive at ``ixpe_qu_slice_pcubes`` as a Python
``str``. Before this fix every explicit numeric keyword crashed with a
``TypeError`` deep inside ``red_ixpe.qu_slice_pcubes`` (str-vs-int/float
comparisons, a str*float multiply), so a caller could only ever run the
tool on its built-in defaults. These tests pin the fix on the real,
already-local Cyg X-1 IXPE ObsID 01002901 (the ObsID the reported defect was
confirmed on), copied into an isolated tmp directory per fixture so the
reduction runs on a scratch copy rather than the shared analyst-env data
directory. ``n_max``/``n_min`` are kept small (2-4), not the tool's
production default (25), to keep the real xpselect/xpbin call count low --
the same convention ``test_ixpe_qu_slice_pcubes.py`` already uses. The two
monkeypatch tests and the fast-fail test share one copy of the real event
files (``shared_event_dir``) since none of the three ever reaches a real
reduction, so no additional real-data cost is added beyond the fixtures
already needed for the regression/capability/no-op/reproducibility checks.

No-op guarantee, at the interface rather than against the committed CSV
(finding, recorded here rather than re-derived by a future reader): the
brief originally framed the no-op check as "reproduce the committed
on-disk ``qu_slices_01002901.csv`` (8 slices) with every numeric kwarg
omitted". That committed artifact is NOT a valid baseline -- its own
provenance header records ``"bkg": "false"`` (a JSON *string*), i.e. it was
produced through the very defect this fix repairs: a raw harness string is
truthy, so ``red_ixpe._decide_background``'s ``bkg if bkg is not None else
...`` returns it unchanged and ``_reduce_window``'s ``if bkg:`` gate
subtracted a background annulus even though the caller asked for none
(``bkg=False``). A defaults run through the FIXED wrapper resolves
``bkg=None`` -> ``_decide_background``'s own auto-decision -> for a
bright, source-dominated pointing like Cyg X-1 01002901, no subtraction --
which legitimately changes ``MDP99_obs``, hence ``C_target``, hence every
slice edge relative to that committed file. Reproducing it would therefore
assert the *defective* numbers as ground truth. The no-op guarantee is
enforced instead where it is exact and unambiguous:
``test_defaults_reach_library_unchanged`` (the 15 defaults reach
``red_ixpe.qu_slice_pcubes`` with identical values AND identical types,
nothing dropped or invented) and ``test_defaults_unchanged_for_native_typed_call``
(a real-data differential: the coercing wrapper reproduces a direct,
non-coercing ``red_ixpe.qu_slice_pcubes`` call's ``n_slice``, whole-sub-obs
PD/PA and every per-slice q/u/PD/PA/t_start/t_stop).
"""

from __future__ import annotations

import csv
import glob
import json
import os
import shutil
import time as _time

import pytest

from obsdb_catalog import neutral
from xspec_utils.lib import red_ixpe
from xspec_utils.lib.toolargs import coerce_tool_args
from xspec_utils.skills import session as ixpe_session

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
OBSID = "01002901"
_PROV_PREFIX = "# qu_slice_pcubes provenance: "


def _copy_real_event_files(dest_dir: str, obsid: str = OBSID) -> list[str]:
    """Copy the real, already-local Level-2 DU1/2/3 event files for
    ``obsid`` into ``<dest_dir>/<obsid>/event_l2/`` and return the sorted
    file list -- never mutates the shared source directory."""
    src_evt_dir = os.path.join(REAL_IXPE_DIR, obsid, "event_l2")
    dst_evt_dir = os.path.join(dest_dir, obsid, "event_l2")
    os.makedirs(dst_evt_dir, exist_ok=True)
    file_list = []
    for det in (1, 2, 3):
        rej = sorted(glob.glob(os.path.join(src_evt_dir, f"ixpe{obsid}_det{det}_evt2_v[0-9][0-9]_rej.fits")))
        plain = sorted(glob.glob(os.path.join(src_evt_dir, f"ixpe{obsid}_det{det}_evt2_v[0-9][0-9].fits")))
        src = rej[0] if rej else (plain[0] if plain else None)
        if src is None:
            pytest.skip(f"real fixture ObsID {obsid} DU{det} not local under {src_evt_dir}")
        dst = os.path.join(dst_evt_dir, os.path.basename(src))
        shutil.copyfile(src, dst)
        file_list.append(dst)
    return sorted(file_list)


def _read_slice_rows(csv_path: str) -> list[dict]:
    with open(csv_path, newline="") as f:
        f.readline()  # provenance comment header
        return list(csv.DictReader(f))


def _read_provenance(csv_path: str) -> dict:
    with open(csv_path, newline="") as f:
        line = f.readline()
    assert line.startswith(_PROV_PREFIX), f"{csv_path}: missing qu_slice_pcubes provenance header"
    return json.loads(line[len(_PROV_PREFIX):])


# --------------------------------------------------------------------------- #
# Fixtures: one real reduction per configuration, shared across the tests
# that only read its result.
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def native_result(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("qu_numeric_native")
    _copy_real_event_files(str(tmp_dir))
    return ixpe_session.ixpe_qu_slice_pcubes(
        str(tmp_dir), OBSID, n_max=3, src_rad=1.0, mdp_frac=0.5, bkg_inner_rad=3.0,
    )


@pytest.fixture(scope="module")
def string_fixture(tmp_path_factory):
    """(result, obs_db_path) for one real reduction called with every
    numeric kwarg string-typed exactly as the sandbox harness delivers
    them. Exposes obs_db_path (unlike a plain result-only fixture) so the
    catalog-record test can read the same obs-db's catalog.json."""
    tmp_dir = tmp_path_factory.mktemp("qu_numeric_string")
    _copy_real_event_files(str(tmp_dir))
    result = ixpe_session.ixpe_qu_slice_pcubes(
        str(tmp_dir), OBSID, n_max="3", src_rad="1.0", mdp_frac="0.5", bkg_inner_rad="3.0",
    )
    return result, str(tmp_dir)


@pytest.fixture(scope="module")
def string_result(string_fixture):
    return string_fixture[0]


@pytest.fixture(scope="module")
def direct_lib_result(tmp_path_factory):
    """Calls red_ixpe.qu_slice_pcubes directly, bypassing the session
    wrapper (and therefore the coerce_tool_args step entirely) with the same
    native-typed arguments as ``native_result`` -- the baseline the wrapper
    must reproduce exactly for a well-typed call."""
    tmp_dir = tmp_path_factory.mktemp("qu_numeric_direct")
    file_list = _copy_real_event_files(str(tmp_dir))
    return red_ixpe.qu_slice_pcubes(
        OBSID, file_list, n_max=3, src_rad=1.0, mdp_frac=0.5, bkg_inner_rad=3.0,
    )


@pytest.fixture(scope="module")
def shared_event_dir(tmp_path_factory):
    """One copy of the real event files, shared by every test below that
    never reaches a real reduction (the two monkeypatch-capture tests, and
    the fast-fail test which is rejected at the coercion step). None of
    them writes anything into this directory."""
    tmp_dir = tmp_path_factory.mktemp("qu_numeric_shared")
    _copy_real_event_files(str(tmp_dir))
    return str(tmp_dir)


# --------------------------------------------------------------------------- #
# Regression: the reported defect -- string-typed numerics, exactly as the
# sandbox harness delivers them, must succeed and agree with native-typed
# equivalents, INCLUDING the on-disk provenance header (the assertion that
# fails if a string ever survived into the library: json would have
# serialised it as "3", not 3).
# --------------------------------------------------------------------------- #

def test_numeric_kwarg_string_vs_native_slicing(native_result, string_result):
    assert native_result["exit_code"] == 0, native_result.get("stderr")
    assert string_result["exit_code"] == 0, string_result.get("stderr")

    info_native = native_result["per_subobs"][OBSID]
    info_string = string_result["per_subobs"][OBSID]
    assert info_native["n_slice"] == info_string["n_slice"]

    rows_native = _read_slice_rows(info_native["csv_path"])
    rows_string = _read_slice_rows(info_string["csv_path"])
    assert len(rows_native) == len(rows_string)
    for row_a, row_b in zip(rows_native, rows_string):
        for field in ("t_start_met", "t_stop_met"):
            assert float(row_a[field]) == pytest.approx(float(row_b[field]), rel=1e-9)

    prov_native = _read_provenance(info_native["csv_path"])
    prov_string = _read_provenance(info_string["csv_path"])

    assert prov_string["n_max"] == 3 and isinstance(prov_string["n_max"], int)
    assert prov_string["src_rad"] == 1.0 and isinstance(prov_string["src_rad"], float)
    assert prov_string["mdp_frac"] == 0.5 and isinstance(prov_string["mdp_frac"], float)
    assert prov_string["bkg_inner_rad"] == 3.0 and isinstance(prov_string["bkg_inner_rad"], float)

    for field in ("n_min", "n_max", "emin_kev", "emax_kev", "src_rad", "mdp_frac",
                  "mdp_floor", "bkg_inner_rad", "bkg_outer_rad", "split_gap_factor"):
        assert prov_native[field] == prov_string[field]
        assert type(prov_native[field]) is type(prov_string[field])


# --------------------------------------------------------------------------- #
# Capability: n_max must genuinely change the slicing, not be silently
# swallowed while the tool falls back to its own default -- and the on-disk
# provenance header must record the real n_max actually used.
# --------------------------------------------------------------------------- #

def test_explicit_n_max_changes_slicing(tmp_path_factory):
    tmp_a = tmp_path_factory.mktemp("qu_numeric_nmax2")
    _copy_real_event_files(str(tmp_a))
    result_2 = ixpe_session.ixpe_qu_slice_pcubes(str(tmp_a), OBSID, n_max=2, n_min=1)
    assert result_2["exit_code"] == 0, result_2.get("stderr")

    tmp_b = tmp_path_factory.mktemp("qu_numeric_nmax4")
    _copy_real_event_files(str(tmp_b))
    result_4 = ixpe_session.ixpe_qu_slice_pcubes(str(tmp_b), OBSID, n_max=4, n_min=1)
    assert result_4["exit_code"] == 0, result_4.get("stderr")

    info_2 = result_2["per_subobs"][OBSID]
    info_4 = result_4["per_subobs"][OBSID]
    assert info_2["n_slice"] <= 2
    assert info_4["n_slice"] <= 4
    assert info_2["n_slice"] != info_4["n_slice"]

    edges_2 = [r["t_start_met"] for r in _read_slice_rows(info_2["csv_path"])]
    edges_4 = [r["t_start_met"] for r in _read_slice_rows(info_4["csv_path"])]
    assert edges_2 != edges_4

    prov_2 = _read_provenance(info_2["csv_path"])
    prov_4 = _read_provenance(info_4["csv_path"])
    assert prov_2["n_max"] == 2 and isinstance(prov_2["n_max"], int)
    assert prov_4["n_max"] == 4 and isinstance(prov_4["n_max"], int)


# --------------------------------------------------------------------------- #
# No-op guarantee: coercion must be an identity for already-native-typed
# input -- the session wrapper (which now coerces) must reproduce the
# underlying lib call (which never coerces) exactly.
# --------------------------------------------------------------------------- #

def test_defaults_unchanged_for_native_typed_call(native_result, direct_lib_result):
    assert native_result["exit_code"] == 0, native_result.get("stderr")

    info_w = native_result["per_subobs"][OBSID]
    info_d = direct_lib_result["per_subobs"][OBSID]
    assert info_w["n_slice"] == info_d["n_slice"]
    assert info_w["whole_pd_pct"] == pytest.approx(info_d["whole_pd_pct"], rel=1e-6)
    assert info_w["whole_pa_deg"] == pytest.approx(info_d["whole_pa_deg"], rel=1e-6)

    rows_w = _read_slice_rows(info_w["csv_path"])
    rows_d = info_d["rows"]
    assert len(rows_w) == len(rows_d)
    for row_w, row_d in zip(rows_w, rows_d):
        for field in ("q", "u", "PD_pct", "PA_deg", "t_start_met", "t_stop_met"):
            assert float(row_w[field]) == pytest.approx(row_d[field], rel=1e-6)


# --------------------------------------------------------------------------- #
# Reproducibility: the newly-working string-typed call path itself gives a
# stable core result across independent calls.
# --------------------------------------------------------------------------- #

def test_repro_string_typed_call_stable_across_calls(string_result, tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("qu_numeric_string_repro")
    _copy_real_event_files(str(tmp_dir))
    second = ixpe_session.ixpe_qu_slice_pcubes(
        str(tmp_dir), OBSID, n_max="3", src_rad="1.0", mdp_frac="0.5", bkg_inner_rad="3.0",
    )
    assert second["exit_code"] == 0, second.get("stderr")

    info_a = string_result["per_subobs"][OBSID]
    info_b = second["per_subobs"][OBSID]
    assert info_a["n_slice"] == info_b["n_slice"]
    for field in ("whole_pd_pct", "whole_pa_deg"):
        assert info_a[field] == pytest.approx(info_b[field], rel=1e-6)

    rows_a = _read_slice_rows(info_a["csv_path"])
    rows_b = _read_slice_rows(info_b["csv_path"])
    for row_a, row_b in zip(rows_a, rows_b):
        for field in ("q", "u", "PD_pct", "PA_deg"):
            assert float(row_a[field]) == pytest.approx(float(row_b[field]), rel=1e-6)


# --------------------------------------------------------------------------- #
# Catalog provenance: the recorded params must carry the COERCED (native)
# values, never the raw string the caller passed -- this fails against the
# pre-fix wrapper, which built the catalog record from the raw un-coerced
# parameters rather than from `coerced`.
# --------------------------------------------------------------------------- #

def test_catalog_record_params_are_native_typed(string_fixture):
    result, obs_db_path = string_fixture
    assert result["exit_code"] == 0, result.get("stderr")

    entry = neutral.read(obs_db_path)["entries"][OBSID]
    products = [p for p in entry["products"] if p["kind"] == "qu_slices"]
    assert products, "no qu_slices catalog product recorded"
    for p in products:
        assert isinstance(p["params"]["src_rad"], float)
        assert p["params"]["src_rad"] == 1.0
        assert isinstance(p["params"]["emin_kev"], float)
        assert isinstance(p["params"]["emax_kev"], float)


# --------------------------------------------------------------------------- #
# Interface boundary, no reduction: the exact set of coerced parameters that
# reaches red_ixpe.qu_slice_pcubes, captured via monkeypatch before the real
# function ever runs.
# --------------------------------------------------------------------------- #

_EXPECTED_DEFAULTS = {
    "ra": None, "dec": None, "emin_kev": 2.0, "emax_kev": 8.0, "src_rad": 1.0,
    "bkg": None, "bkg_inner_rad": 3.0, "bkg_outer_rad": 5.0, "mdp_target": None,
    "mdp_frac": 0.5, "mdp_floor": 0.015, "n_min": 8, "n_max": 25,
    "split_gap_factor": 10.0, "split_multiobs": None,
}


def _assert_exact_type(value, expected):
    if expected is None:
        assert value is None
    elif isinstance(expected, bool):
        assert type(value) is bool
    elif isinstance(expected, int):
        assert type(value) is int
    elif isinstance(expected, float):
        assert type(value) is float


def _capture_and_raise(captured):
    def fake_qu_slice_pcubes(obs, file_list, **kwargs):
        captured.update(kwargs)
        raise RuntimeError("CAPTURED")
    return fake_qu_slice_pcubes


def test_defaults_reach_library_unchanged(shared_event_dir, monkeypatch):
    captured = {}
    monkeypatch.setattr(ixpe_session.red_ixpe, "qu_slice_pcubes", _capture_and_raise(captured))

    result = ixpe_session.ixpe_qu_slice_pcubes(shared_event_dir, OBSID)

    assert result["exit_code"] == 1
    assert "CAPTURED" in result.get("stderr", "")
    assert set(captured) == set(_EXPECTED_DEFAULTS)
    for key, expected in _EXPECTED_DEFAULTS.items():
        assert captured[key] == expected
        _assert_exact_type(captured[key], expected)


def test_string_kwargs_reach_library_as_native_types(shared_event_dir, monkeypatch):
    captured = {}
    monkeypatch.setattr(ixpe_session.red_ixpe, "qu_slice_pcubes", _capture_and_raise(captured))

    result = ixpe_session.ixpe_qu_slice_pcubes(
        shared_event_dir, OBSID,
        ra="299.59", dec="35.2", emin_kev="2.0", emax_kev="8.0", src_rad="1.0",
        bkg="false", bkg_inner_rad="2.5", bkg_outer_rad="4.0", mdp_target="0.06",
        mdp_frac="0.5", mdp_floor="0.015", n_min="8", n_max="40",
        split_gap_factor="10.0", split_multiobs="false",
    )

    assert result["exit_code"] == 1
    assert "CAPTURED" in result.get("stderr", "")

    assert captured["ra"] == pytest.approx(299.59) and type(captured["ra"]) is float
    assert captured["dec"] == pytest.approx(35.2) and type(captured["dec"]) is float
    assert captured["emin_kev"] == 2.0 and type(captured["emin_kev"]) is float
    assert captured["emax_kev"] == 8.0 and type(captured["emax_kev"]) is float
    assert captured["src_rad"] == 1.0 and type(captured["src_rad"]) is float
    assert captured["bkg"] is False
    assert isinstance(captured["bkg"], bool)
    assert captured["bkg_inner_rad"] == 2.5 and type(captured["bkg_inner_rad"]) is float
    assert captured["bkg_outer_rad"] == 4.0 and type(captured["bkg_outer_rad"]) is float
    assert captured["mdp_target"] == 0.06 and type(captured["mdp_target"]) is float
    assert captured["mdp_frac"] == 0.5 and type(captured["mdp_frac"]) is float
    assert captured["mdp_floor"] == 0.015 and type(captured["mdp_floor"]) is float
    assert captured["n_min"] == 8 and type(captured["n_min"]) is int
    assert captured["n_max"] == 40 and type(captured["n_max"]) is int
    assert captured["split_gap_factor"] == 10.0 and type(captured["split_gap_factor"]) is float
    assert captured["split_multiobs"] is False
    assert isinstance(captured["split_multiobs"], bool)


# --------------------------------------------------------------------------- #
# Fast-fail: an un-coercible argument must be rejected before any file I/O.
# --------------------------------------------------------------------------- #

def test_bad_argument_fails_before_io(shared_event_dir):
    evt_dir = os.path.join(shared_event_dir, OBSID, "event_l2")
    before = sorted(glob.glob(os.path.join(evt_dir, "*")))

    t0 = _time.monotonic()
    result = ixpe_session.ixpe_qu_slice_pcubes(shared_event_dir, OBSID, n_min="eight")
    elapsed = _time.monotonic() - t0

    assert result["exit_code"] == 1
    assert "n_min" in result.get("stderr", "")
    assert elapsed < 5.0

    after = sorted(glob.glob(os.path.join(evt_dir, "*")))
    assert after == before  # no wholesrc/pcube product written during coercion
    assert glob.glob(os.path.join(shared_event_dir, OBSID, "qu_slices_*.csv")) == []


# --------------------------------------------------------------------------- #
# Science-level consequence, asserted against the real library code (no
# mock, no reduction): a raw harness string "false" is truthy in Python, so
# red_ixpe._decide_background propagates it unchanged and _reduce_window's
# `if bkg:` gate subtracts a background annulus contrary to the caller's
# explicit bkg=False request -- exactly the path the committed
# qu_slices_01002901.csv went through (its provenance header records
# "bkg": "false", a JSON string, not JSON false). Coercing bkg through
# coerce_tool_args before it ever reaches the library closes this.
# --------------------------------------------------------------------------- #

def test_string_false_bkg_flips_the_background_decision():
    # Pre-fix hazard, against the real (unmocked) library helper: the raw
    # harness string "false" is truthy, so _decide_background's override
    # branch (`bkg if bkg is not None else ...`) returns it unchanged --
    # not the False the caller intended.
    assert red_ixpe._decide_background("false", 5.0) is not False
    assert bool(red_ixpe._decide_background("false", 5.0)) is True

    # Post-fix: coerce_tool_args hands the library a real bool, and the
    # decision the library then makes is the one the caller actually asked for.
    assert coerce_tool_args({"bkg": "bool?"}, {"bkg": "false"})["bkg"] is False
    assert red_ixpe._decide_background(False, 5.0) is False
    assert red_ixpe._decide_background(
        coerce_tool_args({"bkg": "bool?"}, {"bkg": "false"})["bkg"], 5.0
    ) is False
