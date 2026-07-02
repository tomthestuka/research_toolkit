"""Tests for xspec_utils' ixpe_pa_pd_stability_plot (`pa_pd_stability.stability_plot`).

Uses the two real, already-local pcube_timeseg_<obsid>.csv products (01002901,
03006001) and the real pa_pd_stability_catalog.csv already on disk under
/data/fzhou/XSPEC/Cyg_X-1/IXPE -- both written by the already-shipped
ixpe_pcube_timeseg / ixpe_pa_pd_stability tools. stability_plot only reads
these; tmp_path is used solely for save_path_prefix output.
"""

from __future__ import annotations

import csv
import hashlib
import math
import os

import numpy as np
import pytest

from xspec_utils.lib import pa_pd_stability as pps
from xspec_utils.skills import TOOL_MANIFEST

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
REAL_CATALOG_CSV = os.path.join(REAL_IXPE_DIR, "pa_pd_stability_catalog.csv")
OBSID_A = "01002901"   # N_seg=92 -- well above the default sf_min_segments=10
OBSID_B = "03006001"   # N_seg=59


def _require_real_fixtures():
    for obsid in (OBSID_A, OBSID_B):
        path = os.path.join(REAL_IXPE_DIR, obsid, f"pcube_timeseg_{obsid}.csv")
        if not os.path.exists(path):
            pytest.skip(f"real fixture pcube_timeseg CSV not local at {path}")
    if not os.path.exists(REAL_CATALOG_CSV):
        pytest.skip(f"real pa_pd_stability_catalog.csv not local at {REAL_CATALOG_CSV}")


SUFFIXES = ("pa_pd_timeseries", "stokes_plane", "chi2_hist", "structure_function")


# --------------------------------------------------------------------------- #
# Capability: exactly four named, non-empty files, sha256 matches file bytes
# --------------------------------------------------------------------------- #

def test_produces_four_named_files(tmp_path):
    _require_real_fixtures()
    prefix = str(tmp_path / "stability")
    result = pps.stability_plot(REAL_IXPE_DIR, [OBSID_A, OBSID_B], REAL_CATALOG_CSV,
                                 save_path_prefix=prefix)
    assert result["exit_code"] == 0
    assert len(result["figures"]) == 4

    seen_suffixes = set()
    for fig in result["figures"]:
        path = fig["path"]
        assert os.path.exists(path)
        data = open(path, "rb").read()
        assert len(data) > 0
        assert fig["sha256"] == hashlib.sha256(data).hexdigest()
        assert fig["mime_type"] == "image/png"
        assert fig["produced_by"] == "xspec_utils.lib.pa_pd_stability.stability_plot"
        assert fig["caption"]
        matched = [s for s in SUFFIXES if path == f"{prefix}_{s}.png"]
        assert len(matched) == 1, f"unexpected figure path {path}"
        seen_suffixes.add(matched[0])
    assert seen_suffixes == set(SUFFIXES)

    # No stray files beyond the four expected ones at this prefix.
    written = sorted(p for p in os.listdir(tmp_path) if p.startswith("stability_"))
    assert len(written) == 4


def test_pdf_format_written_with_correct_mime(tmp_path):
    _require_real_fixtures()
    prefix = str(tmp_path / "stability_pdf")
    result = pps.stability_plot(REAL_IXPE_DIR, [OBSID_A], REAL_CATALOG_CSV,
                                 save_path_prefix=prefix, fmt="pdf")
    assert result["exit_code"] == 0
    assert len(result["figures"]) == 4
    for fig in result["figures"]:
        assert fig["path"].endswith(".pdf")
        assert fig["mime_type"] == "application/pdf"
        assert os.path.getsize(fig["path"]) > 0


# --------------------------------------------------------------------------- #
# Reproducibility: identical calls -> byte-identical figures
# --------------------------------------------------------------------------- #

def test_repro_stability_plot(tmp_path_factory):
    _require_real_fixtures()
    out_a = str(tmp_path_factory.mktemp("splot_a") / "run")
    out_b = str(tmp_path_factory.mktemp("splot_b") / "run")
    a = pps.stability_plot(REAL_IXPE_DIR, [OBSID_A, OBSID_B], REAL_CATALOG_CSV, save_path_prefix=out_a)
    b = pps.stability_plot(REAL_IXPE_DIR, [OBSID_A, OBSID_B], REAL_CATALOG_CSV, save_path_prefix=out_b)

    def _by_suffix(figs):
        return {s: f["sha256"] for f in figs for s in SUFFIXES if f["path"].endswith(f"_{s}.png")}

    sha_a, sha_b = _by_suffix(a["figures"]), _by_suffix(b["figures"])
    assert set(sha_a) == set(SUFFIXES) == set(sha_b)
    for s in SUFFIXES:
        assert sha_a[s] == sha_b[s], f"figure {s} not byte-identical across repeated calls"


# --------------------------------------------------------------------------- #
# Differential: plotted structure-function values vs an independent recompute
# via the trusted pps.structure_function primitive (never the plotter's own
# helper as the expected side).
# --------------------------------------------------------------------------- #

def test_sf_panel_matches_lib(tmp_path):
    _require_real_fixtures()
    csv_path = os.path.join(REAL_IXPE_DIR, OBSID_A, f"pcube_timeseg_{OBSID_A}.csv")
    provenance, rows_by_segid, _integrated = pps._read_timeseg_csv(csv_path)
    arr = pps._rows_to_arrays(rows_by_segid)

    mid_s = 0.5 * (arr["t_start_mjd"] + arr["t_stop_mjd"]) * 86400.0
    dt = np.diff(np.sort(mid_s))
    half_width_s = float(np.median(dt) * 0.5)
    lags_ks = [5.0, 10.0, 20.0]
    lags_s = [ks * 1000.0 for ks in lags_ks]
    expected_sf_pa = pps.structure_function(mid_s, arr["PA_deg"], lags_s, half_width_s, wrap=True)
    expected_sf_pd = pps.structure_function(mid_s, arr["PD_pct"], lags_s, half_width_s, wrap=False)

    result = pps.stability_plot(REAL_IXPE_DIR, [OBSID_A], REAL_CATALOG_CSV,
                                 save_path_prefix=str(tmp_path / "sf"), sf_lags_ks=lags_ks)
    plotted = result["structure_function"][OBSID_A]
    for lag in lags_s:
        va, vb = plotted["sf_pa"][lag], expected_sf_pa[lag]
        if math.isnan(vb):
            assert math.isnan(va)
        else:
            assert va == pytest.approx(vb, rel=1e-9, abs=1e-9)
        va, vb = plotted["sf_pd"][lag], expected_sf_pd[lag]
        if math.isnan(vb):
            assert math.isnan(va)
        else:
            assert va == pytest.approx(vb, rel=1e-9, abs=1e-9)


def test_sf_panel_omits_obsid_below_min_segments(tmp_path):
    """03006001 has 59 segments; requiring sf_min_segments=1000 must drop it
    from the structure-function data/panel while the other three figures
    still cover it."""
    _require_real_fixtures()
    result = pps.stability_plot(REAL_IXPE_DIR, [OBSID_B], REAL_CATALOG_CSV,
                                 save_path_prefix=str(tmp_path / "sf_gate"),
                                 sf_min_segments=1000)
    assert result["exit_code"] == 0
    assert result["structure_function"] == {}
    assert len(result["figures"]) == 4  # all four still written


# --------------------------------------------------------------------------- #
# Edge cases: missing CSVs
# --------------------------------------------------------------------------- #

def test_missing_csv_skipped(tmp_path):
    _require_real_fixtures()
    result = pps.stability_plot(REAL_IXPE_DIR, [OBSID_A, "99999999"], REAL_CATALOG_CSV,
                                 save_path_prefix=str(tmp_path / "mixed"))
    assert result["exit_code"] == 0
    assert len(result["figures"]) == 4
    assert "99999999" in result["stdout"]
    assert "skipped" in result["stdout"].lower()
    for fig in result["figures"]:
        assert os.path.getsize(fig["path"]) > 0


def test_all_missing_obsids_exit_code_1_no_figures(tmp_path):
    prefix = str(tmp_path / "allmissing")
    result = pps.stability_plot(REAL_IXPE_DIR, ["99999999", "88888888"], REAL_CATALOG_CSV,
                                 save_path_prefix=prefix)
    assert result["exit_code"] == 1
    assert "stderr" in result and result["stderr"]
    assert "figures" not in result
    assert not any(p.startswith(os.path.basename(prefix)) for p in os.listdir(tmp_path))


def test_missing_catalog_csv_fails_loudly(tmp_path):
    """A nonexistent catalog_csv (required upstream product) must fail with
    exit_code=1 and write no figures -- not silently degrade to empty
    overlays / empty chi2 histogram."""
    _require_real_fixtures()
    prefix = str(tmp_path / "nocat")
    missing = str(tmp_path / "does_not_exist_catalog.csv")
    result = pps.stability_plot(REAL_IXPE_DIR, [OBSID_A], missing,
                                 save_path_prefix=prefix)
    assert result["exit_code"] == 1
    assert "stderr" in result and result["stderr"]
    assert "figures" not in result
    assert not any(p.startswith(os.path.basename(prefix)) for p in os.listdir(tmp_path))


def test_unreadable_catalog_csv_fails_loudly(tmp_path):
    """A catalog_csv that exists but has the wrong schema (no 'obsid'
    column) must fail with exit_code=1 and write no figures."""
    _require_real_fixtures()
    bad = tmp_path / "bad_schema_catalog.csv"
    with open(bad, "w", newline="") as f:
        f.write("not_obsid,foo\n1,2\n")
    prefix = str(tmp_path / "badcat")
    result = pps.stability_plot(REAL_IXPE_DIR, [OBSID_A], str(bad),
                                 save_path_prefix=prefix)
    assert result["exit_code"] == 1
    assert "stderr" in result and result["stderr"]
    assert "figures" not in result
    assert not any(p.startswith(os.path.basename(prefix)) for p in os.listdir(tmp_path))


def test_empty_catalog_csv_fails_loudly(tmp_path):
    """A completely empty catalog_csv (zero bytes, no header) has no 'obsid'
    column -> exit_code=1, no figures written."""
    _require_real_fixtures()
    empty = tmp_path / "empty_bytes_catalog.csv"
    empty.write_text("")
    prefix = str(tmp_path / "emptycat")
    result = pps.stability_plot(REAL_IXPE_DIR, [OBSID_A], str(empty),
                                 save_path_prefix=prefix)
    assert result["exit_code"] == 1
    assert "stderr" in result and result["stderr"]
    assert "figures" not in result
    assert not any(p.startswith(os.path.basename(prefix)) for p in os.listdir(tmp_path))


def test_catalog_with_obsid_but_missing_plot_fields_fails_before_writing(tmp_path):
    """A catalog with an 'obsid' column but missing plotter-required fields
    (PA_deg_int, rchi2_PD, ...) must fail up front (exit_code=1) and write
    NO figures -- not raise mid-plot after some files exist."""
    _require_real_fixtures()
    bad = tmp_path / "obsid_only_catalog.csv"
    with open(bad, "w", newline="") as f:
        # header has obsid + N_seg but none of the plotter's overlay/hist fields
        f.write("obsid,N_seg\n")
        f.write(f"{OBSID_A},92\n")
    prefix = str(tmp_path / "partialcat")
    result = pps.stability_plot(REAL_IXPE_DIR, [OBSID_A], str(bad),
                                 save_path_prefix=prefix)
    assert result["exit_code"] == 1
    assert "stderr" in result and result["stderr"]
    assert "figures" not in result
    # Critical: no partial figures on disk despite a valid obsid column.
    assert not any(p.startswith(os.path.basename(prefix)) for p in os.listdir(tmp_path))


def test_obsid_absent_from_catalog_still_plots(tmp_path):
    """An ObsID with a real pcube_timeseg CSV but no row in catalog_csv must
    still get its time-series/Stokes-plane panels -- only the weighted-mean
    overlay is skipped for it."""
    _require_real_fixtures()
    empty_catalog = tmp_path / "empty_catalog.csv"
    with open(empty_catalog, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=pps.CATALOG_FIELDS)
        writer.writeheader()

    result = pps.stability_plot(REAL_IXPE_DIR, [OBSID_A], str(empty_catalog),
                                 save_path_prefix=str(tmp_path / "nocatalog"))
    assert result["exit_code"] == 0
    assert len(result["figures"]) == 4
    for fig in result["figures"]:
        assert os.path.getsize(fig["path"]) > 0


# --------------------------------------------------------------------------- #
# Manifest wiring (AI-facing half) + no science-stack import
# --------------------------------------------------------------------------- #

def test_tool_manifest_wired_to_lib_directly_no_science_stack():
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.ixpe_pa_pd_stability_plot"]
    assert entry["module"] == "xspec_utils.lib.pa_pd_stability"
    assert entry["fn"] == "stability_plot"
    props = entry["args_schema"]["properties"]
    for name in ("obs_db_path", "obsids", "catalog_csv", "save_path_prefix",
                 "sf_lags_ks", "sf_min_segments", "fmt"):
        assert name in props, f"{name} missing from args_schema.properties"
    assert entry["args_schema"]["required"] == ["obs_db_path", "obsids", "catalog_csv"]
    assert props["fmt"]["enum"] == ["png", "pdf"]


def test_invalid_fmt_rejected(tmp_path):
    _require_real_fixtures()
    with pytest.raises(ValueError):
        pps.stability_plot(REAL_IXPE_DIR, [OBSID_A], REAL_CATALOG_CSV,
                            save_path_prefix=str(tmp_path / "bad"), fmt="svg")
