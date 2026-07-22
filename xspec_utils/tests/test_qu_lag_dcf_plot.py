"""Tests for xspec_utils' ixpe_qu_lag_dcf_plot (`lib/qu_lag_dcf.py:qu_lag_dcf_plot`).

Builds real `ccf_master.csv` / `qu_dcf_<subobs_key>.csv` / `population_
correlations.csv` fixtures by running the real `qu_lag_dcf` tool over the
three `qu_slices_<subobs_key>.csv` products already on disk under
/data/fzhou/XSPEC/Cyg_X-1/IXPE (the data catalog confirms all three local):

- 01002901: the real, already-local Cyg X-1 IXPE ObsID.
- synth_qu_lag_evensamp01, synth_qu_lag_null01: the two marked validation
  ObsIDs `test_qu_lag_dcf.py` also uses, carrying their own ground truth in
  the `qu_slice_pcubes provenance` header.

Copied into an isolated tmp obs_db_path per fixture, so the plotter's own
reads never touch the shared analyst-env mission directory. `qu_lag_dcf_plot`
itself writes only PNGs -- no CSV, no catalog entry.
"""

from __future__ import annotations

import ast
import csv
import inspect
import json
import math
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pytest

from xspec_utils.lib import qu_lag_dcf as m
from xspec_utils.skills import TOOL_MANIFEST

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
OBSID_REAL = "01002901"
OBSID_EVENSAMP = "synth_qu_lag_evensamp01"
OBSID_NULL = "synth_qu_lag_null01"
ALL_OBSIDS = [OBSID_REAL, OBSID_EVENSAMP, OBSID_NULL]
SUFFIXES = ("dcf_panels", "population_tau_r_proxy", "population_tau_pd_signif")


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


@pytest.fixture(scope="module")
def real_ccf(tmp_path_factory):
    """One genuine `qu_lag_dcf` run over all three fixture ObsIDs: a real
    `ccf_master.csv`, one `qu_dcf_<subobs_key>.csv` per obsid, and a real
    `population_correlations.csv` -- exactly the products `qu_lag_dcf_plot`
    reads."""
    obs_db = str(tmp_path_factory.mktemp("qu_lag_dcf_plot_fixture"))
    for obsid in ALL_OBSIDS:
        _copy_qu_slices(obs_db, obsid)
    out_dir = os.path.join(obs_db, "out")
    res = m.qu_lag_dcf(obs_db, ALL_OBSIDS, out_dir, n_permutations=50, n_bootstrap=50, seed=0)
    assert res["exit_code"] == 0
    assert all(res["per_subobs"][o]["status"] == "ok" for o in ALL_OBSIDS)
    return {
        "obs_db": obs_db, "out_dir": out_dir,
        "ccf_master_csv": res["ccf_master_csv"], "population_csv": res["population_csv"],
    }


def _by_suffix(figs):
    return {s: f["sha256"] for f in figs for s in SUFFIXES if f["path"].endswith(f"_{s}.png")}


# --------------------------------------------------------------------------- #
# Known-input regression: the exact defect that sank the prior authoring
# attempt -- a bootstrap 68% interval that EXCLUDES tau_peak_s must be drawn
# by its own true endpoints, never a max(0, ...)-clipped span anchored at
# tau_peak_s.
# --------------------------------------------------------------------------- #

def _write_min_ccf_master(path, obsid="00000001", subobs_key="00000001", overrides=None):
    base = {f: "" for f in m.CCF_MASTER_FIELDS}
    base.update(
        obsid=obsid, source="src", subobs_key=subobs_key, n_slices=10, status="ok",
        tau_peak_s=1.0, dcf_peak=0.5, dcf_zero=0.1, n_pairs_peak=20,
        dcf_uncorrected_denominator=False, perm_p=0.005,
        tau_boot_lo=5.0, tau_boot_hi=10.0, dcf_peak_boot_lo=0.3, dcf_peak_boot_hi=0.7,
        coherent_rotation=True, r_proxy=2.0, pd_significance=3.0,
        tau_peak_uniform_s=1.0, median_dt_s=100.0, lag_bin_s_used=100.0,
    )
    base.update(overrides or {})
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=m.CCF_MASTER_FIELDS)
        writer.writeheader()
        writer.writerow(base)


def _write_min_qu_dcf(path, obsid="00000001", subobs_key="00000001"):
    provenance = {"obsid": obsid, "subobs_key": subobs_key, "n_slices": 10}
    with open(path, "w", newline="") as f:
        f.write(m._QU_DCF_HEADER_PREFIX + json.dumps(provenance) + "\n")
        writer = csv.DictWriter(f, fieldnames=m._QU_DCF_FIELDS)
        writer.writeheader()
        for lag, dcf in ((-100.0, 0.1), (0.0, 0.2), (1.0, 0.5), (100.0, 0.3)):
            writer.writerow({"lag_s": lag, "dcf": dcf, "n_pairs": 5,
                              "perm_lo99": -0.4, "perm_hi99": 0.4})


def _write_min_population_csv(path, n=1, rho=0.5, p=0.1):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=m.POPULATION_FIELDS)
        writer.writeheader()
        writer.writerow({"pair": "tau_peak_vs_r_proxy", "n": n, "rho": rho, "p": p})
        writer.writerow({"pair": "tau_peak_vs_pd_significance", "n": n, "rho": rho, "p": p})


def _minimal_fixture_dir(tmp_path):
    """A hand-built ccf_master.csv row exercising the bootstrap-interval-
    excludes-tau_peak edge case, plus its matching qu_dcf CSV and sibling
    population_correlations.csv -- a minimal control-flow fixture, not a
    measurement."""
    obs_db = tmp_path / "obs_db"
    os.makedirs(obs_db / "00000001", exist_ok=True)
    ccf_master_csv = obs_db / "ccf_master.csv"
    _write_min_ccf_master(str(ccf_master_csv))
    _write_min_qu_dcf(str(obs_db / "00000001" / "qu_dcf_00000001.csv"))
    _write_min_population_csv(str(obs_db / "population_correlations.csv"))
    return str(obs_db), str(ccf_master_csv)


def test_bootstrap_interval_drawn_by_true_endpoints_not_distance_clipped(tmp_path):
    obs_db, ccf_master_csv = _minimal_fixture_dir(tmp_path)

    result = m.qu_lag_dcf_plot(obs_db, ccf_master_csv, str(tmp_path / "figs"),
                                save_path_prefix=str(tmp_path / "figs" / "plot"))
    assert result["exit_code"] == 0
    panel = result["plotted"]["dcf_panels"]["00000001"]
    assert panel["tau_boot_lo"] == 5.0
    assert panel["tau_boot_hi"] == 10.0
    assert panel["tau_peak_s"] == 1.0
    assert panel["tau_boot_lo"] > panel["tau_peak_s"]  # interval excludes tau_peak_s

    # The drawn segment's own x-data on the axis must be exactly [5.0, 10.0]
    # -- not derived from max(0, tau_peak-tau_boot_lo)/max(0, tau_boot_hi-tau_peak)
    # anchored at tau_peak_s=1.0 (that defect would draw [1.0, 1.0] or clip a
    # negative span to zero).
    ok_rows = [r for r in m._read_ccf_master_csv(ccf_master_csv) if r["status"] == "ok"]
    fig, axes, _plotted = m._build_dcf_panels_figure(obs_db, ok_rows, max_panels=12)
    try:
        ax = axes[0]
        interval_lines = [ln for ln in ax.lines
                           if ln.get_color() == "tab:red" and len(ln.get_xdata()) == 2]
        assert len(interval_lines) == 1
        xdata = sorted(interval_lines[0].get_xdata())
        assert xdata == pytest.approx([5.0, 10.0])
    finally:
        plt.close(fig)


# --------------------------------------------------------------------------- #
# Capability: dcf_panels + population plotted arrays equal an INDEPENDENT
# direct read of the qu_dcf_*.csv columns and the ccf_master row values.
# --------------------------------------------------------------------------- #

def test_capability_plotted_matches_independent_reads(tmp_path, real_ccf):
    result = m.qu_lag_dcf_plot(real_ccf["obs_db"], real_ccf["ccf_master_csv"],
                                str(tmp_path / "figs"),
                                save_path_prefix=str(tmp_path / "figs" / "plot"))
    assert result["exit_code"] == 0
    assert len(result["figures"]) == 3

    with open(real_ccf["ccf_master_csv"], newline="") as f:
        ccf_rows = {r["subobs_key"]: r for r in csv.DictReader(f) if r["status"] == "ok"}
    assert set(ccf_rows) == set(ALL_OBSIDS)

    plotted = result["plotted"]
    assert set(plotted["dcf_panels"]) == set(ALL_OBSIDS)
    for subobs_key, panel in plotted["dcf_panels"].items():
        row = ccf_rows[subobs_key]
        assert panel["tau_peak_s"] == pytest.approx(float(row["tau_peak_s"]), rel=1e-9, abs=1e-12)
        assert panel["dcf_peak"] == pytest.approx(float(row["dcf_peak"]), rel=1e-9, abs=1e-12)
        assert panel["tau_boot_lo"] == pytest.approx(float(row["tau_boot_lo"]), rel=1e-9, abs=1e-12)
        assert panel["tau_boot_hi"] == pytest.approx(float(row["tau_boot_hi"]), rel=1e-9, abs=1e-12)

        dcf_csv = os.path.join(real_ccf["obs_db"], subobs_key, f"qu_dcf_{subobs_key}.csv")
        with open(dcf_csv, newline="") as f:
            f.readline()  # provenance header
            raw_rows = list(csv.DictReader(f))
        expected_lag = sorted(float(r["lag_s"]) for r in raw_rows)
        assert panel["lag_s"] == pytest.approx(expected_lag, rel=1e-9)
        expected_dcf_by_lag = {round(float(r["lag_s"]), 6): float(r["dcf"]) for r in raw_rows}
        for lag_s, dcf_val in zip(panel["lag_s"], panel["dcf"]):
            assert dcf_val == pytest.approx(expected_dcf_by_lag[round(lag_s, 6)], rel=1e-9, abs=1e-12)

    for suffix, covariate in (("population_tau_r_proxy", "r_proxy"),
                               ("population_tau_pd_signif", "pd_significance")):
        pop = plotted[suffix]
        assert set(pop["subobs_key"]) == set(ALL_OBSIDS)
        for i, subobs_key in enumerate(pop["subobs_key"]):
            row = ccf_rows[subobs_key]
            assert pop["y"][i] == pytest.approx(float(row["tau_peak_s"]), rel=1e-9, abs=1e-12)
            assert pop["x"][i] == pytest.approx(float(row[covariate]), rel=1e-9, abs=1e-12)
            assert pop["tau_boot_lo"][i] == pytest.approx(float(row["tau_boot_lo"]), rel=1e-9, abs=1e-12)
            assert pop["tau_boot_hi"][i] == pytest.approx(float(row["tau_boot_hi"]), rel=1e-9, abs=1e-12)


# --------------------------------------------------------------------------- #
# Differential: each population panel's spearman{n,rho,p} equals an
# independent direct read of population_correlations.csv.
# --------------------------------------------------------------------------- #

def test_differential_spearman_annotation_matches_population_csv(tmp_path, real_ccf):
    result = m.qu_lag_dcf_plot(real_ccf["obs_db"], real_ccf["ccf_master_csv"],
                                str(tmp_path / "figs"),
                                save_path_prefix=str(tmp_path / "figs" / "plot"))
    assert result["exit_code"] == 0

    with open(real_ccf["population_csv"], newline="") as f:
        pop_rows = {r["pair"]: r for r in csv.DictReader(f)}

    for suffix, pair in (("population_tau_r_proxy", "tau_peak_vs_r_proxy"),
                         ("population_tau_pd_signif", "tau_peak_vs_pd_significance")):
        spearman = result["plotted"][suffix]["spearman"]
        expected = pop_rows[pair]
        assert spearman["n"] == int(float(expected["n"]))
        exp_rho = float(expected["rho"])
        if math.isnan(exp_rho):
            assert math.isnan(spearman["rho"])
        else:
            assert spearman["rho"] == pytest.approx(exp_rho, rel=1e-9, abs=1e-12)
            assert spearman["p"] == pytest.approx(float(expected["p"]), rel=1e-9, abs=1e-12)


# --------------------------------------------------------------------------- #
# Reproducibility: two identical calls -> byte-identical PNGs for all three
# suffixes.
# --------------------------------------------------------------------------- #

def test_repro_byte_identical_pngs(tmp_path, real_ccf):
    a = m.qu_lag_dcf_plot(real_ccf["obs_db"], real_ccf["ccf_master_csv"], str(tmp_path / "figs_a"),
                          save_path_prefix=str(tmp_path / "figs_a" / "plot"))
    b = m.qu_lag_dcf_plot(real_ccf["obs_db"], real_ccf["ccf_master_csv"], str(tmp_path / "figs_b"),
                          save_path_prefix=str(tmp_path / "figs_b" / "plot"))
    assert a["exit_code"] == 0 and b["exit_code"] == 0

    sha_a, sha_b = _by_suffix(a["figures"]), _by_suffix(b["figures"])
    assert set(sha_a) == set(SUFFIXES) == set(sha_b)
    for s in SUFFIXES:
        assert sha_a[s] == sha_b[s], f"figure {s} not byte-identical across repeated calls"


# --------------------------------------------------------------------------- #
# Hard-error contract: bad ccf_master_csv / zero ok rows / bad sibling
# population CSV -> exit_code=1, stderr set, no files written.
# --------------------------------------------------------------------------- #

def test_missing_ccf_master_csv_hard_error(tmp_path):
    missing = str(tmp_path / "does_not_exist.csv")
    result = m.qu_lag_dcf_plot(str(tmp_path), missing, str(tmp_path / "figs"))
    assert result["exit_code"] == 1
    assert result["stderr"]
    assert "figures" not in result
    assert not (tmp_path / "figs").exists()


def test_wrong_schema_ccf_master_csv_hard_error(tmp_path):
    bad = tmp_path / "bad_ccf_master.csv"
    with open(bad, "w", newline="") as f:
        f.write("not_obsid,foo\n1,2\n")
    result = m.qu_lag_dcf_plot(str(tmp_path), str(bad), str(tmp_path / "figs"))
    assert result["exit_code"] == 1
    assert result["stderr"]
    assert "figures" not in result
    assert not (tmp_path / "figs").exists()


def test_zero_ok_rows_hard_error(tmp_path):
    obs_db = str(tmp_path / "obs_db")
    os.makedirs(obs_db, exist_ok=True)
    missing_obsid = "09999999"  # no qu_slices_*.csv anywhere -> status='no_qu_slices'
    out_dir = os.path.join(obs_db, "out")
    res = m.qu_lag_dcf(obs_db, [missing_obsid], out_dir, n_permutations=5, n_bootstrap=5)
    assert res["exit_code"] == 0
    assert res["per_subobs"][missing_obsid]["status"] == "no_qu_slices"

    result = m.qu_lag_dcf_plot(obs_db, res["ccf_master_csv"], str(tmp_path / "figs"))
    assert result["exit_code"] == 1
    assert "ok" in result["stderr"]
    assert "figures" not in result
    assert not (tmp_path / "figs").exists()


def test_missing_population_csv_hard_error(tmp_path, real_ccf):
    isolated_dir = tmp_path / "isolated_missing_pop"
    isolated_dir.mkdir()
    isolated_ccf = isolated_dir / "ccf_master.csv"
    shutil.copyfile(real_ccf["ccf_master_csv"], isolated_ccf)
    result = m.qu_lag_dcf_plot(real_ccf["obs_db"], str(isolated_ccf), str(tmp_path / "figs"))
    assert result["exit_code"] == 1
    assert "population_correlations" in result["stderr"]
    assert "figures" not in result
    assert not (tmp_path / "figs").exists()


def test_wrong_schema_population_csv_hard_error(tmp_path, real_ccf):
    isolated_dir = tmp_path / "isolated_bad_pop"
    isolated_dir.mkdir()
    isolated_ccf = isolated_dir / "ccf_master.csv"
    shutil.copyfile(real_ccf["ccf_master_csv"], isolated_ccf)
    with open(isolated_dir / "population_correlations.csv", "w", newline="") as f:
        f.write("not_pair,foo\n1,2\n")
    result = m.qu_lag_dcf_plot(real_ccf["obs_db"], str(isolated_ccf), str(tmp_path / "figs"))
    assert result["exit_code"] == 1
    assert "population_correlations" in result["stderr"]
    assert "figures" not in result
    assert not (tmp_path / "figs").exists()


def test_ok_row_missing_qu_dcf_csv_raises_valueerror(tmp_path, real_ccf):
    isolated_dir = tmp_path / "isolated_no_qu_dcf"
    isolated_dir.mkdir()
    isolated_ccf = isolated_dir / "ccf_master.csv"
    shutil.copyfile(real_ccf["ccf_master_csv"], isolated_ccf)
    shutil.copyfile(real_ccf["population_csv"], isolated_dir / "population_correlations.csv")
    # obs_db_path points at a directory with no <obsid>/qu_dcf_*.csv at all.
    with pytest.raises(ValueError, match="qu_dcf"):
        m.qu_lag_dcf_plot(str(isolated_dir), str(isolated_ccf), str(tmp_path / "figs"))


# --------------------------------------------------------------------------- #
# Mixed-catalog robustness: a genuine 'ok' row plus a genuine non-'ok' nan
# row (both produced by the real qu_lag_dcf) must never crash int-field
# parsing ahead of the status filter.
# --------------------------------------------------------------------------- #

def test_mixed_catalog_robustness_int_fields_never_crash_parsing(tmp_path):
    obs_db = str(tmp_path / "obs_db")
    _copy_qu_slices(obs_db, OBSID_REAL)
    missing_obsid = "09999999"
    out_dir = os.path.join(obs_db, "out")
    res = m.qu_lag_dcf(obs_db, [OBSID_REAL, missing_obsid], out_dir,
                        n_permutations=20, n_bootstrap=20)
    assert res["exit_code"] == 0
    assert res["per_subobs"][OBSID_REAL]["status"] == "ok"
    assert res["per_subobs"][missing_obsid]["status"] == "no_qu_slices"

    result = m.qu_lag_dcf_plot(obs_db, res["ccf_master_csv"], str(tmp_path / "figs"),
                                save_path_prefix=str(tmp_path / "figs" / "plot"))
    assert result["exit_code"] == 0
    assert len(result["figures"]) == 3
    assert list(result["plotted"]["dcf_panels"]) == [OBSID_REAL]
    assert result["plotted"]["population_tau_r_proxy"]["subobs_key"] == [OBSID_REAL]


# --------------------------------------------------------------------------- #
# max_panels bounds dcf_panels only, never the population panels.
# --------------------------------------------------------------------------- #

def test_max_panels_bounds_dcf_panels_not_population(tmp_path, real_ccf):
    result = m.qu_lag_dcf_plot(real_ccf["obs_db"], real_ccf["ccf_master_csv"],
                                str(tmp_path / "figs"),
                                save_path_prefix=str(tmp_path / "figs" / "plot"), max_panels=1)
    assert result["exit_code"] == 0
    assert len(result["plotted"]["dcf_panels"]) == 1
    assert len(result["plotted"]["population_tau_r_proxy"]["subobs_key"]) == len(ALL_OBSIDS)
    assert len(result["plotted"]["population_tau_pd_signif"]["subobs_key"]) == len(ALL_OBSIDS)


# --------------------------------------------------------------------------- #
# Manifest wiring (AI-facing half) + no science-stack import
# --------------------------------------------------------------------------- #

def test_tool_manifest_wired_to_lib_directly_no_science_stack():
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.ixpe_qu_lag_dcf_plot"]
    assert entry["module"] == "xspec_utils.lib.qu_lag_dcf"
    assert entry["fn"] == "qu_lag_dcf_plot"
    props = entry["args_schema"]["properties"]
    for name in ("obs_db_path", "ccf_master_csv", "output_dir", "save_path_prefix", "max_panels"):
        assert name in props, f"{name} missing from args_schema.properties"
    assert entry["args_schema"]["required"] == ["obs_db_path", "ccf_master_csv", "output_dir"]


def test_lib_module_imports_no_xspec_or_ixpeobssim_or_red_ixpe():
    tree = ast.parse(inspect.getsource(m))
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    assert "xspec" not in imported
    assert "ixpeobssim" not in imported
    assert "red_ixpe" not in imported
