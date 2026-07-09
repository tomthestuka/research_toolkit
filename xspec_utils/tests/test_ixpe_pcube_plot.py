"""Tests for xspec_utils' ixpe_pcube_plot (`lib/ixpe_pcube_plot.py`).

Uses the real, already-local Level-2 event files for IXPE ObsIDs 01002901
(Cyg X-1, May 2022) and 01250101 (Cyg X-1, June 2022) -- this brief's own
validation pair -- copied into isolated tmp obs_db_path directories, then
runs the real `ixpe_pcube` skill (xpselect/xpbin) to produce a genuine
`ixpe_pcube_<obsid>.csv` before handing it to the pure-matplotlib plotter
under test. `ixpe_pcube_plot` itself imports neither `xspec` nor
`ixpeobssim` (only `ixpe_pcube`, generating the fixture CSVs, needs them).
"""

from __future__ import annotations

import csv
import hashlib
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.collections import PathCollection

from xspec_utils.lib import ixpe_pcube_plot as plotmod
from xspec_utils.skills import session as ixpe_session
from xspec_utils.skills import TOOL_MANIFEST

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
OBSID_A = "01002901"   # Cyg X-1, May 2022
OBSID_B = "01250101"   # Cyg X-1, June 2022
VERSION = "v01"

RA, DEC = 299.5971468, 35.1922507
IRF_DATE = "20220702"


def _copy_real_event_files(dest_dir: str, obsid: str, dus=(1, 2, 3)) -> list[str]:
    src_evt_dir = os.path.join(REAL_IXPE_DIR, obsid, "event_l2")
    dst_evt_dir = os.path.join(dest_dir, obsid, "event_l2")
    os.makedirs(dst_evt_dir, exist_ok=True)
    file_list = []
    for det in dus:
        src = os.path.join(src_evt_dir, f"ixpe{obsid}_det{det}_evt2_{VERSION}_rej.fits")
        if not os.path.exists(src):
            pytest.skip(f"real fixture ObsID {obsid} not local at {src}")
        dst = os.path.join(dst_evt_dir, os.path.basename(src))
        shutil.copyfile(src, dst)
        file_list.append(dst)
    return sorted(file_list)


@pytest.fixture(scope="module")
def multiband_csv(tmp_path_factory):
    """One real ixpe_pcube reduction of 01002901, default ebinning (4 bands)."""
    tmp_dir = tmp_path_factory.mktemp("ixpe_pcube_plot_multiband")
    obs_db_path = str(tmp_dir)
    _copy_real_event_files(obs_db_path, OBSID_A)
    result = ixpe_session.ixpe_pcube(obs_db_path, OBSID_A, ra=RA, dec=DEC, irf_date=IRF_DATE)
    assert result["exit_code"] == 0, result.get("stderr")
    assert len(result["bands"]) == 4
    return result["csv_path"]


@pytest.fixture(scope="module")
def singleband_csv(tmp_path_factory):
    """One real ixpe_pcube reduction of 01250101, single broadband edge pair."""
    tmp_dir = tmp_path_factory.mktemp("ixpe_pcube_plot_singleband")
    obs_db_path = str(tmp_dir)
    _copy_real_event_files(obs_db_path, OBSID_B)
    result = ixpe_session.ixpe_pcube(obs_db_path, OBSID_B, ra=RA, dec=DEC,
                                      ebinning=[2.0, 8.0], irf_date=IRF_DATE)
    assert result["exit_code"] == 0, result.get("stderr")
    assert len(result["bands"]) == 1
    return result["csv_path"]


def _read_raw_csv_bands(csv_path: str) -> list[dict]:
    with open(csv_path, newline="") as f:
        f.readline()  # provenance header
        return list(csv.DictReader(f))


# --------------------------------------------------------------------------- #
# Capability: exactly one PNG, non-zero size, correct sha256, FigureRef path
# --------------------------------------------------------------------------- #

def test_multiband_produces_one_png(tmp_path, multiband_csv):
    save_path = str(tmp_path / "plot.png")
    result = plotmod.ixpe_pcube_plot(multiband_csv, save_path=save_path)
    assert result["exit_code"] == 0
    assert len(result["figures"]) == 1
    fig = result["figures"][0]
    assert fig["path"] == save_path
    assert os.path.exists(save_path)
    data = open(save_path, "rb").read()
    assert len(data) > 0
    assert fig["sha256"] == hashlib.sha256(data).hexdigest()
    assert fig["mime_type"] == "image/png"
    assert fig["produced_by"] == "xspec_utils.lib.ixpe_pcube_plot.ixpe_pcube_plot"
    assert fig["caption"]
    assert result["plotted"]["n_bands"] == 4
    assert result["plotted"]["n_skipped"] == 0

    # No stray files beyond the requested one.
    assert sorted(os.listdir(tmp_path)) == ["plot.png"]


def test_singleband_renders_without_crash(tmp_path, singleband_csv):
    save_path = str(tmp_path / "single.png")
    result = plotmod.ixpe_pcube_plot(singleband_csv, save_path=save_path)
    assert result["exit_code"] == 0
    assert len(result["figures"]) == 1
    assert os.path.getsize(save_path) > 0
    assert result["plotted"]["n_bands"] == 1
    assert len(result["plotted"]["E_ctr"]) == 1


# --------------------------------------------------------------------------- #
# Differential: plotted arrays equal the CSV columns the ixpe_pcube tool wrote
# --------------------------------------------------------------------------- #

def test_plotted_arrays_match_source_csv(tmp_path, multiband_csv):
    raw_bands = _read_raw_csv_bands(multiband_csv)
    expected_e_ctr = sorted(0.5 * (float(r["E_LO"]) + float(r["E_HI"])) for r in raw_bands)
    expected_pd = {round(0.5 * (float(r["E_LO"]) + float(r["E_HI"])), 9): float(r["PD"]) for r in raw_bands}

    result = plotmod.ixpe_pcube_plot(multiband_csv, save_path=str(tmp_path / "diff.png"))
    plotted = result["plotted"]

    assert len(plotted["E_ctr"]) == len(raw_bands)
    for e_ctr, pd_val in zip(plotted["E_ctr"], plotted["PD"]):
        assert pd_val == pytest.approx(expected_pd[round(e_ctr, 9)], rel=1e-12, abs=1e-12)
    assert sorted(plotted["E_ctr"]) == pytest.approx(expected_e_ctr)

    # QN/UN likewise read back unchanged from the CSV.
    expected_qn = {round(0.5 * (float(r["E_LO"]) + float(r["E_HI"])), 9): float(r["QN"]) for r in raw_bands}
    for e_ctr, qn_val in zip(plotted["E_ctr"], plotted["QN"]):
        assert qn_val == pytest.approx(expected_qn[round(e_ctr, 9)], rel=1e-12, abs=1e-12)


# --------------------------------------------------------------------------- #
# Capability: the drawn artists carry the actual science arrays -- not just
# the `plotted` dict, which is assembled independently of the errorbar/
# scatter calls and would not catch a bug that dropped or corrupted the
# drawing itself.
# --------------------------------------------------------------------------- #

def _vertical_spans(barlinecols):
    """Given an errorbar's `barlinecols` (a tuple of LineCollections, one per
    error-bar direction), return the y-error bars' `(ymin, ymax)` spans --
    the segments whose two endpoints share an x -- sorted by that shared x
    so the result lines up with the energy-sorted plotted arrays. A
    LineCollection carrying x-error bars (endpoints sharing y, differing x)
    is skipped entirely."""
    spans = []
    for coll in barlinecols:
        for (x0, y0), (x1, y1) in coll.get_segments():
            if np.isclose(x0, x1):
                spans.append((x0, min(y0, y1), max(y0, y1)))
    spans.sort(key=lambda s: s[0])
    return [(lo, hi) for _, lo, hi in spans]


def test_axes_carry_the_science_data(multiband_csv):
    fig, axes, plotted, _label = plotmod._build_figure(multiband_csv)
    try:
        n_bands = plotted["n_bands"]

        # PD panel: the errorbar's central data line carries (E_ctr, PD).
        pd_container = axes["pd"].containers[0]
        pd_dataline = pd_container[0]
        assert np.allclose(pd_dataline.get_xdata(), plotted["E_ctr"])
        assert np.allclose(pd_dataline.get_ydata(), plotted["PD"])

        # ... and the 1-sigma yerr bars were actually drawn to +/- PD_err.
        pd_yerr_spans = _vertical_spans(pd_container[2])
        assert len(pd_yerr_spans) == n_bands
        for (lo, hi), y, e in zip(pd_yerr_spans, plotted["PD"], plotted["PD_err"]):
            assert lo == pytest.approx(y - e, abs=1e-9)
            assert hi == pytest.approx(y + e, abs=1e-9)

        # PA panel: central data line carries (E_ctr, PA) and its own yerr.
        pa_container = axes["pa"].containers[0]
        pa_dataline = pa_container[0]
        assert np.allclose(pa_dataline.get_xdata(), plotted["E_ctr"])
        assert np.allclose(pa_dataline.get_ydata(), plotted["PA"])
        pa_yerr_spans = _vertical_spans(pa_container[2])
        assert len(pa_yerr_spans) == n_bands
        for (lo, hi), y, e in zip(pa_yerr_spans, plotted["PA"], plotted["PA_err"]):
            assert lo == pytest.approx(y - e, abs=1e-9)
            assert hi == pytest.approx(y + e, abs=1e-9)

        # Q-U panel: the scatter offsets are exactly (QN, UN).
        path_collections = [c for c in axes["qu"].collections if isinstance(c, PathCollection)]
        assert len(path_collections) == 1
        offsets = path_collections[0].get_offsets()
        assert offsets.shape[0] == n_bands
        assert np.allclose(offsets[:, 0], plotted["QN"])
        assert np.allclose(offsets[:, 1], plotted["UN"])
    finally:
        plt.close(fig)


# --------------------------------------------------------------------------- #
# Edge cases: missing/NaN bands skipped (not zero); empty input raises
# --------------------------------------------------------------------------- #

def test_nan_band_skipped_not_plotted_as_zero(tmp_path, multiband_csv):
    """A pcube CSV fixture with one band row blanked out (control-flow-only,
    not a science measurement) must be skipped, never plotted as PD=0."""
    raw_bands = _read_raw_csv_bands(multiband_csv)
    with open(multiband_csv, newline="") as f:
        header_line = f.readline()

    bad_csv = str(tmp_path / "one_nan_band.csv")
    fieldnames = ["band_id", "E_LO", "E_HI", "PD", "PD_err", "PA", "PA_err",
                  "QN", "QN_err", "UN", "UN_err", "MDP99", "SIGNIF"]
    with open(bad_csv, "w", newline="") as f:
        f.write(header_line)
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i, r in enumerate(raw_bands):
            row = {"band_id": i, **{k: r[k] for k in fieldnames if k != "band_id"}}
            if i == 0:
                row["PD"] = ""  # blank out one band's PD -- must be skipped, not zeroed
            writer.writerow(row)

    result = plotmod.ixpe_pcube_plot(bad_csv, save_path=str(tmp_path / "nanplot.png"))
    assert result["plotted"]["n_bands"] == len(raw_bands) - 1
    assert result["plotted"]["n_skipped"] == 1
    assert 0.0 not in result["plotted"]["PD"]


def test_all_bands_nan_raises(tmp_path, multiband_csv):
    raw_bands = _read_raw_csv_bands(multiband_csv)
    with open(multiband_csv, newline="") as f:
        header_line = f.readline()

    bad_csv = str(tmp_path / "all_nan.csv")
    fieldnames = ["band_id", "E_LO", "E_HI", "PD", "PD_err", "PA", "PA_err",
                  "QN", "QN_err", "UN", "UN_err", "MDP99", "SIGNIF"]
    with open(bad_csv, "w", newline="") as f:
        f.write(header_line)
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i, r in enumerate(raw_bands):
            row = {"band_id": i, **{k: r[k] for k in fieldnames if k != "band_id"}}
            row["PD"] = ""
            writer.writerow(row)

    with pytest.raises(ValueError, match="finite"):
        plotmod.ixpe_pcube_plot(bad_csv)


def test_empty_csv_raises_clear_message(tmp_path):
    empty = str(tmp_path / "empty.csv")
    with open(empty, "w") as f:
        f.write("")
    with pytest.raises(ValueError, match="provenance header"):
        plotmod.ixpe_pcube_plot(empty)


def test_header_only_csv_raises_clear_message(tmp_path, multiband_csv):
    with open(multiband_csv, newline="") as f:
        header_line = f.readline()
    header_only = str(tmp_path / "header_only.csv")
    with open(header_only, "w", newline="") as f:
        f.write(header_line)
        f.write("band_id,E_LO,E_HI,PD,PD_err,PA,PA_err,QN,QN_err,UN,UN_err,MDP99,SIGNIF\n")
    with pytest.raises(ValueError, match="no band rows"):
        plotmod.ixpe_pcube_plot(header_only)


def test_missing_csv_raises_clear_message(tmp_path):
    missing = str(tmp_path / "does_not_exist.csv")
    with pytest.raises(ValueError, match="not found"):
        plotmod.ixpe_pcube_plot(missing)


# --------------------------------------------------------------------------- #
# Reproducibility: identical input -> byte-identical figure
# --------------------------------------------------------------------------- #

def test_repro_byte_identical(tmp_path, multiband_csv):
    out_a = str(tmp_path / "a.png")
    out_b = str(tmp_path / "b.png")
    a = plotmod.ixpe_pcube_plot(multiband_csv, save_path=out_a)
    b = plotmod.ixpe_pcube_plot(multiband_csv, save_path=out_b)
    assert a["figures"][0]["sha256"] == b["figures"][0]["sha256"]


# --------------------------------------------------------------------------- #
# Manifest wiring (AI-facing half), directly to lib/ (no skills/session.py)
# --------------------------------------------------------------------------- #

def test_tool_manifest_wired_to_lib_directly():
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.ixpe_pcube_plot"]
    assert entry["module"] == "xspec_utils.lib.ixpe_pcube_plot"
    assert entry["fn"] == "ixpe_pcube_plot"
    props = entry["args_schema"]["properties"]
    for name in ("csv_path", "save_path", "title"):
        assert name in props, f"{name} missing from args_schema.properties"
    assert entry["args_schema"]["required"] == ["csv_path"]
    assert "png" in entry["description"].lower()
