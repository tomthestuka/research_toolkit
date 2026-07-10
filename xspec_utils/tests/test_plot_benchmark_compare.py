"""Tests for xspec_utils' plot_benchmark_compare (lib/benchmark_compare.py).

`plot_benchmark_compare` is the companion plotter for `compare_to_benchmark`:
a forest/caterpillar plot of SIGNED z per record (direction-preserving,
unlike the CSV's own `abs(z)` 'z' column), so the reviewer can see at a
glance whether a benchmarked quantity's discrepancy from the paper is
systematic (every point on the same side of zero). It is pure matplotlib
(Agg backend) -- no XSPEC/ixpeobssim -- so the majority of this file
exercises it with hand-built records/rows, the same style
`test_benchmark_compare.py` already uses for `compare_to_benchmark` itself.

One test grounds the pipeline in real, already-local data: it extracts a
genuine PD/PA pair from Cyg X-1 IXPE ObsID 01002901's already-produced
`ixpe_pcube` CSV via `pcube_quantity` (the same real fixture
`test_benchmark_compare.py` uses), runs it through the real
`compare_to_benchmark` -> `plot_benchmark_compare` chain, and is skipped if
that fixture is not locally present.
"""

from __future__ import annotations

import math
import os
import shutil

import matplotlib.pyplot as plt
import pytest
from matplotlib.collections import PathCollection

from xspec_utils.lib import benchmark_compare as bc
from xspec_utils.skills import TOOL_MANIFEST

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
OBSID = "01002901"
REAL_PCUBE_CSV = os.path.join(REAL_IXPE_DIR, OBSID, f"ixpe_pcube_{OBSID}.csv")

RECORDS = [
    {"label": "PD_2-8keV", "value_this": 2.34, "err_this": 0.37,
     "value_paper": 2.5, "err_paper": 0.2, "unit": "percent",
     "paper_source": "Krawczynski+2022 Table S2"},
    {"label": "PA_2-8keV", "value_this": -14.79, "err_this": 4.53,
     "value_paper": -20.0, "err_paper": 2.0, "unit": "deg",
     "paper_source": "Krawczynski+2022 Table S2"},
    {"label": "PD_slope", "value_this": 0.614, "err_this": 0.128,
     "value_paper": 0.05, "err_paper": 0.02, "unit": "percent/keV",
     "paper_source": "Krawczynski+2022 Table S2"},
]


def _point_count(ax) -> int:
    return sum(pc.get_offsets().shape[0] for pc in ax.collections if isinstance(pc, PathCollection))


def _band_extent(ax):
    rect = ax.patches[0]
    return float(rect.get_x()), float(rect.get_x() + rect.get_width())


def _all_offsets_x(ax) -> list[float]:
    """Union the x-offsets of every PathCollection on *ax* (the plotter draws
    the unflagged and flagged groups as two separate `ax.scatter` calls, so a
    single collection is not guaranteed to carry every point)."""
    xs = []
    for pc in ax.collections:
        if isinstance(pc, PathCollection):
            xs.extend(float(v) for v in pc.get_offsets()[:, 0])
    return xs


# --------------------------------------------------------------------------- #
# Capability: N points plotted, correct flagging, band spans +/-sigma_flag,
# well-formed FigureRef.
# --------------------------------------------------------------------------- #

def test_capability_n_points_flagged_band_and_figureref(tmp_path):
    save_path = str(tmp_path / "compare.png")
    result = bc.plot_benchmark_compare(rows=RECORDS, save_path=save_path)

    assert result["exit_code"] == 0
    assert len(result["figures"]) == 1
    fig = result["figures"][0]
    assert fig["path"] == save_path
    assert os.path.exists(save_path)
    data = open(save_path, "rb").read()
    assert len(data) > 0
    import hashlib
    assert fig["sha256"] == hashlib.sha256(data).hexdigest()
    assert fig["mime_type"] == "image/png"
    assert fig["produced_by"] == "xspec_utils.lib.benchmark_compare.plot_benchmark_compare"
    assert fig["caption"]

    plotted = result["plotted"]
    assert plotted["n_records"] == len(RECORDS)
    assert plotted["label"] == [r["label"] for r in RECORDS]
    # PD_slope (z~4.35) is flagged at the default sigma_flag=3; the other two are not.
    assert plotted["flagged"] == [False, False, True]
    assert plotted["sigma_flag"] == pytest.approx(3.0)

    # The drawn artists themselves carry exactly N points, split flagged/unflagged.
    fig_obj, ax, plotted2 = bc._build_compare_figure(RECORDS, sigma_flag=3.0)
    try:
        assert _point_count(ax) == len(RECORDS)
        lo, hi = _band_extent(ax)
        assert lo == pytest.approx(-3.0)
        assert hi == pytest.approx(3.0)
        assert plotted2["flagged"] == plotted["flagged"]

        # Each record's marker sits at its own signed z -- read straight off the
        # drawn PathCollection offsets (not the returned `plotted` dict), so this
        # actually proves the figure places points at signed z, including the
        # negative one (PD_2-8keV, z~-0.38) landing left of zero.
        expected_z = [
            (r["value_this"] - r["value_paper"]) / math.sqrt(r["err_this"] ** 2 + r["err_paper"] ** 2)
            for r in RECORDS
        ]
        drawn_x = _all_offsets_x(ax)
        assert sorted(drawn_x) == pytest.approx(sorted(expected_z))
        assert any(x < 0 for x in drawn_x), "expected the negative PD_2-8keV z to be drawn left of zero"
    finally:
        plt.close(fig_obj)


def test_signed_z_matches_independent_recompute():
    for rec in RECORDS:
        expected = (rec["value_this"] - rec["value_paper"]) / math.sqrt(
            rec["err_this"] ** 2 + rec["err_paper"] ** 2)
        row = {"label": rec["label"], "err_this": rec["err_this"], "err_paper": rec["err_paper"],
               "delta": rec["value_this"] - rec["value_paper"]}
        assert bc._row_signed_z(row) == pytest.approx(expected)


# --------------------------------------------------------------------------- #
# Edge cases
# --------------------------------------------------------------------------- #

def test_empty_rows_returns_exit_code_1_no_png(tmp_path):
    save_path = str(tmp_path / "empty.png")
    result = bc.plot_benchmark_compare(rows=[], save_path=save_path)
    assert result["exit_code"] == 1
    assert result["stderr"]
    assert result["figures"] == []
    assert not os.path.exists(save_path)


def test_no_input_raises():
    with pytest.raises(ValueError, match="supply rows or csv_path"):
        bc.plot_benchmark_compare()


def test_infinite_z_clamped_with_annotation(tmp_path):
    rows = [
        {"label": "PD_zero_zero", "err_this": 0.0, "err_paper": 0.0, "delta": 4.0},
        {"label": "PA_normal", "err_this": 3.0, "err_paper": 2.0, "delta": 5.0},
    ]
    fig, ax, plotted = bc._build_compare_figure(rows, sigma_flag=3.0)
    try:
        assert math.isinf(plotted["z"][0]) and plotted["z"][0] > 0
        assert plotted["flagged"][0] is True
        xlim = ax.get_xlim()
        assert all(math.isfinite(v) for v in xlim)
        annotations = [t.get_text() for t in ax.texts]
        assert any("PD_zero_zero" in a and "z=inf" in a for a in annotations)
    finally:
        plt.close(fig)

    # End-to-end: the full call does not raise and still writes a PNG.
    result = bc.plot_benchmark_compare(rows=rows, save_path=str(tmp_path / "inf.png"))
    assert result["exit_code"] == 0
    assert os.path.getsize(result["figures"][0]["path"]) > 0


def test_csv_path_not_found_raises():
    with pytest.raises(ValueError, match="not found"):
        bc.plot_benchmark_compare(csv_path="/nonexistent/path/compare.csv")


def test_csv_path_missing_provenance_header_raises(tmp_path):
    bad_csv = tmp_path / "bad.csv"
    bad_csv.write_text("label,unit,value_this\nPD,percent,1.0\n")
    with pytest.raises(ValueError, match="provenance"):
        bc.plot_benchmark_compare(csv_path=str(bad_csv))


def test_row_missing_required_field_raises():
    bad_rows = [{"label": "PD", "err_this": 0.1}]  # missing err_paper and delta/value_this/value_paper
    with pytest.raises(ValueError, match="missing required field"):
        bc.plot_benchmark_compare(rows=bad_rows)


# --------------------------------------------------------------------------- #
# Differential: rows list vs csv_path (compare_to_benchmark's own product)
# --------------------------------------------------------------------------- #

def test_differential_rows_vs_csv_path_agree(tmp_path):
    out_csv = str(tmp_path / "compare.csv")
    compare_result = bc.compare_to_benchmark(RECORDS, out_csv)

    result_rows = bc.plot_benchmark_compare(rows=compare_result["rows"],
                                             save_path=str(tmp_path / "from_rows.png"))
    result_csv = bc.plot_benchmark_compare(csv_path=compare_result["out_csv"],
                                            save_path=str(tmp_path / "from_csv.png"))

    assert result_rows["exit_code"] == 0
    assert result_csv["exit_code"] == 0
    assert result_rows["plotted"]["n_records"] == result_csv["plotted"]["n_records"]

    flagged_labels_rows = {lbl for lbl, f in zip(result_rows["plotted"]["label"],
                                                  result_rows["plotted"]["flagged"]) if f}
    flagged_labels_csv = {lbl for lbl, f in zip(result_csv["plotted"]["label"],
                                                 result_csv["plotted"]["flagged"]) if f}
    assert flagged_labels_rows == flagged_labels_csv == {"PD_slope"}

    for z_rows, z_csv in zip(result_rows["plotted"]["z"], result_csv["plotted"]["z"]):
        assert z_rows == pytest.approx(z_csv)


def test_sigma_flag_read_from_csv_provenance_when_unset(tmp_path):
    out_csv = str(tmp_path / "compare_sf25.csv")
    bc.compare_to_benchmark(RECORDS, out_csv, sigma_flag=2.5)

    result = bc.plot_benchmark_compare(csv_path=out_csv, save_path=str(tmp_path / "sf25.png"))
    assert result["plotted"]["sigma_flag"] == pytest.approx(2.5)
    # At sigma_flag=2.5, PA_2-8keV (z~1.05) still unflagged, PD_slope (z~4.35) flagged.
    assert result["plotted"]["flagged"] == [False, False, True]


def test_sigma_flag_explicit_override_wins_over_provenance(tmp_path):
    out_csv = str(tmp_path / "compare_sf25b.csv")
    bc.compare_to_benchmark(RECORDS, out_csv, sigma_flag=2.5)

    result = bc.plot_benchmark_compare(csv_path=out_csv, sigma_flag=10.0,
                                        save_path=str(tmp_path / "sf10.png"))
    assert result["plotted"]["sigma_flag"] == pytest.approx(10.0)
    assert result["plotted"]["flagged"] == [False, False, False]


# --------------------------------------------------------------------------- #
# Reproducibility
# --------------------------------------------------------------------------- #

def test_repro_byte_identical(tmp_path):
    out_a = str(tmp_path / "a.png")
    out_b = str(tmp_path / "b.png")
    a = bc.plot_benchmark_compare(rows=RECORDS, save_path=out_a)
    b = bc.plot_benchmark_compare(rows=RECORDS, save_path=out_b)
    assert a["figures"][0]["sha256"] == b["figures"][0]["sha256"]


# --------------------------------------------------------------------------- #
# Manifest wiring (AI-facing half), directly to lib/ (no skills/session.py)
# --------------------------------------------------------------------------- #

def test_tool_manifest_wired_to_lib_directly():
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.plot_benchmark_compare"]
    assert entry["module"] == "xspec_utils.lib.benchmark_compare"
    assert entry["fn"] == "plot_benchmark_compare"
    props = entry["args_schema"]["properties"]
    for name in ("rows", "csv_path", "save_path", "title", "sigma_flag"):
        assert name in props, f"{name} missing from args_schema.properties"
    assert entry["args_schema"]["required"] == []
    assert "png" in entry["description"].lower()
    assert "sigma_flag" in entry["description"]


# --------------------------------------------------------------------------- #
# End-to-end grounded in real, already-local IXPE data (ObsID 01002901)
# --------------------------------------------------------------------------- #

def test_end_to_end_from_real_pcube_quantity(tmp_path):
    if not os.path.exists(REAL_PCUBE_CSV):
        pytest.skip(f"real fixture ixpe_pcube CSV not local at {REAL_PCUBE_CSV}")
    csv_path = str(tmp_path / os.path.basename(REAL_PCUBE_CSV))
    shutil.copyfile(REAL_PCUBE_CSV, csv_path)

    pd_this, pd_this_err = bc.pcube_quantity(csv_path, 0, "PD", unit="percent")
    pa_this, pa_this_err = bc.pcube_quantity(csv_path, 0, "PA")

    records = [
        {"label": "PD_band0", "value_this": pd_this, "err_this": pd_this_err,
         "value_paper": 2.5, "err_paper": 0.2, "unit": "percent",
         "paper_source": "Krawczynski+2022 Table S2"},
        {"label": "PA_band0", "value_this": pa_this, "err_this": pa_this_err,
         "value_paper": -20.0, "err_paper": 2.0, "unit": "deg",
         "paper_source": "Krawczynski+2022 Table S2"},
    ]
    compare_result = bc.compare_to_benchmark(records, str(tmp_path / "real_compare.csv"))

    plot_rows = bc.plot_benchmark_compare(rows=compare_result["rows"],
                                           save_path=str(tmp_path / "real_rows.png"))
    plot_csv = bc.plot_benchmark_compare(csv_path=compare_result["out_csv"],
                                          save_path=str(tmp_path / "real_csv.png"))

    assert plot_rows["exit_code"] == 0
    assert plot_csv["exit_code"] == 0
    assert plot_rows["plotted"]["label"] == plot_csv["plotted"]["label"] == ["PD_band0", "PA_band0"]
    for z_rows, z_csv in zip(plot_rows["plotted"]["z"], plot_csv["plotted"]["z"]):
        assert z_rows == pytest.approx(z_csv)
