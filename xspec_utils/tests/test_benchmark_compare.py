"""Tests for xspec_utils' compare_to_benchmark (lib/benchmark_compare.py) and
its two extractor helpers pcube_quantity/fit_quantity.

`compare_to_benchmark` itself is pure stdlib (no XSPEC/ixpeobssim) -- the
fast majority of this file exercises it with hand-built records, since its
job is the delta/z-score/flag arithmetic and CSV I/O, not any new science.

`pcube_quantity` is checked against the real, already-local `ixpe_pcube`
CSV for Cyg X-1 ObsID 01002901 (`ixpe_pcube_01002901.csv`, produced by a
prior real `ixpe_pcube` reduction) -- copied read-only into tmp_path so this
test never touches the shared analyst-env data directory.

`fit_quantity` is checked against a genuine `fit_tbabs_polconst_powerlaw`
result: one real PyXSPEC polconst spectropolarimetric fit against the
already-local IXPE Level-2 event files for the same ObsID (mirroring
`test_fit_tbabs_polconst_powerlaw.py`'s own fixture pattern), so the PD/PA/
PD_err/PA_err values `fit_quantity` converts are real fitted numbers, not
invented ones.
"""

from __future__ import annotations

import csv
import math
import os
import shutil

import pytest

from xspec_utils.lib import benchmark_compare as bc
from xspec_utils.skills import TOOL_MANIFEST

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
OBSID = "01002901"
REAL_PCUBE_CSV = os.path.join(REAL_IXPE_DIR, OBSID, f"ixpe_pcube_{OBSID}.csv")


def _copy_real_pcube_csv(dest_dir: str) -> str:
    if not os.path.exists(REAL_PCUBE_CSV):
        pytest.skip(f"real fixture ixpe_pcube CSV not local at {REAL_PCUBE_CSV}")
    dst = os.path.join(dest_dir, os.path.basename(REAL_PCUBE_CSV))
    shutil.copyfile(REAL_PCUBE_CSV, dst)
    return dst


def _read_real_band0() -> dict:
    """Independent read of the real CSV's first band row, for comparison
    against pcube_quantity's return values."""
    with open(REAL_PCUBE_CSV, newline='') as f:
        f.readline()  # provenance header
        row = next(csv.DictReader(f))
    return {k: float(v) for k, v in row.items() if k != 'band_id'}


# --------------------------------------------------------------------------- #
# compare_to_benchmark: input validation / edge cases
# --------------------------------------------------------------------------- #

def test_empty_records_raises(tmp_path):
    with pytest.raises(ValueError, match="empty records"):
        bc.compare_to_benchmark([], str(tmp_path / "out.csv"))


def test_empty_paper_source_raises(tmp_path):
    records = [{"label": "PD", "value_this": 1.0, "err_this": 0.1,
                "value_paper": 1.0, "err_paper": 0.1, "unit": "percent",
                "paper_source": ""}]
    with pytest.raises(ValueError, match="paper_source"):
        bc.compare_to_benchmark(records, str(tmp_path / "out.csv"))


def test_missing_paper_source_raises(tmp_path):
    records = [{"label": "PD", "value_this": 1.0, "err_this": 0.1,
                "value_paper": 1.0, "err_paper": 0.1, "unit": "percent"}]
    with pytest.raises(ValueError, match="paper_source"):
        bc.compare_to_benchmark(records, str(tmp_path / "out.csv"))


def test_mismatched_unit_this_paper_raises(tmp_path):
    records = [{"label": "PD", "value_this": 1.0, "err_this": 0.1,
                "value_paper": 1.0, "err_paper": 0.1,
                "unit_this": "keV", "unit_paper": "eV",
                "paper_source": "Krawczynski+2022"}]
    with pytest.raises(ValueError, match="mismatched unit"):
        bc.compare_to_benchmark(records, str(tmp_path / "out.csv"))


def test_unit_this_without_unit_paper_raises(tmp_path):
    records = [{"label": "PD", "value_this": 1.0, "err_this": 0.1,
                "value_paper": 1.0, "err_paper": 0.1,
                "unit_this": "percent",
                "paper_source": "Krawczynski+2022"}]
    with pytest.raises(ValueError, match="unit_this/unit_paper"):
        bc.compare_to_benchmark(records, str(tmp_path / "out.csv"))


def test_unit_and_pair_must_agree(tmp_path):
    records = [{"label": "PD", "value_this": 1.0, "err_this": 0.1,
                "value_paper": 1.0, "err_paper": 0.1,
                "unit": "percent", "unit_this": "percent", "unit_paper": "fraction",
                "paper_source": "Krawczynski+2022"}]
    with pytest.raises(ValueError, match="mismatched unit"):
        bc.compare_to_benchmark(records, str(tmp_path / "out.csv"))


def test_missing_unit_entirely_raises(tmp_path):
    records = [{"label": "PD", "value_this": 1.0, "err_this": 0.1,
                "value_paper": 1.0, "err_paper": 0.1,
                "paper_source": "Krawczynski+2022"}]
    with pytest.raises(ValueError, match="unit"):
        bc.compare_to_benchmark(records, str(tmp_path / "out.csv"))


def test_negative_error_raises(tmp_path):
    records = [{"label": "PD", "value_this": 1.0, "err_this": -0.1,
                "value_paper": 1.0, "err_paper": 0.1, "unit": "percent",
                "paper_source": "Krawczynski+2022"}]
    with pytest.raises(ValueError, match="err_this"):
        bc.compare_to_benchmark(records, str(tmp_path / "out.csv"))


def test_non_numeric_value_raises(tmp_path):
    records = [{"label": "PD", "value_this": "oops", "err_this": 0.1,
                "value_paper": 1.0, "err_paper": 0.1, "unit": "percent",
                "paper_source": "Krawczynski+2022"}]
    with pytest.raises(ValueError, match="value_this"):
        bc.compare_to_benchmark(records, str(tmp_path / "out.csv"))


def test_non_finite_value_raises(tmp_path):
    nan_records = [{"label": "PD", "value_this": float("nan"), "err_this": 0.1,
                     "value_paper": 1.0, "err_paper": 0.1, "unit": "percent",
                     "paper_source": "Krawczynski+2022"}]
    with pytest.raises(ValueError, match="finite"):
        bc.compare_to_benchmark(nan_records, str(tmp_path / "out_nan.csv"))

    inf_err_records = [{"label": "PD", "value_this": 1.0, "err_this": 0.1,
                         "value_paper": 1.0, "err_paper": float("inf"), "unit": "percent",
                         "paper_source": "Krawczynski+2022"}]
    with pytest.raises(ValueError, match="finite"):
        bc.compare_to_benchmark(inf_err_records, str(tmp_path / "out_inf.csv"))


def test_zero_zero_errors_gives_inf_and_flagged(tmp_path):
    records = [{"label": "PD", "value_this": 5.0, "err_this": 0.0,
                "value_paper": 1.0, "err_paper": 0.0, "unit": "percent",
                "paper_source": "Krawczynski+2022"}]
    result = bc.compare_to_benchmark(records, str(tmp_path / "out.csv"))
    row = result["rows"][0]
    assert math.isinf(row["z"]) and row["z"] > 0
    assert row["flagged"] is True
    assert result["flags"] == ["PD"]


# --------------------------------------------------------------------------- #
# Capability: the tool actually computed z, not shaped metadata
# --------------------------------------------------------------------------- #

def test_capability_zero_and_five_sigma_flagging(tmp_path):
    records = [
        {"label": "zero_sigma", "value_this": 2.0, "err_this": 0.5,
         "value_paper": 2.0, "err_paper": 0.5, "unit": "percent",
         "paper_source": "Krawczynski+2022 Table S2"},
        {"label": "five_sigma", "value_this": 5.0, "err_this": 1.0,
         "value_paper": 0.0, "err_paper": 0.0, "unit": "percent",
         "paper_source": "Krawczynski+2022 Table S2"},
        {"label": "neg_delta", "value_this": 1.0, "err_this": 1.0,
         "value_paper": 3.0, "err_paper": 0.0, "unit": "percent",
         "paper_source": "Krawczynski+2022 Table S2"},
    ]
    result = bc.compare_to_benchmark(records, str(tmp_path / "out.csv"))
    assert result["flags"] == ["five_sigma"]

    by_label = {r["label"]: r for r in result["rows"]}
    assert by_label["zero_sigma"]["z"] == pytest.approx(0.0)
    assert by_label["zero_sigma"]["flagged"] is False

    for rec, row in zip(records, result["rows"]):
        independent_z = abs(rec["value_this"] - rec["value_paper"]) / math.sqrt(
            rec["err_this"] ** 2 + rec["err_paper"] ** 2)
        assert row["z"] == pytest.approx(independent_z)
    assert by_label["five_sigma"]["z"] == pytest.approx(5.0)

    # Signed delta = value_this - value_paper -- not reversed, not abs().
    assert by_label["five_sigma"]["delta"] == pytest.approx(5.0)   # 5.0 - 0.0
    assert by_label["neg_delta"]["delta"] == pytest.approx(-2.0)  # 1.0 - 3.0, independently
    assert by_label["neg_delta"]["flagged"] is False


# --------------------------------------------------------------------------- #
# Differential: round-trip the written CSV
# --------------------------------------------------------------------------- #

def test_csv_round_trip_and_provenance_header(tmp_path):
    out_csv = str(tmp_path / "compare.csv")
    records = [
        {"label": "PD_2-8keV", "value_this": 2.34, "err_this": 0.37,
         "value_paper": 2.5, "err_paper": 0.2, "unit": "percent",
         "paper_source": "Krawczynski+2022 Table S2"},
        {"label": "PA_2-8keV", "value_this": -14.79, "err_this": 4.53,
         "value_paper": -20.0, "err_paper": 2.0, "unit": "deg",
         "paper_source": "Krawczynski+2022 Table S2"},
    ]
    sigma_flag = 2.5
    result = bc.compare_to_benchmark(records, out_csv, sigma_flag=sigma_flag)

    with open(result["out_csv"], newline="") as f:
        header_line = f.readline()
        assert header_line.startswith("# compare_to_benchmark provenance: ")
        assert f'"sigma_flag": {sigma_flag}' in header_line
        assert '"created_at":' in header_line
        assert '"produced_by": "xspec_utils.lib.benchmark_compare.compare_to_benchmark"' in header_line
        read_rows = list(csv.DictReader(f))

    assert len(read_rows) == len(result["rows"])
    for got, expected in zip(read_rows, result["rows"]):
        assert got["label"] == expected["label"]
        assert got["unit"] == expected["unit"]
        assert got["paper_source"] == expected["paper_source"]
        assert got["flagged"] == str(expected["flagged"])
        for key in ("value_this", "err_this", "value_paper", "err_paper", "delta", "z"):
            assert float(got[key]) == pytest.approx(expected[key])

    # Signed delta pinned on the round-tripped CSV, computed independently
    # here from the record's own inputs (value_this - value_paper, signed).
    pd_row = next(r for r in read_rows if r["label"] == "PD_2-8keV")
    assert float(pd_row["delta"]) == pytest.approx(2.34 - 2.5)
    assert float(pd_row["delta"]) == pytest.approx(-0.16)


def test_creates_parent_dirs(tmp_path):
    out_csv = str(tmp_path / "nested" / "dir" / "out.csv")
    records = [{"label": "PD", "value_this": 1.0, "err_this": 0.1,
                "value_paper": 1.0, "err_paper": 0.1, "unit": "percent",
                "paper_source": "Krawczynski+2022"}]
    result = bc.compare_to_benchmark(records, out_csv)
    assert os.path.exists(result["out_csv"])


# --------------------------------------------------------------------------- #
# Reproducibility
# --------------------------------------------------------------------------- #

def test_repro_compare_to_benchmark(tmp_path):
    records = [
        {"label": "PD_2-8keV", "value_this": 2.34, "err_this": 0.37,
         "value_paper": 2.5, "err_paper": 0.2, "unit": "percent",
         "paper_source": "Krawczynski+2022 Table S2"},
        {"label": "PA_2-8keV", "value_this": -14.79, "err_this": 4.53,
         "value_paper": -20.0, "err_paper": 2.0, "unit": "deg",
         "paper_source": "Krawczynski+2022 Table S2"},
        {"label": "zero_zero", "value_this": 5.0, "err_this": 0.0,
         "value_paper": 1.0, "err_paper": 0.0, "unit": "percent",
         "paper_source": "Krawczynski+2022 Table S2"},
    ]
    result_a = bc.compare_to_benchmark(records, str(tmp_path / "a.csv"))
    result_b = bc.compare_to_benchmark(records, str(tmp_path / "b.csv"))

    assert result_a["flags"] == result_b["flags"]
    assert len(result_a["rows"]) == len(result_b["rows"])
    for row_a, row_b in zip(result_a["rows"], result_b["rows"]):
        assert row_a["label"] == row_b["label"]
        assert row_a["unit"] == row_b["unit"]
        assert row_a["paper_source"] == row_b["paper_source"]
        assert row_a["flagged"] == row_b["flagged"]
        for key in ("value_this", "err_this", "value_paper", "err_paper", "delta"):
            assert row_a[key] == row_b[key]
        if math.isinf(row_a["z"]) or math.isinf(row_b["z"]):
            assert math.isinf(row_a["z"]) and math.isinf(row_b["z"])
        else:
            assert row_a["z"] == row_b["z"]


# --------------------------------------------------------------------------- #
# Manifest wiring (AI-facing half)
# --------------------------------------------------------------------------- #

def test_tool_manifest_wiring():
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.compare_to_benchmark"]
    assert entry["module"] == "xspec_utils.lib.benchmark_compare"
    assert entry["fn"] == "compare_to_benchmark"
    props = entry["args_schema"]["properties"]
    assert set(entry["args_schema"]["required"]) == {"records", "out_csv"}
    rec_props = props["records"]["items"]["properties"]
    for name in ("label", "value_this", "err_this", "value_paper", "err_paper",
                 "paper_source", "unit", "unit_this", "unit_paper"):
        assert name in rec_props
    assert "sigma_flag" in props
    assert "paper_source" in entry["description"]
    # The extractor helpers keep unit handling out of the caller -- they stay
    # lib-only, never their own AI-facing manifest entries.
    assert "research_toolkit.xspec_utils.pcube_quantity" not in TOOL_MANIFEST
    assert "research_toolkit.xspec_utils.fit_quantity" not in TOOL_MANIFEST


# --------------------------------------------------------------------------- #
# pcube_quantity -- real ixpe_pcube CSV for ObsID 01002901
# --------------------------------------------------------------------------- #

def test_pcube_quantity_pd_percent_is_100x_fraction(tmp_path):
    csv_path = _copy_real_pcube_csv(str(tmp_path))
    real = _read_real_band0()

    value_frac, err_frac = bc.pcube_quantity(csv_path, 0, "PD", unit="fraction")
    value_pct, err_pct = bc.pcube_quantity(csv_path, 0, "PD", unit="percent")

    assert value_frac == pytest.approx(real["PD"])
    assert err_frac == pytest.approx(real["PD_err"])
    assert value_pct == pytest.approx(real["PD"] * 100.0)
    assert err_pct == pytest.approx(real["PD_err"] * 100.0)
    assert value_pct == pytest.approx(value_frac * 100.0)
    assert err_pct == pytest.approx(err_frac * 100.0)


def test_pcube_quantity_qn_un_percent_is_100x_fraction(tmp_path):
    csv_path = _copy_real_pcube_csv(str(tmp_path))
    real = _read_real_band0()
    for name in ("QN", "UN"):
        value_frac, err_frac = bc.pcube_quantity(csv_path, 0, name, unit="fraction")
        value_pct, err_pct = bc.pcube_quantity(csv_path, 0, name, unit="percent")
        assert value_frac == pytest.approx(real[name])
        assert value_pct == pytest.approx(real[name] * 100.0)
        assert err_pct == pytest.approx(err_frac * 100.0)


def test_pcube_quantity_pa_unchanged_across_units(tmp_path):
    csv_path = _copy_real_pcube_csv(str(tmp_path))
    real = _read_real_band0()

    value_frac, err_frac = bc.pcube_quantity(csv_path, 0, "PA", unit="fraction")
    value_pct, err_pct = bc.pcube_quantity(csv_path, 0, "PA", unit="percent")

    assert value_frac == pytest.approx(real["PA"])
    assert err_frac == pytest.approx(real["PA_err"])
    assert value_pct == pytest.approx(value_frac)
    assert err_pct == pytest.approx(err_frac)


def test_pcube_quantity_band_index_out_of_range_raises(tmp_path):
    csv_path = _copy_real_pcube_csv(str(tmp_path))
    with pytest.raises(ValueError, match="band_index"):
        bc.pcube_quantity(csv_path, 999, "PD")


def test_pcube_quantity_unsupported_name_raises(tmp_path):
    csv_path = _copy_real_pcube_csv(str(tmp_path))
    with pytest.raises(ValueError, match="unsupported name"):
        bc.pcube_quantity(csv_path, 0, "SIGNIF")


def test_pcube_quantity_missing_header_raises(tmp_path):
    csv_path = str(tmp_path / "no_header.csv")
    with open(csv_path, "w") as f:
        f.write("band_id,E_LO,E_HI\n0,2.0,3.0\n")
    with pytest.raises(ValueError, match="provenance"):
        bc.pcube_quantity(csv_path, 0, "PD")


def test_pcube_quantity_matches_red_ixpe_reader(tmp_path):
    """Differential: pcube_quantity must delegate to red_ixpe.read_pcube_csv
    (the single source of truth for the ixpe_pcube CSV contract), not a local
    reimplementation -- checked by independently calling read_pcube_csv here
    and asserting pcube_quantity returns exactly the same band-0 values."""
    pytest.importorskip("ixpeobssim")
    csv_path = _copy_real_pcube_csv(str(tmp_path))

    from xspec_utils.lib.red_ixpe import read_pcube_csv
    _, bands = read_pcube_csv(csv_path)

    assert bc.pcube_quantity(csv_path, 0, "PD", unit="fraction") == (bands[0]["PD"], bands[0]["PD_err"])
    assert bc.pcube_quantity(csv_path, 0, "PA", unit="fraction") == (bands[0]["PA"], bands[0]["PA_err"])
    assert bc.pcube_quantity(csv_path, 0, "PD", unit="percent") == (100.0 * bands[0]["PD"], 100.0 * bands[0]["PD_err"])


# --------------------------------------------------------------------------- #
# fit_quantity -- real fit_tbabs_polconst_powerlaw result for ObsID 01002901
# --------------------------------------------------------------------------- #

VERSION = "v01"
SRC_RAD, BKG_INNER_RAD, BKG_OUTER_RAD = 1.333, 2.5, 5.167


def _copy_event_files(dest_dir: str) -> None:
    src_evt_dir = os.path.join(REAL_IXPE_DIR, OBSID, "event_l2")
    dst_evt_dir = os.path.join(dest_dir, OBSID, "event_l2")
    os.makedirs(dst_evt_dir, exist_ok=True)
    for det in (1, 2, 3):
        src = os.path.join(src_evt_dir, f"ixpe{OBSID}_det{det}_evt2_{VERSION}_rej.fits")
        if not os.path.exists(src):
            pytest.skip(f"real fixture ObsID {OBSID} not local at {src}")
        shutil.copyfile(src, os.path.join(dst_evt_dir, os.path.basename(src)))


@pytest.fixture(scope="module")
def real_fit_dict(tmp_path_factory):
    """One genuine PyXSPEC polconst spectropolarimetric fit against the
    already-local IXPE Level-2 event files for ObsID 01002901 -- the actual
    fit_tbabs_polconst_powerlaw result dict fit_quantity is meant to consume,
    not an invented one."""
    xspec = pytest.importorskip("xspec")
    from xspec_utils.lib import xspec_help

    ixpe_dir = str(tmp_path_factory.mktemp("benchmark_compare_ixpe"))
    _copy_event_files(ixpe_dir)
    xspec.AllData.clear()
    xspec.AllModels.clear()
    indices = xspec_help.load_data(
        epoch={"IXPE": [OBSID]},
        ixpe_path=ixpe_dir,
        ixpe_changes=True,
        ixpe_src_rad=SRC_RAD,
        ixpe_bkg_inner_rad=BKG_INNER_RAD,
        ixpe_bkg_outer_rad=BKG_OUTER_RAD,
    )
    return xspec_help.fit_tbabs_polconst_powerlaw(indices)


def test_fit_quantity_pd_percent_is_100x_fraction(real_fit_dict):
    value_frac, err_frac = bc.fit_quantity(real_fit_dict, "PD", unit="fraction")
    value_pct, err_pct = bc.fit_quantity(real_fit_dict, "PD", unit="percent")
    assert value_frac == pytest.approx(real_fit_dict["PD"])
    assert err_frac == pytest.approx(real_fit_dict["PD_err"])
    assert value_pct == pytest.approx(value_frac * 100.0)
    assert err_pct == pytest.approx(err_frac * 100.0)


def test_fit_quantity_pa_unchanged_across_units(real_fit_dict):
    value_frac, err_frac = bc.fit_quantity(real_fit_dict, "PA", unit="fraction")
    value_pct, err_pct = bc.fit_quantity(real_fit_dict, "PA", unit="percent")
    assert value_frac == pytest.approx(real_fit_dict["PA"])
    assert err_frac == pytest.approx(real_fit_dict["PA_err"])
    assert value_pct == pytest.approx(value_frac)
    assert err_pct == pytest.approx(err_frac)


def test_fit_quantity_missing_key_raises(real_fit_dict):
    incomplete = {k: v for k, v in real_fit_dict.items() if k != "PD_err"}
    with pytest.raises(ValueError, match="PD_err"):
        bc.fit_quantity(incomplete, "PD")


def test_fit_quantity_unsupported_name_raises(real_fit_dict):
    with pytest.raises(ValueError, match="unsupported name"):
        bc.fit_quantity(real_fit_dict, "Gamma")
