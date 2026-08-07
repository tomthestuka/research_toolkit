"""Tests for the `tau_axis` argument added to `qu_lag_dcf_plot`
(`lib/qu_lag_dcf.py`).

The reviewer must be able to render the |tau_peak| population panels that
correspond to the abs_tau_peak_* Spearman statistic `ixpe_qu_lag_dcf`
already reports, so the figure and the quoted number never disagree. This
adds `tau_axis: str = "signed"` (`"signed"`/`"abs"`) to `qu_lag_dcf_plot`.
The `dcf_panels` grid is unaffected by `tau_axis` (lag is signed there by
definition); only the two population panels respond.

Uses the same already-local mixed-sign fixture set as
`test_qu_lag_dcf_population_abs.py` (01002901, 02250301,
synth_qu_lag_evensamp01, synth_qu_lag_null01 -- confirmed there to recover
both a negative and a positive peak lag), copied into an isolated tmp
obs_db_path so this suite's own writes never touch the shared analyst-env
mission directory. The differential interval-mapping test instead uses a
minimal hand-built ccf_master.csv/qu_dcf.csv/population_correlations.csv
fixture (mirroring `test_qu_lag_dcf_plot.py`'s own `_write_min_*` fixtures)
to exercise the bracket-zero and non-bracket bootstrap-interval cases under
full control, independent of what the real fixtures happen to recover.
"""

from __future__ import annotations

import csv
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
OBSID_A = "01002901"
OBSID_B = "02250301"
OBSID_EVENSAMP = "synth_qu_lag_evensamp01"
OBSID_NULL = "synth_qu_lag_null01"
MIXED_SIGN_OBSIDS = [OBSID_A, OBSID_B, OBSID_EVENSAMP, OBSID_NULL]

POPULATION_KEYS = ("population_tau_r_proxy", "population_tau_pd_signif")
COVARIATE_OF = {"population_tau_r_proxy": "r_proxy", "population_tau_pd_signif": "pd_significance"}
SIGNED_PAIR_OF = {
    "population_tau_r_proxy": "tau_peak_vs_r_proxy",
    "population_tau_pd_signif": "tau_peak_vs_pd_significance",
}
ABS_PAIR_OF = {
    "population_tau_r_proxy": "abs_tau_peak_vs_r_proxy",
    "population_tau_pd_signif": "abs_tau_peak_vs_pd_significance",
}


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
def real_ccf_mixed(tmp_path_factory):
    """One genuine `qu_lag_dcf` run over the four mixed-sign fixture ObsIDs:
    a real `ccf_master.csv` (with both a negative and a positive
    `tau_peak_s`), one `qu_dcf_<subobs_key>.csv` per obsid, and a real
    `population_correlations.csv` carrying all four population pairs."""
    obs_db = str(tmp_path_factory.mktemp("qu_lag_dcf_plot_tau_axis_fixture"))
    for obsid in MIXED_SIGN_OBSIDS:
        _copy_qu_slices(obs_db, obsid)
    out_dir = os.path.join(obs_db, "out")
    res = m.qu_lag_dcf(obs_db, MIXED_SIGN_OBSIDS, out_dir, n_permutations=50, n_bootstrap=50, seed=0)
    assert res["exit_code"] == 0
    assert all(res["per_subobs"][o]["status"] == "ok" for o in MIXED_SIGN_OBSIDS)
    with open(res["ccf_master_csv"], newline="") as f:
        taus = [float(r["tau_peak_s"]) for r in csv.DictReader(f) if r["status"] == "ok"]
    assert any(t < 0 for t in taus) and any(t > 0 for t in taus), (
        "fixture must carry mixed-sign tau_peak_s or the abs-vs-signed capability "
        "check below is vacuous")
    return {
        "obs_db": obs_db, "out_dir": out_dir,
        "ccf_master_csv": res["ccf_master_csv"], "population_csv": res["population_csv"],
    }


# --------------------------------------------------------------------------- #
# Capability: tau_axis='abs' population panels equal abs(tau_peak_s) from the
# CSV element-by-element (all >= 0), differ from the signed values, and their
# annotation matches the abs_tau_peak_* rows, not the signed ones. Every PNG
# exists/non-empty, and the _abs filenames differ from the default ones.
# --------------------------------------------------------------------------- #

def test_capability_abs_population_matches_csv_and_differs_from_signed(tmp_path, real_ccf_mixed):
    result = m.qu_lag_dcf_plot(real_ccf_mixed["obs_db"], real_ccf_mixed["ccf_master_csv"],
                                str(tmp_path / "figs_abs"),
                                save_path_prefix=str(tmp_path / "figs_abs" / "plot"),
                                tau_axis="abs")
    assert result["exit_code"] == 0
    assert len(result["figures"]) == 3
    for fig in result["figures"]:
        assert os.path.exists(fig["path"])
        assert os.path.getsize(fig["path"]) > 0
        for req_key in ("path", "caption", "mime_type", "produced_by", "created_at", "sha256"):
            assert req_key in fig, f"FigureRef missing {req_key!r}: {fig}"

    # Caption contract: the two abs population figures name |tau_peak|; dcf_panels
    # never does (unaffected by tau_axis).
    abs_pop_captions = [f["caption"] for f in result["figures"]
                         if "population" in os.path.basename(f["path"])]
    assert len(abs_pop_captions) == 2
    for caption in abs_pop_captions:
        assert "|tau_peak|" in caption

    with open(real_ccf_mixed["ccf_master_csv"], newline="") as f:
        ccf_rows = {r["subobs_key"]: r for r in csv.DictReader(f) if r["status"] == "ok"}
    with open(real_ccf_mixed["population_csv"], newline="") as f:
        pop_rows = {r["pair"]: r for r in csv.DictReader(f)}

    for key in POPULATION_KEYS:
        panel = result["plotted"][key]
        assert set(panel["subobs_key"]) == set(MIXED_SIGN_OBSIDS)
        for i, subobs_key in enumerate(panel["subobs_key"]):
            signed_tau = float(ccf_rows[subobs_key]["tau_peak_s"])
            assert panel["y"][i] == pytest.approx(abs(signed_tau), rel=1e-9, abs=1e-12)
            assert panel["y"][i] >= 0.0
            if signed_tau < 0.0:
                assert panel["y"][i] != pytest.approx(signed_tau, abs=1e-9)

        abs_pair, signed_pair = ABS_PAIR_OF[key], SIGNED_PAIR_OF[key]
        exp_abs, exp_signed = pop_rows[abs_pair], pop_rows[signed_pair]
        assert panel["spearman"]["n"] == int(float(exp_abs["n"]))
        exp_abs_rho = float(exp_abs["rho"])
        if math.isnan(exp_abs_rho):
            assert math.isnan(panel["spearman"]["rho"])
        else:
            assert panel["spearman"]["rho"] == pytest.approx(exp_abs_rho, rel=1e-9, abs=1e-12)
            assert panel["spearman"]["p"] == pytest.approx(float(exp_abs["p"]), rel=1e-9, abs=1e-12)
            # Genuinely the abs annotation, not the signed one, for this mixed-sign fixture.
            assert exp_abs_rho != pytest.approx(float(exp_signed["rho"]), abs=1e-9)

    abs_fig_names = {os.path.basename(f["path"]) for f in result["figures"]}
    assert any(n.endswith("_population_tau_r_proxy_abs.png") for n in abs_fig_names)
    assert any(n.endswith("_population_tau_pd_signif_abs.png") for n in abs_fig_names)
    assert any(n.endswith("_dcf_panels.png") for n in abs_fig_names)

    result_signed = m.qu_lag_dcf_plot(real_ccf_mixed["obs_db"], real_ccf_mixed["ccf_master_csv"],
                                       str(tmp_path / "figs_signed"),
                                       save_path_prefix=str(tmp_path / "figs_signed" / "plot"),
                                       tau_axis="signed")
    signed_fig_names = {os.path.basename(f["path"]) for f in result_signed["figures"]}
    assert signed_fig_names != abs_fig_names
    # The two population filenames never collide between axes (dcf_panels is
    # legitimately shared -- it is unaffected by tau_axis).
    signed_population_names = {n for n in signed_fig_names if "population" in n}
    abs_population_names = {n for n in abs_fig_names if "population" in n}
    assert not (signed_population_names & abs_population_names)

    # Signed captions name tau_peak but never the abs form |tau_peak|.
    signed_pop_captions = [f["caption"] for f in result_signed["figures"]
                            if "population" in os.path.basename(f["path"])]
    assert len(signed_pop_captions) == 2
    for caption in signed_pop_captions:
        assert "tau_peak" in caption
        assert "|tau_peak|" not in caption


# --------------------------------------------------------------------------- #
# Differential: an INDEPENDENT (never `_abs_interval`) recompute of the
# bootstrap-interval image under abs(), over a hand-built fixture carrying
# both a bracket-zero and a non-bracket interval.
# --------------------------------------------------------------------------- #

def _make_ok_row(subobs_key, tau_peak_s, tau_boot_lo, tau_boot_hi, r_proxy, pd_significance):
    row = {f: "" for f in m.CCF_MASTER_FIELDS}
    row.update(
        obsid=subobs_key, source=subobs_key, subobs_key=subobs_key, n_slices=10, status="ok",
        tau_peak_s=tau_peak_s, dcf_peak=0.5, dcf_zero=0.1, n_pairs_peak=20,
        dcf_uncorrected_denominator=False, perm_p=0.005,
        tau_boot_lo=tau_boot_lo, tau_boot_hi=tau_boot_hi,
        dcf_peak_boot_lo=0.3, dcf_peak_boot_hi=0.7,
        coherent_rotation=True, r_proxy=r_proxy, pd_significance=pd_significance,
        tau_peak_uniform_s=tau_peak_s, median_dt_s=100.0, lag_bin_s_used=100.0,
    )
    return row


def _write_ccf_master(path, rows):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=m.CCF_MASTER_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def _write_qu_dcf(path, subobs_key):
    provenance = {"obsid": subobs_key, "subobs_key": subobs_key, "n_slices": 10}
    with open(path, "w", newline="") as f:
        f.write(m._QU_DCF_HEADER_PREFIX + json.dumps(provenance) + "\n")
        writer = csv.DictWriter(f, fieldnames=m._QU_DCF_FIELDS)
        writer.writeheader()
        for lag, dcf in ((-100.0, 0.1), (0.0, 0.2), (1.0, 0.5), (100.0, 0.3)):
            writer.writerow({"lag_s": lag, "dcf": dcf, "n_pairs": 5,
                              "perm_lo99": -0.4, "perm_hi99": 0.4})


def _write_population_csv(path, pairs):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=m.POPULATION_FIELDS)
        writer.writeheader()
        for pair, (n, rho, p) in pairs.items():
            writer.writerow({"pair": pair, "n": n, "rho": rho, "p": p})


# (subobs_key, tau_peak_s, tau_boot_lo, tau_boot_hi): one bracket-zero case
# ("aa"), one strictly-positive non-bracket case ("bb"), one strictly-
# negative non-bracket case ("cc").
_ROWS_SPEC = [
    ("aa", -6.0, -5.0, 10.0),
    ("bb", 5.0, 3.0, 8.0),
    ("cc", -7.0, -9.0, -4.0),
]


@pytest.fixture()
def hand_built_ccf(tmp_path):
    obs_db = tmp_path / "obs_db"
    ccf_rows = []
    for i, (subobs_key, tau_peak, lo, hi) in enumerate(_ROWS_SPEC):
        os.makedirs(obs_db / subobs_key, exist_ok=True)
        ccf_rows.append(_make_ok_row(subobs_key, tau_peak, lo, hi,
                                      r_proxy=float(i + 1), pd_significance=float(i + 2)))
        _write_qu_dcf(str(obs_db / subobs_key / f"qu_dcf_{subobs_key}.csv"), subobs_key)
    ccf_master_csv = obs_db / "ccf_master.csv"
    _write_ccf_master(str(ccf_master_csv), ccf_rows)
    population_csv = obs_db / "population_correlations.csv"
    _write_population_csv(str(population_csv), {
        "tau_peak_vs_r_proxy": (3, 0.111, 0.222),
        "tau_peak_vs_pd_significance": (3, 0.111, 0.222),
        "abs_tau_peak_vs_r_proxy": (3, 0.333, 0.444),
        "abs_tau_peak_vs_pd_significance": (3, 0.333, 0.444),
    })
    return str(obs_db), str(ccf_master_csv)


def test_differential_bootstrap_interval_abs_mapping(tmp_path, hand_built_ccf):
    obs_db, ccf_master_csv = hand_built_ccf
    # Non-vacuity: the fixture genuinely carries a bracket-zero interval and a
    # non-bracket one.
    assert any(lo <= 0.0 <= hi for _, _, lo, hi in _ROWS_SPEC)
    assert any(not (lo <= 0.0 <= hi) for _, _, lo, hi in _ROWS_SPEC)

    result = m.qu_lag_dcf_plot(obs_db, ccf_master_csv, str(tmp_path / "figs"),
                                save_path_prefix=str(tmp_path / "figs" / "plot"), tau_axis="abs")
    assert result["exit_code"] == 0

    def expected_interval(lo, hi):
        # Independent recompute of the true image of [lo, hi] under abs(),
        # written directly here -- never via qu_lag_dcf._abs_interval.
        if lo <= 0.0 <= hi:
            return 0.0, max(abs(lo), abs(hi))
        return min(abs(lo), abs(hi)), max(abs(lo), abs(hi))

    by_key = {k: (tau_peak, lo, hi) for k, tau_peak, lo, hi in _ROWS_SPEC}
    for key in POPULATION_KEYS:
        panel = result["plotted"][key]
        for i, subobs_key in enumerate(panel["subobs_key"]):
            tau_peak, lo, hi = by_key[subobs_key]
            exp_lo, exp_hi = expected_interval(lo, hi)
            assert panel["tau_boot_lo"][i] == pytest.approx(exp_lo, abs=1e-12)
            assert panel["tau_boot_hi"][i] == pytest.approx(exp_hi, abs=1e-12)
            assert panel["y"][i] == pytest.approx(abs(tau_peak), abs=1e-12)
        assert panel["spearman"]["rho"] == pytest.approx(0.333, abs=1e-12)
        assert panel["spearman"]["p"] == pytest.approx(0.444, abs=1e-12)


def _expected_abs_interval(lo, hi):
    """Independent recompute of the true image of `[lo, hi]` under `abs()`,
    written directly in this test module -- never via
    `qu_lag_dcf._abs_interval`. NaN in (either endpoint) -> NaN out."""
    if not (math.isfinite(lo) and math.isfinite(hi)):
        return float("nan"), float("nan")
    if lo <= 0.0 <= hi:
        return 0.0, max(abs(lo), abs(hi))
    return min(abs(lo), abs(hi)), max(abs(lo), abs(hi))


# --------------------------------------------------------------------------- #
# Artist-level capability: the FIGURE itself -- not just the returned
# `plotted` dict -- is genuinely re-mapped under tau_axis='abs'. Inspects the
# scatter offsets, error-bar Line2D artists, axis labels/limits and
# annotation text directly off the Axes built by `_build_population_figure`,
# then confirms the public `qu_lag_dcf_plot` entry point forwards `tau_axis`
# into exactly those same drawn values.
# --------------------------------------------------------------------------- #

_PANEL_SPECS = (
    ("r_proxy", "tau_peak_vs_r_proxy", "r_proxy (cps/arcmin^2)", "population_tau_r_proxy"),
    ("pd_significance", "tau_peak_vs_pd_significance", "pd_significance", "population_tau_pd_signif"),
)


def test_capability_abs_axis_artists_really_remapped(tmp_path, real_ccf_mixed):
    ok_rows = [r for r in m._read_ccf_master_csv(real_ccf_mixed["ccf_master_csv"]) if r["status"] == "ok"]
    population = m._read_population_correlations_csv(real_ccf_mixed["population_csv"])

    if not any(r["tau_peak_s"] < 0 for r in ok_rows):
        pytest.fail("fixture must carry a negative tau_peak_s or the abs-vs-signed "
                    "artist-level differencing check below is vacuous")

    result_abs = m.qu_lag_dcf_plot(real_ccf_mixed["obs_db"], real_ccf_mixed["ccf_master_csv"],
                                    str(tmp_path / "figs_public_abs"),
                                    save_path_prefix=str(tmp_path / "figs_public_abs" / "plot"),
                                    tau_axis="abs")
    assert result_abs["exit_code"] == 0

    for covariate_key, pair, xlabel, plotted_key in _PANEL_SPECS:
        fig_a, ax_a, _plotted_a = m._build_population_figure(
            ok_rows, population, covariate_key, pair, xlabel, tau_axis="abs")
        fig_s, ax_s, _plotted_s = m._build_population_figure(
            ok_rows, population, covariate_key, pair, xlabel, tau_axis="signed")
        try:
            # -- Scatter offsets: the actual drawn y (and untouched x). --
            off_a = np.asarray(ax_a.collections[0].get_offsets())
            off_s = np.asarray(ax_s.collections[0].get_offsets())
            expected_y_abs = np.array([abs(r["tau_peak_s"]) for r in ok_rows])
            assert off_a[:, 1] == pytest.approx(expected_y_abs, abs=1e-12)
            assert (off_a[:, 1] >= 0.0).all()
            assert np.any(np.abs(off_a[:, 1] - off_s[:, 1]) > 1e-9)
            assert off_a[:, 0] == pytest.approx(off_s[:, 0], abs=1e-12)

            # -- Error-bar Line2D artists: drawn (x, ylo, yhi) triples. --
            finite_rows = [r for r in ok_rows
                           if math.isfinite(r["tau_boot_lo"]) and math.isfinite(r["tau_boot_hi"])]
            expected_triples = sorted(
                (r[covariate_key],) + _expected_abs_interval(r["tau_boot_lo"], r["tau_boot_hi"])
                for r in finite_rows)
            drawn_triples = []
            for line in ax_a.lines:
                xdata = line.get_xdata()
                ylo, yhi = sorted(line.get_ydata())
                drawn_triples.append((float(xdata[0]), float(ylo), float(yhi)))
            drawn_triples.sort()
            assert len(drawn_triples) == len(expected_triples)
            for (ex, elo, ehi), (dx, dlo, dhi) in zip(expected_triples, drawn_triples):
                assert dx == pytest.approx(ex, abs=1e-12)
                assert dlo == pytest.approx(elo, abs=1e-12)
                assert dhi == pytest.approx(ehi, abs=1e-12)
                assert dlo >= 0.0

            # -- Axis contract. --
            assert ax_a.get_ylabel() == "|tau_peak| (s)"
            assert ax_a.get_ylim()[0] >= 0.0
            assert ax_s.get_ylabel() != ax_a.get_ylabel()

            # -- Annotation text: genuinely the abs Spearman row, not signed. --
            abs_pop, signed_pop = population["abs_" + pair], population[pair]
            assert abs_pop["rho"] != pytest.approx(signed_pop["rho"], abs=1e-3), (
                "fixture's abs and signed rho are indistinguishable at 1e-3 -- the "
                "annotation-text discriminating assertion below would be vacuous")
            texts = [t.get_text() for t in ax_a.texts]
            spearman_texts = [t for t in texts if "Spearman" in t]
            assert len(spearman_texts) == 1
            expected_sub = f"rho={abs_pop['rho']:.3f}"
            not_expected_sub = f"rho={signed_pop['rho']:.3f}"
            assert expected_sub in spearman_texts[0]
            assert not_expected_sub not in spearman_texts[0]

            # -- Close the loop: the public entry forwards tau_axis into the
            # exact artists it then saves. --
            assert result_abs["plotted"][plotted_key]["y"] == pytest.approx(
                list(off_a[:, 1]), abs=1e-12)
        finally:
            plt.close(fig_a)
            plt.close(fig_s)


# --------------------------------------------------------------------------- #
# PNG-byte differential: pixel-level default-path invariance (omitted ==
# tau_axis='signed', byte-identical) and pixel-level proof the abs figure
# content genuinely differs from the signed one (not just its filename).
# dcf_panels is byte-identical across all three calls -- untouched by tau_axis.
# --------------------------------------------------------------------------- #

def _shas_by_suffix(figs):
    return {os.path.basename(f["path"]).rsplit("plot_", 1)[-1]: f["sha256"] for f in figs}


def test_png_bytes_default_path_matches_signed_and_abs_differs(tmp_path, real_ccf_mixed):
    result_omitted = m.qu_lag_dcf_plot(real_ccf_mixed["obs_db"], real_ccf_mixed["ccf_master_csv"],
                                        str(tmp_path / "figs_omit"),
                                        save_path_prefix=str(tmp_path / "figs_omit" / "plot"))
    result_signed = m.qu_lag_dcf_plot(real_ccf_mixed["obs_db"], real_ccf_mixed["ccf_master_csv"],
                                       str(tmp_path / "figs_signed"),
                                       save_path_prefix=str(tmp_path / "figs_signed" / "plot"),
                                       tau_axis="signed")
    result_abs = m.qu_lag_dcf_plot(real_ccf_mixed["obs_db"], real_ccf_mixed["ccf_master_csv"],
                                    str(tmp_path / "figs_abs"),
                                    save_path_prefix=str(tmp_path / "figs_abs" / "plot"),
                                    tau_axis="abs")
    for result in (result_omitted, result_signed, result_abs):
        assert result["exit_code"] == 0

    sha_omit = _shas_by_suffix(result_omitted["figures"])
    sha_signed = _shas_by_suffix(result_signed["figures"])
    sha_abs = _shas_by_suffix(result_abs["figures"])

    # dcf_panels is untouched by tau_axis: byte-identical across all three calls.
    assert sha_omit["dcf_panels.png"] == sha_signed["dcf_panels.png"] == sha_abs["dcf_panels.png"]

    # Pixel-level default-path invariance: omitted == explicit tau_axis='signed'.
    assert sha_omit["population_tau_r_proxy.png"] == sha_signed["population_tau_r_proxy.png"]
    assert sha_omit["population_tau_pd_signif.png"] == sha_signed["population_tau_pd_signif.png"]

    # The abs population figures are genuinely different pixels, not just a
    # different filename.
    assert sha_signed["population_tau_r_proxy.png"] != sha_abs["population_tau_r_proxy_abs.png"]
    assert sha_signed["population_tau_pd_signif.png"] != sha_abs["population_tau_pd_signif_abs.png"]


# --------------------------------------------------------------------------- #
# Default-path invariance: tau_axis='signed' (explicit) and omitted both
# reproduce the pre-existing signed behaviour exactly -- same y values (raw
# signed tau_peak_s), same annotation source (tau_peak_vs_* rows), same
# (un-suffixed) output filenames.
# --------------------------------------------------------------------------- #

def test_default_path_invariance_signed_and_omitted_match(tmp_path, real_ccf_mixed):
    result_omitted = m.qu_lag_dcf_plot(real_ccf_mixed["obs_db"], real_ccf_mixed["ccf_master_csv"],
                                        str(tmp_path / "figs_omit"),
                                        save_path_prefix=str(tmp_path / "figs_omit" / "plot"))
    result_explicit = m.qu_lag_dcf_plot(real_ccf_mixed["obs_db"], real_ccf_mixed["ccf_master_csv"],
                                         str(tmp_path / "figs_explicit"),
                                         save_path_prefix=str(tmp_path / "figs_explicit" / "plot"),
                                         tau_axis="signed")
    assert result_omitted["exit_code"] == 0 and result_explicit["exit_code"] == 0

    with open(real_ccf_mixed["ccf_master_csv"], newline="") as f:
        ccf_rows = {r["subobs_key"]: r for r in csv.DictReader(f) if r["status"] == "ok"}
    with open(real_ccf_mixed["population_csv"], newline="") as f:
        pop_rows = {r["pair"]: r for r in csv.DictReader(f)}

    for result in (result_omitted, result_explicit):
        for key in POPULATION_KEYS:
            panel = result["plotted"][key]
            for i, subobs_key in enumerate(panel["subobs_key"]):
                row = ccf_rows[subobs_key]
                assert panel["y"][i] == pytest.approx(float(row["tau_peak_s"]), rel=1e-9, abs=1e-12)
                assert panel["tau_boot_lo"][i] == pytest.approx(float(row["tau_boot_lo"]), rel=1e-9, abs=1e-12)
                assert panel["tau_boot_hi"][i] == pytest.approx(float(row["tau_boot_hi"]), rel=1e-9, abs=1e-12)
            exp = pop_rows[SIGNED_PAIR_OF[key]]
            assert panel["spearman"]["n"] == int(float(exp["n"]))
            exp_rho = float(exp["rho"])
            if not math.isnan(exp_rho):
                assert panel["spearman"]["rho"] == pytest.approx(exp_rho, rel=1e-9, abs=1e-12)

        fig_suffixes = {os.path.basename(f["path"]).rsplit("plot_", 1)[-1] for f in result["figures"]}
        assert fig_suffixes == {"dcf_panels.png", "population_tau_r_proxy.png",
                                 "population_tau_pd_signif.png"}


# --------------------------------------------------------------------------- #
# Error handling: a bad tau_axis raises ValueError naming the argument and
# the accepted values, before any file is written.
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("bad_value", ["ABS", "absolute", "", "Signed"])
def test_bad_tau_axis_raises_before_any_file_written(tmp_path, real_ccf_mixed, bad_value):
    figs_dir = tmp_path / f"figs_bad_{abs(hash(bad_value))}"
    with pytest.raises(ValueError) as excinfo:
        m.qu_lag_dcf_plot(real_ccf_mixed["obs_db"], real_ccf_mixed["ccf_master_csv"],
                           str(figs_dir), save_path_prefix=str(figs_dir / "plot"),
                           tau_axis=bad_value)
    msg = str(excinfo.value)
    assert "tau_axis" in msg
    assert "signed" in msg and "abs" in msg
    assert not figs_dir.exists()


# --------------------------------------------------------------------------- #
# Error handling: tau_axis='abs' against a population_correlations.csv that
# predates the abs_tau_peak_* population pairs raises ValueError naming the
# missing pair and pointing the caller at ixpe_qu_lag_dcf; writes no PNG.
# --------------------------------------------------------------------------- #

def test_abs_axis_against_older_population_csv_raises_and_writes_nothing(tmp_path):
    obs_db = tmp_path / "obs_db"
    os.makedirs(obs_db / "00000001", exist_ok=True)
    ccf_master_csv = obs_db / "ccf_master.csv"
    _write_ccf_master(str(ccf_master_csv),
                       [_make_ok_row("00000001", 1.0, 5.0, 10.0, 2.0, 3.0)])
    _write_qu_dcf(str(obs_db / "00000001" / "qu_dcf_00000001.csv"), "00000001")
    # Older-format population CSV: only the two signed pairs, no abs_tau_peak_* rows.
    _write_population_csv(str(obs_db / "population_correlations.csv"), {
        "tau_peak_vs_r_proxy": (1, 0.5, 0.1),
        "tau_peak_vs_pd_significance": (1, 0.5, 0.1),
    })

    figs_dir = tmp_path / "figs"
    with pytest.raises(ValueError) as excinfo:
        m.qu_lag_dcf_plot(str(obs_db), str(ccf_master_csv), str(figs_dir),
                           save_path_prefix=str(figs_dir / "plot"), tau_axis="abs")
    msg = str(excinfo.value)
    assert "abs_tau_peak_vs_r_proxy" in msg
    assert "ixpe_qu_lag_dcf" in msg
    assert not figs_dir.exists()


# --------------------------------------------------------------------------- #
# Reproducibility: two successive tau_axis='abs' calls over the same inputs
# return identical plotted-data (y values and interval endpoints) and
# byte-identical PNGs per suffix.
# --------------------------------------------------------------------------- #

def test_repro_qu_lag_dcf_plot_abs(tmp_path, real_ccf_mixed):
    a = m.qu_lag_dcf_plot(real_ccf_mixed["obs_db"], real_ccf_mixed["ccf_master_csv"],
                           str(tmp_path / "figs_a"), save_path_prefix=str(tmp_path / "figs_a" / "plot"),
                           tau_axis="abs")
    b = m.qu_lag_dcf_plot(real_ccf_mixed["obs_db"], real_ccf_mixed["ccf_master_csv"],
                           str(tmp_path / "figs_b"), save_path_prefix=str(tmp_path / "figs_b" / "plot"),
                           tau_axis="abs")
    assert a["exit_code"] == 0 and b["exit_code"] == 0

    for key in POPULATION_KEYS:
        pa, pb = a["plotted"][key], b["plotted"][key]
        assert pa["subobs_key"] == pb["subobs_key"]
        assert pa["y"] == pytest.approx(pb["y"], abs=1e-12)
        assert pa["tau_boot_lo"] == pytest.approx(pb["tau_boot_lo"], abs=1e-12)
        assert pa["tau_boot_hi"] == pytest.approx(pb["tau_boot_hi"], abs=1e-12)
        assert pa["spearman"] == pb["spearman"]

    def by_suffix(figs):
        return {os.path.basename(f["path"]).rsplit("plot_", 1)[-1]: f["sha256"] for f in figs}

    sha_a, sha_b = by_suffix(a["figures"]), by_suffix(b["figures"])
    expected_suffixes = {"dcf_panels.png", "population_tau_r_proxy_abs.png",
                          "population_tau_pd_signif_abs.png"}
    assert set(sha_a) == expected_suffixes == set(sha_b)
    for suffix in expected_suffixes:
        assert sha_a[suffix] == sha_b[suffix], f"figure {suffix} not byte-identical across repeated calls"


# --------------------------------------------------------------------------- #
# Edge cases under tau_axis='abs': a non-finite tau_peak_s (y=NaN, no crash),
# and a non-finite bootstrap endpoint (suppresses only that row's own bar --
# both endpoints reported/drawn as NaN, never a partial/garbage interval).
# --------------------------------------------------------------------------- #

_NONFINITE_ROWS_SPEC = [
    ("ff", -4.0, -1.0, 6.0),
    ("nan_tau", float("nan"), -2.0, 3.0),
    ("nan_lo", 5.0, float("nan"), 3.0),
    ("nan_hi", -5.0, 1.0, float("nan")),
]


@pytest.fixture()
def nonfinite_ccf(tmp_path):
    obs_db = tmp_path / "obs_db"
    ccf_rows = []
    for i, (subobs_key, tau_peak, lo, hi) in enumerate(_NONFINITE_ROWS_SPEC):
        os.makedirs(obs_db / subobs_key, exist_ok=True)
        ccf_rows.append(_make_ok_row(subobs_key, tau_peak, lo, hi,
                                      r_proxy=float(i + 1), pd_significance=float(i + 2)))
        _write_qu_dcf(str(obs_db / subobs_key / f"qu_dcf_{subobs_key}.csv"), subobs_key)
    ccf_master_csv = obs_db / "ccf_master.csv"
    _write_ccf_master(str(ccf_master_csv), ccf_rows)
    population_csv = obs_db / "population_correlations.csv"
    _write_population_csv(str(population_csv), {
        "tau_peak_vs_r_proxy": (4, 0.111, 0.222),
        "tau_peak_vs_pd_significance": (4, 0.111, 0.222),
        "abs_tau_peak_vs_r_proxy": (4, 0.333, 0.444),
        "abs_tau_peak_vs_pd_significance": (4, 0.333, 0.444),
    })
    return str(obs_db), str(ccf_master_csv), str(population_csv)


def test_abs_axis_nonfinite_rows_do_not_crash(tmp_path, nonfinite_ccf):
    obs_db, ccf_master_csv, population_csv = nonfinite_ccf

    result = m.qu_lag_dcf_plot(obs_db, ccf_master_csv, str(tmp_path / "figs"),
                                save_path_prefix=str(tmp_path / "figs" / "plot"), tau_axis="abs")
    assert result["exit_code"] == 0
    for fig in result["figures"]:
        assert os.path.exists(fig["path"])
        assert os.path.getsize(fig["path"]) > 0

    by_key = {k: (tau_peak, lo, hi) for k, tau_peak, lo, hi in _NONFINITE_ROWS_SPEC}
    for key in POPULATION_KEYS:
        panel = result["plotted"][key]
        for i, subobs_key in enumerate(panel["subobs_key"]):
            tau_peak, lo, hi = by_key[subobs_key]
            if subobs_key == "nan_tau":
                assert math.isnan(panel["y"][i])
            else:
                assert panel["y"][i] == pytest.approx(abs(tau_peak), abs=1e-12)
            if subobs_key in ("nan_lo", "nan_hi"):
                assert math.isnan(panel["tau_boot_lo"][i])
                assert math.isnan(panel["tau_boot_hi"][i])
            else:
                assert math.isfinite(panel["tau_boot_lo"][i])
                assert math.isfinite(panel["tau_boot_hi"][i])

    # Artist level: only the two fully-finite-interval rows ("ff", "nan_tau")
    # draw an error bar, identically under either axis.
    ok_rows = [r for r in m._read_ccf_master_csv(ccf_master_csv) if r["status"] == "ok"]
    population = m._read_population_correlations_csv(population_csv)
    for tau_axis in ("abs", "signed"):
        fig, ax, _plotted = m._build_population_figure(
            ok_rows, population, "r_proxy", "tau_peak_vs_r_proxy", "r_proxy (cps/arcmin^2)",
            tau_axis=tau_axis)
        try:
            assert len(ax.lines) == 2
        finally:
            plt.close(fig)


# --------------------------------------------------------------------------- #
# Manifest wiring: tau_axis is documented in args_schema (AI-facing half).
# --------------------------------------------------------------------------- #

def test_tool_manifest_documents_tau_axis():
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.ixpe_qu_lag_dcf_plot"]
    props = entry["args_schema"]["properties"]
    assert "tau_axis" in props
    assert set(props["tau_axis"].get("enum", [])) == {"signed", "abs"}
    desc = props["tau_axis"]["description"]
    assert "abs_tau_peak" in desc
    assert "_abs" in desc
    assert "ixpe_qu_lag_dcf" in desc
    top_desc = entry["args_schema"]["description"]
    assert "tau_axis" in top_desc
