"""What the xspec_utils tools write, and what they return for it.

Each tool writes under the output argument it is given (a folder of the test's own, apart from the
observation database, which stays read-only) and returns every file a report may cite under
``products`` (``{"path", "caption"}``) or ``figures``.
"""

from __future__ import annotations

import hashlib
import os
import runpy
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from xspec_utils.lib import xspec_help
from xspec_utils.skills import session

BOOT = Path(__file__).resolve().parents[2] / "envs" / "xspec" / "boot.py"
IXPE_ROOT = Path("/data/fzhou/XSPEC/Cyg_X-1/IXPE")
IXPE_OBS = "01002901"  # its PHA headers name responses of the human astropy environment's ixpeobssim
NUSTAR = Path("/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/30702017002/products")
SPECTRUM = NUSTAR / "nu30702017002A01_sr_grp.pha"  # RESPFILE absolute; ANCRFILE, BACKFILE relative
EDGES = [4.0, 10.0, 20.0, 50.0]


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def out(tmp_path: Path) -> Path:
    folder = tmp_path / "out"
    folder.mkdir()
    return folder


def _synthetic_pha(folder: Path, **keys: str) -> Path:
    """A ten-channel PHA (1-11 keV, 1 keV channels) and its RMF's EBOUNDS, in ``folder``."""
    lo = np.arange(1.0, 11.0)
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns([
        fits.Column(name="CHANNEL", format="J", array=np.arange(1, 11)),
        fits.Column(name="E_MIN", format="E", array=lo),
        fits.Column(name="E_MAX", format="E", array=lo + 1.0),
    ], name="EBOUNDS")]).writeto(folder / "resp.rmf")
    spectrum = fits.BinTableHDU.from_columns([
        fits.Column(name="CHANNEL", format="J", array=np.arange(1, 11)),
        fits.Column(name="COUNTS", format="J", array=np.full(10, 5)),
    ], name="SPECTRUM")
    for key, value in keys.items():
        spectrum.header[key] = value
    path = folder / "src.pha"
    fits.HDUList([fits.PrimaryHDU(), spectrum]).writeto(path)
    return path


def test_group_pha_groups_and_resolves_only_names_found_next_to_the_input(tmp_path, out):
    data = tmp_path / "data"
    data.mkdir()
    infile = _synthetic_pha(data, RESPFILE="resp.rmf", ANCRFILE="none", BACKFILE="missing_bk.pha")
    before = _sha(infile)

    outfile = out / "grouped.pha"
    result = xspec_help.group_pha(str(infile), str(outfile), [2.0, 5.0, 9.0])

    assert result["exit_code"] == 0
    assert result["products"] == [{"path": str(outfile), "caption": result["products"][0]["caption"]}]
    assert result["products"][0]["caption"]
    assert sorted(p.name for p in out.iterdir()) == ["grouped.pha"]
    assert _sha(infile) == before
    with fits.open(outfile) as hdul:
        spectrum = hdul["SPECTRUM"]
        # Channel midpoints 2.5-4.5 keV and 5.5-8.5 keV: one group per pair of edges.
        assert spectrum.data["GROUPING"].tolist() == [0, 1, -1, -1, 1, -1, -1, -1, 0, 0]
        header = spectrum.header
    assert header["RESPFILE"] == str(data / "resp.rmf")
    assert header["ANCRFILE"] == "none"
    assert header["BACKFILE"] == "missing_bk.pha"  # names no file next to the input: left as written


# A response path ixpeobssim wrote into obs-db PHA files from a pip --user install since removed.
GONE_RMF = ("/grad/fzhou/.local/lib/python3.13/site-packages/ixpeobssim/caldb/ixpe/gpd/cpf/rmf/"
            "ixpe_d1_obssim_alpha075_v012.rmf")


def test_a_response_under_a_removed_ixpeobssim_is_read_from_the_installed_one(tmp_path, out):
    import ixpeobssim

    assert not Path(GONE_RMF).exists()
    installed = Path(ixpeobssim.IXPEOBSSIM_CALDB) / "ixpe/gpd/cpf/rmf/ixpe_d1_obssim_alpha075_v012.rmf"
    assert xspec_help.resolve_response(GONE_RMF) == str(installed)
    # A name under another ixpeobssim CALDB is read from the installed one even when it exists.
    other = tmp_path / "env" / "ixpeobssim" / "caldb" / "ixpe/gpd/cpf/rmf/ixpe_d1_obssim_alpha075_v012.rmf"
    other.parent.mkdir(parents=True)
    other.write_bytes(b"")
    assert xspec_help.resolve_response(str(other)) == str(installed)
    assert xspec_help.resolve_response(str(installed)) == str(installed)
    alias = os.path.realpath(installed)  # the installed file under another spelling: kept as named
    assert xspec_help.resolve_response(alias) == alias
    assert xspec_help.resolve_response(GONE_RMF.replace("caldb/", "caldb//")) == str(installed)
    for kept in ("none", "", "resp.rmf", str(NUSTAR / "nu30702017002A01_sr.rmf")):
        assert xspec_help.resolve_response(kept) == kept
    missing = GONE_RMF.replace("v012", "v999")
    with pytest.raises(FileNotFoundError) as excinfo:
        xspec_help.resolve_response(missing)
    assert missing in str(excinfo.value) and str(installed).replace("v012", "v999") in str(excinfo.value)

    with fits.open(installed) as rmf:
        n = len(rmf["EBOUNDS"].data)
    spectrum = fits.BinTableHDU.from_columns([
        fits.Column(name="CHANNEL", format="J", array=np.arange(n)),
        fits.Column(name="COUNTS", format="J", array=np.ones(n, dtype=int)),
    ], name="SPECTRUM")
    spectrum.header["RESPFILE"] = GONE_RMF
    infile = tmp_path / "src.pha"
    fits.HDUList([fits.PrimaryHDU(), spectrum]).writeto(infile)

    result = xspec_help.group_pha(str(infile), str(out / "grouped.pha"), [2.0, 4.0, 8.0])
    assert fits.getheader(out / "grouped.pha", "SPECTRUM")["RESPFILE"] == str(installed)
    assert "src.pha RESPFILE: ixpe/gpd/cpf/rmf/ixpe_d1_obssim_alpha075_v012.rmf" in result["stdout"]
    assert ixpeobssim.IXPEOBSSIM_CALDB in result["stdout"]
    assert fits.getheader(infile, "SPECTRUM")["RESPFILE"] == GONE_RMF


@pytest.fixture(scope="module")
def xspec_session():
    """XSPEC as the environment's REPL boot file sets it up (session options, local models)."""
    runpy.run_path(str(BOOT))
    from xspec import AllData, AllModels

    yield AllData
    AllData.clear()
    AllModels.clear()


@pytest.fixture
def nustar_session(xspec_session):
    """One NuSTAR FPMA spectrum in a cleared XSPEC session."""
    from xspec import AllData, AllModels

    AllData.clear()
    AllModels.clear()
    AllData(f"1:1 {SPECTRUM}")
    spectrum = AllData(1)
    spectrum.background = str(SPECTRUM.with_name("nu30702017002A01_bk.pha"))
    spectrum.response = str(SPECTRUM.with_name("nu30702017002A01_sr.rmf"))
    spectrum.response.arf = str(SPECTRUM.with_name("nu30702017002A01_sr.arf"))
    spectrum.ignore("**-4.0 50.0-**")
    return {"NuSTAR": {"30702017002": {"FPMA": [1]}}}


@pytest.mark.slow
def test_a_grouped_nustar_copy_loads_from_another_folder(xspec_session, out, tmp_path, monkeypatch):
    from xspec import AllData, AllModels

    outfile = out / "grouped.pha"
    xspec_help.group_pha(str(SPECTRUM), str(outfile), EDGES)
    with fits.open(SPECTRUM) as hdul:
        original = hdul["SPECTRUM"].header
    AllData.clear()
    AllModels.clear()
    monkeypatch.chdir(tmp_path)  # neither the input's folder nor the copy's
    AllData(f"1:1 {outfile}")
    spectrum = AllData(1)
    assert spectrum.background.fileName == str(NUSTAR / original["BACKFILE"])
    assert spectrum.response.arf == str(NUSTAR / original["ANCRFILE"])
    AllData.clear()


@pytest.mark.slow
def test_run_fit_and_plot_fit_return_what_they_write(nustar_session, out):
    from xspec import Model

    Model("tbabs*powerlaw")
    log = out / "run_fit.log"
    fit = session.run_fit(n_iterations=20, log_path=str(log))
    assert fit["exit_code"] == 0
    assert [p["path"] for p in fit["products"]] == [str(log)]
    assert log.stat().st_size > 0
    assert str(fit["fit_stat"]) in fit["stdout"] and str(fit["stat_per_dof"]) in fit["stdout"]

    plot = session.plot_fit(nustar_session, save_path=str(out / "fit"))  # a name with no extension
    png = out / "fit.png"
    assert plot["exit_code"] == 0
    assert [f["path"] for f in plot["figures"]] == [str(png)]
    assert plot["figures"][0]["sha256"] == _sha(png)
    assert sorted(p.name for p in out.iterdir()) == ["fit.png", "run_fit.log"]


@pytest.mark.slow
def test_fit_tbfeo_hss_returns_its_log(nustar_session, out):
    log = out / "fit_tbfeo_hss.log"
    result = xspec_help.fit_tbfeo_plabs_polconst_diskbb_polconst_relxillCp(
        nustar_session, log_path=str(log))
    assert result["exit_code"] == 0
    assert [p["path"] for p in result["products"]] == [str(log)]
    assert str(SPECTRUM) in log.read_text()  # the audit block heads the log
    assert str(result["fit_stat"]) in result["stdout"]
    assert sorted(p.name for p in out.iterdir()) == ["fit_tbfeo_hss.log"]


@pytest.mark.slow
def test_plot_iqu_returns_one_figure_per_du(xspec_session, out, capsys):
    import ixpeobssim
    from xspec import AllData, Model

    indices = xspec_help.load_data({"IXPE": [IXPE_OBS], "STATE": "test"}, str(IXPE_ROOT))
    # Responses that exist where the headers name them are read from the installed CALDB too.
    pha = IXPE_ROOT / IXPE_OBS / "event_l2"
    named = fits.getheader(next(pha.glob("*_det1_*_src_pha1.fits")), "SPECTRUM")["RESPFILE"]
    assert Path(named).is_file(), f"fixture precondition: {named}, which {IXPE_OBS}'s header names, is gone"
    assert not named.startswith(ixpeobssim.IXPEOBSSIM_CALDB)
    assert AllData(1).response.rmf.startswith(ixpeobssim.IXPEOBSSIM_CALDB)
    assert capsys.readouterr().out.count("Response: ") == 18
    Model("polconst*powerlaw")
    session.run_fit(n_iterations=10)
    result = session.plot_iqu(indices, save_path_prefix=str(out / "iqu.png"))
    assert result["exit_code"] == 0
    paths = [out / f"iqu_DU{n}.png" for n in (1, 2, 3)]
    assert [f["path"] for f in result["figures"]] == [str(p) for p in paths]
    assert [f["sha256"] for f in result["figures"]] == [_sha(p) for p in paths]
    assert sorted(p.name for p in out.iterdir()) == [p.name for p in paths]


@pytest.mark.slow
def test_load_data_reads_ixpe_responses_from_the_installed_ixpeobssim(xspec_session, capsys):
    import ixpeobssim
    from xspec import AllData

    obs = "02008201"  # its PHA headers name responses under the removed install
    pha = IXPE_ROOT / obs / "event_l2" / f"ixpe{obs}_det1_evt2_v01_rej_src_pha1.fits"
    assert fits.getheader(pha, "SPECTRUM")["RESPFILE"].startswith("/grad/fzhou/.local/")
    indices = xspec_help.load_data({"IXPE": [obs], "STATE": "test"}, str(IXPE_ROOT))
    printed = capsys.readouterr().out

    assert AllData.nSpectra == 9
    caldb = ixpeobssim.IXPEOBSSIM_CALDB
    for du in ("DU1", "DU2", "DU3"):
        for index in indices["IXPE"][obs][du]:
            response = AllData(index).response
            assert response.rmf.startswith(caldb) and Path(response.rmf).is_file()
            assert response.arf.startswith(caldb) and Path(response.arf).is_file()
    assert printed.count("Response: ") == 18  # RESPFILE and ANCRFILE of 9 spectra
    assert f"{pha.name} RESPFILE: {Path(AllData(1).response.rmf).relative_to(caldb)}" in printed


@pytest.mark.slow
def test_sum_counts_reads_ixpe_responses_from_the_installed_ixpeobssim(xspec_session, capsys):
    import ixpeobssim
    from xspec import AllData

    from xspec_utils.lib import red_ixpe

    obs = "02008201"
    folder = IXPE_ROOT / obs / "event_l2"
    src = [str(p) for p in sorted(folder.glob(f"ixpe{obs}_det*_evt2_v01_rej_src_pha1.fits"))]
    bkg = [str(p) for p in sorted(folder.glob(f"ixpe{obs}_det*_evt2_v01_rej_bkg_pha1.fits"))]
    rate, err = red_ixpe._sum_counts(obs, src, bkg, (2.0, 8.0))
    printed = capsys.readouterr().out
    assert AllData.nSpectra == 3
    for index in (1, 2, 3):
        assert AllData(index).response.rmf.startswith(ixpeobssim.IXPEOBSSIM_CALDB)
    assert printed.count("Response: ") == 6
    # Every loaded spectrum counts, the last one too.
    assert rate == pytest.approx(sum(AllData(index).rate[0] for index in (1, 2, 3)), rel=1e-12)
    # The error from the files themselves: each noticed channel's STAT_ERR, source plus background
    # scaled by BACKSCAL, added in quadrature over the three spectra.
    variance = 0.0
    for index, (s, b) in enumerate(zip(src, bkg), start=1):
        rows = np.asarray(AllData(index).noticed) - 1  # 1-based channel numbers
        with fits.open(s) as hs, fits.open(b) as hb:
            scale = hs["SPECTRUM"].header["BACKSCAL"] / hb["SPECTRUM"].header["BACKSCAL"]
            variance += (hs["SPECTRUM"].data["STAT_ERR"][rows] ** 2).sum()
            variance += scale ** 2 * (hb["SPECTRUM"].data["STAT_ERR"][rows] ** 2).sum()
    assert err == pytest.approx(variance ** 0.5, rel=1e-6)


@pytest.mark.slow
def test_a_spectrum_missing_only_its_arf_keeps_its_rmf(xspec_session, tmp_path):
    import ixpeobssim
    from xspec import AllData

    pha = IXPE_ROOT / "02008201" / "event_l2" / "ixpe02008201_det1_evt2_v01_rej_src_pha1.fits"
    rmf = str(Path(ixpeobssim.IXPEOBSSIM_CALDB) / "ixpe/gpd/cpf/rmf/ixpe_d1_obssim_alpha075_v012.rmf")
    copy = tmp_path / pha.name
    with fits.open(pha) as hdul:
        hdul["SPECTRUM"].header["RESPFILE"] = rmf  # a response that exists; ANCRFILE stays gone
        hdul.writeto(copy)
    AllData.clear()
    AllData(f"1:1 {copy}")
    lines = xspec_help.assign_responses(1, str(copy))
    assert len(lines) == 1 and f"{copy.name} ANCRFILE: " in lines[0]
    assert AllData(1).response.rmf == rmf
    assert AllData(1).response.arf.startswith(ixpeobssim.IXPEOBSSIM_CALDB)
    AllData.clear()


@pytest.mark.slow
def test_audit_paths_labels_each_rate(xspec_session, capsys):
    from xspec import AllData, AllModels, Model

    pha = IXPE_ROOT / "02008201" / "event_l2" / "ixpe02008201_det1_evt2_v01_rej_src_pha1.fits"
    bkg = pha.with_name(pha.name.replace("_src_", "_bkg_"))
    AllData.clear()
    AllModels.clear()
    AllData(f"1:1 {pha}")
    AllData(1).background = str(bkg)
    xspec_help.assign_responses(1, str(pha))
    AllData(1).ignore("**-2.0 8.0-**")
    capsys.readouterr()
    xspec_help.audit_paths()
    before_model = capsys.readouterr().out  # no model defined: no model rate
    Model("powerlaw")
    xspec_help.audit_paths()
    printed = capsys.readouterr().out

    # Each rate from the files and the model themselves, over the noticed channels.
    rows = np.asarray(AllData(1).noticed) - 1
    with fits.open(pha) as hs, fits.open(bkg) as hb:
        src, back = hs["SPECTRUM"], hb["SPECTRUM"]
        scale = src.header["BACKSCAL"] / back.header["BACKSCAL"]
        total = src.data["RATE"][rows].sum()
        net = total - scale * back.data["RATE"][rows].sum()
        sigma = ((src.data["STAT_ERR"][rows] ** 2).sum()
                 + scale ** 2 * (back.data["STAT_ERR"][rows] ** 2).sum()) ** 0.5
    model = float(np.sum(AllModels(1).folded(1)))
    assert f"Rate: net={net:.4g} ± {sigma:.4g} cts/s  (total={total:.4g})\n" in before_model + "\n"
    assert f"Rate: net={net:.4g} ± {sigma:.4g} cts/s  (total={total:.4g}, current model={model:.4g})" in printed
    AllData.clear()
    AllModels.clear()


def _mtimes(folder: Path) -> dict[str, int]:
    return {p.name: p.lstat().st_mtime_ns for p in folder.iterdir()}


def _linked(sources: list[str], folder: Path) -> list[str]:
    """The inputs linked into a folder of the test's own, so a write beside an input lands there."""
    folder.mkdir(parents=True, exist_ok=True)
    for source in sources:
        (folder / Path(source).name).symlink_to(Path(source).resolve())
    return [str(folder / Path(s).name) for s in sources]


@pytest.mark.slow
def test_run_extraction_writes_regions_and_png_to_its_folder(tmp_path, out):
    from xspec_utils.lib import find_src

    source = "/data/fzhou/XSPEC/Cyg_X-1/NuSTAR/30702017002/event_cl/nu30702017002A01_cl.evt"
    evt = Path(_linked([source], tmp_path / "event_cl")[0])
    before = _mtimes(evt.parent)
    result = find_src.run_extraction(str(evt), str(out), 1.0, 2.0, 3.0)

    assert result["exit_code"] == 0
    assert [p["path"] for p in result["products"]] == [str(out / "src.reg"), str(out / "bkg.reg")]
    assert [f["path"] for f in result["figures"]] == [str(out / "verify.png")]
    assert result["figures"][0]["sha256"] == _sha(out / "verify.png")
    assert sorted(p.name for p in out.iterdir()) == ["bkg.reg", "src.reg", "verify.png"]
    assert f"circle({result['ra']:.6f}, {result['dec']:.6f}, 60.00\")" in (out / "src.reg").read_text()
    assert _mtimes(evt.parent) == before


IXPE_PCUBE_OBS = "02008201"


def _rej_events(root: Path, obs: str) -> list[str]:
    return [str(p) for p in sorted((root / obs / "event_l2").glob(f"ixpe{obs}_det*_evt2_v01_rej.fits"))]


@pytest.mark.slow
def test_pcube_writes_every_file_under_its_folder(tmp_path, out):
    from xspec_utils.lib import red_ixpe, src_pos

    obs = IXPE_PCUBE_OBS
    events = _linked(_rej_events(IXPE_ROOT, obs), tmp_path / obs / "event_l2")
    before = _mtimes(tmp_path / obs / "event_l2")
    ra, dec = src_pos.SRC_COORD[obs]
    products: list = []
    src, _bkg = red_ixpe.pcube(obs, events, ebinning=[2.0, 4.0, 8.0], src_rad=1.0, bkg_inner_rad=2.5,
                               bkg_outer_rad=5.0, ra=float(ra), dec=float(dec), irf_date="20230101",
                               outdir=str(out), products=products)

    assert _mtimes(tmp_path / obs / "event_l2") == before
    assert len(products) == 12 and all(set(p) == {"path", "caption"} and p["caption"] for p in products)
    assert {Path(p["path"]).parent for p in products} == {out}
    assert all(Path(p["path"]).is_file() and not Path(p["path"]).is_symlink() for p in products)
    links = [out / Path(e).name for e in events]
    assert all(link.is_symlink() and link.resolve() == Path(e).resolve() for link, e in zip(links, events))
    # Every file in the folder is a link to an input or one the call declares.
    assert {p for p in out.iterdir() if not p.is_symlink()} == {Path(q["path"]) for q in products}

    reread: list = []
    again, _ = red_ixpe.pcube(obs, events, changes=False, outdir=str(out), products=reread)
    assert {q["path"] for q in reread} <= {q["path"] for q in products} and len(reread) == 6
    assert np.allclose(again.PD, src.PD)

    empty = tmp_path / "empty"
    with pytest.raises(FileNotFoundError, match="earlier call"):
        red_ixpe.pcube(obs, events, changes=False, outdir=str(empty))
    assert not empty.exists()


@pytest.mark.slow
def test_pha_writes_into_the_obs_db_and_records_its_spectra(tmp_path):
    from obsdb_catalog import neutral
    from obsdb_catalog.healthcheck import healthcheck

    from xspec_utils.lib import red_ixpe, src_pos

    obs = IXPE_PCUBE_OBS
    db = tmp_path / "IXPE"  # a fixture obs-db named after its mission
    event_l2 = db / obs / "event_l2"
    event_l2.mkdir(parents=True)
    for source in _rej_events(IXPE_ROOT, obs):
        (event_l2 / Path(source).name).symlink_to(Path(source).resolve())
    det1 = f"{obs}/event_l2/ixpe{obs}_det1_evt2_v01_rej.fits"
    neutral.new(str(db), "ixpe", obs, products=[{"kind": "level2_evt", "path": det1, "produced_by": "heasarc"}])
    assert healthcheck(str(db)) is True
    real = _mtimes(IXPE_ROOT / obs / "event_l2")

    ra, dec = src_pos.SRC_COORD[obs]
    src, bkg = red_ixpe.pha(obs, _rej_events(db, obs), src_rad=1.0, bkg_inner_rad=2.5, bkg_outer_rad=5.0,
                            ra=float(ra), dec=float(dec), irf_date="20230101")

    assert all(Path(p).parent == event_l2 and Path(p).is_file() for p in src + bkg)
    assert healthcheck(str(db)) is True
    entry = neutral.read(str(db))["entries"][obs]
    assert entry["max_level"] == 3
    pha = {p["path"]: p["params"] for p in entry["products"] if p["kind"] == "pha"}
    assert sorted(pha) == sorted(str(Path(f).relative_to(db)) for f in src + bkg)  # every Stokes I spectrum
    assert len(pha) == 6
    for params in pha.values():
        assert (params["src_rad_arcmin"], params["bkg_inner_arcmin"], params["bkg_outer_arcmin"]) == (1.0, 2.5, 5.0)
        assert params["irfname"] == "ixpe:obssim20230101_alpha075:v13"
    det1 = pha[f"{obs}/event_l2/ixpe{obs}_det1_evt2_v01_rej_src_pha1.fits"]
    assert (det1["du"], det1["filter"]) == ("DU1", "rej_src")
    assert _mtimes(IXPE_ROOT / obs / "event_l2") == real


def test_pha_refuses_before_writing_what_it_could_not_record(tmp_path):
    from xspec_utils.lib import red_ixpe

    obs = IXPE_PCUBE_OBS
    copy = tmp_path / f"{obs}_picorr" / "event_l2"  # a working copy outside the obs-db layout
    events = _linked(_rej_events(IXPE_ROOT, obs), copy)
    args = dict(src_rad=1.0, bkg_inner_rad=2.5, bkg_outer_rad=5.0, ra=299.59, dec=35.2, irf_date="20230101")
    before = _mtimes(copy)
    with pytest.raises(ValueError, match="event_l2"):
        red_ixpe.pha(obs, events, **args)
    with pytest.raises(ValueError, match="overwrite=True"):
        red_ixpe.pha(obs, events, overwrite=False, **args)
    assert _mtimes(copy) == before
