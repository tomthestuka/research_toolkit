"""What the xspec_utils tools write, and what they return for it.

Each tool writes under the output argument it is given (a folder of the test's own, apart from the
observation database, which stays read-only) and returns every file a report may cite under
``products`` (``{"path", "caption"}``) or ``figures``.
"""

from __future__ import annotations

import hashlib
import runpy
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from xspec_utils.lib import xspec_help
from xspec_utils.skills import session

BOOT = Path(__file__).resolve().parents[2] / "envs" / "xspec" / "boot.py"
IXPE_ROOT = Path("/data/fzhou/XSPEC/Cyg_X-1/IXPE")
IXPE_OBS = "01002901"  # its PHA headers name responses of an installed ixpeobssim
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
    assert GONE_RMF in result["stdout"] and str(installed) in result["stdout"]
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


def _ixpe_responses_readable() -> bool:
    pha = sorted((IXPE_ROOT / IXPE_OBS / "event_l2").glob("*_src_pha1*.fits"))
    names = [fits.getheader(p, "SPECTRUM").get(k) for p in pha for k in ("RESPFILE", "ANCRFILE")]
    return bool(names) and all(Path(n).is_file() for n in names)


@pytest.mark.slow
def test_plot_iqu_returns_one_figure_per_du(xspec_session, out):
    if not _ixpe_responses_readable():
        pytest.skip(f"the responses {IXPE_OBS}'s PHA headers name are not readable here")
    from xspec import Model

    indices = xspec_help.load_data({"IXPE": [IXPE_OBS], "STATE": "test"}, str(IXPE_ROOT))
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
    assert f"{pha.name} RESPFILE: {AllData(1).response.rmf}" in printed


@pytest.mark.slow
def test_sum_counts_reads_ixpe_responses_from_the_installed_ixpeobssim(xspec_session, capsys):
    import ixpeobssim
    from xspec import AllData

    from xspec_utils.lib import red_ixpe

    obs = "02008201"
    folder = IXPE_ROOT / obs / "event_l2"
    src = [str(p) for p in sorted(folder.glob(f"ixpe{obs}_det*_evt2_v01_rej_src_pha1.fits"))]
    bkg = [str(p) for p in sorted(folder.glob(f"ixpe{obs}_det*_evt2_v01_rej_bkg_pha1.fits"))]
    red_ixpe._sum_counts(obs, src, bkg, (2.0, 8.0))
    printed = capsys.readouterr().out
    assert AllData.nSpectra == 3
    for index in (1, 2, 3):
        assert AllData(index).response.rmf.startswith(ixpeobssim.IXPEOBSSIM_CALDB)
    assert printed.count("Response: ") == 6


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
