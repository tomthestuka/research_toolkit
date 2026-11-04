
import os
import glob
import re
import warnings
from pathlib import Path
from astropy.io import fits
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.binning.polarization import xBinnedPolarizationCube
from contextlib import redirect_stdout, redirect_stderr
from xspec import AllData, Spectrum, AllModels
from obsdb_catalog import neutral

def get_all_file_list(ixpe_path: str, chatter=False) -> list[str]:
    """
    Get a list of all files for the given list of observations.

    obs_list: list[str]
        List of 9-digit Observation IDs, as strings.

    ixpe_path: str
        Path to the directory containing IXPE observation files.
    """
    file_list = []
    obs_list = [d.name for d in Path(ixpe_path).iterdir() if d.is_dir() and d.name[0].isdigit()]
    for obs in obs_list:
        obs_dir = Path(ixpe_path) / obs
        if not obs_dir.is_dir():
            warnings.warn(f'Observation directory not found for obs: {obs} at path: {obs_dir}. Skipping this observation.')
            continue
        for fits_file in sorted(glob.glob(f"{obs_dir}/event_l2/*evt2_v[0-9][0-9]_rej.fits", recursive=True)):
            file_list.append(fits_file)
        
        # special cases
        if obs == '01002901_picorr':
            for fits_file in sorted(glob.glob(f"{obs_dir}/event_l2/*evt2_v[0-9][0-9]_picorr.fits", recursive=True)):
                file_list.append(fits_file)
        if obs == '03002599':
            for fits_file in sorted(glob.glob(f"{obs_dir}/event_l2/*evt2_v[0-9][0-9].fits", recursive=True)):
                file_list.append(fits_file)

    if len(file_list) == 0:
        raise ValueError(f'No files found for obs_list: {obs_list} in ixpe_path: {ixpe_path}')
    
    if chatter:
        print(f'Found {len(file_list)} files for obs_list: {obs_list} in ixpe_path: {ixpe_path}')
        for f in file_list:
            print(f)
    
    return file_list

def _get_obs_file_list(obs: str, file_list: list[str], chatter=False) -> list[str]:
    """
    Filter the file list for files matching the given observation ID.

    obs: str
        9-digit Observation ID, as a string.

    file_list: list[str]
        List of file paths to filter.
    """
    obs_file_list = [f for f in file_list if obs in f]
    if len(obs_file_list) == 0:
        raise ValueError(f'No files found in file_list matching obs: {obs}')
    if chatter:
        print(f'Found {len(obs_file_list)} files for observation {obs}:')
        for f in obs_file_list:
            print(f)
    return obs_file_list

def _stage(file_list: list[str], outdir: str) -> list[str]:
    """Symlink each level-2 file into *outdir* and return the links.

    ixpeobssim writes each product beside the file it reads (the input path
    with a suffix, ``ixpeobssim.utils.os_.check_output_file``), so products made
    from the links land in *outdir* and the observation database is only read.
    """
    os.makedirs(outdir, exist_ok=True)
    links = []
    for f in file_list:
        source = Path(f).resolve()
        link = Path(outdir) / Path(f).name
        if link.is_symlink() and link.resolve() == source:
            pass
        elif link.exists() or link.is_symlink():
            raise FileExistsError(f'{link} exists and is not a link to {source}')
        else:
            os.symlink(source, link)
        links.append(str(link))
    return links


def pcube(obs: str,
          file_list: list[str],
          ebinning: list[float] | None=None,
          src_rad: float | None=None,
          bkg_inner_rad: float | None=None,
          bkg_outer_rad: float | None=None,
          ra: float | None=None,
          dec: float | None=None,
          irf_date: str | None=None,
          changes: bool=True,
          overwrite: bool=True,
          chatter: bool=False,
          *,
          outdir: str,
          products: list | None=None,
          ) -> tuple[xBinnedPolarizationCube, xBinnedPolarizationCube]:
    """
    Reduce to IXPE PCUBE data for a given observation and args. returns cleaned src PCUBE, and corresponding bkg PCUBE.

    A PCUBE depends on the chosen energy binning and regions, so it is a
    science artifact: every file this writes (links to the level-2 inputs, the
    selected event lists, the PCUBEs) goes under *outdir*, never beside the
    inputs in the observation database.

    When changes=False: only obs, file_list and outdir are required; the PCUBEs an earlier call wrote in outdir are read.
    When changes=True: ebinning, src_rad, bkg_inner_rad, bkg_outer_rad, ra, dec are required.

    Parameters
    ----------
    obs: str
        9-digit Observation ID, as a string.

    file_list: list[str]
        List of file paths to filter for the given observation.

    Required if changes=True:

    ebinning: list[float]
        Energy binning scheme to use for PCUBE reduction. Required if changes=True.

    src_rad: float
        Radius of the source region in arcmin. Required if changes=True.

    bkg_inner_rad: float
        Inner radius of the background annulus in arcmin. Required if changes=True.

    bkg_outer_rad: float
        Outer radius of the background annulus in arcmin. Required if changes=True.

    ra: float
        Right ascension of the source in degrees. Required if changes=True.

    dec: float
        Declination of the source in degrees. Required if changes=True.

    irf_date: str
        Date string for the IRF version to use, in the format 'YYYYMMDD'. Required if changes=True.
        Avaliable IRF versions at: <installed ixpeobssim>/caldb/ixpe/gpd/cpf/rmf/ (ixpeobssim.IXPEOBSSIM_CALDB)

    outdir: str
        Folder the call writes into (keyword only).

    products: list or None
        When given, extended with one ``{'path', 'caption'}`` per PCUBE file the
        call wrote (changes=True, with the selected event lists) or read
        (changes=False), for a caller that returns its products.

    Returns
    -------
    tuple[xBinnedPolarizationCube, xBinnedPolarizationCube]
        A tuple containing the source PCUBE and background PCUBE as xBinnedPolarizationCube objects.
    """

    # Filter files for the given observation
    obs_file_list = _get_obs_file_list(obs, file_list)
    if len(obs_file_list) == 0:
        raise ValueError(f'No files found in file_list matching obs: {obs}')
    if chatter:
        print(f'Found {len(obs_file_list)} files for observation {obs}:')
        for f in obs_file_list:
            print(f)
    if not changes:
        expected = [Path(outdir) / f"{Path(f).with_suffix('').name}_{role}_pcube{Path(f).suffix}"
                    for f in obs_file_list for role in ('src', 'bkg')]
        absent = [str(p) for p in expected if not p.is_file()]
        if absent:
            raise FileNotFoundError(f'changes=False reads PCUBEs an earlier call wrote in {outdir}; '
                                    f'absent: {", ".join(absent)}')
    obs_file_list = _stage(obs_file_list, outdir)

    if changes:
        # Validate that all required parameters are provided when changes=True
        if any(param is None for param in [ebinning, src_rad, bkg_inner_rad, bkg_outer_rad, ra, dec, irf_date]):
            raise ValueError('When changes=True, all parameters (ebinning, src_rad, bkg_inner_rad, bkg_outer_rad, ra, dec, irf_date) must be provided.')
        '''
        If changes were made to input params, redo reduction
        '''
        # Select source and background regions
        with open(os.devnull, 'w') as fnull:
            with redirect_stdout(fnull):
                selected = {'source': pipeline.xpselect(*obs_file_list , rad= src_rad , suffix='src', ra = ra, dec = dec, overwrite=overwrite),
                            'background': pipeline.xpselect(*obs_file_list, ra = ra, dec = dec , innerrad=bkg_inner_rad , rad=bkg_outer_rad , suffix='bkg', overwrite=overwrite)}
                kwargs = dict(algorithm='PCUBE', ebinalg='LIST',ebinning=ebinning , overwrite=overwrite, irfname=f'ixpe:obssim{irf_date}_alpha075:v13', weights=True)

                # PCUBE results
                src_pcube_file_list = pipeline.xpbin(*selected['source'], **kwargs)
                bkg_pcube_file_list = pipeline.xpbin(*selected['background'], **kwargs)
    else:
        '''
        Use existing reduction results
        '''
        selected = {}
        src_pcube_file_list = []
        bkg_pcube_file_list = []
        for i in range(len(obs_file_list)):
            src_pcube_file_list.append(f"{Path(obs_file_list[i]).with_suffix('')}_src_pcube{Path(obs_file_list[i]).suffix}")
            bkg_pcube_file_list.append(f"{Path(obs_file_list[i]).with_suffix('')}_bkg_pcube{Path(obs_file_list[i]).suffix}")

    if chatter:
            print(f'Selected {len(src_pcube_file_list)} source files and {len(bkg_pcube_file_list)} background files for observation {obs}.')
    if products is not None:
        products.extend({'path': str(f), 'caption': f'IXPE {obs} {role} polarization cube {Path(f).name}'}
                        for role, files in (('source', src_pcube_file_list), ('background', bkg_pcube_file_list))
                        for f in files)
        products.extend({'path': str(f), 'caption': f'IXPE {obs} selected {role} events {Path(f).name}'}
                        for role, files in selected.items() for f in files)

    src_pcube = xBinnedPolarizationCube.from_file_list(src_pcube_file_list)
    bkg_pcube = xBinnedPolarizationCube.from_file_list(bkg_pcube_file_list)

    bkg_pcube *= src_pcube.backscal() / bkg_pcube.backscal()
    src_pcube -= bkg_pcube

    return src_pcube, bkg_pcube

def pha(obs: str,
        file_list: list[str],
        src_rad: float | None=None,
        bkg_inner_rad: float | None=None,
        bkg_outer_rad: float | None=None,
        ra: float | None=None,
        dec: float | None=None,
        irf_date: str | None=None,
        changes: bool=True,
        overwrite: bool=True,
        chatter: bool=False
    ) -> tuple[list[str], list[str]]:
    """
    Reduce to IXPE PHA data for a given observation and args. returns a list of PHA file paths.

    The Stokes I/Q/U spectra are a standard reduction, so they are
    observation-database items: with changes=True they are written beside the
    level-2 inputs (``<obs_db_path>/<obsid>/event_l2/``, checked before anything
    is written) and every Stokes I source and background spectrum is recorded
    through ``obsdb_catalog.neutral.advance`` (kind ``pha``) with the radii,
    position and response set in ``params``. changes=True needs overwrite=True.

    When changes=False: only obs and file_list are required.
    When changes=True: src_rad, bkg_inner_rad, bkg_outer_rad, ra, dec, irf_date are required.

    Parameters
    ----------
    obs: str
        9-digit Observation ID, as a string.

    file_list: list[str]
        List of file paths to filter for the given observation.

    Required if changes=True:

    src_rad: float
        Radius of the source region in arcmin. Required if changes=True.

    bkg_inner_rad: float
        Inner radius of the background annulus in arcmin. Required if changes=True.

    bkg_outer_rad: float
        Outer radius of the background annulus in arcmin. Required if changes=True.

    ra: float
        Right ascension of the source in degrees. Required if changes=True.

    dec: float
        Declination of the source in degrees. Required if changes=True.

    irf_date: str
        Date string for the IRF version to use, in the format 'YYYYMMDD'. Required if changes=True.
        Avaliable IRF versions at: <installed ixpeobssim>/caldb/ixpe/gpd (ixpeobssim.IXPEOBSSIM_CALDB)

    Returns
    -------
    tuple[list[str], list[str]]
        A tuple containing two lists: the first list contains the file paths to the source PHA files, and the second list contains the file paths to the background PHA files.
    """

    # Filter files for the given observation
    obs_file_list = _get_obs_file_list(obs, file_list)
    if len(obs_file_list) == 0:
        raise ValueError(f'No files found in file_list matching obs: {obs}')

    if changes:
        # Validate that all required parameters are provided when changes=True
        if any(param is None for param in [src_rad, bkg_inner_rad, bkg_outer_rad, ra, dec, irf_date]):
            raise ValueError('When changes=True, all parameters (src_rad, bkg_inner_rad, bkg_outer_rad, ra, dec, irf_date) must be provided.')
        if not overwrite:
            # ixpeobssim keeps an existing output under overwrite=False and still returns its path,
            # so the catalog would record settings the file was not made with.
            raise ValueError('changes=True records the spectra with these settings; pass overwrite=True '
                             'so every spectrum is made with them')
        root = _obs_db_of(obs, obs_file_list)
        '''
        If changes were made to input params, redo reduction
        '''

        with open(os.devnull, 'w') as fnull:
            with redirect_stdout(fnull):
                # Select source and background regions (Stokes-agnostic event files)
                src_sel = pipeline.xpselect(*obs_file_list , rad= src_rad , suffix='src', ra = ra, dec = dec, overwrite=overwrite)
                bkg_sel = pipeline.xpselect(*obs_file_list, ra = ra, dec = dec , innerrad=bkg_inner_rad , rad=bkg_outer_rad , suffix='bkg', overwrite=overwrite)

                # Bin the selected events into Stokes I/Q/U PHA spectra.  PHA1 (I)
                # paths are returned; PHA1Q/PHA1U (Q/U) are written to disk with the
                # _pha1q/_pha1u suffixes that get_pha_q/get_pha_u expect.
                bin_kwargs = dict(overwrite=overwrite, irfname=f'ixpe:obssim{irf_date}_alpha075:v13', weights=True)
                src_pha_file_list = pipeline.xpbin(*src_sel, algorithm='PHA1', **bin_kwargs)
                bkg_pha_file_list = pipeline.xpbin(*bkg_sel, algorithm='PHA1', **bin_kwargs)
                pipeline.xpbin(*src_sel, algorithm='PHA1Q', **bin_kwargs)
                pipeline.xpbin(*bkg_sel, algorithm='PHA1Q', **bin_kwargs)
                pipeline.xpbin(*src_sel, algorithm='PHA1U', **bin_kwargs)
                pipeline.xpbin(*bkg_sel, algorithm='PHA1U', **bin_kwargs)
        _record_pha(root, obs, src_pha_file_list + bkg_pha_file_list, {
            'src_rad_arcmin': float(src_rad), 'bkg_inner_arcmin': float(bkg_inner_rad),
            'bkg_outer_arcmin': float(bkg_outer_rad), 'ra': float(ra), 'dec': float(dec),
            'irfname': bin_kwargs['irfname'], 'weights': True,
        })
    else:
        '''
        Use existing reduction results
        '''
        src_pha_file_list = []
        bkg_pha_file_list = []
        for i in range(len(obs_file_list)):
            src_pha_file_list.append(f"{Path(obs_file_list[i]).with_suffix('')}_src_pha1{Path(obs_file_list[i]).suffix}")
            bkg_pha_file_list.append(f"{Path(obs_file_list[i]).with_suffix('')}_bkg_pha1{Path(obs_file_list[i]).suffix}")   

    if chatter:
            print(f'IXPE PHA query results for obs: {obs}')
            print('Source PHA files:')
            for f in src_pha_file_list:
                print(f'Source PHA file: {f}')
            print('Background PHA files:')
            for f in bkg_pha_file_list:
                print(f'Background PHA file: {f}')

    return src_pha_file_list, bkg_pha_file_list

_PHA1 = re.compile(r'_det(?P<du>[123])_evt2_(?P<ver>v\d+)_(?P<filter>(?:rej_)?(?:src|bkg))_pha1\.fits$')
_EVT2 = re.compile(r'_det[123]_evt2_v\d+(?:_rej)?\.fits$')


def _obs_db_of(obs: str, obs_file_list: list[str]) -> Path:
    """The observation database the level-2 inputs lie in, read from the IXPE layout
    ``<obs_db_path>/<obs>/event_l2/ixpe<obs>_det<n>_evt2_v<NN>[_rej].fits``; raises
    ``ValueError`` for an input outside it (a working copy such as ``<obs>_picorr/``), before
    anything is written."""
    for f in obs_file_list:
        path = Path(f)
        if path.parent.name != 'event_l2' or path.parent.parent.name != obs or not _EVT2.search(path.name):
            raise ValueError(f'{f} is not <obs_db_path>/{obs}/event_l2/ixpe{obs}_det<n>_evt2_v<NN>[_rej].fits; '
                             'pha writes its spectra into an observation database and records them in '
                             'its catalog, so it reduces only files in that layout')
    return Path(obs_file_list[0]).parent.parent.parent


def _record_pha(root: Path, obs: str, pha_file_list: list[str], params: dict) -> None:
    """Record every Stokes I spectrum :func:`pha` wrote (source and background per detector
    unit, kind ``pha``) in the catalog of the observation database *root*, each with the
    reduction's *params*; the Q/U spectra follow each I file's name (``_pha1q``/``_pha1u``)."""
    products = []
    for f in pha_file_list:
        match = _PHA1.search(Path(f).name)
        if match is None:
            raise ValueError(f'{f} is not an ixpeobssim PHA1 name; the spectra cannot be recorded')
        products.append({
            'kind': 'pha',
            'path': str(Path(f).relative_to(root)),
            'produced_by': 'xspec_utils.lib.red_ixpe.pha',
            'params': {
                'du': f"DU{match.group('du')}", 'stokes': 'I', 'filter': match.group('filter'),
                'calibration_version': match.group('ver'), **params,
                'note': 'the Stokes Q/U spectra follow the same name with _pha1q/_pha1u',
            },
        })
    neutral.advance(str(root), 'ixpe', obs, level=3, products=products)


def get_pha_q(pha_path: str) -> str:
    """
    Convert PHA1 file path to corresponding PHA1Q file path.
    """
    path = Path(pha_path)
    name = path.name.replace("pha1", "pha1q", 1)
    return str(path.with_name(name))

def get_pha_u(pha_path: str) -> str:
    """
    Convert PHA1 file path to corresponding PHA1U file path.
    """
    path = Path(pha_path)
    name = path.name.replace("pha1", "pha1u", 1)
    return str(path.with_name(name))

def get_mjd(obs, file_list) -> tuple[float, float, float]:
    """
    Get MJD of the observation from the header of the event file. returns (mjd_start, mjd_stop, mjd_mid).

    obs: str
        9-digit Observation ID, as a string.

    file_list: list[str]
        List of file paths to filter for the given observation.
    """

    # Filter files for the given observation
    obs_file_list = _get_obs_file_list(obs, file_list)
    if len(obs_file_list) == 0:
        raise ValueError(f'No files found in file_list matching obs: {obs}')
    
    with fits.open(obs_file_list[0]) as hdul:
        hdr = hdul["EVENTS"].header
        mjdref = hdr['MJDREFI'] + hdr['MJDREFF']
        timezero = hdr.get('TIMEZERO', 0.0)
        tstart, tstop = hdr["TSTART"]+timezero, hdr["TSTOP"]+timezero
        mjd_start = mjdref + tstart/86400.0
        mjd_stop = mjdref + tstop/86400.0
        mjd_mid = mjdref + 0.5*(tstart + tstop)/86400.0
    return mjd_start, mjd_stop, mjd_mid

def _sum_counts(obs, src_pha_file_list, bkg_pha_file_list, e_bin: tuple[float, float]) -> tuple[float, float]:
    """
    sum count rate in energy bin for a given observation, returns [count_rate, error].

    obs: str
        9-digit Observation ID, as a string.

    src_pha_file_list: list[str]
        List of paths to the source PHA files.

    bkg_pha_file_list: list[str]
        List of paths to the background PHA files.

    e_bin: tuple[float, float]
        Energy range (keV) for which to sum the counts.
    """

    # Filter files for the given observation
    src_pha_file_list = _get_obs_file_list(obs, src_pha_file_list)
    bkg_pha_file_list = _get_obs_file_list(obs, bkg_pha_file_list)
    if len(src_pha_file_list) == 0 or len(bkg_pha_file_list) == 0:
        raise ValueError(f'No files found in file_list matching obs: {obs}')
    
    # Clear XSPEC and warn user
    AllData.clear(); AllModels.clear()
    warnings.warn("XSPEC has been reset. All previous data and models have been cleared.")

    # XSPEC loads
    import xspec_help  # resolve_response / assign_responses
    spec_cnt = 0
    resolved = []
    with open(os.devnull, 'w') as fnull:
            with redirect_stdout(fnull):
                for src_pha, bkg_pha in zip(src_pha_file_list, bkg_pha_file_list):
                    spec_cnt += 1
                    AllData(f"{spec_cnt}:{spec_cnt} {src_pha}")
                    AllData(spec_cnt).background = bkg_pha
                    resolved += xspec_help.assign_responses(spec_cnt, src_pha)
    if resolved:
        print(xspec_help.responses_used(resolved))

    # obtain overall rate and error
    AllData.ignore(f"**-{e_bin[0]}  {e_bin[1]}-**")
    rate, var = 0.0, 0.0
    for i in range(1, spec_cnt + 1):
        spec = AllData(i)
        # Spectrum.rate[1] is the net rate's standard deviation (it equals the STAT_ERR sum in
        # quadrature; PyXspec's docstring calls it a variance), so it adds in quadrature.
        this_rate, this_err, *_ = spec.rate
        rate += this_rate
        var += this_err ** 2

    err = var**0.5

    return rate, err

def HR(obs, src_pha_file_list, bkg_pha_file_list, s_bin: tuple[float, float], h_bin: tuple[float, float], chatter=False) -> tuple[float, float]:
    """
    Calculate the hardness ratio for a given observation and PHA file list. returns [HR, HR_error].

    obs: str
        9-digit Observation ID, as a string.

    src_pha_file_list: list[str]
        List of paths to the source PHA files.

    bkg_pha_file_list: list[str]
        List of paths to the background PHA files.

    s_bin: tuple[float, float]
        Energy range for the soft band.

    h_bin: tuple[float, float]
        Energy range for the hard band.
    """

    # Filter files for the given observation
    src_pha_file_list = _get_obs_file_list(obs, src_pha_file_list)
    bkg_pha_file_list = _get_obs_file_list(obs, bkg_pha_file_list)
    if len(src_pha_file_list) == 0 or len(bkg_pha_file_list) == 0:
        raise ValueError(f'No files found in file_list matching obs: {obs}')

    if chatter:
        print(f'Calculating hardness ratio for observation {obs} with {len(src_pha_file_list)} source files and {len(bkg_pha_file_list)} background files.')

    # Get rates and errors for soft and hard bands
    s_rate, s_err = _sum_counts(obs, src_pha_file_list, bkg_pha_file_list, s_bin)
    h_rate, h_err = _sum_counts(obs, src_pha_file_list, bkg_pha_file_list, h_bin)

    if chatter:
        print(f'Soft band rate: {s_rate:.4f} +/- {s_err:.4f}')
        print(f'Hard band rate: {h_rate:.4f} +/- {h_err:.4f}')

    # Calculate hardness ratio and error
    if s_rate <= 0 or h_rate <= 0:
        raise ValueError(f'Non-positive count rate for observation {obs}: soft rate = {s_rate}, hard rate = {h_rate}. Cannot calculate hardness ratio.')
    hr = h_rate/s_rate
    hr_err = (h_err**2 + (hr*s_err)**2)**0.5 / s_rate

    if chatter:
        print(f'Hardness Ratio: {hr:.4f} +/- {hr_err:.4f}')

    return hr, hr_err