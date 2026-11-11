
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


def get_all_file_list(obs_dir: str, obs: str, chatter=False) -> list[str]:
    """
    The level-2 event files of one IXPE observation, in its folder.

    obs_dir: str
        The observation's folder: its product staged from the store (``storelib.stage_in``), laid out
        as HEASARC's ObsID folder (``event_l2/...``).

    obs: str
        The observation's ID (``01002901_picorr`` and ``03002599`` read other files, below).
    """
    file_list = sorted(glob.glob(f"{obs_dir}/event_l2/*evt2_v[0-9][0-9]_rej.fits"))
    # special cases
    if obs == '01002901_picorr':
        file_list += sorted(glob.glob(f"{obs_dir}/event_l2/*evt2_v[0-9][0-9]_picorr.fits"))
    if obs == '03002599':
        file_list += sorted(glob.glob(f"{obs_dir}/event_l2/*evt2_v[0-9][0-9].fits"))

    if len(file_list) == 0:
        raise ValueError(f'No level-2 event files of {obs} in {obs_dir}/event_l2/')

    if chatter:
        print(f'Found {len(file_list)} files of {obs} in {obs_dir}:')
        for f in file_list:
            print(f)

    return file_list


def spectra(folder: str, obs: str) -> tuple[list[str], list[str]]:
    """The Stokes I source and background spectra of ``obs`` below ``folder`` (a staged product
    holding them, or a folder :func:`pha` wrote them into), DU1 to DU3. Where a
    detector unit has spectra of both the background-rejected events (``_rej_``) and the plain ones
    (ObsID ``03002599``), the rejected ones are read, as :func:`get_all_file_list` lists them first.
    Raises ``ValueError`` unless one spectrum is left per role and detector unit."""
    found = {}
    for path in sorted(Path(folder).rglob(f'ixpe{obs}_det*_pha1.fits')):
        match = _PHA1.search(path.name)
        if match:
            found.setdefault((match.group('filter').removeprefix('rej_'), match.group('du')), []).append(str(path))
    out = {}
    for role in ('src', 'bkg'):
        for du in '123':
            files = found.get((role, du), [])
            rejected = [f for f in files if f'_rej_{role}_pha1' in Path(f).name]
            files = rejected if len(files) > 1 and rejected else files
            if len(files) != 1:
                raise ValueError(f'expected one {role} PHA1 spectrum of {obs} DU{du} below {folder}, found {files}')
            out.setdefault(role, []).append(files[0])
    return out['src'], out['bkg']

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
        chatter: bool=False,
        *,
        outdir: str | None=None,
        products: list | None=None,
    ) -> tuple[list[str], list[str]]:
    """
    Reduce to IXPE PHA data for a given observation and args. returns a list of PHA file paths.

    Library code writes no store (only a manifest tool does): with changes=True
    the level-2 inputs are linked into ``outdir`` and the reduction runs there
    (ixpeobssim writes each product beside the file it reads), so the selected
    event lists and the source and background PHA1/PHA1Q/PHA1U spectra all land
    in ``outdir`` and the inputs are only read. A ``products`` list, when given,
    is extended with one ``{'path', 'caption'}`` per spectrum written. A tool
    wrapping this function is what would store the spectra as a product.
    changes=True needs overwrite=True.

    When changes=False: only obs and file_list are required; the spectra are read beside the inputs.
    When changes=True: src_rad, bkg_inner_rad, bkg_outer_rad, ra, dec, irf_date and outdir are required.

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

    outdir: str
        Folder the call writes into (keyword only; required if changes=True).

    products: list or None
        When given, extended with one ``{'path', 'caption'}`` per spectrum written (changes=True).

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
        if outdir is None:
            raise ValueError('changes=True writes the spectra into a folder of its own: pass outdir')
        for f in obs_file_list:
            if not _EVT2.search(Path(f).name):
                raise ValueError(f'{f} is not ixpe{obs}_det<n>_evt2_v<NN>[_rej].fits; pha reduces only level-2 '
                                 'event files of its observation')
        links = _stage(obs_file_list, outdir)
        bin_kwargs = dict(overwrite=overwrite, irfname=f'ixpe:obssim{irf_date}_alpha075:v13', weights=True)
        with open(os.devnull, 'w') as fnull:
            with redirect_stdout(fnull):
                # Select source and background regions (Stokes-agnostic event files)
                src_sel = pipeline.xpselect(*links, rad=src_rad, suffix='src', ra=ra, dec=dec, overwrite=overwrite)
                bkg_sel = pipeline.xpselect(*links, ra=ra, dec=dec, innerrad=bkg_inner_rad, rad=bkg_outer_rad,
                                            suffix='bkg', overwrite=overwrite)
                # Bin the selected events into Stokes I/Q/U PHA spectra.  PHA1 (I)
                # paths are returned; PHA1Q/PHA1U (Q/U) are written to disk with the
                # _pha1q/_pha1u suffixes that get_pha_q/get_pha_u expect.
                written = {}
                for algorithm in ('PHA1', 'PHA1Q', 'PHA1U'):
                    written[algorithm] = (pipeline.xpbin(*src_sel, algorithm=algorithm, **bin_kwargs),
                                          pipeline.xpbin(*bkg_sel, algorithm=algorithm, **bin_kwargs))
        src_pha_file_list, bkg_pha_file_list = (list(files) for files in written['PHA1'])
        if products is not None:
            products.extend({'path': str(f), 'caption': f'IXPE {obs} {role} Stokes {stokes} spectrum {Path(f).name}'}
                            for stokes, algorithm in (('I', 'PHA1'), ('Q', 'PHA1Q'), ('U', 'PHA1U'))
                            for role, files in zip(('source', 'background'), written[algorithm]) for f in files)
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