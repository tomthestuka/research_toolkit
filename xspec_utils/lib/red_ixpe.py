
import os
import csv
import glob
import json
import re
import tempfile
import warnings
from pathlib import Path
import numpy as np
from astropy.io import fits
import ixpeobssim
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.binning.polarization import xBinnedPolarizationCube
from contextlib import redirect_stdout, redirect_stderr

import find_src

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
          chatter: bool=False
          ) -> tuple[xBinnedPolarizationCube, xBinnedPolarizationCube]:
    """
    Reduce to IXPE PCUBE data for a given observation and args. returns cleaned src PCUBE, and corresponding bkg PCUBE.
    
    When changes=False: only obs and file_list are required.
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
        Avaliable IRF versions at: /grad/fzhou/.local/lib/python3.13/site-packages/ixpeobssim/caldb/ixpe/gpd/cpf/rmf/

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
                src_pcube_file_list = pipeline.xpselect(*obs_file_list , rad= src_rad , suffix='src', ra = ra, dec = dec, overwrite=overwrite)
                bkg_pcube_file_list = pipeline.xpselect(*obs_file_list, ra = ra, dec = dec , innerrad=bkg_inner_rad , rad=bkg_outer_rad , suffix='bkg', overwrite=overwrite)
                kwargs = dict(algorithm='PCUBE', ebinalg='LIST',ebinning=ebinning , overwrite=overwrite, irfname=f'ixpe:obssim{irf_date}_alpha075:v13', weights=True)

                # PCUBE results
                src_pcube_file_list = pipeline.xpbin(*src_pcube_file_list, **kwargs)
                bkg_pcube_file_list = pipeline.xpbin(*bkg_pcube_file_list, **kwargs)
    else:
        '''
        Use existing reduction results
        '''
        src_pcube_file_list = []
        bkg_pcube_file_list = []
        for i in range(len(obs_file_list)):
            src_pcube_file_list.append(f"{Path(obs_file_list[i]).with_suffix('')}_src_pcube{Path(obs_file_list[i]).suffix}")
            bkg_pcube_file_list.append(f"{Path(obs_file_list[i]).with_suffix('')}_bkg_pcube{Path(obs_file_list[i]).suffix}")

    if chatter:
            print(f'Selected {len(src_pcube_file_list)} source files and {len(bkg_pcube_file_list)} background files for observation {obs}.')

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
        Avaliable IRF versions at: /grad/fzhou/.local/lib/python3.13/site-packages/ixpeobssim/caldb/ixpe/gpd

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

def reduce_pha(obs: str,
               file_list: list[str],
               src_rad: float = 1.0,
               bkg_inner_rad: float = 3.0,
               bkg_outer_rad: float = 5.0,
               irf_date: str = '20230702',
               overwrite: bool = True,
               chatter: bool = False,
               ) -> dict:
    """Reuse an observation's on-disk per-DU PHA1 (Stokes I) products, or derive them.

    Verifies exactly one Level-2 event file per DU1/DU2/DU3 (`_validate_du_files`),
    then checks whether the src+bkg PHA1 files `pha(..., changes=False)` would
    point at already exist on disk. If so they are reused as-is -- no reduction,
    no source centroiding. Otherwise the source position is centroided from the
    DU1 event file via `find_src.run_extraction` (writing `src.reg`/`bkg.reg`
    into the observation's own directory, one level above `event_l2/`), and
    `pha()` is run (`xpselect` + `xpbin PHA1` per DU) with a circular source
    region of radius `src_rad` and an annular background
    (`bkg_inner_rad`-`bkg_outer_rad`) centred there, at `irf_date`.

    Parameters
    ----------
    obs: str
        9-digit Observation ID, as a string.

    file_list: list[str]
        List of file paths to filter for the given observation.

    src_rad, bkg_inner_rad, bkg_outer_rad: float
        Source region radius / background annulus radii, in arcminutes.
        Used only when a fresh reduction is required.

    irf_date: str
        IRF date string ('YYYYMMDD'); resolves to 'ixpe:obssim{irf_date}_alpha075:v13'.
        Used only when a fresh reduction is required.

    overwrite: bool
        Passed through to `pha()`/`xpselect`/`xpbin` when a fresh reduction is run.

    Returns
    -------
    dict with 'src_pha', 'bkg_pha' (each a 3-element list, DU1/DU2/DU3 Stokes-I
    PHA1 paths), 'ra', 'dec' (None when reused), 'reused' (bool), 'src_reg',
    'bkg_reg' (paths written when a fresh reduction ran; None when reused).

    Raises
    ------
    ValueError
        `file_list` does not resolve to exactly one Level-2 event file per
        DU1/DU2/DU3 for `obs`, or the centroided source region cannot be parsed.
    """
    obs_file_list = _get_obs_file_list(obs, file_list, chatter=chatter)
    obs_file_list = _validate_du_files(obs, obs_file_list, caller='reduce_pha')

    src_pha, bkg_pha = pha(obs, file_list, changes=False, chatter=chatter)
    if all(os.path.exists(f) for f in src_pha + bkg_pha):
        if chatter:
            print(f'reduce_pha: reusing existing PHA1 products for {obs}')
        return {'src_pha': src_pha, 'bkg_pha': bkg_pha, 'ra': None, 'dec': None,
                'reused': True, 'src_reg': None, 'bkg_reg': None}

    obs_dir = Path(obs_file_list[0]).resolve().parent.parent
    src_reg = str(obs_dir / 'src.reg')
    bkg_reg = str(obs_dir / 'bkg.reg')
    with open(os.devnull, 'w') as fnull:
        with redirect_stdout(fnull):
            find_src.run_extraction(obs_file_list[0], src_reg, bkg_reg, None,
                                     src_rad, bkg_inner_rad, bkg_outer_rad)
    with open(src_reg) as f:
        region_line = f.readlines()[1]
    m = re.match(r'circle\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,', region_line)
    if not m:
        raise ValueError(f'reduce_pha: could not parse source region from {src_reg}: {region_line!r}')
    ra, dec = float(m.group(1)), float(m.group(2))

    src_pha, bkg_pha = pha(obs, file_list, src_rad=src_rad, bkg_inner_rad=bkg_inner_rad,
                            bkg_outer_rad=bkg_outer_rad, ra=ra, dec=dec, irf_date=irf_date,
                            changes=True, overwrite=overwrite, chatter=chatter)
    if chatter:
        print(f'reduce_pha: reduced fresh PHA1 products for {obs} at RA={ra:.6f}, Dec={dec:.6f}')
    return {'src_pha': src_pha, 'bkg_pha': bkg_pha, 'ra': ra, 'dec': dec,
            'reused': False, 'src_reg': src_reg, 'bkg_reg': bkg_reg}

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

def _parse_du(filename: str) -> int:
    """Parse the DU number (1/2/3) out of a Level-2 event filename's ``_det[123]_evt2_`` tag."""
    m = re.search(r'_det([123])_evt2_', os.path.basename(filename))
    if not m:
        raise ValueError(f'HR: could not parse a DU number (_det[123]_evt2_) from {filename}')
    return int(m.group(1))

def _installed_du_rmf(du: int) -> str:
    """Path to an installed CALDB RMF for DU ``du`` in the currently-active
    ixpeobssim installation (any IRF date). IXPE's EBOUNDS channel-to-energy
    grid is a fixed instrumental characteristic shared across every IRF date
    and DU -- only the RESPONSE matrix itself varies with date/DU -- so any
    installed RMF for this DU carries the same channel grid a PHA1's own
    (possibly stale) RESPFILE would.
    """
    rmf_dir = os.path.join(os.path.dirname(os.path.abspath(ixpeobssim.__file__)),
                            'caldb', 'ixpe', 'gpd', 'cpf', 'rmf')
    matches = sorted(glob.glob(os.path.join(rmf_dir, f'ixpe_d{du}_obssim*_alpha075_v0*.rmf')))
    if not matches:
        raise ValueError(f'HR: no installed CALDB RMF found for DU{du} under {rmf_dir}')
    return matches[0]

def _resolve_rmf_ebounds(pha_path: str, du: int) -> tuple[np.ndarray, np.ndarray]:
    """Resolve the (E_MIN, E_MAX) EBOUNDS arrays (keV) for ``pha_path``.

    Primary path: ``pha_path``'s own RESPFILE header keyword. Fallback (an
    already-reduced PHA1 whose RESPFILE was baked in under a reduction
    environment that no longer exists on this machine): an installed CALDB
    RMF for the same DU in the currently-active ixpeobssim installation --
    the identical EBOUNDS grid, since it is a fixed instrumental
    characteristic invariant across IRF date/DU (never hardcoded here --
    E_MIN/E_MAX are always read live from whichever RMF resolves). Raises
    ``ValueError`` if the fallback grid's channel count does not match the
    PHA's own channel count (a real mismatch, never silently accepted), or
    if neither the literal RESPFILE nor any installed DU RMF resolves.
    """
    with fits.open(pha_path) as hdul:
        spec = hdul['SPECTRUM']
        rmf = spec.header['RESPFILE']
        if not os.path.isabs(rmf):
            rmf = os.path.join(os.path.dirname(os.path.abspath(pha_path)), rmf)
        n_channels = len(spec.data)

    if not os.path.exists(rmf):
        rmf = _installed_du_rmf(du)

    with fits.open(rmf) as rmf_hdul:
        e_lo = rmf_hdul['EBOUNDS'].data['E_MIN'].astype(float)
        e_hi = rmf_hdul['EBOUNDS'].data['E_MAX'].astype(float)

    if len(e_lo) != n_channels:
        raise ValueError(
            f'HR: fallback RMF {rmf} has {len(e_lo)} EBOUNDS channels but '
            f'{pha_path} (DU{du}) has {n_channels} PHA channels -- refusing a mismatched grid.')

    return e_lo, e_hi

def _ebounds_channel_mask(e_lo: np.ndarray, e_hi: np.ndarray, band: tuple[float, float]) -> np.ndarray:
    """Boolean channel mask selecting EBOUNDS bins whose midpoint falls in ``[lo, hi)`` keV."""
    e_mid = 0.5 * (e_lo + e_hi)
    lo, hi = band
    return (e_mid >= lo) & (e_mid < hi)

def _pha_band_rate(pha_path: str, band: tuple[float, float], du: int) -> tuple[float, float]:
    """Sum one PHA1 file's own RATE/STAT_ERR columns over the channels its
    (resolved, see ``_resolve_rmf_ebounds``) RMF's EBOUNDS places inside
    ``band`` (keV). RATE is already in cts/s and summed directly; STAT_ERR is
    summed in quadrature over the selected channels.
    """
    with fits.open(pha_path) as hdul:
        spec = hdul['SPECTRUM']
        rate = spec.data['RATE'].astype(float)
        stat_err = spec.data['STAT_ERR'].astype(float)
    e_lo, e_hi = _resolve_rmf_ebounds(pha_path, du)
    mask = _ebounds_channel_mask(e_lo, e_hi, band)
    return float(np.nansum(rate[mask])), float(np.sqrt(np.nansum(stat_err[mask] ** 2)))

def _pha_backscal(pha_path: str) -> float:
    with fits.open(pha_path) as hdul:
        return float(hdul['SPECTRUM'].header['BACKSCAL'])

def _du_net_band_rate(src_pha: str, bkg_pha: str, band: tuple[float, float], bkg: bool, du: int) -> tuple[float, float]:
    """Net (optionally background-subtracted) band rate and 1-sigma error for one DU.

    When ``bkg=True``, the background PHA1's rate is scaled by
    ``src_BACKSCAL/bkg_BACKSCAL`` before subtraction, and its (scaled) error is added
    in quadrature. When ``bkg=False`` the background file is never read.
    """
    s_rate, s_err = _pha_band_rate(src_pha, band, du)
    if not bkg:
        return s_rate, s_err
    b_rate, b_err = _pha_band_rate(bkg_pha, band, du)
    scale = _pha_backscal(src_pha) / _pha_backscal(bkg_pha)
    net = s_rate - scale * b_rate
    err = (s_err ** 2 + (scale * b_err) ** 2) ** 0.5
    return net, err

def hr_from_bands(hard_rate: float, hard_err: float, soft_rate: float, soft_err: float) -> tuple[float, float]:
    """Community-normalized hardness ratio ``HR=(H-S)/(H+S)``, bounded in ``[-1, +1]``,
    and its closed-form first-order Gaussian-propagated 1-sigma error
    ``sigma_HR = (2/(H+S)^2) * sqrt(S^2*sigma_H^2 + H^2*sigma_S^2)`` for independent
    net band rates H, S carrying 1-sigma errors sigma_H, sigma_S.

    Returns ``(nan, nan)`` when ``hard_rate + soft_rate <= 0`` -- the ratio is
    undefined for a non-positive net band-rate sum, so this returns an honest NaN
    rather than an invented division result.
    """
    denom = hard_rate + soft_rate
    if denom <= 0:
        return float('nan'), float('nan')
    hr = (hard_rate - soft_rate) / denom
    hr_err = (2.0 / denom ** 2) * (soft_rate ** 2 * hard_err ** 2 + hard_rate ** 2 * soft_err ** 2) ** 0.5
    return hr, hr_err

def HR(obs: str,
       file_list: list[str],
       soft_band: tuple[float, float] = (2.0, 4.0),
       hard_band: tuple[float, float] = (4.0, 8.0),
       bkg: bool = False,
       chatter: bool = False) -> dict:
    """Net Stokes-I hardness ratio ``HR=(H-S)/(H+S)`` for one IXPE observation, read
    directly from its existing per-DU PHA1 (Stokes I) src(+bkg) products -- the same
    PHA1 geometry ``pha()``/``reduce_pha`` produce.

    Attempts DU1, DU2 and DU3 independently via ``pha(..., changes=False)``'s existing
    file-path locators: whichever DUs' Level-2 event files are present in
    ``file_list`` contribute (IXPE can genuinely fly with a DU off -- this never
    hard-fails on a truly absent DU), and every DU that IS present must have its PHA1
    product on disk (a present-but-unreduced DU raises rather than being silently
    dropped). Net per-band rate sums over the contributing DUs; errors add in
    quadrature over those DUs. Each PHA1 file's own RESPFILE RMF EBOUNDS resolves
    ``soft_band``/``hard_band`` to channels -- never a hardcoded channel grid; if a
    PHA1's RESPFILE predates a reduction-environment migration and no longer
    resolves on disk, ``_resolve_rmf_ebounds`` falls back to an installed CALDB RMF
    for the same DU (the EBOUNDS grid is a fixed instrumental characteristic shared
    across every IRF date/DU, so this is the identical grid, not an approximation).

    Background: ``bkg=True`` subtracts each DU's background PHA1 scaled by
    ``src_BACKSCAL/bkg_BACKSCAL``; ``bkg=False`` (default) uses the source PHA1 alone
    and never reads the background file.

    Parameters
    ----------
    obs: str
        9-digit Observation ID, as a string.

    file_list: list[str]
        Level-2 event file paths to filter for ``obs`` (see ``_get_obs_file_list``).
        ``pha(obs, file_list, changes=False)``'s locators derive this observation's
        existing per-DU PHA1 src/bkg paths from these same file paths, so passing
        fewer than 3 DUs' event files here is how a caller expresses "DU off".

    soft_band, hard_band: tuple[float, float]
        Soft/hard energy bands in keV. Defaults (2.0, 4.0) / (4.0, 8.0).

    bkg: bool
        Apply background subtraction per DU. Default False.

    Returns
    -------
    dict with 'HR', 'HR_err' (nan/nan if soft_rate+hard_rate <= 0), 'soft_rate',
    'soft_err', 'hard_rate', 'hard_err' (net, DU-summed cts/s and 1-sigma error), and
    'dus_used' (sorted list of the DU numbers that contributed).

    Raises
    ------
    ValueError
        A DU present in ``file_list`` has no PHA1 product on disk (run
        ``ixpe_reduce_pha``/``reduce_pha`` first); no DU contributes at all; or
        ``file_list`` carries more than one Level-2 event file for the same DU
        (e.g. two processing versions, or a mixed plain+``_rej`` set) -- this is
        never resolved automatically (no "preferred version" is picked), since
        silently picking one would double-count or silently drop data; the
        caller must pass exactly one event file per DU.
    """
    obs_file_list = _get_obs_file_list(obs, file_list, chatter=chatter)
    du_labels = [_parse_du(f) for f in obs_file_list]
    if not du_labels:
        raise ValueError(f'HR: no DU contributes for obs {obs}; file_list resolved to no files.')

    dups = {du for du in du_labels if du_labels.count(du) > 1}
    if dups:
        offending_files = [f for f, du in zip(obs_file_list, du_labels) if du in dups]
        raise ValueError(
            f'HR: obs {obs} has more than one Level-2 event file for DU(s) {sorted(dups)}: '
            f'{offending_files}. Pass exactly one event file per DU (one processing '
            f'version/variant).')

    src_pha_list, bkg_pha_list = pha(obs, file_list, changes=False, chatter=chatter)

    for du, s_pha, b_pha in zip(du_labels, src_pha_list, bkg_pha_list):
        if not os.path.exists(s_pha):
            raise ValueError(
                f'HR: DU{du} is present in file_list for obs {obs} but its source PHA1 '
                f'product is missing: {s_pha}. Run ixpe_reduce_pha/reduce_pha first.')
        if bkg and not os.path.exists(b_pha):
            raise ValueError(
                f'HR: DU{du} is present in file_list for obs {obs} but its background '
                f'PHA1 product is missing: {b_pha} (bkg=True).')

    def _band_sum(band: tuple[float, float]) -> tuple[float, float]:
        total_rate, total_var = 0.0, 0.0
        for du, s_pha, b_pha in zip(du_labels, src_pha_list, bkg_pha_list):
            r, e = _du_net_band_rate(s_pha, b_pha, band, bkg, du)
            total_rate += r
            total_var += e ** 2
        return total_rate, total_var ** 0.5

    soft_rate, soft_err = _band_sum(soft_band)
    hard_rate, hard_err = _band_sum(hard_band)
    hr, hr_err = hr_from_bands(hard_rate, hard_err, soft_rate, soft_err)

    if chatter:
        print(f'HR {obs}: DUs {sorted(du_labels)}  soft={soft_rate:.4f}+/-{soft_err:.4f}  '
              f'hard={hard_rate:.4f}+/-{hard_err:.4f}  HR={hr:.4f}+/-{hr_err:.4f}')

    return {
        'HR': hr, 'HR_err': hr_err,
        'soft_rate': soft_rate, 'soft_err': soft_err,
        'hard_rate': hard_rate, 'hard_err': hard_err,
        'dus_used': sorted(du_labels),
    }

_TIMESEG_FIELDS = ['segment_id', 't_start_mjd', 't_stop_mjd', 'exposure_s',
                   'PD_pct', 'PD_err', 'PA_deg', 'PA_err',
                   'QI', 'QI_err', 'UI', 'UI_err', 'MDP99', 'SIGNIF']

def _segment_edges(tstart: float, tstop: float, duration: float) -> list[float]:
    """Natural elapsed-time segment edges: TSTART + k*duration, plus a final
    partial edge at TSTOP when the span does not divide evenly."""
    n_full = int((tstop - tstart) // duration)
    edges = [tstart + k * duration for k in range(n_full + 1)]
    if tstop - edges[-1] > 0:
        edges.append(tstop)
    return edges

def _gti_overlap(gti_start, gti_stop, t0: float, t1: float) -> tuple[float, float, float] | None:
    """Intersect the observation's GTI intervals with the MET window [t0, t1).

    Returns ``(clipped_start, clipped_stop, overlap_s)`` — the tightest MET
    bounds actually covered by data within the window and the total livetime
    of GTI overlap — or ``None`` if no GTI interval overlaps the window at all
    (zero exposure).
    """
    clipped_start = clipped_stop = None
    overlap_s = 0.0
    for s, e in zip(gti_start, gti_stop):
        lo, hi = max(s, t0), min(e, t1)
        if hi > lo:
            overlap_s += hi - lo
            clipped_start = lo if clipped_start is None else min(clipped_start, lo)
            clipped_stop = hi if clipped_stop is None else max(clipped_stop, hi)
    if clipped_start is None:
        return None
    return clipped_start, clipped_stop, overlap_s

def _validate_du_files(obs: str, obs_file_list: list[str], caller: str = 'pcube_timeseg') -> list[str]:
    """Guarantee exactly one Level-2 event file per DU1/DU2/DU3, in DU order.

    Raises ``ValueError`` naming the offending file(s) on any filename that
    doesn't carry a recognisable DU number, a missing DU, or more than one
    file for the same DU (e.g. two calibration versions, or a ``_rej`` and
    plain variant) — callers must resolve the ambiguity rather than have one
    silently picked for them. ``caller`` is prefixed onto every message so it
    names the actual entry point (e.g. ``'reduce_pha'``), not just its first user.
    """
    by_du: dict[int, list[str]] = {}
    for f in obs_file_list:
        m = re.search(r'_det([123])_evt2_', os.path.basename(f))
        if not m:
            raise ValueError(f'{caller}: could not parse a DU number (_det[123]_evt2_) from {f}')
        by_du.setdefault(int(m.group(1)), []).append(f)

    present = set(by_du)
    if present != {1, 2, 3}:
        raise ValueError(
            f'{caller}: expected Level-2 event files for DU1, DU2 and DU3 for obs {obs}; '
            f'found DUs {sorted(present)} in: {obs_file_list}')

    for n in (1, 2, 3):
        if len(by_du[n]) > 1:
            raise ValueError(
                f'{caller}: multiple Level-2 files for DU{n} (obs {obs}): '
                f'{by_du[n]}; pass exactly one file per DU')

    return [by_du[1][0], by_du[2][0], by_du[3][0]]

def _pcube_row(spc: xBinnedPolarizationCube) -> dict:
    """Extract the single-energy-bin PD/PA/Stokes/MDP99 fields from a PCUBE.

    Does NOT read ``exposure_s`` here: ``xpselect`` only rewrites
    TSTART/TSTOP/ONTIME/LIVETIME/DEADC when ``ltimeupdate=True``, and that
    option raises on Level-2 IXPE event files (no per-event LIVETIME column —
    ixpeobssim only supports an average deadtime correction for L2 data).
    Left alone, every segment's PCUBE would carry the identical, unmodified
    full-observation LIVETIME. The caller computes the real per-segment
    exposure independently from the GTI extension (see ``_gti_overlap``).
    """
    return {
        'PD_pct': float(spc.PD[0]) * 100.0,
        'PD_err': float(spc.PD_ERR[0]) * 100.0,
        'PA_deg': float(spc.PA[0]),
        'PA_err': float(spc.PA_ERR[0]),
        'QI': float(spc.QN[0]),
        'QI_err': float(spc.QN_ERR[0]),
        'UI': float(spc.UN[0]),
        'UI_err': float(spc.UN_ERR[0]),
        'MDP99': float(spc.MDP_99[0]) * 100.0,
        'SIGNIF': float(spc.SIGNIF[0]),
        '_counts': float(spc.COUNTS[0]),
    }

_PROVENANCE_HEADER_PREFIX = '# pcube_timeseg provenance: '

def _provenance_line(provenance: dict) -> str:
    """Serialise the run's reproducibility-relevant parameters as the CSV's first line."""
    return _PROVENANCE_HEADER_PREFIX + json.dumps(provenance, sort_keys=True) + '\n'

def _read_existing_timeseg_csv(csv_path: str) -> tuple[dict, dict, dict | None] | None:
    """Read a previously written ``pcube_timeseg`` CSV back into memory.

    Returns ``(provenance, rows_by_segment_id, integrated_row_or_None)``, or
    ``None`` if ``csv_path`` does not exist or does not carry a recognisable
    provenance header (the latter is treated as absent -- callers override
    rather than resume from an unrecognised file).
    """
    if not os.path.exists(csv_path):
        return None
    with open(csv_path, newline='') as f:
        header_line = f.readline()
        if not header_line.startswith(_PROVENANCE_HEADER_PREFIX):
            return None
        provenance = json.loads(header_line[len(_PROVENANCE_HEADER_PREFIX):])
        rows_by_segid: dict[int, dict] = {}
        integrated = None
        for r in csv.DictReader(f):
            row = {k: (int(v) if k == 'segment_id' else float(v)) for k, v in r.items()}
            if row['segment_id'] == -1:
                integrated = row
            else:
                rows_by_segid[row['segment_id']] = row
    return provenance, rows_by_segid, integrated

def _write_timeseg_csv_header(csv_path: str, provenance: dict) -> None:
    with open(csv_path, 'w', newline='') as f:
        f.write(_provenance_line(provenance))
        csv.DictWriter(f, fieldnames=_TIMESEG_FIELDS).writeheader()

def _append_timeseg_csv_row(csv_path: str, row: dict) -> None:
    """Append one row and durably flush it -- the checkpoint a resumed call reads back."""
    with open(csv_path, 'a', newline='') as f:
        csv.DictWriter(f, fieldnames=_TIMESEG_FIELDS).writerow(row)
        f.flush()
        os.fsync(f.fileno())

def pcube_timeseg(obs: str,
                   file_list: list[str],
                   irf_date: str,
                   ra: float | None = None,
                   dec: float | None = None,
                   emin_kev: float = 2.0,
                   emax_kev: float = 8.0,
                   src_rad: float = 1.0,
                   bkg: bool = False,
                   bkg_inner_rad: float = 3.0,
                   bkg_outer_rad: float = 5.0,
                   segment_duration_s: float = 5000.0,
                   min_seg_exposure_s: float = 1000.0,
                   resume: bool = True,
                   overwrite: bool = True,
                   chatter: bool = False,
                   ) -> dict:
    """Time-resolved IXPE PCUBE polarimetry for one observation, DU1/2/3 merged.

    Splits the observation into natural elapsed-time segments of
    ``segment_duration_s`` (TSTART + k*segment_duration_s, NOT exposure
    equalised) and runs a single-band (``emin_kev``-``emax_kev``) weighted
    PCUBE reduction in each, plus one integrated whole-observation PCUBE.
    Background subtraction (annulus ``bkg_inner_rad``-``bkg_outer_rad``) is
    applied per segment and to the integrated cube only when ``bkg=True`` —
    default off, appropriate for bright source-dominated targets.

    **GTI pre-reduction gate.** Each window's GTI-clipped good exposure is
    computed via ``_gti_overlap`` *before* any ``xpselect``/``xpbin`` call;
    windows whose good exposure is below ``min_seg_exposure_s`` are skipped
    without running the reduction (counted in the returned
    ``n_skipped_empty``). This matters for observations with a long
    inter-pointing gap (e.g. two pointings glued together with a multi-day
    gap between them): without the gate, every window that falls inside the
    gap still runs a full reduction for a spurious near-zero-exposure row.
    A segment that clears the gate but still yields zero counts or a failed
    reduction is warned and skipped (not counted in ``n_skipped_empty``,
    since the gate itself passed it as an expected segment).

    **Incremental write + resume.** Rows are appended to
    ``pcube_timeseg_<obs>.csv`` and flushed as each segment's reduction
    completes, so a call that hits the sandbox's wall-clock cap leaves a
    genuine checkpoint of the segments already done. The CSV's first line is
    a provenance comment recording every parameter that determines the
    segmentation and reduction (``segment_duration_s``, ``emin_kev``,
    ``emax_kev``, ``src_rad``, ``bkg``, ``bkg_inner_rad``, ``bkg_outer_rad``,
    ``ra``, ``dec``, ``irf_date``, ``min_seg_exposure_s``, ``tstart``,
    ``tstop``). When ``resume=True`` and an existing CSV's provenance matches
    this call's exactly, only the missing segments (and the integrated row,
    if absent) are computed; existing rows are re-emitted from disk without
    a repeat reduction. Otherwise (``resume=False``, no existing CSV, or any
    provenance mismatch) the CSV is deleted and rebuilt from scratch — rows
    from a different configuration are never glued onto a fresh run.

    Parameters
    ----------
    obs, file_list : see ``_get_obs_file_list``.
    irf_date : str
        IRF date string ('YYYYMMDD'); resolves to ``ixpe:obssim{irf_date}_alpha075:v13``.
    ra, dec : float or None
        Source position in degrees. If both None, centroided from one Level-2
        event file via ``find_src.run_extraction``.
    emin_kev, emax_kev : float
        Single energy bin edges in keV. Default 2.0-8.0.
    src_rad : float
        Source region radius in arcmin. Default 1.0.
    bkg : bool
        Subtract a background annulus (``bkg_inner_rad``-``bkg_outer_rad``
        arcmin). Default False.
    segment_duration_s : float
        Elapsed-time segment length in seconds. Default 5000.0.
    min_seg_exposure_s : float
        Minimum GTI-clipped good exposure (seconds) a window must have to be
        reduced at all. Default 1000.0 -- above the sub-minute slices a
        window can straddle at a pointing-gap or GTI boundary, below a
        genuine multi-ks segment.
    resume : bool
        Resume from an existing same-provenance CSV, computing only the
        missing segments. Default True.

    Returns
    -------
    dict with 'timeseries' (list of per-segment dicts, including any rows
    re-emitted from a resumed CSV), 'integrated' (dict, same fields,
    segment_id=-1), 'csv_path' (absolute path), 'n_skipped_empty' (windows
    the GTI gate skipped), 'n_expected' (windows the gate passed, i.e.
    attempted), 'complete' (True iff every expected segment plus the
    integrated row is present), 'resumed' (True iff this call resumed from a
    matching on-disk CSV rather than starting fresh).

    Raises
    ------
    ValueError
        ``file_list`` does not resolve to exactly one Level-2 event file per
        DU1/DU2/DU3 for ``obs`` (missing DU, unparseable filename, or more
        than one file for the same DU — never silently picked), or exactly
        one of ``ra``/``dec`` is given without the other.
    """
    obs_file_list = _get_obs_file_list(obs, file_list, chatter=chatter)
    obs_file_list = _validate_du_files(obs, obs_file_list)

    if (ra is None) != (dec is None):
        raise ValueError(
            f'pcube_timeseg: provide both ra and dec, or neither (auto-centroid); '
            f'got ra={ra}, dec={dec}')

    if ra is None:
        with tempfile.TemporaryDirectory(prefix=f'pcube_timeseg_{obs}_') as tmp_dir:
            src_reg = os.path.join(tmp_dir, 'src.reg')
            with open(os.devnull, 'w') as fnull:
                with redirect_stdout(fnull):
                    find_src.run_extraction(
                        obs_file_list[0], src_reg,
                        os.path.join(tmp_dir, 'bkg.reg'), None,
                        src_rad, bkg_inner_rad, bkg_outer_rad)
            with open(src_reg) as f:
                region_line = f.readlines()[1]
            m = re.match(r'circle\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,', region_line)
            if not m:
                raise ValueError(f'Could not parse source region from {src_reg}: {region_line!r}')
            ra, dec = float(m.group(1)), float(m.group(2))

    with fits.open(obs_file_list[0]) as hdul:
        hdr = hdul['EVENTS'].header
        if 'GTI' not in hdul or len(hdul['GTI'].data) == 0:
            raise ValueError(f'{obs_file_list[0]}: no GTI extension data; cannot compute exposure_s')
        gti = hdul['GTI'].data
        gti_start = [float(x) for x in gti['START']]
        gti_stop = [float(x) for x in gti['STOP']]
        tstart, tstop = hdr.get('TSTART'), hdr.get('TSTOP')
        if tstart is None or tstop is None:
            tstart, tstop = gti_start[0], gti_stop[-1]
        tstart, tstop = float(tstart), float(tstop)
        mjdref = float(hdr['MJDREFI']) + float(hdr['MJDREFF'])
        deadc = hdul['PRIMARY'].header.get('DEADC')
        if deadc is None:
            raise ValueError(f'{obs_file_list[0]}: PRIMARY header missing DEADC keyword')
        deadc = float(deadc)

    def met_to_mjd(met: float) -> float:
        return mjdref + met / 86400.0

    bin_kwargs = dict(algorithm='PCUBE', ebinalg='LIST', ebinning=[emin_kev, emax_kev],
                       overwrite=overwrite, irfname=f'ixpe:obssim{irf_date}_alpha075:v13',
                       weights=True)

    def _reduce(suffix: str, tmin: float | None, tmax: float | None) -> xBinnedPolarizationCube:
        sel_kwargs = dict(rad=src_rad, ra=ra, dec=dec, suffix=f'{suffix}src', overwrite=overwrite)
        if tmin is not None:
            sel_kwargs.update(tmin=tmin, tmax=tmax)
        with open(os.devnull, 'w') as fnull:
            with redirect_stdout(fnull):
                src_files = pipeline.xpselect(*obs_file_list, **sel_kwargs)
                spc = xBinnedPolarizationCube.from_file_list(pipeline.xpbin(*src_files, **bin_kwargs))
                if bkg:
                    bkg_sel_kwargs = dict(ra=ra, dec=dec, innerrad=bkg_inner_rad,
                                           rad=bkg_outer_rad, suffix=f'{suffix}bkg', overwrite=overwrite)
                    if tmin is not None:
                        bkg_sel_kwargs.update(tmin=tmin, tmax=tmax)
                    bkg_files = pipeline.xpselect(*obs_file_list, **bkg_sel_kwargs)
                    bpc = xBinnedPolarizationCube.from_file_list(pipeline.xpbin(*bkg_files, **bin_kwargs))
                    bpc *= spc.backscal() / bpc.backscal()
                    spc -= bpc
        return spc

    edges = _segment_edges(tstart, tstop, segment_duration_s)

    provenance = {
        'segment_duration_s': segment_duration_s, 'emin_kev': emin_kev, 'emax_kev': emax_kev,
        'src_rad': src_rad, 'bkg': bkg, 'bkg_inner_rad': bkg_inner_rad,
        'bkg_outer_rad': bkg_outer_rad, 'ra': ra, 'dec': dec, 'irf_date': irf_date,
        'min_seg_exposure_s': min_seg_exposure_s, 'tstart': tstart, 'tstop': tstop,
    }

    obs_dir = Path(obs_file_list[0]).resolve().parent.parent
    csv_path = str(obs_dir / f'pcube_timeseg_{obs}.csv')

    existing = _read_existing_timeseg_csv(csv_path)
    resumed = bool(resume and existing is not None and existing[0] == provenance)
    if resumed:
        _, existing_rows, existing_integrated = existing
    else:
        if os.path.exists(csv_path):
            os.remove(csv_path)
        _write_timeseg_csv_header(csv_path, provenance)
        existing_rows, existing_integrated = {}, None

    timeseries = []
    n_skipped_empty = 0
    n_expected = 0
    for i in range(len(edges) - 1):
        t0, t1 = edges[i], edges[i + 1]
        overlap = _gti_overlap(gti_start, gti_stop, t0, t1)
        good_exposure = overlap[2] * deadc if overlap is not None else 0.0
        if overlap is None or good_exposure < min_seg_exposure_s:
            n_skipped_empty += 1
            continue
        n_expected += 1

        if i in existing_rows:
            timeseries.append(existing_rows[i])
            continue

        try:
            spc = _reduce(f'seg{i:03d}', t0, t1)
        except (Exception, SystemExit) as exc:
            warnings.warn(f'pcube_timeseg: segment {i} [{t0:.1f},{t1:.1f}) failed ({exc}); skipping.')
            continue
        row = _pcube_row(spc)
        counts = row.pop('_counts')
        if counts <= 0:
            warnings.warn(f'pcube_timeseg: segment {i} [{t0:.1f},{t1:.1f}) has zero events; skipping.')
            continue
        clip_start, clip_stop, _overlap_s = overlap
        row['exposure_s'] = good_exposure
        row['segment_id'] = i
        row['t_start_mjd'], row['t_stop_mjd'] = met_to_mjd(clip_start), met_to_mjd(clip_stop)
        timeseries.append(row)
        _append_timeseg_csv_row(csv_path, row)
        if chatter:
            print(f"segment {i}: PD={row['PD_pct']:.3f}% PA={row['PA_deg']:.3f} deg "
                  f"exposure={row['exposure_s']:.1f}s")

    if existing_integrated is not None:
        integrated = existing_integrated
    else:
        spc_int = _reduce('int', None, None)
        integrated = _pcube_row(spc_int)
        del integrated['_counts']
        int_overlap = _gti_overlap(gti_start, gti_stop, tstart, tstop)
        if int_overlap is None:
            raise ValueError(f'{obs_file_list[0]}: no GTI overlap over the full [tstart,tstop) span')
        int_start, int_stop, int_overlap_s = int_overlap
        integrated['exposure_s'] = int_overlap_s * deadc
        integrated['segment_id'] = -1
        integrated['t_start_mjd'], integrated['t_stop_mjd'] = met_to_mjd(int_start), met_to_mjd(int_stop)
        _append_timeseg_csv_row(csv_path, integrated)

    complete = len(timeseries) == n_expected and integrated is not None

    if chatter:
        print(f'pcube_timeseg: {len(timeseries)}/{n_expected} segment(s) + integrated -> {csv_path} '
              f'(skipped {n_skipped_empty}, resumed={resumed}, complete={complete})')

    return {
        'timeseries': timeseries, 'integrated': integrated, 'csv_path': csv_path,
        'n_skipped_empty': n_skipped_empty, 'n_expected': n_expected,
        'complete': complete, 'resumed': resumed,
    }