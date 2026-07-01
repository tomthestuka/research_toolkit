
import os
import csv
import glob
import re
import tempfile
import warnings
from pathlib import Path
from astropy.io import fits
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.binning.polarization import xBinnedPolarizationCube
from contextlib import redirect_stdout, redirect_stderr
from xspec import AllData, Spectrum, AllModels

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
    spec_cnt = 0
    with open(os.devnull, 'w') as fnull:
            with redirect_stdout(fnull):
                for src_pha, bkg_pha in zip(src_pha_file_list, bkg_pha_file_list):
                    spec_cnt += 1
                    AllData(f"{spec_cnt}:{spec_cnt} {src_pha}")
                    AllData(spec_cnt).background = bkg_pha

    # obtain overall rate and error
    AllData.ignore(f"**-{e_bin[0]}  {e_bin[1]}-**")
    rate, var = 0.0, 0.0
    for i in range(1, spec_cnt):
        spec = AllData(i)
        this_rate, this_var, *_ = spec.rate
        rate += this_rate
        var += this_var

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

def _validate_du_files(obs: str, obs_file_list: list[str]) -> list[str]:
    """Guarantee exactly one Level-2 event file per DU1/DU2/DU3, in DU order.

    Raises ``ValueError`` naming the offending file(s) on any filename that
    doesn't carry a recognisable DU number, a missing DU, or more than one
    file for the same DU (e.g. two calibration versions, or a ``_rej`` and
    plain variant) — callers must resolve the ambiguity rather than have one
    silently picked for them.
    """
    by_du: dict[int, list[str]] = {}
    for f in obs_file_list:
        m = re.search(r'_det([123])_evt2_', os.path.basename(f))
        if not m:
            raise ValueError(f'pcube_timeseg: could not parse a DU number (_det[123]_evt2_) from {f}')
        by_du.setdefault(int(m.group(1)), []).append(f)

    present = set(by_du)
    if present != {1, 2, 3}:
        raise ValueError(
            f'pcube_timeseg: expected Level-2 event files for DU1, DU2 and DU3 for obs {obs}; '
            f'found DUs {sorted(present)} in: {obs_file_list}')

    for n in (1, 2, 3):
        if len(by_du[n]) > 1:
            raise ValueError(
                f'pcube_timeseg: multiple Level-2 files for DU{n} (obs {obs}): '
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
    default off, appropriate for bright source-dominated targets. Segments
    with zero events, a failed xpbin, or zero GTI overlap are skipped
    (warned, no row). Writes a CSV of the per-segment + integrated rows next
    to the observation's event files.

    ``exposure_s`` and the segment's ``t_start_mjd``/``t_stop_mjd`` are
    computed from the observation's own GTI extension intersected with each
    segment's ``[t0, t1)`` window, scaled by the (per-observation constant)
    ``DEADC`` deadtime-correction keyword — mirroring the ``exposure = dt *
    deadc`` averaging ixpeobssim's own LC/MDPMAPCUBE binning uses for Level-2
    data. They are deliberately **not** read off the segment PCUBE FITS
    header: ``xpselect`` only rewrites TSTART/TSTOP/ONTIME/LIVETIME/DEADC
    when ``ltimeupdate=True``, and that option raises on Level-2 IXPE event
    files (no per-event LIVETIME column). Left alone, every segment's PCUBE
    header carries the identical, unmodified full-observation values — the
    same defect class as segments silently inheriting the full-observation
    GTI.

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

    Returns
    -------
    dict with 'timeseries' (list of per-segment dicts), 'integrated' (dict,
    same fields, segment_id=-1), 'csv_path' (absolute path written).

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
    timeseries = []
    for i in range(len(edges) - 1):
        t0, t1 = edges[i], edges[i + 1]
        try:
            spc = _reduce(f'seg{i:03d}', t0, t1)
        except (Exception, SystemExit) as exc:
            warnings.warn(f'pcube_timeseg: segment {i} [{t0:.1f},{t1:.1f}) failed ({exc}); skipping.')
            continue
        row = _pcube_row(spc)
        overlap = _gti_overlap(gti_start, gti_stop, t0, t1)
        counts = row.pop('_counts')
        if counts <= 0 or overlap is None:
            warnings.warn(f'pcube_timeseg: segment {i} [{t0:.1f},{t1:.1f}) has zero events/GTI overlap; skipping.')
            continue
        clip_start, clip_stop, overlap_s = overlap
        row['exposure_s'] = overlap_s * deadc
        if row['exposure_s'] <= 0:
            warnings.warn(f'pcube_timeseg: segment {i} [{t0:.1f},{t1:.1f}) has zero exposure; skipping.')
            continue
        row['segment_id'] = i
        row['t_start_mjd'], row['t_stop_mjd'] = met_to_mjd(clip_start), met_to_mjd(clip_stop)
        timeseries.append(row)
        if chatter:
            print(f"segment {i}: PD={row['PD_pct']:.3f}% PA={row['PA_deg']:.3f} deg "
                  f"exposure={row['exposure_s']:.1f}s")

    spc_int = _reduce('int', None, None)
    integrated = _pcube_row(spc_int)
    del integrated['_counts']
    int_start, int_stop, int_overlap_s = _gti_overlap(gti_start, gti_stop, tstart, tstop)
    integrated['exposure_s'] = int_overlap_s * deadc
    integrated['segment_id'] = -1
    integrated['t_start_mjd'], integrated['t_stop_mjd'] = met_to_mjd(int_start), met_to_mjd(int_stop)

    obs_dir = Path(obs_file_list[0]).resolve().parent.parent
    csv_path = str(obs_dir / f'pcube_timeseg_{obs}.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=_TIMESEG_FIELDS)
        writer.writeheader()
        for row in timeseries:
            writer.writerow(row)
        writer.writerow(integrated)

    if chatter:
        print(f'pcube_timeseg: wrote {len(timeseries)} segment(s) + integrated to {csv_path}')

    return {'timeseries': timeseries, 'integrated': integrated, 'csv_path': csv_path}