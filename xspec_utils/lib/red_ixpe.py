
import os
import csv
import glob
import io
import json
import math
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
import irf_epoch

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
    Filter the file list for files belonging to the given observation ID.

    Matches by each file's own observation directory -- the parent of its
    ``event_l2/`` directory, i.e. ``.../<obs>/event_l2/<file>``, the layout
    ``get_all_file_list`` itself builds -- rather than a raw substring test.
    A plain 9-digit ``obs`` must not also pick up a differently-processed
    variant stored under its own suffixed directory (e.g. ``01002901_picorr``,
    whose filenames still literally embed the same 9-digit ObsID): a
    substring filter matches both wherever ``file_list`` spans more than one
    ObsID's directory (as ``get_all_file_list``'s mission-wide file_list
    does), which then reads as a spurious "two files for the same DU"
    ambiguity. Requesting ``obs='01002901_picorr'`` explicitly still resolves
    to that variant's own files, unaffected.

    obs: str
        9-digit Observation ID, as a string (or a special-case directory name
        such as '01002901_picorr').

    file_list: list[str]
        List of file paths to filter.
    """
    obs_file_list = [f for f in file_list if Path(f).parent.parent.name == obs]
    if len(obs_file_list) == 0:
        raise ValueError(f'No files found in file_list matching obs: {obs}')
    if chatter:
        print(f'Found {len(obs_file_list)} files for observation {obs}:')
        for f in obs_file_list:
            print(f)
    return obs_file_list

def pcube(obs: str | list[str],
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
          ) -> tuple[xBinnedPolarizationCube, xBinnedPolarizationCube, str | None]:
    """
    Reduce to IXPE PCUBE data for one observation, or a native multi-ObsID
    combination, and args. Returns cleaned src PCUBE, corresponding bkg PCUBE,
    and the irf_date actually used.

    When changes=False: only obs and file_list are required.
    When changes=True: ebinning, src_rad, bkg_inner_rad, bkg_outer_rad, ra, dec are required.

    Parameters
    ----------
    obs: str | list[str]
        A 9-digit Observation ID, or a list of such IDs. When a list is given,
        every listed ObsID's Level-2 event files are filtered from file_list,
        concatenated (DU-validated per ObsID via _validate_du_files), and
        passed into ONE xpselect/xpbin PCUBE call -- the returned cube sums
        raw Stokes I/Q/U/COUNTS across every DU *and* every listed ObsID
        before PD/PA are derived, i.e. a native directly-combined reduction,
        not a post-hoc merge of separately-computed per-ObsID PD/PA.
        changes=False (reuse on-disk products) is only supported for a single
        ObsID -- see Raises.

    file_list: list[str]
        List of file paths to filter for the given observation(s).

    Required if changes=True:

    ebinning: list[float]
        Energy bin edges (keV) for PCUBE reduction (ebinalg='LIST'), >= 2
        edges (>= 1 band). Required if changes=True.

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

    irf_date: str | None
        IRF calibration epoch 'YYYYMMDD', resolving to 'ixpe:obssim{irf_date}_alpha075:v13'.
        Leave None (default) to auto-select the installed epoch nearest the observation
        date (irf_epoch.resolve_irf_date, midpoint of the earliest start to the latest
        stop across every listed ObsID); pass a value only to override.

    Returns
    -------
    tuple[xBinnedPolarizationCube, xBinnedPolarizationCube, str | None]
        The source PCUBE and background PCUBE as xBinnedPolarizationCube
        objects, and the irf_date used (resolved when changes=True; whatever
        was passed in, possibly None, when changes=False).

    Raises
    ------
    ValueError
        obs resolves to an empty list; any requested ObsID has no matching
        files in file_list; any requested ObsID is missing DU1/DU2/DU3 or
        carries more than one Level-2 file for a DU; changes=True and
        ebinning has fewer than 2 edges; or changes=False with more than one
        ObsID requested (reuse of on-disk products is not defined for a
        combined multi-ObsID request).
    """
    obs_list = [obs] if isinstance(obs, str) else list(obs)
    if not obs_list:
        raise ValueError('pcube: obs resolved to an empty list of ObsIDs')

    # Filter (and DU-validate) files per requested ObsID, then concatenate --
    # this is what makes the subsequent xpselect/xpbin calls below span every
    # listed ObsID's events in one shot.
    obs_file_list = []
    for o in obs_list:
        o_files = _get_obs_file_list(o, file_list, chatter=chatter)
        o_files = _validate_du_files(o, o_files, caller='pcube')
        obs_file_list.extend(o_files)
    label = obs if isinstance(obs, str) else '+'.join(obs_list)
    if chatter:
        print(f'Found {len(obs_file_list)} files for observation(s) {label}:')
        for f in obs_file_list:
            print(f)

    if changes:
        # irf_date=None -> auto-select the installed epoch nearest the observation(s).
        irf_date = _resolve_irf_date(obs, file_list, irf_date, chatter)
        # Validate that all required parameters are provided when changes=True
        if any(param is None for param in [ebinning, src_rad, bkg_inner_rad, bkg_outer_rad, ra, dec]):
            raise ValueError('When changes=True, all parameters (ebinning, src_rad, bkg_inner_rad, bkg_outer_rad, ra, dec) must be provided.')
        if len(ebinning) < 2:
            raise ValueError(f'pcube: ebinning must carry at least 2 edges (>= 1 band), got {ebinning}')
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
        if len(obs_list) > 1:
            raise ValueError('pcube: changes=False (reuse existing products) is not '
                              'supported for a combined multi-ObsID obs list.')
        '''
        Use existing reduction results
        '''
        src_pcube_file_list = []
        bkg_pcube_file_list = []
        for i in range(len(obs_file_list)):
            src_pcube_file_list.append(f"{Path(obs_file_list[i]).with_suffix('')}_src_pcube{Path(obs_file_list[i]).suffix}")
            bkg_pcube_file_list.append(f"{Path(obs_file_list[i]).with_suffix('')}_bkg_pcube{Path(obs_file_list[i]).suffix}")

    if chatter:
            print(f'Selected {len(src_pcube_file_list)} source files and {len(bkg_pcube_file_list)} background files for observation(s) {label}.')

    src_pcube = xBinnedPolarizationCube.from_file_list(src_pcube_file_list)
    bkg_pcube = xBinnedPolarizationCube.from_file_list(bkg_pcube_file_list)

    bkg_pcube *= src_pcube.backscal() / bkg_pcube.backscal()
    src_pcube -= bkg_pcube

    return src_pcube, bkg_pcube, irf_date


def pcube_bands(spc: xBinnedPolarizationCube) -> list[dict]:
    """Per-energy-band PD/PA/Stokes/MDP99/SIGNIF fields from a PCUBE, one dict
    per energy bin in spc's own grid (ENERG_LO/ENERG_HI) -- read via
    polarization() and direct attribute access, never hand-rolled Stokes
    arithmetic. Works for a PCUBE binned with any number of energy bins
    (ebinalg='LIST' with >= 2 edges), including the single-bin (broadband) case.

    Units: PD, PD_err, QN, QN_err, UN, UN_err, MDP99 are dimensionless
    fractions -- the same units xBinnedPolarizationCube itself stores them in
    (not percent). PA, PA_err are degrees (as stored). SIGNIF is dimensionless.
    """
    pd, pd_err, pa, pa_err = spc.polarization()
    return [
        {
            'E_LO': float(spc.ENERG_LO[i]), 'E_HI': float(spc.ENERG_HI[i]),
            'PD': float(pd[i]), 'PD_err': float(pd_err[i]),
            'PA': float(pa[i]), 'PA_err': float(pa_err[i]),
            'QN': float(spc.QN[i]), 'QN_err': float(spc.QN_ERR[i]),
            'UN': float(spc.UN[i]), 'UN_err': float(spc.UN_ERR[i]),
            'MDP99': float(spc.MDP_99[i]), 'SIGNIF': float(spc.SIGNIF[i]),
        }
        for i in range(len(pd))
    ]


_PCUBE_HEADER_PREFIX = '# ixpe_pcube provenance: '
_PCUBE_FIELDS = ['band_id', 'E_LO', 'E_HI', 'PD', 'PD_err', 'PA', 'PA_err',
                  'QN', 'QN_err', 'UN', 'UN_err', 'MDP99', 'SIGNIF']


def write_pcube_csv(csv_path: str, provenance: dict, bands: list[dict]) -> None:
    """Write one provenance-headed CSV row per energy band (pcube_bands' output).

    Mirrors pcube_timeseg's provenance-comment-header convention (first line
    '# ixpe_pcube provenance: {json}'), but writes the whole table in one shot
    -- there is no resume/checkpoint state, since one ixpe_pcube call runs a
    single xpselect/xpbin reduction, not a segment loop. There is no separate
    "integrated" sentinel row (unlike pcube_timeseg's segment_id=-1): a caller
    wanting a single broadband measurement passes a 2-edge ebinning, which
    naturally yields exactly one row spanning the full requested range.
    """
    with open(csv_path, 'w', newline='') as f:
        f.write(_PCUBE_HEADER_PREFIX + json.dumps(provenance, sort_keys=True) + '\n')
        writer = csv.DictWriter(f, fieldnames=_PCUBE_FIELDS)
        writer.writeheader()
        for i, band in enumerate(bands):
            writer.writerow({'band_id': i, **band})


def read_pcube_csv(csv_path: str) -> tuple[dict, list[dict]]:
    """Read a write_pcube_csv CSV back into (provenance, bands).

    Raises ValueError if the provenance header is missing/unrecognised --
    mirrors _read_existing_timeseg_csv's "unrecognised means absent" contract,
    but raises here (there is no resume mode to silently fall back to).
    """
    with open(csv_path, newline='') as f:
        header_line = f.readline()
        if not header_line.startswith(_PCUBE_HEADER_PREFIX):
            raise ValueError(f'{csv_path}: missing ixpe_pcube provenance header')
        provenance = json.loads(header_line[len(_PCUBE_HEADER_PREFIX):])
        rows = list(csv.DictReader(f))
    bands = [{k: (int(v) if k == 'band_id' else float(v)) for k, v in r.items()}
              for r in rows]
    for band in bands:
        del band['band_id']
    return provenance, bands

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
    When changes=True: src_rad, bkg_inner_rad, bkg_outer_rad, ra, dec are required
    (irf_date is auto-selected when None).

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

    irf_date: str | None
        IRF calibration epoch 'YYYYMMDD', resolving to 'ixpe:obssim{irf_date}_alpha075:v13'.
        Leave None (default) to auto-select the installed epoch nearest the observation
        date (irf_epoch.resolve_irf_date); pass a value only to override.

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
        # irf_date=None -> auto-select the installed epoch nearest the observation.
        irf_date = _resolve_irf_date(obs, file_list, irf_date, chatter)
        # Validate that all required parameters are provided when changes=True
        if any(param is None for param in [src_rad, bkg_inner_rad, bkg_outer_rad, ra, dec]):
            raise ValueError('When changes=True, all parameters (src_rad, bkg_inner_rad, bkg_outer_rad, ra, dec) must be provided.')
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
               irf_date: str | None = None,
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

    irf_date: str | None
        IRF calibration epoch 'YYYYMMDD', resolving to 'ixpe:obssim{irf_date}_alpha075:v13'.
        Leave None (default) to auto-select the installed epoch nearest the observation
        date; pass a value only to override. Used only when a fresh reduction is required;
        the resolved epoch is returned under the 'irf_date' key.

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
    # irf_date=None -> auto-select the installed epoch nearest the observation.
    irf_date = _resolve_irf_date(obs, file_list, irf_date, chatter)

    src_pha, bkg_pha = pha(obs, file_list, changes=False, chatter=chatter)
    if all(os.path.exists(f) for f in src_pha + bkg_pha):
        if chatter:
            print(f'reduce_pha: reusing existing PHA1 products for {obs}')
        return {'src_pha': src_pha, 'bkg_pha': bkg_pha, 'ra': None, 'dec': None,
                'reused': True, 'src_reg': None, 'bkg_reg': None, 'irf_date': irf_date}

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
            'reused': False, 'src_reg': src_reg, 'bkg_reg': bkg_reg, 'irf_date': irf_date}

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

def _resolve_irf_date(obs: str | list[str], file_list: list[str], irf_date: str | None,
                      chatter: bool = False) -> str:
    """Return ``irf_date`` when the caller pinned one, else auto-select it.

    Auto-selection reads the observation midpoint MJD (``get_mjd``) and asks
    :func:`irf_epoch.resolve_irf_date` for the installed IXPE calibration epoch
    nearest that date. Centralised here so ``pcube``/``pha``/``reduce_pha``/
    ``pcube_timeseg`` share one policy; the hardcoded ``'YYYYMMDD'`` default is
    gone (a wrong epoch silently biased a live reproduction -- see irf_epoch).

    ``obs`` may be a single ObsID or a list (a native multi-ObsID ``pcube``
    combination): the midpoint used is 0.5*(earliest start + latest stop)
    across every listed ObsID, so one irf_date is resolved for the whole
    combined reduction. For a single ObsID this reduces to that ObsID's own
    midpoint, unchanged.
    """
    if irf_date is not None:
        return irf_date
    obs_list = [obs] if isinstance(obs, str) else list(obs)
    starts, stops = [], []
    for o in obs_list:
        s, e, _ = get_mjd(o, file_list)
        starts.append(s)
        stops.append(e)
    mjd_mid = 0.5 * (min(starts) + max(stops))
    resolved = irf_epoch.resolve_irf_date(mjd_mid)
    if chatter:
        label = obs if isinstance(obs, str) else '+'.join(obs_list)
        print(f'irf_date auto-resolved to {resolved} '
              f'(installed epoch nearest MJD {mjd_mid:.1f}) for {label}')
    return resolved


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
                   irf_date: str | None = None,
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
    irf_date : str or None
        IRF calibration epoch 'YYYYMMDD', resolving to ``ixpe:obssim{irf_date}_alpha075:v13``.
        Leave None (default) to auto-select the installed epoch nearest the observation
        date (``irf_epoch.resolve_irf_date``); pass a value only to override.
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
    # irf_date=None -> auto-select the installed epoch nearest the observation.
    irf_date = _resolve_irf_date(obs, file_list, irf_date, chatter)

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


# --------------------------------------------------------------------------- #
# qu_slice_pcubes: MDP-anchored adaptive q(t)=Q/I, u(t)=U/I intra-observation
# time series, with data-driven N-way splitting of merged observations.
# --------------------------------------------------------------------------- #

# Di Marco+2023 (AJ 165, 143) source-region count-rate surface-brightness
# regimes for the background-subtraction decision: bright (>T_BRIGHT
# cps/arcmin**2) and intermediate (T_FAINT-T_BRIGHT) sources leave the
# background annulus itself source-dominated, so subtraction should be
# avoided; faint (<T_FAINT) sources let rejection work and the residual
# background should be subtracted. Only T_FAINT drives the binary
# subtract_background decision (bright and intermediate both resolve to no
# subtraction); T_BRIGHT is carried in provenance to document the regime
# boundary, not used in the decision itself.
_R_PROXY_T_FAINT = 1.0   # cps/arcmin**2
_R_PROXY_T_BRIGHT = 2.0  # cps/arcmin**2

_QU_SLICE_HEADER_PREFIX = '# qu_slice_pcubes provenance: '
_QU_SLICE_FIELDS = ['slice_id', 't_mean_met', 't_start_met', 't_stop_met',
                     'good_exposure_s', 'counts', 'q', 'q_err', 'u', 'u_err',
                     'PD_pct', 'PA_deg', 'MDP99_pct', 'under_target']


def _common_du_bounds(obs_file_list: list[str]) -> tuple[float, float]:
    """Common physical MET span every file in ``obs_file_list`` actually
    covers: ``(max over files of that file's own EVENTS TSTART, min over
    files of that file's own EVENTS TSTOP)``.

    These are the same per-file TSTART/TSTOP keywords ixpeobssim's own
    ``xEventSelect._validate`` checks a requested ``tmin``/``tmax`` against,
    one file at a time, before running the selection. Real IXPE DU1/DU2/DU3
    Level-2 files for the same ObsID can carry TSTART/TSTOP differing by
    milliseconds up to hundreds of seconds (a merged/re-processed
    observation), so a window edge valid for one DU's own bounds can lie
    outside a sibling DU's -- this is the span every contributing file's own
    bounds cover in common.
    """
    tstarts, tstops = [], []
    for f in obs_file_list:
        with fits.open(f) as hdul:
            hdr = hdul['EVENTS'].header
            tstarts.append(float(hdr['TSTART']))
            tstops.append(float(hdr['TSTOP']))
    return max(tstarts), min(tstops)


def _reduce_window(obs_file_list: list[str], suffix: str, tmin: float | None, tmax: float | None,
                    ra: float, dec: float, emin_kev: float, emax_kev: float, src_rad: float,
                    bkg: bool, bkg_inner_rad: float, bkg_outer_rad: float, irf_date: str,
                    overwrite: bool = True) -> tuple['xBinnedPolarizationCube', float | None, float | None]:
    """One weighted, DU-merged PCUBE over ``[tmin, tmax)`` (or the whole
    observation when ``tmin`` is ``None``), background-subtracted iff
    ``bkg`` -- the same ``xpselect(tmin=, tmax=)``/``xpbin`` pattern
    ``pcube_timeseg``'s own ``_reduce`` closure uses for its fixed
    elapsed-time segments, generalised here to an arbitrary caller-chosen
    window and reused for both a sub-observation's whole-span cube and every
    one of its adaptively-placed slices.

    **Per-DU physical-bounds clipping.** A single ``xpselect(tmin=, tmax=)``
    call spans all three DU files at once, but ``xEventSelect._validate``
    checks that one ``tmin``/``tmax`` pair against EACH file's own EVENTS
    TSTART/TSTOP in turn and aborts the instant either edge falls outside
    any single file's own bounds -- even when the window is entirely valid
    for the other two DUs. When ``tmin`` is not ``None`` this clips the
    requested window to the span every contributing DU file covers in
    common (``_common_du_bounds``): ``tmin_eff = max(tmin, tstart_common)``,
    ``tmax_eff = min(tmax, tstop_common)``. Raises ``ValueError`` -- never
    calling ``xpselect`` at all -- when the clipped window is empty
    (``tmax_eff <= tmin_eff``), i.e. the requested window lies entirely
    outside the common span; this is a genuine, honest failure of that
    window, not a cube to invent. The whole-observation call (``tmin`` is
    ``None``) is never clipped.

    Every ``xpselect``/``xpbin`` call runs with stdout+stderr captured into
    an in-memory buffer, discarded on success. ixpeobssim's own
    ``logging_.abort()`` is a bare ``sys.exit`` -- on that ``SystemExit``
    this raises ``RuntimeError`` embedding the last ~15 captured lines (the
    diagnostic ixpeobssim prints right before "Cannot continue, abort.":
    missing IRF, empty selection, missing extension, ...), so the real cause
    reaches the caller instead of a bare, uninformative abort message.

    Returns
    -------
    tuple[xBinnedPolarizationCube, float | None, float | None]
        The reduced cube, and the effective ``(tmin, tmax)`` actually passed
        to ``xpselect`` -- ``(None, None)`` when ``tmin`` was ``None``, else
        the (possibly per-DU-clipped, and therefore narrower) window. A
        caller computing GTI-derived exposure or a mean event-arrival
        timestamp for this window should recompute it against this
        returned, effective window, not the originally-requested one, so
        the reported quantity matches the events actually selected.
    """
    if tmin is not None:
        tstart_common, tstop_common = _common_du_bounds(obs_file_list)
        tmin_eff = max(tmin, tstart_common)
        tmax_eff = min(tmax, tstop_common)
        if tmax_eff <= tmin_eff:
            raise ValueError(
                f'_reduce_window({suffix}): requested window [{tmin}, {tmax}) does not '
                f'overlap the common per-DU physical span [{tstart_common}, {tstop_common})')
    else:
        tmin_eff = tmax_eff = None

    bin_kwargs = dict(algorithm='PCUBE', ebinalg='LIST', ebinning=[emin_kev, emax_kev],
                       overwrite=overwrite, irfname=f'ixpe:obssim{irf_date}_alpha075:v13', weights=True)
    sel_kwargs = dict(rad=src_rad, ra=ra, dec=dec, suffix=f'{suffix}src', overwrite=overwrite)
    if tmin_eff is not None:
        sel_kwargs.update(tmin=tmin_eff, tmax=tmax_eff)

    buf = io.StringIO()
    try:
        with redirect_stdout(buf), redirect_stderr(buf):
            src_files = pipeline.xpselect(*obs_file_list, **sel_kwargs)
            spc = xBinnedPolarizationCube.from_file_list(pipeline.xpbin(*src_files, **bin_kwargs))
            if bkg:
                bkg_sel_kwargs = dict(ra=ra, dec=dec, innerrad=bkg_inner_rad, rad=bkg_outer_rad,
                                       suffix=f'{suffix}bkg', overwrite=overwrite)
                if tmin_eff is not None:
                    bkg_sel_kwargs.update(tmin=tmin_eff, tmax=tmax_eff)
                bkg_files = pipeline.xpselect(*obs_file_list, **bkg_sel_kwargs)
                bpc = xBinnedPolarizationCube.from_file_list(pipeline.xpbin(*bkg_files, **bin_kwargs))
                bpc *= spc.backscal() / bpc.backscal()
                spc -= bpc
    except SystemExit as exc:
        tail = '\n'.join(buf.getvalue().splitlines()[-15:])
        raise RuntimeError(
            f'_reduce_window({suffix}): ixpeobssim pipeline aborted ({exc}); '
            f'last output:\n{tail}'
        ) from exc
    return spc, tmin_eff, tmax_eff


def _centroid_ra_dec(evt_path: str, src_rad: float, bkg_inner_rad: float, bkg_outer_rad: float) -> tuple[float, float]:
    """Centroid a source position from one Level-2 event file via
    ``find_src.run_extraction`` -- the same auto-centroid ``pcube_timeseg``
    runs when ``ra``/``dec`` are both omitted -- parsed from the DS9 region
    it writes to a scratch directory."""
    with tempfile.TemporaryDirectory(prefix='qu_slice_centroid_') as tmp_dir:
        src_reg = os.path.join(tmp_dir, 'src.reg')
        with open(os.devnull, 'w') as fnull:
            with redirect_stdout(fnull):
                find_src.run_extraction(evt_path, src_reg, os.path.join(tmp_dir, 'bkg.reg'), None,
                                         src_rad, bkg_inner_rad, bkg_outer_rad)
        with open(src_reg) as f:
            region_line = f.readlines()[1]
        m = re.match(r'circle\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,', region_line)
        if not m:
            raise ValueError(f'qu_slice_pcubes: could not parse source region from {src_reg}: {region_line!r}')
        return float(m.group(1)), float(m.group(2))


def _col_num(header, ttype: str) -> int:
    """FITS column number (1-based) whose ``TTYPE`` keyword equals ``ttype``."""
    for key, val in header.items():
        if key.startswith('TTYPE') and val == ttype:
            return int(key[5:])
    raise ValueError(f'qu_slice_pcubes: no {ttype!r} column found in header')


def _sep_arcmin(ra: np.ndarray, dec: np.ndarray, ra0: float, dec0: float) -> np.ndarray:
    """Flat-sky angular separation (arcmin) from ``(ra0, dec0)`` -- adequate
    at the few-arcmin scale of an IXPE source/background region. Used only
    to place the fast slice-edge/``r_proxy`` proxies below, never for the
    PCUBE measurement itself (which always goes through ixpeobssim's own
    ``xpselect`` region selection via ``_reduce_window``)."""
    dra = (ra - ra0) * np.cos(np.radians(dec0))
    ddec = dec - dec0
    return np.sqrt(dra ** 2 + ddec ** 2) * 60.0


def _pi_energy_mask(pi: np.ndarray, du: int, emin_kev: float, emax_kev: float) -> np.ndarray:
    """Boolean mask selecting events whose PI channel's installed-RMF
    EBOUNDS midpoint falls in ``[emin_kev, emax_kev)`` -- the fixed
    instrumental channel grid shared across every IRF date for a given DU
    (see ``_installed_du_rmf``)."""
    rmf = _installed_du_rmf(du)
    with fits.open(rmf) as hdul:
        e_lo = hdul['EBOUNDS'].data['E_MIN'].astype(float)
        e_hi = hdul['EBOUNDS'].data['E_MAX'].astype(float)
    in_band = _ebounds_channel_mask(e_lo, e_hi, (emin_kev, emax_kev))
    pi = np.asarray(pi, dtype=int)
    mask = np.zeros(pi.shape, dtype=bool)
    valid = (pi >= 0) & (pi < len(in_band))
    mask[valid] = in_band[pi[valid]]
    return mask


def _selected_event_times(evt_path: str, du: int, ra: float, dec: float,
                           src_rad: float, emin_kev: float, emax_kev: float) -> np.ndarray:
    """Event ``TIME`` (MET s) values from one DU's Level-2 file that fall
    inside the circular source aperture (``src_rad`` arcmin, via the EVENTS
    table's own linear TAN-projection column WCS for the X/Y columns -- the
    same table-column WCS ``find_src.get_table_wcs`` reads for IXPE) and the
    ``[emin_kev, emax_kev)`` band.

    A fast, un-weighted, un-background-subtracted proxy for two things that
    never touch the actual polarimetric measurement itself: the source-region
    count-rate surface brightness (``r_proxy``) and the adaptive slice-edge
    placement walk (``_greedy_slice_edges``). The measurement always goes
    through ``_reduce_window`` (ixpeobssim's own ``xpselect``/``xpbin``/
    weights machinery, never hand-rolled Stokes arithmetic).
    """
    with fits.open(evt_path) as hdul:
        data = hdul['EVENTS'].data
        header = hdul['EVENTS'].header
        time = np.asarray(data['TIME'], dtype=float)
        pi = data['PI']
        xcol, ycol = _col_num(header, 'X'), _col_num(header, 'Y')
        x = np.asarray(data['X'], dtype=float)
        y = np.asarray(data['Y'], dtype=float)
        ra_evt = header[f'TCRVL{xcol}'] + (x - header[f'TCRPX{xcol}']) * header[f'TCDLT{xcol}']
        dec_evt = header[f'TCRVL{ycol}'] + (y - header[f'TCRPX{ycol}']) * header[f'TCDLT{ycol}']
    region_mask = _sep_arcmin(ra_evt, dec_evt, ra, dec) <= src_rad
    energy_mask = _pi_energy_mask(pi, du, emin_kev, emax_kev)
    return time[region_mask & energy_mask]


def _merged_selected_times(obs_file_list: list[str], ra: float, dec: float,
                            src_rad: float, emin_kev: float, emax_kev: float) -> np.ndarray:
    """Sorted, DU1+DU2+DU3-merged ``_selected_event_times``."""
    times = np.concatenate([
        _selected_event_times(f, _parse_du(f), ra, dec, src_rad, emin_kev, emax_kev)
        for f in obs_file_list
    ])
    return np.sort(times)


def r_proxy(n_events: int, good_exposure_s: float, src_rad: float) -> float:
    """Source-region count-rate surface brightness (cps/arcmin**2) -- the
    Di Marco+2023 (AJ 165, 143) brightness proxy that decides background
    subtraction. Returns ``0.0`` for non-positive exposure (never a
    ``ZeroDivisionError``)."""
    if good_exposure_s <= 0:
        return 0.0
    area = math.pi * src_rad ** 2
    return n_events / (good_exposure_s * area)


def _decide_background(bkg: bool | None, r_proxy_val: float,
                        t_faint: float = _R_PROXY_T_FAINT) -> bool:
    """``bkg`` overrides when given (True/False); otherwise subtract only in
    the faint regime (``r_proxy_val < t_faint``) -- both the intermediate and
    bright regimes resolve to no subtraction (Di Marco+2023's caution that a
    source-dominated background annulus makes subtraction counterproductive)."""
    return bkg if bkg is not None else (r_proxy_val < t_faint)


def _sorted_gap_split_indices(gti_start: list[float], gti_stop: list[float],
                               split_gap_factor: float) -> list[int]:
    """GTI-array indices after which a genuine multi-pointing merge gap
    sits, found by a MAX-RATIO CLIFF on the SORTED gap list -- never against
    the population median/mean of all gaps. Real IXPE GTI structure mixes
    many small on-orbit dropout gaps (fractions of a second to a few ks)
    with, at most, a handful of genuine month-scale merge gaps; comparing
    every gap to the population median flags hundreds of ordinary orbital
    dropouts as "splits".

    Ranks every inter-GTI gap descending, then finds the single LARGEST
    adjacent ratio ``gaps[k] / gaps[k+1]`` in that sorted list -- the cliff
    separating the top tier of genuine (comparable-magnitude) merge gaps
    from the ordinary orbital-dropout gaps below it -- and accepts every gap
    at or above that cliff. This differs from (and replaces) a strictly
    descending chain of outliers, which fails outright on a merge with two
    or more comparable dominant gaps (e.g. sorted gaps
    ``[5.0e6, 4.8e6, 2.8e3, ...]``: the top-vs-second-largest ratio
    ``5.0e6/4.8e6 ~ 1.04`` never clears ``split_gap_factor``, so a naive
    top-down chain returns no split at all -- even though BOTH gaps are
    genuine merge boundaries, comparably huge against every ordinary gap
    below them). The max-ratio cliff instead looks for where the biggest
    proportional DROP occurs, wherever that is in the ranked list:

    - A single dominant gap (e.g. ObsID 03002599's ~5.4e6s literature-merge
      gap against its next-largest ~2.8e3s ordinary gap, ratio ~1900) puts
      the biggest cliff at the very top -- exactly that one gap is accepted.
    - Two (or more) comparable dominant gaps put the biggest cliff below
      ALL of them (their mutual ratio is small; the drop to the ordinary
      gaps below is what is large) -- every one of that top tier is
      accepted together.
    - A tiered case (e.g. sorted gaps ``[5e6, 4e5, 2e3]``, adjacent ratios
      12.5 then 200) prefers the LARGEST cliff: the 200x drop between
      ``4e5`` and ``2e3`` beats the 12.5x drop between ``5e6`` and ``4e5``,
      so the top 2 gaps are accepted, not just the top 1.
    - No adjacent ratio exceeding ``split_gap_factor`` anywhere in the
      ranked list means no genuine cliff was found at all -- no split.

    A ``gap_next <= 0`` (degenerate/overlapping GTI) is treated as an
    infinite ratio, i.e. a definite cliff.
    """
    if len(gti_start) < 2:
        return []
    gaps = [(gti_start[i + 1] - gti_stop[i], i) for i in range(len(gti_start) - 1)]
    ranked = sorted(gaps, key=lambda pair: -pair[0])
    if len(ranked) < 2:
        return []

    best_k, best_ratio = -1, -1.0
    for k in range(len(ranked) - 1):
        gap_k = ranked[k][0]
        gap_next = ranked[k + 1][0]
        ratio = float('inf') if gap_next <= 0 else gap_k / gap_next
        if ratio > best_ratio:
            best_ratio, best_k = ratio, k

    if best_ratio <= split_gap_factor:
        return []

    accepted = [ranked[j][1] for j in range(best_k + 1)]
    return sorted(accepted)


def subobs_groups(gti_start: list[float], gti_stop: list[float],
                   split_gap_factor: float = 10.0) -> list[dict]:
    """Partition one observation's GTI list into N sub-observations at
    genuine multi-pointing merge gaps (see ``_sorted_gap_split_indices``).

    Returns a chronologically-ordered list of dicts, each with ``gti_start``,
    ``gti_stop`` (that sub-observation's own GTI sub-lists) and ``tstart``,
    ``tstop`` (first GTI start / last GTI stop in the group). A single-group
    result (no split detected) is still returned as a 1-element list, so
    callers never special-case "no split" separately from "N-way split".
    """
    cuts = _sorted_gap_split_indices(gti_start, gti_stop, split_gap_factor)
    groups = []
    start = 0
    for cut in cuts + [len(gti_start) - 1]:
        gs = gti_start[start:cut + 1]
        ge = gti_stop[start:cut + 1]
        groups.append({'gti_start': gs, 'gti_stop': ge, 'tstart': gs[0], 'tstop': ge[-1]})
        start = cut + 1
    return groups


def _should_split_subobs(obsid: str, split_multiobs: bool | None = None) -> bool:
    """Whether an N-way GTI-gap split should even be attempted for ``obsid``.

    ``split_multiobs=None`` (default) auto-decides from the IXPE
    literature-merged-ObsID naming convention: ``obsid`` ending in ``'99'``
    (e.g. ``03002599`` -> sub-obs a/b) is judged a merged observation and the
    split machinery (``subobs_groups``/``_sorted_gap_split_indices``, both
    scale-free and otherwise untouched) is allowed to run; any other obsid --
    an ordinary single monitoring pointing, however large its own inter-GTI
    gaps happen to be (e.g. Cyg X-1's 02008201/02008301/...) -- is never
    assumed merged and is always treated as one whole-span sub-observation.
    ``split_multiobs=True``/``False`` forces the decision either way,
    overriding the naming convention.
    """
    if split_multiobs is not None:
        return split_multiobs
    return obsid.endswith('99')


def _resolve_subobs_groups(obs: str, gti_start: list[float], gti_stop: list[float],
                            split_gap_factor: float, split_multiobs: bool | None) -> list[dict]:
    """The exact GTI partition ``qu_slice_pcubes`` reduces per sub-observation:
    ``subobs_groups(...)`` when ``_should_split_subobs(obs, split_multiobs)``
    is True, else a single group spanning the whole GTI list unchanged.
    ``subobs_groups``/``_sorted_gap_split_indices`` themselves stay fully
    general and untouched -- this is the one gate deciding whether they run
    at all.
    """
    if _should_split_subobs(obs, split_multiobs):
        return subobs_groups(gti_start, gti_stop, split_gap_factor)
    return [{'gti_start': list(gti_start), 'gti_stop': list(gti_stop),
             'tstart': gti_start[0], 'tstop': gti_stop[-1]}]


def _letter_suffix(index: int) -> str:
    """0-based ``index`` -> a bijective base-26 lowercase-letter string:
    'a'..'z' for index 0-25 (exactly the original single-letter on-disk
    convention, e.g. ``qu_slices_<obsid>a.csv``), then 'aa','ab',...,'az',
    'ba',... beyond that -- pure ASCII letters at every index, so a >26-way
    split can never emit a non-printable byte (the original ``chr(ord('a')
    + index)`` did, once ``index`` passed 25) into a FITS header value.
    """
    n = index + 1
    letters = []
    while n > 0:
        n, rem = divmod(n - 1, 26)
        letters.append(chr(ord('a') + rem))
    return ''.join(reversed(letters))


def _subobs_key(obs: str, index: int, n_groups: int) -> str:
    """``obs`` unchanged for a single group; ``obs`` + ``_letter_suffix(index)``
    for N>1 groups -- always pure ASCII letters, never a non-printable byte."""
    if n_groups == 1:
        return obs
    return obs + _letter_suffix(index)


def _greedy_slice_edges(times: np.ndarray, tstart: float, tstop: float,
                         target_raw: float, n_min: int, n_max: int) -> tuple[list[float], bool, bool]:
    """Adaptive greedy time-edge placement within one sub-observation.

    Walks ``times`` (sorted MET s, already scoped to the sub-observation's
    own span) accumulating a running raw-event count and closes a slice as
    soon as the count reaches ``target_raw``; a trailing remainder that never
    reaches target is merged into the previous slice (never emitted as its
    own straggler) unless it is the only slice produced. The natural slice
    count is then clipped to ``[n_min, n_max]``: exceeding ``n_max``
    re-partitions ``times`` into exactly ``n_max`` equal-raw-count chunks;
    falling short of ``n_min`` is accepted as-is (this walk never invents
    slices past what the real events support) and signalled via the
    returned ``noise_limited``.

    Returns ``(edges, noise_limited, clipped_to_n_max)`` where ``edges`` has
    ``n_slices + 1`` entries (``edges[0] == tstart``, ``edges[-1] == tstop``).
    """
    n_events = len(times)
    if n_events == 0 or target_raw <= 0:
        return [tstart, tstop], True, False

    edges = [tstart]
    cum = 0
    for t in times:
        cum += 1
        if cum >= target_raw:
            edges.append(float(t))
            cum = 0
    # The walk above always closes a slice the moment it reaches target_raw
    # (even at the very last event), so whatever remains after the final
    # closure is a residual of 0..target_raw-1 events over [last_closure,
    # tstop) -- fold it into the previous slice by dropping its own boundary
    # rather than emit it as its own (possibly EMPTY, if the residual is 0)
    # straggler, unless it is the ONLY slice (nothing to merge into).
    if edges[-1] != tstop:
        edges.append(tstop)
        if len(edges) > 2:
            del edges[-2]

    n_natural = len(edges) - 1
    if n_natural > n_max:
        idxs = sorted({int(round(k * n_events / n_max)) for k in range(1, n_max)})
        idxs = [i for i in idxs if 0 < i < n_events]
        edges = [tstart] + [float(times[i]) for i in idxs] + [tstop]
        return edges, False, True

    return edges, n_natural < n_min, False


def qu_slice_pcubes(obs: str,
                     file_list: list[str],
                     ra: float | None = None,
                     dec: float | None = None,
                     emin_kev: float = 2.0,
                     emax_kev: float = 8.0,
                     src_rad: float = 1.0,
                     bkg: bool | None = None,
                     bkg_inner_rad: float = 3.0,
                     bkg_outer_rad: float = 5.0,
                     mdp_target: float | None = None,
                     mdp_frac: float = 0.5,
                     mdp_floor: float = 0.015,
                     n_min: int = 8,
                     n_max: int = 25,
                     split_gap_factor: float = 10.0,
                     split_multiobs: bool | None = None,
                     irf_date: str | None = None,
                     overwrite: bool = True,
                     chatter: bool = False,
                     ) -> dict:
    """Model-independent, MDP-anchored adaptive q(t)=Q/I, u(t)=U/I intra-
    observation time series, DU1/2/3 merged, with data-driven N-way splitting
    of a merged/multi-pointing observation.

    **N-way split.** The split is attempted at all ONLY when
    ``_should_split_subobs(obs, split_multiobs)`` judges ``obs`` a
    literature-merged observation -- by default, an obsid ending in ``'99'``
    (e.g. 03002599); ``split_multiobs=True``/``False`` overrides that
    naming-convention gate. An ordinary single monitoring pointing (e.g.
    Cyg X-1's 02008201) is never assumed merged, however large its own
    inter-GTI gaps happen to be, and is always reduced as one whole-span
    sub-observation -- ``subobs_groups``/``_sorted_gap_split_indices``
    (the scale-free gap-outlier machinery itself) stay fully general and
    untouched; they simply are not invoked at all when the gate is closed.
    When the gate is open, the observation's GTI list is partitioned into
    one or more sub-observations at genuine multi-pointing merge gaps
    (``subobs_groups``/``split_gap_factor``); each sub-observation is keyed
    ``<obs>`` (no split) or ``<obs>a``, ``<obs>b``, ... (N-way split) and
    processed fully independently below.

    **Per sub-observation:**

    1. A source-region count-rate surface brightness ``r_proxy`` (cps/
       arcmin**2, ``r_proxy()``) decides background subtraction when
       ``bkg=None`` (Di Marco+2023, AJ 165 143): faint (``r_proxy <
       T_FAINT=1.0``) subtracts an annulus, intermediate/bright does not.
       ``bkg=True``/``False`` overrides; the resolved choice is applied to
       both the whole-sub-obs cube and every one of its slices.
    2. One whole-sub-observation weighted, DU-merged PCUBE (``_reduce_window``
       over the sub-obs's own ``[tstart, tstop)``, per-DU-clipped -- see
       ``_reduce_window`` -- when a sibling DU's own physical TSTART/TSTOP is
       narrower) gives ``MDP99_obs`` and total good-time ``COUNTS``
       (``C_obs``). The effective per-slice MDP target is ``mdp_target`` if
       given, else ``max(mdp_frac * PD_obs, mdp_floor)`` (all as
       dimensionless fractions); ``C_target = C_obs * (MDP99_obs / target) ** 2``.
    3. Adaptive greedy slicing (``_greedy_slice_edges``) places slice edges
       at real event arrival times so each slice's raw-proxy count
       approximates ``C_target``, clipped to ``[n_min, n_max]`` slices (see
       that function's docstring for the exact clipping rule).
    4. Each final slice re-runs ``_reduce_window`` over its own
       ``[tmin, tmax)`` (same background choice as the whole sub-obs),
       reading ``q``/``u``/errors/``PD``/``PA``/``MDP99``/``COUNTS`` via
       ``pcube_bands``/direct attribute access -- never hand-rolled Stokes
       arithmetic -- and is timestamped by the MEAN of the real event
       arrival times inside its window (``t_mean_met``), NOT the window
       midpoint -- this deliberately diverges from ``pcube_timeseg``.
       ``good_exposure_s``/``t_mean_met`` are both computed against
       ``_reduce_window``'s returned effective (possibly per-DU-clipped)
       window, so they describe exactly the events the slice's cube was
       built from.

    Writes one provenance-headed CSV per sub-observation,
    ``qu_slices_<key>.csv``, under the observation's own directory (see
    ``_QU_SLICE_FIELDS`` for the row schema). Every sub-observation is fully
    reduced (whole cube + slices) BEFORE any CSV is written -- the write
    pass runs only after the entire per-sub-observation loop has completed
    without raising, so a failure partway through (e.g. sub-observation "b"
    of an N-way split) leaves NO CSV at all on disk, not even for a sub-obs
    ("a") that finished its own reduction successfully first.

    Parameters
    ----------
    obs, file_list : see ``_get_obs_file_list``.
    ra, dec : float or None
        Source position in degrees. If both ``None``, centroided from DU1's
        event file via ``find_src.run_extraction``.
    emin_kev, emax_kev : float
        Single energy bin edges in keV (also the ``r_proxy`` reference band).
        Default 2.0-8.0.
    src_rad, bkg_inner_rad, bkg_outer_rad : float
        Source region radius / background annulus radii, arcmin.
    bkg : bool or None
        Background-subtraction override; ``None`` (default) auto-decides
        from ``r_proxy`` per sub-observation.
    mdp_target : float or None
        Per-slice MDP99 target as a dimensionless fraction (e.g. 0.02 for
        2%). ``None`` (default) derives it from ``mdp_frac``/``mdp_floor``.
    mdp_frac, mdp_floor : float
        ``target = max(mdp_frac * PD_obs, mdp_floor)`` when ``mdp_target``
        is ``None``. Defaults 0.5, 0.015.
    n_min, n_max : int
        Slice-count clipping bounds per sub-observation. Defaults 8, 25.
    split_gap_factor : float
        Outlier-ratio threshold for the N-way GTI-gap split. Default 10.0.
        Only consulted at all when the split is attempted -- see
        ``split_multiobs``.
    split_multiobs : bool or None
        Force (``True``) or forbid (``False``) attempting the N-way GTI-gap
        split at all. ``None`` (default) auto-decides via
        ``_should_split_subobs``: an obsid ending in ``'99'`` (e.g.
        03002599) is judged literature-merged and the split runs; any other
        obsid is treated as one whole-span sub-observation regardless of its
        own GTI gap sizes.
    irf_date : str or None
        IRF calibration epoch; ``None`` auto-selects the installed epoch
        nearest the observation date.

    Returns
    -------
    dict with ``subobs_keys`` (list) and ``per_subobs`` (key -> dict with
    ``mdp99_obs_pct``, ``r_proxy``, ``subtract_background``, ``uses_rej``,
    ``n_slice``, ``noise_limited``, ``whole_pd_pct``, ``whole_pd_err``,
    ``whole_pa_deg``, ``whole_signif``, ``csv_path``, and ``rows`` -- the
    per-slice row dicts written to that sub-obs's CSV).

    Raises
    ------
    ValueError
        ``file_list`` does not resolve to exactly one Level-2 event file per
        DU1/DU2/DU3 for ``obs``; exactly one of ``ra``/``dec`` is given
        without the other; or a sub-observation has no GTI overlap over its
        own span or no usable (region+band-selected) events/counts to derive
        an MDP target or slice from. Any such failure -- even one raised on
        a LATER sub-observation, after an earlier one's reduction already
        succeeded -- leaves no CSV on disk at all, since every CSV write is
        deferred until the whole per-sub-observation loop completes.
    """
    obs_file_list = _get_obs_file_list(obs, file_list, chatter=chatter)
    obs_file_list = _validate_du_files(obs, obs_file_list, caller='qu_slice_pcubes')
    irf_date = _resolve_irf_date(obs, file_list, irf_date, chatter)
    uses_rej = any('_rej' in os.path.basename(f) for f in obs_file_list)

    if (ra is None) != (dec is None):
        raise ValueError(
            f'qu_slice_pcubes: provide both ra and dec, or neither (auto-centroid); '
            f'got ra={ra}, dec={dec}')
    if ra is None:
        ra, dec = _centroid_ra_dec(obs_file_list[0], src_rad, bkg_inner_rad, bkg_outer_rad)

    with fits.open(obs_file_list[0]) as hdul:
        if 'GTI' not in hdul or len(hdul['GTI'].data) == 0:
            raise ValueError(f'{obs_file_list[0]}: no GTI extension data')
        gti = hdul['GTI'].data
        gti_start = [float(x) for x in gti['START']]
        gti_stop = [float(x) for x in gti['STOP']]
        deadc = hdul['PRIMARY'].header.get('DEADC')
        if deadc is None:
            raise ValueError(f'{obs_file_list[0]}: PRIMARY header missing DEADC keyword')
        deadc = float(deadc)

    all_times = _merged_selected_times(obs_file_list, ra, dec, src_rad, emin_kev, emax_kev)
    groups = _resolve_subobs_groups(obs, gti_start, gti_stop, split_gap_factor, split_multiobs)
    obs_dir = Path(obs_file_list[0]).resolve().parent.parent

    subobs_keys = []
    per_subobs = {}
    pending_writes = []  # (csv_path, provenance, rows) -- written only after every sub-obs succeeds
    for i, group in enumerate(groups):
        key = _subobs_key(obs, i, len(groups))
        subobs_keys.append(key)

        sub_tstart, sub_tstop = group['tstart'], group['tstop']
        overlap = _gti_overlap(gti_start, gti_stop, sub_tstart, sub_tstop)
        if overlap is None:
            raise ValueError(f'qu_slice_pcubes: sub-obs {key} has no GTI overlap over its own span')
        good_exposure_whole = overlap[2] * deadc

        sub_times = all_times[(all_times >= sub_tstart) & (all_times < sub_tstop)]
        n_sub_raw = len(sub_times)
        r_proxy_val = r_proxy(n_sub_raw, good_exposure_whole, src_rad)
        subtract_background = _decide_background(bkg, r_proxy_val)

        whole_cube, whole_tmin_eff, whole_tmax_eff = _reduce_window(
            obs_file_list, f'{key}whole', sub_tstart, sub_tstop, ra, dec,
            emin_kev, emax_kev, src_rad, subtract_background,
            bkg_inner_rad, bkg_outer_rad, irf_date, overwrite)
        whole_band = pcube_bands(whole_cube)[0]
        c_obs = float(whole_cube.COUNTS[0])
        if c_obs <= 0 or n_sub_raw == 0:
            raise ValueError(
                f'qu_slice_pcubes: sub-obs {key} has no usable counts '
                f'(C_obs={c_obs}, raw region+band events={n_sub_raw})')

        # Recompute the reported whole-sub-obs exposure against the window
        # _reduce_window actually selected (per-DU-clipped when a sibling
        # DU's own physical TSTART/TSTOP is narrower than sub_tstart/
        # sub_tstop -- see _reduce_window), so good_exposure_whole matches
        # the events the whole-sub-obs cube was actually built from.
        whole_overlap = _gti_overlap(gti_start, gti_stop, whole_tmin_eff, whole_tmax_eff)
        if whole_overlap is None:
            raise ValueError(
                f'qu_slice_pcubes: sub-obs {key} has no GTI overlap over its '
                f'clipped span [{whole_tmin_eff},{whole_tmax_eff})')
        good_exposure_whole = whole_overlap[2] * deadc

        mdp99_obs_frac = whole_band['MDP99']
        pd_obs_frac = whole_band['PD']
        target_frac = mdp_target if mdp_target is not None else max(mdp_frac * pd_obs_frac, mdp_floor)
        if target_frac <= 0:
            raise ValueError(f'qu_slice_pcubes: sub-obs {key} resolved a non-positive MDP target ({target_frac})')
        c_target = c_obs * (mdp99_obs_frac / target_frac) ** 2
        raw_target = c_target * (n_sub_raw / c_obs)

        edges, noise_limited, clipped = _greedy_slice_edges(sub_times, sub_tstart, sub_tstop,
                                                             raw_target, n_min, n_max)
        n_slices = len(edges) - 1
        target_counts = (c_obs / n_slices) if clipped else c_target

        rows = []
        for s_idx in range(n_slices):
            t0, t1 = edges[s_idx], edges[s_idx + 1]
            slice_cube, s_tmin_eff, s_tmax_eff = _reduce_window(
                obs_file_list, f'{key}s{s_idx:03d}', t0, t1, ra, dec,
                emin_kev, emax_kev, src_rad, subtract_background,
                bkg_inner_rad, bkg_outer_rad, irf_date, overwrite)
            band = pcube_bands(slice_cube)[0]
            counts = float(slice_cube.COUNTS[0])

            # GTI overlap and the mean-arrival timestamp both read off the
            # effective (per-DU-clipped) window _reduce_window actually
            # selected, not the originally-requested [t0, t1), so
            # good_exposure_s/t_mean_met describe exactly the events the
            # slice cube was built from.
            slice_overlap = _gti_overlap(gti_start, gti_stop, s_tmin_eff, s_tmax_eff)
            if slice_overlap is None:
                raise ValueError(
                    f'qu_slice_pcubes: slice {s_idx} of {key} [{s_tmin_eff},{s_tmax_eff}) has no GTI overlap')
            clip_start, clip_stop, overlap_s = slice_overlap
            good_exposure_s = overlap_s * deadc

            in_window = sub_times[(sub_times >= s_tmin_eff) & (sub_times < s_tmax_eff)]
            if len(in_window) == 0:
                raise ValueError(
                    f'qu_slice_pcubes: slice {s_idx} of {key} [{s_tmin_eff},{s_tmax_eff}) has no raw '
                    f'region+band events for the mean-arrival timestamp')
            t_mean = float(np.mean(in_window))

            rows.append({
                'slice_id': s_idx, 't_mean_met': t_mean,
                't_start_met': clip_start, 't_stop_met': clip_stop,
                'good_exposure_s': good_exposure_s, 'counts': counts,
                'q': band['QN'], 'q_err': band['QN_err'],
                'u': band['UN'], 'u_err': band['UN_err'],
                'PD_pct': band['PD'] * 100.0, 'PA_deg': band['PA'],
                'MDP99_pct': band['MDP99'] * 100.0,
                'under_target': counts < target_counts,
            })

        csv_path = str(obs_dir / f'qu_slices_{key}.csv')
        provenance = {
            'emin_kev': emin_kev, 'emax_kev': emax_kev, 'src_rad': src_rad,
            'bkg': subtract_background, 'bkg_inner_rad': bkg_inner_rad, 'bkg_outer_rad': bkg_outer_rad,
            'mdp_target': mdp_target, 'mdp_frac': mdp_frac, 'mdp_floor': mdp_floor,
            'n_min': n_min, 'n_max': n_max, 'split_gap_factor': split_gap_factor,
            'split_multiobs': split_multiobs,
            'ra': ra, 'dec': dec, 'irf_date': irf_date, 'uses_rej': uses_rej,
            'subobs_span': [sub_tstart, sub_tstop],
        }
        pending_writes.append((csv_path, provenance, rows))

        per_subobs[key] = {
            'mdp99_obs_pct': mdp99_obs_frac * 100.0,
            'r_proxy': r_proxy_val,
            'subtract_background': subtract_background,
            'uses_rej': uses_rej,
            'n_slice': n_slices,
            'noise_limited': noise_limited,
            'whole_pd_pct': pd_obs_frac * 100.0,
            'whole_pd_err': whole_band['PD_err'] * 100.0,
            'whole_pa_deg': whole_band['PA'],
            'whole_signif': whole_band['SIGNIF'],
            'whole_counts': c_obs,
            'whole_exposure_s': good_exposure_whole,
            'csv_path': csv_path,
            'rows': rows,
        }

        if chatter:
            print(f'qu_slice_pcubes {key}: {n_slices} slice(s), r_proxy={r_proxy_val:.3f} '
                  f'cps/arcmin^2, subtract_background={subtract_background}, '
                  f'noise_limited={noise_limited} -> {csv_path}')

    # Every sub-observation reduced without raising -- only now write the
    # CSVs, so a failure anywhere in the loop above (even on a later
    # sub-observation, after an earlier one already succeeded) leaves no
    # partial CSV at all on disk.
    for csv_path, provenance, rows in pending_writes:
        with open(csv_path, 'w', newline='') as f:
            f.write(_QU_SLICE_HEADER_PREFIX + json.dumps(provenance, sort_keys=True) + '\n')
            writer = csv.DictWriter(f, fieldnames=_QU_SLICE_FIELDS)
            writer.writeheader()
            writer.writerows(rows)

    return {'subobs_keys': subobs_keys, 'per_subobs': per_subobs}