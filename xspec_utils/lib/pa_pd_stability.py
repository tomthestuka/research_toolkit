"""Intra-observation IXPE PA/PD stability statistics -- pure numpy/scipy.

Reads the per-segment PCUBE tables `red_ixpe.pcube_timeseg` writes
(`pcube_timeseg_<obsid>.csv`) and computes weighted constancy statistics,
moment-estimator intrinsic scatter, peak-to-peak ranges, first-order
structure functions, and a cross-source Spearman correlation of MDP99
against the reduced chi-squared statistics. Imports neither `xspec` nor
`ixpeobssim` (unlike `red_ixpe.py`) so this module -- and the
`pa_pd_stability` entry point -- runs standalone, without the science stack.

The CSV reader below (`_read_timeseg_csv`) is an independent reimplementation
of `red_ixpe._read_existing_timeseg_csv`'s contract (same provenance-header
format, same `segment_id=-1` integrated-row convention, same "no header ->
treat as absent" rule) -- not a shared helper, precisely so this module never
imports `red_ixpe` (which pulls in `ixpeobssim`/`xspec` at module load).
"""

from __future__ import annotations

import csv
import json
import os
from typing import Sequence

import numpy as np
from scipy import stats

_TIMESEG_FIELDS = ['segment_id', 't_start_mjd', 't_stop_mjd', 'exposure_s',
                    'PD_pct', 'PD_err', 'PA_deg', 'PA_err',
                    'QI', 'QI_err', 'UI', 'UI_err', 'MDP99', 'SIGNIF']
_PROVENANCE_HEADER_PREFIX = '# pcube_timeseg provenance: '

CATALOG_FIELDS = [
    'obsid', 'source', 'N_seg',
    'PD_pct_int', 'PD_err_int', 'PA_deg_int', 'PA_err_int', 'MDP99_int', 'SIGNIF_int',
    'chi2_PD', 'dof_PD', 'p_PD',
    'chi2_PA', 'dof_PA', 'p_PA',
    'chi2_QI', 'dof_QI', 'p_QI',
    'chi2_UI', 'dof_UI', 'p_UI',
    'chi2_2D', 'dof_2D', 'p_2D', 'variable',
    'sigma_PA_int', 'PA_int_upper_limit',
    'sigma_PD_int', 'PD_int_upper_limit',
    'PD_ptp', 'PA_ptp',
    'rchi2_PD', 'rchi2_PA',
    'mean_MDP99', 'error',
]

_DEFAULT_SF_LAGS_KS = (5.0, 10.0, 20.0)


def _read_timeseg_csv(csv_path: str) -> tuple[dict, dict, dict | None] | None:
    """Read a `pcube_timeseg_<obsid>.csv` into (provenance, rows_by_segid, integrated).

    Mirrors `red_ixpe._read_existing_timeseg_csv` exactly: returns None when
    `csv_path` does not exist, or its first line is not the provenance
    comment header (`# pcube_timeseg provenance: {...}`) -- an unrecognised
    file is treated as absent, never parsed as headerless data.
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


def _wrap180(d):
    """Wrap a PA difference (degrees) into (-90, +90] -- PA is defined mod 180."""
    return np.mod(np.asarray(d, dtype=float) + 90.0, 180.0) - 90.0


def _weighted_stat(x, err) -> dict:
    """Inverse-variance weighted mean and chi-squared constancy test.

    Excludes non-finite values and non-finite/non-positive errors. `chi2`,
    `p` are NaN and `dof = n-1` when fewer than 2 usable points remain
    (n=0 or n=1); `mean` is the single value when n=1, else NaN.
    """
    x = np.asarray(x, dtype=float)
    err = np.asarray(err, dtype=float)
    mask = np.isfinite(x) & np.isfinite(err) & (err > 0)
    n = int(mask.sum())
    if n < 2:
        mean = float(x[mask][0]) if n == 1 else float('nan')
        return {'mean': mean, 'chi2': float('nan'), 'dof': n - 1, 'p': float('nan'), 'n': n, 'mask': mask}
    w = 1.0 / err[mask] ** 2
    mean = float(np.sum(w * x[mask]) / np.sum(w))
    dof = n - 1
    chi2 = float(np.sum(w * (x[mask] - mean) ** 2))
    p = float(stats.chi2.sf(chi2, dof))
    return {'mean': mean, 'chi2': chi2, 'dof': dof, 'p': p, 'n': n, 'mask': mask}


def _weighted_stat_pa(pa_deg, pa_err) -> dict:
    """Circular (mod-180) inverse-variance weighted PA constancy test.

    The weighted mean is the doubled-angle circular mean
    (`atan2(sum(w sin 2th), sum(w cos 2th)) / 2`, wrapped to (-90, +90]);
    residuals are wrapped via `_wrap180` before the chi-squared sum, so a PA
    flipping between e.g. +89 and -89 deg (2 deg apart, physically) is not
    scored as an ~178 deg jump. Same n<2 fallback as `_weighted_stat`.
    """
    pa_deg = np.asarray(pa_deg, dtype=float)
    pa_err = np.asarray(pa_err, dtype=float)
    mask = np.isfinite(pa_deg) & np.isfinite(pa_err) & (pa_err > 0)
    n = int(mask.sum())
    if n < 2:
        mean = float(pa_deg[mask][0]) if n == 1 else float('nan')
        return {'mean': mean, 'chi2': float('nan'), 'dof': n - 1, 'p': float('nan'), 'n': n, 'mask': mask}
    w = 1.0 / pa_err[mask] ** 2
    theta2 = np.deg2rad(2.0 * pa_deg[mask])
    c = np.sum(w * np.cos(theta2))
    s = np.sum(w * np.sin(theta2))
    mean = float(_wrap180(np.rad2deg(np.arctan2(s, c)) * 0.5))
    resid = _wrap180(pa_deg[mask] - mean)
    dof = n - 1
    chi2 = float(np.sum(w * resid ** 2))
    p = float(stats.chi2.sf(chi2, dof))
    return {'mean': mean, 'chi2': chi2, 'dof': dof, 'p': p, 'n': n, 'mask': mask}


def _moment_intrinsic_scatter(values, errs, mask, wrap: bool = False, mean: float | None = None) -> tuple[float, bool]:
    """Moment-estimator intrinsic scatter: sqrt(max(0, sigma_obs^2 - mean_var)).

    `sigma_obs^2` is the sample variance (ddof=1) of the masked values about
    their mean (the wrapped circular residual variance when `wrap=True`);
    `mean_var` is the unweighted mean of the per-point variances (err^2).
    Returns `(sigma_int, upper_limit)`; `upper_limit` is True whenever the
    excess term is <= 0 (no scatter detectable above measurement noise) or
    fewer than 2 points are available (sigma_int is NaN in that case).
    """
    v = np.asarray(values, dtype=float)[mask]
    e = np.asarray(errs, dtype=float)[mask]
    n = v.size
    if n < 2:
        return float('nan'), True
    if wrap:
        resid = _wrap180(v - mean)
        sigma_obs2 = float(np.sum(resid ** 2) / (n - 1))
    else:
        sigma_obs2 = float(np.var(v, ddof=1))
    mean_var = float(np.mean(e ** 2))
    excess = sigma_obs2 - mean_var
    return float(np.sqrt(max(0.0, excess))), bool(excess <= 0)


def _ptp(values, mask) -> float:
    v = np.asarray(values, dtype=float)[mask]
    return float(np.ptp(v)) if v.size else float('nan')


def _ptp_wrapped(values, mask, mean: float) -> float:
    v = np.asarray(values, dtype=float)[mask]
    if v.size == 0:
        return float('nan')
    resid = _wrap180(v - mean)
    return float(np.max(resid) - np.min(resid))


def structure_function(t_sec, x, lags_s: Sequence[float], half_width_s: float, wrap: bool = False) -> dict:
    """First-order structure function SF(tau) = mean((x_j - x_i)^2) over all
    pairs whose time separation falls within `half_width_s` of each `tau` in
    `lags_s`. NaN for a lag with no qualifying pair. `wrap=True` wraps the
    value difference into (-90, +90] before squaring (for PA).
    """
    t = np.asarray(t_sec, dtype=float)
    v = np.asarray(x, dtype=float)
    finite = np.isfinite(t) & np.isfinite(v)
    t, v = t[finite], v[finite]
    order = np.argsort(t)
    t, v = t[order], v[order]
    n = t.size
    if n < 2:
        return {float(tau): float('nan') for tau in lags_s}
    i_idx, j_idx = np.triu_indices(n, k=1)
    dt = t[j_idx] - t[i_idx]
    dx = v[j_idx] - v[i_idx]
    if wrap:
        dx = _wrap180(dx)
    sf = {}
    for tau in lags_s:
        sel = np.abs(dt - tau) <= half_width_s
        sf[float(tau)] = float(np.mean(dx[sel] ** 2)) if np.any(sel) else float('nan')
    return sf


def _rows_to_arrays(rows_by_segid: dict) -> dict:
    seg_ids = sorted(rows_by_segid)

    def col(name):
        return np.array([rows_by_segid[i][name] for i in seg_ids], dtype=float)

    return {
        'segment_id': np.array(seg_ids, dtype=int),
        't_start_mjd': col('t_start_mjd'), 't_stop_mjd': col('t_stop_mjd'),
        'exposure_s': col('exposure_s'),
        'PD_pct': col('PD_pct'), 'PD_err': col('PD_err'),
        'PA_deg': col('PA_deg'), 'PA_err': col('PA_err'),
        'QI': col('QI'), 'QI_err': col('QI_err'),
        'UI': col('UI'), 'UI_err': col('UI_err'),
        'MDP99': col('MDP99'), 'SIGNIF': col('SIGNIF'),
    }


def compute_obsid_stability(
    provenance: dict,
    rows_by_segid: dict,
    integrated: dict | None,
    variability_p_threshold: float = 0.01,
    sf_lags_s: Sequence[float] = (5000.0, 10000.0, 20000.0),
    min_seg_for_sf: int = 10,
) -> dict:
    """PD/PA/QI/UI constancy statistics for one ObsID's pcube_timeseg table.

    Parameters mirror `pa_pd_stability`. Returns a dict with per-quantity
    weighted-mean/chi2/dof/p sub-dicts ('PD', 'PA', 'QI', 'UI'), the
    authoritative 2D Stokes test ('chi2_2D', 'dof_2D', 'p_2D', 'variable'),
    moment-estimator intrinsic scatter ('sigma_PD_int'/'PD_int_upper_limit',
    'sigma_PA_int'/'PA_int_upper_limit'), peak-to-peak ranges ('PD_ptp',
    'PA_ptp'), reduced chi-squared ('rchi2_PD', 'rchi2_PA'), the mean MDP99
    ('mean_MDP99'), structure functions ('sf_pd', 'sf_pa': lag_s -> value),
    'n_seg', and 'integrated' (the PD/PA/MDP99/SIGNIF integrated-row fields,
    or None if absent).
    """
    arr = _rows_to_arrays(rows_by_segid)
    n_seg = arr['segment_id'].size

    pd_stat = _weighted_stat(arr['PD_pct'], arr['PD_err'])
    pa_stat = _weighted_stat_pa(arr['PA_deg'], arr['PA_err'])

    qi_ui_mask = (np.isfinite(arr['QI']) & np.isfinite(arr['QI_err']) & (arr['QI_err'] > 0)
                  & np.isfinite(arr['UI']) & np.isfinite(arr['UI_err']) & (arr['UI_err'] > 0))
    n2 = int(qi_ui_mask.sum())
    if n2 >= 2:
        w_qi = 1.0 / arr['QI_err'][qi_ui_mask] ** 2
        w_ui = 1.0 / arr['UI_err'][qi_ui_mask] ** 2
        qi_mean = float(np.sum(w_qi * arr['QI'][qi_ui_mask]) / np.sum(w_qi))
        ui_mean = float(np.sum(w_ui * arr['UI'][qi_ui_mask]) / np.sum(w_ui))
        chi2_qi = float(np.sum(w_qi * (arr['QI'][qi_ui_mask] - qi_mean) ** 2))
        chi2_ui = float(np.sum(w_ui * (arr['UI'][qi_ui_mask] - ui_mean) ** 2))
        dof_1d = n2 - 1
        qi_stat = {'mean': qi_mean, 'chi2': chi2_qi, 'dof': dof_1d, 'p': float(stats.chi2.sf(chi2_qi, dof_1d)), 'n': n2}
        ui_stat = {'mean': ui_mean, 'chi2': chi2_ui, 'dof': dof_1d, 'p': float(stats.chi2.sf(chi2_ui, dof_1d)), 'n': n2}
        chi2_2d = chi2_qi + chi2_ui
        dof_2d = 2 * dof_1d
        p_2d = float(stats.chi2.sf(chi2_2d, dof_2d))
        variable = bool(p_2d < variability_p_threshold)
    else:
        qi_stat = {'mean': float('nan'), 'chi2': float('nan'), 'dof': n2 - 1, 'p': float('nan'), 'n': n2}
        ui_stat = dict(qi_stat)
        chi2_2d, dof_2d, p_2d, variable = float('nan'), float('nan'), float('nan'), False

    sigma_pd_int, pd_upper = _moment_intrinsic_scatter(arr['PD_pct'], arr['PD_err'], pd_stat['mask'])
    sigma_pa_int, pa_upper = _moment_intrinsic_scatter(
        arr['PA_deg'], arr['PA_err'], pa_stat['mask'], wrap=True, mean=pa_stat['mean'])

    pd_ptp = _ptp(arr['PD_pct'], pd_stat['mask'])
    pa_ptp = _ptp_wrapped(arr['PA_deg'], pa_stat['mask'], pa_stat['mean'])

    rchi2_pd = pd_stat['chi2'] / pd_stat['dof'] if pd_stat['dof'] > 0 else float('nan')
    rchi2_pa = pa_stat['chi2'] / pa_stat['dof'] if pa_stat['dof'] > 0 else float('nan')

    finite_mdp = np.isfinite(arr['MDP99'])
    mean_mdp99 = float(np.mean(arr['MDP99'][finite_mdp])) if finite_mdp.any() else float('nan')

    sf_pd = {float(tau): float('nan') for tau in sf_lags_s}
    sf_pa = dict(sf_pd)
    if n_seg >= min_seg_for_sf:
        half_width_s = provenance.get('segment_duration_s') or None
        half_width_s = 0.5 * half_width_s if half_width_s else None
        if not half_width_s:
            t_sorted = np.sort(arr['t_start_mjd'])
            dt_mjd = np.diff(t_sorted)
            half_width_s = float(np.median(dt_mjd) * 86400.0 * 0.5) if dt_mjd.size else float('nan')
        t_sec = arr['t_start_mjd'] * 86400.0
        sf_pd = structure_function(t_sec, arr['PD_pct'], sf_lags_s, half_width_s)
        sf_pa = structure_function(t_sec, arr['PA_deg'], sf_lags_s, half_width_s, wrap=True)

    integrated_out = None
    if integrated is not None:
        integrated_out = {k: integrated[k] for k in ('PD_pct', 'PD_err', 'PA_deg', 'PA_err', 'MDP99', 'SIGNIF')}

    return {
        'n_seg': n_seg,
        'integrated': integrated_out,
        'PD': {k: v for k, v in pd_stat.items() if k != 'mask'},
        'PA': {k: v for k, v in pa_stat.items() if k != 'mask'},
        'QI': qi_stat, 'UI': ui_stat,
        'chi2_2D': chi2_2d, 'dof_2D': dof_2d, 'p_2D': p_2d, 'variable': variable,
        'sigma_PD_int': sigma_pd_int, 'PD_int_upper_limit': pd_upper,
        'sigma_PA_int': sigma_pa_int, 'PA_int_upper_limit': pa_upper,
        'PD_ptp': pd_ptp, 'PA_ptp': pa_ptp,
        'rchi2_PD': rchi2_pd, 'rchi2_PA': rchi2_pa,
        'mean_MDP99': mean_mdp99,
        'sf_pd': sf_pd, 'sf_pa': sf_pa,
    }


def _spearman(xs: list[float], ys: list[float]) -> dict:
    pairs = [(x, y) for x, y in zip(xs, ys) if np.isfinite(x) and np.isfinite(y)]
    n = len(pairs)
    if n < 3:
        return {'n': n, 'rho': float('nan'), 'p': float('nan'), 'note': 'insufficient data (n<3)'}
    px = [p[0] for p in pairs]
    py = [p[1] for p in pairs]
    res = stats.spearmanr(px, py)
    return {'n': n, 'rho': float(res.statistic), 'p': float(res.pvalue), 'note': None}


def _catalog_row(obsid: str, source: str, s: dict) -> dict:
    integ = s.get('integrated') or {}
    return {
        'obsid': obsid, 'source': source, 'N_seg': s['n_seg'],
        'PD_pct_int': integ.get('PD_pct', float('nan')), 'PD_err_int': integ.get('PD_err', float('nan')),
        'PA_deg_int': integ.get('PA_deg', float('nan')), 'PA_err_int': integ.get('PA_err', float('nan')),
        'MDP99_int': integ.get('MDP99', float('nan')), 'SIGNIF_int': integ.get('SIGNIF', float('nan')),
        'chi2_PD': s['PD']['chi2'], 'dof_PD': s['PD']['dof'], 'p_PD': s['PD']['p'],
        'chi2_PA': s['PA']['chi2'], 'dof_PA': s['PA']['dof'], 'p_PA': s['PA']['p'],
        'chi2_QI': s['QI']['chi2'], 'dof_QI': s['QI']['dof'], 'p_QI': s['QI']['p'],
        'chi2_UI': s['UI']['chi2'], 'dof_UI': s['UI']['dof'], 'p_UI': s['UI']['p'],
        'chi2_2D': s['chi2_2D'], 'dof_2D': s['dof_2D'], 'p_2D': s['p_2D'], 'variable': s['variable'],
        'sigma_PA_int': s['sigma_PA_int'], 'PA_int_upper_limit': s['PA_int_upper_limit'],
        'sigma_PD_int': s['sigma_PD_int'], 'PD_int_upper_limit': s['PD_int_upper_limit'],
        'PD_ptp': s['PD_ptp'], 'PA_ptp': s['PA_ptp'],
        'rchi2_PD': s['rchi2_PD'], 'rchi2_PA': s['rchi2_PA'],
        'mean_MDP99': s['mean_MDP99'], 'error': '',
    }


def pa_pd_stability(
    obs_db_path: str,
    obsids: list[str],
    output_dir: str,
    sources: dict | None = None,
    variability_p_threshold: float = 0.01,
    sf_lags_ks: list[float] | None = None,
    min_seg_for_sf: int = 10,
) -> dict:
    """Cross-ObsID intra-observation PA/PD stability analysis for IXPE.

    Reads `<obs_db_path>/<obsid>/pcube_timeseg_<obsid>.csv` (the CSV
    `red_ixpe.pcube_timeseg` writes) for every ObsID in `obsids`, computes
    per-ObsID constancy statistics via `compute_obsid_stability`, and writes
    `pa_pd_stability_catalog.csv` (one row per ObsID) and
    `pa_pd_variability_summary.txt` (per-ObsID 2D verdict, secondary stats,
    the two cross-source Spearman results, and the cited source CSV paths)
    to `output_dir`. Does not copy the per-segment input CSVs.

    An ObsID whose CSV is missing, or unreadable (see `_read_timeseg_csv`),
    is recorded as an error entry (in `per_obs` and the catalog's `error`
    column) and skipped -- it never raises and never appears in the
    cross-source Spearman correlation.

    Parameters
    ----------
    obs_db_path : str
        The IXPE mission's obs-db directory.
    obsids : list[str]
        ObsIDs to analyse.
    output_dir : str
        Directory the catalog CSV and summary text file are written to
        (created if absent).
    sources : dict[str, str] or None
        Optional obsid -> source-name label for the catalog/summary; falls
        back to the obsid itself when absent.
    variability_p_threshold : float
        p_2D below this flags an ObsID `variable`. Default 0.01.
    sf_lags_ks : list[float] or None
        Structure-function lags in kiloseconds. Default [5, 10, 20].
    min_seg_for_sf : int
        Minimum per-segment row count before structure functions are
        computed at all (else every lag is NaN). Default 10.

    Returns
    -------
    dict with 'per_obs' (obsid -> `compute_obsid_stability` output, or
    `{'error': ...}`), 'spearman' (the two cross-source correlations),
    'catalog_csv', 'summary_txt' (paths), and 'exit_code' (always 0 -- a
    per-ObsID problem is recorded, not raised).
    """
    os.makedirs(output_dir, exist_ok=True)
    sf_lags_s = [ks * 1000.0 for ks in (sf_lags_ks or _DEFAULT_SF_LAGS_KS)]

    per_obs: dict = {}
    catalog_rows = []
    csv_paths_used = []
    for obsid in obsids:
        source = (sources or {}).get(obsid, obsid)
        csv_path = os.path.join(obs_db_path, obsid, f'pcube_timeseg_{obsid}.csv')
        parsed = _read_timeseg_csv(csv_path)
        if parsed is None:
            per_obs[obsid] = {'error': f'no readable pcube_timeseg CSV at {csv_path}'}
            catalog_rows.append({**{f: float('nan') for f in CATALOG_FIELDS},
                                  'obsid': obsid, 'source': source, 'N_seg': 0, 'variable': False,
                                  'PD_int_upper_limit': True, 'PA_int_upper_limit': True,
                                  'error': per_obs[obsid]['error']})
            continue
        provenance, rows_by_segid, integrated = parsed
        stability = compute_obsid_stability(
            provenance, rows_by_segid, integrated,
            variability_p_threshold=variability_p_threshold,
            sf_lags_s=sf_lags_s, min_seg_for_sf=min_seg_for_sf)
        stability['error'] = None
        stability['source'] = source
        stability['csv_path'] = csv_path
        per_obs[obsid] = stability
        catalog_rows.append(_catalog_row(obsid, source, stability))
        csv_paths_used.append(csv_path)

    catalog_csv = os.path.join(output_dir, 'pa_pd_stability_catalog.csv')
    with open(catalog_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=CATALOG_FIELDS, restval='')
        writer.writeheader()
        writer.writerows(catalog_rows)

    ok = {o: s for o, s in per_obs.items() if s.get('error') is None}
    spearman = {
        'mdp99_vs_rchi2_pd': _spearman([s['mean_MDP99'] for s in ok.values()],
                                        [s['rchi2_PD'] for s in ok.values()]),
        'mdp99_vs_rchi2_pa': _spearman([s['mean_MDP99'] for s in ok.values()],
                                        [s['rchi2_PA'] for s in ok.values()]),
    }

    summary_txt = os.path.join(output_dir, 'pa_pd_variability_summary.txt')
    lines = ['IXPE intra-observation PA/PD stability summary', '=' * 48, '']
    for obsid in obsids:
        s = per_obs[obsid]
        if s.get('error'):
            lines.append(f'{obsid}: ERROR -- {s["error"]}')
            continue
        verdict = 'VARIABLE' if s['variable'] else 'constant'
        lines.append(
            f"{obsid} ({s['source']}): N_seg={s['n_seg']}  "
            f"chi2_2D={s['chi2_2D']:.3f} dof_2D={s['dof_2D']} p_2D={s['p_2D']:.3e} -> {verdict}")
        lines.append(
            f"    PD: chi2/dof={s['PD']['chi2']:.3f}/{s['PD']['dof']} p={s['PD']['p']:.3e}  "
            f"sigma_int={s['sigma_PD_int']:.4f}{' (UL)' if s['PD_int_upper_limit'] else ''}  ptp={s['PD_ptp']:.4f}")
        lines.append(
            f"    PA: chi2/dof={s['PA']['chi2']:.3f}/{s['PA']['dof']} p={s['PA']['p']:.3e}  "
            f"sigma_int={s['sigma_PA_int']:.4f}{' (UL)' if s['PA_int_upper_limit'] else ''}  ptp={s['PA_ptp']:.4f}")
        lines.append(f"    mean_MDP99={s['mean_MDP99']:.4f}  source_csv={s['csv_path']}")
        lines.append('')
    lines.append('Cross-source Spearman correlations (mean MDP99 vs reduced chi-squared):')
    for key, res in spearman.items():
        if res['note']:
            lines.append(f"  {key}: {res['note']} (n={res['n']})")
        else:
            lines.append(f"  {key}: rho={res['rho']:.4f} p={res['p']:.4e} (n={res['n']})")
    lines.append('')
    lines.append('Source data:')
    for p in csv_paths_used:
        lines.append(f'  {p}')
    with open(summary_txt, 'w') as f:
        f.write('\n'.join(lines) + '\n')

    return {
        'per_obs': per_obs, 'spearman': spearman,
        'catalog_csv': catalog_csv, 'summary_txt': summary_txt,
        'exit_code': 0,
    }
