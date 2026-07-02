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

`stability_plot` (bottom of this module) is the pure-matplotlib companion
plotter: it renders the four standard diagnostic figures from a
`pa_pd_stability` catalog plus the underlying `pcube_timeseg` CSVs. Only
matplotlib (Agg backend) -- no xspec/ixpeobssim -- so it stays in `lib/` and
is wired directly into `TOOL_MANIFEST`, not through `skills/session.py`.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import time as _time
from pathlib import Path
from typing import Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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


# --------------------------------------------------------------------------- #
# stability_plot -- pure-matplotlib companion plotter (no xspec/ixpeobssim)
# --------------------------------------------------------------------------- #

_MIME_BY_FMT = {'png': 'image/png', 'pdf': 'application/pdf'}

_PLOT_REQUIRED_CATALOG_FIELDS = (
    'obsid',
    'PA_deg_int', 'PA_err_int', 'PD_pct_int', 'PD_err_int',
    'rchi2_PD', 'rchi2_PA',
)


def _safe_float(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return float('nan')


def _load_catalog_rows(catalog_csv: str) -> dict:
    """Read `pa_pd_stability_catalog.csv` keyed on its `obsid` column.

    Validates up front that the header contains every field `stability_plot`
    later reads (`_PLOT_REQUIRED_CATALOG_FIELDS`) and raises KeyError naming
    the missing columns if not -- so a malformed / wrong-schema / header-only
    (or completely empty) catalog is rejected before any figure is drawn,
    rather than degrading to an empty catalog or raising mid-plot after
    files are already written. Fields otherwise follow `CATALOG_FIELDS`; a
    plain header+rows CSV (unlike the provenance-commented `pcube_timeseg`
    CSVs), so a bare `csv.DictReader` suffices. Blank/error cells convert to
    NaN downstream via `_safe_float`.
    """
    rows = {}
    with open(catalog_csv, newline='') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        missing = [c for c in _PLOT_REQUIRED_CATALOG_FIELDS if c not in fieldnames]
        if missing:
            raise KeyError(
                f'catalog_csv missing required column(s) {missing}; '
                f'header was {fieldnames}')
        for r in reader:
            rows[r['obsid']] = r
    return rows


def mid_time_structure_functions(arr: dict, lags_s: Sequence[float]) -> tuple[dict, dict]:
    """PA/PD structure functions vs lag from one ObsID's per-segment arrays.

    The convention `stability_plot`'s structure-function panel draws: times
    are mid-segment (`0.5*(t_start_mjd + t_stop_mjd)`, converted to seconds),
    and the pair-matching tolerance half-width is half the median spacing
    between an ObsID's own sorted mid-times -- independent of any
    `pcube_timeseg` provenance. Delegates the actual pairing/averaging to
    `structure_function` (this module); never reimplements it.
    """
    mid_s = 0.5 * (arr['t_start_mjd'] + arr['t_stop_mjd']) * 86400.0
    dt = np.diff(np.sort(mid_s))
    half_width_s = float(np.median(dt) * 0.5) if dt.size else float('nan')
    sf_pa = structure_function(mid_s, arr['PA_deg'], lags_s, half_width_s, wrap=True)
    sf_pd = structure_function(mid_s, arr['PD_pct'], lags_s, half_width_s, wrap=False)
    return sf_pa, sf_pd


def stability_plot(
    obs_db_path: str,
    obsids: list[str],
    catalog_csv: str,
    save_path_prefix: str | None = None,
    sf_lags_ks: list[float] | None = None,
    sf_min_segments: int = 10,
    fmt: str = 'png',
) -> dict:
    """Render the four standard IXPE PA/PD stability diagnostic figures.

    Reads each requested ObsID's `pcube_timeseg_<obsid>.csv` (written by
    `ixpe_pcube_timeseg`) and the `pa_pd_stability_catalog.csv` (written by
    `pa_pd_stability`, fields per `CATALOG_FIELDS`), and writes exactly four
    files, `{save_path_prefix}_<suffix>.{fmt}`:

    - `_pa_pd_timeseries`: per-ObsID PA (top) / PD (bottom) vs time-elapsed
      (kiloseconds from the earliest segment's `t_start_mjd`), with the
      catalog's weighted-mean +/-1sigma band overlaid when the ObsID is in
      `catalog_csv`.
    - `_stokes_plane`: per-ObsID Q/I vs U/I scatter, color-coded by
      mid-segment elapsed time, with 1sigma error crosses and the origin
      marked.
    - `_chi2_hist`: histograms of `rchi2_PD`/`rchi2_PA` across every row of
      `catalog_csv` (not just `obsids`), with a reference line at 1.0.
    - `_structure_function`: SF_PA (top) / SF_PD (bottom) vs lag, one line
      per ObsID with >= `sf_min_segments` segments, computed via
      `mid_time_structure_functions` (which itself calls the shared
      `structure_function` primitive -- never an independent computation).

    An ObsID with no readable `pcube_timeseg` CSV is skipped (noted in
    `stdout`); if every requested ObsID is unreadable, nothing is written and
    `exit_code=1` with `stderr` set. `catalog_csv` is a required upstream
    product (from `pa_pd_stability`): if it is missing, or exists but cannot
    be read as a catalog (e.g. lacks an `obsid` column), this is a hard input
    error -- `exit_code=1` with `stderr` set and no files written, regardless
    of whether any `pcube_timeseg` CSV is readable. This is distinct from an
    ObsID simply absent from an otherwise-present, readable catalog, which
    still gets its time-series/Stokes-plane panels -- only its weighted-mean
    overlay is skipped. An ObsID short of `sf_min_segments` is omitted from
    the structure-function panel only; if none qualify, that panel carries an
    explicit annotation instead of empty axes.

    Parameters
    ----------
    obs_db_path : str
        The IXPE mission's obs-db directory; each ObsID's per-segment CSV is
        read from `<obs_db_path>/<obsid>/pcube_timeseg_<obsid>.csv`.
    obsids : list[str]
        ObsIDs to draw the time-series/Stokes-plane/structure-function
        panels for.
    catalog_csv : str
        Path to the `pa_pd_stability_catalog.csv` written by
        `pa_pd_stability`.
    save_path_prefix : str or None
        Base path; the four suffixes above plus `.{fmt}` are appended.
        Auto-generated timestamped prefix in `os.getcwd()` when None.
    sf_lags_ks : list[float] or None
        Structure-function lags in kiloseconds. Default [5, 10, 20].
    sf_min_segments : int
        Minimum per-segment row count an ObsID needs to appear in the
        structure-function panel. Default 10.
    fmt : str
        Output image format, 'png' or 'pdf'. Default 'png'.

    Returns
    -------
    dict with 'figures' (list of 4 FigureRef dicts), 'structure_function'
    (obsid -> {'lags_ks', 'sf_pa', 'sf_pd'} for the ObsIDs actually drawn in
    the structure-function panel -- the exact values plotted, for
    differential verification), 'stdout', 'exit_code' (1 with 'stderr' and
    no files written when `catalog_csv` is missing/unreadable, or when every
    requested ObsID is unreadable).
    """
    if fmt not in _MIME_BY_FMT:
        raise ValueError(f"fmt must be one of {sorted(_MIME_BY_FMT)}, got {fmt!r}")

    # catalog_csv is a required upstream product (from pa_pd_stability) -- a
    # missing or malformed file is a hard input error, not an empty catalog,
    # regardless of whether any pcube_timeseg CSV below is readable.
    if not os.path.exists(catalog_csv):
        return {
            'exit_code': 1, 'stdout': '',
            'stderr': (f'catalog_csv not found: {catalog_csv}. It is a required '
                       f'ixpe_pa_pd_stability product; run that tool first.'),
        }
    try:
        cat_rows = _load_catalog_rows(catalog_csv)
    except (OSError, csv.Error, KeyError) as exc:
        return {
            'exit_code': 1, 'stdout': '',
            'stderr': f'catalog_csv unreadable ({catalog_csv}): {exc!r}',
        }

    lags_ks = list(sf_lags_ks) if sf_lags_ks else list(_DEFAULT_SF_LAGS_KS)
    lags_s = [ks * 1000.0 for ks in lags_ks]

    per_obs: dict = {}
    notes = []
    for obsid in obsids:
        csv_path = os.path.join(obs_db_path, obsid, f'pcube_timeseg_{obsid}.csv')
        parsed = _read_timeseg_csv(csv_path)
        if parsed is None:
            notes.append(f'{obsid}: skipped (no readable pcube_timeseg CSV at {csv_path})')
            continue
        _provenance, rows_by_segid, _integrated = parsed
        arr = _rows_to_arrays(rows_by_segid)
        t0_s = float(np.min(arr['t_start_mjd'])) * 86400.0 if arr['segment_id'].size else 0.0
        mid_s = 0.5 * (arr['t_start_mjd'] + arr['t_stop_mjd']) * 86400.0
        per_obs[obsid] = {'arr': arr, 'x_ks': (mid_s - t0_s) / 1000.0}

    if not per_obs:
        return {
            'exit_code': 1, 'stdout': '',
            'stderr': f'No readable pcube_timeseg CSV for any of {obsids} under {obs_db_path}.',
        }

    if save_path_prefix is None:
        ts = _time.strftime('%Y%m%dT%H%M%S')
        save_path_prefix = os.path.join(os.getcwd(), f'stability_plot_{ts}')
    os.makedirs(os.path.dirname(os.path.abspath(save_path_prefix)) or '.', exist_ok=True)

    mime = _MIME_BY_FMT[fmt]
    savefig_kwargs = ({'metadata': {'Software': 'xspec_utils.lib.pa_pd_stability.stability_plot'}}
                       if fmt == 'png' else
                       {'metadata': {'CreationDate': None,
                                     'Producer': 'xspec_utils.lib.pa_pd_stability.stability_plot'}})

    def _save(fig, suffix: str, caption: str) -> dict:
        path = f'{save_path_prefix}_{suffix}.{fmt}'
        fig.savefig(path, dpi=150, bbox_inches='tight', **savefig_kwargs)
        plt.close(fig)
        sha = hashlib.sha256(Path(path).read_bytes()).hexdigest()
        return {
            'path': path, 'caption': caption, 'mime_type': mime,
            'produced_by': 'xspec_utils.lib.pa_pd_stability.stability_plot',
            'created_at': _time.strftime('%Y-%m-%dT%H:%M:%SZ', _time.gmtime()),
            'sha256': sha,
        }

    plotted = [o for o in obsids if o in per_obs]
    n = len(plotted)
    figures = []

    # -- 1. PA/PD time series ------------------------------------------------
    fig1, axes1 = plt.subplots(nrows=2, ncols=n, figsize=(5 * n, 6),
                                sharex='col', squeeze=False)
    for i, obsid in enumerate(plotted):
        arr, x_ks = per_obs[obsid]['arr'], per_obs[obsid]['x_ks']
        ax_pa, ax_pd = axes1[0, i], axes1[1, i]
        ax_pa.errorbar(x_ks, arr['PA_deg'], yerr=arr['PA_err'], fmt='o', ms=3,
                        elinewidth=0.7, capsize=2, color='tab:blue')
        ax_pd.errorbar(x_ks, arr['PD_pct'], yerr=arr['PD_err'], fmt='o', ms=3,
                        elinewidth=0.7, capsize=2, color='tab:orange')
        cat = cat_rows.get(obsid)
        if cat is not None:
            pa_int, pa_err_int = _safe_float(cat['PA_deg_int']), _safe_float(cat['PA_err_int'])
            pd_int, pd_err_int = _safe_float(cat['PD_pct_int']), _safe_float(cat['PD_err_int'])
            if np.isfinite(pa_int):
                ax_pa.axhline(pa_int, color='k', lw=1, ls='--')
                if np.isfinite(pa_err_int):
                    ax_pa.axhspan(pa_int - pa_err_int, pa_int + pa_err_int, color='k', alpha=0.15)
            if np.isfinite(pd_int):
                ax_pd.axhline(pd_int, color='k', lw=1, ls='--')
                if np.isfinite(pd_err_int):
                    ax_pd.axhspan(pd_int - pd_err_int, pd_int + pd_err_int, color='k', alpha=0.15)
        ax_pa.set_title(obsid)
        ax_pa.set_ylabel('PA (deg)' if i == 0 else '')
        ax_pd.set_ylabel('PD (%)' if i == 0 else '')
        ax_pd.set_xlabel('Time elapsed (ks)')
    fig1.tight_layout()
    figures.append(_save(fig1, 'pa_pd_timeseries',
                          f'PA/PD time series for {", ".join(plotted)}'))

    # -- 2. Stokes plane ------------------------------------------------------
    fig2, axes2 = plt.subplots(ncols=n, figsize=(5 * n, 4.5), squeeze=False)
    for i, obsid in enumerate(plotted):
        arr, x_ks = per_obs[obsid]['arr'], per_obs[obsid]['x_ks']
        ax = axes2[0, i]
        ax.errorbar(arr['QI'], arr['UI'], xerr=arr['QI_err'], yerr=arr['UI_err'],
                    fmt='none', ecolor='gray', elinewidth=0.5, alpha=0.5, zorder=1)
        sc = ax.scatter(arr['QI'], arr['UI'], c=x_ks, cmap='viridis', zorder=2, s=20)
        fig2.colorbar(sc, ax=ax, label='Time elapsed (ks)')
        ax.axhline(0.0, color='k', lw=0.5)
        ax.axvline(0.0, color='k', lw=0.5)
        ax.set_xlabel('Q/I')
        ax.set_ylabel('U/I' if i == 0 else '')
        ax.set_title(obsid)
        ax.set_aspect('equal', adjustable='datalim')
    fig2.tight_layout()
    figures.append(_save(fig2, 'stokes_plane',
                          f'Stokes-plane (Q/I vs U/I) evolution for {", ".join(plotted)}'))

    # -- 3. Reduced chi-squared histograms (every catalog row) ----------------
    rchi2_pd_all = np.array([_safe_float(r['rchi2_PD']) for r in cat_rows.values()])
    rchi2_pa_all = np.array([_safe_float(r['rchi2_PA']) for r in cat_rows.values()])
    rchi2_pd_all = rchi2_pd_all[np.isfinite(rchi2_pd_all)]
    rchi2_pa_all = rchi2_pa_all[np.isfinite(rchi2_pa_all)]
    fig3, (ax_hpd, ax_hpa) = plt.subplots(nrows=2, ncols=1, figsize=(6, 8))
    ax_hpd.hist(rchi2_pd_all, bins=20, color='tab:orange', edgecolor='k', alpha=0.8)
    ax_hpd.axvline(1.0, color='k', ls='--')
    ax_hpd.set_xlabel(r'$\chi^2_\nu$(PD)')
    ax_hpd.set_ylabel('N obsids')
    ax_hpa.hist(rchi2_pa_all, bins=20, color='tab:blue', edgecolor='k', alpha=0.8)
    ax_hpa.axvline(1.0, color='k', ls='--')
    ax_hpa.set_xlabel(r'$\chi^2_\nu$(PA)')
    ax_hpa.set_ylabel('N obsids')
    fig3.tight_layout()
    figures.append(_save(fig3, 'chi2_hist',
                          f'Reduced chi-squared histograms across {len(cat_rows)} catalog ObsID(s)'))

    # -- 4. Structure functions -------------------------------------------------
    sf_out: dict = {}
    sf_eligible = [o for o in plotted if per_obs[o]['arr']['segment_id'].size >= sf_min_segments]
    fig4, (ax_sfpa, ax_sfpd) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(6, 8))
    if not sf_eligible:
        ax_sfpa.annotate(f'No ObsID has >= sf_min_segments={sf_min_segments} segments',
                          xy=(0.5, 0.5), xycoords='axes fraction', ha='center', va='center')
        ax_sfpa.set_xticks([]); ax_sfpa.set_yticks([])
        ax_sfpd.set_xticks([]); ax_sfpd.set_yticks([])
    else:
        for obsid in sf_eligible:
            arr = per_obs[obsid]['arr']
            sf_pa, sf_pd = mid_time_structure_functions(arr, lags_s)
            sf_out[obsid] = {
                'lags_ks': lags_ks,
                'sf_pa': {lag: sf_pa[lag] for lag in lags_s},
                'sf_pd': {lag: sf_pd[lag] for lag in lags_s},
            }
            ax_sfpa.plot(lags_ks, [sf_pa[lag] for lag in lags_s], marker='o', label=obsid)
            ax_sfpd.plot(lags_ks, [sf_pd[lag] for lag in lags_s], marker='o', label=obsid)
        ax_sfpa.legend(fontsize='small')
    ax_sfpa.set_ylabel(r'SF$_{PA}$ (deg$^2$)')
    ax_sfpd.set_ylabel(r'SF$_{PD}$ (%$^2$)')
    ax_sfpd.set_xlabel('Lag (ks)')
    fig4.tight_layout()
    figures.append(_save(fig4, 'structure_function',
                          f'PA/PD structure functions for {", ".join(sf_eligible) or "no eligible ObsID"}'))

    stdout_lines = [f'Wrote {len(figures)} figure(s) with prefix {save_path_prefix} ({fmt}).']
    stdout_lines.extend(notes)
    return {
        'figures': figures,
        'structure_function': sf_out,
        'stdout': '\n'.join(stdout_lines),
        'exit_code': 0,
    }
