"""Edelson & Krolik (1988) discrete correlation function between IXPE
q(t)=Q/I and u(t)=U/I slice series -- pure numpy/scipy.

Reads the `qu_slices_<subobs_key>.csv` per-subobs series `red_ixpe.
qu_slice_pcubes`/`ixpe_qu_slice_pcubes` writes, and for each subobs computes a
model-independent DCF between mean-subtracted q(t) and u(t) (never
`scipy.signal.correlate`, which assumes uniform sampling this data does not
have), a Monte-Carlo permutation-null significance, a bootstrap peak-lag/peak-
value uncertainty, and a masked-uniform-grid cross-correlation cross-check.
Also reconstructs the two population covariates the CSV does not itself carry
(a source-brightness proxy and a combined PD significance) and correlates
`tau_peak` against them (Spearman) across every processed subobs.

Same module family as `pa_pd_stability`/`bhxrb_pd_hr_survey`: pure numpy/scipy,
no `xspec`/`ixpeobssim` import, so this module -- and the `qu_lag_dcf` entry
point -- runs standalone, without the science stack.

`_read_qu_slices_csv` is an independent reimplementation of the provenance-
header "no header -> treat as absent" convention `pa_pd_stability.
_read_timeseg_csv` and `bhxrb_pd_hr_survey._read_pcube_integrated` already
follow -- not a shared helper, so this module never imports either of them or
`red_ixpe` (which pulls in `ixpeobssim`/`xspec` at module load).
"""

from __future__ import annotations

import csv
import glob
import json
import math
import os
from typing import Sequence

import numpy as np
from scipy import stats

from obsdb_catalog import neutral

_QU_SLICE_HEADER_PREFIX = '# qu_slice_pcubes provenance: '
_QU_DCF_HEADER_PREFIX = '# qu_lag_dcf provenance: '
_QU_DCF_FIELDS = ['lag_s', 'dcf', 'n_pairs', 'perm_lo99', 'perm_hi99']

CCF_MASTER_FIELDS = [
    'obsid', 'source', 'subobs_key', 'n_slices', 'status',
    'tau_peak_s', 'dcf_peak', 'dcf_zero', 'n_pairs_peak', 'dcf_uncorrected_denominator',
    'perm_p', 'tau_boot_lo', 'tau_boot_hi', 'dcf_peak_boot_lo', 'dcf_peak_boot_hi',
    'coherent_rotation', 'r_proxy', 'pd_significance', 'tau_peak_uniform_s',
    'median_dt_s', 'lag_bin_s_used',
]
POPULATION_FIELDS = ['pair', 'n', 'rho', 'p']


# --------------------------------------------------------------------------- #
# Input: qu_slices_<subobs_key>.csv reader
# --------------------------------------------------------------------------- #

def _read_qu_slices_csv(csv_path: str) -> tuple[dict, list[dict]] | None:
    """Read one `qu_slices_<subobs_key>.csv` into `(provenance, rows)`.

    Mirrors `pa_pd_stability._read_timeseg_csv`'s contract: returns None when
    `csv_path` does not exist, or its first line is not the
    `# qu_slice_pcubes provenance: {json}` comment header -- an unrecognised
    file is treated as absent, never parsed as headerless data.
    """
    if not os.path.exists(csv_path):
        return None
    with open(csv_path, newline='') as f:
        header_line = f.readline()
        if not header_line.startswith(_QU_SLICE_HEADER_PREFIX):
            return None
        provenance = json.loads(header_line[len(_QU_SLICE_HEADER_PREFIX):])
        rows = [{
            'slice_id': int(r['slice_id']),
            't_mean_met': float(r['t_mean_met']),
            't_start_met': float(r['t_start_met']),
            't_stop_met': float(r['t_stop_met']),
            'good_exposure_s': float(r['good_exposure_s']),
            'counts': float(r['counts']),
            'q': float(r['q']), 'q_err': float(r['q_err']),
            'u': float(r['u']), 'u_err': float(r['u_err']),
        } for r in csv.DictReader(f)]
    return provenance, rows


# --------------------------------------------------------------------------- #
# Reconstructed population covariates (never persisted in qu_slices CSV)
# --------------------------------------------------------------------------- #

def r_proxy_from_rows(rows: list[dict], src_rad: float) -> float:
    """Source-region count-rate surface brightness (cps/arcmin**2):
    `sum(counts) / (sum(good_exposure_s) * pi * src_rad**2)`, `src_rad` taken
    exactly from the `qu_slices` provenance header. NaN for non-positive total
    exposure (never a `ZeroDivisionError`)."""
    total_counts = sum(r['counts'] for r in rows)
    total_exposure = sum(r['good_exposure_s'] for r in rows)
    if total_exposure <= 0:
        return float('nan')
    return total_counts / (total_exposure * math.pi * src_rad ** 2)


def _weighted_mean_and_err(x: Sequence[float], err: Sequence[float]) -> tuple[float, float]:
    """Inverse-variance weighted mean and its standard error -- the same
    weighting `pa_pd_stability._weighted_stat` uses, reduced here to just the
    mean/sigma (no chi-squared constancy test needed for this covariate)."""
    x = np.asarray(x, dtype=float)
    err = np.asarray(err, dtype=float)
    mask = np.isfinite(x) & np.isfinite(err) & (err > 0)
    if not mask.any():
        return float('nan'), float('nan')
    w = 1.0 / err[mask] ** 2
    mean = float(np.sum(w * x[mask]) / np.sum(w))
    sigma = float(1.0 / math.sqrt(np.sum(w)))
    return mean, sigma


def pd_significance_from_rows(rows: list[dict]) -> float:
    """Inverse-variance-weighted combined normalized-Stokes significance
    PD/sigma_PD, from the per-slice q/u (a numpy weighted combination, not an
    event-level reduction): weighted-mean q, u give PD=sqrt(q_w**2+u_w**2),
    with sigma_PD from first-order error propagation of the two weighted-mean
    standard errors. NaN if PD or sigma_PD is non-positive."""
    q_w, sigma_q = _weighted_mean_and_err([r['q'] for r in rows], [r['q_err'] for r in rows])
    u_w, sigma_u = _weighted_mean_and_err([r['u'] for r in rows], [r['u_err'] for r in rows])
    if not (np.isfinite(q_w) and np.isfinite(u_w)):
        return float('nan')
    pd = math.sqrt(q_w ** 2 + u_w ** 2)
    if pd <= 0:
        return float('nan')
    sigma_pd = math.sqrt((q_w * sigma_q) ** 2 + (u_w * sigma_u) ** 2) / pd
    return pd / sigma_pd if sigma_pd > 0 else float('nan')


# --------------------------------------------------------------------------- #
# Edelson & Krolik (1988) discrete correlation function
# --------------------------------------------------------------------------- #

def _corrected_var(x: np.ndarray, err: np.ndarray) -> float:
    return float(np.var(x) - np.mean(err ** 2))


def denominator(q: np.ndarray, q_err: np.ndarray, u: np.ndarray, u_err: np.ndarray,
                 variance_floor: float) -> tuple[float, bool]:
    """`D = (var(q)-mean(q_err**2)) * (var(u)-mean(u_err**2))`; falls back to
    the uncorrected `D=var(q)*var(u)` (flagging `uncorrected=True`) when
    either measurement-noise-corrected variance is <= `variance_floor`. A
    still-non-positive `D` (a genuinely constant series) is floored at
    `variance_floor` so the DCF stays finite rather than dividing by zero."""
    var_q_corr, var_u_corr = _corrected_var(q, q_err), _corrected_var(u, u_err)
    if var_q_corr <= variance_floor or var_u_corr <= variance_floor:
        D, uncorrected = float(np.var(q) * np.var(u)), True
    else:
        D, uncorrected = float(var_q_corr * var_u_corr), False
    return (D if D > 0 else variance_floor), uncorrected


def _pairwise_udcf(t: np.ndarray, a: np.ndarray, b: np.ndarray, D: float) -> tuple[np.ndarray, np.ndarray]:
    """`UDCF_ij = a_i*b_j/sqrt(D)`, `dt_ij = t_j-t_i`, over every `(i, j)` pair
    (never just `i<j`: the EK-DCF is directional between the two series)."""
    dt = (t[None, :] - t[:, None]).ravel()
    udcf = (np.outer(a, b) / math.sqrt(D)).ravel()
    return dt, udcf


def dcf_bins(width: float, half_range: float) -> tuple[np.ndarray, np.ndarray]:
    """Symmetric bin edges/centers of `width`, spanning `+/-half_range`, with
    one bin centered exactly on lag 0 (so `dcf_zero` is always well-defined)."""
    n_side = max(1, int(math.floor(half_range / width)))
    centers = np.arange(-n_side, n_side + 1) * width
    edges = np.concatenate([centers - width / 2.0, [centers[-1] + width / 2.0]])
    return edges, centers


def bin_udcf(dt: np.ndarray, udcf: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """`numpy.histogram`-bin `udcf` by `dt`: `DCF(bin)=mean(UDCF)`, NaN for an
    empty bin (never `scipy.signal.correlate`)."""
    n_pairs, _ = np.histogram(dt, bins=edges)
    sum_udcf, _ = np.histogram(dt, bins=edges, weights=udcf)
    with np.errstate(invalid='ignore', divide='ignore'):
        dcf = sum_udcf / n_pairs
    return np.where(n_pairs > 0, dcf, np.nan), n_pairs


def compute_dcf(t: np.ndarray, a: np.ndarray, b: np.ndarray, D: float,
                 lag_bin_s: float | None, max_lag_frac: float) -> dict:
    """The full binned DCF for one subobs: bin width `max(lag_bin_s or
    median_dt, median_dt)` (never finer than the median sampling), range
    `+/-max_lag_frac*span`. Returns `centers`, `dcf`, `n_pairs`, `edges`,
    `tau_peak`/`dcf_peak`/`n_pairs_peak` (bin of max `|DCF|`), `dcf_zero`
    (the bin centered on lag 0, NaN if that bin has no pairs), `median_dt`,
    `width`."""
    dt_sorted = np.diff(np.sort(t))
    median_dt = float(np.median(dt_sorted)) if dt_sorted.size else 0.0
    width = max(lag_bin_s, median_dt) if lag_bin_s is not None else median_dt
    if width <= 0:
        width = 1.0
    span = float(np.max(t) - np.min(t)) if t.size > 1 else 0.0
    half_range = max(max_lag_frac * span, width)

    edges, centers = dcf_bins(width, half_range)
    dt_all, udcf_all = _pairwise_udcf(t, a, b, D)
    dcf, n_pairs = bin_udcf(dt_all, udcf_all, edges)

    finite = np.isfinite(dcf)
    if finite.any():
        idx = int(np.nanargmax(np.abs(dcf)))
        tau_peak, dcf_peak, n_pairs_peak = float(centers[idx]), float(dcf[idx]), int(n_pairs[idx])
    else:
        tau_peak, dcf_peak, n_pairs_peak = float('nan'), float('nan'), 0
    zero_idx = int(np.argmin(np.abs(centers)))
    dcf_zero = float(dcf[zero_idx]) if n_pairs[zero_idx] > 0 else float('nan')

    return {'centers': centers, 'dcf': dcf, 'n_pairs': n_pairs, 'edges': edges,
            'tau_peak': tau_peak, 'dcf_peak': dcf_peak, 'n_pairs_peak': n_pairs_peak,
            'dcf_zero': dcf_zero, 'median_dt': median_dt, 'width': width}


def uniform_grid_crosscheck(t: np.ndarray, a: np.ndarray, b: np.ndarray, D: float,
                             median_dt: float) -> float:
    """Masked-uniform-grid cross-check: resample `a`/`b` onto a grid at
    `median_dt` spacing (nearest sample within half a bin; unmatched grid
    points masked), then a plain masked cross-correlation over every integer
    grid lag (never `scipy.signal.correlate`). Returns the lag (seconds) of
    max `|correlation|`, or NaN if the series is too short/degenerate."""
    if median_dt <= 0 or t.size < 2:
        return float('nan')
    tmin, tmax = float(np.min(t)), float(np.max(t))
    n_grid = int(round((tmax - tmin) / median_dt)) + 1
    if n_grid < 2:
        return float('nan')
    grid = tmin + np.arange(n_grid) * median_dt
    half = median_dt / 2.0

    order = np.argsort(t)
    t_sorted, a_sorted, b_sorted = t[order], a[order], b[order]
    a_grid = np.full(n_grid, np.nan)
    b_grid = np.full(n_grid, np.nan)
    pos = np.searchsorted(t_sorted, grid)
    for gi in range(n_grid):
        cands = [k for k in (pos[gi] - 1, pos[gi]) if 0 <= k < t_sorted.size]
        if not cands:
            continue
        best = min(cands, key=lambda k: abs(t_sorted[k] - grid[gi]))
        if abs(t_sorted[best] - grid[gi]) <= half:
            a_grid[gi], b_grid[gi] = a_sorted[best], b_sorted[best]

    best_abs, best_lag = -1.0, 0
    for lag in range(-(n_grid - 1), n_grid):
        if lag >= 0:
            x, y = a_grid[:n_grid - lag], b_grid[lag:n_grid]
        else:
            x, y = a_grid[-lag:n_grid], b_grid[:n_grid + lag]
        mask = np.isfinite(x) & np.isfinite(y)
        if not mask.any():
            continue
        cc = float(np.mean(x[mask] * y[mask]) / math.sqrt(D))
        if abs(cc) > best_abs:
            best_abs, best_lag = abs(cc), lag
    return float(best_lag * median_dt) if best_abs >= 0 else float('nan')


def permutation_null(t: np.ndarray, a: np.ndarray, b: np.ndarray, D: float,
                      edges: np.ndarray, centers: np.ndarray, observed_dcf: np.ndarray,
                      n_permutations: int, rng: np.random.Generator) -> tuple[float, np.ndarray, np.ndarray]:
    """`n_permutations` time-shuffles of `b` (holding `t`/`a` fixed): an
    empirical p-value on max `|DCF|` (add-one smoothed, so it is never
    exactly 0), plus the per-bin 1st/99th percentile envelope across
    permutations (`perm_lo99`/`perm_hi99`)."""
    observed_abs_max = float(np.nanmax(np.abs(observed_dcf))) if np.isfinite(observed_dcf).any() else float('nan')
    perm_dcfs = np.full((n_permutations, centers.size), np.nan)
    count_ge = 0
    for k in range(n_permutations):
        b_perm = rng.permutation(b)
        dt_all, udcf_all = _pairwise_udcf(t, a, b_perm, D)
        dcf_k, _ = bin_udcf(dt_all, udcf_all, edges)
        perm_dcfs[k] = dcf_k
        perm_abs_max = float(np.nanmax(np.abs(dcf_k))) if np.isfinite(dcf_k).any() else -np.inf
        if np.isfinite(observed_abs_max) and perm_abs_max >= observed_abs_max:
            count_ge += 1
    p = (count_ge + 1) / (n_permutations + 1)
    lo99 = np.nanpercentile(perm_dcfs, 1, axis=0)
    hi99 = np.nanpercentile(perm_dcfs, 99, axis=0)
    return p, lo99, hi99


def bootstrap_ci(t: np.ndarray, q: np.ndarray, q_err: np.ndarray, u: np.ndarray, u_err: np.ndarray,
                  edges: np.ndarray, centers: np.ndarray, variance_floor: float,
                  n_bootstrap: int, rng: np.random.Generator) -> tuple[float, float, float, float]:
    """`n_bootstrap` draws of `q_i+N(0,q_err_i)`, `u_i+N(0,u_err_i)`: the 68%
    (16th/84th percentile) interval on the recovered `tau_peak` and
    `dcf_peak`. `D` is recomputed per draw (the perturbed sample's own
    measurement-noise-corrected variance), same fallback rule as `denominator`."""
    taus = np.full(n_bootstrap, np.nan)
    dcfs = np.full(n_bootstrap, np.nan)
    for k in range(n_bootstrap):
        q_b = q + rng.normal(0.0, q_err)
        u_b = u + rng.normal(0.0, u_err)
        a_b, b_b = q_b - np.mean(q_b), u_b - np.mean(u_b)
        D_b, _ = denominator(q_b, q_err, u_b, u_err, variance_floor)
        dt_all, udcf_all = _pairwise_udcf(t, a_b, b_b, D_b)
        dcf_k, _ = bin_udcf(dt_all, udcf_all, edges)
        if np.isfinite(dcf_k).any():
            idx = int(np.nanargmax(np.abs(dcf_k)))
            taus[k], dcfs[k] = centers[idx], dcf_k[idx]
    tau_lo, tau_hi = (float(v) for v in np.nanpercentile(taus, [16, 84]))
    dcf_lo, dcf_hi = (float(v) for v in np.nanpercentile(dcfs, [16, 84]))
    return tau_lo, tau_hi, dcf_lo, dcf_hi


# --------------------------------------------------------------------------- #
# Per-subobs orchestration
# --------------------------------------------------------------------------- #

def compute_subobs_dcf(rows: list[dict], src_rad: float, lag_bin_s: float | None,
                        n_permutations: int, n_bootstrap: int, variance_floor: float,
                        coherent_p_threshold: float, max_lag_frac: float,
                        rng: np.random.Generator) -> dict:
    """One subobs's full DCF/permutation/bootstrap/covariate analysis, from
    its already-parsed `qu_slices` rows. Assumes `len(rows) >= 2`."""
    t = np.array([r['t_mean_met'] for r in rows], dtype=float)
    q = np.array([r['q'] for r in rows], dtype=float)
    q_err = np.array([r['q_err'] for r in rows], dtype=float)
    u = np.array([r['u'] for r in rows], dtype=float)
    u_err = np.array([r['u_err'] for r in rows], dtype=float)
    a, b = q - np.mean(q), u - np.mean(u)

    D, uncorrected = denominator(q, q_err, u, u_err, variance_floor)
    dcf_result = compute_dcf(t, a, b, D, lag_bin_s, max_lag_frac)
    perm_p, perm_lo99, perm_hi99 = permutation_null(
        t, a, b, D, dcf_result['edges'], dcf_result['centers'], dcf_result['dcf'],
        n_permutations, rng)
    tau_lo, tau_hi, dcf_lo, dcf_hi = bootstrap_ci(
        t, q, q_err, u, u_err, dcf_result['edges'], dcf_result['centers'],
        variance_floor, n_bootstrap, rng)
    tau_peak_uniform = uniform_grid_crosscheck(t, a, b, D, dcf_result['median_dt'])

    coherent = bool(perm_p < coherent_p_threshold and (tau_lo > 0 or tau_hi < 0))

    return {
        'centers': dcf_result['centers'], 'dcf': dcf_result['dcf'], 'n_pairs': dcf_result['n_pairs'],
        'perm_lo99': perm_lo99, 'perm_hi99': perm_hi99,
        'tau_peak_s': dcf_result['tau_peak'], 'dcf_peak': dcf_result['dcf_peak'],
        'dcf_zero': dcf_result['dcf_zero'], 'n_pairs_peak': dcf_result['n_pairs_peak'],
        'dcf_uncorrected_denominator': uncorrected, 'perm_p': perm_p,
        'tau_boot_lo': tau_lo, 'tau_boot_hi': tau_hi,
        'dcf_peak_boot_lo': dcf_lo, 'dcf_peak_boot_hi': dcf_hi,
        'coherent_rotation': coherent, 'tau_peak_uniform_s': tau_peak_uniform,
        'median_dt_s': dcf_result['median_dt'], 'lag_bin_s_used': dcf_result['width'],
        'r_proxy': r_proxy_from_rows(rows, src_rad), 'pd_significance': pd_significance_from_rows(rows),
    }


def _write_csv(path: str, fieldnames: list[str], rows: list[dict]) -> None:
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, restval='')
        writer.writeheader()
        writer.writerows(rows)


def _spearman(xs: Sequence[float], ys: Sequence[float]) -> dict:
    x, y = np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < 3:
        return {'n': n, 'rho': float('nan'), 'p': float('nan')}
    res = stats.spearmanr(x[mask], y[mask])
    return {'n': n, 'rho': float(res.statistic), 'p': float(res.pvalue)}


def qu_lag_dcf(
    obs_db_path: str,
    obsids: list[str],
    output_dir: str,
    sources: dict[str, str] | None = None,
    lag_bin_s: float | None = None,
    n_min_slices: int = 8,
    n_permutations: int = 1000,
    n_bootstrap: int = 1000,
    variance_floor: float = 1e-12,
    coherent_p_threshold: float = 0.01,
    max_lag_frac: float = 0.5,
    exclude_noise_limited: bool = True,
    seed: int = 0,
) -> dict:
    """Edelson & Krolik (1988) q(t)-vs-u(t) discrete correlation function
    across every subobs of every requested ObsID.

    For each `obsid` in `obsids`, globs `<obs_db_path>/<obsid>/qu_slices_*.csv`
    (each file one subobs series, written by `ixpe_qu_slice_pcubes`). An
    `obsid` with NO qu_slices files on disk at all is recorded with
    `subobs_key=obsid`, `status='no_qu_slices'` -- a failed/not-yet-run
    upstream `ixpe_qu_slice_pcubes` step is surfaced explicitly, never
    silently dropped from `ccf_master.csv`/`per_subobs`. A qu_slices file that
    exists but fails the provenance-header check (`_read_qu_slices_csv`
    returning `None`) is recorded `status='unreadable_qu_slices'`, same
    visibility guarantee. For every genuinely-parsed subobs: reconstructs
    `noise_limited` as `n_rows < provenance['n_min']` (the original slicing
    tool's own threshold) and sets `status='insufficient_slices'` (excluded
    from the DCF entirely, never raised) when `n_rows < n_min_slices` (this
    tool's own threshold); else computes the full binned DCF (`compute_dcf`),
    a permutation null (`permutation_null`), a bootstrap peak-lag/peak-value
    interval (`bootstrap_ci`), a masked-uniform-grid cross-check
    (`uniform_grid_crosscheck`), and the two reconstructed covariates
    (`r_proxy_from_rows`, `pd_significance_from_rows`). `coherent_rotation` is
    `perm_p < coherent_p_threshold` AND the bootstrap `tau_peak` 68% interval
    excluding 0.

    Every non-`'ok'` status (`no_qu_slices`, `unreadable_qu_slices`,
    `insufficient_slices`) still gets one row in `ccf_master.csv` and one
    `per_subobs` entry (`csv_path=None`, all-NaN numeric fields) -- nothing
    requested ever vanishes from the output -- and is excluded from the
    population Spearman correlations below (which only ever see `status
    == 'ok'` rows).

    Writes, under `<obs_db_path>/<obsid>/`, one provenance-headed
    `qu_dcf_<subobs_key>.csv` per successfully-processed subobs (`lag_s, dcf,
    n_pairs, perm_lo99, perm_hi99`), recorded at level 3
    (`kind='qu_dcf'`) via `obsdb_catalog.neutral.advance`. Writes, under
    `output_dir`: `ccf_master.csv` (one row per subobs, `CCF_MASTER_FIELDS`)
    and `population_correlations.csv` (Spearman `tau_peak_s` vs `r_proxy` and
    vs `pd_significance`, over `status=='ok'` subobs, dropping
    `noise_limited` ones first when `exclude_noise_limited`).

    All random draws (permutation shuffles, bootstrap noise) are drawn from
    one `numpy.random.default_rng(seed)`, consumed in subobs-glob order, so
    two calls with the same `seed` over the same inputs reproduce identically.

    Parameters
    ----------
    obs_db_path : str
        The IXPE mission's obs-db directory; each ObsID's subobs CSVs are
        read from `<obs_db_path>/<obsid>/qu_slices_*.csv`.
    obsids : list[str]
        ObsIDs to process.
    output_dir : str
        Directory `ccf_master.csv`/`population_correlations.csv` are written
        to (created if absent).
    sources : dict[str, str] or None
        Optional obsid -> source-name label for `ccf_master.csv`; falls back
        to the obsid itself when absent.
    lag_bin_s : float or None
        DCF bin width in seconds. `None` (default) uses each subobs's own
        median inter-slice spacing; an explicit value is still floored at
        that median (never finer than the data actually resolves).
    n_min_slices : int
        Minimum slice count a subobs needs before its DCF is computed at all
        (below this: `status='insufficient_slices'`, excluded). Default 8.
    n_permutations, n_bootstrap : int
        Monte-Carlo draw counts for the permutation null / bootstrap CI.
        Defaults 1000, 1000.
    variance_floor : float
        Threshold below which a measurement-noise-corrected variance is
        treated as degenerate, triggering the uncorrected-denominator
        fallback. Default 1e-12.
    coherent_p_threshold : float
        `perm_p` below this (and a zero-excluding bootstrap `tau_peak` 68%
        interval) flags `coherent_rotation`. Default 0.01.
    max_lag_frac : float
        DCF lag range as a fraction of the subobs's own time span
        (`+/-max_lag_frac*span`). Default 0.5.
    exclude_noise_limited : bool
        Drop `noise_limited` subobs from the population Spearman
        correlations (still appear in `ccf_master.csv`). Default True.
    seed : int
        Seed for every permutation/bootstrap random draw. Default 0.

    Returns
    -------
    dict with `ccf_master_csv`, `population_csv` (paths), `per_subobs`
    (`subobs_key` -> a bounded summary: `status` (`'ok'`, `'no_qu_slices'`,
    `'unreadable_qu_slices'`, or `'insufficient_slices'`), `n_slices`,
    `tau_peak_s`, `dcf_peak`, `perm_p`, `coherent_rotation`, `csv_path` or
    `None` for any non-`'ok'` status), `population_correlations` (the two
    `_spearman` results, computed over `status=='ok'` rows only), `stdout`,
    `exit_code` (always 0 -- a per-obsid/per-subobs problem is recorded, not
    raised or silently dropped).
    """
    os.makedirs(output_dir, exist_ok=True)
    rng = np.random.default_rng(seed)

    master_rows: list[dict] = []
    per_subobs: dict = {}
    stdout_lines = []

    for obsid in obsids:
        source = (sources or {}).get(obsid, obsid)
        csv_paths = sorted(glob.glob(os.path.join(obs_db_path, obsid, 'qu_slices_*.csv')))
        if not csv_paths:
            base_row = {f: float('nan') for f in CCF_MASTER_FIELDS}
            base_row.update(obsid=obsid, source=source, subobs_key=obsid, n_slices=0,
                             status='no_qu_slices', coherent_rotation=False,
                             dcf_uncorrected_denominator=False)
            master_rows.append(base_row)
            per_subobs[obsid] = {'status': 'no_qu_slices', 'n_slices': 0, 'csv_path': None}
            stdout_lines.append(f'{obsid}: no_qu_slices (no qu_slices_*.csv under '
                                 f'{os.path.join(obs_db_path, obsid)})')
            continue
        for slice_csv in csv_paths:
            base = os.path.basename(slice_csv)
            subobs_key = base[len('qu_slices_'):-len('.csv')]
            parsed = _read_qu_slices_csv(slice_csv)
            if parsed is None:
                base_row = {f: float('nan') for f in CCF_MASTER_FIELDS}
                base_row.update(obsid=obsid, source=source, subobs_key=subobs_key, n_slices=0,
                                 status='unreadable_qu_slices', coherent_rotation=False,
                                 dcf_uncorrected_denominator=False)
                master_rows.append(base_row)
                per_subobs[subobs_key] = {'status': 'unreadable_qu_slices', 'n_slices': 0, 'csv_path': None}
                stdout_lines.append(f'{subobs_key}: unreadable_qu_slices (no provenance header at {slice_csv})')
                continue
            provenance, rows = parsed
            n_slices = len(rows)
            noise_limited = n_slices < int(provenance['n_min'])

            base_row = {f: float('nan') for f in CCF_MASTER_FIELDS}
            base_row.update(obsid=obsid, source=source, subobs_key=subobs_key,
                             n_slices=n_slices, coherent_rotation=False,
                             dcf_uncorrected_denominator=False)

            if n_slices < n_min_slices:
                base_row['status'] = 'insufficient_slices'
                master_rows.append(base_row)
                per_subobs[subobs_key] = {'status': 'insufficient_slices', 'n_slices': n_slices,
                                           'csv_path': None}
                continue

            result = compute_subobs_dcf(
                rows, float(provenance['src_rad']), lag_bin_s, n_permutations, n_bootstrap,
                variance_floor, coherent_p_threshold, max_lag_frac, rng)

            dcf_csv = os.path.join(obs_db_path, obsid, f'qu_dcf_{subobs_key}.csv')
            dcf_provenance = {
                'obsid': obsid, 'subobs_key': subobs_key, 'n_slices': n_slices, 'seed': seed,
                'n_permutations': n_permutations, 'n_bootstrap': n_bootstrap,
                'variance_floor': variance_floor, 'coherent_p_threshold': coherent_p_threshold,
                'max_lag_frac': max_lag_frac, 'lag_bin_s': lag_bin_s,
                'lag_bin_s_used': result['lag_bin_s_used'],
                'dcf_uncorrected_denominator': result['dcf_uncorrected_denominator'],
                'source_csv': slice_csv,
            }
            dcf_rows = [{
                'lag_s': float(result['centers'][i]), 'dcf': float(result['dcf'][i]),
                'n_pairs': int(result['n_pairs'][i]),
                'perm_lo99': float(result['perm_lo99'][i]), 'perm_hi99': float(result['perm_hi99'][i]),
            } for i in range(result['centers'].size)]
            with open(dcf_csv, 'w', newline='') as f:
                f.write(_QU_DCF_HEADER_PREFIX + json.dumps(dcf_provenance, sort_keys=True) + '\n')
                writer = csv.DictWriter(f, fieldnames=_QU_DCF_FIELDS)
                writer.writeheader()
                writer.writerows(dcf_rows)

            neutral.advance(obs_db_path, 'ixpe', obsid, level=3, products=[{
                'kind': 'qu_dcf',
                'path': os.path.relpath(dcf_csv, obs_db_path),
                'produced_by': 'xspec_utils.lib.qu_lag_dcf.qu_lag_dcf',
                'params': {'subobs_key': subobs_key, 'n_slices': n_slices,
                           'tau_peak_s': result['tau_peak_s'], 'perm_p': result['perm_p'],
                           'coherent_rotation': result['coherent_rotation']},
            }])

            base_row.update(status='ok', tau_peak_s=result['tau_peak_s'], dcf_peak=result['dcf_peak'],
                             dcf_zero=result['dcf_zero'], n_pairs_peak=result['n_pairs_peak'],
                             dcf_uncorrected_denominator=result['dcf_uncorrected_denominator'],
                             perm_p=result['perm_p'], tau_boot_lo=result['tau_boot_lo'],
                             tau_boot_hi=result['tau_boot_hi'], dcf_peak_boot_lo=result['dcf_peak_boot_lo'],
                             dcf_peak_boot_hi=result['dcf_peak_boot_hi'],
                             coherent_rotation=result['coherent_rotation'], r_proxy=result['r_proxy'],
                             pd_significance=result['pd_significance'],
                             tau_peak_uniform_s=result['tau_peak_uniform_s'],
                             median_dt_s=result['median_dt_s'], lag_bin_s_used=result['lag_bin_s_used'])
            base_row['noise_limited'] = noise_limited
            master_rows.append(base_row)
            per_subobs[subobs_key] = {
                'status': 'ok', 'n_slices': n_slices, 'tau_peak_s': result['tau_peak_s'],
                'dcf_peak': result['dcf_peak'], 'perm_p': result['perm_p'],
                'coherent_rotation': result['coherent_rotation'], 'noise_limited': noise_limited,
                'csv_path': dcf_csv,
            }

    ccf_master_csv = os.path.join(output_dir, 'ccf_master.csv')
    _write_csv(ccf_master_csv, CCF_MASTER_FIELDS,
               [{k: r[k] for k in CCF_MASTER_FIELDS} for r in master_rows])

    ok_rows = [r for r in master_rows if r['status'] == 'ok']
    if exclude_noise_limited:
        ok_rows = [r for r in ok_rows if not r.get('noise_limited')]
    pop = {
        'tau_peak_vs_r_proxy': _spearman([r['tau_peak_s'] for r in ok_rows], [r['r_proxy'] for r in ok_rows]),
        'tau_peak_vs_pd_significance': _spearman([r['tau_peak_s'] for r in ok_rows],
                                                   [r['pd_significance'] for r in ok_rows]),
    }
    population_csv = os.path.join(output_dir, 'population_correlations.csv')
    _write_csv(population_csv, POPULATION_FIELDS,
               [{'pair': k, **v} for k, v in pop.items()])

    return {
        'ccf_master_csv': ccf_master_csv, 'population_csv': population_csv,
        'per_subobs': per_subobs, 'population_correlations': pop,
        'stdout': '\n'.join(stdout_lines), 'exit_code': 0,
    }
