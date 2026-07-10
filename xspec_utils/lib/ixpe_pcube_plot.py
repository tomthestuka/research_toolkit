"""Energy-resolved IXPE PCUBE polarimetry plots -- pure matplotlib.

Companion plotter for `red_ixpe.pcube`/`ixpe_pcube` (`skills/session.py`):
reads the provenance-headed `ixpe_pcube_<obsid(s)>.csv` that tool writes (one
row per energy band -- `band_id, E_LO, E_HI, PD, PD_err, PA, PA_err, QN,
QN_err, UN, UN_err, MDP99, SIGNIF`) and renders the PD(E)/PA(E) panels plus
the Q-U (QN vs UN) plane, in the style of Krawczynski+2022 Fig. S5.

Optionally overlays a linear PD(E)/PA(E) model on the PD/PA panels -- the
reviewer figure for the energy-dependence claim (Krawczynski+2022 Fig. S5
style): pass `pollin` (the `PD_A`/`PD_slope`/`PA_A`/`PA_slope` quartet from
`xspec_help.fit_tbabs_pollin_powerlaw`) to draw that analytic fit, or omit it
to fall back to an inverse-variance-weighted linear fit of the CSV's own
PD(E)/PA(E) points (`scipy.optimize.curve_fit`).

Imports neither `xspec` nor `ixpeobssim` -- like `pa_pd_stability.py` and
`bhxrb_pd_hr_survey.py`, this module reimplements its own small CSV reader
(mirroring `red_ixpe.read_pcube_csv`'s provenance-header contract exactly,
not a shared helper) so it -- and the `ixpe_pcube_plot` entry point -- runs
standalone, without the science stack, and is wired directly into
`TOOL_MANIFEST` rather than through `skills/session.py`.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import time as _time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

_PCUBE_HEADER_PREFIX = '# ixpe_pcube provenance: '
_PCUBE_BAND_FIELDS = ('E_LO', 'E_HI', 'PD', 'PD_err', 'PA', 'PA_err',
                       'QN', 'QN_err', 'UN', 'UN_err', 'MDP99', 'SIGNIF')
_PRODUCED_BY = 'xspec_utils.lib.ixpe_pcube_plot.ixpe_pcube_plot'
_POLLIN_KEYS = ('PD_A', 'PD_slope', 'PA_A', 'PA_slope')


def _pollin_params(pollin: dict) -> dict:
    """Pull the `PD_A`/`PD_slope`/`PA_A`/`PA_slope` quartet out of *pollin*
    (e.g. `fit_tbabs_pollin_powerlaw`'s result dict, or a plain quartet --
    any extra keys such as `Gamma`/`*_err` are ignored). Raises `ValueError`
    naming the missing key(s) if *pollin* is supplied but incomplete."""
    missing = [k for k in _POLLIN_KEYS if k not in pollin]
    if missing:
        raise ValueError(f'ixpe_pcube_plot: pollin block missing {missing}')
    return {k: float(pollin[k]) for k in _POLLIN_KEYS}


def _pollin_line(e: np.ndarray, pd_a: float, pd_slope: float, pa_a: float, pa_slope: float):
    """`pollin`'s own linear parametrisation: PD(E)=pd_a+pd_slope*(E-1 keV),
    PA(E)=pa_a+pa_slope*(E-1 keV) -- the same anchor-at-1-keV convention
    `fit_tbabs_pollin_powerlaw` fits."""
    return pd_a + pd_slope * (e - 1.0), pa_a + pa_slope * (e - 1.0)


def _weighted_linear_fit(x: np.ndarray, y: np.ndarray, yerr: np.ndarray):
    """Inverse-variance-weighted linear fit y=a+b*(x-1 keV) via
    `scipy.optimize.curve_fit` (sigma=yerr, absolute_sigma=True -- the
    standard weighted-least-squares normal equations). Returns
    `(a, b, a_err, b_err)`, or `None` when fewer than two finite points are
    available to constrain a line (e.g. a single-band CSV) -- rather than
    returning a degenerate slope.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    yerr = np.asarray(yerr, dtype=float)
    if len(x) < 2 or not np.all(yerr > 0):
        return None

    def _model(xx, a, b):
        return a + b * (xx - 1.0)

    weights = 1.0 / yerr ** 2
    p0 = (float(np.average(y, weights=weights)), 0.0)
    popt, pcov = curve_fit(_model, x, y, sigma=yerr, absolute_sigma=True, p0=p0)
    a_err, b_err = np.sqrt(np.diag(pcov))
    return float(popt[0]), float(popt[1]), float(a_err), float(b_err)


def _read_pcube_csv(csv_path: str) -> tuple[dict, list[dict]]:
    """Read an `ixpe_pcube_<obsid(s)>.csv` (written by `red_ixpe.write_pcube_csv`)
    into `(provenance, bands)`.

    Raises `ValueError` -- naming the offending path -- when `csv_path` does
    not exist, its first line is not the `ixpe_pcube` provenance header, or
    it carries zero band rows: an absent or empty upstream product is a hard
    input error here, never a silently-empty plot.
    """
    if not os.path.exists(csv_path):
        raise ValueError(f'ixpe_pcube_plot: csv_path not found: {csv_path}')
    with open(csv_path, newline='') as f:
        header_line = f.readline()
        if not header_line.startswith(_PCUBE_HEADER_PREFIX):
            raise ValueError(
                f'ixpe_pcube_plot: {csv_path} is missing the ixpe_pcube provenance header; '
                f'it is not a product of ixpe_pcube.')
        provenance = json.loads(header_line[len(_PCUBE_HEADER_PREFIX):])
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f'ixpe_pcube_plot: {csv_path} carries no band rows to plot.')
    bands = []
    for r in rows:
        band = {}
        for k in _PCUBE_BAND_FIELDS:
            v = r.get(k, '')
            band[k] = float(v) if v not in ('', None) else float('nan')
        bands.append(band)
    return provenance, bands


def _build_figure(csv_path: str, title: str | None = None, pollin: dict | None = None):
    """Read `csv_path`, mask/sort its bands, and build the 3-panel figure.

    Does the full array-preparation + `matplotlib` construction that
    `ixpe_pcube_plot` needs, but neither saves nor closes the figure -- the
    caller owns the figure's lifecycle (so a test can inspect the drawn
    artists directly, and `ixpe_pcube_plot` can `savefig`/`close` it once
    `save_path` is resolved). Raises `ValueError` exactly as `ixpe_pcube_plot`
    documents (see `_read_pcube_csv`, `_pollin_params`, and the all-non-finite
    check below).

    When *pollin* is given, its `PD_A`/`PD_slope`/`PA_A`/`PA_slope` quartet is
    drawn as the analytic PD(E)/PA(E) line on the PD/PA panels. Otherwise an
    inverse-variance-weighted linear fit of the plotted PD(E)/PA(E) points
    themselves (`_weighted_linear_fit`) is drawn instead, unless fewer than
    two bands survive the finite-value mask -- then no line is drawable and
    none is drawn (handled gracefully, not an error).

    Returns
    -------
    (fig, axes, plotted, label) -- `axes` is `{'pd': ax_pd, 'pa': ax_pa,
    'qu': ax_qu}`; `plotted` is the same dict `ixpe_pcube_plot` returns under
    its `'plotted'` key (`E_LO`, `E_HI`, `E_ctr`, `PD`, `PD_err`, `PA`,
    `PA_err`, `QN`, `QN_err`, `UN`, `UN_err`, `n_bands`, `n_skipped`,
    `overlay`); `label` is the obsid(s) label used for the figure
    title/caption.
    """
    provenance, bands = _read_pcube_csv(csv_path)

    e_lo = np.array([b['E_LO'] for b in bands])
    e_hi = np.array([b['E_HI'] for b in bands])
    pd = np.array([b['PD'] for b in bands])
    pd_err = np.array([b['PD_err'] for b in bands])
    pa = np.array([b['PA'] for b in bands])
    pa_err = np.array([b['PA_err'] for b in bands])
    qn = np.array([b['QN'] for b in bands])
    qn_err = np.array([b['QN_err'] for b in bands])
    un = np.array([b['UN'] for b in bands])
    un_err = np.array([b['UN_err'] for b in bands])

    mask = (np.isfinite(e_lo) & np.isfinite(e_hi)
            & np.isfinite(pd) & np.isfinite(pd_err)
            & np.isfinite(pa) & np.isfinite(pa_err)
            & np.isfinite(qn) & np.isfinite(qn_err)
            & np.isfinite(un) & np.isfinite(un_err))
    n_skipped = int(len(bands) - mask.sum())
    if not mask.any():
        raise ValueError(
            f'ixpe_pcube_plot: no band in {csv_path} has finite PD/PA/QN/UN values to plot.')

    order = np.argsort(e_lo[mask])
    e_lo, e_hi = e_lo[mask][order], e_hi[mask][order]
    pd, pd_err = pd[mask][order], pd_err[mask][order]
    pa, pa_err = pa[mask][order], pa_err[mask][order]
    qn, qn_err = qn[mask][order], qn_err[mask][order]
    un, un_err = un[mask][order], un_err[mask][order]
    e_ctr = 0.5 * (e_lo + e_hi)
    e_halfwidth = 0.5 * (e_hi - e_lo)

    obsids = provenance.get('obsids') or []
    label = '+'.join(obsids) if obsids else os.path.basename(csv_path)
    fig_title = title or f'IXPE energy-resolved polarimetry -- {label}'

    fig = plt.figure(figsize=(10, 8))
    gs = fig.add_gridspec(2, 2, width_ratios=(1.3, 1.0), hspace=0.08)
    ax_pd = fig.add_subplot(gs[0, 0])
    ax_pa = fig.add_subplot(gs[1, 0], sharex=ax_pd)
    ax_qu = fig.add_subplot(gs[:, 1])

    ax_pd.errorbar(e_ctr, pd, xerr=e_halfwidth, yerr=pd_err, fmt='o-', ms=5,
                   elinewidth=0.8, capsize=2, color='tab:orange')
    ax_pd.set_ylabel('PD')
    plt.setp(ax_pd.get_xticklabels(), visible=False)

    ax_pa.errorbar(e_ctr, pa, xerr=e_halfwidth, yerr=pa_err, fmt='o-', ms=5,
                   elinewidth=0.8, capsize=2, color='tab:blue')
    ax_pa.set_ylabel('PA (deg)')
    ax_pa.set_xlabel('Energy (keV)')

    e_line = np.linspace(e_lo.min(), e_hi.max(), 100)
    overlay = {'kind': None}
    if pollin is not None:
        params = _pollin_params(pollin)
        pd_line, pa_line = _pollin_line(
            e_line, params['PD_A'], params['PD_slope'], params['PA_A'], params['PA_slope'])
        ax_pd.plot(e_line, pd_line, '--', color='tab:red', lw=1.5, label='pollin fit')
        ax_pa.plot(e_line, pa_line, '--', color='tab:red', lw=1.5, label='pollin fit')
        overlay = {'kind': 'pollin', **params}
    else:
        pd_fit = _weighted_linear_fit(e_ctr, pd, pd_err)
        pa_fit = _weighted_linear_fit(e_ctr, pa, pa_err)
        if pd_fit is not None and pa_fit is not None:
            pd_a, pd_slope, pd_a_err, pd_slope_err = pd_fit
            pa_a, pa_slope, pa_a_err, pa_slope_err = pa_fit
            pd_line = pd_a + pd_slope * (e_line - 1.0)
            pa_line = pa_a + pa_slope * (e_line - 1.0)
            ax_pd.plot(e_line, pd_line, '--', color='tab:green', lw=1.5, label='weighted-linear fit')
            ax_pa.plot(e_line, pa_line, '--', color='tab:green', lw=1.5, label='weighted-linear fit')
            overlay = {
                'kind': 'fallback',
                'PD_A': pd_a, 'PD_slope': pd_slope, 'PD_A_err': pd_a_err, 'PD_slope_err': pd_slope_err,
                'PA_A': pa_a, 'PA_slope': pa_slope, 'PA_A_err': pa_a_err, 'PA_slope_err': pa_slope_err,
            }
    if overlay['kind'] is not None:
        ax_pd.legend(fontsize=8, loc='best')
        ax_pa.legend(fontsize=8, loc='best')

    ax_qu.errorbar(qn, un, xerr=qn_err, yerr=un_err, fmt='none', ecolor='gray',
                   elinewidth=0.6, alpha=0.6, zorder=1)
    sc = ax_qu.scatter(qn, un, c=e_ctr, cmap='viridis', zorder=2, s=30)
    fig.colorbar(sc, ax=ax_qu, label='Energy (keV)')
    ax_qu.axhline(0.0, color='k', lw=0.5)
    ax_qu.axvline(0.0, color='k', lw=0.5)
    ax_qu.set_xlabel('Q/I')
    ax_qu.set_ylabel('U/I')
    ax_qu.set_aspect('equal', adjustable='datalim')

    fig.suptitle(fig_title)

    plotted = {
        'E_LO': e_lo.tolist(), 'E_HI': e_hi.tolist(), 'E_ctr': e_ctr.tolist(),
        'PD': pd.tolist(), 'PD_err': pd_err.tolist(),
        'PA': pa.tolist(), 'PA_err': pa_err.tolist(),
        'QN': qn.tolist(), 'QN_err': qn_err.tolist(),
        'UN': un.tolist(), 'UN_err': un_err.tolist(),
        'n_bands': int(mask.sum()), 'n_skipped': n_skipped,
        'overlay': overlay,
    }
    axes = {'pd': ax_pd, 'pa': ax_pa, 'qu': ax_qu}
    return fig, axes, plotted, label


def ixpe_pcube_plot(csv_path: str, save_path: str | None = None, title: str | None = None,
                     pollin: dict | None = None) -> dict:
    """Render PD(E)/PA(E) and the Q-U plane from an `ixpe_pcube` CSV.

    Reads every band row of `csv_path` (written by `red_ixpe.write_pcube_csv`
    via `ixpe_pcube`) and writes exactly ONE PNG with three panels: PD vs
    band-centre energy (top left), PA vs band-centre energy (bottom left,
    sharing the energy axis), and QN vs UN (the normalized Stokes/Q-U plane,
    right, color-coded by band-centre energy) -- each with 1-sigma error
    bars, in the style of Krawczynski+2022 Fig. S5. `PD`/`PA`/`QN`/`UN` and
    their band-centre energies are plotted exactly as read from `csv_path`
    (no recomputation), and are also returned under `plotted` for
    differential verification against the source CSV.

    The PD/PA panels additionally overlay a linear PD(E)/PA(E) model:
    - `pollin` given: draws PD(E)=PD_A+PD_slope*(E-1 keV) and
      PA(E)=PA_A+PA_slope*(E-1 keV) from the supplied quartet (e.g.
      `fit_tbabs_pollin_powerlaw`'s result dict, passed straight through --
      extra keys are ignored) -- the model actually fit to the spectra.
    - `pollin` omitted: falls back to an inverse-variance-weighted linear fit
      of the CSV's own PD(E)/PA(E) points (`scipy.optimize.curve_fit`), drawn
      the same way. A CSV with fewer than two usable bands has no line to
      fit -- handled gracefully, the points alone are still plotted.

    A band with any non-finite `E_LO`/`E_HI`/`PD`/`PD_err`/`PA`/`PA_err`/
    `QN`/`QN_err`/`UN`/`UN_err` is skipped -- never plotted as zero. A single
    remaining band still renders (a lone marker; there is no second point to
    draw a connecting line to). Raises `ValueError` (see `_read_pcube_csv`,
    `_pollin_params`, and `_build_figure`) if `csv_path` is missing,
    unrecognised, or carries zero band rows, again if every band is skipped
    for non-finite values, or if `pollin` is supplied but missing any of
    `PD_A`/`PD_slope`/`PA_A`/`PA_slope` -- an empty or mislabeled plot is
    never silently produced.

    Parameters
    ----------
    csv_path : str
        Path to the `ixpe_pcube_<obsid(s)>.csv` written by `ixpe_pcube`.
    save_path : str or None
        Full output PNG path. Auto-generated timestamped file in
        `os.getcwd()` when None.
    title : str or None
        Figure suptitle. Defaults to the provenance's `obsids` joined by '+'.
    pollin : dict or None
        `PD_A`/`PD_slope`/`PA_A`/`PA_slope` quartet (fraction, fraction/keV,
        degrees, degrees/keV) to overlay as the analytic PD(E)/PA(E) model,
        e.g. straight from `fit_tbabs_pollin_powerlaw`'s result dict. Omit to
        overlay a `scipy` weighted-linear fit of the CSV's own PD(E)/PA(E)
        points instead.

    Returns
    -------
    dict with 'figures' (one `FigureRef` dict), 'plotted' (the exact E_LO,
    E_HI, E_ctr, PD, PD_err, PA, PA_err, QN, QN_err, UN, UN_err arrays drawn,
    plus 'n_bands'/'n_skipped'/'overlay' -- 'overlay' is `{'kind': None}` when
    no line was drawable, else `{'kind': 'pollin'|'fallback', 'PD_A',
    'PD_slope', 'PA_A', 'PA_slope', ...}`), 'stdout', 'exit_code' (always 0
    -- failures raise instead, see above).
    """
    fig, _axes, plotted, label = _build_figure(csv_path, title=title, pollin=pollin)

    if save_path is None:
        ts = _time.strftime('%Y%m%dT%H%M%S')
        save_path = os.path.join(os.getcwd(), f'ixpe_pcube_plot_{ts}.png')
    os.makedirs(os.path.dirname(os.path.abspath(save_path)) or '.', exist_ok=True)

    savefig_kwargs = {'metadata': {'Software': _PRODUCED_BY}}
    fig.savefig(save_path, dpi=150, bbox_inches='tight', **savefig_kwargs)
    plt.close(fig)
    sha = hashlib.sha256(Path(save_path).read_bytes()).hexdigest()

    return {
        'figures': [{
            'path': save_path,
            'caption': f'PD(E)/PA(E) and Q-U plane -- {label}',
            'mime_type': 'image/png',
            'produced_by': _PRODUCED_BY,
            'created_at': _time.strftime('%Y-%m-%dT%H:%M:%SZ', _time.gmtime()),
            'sha256': sha,
        }],
        'plotted': plotted,
        'stdout': (f'ixpe_pcube_plot: {plotted["n_bands"]} band(s) plotted '
                   f'({plotted["n_skipped"]} skipped), overlay={plotted["overlay"]["kind"]} '
                   f'-> {save_path}'),
        'exit_code': 0,
    }
