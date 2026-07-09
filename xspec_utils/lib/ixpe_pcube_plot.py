"""Energy-resolved IXPE PCUBE polarimetry plots -- pure matplotlib.

Companion plotter for `red_ixpe.pcube`/`ixpe_pcube` (`skills/session.py`):
reads the provenance-headed `ixpe_pcube_<obsid(s)>.csv` that tool writes (one
row per energy band -- `band_id, E_LO, E_HI, PD, PD_err, PA, PA_err, QN,
QN_err, UN, UN_err, MDP99, SIGNIF`) and renders the PD(E)/PA(E) panels plus
the Q-U (QN vs UN) plane, in the style of Krawczynski+2022 Fig. S5.

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

_PCUBE_HEADER_PREFIX = '# ixpe_pcube provenance: '
_PCUBE_BAND_FIELDS = ('E_LO', 'E_HI', 'PD', 'PD_err', 'PA', 'PA_err',
                       'QN', 'QN_err', 'UN', 'UN_err', 'MDP99', 'SIGNIF')
_PRODUCED_BY = 'xspec_utils.lib.ixpe_pcube_plot.ixpe_pcube_plot'


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


def _build_figure(csv_path: str, title: str | None = None):
    """Read `csv_path`, mask/sort its bands, and build the 3-panel figure.

    Does the full array-preparation + `matplotlib` construction that
    `ixpe_pcube_plot` needs, but neither saves nor closes the figure -- the
    caller owns the figure's lifecycle (so a test can inspect the drawn
    artists directly, and `ixpe_pcube_plot` can `savefig`/`close` it once
    `save_path` is resolved). Raises `ValueError` exactly as `ixpe_pcube_plot`
    documents (see `_read_pcube_csv` and the all-non-finite check below).

    Returns
    -------
    (fig, axes, plotted, label) -- `axes` is `{'pd': ax_pd, 'pa': ax_pa,
    'qu': ax_qu}`; `plotted` is the same dict `ixpe_pcube_plot` returns under
    its `'plotted'` key (`E_LO`, `E_HI`, `E_ctr`, `PD`, `PD_err`, `PA`,
    `PA_err`, `QN`, `QN_err`, `UN`, `UN_err`, `n_bands`, `n_skipped`); `label`
    is the obsid(s) label used for the figure title/caption.
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
    }
    axes = {'pd': ax_pd, 'pa': ax_pa, 'qu': ax_qu}
    return fig, axes, plotted, label


def ixpe_pcube_plot(csv_path: str, save_path: str | None = None, title: str | None = None) -> dict:
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

    A band with any non-finite `E_LO`/`E_HI`/`PD`/`PD_err`/`PA`/`PA_err`/
    `QN`/`QN_err`/`UN`/`UN_err` is skipped -- never plotted as zero. A single
    remaining band still renders (a lone marker; there is no second point to
    draw a connecting line to). Raises `ValueError` (see `_read_pcube_csv`
    and `_build_figure`) if `csv_path` is missing, unrecognised, or carries
    zero band rows, and again if every band is skipped for non-finite
    values -- an empty plot is never silently produced.

    Parameters
    ----------
    csv_path : str
        Path to the `ixpe_pcube_<obsid(s)>.csv` written by `ixpe_pcube`.
    save_path : str or None
        Full output PNG path. Auto-generated timestamped file in
        `os.getcwd()` when None.
    title : str or None
        Figure suptitle. Defaults to the provenance's `obsids` joined by '+'.

    Returns
    -------
    dict with 'figures' (one `FigureRef` dict), 'plotted' (the exact E_LO,
    E_HI, E_ctr, PD, PD_err, PA, PA_err, QN, QN_err, UN, UN_err arrays drawn,
    plus 'n_bands'/'n_skipped'), 'stdout', 'exit_code' (always 0 -- failures
    raise instead, see above).
    """
    fig, _axes, plotted, label = _build_figure(csv_path, title=title)

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
                   f'({plotted["n_skipped"]} skipped) -> {save_path}'),
        'exit_code': 0,
    }
