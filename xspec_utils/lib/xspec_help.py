import contextlib
import io
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from xspec import AllData, AllModels, Fit, Model, Plot, Xset

def load_data(
    epoch: dict,
    ixpe_path: str,
    nicer_path: str | None = None,
    nustar_path: str | None = None,
    ixpe_changes: bool = False,
    ixpe_src_rad: float = 1.0,
    ixpe_bkg_inner_rad: float = 3.0,
    ixpe_bkg_outer_rad: float = 5.0,
    ixpe_irf_date: str = '20230702',
    ixpe_erange: tuple[float, float] = (2.0, 8.0),
    nicer_erange: tuple[float, float] = (1.0, 10.0),
    nustar_erange: tuple[float, float] = (4.0, 50.0),
    nustar_use_bkg: bool = False,
    clear_existing: bool = True,
    chatter: bool = False,
) -> dict:
    """Load spectra for one epoch into XSPEC.

    Instruments are loaded in *epoch* key order, so the key order sets the
    data-group numbering.  **Load IXPE first** so IXPE DU1 becomes data group 1:
    the polarized fit models (``plabs`` / ``polconst``) anchor their default
    cross-group links to DG1, which must be a polarization-bearing (IXPE I/Q/U)
    group — anchoring DG1 to a Stokes-I-only spectrum (NuSTAR/NICER) makes the
    fit diverge.  Recommended order: IXPE → NICER → NuSTAR, i.e. IXPE DU1/DU2/DU3
    → groups 1/2/3, NICER → next groups, NuSTAR FPMA/FPMB → last groups.
    Instruments absent from *epoch* are silently skipped (e.g. no IXPE → the
    first present instrument becomes DG1).

    Note: this DG/load order is independent of ``fit_*``'s ``spec_pref``, which
    selects the cross-calibration *reference* group (default NuSTAR) whose
    ``constant.factor`` is frozen to 1 — a different concept from which
    instrument is DG1.

    For IXPE each detector unit contributes three consecutive spectra
    (Stokes I, Q, U).  The Q and U gain slopes are linked to the
    corresponding DU's I gain slope to enforce a shared calibration.

    Parameters
    ----------
    epoch : dict
        Single-epoch sub-dict from ``epoch.epoch``, e.g. ``epoch.epoch['E5']``.
        Expected keys: ``'IXPE'``, ``'NICER'``, ``'NuSTAR'`` (each a list of
        obs-ID strings) and ``'STATE'``.  Only one IXPE observation per epoch
        is supported.
    ixpe_path : str
        Path to the directory containing IXPE observation sub-directories
        (passed to ``red_ixpe.get_all_file_list``).
    nicer_path : str or None
        Path to the NICER data directory.  Required when ``'NICER'`` is in
        *epoch*.
    nustar_path : str or None
        Path to the NuSTAR data directory.  Required when ``'NuSTAR'`` is in
        *epoch*.
    ixpe_changes : bool
        If True, re-run the IXPE reduction pipeline; if False (default), use
        existing pipeline products on disk.
    ixpe_src_rad : float
        IXPE source-region radius in arcminutes (default 1.0).
    ixpe_bkg_inner_rad : float
        IXPE background annulus inner radius in arcminutes (default 3.0).
    ixpe_bkg_outer_rad : float
        IXPE background annulus outer radius in arcminutes (default 5.0).
    ixpe_irf_date : str
        IRF date string in ``'YYYYMMDD'`` format used to select the IXPE
        response (default ``'20230702'``).
    ixpe_erange : tuple[float, float]
        (lo, hi) energy range in keV to *notice* for IXPE (default 2–8 keV).
    nicer_erange : tuple[float, float]
        (lo, hi) energy range in keV to *notice* for NICER (default 1–10 keV).
    nustar_erange : tuple[float, float]
        (lo, hi) energy range in keV to *notice* for NuSTAR (default 4–50 keV).
    nustar_use_bkg : bool
        If True, load the NuSTAR background spectrum.  Disabled by default
        because clean background regions are often unavailable for Cyg X-1.
    clear_existing : bool
        If True (default), call ``AllData.clear()`` and ``AllModels.clear()``
        before loading new spectra.
    chatter : bool
        If True, print file paths and per-instrument loading summaries.

    Returns
    -------
    dict
        Uniform nested dict of XSPEC spectrum indices keyed by instrument →
        observation ID → detector.  Only instruments present in *epoch*
        appear.  Each NuSTAR FPM, each NICER obs, and each IXPE DU occupies
        its own data group; DGs are assigned consecutively from 1 with no
        gaps (recover a DG via ``AllData(spec_idx).dataGroup``)::

            {
                'IXPE': {
                    '<obsid>': {
                        'DU1': [1,  2,  3],
                        'DU2': [4,  5,  6],
                        'DU3': [7,  8,  9],
                    },
                },
                'NICER': {
                    '<obsid>': {'mpu7': [10]},
                },
                'NuSTAR': {
                    '<obsid>': {'FPMA': [11], 'FPMB': [12]},
                },
            }

        IXPE per-DU lists are ordered ``[I, Q, U]``.  NuSTAR FPM labels
        follow alphabetical (sorted) filename order.
    """
    import red_ixpe as _red_ixpe
    import NICER as _NICER
    import NuSTAR as _NuSTAR

    if clear_existing:
        AllData.clear()
        AllModels.clear()

    result = {}
    spec_index = 1  # running 1-based XSPEC spectrum counter
    grp = 1         # running data-group counter; incremented after each instrument block

    for inst in epoch:

        # ------------------------------------------------------------ NuSTAR
        if inst == 'NuSTAR':
            if nustar_path is None:
                raise ValueError(
                    "'nustar_path' is required when the epoch contains NuSTAR observations.")

            FPM_LABELS = ['FPMA', 'FPMB']
            result['NuSTAR'] = {}
            for obs in epoch['NuSTAR']:
                src_lst, bkg_lst = _NuSTAR.pha(obs, nustar_path, chatter=chatter)
                rmf_lst = _NuSTAR.rmf(obs, nustar_path, chatter=chatter)
                arf_lst = _NuSTAR.arf(obs, nustar_path, chatter=chatter)

                for f in src_lst:
                    with fits.open(f, mode='update') as hdul:
                        hdul['SPECTRUM'].header['XFLT0001'] = 'Stokes:0'

                obs_result = {}
                for fpm, f_src, f_bkg, f_rmf, f_arf in zip(
                        FPM_LABELS, src_lst, bkg_lst, rmf_lst, arf_lst):
                    AllData(f'{grp}:{spec_index} {f_src}')
                    if nustar_use_bkg:
                        AllData(spec_index).background = f_bkg
                    AllData(spec_index).response = f_rmf
                    AllData(spec_index).response.arf = f_arf
                    AllData(spec_index).ignore(
                        f'**-{nustar_erange[0]} {nustar_erange[1]}-**')
                    obs_result[fpm] = [spec_index]
                    spec_index += 1
                    grp += 1  # one DG per NuSTAR FPM

                result['NuSTAR'][obs] = obs_result
                if chatter:
                    print(f'Loaded NuSTAR {obs}: {obs_result}')

        # ------------------------------------------------------------- NICER
        elif inst == 'NICER':
            if nicer_path is None:
                raise ValueError(
                    "'nicer_path' is required when the epoch contains NICER observations.")

            result['NICER'] = {}
            for obs in epoch['NICER']:
                src, bkg = _NICER.pha(obs, nicer_path, chatter=chatter)
                f_rmf = _NICER.rmf(obs, nicer_path, chatter=chatter)
                f_arf = _NICER.arf(obs, nicer_path, chatter=chatter)

                for f in (src, bkg):
                    with fits.open(f, mode='update') as hdul:
                        hdul['SPECTRUM'].header['XFLT0001'] = 'Stokes:0'

                AllData(f'{grp}:{spec_index} {src}')
                AllData(spec_index).background = bkg
                AllData(spec_index).response = f_rmf
                AllData(spec_index).response.arf = f_arf
                AllData(spec_index).ignore(
                    f'**-{nicer_erange[0]} {nicer_erange[1]}-**')
                result['NICER'][obs] = {'mpu7': [spec_index]}
                if chatter:
                    print(f'Loaded NICER {obs}: group {grp}, spec {spec_index}')
                spec_index += 1
                grp += 1  # one DG per NICER obs

        # -------------------------------------------------------------- IXPE
        elif inst == 'IXPE':
            if len(epoch['IXPE']) > 1:
                raise NotImplementedError(
                    "load_data currently supports only one IXPE observation per epoch; "
                    f"got {epoch['IXPE']}.")

            ixpe_file_list = _red_ixpe.get_all_file_list(ixpe_path, chatter=chatter)
            obs = epoch['IXPE'][0]
            ixpe_src, ixpe_bkg = _red_ixpe.pha(
                obs=obs,
                file_list=ixpe_file_list,
                src_rad=ixpe_src_rad,
                bkg_inner_rad=ixpe_bkg_inner_rad,
                bkg_outer_rad=ixpe_bkg_outer_rad,
                irf_date=ixpe_irf_date,
                changes=ixpe_changes,
                chatter=chatter,
            )

            DU_LABELS = ['DU1', 'DU2', 'DU3']
            ixpe_obs_result = {}

            for du_idx, du_label in enumerate(DU_LABELS):
                f_src_i = ixpe_src[du_idx]
                f_bkg_i = ixpe_bkg[du_idx]
                du_specs = []

                for sto in ('I', 'Q', 'U'):
                    if sto == 'I':
                        f_src_sto, f_bkg_sto = f_src_i, f_bkg_i
                    elif sto == 'Q':
                        f_src_sto = _red_ixpe.get_pha_q(f_src_i)
                        f_bkg_sto = _red_ixpe.get_pha_q(f_bkg_i)
                    else:  # U
                        f_src_sto = _red_ixpe.get_pha_u(f_src_i)
                        f_bkg_sto = _red_ixpe.get_pha_u(f_bkg_i)

                    AllData(f'{grp}:{spec_index} {f_src_sto}')
                    AllData(spec_index).background = f_bkg_sto
                    AllData(spec_index).ignore(
                        f'**-{ixpe_erange[0]} {ixpe_erange[1]}-**')

                    r = AllData(spec_index).response
                    r.gain.slope = 1.0
                    r.gain.offset = 0.0
                    r.gain.offset.frozen = True
                    if sto in ('Q', 'U'):
                        r.gain.slope.link = AllData(du_specs[0]).response.gain.slope

                    du_specs.append(spec_index)
                    spec_index += 1

                ixpe_obs_result[du_label] = du_specs  # [I, Q, U]
                if chatter:
                    print(f'Loaded IXPE {du_label}: group {grp}, '
                          f'spectra {du_specs} (I, Q, U)')
                grp += 1  # one DG per DU

            result['IXPE'] = {obs: ixpe_obs_result}

    return result


def freeze_all(thaw=False):
    """Freeze (or thaw) every parameter across all data groups.

    Parameters
    ----------
    thaw : bool
        If True, thaw instead of freeze.
    """
    for dg in range(1, AllData.nGroups + 1):
        m = AllModels(dg)
        for i in range(1, m.nParameters + 1):
            p = m(i)
            if not p.link:
                p.frozen = not thaw


def robust_fit(verbose=True, max_iter_simplex=10, max_iter_leven=100):
    """Replace Fit.perform() with a Simplex → Levenberg-Marquardt strategy.

    Simplex provides derivative-free global exploration; Levenberg-Marquardt
    polishes the result.  The fit method is restored to Levenberg-Marquardt
    after the call.

    Parameters
    ----------
    verbose : bool
        Print per-step and final fit statistics.
    max_iter_simplex : int
        Maximum iterations for the Simplex fit (default 10).
    max_iter_leven : int
        Maximum iterations for the Levenberg-Marquardt fit (default 100).

    Returns
    -------
    stat : float
        Final fit statistic.
    dof : int
        Degrees of freedom.
    """
    methods = [
        ("simplex",       "Simplex",              max_iter_simplex),
        ("leven 10 0.01", "Levenberg-Marquardt",  max_iter_leven),
    ]

    for method_str, label, max_iter in methods:
        Fit.nIterations = max_iter
        Fit.method = method_str
        try:
            Fit.perform()
        except Exception as e:
            if verbose:
                print(f"  [{label}] {e}")
        if verbose:
            print(f"  [{label}]  stat = {Fit.statistic:.4f} / {Fit.dof} dof")

    if verbose:
        print(f"  Final: stat/dof = {Fit.statistic / Fit.dof:.4f}")

    return Fit.statistic, Fit.dof


def inspect_fit(
    spec_groups,
    plot_type='ufspec',
    residuals=('delchi',),
    title='',
    xlim=None,
    xscale='log',
    yscale='log',
    figsize=None,
    legend=True,
    show_model=True,
    show_components=True,
):
    """Plot XSPEC fit results for grouped spectra.

    Parameters
    ----------
    spec_groups : dict
        Mapping of label prefix to list of XSPEC spectrum indices,
        e.g. {'IXPE Det': [5, 6, 7], 'NICER': [3, 4], 'NuSTAR': [1, 2]}.
    plot_type : str
        XSPEC plot type for the main panel: 'ufspec' or 'data'.
    residuals : sequence of str
        Ordered list of residual panel types.  Each entry must be one of
        'delchi' or 'ratio'.  One panel is created per entry.
    title : str
        Figure suptitle.
    xlim : tuple or None
        (xmin, xmax) for the shared x-axis.
    xscale : str
        'log' or 'linear' for the x-axis.
    yscale : str
        'log' or 'linear' for the main spectrum y-axis.
    figsize : tuple or None
        Figure size.  Defaults to (10, 4 + 2*len(residuals)).
    legend : bool
        Whether to call fig.legend() with an outside-right legend.
    show_model : bool
        Whether to overplot the total folded model.
    show_components : bool
        Whether to overplot individual additive model components as dashed
        lines.  Only one copy per component is drawn (from the first
        spectrum of the first group).  Requires Plot.add = True, which is
        set temporarily inside this function.

    Returns
    -------
    fig, axes : matplotlib Figure and list of Axes (spec first, then residuals).
    """
    n_resid = len(residuals)
    if figsize is None:
        figsize = (10, 4 + 2 * n_resid)

    height_ratios = [2] + [1] * n_resid
    fig, axes = plt.subplots(
        nrows=1 + n_resid,
        ncols=1,
        figsize=figsize,
        sharex=True,
        gridspec_kw={'height_ratios': height_ratios},
    )
    if n_resid == 0:
        axes = [axes]

    ax_spec = axes[0]

    # ---- main spectrum panel ----
    # Enable add-component extraction when requested
    _prev_add = Plot.add
    if show_model and show_components:
        Plot.add = True
    Plot(plot_type)

    for label, indices in spec_groups.items():
        for idx in indices:
            E      = Plot.x(idx)
            xerr   = Plot.xErr(idx)
            y      = Plot.y(idx)
            yerr   = Plot.yErr(idx)
            ax_spec.errorbar(
                E, y, xerr=xerr, yerr=yerr,
                fmt='o', ms=2, mew=0.5, elinewidth=0.5,
                capsize=1.5, capthick=0.5, alpha=0.7,
                label=f'{label} {idx}',
            )
            if show_model:
                model_y = Plot.model(idx)
                ax_spec.plot(E, model_y, lw=1, zorder=99, color='r')

    # ---- additive components (one copy each, from first spectrum) ----
    if show_model and show_components:
        _ref_idx = next(iter(next(iter(spec_groups.values()))))
        n_comp = Plot.nAddComps(_ref_idx)
        if n_comp > 0:
            # Retrieve component names from model
            try:
                from xspec import AllData as _AD, AllModels as _AM
                _dg = _AD(_ref_idx).dataGroup
                _m  = _AM(_dg)
                _add_names = [cn for cn in _m.componentNames
                              if getattr(_m, cn).type == 'add']
            except Exception:
                _add_names = []
            _E_ref = np.array(Plot.x(_ref_idx))
            _comp_ls     = ['--', ':', '-.', (0, (3, 1, 1, 1)), (0, (5, 1))]
            _comp_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
            for ci in range(1, n_comp + 1):
                cname  = _add_names[ci - 1] if ci <= len(_add_names) else f'comp{ci}'
                comp_y = np.array(Plot.addComp(ci, _ref_idx))
                ax_spec.plot(
                    _E_ref, comp_y,
                    ls=_comp_ls[(ci - 1) % len(_comp_ls)],
                    color=_comp_colors[ci % len(_comp_colors)],
                    lw=1, alpha=0.8,
                    label=cname,
                    zorder=97,
                )
        Plot.add = _prev_add

    ax_spec.set_ylabel('Counts s$^{-1}$ keV$^{-1}$')
    ax_spec.set_xscale(xscale)
    ax_spec.set_yscale(yscale)
    if xlim is not None:
        ax_spec.set_xlim(xlim)

    # ---- residual panels ----
    for ax_r, resid_type in zip(axes[1:], residuals):
        Plot(resid_type)
        for label, indices in spec_groups.items():
            for idx in indices:
                E = Plot.x(idx)
                y = Plot.y(idx)
                if resid_type == 'ratio':
                    yerr   = Plot.yErr(idx)
                    xerr   = Plot.xErr(idx)
                    ax_r.errorbar(
                        E, y, xerr=xerr, yerr=yerr,
                        fmt='o', ms=2, mew=0.5, elinewidth=0.5,
                        capsize=1.5, capthick=0.5, alpha=0.7,
                        label=f'{label} {idx}',
                    )
                    ax_r.axhline(1, lw=0.8, color='k')
                    ax_r.set_ylabel('ratio')
                else:  # delchi
                    ax_r.scatter(E, y, s=5, label=f'{label} {idx}')
                    ax_r.axhline(0, lw=0.8, color='k')
                    ax_r.set_ylabel('$\\Delta\\chi$')

        ax_r.set_xlabel('Energy (keV)')

    if title:
        fig.suptitle(title)

    if legend:
        fig.legend(loc='upper right', fontsize='small',
                   bbox_to_anchor=(1.0, 0.95))
        fig.tight_layout(rect=[0, 0, 0.82, 0.95])
    else:
        fig.tight_layout()

    return fig, axes


def group_pha(infile, outfile, e_edges):
    """Write a grouped copy of a PHA file with custom energy bin edges.

    Reads the RMF path from the PHA RESPFILE header, maps e_edges to
    channel indices via the RMF EBOUNDS extension, and writes the result
    to outfile.  Load outfile into XSPEC instead of infile.

    Parameters
    ----------
    infile : str
        Input PHA file path.
    outfile : str
        Output grouped PHA file path.
    e_edges : sequence of float
        Bin edges in keV, e.g. [2., 3., 4., 5., 6., 8.].

    Example
    -------
    >>> e_bins = [2., 3., 4., 5., 6., 8.]
    >>> q_src = red_ixpe.get_pha_q(f_src)
    >>> q_grp = q_src.replace('.fits', '_grp.fits')
    >>> xspec_help.group_pha(q_src, q_grp, e_bins)
    >>> AllData(f'{grp}:{idx} {q_grp}')
    """
    with fits.open(infile) as hdul:
        rmf = hdul['SPECTRUM'].header['RESPFILE']
        if not os.path.isabs(rmf):
            rmf = os.path.join(os.path.dirname(os.path.abspath(infile)), rmf)
        with fits.open(rmf) as rmf_hdul:
            e_lo = rmf_hdul['EBOUNDS'].data['E_MIN'].astype(float)
            e_hi = rmf_hdul['EBOUNDS'].data['E_MAX'].astype(float)
        e_mid = 0.5 * (e_lo + e_hi)
        group = np.zeros(len(e_mid), dtype=np.int16)
        for lo, hi in zip(e_edges[:-1], e_edges[1:]):
            idx = np.where((e_mid >= lo) & (e_mid < hi))[0]
            if len(idx):
                group[idx[0]]  = 1
                group[idx[1:]] = -1
        spec = hdul['SPECTRUM']
        if 'GROUPING' in spec.columns.names:
            spec.data['GROUPING'] = group
        else:
            new_col = fits.Column(name='GROUPING', format='I', array=group)
            new_hdu = fits.BinTableHDU.from_columns(
                spec.columns + fits.ColDefs([new_col]), header=spec.header)
            hdul['SPECTRUM'] = new_hdu
        hdul.writeto(outfile, overwrite=True)


def rebin_plot(spec_idx, e_edges, plot_type='ufspec'):
    """Extract XSPEC plot data and rebin to custom energy edges.

    For visualisation only — does not affect the fit.  Useful when the
    fit is run on fine channels but you want to display coarse-band points.

    Parameters
    ----------
    spec_idx : int
        XSPEC spectrum index (1-based).
    e_edges : sequence of float
        Bin edges in keV.
    plot_type : str
        XSPEC plot type to call before extraction (default 'ufspec').

    Returns
    -------
    e_cen, y, yerr : np.ndarray
        Band-centre energies, inverse-variance-weighted mean y values,
        and propagated 1-sigma errors.
    """
    Plot(plot_type)
    e    = np.array(Plot.x(spec_idx))
    y    = np.array(Plot.y(spec_idx))
    yerr = np.array(Plot.yErr(spec_idx))
    de   = np.array(Plot.xErr(spec_idx)) * 2   # full channel width

    e_lo_arr = e - de / 2
    e_hi_arr = e + de / 2

    e_cen_out, y_out, yerr_out = [], [], []
    for lo, hi in zip(e_edges[:-1], e_edges[1:]):
        mask = (e_lo_arr >= lo - 1e-6) & (e_hi_arr <= hi + 1e-6)
        if not mask.any():
            continue
        w = 1.0 / yerr[mask] ** 2
        e_cen_out.append(0.5 * (lo + hi))
        y_out.append(np.average(y[mask], weights=w))
        yerr_out.append(1.0 / np.sqrt(w.sum()))
    return np.array(e_cen_out), np.array(y_out), np.array(yerr_out)


def plot_iqu_decomp(
    iqu_groups,
    e_bins,
    show_components=True,
    title_prefix='',
    figsize=(14, 4),
):
    """Plot per-DU Stokes I / Q / U panels with all additive component overlays.

    Creates one figure per DU.  All three Stokes spectra are rebinned to
    ``e_bins`` using inverse-variance weighting.  When ``show_components=True``,
    every additive model component is extracted via ``Plot.addComp`` and drawn
    as a step line (one colour/style per component, consistent across panels).

    Parameters
    ----------
    iqu_groups : dict
        ``{label: [I_idx, Q_idx, U_idx]}`` — 1-based XSPEC spectrum indices
        for each DU, e.g. ``{'DU1': [1, 4, 7], 'DU2': [2, 5, 8], ...}``.
    e_bins : sequence of float
        Energy bin edges in keV, e.g. ``[2., 3., 4., 5., 6., 7., 8.]``.
    show_components : bool
        Overlay all additive model components as step lines (default True).
    title_prefix : str
        Prepended to each figure suptitle, e.g. the observation ID.
    figsize : tuple
        Size of each per-DU figure.

    Returns
    -------
    figs : list of ``(fig, (ax_i, ax_q, ax_u))`` tuples, one per DU label.

    Example
    -------
    >>> figs = xspec_help.plot_iqu_decomp(
    ...     iqu_groups={'DU1': [1, 4, 7], 'DU2': [2, 5, 8], 'DU3': [3, 6, 9]},
    ...     e_bins=np.array([2., 3., 4., 5., 6., 7., 8.]),
    ...     title_prefix='02008401',
    ... )
    """
    Plot.xAxis = 'keV'
    _prev_add = Plot.add
    if show_components:
        Plot.add = True
    Plot('data')

    e_bins  = np.asarray(e_bins)
    e_cen   = 0.5 * (e_bins[:-1] + e_bins[1:])
    e_hw    = 0.5 * (e_bins[1:]  - e_bins[:-1])
    n_bands = len(e_bins) - 1

    # component names from the model of the first I spectrum
    _first_i = int(next(iter(iqu_groups.values()))[0])
    _dg      = AllData(_first_i).dataGroup
    _m       = AllModels(_dg)
    try:
        _add_names = [cn for cn in _m.componentNames
                      if getattr(_m, cn).type == 'add']
    except Exception:
        _add_names = []
    n_comp = Plot.nAddComps(_first_i) if show_components else 0

    _comp_ls     = ['--', ':', '-.', (0, (3, 1, 1, 1)), (0, (5, 1))]
    _comp_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    def _rebin(spec_idx):
        """Rebin Plot.y for spec_idx to e_bins; return e_cen, y_rb, yerr_rb."""
        E    = np.array(Plot.x(spec_idx))
        dE   = np.array(Plot.xErr(spec_idx)) * 2
        y    = np.array(Plot.y(spec_idx))
        yerr = np.array(Plot.yErr(spec_idx))
        e_lo = E - dE / 2
        e_hi = E + dE / 2
        y_rb, yerr_rb, x_rb, xerr_rb = [], [], [], []
        for j, (lo, hi) in enumerate(zip(e_bins[:-1], e_bins[1:])):
            mask = (e_lo >= lo - 1e-6) & (e_hi <= hi + 1e-6)
            if mask.any():
                w = 1.0 / yerr[mask] ** 2
                y_rb.append(np.average(y[mask], weights=w))
                yerr_rb.append(1.0 / np.sqrt(w.sum()))
                x_rb.append(e_cen[j])
                xerr_rb.append(e_hw[j])
        return (np.array(x_rb), np.array(xerr_rb),
                np.array(y_rb), np.array(yerr_rb))

    def _rebin_comp(ci, spec_idx):
        """Rebin additive component ci for spec_idx to e_bins."""
        E    = np.array(Plot.x(spec_idx))
        dE   = np.array(Plot.xErr(spec_idx)) * 2
        comp = np.array(Plot.addComp(ci, spec_idx))
        e_lo = E - dE / 2
        e_hi = E + dE / 2
        rb = np.full(n_bands, np.nan)
        for j, (lo, hi) in enumerate(zip(e_bins[:-1], e_bins[1:])):
            mask = (e_lo >= lo - 1e-6) & (e_hi <= hi + 1e-6)
            if mask.any():
                rb[j] = np.average(comp[mask], weights=dE[mask])
        return rb

    figs = []
    for label, (i_idx, q_idx, u_idx) in iqu_groups.items():
        i_idx, q_idx, u_idx = int(i_idx), int(q_idx), int(u_idx)

        fig, (ax_i, ax_q, ax_u) = plt.subplots(1, 3, figsize=figsize,
                                                sharex=True)
        fig.suptitle(f'{title_prefix} {label}'.strip())

        for ax, sidx, stitle in [
            (ax_i, i_idx, 'Stokes I'),
            (ax_q, q_idx, 'Stokes Q'),
            (ax_u, u_idx, 'Stokes U'),
        ]:
            xc, xe, yc, ye = _rebin(sidx)
            ax.errorbar(xc, yc, xerr=xe, yerr=ye,
                        fmt='o', ms=4, color='k', capsize=2)
            ax.set_title(stitle)
            ax.set_xlabel('Energy (keV)')

            if show_components:
                for ci in range(1, n_comp + 1):
                    cname = (_add_names[ci - 1] if ci <= len(_add_names)
                             else f'comp{ci}')
                    rb = _rebin_comp(ci, sidx)
                    ax.step(e_bins[:-1], rb, where='post',
                            ls=_comp_ls[(ci - 1) % len(_comp_ls)],
                            color=_comp_colors[ci % len(_comp_colors)],
                            lw=1.5, label=cname)
            ax.legend(fontsize='small')

        ax_i.set_ylabel('counts s$^{-1}$ keV$^{-1}$')
        fig.tight_layout()
        figs.append((fig, (ax_i, ax_q, ax_u)))

    Plot.add = _prev_add
    return figs


def fit_tbabs_bknpower_polconst_diskbb_polconst_relxillCp(
    indices: dict,
    spec_pref: list[str] = ['NuSTAR', 'NICER', 'IXPE'],
    audit_data: bool = True,
    log_path: str | None = None,
    chatter_level: int = 10,
) -> None:
    """Define and fit the model ``constant * tbabs * bknpower * (polconst*diskbb + polconst*relxillCp)``.

    Sets up the model, assigns initial parameter values, performs the fit, and
    optionally logs the verbose model/fit summary to a file.  The
    cross-calibration constant of the preferred reference instrument is
    frozen to 1; all other instrument constants are freed.

    Parameters
    ----------
    indices : dict
        Return value of :func:`load_data`.  Used to look up data-group
        numbers for each instrument.
    spec_pref : list[str]
        Instrument preference list for the reference constant (frozen to 1).
        The first instrument in this list that appears in *indices* is
        selected.  For example, ``['NuSTAR', 'NICER', 'IXPE']`` tries
        NuSTAR first; if absent, NICER; if absent, IXPE.
    audit_data : bool
        If True (default), print the loaded file paths and count rates for each spectrum
        at the beginning of the log file.
    log_path : str or None
        If given, the verbose ``AllModels.show`` / ``Fit.show`` output (at
        ``chatter_level``) is written to this file via ``Xset.openLog``.
        If None, the output is printed to the notebook.
    chatter_level : int
        XSPEC chatter level used for the logged model/fit summary
        (default 10).  Restored to 5 afterwards.
    """
    if audit_data:
        if log_path is not None:
            _audit_buf = io.StringIO()
            with contextlib.redirect_stdout(_audit_buf):
                audit_paths()
        else:
            audit_paths()

    m1 = Model("constant * tbabs * bknpower * (polconst*diskbb + polconst*relxillCp)")

    # Cross-calibration constants: preferred instrument frozen to 1, all others free
    ref_grp = None
    for inst in spec_pref:
        if inst in indices:
            first_obs = next(iter(indices[inst]))
            first_det = next(iter(indices[inst][first_obs]))
            first_idx = indices[inst][first_obs][first_det][0]
            ref_grp = AllData(first_idx).dataGroup
            break
    for dg in range(1, AllData.nGroups + 1):
        c = AllModels(dg).constant.factor
        if dg == ref_grp:
            c.values = 1.0
            c.frozen = True
        else:
            c.link = ""
            c.frozen = False

    # Initial parameter values
    m1.polconst.psi        = 63.0
    m1.polconst.psi.frozen = False
    m1.polconst_5.psi        = -27.0
    m1.polconst_5.psi.frozen = True
    m1.polconst_5.A          = 0.03
    m1.TBabs.nH              = 0.45
    m1.diskbb.Tin            = 0.43
    m1.diskbb.norm           = 1.3e5
    m1.relxillCp.Incl        = 37
    m1.relxillCp.Incl.frozen = True
    m1.relxillCp.a           = 0.998
    m1.relxillCp.a.frozen    = True
    m1.relxillCp.Rin.values    = [-2.5, 0.05, -6.0, -6.0, -1.0, -1.0]
    m1.relxillCp.Rin.frozen    = False
    m1.relxillCp.Index1.values = [3.0, 0.1, 3.0, 3.0, 5.0, 9.0]
    m1.relxillCp.Index1.frozen = True
    m1.relxillCp.gamma         = 2.0
    m1.relxillCp.gamma.frozen  = False
    m1.relxillCp.logxi         = 3.0
    m1.relxillCp.logN          = 15.0
    m1.relxillCp.logN.frozen   = True
    m1.relxillCp.Afe.values    = [1.0, 0.1, 0.5, 0.5, 2.0, 3.0]
    m1.relxillCp.Afe.frozen    = True
    m1.relxillCp.kTe           = 300.0
    m1.relxillCp.kTe.frozen    = True
    m1.relxillCp.refl_frac     = 1.0

    AllModels.show()
    Fit.perform()
    print(f"Fit statistic: {Fit.statistic}, dof: {Fit.dof}")

    # Verbose model+fit summary, redirected to log if requested
    Xset.chatter = chatter_level
    try:
        if log_path is not None:
            Xset.openLog(log_path)
        AllModels.show()
        Fit.show()
    finally:
        if log_path is not None:
            Xset.closeLog()
        Xset.chatter = 5

    # Prepend audit output so it appears first in the log
    if audit_data and log_path is not None:
        with open(log_path, 'r') as f:
            xspec_content = f.read()
        with open(log_path, 'w') as f:
            f.write(_audit_buf.getvalue())
            f.write(xspec_content)


def fit_tbfeo_plabs_polconst_diskbb_polconst_relxillCp(
    indices: dict,
    spec_pref: list[str] = ['NuSTAR', 'NICER', 'IXPE'],
    audit_data: bool = True,
    log_path: str | None = None,
    chatter_level: int = 10,
) -> None:
    """Define and fit the model ``constant * tbfeo * plabs * (polconst*diskbb + polconst*relxillCp)``.

    Sets up the model, assigns initial parameter values, performs the fit, and
    optionally logs the verbose model/fit summary to a file.  The reference
    instrument selected via *spec_pref* has both its cross-calibration
    ``constant.factor`` **and** its ``plabs.index`` frozen; all other data
    groups have both freed.  ``plabs.coef`` is frozen at 1.0 on DG 1 (the
    default-link master).

    Parameters
    ----------
    indices : dict
        Return value of :func:`load_data`.  Used to look up data-group
        numbers for each instrument.
    spec_pref : list[str]
        Instrument preference list for the reference data group.  The first
        instrument in this list that appears in *indices* is selected; its
        ``constant.factor`` is frozen to 1 and its ``plabs.index`` is frozen
        at the model default.  For example, ``['NuSTAR', 'NICER', 'IXPE']``
        tries NuSTAR first; if absent, NICER; if absent, IXPE.
    audit_data : bool
        If True (default), print the loaded file paths and count rates for
        each spectrum at the beginning of the log file.
    log_path : str or None
        If given, the verbose ``AllModels.show`` / ``Fit.show`` output (at
        ``chatter_level``) is written to this file via ``Xset.openLog``.
        If None, the output is printed to the notebook.
    chatter_level : int
        XSPEC chatter level used for the logged model/fit summary
        (default 10).  Restored to 5 afterwards.
    """
    if audit_data:
        if log_path is not None:
            _audit_buf = io.StringIO()
            with contextlib.redirect_stdout(_audit_buf):
                audit_paths()
        else:
            audit_paths()

    m1 = Model("constant * tbfeo * plabs * (polconst*diskbb + polconst*relxillCp)")

    # Pick reference DG from spec_pref
    ref_grp = None
    for inst in spec_pref:
        if inst in indices:
            first_obs = next(iter(indices[inst]))
            first_det = next(iter(indices[inst][first_obs]))
            first_idx = indices[inst][first_obs][first_det][0]
            ref_grp = AllData(first_idx).dataGroup
            break

    # Per-DG cross-calibration constant and plabs.index controls
    for dg in range(1, AllData.nGroups + 1):
        m = AllModels(dg)
        if dg == ref_grp:
            m.constant.factor = 1.0
            m.constant.factor.frozen = True
            m.plabs.index.frozen = True
        else:
            m.constant.factor.link = ""
            m.constant.factor.frozen = False
            m.plabs.index.link = ""
            m.plabs.index.values = [0.0, 0.1, -3.0, -3.0, 3.0, 3.0]
            m.plabs.index.frozen = False

    # plabs.coef: leave the default DG-1 link in place; freeze the master at 1.0
    AllModels(1).plabs.coef = 1.0
    AllModels(1).plabs.coef.frozen = True

    # Initial parameter values (match Fit 5 in thermal_emission.ipynb)
    m1.polconst.psi          = 63.0
    m1.polconst.psi.frozen   = False
    m1.polconst_6.psi        = -27.0
    m1.polconst_6.psi.frozen = True
    m1.polconst_6.A          = 0.03
    m1.TBfeo.nH              = 0.40
    m1.TBfeo.nH.frozen       = True
    m1.diskbb.Tin            = 0.43
    m1.diskbb.norm           = 1.3e5
    m1.relxillCp.Incl        = 37
    m1.relxillCp.Incl.frozen = True
    m1.relxillCp.a           = 0.998
    m1.relxillCp.a.frozen    = True
    m1.relxillCp.Rin.values    = [-2.5, 0.05, -6.0, -6.0, -1.0, -1.0]
    m1.relxillCp.Rin.frozen    = False
    m1.relxillCp.Index1.values = [3.0, 0.1, 3.0, 3.0, 5.0, 9.0]
    m1.relxillCp.Index1.frozen = True
    m1.relxillCp.gamma         = 2.05
    m1.relxillCp.gamma.frozen  = False
    m1.relxillCp.logxi         = 3.0
    m1.relxillCp.logN          = 15.0
    m1.relxillCp.logN.frozen   = True
    m1.relxillCp.Afe           = 7.0
    m1.relxillCp.Afe.frozen    = True
    m1.relxillCp.kTe           = 100.0
    m1.relxillCp.kTe.frozen    = True
    m1.relxillCp.refl_frac     = 1.0

    AllModels.show()
    Fit.perform()
    print(f"Fit statistic: {Fit.statistic}, dof: {Fit.dof}")

    Xset.chatter = chatter_level
    try:
        if log_path is not None:
            Xset.openLog(log_path)
        AllModels.show()
        Fit.show()
    finally:
        if log_path is not None:
            Xset.closeLog()
        Xset.chatter = 5

    if audit_data and log_path is not None:
        with open(log_path, 'r') as f:
            xspec_content = f.read()
        with open(log_path, 'w') as f:
            f.write(_audit_buf.getvalue())
            f.write(xspec_content)


def inspect_epoch(
    indices: dict,
    e_bins,
    title: str = '',
    xlim: tuple[float, float] = (1.0, 79.0),
    xscale: str = 'log',
    yscale: str = 'log',
    plot_type: str = 'ufspec',
    residuals: tuple = ('delchi', 'ratio'),
):
    """Plot the fit result for one epoch: spectrum panel + per-DU Stokes I/Q/U.

    Combines :func:`inspect_fit` and :func:`plot_iqu_decomp` into a single
    call.  ``spec_groups`` and ``iqu_groups`` are built automatically from
    *indices* (the return value of :func:`load_data`), so only instruments
    present in the epoch appear.

    Parameters
    ----------
    indices : dict
        Return value of :func:`load_data`.
    e_bins : array-like
        Energy bin edges in keV passed to :func:`plot_iqu_decomp`.
    title : str
        Figure suptitle for the spectrum panel; also used as the
        ``title_prefix`` for the per-DU IQU figures.
    xlim : tuple[float, float]
        (xmin, xmax) for the shared x-axis of the spectrum panel
        (default 1–79 keV).
    xscale : str
        X-axis scale for the spectrum panel (default ``'log'``).
    yscale : str
        Y-axis scale for the spectrum panel (default ``'log'``).
    plot_type : str
        XSPEC plot type for the main panel (default ``'ufspec'``).
    residuals : tuple of str
        Residual panel types, each ``'delchi'`` or ``'ratio'``
        (default ``('delchi', 'ratio')``).

    Returns
    -------
    fig : matplotlib Figure
        The spectrum + residuals figure from :func:`inspect_fit`.
    axes : list of Axes
        Axes from *fig* (spectrum panel first, then residual panels).
    figs : list of tuples
        Per-DU ``(fig, (ax_i, ax_q, ax_u))`` tuples from
        :func:`plot_iqu_decomp`.  Empty list if IXPE is not in *indices*.
    """
    Plot.xAxis = 'keV'

    # Build spec_groups: IXPE I-only per DU, then NICER, then NuSTAR
    spec_groups = {}
    if 'IXPE' in indices:
        ixpe_specs = []
        for dus in indices['IXPE'].values():
            for du in ('DU1', 'DU2', 'DU3'):
                if du in dus:
                    ixpe_specs.append(dus[du][0])  # I only
        spec_groups['IXPE Det'] = ixpe_specs
    if 'NICER' in indices:
        nicer_specs = []
        for dets in indices['NICER'].values():
            for specs in dets.values():
                nicer_specs.extend(specs)
        spec_groups['NICER'] = nicer_specs
    if 'NuSTAR' in indices:
        nustar_specs = []
        for dets in indices['NuSTAR'].values():
            for det in ('FPMA', 'FPMB'):
                if det in dets:
                    nustar_specs.extend(dets[det])
        spec_groups['NuSTAR'] = nustar_specs

    fig, axes = inspect_fit(
        spec_groups,
        plot_type=plot_type,
        residuals=residuals,
        title=title,
        xlim=xlim,
        xscale=xscale,
        yscale=yscale,
    )

    # Build iqu_groups and plot per-DU Stokes decomposition (one IXPE obs enforced upstream)
    figs = []
    if 'IXPE' in indices:
        iqu_groups = {}
        for dus in indices['IXPE'].values():
            for du in ('DU1', 'DU2', 'DU3'):
                if du in dus:
                    iqu_groups[du] = dus[du]  # [I, Q, U]
        figs = plot_iqu_decomp(iqu_groups, e_bins, title_prefix=title)

    return fig, axes, figs


def audit_paths():
    rows = []
    for i in range(1, AllData.nSpectra + 1):
        s = AllData(i)

        try: bkg = s.background.fileName
        except Exception: bkg = None

        try:
            rmf = s.response.rmf
            arf = s.response.arf
        except Exception:
            rmf, arf = None, None

        # s.rate -> (total, background, net, net_error) cts/s over noticed channels
        try: rate = s.rate
        except Exception: rate = None

        rows.append((i, s.dataGroup, s.fileName, bkg, rmf, arf, rate))

    for r in rows:
        i, dg, pha, bkg, rmf, arf, rate = r
        if rate is not None:
            rate_str = (f"net={rate[2]:.4g} ± {rate[3]:.4g} cts/s  "
                        f"(total={rate[0]:.4g}, bkg={rate[1]:.4g})")
        else:
            rate_str = "N/A"
        print(f"\n[{i}] DG={dg}\n  PHA: {pha}\n  BKG: {bkg}\n  RMF: {rmf}\n  ARF/MRF: {arf}\n  Rate: {rate_str}")