"""XSPEC session wrappers for the auto_research sandbox REPL subprocess.

Imported ONLY inside the persistent REPL subprocess (where xspec and HEASoft libs
are on the Python path). Never imported in the orchestration environment.

Wraps PyXSPEC calls that return serializable dicts for the sandbox JSON-line
protocol. Each plotting function produces exactly one PNG per call (or one per DU
for plot_iqu) so figure refs map cleanly to single paths.
matplotlib objects never cross the subprocess boundary.
"""

from __future__ import annotations

import glob
import hashlib
import os
import time as _time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from xspec import AllData, AllModels, Fit, Xset

import xspec_help  # human-authored lib; on sys.path via xspec_utils/__init__.py
import red_ixpe

from obsdb_catalog import neutral


def run_fit(n_iterations: int = 100, log_path: str | None = None) -> dict:
    """Run XSPEC Fit.perform() on the current session and return fit statistics.

    Call after load_data + fit_tbabs_hss (or fit_tbfeo_hss) to refine
    convergence. May be called multiple times until stat/dof stabilises.

    Parameters
    ----------
    n_iterations : int
        XSPEC Fit.nIterations. Default 100.
    log_path : str or None
        If given, write verbose XSPEC fit output to this file.

    Returns
    -------
    dict with fit_stat, dof, stat_per_dof, stdout, exit_code.
    """
    if log_path:
        Xset.openLog(log_path)
    try:
        Fit.nIterations = n_iterations
        Fit.perform()
    finally:
        if log_path:
            Xset.closeLog()

    dof = Fit.dof
    stat = Fit.statistic
    return {
        "stdout": f"Fit complete: stat={stat:.4f} dof={dof} stat/dof={stat/dof:.4f}",
        "fit_stat": stat,
        "dof": dof,
        "stat_per_dof": stat / dof if dof > 0 else None,
        "exit_code": 0,
    }


def xspec_status() -> dict:
    """Return the current XSPEC session state: loaded spectra, model, fit statistics.

    Call at any point to verify session progress. Returns a structured dict.
    Use xspec_help.audit_paths() for the full file-path / count-rate table.

    Returns
    -------
    dict with n_spectra, model_expression, fit_stat, dof, stat_per_dof,
    stdout, exit_code.
    """
    n_groups = AllData.nGroups
    if n_groups == 0:
        return {
            "stdout": "No spectra loaded.",
            "n_spectra": 0,
            "model_expression": None,
            "fit_stat": None,
            "dof": None,
            "exit_code": 0,
        }

    try:
        model_expr = AllModels(1).expression
    except Exception:
        model_expr = None

    try:
        stat = Fit.statistic
        dof = Fit.dof
    except Exception:
        stat = None
        dof = None

    lines = [f"Loaded spectra: {n_groups} group(s)"]
    if model_expr:
        lines.append(f"Model: {model_expr}")
    else:
        lines.append("No model set yet.")
    if stat is not None and dof and dof > 0:
        lines.append(f"Fit: stat={stat:.4f}  dof={dof}  stat/dof={stat/dof:.4f}")

    return {
        "stdout": "\n".join(lines),
        "n_spectra": n_groups,
        "model_expression": model_expr,
        "fit_stat": stat,
        "dof": dof,
        "stat_per_dof": stat / dof if (stat is not None and dof and dof > 0) else None,
        "exit_code": 0,
    }


def plot_fit(
    indices: dict,
    save_path: str | None = None,
    title: str | None = None,
    plot_type: str = "ufspec",
    residuals: list | None = None,
    xlim: list | None = None,
) -> dict:
    """Save a broadband spectrum + residuals panel for the current XSPEC session.

    Wraps xspec_help.inspect_fit. Accepts the indices dict returned by load_data
    to build the per-instrument spectrum groups (IXPE I-only, NICER, NuSTAR).
    Produces exactly one PNG. Call after load_data and at least one fit step.

    Parameters
    ----------
    indices : dict
        Return value of load_data.
    save_path : str or None
        Full output PNG path. Auto-generated timestamped file if None.
    title : str or None
        Figure suptitle.
    plot_type : str
        XSPEC plot type: 'ufspec' or 'data'. Default 'ufspec'.
    residuals : list of str or None
        Residual panel types ('delchi' or 'ratio'). Default ['delchi'].
    xlim : list [xmin, xmax] or None
        Energy range in keV for the shared x-axis.

    Returns
    -------
    dict with figures (one entry), stdout, exit_code.
    """
    if AllData.nGroups == 0:
        return {"exit_code": 1, "stderr": "No spectra loaded. Call load_data first.", "stdout": ""}

    # Build spec_groups from indices (IXPE I-only, then NICER, then NuSTAR)
    spec_groups: dict = {}
    if "IXPE" in indices:
        ixpe = [dus[du][0] for dus in indices["IXPE"].values()
                for du in ("DU1", "DU2", "DU3") if du in dus]
        if ixpe:
            spec_groups["IXPE Det"] = ixpe
    if "NICER" in indices:
        nicer = [s for dets in indices["NICER"].values() for specs in dets.values() for s in specs]
        if nicer:
            spec_groups["NICER"] = nicer
    if "NuSTAR" in indices:
        nustar = [s for dets in indices["NuSTAR"].values()
                  for det in ("FPMA", "FPMB") if det in dets for s in dets[det]]
        if nustar:
            spec_groups["NuSTAR"] = nustar

    if not spec_groups:
        return {"exit_code": 1, "stderr": "No instruments found in indices.", "stdout": ""}

    if save_path is None:
        ts = _time.strftime("%Y%m%dT%H%M%S")
        save_path = os.path.join(os.getcwd(), f"fit_{ts}.png")
    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)

    fig, _axes = xspec_help.inspect_fit(
        spec_groups,
        plot_type=plot_type,
        residuals=tuple(residuals) if residuals else ("delchi",),
        title=title or "",
        xlim=tuple(xlim) if xlim else None,
    )
    fig.savefig(save_path, bbox_inches="tight", dpi=150)
    plt.close(fig)

    sha = hashlib.sha256(Path(save_path).read_bytes()).hexdigest()
    return {
        "stdout": f"Figure saved: {save_path}",
        "figures": [{
            "path": save_path,
            "caption": title or "Broadband spectrum and residuals",
            "mime_type": "image/png",
            "produced_by": "xspec_utils.skills.session.plot_fit",
            "created_at": _time.strftime("%Y-%m-%dT%H:%M:%SZ", _time.gmtime()),
            "sha256": sha,
        }],
        "exit_code": 0,
    }


def plot_iqu(
    indices: dict,
    e_bins: list | None = None,
    save_path_prefix: str | None = None,
    title_prefix: str | None = None,
) -> dict:
    """Save per-DU Stokes I/Q/U decomposition panels for the current XSPEC session.

    Wraps xspec_help.plot_iqu_decomp. Accepts the indices dict returned by
    load_data to build per-DU IQU groups. Produces one PNG per IXPE DU present
    (DU1, DU2, DU3). Call after load_data and at least one fit step.

    Parameters
    ----------
    indices : dict
        Return value of load_data (must contain 'IXPE').
    e_bins : list of float or None
        Energy bin edges in keV. Default [2, 3, 4, 5, 6, 8].
    save_path_prefix : str or None
        Base path; '_DU1.png', '_DU2.png', etc. are appended per DU.
        Auto-generated timestamped prefix if None.
    title_prefix : str or None
        Prepended to each figure suptitle.

    Returns
    -------
    dict with figures (one entry per DU), stdout, exit_code.
    """
    if AllData.nGroups == 0:
        return {"exit_code": 1, "stderr": "No spectra loaded. Call load_data first.", "stdout": ""}

    # Build iqu_groups from indices (one [I, Q, U] list per DU)
    iqu_groups: dict = {}
    if "IXPE" in indices:
        for dus in indices["IXPE"].values():
            for du in ("DU1", "DU2", "DU3"):
                if du in dus:
                    iqu_groups[du] = dus[du]

    if not iqu_groups:
        return {"exit_code": 1, "stderr": "No IXPE data found in indices.", "stdout": ""}

    _e_bins = e_bins or [2.0, 3.0, 4.0, 5.0, 6.0, 8.0]
    _title_prefix = title_prefix or ""

    if save_path_prefix is None:
        ts = _time.strftime("%Y%m%dT%H%M%S")
        save_path_prefix = os.path.join(os.getcwd(), f"iqu_{ts}")
    os.makedirs(os.path.dirname(os.path.abspath(save_path_prefix)), exist_ok=True)

    du_figs = xspec_help.plot_iqu_decomp(iqu_groups, _e_bins, title_prefix=_title_prefix)
    du_names = list(iqu_groups.keys())

    figures = []
    saved = []
    for i, (du_fig, _axes) in enumerate(du_figs):
        du_label = du_names[i] if i < len(du_names) else f"DU{i + 1}"
        path = f"{save_path_prefix}_{du_label}.png"
        du_fig.savefig(path, bbox_inches="tight", dpi=150)
        plt.close(du_fig)
        sha = hashlib.sha256(Path(path).read_bytes()).hexdigest()
        caption = f"Stokes IQU — {du_label}" + (f" ({_title_prefix})" if _title_prefix else "")
        figures.append({
            "path": path,
            "caption": caption,
            "mime_type": "image/png",
            "produced_by": "xspec_utils.skills.session.plot_iqu",
            "created_at": _time.strftime("%Y-%m-%dT%H:%M:%SZ", _time.gmtime()),
            "sha256": sha,
        })
        saved.append(path)

    return {
        "stdout": f"Saved {len(saved)} IQU figure(s): {', '.join(saved)}",
        "figures": figures,
        "exit_code": 0,
    }


def ixpe_pcube_timeseg(
    obs_db_path: str,
    obsid: str,
    irf_date: str = "20230702",
    ra: float | None = None,
    dec: float | None = None,
    emin_kev: float = 2.0,
    emax_kev: float = 8.0,
    src_rad: float = 1.0,
    bkg: bool = False,
    bkg_inner_rad: float = 3.0,
    bkg_outer_rad: float = 5.0,
    segment_duration_s: float = 5000.0,
) -> dict:
    """Time-resolved 2-8 keV IXPE PCUBE polarimetry for one local ObsID.

    Locates the ObsID's Level-2 DU1/2/3 event files under
    ``<obs_db_path>/<obsid>/event_l2/`` (preferring the instrumental-background-
    rejected ``_rej`` files; falling back to plain ``evt2_v*.fits`` when no
    ``_rej`` variant exists), runs ``red_ixpe.pcube_timeseg``, and records the
    written CSV in the mission catalog at level 3.

    Returns
    -------
    dict with 'timeseries', 'integrated', 'csv_path' (the three ``lib`` keys),
    plus 'stdout' and 'exit_code' (0 on success; 1 with 'stderr' on failure —
    nothing is recorded on failure).
    """
    evt_dir = os.path.join(obs_db_path, obsid, "event_l2")
    # v[0-9][0-9] anchored immediately against .fits (or _rej.fits) -- NOT
    # v*.fits -- so this matches only original Level-2 event files, never this
    # tool's own generated products (e.g. ..._v03_seg000src.fits,
    # ..._v03_intsrc_pcube.fits), which have extra text between the version
    # digits and .fits and would otherwise be rediscovered as inputs on a
    # second run. Mirrors red_ixpe.get_all_file_list's version convention.
    file_list = sorted(glob.glob(os.path.join(evt_dir, f"ixpe{obsid}_det[123]_evt2_v[0-9][0-9]_rej.fits")))
    if not file_list:
        file_list = sorted(glob.glob(os.path.join(evt_dir, f"ixpe{obsid}_det[123]_evt2_v[0-9][0-9].fits")))
    if not file_list:
        return {"exit_code": 1, "stdout": "",
                "stderr": f"No Level-2 event files found for {obsid} under {evt_dir}"}

    try:
        result = red_ixpe.pcube_timeseg(
            obsid, file_list, irf_date, ra=ra, dec=dec,
            emin_kev=emin_kev, emax_kev=emax_kev, src_rad=src_rad, bkg=bkg,
            bkg_inner_rad=bkg_inner_rad, bkg_outer_rad=bkg_outer_rad,
            segment_duration_s=segment_duration_s,
        )
    except (Exception, SystemExit) as exc:
        return {"exit_code": 1, "stdout": "", "stderr": str(exc)}

    csv_path = result["csv_path"]
    if not os.path.exists(csv_path):
        return {"exit_code": 1, "stdout": "",
                "stderr": f"pcube_timeseg did not produce {csv_path}"}

    neutral.advance(obs_db_path, "ixpe", obsid, level=3, products=[{
        "kind": "pcube_timeseg",
        "path": os.path.relpath(csv_path, obs_db_path),
        "produced_by": "xspec_utils.red_ixpe.pcube_timeseg",
        "params": {"emin_kev": emin_kev, "emax_kev": emax_kev, "src_rad": src_rad,
                   "bkg": bkg, "segment_duration_s": segment_duration_s},
    }])

    return {
        "timeseries": result["timeseries"],
        "integrated": result["integrated"],
        "csv_path": csv_path,
        "stdout": f"pcube_timeseg: {len(result['timeseries'])} segment(s) + integrated -> {csv_path}",
        "exit_code": 0,
    }
