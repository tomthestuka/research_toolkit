"""XSPEC session wrappers for the auto_research sandbox REPL subprocess.

Imported ONLY inside the persistent REPL subprocess (where xspec and HEASoft libs
are on the Python path). Never imported in the orchestration environment.

Wraps PyXSPEC calls that return serializable dicts for the sandbox JSON-line
protocol. Each plotting function produces exactly one PNG per call (or one per DU
for plot_iqu) so figure refs map cleanly to single paths.
matplotlib objects never cross the subprocess boundary.
"""

from __future__ import annotations

import csv
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
) -> dict:
    """Time-resolved 2-8 keV IXPE PCUBE polarimetry for one local ObsID.

    Locates the ObsID's Level-2 DU1/2/3 event files under
    ``<obs_db_path>/<obsid>/event_l2/`` (preferring the instrumental-background-
    rejected ``_rej`` files; falling back to plain ``evt2_v*.fits`` when no
    ``_rej`` variant exists), runs ``red_ixpe.pcube_timeseg`` (which gates
    each window on GTI-clipped exposure before reducing it, and checkpoints
    to CSV incrementally so a call near the sandbox time cap can resume next
    call), and records the CSV in the mission catalog at level 3.

    The IRF calibration epoch is auto-selected as the installed epoch nearest
    the observation date -- no epoch argument is needed or accepted here.

    The full per-segment timeseries is never returned inline -- for a long
    observation (O(100) segments) it can exceed the sandbox transport chunk
    limit even though the CSV itself is written fine. Read the CSV at
    ``csv_path`` for the complete table; this wrapper returns only a bounded
    summary (integrated row plus the first/last 3 segment rows).

    Returns
    -------
    dict with 'csv_path', 'n_segments' (rows in the in-memory timeseries,
    resumed + newly computed), 'n_skipped_empty', 'n_expected', 'complete'
    (bool), 'resumed' (bool), 'integrated' (dict, segment_id=-1), 'head'
    and 'tail' (first/last 3 segment-row dicts), 'stdout', 'exit_code'
    (0 on success; 1 with 'stderr' on failure -- nothing is recorded on
    failure).
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
            obsid, file_list, ra=ra, dec=dec,
            emin_kev=emin_kev, emax_kev=emax_kev, src_rad=src_rad, bkg=bkg,
            bkg_inner_rad=bkg_inner_rad, bkg_outer_rad=bkg_outer_rad,
            segment_duration_s=segment_duration_s, min_seg_exposure_s=min_seg_exposure_s,
            resume=resume,
        )
    except (Exception, SystemExit) as exc:
        return {"exit_code": 1, "stdout": "", "stderr": str(exc)}

    csv_path = result["csv_path"]
    if not os.path.exists(csv_path):
        return {"exit_code": 1, "stdout": "",
                "stderr": f"pcube_timeseg did not produce {csv_path}"}

    timeseries = result["timeseries"]
    complete = result["complete"]

    neutral.advance(obs_db_path, "ixpe", obsid, level=3, products=[{
        "kind": "pcube_timeseg",
        "path": os.path.relpath(csv_path, obs_db_path),
        "produced_by": "xspec_utils.red_ixpe.pcube_timeseg",
        "params": {"emin_kev": emin_kev, "emax_kev": emax_kev, "src_rad": src_rad,
                   "bkg": bkg, "segment_duration_s": segment_duration_s,
                   "min_seg_exposure_s": min_seg_exposure_s, "complete": complete},
    }])

    return {
        "csv_path": csv_path,
        "n_segments": len(timeseries),
        "n_skipped_empty": result["n_skipped_empty"],
        "n_expected": result["n_expected"],
        "complete": complete,
        "resumed": result["resumed"],
        "integrated": result["integrated"],
        "head": timeseries[:3],
        "tail": timeseries[-3:],
        "stdout": (f"pcube_timeseg: {len(timeseries)}/{result['n_expected']} segment(s) "
                   f"(+integrated), {result['n_skipped_empty']} skipped, complete={complete}, "
                   f"resumed={result['resumed']} -> {csv_path}"),
        "exit_code": 0,
    }


def ixpe_reduce_pha(
    obs_db_path: str,
    obsid: str,
    src_rad: float = 1.0,
    bkg_inner_rad: float = 3.0,
    bkg_outer_rad: float = 5.0,
) -> dict:
    """Reuse or reduce per-DU IXPE PHA1 (Stokes I) spectra for one local ObsID.

    Locates the ObsID's Level-2 DU1/2/3 event files under
    ``<obs_db_path>/<obsid>/event_l2/`` (preferring the instrumental-background-
    rejected ``_rej`` files; falling back to plain ``evt2_v*.fits``). If the
    per-DU src+bkg PHA1 files already exist on disk, they are reused as-is --
    no reduction, no region files written. Otherwise the source position is
    centroided from the DU1 event file via ``find_src.run_extraction``
    (writing ``src.reg``/``bkg.reg`` under ``<obs_db_path>/<obsid>/``), and
    ``red_ixpe.pha`` is run (``xpselect`` + ``xpbin PHA1`` per DU) with a
    circular source region (``src_rad``) and annular background
    (``bkg_inner_rad``-``bkg_outer_rad``) at that position. The IRF calibration
    epoch is auto-selected as the installed epoch nearest the observation date --
    no epoch argument is needed or accepted here.

    Every PHA1 product on disk after this call (3 DUs x src/bkg, Stokes I) is
    recorded in the mission catalog at level 3 via
    ``obsdb_catalog.neutral.advance``, whether reused or freshly reduced.

    Returns
    -------
    dict with 'src_pha', 'bkg_pha' (each a 3-element list, DU1/DU2/DU3
    Stokes-I PHA1 paths), 'ra', 'dec' (None when reused), 'reused' (bool),
    'src_reg', 'bkg_reg' (None when reused), 'stdout', 'exit_code' (1 with
    'stderr' on failure -- nothing is recorded on failure).
    """
    evt_dir = os.path.join(obs_db_path, obsid, "event_l2")
    file_list = sorted(glob.glob(os.path.join(evt_dir, f"ixpe{obsid}_det[123]_evt2_v[0-9][0-9]_rej.fits")))
    if not file_list:
        file_list = sorted(glob.glob(os.path.join(evt_dir, f"ixpe{obsid}_det[123]_evt2_v[0-9][0-9].fits")))
    if not file_list:
        return {"exit_code": 1, "stdout": "",
                "stderr": f"No Level-2 event files found for {obsid} under {evt_dir}"}

    try:
        result = red_ixpe.reduce_pha(
            obsid, file_list, src_rad=src_rad, bkg_inner_rad=bkg_inner_rad,
            bkg_outer_rad=bkg_outer_rad,
        )
    except (Exception, SystemExit) as exc:
        return {"exit_code": 1, "stdout": "", "stderr": str(exc)}

    src_pha, bkg_pha = result["src_pha"], result["bkg_pha"]
    missing = [f for f in src_pha + bkg_pha if not os.path.exists(f)]
    if missing:
        return {"exit_code": 1, "stdout": "",
                "stderr": f"reduce_pha did not produce: {missing}"}

    du_labels = ("DU1", "DU2", "DU3")
    products = [
        {
            "kind": "pha",
            "path": os.path.relpath(path, obs_db_path),
            "produced_by": "xspec_utils.red_ixpe.pha",
            "params": {"du": du, "stokes": "I", "region": region,
                       "src_rad": src_rad, "bkg_inner_rad": bkg_inner_rad,
                       "bkg_outer_rad": bkg_outer_rad,
                       "irf_date": result.get("irf_date"),
                       "reused": result["reused"]},
        }
        for du, s, b in zip(du_labels, src_pha, bkg_pha)
        for region, path in (("src", s), ("bkg", b))
    ]
    neutral.advance(obs_db_path, "ixpe", obsid, level=3, products=products)

    return {
        "src_pha": src_pha,
        "bkg_pha": bkg_pha,
        "ra": result["ra"],
        "dec": result["dec"],
        "reused": result["reused"],
        "src_reg": result["src_reg"],
        "bkg_reg": result["bkg_reg"],
        "stdout": (f"ixpe_reduce_pha {obsid}: reused={result['reused']} "
                   f"src_pha={src_pha} bkg_pha={bkg_pha}"),
        "exit_code": 0,
    }


def ixpe_pcube(
    obs_db_path: str,
    obs,
    ra: float,
    dec: float,
    ebinning: list | None = None,
    src_rad: float = 1.0,
    bkg_inner_rad: float = 3.0,
    bkg_outer_rad: float = 5.0,
    irf_date: str | None = None,
) -> dict:
    """Energy-resolved IXPE PCUBE polarimetry for one ObsID, or a native
    multi-ObsID combination, background-subtracted.

    Resolves the mission's full Level-2 event file list via
    ``red_ixpe.get_all_file_list(obs_db_path)``, then calls ``red_ixpe.pcube``
    -- ``obs`` a single 9-digit ObsID string reduces just that observation;
    ``obs`` a list of ObsID strings filters and concatenates every listed
    ObsID's Level-2 files and passes them into ONE xpselect/xpbin PCUBE call,
    so the returned cube is the native directly-combined reduction (raw
    Stokes I/Q/U/COUNTS summed across every DU *and* every listed ObsID before
    PD/PA are derived), not a post-hoc merge of separately-computed per-ObsID
    PD/PA. ``ra``/``dec`` are required (no auto-centroid here -- unlike
    ``ixpe_pcube_timeseg``, this tool always needs an explicit region centre,
    especially for a multi-ObsID combination). The IRF calibration epoch is
    auto-selected as the installed epoch nearest the midpoint across every
    requested ObsID -- no epoch argument is AI-facing.

    Reads PD, PD_err, PA, PA_err, QN, QN_err, UN, UN_err, MDP99, SIGNIF,
    E_LO, E_HI straight off the background-subtracted source cube's
    ``polarization()``/attribute accessors (``red_ixpe.pcube_bands`` -- never
    hand-rolled Stokes arithmetic). PD/PD_err/QN/QN_err/UN/UN_err/MDP99 are
    dimensionless fractions (not percent); PA/PA_err are degrees.

    Writes a provenance-headed CSV, one row per energy band (band_id, E_LO,
    E_HI, PD, PD_err, PA, PA_err, QN, QN_err, UN, UN_err, MDP99, SIGNIF) --
    there is no separate "integrated" row (unlike ``ixpe_pcube_timeseg``'s
    ``segment_id=-1``); a single-edge-pair ``ebinning`` (e.g. ``[2, 8]``)
    naturally yields exactly one row spanning the full requested band. A
    single ObsID's CSV is written to
    ``<obs_db_path>/<obsid>/ixpe_pcube_<obsid>.csv``; a combined request's CSV
    is written to ``<obs_db_path>/ixpe_pcube_<obsid1>_<obsid2>_....csv`` (no
    single ObsID directory owns a combined product, mirroring
    ``bhxrb_hardness_ratios.csv``'s placement at the mission root). Records
    one ``kind='pcube'`` product per contributing ObsID via
    ``obsdb_catalog.neutral.advance`` (level 3), all pointing at the same CSV
    for a combined request.

    Parameters
    ----------
    obs_db_path : str
        The IXPE mission's obs-db directory (e.g.
        '/data/fzhou/XSPEC/Cyg_X-1/IXPE'). Passed directly to
        ``red_ixpe.get_all_file_list``.
    obs : str or list[str]
        A single 9-digit IXPE ObsID, or a list of such ObsIDs to natively
        combine.
    ra, dec : float
        Source position in decimal degrees, shared across every listed ObsID.
    ebinning : list[float] or None
        Energy bin edges in keV, ebinalg='LIST'. Default [2, 3, 4, 6, 8]
        (4 bands). Must carry >= 2 edges; a 2-edge list yields one broadband
        row.
    src_rad, bkg_inner_rad, bkg_outer_rad : float
        Source region radius / background annulus radii, in arcminutes.
        Defaults 1.0 / 3.0 / 5.0.

    Returns
    -------
    dict with 'bands' (list of per-band dicts, one per energy bin),
    'provenance' (obsids, ebinning, src_rad, bkg_inner_rad, bkg_outer_rad, ra,
    dec, irf_date, combined), 'csv_path', 'stdout', 'exit_code' (1 with
    'stderr' on failure -- e.g. a requested ObsID has no local event files, is
    missing a DU, or ebinning has fewer than 2 edges -- nothing is written or
    recorded on failure).
    """
    obs_list = [obs] if isinstance(obs, str) else list(obs)
    if not obs_list:
        return {"exit_code": 1, "stdout": "",
                "stderr": "obs resolved to an empty list of ObsIDs"}

    _ebinning = [2.0, 3.0, 4.0, 6.0, 8.0] if ebinning is None else list(ebinning)

    try:
        file_list = red_ixpe.get_all_file_list(obs_db_path)
    except (Exception, SystemExit) as exc:
        return {"exit_code": 1, "stdout": "", "stderr": str(exc)}

    try:
        src_pcube, _bkg_pcube, resolved_irf_date = red_ixpe.pcube(
            obs_list[0] if len(obs_list) == 1 else obs_list, file_list,
            ebinning=_ebinning, src_rad=src_rad, bkg_inner_rad=bkg_inner_rad,
            bkg_outer_rad=bkg_outer_rad, ra=ra, dec=dec, irf_date=irf_date,
        )
    except (Exception, SystemExit) as exc:
        return {"exit_code": 1, "stdout": "", "stderr": str(exc)}

    bands = red_ixpe.pcube_bands(src_pcube)
    combined = len(obs_list) > 1
    provenance = {
        "obsids": sorted(obs_list), "ebinning": _ebinning, "src_rad": src_rad,
        "bkg_inner_rad": bkg_inner_rad, "bkg_outer_rad": bkg_outer_rad,
        "ra": ra, "dec": dec, "irf_date": resolved_irf_date, "combined": combined,
    }

    if combined:
        csv_path = os.path.join(obs_db_path, f"ixpe_pcube_{'_'.join(sorted(obs_list))}.csv")
    else:
        csv_dir = os.path.join(obs_db_path, obs_list[0])
        os.makedirs(csv_dir, exist_ok=True)
        csv_path = os.path.join(csv_dir, f"ixpe_pcube_{obs_list[0]}.csv")

    red_ixpe.write_pcube_csv(csv_path, provenance, bands)

    for obsid in obs_list:
        neutral.advance(obs_db_path, "ixpe", obsid, level=3, products=[{
            "kind": "pcube",
            "path": os.path.relpath(csv_path, obs_db_path),
            "produced_by": "xspec_utils.red_ixpe.pcube",
            "params": dict(provenance),
        }])

    return {
        "bands": bands,
        "provenance": provenance,
        "csv_path": csv_path,
        "stdout": (f"ixpe_pcube: {len(bands)} band(s) for "
                   f"{'+'.join(sorted(obs_list))} -> {csv_path}"),
        "exit_code": 0,
    }


def ixpe_qu_slice_pcubes(
    obs_db_path: str,
    obsid: str,
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
) -> dict:
    """MDP-anchored adaptive q(t)=Q/I, u(t)=U/I intra-observation time series
    for one local IXPE ObsID, with data-driven N-way splitting of a
    merged/multi-pointing observation.

    Locates the ObsID's Level-2 DU1/2/3 event files under
    ``<obs_db_path>/<obsid>/event_l2/`` (preferring the instrumental-
    background-rejected ``_rej`` files; falling back to plain
    ``evt2_v*.fits`` -- the same convention as ``ixpe_pcube_timeseg``/
    ``ixpe_reduce_pha``, which naturally resolves a literature-merged ObsID
    such as 03002599, whose own reduction script skips ``_rej``, to its
    plain files), then calls ``red_ixpe.qu_slice_pcubes`` -- see that
    function's docstring for the full N-way-split / background-decision /
    MDP-anchored-slicing / mean-arrival-timestamp procedure.

    The N-way split is attempted at all only when ``split_multiobs`` (or,
    when omitted, the ``'...99'`` literature-merged-ObsID naming convention)
    judges this obsid merged -- an ordinary single monitoring pointing
    (e.g. Cyg X-1's 02008201) is never assumed merged, however large its own
    GTI gaps happen to be, and is always reduced as one whole-span
    sub-observation.

    Records one ``kind='qu_slices'`` product per sub-observation CSV in the
    mission catalog at level 3 via ``obsdb_catalog.neutral.advance`` (all
    against the real 9-digit ``obsid`` -- a sub-observation key such as
    ``<obsid>a`` is a virtual partition of that same ObsID, not a separate
    catalog entry).

    Returns
    -------
    dict with 'obsid', 'subobs_keys' (list, chronological order), and
    'per_subobs' (sub-obs key -> dict with 'mdp99_obs_pct', 'r_proxy',
    'subtract_background', 'uses_rej', 'n_slice', 'noise_limited',
    'whole_pd_pct', 'whole_pd_err', 'whole_pa_deg', 'whole_signif',
    'csv_path' -- read that CSV for the full per-slice table), plus
    'stdout'/'exit_code' (1 with 'stderr' on failure -- nothing written or
    recorded on failure).
    """
    evt_dir = os.path.join(obs_db_path, obsid, "event_l2")
    file_list = sorted(glob.glob(os.path.join(evt_dir, f"ixpe{obsid}_det[123]_evt2_v[0-9][0-9]_rej.fits")))
    if not file_list:
        file_list = sorted(glob.glob(os.path.join(evt_dir, f"ixpe{obsid}_det[123]_evt2_v[0-9][0-9].fits")))
    if not file_list:
        return {"exit_code": 1, "stdout": "",
                "stderr": f"No Level-2 event files found for {obsid} under {evt_dir}"}

    try:
        result = red_ixpe.qu_slice_pcubes(
            obsid, file_list, ra=ra, dec=dec, emin_kev=emin_kev, emax_kev=emax_kev,
            src_rad=src_rad, bkg=bkg, bkg_inner_rad=bkg_inner_rad, bkg_outer_rad=bkg_outer_rad,
            mdp_target=mdp_target, mdp_frac=mdp_frac, mdp_floor=mdp_floor,
            n_min=n_min, n_max=n_max, split_gap_factor=split_gap_factor,
            split_multiobs=split_multiobs,
        )
    except (Exception, SystemExit) as exc:
        return {"exit_code": 1, "stdout": "", "stderr": str(exc)}

    per_subobs = result["per_subobs"]
    missing = [k for k, v in per_subobs.items() if not os.path.exists(v["csv_path"])]
    if missing:
        return {"exit_code": 1, "stdout": "",
                "stderr": f"qu_slice_pcubes did not produce CSV(s) for sub-obs: {missing}"}

    for key, info in per_subobs.items():
        neutral.advance(obs_db_path, "ixpe", obsid, level=3, products=[{
            "kind": "qu_slices",
            "path": os.path.relpath(info["csv_path"], obs_db_path),
            "produced_by": "xspec_utils.red_ixpe.qu_slice_pcubes",
            "params": {"subobs_key": key, "emin_kev": emin_kev, "emax_kev": emax_kev,
                       "src_rad": src_rad, "subtract_background": info["subtract_background"],
                       "n_slice": info["n_slice"], "noise_limited": info["noise_limited"]},
        }])

    summary = {k: {kk: vv for kk, vv in v.items() if kk != "rows"} for k, v in per_subobs.items()}
    return {
        "obsid": obsid,
        "subobs_keys": result["subobs_keys"],
        "per_subobs": summary,
        "stdout": (f"ixpe_qu_slice_pcubes {obsid}: {len(result['subobs_keys'])} sub-obs, "
                   f"{sum(v['n_slice'] for v in per_subobs.values())} slice(s) total -> "
                   f"{[v['csv_path'] for v in per_subobs.values()]}"),
        "exit_code": 0,
    }


def ixpe_HR(
    obs_db_path: str,
    obsids: list,
    soft_band: list | None = None,
    hard_band: list | None = None,
    bkg: bool = False,
) -> dict:
    """Net Stokes-I hardness ratio HR=(H-S)/(H+S) for a survey of local IXPE ObsIDs.

    Reads directly from each ObsID's existing per-DU PHA1 (Stokes I) src(+bkg)
    products under ``<obs_db_path>/<obsid>/event_l2/`` (preferring the
    instrumental-background-rejected ``_rej`` files, falling back to plain
    ``evt2_v*.fits`` -- the same convention as ``ixpe_reduce_pha``/
    ``ixpe_pcube_timeseg``) via ``red_ixpe.HR``, so PD (from ``ixpe_pcube_timeseg``)
    and HR share identical extraction regions. This never re-reduces PHA -- run
    ``ixpe_reduce_pha`` first for any ObsID missing its PHA1 products.

    An ObsID with no local Level-2 event files, or whose ``red_ixpe.HR`` call raises
    (e.g. a DU present in the local event files has no reduced PHA1 product), is
    recorded in the returned ``errors`` dict and excluded from the CSV -- never a
    dummy row, and never a crash of the rest of the survey.

    Writes ``bhxrb_hardness_ratios.csv`` (columns obsid, soft_rate, soft_err,
    hard_rate, hard_err, HR, HR_err) at ``<obs_db_path>/bhxrb_hardness_ratios.csv``,
    overwriting any previous run, and records one ``kind='hardness_ratio'`` product
    per successful ObsID via ``obsdb_catalog.neutral.advance`` (level 3), with
    ``params`` carrying the real ``soft_band``/``hard_band``/``bkg``/``dus_used``
    used for that ObsID.

    Parameters
    ----------
    obs_db_path : str
        The IXPE mission's obs-db directory.
    obsids : list[str]
        9-digit IXPE Observation IDs to include in the survey.
    soft_band, hard_band : list[float] or None
        [lo, hi] energy band edges in keV. Default [2.0, 4.0] / [4.0, 8.0].
    bkg : bool
        Apply background subtraction (per DU, BACKSCAL-scaled) for every ObsID in
        the survey. Default False (bright source-dominated targets); set True only
        where background is a non-negligible fraction of the source aperture.

    Returns
    -------
    dict with 'csv_path', 'rows' (the successful per-obsid row dicts written to the
    CSV), 'errors' (obsid -> error message, for ObsIDs that failed), 'stdout',
    'exit_code' (1 with 'stderr' set only if EVERY ObsID failed -- nothing is
    written or recorded in that case; 0 otherwise, with any per-ObsID failures
    surfaced in 'errors').
    """
    soft = tuple(soft_band) if soft_band else (2.0, 4.0)
    hard = tuple(hard_band) if hard_band else (4.0, 8.0)

    rows, errors = [], {}
    for obsid in obsids:
        evt_dir = os.path.join(obs_db_path, obsid, "event_l2")
        file_list = sorted(glob.glob(os.path.join(evt_dir, f"ixpe{obsid}_det[123]_evt2_v[0-9][0-9]_rej.fits")))
        if not file_list:
            file_list = sorted(glob.glob(os.path.join(evt_dir, f"ixpe{obsid}_det[123]_evt2_v[0-9][0-9].fits")))
        if not file_list:
            errors[obsid] = f"No Level-2 event files found for {obsid} under {evt_dir}"
            continue
        try:
            result = red_ixpe.HR(obsid, file_list, soft_band=soft, hard_band=hard, bkg=bkg)
        except (Exception, SystemExit) as exc:
            errors[obsid] = str(exc)
            continue
        rows.append({
            "obsid": obsid,
            "soft_rate": result["soft_rate"], "soft_err": result["soft_err"],
            "hard_rate": result["hard_rate"], "hard_err": result["hard_err"],
            "HR": result["HR"], "HR_err": result["HR_err"],
            "_dus_used": result["dus_used"],
        })

    if not rows:
        return {"exit_code": 1, "stdout": "",
                "stderr": f"No ObsID produced a hardness ratio; errors={errors}"}

    csv_path = os.path.join(obs_db_path, "bhxrb_hardness_ratios.csv")
    fieldnames = ["obsid", "soft_rate", "soft_err", "hard_rate", "hard_err", "HR", "HR_err"]
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in fieldnames})

    for row in rows:
        neutral.advance(obs_db_path, "ixpe", row["obsid"], level=3, products=[{
            "kind": "hardness_ratio",
            "path": os.path.relpath(csv_path, obs_db_path),
            "produced_by": "xspec_utils.red_ixpe.HR",
            "params": {"soft_band": list(soft), "hard_band": list(hard), "bkg": bkg,
                       "dus_used": row["_dus_used"], "HR": row["HR"], "HR_err": row["HR_err"]},
        }])

    return {
        "csv_path": csv_path,
        "rows": [{k: v for k, v in row.items() if k != "_dus_used"} for row in rows],
        "errors": errors,
        "stdout": (f"ixpe_HR: {len(rows)}/{len(obsids)} ObsID(s) -> {csv_path} "
                   f"({len(errors)} error(s))"),
        "exit_code": 0,
    }
