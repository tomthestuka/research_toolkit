"""xspec_utils — XSPEC spectral analysis tools for the auto_research sandbox.

Exposes TOOL_MANIFEST for ToolkitRegistry discovery.  All XSPEC-dependent
imports are deferred to inside tool functions so this module can be imported
even when XSPEC is not on the Python path (e.g. during registry loading in the
orchestration env).
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add lib/ to sys.path at import time so bare imports in xspec_help.py work
# (import NICER, import red_ixpe, etc.) when this package is imported in the
# analysis subprocess.
_lib = Path(__file__).parent / "lib"
if str(_lib) not in sys.path:
    sys.path.insert(0, str(_lib))


def run_joint_spectral_fit(
    epoch: dict,
    ixpe_path: str,
    nicer_path: str | None = None,
    nustar_path: str | None = None,
    model: str = "tbabs_bknpower",
    log_path: str | None = None,
    figure_out: str | None = None,
    e_bins: list | None = None,
    relxill_path: str | None = None,
    ixpe_erange: list | None = None,
    nicer_erange: list | None = None,
    nustar_erange: list | None = None,
) -> dict:
    """Load multi-instrument spectra for one epoch and run a joint spectral fit.

    All XSPEC work (load_data → model setup → fit → optional figure) happens
    inside a single call so the full XSPEC session state is self-contained in
    one sandbox subprocess.

    Parameters
    ----------
    epoch : dict
        Epoch definition: {'IXPE': ['obsid'], 'NICER': ['obsid', ...],
        'NuSTAR': ['obsid', ...], 'STATE': 'HSS'}.  Key order sets XSPEC
        data-group numbering.  Absent instruments are silently skipped.
    ixpe_path : str
        Path to directory containing IXPE observation sub-directories.
    nicer_path : str or None
        Path to NICER data directory (required if epoch has NICER).
    nustar_path : str or None
        Path to NuSTAR data directory (required if epoch has NuSTAR).
    model : str
        Spectral model.  One of:
        - "tbabs_bknpower"  →  tbabs*bknpower*(polconst*diskbb+polconst*relxillCp)
        - "tbfeo_plabs"     →  tbfeo*plabs*(polconst*diskbb+polconst*relxillCp)
    log_path : str or None
        If given, write verbose XSPEC model/fit output to this file.
    figure_out : str or None
        If given, save the broadband spectrum+residuals figure to this PNG path.
    e_bins : list or None
        Energy bin edges in keV for the Stokes IQU decomposition panels.
        Defaults to [2, 3, 4, 5, 6, 8] when IXPE is in the epoch.
    relxill_path : str or None
        Path to the relxill local-model directory.  Passed to
        AllModels.lmod("relxill", relxill_path).  Skip if relxill is already
        in the XSPEC model path.
    ixpe_erange : list or None
        [lo, hi] energy range in keV for IXPE (default [2.0, 8.0]).
    nicer_erange : list or None
        [lo, hi] energy range in keV for NICER (default [1.0, 10.0]).
    nustar_erange : list or None
        [lo, hi] energy range in keV for NuSTAR (default [4.0, 50.0]).

    Returns
    -------
    dict
        JSON-serialisable result with keys:
        - fit_stat : float       — XSPEC fit statistic
        - dof : int              — degrees of freedom
        - stat_per_dof : float   — fit_stat / dof
        - loaded_instruments : list[str]
        - figures : list of FigureRef-compatible dicts (path, caption,
          mime_type, produced_by, created_at, sha256)
    """
    import hashlib
    import os
    import time as _time

    import matplotlib
    matplotlib.use("Agg")  # non-interactive backend — required in subprocess
    import matplotlib.pyplot as plt
    from xspec import AllModels, Fit, Xset

    # Deferred import of lib functions (requires lib/ on sys.path, set at module level).
    from xspec_help import (
        fit_tbabs_bknpower_polconst_diskbb_polconst_relxillCp,
        fit_tbfeo_plabs_polconst_diskbb_polconst_relxillCp,
        inspect_epoch,
        load_data,
    )

    Fit.query = "no"
    Xset.allowPrompting = False

    if relxill_path is not None:
        AllModels.lmod("relxill", relxill_path)

    load_kwargs: dict = {}
    if ixpe_erange:
        load_kwargs["ixpe_erange"] = tuple(ixpe_erange)
    if nicer_erange:
        load_kwargs["nicer_erange"] = tuple(nicer_erange)
    if nustar_erange:
        load_kwargs["nustar_erange"] = tuple(nustar_erange)

    indices = load_data(
        epoch=epoch,
        ixpe_path=ixpe_path,
        nicer_path=nicer_path,
        nustar_path=nustar_path,
        **load_kwargs,
    )

    _FIT_FNS = {
        "tbabs_bknpower": fit_tbabs_bknpower_polconst_diskbb_polconst_relxillCp,
        "tbfeo_plabs": fit_tbfeo_plabs_polconst_diskbb_polconst_relxillCp,
    }
    fit_fn = _FIT_FNS.get(model, fit_tbabs_bknpower_polconst_diskbb_polconst_relxillCp)
    fit_fn(indices, log_path=log_path)

    dof = Fit.dof
    stat = Fit.statistic
    result: dict = {
        "fit_stat": stat,
        "dof": dof,
        "stat_per_dof": stat / dof if dof > 0 else None,
        "loaded_instruments": list(indices.keys()),
        "figures": [],
    }

    if figure_out:
        _e_bins = e_bins or ([2.0, 3.0, 4.0, 5.0, 6.0, 8.0] if "IXPE" in indices else [1.0, 3.0, 6.0, 10.0, 50.0])
        fig, _axes, du_figs = inspect_epoch(indices, _e_bins, title=os.path.basename(figure_out))
        fig.savefig(figure_out, bbox_inches="tight", dpi=150)
        plt.close(fig)
        for du_fig, _ in du_figs:
            plt.close(du_fig)
        sha = hashlib.sha256(Path(figure_out).read_bytes()).hexdigest()
        result["figures"].append({
            "path": figure_out,
            "caption": f"Joint spectral fit ({model}): broadband spectrum + residuals",
            "mime_type": "image/png",
            "produced_by": "xspec_utils.run_joint_spectral_fit",
            "created_at": _time.strftime("%Y-%m-%dT%H:%M:%SZ", _time.gmtime()),
            "sha256": sha,
        })

    return result


def group_pha(infile: str, outfile: str, e_edges: list) -> None:
    """Write a grouped PHA file with custom energy bin edges.

    Reads the RMF path from the PHA RESPFILE header, maps e_edges to channel
    indices via the RMF EBOUNDS extension, and writes the grouped file to
    outfile.  Does not require XSPEC — only astropy and numpy.

    Parameters
    ----------
    infile : str
        Input PHA file path.
    outfile : str
        Output grouped PHA file path.
    e_edges : list of float
        Bin edges in keV, e.g. [2.0, 3.0, 4.0, 5.0, 6.0, 8.0].
    """
    import os

    import numpy as np
    from astropy.io import fits

    with fits.open(infile) as hdul:
        rmf = hdul["SPECTRUM"].header["RESPFILE"]
        if not os.path.isabs(rmf):
            rmf = os.path.join(os.path.dirname(os.path.abspath(infile)), rmf)
        with fits.open(rmf) as rmf_hdul:
            e_lo = rmf_hdul["EBOUNDS"].data["E_MIN"].astype(float)
            e_hi = rmf_hdul["EBOUNDS"].data["E_MAX"].astype(float)
        e_mid = 0.5 * (e_lo + e_hi)
        group = np.zeros(len(e_mid), dtype=np.int16)
        for lo, hi in zip(e_edges[:-1], e_edges[1:]):
            idx = np.where((e_mid >= lo) & (e_mid < hi))[0]
            if len(idx):
                group[idx[0]] = 1
                group[idx[1:]] = -1
        spec = hdul["SPECTRUM"]
        if "GROUPING" in spec.columns.names:
            spec.data["GROUPING"] = group
        else:
            new_col = fits.Column(name="GROUPING", format="I", array=group)
            new_hdu = fits.BinTableHDU.from_columns(
                spec.columns + fits.ColDefs([new_col]), header=spec.header
            )
            hdul["SPECTRUM"] = new_hdu
        hdul.writeto(outfile, overwrite=True)


TOOL_MANIFEST: dict = {
    "xspec_joint_fit": {
        "module": "xspec_utils",
        "fn": "run_joint_spectral_fit",
        "description": (
            "Load multi-instrument X-ray spectra (IXPE / NICER / NuSTAR) for one epoch "
            "and run a joint spectro-polarimetric fit in XSPEC. "
            "Supports model 'tbabs_bknpower' "
            "(constant*tbabs*bknpower*(polconst*diskbb+polconst*relxillCp), tuned for Cyg X-1 HSS) "
            "and 'tbfeo_plabs' "
            "(constant*tbfeo*plabs*(polconst*diskbb+polconst*relxillCp)). "
            "Returns fit_stat, dof, stat_per_dof, loaded_instruments, and optional figures. "
            "The epoch dict must be constructed by the analyst from the local obs DB — "
            "keys 'IXPE', 'NICER', 'NuSTAR' (each a list of obs-ID strings) plus 'STATE'."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "epoch": {
                    "type": "object",
                    "description": (
                        "Epoch definition: "
                        "{'IXPE': ['obsid'], 'NICER': ['obsid', ...], "
                        "'NuSTAR': ['obsid', ...], 'STATE': 'HSS'}. "
                        "Absent instruments are silently skipped. "
                        "Key order determines XSPEC data-group numbering "
                        "(NuSTAR first is conventional)."
                    ),
                },
                "ixpe_path": {
                    "type": "string",
                    "description": "Path to directory containing IXPE observation sub-directories.",
                },
                "nicer_path": {
                    "type": ["string", "null"],
                    "description": "Path to NICER data directory. Required if epoch has 'NICER'.",
                },
                "nustar_path": {
                    "type": ["string", "null"],
                    "description": "Path to NuSTAR data directory. Required if epoch has 'NuSTAR'.",
                },
                "model": {
                    "type": "string",
                    "enum": ["tbabs_bknpower", "tbfeo_plabs"],
                    "description": "Spectral model to fit. Default 'tbabs_bknpower'.",
                    "default": "tbabs_bknpower",
                },
                "log_path": {
                    "type": ["string", "null"],
                    "description": "If given, write verbose XSPEC model+fit output to this file.",
                },
                "figure_out": {
                    "type": ["string", "null"],
                    "description": "If given, save broadband spectrum+residuals figure to this PNG path.",
                },
                "e_bins": {
                    "type": ["array", "null"],
                    "items": {"type": "number"},
                    "description": (
                        "Energy bin edges in keV for Stokes IQU decomposition panels. "
                        "Defaults to [2, 3, 4, 5, 6, 8] when IXPE is present."
                    ),
                },
                "relxill_path": {
                    "type": ["string", "null"],
                    "description": (
                        "Path to relxill local-model directory. "
                        "Required unless relxill is already loaded in XSPEC model path."
                    ),
                },
                "ixpe_erange": {
                    "type": ["array", "null"],
                    "items": {"type": "number"},
                    "description": "[lo, hi] energy range in keV for IXPE (default [2.0, 8.0]).",
                },
                "nicer_erange": {
                    "type": ["array", "null"],
                    "items": {"type": "number"},
                    "description": "[lo, hi] energy range in keV for NICER (default [1.0, 10.0]).",
                },
                "nustar_erange": {
                    "type": ["array", "null"],
                    "items": {"type": "number"},
                    "description": "[lo, hi] energy range in keV for NuSTAR (default [4.0, 50.0]).",
                },
            },
            "required": ["epoch", "ixpe_path"],
        },
    },
    "xspec_group_pha": {
        "module": "xspec_utils",
        "fn": "group_pha",
        "description": (
            "Write a grouped copy of a PHA spectrum file with custom energy bin edges "
            "derived from the file's own RMF EBOUNDS extension. "
            "Does not require XSPEC — uses only astropy and numpy. "
            "Load the output file into XSPEC instead of the original when coarser binning is needed."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "infile": {"type": "string", "description": "Input PHA file path."},
                "outfile": {"type": "string", "description": "Output grouped PHA file path."},
                "e_edges": {
                    "type": "array",
                    "items": {"type": "number"},
                    "description": "Bin edges in keV, e.g. [2.0, 3.0, 4.0, 5.0, 6.0, 8.0].",
                },
            },
            "required": ["infile", "outfile", "e_edges"],
        },
    },
}
