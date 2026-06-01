"""xspec_utils.skills — AI-dedicated tool interface for the auto_research sandbox.

Exposes TOOL_MANIFEST for ToolkitRegistry discovery.  Pure data module: no xspec
imports here so the registry can load this in the orchestration environment (which
has no XSPEC/HEASoft).  All XSPEC-dependent code lives in session.py, which is
only imported inside the persistent REPL subprocess.

Human-authored domain code in lib/ (xspec_help, NICER, NuSTAR, red_ixpe, …) is
never duplicated here — tool entries that wrap lib/ functions point directly at
"xspec_utils.lib.xspec_help" as their module.

Workflow (all steps share the persistent XSPEC session):
  1. load_data      — load multi-instrument spectra; returns indices dict
  2. fit_tbabs_hss  — set Cyg X-1 HSS tbabs model and fit
     fit_tbfeo_hss  — alternative tbfeo*plabs model
  3. run_fit        — refine convergence with Fit.perform()
  4. xspec_status   — query session state (spectra, model, fit stats)
  5. plot_fit       — save broadband spectrum+residuals PNG (one file)
     plot_iqu       — save per-DU Stokes IQU panels (one file per DU)
  6. group_pha      — pre-process PHA grouping (no XSPEC required)
"""

from __future__ import annotations

TOOL_MANIFEST: dict = {
    # ── Step 1: load data ─────────────────────────────────────────────────
    "research_toolkit.xspec_utils.load_data": {
        "module": "xspec_utils.lib.xspec_help",
        "fn": "load_data",
        "description": (
            "Load all instruments for one epoch into XSPEC and return a nested indices dict. "
            "Accepts an epoch dict with keys 'IXPE', 'NICER', 'NuSTAR' (each a list of obs-ID "
            "strings) and 'STATE'. Absent instruments are silently skipped. "
            "Data-group order follows key insertion order (convention: NuSTAR → NICER → IXPE). "
            "Returns dict: instrument → obsid → detector → [XSPEC spectrum indices]. "
            "Pass this indices dict unchanged to fit_*, xspec_status, plot_fit, and plot_iqu."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "epoch": {
                    "type": "string",
                    "description": (
                        "JSON-encoded epoch dict, e.g. "
                        "'{\"NuSTAR\":[\"80902318002\"],\"NICER\":[\"6643010101\",\"6643010102\"],"
                        "\"IXPE\":[\"02008401\"],\"STATE\":\"HSS\"}'. "
                        "Key order determines data-group numbering."
                    ),
                },
                "ixpe_path": {
                    "type": "string",
                    "description": (
                        "Directory containing IXPE observation sub-directories, "
                        "e.g. '/data/fzhou/XSPEC/Cyg_X-1/IXPE'."
                    ),
                },
                "nicer_path": {
                    "type": "string",
                    "description": (
                        "NICER data directory, e.g. '/data/fzhou/XSPEC/Cyg_X-1/NICER'. "
                        "Required if epoch has 'NICER'."
                    ),
                },
                "nustar_path": {
                    "type": "string",
                    "description": (
                        "NuSTAR data directory, e.g. '/data/fzhou/XSPEC/Cyg_X-1/NuSTAR'. "
                        "Required if epoch has 'NuSTAR'."
                    ),
                },
                "clear_existing": {
                    "type": "string",
                    "description": (
                        "JSON bool. Call AllData.clear() + AllModels.clear() before loading. "
                        "Default 'true'."
                    ),
                },
                "chatter": {
                    "type": "string",
                    "description": "JSON bool. Print file paths and loading summaries. Default 'false'.",
                },
            },
            "required": ["epoch", "ixpe_path"],
        },
    },

    # ── Step 2a: Cyg X-1 HSS tbabs model ──────────────────────────────────
    "research_toolkit.xspec_utils.fit_tbabs_hss": {
        "module": "xspec_utils.lib.xspec_help",
        "fn": "fit_tbabs_bknpower_polconst_diskbb_polconst_relxillCp",
        "description": (
            "Define and fit constant*tbabs*(polconst*diskbb + polconst*relxillCp) "
            "with Cyg X-1 HSS initial values. "
            "Must be called after load_data. Pass the indices dict returned by load_data. "
            "Sets nH=0.45, Tin=0.43 keV, diskbb psi=63 deg, relxillCp i=37 deg (frozen), "
            "a=0.998 (frozen). NuSTAR cross-calibration constant frozen to 1; all others freed. "
            "Returns None; call xspec_status or run_fit after."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "indices": {
                    "type": "string",
                    "description": "JSON-encoded return value of load_data.",
                },
                "log_path": {
                    "type": "string",
                    "description": "Write verbose XSPEC model/fit output to this file path.",
                },
                "audit_data": {
                    "type": "string",
                    "description": "JSON bool. Print file paths and count rates at start. Default 'true'.",
                },
            },
            "required": ["indices"],
        },
    },

    # ── Step 2b: Cyg X-1 HSS tbfeo model ──────────────────────────────────
    "research_toolkit.xspec_utils.fit_tbfeo_hss": {
        "module": "xspec_utils.lib.xspec_help",
        "fn": "fit_tbfeo_plabs_polconst_diskbb_polconst_relxillCp",
        "description": (
            "Define and fit constant*tbfeo*plabs*(polconst*diskbb + polconst*relxillCp) "
            "for Cyg X-1 HSS — alternative to fit_tbabs_hss with polarized absorption. "
            "Same calling convention as fit_tbabs_hss. Pass indices from load_data."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "indices": {
                    "type": "string",
                    "description": "JSON-encoded return value of load_data.",
                },
                "log_path": {"type": "string"},
                "audit_data": {"type": "string"},
            },
            "required": ["indices"],
        },
    },

    # ── Step 3: refine fit ─────────────────────────────────────────────────
    "research_toolkit.xspec_utils.run_fit": {
        "module": "xspec_utils.skills.session",
        "fn": "run_fit",
        "description": (
            "Run XSPEC Fit.perform() to refine convergence after fit_tbabs_hss or fit_tbfeo_hss. "
            "Call multiple times until stat/dof stabilises. "
            "Returns dict with fit_stat, dof, stat_per_dof."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "n_iterations": {
                    "type": "string",
                    "description": "JSON int. XSPEC Fit.nIterations. Default 100.",
                },
                "log_path": {
                    "type": "string",
                    "description": "Write verbose XSPEC fit output to this file path.",
                },
            },
            "required": [],
        },
    },

    # ── Session status ─────────────────────────────────────────────────────
    "research_toolkit.xspec_utils.xspec_status": {
        "module": "xspec_utils.skills.session",
        "fn": "xspec_status",
        "description": (
            "Return current XSPEC session state: n_spectra, model_expression, "
            "fit_stat, dof, stat_per_dof. Call at any point to verify session progress."
        ),
        "args_schema": {"type": "object", "properties": {}, "required": []},
    },

    # ── Plots ──────────────────────────────────────────────────────────────
    "research_toolkit.xspec_utils.plot_fit": {
        "module": "xspec_utils.skills.session",
        "fn": "plot_fit",
        "description": (
            "Save a broadband spectrum+residuals panel for the current XSPEC session. "
            "Accepts the indices dict from load_data to build per-instrument groups. "
            "Produces exactly one PNG. "
            "Call after load_data and at least one fit step."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "indices": {
                    "type": "string",
                    "description": "JSON-encoded return value of load_data.",
                },
                "save_path": {
                    "type": "string",
                    "description": "Full output PNG path. Auto-generated if omitted.",
                },
                "title": {"type": "string"},
                "plot_type": {
                    "type": "string",
                    "description": "'ufspec' or 'data'. Default 'ufspec'.",
                },
                "residuals": {
                    "type": "string",
                    "description": "JSON list of residual types, e.g. '[\"delchi\"]'. Default ['delchi'].",
                },
                "xlim": {
                    "type": "string",
                    "description": "JSON [xmin, xmax] in keV, e.g. '[1, 79]'.",
                },
            },
            "required": ["indices"],
        },
    },

    "research_toolkit.xspec_utils.plot_iqu": {
        "module": "xspec_utils.skills.session",
        "fn": "plot_iqu",
        "description": (
            "Save per-DU Stokes I/Q/U decomposition panels for the current XSPEC session. "
            "Accepts the indices dict from load_data; requires IXPE data. "
            "Produces one PNG per DU (DU1, DU2, DU3). "
            "Returns a figures list with one entry per DU. "
            "Call after load_data and at least one fit step."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "indices": {
                    "type": "string",
                    "description": "JSON-encoded return value of load_data (must contain 'IXPE').",
                },
                "e_bins": {
                    "type": "string",
                    "description": "JSON list of energy bin edges in keV. Default '[2,3,4,5,6,8]'.",
                },
                "save_path_prefix": {
                    "type": "string",
                    "description": "Base path; '_DU1.png', '_DU2.png', etc. appended per DU.",
                },
                "title_prefix": {
                    "type": "string",
                    "description": "Prepended to each figure suptitle.",
                },
            },
            "required": ["indices"],
        },
    },

    # ── PHA grouping (no XSPEC required) ──────────────────────────────────
    "research_toolkit.xspec_utils.group_pha": {
        "module": "xspec_utils.lib.xspec_help",
        "fn": "group_pha",
        "description": (
            "Write a grouped copy of a PHA file with custom energy bin edges "
            "derived from the file's own RMF EBOUNDS. No XSPEC required. "
            "Load outfile into XSPEC instead of infile when coarser binning is needed."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "infile": {"type": "string", "description": "Input PHA file path."},
                "outfile": {"type": "string", "description": "Output grouped PHA file path."},
                "e_edges": {
                    "type": "string",
                    "description": "JSON list of bin edges in keV, e.g. '[2.0,3.0,4.0,5.0,6.0,8.0]'.",
                },
            },
            "required": ["infile", "outfile", "e_edges"],
        },
    },
}
