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
  2. fit_tbfeo_hss  — set Cyg X-1 HSS tbfeo*plabs model and fit
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
            "Data-group order follows the epoch key order. LIST IXPE FIRST so IXPE DU1 is "
            "data group 1 — the polarized fit models require DG1 to be an IXPE (I/Q/U) group; "
            "a non-IXPE DG1 makes the fit diverge. Recommended key order: IXPE, NICER, NuSTAR "
            "(skip any not present). This is NOT the same as the cross-calibration reference, "
            "which fit_* selects internally (NuSTAR by default). "
            "Returns dict: instrument → obsid → detector → [XSPEC spectrum indices]. "
            "Pass this indices dict unchanged to fit_*, xspec_status, plot_fit, and plot_iqu. "
            "On first load XSPEC prints harmless warnings about unresolved response-file "
            "paths (arf/rmf are loaded explicitly) — ignore them. For bright sources like "
            "Cyg X-1, NuSTAR/NICER background is None by design — not an error or missing data."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "epoch": {
                    "type": "object",
                    "additionalProperties": True,
                    "description": (
                        "Epoch dict (pass the object directly — do NOT JSON-encode it). "
                        "List IXPE FIRST (IXPE = data group 1), then NICER, then NuSTAR; "
                        "key order determines data-group numbering. Keys 'IXPE'/'NICER'/'NuSTAR' "
                        "map to lists of obs-ID strings; 'STATE' is a label. E.g. "
                        "{\"IXPE\":[\"02008401\"],\"NICER\":[\"6643010101\",\"6643010102\"],"
                        "\"NuSTAR\":[\"80902318002\"],\"STATE\":\"HSS\"}. "
                        "Omit any instrument not present in the epoch."
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
                    "type": "boolean",
                    "description": (
                        "Call AllData.clear() + AllModels.clear() before loading. Default true."
                    ),
                },
                "chatter": {
                    "type": "boolean",
                    "description": "Print file paths and loading summaries. Default false.",
                },
            },
            "required": ["epoch", "ixpe_path"],
        },
    },

    # ── Step 2: Cyg X-1 HSS tbfeo model ───────────────────────────────────
    "research_toolkit.xspec_utils.fit_tbfeo_hss": {
        "module": "xspec_utils.lib.xspec_help",
        "fn": "fit_tbfeo_plabs_polconst_diskbb_polconst_relxillCp",
        "description": (
            "Define and fit constant*tbfeo*plabs*(polconst*diskbb + polconst*relxillCp) "
            "for Cyg X-1 HSS with polarized absorption. "
            "Pass indices from load_data."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "indices": {
                    "type": "object",
                    "additionalProperties": True,
                    "description": "Return value of load_data — pass the object directly, do NOT JSON-encode it.",
                },
                "log_path": {"type": "string"},
                "audit_data": {"type": "boolean", "description": "Audit loaded data before fitting. Default true."},
            },
            "required": ["indices"],
        },
    },

    # ── Step 3: refine fit ─────────────────────────────────────────────────
    "research_toolkit.xspec_utils.run_fit": {
        "module": "xspec_utils.skills.session",
        "fn": "run_fit",
        "description": (
            "Run XSPEC Fit.perform() to refine convergence after fit_tbfeo_hss. "
            "Call multiple times until stat/dof stabilises. "
            "Returns dict with fit_stat, dof, stat_per_dof."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "n_iterations": {
                    "type": "integer",
                    "description": "XSPEC Fit.nIterations. Default 100.",
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
                    "type": "object",
                    "additionalProperties": True,
                    "description": "Return value of load_data — pass the object directly, do NOT JSON-encode it.",
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
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Residual types, e.g. [\"delchi\"]. Default ['delchi'].",
                },
                "xlim": {
                    "type": "array",
                    "items": {"type": "number"},
                    "description": "[xmin, xmax] in keV, e.g. [1, 79].",
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
                    "type": "object",
                    "additionalProperties": True,
                    "description": (
                        "Return value of load_data (must contain 'IXPE') — pass the object "
                        "directly, do NOT JSON-encode it."
                    ),
                },
                "e_bins": {
                    "type": "array",
                    "items": {"type": "number"},
                    "description": "Energy bin edges in keV. Default [2,3,4,5,6,8].",
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
                    "type": "array",
                    "items": {"type": "number"},
                    "description": "Bin edges in keV, e.g. [2.0,3.0,4.0,5.0,6.0,8.0]. This argument is named 'e_edges' — not 'bin_edges' or 'edges'.",
                },
            },
            "required": ["infile", "outfile", "e_edges"],
        },
    },
}
