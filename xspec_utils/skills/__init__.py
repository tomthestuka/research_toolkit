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
  7. ixpe_pa_pd_stability — constancy/variability stats on an existing
     pcube_timeseg CSV (no XSPEC/ixpeobssim required)
  8. ixpe_pa_pd_stability_plot — the four standard PA/PD stability diagnostic
     figures from an ixpe_pa_pd_stability catalog (no XSPEC/ixpeobssim required)
  9. ixpe_reduce_pha — reuse or reduce per-DU IXPE PHA1 (Stokes I) src+bkg
     spectra for one local ObsID, the input HR() and load_data need
  10. ixpe_HR — net Stokes-I hardness ratio HR=(H-S)/(H+S) survey across local
      IXPE ObsIDs, read directly from each ObsID's existing per-DU PHA1
      products (no XSPEC AllData loading -- pure PHA1/RMF FITS reads)
  11. bhxrb_pd_hr_survey — correlates PD (from ixpe_pcube_timeseg's integrated
      row) against HR (from ixpe_HR's CSV) across a survey of BH-XRB IXPE
      ObsIDs: detection classification, per-source/combined Spearman, a
      soft-vs-hard-state KS test, and a detection-bias Mann-Whitney U test
      (no XSPEC/ixpeobssim required)
  12. bhxrb_pd_hr_survey_plot — the five standard PD-hardness survey
      reviewer figures from a bhxrb_pd_hr_survey CSV (no XSPEC/ixpeobssim
      required)
  13. ixpe_pcube — energy-resolved IXPE PCUBE polarimetry for one local
      ObsID, or a native multi-ObsID combination (one xpselect/xpbin call
      across every listed ObsID's Level-2 files), background-subtracted,
      per caller-supplied energy bin edges
  14. ixpe_pcube_plot — PD(E)/PA(E) and Q-U plane companion figure for an
      ixpe_pcube CSV (no XSPEC/ixpeobssim required)
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

    # ── Time-resolved IXPE polarimetry (intra-observation stability) ───────
    "research_toolkit.xspec_utils.ixpe_pcube_timeseg": {
        "module": "xspec_utils.skills.session",
        "fn": "ixpe_pcube_timeseg",
        "description": (
            "Time-resolved 2-8 keV IXPE PCUBE polarimetry for one local ObsID, via "
            "ixpeobssim (xpselect -> xpbin PCUBE, weights=True), DU1/2/3 merged. "
            "Splits the observation into natural elapsed-time segments of "
            "segment_duration_s (default 5000 s; TSTART + k*segment_duration_s, NOT "
            "exposure-equalised) and additionally runs one integrated whole-observation "
            "reduction. Use this to test intra-observation polarization constancy "
            "(PD/PA/normalized-Stokes/MDP99 vs time) for a bright, source-dominated "
            "target already local at level >= 2 (check with a data-catalog tool first; "
            "this tool does not download data). Background subtraction is OFF by default "
            "(bkg=False) — appropriate for bright source-dominated targets (e.g. Cyg X-1, "
            "4U 1630-47); set bkg=True only when the source is faint enough that "
            "instrumental/diffuse background is a non-negligible fraction of the source "
            "aperture counts. When ra/dec are omitted, the source position is centroided "
            "automatically from one Level-2 event file. Before reducing a window, its "
            "GTI-clipped good exposure is checked against min_seg_exposure_s (default "
            "1000 s) and windows below it are skipped without running xpselect/xpbin -- "
            "this matters for observations with a long inter-pointing gap, where whole "
            "windows can fall inside the gap. Rows are appended to the CSV and flushed "
            "as each segment finishes, so a call that hits the sandbox time cap leaves a "
            "genuine checkpoint; call again with the same arguments (resume=True, the "
            "default) to compute only the remaining segments, or resume=False (or any "
            "changed argument) to force a full rebuild. Writes/updates a CSV of "
            "per-segment + integrated rows (segment_id, t_start_mjd, t_stop_mjd, "
            "exposure_s, PD_pct, PD_err, PA_deg, PA_err, QI, QI_err, UI, UI_err, MDP99, "
            "SIGNIF; the integrated row has segment_id=-1, plus a provenance comment "
            "header) to <obs_db_path>/<obsid>/pcube_timeseg_<obsid>.csv and records it "
            "in the mission catalog at level 3 with a 'complete' flag. This tool's "
            "return value NEVER carries the full per-segment table (it can exceed the "
            "sandbox transport limit for long observations) -- read csv_path for every "
            "row. Returns 'csv_path', 'n_segments' (rows resumed + newly computed), "
            "'n_skipped_empty', 'n_expected', 'complete' (bool: every expected segment "
            "plus the integrated row is present), 'resumed' (bool), 'integrated' (dict), "
            "'head'/'tail' (first/last 3 segment-row dicts), plus 'stdout'/'exit_code' "
            "(1 with 'stderr' on failure; nothing is recorded on failure)."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "obs_db_path": {
                    "type": "string",
                    "description": (
                        "The IXPE mission's obs-db directory (e.g. "
                        "'/data/fzhou/XSPEC/Cyg_X-1/IXPE'). Event files are expected at "
                        "<obs_db_path>/<obsid>/event_l2/."
                    ),
                },
                "obsid": {
                    "type": "string",
                    "description": "9-digit IXPE Observation ID, as a string, e.g. '01250401'.",
                },
                "ra": {
                    "type": "number",
                    "description": (
                        "Source RA in decimal degrees. Omit together with dec to "
                        "auto-centroid from one Level-2 event file."
                    ),
                },
                "dec": {
                    "type": "number",
                    "description": (
                        "Source Dec in decimal degrees. Omit together with ra to "
                        "auto-centroid from one Level-2 event file."
                    ),
                },
                "emin_kev": {
                    "type": "number",
                    "description": "Lower edge of the single energy bin, in keV. Default 2.0.",
                },
                "emax_kev": {
                    "type": "number",
                    "description": "Upper edge of the single energy bin, in keV. Default 8.0.",
                },
                "src_rad": {
                    "type": "number",
                    "description": "Source region radius in arcminutes. Default 1.0.",
                },
                "bkg": {
                    "type": "boolean",
                    "description": (
                        "Subtract a background annulus (bkg_inner_rad-bkg_outer_rad) from "
                        "each segment and the integrated cube. Default false (bright "
                        "source-dominated targets); set true for fainter sources."
                    ),
                },
                "bkg_inner_rad": {
                    "type": "number",
                    "description": "Background annulus inner radius in arcminutes (used only if bkg=true). Default 3.0.",
                },
                "bkg_outer_rad": {
                    "type": "number",
                    "description": "Background annulus outer radius in arcminutes (used only if bkg=true). Default 5.0.",
                },
                "segment_duration_s": {
                    "type": "number",
                    "description": (
                        "Elapsed-time segment length in seconds (natural TSTART + "
                        "k*segment_duration_s edges, not exposure-equalised). Default 5000.0. "
                        "Increase for a coarser, faster time series; decrease for finer time "
                        "resolution at the cost of larger per-segment MDP99/uncertainties."
                    ),
                },
                "min_seg_exposure_s": {
                    "type": "number",
                    "description": (
                        "Minimum GTI-clipped good exposure (seconds) a window must have "
                        "before it is reduced at all; windows below this are skipped "
                        "before any xpselect/xpbin call and counted in n_skipped_empty. "
                        "Default 1000.0. Lower it only if a target observation's real "
                        "GTI structure has genuine short segments worth keeping; raise it "
                        "to be more conservative near pointing-gap/GTI boundaries."
                    ),
                },
                "resume": {
                    "type": "boolean",
                    "description": (
                        "Resume from an existing pcube_timeseg_<obsid>.csv whose "
                        "provenance (all reduction parameters plus the derived segment "
                        "edges) exactly matches this call, computing only the segments "
                        "still missing. Default true. Any mismatch (including a prior "
                        "call with different arguments) deletes the stale CSV and "
                        "rebuilds from scratch rather than mixing rows from different "
                        "configurations; set false to always force a full rebuild."
                    ),
                },
            },
            "required": ["obs_db_path", "obsid"],
        },
    },

    # ── Per-DU IXPE PHA1 reduction (reuse-or-reduce) ───────────────────────
    "research_toolkit.xspec_utils.ixpe_reduce_pha": {
        "module": "xspec_utils.skills.session",
        "fn": "ixpe_reduce_pha",
        "description": (
            "Reuse or reduce per-DU IXPE PHA1 (Stokes I) source+background spectra "
            "for one local ObsID -- the input HR() (hardness ratio) and load_data need, "
            "surfaced as its own tool rather than folded into an HR tool since HR reads "
            "PHA. Locates Level-2 DU1/DU2/DU3 event files under "
            "<obs_db_path>/<obsid>/event_l2/ (preferring '_rej' instrumental-background- "
            "rejected files, falling back to plain evt2_v*.fits); this tool does not "
            "download data, so check ObsID locality with a data-catalog tool first. "
            "If the per-DU src+bkg PHA1 files already exist on disk, they are reused "
            "as-is -- no re-reduction, no region files written, 'reused'=true, 'ra'/'dec' "
            "null. Otherwise the source RA/Dec is centroided from the DU1 event file via "
            "find_src.run_extraction (writing DS9 'src.reg'/'bkg.reg' under "
            "<obs_db_path>/<obsid>/ -- NOT read from FITS headers, NOT from a fixed "
            "source-coordinate table), and xpselect+xpbin PHA1 is run per DU with a "
            "circular source region (src_rad) and annular background (bkg_inner_rad- "
            "bkg_outer_rad) centred there, at the IRF epoch auto-selected as the installed "
            "epoch nearest the observation date. Every PHA1 product on disk after "
            "the call (3 DUs x src/bkg) is recorded in the mission catalog at level 3 via "
            "obsdb_catalog.neutral.advance, whether reused or freshly reduced -- so a "
            "downstream HR/load_data call always has real, catalogued products to read. "
            "Returns 'src_pha'/'bkg_pha' (each a 3-element DU1/DU2/DU3 list of Stokes-I "
            "PHA1 paths), 'ra'/'dec' (null when reused), 'reused' (bool), 'src_reg'/ "
            "'bkg_reg' (null when reused), 'stdout', 'exit_code' (1 with 'stderr' on "
            "failure -- e.g. no local Level-2 files, or fewer/more than one file per DU -- "
            "nothing is recorded on failure)."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "obs_db_path": {
                    "type": "string",
                    "description": (
                        "The IXPE mission's obs-db directory (e.g. "
                        "'/data/fzhou/XSPEC/Cyg_X-1/IXPE'). Event files are expected at "
                        "<obs_db_path>/<obsid>/event_l2/; src.reg/bkg.reg and the PHA1 "
                        "products (when freshly reduced) land under "
                        "<obs_db_path>/<obsid>/."
                    ),
                },
                "obsid": {
                    "type": "string",
                    "description": "9-digit IXPE Observation ID, as a string, e.g. '01002901'.",
                },
                "src_rad": {
                    "type": "number",
                    "description": (
                        "Source region radius in arcminutes. Default 1.0. Used only "
                        "when a fresh reduction is required (ignored when PHA1 products "
                        "are reused)."
                    ),
                },
                "bkg_inner_rad": {
                    "type": "number",
                    "description": (
                        "Background annulus inner radius in arcminutes. Default 3.0. "
                        "Used only when a fresh reduction is required."
                    ),
                },
                "bkg_outer_rad": {
                    "type": "number",
                    "description": (
                        "Background annulus outer radius in arcminutes. Default 5.0. "
                        "Used only when a fresh reduction is required."
                    ),
                },
            },
            "required": ["obs_db_path", "obsid"],
        },
    },

    # ── Constancy/variability statistics on an existing pcube_timeseg CSV ──
    "research_toolkit.xspec_utils.ixpe_pa_pd_stability": {
        "module": "xspec_utils.lib.pa_pd_stability",
        "fn": "pa_pd_stability",
        "description": (
            "Intra-observation PA/PD/Stokes constancy analysis across one or more local "
            "IXPE ObsIDs, from each ObsID's existing pcube_timeseg_<obsid>.csv (produced "
            "by ixpe_pcube_timeseg -- run that first; this tool never reduces data itself, "
            "so it needs no XSPEC or ixpeobssim and runs fast). For each ObsID: inverse-"
            "variance-weighted chi-squared constancy tests on PD, PA (mod-180 circular "
            "handling -- PA residuals wrapped to (-90,+90] so e.g. a +89/-89 deg flip "
            "scores as ~2 deg, not ~178), QI, and UI; the AUTHORITATIVE variability call "
            "is the 2D Stokes test chi2_2D=chi2_QI+chi2_UI, dof_2D=2*(N-1), flagging "
            "'variable' when p_2D < variability_p_threshold (PD/PA chi-squared are "
            "secondary/descriptive, not the variability call); moment-estimator intrinsic "
            "scatter (sigma_int=sqrt(max(0, sample_variance - mean_measurement_variance)), "
            "upper-limit flagged when non-positive) for PD and PA; peak-to-peak PD/PA "
            "ranges; first-order structure functions SF_PD/SF_PA at sf_lags_ks (NaN for a "
            "lag with no segment pair inside its tolerance window, and only computed at "
            "all when an ObsID has >= min_seg_for_sf segments); and, across all requested "
            "ObsIDs together, a Spearman rank correlation of each ObsID's mean MDP99 "
            "against its reduced chi-squared for PD and for PA (NaN with a note when fewer "
            "than 3 ObsIDs have finite values). An ObsID with no pcube_timeseg CSV (or an "
            "unrecognised one -- missing the provenance header ixpe_pcube_timeseg writes) "
            "is recorded as an error entry and skipped, never raised. Writes "
            "pa_pd_stability_catalog.csv (one row per ObsID) and "
            "pa_pd_variability_summary.txt (per-ObsID verdict, secondary stats, the two "
            "Spearman results, and the source CSV paths) under output_dir. Returns "
            "'per_obs' (obsid -> stats dict, including 'sf_pd'/'sf_pa'), 'spearman', "
            "'catalog_csv', 'summary_txt', 'exit_code' (always 0 -- per-ObsID problems are "
            "recorded, not raised)."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "obs_db_path": {
                    "type": "string",
                    "description": (
                        "The IXPE mission's obs-db directory (e.g. "
                        "'/data/fzhou/XSPEC/Cyg_X-1/IXPE'). Each ObsID's CSV is read from "
                        "<obs_db_path>/<obsid>/pcube_timeseg_<obsid>.csv."
                    ),
                },
                "obsids": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "ObsIDs to analyse (each must already carry a pcube_timeseg CSV).",
                },
                "output_dir": {
                    "type": "string",
                    "description": "Directory the catalog CSV and summary text file are written to (created if absent).",
                },
                "sources": {
                    "type": "object",
                    "additionalProperties": {"type": "string"},
                    "description": (
                        "Optional {obsid: source_name} label map for the catalog/summary "
                        "(e.g. distinguishing which BH-XRB each ObsID belongs to). Falls "
                        "back to the obsid itself when omitted or an obsid is absent from it."
                    ),
                },
                "variability_p_threshold": {
                    "type": "number",
                    "description": (
                        "p_2D below this flags an ObsID 'variable' in the authoritative "
                        "2D Stokes (QI+UI) test. Default 0.01."
                    ),
                },
                "sf_lags_ks": {
                    "type": "array",
                    "items": {"type": "number"},
                    "description": (
                        "Structure-function lags to evaluate, in kiloseconds. Default "
                        "[5, 10, 20]. The tolerance half-width around each lag is half the "
                        "ObsID's own segment_duration_s (from the CSV's provenance header)."
                    ),
                },
                "min_seg_for_sf": {
                    "type": "integer",
                    "description": (
                        "Minimum per-segment row count an ObsID must have before its "
                        "structure functions are computed at all (else every lag is NaN "
                        "for that ObsID). Default 10."
                    ),
                },
            },
            "required": ["obs_db_path", "obsids", "output_dir"],
        },
    },

    # ── PA/PD stability diagnostic figures (companion plotter) ─────────────
    "research_toolkit.xspec_utils.ixpe_pa_pd_stability_plot": {
        "module": "xspec_utils.lib.pa_pd_stability",
        "fn": "stability_plot",
        "description": (
            "Renders the four standard IXPE PA/PD stability diagnostic figures from the "
            "products of ixpe_pa_pd_stability -- run that tool first to get catalog_csv "
            "(and ixpe_pcube_timeseg before it, to get each ObsID's pcube_timeseg CSV). "
            "Pure matplotlib (Agg backend); imports no XSPEC/ixpeobssim, so it runs fast "
            "and needs no science stack. Writes exactly four files, named by appending "
            "'_pa_pd_timeseries', '_stokes_plane', '_chi2_hist', '_structure_function' "
            "(plus '.png'/'.pdf') to save_path_prefix: (1) per-ObsID PA (top) / PD (bottom) "
            "panels vs time elapsed in kiloseconds from the earliest segment's t_start, "
            "with the catalog's weighted-mean +/-1sigma band overlaid when the ObsID is in "
            "catalog_csv; (2) per-ObsID Q/I vs U/I Stokes-plane scatter color-coded by "
            "elapsed time, with 1sigma error crosses and the origin marked; (3) histograms "
            "of rchi2_PD/rchi2_PA across EVERY row of catalog_csv (not just obsids), with a "
            "reference line at 1.0; (4) PA/PD structure functions vs lag, one line per "
            "ObsID with >= sf_min_segments segments (an ObsID short of that is omitted from "
            "this panel only, still appears in the other three; if none qualify the panel "
            "carries an explicit annotation instead of empty axes). An ObsID with no "
            "readable pcube_timeseg CSV is skipped and noted in stdout; if every requested "
            "ObsID is unreadable, nothing is written and exit_code=1 with stderr set. "
            "catalog_csv is a REQUIRED upstream product: a missing or unreadable (e.g. "
            "wrong-schema) catalog_csv is a hard error -- exit_code=1 with stderr set and "
            "no files written -- it never silently degrades to an empty catalog. This is "
            "distinct from an ObsID simply absent from an otherwise-present, readable "
            "catalog_csv, which still gets its time-series/Stokes-plane panels -- only its "
            "weighted-mean overlay is skipped. Returns 'figures' (4 FigureRef dicts), "
            "'structure_function' (obsid -> the lags_ks/sf_pa/sf_pd values actually "
            "plotted, for verification), 'stdout', 'exit_code'."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "obs_db_path": {
                    "type": "string",
                    "description": (
                        "The IXPE mission's obs-db directory (e.g. "
                        "'/data/fzhou/XSPEC/Cyg_X-1/IXPE'). Each ObsID's per-segment CSV is "
                        "read from <obs_db_path>/<obsid>/pcube_timeseg_<obsid>.csv."
                    ),
                },
                "obsids": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": (
                        "ObsIDs to draw the time-series/Stokes-plane/structure-function "
                        "panels for. The chi2 histogram panel instead covers every row "
                        "already in catalog_csv, regardless of this list."
                    ),
                },
                "catalog_csv": {
                    "type": "string",
                    "description": (
                        "Path to the pa_pd_stability_catalog.csv written by "
                        "ixpe_pa_pd_stability."
                    ),
                },
                "save_path_prefix": {
                    "type": "string",
                    "description": (
                        "Base path; '_pa_pd_timeseries', '_stokes_plane', '_chi2_hist', "
                        "'_structure_function' plus '.png' are appended per figure. "
                        "Auto-generated timestamped prefix in the current directory if "
                        "omitted."
                    ),
                },
                "sf_lags_ks": {
                    "type": "array",
                    "items": {"type": "number"},
                    "description": "Structure-function lags to plot, in kiloseconds. Default [5, 10, 20].",
                },
                "sf_min_segments": {
                    "type": "integer",
                    "description": (
                        "Minimum per-segment row count an ObsID needs to appear in the "
                        "structure-function panel. Default 10."
                    ),
                },
            },
            "required": ["obs_db_path", "obsids", "catalog_csv"],
        },
    },

    # ── Energy-resolved / multi-ObsID-combined IXPE polarimetry ────────────
    "research_toolkit.xspec_utils.ixpe_pcube": {
        "module": "xspec_utils.skills.session",
        "fn": "ixpe_pcube",
        "description": (
            "Energy-resolved IXPE PCUBE polarimetry for one local ObsID, OR a native "
            "multi-ObsID combination, background-subtracted -- e.g. reproducing "
            "Krawczynski+2022 Cyg X-1 energy-resolved polarimetry (Fig. S5/Table S2) or "
            "its native multi-ObsID (e.g. May 01002901 + June 01250101) combined cube. "
            "Locates the mission's full Level-2 event file list via "
            "red_ixpe.get_all_file_list(obs_db_path); this tool does not download data, so "
            "check ObsID locality with a data-catalog tool first. "
            "'obs' as a single 9-digit ObsID string reduces just that observation. 'obs' as "
            "a LIST of ObsID strings filters and concatenates every listed ObsID's Level-2 "
            "files and passes them into ONE xpselect/xpbin PCUBE call, so the returned cube "
            "sums raw Stokes I/Q/U/COUNTS across every DU *and* every listed ObsID before "
            "PD/PA are derived -- the native directly-combined reduction, NOT an inverse-"
            "variance approximation built by merging separately-computed per-ObsID PD/PA. "
            "'ebinning' is the caller's energy bin edges in keV (ebinalg='LIST' unchanged), "
            "default [2,3,4,6,8] (4 bands); pass e.g. 11 equal-width edges over 2-8 keV for "
            "a finer energy grid; a 2-edge list yields a single broadband row. 'ra'/'dec' "
            "are REQUIRED (no auto-centroid here -- unlike ixpe_pcube_timeseg -- since a "
            "multi-ObsID combination especially needs one explicit, shared region centre). "
            "The IRF calibration epoch is auto-selected as the installed epoch nearest the "
            "midpoint across every requested ObsID -- no epoch argument is accepted here. "
            "Reads PD, PD_err, PA, PA_err, QN, QN_err, UN, UN_err, MDP99, SIGNIF, E_LO, E_HI "
            "straight off the background-subtracted source cube's polarization()/attribute "
            "accessors -- never hand-rolled Stokes arithmetic. PD, PD_err, QN, QN_err, UN, "
            "UN_err, MDP99 are DIMENSIONLESS FRACTIONS (not percent, unlike "
            "ixpe_pcube_timeseg's PD_pct/MDP99) -- the same units the cube itself stores "
            "them in; PA, PA_err are degrees. Raises (returned as exit_code=1, nothing "
            "written or recorded) if a requested ObsID has no local Level-2 event files, is "
            "missing DU1/DU2/DU3 (or carries more than one file for a DU), or ebinning has "
            "fewer than 2 edges. Writes a provenance-headed CSV, one row per energy band -- "
            "there is no separate 'integrated' sentinel row (unlike ixpe_pcube_timeseg's "
            "segment_id=-1); a single-edge-pair ebinning naturally yields that one broadband "
            "row instead. A single ObsID's CSV lands at "
            "<obs_db_path>/<obsid>/ixpe_pcube_<obsid>.csv; a combined request's CSV lands at "
            "<obs_db_path>/ixpe_pcube_<obsid1>_<obsid2>_....csv (mission root -- no single "
            "ObsID directory owns a combined product, mirroring bhxrb_hardness_ratios.csv's "
            "placement). Records one kind='pcube' product per contributing ObsID via "
            "obsdb_catalog.neutral.advance (level 3), all pointing at the same CSV for a "
            "combined request. Returns 'bands' (list of per-band dicts), 'provenance' "
            "(obsids, ebinning, src_rad, bkg_inner_rad, bkg_outer_rad, ra, dec, irf_date, "
            "combined:bool), 'csv_path', 'stdout', 'exit_code'. This is a data tool -- it "
            "writes no plot; pair it with a separate plotter tool for figures."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "obs_db_path": {
                    "type": "string",
                    "description": (
                        "The IXPE mission's obs-db directory (e.g. "
                        "'/data/fzhou/XSPEC/Cyg_X-1/IXPE'). Event files are expected at "
                        "<obs_db_path>/<obsid>/event_l2/ for every requested ObsID."
                    ),
                },
                "obs": {
                    "oneOf": [
                        {"type": "string"},
                        {"type": "array", "items": {"type": "string"}, "minItems": 1},
                    ],
                    "description": (
                        "A single 9-digit IXPE ObsID string, OR a list of ObsID strings to "
                        "natively combine into one xpselect/xpbin PCUBE call (e.g. "
                        "[\"01002901\", \"01250101\"]). A list of one ObsID behaves like "
                        "passing that ObsID as a plain string."
                    ),
                },
                "ra": {
                    "type": "number",
                    "description": (
                        "Source RA in decimal degrees, shared across every listed ObsID. "
                        "Required -- this tool never auto-centroids."
                    ),
                },
                "dec": {
                    "type": "number",
                    "description": (
                        "Source Dec in decimal degrees, shared across every listed ObsID. "
                        "Required -- this tool never auto-centroids."
                    ),
                },
                "ebinning": {
                    "type": "array",
                    "items": {"type": "number"},
                    "minItems": 2,
                    "description": (
                        "Energy bin edges in keV, e.g. [2,3,4,6,8]. Default [2,3,4,6,8] "
                        "(4 bands, the Krawczynski+2022 Table S2 grid). Must carry >= 2 "
                        "edges; a 2-edge list (e.g. [2,8]) yields a single broadband row."
                    ),
                },
                "src_rad": {
                    "type": "number",
                    "description": "Source region radius in arcminutes. Default 1.0.",
                },
                "bkg_inner_rad": {
                    "type": "number",
                    "description": "Background annulus inner radius in arcminutes. Default 3.0.",
                },
                "bkg_outer_rad": {
                    "type": "number",
                    "description": "Background annulus outer radius in arcminutes. Default 5.0.",
                },
            },
            "required": ["obs_db_path", "obs", "ra", "dec"],
        },
    },

    # ── Net Stokes-I hardness ratio survey (fix + surface) ─────────────────
    "research_toolkit.xspec_utils.ixpe_HR": {
        "module": "xspec_utils.skills.session",
        "fn": "ixpe_HR",
        "description": (
            "Net Stokes-I hardness ratio HR=(H-S)/(H+S) for a survey of local IXPE "
            "ObsIDs, read directly from each ObsID's existing per-DU PHA1 (Stokes I) "
            "src(+bkg) products -- the SAME PHA1 geometry ixpe_reduce_pha produces and "
            "that feeds ixpe_pcube_timeseg's PD, so PD and HR share identical "
            "extraction regions. This never re-reduces PHA -- run ixpe_reduce_pha "
            "first for any ObsID missing its PHA1 products. "
            "Fixes three defects in the prior red_ixpe.HR: (1) DU3 was silently "
            "dropped from the summed rate by an off-by-one loop bound over the loaded "
            "XSPEC spectra -- this version attempts DU1/DU2/DU3 independently via "
            "pha()'s existing file locators, sums net rate over every DU actually "
            "present (never silently dropping a present DU -- a present DU whose PHA1 "
            "is missing on disk is reported as this ObsID's error instead), and drops "
            "only a genuinely-missing DU (IXPE can fly with a DU off) without raising; "
            "(2) HR was H/S -- now the community-normalized (H-S)/(H+S), bounded "
            "[-1,+1]; (3) the returned error came from misreading XSPEC's "
            "Spectrum.rate tuple as (rate, variance) when it is actually (total, "
            "background, net, net_error) -- now a closed-form first-order Gaussian "
            "propagation sigma_HR=(2/(H+S)^2)*sqrt(S^2*sigma_H^2+H^2*sigma_S^2), with "
            "each band's net rate and error read straight from the PHA RATE/STAT_ERR "
            "columns (summed over in-band channels and over contributing DUs) -- no "
            "XSPEC AllData loading at all. soft_band/hard_band channels are resolved "
            "per PHA file from its own RESPFILE RMF EBOUNDS -- never a hardcoded "
            "channel grid -- so any band works, not just the survey's default 2-4/"
            "4-8 keV; when a PHA1's RESPFILE predates a reduction-environment "
            "migration and no longer resolves on disk (the real state of most "
            "already-reduced Cyg X-1 IXPE PHA1 products), the same DU's installed "
            "CALDB RMF is used instead -- the identical EBOUNDS grid, since it is a "
            "fixed instrumental characteristic shared across every IRF date/DU. "
            "bkg=True subtracts each DU's background PHA1 scaled by "
            "(src_BACKSCAL/bkg_BACKSCAL) before summing; bkg=False (default, right "
            "for source-dominated bright BH-XRBs) uses the source PHA1 alone and "
            "never reads the background file. An ObsID with no local Level-2 event "
            "files, or that red_ixpe.HR raises on, is recorded in the returned "
            "'errors' dict and excluded from the CSV -- never a dummy row, "
            "never a crash of the rest of the survey. Writes bhxrb_hardness_ratios.csv "
            "(columns obsid, soft_rate, soft_err, hard_rate, hard_err, HR, HR_err) at "
            "<obs_db_path>/bhxrb_hardness_ratios.csv (overwriting any previous run), "
            "and records one 'hardness_ratio' product per successful ObsID in the "
            "mission catalog at level 3 via obsdb_catalog.neutral.advance, with "
            "'params' carrying the real soft_band/hard_band/bkg/dus_used used for "
            "that ObsID. Returns 'csv_path', 'rows' (the successful per-obsid rows), "
            "'errors' (obsid -> message), 'stdout', 'exit_code' (1 with 'stderr' only "
            "if EVERY ObsID failed -- nothing written or recorded in that case; 0 "
            "otherwise, with per-ObsID failures surfaced in 'errors')."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "obs_db_path": {
                    "type": "string",
                    "description": (
                        "The IXPE mission's obs-db directory (e.g. "
                        "'/data/fzhou/XSPEC/Cyg_X-1/IXPE'). Event files are expected "
                        "at <obs_db_path>/<obsid>/event_l2/; bhxrb_hardness_ratios.csv "
                        "is written at <obs_db_path>/."
                    ),
                },
                "obsids": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "9-digit IXPE Observation IDs to include in the survey.",
                },
                "soft_band": {
                    "type": "array",
                    "items": {"type": "number"},
                    "description": "[lo, hi] soft-band energy edges in keV. Default [2.0, 4.0].",
                },
                "hard_band": {
                    "type": "array",
                    "items": {"type": "number"},
                    "description": "[lo, hi] hard-band energy edges in keV. Default [4.0, 8.0].",
                },
                "bkg": {
                    "type": "boolean",
                    "description": (
                        "Subtract each DU's background PHA1 (scaled by "
                        "src_BACKSCAL/bkg_BACKSCAL) before summing band rates. Default "
                        "false, appropriate for source-dominated bright BH-XRBs; set "
                        "true only where background is a non-negligible aperture "
                        "fraction."
                    ),
                },
            },
            "required": ["obs_db_path", "obsids"],
        },
    },

    # ── PD-vs-HR correlation/statistics survey (data only; no figures) ─────
    "research_toolkit.xspec_utils.bhxrb_pd_hr_survey": {
        "module": "xspec_utils.lib.bhxrb_pd_hr_survey",
        "fn": "pd_hr_survey",
        "description": (
            "Correlates IXPE polarization degree (PD) against hardness ratio (HR) "
            "across a survey of BH-XRB ObsIDs. For each obsid in obsids_with_sources, "
            "reads PD_pct/PD_err/MDP99 from that obsid's existing "
            "<pcube_csv_dir>/<obsid>/pcube_timeseg_<obsid>.csv integrated (segment_id=-1) "
            "row -- run ixpe_pcube_timeseg first -- and HR/HR_err from hr_csv (the "
            "bhxrb_hardness_ratios.csv ixpe_HR writes -- run that first). Classifies "
            "detected=PD_pct>MDP99 and flags mdp_flagged=MDP99>mdp_flag_pct (a data-"
            "quality note, not an exclusion). An obsid with no readable pcube CSV, or no "
            "row in hr_csv, gets a recorded note and is excluded from every statistic -- "
            "never raised, never a dummy value. Over the detected-only rows (a non-"
            "detection's PD is MDP99-noise-dominated, not a physical measurement): "
            "per-source Spearman rank correlation (PD_pct vs HR) for every source with "
            ">= min_n_intra detections (a source short of that is omitted from "
            "intra_source_correlations.csv, not crashed on); the combined Spearman "
            "across every detected point regardless of source; a two-sample "
            "Kolmogorov-Smirnov test between the detected PD_pct of soft_state_sources "
            "vs hard_state_sources (both passed by the caller -- never hardcoded -- so "
            "this serves any source-group split). Over every obsid with both a readable "
            "PD and HR: a two-sided Mann-Whitney U test on HR between the detected and "
            "non-detected groups (does HR predict whether a detection was made). Every "
            "statistic is NaN, never an exception, when its input group is empty or -- "
            "for the two Spearman correlations -- has fewer than 3 finite pairs. Writes "
            "bhxrb_pd_hardness_survey.csv, intra_source_correlations.csv, "
            "combined_correlation.csv, ks_test_soft_vs_hard.csv, "
            "detection_bias_assessment.csv under obs_db_path, and records one "
            "'pd_hardness_survey' product per obsid with both PD and HR in the mission "
            "catalog at level 3 via obsdb_catalog.neutral.advance. Returns 'rows' (one "
            "per requested obsid), 'skipped' (obsids with a note), "
            "'per_source_correlations', 'spearman_combined', 'ks_result', 'mwu_result', "
            "the five output CSV paths, 'exit_code' (always 0 -- a per-obsid problem is "
            "recorded, not raised)."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "obs_db_path": {
                    "type": "string",
                    "description": (
                        "Directory the five output CSVs are written to and the mission "
                        "catalog (catalog.json) is recorded in, e.g. "
                        "'/data/fzhou/XSPEC/Cyg_X-1/IXPE'."
                    ),
                },
                "obsids_with_sources": {
                    "type": "object",
                    "additionalProperties": {"type": "string"},
                    "description": (
                        "Map of obsid -> source-name label (pass the object directly, "
                        "do NOT JSON-encode it), e.g. {\"01002901\": \"Cyg X-1\", "
                        "\"01250401\": \"4U 1630-47\"}. Every obsid appears in "
                        "bhxrb_pd_hardness_survey.csv (with a note if skipped); the "
                        "label groups obsids for the intra-source correlations and the "
                        "soft/hard KS split."
                    ),
                },
                "hr_csv": {
                    "type": "string",
                    "description": (
                        "Path to bhxrb_hardness_ratios.csv written by ixpe_HR. A missing "
                        "file is a hard error; a missing per-obsid row is a per-obsid note."
                    ),
                },
                "pcube_csv_dir": {
                    "type": "string",
                    "description": (
                        "Directory holding each obsid's "
                        "<obsid>/pcube_timeseg_<obsid>.csv, written by ixpe_pcube_timeseg "
                        "(often the same directory as obs_db_path)."
                    ),
                },
                "mdp_flag_pct": {
                    "type": "number",
                    "description": (
                        "MDP99 (percent) above this flags the obsid mdp_flagged in the "
                        "survey CSV -- a data-quality note, not an exclusion. Default 15.0."
                    ),
                },
                "min_n_intra": {
                    "type": "integer",
                    "description": (
                        "Minimum detected-obsid count a source needs before its "
                        "intra-source Spearman correlation is computed at all. Default 4."
                    ),
                },
                "soft_state_sources": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": (
                        "Source-name labels (matched against obsids_with_sources' "
                        "values) forming the KS test's soft-state group. Default empty "
                        "-- the KS test is NaN until both groups are supplied."
                    ),
                },
                "hard_state_sources": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": (
                        "Source-name labels forming the KS test's hard/intermediate-"
                        "state group. Default empty."
                    ),
                },
            },
            "required": ["obs_db_path", "obsids_with_sources", "hr_csv", "pcube_csv_dir"],
        },
    },

    # ── PD-vs-HR survey diagnostic figures (companion plotter) ─────────────
    "research_toolkit.xspec_utils.bhxrb_pd_hr_survey_plot": {
        "module": "xspec_utils.lib.bhxrb_pd_hr_survey",
        "fn": "pd_hr_survey_plot",
        "description": (
            "Renders the five standard PD-hardness survey reviewer figures from "
            "bhxrb_pd_hardness_survey.csv -- run bhxrb_pd_hr_survey first to get "
            "survey_csv. Pure matplotlib (Agg backend); imports no XSPEC/ixpeobssim, "
            "so it runs fast and needs no science stack. Writes exactly five files, "
            "named by appending 'pd_vs_hr_combined', 'pd_vs_hr_per_source', "
            "'pd_distribution_comparison', 'hr_detection_bias', "
            "'feasibility_mdp99_exposure' (plus '.png') to save_path_prefix: "
            "(1) PD_pct vs HR scatter over every detected row, color-coded by source, "
            "annotated with the combined Spearman rho/p; (2) the same scatter split "
            "one panel per source (every source with a usable row), each annotated "
            "with its own Spearman rho/p when it has >= 3 detections; (3) detected "
            "PD_pct histograms for soft_state_sources vs hard_state_sources (both "
            "caller-supplied, never hardcoded -- the same grouping "
            "bhxrb_pd_hr_survey's KS test uses), annotated with the KS statistic/p; "
            "(4) HR histograms for the detected vs non-detected groups (every usable "
            "row), annotated with the Mann-Whitney U/p; (5) MDP99 (percent) vs "
            "exposure (kiloseconds) scatter, color-coded by detection, with a "
            "mdp_flag_pct reference line -- exposure is read per obsid from its own "
            "<obs_db_path>/<obsid>/pcube_timeseg_<obsid>.csv integrated row (not "
            "carried by survey_csv itself), so an obsid whose pcube_timeseg CSV "
            "isn't locally readable is skipped from this panel only and noted in "
            "stdout. An empty group (no detections, no soft/hard match, no readable "
            "exposure data, ...) never degrades to a misleading empty plot -- the "
            "affected panel carries an explicit annotation instead. survey_csv is a "
            "REQUIRED upstream product: missing, wrong-schema, or entirely-skipped "
            "(every row carries a note, i.e. every requested obsid failed upstream) "
            "is a hard error -- exit_code=1 with stderr set and no files written. "
            "Returns 'figures' (5 FigureRef dicts), 'plotted' (the exact per-panel "
            "values drawn -- 'combined_points', 'per_source_points', "
            "'distribution_groups', 'detection_bias_hr', 'mdp99_exposure' -- for "
            "verification), 'stdout', 'exit_code'."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "obs_db_path": {
                    "type": "string",
                    "description": (
                        "The IXPE mission's obs-db directory (e.g. "
                        "'/data/fzhou/XSPEC/Cyg_X-1/IXPE'). Only used for the "
                        "feasibility panel: each usable obsid's exposure is read from "
                        "<obs_db_path>/<obsid>/pcube_timeseg_<obsid>.csv."
                    ),
                },
                "survey_csv": {
                    "type": "string",
                    "description": "Path to bhxrb_pd_hardness_survey.csv written by bhxrb_pd_hr_survey.",
                },
                "save_path_prefix": {
                    "type": "string",
                    "description": (
                        "Base path; 'pd_vs_hr_combined', 'pd_vs_hr_per_source', "
                        "'pd_distribution_comparison', 'hr_detection_bias', "
                        "'feasibility_mdp99_exposure' plus '.png' are appended per "
                        "figure. Auto-generated timestamped prefix in the current "
                        "directory if omitted."
                    ),
                },
                "soft_state_sources": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": (
                        "Source-name labels (matched against survey_csv's 'source' "
                        "column) forming the pd_distribution_comparison panel's "
                        "soft-state group -- the same split passed as "
                        "bhxrb_pd_hr_survey's own soft_state_sources. Default empty "
                        "-- the panel is annotated 'no groups supplied' until both "
                        "this and hard_state_sources are given."
                    ),
                },
                "hard_state_sources": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Source-name labels forming the panel's hard/intermediate-state group. Default empty.",
                },
                "mdp_flag_pct": {
                    "type": "number",
                    "description": (
                        "MDP99 (percent) reference line drawn on "
                        "feasibility_mdp99_exposure -- the same data-quality "
                        "threshold bhxrb_pd_hr_survey flags rows with. Default 15.0."
                    ),
                },
            },
            "required": ["obs_db_path", "survey_csv"],
        },
    },

    # ── Energy-resolved IXPE polarimetry figures (companion plotter) ───────
    "research_toolkit.xspec_utils.ixpe_pcube_plot": {
        "module": "xspec_utils.lib.ixpe_pcube_plot",
        "fn": "ixpe_pcube_plot",
        "description": (
            "Renders the PD(E)/PA(E) and Q-U-plane reviewer figure for an ixpe_pcube "
            "CSV -- run ixpe_pcube first to get csv_path (Krawczynski+2022 Fig. S5 "
            "style). Pure matplotlib (Agg backend); imports no XSPEC/ixpeobssim, so it "
            "runs fast and needs no science stack. Writes exactly ONE PNG with three "
            "panels: PD vs band-centre energy (top left), PA vs band-centre energy "
            "(bottom left, sharing the energy axis), and QN vs UN -- the normalized "
            "Stokes/Q-U plane (right, color-coded by band-centre energy) -- each with "
            "1-sigma error bars read straight from csv_path's own columns (PD, PD_err, "
            "PA, PA_err, QN, QN_err, UN, UN_err are plotted exactly as ixpe_pcube wrote "
            "them -- dimensionless fractions for PD/QN/UN, degrees for PA -- never "
            "recomputed). A single-band CSV still renders (one marker; no line to "
            "connect). A band with any non-finite value in the plotted columns is "
            "skipped -- never drawn as zero. Raises (no PNG written) if csv_path does "
            "not exist, is missing the ixpe_pcube provenance header, carries zero band "
            "rows, or every band is skipped for non-finite values. Returns 'figures' "
            "(one FigureRef dict), 'plotted' (the exact E_LO/E_HI/E_ctr/PD/PD_err/PA/"
            "PA_err/QN/QN_err/UN/UN_err arrays drawn, plus n_bands/n_skipped -- for "
            "differential verification against csv_path), 'stdout', 'exit_code' "
            "(always 0 -- failures raise instead of degrading to an empty plot)."
        ),
        "args_schema": {
            "type": "object",
            "properties": {
                "csv_path": {
                    "type": "string",
                    "description": (
                        "Path to the ixpe_pcube_<obsid(s)>.csv written by ixpe_pcube "
                        "(one row per energy band, with the '# ixpe_pcube provenance: "
                        "...' header line)."
                    ),
                },
                "save_path": {
                    "type": "string",
                    "description": "Full output PNG path. Auto-generated timestamped file if omitted.",
                },
                "title": {
                    "type": "string",
                    "description": (
                        "Figure suptitle. Defaults to the CSV provenance's obsids "
                        "joined by '+'."
                    ),
                },
            },
            "required": ["csv_path"],
        },
    },
}
