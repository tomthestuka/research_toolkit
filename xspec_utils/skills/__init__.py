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
                "irf_date": {
                    "type": "string",
                    "description": (
                        "IRF calibration date 'YYYYMMDD', resolved to "
                        "'ixpe:obssim{irf_date}_alpha075:v13'. Default '20230702'."
                    ),
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
            "bkg_outer_rad) centred there, at irf_date. Every PHA1 product on disk after "
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
                "irf_date": {
                    "type": "string",
                    "description": (
                        "IRF calibration date 'YYYYMMDD', resolved to "
                        "'ixpe:obssim{irf_date}_alpha075:v13'. Default '20230702'. "
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
                        "'_structure_function' plus '.{fmt}' are appended per figure. "
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
                "fmt": {
                    "type": "string",
                    "enum": ["png", "pdf"],
                    "description": "Output image format for all four figures. Default 'png'.",
                },
            },
            "required": ["obs_db_path", "obsids", "catalog_csv"],
        },
    },
}
