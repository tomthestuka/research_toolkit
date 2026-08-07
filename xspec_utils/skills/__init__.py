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
      ixpe_pcube CSV, with an optional fit_tbabs_pollin_powerlaw PD(E)/PA(E)
      overlay (falling back to an internal scipy weighted-linear fit of the
      PCUBE points when no pollin quartet is supplied) (no XSPEC/ixpeobssim
      required)
  15. fit_tbabs_polconst_powerlaw — polconst spectropolarimetric cross-check
      of the standard PCUBE polarization: fits tbabs*polconst*powerlaw to
      every IXPE Stokes I/Q/U triplet in a load_data indices dict, tying
      polconst's PD/PA into one free pair for the whole session
  16. fit_tbabs_pollin_powerlaw — the paper's actual linear energy-dependence
      route for Cyg X-1 polarization: fits tbabs*pollin*powerlaw to every
      IXPE Stokes I/Q/U triplet in a load_data indices dict, tying pollin's
      linear PD(E)/PA(E) (a constant plus an energy-slope term each) into
      one free quartet for the whole session -- not a scipy fit on
      ixpe_pcube's PD(E) points
  17. compare_to_benchmark — this-work-vs-paper benchmark comparison: takes
      a caller-built list of records (this-work value+error, paper
      value+error, a shared unit, and a mandatory paper citation),
      computes delta/z-score/>sigma_flag flags, and writes a
      provenance-headed comparison CSV (no XSPEC/ixpeobssim required, no
      paper numbers hardcoded)
  18. plot_benchmark_compare — the compare_to_benchmark reviewer figure: a
      forest/caterpillar plot of SIGNED z per record (direction-preserving,
      unlike the CSV's own abs(z) column), so a systematic this-work-vs-
      paper offset is visible at a glance; reads either compare_to_
      benchmark's returned rows list or the CSV it wrote (no XSPEC/
      ixpeobssim required)
  19. ixpe_qu_slice_pcubes — MDP-anchored adaptive q(t)=Q/I, u(t)=U/I
      intra-observation time series for one local ObsID: data-driven N-way
      splitting of a merged/multi-pointing observation at genuine GTI merge
      gaps, a Di Marco+2023 source-brightness-proxy background decision,
      adaptive greedy slice-edge placement anchored to a per-slice MDP99
      target, and mean-event-arrival-time (not window-midpoint) slice
      timestamps
  20. ixpe_qu_lag_dcf — model-independent Edelson & Krolik (1988) discrete
      correlation function between mean-subtracted q(t) and u(t) for every
      subobs of a set of local IXPE ObsIDs, from each ObsID's existing
      qu_slices_<subobs_key>.csv series (produced by ixpe_qu_slice_pcubes),
      with Monte-Carlo permutation-null significance, bootstrap peak-lag/
      peak-value uncertainty, a masked-uniform-grid cross-correlation
      cross-check, and a cross-observation population Spearman correlation
      of the recovered lag against two reconstructed covariates (a source-
      brightness proxy and a combined PD significance) (no XSPEC/ixpeobssim
      required)
  21. ixpe_qu_lag_dcf_plot — the three standard reviewer figures for an
      ixpe_qu_lag_dcf run: per-subobs DCF-vs-lag panels (permutation 99%
      band, tau_peak/dcf_peak marker, bootstrap 68% interval drawn by its
      own true endpoints) plus the two cross-observation population panels
      of tau_peak vs the two reconstructed covariates, annotated with the
      Spearman rho/p read straight from population_correlations.csv (no
      XSPEC/ixpeobssim required)
"""

from __future__ import annotations

TOOL_MANIFEST: dict = {
    # ── Step 1: load data ─────────────────────────────────────────────────
    "research_toolkit.xspec_utils.load_data": {
        "module": "xspec_utils.lib.xspec_help",
        "fn": "load_data",
        "description": (
            "Load every instrument of one epoch (IXPE, NICER, NuSTAR spectra) into the XSPEC "
            "session and return the nested spectrum-indices dict that fit_*, xspec_status, "
            "plot_fit and plot_iqu consume. Start here before any joint spectral fit."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "DATA-GROUP ORDER: data groups follow the epoch dict's key order, so LIST IXPE "
                "FIRST -- IXPE DU1 must be data group 1, because the polarized fit models "
                "require DG1 to be an IXPE (I/Q/U) group and any other DG1 makes the fit "
                "diverge. Recommended key order: IXPE, NICER, NuSTAR, skipping any instrument "
                "the epoch lacks (absent instruments are simply skipped). Data-group order is a "
                "separate matter from the cross-calibration reference, which each fit_* tool "
                "selects internally (NuSTAR by default).\n"
                "MULTI-OBSID IXPE: epoch['IXPE'] may list any number of ObsIDs (e.g. combining "
                "a validation pair for one tied joint polconst/pollin fit). Each ObsID resolves "
                "its own centroid, response files/BACKSCAL and IRF calibration epoch "
                "independently, occupying consecutive gapless data groups (3 DGs per ObsID); "
                "Q/U gain slopes tie only to that same ObsID+DU's own I gain slope, staying "
                "within one ObsID and one DU.\n"
                "RETURNS: a nested dict, instrument -> obsid -> detector -> [XSPEC spectrum "
                "indices]. Pass it unchanged to fit_*, xspec_status, plot_fit and plot_iqu.\n"
                "EXPECTED OUTPUT NOISE: on first load XSPEC prints harmless warnings about "
                "unresolved response-file paths (arf/rmf are loaded explicitly, so these are "
                "safe to ignore). For bright sources such as Cyg X-1 the NuSTAR/NICER "
                "background is None by design -- that is a deliberate choice, not missing data."
            ),
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
                "ixpe_changes": {
                    "type": "boolean",
                    "description": (
                        "If true, force a fresh IXPE reduction: centroid the source position "
                        "(find_src.run_extraction -- measured per ObsID, rather than read from "
                        "a fixed coordinate table) and re-run red_ixpe.pha at ixpe_src_rad/"
                        "ixpe_bkg_inner_rad/ixpe_bkg_outer_rad, so each spectrum's own "
                        "RESPFILE/ANCRFILE (arf for I, mrf for Q/U) and BACKSCAL resolve for "
                        "that region. Set this to reproduce a paper's specific region radii. "
                        "If false (default), reuse whatever per-DU PHA1 products are already "
                        "on disk as-is."
                    ),
                },
                "ixpe_src_rad": {
                    "type": "number",
                    "description": (
                        "IXPE source-region radius in arcminutes. Default 1.0. Used only "
                        "when ixpe_changes=true."
                    ),
                },
                "ixpe_bkg_inner_rad": {
                    "type": "number",
                    "description": (
                        "IXPE background annulus inner radius in arcminutes. Default 3.0. "
                        "Used only when ixpe_changes=true."
                    ),
                },
                "ixpe_bkg_outer_rad": {
                    "type": "number",
                    "description": (
                        "IXPE background annulus outer radius in arcminutes. Default 5.0. "
                        "Used only when ixpe_changes=true."
                    ),
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

    # ── polconst spectropolarimetric cross-check (IXPE-only) ──────────────
    "research_toolkit.xspec_utils.fit_tbabs_polconst_powerlaw": {
        "module": "xspec_utils.lib.xspec_help",
        "fn": "fit_tbabs_polconst_powerlaw",
        "description": (
            "Fit tbabs*polconst*powerlaw to every IXPE Stokes I/Q/U triplet, tying one "
            "session-wide polarization degree (PD) and angle (PA) pair -- the polconst "
            "spectropolarimetric cross-check of a PCUBE polarization measurement."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "TIES: polconst's polarization degree (PD, fraction 0-1) and angle (PA, "
                "degrees) are tied into ONE free pair for the whole session -- a single pair "
                "shared across every detector unit and ObsID. This is the XSPEC column of the "
                "Krawczynski+2022 Cyg X-1 measurement.\n"
                "PREREQUISITE: call load_data first with ixpe_changes=true (a fresh, "
                "freshly-centroided reduction) so each spectrum's own RESPFILE/ANCRFILE and "
                "BACKSCAL resolve for the region actually used; a stale on-disk PHA1 "
                "(ixpe_changes=false) can carry dead response paths and another region's "
                "BACKSCAL, producing an unphysically poor fit. Pass an IXPE-only indices dict "
                "-- Model() here applies to every currently active data group, so this is a "
                "single-instrument fit rather than a joint multi-instrument one. One ObsID in "
                "indices['IXPE'] gives a single-ObsID fit; several give a joint combined fit, "
                "still with exactly one tied PD/PA pair.\n"
                "PROCEDURE: notices 2-8 keV ('erange') on every IXPE spectrum, then runs "
                "Fit.perform() to convergence.\n"
                "ERRORS: PD_err/PA_err/Gamma_err are the covariance-matrix parameter-table "
                "1-sigma (Parameter.sigma), scaled to 'conf_level' via the chi-squared "
                "quantile for 1 degree of freedom. Raises RuntimeError naming the parameter "
                "when any of the three has no valid covariance sigma, in place of returning a "
                "meaningless 0 or negative error.\n"
                "EXPECTED CHI2: for Cyg X-1 the global reduced chi2 typically exceeds 2, "
                "because a single power law does not fully describe the Stokes-I continuum "
                "over 2-8 keV (polconst leaves I unpolarised). polconst multiplies only Q/U, "
                "so PD/PA are set by the well-fit polarization sector (reduced chi2 ~ 1) and "
                "stay robust despite the elevated global chi2.\n"
                "OUTPUT: every fitted number is returned in the result dict (XSPEC's terminal "
                "output is muted, so downstream steps need nothing from stdout). Companion "
                "figure: plot_fit."
            ),
            "properties": {
                "indices": {
                    "type": "object",
                    "additionalProperties": True,
                    "description": (
                        "Return value of load_data (must contain 'IXPE') — pass the "
                        "object directly, do NOT JSON-encode it."
                    ),
                },
                "gamma_init": {
                    "type": "number",
                    "description": "Initial powerlaw photon index. Default 2.0 (Cyg X-1 starting point).",
                },
                "pd_init": {
                    "type": "number",
                    "description": (
                        "Initial polconst polarization fraction (0-1, NOT percent). "
                        "Default 0.04 (Cyg X-1 starting point)."
                    ),
                },
                "pa_init": {
                    "type": "number",
                    "description": "Initial polconst polarization angle in degrees. Default -20.0 (Cyg X-1 starting point).",
                },
                "nh": {
                    "type": "number",
                    "description": (
                        "TBabs equivalent hydrogen column, in units of 1e22 cm^-2. "
                        "Default 0.6 (Cyg X-1 ISM value)."
                    ),
                },
                "freeze_nh": {
                    "type": "boolean",
                    "description": "Freeze nH at the supplied value. Default true; set false to let nH float.",
                },
                "erange": {
                    "type": "array",
                    "items": {"type": "number"},
                    "description": "[lo, hi] keV range noticed on every IXPE spectrum before fitting. Default [2.0, 8.0].",
                },
                "conf_level": {
                    "type": "number",
                    "description": (
                        "Confidence level the reported covariance-matrix 1-sigma "
                        "errors are scaled to. Default 0.683 (1 sigma, Delta-chi2=1.0 "
                        "for one interesting parameter); pass e.g. 0.90 for a 90% "
                        "interval."
                    ),
                },
                "n_iterations": {
                    "type": "integer",
                    "description": "Fit.nIterations per convergence round. Default 100.",
                },
                "log_path": {
                    "type": "string",
                    "description": "If given, write the verbose AllModels.show()/Fit.show() summary here.",
                },
            },
            "required": ["indices"],
        },
    },

    # ── pollin linear-energy-dependence fit (IXPE-only) ────────────────────
    "research_toolkit.xspec_utils.fit_tbabs_pollin_powerlaw": {
        "module": "xspec_utils.lib.xspec_help",
        "fn": "fit_tbabs_pollin_powerlaw",
        "description": (
            "Fit tbabs*pollin*powerlaw to every IXPE Stokes I/Q/U triplet, tying one "
            "session-wide quartet for linear energy-dependent polarization PD(E) and PA(E) "
            "-- the XSPEC route to a polarization energy slope."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "TIES: pollin's PD(E)=A1+Aslope*(E-1keV) (A1 fraction 0-1, Aslope fraction "
                "per keV) and PA(E)=psi1+psislope*(E-1keV) (psi1 degrees, psislope degrees "
                "per keV) are tied into ONE free (A1, Aslope, psi1, psislope) quartet for the "
                "whole session -- a single quartet shared across every detector unit and "
                "ObsID. This is the measurement route for a linear polarization "
                "energy-dependence: fitting the folded model in XSPEC keeps the result "
                "toolkit-grounded, where a regression fitted inline over ixpe_pcube's PD(E) "
                "points would be an analysis step performed outside the toolkit.\n"
                "PREREQUISITE: same load path as fit_tbabs_polconst_powerlaw -- call load_data "
                "first with ixpe_changes=true (a fresh, freshly-centroided reduction) so each "
                "spectrum's own RESPFILE/ANCRFILE and BACKSCAL resolve for the region actually "
                "used; a stale on-disk PHA1 (ixpe_changes=false) can carry dead response paths "
                "and another region's BACKSCAL, producing an unphysically poor fit. Pass an "
                "IXPE-only indices dict -- Model() here applies to every currently active data "
                "group, so this is a single-instrument fit. One ObsID in indices['IXPE'] gives "
                "a single-ObsID fit; several give a joint combined fit, still with exactly one "
                "tied pollin quartet.\n"
                "PROCEDURE: notices 2-8 keV ('erange') on every IXPE spectrum, then runs "
                "Fit.perform() to convergence.\n"
                "ERRORS: every _err is the covariance-matrix parameter-table 1-sigma "
                "(Parameter.sigma), scaled to 'conf_level' via the chi-squared quantile for 1 "
                "degree of freedom. Raises RuntimeError naming the parameter when any of "
                "PD_A/PD_slope/PA_A/PA_slope/Gamma has no valid covariance sigma, in place of "
                "returning a meaningless 0 or negative error.\n"
                "OUTPUT: every fitted number is returned in the result dict (XSPEC's terminal "
                "output is muted). Companion figure: plot_fit; the PD(E)/PA(E) "
                "energy-dependence overlay across ObsIDs is a separate tool."
            ),
            "properties": {
                "indices": {
                    "type": "object",
                    "additionalProperties": True,
                    "description": (
                        "Return value of load_data (must contain 'IXPE') — pass the "
                        "object directly, do NOT JSON-encode it."
                    ),
                },
                "gamma_init": {
                    "type": "number",
                    "description": "Initial powerlaw photon index. Default 2.0 (Cyg X-1 starting point).",
                },
                "pd_a_init": {
                    "type": "number",
                    "description": (
                        "Initial pollin PD at 1 keV (A1, fraction 0-1, NOT "
                        "percent). Default 0.02."
                    ),
                },
                "pd_slope_init": {
                    "type": "number",
                    "description": (
                        "Initial pollin PD energy slope (Aslope, fraction per keV). "
                        "Default 0.004 (Cyg X-1's PD rising roughly 2.6%->5.0% over "
                        "2-8 keV)."
                    ),
                },
                "pa_a_init": {
                    "type": "number",
                    "description": "Initial pollin PA at 1 keV (psi1, degrees). Default -20.0.",
                },
                "pa_slope_init": {
                    "type": "number",
                    "description": "Initial pollin PA energy slope (psislope, degrees per keV). Default 0.0.",
                },
                "nh": {
                    "type": "number",
                    "description": (
                        "TBabs equivalent hydrogen column, in units of 1e22 cm^-2. "
                        "Default 0.6 (Cyg X-1 ISM value)."
                    ),
                },
                "freeze_nh": {
                    "type": "boolean",
                    "description": "Freeze nH at the supplied value. Default true; set false to let nH float.",
                },
                "erange": {
                    "type": "array",
                    "items": {"type": "number"},
                    "description": "[lo, hi] keV range noticed on every IXPE spectrum before fitting. Default [2.0, 8.0].",
                },
                "conf_level": {
                    "type": "number",
                    "description": (
                        "Confidence level the reported covariance-matrix 1-sigma "
                        "errors are scaled to. Default 0.683 (1 sigma, Delta-chi2=1.0 "
                        "for one interesting parameter); pass e.g. 0.90 for a 90% "
                        "interval."
                    ),
                },
                "n_iterations": {
                    "type": "integer",
                    "description": "Fit.nIterations per convergence round. Default 100.",
                },
                "log_path": {
                    "type": "string",
                    "description": "If given, write the verbose AllModels.show()/Fit.show() summary here.",
                },
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
            "Time-resolved 2-8 keV IXPE PCUBE polarimetry for one local ObsID: splits the "
            "observation into elapsed-time segments and measures PD, PA, normalized Stokes "
            "and MDP99 per segment, for testing intra-observation polarization constancy. "
            "Data only -- pair it with a plotter tool for figures."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "METHOD: ixpeobssim (xpselect -> xpbin PCUBE, weights=True) with DU1/2/3 "
                "merged. Segments are natural elapsed-time windows of 'segment_duration_s' "
                "(TSTART + k*segment_duration_s, i.e. elapsed time rather than equalised "
                "exposure), plus one integrated whole-observation reduction. Suited to a "
                "bright, source-dominated target already local at level >= 2 -- confirm "
                "locality with a data-catalog tool first, and acquire data with the "
                "archive-acquisition tool.\n"
                "EMPTY-WINDOW SKIPPING: each window's GTI-clipped good exposure is checked "
                "against 'min_seg_exposure_s' before reduction, and thinner windows are "
                "skipped without running xpselect/xpbin -- this matters for observations with "
                "a long inter-pointing gap, where whole windows can fall inside the gap.\n"
                "RESUMABLE: rows are appended and flushed as each segment finishes, so a call "
                "that hits the sandbox time cap leaves a genuine checkpoint. Call again with "
                "the same arguments (resume=True, the default) to compute only the remaining "
                "segments, or pass resume=False (or change any argument) to force a full "
                "rebuild.\n"
                "OUTPUT: writes/updates <obs_db_path>/<obsid>/pcube_timeseg_<obsid>.csv with "
                "one row per segment plus the integrated row (segment_id, t_start_mjd, "
                "t_stop_mjd, exposure_s, PD_pct, PD_err, PA_deg, PA_err, QI, QI_err, UI, "
                "UI_err, MDP99, SIGNIF; the integrated row carries segment_id=-1) under a "
                "provenance comment header, and records it in the mission catalog at level 3 "
                "with a 'complete' flag. PD_pct and MDP99 are PERCENT here (ixpe_pcube reports "
                "dimensionless fractions instead); PA_deg is degrees.\n"
                "RETURNS: 'csv_path', 'n_segments' (resumed + newly computed), "
                "'n_skipped_empty', 'n_expected', 'complete' (bool: every expected segment "
                "plus the integrated row is present), 'resumed' (bool), 'integrated' (dict), "
                "'head'/'tail' (first/last 3 segment-row dicts), 'stdout' and 'exit_code' (1 "
                "with 'stderr' on failure, recording nothing). Read csv_path for the full "
                "table -- the return value carries only head/tail, since a long observation's "
                "per-segment table can exceed the sandbox transport limit."
            ),
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
                        "each segment and the integrated cube. Default false, which suits "
                        "bright source-dominated targets (e.g. Cyg X-1, 4U 1630-47). Set true "
                        "when the source is faint enough that instrumental/diffuse background "
                        "is a non-negligible fraction of the source-aperture counts."
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
            "Reuse or reduce the per-DU IXPE PHA1 (Stokes I) source+background spectra for "
            "one local ObsID -- the PHA products that the hardness-ratio and load_data tools "
            "read. Run this first for any ObsID lacking PHA1 products."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "INPUT: locates Level-2 DU1/DU2/DU3 event files under "
                "<obs_db_path>/<obsid>/event_l2/, preferring '_rej' "
                "instrumental-background-rejected files and falling back to plain "
                "evt2_v*.fits. Run this on ObsIDs already held locally -- confirm locality "
                "with a data-catalog tool first, and acquire missing data with the "
                "archive-acquisition tool.\n"
                "REUSE PATH: when the per-DU src+bkg PHA1 files already exist on disk they "
                "are reused as-is, leaving the reduction and region files untouched and "
                "returning 'reused'=true with null 'ra'/'dec'.\n"
                "FRESH-REDUCTION PATH: otherwise the source RA/Dec is centroided from the DU1 "
                "event file via find_src.run_extraction (measured from the events themselves, "
                "rather than read from FITS headers or a fixed source-coordinate table), "
                "writing DS9 'src.reg'/'bkg.reg' under <obs_db_path>/<obsid>/; xpselect+xpbin "
                "PHA1 then runs per DU with a circular source region ('src_rad') and annular "
                "background ('bkg_inner_rad'-'bkg_outer_rad') centred there, at the IRF epoch "
                "auto-selected as the installed epoch nearest the observation date.\n"
                "CATALOG: every PHA1 product on disk after the call (3 DUs x src/bkg) is "
                "recorded in the mission catalog at level 3 via obsdb_catalog.neutral.advance, "
                "whether reused or freshly reduced, so a downstream hardness-ratio or "
                "load_data call always finds real, catalogued products.\n"
                "RETURNS: 'src_pha'/'bkg_pha' (each a 3-element DU1/DU2/DU3 list of Stokes-I "
                "PHA1 paths), 'ra'/'dec' (null when reused), 'reused' (bool), "
                "'src_reg'/'bkg_reg' (null when reused), 'stdout', 'exit_code' (1 with "
                "'stderr' on failure -- e.g. no local Level-2 files, or a DU carrying other "
                "than exactly one file -- recording nothing)."
            ),
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
            "Intra-observation PA/PD/Stokes constancy statistics across one or more local "
            "IXPE ObsIDs: chi-squared variability tests, intrinsic scatter, peak-to-peak "
            "ranges and structure functions, computed from existing pcube_timeseg CSVs. "
            "Data only -- pair it with a plotter tool for figures."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "PREREQUISITE: reads each ObsID's existing "
                "<obs_db_path>/<obsid>/pcube_timeseg_<obsid>.csv, so run ixpe_pcube_timeseg "
                "first. This tool computes statistics on those tables and reduces no data "
                "itself, so it needs neither XSPEC nor ixpeobssim and runs fast.\n"
                "VARIABILITY CALL: the AUTHORITATIVE test is the 2D Stokes one, "
                "chi2_2D=chi2_QI+chi2_UI with dof_2D=2*(N-1), flagging 'variable' when p_2D < "
                "'variability_p_threshold'. The PD and PA chi-squared values are secondary "
                "and descriptive. All constancy tests are inverse-variance weighted, and PA "
                "uses mod-180 circular handling -- PA residuals wrap to (-90,+90], so a "
                "+89/-89 deg flip scores as ~2 deg rather than ~178.\n"
                "SECONDARY STATISTICS: moment-estimator intrinsic scatter for PD and PA "
                "(sigma_int=sqrt(max(0, sample_variance - mean_measurement_variance)), flagged "
                "as an upper limit when non-positive); peak-to-peak PD/PA ranges; first-order "
                "structure functions SF_PD/SF_PA at 'sf_lags_ks' (NaN for a lag whose "
                "tolerance window holds no segment pair, and computed only for an ObsID with "
                ">= 'min_seg_for_sf' segments); and, across all requested ObsIDs together, a "
                "Spearman rank correlation of each ObsID's mean MDP99 against its reduced "
                "chi-squared for PD and for PA (NaN with a note when fewer than 3 ObsIDs have "
                "finite values).\n"
                "MISSING INPUTS: an ObsID with no pcube_timeseg CSV, or one lacking the "
                "provenance header ixpe_pcube_timeseg writes, is recorded as an error entry "
                "and skipped, keeping the run alive and the omission visible.\n"
                "OUTPUT: writes pa_pd_stability_catalog.csv (one row per ObsID) and "
                "pa_pd_variability_summary.txt (per-ObsID verdict, secondary stats, the two "
                "Spearman results, and the source CSV paths) under 'output_dir'.\n"
                "RETURNS: 'per_obs' (obsid -> stats dict, including 'sf_pd'/'sf_pa'), "
                "'spearman', 'catalog_csv', 'summary_txt', 'exit_code' (always 0 -- per-ObsID "
                "problems are recorded rather than raised)."
            ),
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
            "Renders the four standard IXPE PA/PD stability reviewer figures (time series, "
            "Stokes plane, chi2 histogram, structure function) from the products of "
            "ixpe_pa_pd_stability."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "PREREQUISITE: run ixpe_pa_pd_stability first for 'catalog_csv', and "
                "ixpe_pcube_timeseg before that for each ObsID's pcube_timeseg CSV. Pure "
                "matplotlib (Agg backend) importing neither XSPEC nor ixpeobssim, so it runs "
                "fast without the science stack.\n"
                "FIGURES: writes exactly four files, appending '_pa_pd_timeseries', "
                "'_stokes_plane', '_chi2_hist', '_structure_function' (plus '.png'/'.pdf') to "
                "'save_path_prefix'. (1) Per-ObsID PA (top) / PD (bottom) panels against time "
                "elapsed in kiloseconds from the earliest segment's t_start, overlaying the "
                "catalog's weighted-mean +/-1sigma band for any ObsID present in catalog_csv. "
                "(2) Per-ObsID Q/I vs U/I Stokes-plane scatter colour-coded by elapsed time, "
                "with 1sigma error crosses and the origin marked. (3) Histograms of "
                "rchi2_PD/rchi2_PA over EVERY row of catalog_csv (the full catalog, beyond the "
                "requested 'obsids'), with a reference line at 1.0. (4) PA/PD structure "
                "functions against lag, one line per ObsID holding >= 'sf_min_segments' "
                "segments; a shorter ObsID is omitted from this panel alone and still appears "
                "in the other three, and when none qualify the panel carries an explicit "
                "annotation in place of empty axes.\n"
                "FAILURE MODES: an ObsID with no readable pcube_timeseg CSV is skipped and "
                "noted in stdout; when every requested ObsID is unreadable, exit_code=1 with "
                "stderr set and nothing written. 'catalog_csv' is a required upstream product "
                "-- a missing or wrong-schema catalog_csv is a hard error (exit_code=1, stderr "
                "set, no files written) rather than a silent degrade to an empty catalog. That "
                "is a separate case from an ObsID merely absent from an otherwise-readable "
                "catalog_csv, which still gets its time-series and Stokes-plane panels and "
                "loses only its weighted-mean overlay.\n"
                "RETURNS: 'figures' (4 FigureRef dicts), 'structure_function' (obsid -> the "
                "lags_ks/sf_pa/sf_pd values actually plotted, for verification), 'stdout', "
                "'exit_code'."
            ),
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
            "Energy-resolved IXPE PCUBE polarimetry (PD, PA, normalized Stokes Q/U, MDP99) "
            "for one local ObsID or a native multi-ObsID combination, background-subtracted; "
            "writes one CSV row per energy band. Data only -- pair it with a plotter tool "
            "for figures."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "Reduces Level-2 event files already on disk -- e.g. reproducing "
                "Krawczynski+2022 Cyg X-1 energy-resolved polarimetry (Fig. S5/Table S2), or "
                "its native multi-ObsID (May 01002901 + June 01250101) combined cube.\n"
                "UNITS: PD, PD_err, QN, QN_err, UN, UN_err and MDP99 are DIMENSIONLESS "
                "FRACTIONS -- the units the cube itself stores (ixpe_pcube_timeseg reports "
                "PD_pct/MDP99 in percent instead); PA and PA_err are degrees. PD, PD_err, PA, "
                "PA_err, QN, QN_err, UN, UN_err, MDP99, SIGNIF, E_LO and E_HI are each read "
                "straight off the background-subtracted source cube's polarization()/attribute "
                "accessors.\n"
                "IRF CALIBRATION: the epoch is chosen automatically as the installed epoch "
                "nearest the midpoint across every requested ObsID; this tool takes no epoch "
                "argument.\n"
                "OUTPUT: a provenance-headed CSV, one row per energy band; a single-edge-pair "
                "'ebinning' yields exactly one broadband row (ixpe_pcube_timeseg's "
                "segment_id=-1 sentinel row has no counterpart here). A single ObsID's CSV "
                "lands at <obs_db_path>/<obsid>/ixpe_pcube_<obsid>.csv; a combined request's "
                "CSV lands at <obs_db_path>/ixpe_pcube_<obsid1>_<obsid2>_....csv (mission "
                "root, since a combined product belongs to no single ObsID directory -- "
                "mirroring bhxrb_hardness_ratios.csv's placement). Records one kind='pcube' "
                "product per contributing ObsID via obsdb_catalog.neutral.advance (level 3), "
                "all pointing at the same CSV for a combined request.\n"
                "RETURNS: 'bands' (list of per-band dicts), 'provenance' (obsids, ebinning, "
                "src_rad, bkg_inner_rad, bkg_outer_rad, ra, dec, irf_date, combined:bool), "
                "'csv_path', 'stdout', 'exit_code'.\n"
                "FAILS with exit_code=1, writing and recording nothing, when a requested "
                "ObsID has no local Level-2 event files, is missing DU1/DU2/DU3 (or carries "
                "more than one file for a DU), or 'ebinning' carries fewer than 2 edges."
            ),
            "properties": {
                "obs_db_path": {
                    "type": "string",
                    "description": (
                        "The IXPE mission's obs-db directory (e.g. "
                        "'/data/fzhou/XSPEC/Cyg_X-1/IXPE'). Event files are expected at "
                        "<obs_db_path>/<obsid>/event_l2/ for every requested ObsID; the full "
                        "Level-2 file list is located via "
                        "red_ixpe.get_all_file_list(obs_db_path). Run this on ObsIDs already "
                        "held locally -- confirm locality with a data-catalog tool first, and "
                        "acquire missing data with the archive-acquisition tool."
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
                        "passing that ObsID as a plain string. A multi-ObsID list filters and "
                        "concatenates every listed ObsID's Level-2 files into ONE PCUBE call, "
                        "so the returned cube sums raw Stokes I/Q/U/COUNTS across every DU "
                        "*and* every listed ObsID before PD/PA are derived -- the native "
                        "directly-combined reduction. (Merging separately-computed per-ObsID "
                        "PD/PA by inverse variance is a different, approximate operation this "
                        "argument deliberately replaces.)"
                    ),
                },
                "ra": {
                    "type": "number",
                    "description": (
                        "Source RA in decimal degrees, shared across every listed ObsID. "
                        "Required: supply an explicit region centre (ixpe_pcube_timeseg "
                        "auto-centroids, this tool does not, because a multi-ObsID "
                        "combination needs one shared centre)."
                    ),
                },
                "dec": {
                    "type": "number",
                    "description": (
                        "Source Dec in decimal degrees, shared across every listed ObsID. "
                        "Required alongside 'ra' as the explicit shared region centre."
                    ),
                },
                "ebinning": {
                    "type": "array",
                    "items": {"type": "number"},
                    "minItems": 2,
                    "description": (
                        "Energy bin edges in keV (ebinalg='LIST'), e.g. [2,3,4,6,8]. Default "
                        "[2,3,4,6,8] (4 bands, the Krawczynski+2022 Table S2 grid). Supply >= "
                        "2 edges; a 2-edge list (e.g. [2,8]) yields a single broadband row, "
                        "and e.g. 11 equal-width edges over 2-8 keV give a finer energy grid."
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
            "Net Stokes-I hardness ratio HR=(H-S)/(H+S) for a survey of local IXPE ObsIDs, "
            "read from each ObsID's existing per-DU PHA1 products. Data only -- pair it with "
            "a plotter tool for figures."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "PREREQUISITE: reads each ObsID's existing per-DU PHA1 (Stokes I) src(+bkg) "
                "products, the same PHA1 geometry ixpe_reduce_pha produces and that feeds "
                "ixpe_pcube_timeseg's PD -- so PD and HR share identical extraction regions. "
                "Run ixpe_reduce_pha first for any ObsID missing its PHA1 products; this tool "
                "reads PHA rather than re-reducing it.\n"
                "DEFINITION: HR is the community-normalized (H-S)/(H+S), bounded [-1,+1].\n"
                "PER-DU HANDLING: DU1/DU2/DU3 are located independently via pha()'s file "
                "locators and the net rate is summed over every DU actually present. IXPE can "
                "fly with a DU off, so a genuinely-absent DU is dropped and the survey "
                "continues; a DU that is present but whose PHA1 is missing on disk is reported "
                "as that ObsID's error instead.\n"
                "ERRORS: closed-form first-order Gaussian propagation "
                "sigma_HR=(2/(H+S)^2)*sqrt(S^2*sigma_H^2+H^2*sigma_S^2), with each band's net "
                "rate and error read straight from the PHA RATE/STAT_ERR columns (summed over "
                "in-band channels and contributing DUs). The computation reads the PHA files "
                "directly and loads nothing into XSPEC AllData.\n"
                "BANDS: 'soft_band'/'hard_band' channels resolve per PHA file from its own "
                "RESPFILE RMF EBOUNDS, so any band works beyond the survey default 2-4/4-8 "
                "keV. When a PHA1's RESPFILE predates a reduction-environment migration and no "
                "longer resolves on disk (the real state of most already-reduced Cyg X-1 IXPE "
                "PHA1 products), the same DU's installed CALDB RMF supplies the grid instead -- "
                "an identical EBOUNDS grid, since it is a fixed instrumental characteristic "
                "shared across every IRF date and DU.\n"
                "PER-OBSID FAILURES: an ObsID with no local Level-2 event files, or one the "
                "reduction raises on, is recorded in the returned 'errors' dict and excluded "
                "from the CSV, keeping the rest of the survey alive and the omission visible "
                "in place of a dummy row.\n"
                "OUTPUT: writes <obs_db_path>/bhxrb_hardness_ratios.csv (columns obsid, "
                "soft_rate, soft_err, hard_rate, hard_err, HR, HR_err), replacing any previous "
                "run, and records one 'hardness_ratio' product per successful ObsID in the "
                "mission catalog at level 3 via obsdb_catalog.neutral.advance, with 'params' "
                "carrying the real soft_band/hard_band/bkg/dus_used for that ObsID.\n"
                "RETURNS: 'csv_path', 'rows' (successful per-obsid rows), 'errors' (obsid -> "
                "message), 'stdout', 'exit_code' (1 with 'stderr' only when EVERY ObsID "
                "failed, writing and recording nothing in that case; 0 otherwise, with "
                "per-ObsID failures surfaced in 'errors')."
            ),
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
            "Correlates IXPE polarization degree (PD) against hardness ratio (HR) across a "
            "survey of BH-XRB ObsIDs: per-source and combined Spearman correlations, a "
            "soft-vs-hard KS test and a detection-bias test. Data only -- pair it with a "
            "plotter tool for figures."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "PREREQUISITES: for each obsid in 'obsids_with_sources', reads "
                "PD_pct/PD_err/MDP99 from that obsid's existing "
                "<pcube_csv_dir>/<obsid>/pcube_timeseg_<obsid>.csv integrated (segment_id=-1) "
                "row -- run ixpe_pcube_timeseg first -- and HR/HR_err from 'hr_csv', the "
                "bhxrb_hardness_ratios.csv that the hardness-ratio tool writes -- run that "
                "first too.\n"
                "CLASSIFICATION: detected = PD_pct > MDP99; mdp_flagged = MDP99 > "
                "'mdp_flag_pct' is a data-quality note that leaves the row in the sample. An "
                "obsid with no readable pcube CSV, or no row in hr_csv, carries a recorded "
                "note and is excluded from every statistic, keeping the survey alive and the "
                "omission visible in place of a dummy value.\n"
                "STATISTICS OVER DETECTED ROWS ONLY (a non-detection's PD is "
                "MDP99-noise-dominated rather than a physical measurement): per-source "
                "Spearman rank correlation (PD_pct vs HR) for every source holding >= "
                "'min_n_intra' detections, with shorter sources omitted from "
                "intra_source_correlations.csv; the combined Spearman across every detected "
                "point regardless of source; and a two-sample Kolmogorov-Smirnov test between "
                "the detected PD_pct of 'soft_state_sources' vs 'hard_state_sources', both "
                "supplied by the caller so the split serves any source grouping.\n"
                "STATISTIC OVER EVERY OBSID WITH BOTH PD AND HR: a two-sided Mann-Whitney U "
                "test on HR between the detected and non-detected groups, i.e. whether HR "
                "predicts that a detection was made.\n"
                "EMPTY INPUTS: every statistic returns NaN when its input group is empty, or "
                "-- for the two Spearman correlations -- holds fewer than 3 finite pairs.\n"
                "OUTPUT: writes bhxrb_pd_hardness_survey.csv, intra_source_correlations.csv, "
                "combined_correlation.csv, ks_test_soft_vs_hard.csv and "
                "detection_bias_assessment.csv under 'obs_db_path', and records one "
                "'pd_hardness_survey' product per obsid with both PD and HR in the mission "
                "catalog at level 3 via obsdb_catalog.neutral.advance.\n"
                "RETURNS: 'rows' (one per requested obsid), 'skipped' (obsids with a note), "
                "'per_source_correlations', 'spearman_combined', 'ks_result', 'mwu_result', "
                "the five output CSV paths, and 'exit_code' (always 0 -- a per-obsid problem "
                "is recorded rather than raised)."
            ),
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
            "Renders the five standard PD-hardness survey reviewer figures (combined and "
            "per-source PD vs HR scatters, PD distribution comparison, HR detection bias, "
            "MDP99-exposure feasibility) from bhxrb_pd_hardness_survey.csv."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "PREREQUISITE: run bhxrb_pd_hr_survey first for 'survey_csv'. Pure matplotlib "
                "(Agg backend) importing neither XSPEC nor ixpeobssim, so it runs fast without "
                "the science stack.\n"
                "FIGURES: writes exactly five files, appending 'pd_vs_hr_combined', "
                "'pd_vs_hr_per_source', 'pd_distribution_comparison', 'hr_detection_bias', "
                "'feasibility_mdp99_exposure' (plus '.png') to 'save_path_prefix'. (1) PD_pct "
                "vs HR scatter over every detected row, colour-coded by source, annotated with "
                "the combined Spearman rho/p. (2) The same scatter split one panel per source "
                "(every source with a usable row), each annotated with its own Spearman rho/p "
                "once it holds >= 3 detections. (3) Detected PD_pct histograms for "
                "'soft_state_sources' vs 'hard_state_sources' -- both caller-supplied, the same "
                "grouping bhxrb_pd_hr_survey's KS test uses -- annotated with the KS "
                "statistic/p. (4) HR histograms for the detected vs non-detected groups (every "
                "usable row), annotated with the Mann-Whitney U/p. (5) MDP99 (percent) vs "
                "exposure (kiloseconds) scatter, colour-coded by detection, with an "
                "'mdp_flag_pct' reference line; exposure is read per obsid from its own "
                "<obs_db_path>/<obsid>/pcube_timeseg_<obsid>.csv integrated row rather than "
                "from survey_csv, so an obsid whose pcube_timeseg CSV is locally unreadable is "
                "skipped from this panel alone and noted in stdout.\n"
                "EMPTY GROUPS: when a group is empty (no detections, no soft/hard match, no "
                "readable exposure data, ...) the affected panel carries an explicit "
                "annotation, keeping an empty plot from reading as a real null result.\n"
                "REQUIRED UPSTREAM PRODUCT: a 'survey_csv' that is missing, wrong-schema, or "
                "entirely skipped (every row carrying a note, i.e. every requested obsid "
                "failed upstream) is a hard error -- exit_code=1 with stderr set and no files "
                "written.\n"
                "RETURNS: 'figures' (5 FigureRef dicts), 'plotted' (the exact per-panel values "
                "drawn -- 'combined_points', 'per_source_points', 'distribution_groups', "
                "'detection_bias_hr', 'mdp99_exposure' -- for verification), 'stdout', "
                "'exit_code'."
            ),
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
            "Renders the PD(E)/PA(E) and normalized-Stokes Q-U-plane reviewer figure for an "
            "ixpe_pcube CSV, with an optional fitted linear energy-dependence overlay."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "PREREQUISITE: run ixpe_pcube first for 'csv_path' (Krawczynski+2022 Fig. S5 "
                "style). Pure matplotlib (Agg backend) importing neither XSPEC nor ixpeobssim, "
                "so it runs fast without the science stack.\n"
                "FIGURE: exactly ONE PNG with three panels -- PD vs band-centre energy (top "
                "left), PA vs band-centre energy (bottom left, sharing the energy axis), and "
                "QN vs UN, the normalized Stokes/Q-U plane (right, colour-coded by band-centre "
                "energy). Each carries 1-sigma error bars taken straight from csv_path's own "
                "columns: PD, PD_err, PA, PA_err, QN, QN_err, UN, UN_err are plotted exactly as "
                "ixpe_pcube wrote them (dimensionless fractions for PD/QN/UN, degrees for PA), "
                "with the values reproduced as stored.\n"
                "ENERGY-DEPENDENCE OVERLAY (the reviewer check for an energy-dependence "
                "claim): pass 'pollin' -- the PD_A/PD_slope/PA_A/PA_slope quartet from "
                "fit_tbabs_pollin_powerlaw's result dict, extra keys ignored -- to draw "
                "PD(E)=PD_A+PD_slope*(E-1 keV) and PA(E)=PA_A+PA_slope*(E-1 keV), the model "
                "actually fit to the spectra. Omit 'pollin' to fall back to an "
                "inverse-variance-weighted linear fit of the CSV's own PD(E)/PA(E) points "
                "(scipy.optimize.curve_fit), drawn the same way.\n"
                "EDGE CASES: a single-band CSV renders one marker and omits the overlay line, "
                "since a single point supports no line. A band carrying any non-finite value "
                "in the plotted columns is skipped, keeping it out of the figure rather than "
                "drawing it as zero.\n"
                "FAILS (no PNG written) when csv_path is absent, lacks the ixpe_pcube "
                "provenance header, carries zero band rows, has every band skipped for "
                "non-finite values, or when 'pollin' is supplied while missing any of "
                "PD_A/PD_slope/PA_A/PA_slope.\n"
                "RETURNS: 'figures' (one FigureRef dict), 'plotted' (the exact "
                "E_LO/E_HI/E_ctr/PD/PD_err/PA/PA_err/QN/QN_err/UN/UN_err arrays drawn, plus "
                "n_bands/n_skipped and 'overlay' -- {'kind': None} when no line was drawable, "
                "else {'kind': 'pollin'|'fallback', PD_A, PD_slope, PA_A, PA_slope, ...} -- for "
                "differential verification against csv_path), 'stdout', 'exit_code' (always 0; "
                "the failure cases above raise instead of degrading to an empty plot)."
            ),
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
                "pollin": {
                    "type": "object",
                    "additionalProperties": True,
                    "description": (
                        "Optional PD_A/PD_slope/PA_A/PA_slope quartet (fraction, "
                        "fraction/keV, degrees, degrees/keV) to overlay as the "
                        "analytic PD(E)/PA(E) model -- pass fit_tbabs_pollin_powerlaw's "
                        "result dict directly (do NOT JSON-encode it; extra keys such "
                        "as Gamma/*_err are ignored). Omit to overlay a scipy "
                        "weighted-linear fit of the CSV's own PD(E)/PA(E) points "
                        "instead."
                    ),
                },
            },
            "required": ["csv_path"],
        },
    },

    # ── This-work-vs-paper benchmark comparison (data product) ─────────────
    "research_toolkit.xspec_utils.compare_to_benchmark": {
        "module": "xspec_utils.lib.benchmark_compare",
        "fn": "compare_to_benchmark",
        "description": (
            "Compares this-work measured quantities against cited paper values, computing "
            "signed deltas and z-scores and flagging significant discrepancies -- the "
            "this-work-vs-paper comparison table a reproduction study needs. Data only -- "
            "pair it with a plotter tool for the forest figure."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "SOURCE OF VALUES: every value in 'records' comes from the caller. The "
                "this-work side typically comes from ixpe_pcube / "
                "fit_tbabs_polconst_powerlaw / fit_tbabs_pollin_powerlaw results via the "
                "lib-only pcube_quantity/fit_quantity extractors (which apply the "
                "fraction<->percent unit contract, so this tool receives values already in "
                "their final units), and the paper side from a literature value the caller "
                "cites. This tool holds no paper numbers of its own and derives no science; it "
                "is pure stdlib and runs standalone without XSPEC or ixpeobssim. "
                "(pcube_quantity delegates to red_ixpe.read_pcube_csv -- the single source of "
                "truth for the ixpe_pcube CSV's provenance-header and band-field contract -- "
                "via a lazy function-local import, so only a pcube_quantity call pulls in "
                "ixpeobssim.)\n"
                "RECORD CONTRACT: each record carries 'label' (non-empty str), 'value_this' "
                "(finite number), 'err_this' (finite number >= 0), 'value_paper' (finite "
                "number), 'err_paper' (finite number >= 0), 'paper_source' (non-empty str, the "
                "literature citation), and EITHER a shared 'unit' (non-empty str) OR the "
                "'unit_this'/'unit_paper' pair together (when 'unit' is also given, all three "
                "must agree). NaN/inf on any of the four numeric fields raises ValueError, so "
                "a failed upstream measurement surfaces instead of propagating into a "
                "z=nan row; a missing or empty 'paper_source' likewise raises, since a value "
                "can only be benchmarked against an attributed source.\n"
                "UNITS: supply both sides already expressed in the same unit -- this tool "
                "performs no conversion. A mismatch between unit_this and unit_paper (or "
                "between either and 'unit') raises ValueError naming the record's label and "
                "the conflicting units, surfacing the conflict rather than resolving it "
                "silently.\n"
                "STATISTICS: delta = value_this - value_paper (signed, positive when this-work "
                "exceeds the paper value); z = abs(delta) / sqrt(err_this**2 + err_paper**2); "
                "flagged = abs(z) > 'sigma_flag'. When err_this and err_paper are both 0 the "
                "z is +inf (a computed zero denominator, distinct from the input-level NaN/inf "
                "refusal that happens before this arithmetic) and is always flagged.\n"
                "OUTPUT: writes a provenance-headed CSV to 'out_csv' (columns label, unit, "
                "value_this, err_this, value_paper, err_paper, paper_source, delta, z, flagged, "
                "under a '# compare_to_benchmark provenance: {...}' JSON comment header "
                "carrying sigma_flag/created_at/produced_by), creating parent directories as "
                "needed. An empty records list raises ValueError('empty records'), and any "
                "record failing validation raises, writing nothing.\n"
                "RETURNS: 'rows' (one dict per record, the exact fields written to the CSV), "
                "'flags' (labels with abs(z) > sigma_flag), 'out_csv' (absolute path written), "
                "'n_records', 'stdout', 'exit_code' (always 0; invalid input raises instead)."
            ),
            "properties": {
                "records": {
                    "type": "array",
                    "minItems": 1,
                    "items": {
                        "type": "object",
                        "additionalProperties": True,
                        "properties": {
                            "label": {"type": "string", "description": "Non-empty row label, e.g. 'PD_2-8keV'."},
                            "value_this": {"type": "number", "description": "This-work measured value. Must be finite -- NaN/inf are refused."},
                            "err_this": {"type": "number", "description": "This-work 1-sigma error. Must be finite and >= 0 -- NaN/inf are refused."},
                            "value_paper": {"type": "number", "description": "Cited paper value, in the SAME unit as value_this. Must be finite -- NaN/inf are refused."},
                            "err_paper": {"type": "number", "description": "Cited paper 1-sigma error. Must be finite and >= 0 -- NaN/inf are refused."},
                            "unit": {
                                "type": "string",
                                "description": (
                                    "Shared unit for both value_this and value_paper, e.g. "
                                    "'percent' or 'deg'. Omit only if supplying the "
                                    "unit_this/unit_paper pair instead."
                                ),
                            },
                            "unit_this": {
                                "type": "string",
                                "description": (
                                    "This-work side's unit. Must be supplied together with "
                                    "unit_paper (never alone); both must agree with each other "
                                    "(and with 'unit' if also given) or the record is refused."
                                ),
                            },
                            "unit_paper": {
                                "type": "string",
                                "description": "Paper side's unit -- see unit_this.",
                            },
                            "paper_source": {
                                "type": "string",
                                "description": "Non-empty literature citation this record is benchmarked against, e.g. 'Krawczynski+2022 Table S2'.",
                            },
                        },
                        "required": ["label", "value_this", "err_this", "value_paper", "err_paper", "paper_source"],
                    },
                    "description": (
                        "Non-empty list of comparison records (pass the list of objects "
                        "directly -- do NOT JSON-encode it). Build value_this/err_this via "
                        "pcube_quantity/fit_quantity (lib helpers) so the fraction<->percent "
                        "unit choice is made once, correctly, and value_paper/err_paper/unit "
                        "match it exactly."
                    ),
                },
                "out_csv": {
                    "type": "string",
                    "description": "Output CSV path. Parent directories are created as needed.",
                },
                "sigma_flag": {
                    "type": "number",
                    "description": "z-score threshold above which a record is flagged. Default 3.0.",
                },
            },
            "required": ["records", "out_csv"],
        },
    },

    # ── This-work-vs-paper benchmark comparison figure (companion plotter) ─
    "research_toolkit.xspec_utils.plot_benchmark_compare": {
        "module": "xspec_utils.lib.benchmark_compare",
        "fn": "plot_benchmark_compare",
        "description": (
            "Renders the compare_to_benchmark reviewer figure: a forest/caterpillar plot of "
            "signed z-scores, one tick per record, with the flagging band shaded."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "AXIS: one horizontal tick per record in the order supplied, with the x-axis "
                "the SIGNED z = (value_this - value_paper) / sqrt(err_this**2 + err_paper**2). "
                "The sign is recomputed here because the CSV's own 'z' column stores abs(z) "
                "for flagging and drops direction -- and the sign is what makes a systematic "
                "offset visible (every PD point pulled to the same side of zero reveals a "
                "consistent this-work-vs-paper discrepancy that an abs(z) table hides). "
                "Heterogeneous quantities (PD in percent, PA in degrees, a slope in "
                "percent/keV, per-DU values, ...) legitimately share this one axis because z "
                "is dimensionless.\n"
                "MARKINGS: a vertical reference line at z=0 and a shaded band at "
                "+/-'sigma_flag'; any record with abs(z) > sigma_flag is drawn in a distinct "
                "colour and annotated with its label. A record whose value_this and "
                "value_paper both carry zero error has an infinite z (signed by delta) and is "
                "clamped to the finite-z axis edge with an explicit '[z=inf]' annotation, "
                "keeping the axis autoscale usable.\n"
                "INPUT: supply EITHER 'rows' (e.g. compare_to_benchmark's own returned 'rows' "
                "list) OR 'csv_path' (the CSV it wrote) -- exactly one; 'rows' wins when both "
                "are given. Leave 'sigma_flag' unset to read the threshold from the CSV's own "
                "provenance header when 'csv_path' is used (falling back to 3.0), so the "
                "plotted band matches the threshold compare_to_benchmark actually applied; "
                "pass an explicit number to override.\n"
                "OUTPUT: exactly ONE PNG.\n"
                "EMPTY AND INVALID INPUT: an empty 'rows' list, or a 'csv_path' with zero data "
                "rows, writes no PNG and returns exit_code=1 with stderr set. A missing "
                "csv_path, a malformed CSV (one that is not a compare_to_benchmark product), "
                "or a row missing a required field ('label', 'err_this', 'err_paper', and "
                "either 'delta' or 'value_this'/'value_paper') raises ValueError.\n"
                "RETURNS: 'figures' (one FigureRef dict, or [] in the empty-input case), "
                "'plotted' (the exact label/z/flagged arrays drawn, plus "
                "sigma_flag/n_records, for differential verification), 'stdout', 'exit_code'."
            ),
            "properties": {
                "rows": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "additionalProperties": True,
                        "properties": {
                            "label": {"type": "string"},
                            "err_this": {"type": "number"},
                            "err_paper": {"type": "number"},
                            "delta": {"type": "number"},
                            "value_this": {"type": "number"},
                            "value_paper": {"type": "number"},
                        },
                    },
                    "description": (
                        "compare_to_benchmark's own returned 'rows' list (pass the object "
                        "directly -- do NOT JSON-encode it) -- or any equivalent list of "
                        "dicts carrying 'label'/'err_this'/'err_paper' and either 'delta' "
                        "or 'value_this'/'value_paper'. Omit to read from csv_path instead."
                    ),
                },
                "csv_path": {
                    "type": "string",
                    "description": (
                        "Path to the CSV written by compare_to_benchmark (its 'out_csv'). "
                        "Ignored if 'rows' is also given. Must carry the "
                        "'# compare_to_benchmark provenance: ...' header line."
                    ),
                },
                "save_path": {
                    "type": "string",
                    "description": "Full output PNG path. Auto-generated if omitted.",
                },
                "title": {
                    "type": "string",
                    "description": "Figure title. Defaults to a generic caption.",
                },
                "sigma_flag": {
                    "type": "number",
                    "description": (
                        "z-score threshold for the shaded band / flag colour. Left "
                        "unset (the default): read from the CSV's own provenance header "
                        "when csv_path is used, else 3.0. Pass an explicit value to "
                        "override either way."
                    ),
                },
            },
            "required": [],
        },
    },

    # ── MDP-anchored adaptive intra-observation q(t)/u(t) slicing ──────────
    "research_toolkit.xspec_utils.ixpe_qu_slice_pcubes": {
        "module": "xspec_utils.skills.session",
        "fn": "ixpe_qu_slice_pcubes",
        "description": (
            "Model-independent intra-observation q(t)=Q/I, u(t)=U/I normalized-Stokes time "
            "series for one local IXPE ObsID, sliced adaptively so each slice reaches a target "
            "MDP99. Data only -- pair it with a plotter tool for figures."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "SLICING CONTRAST: slices are adaptive to a TARGET MDP99, where "
                "ixpe_pcube_timeseg uses fixed-duration, GTI-window-timed segments. DU1/2/3 "
                "are merged.\n"
                "INPUT: locates Level-2 DU1/2/3 event files under "
                "<obs_db_path>/<obsid>/event_l2/, preferring '_rej' "
                "instrumental-background-rejected files and falling back to plain "
                "evt2_v*.fits -- the same convention as ixpe_pcube_timeseg/ixpe_reduce_pha, "
                "which resolves a literature-merged ObsID such as 03002599 to its plain files. "
                "Run this on ObsIDs already held locally -- confirm locality with a "
                "data-catalog tool first.\n"
                "N-WAY SPLIT: attempted only when 'split_multiobs' (or, when omitted, the IXPE "
                "8-digit ObsID convention PPnnnnpp -- 2-digit mission year + 4-digit target ID "
                "+ 2-digit segment number, per the HEASARC IXPE Quick Start guide -- where "
                "segment '99' marks several pointings MERGED into one ObsID, e.g. 03002599) "
                "judges this obsid a merged observation. An ordinary single monitoring pointing "
                "(e.g. Cyg X-1's 02008201/02008301, segment '01') is treated as one whole-span "
                "sub-observation regardless of how large its own GTI gaps are. With the gate "
                "open, the GTI list is partitioned at genuine multi-pointing merge gaps: a gap "
                "must clear an outlier-ratio test on the SORTED inter-GTI gap list (gated by "
                "'split_gap_factor' -- sorted rather than compared against the population "
                "median, which misfires on real IXPE GTI structure where on-orbit dropout gaps "
                "sit orders of magnitude below a genuine month-scale merge gap) AND BOTH scale "
                "floors 'split_gap_min_s' and 'split_gap_min_frac' -- without the floors, a "
                "long campaign's own routine orbital-visibility gaps (IXPE flies a ~96.6 min "
                "orbit, so a routine occultation/SAA gap is at most a few ks) sit at the same "
                "ratio-cliff tier as the real merge gap and get accepted as spurious 'merge "
                "boundaries' (the reported defect: one 386 ks '...99' campaign over-split into "
                "~230 single-slice fragments). At most 'max_subobs' sub-observations are ever "
                "emitted -- when more boundaries survive the ratio test and both floors, only "
                "the largest 'max_subobs'-1 are kept (ties broken by earlier GTI index). No "
                "split gives one sub-obs keyed <obsid>; an N-way split gives <obsid>a, "
                "<obsid>b, ... in time order (bijective base-26 letter suffixes -- 'a'..'z' "
                "then 'aa','ab',... past 26 sub-obs, pure ASCII so a FITS header stays "
                "printable), each analysed independently, with every slice contained inside "
                "one sub-observation.\n"
                "PER SUB-OBSERVATION: (1) BACKGROUND -- bkg=None auto-decides from a "
                "source-region count-rate surface-brightness proxy r_proxy (cps/arcmin^2, Di "
                "Marco+2023 AJ 165 143): faint (r_proxy<1.0) subtracts a background annulus, "
                "intermediate/bright (r_proxy>=1.0) skips it, since a source-dominated "
                "background annulus makes subtraction counterproductive. bkg=True/False "
                "overrides that per-sub-obs decision, and the resolved choice applies "
                "identically to the whole-sub-obs cube and every one of its slices. (2) One "
                "whole-sub-obs weighted, DU-merged PCUBE gives MDP99_obs and total good-time "
                "COUNTS (C_obs); the per-slice MDP target is 'mdp_target' when given, else "
                "max(mdp_frac*PD_obs, mdp_floor) (dimensionless fractions), giving "
                "C_target=C_obs*(MDP99_obs/target)^2. (3) ADAPTIVE GREEDY SLICING places edges "
                "at real event arrival times so each slice's counts approximate C_target, "
                "clipped to ['n_min', 'n_max'] slices per sub-obs: exceeding n_max "
                "re-partitions into exactly n_max equal-count slices, raising the effective "
                "target accordingly, while falling short of n_min is accepted as-is and "
                "flagged noise_limited=true for that sub-obs, so the slice count stays within "
                "what the real events support. A trailing under-target remainder merges into "
                "the previous slice in place of an emitted straggler. (4) Each final slice "
                "re-runs the same weighted DU-merged PCUBE over its own [tmin,tmax) window and "
                "is TIMESTAMPED by the MEAN of the real event arrival times inside that window "
                "(t_mean_met) rather than the window midpoint -- a deliberate divergence from "
                "ixpe_pcube_timeseg that matters across a real intra-sub-obs GTI gap.\n"
                "CENTROID: auto-centroids ra/dec via source extraction when both are omitted "
                "(same as ixpe_pcube_timeseg); supply both or neither.\n"
                "OUTPUT: one provenance-headed qu_slices_<subobs_key>.csv per sub-observation "
                "under <obs_db_path>/<obsid>/ (columns slice_id, t_mean_met, t_start_met, "
                "t_stop_met, good_exposure_s, counts, q, q_err, u, u_err, PD_pct, PA_deg, "
                "MDP99_pct, under_target; MET times in seconds; q/u are dimensionless "
                "normalized Stokes fractions Q/I, U/I; PD_pct and MDP99_pct are percent, "
                "matching ixpe_pcube_timeseg), plus one 'qu_slices' product per sub-obs CSV "
                "recorded in the mission catalog at level 3 via obsdb_catalog.neutral.advance, "
                "all against the real 9-digit obsid -- every sub-obs CSV is written and "
                "recorded regardless of the bounded-return truncation below.\n"
                "FAILS with exit_code=1/stderr, writing and recording nothing, on a "
                "missing/incomplete DU set, a mismatched ra/dec pair, or a sub-observation "
                "with no GTI overlap or no usable counts over its own span.\n"
                "BOUNDED RETURN: 'stdout' is always a short, fixed-size summary line (obsid, "
                "sub-obs count, total slice count, the directory the CSVs landed in) -- it "
                "never enumerates per-sub-obs paths, which is what let a many-way split's "
                "response exceed the sandbox transport chunk limit before this fix. A sidecar "
                "'qu_slice_subobs_listing_<obsid>.csv' (subobs_key, csv_path, n_slice; one row "
                "per sub-observation, chronological) is ALWAYS written next to the sub-obs CSVs "
                "-- not catalogued, a listing index only -- and 'per_subobs_truncated' (bool), "
                "'n_subobs_total' (the true count) and 'subobs_listing_path' are ALWAYS present "
                "in the return value. When the true sub-observation count exceeds "
                "'max_return_subobs', 'subobs_keys' and 'per_subobs' in the RETURN VALUE are "
                "truncated to the first 'max_return_subobs' (chronological order) and "
                "'per_subobs_truncated'=true; otherwise every sub-obs is returned and "
                "'per_subobs_truncated'=false. Truncation affects only this returned payload, "
                "never which CSVs are written or catalogued.\n"
                "RETURNS: 'obsid', 'subobs_keys' (list, time order, possibly truncated), "
                "'per_subobs' (key -> 'mdp99_obs_pct', 'r_proxy', 'subtract_background', "
                "'uses_rej', 'n_slice', 'noise_limited', 'whole_pd_pct', 'whole_pd_err', "
                "'whole_pa_deg', 'whole_signif', 'whole_counts', 'whole_exposure_s', "
                "'csv_path', plus the diagnostics 'mdp_target_frac_used' (the resolved "
                "per-slice MDP99 target fraction actually used), 'c_target_counts' (the "
                "derived target raw-proxy count), 'n_natural_slices' (the greedy walk's slice "
                "count before any n_max clipping) and 'clipped_to_n_max' (bool)), the "
                "whole-observation GTI diagnostics 'n_gti', 'median_gti_gap_s' (0.0 for a "
                "single-GTI observation), 'total_elapsed_span_s' and 'split_gaps_s' (one entry "
                "per accepted split boundary gap duration in seconds, chronological order; "
                "len(split_gaps_s) == len(subobs_keys)-1; empty when no split occurred), plus "
                "'per_subobs_truncated' (bool), 'n_subobs_total' (int, the true count "
                "irrespective of truncation) and 'subobs_listing_path' (str) -- these three are "
                "ALWAYS present, plus 'stdout', 'exit_code'. Read the CSV for the full "
                "per-slice table -- the return value carries the summary only."
            ),
            "properties": {
                "obs_db_path": {
                    "type": "string",
                    "description": (
                        "The IXPE mission's obs-db directory (e.g. "
                        "'/data/fzhou/XSPEC/Cyg_X-1/IXPE'). Event files are expected at "
                        "<obs_db_path>/<obsid>/event_l2/; qu_slices_<subobs_key>.csv lands "
                        "under <obs_db_path>/<obsid>/."
                    ),
                },
                "obsid": {
                    "type": "string",
                    "description": "9-digit IXPE Observation ID, as a string, e.g. '01002901'.",
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
                    "description": (
                        "Lower edge of the single energy bin, in keV -- also the r_proxy "
                        "background-decision reference band. Default 2.0."
                    ),
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
                        "Background-subtraction override, applied identically to the "
                        "whole-sub-obs cube and every slice within each sub-observation. "
                        "Omit (null/default) to auto-decide per sub-observation from the "
                        "r_proxy source-brightness proxy (Di Marco+2023): faint "
                        "(r_proxy<1.0 cps/arcmin^2) subtracts, intermediate/bright does not."
                    ),
                },
                "bkg_inner_rad": {
                    "type": "number",
                    "description": "Background annulus inner radius in arcminutes. Default 3.0.",
                },
                "bkg_outer_rad": {
                    "type": "number",
                    "description": "Background annulus outer radius in arcminutes. Default 5.0.",
                },
                "mdp_target": {
                    "type": "number",
                    "description": (
                        "Per-slice MDP99 target as a dimensionless fraction (e.g. 0.02 for "
                        "2%, NOT percent). Omit (null/default) to derive it per "
                        "sub-observation from mdp_frac/mdp_floor instead."
                    ),
                },
                "mdp_frac": {
                    "type": "number",
                    "description": (
                        "When mdp_target is omitted: the per-slice MDP99 target is "
                        "max(mdp_frac * whole-sub-obs PD, mdp_floor), i.e. mdp_frac is the "
                        "fraction of the sub-observation's own detected polarization degree "
                        "considered an adequate per-slice sensitivity. Default 0.5."
                    ),
                },
                "mdp_floor": {
                    "type": "number",
                    "description": (
                        "Floor (dimensionless fraction, e.g. 0.015 for 1.5%) on the derived "
                        "per-slice MDP99 target when mdp_target is omitted -- prevents an "
                        "unreachably tight target for a weakly-polarized or short "
                        "sub-observation. Default 0.015."
                    ),
                },
                "n_min": {
                    "type": "integer",
                    "description": (
                        "Minimum slice count per sub-observation. Fewer natural slices than "
                        "this are still emitted as-is (never invented); that sub-obs is "
                        "flagged noise_limited=true instead. Default 8."
                    ),
                },
                "n_max": {
                    "type": "integer",
                    "description": (
                        "Maximum slice count per sub-observation. More natural slices than "
                        "this are re-partitioned into exactly n_max equal-count slices "
                        "(raising the effective per-slice target accordingly). Default 25."
                    ),
                },
                "split_gap_factor": {
                    "type": "number",
                    "description": (
                        "N-way split threshold: an inter-GTI gap is accepted as a genuine "
                        "multi-pointing merge boundary only when it exceeds this many times "
                        "its own next-largest neighbour in the sorted gap list (an outlier "
                        "test on the sorted gaps, never against the population median/mean "
                        "of all gaps). Default 10.0. Only consulted when the split is "
                        "attempted at all -- see split_multiobs."
                    ),
                },
                "split_multiobs": {
                    "type": "boolean",
                    "description": (
                        "Force (true) or forbid (false) attempting the N-way GTI-gap split "
                        "at all. Omit (null/default) to auto-decide from the IXPE "
                        "literature-merged-ObsID naming convention: an obsid ending in "
                        "'99' (e.g. 03002599) is judged merged and the split criterion "
                        "(split_gap_factor AND split_gap_min_s AND split_gap_min_frac, "
                        "capped by max_subobs) is applied; any other obsid -- an ordinary "
                        "single monitoring pointing, however large its own GTI gaps happen "
                        "to be (e.g. 02008201) -- is treated as one whole-span "
                        "sub-observation and never split."
                    ),
                },
                "split_gap_min_s": {
                    "type": "number",
                    "description": (
                        "Absolute scale floor (seconds): a gap that clears split_gap_factor "
                        "must ALSO be at least this many seconds to become a split boundary. "
                        "Default 86400.0 (1 day) -- far above any routine IXPE "
                        "occultation/SAA visibility gap (~96.6 min orbit caps a routine gap "
                        "at a few ks), far below a genuine merged-pointing gap (days to "
                        "months). Only consulted when the split is attempted at all -- see "
                        "split_multiobs."
                    ),
                },
                "split_gap_min_frac": {
                    "type": "number",
                    "description": (
                        "Fractional scale floor: a gap that clears split_gap_factor must "
                        "ALSO be at least this fraction of the observation's own total "
                        "elapsed span (last GTI stop minus first GTI start) to become a "
                        "split boundary. Default 0.10 (10%). This is what separates 'two "
                        "distinct pointings merged' (a dominant fraction of the elapsed "
                        "span is empty) from 'one campaign with many day-scale visit gaps' "
                        "(each gap a small fraction of the whole span). Only consulted when "
                        "the split is attempted at all -- see split_multiobs."
                    ),
                },
                "max_subobs": {
                    "type": "integer",
                    "description": (
                        "Hard cap on how many sub-observations the N-way split can ever "
                        "produce. After split_gap_factor AND both scale floors, only the "
                        "(max_subobs - 1) LARGEST surviving boundary gaps are kept (ties "
                        "broken by earlier GTI index), so at most max_subobs sub-obs are "
                        "emitted. Default 4. Omit (null) to disable the cap; <= 1 forbids "
                        "any split outright. Only consulted when the split is attempted at "
                        "all -- see split_multiobs."
                    ),
                },
                "max_return_subobs": {
                    "type": "integer",
                    "description": (
                        "Bounds the RETURNED 'subobs_keys'/'per_subobs' payload to at most "
                        "this many sub-observations (chronological order) -- a caller-side "
                        "transport limit, distinct from max_subobs (which bounds how many "
                        "sub-observations are produced at all). Default 25. Every sub-obs "
                        "CSV is still written and catalogued regardless; when the true "
                        "count exceeds this, 'per_subobs_truncated'/'n_subobs_total'/"
                        "'subobs_listing_path' are added to the return value -- see RETURNS."
                    ),
                },
            },
            "required": ["obs_db_path", "obsid"],
        },
    },

    # ── Model-independent q(t)-u(t) discrete correlation function ──────────
    "research_toolkit.xspec_utils.ixpe_qu_lag_dcf": {
        "module": "xspec_utils.lib.qu_lag_dcf",
        "fn": "qu_lag_dcf",
        "description": (
            "Model-independent Edelson & Krolik (1988) discrete correlation function (DCF) "
            "between q(t)=Q/I and u(t)=U/I for every sub-observation of a set of local IXPE "
            "ObsIDs, with permutation significance and bootstrap lag intervals. Data only -- "
            "pair it with a plotter tool for figures."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "PREREQUISITE: reads each ObsID's existing qu_slices_<subobs_key>.csv series "
                "(globbing <obs_db_path>/<obsid>/qu_slices_*.csv, one file per subobs), so run "
                "ixpe_qu_slice_pcubes first. This tool computes statistics on those series and "
                "reduces no data itself, so it needs neither XSPEC nor ixpeobssim and runs "
                "fast.\n"
                "EVERY REQUESTED SUBOBS IS RECORDED: an obsid with no qu_slices files on disk "
                "gets status='no_qu_slices' (subobs_key=obsid), surfacing a failed or "
                "not-yet-run slicing step in place of a plausible-looking partial result; a "
                "qu_slices file that fails the provenance-header check gets "
                "status='unreadable_qu_slices'; a subobs with fewer than 'n_min_slices' rows "
                "gets status='insufficient_slices'. All three are excluded from every "
                "statistic and from the population correlations, while still getting one row "
                "in ccf_master.csv and one per_subobs entry (csv_path=None), so each one stays "
                "visible.\n"
                "DCF: mean-subtracts q/u, forms UDCF_ij=a_i*b_j/sqrt(D) over every (i,j) pair "
                "(dt_ij=t_mean_j-t_mean_i in MET seconds; D is the measurement-noise-corrected "
                "(var(q)-mean(q_err^2))*(var(u)-mean(u_err^2)), falling back to the "
                "uncorrected var(q)*var(u) with a dcf_uncorrected_denominator flag when either "
                "corrected variance drops to 'variance_floor' or below), then bins UDCF by dt "
                "at width max(lag_bin_s or median_dt, median_dt) over "
                "+/-'max_lag_frac'*span. Binning the irregular pair set directly is what makes "
                "this valid for a slice series with non-uniform sampling, which a uniform-grid "
                "correlate would misrepresent. Records tau_peak_s/dcf_peak (the bin of max "
                "|DCF|), dcf_zero (the bin centred on lag 0, NaN when it holds no pairs) and "
                "n_pairs_peak, and cross-checks against a masked-uniform-grid "
                "cross-correlation (resampled onto a median-spacing grid, nearest sample "
                "within half a bin, gaps masked) as tau_peak_uniform_s.\n"
                "SIGNIFICANCE: 'n_permutations' time-shuffles of u(t) give an empirical perm_p "
                "on max |DCF| (add-one smoothed) plus a per-bin 1st/99th-percentile envelope "
                "(perm_lo99/perm_hi99 for the plotter's band); 'n_bootstrap' draws of "
                "q_i+N(0,q_err_i) / u_i+N(0,u_err_i) give a 68% interval on the recovered "
                "tau_peak/dcf_peak (tau_boot_lo/hi, dcf_peak_boot_lo/hi). coherent_rotation = "
                "(perm_p < 'coherent_p_threshold') AND (the bootstrap tau_peak 68% interval "
                "excludes 0).\n"
                "REPRODUCIBILITY: all random draws (permutation shuffles, bootstrap noise) "
                "come from one numpy.random.default_rng('seed'), so two calls with the same "
                "seed over the same inputs reproduce identically.\n"
                "POPULATION COVARIATES: reconstructs the two covariates the qu_slices CSV does "
                "not carry, from that CSV rather than by re-slicing every ObsID -- r_proxy = "
                "sum(counts)/(sum(good_exposure_s)*pi*src_rad^2), with src_rad taken from the "
                "qu_slices provenance header, and pd_significance = an "
                "inverse-variance-weighted combined normalized-Stokes PD/sigma_PD from the "
                "per-slice q/u (a weighted combination of existing values, not an event-level "
                "reduction). Runs a cross-observation Spearman correlation of tau_peak_s "
                "against each covariate over status=='ok' subobs, dropping noise_limited ones "
                "when 'exclude_noise_limited' is true (reconstructed as n_rows < the qu_slices "
                "provenance's own n_min -- a different threshold from this tool's "
                "'n_min_slices').\n"
                "OUTPUT: one provenance-headed qu_dcf_<subobs_key>.csv per "
                "successfully-processed subobs (lag_s, dcf, n_pairs, perm_lo99, perm_hi99) "
                "under <obs_db_path>/<obsid>/, recorded at level 3 (kind='qu_dcf') via "
                "obsdb_catalog.neutral.advance; plus ccf_master.csv (one row per subobs) and "
                "population_correlations.csv under 'output_dir'.\n"
                "RETURNS: 'ccf_master_csv', 'population_csv', 'per_subobs' (subobs_key -> a "
                "bounded summary: 'status' in 'ok'/'no_qu_slices'/'unreadable_qu_slices'/"
                "'insufficient_slices', 'n_slices', 'tau_peak_s', 'dcf_peak', 'perm_p', "
                "'coherent_rotation', and 'csv_path' or None for any non-'ok' status), "
                "'population_correlations' (over status=='ok' rows only), 'stdout', "
                "'exit_code' (always 0 -- a per-obsid/per-subobs problem is recorded rather "
                "than raised)."
            ),
            "properties": {
                "obs_db_path": {
                    "type": "string",
                    "description": (
                        "The IXPE mission's obs-db directory (e.g. "
                        "'/data/fzhou/XSPEC/Cyg_X-1/IXPE'). Each ObsID's subobs series "
                        "are read from <obs_db_path>/<obsid>/qu_slices_*.csv; each "
                        "subobs's qu_dcf_<subobs_key>.csv lands alongside it."
                    ),
                },
                "obsids": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "IXPE ObsIDs to process (each must already carry qu_slices CSV(s)).",
                },
                "output_dir": {
                    "type": "string",
                    "description": (
                        "Directory ccf_master.csv and population_correlations.csv are "
                        "written to (created if absent)."
                    ),
                },
                "sources": {
                    "type": "object",
                    "additionalProperties": {"type": "string"},
                    "description": (
                        "Optional {obsid: source_name} label map for ccf_master.csv. "
                        "Falls back to the obsid itself when omitted or an obsid is "
                        "absent from it."
                    ),
                },
                "lag_bin_s": {
                    "type": "number",
                    "description": (
                        "DCF bin width in seconds. Omit (null/default) to use each "
                        "subobs's own median inter-slice spacing; an explicit value is "
                        "still floored at that median -- never finer than the data "
                        "actually resolves."
                    ),
                },
                "n_min_slices": {
                    "type": "integer",
                    "description": (
                        "Minimum slice-row count a subobs needs before its DCF is "
                        "computed at all; fewer gets status='insufficient_slices' and is "
                        "excluded. Default 8. This is this tool's OWN threshold, "
                        "distinct from the qu_slices provenance's own n_min (used only "
                        "to reconstruct noise_limited)."
                    ),
                },
                "n_permutations": {
                    "type": "integer",
                    "description": "Time-shuffle draw count for the permutation-null significance test. Default 1000.",
                },
                "n_bootstrap": {
                    "type": "integer",
                    "description": "Noise-draw count for the bootstrap peak-lag/peak-value uncertainty. Default 1000.",
                },
                "variance_floor": {
                    "type": "number",
                    "description": (
                        "Threshold below which a measurement-noise-corrected variance is "
                        "treated as degenerate, triggering the uncorrected-denominator "
                        "fallback. Default 1e-12."
                    ),
                },
                "coherent_p_threshold": {
                    "type": "number",
                    "description": (
                        "perm_p below this (together with a zero-excluding bootstrap "
                        "tau_peak 68% interval) flags coherent_rotation. Default 0.01."
                    ),
                },
                "max_lag_frac": {
                    "type": "number",
                    "description": (
                        "DCF lag range as a fraction of the subobs's own time span "
                        "(+/-max_lag_frac*span). Default 0.5."
                    ),
                },
                "exclude_noise_limited": {
                    "type": "boolean",
                    "description": (
                        "Drop noise_limited subobs from the population Spearman "
                        "correlations (they still appear in ccf_master.csv). Default true."
                    ),
                },
                "seed": {
                    "type": "integer",
                    "description": (
                        "Seed for every permutation/bootstrap random draw -- drives exact "
                        "reproducibility across repeated calls. Default 0."
                    ),
                },
            },
            "required": ["obs_db_path", "obsids", "output_dir"],
        },
    },

    # ── qu(t)-u(t) DCF reviewer figures (companion plotter) ─────────────────
    "research_toolkit.xspec_utils.ixpe_qu_lag_dcf_plot": {
        "module": "xspec_utils.lib.qu_lag_dcf",
        "fn": "qu_lag_dcf_plot",
        "description": (
            "Renders the three standard reviewer figures for an ixpe_qu_lag_dcf run: a grid of "
            "per-subobs DCF panels and two population tau-vs-covariate scatters."
        ),
        "args_schema": {
            "type": "object",
            "description": (
                "PREREQUISITE: run ixpe_qu_lag_dcf first for 'ccf_master_csv'; its sibling "
                "population_correlations.csv in the same directory is read automatically. Pure "
                "matplotlib (Agg backend) importing neither XSPEC, ixpeobssim nor red_ixpe, so "
                "it runs fast without the science stack.\n"
                "FIGURES: writes exactly three files, appending 'dcf_panels', "
                "'population_tau_r_proxy', 'population_tau_pd_signif' (plus '.png') to "
                "'save_path_prefix'. (1) dcf_panels -- a grid of up to 'max_panels' "
                "status=='ok' subobs in ccf_master.csv row order, each panel reading its own "
                "<obs_db_path>/<obsid>/qu_dcf_<subobs_key>.csv for DCF vs lag, with the "
                "permutation 99% band (perm_lo99/perm_hi99) shaded, the tau_peak_s/dcf_peak "
                "marker, and the bootstrap 68% interval drawn as the literal segment between "
                "its own [tau_boot_lo, tau_boot_hi] endpoints read straight from the "
                "ccf_master row -- reproduced faithfully even when that interval does not "
                "bracket tau_peak_s, so an asymmetric or offset interval stays visible as "
                "measured. Each overlay is drawn once both its values are finite. (2) "
                "population_tau_r_proxy -- tau_peak_s (its bootstrap 68% interval as a "
                "vertical bar) vs r_proxy over EVERY status=='ok' subobs, unbounded by "
                "max_panels, annotated with the Spearman rho/p read directly from "
                "population_correlations.csv's tau_peak_vs_r_proxy row. That statistic is "
                "reused rather than recomputed because ixpe_qu_lag_dcf drops noise_limited "
                "subobs from the correlation by default, and noise_limited is not a "
                "ccf_master column this plotter could reconstruct. (3) "
                "population_tau_pd_signif -- the same for pd_significance and the "
                "tau_peak_vs_pd_significance row.\n"
                "UNITS: tau_peak_s, lag_s, tau_boot_lo and tau_boot_hi are MET seconds; "
                "r_proxy is cps/arcmin^2; pd_significance is dimensionless.\n"
                "INPUT ERRORS (checked before any figure is drawn): a missing or wrong-schema "
                "ccf_master_csv, zero status=='ok' rows, or a missing/wrong-schema sibling "
                "population_correlations.csv gives exit_code=1 with stderr set and no files "
                "written. A status=='ok' row whose own qu_dcf_<subobs_key>.csv is absent or "
                "unreadable raises ValueError in place of an empty panel, since an 'ok' "
                "ccf_master row promises that file exists.\n"
                "OUTPUT LOCATION: 'output_dir' is created as needed; omitting "
                "'save_path_prefix' builds an auto-generated timestamped prefix under it.\n"
                "RETURNS: 'figures' (3 FigureRef dicts), 'plotted' (dcf_panels: subobs_key -> "
                "lag_s/dcf/perm_lo99/perm_hi99/tau_peak_s/dcf_peak/tau_boot_lo/tau_boot_hi; "
                "each population panel: subobs_key/x/y/tau_boot_lo/tau_boot_hi/spearman "
                "{n, rho, p}), 'stdout', 'exit_code'."
            ),
            "properties": {
                "obs_db_path": {
                    "type": "string",
                    "description": (
                        "The IXPE mission's obs-db directory (e.g. "
                        "'/data/fzhou/XSPEC/Cyg_X-1/IXPE'). Each status=='ok' subobs's own "
                        "qu_dcf_<subobs_key>.csv is read from "
                        "<obs_db_path>/<obsid>/qu_dcf_<subobs_key>.csv."
                    ),
                },
                "ccf_master_csv": {
                    "type": "string",
                    "description": (
                        "Path to ccf_master.csv written by ixpe_qu_lag_dcf. Its sibling "
                        "population_correlations.csv is located automatically at "
                        "os.path.dirname(os.path.abspath(ccf_master_csv))/"
                        "population_correlations.csv -- both must be present."
                    ),
                },
                "output_dir": {
                    "type": "string",
                    "description": (
                        "Directory created for this call's output; the auto-generated "
                        "save_path_prefix (when one is not supplied) is built under it."
                    ),
                },
                "save_path_prefix": {
                    "type": "string",
                    "description": (
                        "Base path; 'dcf_panels', 'population_tau_r_proxy', "
                        "'population_tau_pd_signif' plus '.png' are appended per figure. "
                        "Auto-generated timestamped prefix under output_dir if omitted."
                    ),
                },
                "max_panels": {
                    "type": "integer",
                    "description": (
                        "Maximum status=='ok' subobs drawn in the dcf_panels grid (ccf_master "
                        "row order). Default 12. Never bounds the two population panels, "
                        "which always cover every status=='ok' subobs."
                    ),
                },
            },
            "required": ["obs_db_path", "ccf_master_csv", "output_dir"],
        },
    },
}
