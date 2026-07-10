# xspec_utils

A Python/Bash toolkit for joint spectro-polarimetric analysis of X-ray binary data from **IXPE**, **NICER**, and **NuSTAR**, built on top of [HEASoft/XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/) and [ixpeobssim](https://ixpeobssim.readthedocs.io/). Developed for the Cyg X-1 campaign but designed to be reusable for other bright X-ray binaries. The current instrument set (IXPE + NICER + NuSTAR) is a starting point — the loader, reduction scripts, and epoch definitions are all designed to be extended modularly to additional missions (e.g. Chandra, XMM, Swift) by adding a new file locator module and a matching branch in `xspec_help.load_data`.

## Features

- **Batch data reduction** shell scripts for IXPE, NICER, and NuSTAR pipelines
- **Automated point-source detection** via image centroiding (`find_src.py`) with DS9 region output
- **Unified data loader** that assigns XSPEC data groups across instruments in a single call (`xspec_help.load_data`), extensible to new instruments
- **Stokes I/Q/U spectral extraction** for IXPE polarization cubes and PHA spectra
- **Time-resolved IXPE polarimetry** (`red_ixpe.pcube_timeseg`) for intra-observation PD/PA/Stokes/MDP99 stability, DU1/2/3 merged, exposed as the `ixpe_pcube_timeseg` AI tool
- **Energy-resolved / multi-ObsID-combined IXPE polarimetry** (`red_ixpe.pcube`, extended): a caller-supplied energy `ebinning` per band, and `obs` as either one ObsID or a list of ObsIDs natively combined into a single `xpselect`/`xpbin PCUBE` call (raw Stokes I/Q/U/COUNTS summed across every DU *and* ObsID, not a post-hoc PD/PA merge) — exposed as the `ixpe_pcube` AI tool
- **Energy-resolved IXPE polarimetry figures** (`ixpe_pcube_plot.ixpe_pcube_plot`): the PD(E)/PA(E) + Q-U-plane reviewer figure (Krawczynski+2022 Fig. S5 style) from an `ixpe_pcube` CSV, with an optional `fit_tbabs_pollin_powerlaw` PD(E)/PA(E) model overlay (falling back to an internal `scipy` weighted-linear fit of the PCUBE points when no `pollin` quartet is supplied) — pure matplotlib, no XSPEC/ixpeobssim, exposed as the `ixpe_pcube_plot` AI tool
- **Per-DU IXPE PHA1 reduction, reuse-or-reduce** (`red_ixpe.reduce_pha`): reuses an ObsID's existing per-DU Stokes-I PHA1 src+bkg products when already on disk, otherwise centroids the source position from the Level-2 event file via `find_src.run_extraction` (writing DS9 `src.reg`/`bkg.reg`) and runs `red_ixpe.pha`; exposed as the `ixpe_reduce_pha` AI tool -- the input `red_ixpe.HR` (hardness ratio) and `xspec_help.load_data` need, kept as its own tool rather than folded into an HR tool since HR reads PHA
- **Net Stokes-I hardness ratio survey** (`red_ixpe.HR`): community-normalized `HR=(H-S)/(H+S)` read directly from each ObsID's existing per-DU PHA1 products (no XSPEC `AllData` loading), attempting DU1/DU2/DU3 independently so a genuinely-absent DU is dropped gracefully while a present-but-unreduced DU is never silently skipped; closed-form Gaussian error propagation from PHA `RATE`/`STAT_ERR`; soft/hard bands resolved per PHA from its own RMF EBOUNDS (falling back to an installed CALDB RMF for the same DU when a PHA's `RESPFILE` predates an environment migration); exposed as the `ixpe_HR` AI tool, which surveys a list of ObsIDs and writes `bhxrb_hardness_ratios.csv`
- **PA/PD constancy analysis** (`pa_pd_stability.pa_pd_stability`) on an existing `pcube_timeseg` CSV: weighted chi-squared constancy tests (PD/PA/QI/UI), the authoritative 2D Stokes variability test, moment-estimator intrinsic scatter, peak-to-peak ranges, structure functions, and cross-source Spearman correlations — pure numpy/scipy, no XSPEC/ixpeobssim, exposed as the `ixpe_pa_pd_stability` AI tool
- **PA/PD stability diagnostic figures** (`pa_pd_stability.stability_plot`): the four standard reviewer figures (PA/PD time series, Stokes-plane Q/I vs U/I evolution, reduced-chi2 histograms, structure functions) from an `ixpe_pa_pd_stability` catalog — pure matplotlib, no XSPEC/ixpeobssim, exposed as the `ixpe_pa_pd_stability_plot` AI tool
- **PD-vs-HR correlation survey** (`bhxrb_pd_hr_survey.pd_hr_survey`): correlates PD (an `ixpe_pcube_timeseg` integrated row) against HR (an `ixpe_HR` CSV) across a BH-XRB ObsID survey — detection classification (`PD>MDP99`), per-source/combined Spearman, a caller-defined soft-vs-hard-state KS test, and a detection-bias Mann-Whitney U test — pure numpy/scipy, no XSPEC/ixpeobssim, exposed as the `bhxrb_pd_hr_survey` AI tool
- **PD-vs-HR survey diagnostic figures** (`bhxrb_pd_hr_survey.pd_hr_survey_plot`): the five standard reviewer figures (combined and per-source PD-vs-HR scatter, soft-vs-hard PD distribution comparison, HR detection-bias histogram, MDP99-vs-exposure feasibility scatter) from a `bhxrb_pd_hr_survey` CSV — pure matplotlib, no XSPEC/ixpeobssim, exposed as the `bhxrb_pd_hr_survey_plot` AI tool
- **Joint fitting helpers**: multi-panel plotting with model component decomposition, fit-state freezing for snapshotting
- **polconst spectropolarimetric cross-check** (`xspec_help.fit_tbabs_polconst_powerlaw`): fits `tbabs*polconst*powerlaw` to every IXPE Stokes I/Q/U triplet in a `load_data` indices dict, tying polarization degree/angle into one free pair for the whole session -- the polconst XSPEC cross-check of the standard PCUBE polarization (Krawczynski+2022 Cyg X-1), exposed as the `fit_tbabs_polconst_powerlaw` AI tool
- **pollin linear energy-dependence fit** (`xspec_help.fit_tbabs_pollin_powerlaw`): fits `tbabs*pollin*powerlaw` to every IXPE Stokes I/Q/U triplet in a `load_data` indices dict, tying `pollin`'s linear PD(E)/PA(E) (a constant plus an energy-slope term each) into one free quartet for the whole session -- the paper's actual linear energy-dependence route for Cyg X-1 polarization, replacing a non-toolkit-grounded `scipy` fit on `ixpe_pcube`'s PD(E) points, exposed as the `fit_tbabs_pollin_powerlaw` AI tool
- **This-work-vs-paper benchmark comparison** (`benchmark_compare.compare_to_benchmark`): a general, record-based comparison of this-work measured quantities against cited paper values -- delta, z-score significance (`z=abs(delta)/sqrt(err_this**2+err_paper**2)`, `z=inf` when both errors are zero), and a `>sigma_flag` flag list -- written to a provenance-headed CSV; hardcodes no paper numbers (every record, including the mandatory `paper_source` citation, comes from the caller). Two lib-only extractor helpers, `pcube_quantity`/`fit_quantity`, keep the fraction<->percent unit choice out of the caller by reading straight from an `ixpe_pcube` CSV or a `fit_tbabs_polconst_powerlaw`/`fit_tbabs_pollin_powerlaw` result dict. Pure stdlib -- no XSPEC/ixpeobssim import -- exposed as the `compare_to_benchmark` AI tool
- **Fully customisable epoch management** for multi-observation campaigns, applicable to any source

## Repository layout

```
lib/                        # Python library modules
  xspec_help.py             # Core: data loading, fitting, inspection, rebinning
  red_ixpe.py               # IXPE reduction: PCUBE, PHA1/Q/U extraction, MJD, hardness ratio
  pa_pd_stability.py        # PA/PD/Stokes constancy statistics + stability_plot diagnostic figures (pure numpy/scipy/matplotlib)
  bhxrb_pd_hr_survey.py     # PD-vs-HR correlation/statistics survey + pd_hr_survey_plot diagnostic figures (pure numpy/scipy/matplotlib)
  ixpe_pcube_plot.py        # PD(E)/PA(E) + Q-U-plane figure for an ixpe_pcube CSV (pure matplotlib)
  benchmark_compare.py      # This-work-vs-paper delta/z-score/flag comparison + pcube_quantity/fit_quantity unit extractors (pure stdlib)
  NICER.py                  # NICER file locators (PHA, RMF, ARF)
  NuSTAR.py                 # NuSTAR file locators (PHA grouped, RMF, ARF)
  find_src.py               # Automated source centroiding and region file generation
  src_pos.py                # Catalogue of IXPE source coordinates (RA, Dec)
  epoch.py                  # Epoch definitions: obs-ID groupings and spectral states

IXPE/
  filter_inst_bkg_IXPE.sh   # Batch unzip + instrumental background rejection

NICER/
  nicer_red.sh              # Batch nicerl3-spect reduction (3C50 background)

NuSTAR/
  nustar_red.sh             # Batch nuproducts (FPMA/FPMB) + source detection
  nustar_grp.sh             # Batch optimal grouping via ftgrouppha

Joint_data_ana/
  spec_decomp.ipynb         # Example: joint spectro-polarimetric fitting across epochs
```

## Prerequisites

| Dependency | Version tested | Purpose |
|---|---|---|
| [HEASoft](https://heasarc.gsfc.nasa.gov/lheasoft/) | 6.34 | XSPEC, `nicerl2`/`nicerl3-spect`, `nuproducts`, `ftgrouppha` |
| [CALDB](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/) | current | Calibration files for NICER, NuSTAR, IXPE |
| Python | >= 3.10 | Type hint syntax (`X | None`, `tuple[...]`) |
| [ixpeobssim](https://ixpeobssim.readthedocs.io/) | latest | IXPE pipeline (`xpselect`, `xpbin`), IRFs, polarization cubes |
| [PyXSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/) | (bundled with HEASoft) | Scripted spectral fitting |
| [relxill](http://www.sternwarte.uni-erlangen.de/~dauser/research/relxill/) | latest | Relativistic reflection model (loaded via `AllModels.lmod`) |
| astropy, numpy, scipy, matplotlib | standard | FITS I/O, numerics, plotting |

### Environment setup

Source the HEASoft and CALDB environment before launching Python/Jupyter:

```bash
# HEASoft
export HEADAS=/opt/heasoft/heasoft-6.34/x86_64-pc-linux-gnu-libc2.34
. "$HEADAS/headas-init.sh"

# Conda (or your preferred Python env)
source /opt/Anaconda/anaconda3/bin/activate
conda activate astropy_env

# CALDB
export CALDB=/path/to/CALDB
export CALDBCONFIG="$CALDB/software/tools/caldb.config"
export CALDBALIAS="$CALDB/software/tools/alias_config.fits"
export GEOMAG_PATH="$CALDB/data/gen/pcf/geomag"
```

## Quick start

### 1. Data reduction (shell)

Run from within each instrument's data directory:

```bash
# IXPE: unzip + instrumental background filtering
cd IXPE/ && nohup ./filter_inst_bkg_IXPE.sh > filter_inst_bkg_IXPE.out 2>&1 &

# NICER: Level 3 spectral products
cd NICER/ && nohup ./nicer_red.sh > nicer_red.out 2>&1 &

# NuSTAR: source detection + nuproducts + optimal grouping
cd NuSTAR/ && nohup ./nustar_red.sh > nustar_red.out 2>&1 &
cd NuSTAR/ && nohup ./nustar_grp.sh > nustar_grp.out 2>&1 &
```

> **Important — do not parallelise.** The reduction shell scripts loop over observations sequentially by design. Running multiple obs-IDs of the same instrument in parallel causes race conditions and CALDB cache corruption. Launch one script per instrument at a time, and let each complete before starting another pass.
>
> **CALDB must be reliable and local.** These pipelines read calibration files intensively; remote / FTP-mounted CALDB is very likely to break the run mid-way. Keep CALDB on local disk and do **not** schedule cron jobs that update or sync CALDB while a reduction is in flight.

### 2. Joint spectral fitting (Python / Jupyter)

```python
import sys
sys.path.append("/path/to/lib")

from xspec import AllModels, Fit, Xset
AllModels.lmod("relxill", "/path/to/relxill/")
Fit.query = "no"
Xset.allowPrompting = False

import xspec_help, epoch
from src_pos import SRC_COORD

# Load one epoch (NuSTAR + NICER + IXPE Stokes I/Q/U)
idx = xspec_help.load_data(
    epoch=epoch.epoch["E5"],
    ixpe_path="/path/to/IXPE",
    nicer_path="/path/to/NICER",
    nustar_path="/path/to/NuSTAR",
)

# Fit a model
xspec_help.fit_tbabs_polconst_diskbb_polconst_relxillCp(
    indices=idx,
    log_path="fit_E5.log",
)

# Inspect: broadband spectrum + per-DU Stokes decomposition
fig, ax, figs = xspec_help.inspect_epoch(
    indices=idx,
    e_bins=[2., 3., 4., 5., 6., 8.],
    title="Epoch E5",
)
```

### 3. Source detection (CLI)

```bash
python3 lib/find_src.py event_file.evt \
    --src_rad 1.0 --bkg_inner_rad 2.0 --bkg_outer_rad 3.0
# Outputs: src.reg, bkg.reg, verify.png in the event file directory
```

## API reference

---

### `xspec_help` — core fitting and visualisation

---

#### `xspec_help.load_data(epoch, ixpe_path, ...)`

Loads all instruments for one epoch into XSPEC and returns a nested index dict.

```python
def load_data(
    epoch: dict,
    ixpe_path: str,
    nicer_path: str | None = None,
    nustar_path: str | None = None,
    ixpe_changes: bool = False,
    ixpe_src_rad: float = 1.0,
    ixpe_bkg_inner_rad: float = 3.0,
    ixpe_bkg_outer_rad: float = 5.0,
    ixpe_irf_date: str | None = None,
    ixpe_erange: tuple[float, float] = (2.0, 8.0),
    nicer_erange: tuple[float, float] = (1.0, 10.0),
    nustar_erange: tuple[float, float] = (4.0, 50.0),
    nustar_use_bkg: bool = False,
    clear_existing: bool = True,
    chatter: bool = False,
) -> dict
```

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `epoch` | `dict` | — | Single-epoch sub-dict from `epoch.epoch`, e.g. `epoch.epoch['E5']`. Expected keys: `'IXPE'`, `'NICER'`, `'NuSTAR'` (each a list of obs-ID strings) and `'STATE'`. Instruments absent from the dict are silently skipped. `epoch['IXPE']` may list any number of ObsIDs, exactly like the `'NICER'`/`'NuSTAR'` lists — each ObsID independently resolves its own centroid, response files/`BACKSCAL`, and IRF epoch (see **IXPE Stokes gain linking** below).|
| `ixpe_path` | `str` | — | Path to the directory containing IXPE observation sub-directories. |
| `nicer_path` | `str \| None` | `None` | Path to the NICER data directory. Required when `'NICER'` is in `epoch`. |
| `nustar_path` | `str \| None` | `None` | Path to the NuSTAR data directory. Required when `'NuSTAR'` is in `epoch`. |
| `ixpe_changes` | `bool` | `False` | Force a fresh IXPE reduction if `True`: the source position is centroided from the DU1 Level-2 event file via `find_src.run_extraction` (mirroring `red_ixpe.reduce_pha` -- never `src_pos.SRC_COORD`) and `red_ixpe.pha` is re-run at that position and `ixpe_src_rad`/`ixpe_bkg_inner_rad`/`ixpe_bkg_outer_rad`, so a stale on-disk PHA1 (wrong region, dead response paths) is never silently reused. If `False` (default), reuse products already on disk as-is, no centroiding. |
| `ixpe_src_rad` | `float` | `1.0` | IXPE source-region radius in arcminutes. |
| `ixpe_bkg_inner_rad` | `float` | `3.0` | IXPE background annulus inner radius in arcminutes. |
| `ixpe_bkg_outer_rad` | `float` | `5.0` | IXPE background annulus outer radius in arcminutes. |
| `ixpe_irf_date` | `str \| None` | `None` | IRF calibration epoch (`'YYYYMMDD'`) for the IXPE response. `None` auto-selects the installed epoch nearest the observation date (see **IRF calibration epoch**); pass a value only to override. Used only when `ixpe_changes=True`. |
| `ixpe_erange` | `tuple[float, float]` | `(2.0, 8.0)` | Energy range in keV to *notice* for IXPE spectra. |
| `nicer_erange` | `tuple[float, float]` | `(1.0, 10.0)` | Energy range in keV to *notice* for NICER spectra. |
| `nustar_erange` | `tuple[float, float]` | `(4.0, 50.0)` | Energy range in keV to *notice* for NuSTAR spectra. |
| `nustar_use_bkg` | `bool` | `False` | Load NuSTAR background spectrum. Disabled by default because clean background regions are often unavailable for bright sources like Cyg X-1. |
| `clear_existing` | `bool` | `True` | Call `AllData.clear()` and `AllModels.clear()` before loading. |
| `chatter` | `bool` | `False` | Print file paths and per-instrument loading summaries. |

**Data-group assignment**

Data groups are assigned **first-come-first-served with respect to the epoch definition**: the loader iterates over the keys of the epoch dict in insertion order, and the first instrument encountered starts at group 1. The `epoch.py` convention lists `NuSTAR → NICER → IXPE`. Reorder the keys in `epoch.py` and the data-group numbering changes accordingly. Each NuSTAR FPM (FPMA and FPMB) occupies its own data group, so one NuSTAR observation contributes two groups and two observations contribute four groups; each NICER observation occupies one group; each IXPE detector unit gets one group.

Within each instrument, ObsIDs are processed in `epoch[<inst>]` list order via an inner loop (`for obs in epoch[<inst>]`) — the NuSTAR, NICER, and IXPE branches all share this shape. For IXPE, each ObsID contributes 3 consecutive data groups (one per DU), so `epoch['IXPE'] = [obsA, obsB]` yields `3*2 = 6` IXPE data groups, consecutive and gapless, before the next instrument's groups begin.

**IXPE Stokes gain linking**

For each IXPE DU the three spectra (Stokes I, Q, U) are loaded into consecutive XSPEC indices within the same data group. The gain slope of Q and U are linked to the I gain slope of **that same (ObsID, DU)** (`r.gain.slope.link = AllData(I_idx).response.gain.slope`, where `I_idx` is that triplet's own running spectrum index) to enforce a shared per-DU energy-scale calibration. The gain offset is frozen to zero for all three. This tying is scoped strictly to one (ObsID, DU) at a time — a later ObsID's Q/U in a multi-ObsID epoch never links back to an earlier ObsID's I, and no gain parameter is ever tied *across* ObsIDs or DUs (the same convention the NuSTAR/NICER branches implicitly follow, since they carry no gain-linking at all).

**Returns**

A uniform nested dict keyed by `instrument → obsid → detector → list of XSPEC spectrum indices`. Only instruments present in `epoch` appear. Data groups are not returned explicitly; recover a spectrum's data group via `AllData(spec_idx).dataGroup`.

```python
# example: 2 NuSTAR obs (FPMA+FPMB each), 2 NICER obs, 1 IXPE obs
{
    'NuSTAR': {
        '<nustar_obsid_1>': {'FPMA': [1], 'FPMB': [2]},
        '<nustar_obsid_2>': {'FPMA': [3], 'FPMB': [4]},
    },
    'NICER': {
        '<nicer_obsid_1>': {'mpu7': [5]},
        '<nicer_obsid_2>': {'mpu7': [6]},
    },
    'IXPE': {
        '<ixpe_obsid>': {
            'DU1': [7,  8,  9],   # I, Q, U
            'DU2': [10, 11, 12],
            'DU3': [13, 14, 15],
        },
    },
}
```

Per-detector lists are always ordered by spectrum index; IXPE DU lists are ordered `[I, Q, U]`. NuSTAR FPM labels follow sorted filename order. A multi-ObsID `epoch['IXPE']` (e.g. `['01002901', '01250101']`, the Krawczynski+2022 May+June validation pair) simply adds another top-level key under `indices['IXPE']`, each with its own 3 DU sub-dicts — `indices['IXPE']['01002901']['DU1']`, `indices['IXPE']['01250101']['DU1']`, etc.

**Allowed keys**

| Level | Instrument | Allowed keys |
|---|---|---|
| `instrument` | — | `'NuSTAR'`, `'NICER'`, `'IXPE'` |
| `detector`   | `NuSTAR` | `'FPMA'`, `'FPMB'` |
| `detector`   | `NICER`  | `'mpu7'` |
| `detector`   | `IXPE`   | `'DU1'`, `'DU2'`, `'DU3'` |

The NICER `'mpu7'` key is not a physical detector — it is the conventional label for the co-added, MPU-averaged spectrum produced by the standard NICER pipeline (`nicerl3-spect`) after good-detector screening. It is the **only** allowed detector key under NICER, so `indices['NICER'][obs]['mpu7']` is always a single-element list of one XSPEC spectrum index. The key is kept to preserve the uniform `instrument → obsid → detector → [spec_idx, ...]` shape across all instruments.

---

#### `xspec_help.fit_tbabs_polconst_diskbb_polconst_relxillCp(indices, ...)`

Defines and fits the model `constant * tbabs * (polconst*diskbb + polconst*relxillCp)` with physically motivated initial values for Cyg X-1.

> **Scope warning.** This helper is **specifically tuned for Cyg X-1 in its high/soft state (HSS)**. The initial values (nH, `diskbb` `Tin`/`norm`, `relxillCp` geometry, polarisation angles) encode HSS-appropriate priors and will not converge sensibly on other sources or other spectral states. Use it as a template, not as a black box: for new sources, states, or model compositions, prototype interactively first, validate the fit, and only then wrap the recipe into its own `fit_*` function following the same pattern.

```python
def fit_tbabs_polconst_diskbb_polconst_relxillCp(
    indices: dict,
    spec_pref: list[str] = ['NuSTAR', 'NICER', 'IXPE'],
    audit_data: bool = True,
    log_path: str | None = None,
    chatter_level: int = 10,
) -> None
```

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `indices` | `dict` | — | Return value of `load_data`. Used to look up data-group numbers. |
| `spec_pref` | `list[str]` | `['NuSTAR', 'NICER', 'IXPE']` | Instrument preference list for the reference cross-calibration constant. See below. |
| `audit_data` | `bool` | `True` | If `True`, print file paths and count rates for every loaded spectrum at the top of the log. |
| `log_path` | `str \| None` | `None` | If given, the verbose `AllModels.show` / `Fit.show` output is written to this file. If `None`, output goes to the notebook/terminal. |
| `chatter_level` | `int` | `10` | XSPEC chatter level used for the logged model/fit summary. Restored to 5 afterwards. |

**The `spec_pref` parameter in detail**

`spec_pref` controls which instrument's cross-calibration `constant` factor is designated as the reference (frozen to 1.0). All other data groups have their constants freed to float. This is the standard approach when jointly fitting multi-instrument spectra: one instrument anchors the flux scale and the others absorb residual calibration offsets.

The selection algorithm walks through `spec_pref` left-to-right and picks the **first entry that is a key in `indices`**, then uses the first obs / first detector / first spectrum under that instrument as the reference:

```python
for inst in spec_pref:
    if inst in indices:
        first_obs = next(iter(indices[inst]))
        first_det = next(iter(indices[inst][first_obs]))
        first_idx = indices[inst][first_obs][first_det][0]
        ref_grp  = AllData(first_idx).dataGroup
        break
```

Then every data group is iterated; the reference group gets `constant.frozen = True, value = 1.0` and all others get `constant.link = "", constant.frozen = False`.

**Practical examples**

- Default `['NuSTAR', 'NICER', 'IXPE']`: NuSTAR is preferred. If the epoch has NuSTAR, the constant on its first FPM is frozen to 1; every other NuSTAR FPM, NICER obs, and IXPE DU constant is free.
- IXPE-only epoch (e.g. E3, E4): NuSTAR and NICER are absent from `indices`, so the search falls through to IXPE and DU1's group becomes the reference.
- NuSTAR + IXPE epoch (no NICER, e.g. E6): NuSTAR is found first; the remaining NuSTAR FPMs and all IXPE DUs are freed relative to the first NuSTAR FPM.
- To override and use NICER as reference: pass `spec_pref=['NICER', 'NuSTAR', 'IXPE']`.

**Initial parameter values set by this function**

| Component | Parameter | Value | Frozen? |
|---|---|---|---|
| `tbabs` | `nH` | 0.45 | no |
| `polconst` (diskbb) | `psi` | 63.0° | no |
| `polconst` (relxillCp) | `psi` | −27.0° | yes |
| `polconst` (relxillCp) | `A` | 0.03 | no |
| `diskbb` | `Tin` | 0.43 keV | no |
| `diskbb` | `norm` | 1.3×10⁵ | no |
| `relxillCp` | `Incl` | 37° | yes |
| `relxillCp` | `a` | 0.998 | yes |
| `relxillCp` | `Rin` | −2.5 (RISCO units) | no |
| `relxillCp` | `Index1` | 3.0 | yes |
| `relxillCp` | `gamma` | 2.4 | no |
| `relxillCp` | `logxi` | 3.0 | no |
| `relxillCp` | `logN` | 15.0 | yes |
| `relxillCp` | `Afe` | 1.0 | yes |
| `relxillCp` | `kTe` | 300.0 keV | yes |
| `relxillCp` | `refl_frac` | 1.0 | no |

---

#### `xspec_help.fit_tbabs_polconst_powerlaw(indices, ...)`

Fits `tbabs*polconst*powerlaw` to every IXPE Stokes I/Q/U triplet in an indices dict, tying `polconst`'s polarization degree (`A`) and angle (`psi`) into **one free pair for the whole session** -- the polconst spectropolarimetric cross-check of the standard PCUBE polarization measurement (the Krawczynski+2022 Cyg X-1 XSPEC column).

```python
def fit_tbabs_polconst_powerlaw(
    indices: dict,
    gamma_init: float = 2.0,
    pd_init: float = 0.04,
    pa_init: float = -20.0,
    nh: float = 0.6,
    freeze_nh: bool = True,
    erange: tuple[float, float] = (2.0, 8.0),
    conf_level: float = 0.683,
    n_iterations: int = 100,
    log_path: str | None = None,
) -> dict
```

Call after `load_data(..., ixpe_changes=True, ixpe_src_rad=..., ixpe_bkg_inner_rad=..., ixpe_bkg_outer_rad=...)` -- a fresh, freshly-centroided reduction is required so each spectrum's own `RESPFILE`/`ANCRFILE` (arf for I, mrf for Q/U) and `BACKSCAL` resolve for the region actually used. Reusing a stale on-disk PHA1 (`ixpe_changes=False`) can carry dead response paths and a different region's `BACKSCAL` -- the root cause of a manual `tbabs*polconst*powerlaw` attempt (see Changelog) that returned `chi2/dof=64.1` and an unphysical `PD=0.85%`. Pass an IXPE-only `indices` dict: `Model()` here applies to every currently active data group, so a mixed epoch (NICER/NuSTAR also loaded) would incorrectly fold those groups into this fit too.

**Why no explicit tie/untie call is needed.** XSPEC ties every parameter of a model defined once across several already-active data groups to the first group's value by default (confirmed empirically: `AllModels(dg).polconst.A.link == '= p2'` for `dg > 1`). `fit_tbfeo_hss`'s cross-calibration constant explicitly *unties* per instrument (`constant.link = ""`) because that fit wants per-instrument freedom; this fit wants the opposite (one PD/PA for the whole session) and gets it by doing nothing extra.

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `indices` | `dict` | -- | Return value of `load_data`; must contain `'IXPE'`. Every obsid/DU I/Q/U triplet present is loaded into the joint fit -- one ObsID gives a single-ObsID fit, more than one a joint combined fit, always with exactly one tied PD/PA pair. |
| `gamma_init`, `pd_init`, `pa_init` | `float` | `2.0`, `0.04`, `-20.0` | Initial powerlaw photon index, `polconst` polarization fraction (0-1, not percent) and polarization angle (degrees) -- the physically-motivated Cyg X-1 starting point. |
| `nh` | `float` | `0.6` | `TBabs` equivalent hydrogen column, in units of 1e22 cm^-2 (Cyg X-1 ISM value). |
| `freeze_nh` | `bool` | `True` | Freeze `nH` at `nh`. Set `False` to let it float. |
| `erange` | `tuple[float, float]` | `(2.0, 8.0)` | keV range noticed on every IXPE spectrum before fitting. |
| `conf_level` | `float` | `0.683` | Confidence level the reported errors are scaled to (0.683 = 1 sigma, i.e. Delta-chi2=1.0 for one interesting parameter); converted via the chi-squared quantile for 1 degree of freedom (`scipy.stats.chi2.ppf`). Pass e.g. `0.90` for a 90% interval. |
| `n_iterations` | `int` | `100` | `Fit.nIterations` per convergence round (up to 5 rounds, stopping early once the fit statistic stabilises). |
| `log_path` | `str \| None` | `None` | If given, the verbose `AllModels.show()`/`Fit.show()` summary is logged here (chatter 10, restored to 5 afterwards). |

**Returns** a dict with `PD`, `PA` (fraction, degrees), `PD_err`, `PA_err`, `Gamma`, `Gamma_err`, `nH`, `chi2`, `dof`, `stat_per_dof`, `conf_level`, `indices` (the input, unchanged, for downstream `plot_fit`), `stdout`, `exit_code` -- every fitted number the caller needs, since XSPEC's own terminal output is muted in the sandbox.

**Errors are the covariance-matrix parameter-table 1-sigma.** `PD_err`/`PA_err`/`Gamma_err` come from `Parameter.sigma` (populated by `Fit.perform()`), scaled by `sqrt(chi2.ppf(conf_level, 1))`. If any of the three has no valid covariance sigma (frozen, or XSPEC returns a non-positive/unavailable value), the function raises `RuntimeError` naming the parameter instead of returning a meaningless 0 or negative error.

**On the elevated global chi2.** For Cyg X-1 the global reduced chi2 is typically >2: `tbabs*powerlaw` (I is unpolarised by `polconst`) does not fully describe the Stokes-I continuum over 2-8 keV. But `polconst` multiplies only Q/U, so PD/PA are set by the well-fit (reduced chi2 ~ 1) polarization sector and are robust despite the elevated global chi2 -- confirmed by a per-spectrum chi2 decomposition (I: reduced chi2 ~5-11; Q/U: reduced chi2 ~0.9-1.15) and by agreement with the independent `ixpe_pcube` broadband PD/PA to <0.4 sigma (see Changelog).

Companion figure: the existing `plot_fit` (I-spectrum + residuals) -- no separate plotter for this tool.

---

#### `xspec_help.fit_tbabs_pollin_powerlaw(indices, ...)`

Fits `tbabs*pollin*powerlaw` to every IXPE Stokes I/Q/U triplet in an indices dict, tying `pollin`'s linear PD(E)=`A1`+`Aslope`*(E-1 keV) and PA(E)=`psi1`+`psislope`*(E-1 keV) into **one free (A1, Aslope, psi1, psislope) quartet for the whole session** -- the paper's actual linear energy-dependence route for Cyg X-1 polarization (an XSPEC `pollin` fit), replacing worker 0's inline `runfit.py` `scipy` fit on `ixpe_pcube`'s PD(E) points, which was not toolkit-grounded.

```python
def fit_tbabs_pollin_powerlaw(
    indices: dict,
    gamma_init: float = 2.0,
    pd_a_init: float = 0.02,
    pd_slope_init: float = 0.004,
    pa_a_init: float = -20.0,
    pa_slope_init: float = 0.0,
    nh: float = 0.6,
    freeze_nh: bool = True,
    erange: tuple[float, float] = (2.0, 8.0),
    conf_level: float = 0.683,
    n_iterations: int = 100,
    log_path: str | None = None,
) -> dict
```

Same load path as `fit_tbabs_polconst_powerlaw`: call after `load_data(..., ixpe_changes=True, ixpe_src_rad=..., ixpe_bkg_inner_rad=..., ixpe_bkg_outer_rad=...)`, and pass an IXPE-only `indices` dict (`Model()` here applies to every currently active data group). The same default cross-data-group parameter tying `fit_tbabs_polconst_powerlaw` relies on ties `pollin`'s four parameters across every DU/ObsID with no explicit tie/untie call.

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `indices` | `dict` | -- | Return value of `load_data`; must contain `'IXPE'`. One ObsID gives a single-ObsID fit, more than one a joint combined fit, always with exactly one tied pollin quartet. |
| `gamma_init` | `float` | `2.0` | Initial powerlaw photon index (Cyg X-1 starting point). |
| `pd_a_init` | `float` | `0.02` | Initial pollin PD at 1 keV (`A1`, fraction 0-1, not percent). |
| `pd_slope_init` | `float` | `0.004` | Initial pollin PD energy slope (`Aslope`, fraction per keV) -- Cyg X-1's PD rises roughly 2.6%->5.0% over 2-8 keV. |
| `pa_a_init` | `float` | `-20.0` | Initial pollin PA at 1 keV (`psi1`, degrees). |
| `pa_slope_init` | `float` | `0.0` | Initial pollin PA energy slope (`psislope`, degrees per keV). |
| `nh` | `float` | `0.6` | `TBabs` equivalent hydrogen column, in units of 1e22 cm^-2 (Cyg X-1 ISM value). |
| `freeze_nh` | `bool` | `True` | Freeze `nH` at `nh`. Set `False` to let it float. |
| `erange` | `tuple[float, float]` | `(2.0, 8.0)` | keV range noticed on every IXPE spectrum before fitting. |
| `conf_level` | `float` | `0.683` | Confidence level the reported errors are scaled to, via the chi-squared quantile for 1 degree of freedom (`scipy.stats.chi2.ppf`). |
| `n_iterations` | `int` | `100` | `Fit.nIterations` per convergence round (up to 5 rounds, stopping early once the fit statistic stabilises). |
| `log_path` | `str \| None` | `None` | If given, the verbose `AllModels.show()`/`Fit.show()` summary is logged here (chatter 10, restored to 5 afterwards). |

**Returns** a dict with `PD_A`, `PD_A_err` (fraction), `PD_slope`, `PD_slope_err` (fraction per keV), `PA_A`, `PA_A_err` (degrees), `PA_slope`, `PA_slope_err` (degrees per keV), `Gamma`, `Gamma_err`, `nH`, `chi2`, `dof`, `stat_per_dof`, `conf_level`, `indices` (unchanged, for downstream `plot_fit`), `stdout`, `exit_code`.

**Errors are the covariance-matrix parameter-table 1-sigma**, exactly as in `fit_tbabs_polconst_powerlaw`: `Parameter.sigma` scaled by `sqrt(chi2.ppf(conf_level, 1))`. If any of PD_A/PD_slope/PA_A/PA_slope/Gamma has no valid covariance sigma, the function raises `RuntimeError` naming the parameter instead of returning a meaningless 0 or negative error.

Companion figure: the existing `plot_fit` -- no separate plotter for this tool. A separate PD(E)/PA(E) energy-dependence overlay across ObsIDs is a distinct tool.

---

#### `xspec_help.robust_fit(verbose=True, max_iter_simplex=10, max_iter_leven=100)`

> ⚠ **Not currently working.** This helper is presently broken and should not be relied on — the Simplex → Levenberg-Marquardt chain does not behave as intended. Use XSPEC's `Fit.perform()` directly until this is fixed.

Intended design: two-stage fitting — Simplex (derivative-free global exploration) followed by Levenberg-Marquardt (gradient polish) — with the fit method restored to Levenberg-Marquardt on exit.

**Returns** `(stat, dof)` — final fit statistic and degrees of freedom.

---

#### `xspec_help.inspect_fit(spec_groups, ...)`

Multi-instrument spectrum + residual plot.

```python
def inspect_fit(
    spec_groups: dict,          # {'label': [spec_idx, ...], ...}
    plot_type: str = 'ufspec',
    residuals: tuple = ('delchi',),
    title: str = '',
    xlim: tuple | None = None,
    xscale: str = 'log',
    yscale: str = 'log',
    figsize: tuple | None = None,
    legend: bool = True,
    show_model: bool = True,
    show_components: bool = True,
) -> tuple[Figure, list[Axes]]
```

`spec_groups` maps a label string to a list of XSPEC spectrum indices, e.g. `{'IXPE Det': [7, 10, 13], 'NICER': [5, 6], 'NuSTAR': [1, 2, 3, 4]}`. Each spectrum index produces one error-bar series. When `show_components=True`, `Plot.add = True` is set temporarily and every additive model component is overlaid as a dashed line (one copy, from the first spectrum of the first group). Returns `(fig, axes)`.

---

#### `xspec_help.inspect_epoch(indices, e_bins, ...)`

One-call wrapper: broadband spectrum + residuals panel **and** per-DU Stokes I/Q/U decomposition figures. Builds `spec_groups` and `iqu_groups` automatically from the `load_data` return dict, so only instruments present in the epoch appear.

```python
def inspect_epoch(
    indices: dict,
    e_bins,
    title: str = '',
    xlim: tuple[float, float] = (1.0, 79.0),
    xscale: str = 'log',
    yscale: str = 'log',
    plot_type: str = 'ufspec',
    residuals: tuple = ('delchi', 'ratio'),
) -> tuple[Figure, list[Axes], list[tuple]]
```

**Returns** `(fig, axes, figs)` where `figs` is the list of per-DU `(fig, (ax_i, ax_q, ax_u))` tuples from `plot_iqu_decomp` (empty if IXPE is absent).

---

#### `xspec_help.plot_iqu_decomp(iqu_groups, e_bins, ...)`

Per-DU Stokes I/Q/U three-panel figures with additive component overlays, rebinned to `e_bins` using inverse-variance weighting.

```python
def plot_iqu_decomp(
    iqu_groups: dict,       # {'DU1': [I_idx, Q_idx, U_idx], ...}
    e_bins,
    show_components: bool = True,
    title_prefix: str = '',
    figsize: tuple = (14, 4),
) -> list[tuple[Figure, tuple[Axes, Axes, Axes]]]
```

Returns one `(fig, (ax_i, ax_q, ax_u))` tuple per DU label.

---

#### `xspec_help.rebin_plot(spec_idx, e_edges, plot_type='ufspec')`

Extracts XSPEC plot data for one spectrum and rebins to coarser energy edges using inverse-variance weighting. For visualisation only — does not affect the fit.

**Returns** `(e_cen, y, yerr)` — band-centre energies, weighted mean y values, and propagated 1σ errors.

---

#### `xspec_help.group_pha(infile, outfile, e_edges)`

Writes a grouped copy of a PHA file with custom energy bin edges derived from the file's own RMF EBOUNDS. Load `outfile` into XSPEC instead of `infile` when coarser binning is needed.

```python
group_pha(infile: str, outfile: str, e_edges: sequence[float]) -> None
```

Example:
```python
e_bins = [2., 3., 4., 5., 6., 8.]
q_src = red_ixpe.get_pha_q(f_src)
q_grp = q_src.replace('.fits', '_grp.fits')
xspec_help.group_pha(q_src, q_grp, e_bins)
AllData(f'{grp}:{idx} {q_grp}')
```

---

#### `xspec_help.freeze_all(thaw=False)`

Freezes (or thaws when `thaw=True`) every non-linked parameter across all data groups. Handy for **snapshotting the current fit**: call `freeze_all()` once you're happy with a set of best-fit values so you can save them. Unfreezing later with `freeze_all(thaw=True)` is not preferred but workable: it will break the frozen const=1 on the reference spectrum.

---

#### `xspec_help.audit_paths()`

Prints a formatted table of file paths (PHA, background, RMF, ARF) and net count rates for every currently loaded XSPEC spectrum. Called automatically by `fit_tbabs_polconst_diskbb_polconst_relxillCp` when `audit_data=True`.

---

### `red_ixpe` — IXPE reduction

---

#### `red_ixpe.get_all_file_list(ixpe_path, chatter=False)`

Scans `ixpe_path` for numeric observation sub-directories and returns a flat list of all Level-2 event files (`*evt2_v??_rej.fits`). Handles the special cases `01002901_picorr` (pi-corrected files) and `03002599` (non-`_rej` files). Raises `ValueError` if no files are found.

**Returns** `list[str]` — absolute paths to all event files, sorted within each observation.

---

#### IRF calibration epoch (auto-selected)

ixpeobssim has no concept of IRF validity dates (its own `irf/caldb.py`: *"we don't really have a concept of validity date"*), so the reduction must choose which `obssimYYYYMMDD` calibration epoch to use. Every IXPE reduction entry (`pcube`, `pha`, `reduce_pha`, `pcube_timeseg`, and `xspec_help.load_data`) takes `irf_date` as `None` by default and delegates to `lib/irf_epoch.py`, which reads the observation's midpoint date and selects the **installed epoch nearest in time** (`irf_epoch.resolve_irf_date`) — so omitting `irf_date` yields the temporally-appropriate response and a wrong-epoch bias is not something the caller can introduce by accident. Passing an explicit `'YYYYMMDD'` overrides it. The AI-facing `skills/` tools do **not** expose `irf_date` at all — the epoch is chosen for them.

#### `red_ixpe.pcube(obs, file_list, ...)`

Extracts an IXPE polarization cube (PCUBE) for one observation, **or a native multi-ObsID combination**, background-subtracts it, and returns both plus the `irf_date` used.

```python
def pcube(
    obs: str | list[str],
    file_list: list[str],
    ebinning: list[float] | None = None,
    src_rad: float | None = None,
    bkg_inner_rad: float | None = None,
    bkg_outer_rad: float | None = None,
    ra: float | None = None,
    dec: float | None = None,
    irf_date: str | None = None,
    changes: bool = True,
    overwrite: bool = True,
    chatter: bool = False,
) -> tuple[xBinnedPolarizationCube, xBinnedPolarizationCube, str | None]
```

`obs` is a single 9-digit ObsID, or a **list** of ObsIDs to natively combine: every listed ObsID's Level-2 files are filtered from `file_list` (DU-validated per ObsID via `_validate_du_files` -- each must resolve to exactly one file per DU1/DU2/DU3), concatenated, and passed into **one** `xpselect`/`xpbin PCUBE` call, so the returned cube sums raw Stokes I/Q/U/COUNTS across every DU *and* every listed ObsID before PD/PA are derived -- a native directly-combined reduction, not a post-hoc merge of separately-computed per-ObsID PD/PA (this falls out of `xBinnedPolarizationCube.from_file_list` itself: it `__iadd__`s every file `xpbin` returns, regardless of which DU/ObsID produced it).

When `changes=False` only `obs` and `file_list` are required (and `obs` must be a single ObsID -- reuse of on-disk products is not defined for a combined multi-ObsID request); existing `_src_pcube` / `_bkg_pcube` FITS files on disk are read directly. When `changes=True` the spatial parameters (`ebinning`, `src_rad`, `bkg_inner_rad`, `bkg_outer_rad`, `ra`, `dec`) are required, `ebinning` must carry >= 2 edges, and the full `ixpeobssim` pipeline (`xpselect` → `xpbin PCUBE`) is re-run; `irf_date` is auto-selected when `None` -- the midpoint of the earliest start to the latest stop across every listed ObsID (see **IRF calibration epoch** above; for a single ObsID this is unchanged). Background scaling is applied as `bkg_pcube *= src_pcube.backscal() / bkg_pcube.backscal()` before subtraction.

**Returns** `(src_pcube, bkg_pcube, irf_date)` -- the last element is the epoch actually used (resolved when `changes=True`; whatever was passed, possibly `None`, when `changes=False`).

**`red_ixpe.pcube_bands(spc)`** reads the per-energy-band PD/PA/Stokes/MDP99/SIGNIF fields off a PCUBE (one dict per bin in `spc`'s own `ENERG_LO`/`ENERG_HI` grid) via `polarization()` and direct attribute access -- never hand-rolled Stokes arithmetic. `PD`, `PD_err`, `QN`, `QN_err`, `UN`, `UN_err`, `MDP99` are dimensionless fractions (the units `xBinnedPolarizationCube` itself stores them in, not percent); `PA`, `PA_err` are degrees.

**`red_ixpe.write_pcube_csv`/`read_pcube_csv`** write/read a provenance-headed CSV (`# ixpe_pcube provenance: {...}` first line, mirroring `pcube_timeseg`'s convention), one row per energy band. There is no separate "integrated" sentinel row (unlike `pcube_timeseg`'s `segment_id=-1`): a 2-edge `ebinning` naturally yields exactly one row spanning the full requested band.

---

#### `_get_obs_file_list` matches by observation directory, not substring

`_get_obs_file_list(obs, file_list)` (used internally by `pcube`/`pha`/`reduce_pha`/`HR`/`pcube_timeseg`) filters `file_list` by each file's own observation directory (the parent of its `event_l2/`, i.e. `.../<obs>/event_l2/<file>` -- the layout `get_all_file_list` itself builds), not a raw substring test. This matters once `file_list` spans more than one ObsID's directory (as `get_all_file_list`'s mission-wide file_list does, e.g. for `ixpe_pcube`'s combined-ObsID case): a substring filter for `obs='01002901'` also matched `01002901_picorr`'s files (a differently-processed variant stored under its own suffixed directory, whose filenames still literally embed the same 9-digit ObsID), which then read as a spurious "two files for the same DU" ambiguity and made `_validate_du_files` raise on a perfectly good request. Requesting `obs='01002901_picorr'` explicitly still resolves to that variant's own files, unaffected.

---

#### `red_ixpe.pcube`/`ixpe_pcube` -- energy-resolved and multi-ObsID-combined polarimetry

The `ixpe_pcube` AI tool (`skills/session.py`) surfaces `red_ixpe.pcube` for **energy-resolved** polarimetry (a caller-supplied `ebinning` list, default `[2, 3, 4, 6, 8]` keV -- 4 bands, the Krawczynski+2022 Table S2 grid) of one local IXPE ObsID, or a **native multi-ObsID combination**, e.g. reproducing Krawczynski+2022's Cyg X-1 Fig. S5/Table S2 and its combined May (`01002901`) + June (`01250101`) cube.

```python
def ixpe_pcube(
    obs_db_path: str,
    obs: str | list,
    ra: float,
    dec: float,
    ebinning: list | None = None,
    src_rad: float = 1.0,
    bkg_inner_rad: float = 3.0,
    bkg_outer_rad: float = 5.0,
    irf_date: str | None = None,
) -> dict
```

Resolves the mission's full Level-2 event file list via `red_ixpe.get_all_file_list(obs_db_path)`, then calls `red_ixpe.pcube` (`obs` a string reduces just that ObsID; `obs` a list combines every listed ObsID's files into one `xpselect`/`xpbin` call, per `pcube`'s own contract above). `ra`/`dec` are **required** -- unlike `ixpe_pcube_timeseg`, this tool never auto-centroids, since a multi-ObsID combination especially needs one explicit, shared region centre. The IRF epoch is auto-selected (not AI-facing) as the installed epoch nearest the midpoint across every requested ObsID.

Reads `PD`, `PD_err`, `PA`, `PA_err`, `QN`, `QN_err`, `UN`, `UN_err`, `MDP99`, `SIGNIF`, `E_LO`, `E_HI` per band via `red_ixpe.pcube_bands` on the background-subtracted source cube. `PD`/`PD_err`/`QN`/`QN_err`/`UN`/`UN_err`/`MDP99` are dimensionless fractions (not percent, unlike `ixpe_pcube_timeseg`'s `PD_pct`/`MDP99`); `PA`/`PA_err` are degrees.

**Raises** (surfaced by the wrapper as `exit_code=1`/`stderr`, nothing written or recorded) if a requested ObsID has no local Level-2 event files, is missing DU1/DU2/DU3 (or carries more than one file for a DU), or `ebinning` has fewer than 2 edges.

**Writes** a provenance-headed CSV, one row per energy band -- a single ObsID's CSV at `<obs_db_path>/<obsid>/ixpe_pcube_<obsid>.csv`; a combined request's CSV at `<obs_db_path>/ixpe_pcube_<obsid1>_<obsid2>_....csv` (mission root -- no single ObsID directory owns a combined product, mirroring `bhxrb_hardness_ratios.csv`'s placement). Records one `kind='pcube'` product per contributing ObsID via `obsdb_catalog.neutral.advance` (level 3), all pointing at the same CSV for a combined request.

**Returns** a dict with `bands` (list of per-band dicts), `provenance` (`obsids`, `ebinning`, `src_rad`, `bkg_inner_rad`, `bkg_outer_rad`, `ra`, `dec`, `irf_date`, `combined`:bool), `csv_path`, `stdout`, `exit_code`.

**No figures.** This is a data-only tool; `ixpe_pcube_plot` (below) is the companion plotter, mirroring how `pa_pd_stability.stability_plot` sits alongside `pa_pd_stability.pa_pd_stability`.

---

#### `ixpe_pcube_plot.ixpe_pcube_plot(csv_path, ...)`

Pure-matplotlib (Agg backend) companion plotter for `ixpe_pcube`: renders the PD(E)/PA(E) + Q-U-plane reviewer figure (Krawczynski+2022 Fig. S5 style) from an `ixpe_pcube_<obsid(s)>.csv`. Imports neither `xspec` nor `ixpeobssim` (its own small, independent reimplementation of `red_ixpe.read_pcube_csv`'s provenance-header contract avoids importing `red_ixpe`, which pulls in `ixpeobssim` at module load) -- so it stays in `lib/` and is wired directly into `TOOL_MANIFEST` (`ixpe_pcube_plot`), never through `skills/session.py`.

```python
def ixpe_pcube_plot(
    csv_path: str,
    save_path: str | None = None,
    title: str | None = None,
    pollin: dict | None = None,
) -> dict
```

**Writes exactly one PNG** with three panels: PD vs band-centre energy (top left), PA vs band-centre energy (bottom left, sharing the energy axis), and QN vs UN -- the normalized Stokes/Q-U plane (right, color-coded by band-centre energy) -- each with 1-sigma error bars. `PD`/`PD_err`/`PA`/`PA_err`/`QN`/`QN_err`/`UN`/`UN_err` and the band-centre energies (`0.5*(E_LO+E_HI)`) are plotted exactly as read from `csv_path`'s own columns -- never recomputed -- and are also returned under `plotted` for differential verification against the source CSV.

**PD(E)/PA(E) model overlay.** The reviewer figure for the energy-dependence claim (Krawczynski+2022 Fig. S5 style): the PD/PA panels overlay a linear model line on top of the PCUBE points.
- `pollin` given -- a `PD_A`/`PD_slope`/`PA_A`/`PA_slope` quartet (e.g. `fit_tbabs_pollin_powerlaw`'s result dict, passed straight through; extra keys such as `Gamma`/`*_err` are ignored) -- draws the analytic `pollin` fit PD(E)=`PD_A`+`PD_slope`*(E-1 keV), PA(E)=`PA_A`+`PA_slope`*(E-1 keV): the model actually fit to the spectra, the A/B-else-fallback design's "A" path.
- `pollin` omitted (the "B" fallback) -- an inverse-variance-weighted linear fit of the CSV's own PD(E)/PA(E) points (`_weighted_linear_fit`, `scipy.optimize.curve_fit` with `sigma=`the column's own `*_err`, `absolute_sigma=True`) is drawn instead.

Both paths draw the same way (dashed line spanning `[min(E_LO), max(E_HI)]`) and are returned under `plotted['overlay']`.

**Edge cases.** A band with any non-finite value among `E_LO`/`E_HI`/`PD`/`PD_err`/`PA`/`PA_err`/`QN`/`QN_err`/`UN`/`UN_err` is skipped -- never plotted as zero (counted in `plotted['n_skipped']`). A single remaining band still renders as a lone marker (no second point to draw a connecting line to); in the fallback path, fewer than two usable bands leaves no line drawable, so none is drawn (`plotted['overlay'] == {'kind': None}`) -- handled gracefully, not an error. `csv_path` missing, missing the `ixpe_pcube` provenance header, carrying zero band rows, having every band skipped for non-finite values, or `pollin` supplied but missing any of `PD_A`/`PD_slope`/`PA_A`/`PA_slope` all **raise** `ValueError` with a message naming the problem -- an empty or mislabeled input never degrades to a silently-wrong plot.

**Determinism.** Same convention as `stability_plot`/`pd_hr_survey_plot`: fixed `dpi=150`, `bbox_inches='tight'`, and explicit `savefig` metadata (a fixed `Software` string, no embedded timestamp), so repeated calls on unchanged input are byte-identical -- including the fallback fit, which is itself deterministic (closed-form weighted least squares) on unchanged input.

**Returns** a dict with `figures` (one `FigureRef` dict -- `path`, `caption`, `mime_type` (`image/png`), `produced_by`, `created_at`, `sha256` of the written bytes), `plotted` (the exact `E_LO`/`E_HI`/`E_ctr`/`PD`/`PD_err`/`PA`/`PA_err`/`QN`/`QN_err`/`UN`/`UN_err` arrays drawn, plus `n_bands`/`n_skipped`, and `overlay` -- `{'kind': None}` when no line was drawn, else `{'kind': 'pollin'|'fallback', 'PD_A', 'PD_slope', 'PA_A', 'PA_slope', ...}` (the fallback also carries each parameter's fit error)), `stdout`, `exit_code` (always 0 -- failures raise instead of degrading to an empty plot).

---

#### `red_ixpe.pha(obs, file_list, ...)`

Extracts IXPE PHA1 Stokes I spectra (plus matching PHA1Q and PHA1U files) for one observation via `ixpeobssim`.

```python
def pha(
    obs: str,
    file_list: list[str],
    src_rad: float | None = None,
    bkg_inner_rad: float | None = None,
    bkg_outer_rad: float | None = None,
    ra: float | None = None,
    dec: float | None = None,
    irf_date: str | None = None,
    changes: bool = True,
    overwrite: bool = True,
    chatter: bool = False,
) -> tuple[list[str], list[str]]
```

Same `changes` toggle as `pcube`. When `changes=True`, runs `xpselect` + `xpbin PHA1` with IRF name `ixpe:obssim{irf_date}_alpha075:v13` (`irf_date` auto-selected when `None` — see **IRF calibration epoch** above). One source and one background PHA1 file is produced per detector unit (DU1, DU2, DU3). The Q and U files are derived from the I paths by `red_ixpe.get_pha_q` / `red_ixpe.get_pha_u`.

**Returns** `(src_pha_list, bkg_pha_list)` — each a list of three paths ordered DU1, DU2, DU3.

---

#### `red_ixpe.reduce_pha(obs, file_list, ...)`

Reuses an observation's on-disk per-DU PHA1 (Stokes I) src+bkg products, or derives them. This is the tool wrapped as the `ixpe_reduce_pha` AI-facing entry (`skills/session.py`).

```python
def reduce_pha(
    obs: str,
    file_list: list[str],
    src_rad: float = 1.0,
    bkg_inner_rad: float = 3.0,
    bkg_outer_rad: float = 5.0,
    irf_date: str | None = None,
    overwrite: bool = True,
    chatter: bool = False,
) -> dict
```

Verifies exactly one Level-2 event file per DU1/DU2/DU3 (`_validate_du_files`), then checks whether the src+bkg PHA1 paths `pha(..., changes=False)` derives already exist on disk. If so, they are reused as-is: no reduction, no source centroiding, `ra`/`dec`/`src_reg`/`bkg_reg` are `None`. Otherwise the source position is centroided from the DU1 event file via `find_src.run_extraction` (writing `src.reg`/`bkg.reg` into the observation's own directory, one level above `event_l2/` -- **not** read from a FITS header, and **not** from the fixed `src_pos.SRC_COORD` table), and `pha()` is run (`xpselect` + `xpbin PHA1`) with a circular source region of radius `src_rad` and an annular background (`bkg_inner_rad`-`bkg_outer_rad`) centred there, at the auto-selected IRF epoch (or `irf_date` when given). The resolved epoch is returned under the `irf_date` key.

**Returns** a dict with `src_pha`, `bkg_pha` (each a 3-element DU1/DU2/DU3 list of Stokes-I PHA1 paths), `ra`, `dec` (`None` when reused), `reused` (bool), `src_reg`, `bkg_reg` (paths when freshly reduced, `None` when reused).

**`skills/session.py`'s `ixpe_reduce_pha` wrapper** locates the ObsID's Level-2 event files under `<obs_db_path>/<obsid>/event_l2/` (preferring `_rej` files, falling back to plain `evt2_v*.fits`, the same convention as `ixpe_pcube_timeseg`), calls `reduce_pha`, and records every PHA1 product on disk after the call (3 DUs x src/bkg) in the mission catalog at level 3 via `obsdb_catalog.neutral.advance` -- one `kind='pha'` product per file, `params` carrying `du`, `stokes='I'`, `region` (`'src'`/`'bkg'`), the region radii, `irf_date`, and `reused` -- whether the call reused existing products or ran a fresh reduction. Nothing is recorded on failure (e.g. no local Level-2 files, or other than one file per DU).

---

#### `red_ixpe.get_pha_q(pha_path)` / `red_ixpe.get_pha_u(pha_path)`

Derive the Stokes Q or U PHA path from the Stokes I PHA path by replacing the first occurrence of `pha1` with `pha1q` or `pha1u` in the filename.

---

#### `red_ixpe.get_mjd(obs, file_list)`

Reads `TSTART`, `TSTOP`, `MJDREFI`, `MJDREFF`, and `TIMEZERO` from the EVENTS header of the first event file for `obs`.

**Returns** `(mjd_start, mjd_stop, mjd_mid)` as floats.

---

#### `red_ixpe.HR(obs, file_list, soft_band=(2.0, 4.0), hard_band=(4.0, 8.0), bkg=False, chatter=False)`

Net Stokes-I hardness ratio `HR=(H-S)/(H+S)` for one IXPE observation, read directly from its existing per-DU PHA1 (Stokes I) src(+bkg) products -- no XSPEC `AllData` loading at all. This is the tool wrapped as the `ixpe_HR` AI-facing entry (`skills/session.py`).

```python
def HR(
    obs: str,
    file_list: list[str],
    soft_band: tuple[float, float] = (2.0, 4.0),
    hard_band: tuple[float, float] = (4.0, 8.0),
    bkg: bool = False,
    chatter: bool = False,
) -> dict
```

**DU handling.** Attempts DU1, DU2 and DU3 independently via `pha(..., changes=False)`'s existing file-path locators: whichever DUs' Level-2 event files are present in `file_list` contribute -- passing fewer than 3 DUs' event files is how a caller expresses "DU off" (IXPE can genuinely fly with a DU disabled), and this never hard-fails on a truly absent DU. A DU that IS present in `file_list` but has no PHA1 product on disk raises (run `reduce_pha`/`ixpe_reduce_pha` first) rather than being silently dropped -- the defect this replaces (the prior implementation dropped DU3 from the summed rate via an off-by-one loop bound over the XSPEC spectra it had loaded). `file_list` carrying more than one Level-2 event file for the same DU (two processing versions, or a mixed plain+`_rej` set) also raises, naming the duplicated DU(s) and the offending files, before any PHA read -- `HR` never auto-picks a "preferred" version, since doing so silently would either double-count or silently drop real data; the caller/wrapper owns variant selection and must pass exactly one file per DU. Net per-band rate sums over the contributing DUs; errors add in quadrature over those DUs.

**Convention.** `HR=(H-S)/(H+S)`, bounded in `[-1, +1]` -- not `H/S`, the prior (non-community-standard) convention.

**Error.** `sigma_HR = (2/(H+S)^2) * sqrt(S^2*sigma_H^2 + H^2*sigma_S^2)`, the closed-form first-order Gaussian propagation for independent net band rates H, S with 1-sigma errors sigma_H, sigma_S (`hr_from_bands`). Per-band rate and error are read straight from each PHA1's own `RATE` (cts/s, summed directly) and `STAT_ERR` (quadrature-summed) columns over the in-band channels, then over contributing DUs -- replacing the prior implementation's error, which came from misreading PyXSPEC's `Spectrum.rate` tuple (`total, background, net, net_error`) as `(rate, variance)`.

**Band -> channel.** Each PHA1's own `RESPFILE` header keyword names its RMF; `E_MIN`/`E_MAX` from that RMF's `EBOUNDS` extension resolve `soft_band`/`hard_band` to channels (a channel's bin midpoint must fall inside `[lo, hi)`) -- never a hardcoded channel grid. When a PHA1's `RESPFILE` predates a reduction-environment migration and no longer resolves on disk (the real state of most already-reduced Cyg X-1 IXPE PHA1 products on this machine -- their `RESPFILE` points at a now-removed `~/.local` ixpeobssim install), `_resolve_rmf_ebounds` falls back to an installed CALDB RMF for the same DU in the currently-active ixpeobssim installation: verified byte-identical to every other IXPE IRF date/DU's `EBOUNDS` (375 channels, 0.04 keV wide -- a fixed instrumental characteristic; only the RESPONSE matrix varies with date/DU), so this is the same grid, not an approximation. Guards the fallback by raising if the resolved grid's channel count disagrees with the PHA's own channel count.

**Background.** `bkg=True` subtracts each DU's background PHA1, scaled by `src_BACKSCAL/bkg_BACKSCAL`, before summing; `bkg=False` (default) uses the source PHA1 alone and never reads the background file.

**Returns** a dict with `HR`, `HR_err` (`nan`/`nan` if `soft_rate+hard_rate <= 0` -- an honest, genuinely undefined ratio, not a crash), `soft_rate`, `soft_err`, `hard_rate`, `hard_err` (net, DU-summed cts/s and 1-sigma error), and `dus_used` (sorted list of the DU numbers that contributed).

**`skills/session.py`'s `ixpe_HR` wrapper** surveys a list of ObsIDs: for each, locates Level-2 event files under `<obs_db_path>/<obsid>/event_l2/` (preferring `_rej` files, falling back to plain `evt2_v*.fits`, the same convention as `ixpe_reduce_pha`/`ixpe_pcube_timeseg`), calls `HR`, and either appends a row or -- on any failure (no local event files, or `HR` raising) -- records the ObsID in the returned `errors` dict and excludes it from the CSV, never writing a dummy row for it or crashing the rest of the survey. Writes `bhxrb_hardness_ratios.csv` (columns `obsid, soft_rate, soft_err, hard_rate, hard_err, HR, HR_err`) at `<obs_db_path>/bhxrb_hardness_ratios.csv`, overwriting any previous run, and records one `kind='hardness_ratio'` product per successful ObsID via `obsdb_catalog.neutral.advance` (level 3), with `params` carrying that ObsID's real `soft_band`/`hard_band`/`bkg`/`dus_used`. If every requested ObsID fails, nothing is written or recorded and `exit_code=1` with `stderr` set.

---

#### `red_ixpe.pcube_timeseg(obs, file_list, irf_date, ...)`

Time-resolved PCUBE polarimetry for one observation: splits it into natural elapsed-time segments and runs one single-band, DU-merged PCUBE reduction per segment plus one integrated whole-observation reduction, to test intra-observation polarization constancy (e.g. Cyg X-1, 4U 1630-47). This is the tool wrapped as the `ixpe_pcube_timeseg` AI-facing entry (`skills/session.py`).

```python
def pcube_timeseg(
    obs: str,
    file_list: list[str],
    irf_date: str | None = None,
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
    overwrite: bool = True,
    chatter: bool = False,
) -> dict
```

**Segmentation.** Edges are `TSTART + k*segment_duration_s` for `k = 0..n_full`, plus a final partial edge at `TSTOP` when the span does not divide evenly — natural elapsed time, **not** exposure-equalised.

**GTI pre-reduction gate.** Each window's GTI-clipped good exposure is computed via `_gti_overlap` *before* any `xpselect`/`xpbin` call; a window whose good exposure is below `min_seg_exposure_s` (default 1000.0s) is skipped without running the reduction at all, and counted in the returned `n_skipped_empty`. This matters for an observation with a long inter-pointing gap (two pointings glued together with a multi-day gap between them, e.g. 03002599): without the gate, every window landing inside the gap still ran a full reduction for a spurious near-zero-exposure row, and for a long observation this alone could exhaust the whole per-call time budget. The 1000s default was set from the real GTI structure of the reported pathological windows (a 1.7s window at 01250401's default-segmentation segment 7, a 241.0s window at 02001901's segment 148) — comfortably above both, comfortably below every other real segment measured in either observation (see `tests/test_ixpe_pcube_timeseg.py::test_empty_slice_gate_default_separates_known_pathological_segments`). A window that clears the gate but still yields zero counts or a failed reduction is warned and skipped without incrementing `n_skipped_empty` (the gate already counted it as expected).

**Incremental write + resume/override.** Rows are appended to `pcube_timeseg_<obs>.csv` and flushed (`flush()` + `fsync()`) as each segment's reduction completes, so a call that hits the sandbox's 5400s wall-clock cap (`tool_config.yaml`, held fixed rather than raised) leaves a genuine checkpoint instead of losing all progress. The CSV's first line is a provenance comment (`# pcube_timeseg provenance: {...}`) recording every parameter that determines the segmentation and reduction: `segment_duration_s, emin_kev, emax_kev, src_rad, bkg, bkg_inner_rad, bkg_outer_rad, ra, dec, irf_date, min_seg_exposure_s, tstart, tstop`. When `resume=True` (default) and an existing CSV's provenance matches this call's exactly, only the segments still missing (and the integrated row, if absent) are computed; matching rows already on disk are re-emitted without a repeat reduction. On any mismatch — `resume=False`, no existing CSV, or any provenance field differing (including a prior call with different arguments) — the stale CSV is deleted and rebuilt from scratch; rows from a different configuration are never glued onto a fresh run.

**Source position.** If both `ra` and `dec` are `None`, the source is centroided from one Level-2 event file via `find_src.run_extraction` (validated against this tool's accuracy anchor — see Changelog); the resolved value (not `None`) is what is recorded in the CSV's provenance.

**Background.** Off by default (`bkg=False`), matching `pcube`'s convention and appropriate for bright source-dominated targets. When `bkg=True`, an annulus (`bkg_inner_rad`-`bkg_outer_rad`) is reduced the same way and subtracted exactly as `pcube` does (`bpc *= spc.backscal()/bpc.backscal(); spc -= bpc`), per segment and for the integrated cube.

**`exposure_s` is *not* read off the PCUBE FITS header.** `pipeline.xpselect` only rewrites `TSTART`/`TSTOP`/`ONTIME`/`LIVETIME`/`DEADC` when `ltimeupdate=True`, and that option raises (`KeyError: 'LIVETIME'`) on Level-2 IXPE event files, which carry no per-event `LIVETIME` column (ixpeobssim supports only an average deadtime correction for L2 data). Left alone, every segment's PCUBE would silently inherit the identical, unmodified full-observation header values — the same defect class as segments inheriting the full-observation GTI. Instead, `exposure_s` and the segment's `t_start_mjd`/`t_stop_mjd` are computed from the observation's own GTI extension intersected with `[t0, t1)`, scaled by the (per-observation constant) `DEADC` keyword — mirroring the `exposure = dt * deadc` averaging ixpeobssim's own `LC`/`MDPMAPCUBE` binning uses for Level-2 data. The gate reuses this same GTI-clipped, DEADC-scaled value as `exposure_s`, so a surviving segment's `exposure_s` is exactly the quantity that was checked against `min_seg_exposure_s`.

**Returns** a dict with:
- `timeseries` — list of per-segment dicts, including any rows re-emitted from a resumed CSV
- `integrated` — dict, same fields, `segment_id=-1`
- `csv_path` — absolute path written
- `n_skipped_empty` — windows the GTI gate skipped before any reduction
- `n_expected` — windows the gate passed (i.e. attempted)
- `complete` — `True` iff every expected segment plus the integrated row is present
- `resumed` — `True` iff this call resumed from a matching on-disk CSV rather than starting fresh

Each row (segment and integrated) has keys `segment_id, t_start_mjd, t_stop_mjd, exposure_s, PD_pct, PD_err, PA_deg, PA_err, QI, QI_err, UI, UI_err, MDP99, SIGNIF` (the integrated row uses `segment_id=-1`). `PD_pct`/`PD_err`/`MDP99` are in percent; `PA_deg`/`PA_err` in degrees; `QI`/`UI` are the normalized Stokes parameters Q/I, U/I. The CSV is written to `<obsid_dir>/pcube_timeseg_<obs>.csv` (the observation directory inferred from `file_list`) with a provenance comment line followed by these same columns, in the same order.

**`skills/session.py`'s `ixpe_pcube_timeseg` wrapper never returns the full `timeseries` inline** — for a long observation (O(100-170) segments) the serialized per-segment array can exceed the sandbox's transport chunk limit even though the CSV itself is written fine. The wrapper instead returns a bounded summary: `csv_path`, `n_segments`, `n_skipped_empty`, `n_expected`, `complete`, `resumed`, `integrated`, `head`/`tail` (first/last 3 segment rows), plus `stdout`/`exit_code`. Read the CSV at `csv_path` for the complete per-segment table. The wrapper also records the CSV in the mission catalog at level 3 via `obsdb_catalog.neutral.advance`, with `complete` included in the recorded `params` so a partial (in-progress) reduction is never mistaken for a finished one.

---

### `pa_pd_stability` — PA/PD/Stokes constancy analysis

---

#### `pa_pd_stability.pa_pd_stability(obs_db_path, obsids, output_dir, ...)`

Cross-ObsID intra-observation PA/PD stability analysis for IXPE, computed entirely from the `pcube_timeseg_<obsid>.csv` files `red_ixpe.pcube_timeseg` already wrote (run `ixpe_pcube_timeseg` first — this tool never reduces data itself). Pure numpy/scipy: imports neither `xspec` nor `ixpeobssim`, unlike every other module in this package, so it runs standalone and fast, and is wired directly in `TOOL_MANIFEST` (module `xspec_utils.lib.pa_pd_stability`, no `skills/session.py` involvement — the same pattern as `group_pha`).

```python
def pa_pd_stability(
    obs_db_path: str,
    obsids: list[str],
    output_dir: str,
    sources: dict | None = None,
    variability_p_threshold: float = 0.01,
    sf_lags_ks: list[float] | None = None,
    min_seg_for_sf: int = 10,
) -> dict
```

**Per-ObsID statistics** (`compute_obsid_stability`, one call per `obsid`):

1. **Weighted constancy tests.** Inverse-variance weighted mean `x_bar = sum(w_i x_i)/sum(w_i)`, `w_i = 1/sigma_i^2`, and `chi2 = sum(w_i (x_i - x_bar)^2)`, `dof = N-1`, for PD, QI, and UI. PA uses mod-180 circular handling: the weighted mean is the doubled-angle circular mean (`atan2(sum(w sin 2θ), sum(w cos 2θ))/2`, wrapped to `(-90, +90]`), and residuals are wrapped via `((d+90) % 180) - 90` before the chi-squared sum, so a PA flip between e.g. +89° and -89° (physically ~2° apart) is never scored as a ~178° jump. A quantity's non-finite values and non-finite/non-positive errors are excluded per quantity (independently) before the fit; `chi2`/`p` are NaN and `variable=False` when fewer than 2 usable points remain.
2. **Authoritative 2D Stokes test.** `chi2_2D = chi2_QI + chi2_UI` (QI and UI computed on their shared valid-data mask so `dof_2D = 2*(N-1)` is well defined), `p_2D = scipy.stats.chi2.sf(chi2_2D, dof_2D)`, and `variable = p_2D < variability_p_threshold`. This is the authoritative variability call; the PD/PA chi-squared tests above are secondary/descriptive.
3. **Moment-estimator intrinsic scatter**, for PD and PA (circular residuals): `sigma_int = sqrt(max(0, sigma_obs^2 - mean_var))`, where `sigma_obs^2` is the sample variance (ddof=1) of the values about their mean and `mean_var` is the unweighted mean of the per-point `err^2`; `*_upper_limit=True` when the excess term is `<= 0`.
4. **Peak-to-peak ranges** `PD_ptp`, `PA_ptp` (PA computed from the wrapped residuals about its circular mean, so a branch-cut crossing does not inflate the range).
5. **First-order structure functions** `SF_PD(tau)`, `SF_PA(tau)` (PA differences wrapped) at each lag in `sf_lags_ks` (default `[5, 10, 20]` ks), only when an ObsID has `>= min_seg_for_sf` (default 10) segments; the tolerance half-width is half the ObsID's own `segment_duration_s` (read from the CSV's provenance header), falling back to half the median inter-segment spacing if absent. NaN for a lag with no segment pair inside its tolerance window.
6. **Cross-source Spearman correlation** (across all successfully-processed ObsIDs together, not per-ObsID) of each ObsID's mean MDP99 against its reduced chi-squared for PD and for PA (`scipy.stats.spearmanr`); NaN with a `note` when fewer than 3 ObsIDs have finite values.

**Edge cases handled without raising.** An ObsID whose CSV is missing, or fails the same provenance-header check `red_ixpe._read_existing_timeseg_csv` uses (see below), is recorded as an error entry in `per_obs` and the catalog's `error` column, and skipped — it never appears in the Spearman correlation. `N < 2` valid points for any quantity yields `chi2=NaN`, `variable=False`, never a crash.

**The CSV reader (`_read_timeseg_csv`) is an independent reimplementation**, not a shared helper: it never imports `red_ixpe` (which pulls in `ixpeobssim`/`xspec` at module load), but mirrors its contract exactly — same provenance-comment-header format, same `segment_id=-1` integrated-row convention, same "missing/unrecognised header → treat as absent" rule (verified by a differential test against `red_ixpe._read_existing_timeseg_csv` on real data).

**Writes** (to `output_dir`, created if absent): `pa_pd_stability_catalog.csv` (one row per ObsID: `obsid, source, N_seg`, integrated PD/PA/MDP99/SIGNIF, `chi2/dof/p` for PD/PA/QI/UI, `chi2_2D/dof_2D/p_2D/variable`, `sigma_PD_int/PD_int_upper_limit`, `sigma_PA_int/PA_int_upper_limit`, `PD_ptp/PA_ptp`, `rchi2_PD/rchi2_PA`, `mean_MDP99`, `error`) and `pa_pd_variability_summary.txt` (per-ObsID 2D verdict + secondary stats, the two Spearman results, and the cited source CSV paths). Does not copy the per-segment input CSVs.

**Returns** a dict with `per_obs` (obsid → the stats above, including `sf_pd`/`sf_pa`, or `{'error': ...}`), `spearman`, `catalog_csv`, `summary_txt`, `exit_code` (always 0 — a per-ObsID problem is recorded, not raised).

---

#### `pa_pd_stability.stability_plot(obs_db_path, obsids, catalog_csv, ...)`

Pure-matplotlib (Agg backend) companion plotter for `pa_pd_stability`: renders the four standard reviewer diagnostic figures from its outputs. Imports neither `xspec` nor `ixpeobssim`, so — like `pa_pd_stability` itself — it stays in `lib/` and is wired directly into `TOOL_MANIFEST` (`ixpe_pa_pd_stability_plot`), never through `skills/session.py`.

```python
def stability_plot(
    obs_db_path: str,
    obsids: list[str],
    catalog_csv: str,
    save_path_prefix: str | None = None,
    sf_lags_ks: list[float] | None = None,
    sf_min_segments: int = 10,
) -> dict
```

**Writes exactly four files**, `{save_path_prefix}_<suffix>.png`:

1. `_pa_pd_timeseries` — per-ObsID stacked panels (PA top, PD bottom) vs time elapsed in kiloseconds from the earliest segment's `t_start_mjd` (mid-segment time `0.5*(t_start_mjd+t_stop_mjd)`, MJD→s via `×86400`), with per-segment errorbars and the catalog's weighted-mean line (`PA_deg_int`/`PD_pct_int`) plus a shaded ±1σ band (`PA_err_int`/`PD_err_int`) overlaid when the ObsID has a row in `catalog_csv`.
2. `_stokes_plane` — per-ObsID Q/I vs U/I scatter, color-coded by mid-segment elapsed time (kiloseconds, own colorbar per ObsID), 1σ error crosses from `QI_err`/`UI_err`, origin marked.
3. `_chi2_hist` — histograms of `rchi2_PD`/`rchi2_PA` across **every** row already in `catalog_csv` (not filtered to `obsids`), reference line at 1.0.
4. `_structure_function` — SF_PA (top) / SF_PD (bottom) vs lag (`sf_lags_ks`, default `[5,10,20]`), one line per ObsID with `>= sf_min_segments` segments. Recomputed via the module's own `mid_time_structure_functions` helper, which itself calls the shared `structure_function` primitive on mid-segment times (seconds) with `wrap=True` for PA and `half_width_s = 0.5*median(diff(sorted mid-times))` — never an independently-invented calculation. An ObsID short of `sf_min_segments` is omitted from this panel only (still drawn in the other three); if none qualify, the panel carries an explicit annotation instead of empty axes.

**Edge cases.** An ObsID with no readable `pcube_timeseg_<obsid>.csv` is skipped and noted in `stdout`; if every requested ObsID is unreadable, nothing is written and `exit_code=1` with `stderr` set (no `figures` key). `catalog_csv` is a **required** upstream product: if it does not exist, or exists but cannot be read as a catalog (e.g. missing the `obsid` column), this is a hard input error — `exit_code=1` with `stderr` set and no files written, checked before any figure is drawn regardless of whether the per-ObsID `pcube_timeseg` CSVs are readable. This is distinct from an ObsID simply *absent from* an otherwise-present, readable `catalog_csv`, which still gets its time-series/Stokes-plane panels — only its weighted-mean overlay is skipped.

**Determinism.** Fixed `dpi=150`, `bbox_inches='tight'`, and explicit `savefig` metadata (a fixed `Software` string) so repeated calls on unchanged input are byte-identical — no timestamps are baked into the pixels.

**Returns** a dict with `figures` (list of 4 `FigureRef` dicts — `path`, `caption`, `mime_type` (`image/png`), `produced_by`, `created_at`, `sha256` of the written bytes), `structure_function` (obsid → `{lags_ks, sf_pa, sf_pd}` for the ObsIDs actually drawn in panel 4 — the exact plotted values, exposed for differential verification), `stdout`, `exit_code`.

---

### `bhxrb_pd_hr_survey` — PD-vs-HR correlation/statistics survey

---

#### `bhxrb_pd_hr_survey.pd_hr_survey(obs_db_path, obsids_with_sources, hr_csv, pcube_csv_dir, ...)`

Correlates IXPE polarization degree (PD) against hardness ratio (HR) across a survey of BH-XRB ObsIDs, from two already-written upstream products — each ObsID's `pcube_timeseg_<obsid>.csv` integrated row (`ixpe_pcube_timeseg`) and the survey's `bhxrb_hardness_ratios.csv` (`ixpe_HR`). Same module family as `red_ixpe.HR`/`pa_pd_stability.pa_pd_stability`: pure numpy/scipy, imports neither `xspec` nor `ixpeobssim`, so `TOOL_MANIFEST` wires it directly to `xspec_utils.lib.bhxrb_pd_hr_survey` (no `skills/session.py` involvement).

```python
def pd_hr_survey(
    obs_db_path: str,
    obsids_with_sources: dict[str, str],
    hr_csv: str,
    pcube_csv_dir: str,
    mdp_flag_pct: float = 15.0,
    min_n_intra: int = 4,
    soft_state_sources: Sequence[str] = (),
    hard_state_sources: Sequence[str] = (),
) -> dict
```

**Per-obsid assembly.** For each `obsid` in `obsids_with_sources` (obsid → source-name label), `PD_pct`/`PD_err`/`MDP99` are read from `<pcube_csv_dir>/<obsid>/pcube_timeseg_<obsid>.csv`'s integrated (`segment_id=-1`) row via `_read_pcube_integrated` — an independent reimplementation of the provenance-header/`segment_id=-1` convention `pa_pd_stability._read_timeseg_csv`/`red_ixpe._read_existing_timeseg_csv` already follow (not a shared helper, so this module never imports either), and `HR`/`HR_err` from `hr_csv` via `_read_hr_rows`. `detected = PD_pct > MDP99`; `mdp_flagged = MDP99 > mdp_flag_pct` (a data-quality note, never an exclusion). An obsid with no readable pcube CSV, or absent from `hr_csv`, gets a recorded `note` and is excluded from every statistic below — never raised, never a dummy value. `hr_csv` itself missing entirely is a hard `ValueError` (a required upstream product); a single obsid's row being absent from an otherwise-present file is only a per-obsid note.

**Statistics, computed over the detected-only rows** (a non-detection's PD is MDP99-noise-dominated, not a physical measurement):
1. **Per-source Spearman** (`PD_pct` vs `HR`, `scipy.stats.spearmanr`) for every source with `>= min_n_intra` detections — a source short of that is omitted from `intra_source_correlations.csv` and `per_source_correlations`, never crashed on.
2. **Combined Spearman** across every detected point, regardless of source.
3. **Two-sample Kolmogorov-Smirnov** (`scipy.stats.ks_2samp`) between the detected `PD_pct` of `soft_state_sources` vs `hard_state_sources` — both passed by the caller, never hardcoded, so the same tool serves any source-group split.

**Statistics, computed over every obsid with both a readable PD and HR** (detected and non-detected together):

4. **Two-sided Mann-Whitney U** (`scipy.stats.mannwhitneyu`) on `HR` between the detected and non-detected groups — a detection-bias check (does HR predict whether a polarization detection was made).

**NaN, never an exception**, when a Spearman correlation's finite-pair count is `< 3` (`_spearman`), or when the KS/Mann-Whitney input group is empty (`_ks_two_sample`/`_mann_whitney`) — computed by hand rather than letting the underlying `scipy.stats` call raise on an empty array.

**Writes** (to `obs_db_path`, created if absent): `bhxrb_pd_hardness_survey.csv` (one row per requested obsid: `obsid, source, PD_pct, PD_err, MDP99, detected, mdp_flagged, HR, HR_err, note`), `intra_source_correlations.csv` (`source, n, rho, p`), `combined_correlation.csv` (`n, rho, p`), `ks_test_soft_vs_hard.csv` (`n_soft, n_hard, statistic, p`), `detection_bias_assessment.csv` (`n_detected, n_non_detected, statistic, p`). Records one `kind='pd_hardness_survey'` product per obsid with both a readable PD and HR via `obsdb_catalog.neutral.advance` (mission `'ixpe'`, level 3), pointing at the shared survey CSV with that obsid's own `source`/`PD_pct`/`MDP99`/`detected`/`mdp_flagged`/`HR`/`mdp_flag_pct` in `params`.

**Returns** a dict with `rows` (one per requested obsid, per `SURVEY_FIELDS`), `skipped` (obsids with a note), `per_source_correlations` (source → `_spearman` result, only sources meeting `min_n_intra`), `spearman_combined`, `ks_result`, `mwu_result` (each a dict per the statistics above), `survey_csv`/`intra_csv`/`combined_csv`/`ks_csv`/`bias_csv` (output paths), `exit_code` (always 0 — a per-obsid problem is recorded, not raised).

**No figures.** This tool is data-only, per its brief; `pd_hr_survey_plot` (below) is the companion plotter, mirroring how `pa_pd_stability.stability_plot` sits alongside `pa_pd_stability.pa_pd_stability`.

---

#### `bhxrb_pd_hr_survey.pd_hr_survey_plot(obs_db_path, survey_csv, ...)`

Pure-matplotlib (Agg backend) companion plotter for `pd_hr_survey`: renders the five standard reviewer diagnostic figures from `bhxrb_pd_hardness_survey.csv`. Imports neither `xspec` nor `ixpeobssim`, so — like `pd_hr_survey` itself — it stays in `lib/` and is wired directly into `TOOL_MANIFEST` (`bhxrb_pd_hr_survey_plot`), never through `skills/session.py`.

```python
def pd_hr_survey_plot(
    obs_db_path: str,
    survey_csv: str,
    save_path_prefix: str | None = None,
    soft_state_sources: Sequence[str] = (),
    hard_state_sources: Sequence[str] = (),
    mdp_flag_pct: float = 15.0,
) -> dict
```

**Writes exactly five files**, `{save_path_prefix}_<suffix>.png`:

1. `pd_vs_hr_combined` — PD_pct-vs-HR scatter (with PD_err/HR_err crosses) over every **detected** row, color-coded by source, annotated with the combined Spearman rho/p (`_spearman`, same helper `pd_hr_survey` uses — same module, reused directly).
2. `pd_vs_hr_per_source` — the same scatter split one panel per source (every source with a usable row, not just ones with detections); each panel annotated with its own Spearman rho/p when it has `>= 3` detections, else an explicit "No detected obsids" annotation.
3. `pd_distribution_comparison` — detected `PD_pct` histograms for `soft_state_sources` vs `hard_state_sources` (both caller-supplied, never hardcoded — the identical grouping argument `pd_hr_survey`'s own KS test takes), annotated with the KS statistic/p (`_ks_two_sample`). Annotated "no groups supplied or matched" instead of empty axes when neither group has a match (the default, empty-tuple case).
4. `hr_detection_bias` — HR histograms for the detected vs non-detected groups (every usable row), annotated with the Mann-Whitney U/p (`_mann_whitney`).
5. `feasibility_mdp99_exposure` — MDP99 (%) vs exposure (kiloseconds) scatter, color-coded by detection, with a `mdp_flag_pct` reference line. `survey_csv` does not carry `exposure_s`, so this panel alone reads each usable obsid's own `<obs_db_path>/<obsid>/pcube_timeseg_<obsid>.csv` integrated (`segment_id=-1`) row via a new `_read_integrated_row` helper (an independent reader returning the *whole* row, kept separate from `_read_pcube_integrated` — which `test_read_pcube_integrated_round_trip` pins to exactly `PD_pct`/`PD_err`/`MDP99` — precisely so that existing contract is untouched). An obsid whose `pcube_timeseg` CSV isn't locally readable is skipped from this panel only (noted in `stdout`), never from the other four.

**Edge cases.** `survey_csv` is a **required** upstream product (from `pd_hr_survey`): missing, or present but failing `_load_survey_rows`'s `SURVEY_FIELDS` schema check, is a hard error — `exit_code=1` with `stderr` set and no files written, checked before any figure is drawn. A `survey_csv` with **zero usable rows** (every row's `note` column is set — i.e. every requested obsid was skipped upstream, or the CSV is header-only) is the same hard error, since nothing is left to plot. Within an otherwise-usable `survey_csv`, an empty per-panel group (no detections, no soft/hard match, no readable exposure data for any obsid, ...) never degrades to a misleadingly empty plot — the affected panel(s) carry an explicit annotation instead, and `exit_code` stays `0`.

**Determinism.** Same convention as `stability_plot`: fixed `dpi=150`, `bbox_inches='tight'`, and explicit `savefig` metadata (no embedded timestamp), so repeated calls on unchanged input are byte-identical.

**Returns** a dict with `figures` (list of 5 `FigureRef` dicts), `plotted` (`combined_points`, `per_source_points`, `distribution_groups`, `detection_bias_hr`, `mdp99_exposure` — the exact values drawn in each panel, exposed for differential verification against an independent `survey_csv`/`pcube_timeseg` read), `stdout`, `exit_code`.

---

### `NICER` — NICER file locators

All three functions share the same signature pattern: `(obs: str, nicer_path: str, chatter=False)`. They glob under `<nicer_path>/<obs>/xti/event_cl/` and raise `ValueError` if exactly one match is not found.

#### `NICER.pha(obs, nicer_path)`
Locates source (`*sr.pha`) and background (`*bg.pha`) PHA files.
**Returns** `tuple[str, str]` — `(src_pha, bkg_pha)`.

#### `NICER.rmf(obs, nicer_path)`
Locates the RMF file (`*mpu7.rmf`).
**Returns** `str`.

#### `NICER.arf(obs, nicer_path)`
Locates the ARF file (`*mpu7.arf`).
**Returns** `str`.

---

### `NuSTAR` — NuSTAR file locators

All three functions share the same signature pattern: `(obs: str, nustar_path: str, chatter=False)`. They glob under `<nustar_path>/<obs>/products/` and raise `ValueError` if exactly two matches are not found (one per focal-plane module).

#### `NuSTAR.pha(obs, nustar_path)`
Locates optimally-grouped source (`*sr_grp.pha`) and background (`*bk.pha`) PHA files for FPMA and FPMB.
**Returns** `tuple[list[str], list[str]]` — `(src_phas, bkg_phas)`, each a two-element list.

#### `NuSTAR.rmf(obs, nustar_path)`
Locates RMF files for FPMA and FPMB (sorted alphabetically, giving FPMA first).
**Returns** `list[str]`.

#### `NuSTAR.arf(obs, nustar_path)`
Locates ARF files for FPMA and FPMB (sorted alphabetically).
**Returns** `list[str]`.

---

> **Adding a new instrument.** `NICER.py` and `NuSTAR.py` are deliberately thin: each is just a small set of `pha` / `rmf` / `arf` glob helpers with raise-on-ambiguity behaviour. To support a new mission (Chandra, XMM, Swift XRT, ...), drop in an analogous `<INST>.py` module exposing the same trio, then add a matching `elif inst == '<INST>':` branch in `xspec_help.load_data` that calls those helpers, sets `XFLT0001` if needed, handles energy noticing, and appends to the returned index dict. No other code needs to change.

---

### `find_src` — automated point-source centroiding

---

#### `find_src.run_extraction(evt_path, src_reg_out, bkg_reg_out, png_out, src_rad, bkg_inner_rad, bkg_outer_rad)`

Centroiding pipeline for a single FITS event file. Constructs a 1000×1000 pixel histogram of the X/Y event coordinates, identifies the brightest pixel, and refines the centroid within a 20-pixel window using `scipy.ndimage.center_of_mass`. Converts the pixel centroid to sky coordinates using either table-column WCS (NuSTAR) or standard image WCS (IXPE).

> **Point sources only.** The centroiding assumes a single dominant peak in the event image. It is not appropriate for extended sources, crowded fields, or targets that are not the brightest emitter in the FOV — in those cases, build regions by hand in DS9 instead.

```python
def run_extraction(
    evt_path: str,
    src_reg_out: str,
    bkg_reg_out: str,
    png_out: str,
    src_rad: float,           # arcminutes
    bkg_inner_rad: float,     # arcminutes
    bkg_outer_rad: float,     # arcminutes
) -> None
```

Writes:
- `src_reg_out` — DS9 `fk5; circle(RA, Dec, radius")` region
- `bkg_reg_out` — DS9 `fk5; annulus(RA, Dec, r_in", r_out")` region
- `png_out` — log-scale event image with source and background overlays

**CLI usage:**
```bash
python3 lib/find_src.py event_file.evt \
    --src_rad 1.0 --bkg_inner_rad 2.0 --bkg_outer_rad 3.0
# outputs: src.reg, bkg.reg, verify.png in the event file's directory
```

---

### `epoch` — observing epoch definitions

`epoch.epoch` is a `dict[str, dict]`. The shipped version holds 12 entries (E1–E12) for the Cyg X-1 campaign, but the structure is **fully customisable and source-agnostic**: you can rename the epoch keys, change the number of epochs, drop or add instruments, and define whatever `STATE` labels fit your target (e.g. orbital phases, outburst stages, flare/quiescent, spectral hardness bins). The loader only reads the instrument keys present in each entry, so adapting the file is the main integration point for a new source.

Each entry has the structure:

```python
{
    'IXPE':    ['XXXXXXXXX'],          # list of 9-digit obs IDs (optional)
    'NICER':   ['XXXXXXXXXX', ...],    # one or more obs IDs (optional)
    'NuSTAR':  ['XXXXXXXXXXX', ...],   # one or more obs IDs (optional)
    'STATE':   'LHS' | 'HSS' | ...,    # arbitrary user-defined label
}
```

Any subset of the instrument keys may be present; the loader silently skips absent ones. The order of keys is also the data-group ordering — see [Data-group assignment](#xspec_helpload_dataepoch-ixpe_path-).

| Epoch | State | IXPE obs | NICER | NuSTAR |
|---|---|---|---|---|
| E1 | LHS | 01002901 | yes | yes |
| E2 | LHS | 01250101 | yes | yes |
| E3 | HSS | 02008201 | — | — |
| E4 | HSS | 02008301 | — | — |
| E5 | HSS | 02008401 | yes (×2) | yes |
| E6 | HSS | 02008501 | — | yes (×2) |
| E7 | HSS | 02008601 | yes | yes |
| E8 | LHS | 03002201 | yes | yes |
| E9 | LHS | 03003101 | — | — |
| E10 | LHS | 03010001 | — | — |
| E11 | LHS | 03010101 | yes (×3) | — |
| E12 | HSS | 03002599 | — | — |

---

### `src_pos` — IXPE source coordinates

`SRC_COORD: dict[str, tuple[float, float]]` maps IXPE observation IDs to manually determined `[RA, Dec]` in degrees (J2000 / FK5). Currently covers all 12 Cyg X-1 observations. Coordinates are determined from the IXPE event data; the pi-corrected variant `01002901_picorr` has a slightly different centroid.

This catalogue is superseded by the on-the-fly centroiding in `find_src.py`: `find_src` is validated on IXPE event files by `red_ixpe.pcube_timeseg` (its centroid on 01250401 matches the source's published position and polarization within tolerance — see Changelog), and as of the `ixpe_changes=True` fix (see Changelog) `xspec_help.load_data`'s fresh-reduction path also centroids fresh via `find_src.run_extraction` rather than reading `src_pos.SRC_COORD` — the same convention `red_ixpe.reduce_pha` already used. `src_pos` remains only as a historical/manual record.

## IXPE background filtering

Instrumental background rejection follows the method described in [Di Marco et al. (IXPE background)](https://github.com/aledimarco/IXPE-background). The script `filter_inst_bkg_IXPE.sh` applies `filter_background.py` from that repository to all observations.

## Notes

- The `epoch.py` file defines 12 observing epochs (E1--E12) for Cyg X-1, spanning both the low/hard state (LHS) and the high/soft state (HSS). Adapt this file for other sources.
- `src_pos.py` stores manually determined IXPE source coordinates, kept only as a historical/manual record. `find_src.py` centroiding is validated on IXPE event files and is what `xspec_help.load_data`'s `ixpe_changes=True` path, `red_ixpe.reduce_pha`, and `red_ixpe.pcube_timeseg` (when `ra`/`dec` are omitted) all actually use. For NuSTAR, centroids are computed on-the-fly by `find_src.py`.
- IXPE Stokes Q and U PHA files are derived from the I file path by replacing `pha1` with `pha1q`/`pha1u`.
- NuSTAR background loading is disabled by default (`nustar_use_bkg=False`) because clean background regions are often unavailable for bright sources like Cyg X-1.

## Changelog

- **2026-07-10** — Added `lib/benchmark_compare.py` (`compare_to_benchmark`, `pcube_quantity`, `fit_quantity`) and the `compare_to_benchmark` AI tool, per the reproduction-study brief: worker 0's this-work-vs-paper significance table (Cyg X-1 IXPE polarimetry reproduction of Krawczynski+2022) had been computed by hand arithmetic, violating the derived-quantity/data-product workflow, and worker 1 had no tool to produce it at all. `compare_to_benchmark(records, out_csv, sigma_flag=3.0)` takes a caller-built list of records (`label`, `value_this`/`err_this`, `value_paper`/`err_paper`, a shared `unit` or the `unit_this`/`unit_paper` pair, and a mandatory non-empty `paper_source`), computes `delta=value_this-value_paper` and `z=abs(delta)/sqrt(err_this**2+err_paper**2)` (`z=inf`, never a `ZeroDivisionError`, when both errors are zero), flags every record with `abs(z)>sigma_flag`, and writes a provenance-headed CSV (`sigma_flag`/`created_at`/`produced_by` in a JSON comment header line). It hardcodes no paper numbers and re-derives no science -- every value comes from the caller. Refuses (`ValueError`) an empty `records` list, an empty/missing `paper_source`, or a genuine unit mismatch between the paired sides (`unit_this`/`unit_paper` disagreeing, or disagreeing with a supplied `unit`) -- this tool does no unit conversion itself, only verifies the caller already matched units. Two lib-only extractor helpers (never their own `TOOL_MANIFEST` entries, so the error-prone unit choice stays out of the caller, per the lead's steer) apply the identical fraction<->percent contract: `pcube_quantity(csv_path, band_index, name, unit='percent')` for `PD`/`QN`/`UN`/`PA` off an `ixpe_pcube` CSV (`PA` always returned in degrees, unchanged by `unit`), and `fit_quantity(fit_dict, name, unit='percent')` for the same contract against `fit_tbabs_polconst_powerlaw`'s (`PD`/`PA`) and `fit_tbabs_pollin_powerlaw`'s (`PD_A`/`PD_slope`/`PA_A`/`PA_slope`) result dicts. `pcube_quantity` reimplements `red_ixpe.read_pcube_csv`'s provenance-header/band-field contract locally (same header prefix, field names, and hard-`ValueError` behaviour) rather than importing `red_ixpe` (which pulls in `ixpeobssim` at module top level), matching `ixpe_pcube_plot.py`/`bhxrb_pd_hr_survey.py`/`pa_pd_stability.py`'s established convention -- the whole module stays import- and call-clean without the science stack, wired directly into `TOOL_MANIFEST`. Validated on the real, already-local `ixpe_pcube_01002901.csv` (`pcube_quantity`: `PD`/`QN`/`UN` percent values are exactly `100x` their fraction counterparts, `PA` is unchanged across units, both checked against an independent direct CSV read) and a genuine `fit_tbabs_polconst_powerlaw` fit against real Cyg X-1 IXPE ObsID `01002901` Level-2 data (`fit_quantity`: same contract on real fitted `PD`/`PA`). `compare_to_benchmark` itself is tested with hand-built records (no new science to validate): a capability test confirms a zero-sigma record is never flagged and an exactly-5-sigma record is, with every `z` cross-checked against an independent in-test `abs(delta)/sqrt(...)` recompute; a differential test round-trips the written CSV byte-for-byte against the returned rows; a reproducibility test confirms two calls on identical records agree exactly, including the `z=inf`/flagged edge case; edge-case tests cover empty records, empty `paper_source`, and mismatched `unit_this`/`unit_paper`.

- **2026-07-10** — Generalized `xspec_help.load_data`'s IXPE branch from a single ObsID to an inner `for obs in epoch['IXPE']` loop of arbitrary length N (>= 1), mirroring the NuSTAR/NICER branches' own inner ObsID loop exactly (previously the lone exception, with an explicit `NotImplementedError` guard for `len(epoch['IXPE']) > 1`, now removed). Needed for Krawczynski+2022's Cyg X-1 single-session tied-`polconst` cross-check, which ties PD/PA across *all* DUs of *both* the May (`01002901`) and June (`01250101`) IXPE ObsIDs in one `tbabs*polconst*powerlaw` fit — previously blocked since `load_data` could never produce a multi-ObsID `indices['IXPE']` for `fit_tbabs_polconst_powerlaw`/`fit_tbabs_pollin_powerlaw` to consume. Per ObsID in the loop: the source centroid (`find_src.run_extraction`, when `ixpe_changes=True`), per-spectrum `RESPFILE`/`ANCRFILE`/`BACKSCAL`, and the auto-resolved IRF epoch (`irf_epoch.resolve_irf_date`, via `red_ixpe.pha`) all resolve independently — nothing is shared, cached, or reused across ObsIDs (confirmed: 01002901's May midpoint and 01250101's June midpoint are each fed to `red_ixpe._resolve_irf_date` on their own, even though both happen to resolve to the same nearest installed epoch, `20220702`, for this particular pair). Data groups are assigned consecutively and gaplessly across every ObsID (3 DGs each; `epoch['IXPE']=[A,B]` yields DGs 1-3 for A, 4-6 for B); DG1 remains an IXPE I/Q/U group, unchanged. **Gain convention (settled with the human, Option A):** within each (ObsID, DU) triplet, Q/U's gain slope link to that same (ObsID, DU)'s own I gain-slope parameter, computed strictly from the current triplet's own running spectrum index (`AllData(du_specs[0])`, `du_specs` scoped per DU-loop iteration) — never a fixed offset, so a later ObsID's Q/U can never mislink back to an earlier ObsID's I; no gain parameter is ever tied *across* ObsIDs, matching how NuSTAR/NICER already treat a multi-ObsID epoch (no cross-obs calibration tie at all). N=1 is unchanged (byte-for-byte the same code path, just wrapped in a single-iteration loop). `TOOL_MANIFEST`'s `load_data` entry description updated to document the multi-ObsID IXPE capability and gain convention; no `args_schema` change (the existing `epoch` argument already accepted an obs-ID list). `tool_config.yaml`'s `load_data` `call_timeout_s` raised 900s -> 2700s, since the real runtime now scales with ObsID count rather than being capped at one fresh reduction. Validated on the real, already-local validation pair (`01002901`, `01250101`), copied into isolated tmp directories per fixture: a structure test confirms `AllData.nGroups == 3*len(epoch['IXPE'])`, gapless consecutive DG numbering, and both ObsIDs present in `indices['IXPE']`, parameterized on N rather than the literal 6; a capability test walks the actual `AllModels`/`AllData` gain links and confirms every (ObsID, DU)'s Q/U resolves to that same (ObsID, DU)'s own I gain-slope parameter and never to a spectrum owned by a different ObsID (decoding XSPEC's global `p<N>` gain-parameter numbering, since gain is set with no `Model()` ever defined); a differential test confirms the multi-ObsID load's per-DU net Stokes-I rate and resolved IRF epoch match a separate, independent single-ObsID `load_data` call for each ObsID (to `rel=1e-6` on the rate; exact match on the epoch string read off each spectrum's own `RESPFILE` path); a second differential test confirms each resolved epoch independently against `red_ixpe.get_mjd`+`irf_epoch.resolve_irf_date` (the trusted primitives, not `load_data`'s own code path); a reproducibility test confirms two independent multi-ObsID `load_data` calls agree on every per-DU net rate (`rel=1e-6`) and resolved epoch; a regression test confirms each single-ObsID (N=1) load still gets the pre-existing DG/gain-link behaviour unchanged. The existing `test_guardrail_multi_obsid_raises_cleanly` test (asserting the now-removed guard fired) was rewritten to `test_guardrail_multi_obsid_no_longer_raises_not_implemented`, confirming instead that a 2-ObsID IXPE list reaches (and is stopped only by) `red_ixpe.get_all_file_list`'s own "no files found" check on a deliberately-empty `ixpe_path` — proving the `NotImplementedError` guard is gone without touching the shared XSPEC session `test_fit_tbabs_polconst_powerlaw.py`'s other tests rely on.

- **2026-07-09** — Extended `ixpe_pcube_plot.ixpe_pcube_plot`/the `ixpe_pcube_plot` AI tool with an optional `pollin` PD(E)/PA(E) model overlay, per the reviewer-figure brief for the energy-dependence claim (Krawczynski+2022 Fig. S5 style: the fitted model line drawn over the energy-binned PCUBE points). Added `pollin: dict | None = None`: when given (e.g. `fit_tbabs_pollin_powerlaw`'s result dict, extra keys ignored), draws the analytic PD(E)=`PD_A`+`PD_slope`*(E-1 keV)/PA(E)=`PA_A`+`PA_slope`*(E-1 keV) line on the PD/PA panels; when omitted, falls back to an inverse-variance-weighted linear fit of the CSV's own PD(E)/PA(E) points via a new `_weighted_linear_fit` helper (`scipy.optimize.curve_fit`, `sigma=`the column's own `*_err`, `absolute_sigma=True`) — the human's A/B-else-fallback design, so the figure always carries a line unless too few bands survive to fit one. Both paths draw the same way and are reported under a new `plotted['overlay']` key (`{'kind': None}` when nothing was drawable). A single-band CSV still renders the bare point with no fallback line (handled gracefully, not an error); `pollin` missing any of `PD_A`/`PD_slope`/`PA_A`/`PA_slope` raises `ValueError` naming the missing keys. No change to the existing PD(E)/PA(E)/Q-U-plane point-plotting behaviour, its `ValueError` contract, or its byte-determinism convention. Validated on the real, already-local Cyg X-1 IXPE ObsID `01002901` (this brief's validation ObsID) via the real `ixpe_pcube` skill: a capability test confirms the analytic `pollin` line is actually drawn (line data extracted from the PD/PA axes' `Line2D` artists and checked against `PD_A+PD_slope*(E-1)`/`PA_A+PA_slope*(E-1)` at the drawn energies) and that the fallback line is drawn instead when `pollin` is omitted, plus the existing exactly-one-PNG/sha256 checks; a differential test confirms `_weighted_linear_fit`'s slope/intercept (and their errors) against an independent closed-form weighted-least-squares matrix computation (`(XᵀWX)⁻¹XᵀWy`, `W=diag(1/σ²)`) on the same PD(E) points, not the plotter's own helper; a reproducibility test confirms two calls on the same CSV give byte-identical fallback-fit parameters (deterministic `curve_fit` on unchanged input). No new `tool_config.yaml` entries needed — the added work is pure numpy/scipy on already-read CSV columns, well inside the existing test-timeout budget.

- **2026-07-09** — Added `xspec_help.fit_tbabs_pollin_powerlaw` and the `fit_tbabs_pollin_powerlaw` AI tool: the paper's *actual* linear energy-dependence route for Cyg X-1 polarization (an XSPEC `pollin` fit), per the brief flagging that worker 0's reported linear PD(E) numbers (constant `A=3.187%±0.140%`, linear `A=1.724%±0.341%`, `B=0.534%/keV±0.114%/keV`) came from an inline `runfit.py` `scipy` fit on `ixpe_pcube`'s PD(E) points — a data-analysis boundary violation, not toolkit-grounded. Fits `tbabs*pollin*powerlaw` to every IXPE Stokes I/Q/U triplet in a `load_data` indices dict, tying `pollin`'s PD(E)=`A1`+`Aslope`*(E-1 keV) (fraction, fraction/keV) and PA(E)=`psi1`+`psislope`*(E-1 keV) (degrees, degrees/keV) into **one free (A1, Aslope, psi1, psislope) quartet for the whole session**, via the exact same default cross-data-group parameter tying `fit_tbabs_polconst_powerlaw` already relies on (confirmed empirically here too — no explicit tie/untie call needed). Same load path, same guard-rails, and the same covariance-matrix `Parameter.sigma`-scaled-to-`conf_level` error convention as `fit_tbabs_polconst_powerlaw` (`RuntimeError` naming the parameter if PD_A/PD_slope/PA_A/PA_slope/Gamma has no valid sigma). Validated on the real, already-local Cyg X-1 IXPE ObsID `01002901` via a fresh `load_data(ixpe_changes=True)` reduction at the paper's region radii (`src_rad=1.333`, `bkg_inner_rad=2.5`, `bkg_outer_rad=5.167`): the fit converges to `chi2/dof=4482.47/1332=3.37` (elevated for the same reason `fit_tbabs_polconst_powerlaw`'s global chi2 is — `pollin` multiplies only Q/U, so the Stokes-I continuum's poor description by a bare `powerlaw` does not compromise PD/PA), `PD_A=1.75%±0.37%`, `PD_slope=0.614%/keV±0.128%/keV` (positive, consistent in sign and order of magnitude with the paper's reported rise), `PA_A=-14.7°±3.2°`, `PA_slope=-1.94°/keV±0.94°/keV`. A capability test confirms exactly one free parameter of each of A1/Aslope/psi1/psislope across the whole session and that background was genuinely subtracted (net Stokes-I rate < gross); a second capability test confirms the linear term is genuinely exercised — freezing `Aslope` at 0 and refitting gives a nonzero Delta-chi2 relative to the reported slope-free fit; a differential test confirms `PD_slope`'s sign and rough (factor-of-4) magnitude agree with an independent inverse-variance-weighted linear regression of `ixpe_pcube`'s own PD(E) points (not this helper's own fit code); a reproducibility test confirms two calls on identical already-loaded fresh-reduced spectra agree on PD_A/PD_slope/PA_A/PA_slope/chi2 within 0.1%. No new `tool_config.yaml` entries needed — like `fit_tbabs_polconst_powerlaw`, this is a fit-only call on already-loaded spectra, not a fresh xpselect/xpbin reduction, and comfortably finishes well inside the existing package-wide `gate_test_timeout_s=1500`.

- **2026-07-09** — Simplified `xspec_help.fit_tbabs_polconst_powerlaw`'s error interface: removed the XSPEC `error`-command machinery entirely (no more `Fit.error(...)` call, `try/except`, or `error_method` bookkeeping) in favor of a single, always-taken path -- the covariance-matrix parameter-table 1-sigma (`Parameter.sigma`, scaled to `conf_level` via `sqrt(chi2.ppf(conf_level, 1))`). Diagnosed first (no code changes) why `error_method` always came back `'covariance'` on real Cyg X-1 data: `Fit.error` refuses outright via its documented `maximum`-reduced-chi2 guard (default 2.0; [HEASARC manual](https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSerror.html)) whenever reduced chi2 exceeds it -- confirmed a legitimate, hardcoded XSPEC safeguard rather than a mis-call (`Fit.query="no"` made no difference; raising `maximum` unblocked the command). A follow-up diagnostic pass (gamma_init swept over {2.5, 2.0, 1.6, 1.2, 0.8}, extra `Fit.perform()` iterations at `nIterations=1000`) confirmed Gamma=1.609 is a true, unique global minimum (not an early-stopping or local-minimum artifact), and a per-spectrum chi2 decomposition showed the elevated global reduced chi2 (~3.4) sits almost entirely in the three Stokes-I spectra (reduced chi2 ~5-11) while every Q/U spectrum fits at reduced chi2 ~0.9-1.15 -- since `polconst` multiplies only Q/U, PD/PA are set by the well-fit polarization sector regardless of the I continuum's poor description by a bare `powerlaw`. Given that, the `error` command's extra rigor (asymmetric bounds, non-linear parameter search) was buying nothing our fixed Cyg X-1 chi2 regime could reliably use, while its refuse/fallback dance added a code path with no test coverage of its own status-string edge cases; the lead settled on trimming to the single covariance path, kept honest by raising `RuntimeError` (naming the parameter) if any of PD/PA/Gamma has no valid covariance sigma (`Parameter.sigma <= 0` or unavailable), instead of ever returning a meaningless 0/negative error. Output contract trimmed to drop `error_method` (kept `conf_level`, since it still documents what the returned errors mean). `test_capability_output_contract` updated to match; `test_differential_polconst_matches_pcube_broadband`'s tolerance tightened from 3-sigma to 2-sigma (observed agreement is <0.4 sigma, so 2-sigma stays safe while being a stricter regression guard).

- **2026-07-09** — Added `xspec_help.fit_tbabs_polconst_powerlaw` and the `fit_tbabs_polconst_powerlaw` AI tool: the polconst spectropolarimetric cross-check of the standard PCUBE polarization (the Krawczynski+2022 Cyg X-1 XSPEC column), fitting `tbabs*polconst*powerlaw` to every IXPE Stokes I/Q/U triplet in a `load_data` indices dict with `polconst`'s PD (`A`)/PA (`psi`) tied into one free pair for the whole session — never per-DU or per-ObsID — which XSPEC does by default for a model defined once across several already-active data groups, so no explicit tie/untie call was needed (unlike `fit_tbfeo_hss`'s cross-calibration constant, which explicitly unties). A prior manual `tbabs*polconst*powerlaw` attempt returned `chi2/dof=64.1` and an unphysical `PD=0.85%`, traced to `load_data(ixpe_changes=False)` reusing a stale on-disk PHA1 with dead `RESPFILE` paths and a non-paper-region `BACKSCAL`; this helper instead requires a fresh reduction (`load_data(ixpe_changes=True, ixpe_src_rad=1.333, ixpe_bkg_inner_rad=2.5, ixpe_bkg_outer_rad=5.167)` for the paper's regions), which resolves both correctly. Runs `Fit.perform()` to a stable minimum, then the XSPEC `error` command at `conf_level` (default 0.683 = 1 sigma, converted via the chi-squared quantile for 1 dof); when `error` refuses on a statistically poor fit (`Reduced Chi^2 > 2`, a real XSPEC safeguard — real Cyg X-1 IXPE data does not always reach `chi2/dof<2` for this simple a model), `PD_err`/`PA_err`/`Gamma_err` fall back to the covariance-matrix `Parameter.sigma` scaled by `sqrt(delta_chi2)`, so a real number is always returned and `error_method` records which path was used (the manual worker script's own fallback used `Parameter.values[1]` — the fit step size, not the uncertainty — an error this helper does not repeat). Every fitted number is returned in the result dict, since XSPEC's own terminal output is muted in the sandbox.

  **Companion fix (uncovered while authoring, confirmed with the lead before making):** `xspec_help.load_data`'s IXPE `ixpe_changes=True` branch had never actually been reachable — it called `red_ixpe.pha(..., ra=None, dec=None, changes=True)`, and `pha` raises `ValueError` whenever `changes=True` and `ra`/`dec` are `None`. Fixed by centroiding the source position via `find_src.run_extraction` on the DU1 Level-2 event file first (mirroring `red_ixpe.reduce_pha`'s existing fresh-reduction path exactly — never `src_pos.SRC_COORD`, a legacy manual table `ixpe_reduce_pha` already avoids), then passing the resolved `ra`/`dec` into `pha`. `ixpe_changes=False` (reuse on-disk products) is unchanged. Also extended the AI-facing `research_toolkit.xspec_utils.load_data` `TOOL_MANIFEST` entry's `args_schema` to expose `ixpe_changes`/`ixpe_src_rad`/`ixpe_bkg_inner_rad`/`ixpe_bkg_outer_rad` (previously absent — since `load_data`/`fit_*`/`plot_fit` are driven entirely through sequential `TOOL_MANIFEST` tool calls sharing the persistent XSPEC session, with no separate REPL/script path exposed to a tool-calling analyst, this omission made the fresh-reduction path unreachable by any analyst regardless of the lib-level fix). `tool_config.yaml` gained a `load_data: 900` `call_timeout_s` entry (matching `ixpe_reduce_pha`'s budget for the same find_src + xpselect/xpbin PHA1/Q/U cost), since this path previously failed instantly and now does real work.

  Validated end-to-end on real, already-local Cyg X-1 IXPE ObsIDs `01002901`, `01250101`, `02008401`: a fresh `load_data(ixpe_changes=True)` reduction at the paper's region radii succeeds (previously impossible) and assigns the correct per-Stokes response (arf for I, mrf for Q/U, confirmed in the XSPEC log); a differential test confirms the auto-resolved centroid matches an independent `find_src.run_extraction` call on the same DU1 file to `1e-6` degrees; a capability test confirms background was genuinely subtracted (net Stokes-I rate < gross, read via PyXspec's actual `Spectrum.rate` tuple order — `(net, variance, total, model)`, not the order `xspec_help.audit_paths`'s own inline comment claims) and that exactly one PD and one PA parameter are free across the session (`AllModels(dg).polconst.A.link`/`.frozen` inspected per group); a differential test confirms the fitted broadband (2-8 keV) PD/PA agree with an independently-computed `ixpe_pcube` broadband PD/PA for the same ObsID/region within combined error (for 01002901: PD=3.36%±0.14% vs pcube's 3.45%±0.19%; PA=-20.30±1.21 deg vs pcube's -20.66±1.61 deg); a reproducibility test confirms two calls on identical already-loaded fresh-reduced spectra agree on PD/PA/Gamma/chi2 within 0.1%. Manually validated for robustness on all three ObsIDs outside the checked-in suite (to keep the gate fast): 01250101 gave `chi2/dof=1.75` (the `error` command succeeded directly); 02008401 gave `chi2/dof=3.93` (covariance fallback used); both PD (2-4%) and PA (-15 to -25 deg) are consistent with Cyg X-1's known polarization. `plot_fit` visually confirmed a converged fit (folded model tracking the data, residuals bounded within the noticed band) before declaring the tool done.

- **2026-07-09** — Added `lib/ixpe_pcube_plot.py` (`ixpe_pcube_plot`) and the `ixpe_pcube_plot` AI tool: the companion figure `ixpe_pcube`'s own changelog entry (below) flagged as the natural next chain link, so the reviewer can see the energy-resolved polarimetry `ixpe_pcube` produces (Krawczynski+2022 Fig. S5 PD(E)/PA(E) panels plus the Q-U plane). Pure matplotlib (Agg backend); imports neither `xspec` nor `ixpeobssim` -- a small independent reimplementation of `red_ixpe.read_pcube_csv`'s provenance-header contract (mirroring the same pattern `pa_pd_stability`/`bhxrb_pd_hr_survey` already use to avoid pulling in `red_ixpe`) keeps this module and its entry point standalone, wired directly into `TOOL_MANIFEST`, never through `skills/session.py`. Writes exactly one PNG (three panels: PD(E) top-left, PA(E) bottom-left sharing the energy axis, QN-vs-UN Q-U plane on the right color-coded by band-centre energy, 1-sigma error bars throughout), with every plotted value read straight from the CSV's own columns -- never recomputed -- and also returned under `plotted` for differential verification. A band with any non-finite plotted value is skipped, never drawn as zero; a single remaining band still renders as a lone marker (no line to connect); an empty, headerless, or all-NaN CSV **raises** `ValueError` naming the problem rather than silently producing an empty figure -- the brief's own edge-case contract, distinct from `stability_plot`/`pd_hr_survey_plot`'s exit_code=1 convention for a missing/malformed upstream CSV. Same byte-determinism convention as `stability_plot`/`pd_hr_survey_plot` (fixed `dpi`, `bbox_inches`, `savefig` metadata with no embedded timestamp). Validated on the real, already-local Cyg X-1 IXPE ObsIDs `01002901` (May 2022, default 4-band `ebinning`) and `01250101` (June 2022, single-edge-pair broadband case) via the real `ixpe_pcube` skill (genuine `xpselect`/`xpbin` reductions -- every CSV under test comes from the real pipeline, none hand-written): capability tests confirm exactly one non-empty, correctly-`sha256`'d PNG for both the multi-band and single-band case; a differential test confirms the plotted band-centre energies/PD/QN arrays equal `ixpe_pcube`'s own CSV columns read back independently; edge-case tests confirm a blanked-out band is skipped (not zeroed), an all-NaN CSV/empty CSV/header-only CSV/missing path each raise a clear `ValueError`; a reproducibility test confirms two calls on the same CSV are byte-identical.

- **2026-07-09** — Surfaced `red_ixpe.pcube` as the `ixpe_pcube` AI tool, per the Krawczynski+2022 Cyg X-1 energy-resolved-polarimetry (Fig. S5/Table S2) and native multi-ObsID-combination brief: the prior worker had only `ixpe_pcube_timeseg` (a fixed single 2-8 keV band) to mis-call for both the energy-resolved and combined-ObsID cases. Extended `pcube` itself (no other caller existed) so `obs` accepts either a single ObsID or a **list** of ObsIDs -- a list is filtered/DU-validated per ObsID via the existing `_validate_du_files`, concatenated, and passed into **one** `xpselect`/`xpbin PCUBE` call, so the returned cube is the native directly-combined reduction (raw Stokes I/Q/U/COUNTS summed across every DU *and* ObsID before PD/PA are derived), never an inverse-variance approximation built from separately-computed per-ObsID PD/PA -- this falls out for free from `xBinnedPolarizationCube.from_file_list`'s own `__iadd__` reduction over whatever files `xpbin` returns. `pcube` also gained an `ebinning`-length guard (`>= 2` edges) and now returns the `irf_date` actually used (`_resolve_irf_date` extended to average the earliest-start/latest-stop MJD span across every listed ObsID, reducing to the prior single-ObsID behaviour unchanged). Added `pcube_bands` (per-band PD/PA/Stokes/MDP99/SIGNIF via `polarization()`/attribute access -- never hand-rolled Stokes arithmetic) and `write_pcube_csv`/`read_pcube_csv` (a provenance-headed CSV, one row per band, no separate "integrated" sentinel row -- a 2-edge `ebinning` naturally yields that single broadband row instead). While building the wrapper (which the brief specifies must resolve `file_list` via `get_all_file_list(obs_db_path)` rather than a narrow per-ObsID glob), found and fixed a latent `_get_obs_file_list` defect this exposed: its raw substring match on the mission-wide file list collided plain `01002901` with `01002901_picorr` (a differently-processed variant whose filenames still embed the same 9-digit ObsID), making `_validate_du_files` see two files per DU and raise on a perfectly good request -- `_get_obs_file_list` now matches by each file's own observation directory instead. `ixpe_pcube` (`skills/session.py`) requires `ra`/`dec` explicitly (no auto-centroid, unlike `ixpe_pcube_timeseg` -- a combined request especially needs one shared region centre), auto-selects the IRF epoch (not AI-facing), returns `PD`/`PD_err`/`QN`/`QN_err`/`UN`/`UN_err`/`MDP99` as dimensionless fractions (not percent, documented explicitly since this differs from `ixpe_pcube_timeseg`'s `PD_pct`/`MDP99` convention) and `PA`/`PA_err` in degrees, writes a single ObsID's CSV under `<obs_db_path>/<obsid>/` and a combined request's CSV at the mission root (mirroring `bhxrb_hardness_ratios.csv`'s placement), and records one `kind='pcube'` product per contributing ObsID via `obsdb_catalog.neutral.advance` (level 3). Data-only tool, no figures, per the brief. `tool_config.yaml` gained a `call_timeout_s` entry (1800s, above `ixpe_reduce_pha`'s single-ObsID 900s, since a combined request's single `xpselect`/`xpbin` pair scales with total combined event volume). Validated on the real, already-local Cyg X-1 IXPE ObsIDs `01002901` (May 2022) and `01250101` (June 2022) -- the brief's own example combination: a capability test confirms the combined cube's raw per-DU-per-ObsID FITS products (independently read straight off disk, bypassing `pcube_bands`/`xBinnedPolarizationCube`) sum to the combined total, proving the merge used both ObsIDs' real events; a differential test recomputes PD/PA from the tool's own returned QN/UN (`PD=sqrt(QN^2+UN^2)`, `PA=0.5*atan2(UN,QN)`) against the tool's separately-returned PD/PA; a reproducibility test confirms two independent single-ObsID calls agree on every band's PD/PA within tolerance; a dedicated regression test confirms the `_get_obs_file_list`/`01002901_picorr` fix directly against the real, unmodified mission directory (read-only).

- **2026-07-07** — Added `bhxrb_pd_hr_survey.pd_hr_survey_plot` and the `bhxrb_pd_hr_survey_plot` AI tool: the five prescription figures the PD-hardness survey brief calls for (`pd_vs_hr_combined`, `pd_vs_hr_per_source`, `pd_distribution_comparison`, `hr_detection_bias`, `feasibility_mdp99_exposure`), consuming `bhxrb_pd_hardness_survey.csv` (written by `pd_hr_survey`, added earlier the same day) — the companion plotter that tool's own changelog entry flagged as the natural next chain link, following the established `pa_pd_stability`/`stability_plot` data/plot split exactly: pure matplotlib (Agg backend), no `xspec`/`ixpeobssim` import, wired directly into `TOOL_MANIFEST` (never through `skills/session.py`), one `FigureRef` per file with `path`/`caption`/`mime_type`/`produced_by`/`created_at`/`sha256` matching the written bytes. Lives in the same module as `pd_hr_survey` (mirroring how `stability_plot` lives alongside `pa_pd_stability` in one file) and reuses its `_spearman`/`_ks_two_sample`/`_mann_whitney` helpers directly rather than reimplementing them. The feasibility panel needs `exposure_s`, which `survey_csv` does not carry, so it reads each usable obsid's own `pcube_timeseg_<obsid>.csv` integrated row via a new `_read_integrated_row` helper — kept independent from the existing `_read_pcube_integrated` (whose return value a regression test pins to exactly `PD_pct`/`PD_err`/`MDP99`) rather than widening that contract. `soft_state_sources`/`hard_state_sources` (the distribution-comparison grouping) and `mdp_flag_pct` (the feasibility reference line) are caller arguments, not hardcoded, matching `pd_hr_survey`'s own overridable knobs, per the brief's requirement that every survey-CSV-exposed override stay overridable in the plotter too. `survey_csv` missing, wrong-schema, or carrying zero usable rows (every row's `note` set) is a hard `exit_code=1` with no files written — mirroring `stability_plot`'s `catalog_csv`/all-obsids-unreadable contract exactly; within an otherwise-usable CSV, an empty per-panel group (no detections, no soft/hard match, no obsid with a readable `pcube_timeseg` CSV) is annotated rather than rendered as a misleading empty plot. Validated on the real, already-local four-ObsID Cyg X-1 survey CSV built by running the real `ixpe_HR` and `pd_hr_survey` chain (same fixture pattern as `test_bhxrb_pd_hr_survey.py`): capability tests confirm exactly five named, non-empty, correctly-`sha256`'d files; differential tests confirm the plotted PD/HR scatter points, the detected/non-detected HR split (cross-checked against a direct `scipy.stats.mannwhitneyu` recompute), and the MDP99/exposure feasibility points all match an independent direct CSV read — never the plotter's own helpers; a reproducibility test confirms two calls are byte-identical.

- **2026-07-07** — Added `lib/bhxrb_pd_hr_survey.py` (`pd_hr_survey`) and the `bhxrb_pd_hr_survey` AI tool: correlates PD (each ObsID's `pcube_timeseg_<obsid>.csv` integrated `segment_id=-1` row, from `ixpe_pcube_timeseg`) against HR (`bhxrb_hardness_ratios.csv`, from `ixpe_HR`) across a survey of BH-XRB IXPE ObsIDs — the survey step downstream of both tools per the BH-XRB polarization-hardness campaign brief. Classifies `detected=PD_pct>MDP99` and flags `mdp_flagged=MDP99>mdp_flag_pct` (default 15.0); over the detected-only rows, computes per-source Spearman (`>= min_n_intra` detections, default 4, else omitted rather than crashed on), the combined Spearman across all detected points, and a two-sample KS test between caller-supplied `soft_state_sources`/`hard_state_sources` groups (never hardcoded); over every obsid with both a readable PD and HR, a two-sided Mann-Whitney U test on HR between detected and non-detected groups. Same module family as `red_ixpe.HR`/`pa_pd_stability.pa_pd_stability`: pure numpy/scipy, no `xspec`/`ixpeobssim` import, wired directly into `TOOL_MANIFEST`. `_read_pcube_integrated` is an independent reimplementation of the provenance-header/`segment_id=-1` reader convention (never imports `pa_pd_stability`/`red_ixpe`). An obsid with no readable pcube CSV, or missing from `hr_csv`, is recorded with a `note` and excluded from every statistic rather than raised; every statistic is NaN (never an exception) on an empty input group or a Spearman pair-count `< 3`. Writes five CSVs (`bhxrb_pd_hardness_survey.csv`, `intra_source_correlations.csv`, `combined_correlation.csv`, `ks_test_soft_vs_hard.csv`, `detection_bias_assessment.csv`) and records one `pd_hardness_survey` product per successful obsid via `obsdb_catalog.neutral.advance` (level 3). Validated on four real, already-local Cyg X-1 IXPE ObsIDs (01002901, 01250401, 02250901, 03005101): their real integrated rows give a genuine 3-detected/1-non-detected split (PD>MDP99 for the first three, PD<MDP99 for 03005101), independently confirmed against a direct CSV read before writing the tests; HR came from running the real `ixpe_HR` skill against each ObsID's real, already-local per-DU PHA1 products, copied into an isolated scratch mission directory so the shared mission `catalog.json` is never touched by the test suite. Differential tests confirm the combined Spearman/KS/Mann-Whitney statistics against direct, independent `scipy.stats` recomputation (not the tool's own helpers); a reproducibility test confirms two calls agree exactly (pure arithmetic, no fitting). No figures — a companion plotter is the natural next addition, per the brief's own "figures are the next chain link" note.

- **2026-07-07** — Fixed a second `red_ixpe.HR` silent-wrong-result mode found in co-lead review: `file_list` carrying more than one Level-2 event file for the same DU (two processing versions, e.g. `..._det1_evt2_v02.fits` and `..._det1_evt2_v03.fits`, or a mixed plain+`_rej` set) produced duplicate `du_labels`, so `HR` zipped over the duplicated PHA paths and double-counted that DU into an inflated `soft_rate`/`hard_rate` (and a `dus_used` list carrying a repeated label) -- silently, no error. `HR` now detects duplicate DU labels immediately after parsing `du_labels` (before any PHA existence check or read, regardless of `bkg`) and raises `ValueError` naming the duplicated DU(s) and the offending event files; it never auto-selects a "preferred" version -- that choice belongs to the caller/wrapper. Two regression tests added: `test_hr_raises_on_duplicate_du_versions` (a manufactured second DU1 version fed directly to `HR` raises) and `test_skill_wrapper_glob_duplicate_du_versions_recorded_as_error` (a real event_l2 directory carrying two DU1 versions, with no `_rej` variant so the `ixpe_HR` wrapper's plain-glob fallback matches both, is caught by the new guard and recorded in the survey's `errors` dict -- no inflated row, and the rest of the survey still succeeds).

- **2026-07-07** — Fixed `red_ixpe.HR` and added the `ixpe_HR` AI tool, per three defects the analyst-env owner reported against the prior implementation: (1) DU3 was silently dropped from the summed rate — the prior `_sum_counts` loaded each DU's PHA into XSPEC, incrementing `spec_cnt` to 3, then summed over `range(1, spec_cnt)` = `range(1, 3)` = spectra 1-2 only, always omitting the last-loaded DU; the rewrite reads PHA1 `RATE`/`STAT_ERR` directly (no XSPEC `AllData` at all) and attempts DU1/DU2/DU3 independently via `pha(..., changes=False)`'s existing file locators, summing over every DU actually present in `file_list` (a DU absent from `file_list` is dropped gracefully — IXPE can fly with a DU off — but a DU present in `file_list` with no PHA1 product on disk now raises, never silently drops); (2) `HR` was `H/S` — now the community-normalized `HR=(H-S)/(H+S)`, bounded `[-1,+1]` (`hr_from_bands`); (3) the returned error was computed from PyXSPEC's `Spectrum.rate` tuple misread as `(rate, variance)` when it is actually `(total, background, net, net_error)` — replaced with the closed-form first-order Gaussian propagation `sigma_HR=(2/(H+S)^2)*sqrt(S^2*sigma_H^2+H^2*sigma_S^2)` from each band's PHA-native `RATE`/`STAT_ERR`. `soft_band`/`hard_band` resolve to channels via each PHA's own `RESPFILE` RMF `EBOUNDS` (never hardcoded); discovered while validating this against the real, already-local Cyg X-1 IXPE PHA1 products that 11 of 13 already-reduced ObsIDs (e.g. 01002901) carry a `RESPFILE` baked in under a now-removed `~/.local` ixpeobssim install from an earlier reduction environment, so `_resolve_rmf_ebounds` falls back to an installed CALDB RMF for the same DU when the literal path doesn't resolve — confirmed byte-identical to every other IXPE IRF date/DU's EBOUNDS (a fixed instrumental grid; only the RESPONSE matrix varies), with a guard that raises on any channel-count mismatch rather than silently accepting a different grid. Validated on real, already-local IXPE ObsIDs 01250401 (all 3 DUs, confirms dropping DU3 changes the summed rate) and 01002901 (the stale-`RESPFILE` fallback case): differential tests confirm channel selection and the propagated error against independent astropy/quadrature recomputation; a capability test confirms the fallback grid and resulting band rates match the working-`RESPFILE` case; a reproducibility test confirms two calls agree exactly (pure FITS reads, no fitting). The `ixpe_HR` wrapper surveys a list of ObsIDs, writes `bhxrb_hardness_ratios.csv` (columns `obsid, soft_rate, soft_err, hard_rate, hard_err, HR, HR_err`) at the mission's `obs_db_path` root, and records one `hardness_ratio` product per successful ObsID via `obsdb_catalog.neutral.advance` (level 3); an ObsID that fails is recorded in the returned `errors` dict and excluded from the CSV, never a dummy row.

- **2026-07-06** — Added `red_ixpe.reduce_pha` and the `ixpe_reduce_pha` AI tool: surfaces the existing `red_ixpe.pha()` reduction (previously used only internally by `xspec_help.load_data`) as its own standalone tool, per the BH-XRB polarization-hardness survey brief's design (PHA reduction kept separate from an HR tool, since HR reads PHA). For each ObsID, reuses on-disk per-DU Stokes-I PHA1 src+bkg products when present (no re-reduction); otherwise centroids the source position from the Level-2 DU1 event file via `find_src.run_extraction` (writing DS9 `src.reg`/`bkg.reg` -- per the human's steer, coordinates come from `find_src`, never from FITS headers or the Cyg-X-1-only `src_pos.SRC_COORD` table) and runs `xpselect`+`xpbin PHA1` with a circular source region (`src_rad`, default 1.0 arcmin) and annular background (`bkg_inner_rad`-`bkg_outer_rad`, default 3.0-5.0 arcmin) at the caller's `irf_date`. Every PHA1 product on disk after the call is recorded in the mission catalog at level 3 via `obsdb_catalog.neutral.advance`, whether reused or freshly reduced. `_validate_du_files` gained an optional `caller` parameter (default `'pcube_timeseg'`, unchanged for existing callers) so its raised messages name the actual entry point instead of always crediting `pcube_timeseg`. Validated on real, already-local IXPE ObsID 01002901 (one of the twelve Cyg X-1 ObsIDs with PHA already on disk): the reuse path resolves the real on-disk PHA1 products (non-zero net `RATE`, not just well-shaped metadata); a differential test confirms the fresh-reduction path's centroided RA/Dec against an independent direct call to `find_src.run_extraction` on the same event file; a reproducibility test confirms two independent fresh reductions agree on RA/Dec and net PHA1 rates. `tool_config.yaml` gained a `call_timeout_s` entry (900s) for the tool's own sandbox backstop.

- **2026-07-02** — Tightened `pa_pd_stability._load_catalog_rows`'s validation: previously a `KeyError` was only raised once a data *row* lacking `obsid` was iterated, so an empty or header-only `catalog_csv` (or one whose header dropped `obsid` but had zero data rows) yielded `{}` silently, and a catalog with `obsid` but missing a field the plotter reads later (`PA_deg_int`/`PA_err_int`/`PD_pct_int`/`PD_err_int`/`rchi2_PD`/`rchi2_PA`) raised mid-plot -- potentially after some of the four figures were already written. `_load_catalog_rows` now validates `DictReader.fieldnames` against the new `_PLOT_REQUIRED_CATALOG_FIELDS` up front and raises `KeyError` naming the missing columns before building any row dict, so every malformed-catalog case (empty, header-only, wrong-schema, or `obsid`-only) is caught by `stability_plot`'s existing `except (OSError, csv.Error, KeyError)` guard *before* any figure is drawn. `CATALOG_FIELDS` is a superset of `_PLOT_REQUIRED_CATALOG_FIELDS`, so a correctly-schema'd catalog with zero rows (as `test_obsid_absent_from_catalog_still_plots` writes) is still accepted. Two more regression tests added (`test_empty_catalog_csv_fails_loudly`, `test_catalog_with_obsid_but_missing_plot_fields_fails_before_writing`, the latter asserting explicitly that no partial figures land on disk).

- **2026-07-02** — Fixed `pa_pd_stability.stability_plot`: a missing or unreadable `catalog_csv` (a required upstream `pa_pd_stability` product) previously fell back silently to an empty catalog (`cat_rows = {}`) and returned `exit_code=0` with four plausible-looking figures — no weighted-mean overlays and an empty chi2-histogram panel — masking a typo'd path or a failed upstream run. `catalog_csv` existence and readability (parseable with an `obsid` column) are now checked before any figure is drawn; either failure now returns `exit_code=1` with `stderr` set and writes no files, the same contract as the all-ObsIDs-unreadable case. This is distinct from — and does not change — the already-correct handling of an ObsID simply *absent from* an otherwise-present, readable catalog, which still gets its time-series/Stokes-plane panels with only the weighted-mean overlay skipped. Two regression tests added (`test_missing_catalog_csv_fails_loudly`, `test_unreadable_catalog_csv_fails_loudly`).

- **2026-07-02** — Added `pa_pd_stability.stability_plot` and the `ixpe_pa_pd_stability_plot` AI tool: the four standard reviewer diagnostic figures for the six-ObsID BH-XRB PA/PD stability workflow (Cyg X-1 01002901/03002599, Swift J1727.8-1613 03005801/03006001, 4U 1630-47 01250401, LMC X-1 02001901), which `ixpe_pa_pd_stability` computes statistics for but does not plot. Pure matplotlib (Agg backend) — no XSPEC/ixpeobssim — so it lives alongside `pa_pd_stability` in `lib/` and is wired directly into `TOOL_MANIFEST`, never through `skills/session.py` (which imports XSPEC). Writes exactly four files per call (`_pa_pd_timeseries`, `_stokes_plane`, `_chi2_hist`, `_structure_function`): per-ObsID PA/PD time series (kiloseconds elapsed from each ObsID's own earliest segment, catalog weighted-mean ±1σ band overlaid), per-ObsID Q/I vs U/I Stokes-plane scatter color-coded by elapsed time, reduced-chi2 histograms across every row of the supplied catalog, and structure functions vs lag (one line per ObsID with enough segments) computed via a `mid_time_structure_functions` helper that itself calls the already-trusted `structure_function` primitive — never an independently-invented calculation. Handles a missing `pcube_timeseg` CSV per ObsID (skipped, noted in `stdout`; all-missing raises `exit_code=1` with no files written) and an ObsID absent from the catalog (still gets its time-series/Stokes-plane panels, only the weighted-mean overlay is skipped) without raising. Deterministic byte output (fixed `dpi`, `bbox_inches`, and `savefig` metadata with no embedded timestamp) verified by a reproducibility test; a differential test confirms the plotted structure-function values against an independent recompute via `structure_function` from the same CSV. Validated on real, already-local data (01002901, 03006001) with the real `pa_pd_stability_catalog.csv` already on disk.

- **2026-07-02** — Added `lib/pa_pd_stability.py` (`pa_pd_stability`) and the `ixpe_pa_pd_stability` AI tool: intra-observation PA/PD/Stokes constancy analysis across one or more local IXPE ObsIDs, computed purely from each ObsID's existing `pcube_timeseg_<obsid>.csv` (all six BH-XRB ObsIDs — 01002901, 01250401, 02001901, 03002599, 03005801, 03006001 — already carry one). Weighted chi-squared constancy tests on PD/PA (mod-180 circular handling)/QI/UI; the authoritative 2D Stokes test (`chi2_2D = chi2_QI + chi2_UI`, `dof_2D = 2(N-1)`, `variable = p_2D < variability_p_threshold`); moment-estimator intrinsic scatter (upper-limit flagged); peak-to-peak ranges; first-order structure functions at configurable lags (`sf_lags_ks`, default `[5,10,20]`, gated on `min_seg_for_sf`); and a cross-source Spearman correlation of mean MDP99 against reduced chi-squared. Pure numpy/scipy — imports neither `xspec` nor `ixpeobssim` (its own small, independently-implemented CSV reader avoids importing `red_ixpe`, which would otherwise pull in the whole science stack) — so `TOOL_MANIFEST` wires it directly to `xspec_utils.lib.pa_pd_stability` (no `skills/session.py` involvement, same pattern as `group_pha`). Validated on real, already-local data: the 2D test correctly separates a known-constant real ObsID (01002901, p_2D ≈ 0.96) from a known-variable one at the default threshold (03006001, p_2D ≈ 0.0079); an independent differential recompute of `chi2_PD` and a differential CSV-parse check against `red_ixpe._read_existing_timeseg_csv` both agree; two independent calls reproduce every catalog number near-exactly (pure arithmetic, no fitting). No `obsdb_catalog` write (a derived summary, not observation data). A companion `plot_pa_pd_stability` plotting tool (consuming this tool's catalog/summary) is the natural next addition.

- **2026-07-01** — Fixed `red_ixpe.pcube_timeseg`/`ixpe_pcube_timeseg` for long observations (03002599, 01250401, 02001901). Three coupled fixes: (1) a GTI pre-reduction gate — each window's GTI-clipped exposure is now checked (new `min_seg_exposure_s`, default 1000.0s) *before* any `xpselect`/`xpbin` call, so windows falling inside a long inter-pointing gap (e.g. 03002599's two glued pointings) are skipped (`n_skipped_empty`) instead of wasting the call budget on spurious near-zero-exposure rows; the 1000s default was set from the real GTI structure of the reported pathological windows (1.7s at 01250401 segment 7, 241.0s at 02001901 segment 148) and validated to sit below every other real segment measured in both observations; (2) incremental write + resume/override — rows are now appended and flushed to the CSV as each segment completes (a provenance comment header records every parameter plus the derived `tstart`/`tstop`), so a call hitting the 5400s cap leaves a real checkpoint; a rerun with `resume=True` (default) and matching provenance computes only the missing segments, while any mismatch deletes and rebuilds the CSV from scratch rather than gluing rows from different configurations; (3) the `ixpe_pcube_timeseg` wrapper no longer returns the full per-segment `timeseries` inline (it could exceed the sandbox transport chunk limit at the ~150-170 segment production default) — it now returns a bounded summary (`csv_path`, `n_segments`, `n_skipped_empty`, `n_expected`, `complete`, `resumed`, `integrated`, `head`/`tail`) and records `complete` in the catalog params. The 5400s per-call cap is held fixed; the gate and checkpointing are what let a long observation finish, across repeated calls if needed. Held-timeout tests (`tests/test_ixpe_pcube_timeseg.py`) run ~30 real `xpselect`/`xpbin` reductions on real, already-local 01250401/02001901 event data to validate the gate, resume, override, and bounded-return behaviour against genuine data — `tool_config.yaml`'s `gate_test_timeout_s` was raised to 1500s accordingly.

- **2026-06-30** — Added `red_ixpe.pcube_timeseg` (time-resolved 2-8 keV IXPE PCUBE polarimetry, DU1/2/3 merged, natural elapsed-time segmentation, background subtraction off by default) and the `ixpe_pcube_timeseg` AI tool (`skills/session.py`, records the written CSV in the mission catalog at level 3 via `obsdb_catalog.neutral.advance`). Validated against a published 3-DU-merged accuracy anchor on IXPE ObsID 01250401 (PD, PA within tolerance); this run also validated `find_src.py` centroiding on IXPE event files for the first time (previously untested — see the `src_pos` notes above).