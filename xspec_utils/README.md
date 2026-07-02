# xspec_utils

A Python/Bash toolkit for joint spectro-polarimetric analysis of X-ray binary data from **IXPE**, **NICER**, and **NuSTAR**, built on top of [HEASoft/XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/) and [ixpeobssim](https://ixpeobssim.readthedocs.io/). Developed for the Cyg X-1 campaign but designed to be reusable for other bright X-ray binaries. The current instrument set (IXPE + NICER + NuSTAR) is a starting point — the loader, reduction scripts, and epoch definitions are all designed to be extended modularly to additional missions (e.g. Chandra, XMM, Swift) by adding a new file locator module and a matching branch in `xspec_help.load_data`.

## Features

- **Batch data reduction** shell scripts for IXPE, NICER, and NuSTAR pipelines
- **Automated point-source detection** via image centroiding (`find_src.py`) with DS9 region output
- **Unified data loader** that assigns XSPEC data groups across instruments in a single call (`xspec_help.load_data`), extensible to new instruments
- **Stokes I/Q/U spectral extraction** for IXPE polarization cubes and PHA spectra
- **Time-resolved IXPE polarimetry** (`red_ixpe.pcube_timeseg`) for intra-observation PD/PA/Stokes/MDP99 stability, DU1/2/3 merged, exposed as the `ixpe_pcube_timeseg` AI tool
- **Joint fitting helpers**: multi-panel plotting with model component decomposition, fit-state freezing for snapshotting
- **Fully customisable epoch management** for multi-observation campaigns, applicable to any source

## Repository layout

```
lib/                        # Python library modules
  xspec_help.py             # Core: data loading, fitting, inspection, rebinning
  red_ixpe.py               # IXPE reduction: PCUBE, PHA1/Q/U extraction, MJD, hardness ratio
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
    ixpe_irf_date: str = '20230702',
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
| `epoch` | `dict` | — | Single-epoch sub-dict from `epoch.epoch`, e.g. `epoch.epoch['E5']`. Expected keys: `'IXPE'`, `'NICER'`, `'NuSTAR'` (each a list of obs-ID strings) and `'STATE'`. Instruments absent from the dict are silently skipped.|
| `ixpe_path` | `str` | — | Path to the directory containing IXPE observation sub-directories. |
| `nicer_path` | `str \| None` | `None` | Path to the NICER data directory. Required when `'NICER'` is in `epoch`. |
| `nustar_path` | `str \| None` | `None` | Path to the NuSTAR data directory. Required when `'NuSTAR'` is in `epoch`. |
| `ixpe_changes` | `bool` | `False` | Re-run the IXPE `ixpeobssim` pipeline if `True`; reuse products already on disk if `False`. |
| `ixpe_src_rad` | `float` | `1.0` | IXPE source-region radius in arcminutes. |
| `ixpe_bkg_inner_rad` | `float` | `3.0` | IXPE background annulus inner radius in arcminutes. |
| `ixpe_bkg_outer_rad` | `float` | `5.0` | IXPE background annulus outer radius in arcminutes. |
| `ixpe_irf_date` | `str` | `'20230702'` | IRF date string (`'YYYYMMDD'`) used to select the IXPE response matrix. Available IRFs live under `/grad/fzhou/.local/lib/python3.13/site-packages/ixpeobssim/caldb/ixpe/gpd/cpf/rmf/`. |
| `ixpe_erange` | `tuple[float, float]` | `(2.0, 8.0)` | Energy range in keV to *notice* for IXPE spectra. |
| `nicer_erange` | `tuple[float, float]` | `(1.0, 10.0)` | Energy range in keV to *notice* for NICER spectra. |
| `nustar_erange` | `tuple[float, float]` | `(4.0, 50.0)` | Energy range in keV to *notice* for NuSTAR spectra. |
| `nustar_use_bkg` | `bool` | `False` | Load NuSTAR background spectrum. Disabled by default because clean background regions are often unavailable for bright sources like Cyg X-1. |
| `clear_existing` | `bool` | `True` | Call `AllData.clear()` and `AllModels.clear()` before loading. |
| `chatter` | `bool` | `False` | Print file paths and per-instrument loading summaries. |

**Data-group assignment**

Data groups are assigned **first-come-first-served with respect to the epoch definition**: the loader iterates over the keys of the epoch dict in insertion order, and the first instrument encountered starts at group 1. The `epoch.py` convention lists `NuSTAR → NICER → IXPE`. Reorder the keys in `epoch.py` and the data-group numbering changes accordingly. Each NuSTAR FPM (FPMA and FPMB) occupies its own data group, so one NuSTAR observation contributes two groups and two observations contribute four groups; each NICER observation occupies one group; each IXPE detector unit gets one group.

**IXPE Stokes gain linking**

For each IXPE DU the three spectra (Stokes I, Q, U) are loaded into consecutive XSPEC indices within the same data group. The gain slope of Q and U are linked to the I gain slope of the same DU (`r.gain.slope.link = AllData(I_idx).response.gain.slope`) to enforce a shared energy-scale calibration. The gain offset is frozen to zero for all three.

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

Per-detector lists are always ordered by spectrum index; IXPE DU lists are ordered `[I, Q, U]`. NuSTAR FPM labels follow sorted filename order.

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

#### `red_ixpe.pcube(obs, file_list, ...)`

Extracts an IXPE polarization cube (PCUBE) for one observation, background-subtracts it, and returns both.

```python
def pcube(
    obs: str,
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
) -> tuple[xBinnedPolarizationCube, xBinnedPolarizationCube]
```

When `changes=False` only `obs` and `file_list` are required; existing `_src_pcube` / `_bkg_pcube` FITS files on disk are read directly. When `changes=True` all spatial and IRF parameters are required and the full `ixpeobssim` pipeline (`xpselect` → `xpbin PCUBE`) is re-run. Background scaling is applied as `bkg_pcube *= src_pcube.backscal() / bkg_pcube.backscal()` before subtraction.

**Returns** `(src_pcube, bkg_pcube)` as `xBinnedPolarizationCube` objects.

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

Same `changes` toggle as `pcube`. When `changes=True`, runs `xpselect` + `xpbin PHA1` with IRF name `ixpe:obssim{irf_date}_alpha075:v13`. One source and one background PHA1 file is produced per detector unit (DU1, DU2, DU3). The Q and U files are derived from the I paths by `red_ixpe.get_pha_q` / `red_ixpe.get_pha_u`.

**Returns** `(src_pha_list, bkg_pha_list)` — each a list of three paths ordered DU1, DU2, DU3.

---

#### `red_ixpe.get_pha_q(pha_path)` / `red_ixpe.get_pha_u(pha_path)`

Derive the Stokes Q or U PHA path from the Stokes I PHA path by replacing the first occurrence of `pha1` with `pha1q` or `pha1u` in the filename.

---

#### `red_ixpe.get_mjd(obs, file_list)`

Reads `TSTART`, `TSTOP`, `MJDREFI`, `MJDREFF`, and `TIMEZERO` from the EVENTS header of the first event file for `obs`.

**Returns** `(mjd_start, mjd_stop, mjd_mid)` as floats.

---

#### `red_ixpe.HR(obs, src_pha_file_list, bkg_pha_file_list, s_bin, h_bin, chatter=False)`

Computes the hardness ratio `H/S` for a given observation by loading PHA files into XSPEC, noticing the soft (`s_bin`) and hard (`h_bin`) energy bands, and summing net count rates.

**Returns** `(hr, hr_err)`. Raises `ValueError` if either band has a non-positive net rate.

---

#### `red_ixpe.pcube_timeseg(obs, file_list, irf_date, ...)`

Time-resolved PCUBE polarimetry for one observation: splits it into natural elapsed-time segments and runs one single-band, DU-merged PCUBE reduction per segment plus one integrated whole-observation reduction, to test intra-observation polarization constancy (e.g. Cyg X-1, 4U 1630-47). This is the tool wrapped as the `ixpe_pcube_timeseg` AI-facing entry (`skills/session.py`).

```python
def pcube_timeseg(
    obs: str,
    file_list: list[str],
    irf_date: str,
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

In principle this catalogue is superseded by the on-the-fly centroiding in `find_src.py`; `find_src` is now validated on IXPE event files by `red_ixpe.pcube_timeseg` (its centroid on 01250401 matches the source's published position and polarization within tolerance — see Changelog), but `src_pos` remains the operational source for the fixed-epoch `xspec_help.load_data` reductions.

## IXPE background filtering

Instrumental background rejection follows the method described in [Di Marco et al. (IXPE background)](https://github.com/aledimarco/IXPE-background). The script `filter_inst_bkg_IXPE.sh` applies `filter_background.py` from that repository to all observations.

## Notes

- The `epoch.py` file defines 12 observing epochs (E1--E12) for Cyg X-1, spanning both the low/hard state (LHS) and the high/soft state (HSS). Adapt this file for other sources.
- `src_pos.py` stores manually determined IXPE source coordinates, used by the fixed-epoch `xspec_help.load_data` path. `find_src.py` centroiding is now also validated on IXPE event files (used automatically by `red_ixpe.pcube_timeseg` when `ra`/`dec` are omitted). For NuSTAR, centroids are computed on-the-fly by `find_src.py`.
- IXPE Stokes Q and U PHA files are derived from the I file path by replacing `pha1` with `pha1q`/`pha1u`.
- NuSTAR background loading is disabled by default (`nustar_use_bkg=False`) because clean background regions are often unavailable for bright sources like Cyg X-1.

## Changelog

- **2026-07-01** — Fixed `red_ixpe.pcube_timeseg`/`ixpe_pcube_timeseg` for long observations (03002599, 01250401, 02001901). Three coupled fixes: (1) a GTI pre-reduction gate — each window's GTI-clipped exposure is now checked (new `min_seg_exposure_s`, default 1000.0s) *before* any `xpselect`/`xpbin` call, so windows falling inside a long inter-pointing gap (e.g. 03002599's two glued pointings) are skipped (`n_skipped_empty`) instead of wasting the call budget on spurious near-zero-exposure rows; the 1000s default was set from the real GTI structure of the reported pathological windows (1.7s at 01250401 segment 7, 241.0s at 02001901 segment 148) and validated to sit below every other real segment measured in both observations; (2) incremental write + resume/override — rows are now appended and flushed to the CSV as each segment completes (a provenance comment header records every parameter plus the derived `tstart`/`tstop`), so a call hitting the 5400s cap leaves a real checkpoint; a rerun with `resume=True` (default) and matching provenance computes only the missing segments, while any mismatch deletes and rebuilds the CSV from scratch rather than gluing rows from different configurations; (3) the `ixpe_pcube_timeseg` wrapper no longer returns the full per-segment `timeseries` inline (it could exceed the sandbox transport chunk limit at the ~150-170 segment production default) — it now returns a bounded summary (`csv_path`, `n_segments`, `n_skipped_empty`, `n_expected`, `complete`, `resumed`, `integrated`, `head`/`tail`) and records `complete` in the catalog params. The 5400s per-call cap is held fixed; the gate and checkpointing are what let a long observation finish, across repeated calls if needed. Held-timeout tests (`tests/test_ixpe_pcube_timeseg.py`) run ~30 real `xpselect`/`xpbin` reductions on real, already-local 01250401/02001901 event data to validate the gate, resume, override, and bounded-return behaviour against genuine data — `tool_config.yaml`'s `gate_test_timeout_s` was raised to 1500s accordingly.

- **2026-06-30** — Added `red_ixpe.pcube_timeseg` (time-resolved 2-8 keV IXPE PCUBE polarimetry, DU1/2/3 merged, natural elapsed-time segmentation, background subtraction off by default) and the `ixpe_pcube_timeseg` AI tool (`skills/session.py`, records the written CSV in the mission catalog at level 3 via `obsdb_catalog.neutral.advance`). Validated against a published 3-DU-merged accuracy anchor on IXPE ObsID 01250401 (PD, PA within tolerance); this run also validated `find_src.py` centroiding on IXPE event files for the first time (previously untested — see the `src_pos` notes above).