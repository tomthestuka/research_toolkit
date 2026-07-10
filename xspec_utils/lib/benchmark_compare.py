"""This-work-vs-paper benchmark comparison -- pure stdlib (no XSPEC/ixpeobssim).

Reproduction studies (e.g. the Cyg X-1 IXPE polarimetry reproduction of
Krawczynski+2022) need a final, reproducible data product that compares
this-work measured quantities against published benchmark values and flags
statistically significant discrepancies. `compare_to_benchmark` is that
capability: it takes a caller-built list of already-measured/cited records
(this-work value+error, paper value+error, a shared unit, and the paper's
citation) and computes delta/z-score/flag arithmetic, writing a
provenance-headed CSV. It hardcodes no paper numbers and re-derives no
science -- the numbers come from the caller (typically via `pcube_quantity`/
`fit_quantity` for the this-work side, and a literature value for the paper
side).

`pcube_quantity`/`fit_quantity` are the two extractor helpers that keep the
error-prone fraction<->percent unit choice out of the caller: they read a
`PD`/`PA`/`QN`/`UN` fraction straight from an `ixpe_pcube` CSV, or a
`PD`/`PA`/`PD_A`/`PD_slope`/`PA_A`/`PA_slope` value straight from a
`fit_tbabs_polconst_powerlaw`/`fit_tbabs_pollin_powerlaw` result dict, and
apply the identical fraction->percent contract on request. `pcube_quantity`
delegates to `red_ixpe.read_pcube_csv` -- the single source of truth for the
`ixpe_pcube` CSV's provenance-header + band-field contract -- via a lazy,
function-local import, so that contract is never duplicated and can't drift.
`red_ixpe.py` imports `ixpeobssim` at module top level, so this module keeps
that import inside `pcube_quantity` itself: `compare_to_benchmark` (pure
arithmetic) and the module's own top level stay import-clean and callable
without the science stack, wired directly into `TOOL_MANIFEST`; only a call
to `pcube_quantity` pulls in `ixpeobssim`, and any caller with a real pcube
CSV in hand already has it installed.

`plot_benchmark_compare` (bottom of this module) is the pure-matplotlib
companion plotter -- the reviewer figure for a `compare_to_benchmark` CSV or
its returned `rows` list: a forest/caterpillar plot of signed z per record,
so a systematic offset (every point pulled the same direction) is visible at
a glance, which a CSV of `abs(z)` values alone cannot show. Only `matplotlib`/
`numpy` beyond stdlib (already used by `pa_pd_stability.py`/
`bhxrb_pd_hr_survey.py`/`ixpe_pcube_plot.py`), so it stays in `lib/` and is
wired directly into `TOOL_MANIFEST`, never through `skills/session.py`.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import time as _time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_PCUBE_QUANTITY_NAMES = ('PD', 'PA', 'QN', 'UN')

_FIT_ANGLE_NAMES = ('PA', 'PA_A', 'PA_slope')
_FIT_FRACTION_NAMES = ('PD', 'PD_A', 'PD_slope')

_HEADER_PREFIX = '# compare_to_benchmark provenance: '
_PRODUCED_BY = 'xspec_utils.lib.benchmark_compare.compare_to_benchmark'
_ROW_FIELDS = ['label', 'unit', 'value_this', 'err_this', 'value_paper',
               'err_paper', 'paper_source', 'delta', 'z', 'flagged']


def _apply_unit(value: float, err: float, unit: str) -> tuple[float, float]:
    """Fraction->percent contract shared by `pcube_quantity`/`fit_quantity`:
    multiply by 100 when `unit=='percent'`, pass the raw fraction through
    when `unit=='fraction'`. Callers of this helper have already screened
    out angle-like quantities (PA and friends), which never convert."""
    if unit == 'percent':
        return value * 100.0, err * 100.0
    return value, err


def pcube_quantity(csv_path: str, band_index: int, name: str, unit: str = 'percent') -> tuple[float, float]:
    """Read one energy band's `(value, err)` for `name` out of an `ixpe_pcube`
    CSV, applying the fraction->percent contract so a caller never has to
    remember it.

    Delegates the actual CSV read to `red_ixpe.read_pcube_csv` -- the single
    source of truth for the `ixpe_pcube` provenance-header + band-field
    contract -- via a lazy, function-local import, so this contract is never
    duplicated. The import is local to this function (not at module top
    level) so that importing `benchmark_compare` itself, and calling
    `compare_to_benchmark`/`fit_quantity`, never requires `ixpeobssim`; only
    a call to `pcube_quantity` does.

    Parameters
    ----------
    csv_path : str
        Path to the `ixpe_pcube_<obsid(s)>.csv` written by `ixpe_pcube`.
    band_index : int
        0-based index into the CSV's own band rows (in file order).
    name : str
        One of 'PD', 'PA', 'QN', 'UN'. For 'PD'/'QN'/'UN' (stored as
        dimensionless fractions), the returned value/err is multiplied by
        100 when `unit=='percent'` and left as the raw fraction when
        `unit=='fraction'`. For 'PA' the stored degrees are returned
        unchanged regardless of `unit`.
    unit : str
        'percent' (default) or 'fraction'.

    Returns
    -------
    (value, err) : tuple[float, float]

    Raises
    ------
    ValueError
        `name`/`unit` unrecognised; `csv_path` missing the `ixpe_pcube`
        provenance header (raised by `red_ixpe.read_pcube_csv` itself and
        propagated unchanged); the CSV carries zero band rows; or
        `band_index` is out of range.
    """
    if name not in _PCUBE_QUANTITY_NAMES:
        raise ValueError(f"pcube_quantity: unsupported name {name!r}; expected one of {_PCUBE_QUANTITY_NAMES}")
    if unit not in ('percent', 'fraction'):
        raise ValueError(f"pcube_quantity: unsupported unit {unit!r}; expected 'percent' or 'fraction'")

    from xspec_utils.lib.red_ixpe import read_pcube_csv
    _, bands = read_pcube_csv(csv_path)

    if not bands:
        raise ValueError(f'pcube_quantity: {csv_path} carries no band rows.')
    if not (0 <= band_index < len(bands)):
        raise ValueError(
            f"pcube_quantity: band_index {band_index} out of range for "
            f"{len(bands)} band(s) in {csv_path}")
    band = bands[band_index]
    value, err = band[name], band[f'{name}_err']
    if name == 'PA':
        return value, err
    return _apply_unit(value, err, unit)


def fit_quantity(fit_dict: dict, name: str, unit: str = 'percent') -> tuple[float, float]:
    """Read `(value, err)` for `name` out of an XSPEC fit result dict (e.g.
    `fit_tbabs_polconst_powerlaw`'s or `fit_tbabs_pollin_powerlaw`'s return
    value), applying the identical fraction->percent contract `pcube_quantity`
    applies.

    Parameters
    ----------
    fit_dict : dict
        A fit result dict carrying `name` and `f'{name}_err'`.
    name : str
        One of 'PD', 'PA' (`fit_tbabs_polconst_powerlaw`) or 'PD_A',
        'PD_slope', 'PA_A', 'PA_slope' (`fit_tbabs_pollin_powerlaw`). 'PD'/
        'PD_A'/'PD_slope' are fractions (or fraction/keV for the slope):
        multiplied by 100 when `unit=='percent'`, left raw when
        `unit=='fraction'`. 'PA'/'PA_A'/'PA_slope' are angles (degrees or
        degrees/keV) and are returned unchanged regardless of `unit`.
    unit : str
        'percent' (default) or 'fraction'.

    Returns
    -------
    (value, err) : tuple[float, float]

    Raises
    ------
    ValueError
        `name`/`unit` unrecognised, or `fit_dict` is missing `name` or
        `f'{name}_err'`.
    """
    if name not in _FIT_ANGLE_NAMES + _FIT_FRACTION_NAMES:
        raise ValueError(
            f"fit_quantity: unsupported name {name!r}; expected one of "
            f"{_FIT_ANGLE_NAMES + _FIT_FRACTION_NAMES}")
    if unit not in ('percent', 'fraction'):
        raise ValueError(f"fit_quantity: unsupported unit {unit!r}; expected 'percent' or 'fraction'")
    err_key = f'{name}_err'
    if name not in fit_dict or err_key not in fit_dict:
        raise ValueError(f"fit_quantity: fit_dict is missing {name!r}/{err_key!r}")
    value, err = float(fit_dict[name]), float(fit_dict[err_key])
    if name in _FIT_ANGLE_NAMES:
        return value, err
    return _apply_unit(value, err, unit)


def _resolve_unit(rec: dict, label: str) -> str:
    """Resolve one record's effective shared unit.

    A record supplies either a single shared `unit`, or the pair
    `unit_this`/`unit_paper` (both required together), or both -- in which
    case every supplied value must agree. Whitespace is stripped before
    comparison; an empty unit on any side, or a genuine mismatch between the
    paired sides, is a `ValueError` naming `label` and the conflicting units
    (this tool does no unit conversion -- it only verifies the caller already
    put both sides in the same unit).
    """
    has_this = 'unit_this' in rec
    has_paper = 'unit_paper' in rec
    if has_this != has_paper:
        raise ValueError(
            f"compare_to_benchmark: record {label!r} supplies only one of "
            f"unit_this/unit_paper -- provide both or neither")

    pair_unit = None
    if has_this and has_paper:
        unit_this = str(rec['unit_this']).strip()
        unit_paper = str(rec['unit_paper']).strip()
        if not unit_this or not unit_paper:
            raise ValueError(f"compare_to_benchmark: record {label!r} has an empty unit_this/unit_paper")
        if unit_this != unit_paper:
            raise ValueError(
                f"compare_to_benchmark: record {label!r} has mismatched unit between "
                f"the paired sides: unit_this={unit_this!r} vs unit_paper={unit_paper!r}")
        pair_unit = unit_this

    base_unit = None
    if 'unit' in rec:
        base_unit = str(rec['unit']).strip()
        if not base_unit:
            raise ValueError(f"compare_to_benchmark: record {label!r} has an empty unit")

    if pair_unit is not None and base_unit is not None and pair_unit != base_unit:
        raise ValueError(
            f"compare_to_benchmark: record {label!r} has mismatched unit between "
            f"the paired sides: unit={base_unit!r} vs unit_this/unit_paper={pair_unit!r}")

    resolved = pair_unit if pair_unit is not None else base_unit
    if resolved is None:
        raise ValueError(
            f"compare_to_benchmark: record {label!r} is missing 'unit' (or the "
            f"unit_this/unit_paper pair)")
    return resolved


def _validate_record(rec: dict) -> dict:
    """Validate one input record and return its normalized fields (see
    `compare_to_benchmark`'s docstring for the required keys)."""
    label = rec.get('label')
    if not isinstance(label, str) or not label.strip():
        raise ValueError(f"compare_to_benchmark: record missing a non-empty 'label': {rec!r}")
    label = label.strip()

    paper_source = rec.get('paper_source')
    if not isinstance(paper_source, str) or not paper_source.strip():
        raise ValueError(
            f"compare_to_benchmark: record {label!r} has an empty/missing paper_source "
            f"-- refusing to benchmark against an unattributed value")

    for key in ('value_this', 'value_paper'):
        val = rec.get(key)
        if isinstance(val, bool) or not isinstance(val, (int, float)):
            raise ValueError(f"compare_to_benchmark: record {label!r} needs a numeric {key!r}")
        if not math.isfinite(val):
            raise ValueError(f"compare_to_benchmark: record {label!r} needs a finite {key!r}, got {val!r}")

    for key in ('err_this', 'err_paper'):
        val = rec.get(key)
        if (isinstance(val, bool) or not isinstance(val, (int, float))
                or not math.isfinite(val) or val < 0):
            raise ValueError(
                f"compare_to_benchmark: record {label!r} needs a finite non-negative "
                f"numeric {key!r}, got {val!r}")

    unit = _resolve_unit(rec, label)

    return {
        'label': label,
        'unit': unit,
        'value_this': float(rec['value_this']),
        'err_this': float(rec['err_this']),
        'value_paper': float(rec['value_paper']),
        'err_paper': float(rec['err_paper']),
        'paper_source': paper_source.strip(),
    }


def _write_csv(out_csv: str, rows: list[dict], sigma_flag: float) -> None:
    out_dir = os.path.dirname(os.path.abspath(out_csv))
    os.makedirs(out_dir, exist_ok=True)
    provenance = {
        'sigma_flag': sigma_flag,
        'created_at': _time.strftime('%Y-%m-%dT%H:%M:%SZ', _time.gmtime()),
        'produced_by': _PRODUCED_BY,
    }
    with open(out_csv, 'w', newline='') as f:
        f.write(_HEADER_PREFIX + json.dumps(provenance, sort_keys=True) + '\n')
        writer = csv.DictWriter(f, fieldnames=_ROW_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def compare_to_benchmark(records: list[dict], out_csv: str, sigma_flag: float = 3.0) -> dict:
    """Compare this-work measured quantities against cited paper values and
    flag statistically significant discrepancies.

    Each element of `records` must carry: `label` (non-empty str),
    `value_this` (finite float), `err_this` (finite float >= 0),
    `value_paper` (finite float), `err_paper` (finite float >= 0) --
    `NaN`/`inf` are refused on all four, since a failed upstream measurement
    must not silently propagate into a `z=nan`/unreliable-`flagged` row --
    `paper_source` (non-empty str, the literature citation this record is
    benchmarked against), and either a shared `unit` (non-empty str) or the
    `unit_this`/`unit_paper` pair (both required together, and must agree
    with each other and with `unit` if all are given) -- see
    `_resolve_unit`. This tool does NO unit conversion: both sides must
    already be expressed in the same unit (e.g. via `pcube_quantity`/
    `fit_quantity`); a genuine physical-unit mismatch is a caller error this
    tool refuses rather than silently reconciles.

    For each record: `delta = value_this - value_paper` (signed: positive
    when this-work exceeds the paper value),
    `z = abs(delta) / sqrt(err_this**2 + err_paper**2)` (`z = inf` when both
    errors are legitimately zero -- a computed denominator of zero, never a
    `ZeroDivisionError` -- distinct from the input-level `NaN`/`inf` refusal
    above, which is refused before this arithmetic ever runs), and
    `flagged = abs(z) > sigma_flag`. Writes a provenance-headed CSV
    (`sigma_flag`, `created_at` UTC, tool provenance in a
    `# compare_to_benchmark provenance: {...}` JSON comment line) to
    `out_csv`, creating parent directories as needed.

    Parameters
    ----------
    records : list[dict]
        Non-empty list of records as described above.
    out_csv : str
        Output CSV path.
    sigma_flag : float
        z-score threshold above which a record is flagged. Default 3.0.

    Returns
    -------
    dict with 'rows' (one dict per record: label, unit, value_this, err_this,
    value_paper, err_paper, paper_source, delta, z, flagged), 'flags' (list
    of flagged labels), 'out_csv' (absolute path written), 'n_records',
    'stdout', 'exit_code' (always 0 -- invalid input raises instead).

    Raises
    ------
    ValueError
        `records` is empty, or any record fails validation (see above).
    """
    if not records:
        raise ValueError('empty records')

    rows = []
    flags = []
    for rec in records:
        v = _validate_record(rec)
        delta = v['value_this'] - v['value_paper']
        denom = math.hypot(v['err_this'], v['err_paper'])
        z = float('inf') if denom == 0 else abs(delta) / denom
        flagged = bool(abs(z) > sigma_flag)
        row = {**v, 'delta': delta, 'z': z, 'flagged': flagged}
        rows.append(row)
        if flagged:
            flags.append(v['label'])

    out_csv = os.path.abspath(out_csv)
    _write_csv(out_csv, rows, sigma_flag)

    return {
        'rows': rows,
        'flags': flags,
        'out_csv': out_csv,
        'n_records': len(rows),
        'stdout': (f'compare_to_benchmark: {len(rows)} record(s), {len(flags)} '
                   f'flagged at >{sigma_flag} sigma -> {out_csv}'),
        'exit_code': 0,
    }


# --------------------------------------------------------------------------- #
# plot_benchmark_compare -- pure-matplotlib companion plotter
# --------------------------------------------------------------------------- #

_PLOT_PRODUCED_BY = 'xspec_utils.lib.benchmark_compare.plot_benchmark_compare'
_PLOT_ROW_KEYS = ('label', 'err_this', 'err_paper')


def _read_compare_csv(csv_path: str) -> tuple[dict, list[dict]]:
    """Read a `compare_to_benchmark` CSV into `(provenance, rows)`, each row
    reduced to the fields `plot_benchmark_compare` needs (`label`, `err_this`,
    `err_paper`, `delta`). Raises `ValueError` -- naming the offending path --
    when `csv_path` does not exist or its first line is not the
    `compare_to_benchmark` provenance header."""
    if not os.path.exists(csv_path):
        raise ValueError(f'plot_benchmark_compare: csv_path not found: {csv_path}')
    with open(csv_path, newline='') as f:
        header_line = f.readline()
        if not header_line.startswith(_HEADER_PREFIX):
            raise ValueError(
                f'plot_benchmark_compare: {csv_path} is missing the compare_to_benchmark '
                f'provenance header; it is not a product of compare_to_benchmark.')
        provenance = json.loads(header_line[len(_HEADER_PREFIX):])
        raw_rows = list(csv.DictReader(f))
    rows = [{
        'label': r['label'],
        'err_this': float(r['err_this']),
        'err_paper': float(r['err_paper']),
        'delta': float(r['delta']),
    } for r in raw_rows]
    return provenance, rows


def _row_signed_z(row: dict) -> float:
    """Signed z = delta / sqrt(err_this**2 + err_paper**2) for one
    `compare_to_benchmark` row (direct or CSV-parsed) -- the direction-
    preserving companion to the CSV's own `abs(z)` 'z' column. `delta`
    (`value_this - value_paper`) is read directly if present, else
    recomputed from `value_this`/`value_paper`. +/-inf, signed by `delta`
    (never a `ZeroDivisionError`), when both errors are zero -- `delta==0`
    resolves to +inf, mirroring `compare_to_benchmark`'s own
    always-flagged treatment of a zero-zero-error record."""
    err_this = float(row['err_this'])
    err_paper = float(row['err_paper'])
    if row.get('delta') not in (None, ''):
        delta = float(row['delta'])
    else:
        delta = float(row['value_this']) - float(row['value_paper'])
    denom = math.hypot(err_this, err_paper)
    if denom == 0:
        return math.inf if delta >= 0 else -math.inf
    return delta / denom


def _build_compare_figure(rows: list[dict], sigma_flag: float, title: str | None = None):
    """Build the forest/caterpillar plot of signed z, one tick per record in
    `rows`' own order. Returns `(fig, ax, plotted)`; the caller owns the
    figure's lifecycle. `plotted` carries `label`/`z`/`flagged` (the exact
    per-record values drawn) plus `sigma_flag`/`n_records`."""
    labels = [str(r['label']) for r in rows]
    zs = [_row_signed_z(r) for r in rows]
    n = len(rows)

    finite_abs = [abs(z) for z in zs if math.isfinite(z)]
    half_range = max(finite_abs) * 1.15 if finite_abs else 0.0
    half_range = max(half_range, sigma_flag * 1.3, 1.0)
    clamp = 0.95 * half_range

    y = np.arange(n)
    flagged = np.array([abs(z) > sigma_flag for z in zs])
    x = np.array([math.copysign(clamp, z) if math.isinf(z) else z for z in zs])

    fig, ax = plt.subplots(figsize=(7.5, max(2.5, 0.4 * n + 1.2)))
    ax.axvspan(-sigma_flag, sigma_flag, color='0.5', alpha=0.15, zorder=0,
               label=f'|z| <= {sigma_flag}')
    ax.axvline(0.0, color='k', lw=0.8, zorder=1)

    if (~flagged).any():
        ax.scatter(x[~flagged], y[~flagged], color='tab:blue', zorder=3,
                   label='within threshold')
    if flagged.any():
        ax.scatter(x[flagged], y[flagged], color='tab:red', zorder=3,
                   label=f'|z| > {sigma_flag}')
        for xi, yi, lbl, zi in zip(x[flagged], y[flagged],
                                    np.array(labels)[flagged], np.array(zs)[flagged]):
            note = lbl + (' [z=inf]' if math.isinf(zi) else '')
            ax.annotate(note, (xi, yi), textcoords='offset points', xytext=(6, 0),
                        va='center', fontsize=8, color='tab:red')

    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlim(-half_range, half_range)
    ax.set_xlabel('signed z = (this-work - paper) / sqrt(err_this^2 + err_paper^2)')
    ax.set_title(title or 'This-work vs. paper benchmark comparison')
    ax.legend(fontsize=8, loc='best')
    fig.tight_layout()

    plotted = {
        'label': labels,
        'z': [float(z) for z in zs],
        'flagged': [bool(f) for f in flagged],
        'sigma_flag': float(sigma_flag),
        'n_records': n,
    }
    return fig, ax, plotted


def plot_benchmark_compare(rows: list[dict] | None = None, csv_path: str | None = None,
                            save_path: str | None = None, title: str | None = None,
                            sigma_flag: float | None = None) -> dict:
    """Render the forest/caterpillar reviewer figure for a `compare_to_benchmark`
    result: one tick per record (in the supplied order), x-axis the SIGNED
    z = (value_this-value_paper)/sqrt(err_this**2+err_paper**2) -- so a
    systematic offset (every point on the same side of zero) is visible at a
    glance, which the CSV's own `abs(z)` 'z' column cannot show. A vertical
    line marks z=0 and a shaded band marks +/-`sigma_flag`; a record with
    `abs(z) > sigma_flag` is drawn in a distinct colour and annotated with its
    label. Writes exactly ONE PNG.

    Takes EITHER `rows` (e.g. `compare_to_benchmark`'s own returned `rows`
    list) OR `csv_path` (the CSV it wrote) -- supply exactly one; `rows`
    takes precedence if both are given. `sigma_flag=None` (default) reads the
    threshold from the CSV's own provenance header when `csv_path` is used,
    else falls back to 3.0; pass an explicit float to override either way.

    A record with both errors zero has an infinite z (signed by `delta`);
    it is clamped to the finite-z axis edge (with an explicit "[z=inf]"
    annotation) rather than leaving `set_xlim` an infinite/NaN range. An
    empty `rows` (or a `csv_path` with zero data rows) is not an error: no
    PNG is written and `exit_code=1` with `stderr` set.

    Parameters
    ----------
    rows : list[dict] or None
        Records carrying `label`, `err_this`, `err_paper`, and either `delta`
        or `value_this`/`value_paper` (`compare_to_benchmark`'s own row shape
        satisfies this).
    csv_path : str or None
        Path to a `compare_to_benchmark` output CSV.
    save_path : str or None
        Full output PNG path. Auto-generated timestamped file in
        `os.getcwd()` when None.
    title : str or None
        Figure title. Defaults to a generic caption.
    sigma_flag : float or None
        z-score threshold for the shaded band / flag colour. See above.

    Returns
    -------
    dict with 'figures' (one `FigureRef` dict, or `[]` on the empty-rows
    edge case), 'plotted' (`label`/`z`/`flagged`/`sigma_flag`/`n_records` --
    the exact per-record values drawn), 'stdout', 'exit_code' (1 with
    'stderr' set only for empty rows; raises `ValueError` for a missing
    csv_path/malformed CSV or a row missing required fields).
    """
    if rows is None and csv_path is None:
        raise ValueError('plot_benchmark_compare: supply rows or csv_path')

    provenance: dict = {}
    if rows is None:
        provenance, rows = _read_compare_csv(csv_path)

    if sigma_flag is None:
        sigma_flag = float(provenance['sigma_flag']) if 'sigma_flag' in provenance else 3.0

    if not rows:
        return {
            'figures': [],
            'stdout': '',
            'stderr': 'plot_benchmark_compare: no records to plot (rows is empty).',
            'exit_code': 1,
        }

    for i, row in enumerate(rows):
        missing = [k for k in _PLOT_ROW_KEYS if k not in row]
        if 'delta' not in row and not {'value_this', 'value_paper'} <= row.keys():
            missing.append("'delta' (or 'value_this'/'value_paper')")
        if missing:
            raise ValueError(f'plot_benchmark_compare: row {i} missing required field(s) {missing}')

    fig, _ax, plotted = _build_compare_figure(rows, sigma_flag, title=title)

    if save_path is None:
        ts = _time.strftime('%Y%m%dT%H%M%S')
        save_path = os.path.join(os.getcwd(), f'plot_benchmark_compare_{ts}.png')
    os.makedirs(os.path.dirname(os.path.abspath(save_path)) or '.', exist_ok=True)

    fig.savefig(save_path, dpi=150, bbox_inches='tight', metadata={'Software': _PLOT_PRODUCED_BY})
    plt.close(fig)
    sha = hashlib.sha256(Path(save_path).read_bytes()).hexdigest()

    n_flagged = sum(plotted['flagged'])
    return {
        'figures': [{
            'path': save_path,
            'caption': (f'This-work vs. paper benchmark comparison: {len(rows)} record(s), '
                        f'{n_flagged} flagged at |z| > {sigma_flag}'),
            'mime_type': 'image/png',
            'produced_by': _PLOT_PRODUCED_BY,
            'created_at': _time.strftime('%Y-%m-%dT%H:%M:%SZ', _time.gmtime()),
            'sha256': sha,
        }],
        'plotted': plotted,
        'stdout': (f'plot_benchmark_compare: {len(rows)} record(s), {n_flagged} flagged '
                   f'at |z| > {sigma_flag} -> {save_path}'),
        'exit_code': 0,
    }
