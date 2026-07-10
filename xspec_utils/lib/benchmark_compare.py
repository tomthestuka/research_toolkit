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
"""

from __future__ import annotations

import csv
import json
import math
import os
import time as _time

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
