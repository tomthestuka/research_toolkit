"""BH-XRB polarization-degree vs hardness-ratio survey statistics -- pure numpy/scipy.

Correlates PD (from the existing `pcube_timeseg_<obsid>.csv` integrated
`segment_id=-1` row: `PD_pct`, `PD_err`, `MDP99`) against HR (from
`bhxrb_hardness_ratios.csv`, written by `red_ixpe.HR`/`ixpe_HR`) across a
survey of IXPE ObsIDs spanning one or more BH-XRB sources, classifies
detections (`PD_pct > MDP99`), and runs the standard PD-HR correlation and
group-comparison tests: per-source and combined Spearman rank correlation,
a two-sample Kolmogorov-Smirnov test between caller-defined source groups,
and a Mann-Whitney U test on HR between detected and non-detected ObsIDs.

Same module family as `red_ixpe.HR` and `pa_pd_stability.pa_pd_stability`:
pure numpy/scipy, no `xspec`/`ixpeobssim` import, so this module -- and the
`pd_hr_survey` entry point -- runs standalone, without the science stack.

`_read_pcube_integrated` is an independent reimplementation of the
provenance-header / `segment_id=-1` integrated-row convention
`pa_pd_stability._read_timeseg_csv` and `red_ixpe._read_existing_timeseg_csv`
already follow -- not a shared helper, so this module never imports either
of them.
"""

from __future__ import annotations

import csv
import hashlib
import os
import time as _time
from pathlib import Path
from typing import Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from obsdb_catalog import neutral

_PCUBE_PROVENANCE_PREFIX = '# pcube_timeseg provenance: '

SURVEY_FIELDS = ['obsid', 'source', 'PD_pct', 'PD_err', 'MDP99',
                 'detected', 'mdp_flagged', 'HR', 'HR_err', 'note']
INTRA_FIELDS = ['source', 'n', 'rho', 'p']
COMBINED_FIELDS = ['n', 'rho', 'p']
KS_FIELDS = ['n_soft', 'n_hard', 'statistic', 'p']
BIAS_FIELDS = ['n_detected', 'n_non_detected', 'statistic', 'p']


def _read_pcube_integrated(csv_path: str) -> dict | None:
    """Read the `segment_id=-1` integrated row's `PD_pct`/`PD_err`/`MDP99`
    from a `pcube_timeseg_<obsid>.csv`.

    Mirrors `pa_pd_stability._read_timeseg_csv`'s provenance-header /
    `segment_id=-1` convention: returns None when `csv_path` does not exist,
    its first line is not the provenance comment header, or no integrated
    row is present -- an unrecognised or incomplete file is treated as
    absent, never parsed as headerless or partial data.
    """
    if not os.path.exists(csv_path):
        return None
    with open(csv_path, newline='') as f:
        header_line = f.readline()
        if not header_line.startswith(_PCUBE_PROVENANCE_PREFIX):
            return None
        for row in csv.DictReader(f):
            if int(row['segment_id']) == -1:
                return {'PD_pct': float(row['PD_pct']), 'PD_err': float(row['PD_err']),
                        'MDP99': float(row['MDP99'])}
    return None


def _read_hr_rows(hr_csv: str) -> dict[str, dict]:
    """Read `bhxrb_hardness_ratios.csv` (written by `red_ixpe.HR`/`ixpe_HR`)
    into `obsid -> {'HR': ..., 'HR_err': ...}`.

    `hr_csv` is a required upstream product -- a missing file is a hard
    input error (raises `ValueError`), distinct from a single obsid's row
    being absent from an otherwise-present file (recorded as a per-obsid
    note, never raised).
    """
    if not os.path.exists(hr_csv):
        raise ValueError(f'pd_hr_survey: hr_csv not found: {hr_csv}')
    with open(hr_csv, newline='') as f:
        return {r['obsid']: {'HR': float(r['HR']), 'HR_err': float(r['HR_err'])}
                for r in csv.DictReader(f)}


def _spearman(xs: Sequence[float], ys: Sequence[float]) -> dict:
    """Spearman rank correlation over the finite (x, y) pairs. NaN with n<3
    (matches `pa_pd_stability._spearman`'s convention)."""
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < 3:
        return {'n': n, 'rho': float('nan'), 'p': float('nan')}
    res = stats.spearmanr(x[mask], y[mask])
    return {'n': n, 'rho': float(res.statistic), 'p': float(res.pvalue)}


def _ks_two_sample(soft: Sequence[float], hard: Sequence[float]) -> dict:
    """Two-sample KS test between the two PD distributions. NaN if either
    group is empty (rather than letting `scipy.stats.ks_2samp` raise)."""
    a = np.asarray(soft, dtype=float)
    a = a[np.isfinite(a)]
    b = np.asarray(hard, dtype=float)
    b = b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return {'n_soft': int(a.size), 'n_hard': int(b.size),
                'statistic': float('nan'), 'p': float('nan')}
    res = stats.ks_2samp(a, b)
    return {'n_soft': int(a.size), 'n_hard': int(b.size),
            'statistic': float(res.statistic), 'p': float(res.pvalue)}


def _mann_whitney(detected: Sequence[float], non_detected: Sequence[float]) -> dict:
    """Two-sided Mann-Whitney U test on HR between the two groups. NaN if
    either group is empty (rather than letting `scipy.stats.mannwhitneyu`
    raise)."""
    a = np.asarray(detected, dtype=float)
    a = a[np.isfinite(a)]
    b = np.asarray(non_detected, dtype=float)
    b = b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return {'n_detected': int(a.size), 'n_non_detected': int(b.size),
                'statistic': float('nan'), 'p': float('nan')}
    res = stats.mannwhitneyu(a, b, alternative='two-sided')
    return {'n_detected': int(a.size), 'n_non_detected': int(b.size),
            'statistic': float(res.statistic), 'p': float(res.pvalue)}


def _write_csv(path: str, fieldnames: list[str], rows: list[dict]) -> None:
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, restval='')
        writer.writeheader()
        writer.writerows(rows)


def pd_hr_survey(
    obs_db_path: str,
    obsids_with_sources: dict[str, str],
    hr_csv: str,
    pcube_csv_dir: str,
    mdp_flag_pct: float = 15.0,
    min_n_intra: int = 4,
    soft_state_sources: Sequence[str] = (),
    hard_state_sources: Sequence[str] = (),
) -> dict:
    """PD-vs-HR survey correlation/statistics across a set of IXPE BH-XRB ObsIDs.

    For each `obsid` in `obsids_with_sources` (obsid -> source-name label),
    reads its PD from `<pcube_csv_dir>/<obsid>/pcube_timeseg_<obsid>.csv`'s
    integrated (`segment_id=-1`) row (`PD_pct`, `PD_err`, `MDP99`) and its HR
    from `hr_csv` (the `bhxrb_hardness_ratios.csv` `red_ixpe.HR`/`ixpe_HR`
    writes), classifies `detected = PD_pct > MDP99`, and flags
    `mdp_flagged = MDP99 > mdp_flag_pct`.

    An obsid with no readable pcube CSV, or with no row in `hr_csv`, is
    recorded with a non-empty `note` and excluded from every statistic below
    -- never raised, never a dummy detected/HR value.

    Computes, over the detected-only rows (a non-detection's PD is
    MDP99-dominated, not a physically meaningful polarization measurement):

    - per-source Spearman rank correlation (PD_pct vs HR) for every source
      with >= `min_n_intra` detections (a source short of that is omitted
      from `intra_source_correlations.csv`, not crashed on);
    - the combined Spearman correlation across every detected point,
      regardless of source;
    - a two-sample Kolmogorov-Smirnov test between the detected PD_pct
      values of `soft_state_sources` vs `hard_state_sources` (both passed by
      the caller -- never hardcoded -- so this serves any source-group
      split, not just one survey's own state assignments).

    And, over every obsid with both a readable PD and HR (detected and
    non-detected together):

    - a two-sided Mann-Whitney U test on HR between the detected and
      non-detected groups (a detection-bias check: does HR predict whether a
      polarization detection was made).

    Every statistic above is NaN (never an exception) when its input group
    is empty or -- for the two Spearman correlations -- has fewer than 3
    finite pairs.

    Writes, under `obs_db_path`: `bhxrb_pd_hardness_survey.csv` (one row per
    requested obsid, `SURVEY_FIELDS`), `intra_source_correlations.csv`
    (`INTRA_FIELDS`), `combined_correlation.csv` (`COMBINED_FIELDS`),
    `ks_test_soft_vs_hard.csv` (`KS_FIELDS`), `detection_bias_assessment.csv`
    (`BIAS_FIELDS`). Records one `kind='pd_hardness_survey'` product (via
    `obsdb_catalog.neutral.advance`, mission `'ixpe'`, level 3) per obsid
    that has both a readable PD and HR, pointing at the shared survey CSV
    with that obsid's own PD/HR/detection parameters.

    Parameters
    ----------
    obs_db_path : str
        Directory the five output CSVs are written to, and the mission
        catalog (`catalog.json`) is recorded in.
    obsids_with_sources : dict[str, str]
        obsid -> source-name label, e.g. `{"01002901": "Cyg X-1", ...}`.
        Every requested obsid appears in `bhxrb_pd_hardness_survey.csv`
        (with a `note` if skipped), and the label groups obsids for the
        intra-source correlations and the soft/hard KS split.
    hr_csv : str
        Path to `bhxrb_hardness_ratios.csv` (written by `red_ixpe.HR`/
        `ixpe_HR`). A missing file is a hard error; a missing per-obsid row
        is a per-obsid note.
    pcube_csv_dir : str
        Directory holding each obsid's `<obsid>/pcube_timeseg_<obsid>.csv`
        (written by `red_ixpe.pcube_timeseg`/`ixpe_pcube_timeseg`).
    mdp_flag_pct : float
        `MDP99` above this (percent) flags the obsid `mdp_flagged` in the
        survey CSV -- a data-quality note, not an exclusion. Default 15.0.
    min_n_intra : int
        Minimum detected-obsid count a source needs before its intra-source
        Spearman correlation is computed at all. Default 4.
    soft_state_sources, hard_state_sources : Sequence[str]
        Source-name labels (matched against `obsids_with_sources`' values)
        forming the two KS-test groups. Default empty -- the KS test is
        NaN until the caller supplies both groups.

    Returns
    -------
    dict with 'rows' (one dict per requested obsid, per `SURVEY_FIELDS`),
    'skipped' (obsids with a non-empty note), 'per_source_correlations'
    (source -> `_spearman` result, only sources with >= `min_n_intra`
    detections), 'spearman_combined', 'ks_result', 'mwu_result' (each a
    dict per the helpers above), 'survey_csv', 'intra_csv', 'combined_csv',
    'ks_csv', 'bias_csv' (output paths), 'exit_code' (always 0 -- a
    per-obsid problem is recorded, not raised).
    """
    hr_rows = _read_hr_rows(hr_csv)
    soft_set, hard_set = set(soft_state_sources), set(hard_state_sources)

    os.makedirs(obs_db_path, exist_ok=True)

    rows = []
    for obsid, source in obsids_with_sources.items():
        pcube_path = os.path.join(pcube_csv_dir, obsid, f'pcube_timeseg_{obsid}.csv')
        pd = _read_pcube_integrated(pcube_path)
        base = {'obsid': obsid, 'source': source, 'PD_pct': float('nan'),
                'PD_err': float('nan'), 'MDP99': float('nan'), 'detected': '',
                'mdp_flagged': '', 'HR': float('nan'), 'HR_err': float('nan'), 'note': ''}
        if pd is None:
            base['note'] = f'no readable pcube_timeseg CSV at {pcube_path}'
            rows.append(base)
            continue
        hr = hr_rows.get(obsid)
        if hr is None:
            base.update(pd)
            base['note'] = f'no HR row for {obsid} in {hr_csv}'
            rows.append(base)
            continue
        detected = bool(pd['PD_pct'] > pd['MDP99'])
        base.update(pd, HR=hr['HR'], HR_err=hr['HR_err'],
                    detected=detected, mdp_flagged=bool(pd['MDP99'] > mdp_flag_pct))
        rows.append(base)

    survey_csv = os.path.join(obs_db_path, 'bhxrb_pd_hardness_survey.csv')
    _write_csv(survey_csv, SURVEY_FIELDS, rows)

    ok = [r for r in rows if not r['note']]
    detected_rows = [r for r in ok if r['detected']]
    non_detected_rows = [r for r in ok if not r['detected']]

    per_source: dict[str, dict] = {}
    intra_rows = []
    by_source: dict[str, list[dict]] = {}
    for r in detected_rows:
        by_source.setdefault(r['source'], []).append(r)
    for source, srows in by_source.items():
        if len(srows) < min_n_intra:
            continue
        res = _spearman([r['PD_pct'] for r in srows], [r['HR'] for r in srows])
        per_source[source] = res
        intra_rows.append({'source': source, **res})
    intra_csv = os.path.join(obs_db_path, 'intra_source_correlations.csv')
    _write_csv(intra_csv, INTRA_FIELDS, intra_rows)

    spearman_combined = _spearman([r['PD_pct'] for r in detected_rows],
                                   [r['HR'] for r in detected_rows])
    combined_csv = os.path.join(obs_db_path, 'combined_correlation.csv')
    _write_csv(combined_csv, COMBINED_FIELDS, [spearman_combined])

    soft_pd = [r['PD_pct'] for r in detected_rows if r['source'] in soft_set]
    hard_pd = [r['PD_pct'] for r in detected_rows if r['source'] in hard_set]
    ks_result = _ks_two_sample(soft_pd, hard_pd)
    ks_csv = os.path.join(obs_db_path, 'ks_test_soft_vs_hard.csv')
    _write_csv(ks_csv, KS_FIELDS, [ks_result])

    mwu_result = _mann_whitney([r['HR'] for r in detected_rows],
                                [r['HR'] for r in non_detected_rows])
    bias_csv = os.path.join(obs_db_path, 'detection_bias_assessment.csv')
    _write_csv(bias_csv, BIAS_FIELDS, [mwu_result])

    for r in ok:
        neutral.advance(obs_db_path, 'ixpe', r['obsid'], level=3, products=[{
            'kind': 'pd_hardness_survey',
            'path': os.path.relpath(survey_csv, obs_db_path),
            'produced_by': 'xspec_utils.lib.bhxrb_pd_hr_survey.pd_hr_survey',
            'params': {'source': r['source'], 'PD_pct': r['PD_pct'], 'MDP99': r['MDP99'],
                       'detected': r['detected'], 'mdp_flagged': r['mdp_flagged'],
                       'HR': r['HR'], 'mdp_flag_pct': mdp_flag_pct},
        }])

    return {
        'rows': rows,
        'skipped': [r['obsid'] for r in rows if r['note']],
        'per_source_correlations': per_source,
        'spearman_combined': spearman_combined,
        'ks_result': ks_result,
        'mwu_result': mwu_result,
        'survey_csv': survey_csv, 'intra_csv': intra_csv, 'combined_csv': combined_csv,
        'ks_csv': ks_csv, 'bias_csv': bias_csv,
        'exit_code': 0,
    }


# --------------------------------------------------------------------------- #
# pd_hr_survey_plot -- pure-matplotlib companion plotter (no xspec/ixpeobssim)
# --------------------------------------------------------------------------- #

_PRODUCED_BY = 'xspec_utils.lib.bhxrb_pd_hr_survey.pd_hr_survey_plot'


def _safe_float(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return float('nan')


def _load_survey_rows(survey_csv: str) -> list[dict]:
    """Read `bhxrb_pd_hardness_survey.csv` (written by `pd_hr_survey`) into a
    list of typed row dicts. Raises `KeyError` naming any missing `SURVEY_FIELDS`
    column up front, before any row is parsed -- a malformed/wrong-schema (or
    header-only) CSV is rejected before any figure is drawn."""
    with open(survey_csv, newline='') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        missing = [c for c in SURVEY_FIELDS if c not in fieldnames]
        if missing:
            raise KeyError(f'survey_csv missing required column(s) {missing}; header was {fieldnames}')
        return [{
            'obsid': r['obsid'], 'source': r['source'], 'note': r['note'],
            'PD_pct': _safe_float(r['PD_pct']), 'PD_err': _safe_float(r['PD_err']),
            'MDP99': _safe_float(r['MDP99']),
            'detected': r['detected'] == 'True', 'mdp_flagged': r['mdp_flagged'] == 'True',
            'HR': _safe_float(r['HR']), 'HR_err': _safe_float(r['HR_err']),
        } for r in reader]


def _read_integrated_row(csv_path: str) -> dict | None:
    """Read the full `segment_id=-1` integrated row (every `pcube_timeseg`
    column, e.g. `exposure_s`) a `pcube_timeseg_<obsid>.csv` carries.

    Same provenance-header / "unrecognised file -> absent" contract as
    `_read_pcube_integrated`, kept as an independent function (rather than
    widening `_read_pcube_integrated`'s return, which
    `test_read_pcube_integrated_round_trip` pins to exactly `PD_pct`/`PD_err`/
    `MDP99`) since the feasibility panel below needs `exposure_s` too.
    """
    if not os.path.exists(csv_path):
        return None
    with open(csv_path, newline='') as f:
        header_line = f.readline()
        if not header_line.startswith(_PCUBE_PROVENANCE_PREFIX):
            return None
        for row in csv.DictReader(f):
            if int(row['segment_id']) == -1:
                return row
    return None


def pd_hr_survey_plot(
    obs_db_path: str,
    survey_csv: str,
    save_path_prefix: str | None = None,
    soft_state_sources: Sequence[str] = (),
    hard_state_sources: Sequence[str] = (),
    mdp_flag_pct: float = 15.0,
) -> dict:
    """Render the five standard PD-vs-HR survey reviewer figures.

    Reads `bhxrb_pd_hardness_survey.csv` (`survey_csv`, written by
    `pd_hr_survey`) and, for the feasibility panel only, each usable obsid's
    own `<obs_db_path>/<obsid>/pcube_timeseg_<obsid>.csv` integrated row (for
    `exposure_s`, which `survey_csv` does not carry). Writes exactly five
    files, `{save_path_prefix}_<suffix>.png`:

    - `pd_vs_hr_combined`: PD_pct vs HR scatter (with PD_err/HR_err
      crosses), color-coded by source, over every detected row; annotated
      with the combined Spearman rho/p.
    - `pd_vs_hr_per_source`: the same scatter split one panel per source
      (every source with a usable row, not just detected ones), each
      annotated with its own Spearman rho/p when it has >= 3 detections.
    - `pd_distribution_comparison`: detected PD_pct histograms for
      `soft_state_sources` vs `hard_state_sources` (both caller-supplied,
      never hardcoded -- matching `pd_hr_survey`'s own KS-test grouping),
      annotated with the KS statistic/p.
    - `hr_detection_bias`: HR histograms for the detected vs non-detected
      groups (every usable row), annotated with the Mann-Whitney U/p.
    - `feasibility_mdp99_exposure`: MDP99 (%) vs exposure (ks) scatter for
      every usable obsid with a readable `pcube_timeseg` CSV, color-coded by
      detection, with a `mdp_flag_pct` reference line.

    An empty group (no detections, no soft/hard match, an obsid whose
    `pcube_timeseg` CSV is unreadable for the feasibility panel, ...) is
    never degraded to a misleading empty plot: the affected panel carries an
    explicit annotation instead, noted in `stdout` where relevant. `survey_csv`
    is a required upstream product: missing, or present but failing
    `_load_survey_rows`'s schema check, is a hard error -- `exit_code=1` with
    `stderr` set and no files written. A `survey_csv` with zero usable rows
    (every row's `note` is set, i.e. every requested obsid was skipped
    upstream) is the same hard error, since there is nothing left to plot.

    Parameters
    ----------
    obs_db_path : str
        The IXPE mission's obs-db directory; each obsid's per-segment CSV
        (for the feasibility panel's `exposure_s`) is read from
        `<obs_db_path>/<obsid>/pcube_timeseg_<obsid>.csv`.
    survey_csv : str
        Path to `bhxrb_pd_hardness_survey.csv`, written by `pd_hr_survey`.
    save_path_prefix : str or None
        Base path; the five suffixes above plus `.png` are appended.
        Auto-generated timestamped prefix in `os.getcwd()` when None.
    soft_state_sources, hard_state_sources : Sequence[str]
        Source-name labels (matched against `survey_csv`'s `source` column)
        forming the `pd_distribution_comparison` groups -- the same
        caller-supplied split `pd_hr_survey`'s KS test uses. Default empty.
    mdp_flag_pct : float
        MDP99 (percent) reference line drawn on `feasibility_mdp99_exposure`
        -- the same data-quality threshold `pd_hr_survey` flags rows with.
        Default 15.0.

    Returns
    -------
    dict with 'figures' (list of 5 `FigureRef` dicts), 'plotted' (the exact
    values drawn in each panel -- 'combined_points', 'per_source_points',
    'distribution_groups', 'detection_bias_hr', 'mdp99_exposure' -- for
    differential verification), 'stdout', 'exit_code' (1 with 'stderr' and
    no files written when `survey_csv` is missing/unreadable/all-skipped).
    """
    if not os.path.exists(survey_csv):
        return {
            'exit_code': 1, 'stdout': '',
            'stderr': (f'survey_csv not found: {survey_csv}. It is a required '
                       f'bhxrb_pd_hr_survey product; run that tool first.'),
        }
    try:
        rows = _load_survey_rows(survey_csv)
    except (OSError, csv.Error, KeyError) as exc:
        return {
            'exit_code': 1, 'stdout': '',
            'stderr': f'survey_csv unreadable ({survey_csv}): {exc!r}',
        }

    usable = [r for r in rows if not r['note']]
    if not usable:
        return {
            'exit_code': 1, 'stdout': '',
            'stderr': f'No usable (PD+HR) rows in survey_csv: {survey_csv}.',
        }

    if save_path_prefix is None:
        ts = _time.strftime('%Y%m%dT%H%M%S')
        save_path_prefix = os.path.join(os.getcwd(), f'pd_hr_survey_plot_{ts}')
    os.makedirs(os.path.dirname(os.path.abspath(save_path_prefix)) or '.', exist_ok=True)

    mime = 'image/png'
    savefig_kwargs = {'metadata': {'Software': _PRODUCED_BY}}

    def _save(fig, suffix: str, caption: str) -> dict:
        path = f'{save_path_prefix}_{suffix}.png'
        fig.savefig(path, dpi=150, bbox_inches='tight', **savefig_kwargs)
        plt.close(fig)
        sha = hashlib.sha256(Path(path).read_bytes()).hexdigest()
        return {
            'path': path, 'caption': caption, 'mime_type': mime, 'produced_by': _PRODUCED_BY,
            'created_at': _time.strftime('%Y-%m-%dT%H:%M:%SZ', _time.gmtime()),
            'sha256': sha,
        }

    figures = []
    notes = []
    detected = [r for r in usable if r['detected']]

    # -- 1. Combined PD vs HR (detected only) --------------------------------
    fig1, ax1 = plt.subplots(figsize=(6, 5))
    if detected:
        sources = sorted(set(r['source'] for r in detected))
        cmap = plt.get_cmap('tab10')
        colors = {s: cmap(i % 10) for i, s in enumerate(sources)}
        for s in sources:
            srows = [r for r in detected if r['source'] == s]
            ax1.errorbar([r['HR'] for r in srows], [r['PD_pct'] for r in srows],
                          xerr=[r['HR_err'] for r in srows], yerr=[r['PD_err'] for r in srows],
                          fmt='o', ms=5, elinewidth=0.7, capsize=2, color=colors[s], label=s)
        combined = _spearman([r['PD_pct'] for r in detected], [r['HR'] for r in detected])
        if not np.isnan(combined['rho']):
            ax1.annotate(f"combined Spearman rho={combined['rho']:.3f} p={combined['p']:.3e} (n={combined['n']})",
                         xy=(0.02, 0.98), xycoords='axes fraction', va='top', fontsize=8)
        ax1.legend(fontsize='small')
    else:
        ax1.annotate('No detected obsids', xy=(0.5, 0.5), xycoords='axes fraction', ha='center', va='center')
        ax1.set_xticks([]); ax1.set_yticks([])
    ax1.set_xlabel('HR'); ax1.set_ylabel('PD (%)')
    fig1.tight_layout()
    figures.append(_save(fig1, 'pd_vs_hr_combined', f'Combined PD vs HR for {len(detected)} detected obsid(s)'))
    combined_points = [{'obsid': r['obsid'], 'source': r['source'], 'PD_pct': r['PD_pct'], 'HR': r['HR']}
                       for r in detected]

    # -- 2. Per-source PD vs HR -----------------------------------------------
    sources_all = sorted(set(r['source'] for r in usable))
    n = len(sources_all)
    fig2, axes2 = plt.subplots(ncols=n, figsize=(4 * n, 4), squeeze=False)
    per_source_points: dict = {}
    for i, s in enumerate(sources_all):
        ax = axes2[0, i]
        srows = [r for r in detected if r['source'] == s]
        per_source_points[s] = [{'obsid': r['obsid'], 'PD_pct': r['PD_pct'], 'HR': r['HR']} for r in srows]
        if srows:
            ax.errorbar([r['HR'] for r in srows], [r['PD_pct'] for r in srows],
                        xerr=[r['HR_err'] for r in srows], yerr=[r['PD_err'] for r in srows],
                        fmt='o', ms=5, elinewidth=0.7, capsize=2, color='tab:blue')
            res = _spearman([r['PD_pct'] for r in srows], [r['HR'] for r in srows])
            if not np.isnan(res['rho']):
                ax.annotate(f"rho={res['rho']:.3f}\np={res['p']:.3e}\n(n={res['n']})",
                            xy=(0.02, 0.98), xycoords='axes fraction', va='top', fontsize=7)
        else:
            ax.annotate('No detected obsids', xy=(0.5, 0.5), xycoords='axes fraction', ha='center', va='center')
            ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(s)
        ax.set_xlabel('HR')
        ax.set_ylabel('PD (%)' if i == 0 else '')
    fig2.tight_layout()
    figures.append(_save(fig2, 'pd_vs_hr_per_source', f'Per-source PD vs HR for {", ".join(sources_all)}'))

    # -- 3. PD distribution: soft vs hard state sources ------------------------
    soft_set, hard_set = set(soft_state_sources), set(hard_state_sources)
    soft_pd = [r['PD_pct'] for r in detected if r['source'] in soft_set]
    hard_pd = [r['PD_pct'] for r in detected if r['source'] in hard_set]
    fig3, ax3 = plt.subplots(figsize=(6, 4.5))
    if soft_pd or hard_pd:
        if soft_pd:
            ax3.hist(soft_pd, bins=10, alpha=0.6, color='tab:blue', edgecolor='k', label=f'soft (n={len(soft_pd)})')
        if hard_pd:
            ax3.hist(hard_pd, bins=10, alpha=0.6, color='tab:red', edgecolor='k', label=f'hard (n={len(hard_pd)})')
        ks = _ks_two_sample(soft_pd, hard_pd)
        if not np.isnan(ks['statistic']):
            ax3.annotate(f"KS D={ks['statistic']:.3f} p={ks['p']:.3e}",
                         xy=(0.02, 0.98), xycoords='axes fraction', va='top', fontsize=8)
        ax3.legend(fontsize='small')
    else:
        ax3.annotate('No soft/hard source groups supplied or matched',
                     xy=(0.5, 0.5), xycoords='axes fraction', ha='center', va='center')
        ax3.set_xticks([]); ax3.set_yticks([])
    ax3.set_xlabel('PD (%)'); ax3.set_ylabel('N obsids')
    fig3.tight_layout()
    figures.append(_save(fig3, 'pd_distribution_comparison',
                          f'Detected PD distribution: soft (n={len(soft_pd)}) vs hard (n={len(hard_pd)})'))

    # -- 4. HR detection bias ---------------------------------------------------
    det_hr = [r['HR'] for r in usable if r['detected']]
    nondet_hr = [r['HR'] for r in usable if not r['detected']]
    fig4, ax4 = plt.subplots(figsize=(6, 4.5))
    if det_hr or nondet_hr:
        if det_hr:
            ax4.hist(det_hr, bins=10, alpha=0.6, color='tab:green', edgecolor='k', label=f'detected (n={len(det_hr)})')
        if nondet_hr:
            ax4.hist(nondet_hr, bins=10, alpha=0.6, color='tab:gray', edgecolor='k',
                     label=f'non-detected (n={len(nondet_hr)})')
        mwu = _mann_whitney(det_hr, nondet_hr)
        if not np.isnan(mwu['statistic']):
            ax4.annotate(f"Mann-Whitney U={mwu['statistic']:.3f} p={mwu['p']:.3e}",
                         xy=(0.02, 0.98), xycoords='axes fraction', va='top', fontsize=8)
        ax4.legend(fontsize='small')
    else:
        ax4.annotate('No obsids with both PD and HR', xy=(0.5, 0.5), xycoords='axes fraction', ha='center', va='center')
        ax4.set_xticks([]); ax4.set_yticks([])
    ax4.set_xlabel('HR'); ax4.set_ylabel('N obsids')
    fig4.tight_layout()
    figures.append(_save(fig4, 'hr_detection_bias',
                          f'HR distribution: detected (n={len(det_hr)}) vs non-detected (n={len(nondet_hr)})'))

    # -- 5. Feasibility: MDP99 vs exposure ---------------------------------------
    mdp_exp_points = []
    for r in usable:
        csv_path = os.path.join(obs_db_path, r['obsid'], f"pcube_timeseg_{r['obsid']}.csv")
        integ = _read_integrated_row(csv_path)
        if integ is None:
            notes.append(f"{r['obsid']}: skipped from feasibility panel (no readable pcube_timeseg CSV at {csv_path})")
            continue
        mdp_exp_points.append({'obsid': r['obsid'], 'MDP99': r['MDP99'],
                               'exposure_s': float(integ['exposure_s']), 'detected': r['detected']})
    fig5, ax5 = plt.subplots(figsize=(6, 4.5))
    if mdp_exp_points:
        det_pts = [p for p in mdp_exp_points if p['detected']]
        nondet_pts = [p for p in mdp_exp_points if not p['detected']]
        if det_pts:
            ax5.scatter([p['exposure_s'] / 1000.0 for p in det_pts], [p['MDP99'] for p in det_pts],
                       color='tab:green', s=25, label=f'detected (n={len(det_pts)})')
        if nondet_pts:
            ax5.scatter([p['exposure_s'] / 1000.0 for p in nondet_pts], [p['MDP99'] for p in nondet_pts],
                       color='tab:gray', s=25, label=f'non-detected (n={len(nondet_pts)})')
        ax5.axhline(mdp_flag_pct, color='k', ls='--', lw=1, label=f'mdp_flag_pct={mdp_flag_pct:g}%')
        ax5.legend(fontsize='small')
    else:
        ax5.annotate('No obsid has a readable pcube_timeseg CSV for exposure',
                     xy=(0.5, 0.5), xycoords='axes fraction', ha='center', va='center')
        ax5.set_xticks([]); ax5.set_yticks([])
    ax5.set_xlabel('Exposure (ks)'); ax5.set_ylabel('MDP99 (%)')
    fig5.tight_layout()
    figures.append(_save(fig5, 'feasibility_mdp99_exposure', f'MDP99 vs exposure for {len(mdp_exp_points)} obsid(s)'))

    stdout_lines = [f'Wrote {len(figures)} figure(s) with prefix {save_path_prefix} (png).']
    stdout_lines.extend(notes)
    return {
        'figures': figures,
        'plotted': {
            'combined_points': combined_points,
            'per_source_points': per_source_points,
            'distribution_groups': {'soft': soft_pd, 'hard': hard_pd},
            'detection_bias_hr': {'detected': det_hr, 'non_detected': nondet_hr},
            'mdp99_exposure': mdp_exp_points,
        },
        'stdout': '\n'.join(stdout_lines),
        'exit_code': 0,
    }
