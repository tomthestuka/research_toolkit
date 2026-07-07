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
import os
from typing import Sequence

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
