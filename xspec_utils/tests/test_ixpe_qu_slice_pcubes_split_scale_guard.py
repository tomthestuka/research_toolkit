"""Tests for the ixpe-qu-slice-split-scale-guard fix (xspec_utils/lib/red_ixpe.py:
``_sorted_gap_split_indices``, ``subobs_groups``, ``_resolve_subobs_groups``,
``qu_slice_pcubes``, plus the ``ixpe_qu_slice_pcubes`` wrapper in
xspec_utils/skills/session.py).

DEFECT: the N-way GTI-gap splitter accepted a merge boundary on a pure
sorted-gap ratio cliff, with no absolute scale guard and no cap. In a long,
gap-riddled '...99' campaign the biggest proportional drop in the ranked gap
list sits between the routine orbital-visibility-gap tier (order of IXPE's
~96.6 min orbit, at most a few ks) and the sub-second dropout tier, so every
orbital gap got accepted as a "merge boundary" -- one real 386 ks '...99'
ObsID (04002899) produced ~230 single-slice sub-observations instead of a
handful.

FIX: a cliff-accepted gap must ALSO clear two scale floors
(``split_gap_min_s``, ``split_gap_min_frac``) before becoming a boundary,
and at most ``max_subobs`` sub-observations are ever emitted (the largest
surviving boundaries are kept when more clear the floors than the cap
allows). The ratio-cliff computation itself, and the '...99'-only
``_should_split_subobs`` gate, are UNCHANGED.

Real, already-local Cyg X-1 IXPE ObsIDs used as fixtures: 04002899 (the
386 ks '...99' ObsID the defect was reported on -- read here at the GTI
level only, never reduced: a full reduction is far too expensive for the
gate budget), 03002599 (a smaller, well-understood '...99' merge: a single
~5.4e6s literature-merge gap against ~2.8e3s ordinary gaps), 02008201 (an
ordinary '...01' monitoring pointing), and 01002901 (single continuous
pointing, used for the reproducibility check).

Measured real-GTI numbers (read-only, no reduction -- see the module
docstrings of each real-data test below for where each number is used):

04002899: n_gti=260, total_elapsed_span_s=1032987.9369999766,
  total_good_time_s=387109.0391755104 (~387 ks, matches the reported "386 ks"
  ObsID). Top gaps (idx, gap_s): (121, 156810.0), (198, 151544.51999998093),
  (108, 3010.36), (80, 2989.42), (50, 2965.49), (15, 2941.56), (33, 2907.13),
  (61, 2883.21). median_gap_s=736.68. threshold=max(86400,0.10*span)=
  103298.79. n_eligible=2. _sorted_gap_split_indices DEFAULT=[121, 198]
  (len=2). _resolve_subobs_groups DEFAULT groups len=3. Floors/cap disabled
  (pre-fix ratio-cliff-only behaviour): len=182 (183 groups) -- the reported
  over-split, reproduced exactly against this real ObsID.

03002599: n_gti=207, total_elapsed_span_s=5612385.0,
  total_good_time_s=119436.03970080614. Dominant gap (idx 106)=
  5395962.576999992s; next-largest=2798.38s. median_gap_s=3.92.
  threshold=561238.5. n_eligible=1. DEFAULT=[106] (len=1), groups len=2.
  Floors/cap disabled: identical, [106] (len=1) -- only one gap ever clears
  the ratio cliff here regardless of the floors.

02008201: n_gti=54, total_elapsed_span_s=43964.49700000882,
  total_good_time_s=23121.368685036898. Top gaps: 2682.67s, 2435.13s,
  2126.0s, 2126.0s, 2125.0s, 2125.0s, 2125.0s, 876.39s. median_gap_s=1.51.
  threshold=86400.0 (absolute floor dominates; 10% of span is only
  4396.45). n_eligible=0. DEFAULT=[] (len=0), groups len=1. Floors/cap
  disabled: len=51 (52 groups) -- confirms the pure ratio cliff over-splits
  even an ordinary, non-'...99' pointing's own real GTI structure; only the
  '...99' naming gate (unchanged by this fix) keeps 02008201 itself from
  ever reaching the split machinery at all.
"""

from __future__ import annotations

import csv
import glob
import os
import shutil

import numpy as np
import pytest
from astropy.io import fits

from obsdb_catalog import neutral
from xspec_utils.lib import red_ixpe
from xspec_utils.skills import session as ixpe_session

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
OBSID_386KS_MERGED = "04002899"   # the reported ~230-fragment over-split ObsID
OBSID_MERGED = "03002599"        # smaller, well-understood '...99' merge -- splits a/b
OBSID_ORDINARY = "02008201"      # ordinary '...01' monitoring pointing -- never split
OBSID_SINGLE = "01002901"        # single continuous pointing -- reproducibility fixture


def _real_gti(obsid: str) -> tuple[list[float], list[float]]:
    """Read one DU's real GTI START/STOP (no reduction) for ``obsid``, preferring
    the '_rej' variant, falling back to plain evt2_v*.fits (mirrors the tool's
    own file-discovery convention)."""
    evt_dir = os.path.join(REAL_IXPE_DIR, obsid, "event_l2")
    files = sorted(glob.glob(os.path.join(evt_dir, f"ixpe{obsid}_det1_evt2_v[0-9][0-9]_rej.fits")))
    if not files:
        files = sorted(glob.glob(os.path.join(evt_dir, f"ixpe{obsid}_det1_evt2_v[0-9][0-9].fits")))
    if not files:
        pytest.skip(f"real fixture ObsID {obsid} not local under {evt_dir}")
    with fits.open(files[0]) as hdul:
        gti = hdul["GTI"].data
        return [float(x) for x in gti["START"]], [float(x) for x in gti["STOP"]]


def _copy_real_event_files(dest_dir: str, obsid: str) -> list[str]:
    src_evt_dir = os.path.join(REAL_IXPE_DIR, obsid, "event_l2")
    dst_evt_dir = os.path.join(dest_dir, obsid, "event_l2")
    os.makedirs(dst_evt_dir, exist_ok=True)
    file_list = []
    for det in (1, 2, 3):
        rej = sorted(glob.glob(os.path.join(src_evt_dir, f"ixpe{obsid}_det{det}_evt2_v[0-9][0-9]_rej.fits")))
        plain = sorted(glob.glob(os.path.join(src_evt_dir, f"ixpe{obsid}_det{det}_evt2_v[0-9][0-9].fits")))
        src = rej[0] if rej else (plain[0] if plain else None)
        if src is None:
            pytest.skip(f"real fixture ObsID {obsid} DU{det} not local under {src_evt_dir}")
        dst = os.path.join(dst_evt_dir, os.path.basename(src))
        shutil.copyfile(src, dst)
        file_list.append(dst)
    return sorted(file_list)


def _gti_from_gaps(gaps: list[float], width: float = 1.0, t0: float = 0.0) -> tuple[list[float], list[float]]:
    """Build a synthetic GTI start/stop pair with exactly ``gaps`` (in order)
    as the inter-GTI gaps, each GTI itself ``width`` seconds wide."""
    starts = [t0]
    stops = [t0 + width]
    for g in gaps:
        t0 = stops[-1] + g
        starts.append(t0)
        stops.append(t0 + width)
    return starts, stops


def _ten_day_scale_gaps_fixture() -> tuple[list[float], list[float], list[float]]:
    """Ten day-scale gaps (90000s..990000s in 100000s steps) plus two
    sub-minute gaps (10s, 15s). All ten mutually-mild-ratio big gaps clear
    the ratio cliff together (the real cliff sits below all ten, against the
    trailing sub-minute tier: 90000/15 = 6000, decisively above
    split_gap_factor=10). Total elapsed span is ~5.400038e6s, so the 10%
    fractional floor is ~5.40004e5s -- verified (module-level, real
    arithmetic, not asserted against): exactly FIVE of the ten (590000,
    690000, 790000, 890000, 990000) clear BOTH the 1-day absolute floor and
    the 10% fractional floor at their defaults; the other five (90000..
    490000) clear only the absolute floor. Returns (gti_start, gti_stop,
    big_gaps).
    """
    big_gaps = [90000.0 + i * 100000.0 for i in range(10)]  # 90000..990000
    small_gaps = [10.0, 15.0]
    gti_start, gti_stop = _gti_from_gaps(big_gaps + small_gaps)
    return gti_start, gti_stop, big_gaps


# --------------------------------------------------------------------------- #
# Capability (real GTI, no reduction): the 386 ks '...99' ObsID that produced
# ~230 fragments must now resolve to exactly the measured group count.
# --------------------------------------------------------------------------- #

def test_capability_04002899_real_gti_bounded_by_max_subobs():
    gti_start, gti_stop = _real_gti(OBSID_386KS_MERGED)
    groups = red_ixpe._resolve_subobs_groups(OBSID_386KS_MERGED, gti_start, gti_stop, 10.0, None)

    assert len(groups) <= 4  # max_subobs default -- fails loudly against the ~230-fragment defect
    # Measured (see module docstring): 2 accepted boundaries -> 3 groups. A
    # degenerate always-one-group implementation would also pass "<= 4"; the
    # exact count, together with the 03002599 "== 2" regression below, pins
    # the real behaviour from both sides.
    assert len(groups) == 3

    total_span = gti_stop[-1] - gti_start[0]
    scale_floor = max(86400.0, 0.10 * total_span)
    for prev, nxt in zip(groups, groups[1:]):
        gap = nxt["tstart"] - prev["tstop"]
        assert gap >= 86400.0
        assert gap >= 0.10 * total_span
        assert gap >= scale_floor

    # No GTI dropped or duplicated by the partition.
    recon_start = [s for g in groups for s in g["gti_start"]]
    recon_stop = [s for g in groups for s in g["gti_stop"]]
    assert recon_start == gti_start
    assert recon_stop == gti_stop


# --------------------------------------------------------------------------- #
# Differential (real GTI, independent numpy recomputation -- never calls the
# tool's own helper for the eligibility computation itself).
# --------------------------------------------------------------------------- #

def test_differential_04002899_boundaries_subset_of_independent_eligible_set():
    gti_start, gti_stop = _real_gti(OBSID_386KS_MERGED)
    starts = np.array(gti_start)
    stops = np.array(gti_stop)
    gaps = starts[1:] - stops[:-1]
    threshold = max(86400.0, 0.10 * (stops[-1] - starts[0]))
    independent_eligible = set(np.nonzero(gaps >= threshold)[0].tolist())

    returned = set(red_ixpe._sorted_gap_split_indices(gti_start, gti_stop, 10.0))
    # Measured (see module docstring): the two dominant real gaps at indices
    # 121 (156810.0s) and 198 (151544.52s) are the only ones that clear both
    # the ratio cliff and both scale floors.
    assert len(returned) == 2
    assert returned <= independent_eligible  # floors respected
    assert len(returned) <= 3                # max_subobs=4 -> at most 3 boundaries

    if returned:
        smallest_accepted = min(gaps[i] for i in returned)
        for i in independent_eligible - returned:
            # The cliff accepts a PREFIX of the descending-sorted gaps and the
            # cap keeps the LARGEST survivors, so nothing eligible-but-dropped
            # can exceed an accepted gap -- exactly the invariant the
            # ~230-fragment defect violated (it accepted gaps far below the
            # genuine merge-gap tier).
            assert gaps[i] <= smallest_accepted


def test_differential_03002599_single_boundary_equals_independent_argmax():
    gti_start, gti_stop = _real_gti(OBSID_MERGED)
    starts = np.array(gti_start)
    stops = np.array(gti_stop)
    gaps = starts[1:] - stops[:-1]
    independent_argmax = int(np.argmax(gaps))

    returned = red_ixpe._sorted_gap_split_indices(gti_start, gti_stop, 10.0)
    assert returned == [independent_argmax]


# --------------------------------------------------------------------------- #
# Regression (real GTI, no reduction): the two documented real-data outcomes
# must hold at the fix's own defaults.
# --------------------------------------------------------------------------- #

def test_regression_03002599_still_splits_into_exactly_two():
    gti_start, gti_stop = _real_gti(OBSID_MERGED)
    groups = red_ixpe._resolve_subobs_groups(OBSID_MERGED, gti_start, gti_stop, 10.0, None)
    assert len(groups) == 2


def test_regression_02008201_still_yields_exactly_one_whole_span():
    gti_start, gti_stop = _real_gti(OBSID_ORDINARY)
    groups = red_ixpe._resolve_subobs_groups(OBSID_ORDINARY, gti_start, gti_stop, 10.0, None)
    assert len(groups) == 1
    assert groups[0]["tstart"] == gti_start[0] and groups[0]["tstop"] == gti_stop[-1]


# --------------------------------------------------------------------------- #
# Unit tests: the two new floors and the cap, exercised at their DEFAULTS on
# synthetic (unphysical-magnitude but structurally realistic) GTI lists.
# --------------------------------------------------------------------------- #

def test_floors_many_orbital_gaps_in_multiweek_span_no_split():
    """300 routine ~3ks orbital gaps plus one 60 ks outlier (comfortably
    clearing the ratio cliff, ratio=20 > default split_gap_factor=10) in a
    multi-week (~11 day) span: the 60 ks outlier is below BOTH the 1-day
    absolute floor and the 10%-of-span fractional floor, so it must be
    rejected -- no split, at the fix's own defaults."""
    gaps = [3000.0] * 150 + [60000.0] + [3000.0] * 150
    gti_start, gti_stop = _gti_from_gaps(gaps)
    assert red_ixpe._sorted_gap_split_indices(gti_start, gti_stop, 10.0) == []
    assert len(red_ixpe.subobs_groups(gti_start, gti_stop, 10.0)) == 1


def test_defect_reproduction_orbital_plus_subsecond_dropout_tiers():
    """Faithful reproduction of the REPORTED failure mode -- a two-tier gap
    population, not the uniform-gap approximation above. ~150 routine
    orbital-visibility gaps (~3 ks; IXPE's ~96.6 min orbit) interleaved with
    ~150 sub-second telemetry dropouts (~0.05 s), over a ~3-week span.
    PRE-FIX, the max ratio in the sorted gap list sits at the 3000/0.05 =
    6e4 cliff BETWEEN the two tiers, so every orbital gap was accepted as a
    merge boundary -- the reported ~230-fragment over-split."""
    gaps = []
    for _ in range(150):
        gaps.extend([3000.0, 0.05])
    # width=4534.0 makes the 301-GTI total elapsed span come out to ~1.81e6s
    # (~21 days, ~3 weeks) -- the width itself carries no orbital-gap
    # information, it only sets the campaign's overall duration.
    gti_start, gti_stop = _gti_from_gaps(gaps, width=4534.0)

    # At the fix's defaults: nothing clears the 1-day absolute floor -> no split.
    assert red_ixpe._sorted_gap_split_indices(gti_start, gti_stop, 10.0) == []
    assert len(red_ixpe.subobs_groups(gti_start, gti_stop, 10.0)) == 1

    # Floors/cap disabled == the pre-fix code path: the fixture genuinely
    # reproduces the defect, so the FLOORS (not some incidental fixture
    # property) are what suppress it -- without this half, the assertion
    # above would be unfalsifiable.
    prefix_only = red_ixpe._sorted_gap_split_indices(
        gti_start, gti_stop, 10.0,
        split_gap_min_s=0.0, split_gap_min_frac=0.0, max_subobs=None,
    )
    assert len(prefix_only) > 50


def test_floors_one_month_scale_gap_exactly_one_boundary():
    """The same routine-gap background, but the outlier is month-scale
    (~31 days): clears both floors comfortably -- exactly one boundary."""
    gaps = [3000.0] * 150 + [2700000.0] + [3000.0] * 150
    gti_start, gti_stop = _gti_from_gaps(gaps)
    idxs = red_ixpe._sorted_gap_split_indices(gti_start, gti_stop, 10.0)
    assert len(idxs) == 1
    total_span = gti_stop[-1] - gti_start[0]
    gap = gti_start[idxs[0] + 1] - gti_stop[idxs[0]]
    assert gap == pytest.approx(2700000.0)
    assert gap >= 86400.0
    assert gap >= 0.10 * total_span


def test_cap_ten_eligible_day_scale_gaps_at_defaults_keeps_three_largest():
    """At the fix's own DEFAULTS (both floors on): of the ten day-scale
    gaps, exactly five (590000..990000) clear both the 1-day absolute floor
    and the 10%-of-span (~5.4e5s) fractional floor -- see
    ``_ten_day_scale_gaps_fixture``'s docstring for the arithmetic. The cap
    (max_subobs=4) then keeps only the THREE largest of those five."""
    gti_start, gti_stop, big_gaps = _ten_day_scale_gaps_fixture()
    idxs = red_ixpe._sorted_gap_split_indices(gti_start, gti_stop, 10.0)
    assert len(idxs) == 3

    selected_gaps = sorted((gti_start[i + 1] - gti_stop[i] for i in idxs), reverse=True)
    three_largest = sorted(big_gaps, reverse=True)[:3]  # 990000, 890000, 790000
    assert selected_gaps == pytest.approx(three_largest)


def test_cap_disabled_keeps_all_five_floor_eligible_gaps():
    """Same fixture, ``max_subobs=None`` (cap off), both floors still at
    their defaults: all FIVE floor-eligible gaps survive (the cap is the
    only thing STEP 5's default-cap test disables)."""
    gti_start, gti_stop, big_gaps = _ten_day_scale_gaps_fixture()
    idxs = red_ixpe._sorted_gap_split_indices(gti_start, gti_stop, 10.0, max_subobs=None)
    assert len(idxs) == 5

    selected_gaps = sorted((gti_start[i + 1] - gti_stop[i] for i in idxs), reverse=True)
    five_largest = sorted(big_gaps, reverse=True)[:5]  # 990000..590000
    assert selected_gaps == pytest.approx(five_largest)


def test_floors_and_cap_disabled_keeps_all_ten_ratio_cliff_only():
    """Same fixture, both floors AND the cap disabled: the ratio cliff run
    in isolation accepts all ten big gaps (its mutual-ratio structure is
    unaffected by either floor or the cap)."""
    gti_start, gti_stop, big_gaps = _ten_day_scale_gaps_fixture()
    idxs = red_ixpe._sorted_gap_split_indices(
        gti_start, gti_stop, 10.0,
        split_gap_min_s=0.0, split_gap_min_frac=0.0, max_subobs=None,
    )
    assert len(idxs) == 10


def test_max_subobs_le_one_forbids_any_split():
    gti_start, gti_stop = _gti_from_gaps([2700000.0])
    assert red_ixpe._sorted_gap_split_indices(gti_start, gti_stop, 10.0, max_subobs=1) == []
    assert red_ixpe._sorted_gap_split_indices(gti_start, gti_stop, 10.0, max_subobs=0) == []


def test_zero_elapsed_span_never_divides_by_zero_and_never_splits():
    # >= 2 gaps are required to reach the len(ranked) < 2 guard's else-branch
    # at all (a single-gap GTI list short-circuits there and never touches
    # the span/floor code this test targets) -- 3 GTIs give exactly 2 gaps.
    # Degenerate (non-positive) gaps keep their infinite-ratio treatment for
    # the CLIFF (both gaps tie at ratio=inf/inf, so the cliff would accept
    # everything), so the FLOORS must be what rejects them here: a
    # non-positive total_elapsed_span makes the fractional floor infinite
    # (never a ZeroDivisionError), and no gap can ever clear an infinite
    # floor.
    gti_start = [0.0, 0.0, 0.0]
    gti_stop = [0.0, 0.0, 0.0]
    assert red_ixpe._sorted_gap_split_indices(gti_start, gti_stop, 10.0) == []


def test_negative_elapsed_span_never_splits():
    # Same >= 2-gaps requirement as above; a decreasing GTI list is
    # degenerate but must not raise or divide by a negative number.
    gti_start = [10.0, 0.0, 5.0]
    gti_stop = [10.0, 0.0, 5.0]
    assert red_ixpe._sorted_gap_split_indices(gti_start, gti_stop, 10.0) == []


# --------------------------------------------------------------------------- #
# Payload-bound: forced truncation on a real multi-sub-obs reduction. The
# un-truncated branch is covered by
# test_ixpe_qu_slice_pcubes.py::test_skill_wrapper_falls_back_to_plain_evt2_and_records_catalog
# (which already runs a 03002599 reduction at the default max_return_subobs)
# -- no need to spend a second real 03002599 reduction here on it.
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def merged_truncated_result(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("qu_split_guard_truncate")
    _copy_real_event_files(str(tmp_dir), OBSID_MERGED)
    result = ixpe_session.ixpe_qu_slice_pcubes(
        str(tmp_dir), OBSID_MERGED, n_max=1, n_min=1, max_return_subobs=1,
    )
    return result, str(tmp_dir)


def test_payload_bound_truncates_return_but_not_products(merged_truncated_result):
    result, obs_db_path = merged_truncated_result
    assert result["exit_code"] == 0, result.get("stderr")

    n_total = result["n_subobs_total"]
    assert n_total >= 2  # 03002599 splits into >= 2 sub-obs at defaults
    assert result["per_subobs_truncated"] is True
    assert len(result["per_subobs"]) <= 1
    assert len(result["subobs_keys"]) <= 1

    listing_path = result["subobs_listing_path"]
    assert os.path.exists(listing_path)
    with open(listing_path, newline="") as f:
        rows = list(csv.DictReader(f))
    assert {r["subobs_key"] for r in rows} == set(
        p["params"]["subobs_key"]
        for p in neutral.read(obs_db_path)["entries"][OBSID_MERGED]["products"]
        if p["kind"] == "qu_slices"
    )
    for r in rows:
        assert os.path.exists(r["csv_path"]) or os.path.exists(os.path.join(obs_db_path, r["csv_path"]))
        assert int(r["n_slice"]) >= 1

    # stdout is a short, fixed-size summary -- never a per-sub-obs path enumeration.
    assert len(result["stdout"]) < 400
    assert ".csv" not in result["stdout"]
    assert obs_db_path in result["stdout"] or OBSID_MERGED in result["stdout"]

    # Every sub-obs CSV is on disk and catalogued regardless of truncation --
    # the catalog count, NOT len(per_subobs), proves truncation touched only
    # the returned payload.
    csvs_on_disk = glob.glob(os.path.join(obs_db_path, OBSID_MERGED, "qu_slices_*.csv"))
    assert len(csvs_on_disk) == n_total
    entry = neutral.read(obs_db_path)["entries"][OBSID_MERGED]
    products = [p for p in entry["products"] if p["kind"] == "qu_slices"]
    assert len(products) == n_total



# --------------------------------------------------------------------------- #
# Reproducibility: run the TOOL (the wrapper, not the bare lib function)
# twice on real ObsID 01002901, into two independent obs-db directories, and
# compare the written CSVs -- the actual downstream product.
# --------------------------------------------------------------------------- #

def test_repro_qu_slice_pcubes(tmp_path_factory):
    results = []
    for _ in range(2):
        tmp_dir = tmp_path_factory.mktemp("qu_split_guard_repro")
        _copy_real_event_files(str(tmp_dir), OBSID_SINGLE)
        results.append(ixpe_session.ixpe_qu_slice_pcubes(str(tmp_dir), OBSID_SINGLE, n_max=2, n_min=2))

    a, b = results
    assert a["exit_code"] == 0, a.get("stderr")
    assert b["exit_code"] == 0, b.get("stderr")
    assert a["subobs_keys"] == b["subobs_keys"] == [OBSID_SINGLE]
    assert a["per_subobs_truncated"] == b["per_subobs_truncated"] is False
    assert a["n_subobs_total"] == b["n_subobs_total"] == 1

    info_a = a["per_subobs"][OBSID_SINGLE]
    info_b = b["per_subobs"][OBSID_SINGLE]
    assert info_a["n_slice"] == info_b["n_slice"]

    def _read_slice_rows(csv_path):
        with open(csv_path, newline="") as f:
            f.readline()  # provenance comment header
            return list(csv.DictReader(f))

    rows_a = _read_slice_rows(info_a["csv_path"])
    rows_b = _read_slice_rows(info_b["csv_path"])
    assert len(rows_a) == len(rows_b) == info_a["n_slice"]
    for row_a, row_b in zip(rows_a, rows_b):
        assert row_a["slice_id"] == row_b["slice_id"]
        assert float(row_a["t_mean_met"]) == pytest.approx(float(row_b["t_mean_met"]), abs=1e-6)
        for field in ("q", "u", "q_err", "u_err"):
            assert float(row_a[field]) == pytest.approx(float(row_b[field]), abs=1e-9)

    # New diagnostics are present and internally consistent, on both runs.
    for info in (info_a, info_b):
        assert "mdp_target_frac_used" in info and info["mdp_target_frac_used"] > 0
        assert "c_target_counts" in info and info["c_target_counts"] > 0
        assert "n_natural_slices" in info and info["n_natural_slices"] >= 1
        assert isinstance(info["clipped_to_n_max"], bool)
    for res in (a, b):
        assert res["n_gti"] >= 1
        assert res["median_gti_gap_s"] >= 0.0
        assert res["total_elapsed_span_s"] > 0.0
        assert res["split_gaps_s"] == []  # OBSID_SINGLE is not a '...99' ObsID -- never split
