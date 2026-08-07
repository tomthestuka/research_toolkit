"""Tests for xspec_utils' ixpe_qu_slice_pcubes (lib/red_ixpe.py + skills/session.py).

Uses the real, already-local Level-2 event files for IXPE ObsIDs 01002901 and
03002599 (both real validation ObsIDs named in the tool's brief) as fixture
data, copied into an isolated tmp directory per fixture so the reduction runs
on a scratch copy rather than mutating the shared analyst-env data directory.
n_max is kept small (2-3) across the real-data tests -- not the tool's
production default (25) -- to keep the real ixpeobssim pipeline call count
low while still exercising genuine multi-slice, N-way-split, and mean-
arrival-timestamp behaviour on real data.
"""

from __future__ import annotations

import csv
import glob
import json
import os
import shutil

import numpy as np
import pytest
from astropy.io import fits

from obsdb_catalog import neutral
from xspec_utils.lib import red_ixpe
from xspec_utils.skills import session as ixpe_session
from xspec_utils.skills import TOOL_MANIFEST

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
OBSID_SINGLE = "01002901"   # single continuous pointing -- no N-way split
OBSID_MERGED = "03002599"   # literature-merged pointings -- N-way split into a/b
OBSID_ORDINARY_GAPPY = "02008201"  # ordinary monitoring pointing that over-split pre-fix


def _copy_real_event_files(dest_dir: str, obsid: str) -> list[str]:
    """Copy the real, already-local Level-2 DU1/2/3 event files for ``obsid``
    into ``<dest_dir>/<obsid>/event_l2/`` and return the sorted file list.
    Prefers the ``_rej`` variant per DU, falling back to the plain
    ``evt2_v??.fits`` file -- the same discovery convention the tool itself
    uses (and exactly what 03002599's real on-disk layout requires, since it
    carries no ``_rej`` files)."""
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


# --------------------------------------------------------------------------- #
# Fast, pure-logic unit tests: N-way GTI-gap split (outlier-on-sorted-gaps)
# --------------------------------------------------------------------------- #

def test_sorted_gap_split_indices_no_outlier_no_split():
    # All gaps within a factor of ~2 of each other -- no genuine outlier.
    gti_start = [0.0, 10.0, 22.0, 35.0]
    gti_stop = [5.0, 15.0, 27.0, 40.0]
    assert red_ixpe._sorted_gap_split_indices(gti_start, gti_stop, 10.0) == []


def test_sorted_gap_split_indices_single_dominant_gap():
    # One huge merge gap (idx 1) against many small comparable orbital gaps.
    gti_start = [0.0, 10.0, 1000010.0, 1000020.0, 1000032.0]
    gti_stop = [5.0, 15.0, 1000015.0, 1000025.0, 1000037.0]
    assert red_ixpe._sorted_gap_split_indices(gti_start, gti_stop, 10.0) == [1]


def test_sorted_gap_split_indices_two_dominant_gaps():
    # ~2e9s total span: with the new scale floors at their DEFAULTS the
    # first (~1e6s) dominant gap would fail split_gap_min_frac (10% of a
    # ~2e9s span is ~2e8s) even though it clears the ratio cliff -- disable
    # both floors and the cap explicitly to keep testing the ratio cliff
    # (Step 1) in isolation, per split-scale-guard fix.
    gti_start = [0.0, 10.0, 1000010.0, 1000020.0, 2000000030.0, 2000000042.0]
    gti_stop = [5.0, 15.0, 1000015.0, 1000025.0, 2000000035.0, 2000000047.0]
    assert red_ixpe._sorted_gap_split_indices(
        gti_start, gti_stop, 10.0,
        split_gap_min_s=0.0, split_gap_min_frac=0.0, max_subobs=None,
    ) == [1, 3]


def test_sorted_gap_split_indices_two_comparable_dominant_gaps():
    """Regression: a merge with TWO (or more) comparable dominant gaps must
    accept BOTH, even though their mutual ratio never clears
    split_gap_factor -- the max-ratio cliff sits below them (against the
    ordinary orbital gaps), not between them. A naive strictly-descending
    outlier chain fails this outright (the very first comparison,
    ~5e6/~4.9e6, never clears 10x, so it would return no split at all)."""
    gti_start = [0.0, 10.0, 5000000.0, 5000010.0, 9900000.0, 9900012.0]
    gti_stop = [5.0, 15.0, 5000005.0, 5000015.0, 9900005.0, 9900017.0]
    idxs = red_ixpe._sorted_gap_split_indices(gti_start, gti_stop, 10.0)
    assert len(idxs) == 2
    assert len(red_ixpe.subobs_groups(gti_start, gti_stop, 10.0)) == 3


def test_sorted_gap_split_indices_tiered_prefers_largest_cliff():
    """A tiered gap structure (sorted gaps 5e6, 4e5, 2e3 -- adjacent ratios
    12.5 then 200) must prefer the LARGEST cliff (the 200x drop), accepting
    the top TWO gaps, not just the single largest gap (which only clears the
    smaller 12.5x cliff). The second gap (4e5s) is below split_gap_min_frac's
    DEFAULT threshold for this ~5.4e6s total span, so both floors and the
    cap are explicitly disabled to keep testing the ratio cliff (Step 1) in
    isolation, per split-scale-guard fix."""
    gti_start = [0.0, 5000000.0, 5400000.0, 5402000.0]
    gti_stop = [0.0, 4999999.0, 5399999.0, 5402001.0]
    idxs = red_ixpe._sorted_gap_split_indices(
        gti_start, gti_stop, 10.0,
        split_gap_min_s=0.0, split_gap_min_frac=0.0, max_subobs=None,
    )
    assert idxs == [0, 1]
    assert len(red_ixpe.subobs_groups(
        gti_start, gti_stop, 10.0,
        split_gap_min_s=0.0, split_gap_min_frac=0.0, max_subobs=None,
    )) == 3


def test_sorted_gap_split_indices_fewer_than_two_gtis():
    assert red_ixpe._sorted_gap_split_indices([0.0], [5.0], 10.0) == []
    assert red_ixpe._sorted_gap_split_indices([], [], 10.0) == []


def test_subobs_groups_partitions_at_split_indices():
    gti_start = [0.0, 10.0, 1000010.0, 1000020.0]
    gti_stop = [5.0, 15.0, 1000015.0, 1000025.0]
    groups = red_ixpe.subobs_groups(gti_start, gti_stop, 10.0)
    assert len(groups) == 2
    assert groups[0]["tstart"] == 0.0 and groups[0]["tstop"] == 15.0
    assert groups[1]["tstart"] == 1000010.0 and groups[1]["tstop"] == 1000025.0


def test_subobs_groups_single_group_when_no_split():
    gti_start = [0.0, 10.0, 22.0]
    gti_stop = [5.0, 15.0, 27.0]
    groups = red_ixpe.subobs_groups(gti_start, gti_stop, 10.0)
    assert len(groups) == 1
    assert groups[0]["tstart"] == 0.0 and groups[0]["tstop"] == 27.0


def test_subobs_key_convention():
    assert red_ixpe._subobs_key("X", 0, 1) == "X"
    assert red_ixpe._subobs_key("X", 0, 2) == "Xa"
    assert red_ixpe._subobs_key("X", 1, 2) == "Xb"
    assert red_ixpe._subobs_key("X", 2, 3) == "Xc"


def test_real_03002599_gti_split_matches_known_merge_gap():
    """Grounds split_gap_factor=10.0's default against 03002599's real GTI
    structure directly (no reduction): a single ~5.4e6s literature-merge gap
    against every other real inter-GTI gap (<=2.8e3s) -- decisively over the
    default outlier ratio -- must isolate exactly one split."""
    evt_dir = os.path.join(REAL_IXPE_DIR, OBSID_MERGED, "event_l2")
    files = sorted(glob.glob(os.path.join(evt_dir, f"ixpe{OBSID_MERGED}_det1_evt2_v[0-9][0-9].fits")))
    if not files:
        pytest.skip(f"real fixture ObsID {OBSID_MERGED} not local under {evt_dir}")
    with fits.open(files[0]) as hdul:
        gti = hdul["GTI"].data
        gti_start = [float(x) for x in gti["START"]]
        gti_stop = [float(x) for x in gti["STOP"]]
    cuts = red_ixpe._sorted_gap_split_indices(gti_start, gti_stop, 10.0)
    assert len(cuts) == 1
    groups = red_ixpe.subobs_groups(gti_start, gti_stop, 10.0)
    assert len(groups) == 2
    gap = groups[1]["tstart"] - groups[0]["tstop"]
    assert gap > 1e6  # the real ~5.4e6s literature-merge gap


# --------------------------------------------------------------------------- #
# Fix-under-test: the split-gate must fire ONLY on '...99' merged ObsIDs,
# never on an ordinary monitoring pointing's own (possibly large) GTI gaps.
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("obsid,split_multiobs,expected", [
    ("03002599", None, True),
    ("02001899", None, True),
    ("02008201", None, False),
    ("02008401", None, False),
    ("01002901", None, False),
    ("02008201", True, True),
    ("03002599", False, False),
])
def test_should_split_subobs_capability_gate(obsid, split_multiobs, expected):
    assert red_ixpe._should_split_subobs(obsid, split_multiobs) is expected


def test_real_02008201_gate_prevents_over_split_ungated_path_still_over_splits():
    """Differential, real-data regression for the reported crash: a cheap GTI-
    header read (no ixpeobssim reduction) on the real 02008201 fixture shows
    the UNGATED subobs_groups(...) call with the ratio cliff run in
    isolation (both scale floors and the cap explicitly disabled --
    reproducing what qu_slice_pcubes did before either the '...99' gate fix
    or the split-scale-guard floors) over-splits an ordinary monitoring
    pointing into dozens of groups -- enough to have driven the old
    chr(ord('a')+index) past 26 and emitted a non-ASCII byte -- while the new
    GATED path (_resolve_subobs_groups, which _should_split_subobs judges
    this obsid as NOT merged) yields exactly one whole-span sub-observation.
    """
    evt_dir = os.path.join(REAL_IXPE_DIR, OBSID_ORDINARY_GAPPY, "event_l2")
    files = sorted(glob.glob(os.path.join(evt_dir, f"ixpe{OBSID_ORDINARY_GAPPY}_det1_evt2_v[0-9][0-9]_rej.fits")))
    if not files:
        pytest.skip(f"real fixture ObsID {OBSID_ORDINARY_GAPPY} not local under {evt_dir}")
    with fits.open(files[0]) as hdul:
        gti = hdul["GTI"].data
        gti_start = [float(x) for x in gti["START"]]
        gti_stop = [float(x) for x in gti["STOP"]]

    ungated = red_ixpe.subobs_groups(
        gti_start, gti_stop, 10.0,
        split_gap_min_s=0.0, split_gap_min_frac=0.0, max_subobs=None,
    )
    assert len(ungated) > 1  # the pre-fix behaviour: a real over-split

    gated = red_ixpe._resolve_subobs_groups(OBSID_ORDINARY_GAPPY, gti_start, gti_stop, 10.0, None)
    assert len(gated) == 1
    assert gated[0]["tstart"] == gti_start[0] and gated[0]["tstop"] == gti_stop[-1]


def test_real_02008201_ungated_default_floors_also_prevent_over_split():
    """Belt-and-suspenders regression at the split-scale-guard fix's own
    DEFAULTS: on this same real, ordinary monitoring pointing, even the
    UNGATED subobs_groups(...) call (no '...99' naming-convention gate at
    all) no longer over-splits once the new absolute/fractional scale floors
    are applied at their defaults -- every one of 02008201's real inter-GTI
    gaps is order ~1-3 ks, far below the 1-day absolute floor."""
    evt_dir = os.path.join(REAL_IXPE_DIR, OBSID_ORDINARY_GAPPY, "event_l2")
    files = sorted(glob.glob(os.path.join(evt_dir, f"ixpe{OBSID_ORDINARY_GAPPY}_det1_evt2_v[0-9][0-9]_rej.fits")))
    if not files:
        pytest.skip(f"real fixture ObsID {OBSID_ORDINARY_GAPPY} not local under {evt_dir}")
    with fits.open(files[0]) as hdul:
        gti = hdul["GTI"].data
        gti_start = [float(x) for x in gti["START"]]
        gti_stop = [float(x) for x in gti["STOP"]]

    ungated_default_floors = red_ixpe.subobs_groups(gti_start, gti_stop, 10.0)
    assert len(ungated_default_floors) == 1


# --------------------------------------------------------------------------- #
# Fix-under-test: _subobs_key must stay pure-ASCII past a 26-way split.
# --------------------------------------------------------------------------- #

def test_subobs_key_ascii_safe_past_26_groups():
    n_groups = 45
    keys = [red_ixpe._subobs_key("X", i, n_groups) for i in range(41)]
    for key in keys:
        assert key.isascii()
        assert all(31 < ord(c) < 127 for c in key[1:])  # no control characters
    assert keys[25] == "Xz"
    assert keys[26] == "Xaa"
    assert keys[27] == "Xab"


# --------------------------------------------------------------------------- #
# Fast, pure-logic unit tests: r_proxy / background decision
# --------------------------------------------------------------------------- #

def test_r_proxy_zero_exposure_is_zero_not_a_crash():
    assert red_ixpe.r_proxy(100, 0.0, 1.0) == 0.0


def test_r_proxy_known_value():
    # 100 events / (1000s * pi*1^2 arcmin^2) -- a plain division check.
    import math
    assert red_ixpe.r_proxy(100, 1000.0, 1.0) == pytest.approx(100.0 / (1000.0 * math.pi))


def test_decide_background_override_wins_over_r_proxy():
    assert red_ixpe._decide_background(True, 5.0) is True
    assert red_ixpe._decide_background(False, 0.001) is False


def test_decide_background_faint_subtracts_bright_does_not():
    assert red_ixpe._decide_background(None, 0.5) is True    # faint: r_proxy < 1.0
    assert red_ixpe._decide_background(None, 1.5) is False   # intermediate
    assert red_ixpe._decide_background(None, 5.0) is False   # bright


# --------------------------------------------------------------------------- #
# Fast, pure-logic unit tests: adaptive greedy slice-edge placement
# --------------------------------------------------------------------------- #

def test_greedy_slice_edges_exact_multiple_no_empty_trailing_slice():
    """Regression: closing exactly on the final event must not leave a
    spurious empty trailing slice -- the residual (0..target-1 events) is
    always folded into the previous slice."""
    times = np.arange(0.0, 100.0, 1.0)
    edges, noise_limited, clipped, n_natural = red_ixpe._greedy_slice_edges(times, 0.0, 100.0, 10.0, n_min=2, n_max=25)
    assert edges[0] == 0.0 and edges[-1] == 100.0
    assert all(edges[i] < edges[i + 1] for i in range(len(edges) - 1))
    assert not noise_limited and not clipped
    assert len(edges) - 1 == 10
    assert n_natural == 10


def test_greedy_slice_edges_non_exact_remainder_merged():
    times = np.arange(0.0, 103.0, 1.0)
    edges, noise_limited, clipped, n_natural = red_ixpe._greedy_slice_edges(times, 0.0, 103.0, 10.0, n_min=2, n_max=25)
    assert edges[0] == 0.0 and edges[-1] == 103.0
    assert all(edges[i] < edges[i + 1] for i in range(len(edges) - 1))
    # The 3-event remainder folds into the last full slice -- still 10 slices.
    assert len(edges) - 1 == 10
    assert n_natural == 10


def test_greedy_slice_edges_clips_to_n_max():
    times = np.arange(0.0, 100.0, 1.0)
    edges, noise_limited, clipped, n_natural = red_ixpe._greedy_slice_edges(times, 0.0, 100.0, 10.0, n_min=2, n_max=3)
    assert clipped is True
    assert noise_limited is False
    assert len(edges) - 1 == 3
    assert edges[0] == 0.0 and edges[-1] == 100.0
    assert n_natural == 10  # the natural (pre-clip) slice count, not the clipped 3


def test_greedy_slice_edges_below_n_min_is_noise_limited_not_invented():
    times = np.arange(0.0, 100.0, 1.0)
    edges, noise_limited, clipped, n_natural = red_ixpe._greedy_slice_edges(times, 0.0, 100.0, 1000.0, n_min=2, n_max=25)
    assert noise_limited is True
    assert clipped is False
    assert len(edges) - 1 == 1  # never invents slices the real events don't support
    assert n_natural == 1


def test_greedy_slice_edges_empty_times_is_single_noise_limited_slice():
    edges, noise_limited, clipped, n_natural = red_ixpe._greedy_slice_edges(np.array([]), 0.0, 100.0, 10.0, n_min=2, n_max=25)
    assert edges == [0.0, 100.0]
    assert noise_limited is True and clipped is False
    assert n_natural == 1


# --------------------------------------------------------------------------- #
# Fix-under-test: _reduce_window must surface the real ixpeobssim abort
# diagnostic, not a bare 'Cannot continue, abort.'
# --------------------------------------------------------------------------- #

def test_reduce_window_surfaces_captured_diagnostic_on_abort(monkeypatch):
    """Monkeypatches the ixpeobssim pipeline call to print a sentinel line
    then raise SystemExit -- exactly ixpeobssim's own logging_.abort() bare
    sys.exit pattern (missing IRF, empty selection, missing extension, ...,
    per ixpeobssim tracker #699) -- and asserts the surfaced RuntimeError
    embeds that sentinel line, proving the captured stdout/stderr buffer
    (not os.devnull) is what backs the raised error."""
    sentinel = "SENTINEL: could not find IRF ixpe:obssimBOGUS_alpha075:v13"

    def boom(*args, **kwargs):
        print(sentinel)
        raise SystemExit(1)

    monkeypatch.setattr(red_ixpe.pipeline, "xpselect", boom)
    with pytest.raises(RuntimeError) as excinfo:
        red_ixpe._reduce_window(
            ["dummy_det1_evt2_v01.fits"], "test", None, None, 299.59, 35.2,
            2.0, 8.0, 1.0, False, 3.0, 5.0, "20220101",
        )
    assert sentinel in str(excinfo.value)


# --------------------------------------------------------------------------- #
# Fix-under-test: _reduce_window must clip a windowed xpselect(tmin=,tmax=)
# call to the common per-DU physical span rather than abort on a per-DU
# TSTART/TSTOP mismatch (toolspec-reduce-window-perdu-tstop). OBSID_MISMATCH
# (02008501) is one of the five real ObsIDs the brief reports aborting at
# _reduce_window(<obsid>whole): its DU1 TSTART (203469478.24) precedes
# DU2/DU3's own TSTART (203470050.77) by ~572s, so the DU1-derived
# whole-sub-obs window used to fail ixpeobssim's per-file physical-bounds
# check for DU2/DU3.
# --------------------------------------------------------------------------- #

OBSID_MISMATCH = "02008501"  # real per-DU TSTART mismatch; aborted pre-fix


def _events_tstart_tstop(evt_path: str) -> tuple[float, float]:
    with fits.open(evt_path) as hdul:
        hdr = hdul["EVENTS"].header
        return float(hdr["TSTART"]), float(hdr["TSTOP"])


def test_common_du_bounds_real_mismatched_obsid_is_narrower_than_du1(tmp_path):
    """Grounds _common_du_bounds against the real, reported mismatch: the
    common per-DU intersection must be strictly narrower than DU1's own
    individual span, exactly matching the mismatch the brief reports
    (DU1 TSTART precedes DU2/DU3's own TSTART by ~572s; TSTOP agrees)."""
    file_list = _copy_real_event_files(str(tmp_path), OBSID_MISMATCH)
    du1_tstart, du1_tstop = _events_tstart_tstop(file_list[0])
    tstart_common, tstop_common = red_ixpe._common_du_bounds(file_list)
    assert tstart_common > du1_tstart
    assert tstart_common == pytest.approx(203470050.77, abs=0.1)
    assert tstop_common == pytest.approx(du1_tstop, abs=1e-6)


def test_common_du_bounds_matches_du_span_when_all_dus_agree(tmp_path):
    """Differential: on 01002901 -- whose three DU files share identical
    TSTART/TSTOP -- the per-DU-clipped span must equal each DU's own
    (unclipped) span exactly, so the fix is a genuine no-op for the
    already-working case."""
    file_list = _copy_real_event_files(str(tmp_path), OBSID_SINGLE)
    du_tstart, du_tstop = _events_tstart_tstop(file_list[0])
    tstart_common, tstop_common = red_ixpe._common_du_bounds(file_list)
    assert tstart_common == pytest.approx(du_tstart, abs=1e-6)
    assert tstop_common == pytest.approx(du_tstop, abs=1e-6)


def test_reduce_window_whole_observation_call_never_clipped(monkeypatch):
    """The whole-observation call (tmin=None) must stay exactly as before:
    _common_du_bounds is never even consulted (so it never opens a file --
    a dummy, non-existent path is fine here)."""

    def bounds_boom(*args, **kwargs):
        raise AssertionError("_common_du_bounds must not be called when tmin is None")

    def xpselect_exit(*args, **kwargs):
        raise SystemExit(1)

    monkeypatch.setattr(red_ixpe, "_common_du_bounds", bounds_boom)
    monkeypatch.setattr(red_ixpe.pipeline, "xpselect", xpselect_exit)
    with pytest.raises(RuntimeError):
        red_ixpe._reduce_window(
            ["dummy_det1_evt2_v01.fits"], "wholecall", None, None, 299.59, 35.20,
            2.0, 8.0, 1.0, False, 3.0, 5.0, "20220101",
        )


def test_reduce_window_edge_case_empty_intersection_raises_before_xpselect(tmp_path, monkeypatch):
    """When the requested window lies entirely outside the common per-DU
    span, _reduce_window must raise ValueError WITHOUT ever calling
    xpselect -- an honest failure of that window, never a silent
    empty/invented cube."""
    file_list = _copy_real_event_files(str(tmp_path), OBSID_MISMATCH)
    _, tstop_common = red_ixpe._common_du_bounds(file_list)

    def boom(*args, **kwargs):
        raise AssertionError("xpselect must not be called for an empty clipped window")

    monkeypatch.setattr(red_ixpe.pipeline, "xpselect", boom)
    with pytest.raises(ValueError, match="does not overlap"):
        red_ixpe._reduce_window(
            file_list, "edge", tstop_common + 100.0, tstop_common + 200.0, 299.59, 35.20,
            2.0, 8.0, 1.0, False, 3.0, 5.0, "20220101",
        )


def test_reduce_window_passes_clipped_window_to_xpselect(tmp_path, monkeypatch):
    """Direct capability proof (toolspec check a): when the requested window
    extends past the common per-DU physical span on BOTH edges,
    _reduce_window must pass the CLIPPED window -- exactly
    [max-TSTART, min-TSTOP] from _common_du_bounds -- to xpselect, not the
    wider requested window. Captures the tmin/tmax kwargs xpselect actually
    receives (short-circuiting with SystemExit before xpbin runs) and asserts
    they equal the common per-DU bounds, and are strictly inside the requested
    window."""
    file_list = _copy_real_event_files(str(tmp_path), OBSID_MISMATCH)
    tstart_common, tstop_common = red_ixpe._common_du_bounds(file_list)

    captured = {}

    def fake_xpselect(*args, **kwargs):
        captured["tmin"] = kwargs.get("tmin")
        captured["tmax"] = kwargs.get("tmax")
        raise SystemExit(1)  # short-circuit before xpbin/from_file_list

    monkeypatch.setattr(red_ixpe.pipeline, "xpselect", fake_xpselect)

    req_tmin = tstart_common - 100.0
    req_tmax = tstop_common + 100.0
    with pytest.raises(RuntimeError):
        red_ixpe._reduce_window(
            file_list, "clip", req_tmin, req_tmax, 299.59, 35.20,
            2.0, 8.0, 1.0, False, 3.0, 5.0, "20220101",
        )
    assert captured["tmin"] == pytest.approx(tstart_common, abs=1e-6)
    assert captured["tmax"] == pytest.approx(tstop_common, abs=1e-6)
    # The clip genuinely narrowed the requested window on both edges.
    assert captured["tmin"] > req_tmin
    assert captured["tmax"] < req_tmax


@pytest.fixture(scope="module")
def mismatch_result(tmp_path_factory):
    """One real qu_slice_pcubes call on a scratch copy of OBSID_MISMATCH --
    the ObsID that aborted at _reduce_window(<obsid>whole) pre-fix. n_max
    kept small to bound the real xpselect/xpbin call count."""
    tmp_dir = tmp_path_factory.mktemp("qu_slice_mismatch")
    file_list = _copy_real_event_files(str(tmp_dir), OBSID_MISMATCH)
    result = red_ixpe.qu_slice_pcubes(OBSID_MISMATCH, file_list, n_max=2, n_min=1)
    return result, file_list


def test_capability_mismatch_obsid_completes_and_window_was_clipped(mismatch_result):
    """Capability regression for the reported abort: qu_slice_pcubes must
    now complete on OBSID_MISMATCH and write a real qu_slices CSV, and the
    window it actually needed to reduce (DU1's own GTI-derived whole-sub-obs
    span) must be genuinely narrower once clipped to the common per-DU
    intersection -- proving the clip is what let the reduction succeed,
    rather than it succeeding some other way with the clip a no-op."""
    result, file_list = mismatch_result
    assert result["subobs_keys"] == [OBSID_MISMATCH]  # not a '...99' ObsID -- no N-way split
    info = result["per_subobs"][OBSID_MISMATCH]
    assert os.path.exists(info["csv_path"])
    assert info["whole_counts"] > 0
    assert info["n_slice"] >= 1
    assert len(info["rows"]) == info["n_slice"]

    with fits.open(file_list[0]) as hdul:
        gti = hdul["GTI"].data
        du1_tstart, du1_tstop = float(gti["START"][0]), float(gti["STOP"][-1])
    tstart_common, tstop_common = red_ixpe._common_du_bounds(file_list)
    assert tstart_common > du1_tstart  # the clip genuinely narrowed the requested window
    assert tstop_common <= du1_tstop


def test_repro_mismatch_obsid_stable_across_calls(mismatch_result, tmp_path_factory):
    """Reproducibility, on the previously-aborting ObsID itself: an
    independent second qu_slice_pcubes call on a fresh scratch copy with the
    same arguments gives the same whole-sub-obs PD and the same per-slice
    q/u within tolerance."""
    first, _ = mismatch_result
    tmp_dir = tmp_path_factory.mktemp("qu_slice_mismatch_repro")
    file_list = _copy_real_event_files(str(tmp_dir), OBSID_MISMATCH)
    second = red_ixpe.qu_slice_pcubes(OBSID_MISMATCH, file_list, n_max=2, n_min=1)

    info_a = first["per_subobs"][OBSID_MISMATCH]
    info_b = second["per_subobs"][OBSID_MISMATCH]
    assert info_a["n_slice"] == info_b["n_slice"]
    assert info_a["whole_pd_pct"] == pytest.approx(info_b["whole_pd_pct"], rel=1e-6)
    assert info_a["whole_pa_deg"] == pytest.approx(info_b["whole_pa_deg"], rel=1e-6)
    for row_a, row_b in zip(info_a["rows"], info_b["rows"]):
        for field in ("q", "u", "PD_pct", "PA_deg"):
            assert row_a[field] == pytest.approx(row_b[field], rel=1e-6)


# --------------------------------------------------------------------------- #
# Manifest wiring (AI-facing half)
# --------------------------------------------------------------------------- #

def test_tool_manifest_exposes_class_parameters():
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.ixpe_qu_slice_pcubes"]
    assert entry["module"] == "xspec_utils.skills.session"
    assert entry["fn"] == "ixpe_qu_slice_pcubes"
    props = entry["args_schema"]["properties"]
    for name in ("ra", "dec", "emin_kev", "emax_kev", "src_rad", "bkg",
                 "bkg_inner_rad", "bkg_outer_rad", "mdp_target", "mdp_frac",
                 "mdp_floor", "n_min", "n_max", "split_gap_factor",
                 "split_multiobs", "split_gap_min_s", "split_gap_min_frac",
                 "max_subobs", "max_return_subobs"):
        assert name in props, f"{name} missing from args_schema.properties"
    assert entry["args_schema"]["required"] == ["obs_db_path", "obsid"]
    assert props["n_min"]["type"] == "integer"
    assert props["n_max"]["type"] == "integer"
    assert props["bkg"]["type"] == "boolean"
    assert props["split_multiobs"]["type"] == "boolean"
    assert props["max_subobs"]["type"] == "integer"
    assert props["max_return_subobs"]["type"] == "integer"
    assert props["split_gap_min_s"]["type"] == "number"
    assert props["split_gap_min_frac"]["type"] == "number"
    assert "irf_date" not in props  # auto-selected internally, not AI-facing


# --------------------------------------------------------------------------- #
# Silent-misuse guards -- fail loudly before any reduction is attempted
# --------------------------------------------------------------------------- #

def test_missing_du_raises(tmp_path):
    file_list = _copy_real_event_files(str(tmp_path), OBSID_SINGLE)
    two_du_list = [f for f in file_list if "_det3_" not in f]
    with pytest.raises(ValueError, match="DU3"):
        red_ixpe.qu_slice_pcubes(OBSID_SINGLE, two_du_list)


def test_partial_radec_raises(tmp_path):
    file_list = _copy_real_event_files(str(tmp_path), OBSID_SINGLE)
    with pytest.raises(ValueError, match="both ra and dec"):
        red_ixpe.qu_slice_pcubes(OBSID_SINGLE, file_list, ra=299.59, dec=None)
    with pytest.raises(ValueError, match="both ra and dec"):
        red_ixpe.qu_slice_pcubes(OBSID_SINGLE, file_list, ra=None, dec=35.2)


# --------------------------------------------------------------------------- #
# Real-data reduction: single-pointing ObsID (01002901) -- no N-way split
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def single_result(tmp_path_factory):
    """One real qu_slice_pcubes call on a scratch copy of 01002901, n_max
    kept small to bound the real xpselect/xpbin call count."""
    tmp_dir = tmp_path_factory.mktemp("qu_slice_single")
    file_list = _copy_real_event_files(str(tmp_dir), OBSID_SINGLE)
    result = red_ixpe.qu_slice_pcubes(OBSID_SINGLE, file_list, n_max=3, n_min=2)
    return result


def test_single_pointing_yields_one_subobs(single_result):
    assert single_result["subobs_keys"] == [OBSID_SINGLE]
    assert set(single_result["per_subobs"]) == {OBSID_SINGLE}


def test_capability_slices_narrower_than_whole_and_partition_exposure(single_result):
    """Pins the tool's core defect-to-avoid: a slice must genuinely reduce a
    NARROWER time window than the whole (sub-)observation, not re-run the
    whole thing N times. Every slice's exposure/counts is strictly less than
    the whole-sub-obs value, and slice exposures sum back to (approximately)
    the whole-sub-obs exposure."""
    info = single_result["per_subobs"][OBSID_SINGLE]
    rows = info["rows"]
    assert len(rows) == info["n_slice"] >= 2

    for row in rows:
        assert row["good_exposure_s"] < info["whole_exposure_s"]
        assert row["counts"] < info["whole_counts"]

    total_exposure = sum(row["good_exposure_s"] for row in rows)
    assert total_exposure == pytest.approx(info["whole_exposure_s"], rel=0.02)


def test_known_input_mean_arrival_time_within_window(single_result):
    """t_mean_met must lie within [t_start_met, t_stop_met) for every slice
    -- proves the mean was computed from that slice's own selected events."""
    rows = single_result["per_subobs"][OBSID_SINGLE]["rows"]
    for row in rows:
        assert row["t_start_met"] <= row["t_mean_met"] < row["t_stop_met"] + 1e-6


def test_csv_matches_returned_rows(single_result):
    info = single_result["per_subobs"][OBSID_SINGLE]
    with open(info["csv_path"], newline="") as f:
        provenance_line = f.readline()
        csv_rows = list(csv.DictReader(f))
    assert provenance_line.startswith("# qu_slice_pcubes provenance: ")
    provenance = json.loads(provenance_line[len("# qu_slice_pcubes provenance: "):])
    assert provenance["n_max"] == 3
    assert list(csv_rows[0].keys()) == red_ixpe._QU_SLICE_FIELDS
    assert len(csv_rows) == len(info["rows"])
    assert float(csv_rows[0]["PD_pct"]) == pytest.approx(info["rows"][0]["PD_pct"])


def test_independent_differential_whole_cube_matches_pcube_timeseg(single_result, tmp_path_factory):
    """Independent differential: qu_slice_pcubes' whole-sub-obs PD/PA (via
    _reduce_window) must reproduce the already-trusted pcube_timeseg tool's
    integrated whole-observation PD/PA (via its own separate _reduce path) on
    the same single-pointing ObsID, same energy band / src region / background
    choice. A bug shared only inside qu_slice (e.g. wrong weights, wrong
    energy binning, wrong Stokes attribute) is invisible to the tool's own
    internal capability/repro checks but is caught here against independent code."""
    info = single_result["per_subobs"][OBSID_SINGLE]
    # 01002901 is a bright source-dominated single pointing: qu_slice's
    # bkg=None auto-decision must resolve to NO subtraction, matching
    # pcube_timeseg's own default bkg=False -- assert that so the two
    # reductions are genuinely comparing the same (unsubtracted) measurement.
    assert info["subtract_background"] is False
    assert single_result["subobs_keys"] == [OBSID_SINGLE]  # no split

    tmp_dir = tmp_path_factory.mktemp("qu_slice_diff_timeseg")
    file_list = _copy_real_event_files(str(tmp_dir), OBSID_SINGLE)
    # A segment_duration_s far larger than the observation span -> exactly one
    # elapsed-time window, so the pipeline runs the whole-obs integrated
    # reduction (plus at most that single window) rather than ~170 segments.
    ts = red_ixpe.pcube_timeseg(
        OBSID_SINGLE, file_list, emin_kev=2.0, emax_kev=8.0, src_rad=1.0,
        bkg=False, segment_duration_s=1.0e9, resume=False,
    )
    integrated = ts["integrated"]
    # Same physical quantity, two independent reduction paths -> agree tightly.
    assert integrated["PD_pct"] == pytest.approx(info["whole_pd_pct"], rel=0.02)
    # PA is an angle in degrees; compare with a mod-180 tolerance so a wrap
    # near +/-90 doesn't spuriously fail. They should agree to well under a deg.
    dpa = (integrated["PA_deg"] - info["whole_pa_deg"] + 90.0) % 180.0 - 90.0
    assert abs(dpa) < 1.0


def test_repro_qu_slice_pcubes_stable_across_calls(tmp_path_factory):
    """Reproducibility: two independent qu_slice_pcubes calls on independent
    scratch copies with identical arguments give the same whole-sub-obs PD
    and the same (single, n_max=1) slice's q/u within tolerance."""
    results = []
    for _ in range(2):
        tmp_dir = tmp_path_factory.mktemp("qu_slice_repro")
        file_list = _copy_real_event_files(str(tmp_dir), OBSID_SINGLE)
        results.append(red_ixpe.qu_slice_pcubes(OBSID_SINGLE, file_list, n_max=1, n_min=1))

    a, b = results
    info_a = a["per_subobs"][OBSID_SINGLE]
    info_b = b["per_subobs"][OBSID_SINGLE]
    assert info_a["n_slice"] == info_b["n_slice"] == 1

    for field in ("whole_pd_pct", "whole_pa_deg"):
        assert info_a[field] == pytest.approx(info_b[field], rel=1e-6)
    for field in ("q", "u", "PD_pct", "PA_deg"):
        assert info_a["rows"][0][field] == pytest.approx(info_b["rows"][0][field], rel=1e-6)


# --------------------------------------------------------------------------- #
# Real-data reduction: literature-merged ObsID (03002599) -- N-way split
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def merged_result(tmp_path_factory):
    """One real qu_slice_pcubes call on a scratch copy of 03002599, n_max
    kept small (2) to bound the real xpselect/xpbin call count across both
    resulting sub-observations."""
    tmp_dir = tmp_path_factory.mktemp("qu_slice_merged")
    file_list = _copy_real_event_files(str(tmp_dir), OBSID_MERGED)
    result = red_ixpe.qu_slice_pcubes(OBSID_MERGED, file_list, n_max=2, n_min=1)
    return result


def test_split_structure_yields_at_least_two_subobs(merged_result):
    assert len(merged_result["subobs_keys"]) >= 2
    assert merged_result["subobs_keys"] == sorted(merged_result["subobs_keys"])
    assert set(merged_result["subobs_keys"]) == set(merged_result["per_subobs"])


def test_no_slice_straddles_the_dominant_central_gti_gap(merged_result):
    """Every slice's [t_start_met, t_stop_met] must lie within exactly one
    GTI-contiguous sub-observation stretch -- no slice spans the dominant
    ~5.4e6s literature-merge gap between the two sub-observations."""
    spans = []
    for key in merged_result["subobs_keys"]:
        info = merged_result["per_subobs"][key]
        rows = info["rows"]
        assert rows, f"sub-obs {key} produced no slices"
        span_lo = min(r["t_start_met"] for r in rows)
        span_hi = max(r["t_stop_met"] for r in rows)
        spans.append((span_lo, span_hi))
        for row in rows:
            assert span_lo <= row["t_start_met"] < row["t_stop_met"] <= span_hi

    spans.sort()
    for (_, hi_prev), (lo_next, _) in zip(spans, spans[1:]):
        assert lo_next - hi_prev > 1e6  # the real dominant merge gap, not an ordinary orbital gap


def test_mean_arrival_time_differs_from_window_midpoint_for_at_least_one_slice(merged_result):
    """Proves mean-arrival timing was actually used, not the window midpoint
    -- expected on a real, gap-riddled observation like 03002599."""
    any_differs = False
    for key in merged_result["subobs_keys"]:
        for row in merged_result["per_subobs"][key]["rows"]:
            assert row["t_start_met"] <= row["t_mean_met"] < row["t_stop_met"] + 1e-6
            midpoint = 0.5 * (row["t_start_met"] + row["t_stop_met"])
            if abs(row["t_mean_met"] - midpoint) > 1.0:
                any_differs = True
    assert any_differs


# --------------------------------------------------------------------------- #
# Skill wrapper: file discovery (rej/plain fallback), catalog recording,
# bounded return contract
# --------------------------------------------------------------------------- #

def test_skill_wrapper_falls_back_to_plain_evt2_and_records_catalog(tmp_path):
    """03002599 carries no '_rej' variant on disk -- the skill wrapper's glob
    must fall back to plain evt2_v*.fits (mirroring get_all_file_list's own
    03002599 special case). Confirms end-to-end: file discovery, reduction,
    CSV(s) on disk, level-3 catalog recording (one product per sub-obs, all
    against the real 9-digit obsid), a bounded return (no 'rows' key), and
    -- since this call's sub-obs count is well under the default
    max_return_subobs -- the split-scale-guard fix's UN-TRUNCATED bounded-
    return contract: 'per_subobs_truncated'/'n_subobs_total'/
    'subobs_listing_path' are always present, the listing CSV is always
    written and lists every sub-obs, and stdout never enumerates paths."""
    obs_db_path = str(tmp_path)
    _copy_real_event_files(obs_db_path, OBSID_MERGED)

    result = ixpe_session.ixpe_qu_slice_pcubes(obs_db_path, OBSID_MERGED, n_max=2, n_min=1)

    assert result["exit_code"] == 0, result.get("stderr")
    assert len(result["subobs_keys"]) >= 2
    for key in result["subobs_keys"]:
        info = result["per_subobs"][key]
        assert "rows" not in info
        assert os.path.exists(info["csv_path"])

    entry = neutral.read(obs_db_path)["entries"][OBSID_MERGED]
    assert entry["max_level"] == 3
    products = [p for p in entry["products"] if p["kind"] == "qu_slices"]
    assert len(products) == len(result["subobs_keys"])
    recorded_keys = {p["params"]["subobs_key"] for p in products}
    assert recorded_keys == set(result["subobs_keys"])
    for p in products:
        assert p["produced_by"] == "xspec_utils.red_ixpe.qu_slice_pcubes"

    # Un-truncated branch: the three bounded-return keys are always present.
    assert result["per_subobs_truncated"] is False
    n_total = result["n_subobs_total"]
    assert n_total == len(result["subobs_keys"]) == len(result["per_subobs"])
    listing_path = result["subobs_listing_path"]
    assert os.path.exists(listing_path)
    with open(listing_path, newline="") as f:
        listing_rows = list(csv.DictReader(f))
    assert len(listing_rows) == n_total
    assert {r["subobs_key"] for r in listing_rows} == set(result["subobs_keys"])
    assert ".csv" not in result["stdout"] and len(result["stdout"]) < 400


def test_skill_wrapper_missing_event_files_fails_loudly_no_record(tmp_path):
    result = ixpe_session.ixpe_qu_slice_pcubes(str(tmp_path), "99999999")
    assert result["exit_code"] == 1
    assert "stderr" in result
    assert neutral.read(str(tmp_path))["entries"] == {}


def test_failure_in_later_subobs_leaves_no_csv(tmp_path, monkeypatch):
    """A later sub-obs raising must leave NO qu_slices CSV on disk and record
    nothing -- the 'nothing written on failure' contract, even though the
    earlier sub-obs's reduction succeeded. Pins the deferred-write fix: CSVs
    are written only after every sub-observation in the loop has completed
    without raising, so sub-obs a's already-successful reduction must not
    leave its CSV behind when sub-obs b subsequently fails."""
    _copy_real_event_files(str(tmp_path), OBSID_MERGED)
    # The skill wrapper calls qu_slice_pcubes through session.py's own
    # `import red_ixpe` (the top-level module alias), which is a DIFFERENT
    # module object than `xspec_utils.lib.red_ixpe`. Patch the object the
    # wrapper actually executes against (ixpe_session.red_ixpe), else the
    # patch is a no-op and the tool runs to completion.
    target = ixpe_session.red_ixpe
    real_reduce = target._reduce_window

    def boom(obs_file_list, suffix, *args, **kwargs):
        if suffix.startswith(OBSID_MERGED + "b"):
            raise RuntimeError("boom")
        return real_reduce(obs_file_list, suffix, *args, **kwargs)

    monkeypatch.setattr(target, "_reduce_window", boom)
    result = ixpe_session.ixpe_qu_slice_pcubes(str(tmp_path), OBSID_MERGED, n_max=2, n_min=1)
    assert result["exit_code"] == 1
    assert result.get("stderr")
    assert glob.glob(os.path.join(str(tmp_path), OBSID_MERGED, "qu_slices_*.csv")) == []
    assert neutral.read(str(tmp_path))["entries"] == {}
