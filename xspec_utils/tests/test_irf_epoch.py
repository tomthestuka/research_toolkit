"""Unit tests for the nearest-epoch IXPE IRF resolver (lib/irf_epoch.py).

Pure date/selection math -- no ixpeobssim / science stack (the CALDB scan is
exercised against a temp directory of fake MRF filenames), so these run in any
environment with the toolkit on the path.
"""

from __future__ import annotations

import pytest

from xspec_utils.lib import irf_epoch

# The installed epochs on this machine, as of writing -- used only as a realistic
# fixture; the resolver takes the epoch list by argument so nothing here touches
# the actual CALDB.
EPOCHS = [
    "20211209", "20220702", "20230101", "20230702",
    "20240101", "20240701", "20250101", "20250701",
]


def test_epoch_to_mjd_matches_known_anchor():
    # 2000-01-01 is MJD 51544; 1858-11-17 is MJD 0.
    assert irf_epoch._epoch_to_mjd("20000101") == 51544.0
    assert irf_epoch._epoch_to_mjd("18581117") == 0.0


def test_resolve_picks_nearest_even_when_later_than_obs():
    # Cyg X-1 May 2022 obs (MJD ~59714, 2022-05-15) sits between 20211209
    # (MJD 59557) and 20220702 (MJD 59762). Nearest by |distance| is 20220702
    # (~47 d away) over 20211209 (~158 d), though it begins after the obs.
    assert irf_epoch.resolve_irf_date(59714.0, epochs=EPOCHS) == "20220702"


def test_resolve_picks_earlier_when_closer():
    # An obs a few days after 20220702 stays on that epoch (nearest by distance).
    assert irf_epoch.resolve_irf_date(59765.0, epochs=EPOCHS) == "20220702"


def test_resolve_obs_before_earliest_returns_earliest():
    # Nothing to fall back to is needed -- nearest simply returns the closest,
    # which for a pre-launch date is the earliest installed epoch.
    assert irf_epoch.resolve_irf_date(58000.0, epochs=EPOCHS) == "20211209"


def test_resolve_obs_after_latest_returns_latest():
    assert irf_epoch.resolve_irf_date(70000.0, epochs=EPOCHS) == "20250701"


def test_resolve_tie_breaks_to_earlier():
    # Obs exactly midway between two adjacent epochs -> the earlier one wins
    # (its calibration window has already opened).
    a, b = "20230101", "20230701"
    mid = 0.5 * (irf_epoch._epoch_to_mjd(a) + irf_epoch._epoch_to_mjd(b))
    assert irf_epoch.resolve_irf_date(mid, epochs=[a, b]) == a


def test_resolve_single_epoch():
    assert irf_epoch.resolve_irf_date(59714.0, epochs=["20220702"]) == "20220702"


def test_resolve_empty_raises():
    with pytest.raises(RuntimeError, match="No installed IXPE IRF epochs"):
        irf_epoch.resolve_irf_date(59714.0, epochs=[])


def test_installed_irf_epochs_parses_and_dedups(tmp_path):
    # Fake CALDB MRF dir: several DUs and calibration flavours repeat each epoch;
    # the scan must return one sorted entry per epoch and ignore non-matching files.
    names = [
        "ixpe_d1_obssim20211209_alpha075_v013.mrf",
        "ixpe_d2_obssim20211209_alpha075_v013.mrf",
        "ixpe_d1_obssim20211209_gray_v013.mrf",
        "ixpe_d3_obssim20220702_alpha075_v013.mrf",
        "ixpe_d1_obssim20230101_alpha075simple_v013.mrf",
        "readme.txt",  # ignored
        "ixpe_d1_obssim_v013.mrf",  # no date token -> ignored
    ]
    for n in names:
        (tmp_path / n).write_text("")
    assert irf_epoch.installed_irf_epochs(str(tmp_path)) == [
        "20211209", "20220702", "20230101",
    ]


def test_installed_irf_epochs_empty_dir(tmp_path):
    assert irf_epoch.installed_irf_epochs(str(tmp_path)) == []
