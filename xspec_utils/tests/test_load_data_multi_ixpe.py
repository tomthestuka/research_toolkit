"""Tests for xspec_utils' load_data multi-ObsID IXPE capability
(lib/xspec_help.py's IXPE branch).

Reproducing Krawczynski+2022's Cyg X-1 single-session tied-polconst
cross-check needs one tbabs*polconst*powerlaw fit tying PD/PA across every
DU of BOTH the May (01002901) and June (01250101) IXPE ObsIDs -- but
load_data's IXPE branch raised NotImplementedError for any epoch['IXPE']
list longer than one ObsID. This generalizes it to an inner
``for obs in epoch['IXPE']`` loop, mirroring the NuSTAR/NICER branches
exactly: each ObsID independently resolves its own source centroid,
region, per-spectrum RESPFILE/ANCRFILE/BACKSCAL, and IRF calibration
epoch; data groups are assigned consecutively and gaplessly (3 per
ObsID); and each (ObsID, DU)'s Q/U gain slope ties only to that same
(ObsID, DU)'s own I gain slope -- never across ObsIDs or DUs.

Runs real ixpeobssim fresh reductions (ixpe_changes=True) against the
already-local Level-2 event files for both validation ObsIDs, copied into
isolated tmp directories per fixture so this never mutates the shared
analyst-env data directory.
"""

from __future__ import annotations

import os
import re
import shutil

import pytest
from xspec import AllData, AllModels

from xspec_utils.lib import red_ixpe, irf_epoch, xspec_help
from xspec_utils.skills import TOOL_MANIFEST

REAL_IXPE_DIR = "/data/fzhou/XSPEC/Cyg_X-1/IXPE"
OBSID_A = "01002901"  # 2022 May
OBSID_B = "01250101"  # 2022 June
OBSIDS = [OBSID_A, OBSID_B]
VERSION = "v01"
SRC_RAD, BKG_INNER_RAD, BKG_OUTER_RAD = 1.333, 2.5, 5.167

_EPOCH_RE = re.compile(r"obssim(\d{8})")
_LINK_RE = re.compile(r"^\s*=\s*p(\d+)\s*$")


def _copy_event_files(dest_dir: str, obsid: str) -> None:
    src_evt_dir = os.path.join(REAL_IXPE_DIR, obsid, "event_l2")
    dst_evt_dir = os.path.join(dest_dir, obsid, "event_l2")
    os.makedirs(dst_evt_dir, exist_ok=True)
    for det in (1, 2, 3):
        src = os.path.join(src_evt_dir, f"ixpe{obsid}_det{det}_evt2_{VERSION}_rej.fits")
        if not os.path.exists(src):
            pytest.skip(f"real fixture ObsID {obsid} not local at {src}")
        shutil.copyfile(src, os.path.join(dst_evt_dir, os.path.basename(src)))


def _linked_i_spec_index(link_str: str) -> int:
    """Decode a gain.slope.link string (e.g. '= p7') to the 1-based XSPEC
    spectrum index it ties to. Gain carries exactly 2 sub-parameters
    (slope, offset) per response, assigned sequentially in load order with
    no Model() ever defined for these spectra, so response n's slope is the
    global parameter number 2n-1; inverting recovers n."""
    m = _LINK_RE.match(link_str)
    assert m, f"unexpected gain.slope.link string: {link_str!r}"
    p = int(m.group(1))
    assert p % 2 == 1, f"expected an odd global gain-slope parameter number, got p{p}"
    return (p + 1) // 2


def _spec_owner_map(indices: dict) -> dict[int, tuple[str, str]]:
    """spec_index -> (obsid, du) for every IXPE spectrum in *indices*."""
    owner = {}
    for obsid, dus in indices["IXPE"].items():
        for du, triplet in dus.items():
            for si in triplet:
                owner[si] = (obsid, du)
    return owner


def _rate_and_epoch(i_idx: int) -> tuple[float, str]:
    """Net Stokes-I rate and the obssimYYYYMMDD IRF epoch actually resolved
    for spectrum *i_idx*, read from its own (real, on-disk) RESPFILE path --
    never a hand-carried value."""
    net, _var, _total, _model = AllData(i_idx).rate
    rmf_path = AllData(i_idx).response.rmf
    m = _EPOCH_RE.search(rmf_path)
    assert m, f"no obssimYYYYMMDD epoch token found in RESPFILE path {rmf_path!r}"
    return net, m.group(1)


def _expected_irf_epoch(obsid: str, file_list: list[str]) -> str:
    """The IRF epoch load_data's IXPE branch *should* auto-resolve for
    *obsid*, computed independently via the same trusted primitives
    red_ixpe._resolve_irf_date itself calls (get_mjd + irf_epoch.resolve_irf_date)
    -- not by re-invoking red_ixpe._resolve_irf_date itself."""
    _start, _stop, mjd_mid = red_ixpe.get_mjd(obsid, file_list)
    return irf_epoch.resolve_irf_date(mjd_mid)


def _extract_multi_load(ixpe_dir: str) -> dict:
    """Run one real fresh multi-ObsID load_data call for OBSIDS and
    materialize every plain-value fact the tests below need, since a
    later load_data call (single-ObsID comparisons, the repro call) clears
    and replaces the shared global XSPEC session."""
    AllData.clear()
    AllModels.clear()
    indices = xspec_help.load_data(
        epoch={"IXPE": OBSIDS},
        ixpe_path=ixpe_dir,
        ixpe_changes=True,
        ixpe_src_rad=SRC_RAD,
        ixpe_bkg_inner_rad=BKG_INNER_RAD,
        ixpe_bkg_outer_rad=BKG_OUTER_RAD,
    )

    owner = _spec_owner_map(indices)
    n_groups = AllData.nGroups

    dg_by_obsid_du = {}
    rates_epochs = {}
    gain_links = {}
    for obsid, dus in indices["IXPE"].items():
        rates_epochs[obsid] = {}
        for du, (i_idx, q_idx, u_idx) in dus.items():
            dg_by_obsid_du[(obsid, du)] = AllData(i_idx).dataGroup
            net, epoch = _rate_and_epoch(i_idx)
            rates_epochs[obsid][du] = {"net": net, "epoch": epoch}

            i_link = AllData(i_idx).response.gain.slope.link
            q_link = AllData(q_idx).response.gain.slope.link
            u_link = AllData(u_idx).response.gain.slope.link
            gain_links[(obsid, du)] = {
                "i_idx": i_idx, "i_link": i_link,
                "q_target": _linked_i_spec_index(q_link) if q_link else None,
                "u_target": _linked_i_spec_index(u_link) if u_link else None,
            }

    return {
        "indices": indices,
        "owner": owner,
        "n_groups": n_groups,
        "dg_by_obsid_du": dg_by_obsid_du,
        "rates_epochs": rates_epochs,
        "gain_links": gain_links,
    }


@pytest.fixture(scope="module")
def multi_load(tmp_path_factory):
    """One real fresh multi-ObsID (May + June) reduction + load_data call,
    shared read-only across the structure/gain-linking/differential tests
    below to avoid repeating the expensive xpselect/xpbin pipeline."""
    ixpe_dir = str(tmp_path_factory.mktemp("multi_ixpe"))
    for obsid in OBSIDS:
        _copy_event_files(ixpe_dir, obsid)
    return _extract_multi_load(ixpe_dir), ixpe_dir


@pytest.fixture(scope="module")
def single_loads(tmp_path_factory):
    """Independent single-ObsID load_data calls (real fresh reductions), one
    per validation ObsID, each in its own isolated tmp directory -- the
    baseline the multi-ObsID load must reproduce per-ObsID."""
    results = {}
    for obsid in OBSIDS:
        ixpe_dir = str(tmp_path_factory.mktemp(f"single_ixpe_{obsid}"))
        _copy_event_files(ixpe_dir, obsid)
        AllData.clear()
        AllModels.clear()
        indices = xspec_help.load_data(
            epoch={"IXPE": [obsid]},
            ixpe_path=ixpe_dir,
            ixpe_changes=True,
            ixpe_src_rad=SRC_RAD,
            ixpe_bkg_inner_rad=BKG_INNER_RAD,
            ixpe_bkg_outer_rad=BKG_OUTER_RAD,
        )
        owner = _spec_owner_map(indices)
        dg_by_du, rates_epochs, gain_links = {}, {}, {}
        for du, (i_idx, q_idx, u_idx) in indices["IXPE"][obsid].items():
            dg_by_du[du] = AllData(i_idx).dataGroup
            net, epoch = _rate_and_epoch(i_idx)
            rates_epochs[du] = {"net": net, "epoch": epoch}
            q_link = AllData(q_idx).response.gain.slope.link
            u_link = AllData(u_idx).response.gain.slope.link
            gain_links[du] = {
                "i_idx": i_idx,
                "q_target": _linked_i_spec_index(q_link),
                "u_target": _linked_i_spec_index(u_link),
            }
        results[obsid] = {
            "indices": indices, "owner": owner, "n_groups": AllData.nGroups,
            "dg_by_du": dg_by_du, "rates_epochs": rates_epochs,
            "gain_links": gain_links,
        }
    return results


@pytest.fixture(scope="module")
def second_multi_load(tmp_path_factory):
    """A second, fully independent multi-ObsID load_data call (own tmp
    copy, own reduction) -- the reproducibility check's second call."""
    ixpe_dir = str(tmp_path_factory.mktemp("multi_ixpe_repro"))
    for obsid in OBSIDS:
        _copy_event_files(ixpe_dir, obsid)
    return _extract_multi_load(ixpe_dir)


# --------------------------------------------------------------------------- #
# Structure: parameterized on N, not the literal 6
# --------------------------------------------------------------------------- #

def test_structure_dg_gapless_and_indices_present(multi_load):
    ctx, _ixpe_dir = multi_load
    n = len(OBSIDS)

    assert ctx["n_groups"] == 3 * n
    assert set(ctx["indices"]["IXPE"].keys()) == set(OBSIDS)

    for obsid in OBSIDS:
        assert set(ctx["indices"]["IXPE"][obsid].keys()) == {"DU1", "DU2", "DU3"}
        for triplet in ctx["indices"]["IXPE"][obsid].values():
            assert len(triplet) == 3

    # Every spectrum index 1..9*n is used exactly once, and every data
    # group 1..3*n is used by exactly 3 spectra (I, Q, U of one DU) --
    # consecutive and gapless across both ObsIDs.
    all_specs = sorted(ctx["owner"].keys())
    assert all_specs == list(range(1, 9 * n + 1))

    dgs = sorted(ctx["dg_by_obsid_du"].values())
    assert dgs == list(range(1, 3 * n + 1))

    # Data-group order follows epoch key order: obsA's 3 DUs occupy the
    # first 3n/2... i.e. the lower-numbered data groups than obsB's.
    obsA_dgs = [ctx["dg_by_obsid_du"][(OBSID_A, du)] for du in ("DU1", "DU2", "DU3")]
    obsB_dgs = [ctx["dg_by_obsid_du"][(OBSID_B, du)] for du in ("DU1", "DU2", "DU3")]
    assert max(obsA_dgs) < min(obsB_dgs)


# --------------------------------------------------------------------------- #
# Capability: gain linking proves the fix, not just well-shaped metadata
# --------------------------------------------------------------------------- #

def test_capability_gain_linking_per_obsid_du_no_cross_mislink(multi_load):
    """Every (ObsID, DU)'s Q and U gain-slope link resolves to that same
    (ObsID, DU)'s own I gain-slope parameter, and never to a spectrum
    belonging to a different ObsID (or a different DU) -- the actual defect
    a naive generalization (e.g. always linking back to spectrum 1) would
    produce."""
    ctx, _ixpe_dir = multi_load
    owner = ctx["owner"]

    for obsid in OBSIDS:
        for du in ("DU1", "DU2", "DU3"):
            link = ctx["gain_links"][(obsid, du)]
            assert link["i_link"] == "", (
                f"{obsid}/{du}: I gain slope must be free (unlinked), got {link['i_link']!r}")
            for sto in ("q_target", "u_target"):
                target_idx = link[sto]
                assert target_idx == link["i_idx"], (
                    f"{obsid}/{du}: {sto} does not resolve to this DU's own I "
                    f"spectrum {link['i_idx']} (resolved to {target_idx})")
                assert owner[target_idx] == (obsid, du), (
                    f"{obsid}/{du}: {sto} gain link resolves to a spectrum "
                    f"owned by {owner[target_idx]} -- cross-ObsID/DU mislink")


# --------------------------------------------------------------------------- #
# Differential: multi-ObsID load reproduces a separate single-ObsID load,
# per ObsID -- proving independent response/BACKSCAL/IRF-epoch resolution
# --------------------------------------------------------------------------- #

def test_differential_per_obsid_matches_independent_single_load(multi_load, single_loads):
    ctx, _ixpe_dir = multi_load
    for obsid in OBSIDS:
        multi_re = ctx["rates_epochs"][obsid]
        single_re = single_loads[obsid]["rates_epochs"]
        for du in ("DU1", "DU2", "DU3"):
            assert multi_re[du]["net"] == pytest.approx(single_re[du]["net"], rel=1e-6), (
                f"{obsid}/{du}: multi-ObsID net Stokes-I rate does not match an "
                f"independent single-ObsID load_data call")
            assert multi_re[du]["epoch"] == single_re[du]["epoch"], (
                f"{obsid}/{du}: multi-ObsID resolved IRF epoch "
                f"{multi_re[du]['epoch']} != independent single-ObsID epoch "
                f"{single_re[du]['epoch']}")


def test_differential_irf_epoch_matches_independent_resolution(multi_load):
    """Each ObsID's resolved IRF epoch (read off its own RESPFILE) must
    match an independent computation via get_mjd + irf_epoch.resolve_irf_date
    -- the trusted primitives red_ixpe._resolve_irf_date itself calls, not
    load_data's own code path -- confirming per-ObsID (not shared/cached)
    resolution. May (01002901) and June (01250101) 2022 happen to both
    resolve nearest the same installed 20220702 epoch, but each is computed
    from that ObsID's own midpoint MJD, independently."""
    ctx, ixpe_dir = multi_load
    file_list = red_ixpe.get_all_file_list(ixpe_dir)
    for obsid in OBSIDS:
        expected = _expected_irf_epoch(obsid, file_list)
        for du in ("DU1", "DU2", "DU3"):
            assert ctx["rates_epochs"][obsid][du]["epoch"] == expected, (
                f"{obsid}/{du}: resolved epoch does not match independent "
                f"get_mjd+resolve_irf_date computation")


# --------------------------------------------------------------------------- #
# Reproducibility: two independent multi-ObsID load_data calls agree
# --------------------------------------------------------------------------- #

def test_repro_stable_across_independent_multi_obsid_calls(multi_load, second_multi_load):
    first, _ixpe_dir = multi_load
    second = second_multi_load
    for obsid in OBSIDS:
        for du in ("DU1", "DU2", "DU3"):
            a = first["rates_epochs"][obsid][du]
            b = second["rates_epochs"][obsid][du]
            assert a["net"] == pytest.approx(b["net"], rel=1e-6), (
                f"{obsid}/{du}: net Stokes-I rate not stable across two "
                f"independent multi-ObsID load_data calls")
            assert a["epoch"] == b["epoch"], (
                f"{obsid}/{du}: resolved IRF epoch not stable across two "
                f"independent multi-ObsID load_data calls")


# --------------------------------------------------------------------------- #
# Regression: N=1 must behave exactly as before this change
# --------------------------------------------------------------------------- #

def test_regression_single_obsid_unchanged_schema_and_gain_link(single_loads):
    for obsid in OBSIDS:
        result = single_loads[obsid]
        assert result["n_groups"] == 3
        assert set(result["indices"]["IXPE"].keys()) == {obsid}
        assert set(result["indices"]["IXPE"][obsid].keys()) == {"DU1", "DU2", "DU3"}

        dgs = sorted(result["dg_by_du"].values())
        assert dgs == [1, 2, 3]

        for du in ("DU1", "DU2", "DU3"):
            i_idx = result["gain_links"][du]["i_idx"]
            assert result["gain_links"][du]["q_target"] == i_idx
            assert result["gain_links"][du]["u_target"] == i_idx


# --------------------------------------------------------------------------- #
# Manifest wiring (AI-facing half)
# --------------------------------------------------------------------------- #

def test_tool_manifest_documents_multi_obsid_ixpe():
    entry = TOOL_MANIFEST["research_toolkit.xspec_utils.load_data"]
    assert "any number of ObsIDs" in entry["description"]
    assert "independently resolves" in entry["description"]
