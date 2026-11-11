"""High-level acquire entry point — download an ObsID and store it as a product.

This is the abstracted, human-callable API (and the function the LLM tool wraps). An archive download
is a product like any other (``storelib.produce``): bucket = the mission (HEASARC's name for it),
kind = ``raw``, subject = the ObsID, settings = what identifies the archive's copy (the master-table
row's processing date, software and CALDB versions, and whether the observation is public yet). The
same recipe is downloaded once: a repeat acquire is answered from the store. A reprocessed
observation, or one that has since become public, is a new recipe and a fresh download.

**Encrypted (proprietary) observations.** HEASARC serves data still inside its proprietary period
GPG-encrypted. Such a download is complete and correct — it is simply not readable without the key —
so it is stored rather than treated as a failed download, and the result says so (``locked``). Asking
again while the observation is still proprietary is the same recipe, so the stored copy is returned
without a transfer; once it is public the recipe changes and the readable copy is downloaded.
"""

from __future__ import annotations

import posixpath

import storelib

from . import heasarc

__all__ = ["acquire", "TOOL", "ALGORITHM_VERSION", "KIND"]

#: The manifest key of the tool this function backs: the recipe names it.
TOOL = "research_toolkit.heasarc_utils.acquire"
#: Bumped whenever what a download stores changes (the manifest's ``algorithm_version``).
ALGORITHM_VERSION = "1"
#: The product kind of an archive download.
KIND = "raw"
#: How long a download may hold its claim (seconds): above the tool's call timeout ceiling (6 h).
LEASE_S = 7 * 3600


def acquire(obsid: str, mission: str, *, host: str | None = None, chatter: bool = False,
            timeout_s: float = 1800.0) -> dict:
    """Acquire archive data for ``obsid`` of ``mission`` into the store and return its record.

    Resolves the ObsID authoritatively via the mission's HEASARC master table (one row, which also
    states the archive's processing and public date), then stores the observation as the product
    ``raw/<obsid>/<short hash>/`` in the mission's bucket of the store the process's keys reach (``storelib.store_name``):
    downloaded (astroquery only — no fallback) into a temporary folder in the working folder,
    gunzipped, checked to hold the representative level-2 event, and uploaded file by file with the
    ObsID folder's own layout (``event_l2/...``). When the store already holds this recipe, nothing is
    downloaded.

    Parameters
    ----------
    obsid : str
        Observation ID (mission-native format, e.g. ``"01002901"`` for IXPE).
    mission : str
        ``"ixpe"`` | ``"nicer"`` | ``"nustar"`` (or any registered mission).
    host : str, optional
        astroquery download host (``"heasarc"`` | ``"sciserver"`` | ``"aws"``). Default lets
        astroquery choose based on the environment.
    timeout_s : float, optional
        HTTP transfer timeout in seconds for the archive download. Default 1800 s (30 min). The
        astroquery factory default (30 s) is too short for large IXPE observations (> 1 GB).

    Returns
    -------
    dict
        ``{"obsid", "mission", "record", "level2_evt", "downloaded", "locked", "archive", "files"}``:
        ``record`` the product's record (pass it to tools that read the observation), ``level2_evt``
        the record of its level-2 event, ``downloaded`` false when the store already held it,
        ``locked`` true when that event is GPG-encrypted, ``archive`` the master-table state the
        recipe holds, ``files`` how many files the product holds. A locked observation adds
        ``days_until_public`` and a ``note``.

    Raises
    ------
    heasarc.HeasarcError
        Unknown mission, ObsID not found, an archive folder other than the mission's, or a download
        that does not yield the expected level-2 file. Nothing is stored on failure (the claim is
        given up with the reason).
    storelib.StoreError
        The store refused, or this process holds no store keys; ``storelib.ProductBusy`` while another
        unit downloads the same recipe.
    """
    mission = mission.lower()
    profile = heasarc.profile(mission)
    tab, links = heasarc.resolve(profile.catalog, obsid)
    folder, _tree = heasarc.archive_path(links, obsid)
    if folder.lower() != mission:
        raise heasarc.HeasarcError(f"HEASARC keeps {mission} ObsID {obsid} under FTP/{folder}/; register the "
                                   f"mission under its archive name {folder!r}")
    state = heasarc.archive_state(tab)
    if state["processing_date"] is None:
        raise heasarc.HeasarcError(
            f"{profile.catalog} states no processing_date for ObsID {obsid!r}: the archive's processing state "
            "identifies a download, so there is nothing to store until the archive has processed it")
    public = state["public_date"] is None or state["public_date"] <= heasarc.mjd_now()
    settings = {"mission": mission, "processing_date": state["processing_date"],
                "software_version": state["software_version"], "caldb_version": state["caldb_version"],
                "public": public}

    def make(work):
        paths = heasarc.download_obsid(links, obsid, str(work), host=host, timeout_s=timeout_s)
        level2 = profile.level2(paths, obsid)  # raises when the download holds no level-2 event
        if not heasarc.is_locked(level2):
            # The recipe names the row's processing; the bytes must be that processing.
            lag = heasarc.state_mismatch(state, str(work / level2))
            if lag:
                raise heasarc.HeasarcError(
                    f"{obsid}: the archive's level-2 file {level2} holds {lag}: the archive still serves an "
                    "older processing than its catalog names (it lags its catalog); nothing is stored. "
                    "Acquire it again later")
        if public and heasarc.is_locked(level2):
            # The archive still serves the encrypted copy after the public date (a lag): storing it
            # under a public recipe would answer every later acquire of this state with it.
            raise heasarc.HeasarcError(
                f"{obsid} is public since MJD {state['public_date']:.0f}, but the archive still serves it "
                "GPG-encrypted; nothing is stored. Acquire it again later")
        if chatter:
            print(f"[heasarc_utils] {mission} {obsid}: {len(paths)} files downloaded (timeout_s={timeout_s})")
        return paths

    out = storelib.produce(bucket=mission, kind=KIND, subject=str(obsid), tool=TOOL,
                           algorithm_version=ALGORITHM_VERSION, settings=settings, inputs=[], make=make,
                           lease_s=LEASE_S)
    folder_key = storelib.product_folder(out["record"])
    files = {posixpath.relpath(o["key"], folder_key): o for o in out["objects"]
             if o["key"] != out["record"]["key"]}
    level2 = files[profile.level2(sorted(files), obsid)]
    locked = heasarc.is_locked(level2["key"])
    summary = {
        "obsid": str(obsid),
        "mission": mission,
        "record": out["record"],
        "level2_evt": {k: level2[k] for k in ("store", "bucket", "key", "version_id")},
        "downloaded": out["status"] == "stored",
        "locked": locked,
        "archive": {k: v for k, v in settings.items() if k != "mission"},
        "files": len(files),
    }
    if locked:
        days = max(0.0, state["public_date"] - heasarc.mjd_now()) if state["public_date"] is not None else None
        summary["days_until_public"] = days
        summary["note"] = (
            f"{obsid} was served GPG-encrypted (proprietary period"
            f"{'' if days is None else f', public in {days:.0f} days'}): it is stored but not directly usable. "
            "Acquire again once the observation is public; until then the store answers with this copy.")
    if chatter:
        print(f"[heasarc_utils] {mission} {obsid}: {out['status']} as {out['record']['key']}")
    return summary
