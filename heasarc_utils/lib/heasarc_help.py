"""High-level acquire entry point — download an ObsID and record it.

This is the abstracted, human-callable API (and the function the LLM tool wraps).
It is the first producer wired into ``obsdb_catalog``: after the archive data is
on disk it records the level-2 event through the single trusted writer
(``obsdb_catalog.neutral.new``) — the only path that may mutate ``catalog.json``.

**Encrypted (proprietary) observations.** HEASARC serves data still inside its
proprietary period GPG-encrypted. Such a download is complete and correct — it is
simply not readable without the key — so it is *recorded* rather than treated as a
failed download, and re-running acquire consults the archive's ``public_date``
instead of transferring the payload again. The lifecycle is: download once, then on
a repeat request ask whether the observation has become public, and download again
only once it has. Before this, an encrypted download recorded nothing, so the
catalog could never learn the ObsID existed and every later request re-downloaded
gigabytes that could never be read.
"""

from __future__ import annotations

import os

from obsdb_catalog import neutral

from . import heasarc

__all__ = ["acquire"]


def _locked_copy_on_disk(obs_db_path: str, obsid: str) -> str | None:
    """Relative path of an already-recorded encrypted product still on disk, else None.

    Checked live rather than trusted from the catalog alone: a human may decrypt the
    archive out of band, at which point the recorded ``.gpg`` path stops existing and
    this correctly reports nothing to hold back for.
    """
    try:
        entries = neutral.read(obs_db_path).get("entries") or {}
    except neutral.CatalogError:
        return None  # a broken catalog is repair's problem, never a reason to skip work
    for product in (entries.get(str(obsid)) or {}).get("products") or []:
        path = (product or {}).get("path") or ""
        if neutral.is_locked(path) and os.path.exists(os.path.join(obs_db_path, path)):
            return path
    return None


def acquire(obsid: str, mission: str, obs_db_path: str,
            *, host: str | None = None, chatter: bool = False,
            timeout_s: float = 1800.0) -> dict:
    """Acquire archive data for ``obsid`` into ``obs_db_path`` and record it.

    Resolves the ObsID authoritatively via the mission's HEASARC master table,
    downloads it (astroquery only — no fallback), confirms the representative
    level-2 event is on disk, and records it via ``obsdb_catalog.neutral.new``
    (level 2). **Re-downloads on every call**, so the caller decides when to
    acquire; a deliberate re-acquire regresses ``max_level`` to 2, flagging stale
    downstream products for re-reduction.

    **One exception, and only one:** when this ObsID was already acquired and its
    data is on disk **GPG-encrypted** (HEASARC serves an observation inside its
    proprietary period that way), the archive is asked for its ``public_date``
    first. While the observation is still proprietary the download is skipped and
    the encrypted state reported, because re-transferring gigabytes would return
    the same unreadable bytes. Once the observation is public, acquire downloads
    normally and the readable copy replaces the encrypted record.

    Parameters
    ----------
    obsid : str
        Observation ID (mission-native format, e.g. ``"01002901"`` for IXPE).
    mission : str
        ``"ixpe"`` | ``"nicer"`` | ``"nustar"`` (or any registered mission).
    obs_db_path : str
        The mission's obs-db directory; data lands at ``<obs_db_path>/<obsid>/``.
    host : str, optional
        astroquery download host (``"heasarc"`` | ``"sciserver"`` | ``"aws"``).
        Default lets astroquery choose based on the environment.
    timeout_s : float, optional
        HTTP transfer timeout in seconds for the archive download. Default
        1800 s (30 min). The astroquery factory default (30 s) is too short
        for large IXPE observations (> 1 GB). Increase further if transfers
        are still being cut off on very slow connections.

    Returns
    -------
    dict
        ``{"obsid", "mission", "level2_evt": <relative path>, "catalog_entry",
        "downloaded": bool, "locked": bool}``. ``locked`` is true when the
        level-2 product on disk is GPG-encrypted and therefore not directly
        usable; a skipped re-download additionally carries ``"public_date_mjd"``,
        ``"days_until_public"`` and a human-readable ``"note"``.

    Raises
    ------
    heasarc.HeasarcError
        Unknown mission, ObsID not found, or the download did not yield the
        expected level-2 file. Nothing is recorded on failure.
    """
    mission = mission.lower()

    locked_rel = _locked_copy_on_disk(obs_db_path, obsid)
    if locked_rel is not None:
        profile = heasarc.MISSIONS.get(mission)
        status = heasarc.availability(profile.catalog, obsid) if profile else {}
        if status.get("is_public") is False:
            days = status.get("days_until_public")
            note = (
                f"{obsid} is already on disk GPG-encrypted and is still inside its "
                f"proprietary period (public in {days:.0f} days, MJD "
                f"{status['public_date_mjd']:.0f}), so the archive would serve the "
                f"same encrypted bytes — the download was skipped. Decrypt the local "
                f"copy out of band, or acquire again once the observation is public."
            )
            if chatter:
                print(f"[heasarc_utils] {mission} {obsid}: {note}")
            return {
                "obsid": obsid,
                "mission": mission,
                "level2_evt": locked_rel,
                "catalog_entry": neutral.read(obs_db_path)["entries"].get(str(obsid)),
                "downloaded": False,
                "locked": True,
                "public_date_mjd": status.get("public_date_mjd"),
                "days_until_public": days,
                "note": note,
            }

    product = heasarc.acquire_obsid(mission, obsid, obs_db_path,
                                    host=host, timeout_s=timeout_s)
    locked = neutral.is_locked(product["path"])
    if chatter:
        print(f"[heasarc_utils] {mission} {obsid}: recording {product['path']}"
              f"{' (GPG-encrypted)' if locked else ''} (timeout_s={timeout_s})")
    entry = neutral.new(obs_db_path, mission, obsid, products=[product])
    summary = {
        "obsid": obsid,
        "mission": mission,
        "level2_evt": product["path"],
        "catalog_entry": entry,
        "downloaded": True,
        "locked": locked,
    }
    if locked:
        summary["note"] = (
            f"{obsid} was served GPG-encrypted (proprietary period) and is recorded "
            f"as present but not directly usable. Decrypt the local copy out of band, "
            f"or acquire again once the observation is public."
        )
    return summary
