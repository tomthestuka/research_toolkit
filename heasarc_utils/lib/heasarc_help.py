"""High-level acquire entry point — download an ObsID and record it.

This is the abstracted, human-callable API (and the function the LLM tool wraps).
It is the first producer wired into ``obsdb_catalog``: after the archive data is
on disk it records the level-2 event through the single trusted writer
(``obsdb_catalog.neutral.new``) — the only path that may mutate ``catalog.json``.
"""

from __future__ import annotations

from obsdb_catalog import neutral

from . import heasarc

__all__ = ["acquire"]


def acquire(obsid: str, mission: str, obs_db_path: str,
            *, host: str | None = None, chatter: bool = False) -> dict:
    """Acquire archive data for ``obsid`` into ``obs_db_path`` and record it.

    Resolves the ObsID authoritatively via the mission's HEASARC master table,
    downloads it (astroquery only — no fallback), confirms the representative
    level-2 event is on disk, and records it via ``obsdb_catalog.neutral.new``
    (level 2). **Always (re)downloads** — the caller decides when to acquire; a
    deliberate re-acquire regresses ``max_level`` to 2, flagging stale downstream
    products for re-reduction.

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

    Returns
    -------
    dict
        ``{"obsid", "mission", "level2_evt": <relative path>, "catalog_entry"}``.

    Raises
    ------
    heasarc.HeasarcError
        Unknown mission, ObsID not found, or the download did not yield the
        expected level-2 file. Nothing is recorded on failure.
    """
    mission = mission.lower()
    product = heasarc.acquire_obsid(mission, obsid, obs_db_path, host=host)
    if chatter:
        print(f"[heasarc_utils] {mission} {obsid}: recording {product['path']}")
    entry = neutral.new(obs_db_path, mission, obsid, products=[product])
    return {
        "obsid": obsid,
        "mission": mission,
        "level2_evt": product["path"],
        "catalog_entry": entry,
    }
