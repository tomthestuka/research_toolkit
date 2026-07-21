"""Local data-catalog read + listing (pure standard library).

Parses a mission's ``<obs_db_path>/catalog.json`` via the single trusted reader
(``obsdb_catalog.neutral.read``) and renders the "what's already local" listing an
analyst consults before requesting acquisition. No science-stack import, so this
loads and runs in the orchestration environment as well as the sandbox.

The caller supplies ``obs_db_path`` (the orchestration layer resolves the mission's
directory from its own obs-db-location config); this module never reads that config.

Beyond "is it catalogued", the listing flags an ObsID whose data is present on disk but
only **GPG-encrypted** (a local at-rest convention) — a live filesystem check at read
time, since encryption state is mutable (a human may decrypt out of band). Such data is
not directly usable: decrypt it, or re-acquire it if the remote archive is unlocked.
"""

from __future__ import annotations

import os

from obsdb_catalog import neutral

_LOCKED_NOTE = (
    "  [LOCKED ON DISK: level-2 data present but GPG-encrypted — not directly usable; "
    "decrypt the local archive out of band, or re-acquire if the remote archive is unlocked]"
)


def _locked_on_disk(obs_db_path: str, obsid: str, entry: dict) -> bool:
    """True when the ObsID's recorded products are absent on disk but a ``.gpg`` copy exists.

    The data is present but encrypted (locked): unusable until decrypted. Detected live,
    not stored, because encryption state changes out of band. Cheap: only walks the ObsID
    directory when no recorded product file is actually present.
    """
    products = entry.get("products") or []
    if not products:
        return False
    for p in products:
        rel = (p or {}).get("path")
        if rel and os.path.exists(os.path.join(obs_db_path, rel)):
            return False  # at least one recorded product is present & usable
    obsid_dir = os.path.join(obs_db_path, str(obsid))
    if not os.path.isdir(obsid_dir):
        return False
    for _root, _dirs, files in os.walk(obsid_dir):
        if any(f.endswith(".gpg") for f in files):
            return True
    return False


def data_catalog(obs_db_path, min_level=None) -> str:
    """List locally-available ObsIDs for the mission stored at ``obs_db_path``.

    Reads ``<obs_db_path>/catalog.json`` through the trusted ``obsdb_catalog`` reader
    and returns a text listing: one line per catalogued ObsID with its highest
    processing level (HEASARC supplies up to level 2) and the product kinds derived
    from it. An ObsID whose products are present but only GPG-encrypted is flagged
    LOCKED ON DISK (not directly usable). ``min_level`` (optional) restricts the listing
    to ObsIDs available at or above that processing level. An absent catalog lists zero
    ObsIDs. Raises ``obsdb_catalog.neutral.CatalogError`` on a present-but-malformed
    catalog (a torn inventory fails loudly rather than parsing to junk).
    """
    obs_db_path = str(obs_db_path)
    catalog = neutral.read(obs_db_path)
    mission = catalog.get("mission", "") or ""
    entries = catalog.get("entries") or {}
    lines: list[str] = []
    for obsid in sorted(entries):
        entry = entries[obsid] or {}
        max_level = entry.get("max_level")
        if min_level is not None and not (isinstance(max_level, int) and max_level >= int(min_level)):
            continue
        kinds = sorted({(p or {}).get("kind", "?") for p in (entry.get("products") or [])})
        products = ", ".join(kinds) or "(none)"
        line = f"- {obsid}: max_level={max_level}, products=[{products}]"
        if _locked_on_disk(obs_db_path, obsid, entry):
            line += "\n" + _LOCKED_NOTE
        lines.append(line)
    body = "\n".join(lines) or "(no local data catalogued)"
    return f"Local {mission} data ({len(lines)} ObsIDs):\n{body}"
