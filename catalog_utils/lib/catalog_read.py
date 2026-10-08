"""Local data-catalog read + listing (pure standard library).

Parses a mission's ``<obs_db_path>/catalog.json`` via the single trusted reader
(``obsdb_catalog.neutral.read``) and renders the "what's already local" listing an
analyst consults before requesting acquisition. No science-stack import, so this
loads and runs in the orchestration environment as well as the sandbox.

The caller supplies ``obs_db_path`` (the orchestration layer resolves the mission's
directory from its own obs-db-location config); this module never reads that config.

Beyond "is it catalogued", the listing flags an ObsID whose data is present but
**GPG-encrypted** — how HEASARC serves an observation still inside its proprietary
period. Such data is not directly usable: decrypt it, or acquire it again once the
observation is public. The flag is derived live from the recorded product paths, since
encryption state is mutable (a human may decrypt out of band).
"""

from __future__ import annotations

import os

from obsdb_catalog import neutral

_LOCKED_MARK = " -- LOCKED, NOT USABLE"
_LOCKED_NOTE = (
    "  [GPG-encrypted on disk (proprietary period): the data is present but cannot be "
    "read by any analysis tool as-is. Decrypt the local copy out of band, or acquire it "
    "again once the observation becomes public -- re-acquiring while it is still "
    "proprietary returns the same encrypted bytes]"
)


def _usable_and_locked(obs_db_path: str, entry: dict) -> tuple[bool, bool]:
    """Split an entry's recorded products into ``(has_usable, has_locked)`` — live.

    Encryption reaches a catalogued product from two directions, and both count as
    locked:

    * the recorded path **is** the encrypted copy and it is on disk — the archive
      served the observation that way, inside its proprietary period;
    * the recorded path is a readable file that has since **gone**, with its ``.gpg``
      twin sitting where it was — someone encrypted it locally after it was recorded.

    Either way the product is present and unreadable. Both are read from the paths
    on disk at call time rather than from a stored flag, because encryption state
    changes out of band and a stored boolean would drift out of agreement with the
    path beside it. An entry can be usable *and* locked at once; saying so beats
    collapsing the ObsID to one verdict.
    """
    has_usable = has_locked = False
    for product in entry.get("products") or []:
        rel = (product or {}).get("path")
        if not rel:
            continue
        if os.path.exists(os.path.join(obs_db_path, rel)):
            if neutral.is_locked(rel):
                has_locked = True
            else:
                has_usable = True
        elif os.path.exists(os.path.join(obs_db_path, rel + neutral.LOCK_SUFFIX)):
            has_locked = True
    return has_usable, has_locked


def data_catalog(obs_db_path, min_level=None) -> str:
    """List locally-available ObsIDs for the mission stored at ``obs_db_path``.

    Reads ``<obs_db_path>/catalog.json`` through the trusted ``obsdb_catalog`` reader
    and returns a text listing: one line per catalogued ObsID with its highest
    processing level (HEASARC supplies up to level 2) and the product kinds derived
    from it. An ObsID whose data is present but GPG-encrypted carries LOCKED, NOT
    USABLE on its own line: it is catalogued at the level its products genuinely are,
    and the mark is what says the bytes cannot be read yet. ``min_level`` (optional)
    restricts the listing to ObsIDs available at or above that processing level — a
    locked ObsID is levelled honestly and so still appears, which is why the mark, not
    the level, is what a caller must read. An absent catalog lists zero ObsIDs. Raises
    ``obsdb_catalog.neutral.CatalogError`` on a present-but-malformed catalog (a torn
    inventory fails loudly rather than parsing to junk).
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
        has_usable, has_locked = _usable_and_locked(obs_db_path, entry)
        line = f"- {obsid}: max_level={max_level}, products=[{products}]"
        if has_locked:
            # On the ObsID's own line, not only as a trailing note: the level is
            # honest about what the products are, so this is the part that says
            # they cannot be read.
            line += _LOCKED_MARK + (" (some products)" if has_usable else "")
            line += "\n" + _LOCKED_NOTE
        lines.append(line)
    body = "\n".join(lines) or "(no local data catalogued)"
    return f"Local {mission} data ({len(lines)} ObsIDs):\n{body}"
