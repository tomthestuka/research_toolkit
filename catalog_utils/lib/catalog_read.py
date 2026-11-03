"""Local data-catalog read + listing (pure standard library).

Parses a mission's ``<obs_db_path>/catalog.json`` via the single trusted reader
(``obsdb_catalog.neutral.read``) and renders the "what's already local" listing an
analyst consults before requesting acquisition. No science-stack import, so this
loads and runs in the orchestration environment as well as the sandbox.

The caller supplies ``obs_db_path`` (the orchestration layer resolves the mission's
directory from its own obs-db-location config); this module never reads that config.

Each line names the observed source when its event file's header records one (the FITS
``OBJECT`` keyword), so a source can be found by name as well as by ObsID.

Beyond "is it catalogued", the listing flags an ObsID whose data is present but
**GPG-encrypted** — how HEASARC serves an observation still inside its proprietary
period. Such data is not directly usable: decrypt it, or acquire it again once the
observation is public. The flag is derived live from the recorded product paths, since
encryption state is mutable (a human may decrypt out of band).
"""

from __future__ import annotations

import gzip
import os
import re
import stat

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


#: A source name as the listing shows it: printable, and free of the characters the listing itself
#: uses (``:``, ``(``, ``)``, ``[``, ``]``), so header bytes can never forge a line or a mark.
_NAME_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9 +\-.'*_/]{0,47}")
#: Headers read per event file: the primary and the first extension (IXPE keeps OBJECT in EVENTS).
_HEADERS_READ = 2
#: A header longer than this many blocks is not a header this listing reads.
_MAX_HEADER_BLOCKS = 36


def _string_value(card: str) -> str | None:
    """A card's quoted string value per FITS (``''`` is one quote), or ``None`` when unquoted."""
    rest = card[10:].lstrip()
    if not rest.startswith("'"):
        return None
    out, i = [], 1
    while i < len(rest):
        if rest[i] == "'":
            if rest[i + 1:i + 2] == "'":
                out.append("'")
                i += 2
                continue
            return "".join(out).rstrip()
        out.append(rest[i])
        i += 1
    return None


def _header_cards(fh) -> dict[str, str]:
    """The cards of the first ``_HEADERS_READ`` headers, by keyword (the first one wins), read block
    by block up to each END and stepping over the primary's data, if any."""
    cards: dict[str, str] = {}
    for _ in range(_HEADERS_READ):
        header: list[str] = []
        for _block in range(_MAX_HEADER_BLOCKS):
            block = fh.read(2880)
            if len(block) < 2880:
                return cards
            text = block.decode("ascii", "replace")
            header.extend(text[k:k + 80] for k in range(0, 2880, 80))
            if any(card[:8] == "END     " for card in header[-36:]):
                break
        else:
            return cards
        for card in header:
            cards.setdefault(card[:8].rstrip(), card)
        here = {card[:8].rstrip(): card[10:30].strip() for card in header if card[8:10] == "= "}
        try:
            naxis = int(here.get("NAXIS", "0"))
            size = abs(int(here.get("BITPIX", "8"))) // 8 if naxis else 0
            for n in range(1, naxis + 1):
                size *= int(here.get(f"NAXIS{n}", "0"))
            size *= int(here.get("GCOUNT", "1") or 1) if "XTENSION" in here else 1
        except ValueError:
            return cards
        fh.seek(-(-size // 2880) * 2880, 1)  # the HDU's data, padded to whole blocks
    return cards


def _source_name(obs_db_path: str, entry: dict) -> str | None:
    """The source an entry observed, as its first readable event file's header names it (FITS
    ``OBJECT``; without one, its ``RA_OBJ``/``DEC_OBJ``), or ``None``. Only a regular file inside
    ``obs_db_path`` is read, and a name outside :data:`_NAME_RE` is not shown."""
    root = os.path.realpath(obs_db_path)
    for product in entry.get("products") or []:
        rel = (product or {}).get("path")
        if (product or {}).get("kind") != "level2_evt" or not rel or neutral.is_locked(rel):
            continue
        path = os.path.realpath(os.path.join(root, rel))
        if os.path.commonpath([root, path]) != root or not stat.S_ISREG(_mode(path)):
            continue
        try:
            with (gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")) as fh:
                cards = _header_cards(fh)
        except (OSError, EOFError):
            continue
        name = _string_value(cards["OBJECT"]) if "OBJECT" in cards else None
        if name and _NAME_RE.fullmatch(name.strip()):
            return name.strip()
        try:
            ra, dec = (float(cards[k][10:30].split("/")[0]) for k in ("RA_OBJ", "DEC_OBJ"))
        except (KeyError, ValueError):
            continue
        return f"unnamed in its header; RA {ra:.3f}, Dec {dec:.3f}"
    return None


def _mode(path: str) -> int:
    try:
        return os.stat(path).st_mode
    except OSError:
        return 0


def data_catalog(obs_db_path, min_level=None) -> str:
    """List locally-available ObsIDs for the mission stored at ``obs_db_path``.

    Reads ``<obs_db_path>/catalog.json`` through the trusted ``obsdb_catalog`` reader
    and returns a text listing: one line per catalogued ObsID with its highest
    processing level (HEASARC supplies up to level 2) and the product kinds derived
    from it, and the observed source as its event file's header names it (or that header's
    coordinates). An ObsID whose data is present but GPG-encrypted carries LOCKED, NOT
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
        name = _source_name(obs_db_path, entry)
        line = f"- {obsid}{f' ({name})' if name else ''}: max_level={max_level}, products=[{products}]"
        if has_locked:
            # On the ObsID's own line, not only as a trailing note: the level is
            # honest about what the products are, so this is the part that says
            # they cannot be read.
            line += _LOCKED_MARK + (" (some products)" if has_usable else "")
            line += "\n" + _LOCKED_NOTE
        lines.append(line)
    body = "\n".join(lines) or "(no local data catalogued)"
    return f"Local {mission} data ({len(lines)} ObsIDs):\n{body}"
