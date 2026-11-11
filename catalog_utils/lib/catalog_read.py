"""The data catalog as an analyst reads it: which products the store holds (standard library,
``storelib``, and ``heasarc_utils``' mission profiles to find a raw download's level-2 event).

One line per product the catalog holds as committed or verified: the current product of each
(bucket, kind, subject) and any made after it, newest first, with its record (what a tool that reads it takes).
A raw download names the observed source as its level-2 event's header does (the FITS ``OBJECT``
keyword, read from the first blocks of the file by one version-pinned ranged GET), and is marked
**LOCKED, NOT USABLE** when that event is GPG-encrypted (an observation still in its proprietary
period).
"""

from __future__ import annotations

import io
import re
import json
import posixpath

import storelib

_LOCKED_MARK = " -- LOCKED, NOT USABLE"
_LOCKED_NOTE = (
    "  [GPG-encrypted (proprietary period): the data is stored but no analysis tool can read it as is. "
    "Acquire it again once the observation is public; until then acquiring returns this copy]"
)
#: Bytes read from the head of an event file: two headers of at most ``_MAX_HEADER_BLOCKS`` blocks.
_HEAD_BYTES = 2 * 36 * 2880

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


def _source_name(s3, record: dict) -> str | None:
    """The source a level-2 event file names in its header (FITS ``OBJECT``; without one, its
    ``RA_OBJ``/``DEC_OBJ``), or ``None``. A name outside :data:`_NAME_RE` is not shown."""
    try:
        head = s3.get_object(Bucket=record["bucket"], Key=record["key"], VersionId=record["version_id"],
                             Range=f"bytes=0-{_HEAD_BYTES - 1}")["Body"].read()
    except Exception:  # noqa: BLE001 - a header that cannot be read leaves the line without a name
        return None
    cards = _header_cards(io.BytesIO(head))
    name = _string_value(cards["OBJECT"]) if "OBJECT" in cards else None
    if name and _NAME_RE.fullmatch(name.strip()):
        return name.strip()
    try:
        ra, dec = (float(cards[k][10:30].split("/")[0]) for k in ("RA_OBJ", "DEC_OBJ"))
    except (KeyError, ValueError):
        return None
    return f"unnamed in its header; RA {ra:.3f}, Dec {dec:.3f}"


_QUERY = """
SELECT o.store, o.bucket, o.key, o.version_id, p.recipe_hash, p.kind, p.subject, p.state,
       p.recipe_hash = c.recipe_hash AS is_current,
       (SELECT jsonb_agg(jsonb_build_array(x.key, x.version_id) ORDER BY x.key) FROM catalog.object x
         WHERE x.recipe_hash = p.recipe_hash) AS files
FROM catalog.product p
JOIN catalog.object o ON o.recipe_hash = p.recipe_hash AND o.key ~ '^[^/]+/[^/]+/[^/]+/recipe\\.json$'
LEFT JOIN catalog.current c ON c.bucket = p.bucket AND c.kind = p.kind AND c.subject = p.subject
LEFT JOIN catalog.product cp ON cp.recipe_hash = c.recipe_hash
WHERE p.state IN ('committed', 'verified')
  AND (cp.recipe_hash IS NULL OR p.recipe_hash = cp.recipe_hash OR p.claimed_at > cp.claimed_at)
  AND (%(bucket)s::text IS NULL OR o.bucket = %(bucket)s)
  AND (%(kind)s::text IS NULL OR p.kind = %(kind)s)
  AND (%(subject)s::text IS NULL OR p.subject = %(subject)s)
ORDER BY o.bucket, p.subject, p.kind, p.claimed_at DESC
"""


def data_catalog(mission=None, kind=None, subject=None) -> str:
    """List the products the store's catalog holds, as text.

    One line per product the catalog holds as committed or verified, grouped by bucket (the mission):
    for each (bucket, kind, subject) the current product and any made after it, newest first. Each line names
    the subject (the ObsID), the kind, its state (``current`` when the catalog's pointer names it;
    ``verified`` once the runtime has checked its bytes; ``committed`` before that), and its record as
    JSON, which a tool reading the product takes. A raw download adds the observed source as its
    level-2 event's header names it, and the mark LOCKED, NOT USABLE when that event is encrypted.
    ``mission``, ``kind`` and ``subject`` keep the products of that bucket, kind or subject.
    """
    from heasarc_utils.lib import heasarc

    with storelib.connect() as conn:
        rows = conn.execute(_QUERY, {"bucket": mission, "kind": kind, "subject": subject}).fetchall()
    s3 = storelib.client() if any(r[5] == "raw" for r in rows) else None
    lines, bucket = [], None
    for store, bkt, key, version_id, _hash, kind_, subject_, state, is_current, pairs in rows:
        if bkt != bucket:
            bucket = bkt
            lines.append(f"{bkt}:")
        record = {"store": store, "bucket": bkt, "key": key, "version_id": version_id}
        folder = posixpath.dirname(key)
        files = {posixpath.relpath(k, folder): v for k, v in (pairs or []) if k != key}
        name = locked = None
        profile = heasarc.MISSIONS.get(bkt)
        if kind_ == "raw" and profile is not None:
            try:
                level2 = profile.level2(sorted(files), subject_)
            except heasarc.HeasarcError:
                level2 = None
            if level2 is not None:
                locked = heasarc.is_locked(level2)
                if not locked:
                    name = _source_name(s3, {"bucket": bkt, "key": f"{folder}/{level2}",
                                             "version_id": files[level2]})
        status = "current" if is_current else state
        line = (f"- {subject_}{f' ({name})' if name else ''} {kind_} [{status}, {len(files)} files]: "
                f"record {json.dumps(record)}")
        if locked:
            line += _LOCKED_MARK + "\n" + _LOCKED_NOTE
        lines.append(line)
    count = sum(1 for line in lines if line.startswith("- "))
    body = "\n".join(lines) or "(no product stored)"
    return f"Stored products ({count}):\n{body}"
