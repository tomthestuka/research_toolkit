"""Numeric/boolean tool-argument coercion -- pure stdlib.

The sandbox harness invokes a tool as a bare ``result = fn(**args)`` with no
schema coercion, so a numeric ``TOOL_MANIFEST`` argument declared ``integer``/
``number`` can arrive at the wrapped function as a Python ``str`` (e.g.
``n_max="3"``). A function whose body compares or multiplies that value
against a native ``int``/``float`` then raises a ``TypeError`` deep inside its
own logic, well after any early, attributable check.

``coerce_tool_args`` normalises a declared subset of a call's keyword
arguments to their real Python types (or raises a ``ValueError`` naming the
offending argument) before the wrapped function's own body runs. Stdlib-only
-- no numpy, no xspec, no ixpeobssim, no ixpeobssim-importing sibling module
-- so both an ``xspec``-dependent ``skills/session.py`` wrapper and a pure-
numpy ``lib/`` module can import it without pulling in the science stack.
"""

from __future__ import annotations

import math

_BOOL_TRUE = {"true", "t", "yes", "y", "1"}
_BOOL_FALSE = {"false", "f", "no", "n", "0"}
_NONE_LIKE = {"", "none", "null"}


def _prefix(context: str | None) -> str:
    return f"{context}: " if context is not None else ""


def _fail(name: str, value, expects: str, context: str | None = None) -> None:
    raise ValueError(
        f"{_prefix(context)}argument {name!r} expects {expects}, "
        f"got {value!r} ({type(value).__name__})"
    )


def _fail_bad_spec(name: str, base: str, context: str | None = None) -> None:
    raise ValueError(
        f"{_prefix(context)}argument {name!r} has unsupported type code "
        f"{base!r} in coerce_tool_args spec"
    )


def _coerce_bool(name: str, value, context: str | None = None):
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        s = value.strip().lower()
        if s in _BOOL_TRUE:
            return True
        if s in _BOOL_FALSE:
            return False
        _fail(name, value, "a boolean", context)
    if isinstance(value, int) and value in (0, 1):
        return bool(value)
    _fail(name, value, "a boolean", context)


def _coerce_int(name: str, value, context: str | None = None):
    if isinstance(value, bool):
        _fail(name, value, "an integer", context)
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            _fail(name, value, "an integer", context)
        if value.is_integer():
            return int(value)
        _fail(name, value, "an integer", context)
    if isinstance(value, str):
        try:
            f = float(value.strip())
        except ValueError:
            _fail(name, value, "an integer", context)
        if not math.isfinite(f):
            _fail(name, value, "an integer", context)
        if f.is_integer():
            return int(f)
        _fail(name, value, "an integer", context)
    _fail(name, value, "an integer", context)


def _coerce_float(name: str, value, context: str | None = None):
    if isinstance(value, bool):
        _fail(name, value, "a floating-point number", context)
    if isinstance(value, (int, float)):
        f = float(value)
        if not math.isfinite(f):
            _fail(name, value, "a floating-point number", context)
        return f
    if isinstance(value, str):
        try:
            f = float(value.strip())
        except ValueError:
            _fail(name, value, "a floating-point number", context)
        if not math.isfinite(f):
            _fail(name, value, "a floating-point number", context)
        return f
    _fail(name, value, "a floating-point number", context)


_COERCERS = {"bool": _coerce_bool, "int": _coerce_int, "float": _coerce_float}


def coerce_tool_args(spec: dict, values: dict, context: str | None = None) -> dict:
    """Normalise ``values`` to native Python types per ``spec``.

    ``spec`` maps argument name -> type code (``"int"``, ``"float"``,
    ``"bool"``, each optionally suffixed with ``"?"`` to mark the argument
    optional). Returns a NEW dict; ``values`` is never mutated.

    ``context``, when given (typically the tool's own name, e.g.
    ``"ixpe_qu_slice_pcubes"``), is prefixed to every raised message as
    ``f"{context}: ..."``, so a caller wrapping several coerced arguments
    for one tool gets an attributable error. Omitting it (``None``, the
    default) leaves the message exactly as it would read with no tool
    context at all.

    Per key present in ``spec``:
      - A type code (after stripping any trailing ``"?"``) other than
        ``"int"``/``"float"``/``"bool"`` raises ``ValueError`` naming the
        argument and the offending code -- a malformed spec is a
        programmer error, but it must be a loud, attributable one, never a
        bare ``KeyError`` a caller's broad ``except`` swallows into a
        mystery.
      - Absent from ``values``: omitted from the result (the caller's own
        signature default then applies untouched); never invented.
      - Optional (``"?"``): ``None``, or a str whose ``.strip().lower()`` is
        ``""``/``"none"``/``"null"``, coerces to ``None``.
      - ``"bool"``: native ``bool`` passes through; a str
        ``true``/``t``/``yes``/``y``/``1`` (case-insensitive, stripped) ->
        ``True``, ``false``/``f``/``no``/``n``/``0`` -> ``False``; native
        ``0``/``1`` (not ``bool``) -> ``False``/``True``.
      - ``"int"``: native ``int`` (not ``bool``) passes through; a native or
        str-parsed ``float`` coerces only when integral (``8.0`` -> ``8``,
        ``8.4`` raises); a str is ``float()``-parsed first, so ``"8"`` and
        ``"8.0"`` both give ``8``.
      - ``"float"``: native ``int``/``float`` (not ``bool``) coerces to
        ``float``; a str is ``float()``-parsed after stripping. Non-finite
        results (``nan``/``inf``) raise.
      - Any failure raises ``ValueError`` naming the argument, ``repr()`` of
        the received value, and ``type(value).__name__``. Never truncates,
        never silently falls back to a default.

    Keys present in ``values`` but absent from ``spec`` are copied through
    unchanged -- this is a filter over the declared numeric/boolean
    arguments, not a whitelist of the whole call.
    """
    out: dict = {}
    for name, typecode in spec.items():
        if name not in values:
            continue
        optional = typecode.endswith("?")
        base = typecode[:-1] if optional else typecode
        if base not in _COERCERS:
            _fail_bad_spec(name, base, context)
        value = values[name]
        if optional and (
            value is None
            or (isinstance(value, str) and value.strip().lower() in _NONE_LIKE)
        ):
            out[name] = None
            continue
        out[name] = _COERCERS[base](name, value, context)

    for name, value in values.items():
        if name not in spec:
            out[name] = value
    return out
