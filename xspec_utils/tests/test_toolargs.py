"""Tests for xspec_utils.lib.toolargs.coerce_tool_args -- the sandbox
numeric/boolean argument-marshalling fix (fix-ixpe-qu-slice-pcubes-numeric-
args). Pure logic, no real data required: pins the coercion rules against a
hand-written table computed independently of the module's own coercers
(plain int()/float()/bool literals), plus the rejection cases the sandbox
harness's bare ``fn(**args)`` call can actually produce.
"""

from __future__ import annotations

import pytest

from xspec_utils.lib.toolargs import coerce_tool_args

# --------------------------------------------------------------------------- #
# Differential: independently-computed expected values (not the module's own
# helpers) for every accept path the spec documents.
# --------------------------------------------------------------------------- #

PASS_CASES = [
    # (spec, values, key, expected)
    ({"n": "int"}, {"n": 8}, "n", 8),                        # native int, identity
    ({"n": "int"}, {"n": 8.0}, "n", int(8.0)),                # integral float -> int
    ({"n": "int"}, {"n": "8"}, "n", int(float("8"))),         # '8' -> 8
    ({"n": "int"}, {"n": "8.0"}, "n", int(float("8.0"))),     # '8.0' -> 8
    ({"x": "float"}, {"x": 0.5}, "x", float(0.5)),            # native float, identity
    ({"x": "float"}, {"x": 5}, "x", float(5)),                # native int -> float
    ({"x": "float"}, {"x": "1e-12"}, "x", float("1e-12")),
    ({"x": "float"}, {"x": " 0.5 "}, "x", float(" 0.5 ".strip())),
    ({"x": "float"}, {"x": "-1.0"}, "x", float("-1.0")),
    ({"b": "bool"}, {"b": True}, "b", True),                  # native bool, identity
    ({"b": "bool"}, {"b": False}, "b", False),
    ({"b": "bool"}, {"b": "true"}, "b", True),
    ({"b": "bool"}, {"b": "False"}, "b", False),
    ({"b": "bool"}, {"b": 1}, "b", True),
    ({"b": "bool"}, {"b": 0}, "b", False),
    ({"o": "float?"}, {"o": None}, "o", None),
    ({"o": "float?"}, {"o": ""}, "o", None),
    ({"o": "float?"}, {"o": "  NoNe  "}, "o", None),
    ({"o": "int?"}, {"o": "null"}, "o", None),
    ({"o": "bool?"}, {"o": None}, "o", None),
]


@pytest.mark.parametrize("spec, values, key, expected", PASS_CASES)
def test_coerce_tool_args_matches_stdlib(spec, values, key, expected):
    out = coerce_tool_args(spec, values)
    assert out[key] == expected
    assert type(out[key]) is type(expected)


# --------------------------------------------------------------------------- #
# Rejection cases: lossy or nonsensical input must raise, never truncate or
# silently fall back to a default.
# --------------------------------------------------------------------------- #

FAIL_CASES = [
    ({"x": "float"}, {"x": "abc"}, "x"),          # not parseable
    ({"n": "int"}, {"n": "8.4"}, "n"),            # non-integral string
    ({"n": "int"}, {"n": 8.4}, "n"),              # non-integral float
    ({"n": "int"}, {"n": True}, "n"),             # bool rejected for an int slot
    ({"n": "int"}, {"n": None}, "n"),             # None for a required arg
    ({"x": "float"}, {"x": None}, "x"),           # None for a required arg
    ({"x": "float"}, {"x": "nan"}, "x"),          # non-finite
    ({"x": "float"}, {"x": "inf"}, "x"),          # non-finite
    ({"b": "bool"}, {"b": "maybe"}, "b"),         # not a recognised boolean token
    ({"b": "bool"}, {"b": 2}, "b"),               # int outside {0, 1}
]


@pytest.mark.parametrize("spec, values, key", FAIL_CASES)
def test_coerce_tool_args_rejects_bad_input(spec, values, key):
    with pytest.raises(ValueError) as excinfo:
        coerce_tool_args(spec, values)
    msg = str(excinfo.value)
    assert key in msg
    assert repr(values[key]) in msg


# --------------------------------------------------------------------------- #
# Contract: never mutate the input, always return a new dict, omit absent
# keys rather than inventing them, and pass unspecified keys through as-is.
# --------------------------------------------------------------------------- #

def test_does_not_mutate_input_and_returns_new_dict():
    values = {"n": 8, "x": "1.5"}
    spec = {"n": "int", "x": "float"}
    out = coerce_tool_args(spec, values)
    assert out is not values
    assert values == {"n": 8, "x": "1.5"}
    assert out == {"n": 8, "x": 1.5}


def test_absent_key_is_omitted_not_invented():
    out = coerce_tool_args({"n": "int"}, {})
    assert "n" not in out


def test_unspecced_keys_pass_through_unchanged():
    payload = {"a": 1}
    out = coerce_tool_args({"n": "int"}, {"n": "5", "other": payload, "path": "/x"})
    assert out["n"] == 5
    assert out["other"] is payload
    assert out["path"] == "/x"


# --------------------------------------------------------------------------- #
# `context`: prefixes every raised message with the caller's own tool name,
# reproducing the brief's exact example string.
# --------------------------------------------------------------------------- #

def test_error_message_includes_tool_context():
    with pytest.raises(ValueError) as excinfo:
        coerce_tool_args({"n_min": "int"}, {"n_min": "eight"}, context="ixpe_qu_slice_pcubes")
    assert str(excinfo.value) == "ixpe_qu_slice_pcubes: argument 'n_min' expects an integer, got 'eight' (str)"


def test_error_message_without_context_unchanged():
    with pytest.raises(ValueError) as excinfo:
        coerce_tool_args({"n_min": "int"}, {"n_min": "eight"})
    assert str(excinfo.value) == "argument 'n_min' expects an integer, got 'eight' (str)"


# --------------------------------------------------------------------------- #
# Malformed spec: a programmer error must be a loud, attributable ValueError,
# never a bare KeyError a caller's broad `except` swallows into a mystery.
# --------------------------------------------------------------------------- #

def test_unsupported_type_code_raises_value_error():
    with pytest.raises(ValueError) as excinfo:
        coerce_tool_args({"n": "integer"}, {"n": 8})
    msg = str(excinfo.value)
    assert "'n'" in msg
    assert "'integer'" in msg
