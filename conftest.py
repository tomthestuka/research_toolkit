"""Test classes for the whole toolkit, so a slow suite can be told from a broken one.

Three classes, by marker:

* unmarked — fast and deterministic: synthetic or tiny inputs, no network, seconds per test.
  The default for every test.
* ``slow`` — runs on real observation data or a real model and takes minutes. Runs by default;
  ``pytest -m "not slow"`` gives the quick suite.
* ``live`` — reaches the network or an external service (an archive, a web API). Skipped unless
  asked for with ``--live`` or ``TOOLKIT_LIVE=1``; the authoring gate deselects it (no network).
* ``data("<name>")`` — reads the named data set (real observations), whose folder the variable
  ``TOOLKIT_DATA_<NAME>`` names; the test reads its folder from there. Skipped, with the reason,
  when the variable is unset or its folder absent. The authoring gate sets it for every data set
  its host's ``test_data`` names and binds the folder read-only; a data test it skips fails the
  change. Combine with ``slow`` when it takes minutes.

Mark a module with ``pytestmark = pytest.mark.slow`` (or ``live``), or one test with
``@pytest.mark.slow`` (``@pytest.mark.data("cygx1_ixpe")``).
"""

from __future__ import annotations

import os

import pytest


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption("--live", action="store_true", default=False,
                     help="also run tests marked live (network or external services)")


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line("markers", "slow: runs on real data or a real model and takes minutes")
    config.addinivalue_line("markers", "live: reaches the network or an external service; opt-in")
    config.addinivalue_line("markers", "data(name): reads the named data set from the folder TOOLKIT_DATA_<NAME> names")


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    live = config.getoption("--live") or os.environ.get("TOOLKIT_LIVE") == "1"
    for item in items:
        if item.get_closest_marker("live") and not live:
            item.add_marker(pytest.mark.skip(reason="live: run with --live or TOOLKIT_LIVE=1"))
        for mark in item.iter_markers("data"):
            name = str(mark.args[0]) if mark.args else ""
            variable = f"TOOLKIT_DATA_{name.upper()}"
            folder = os.environ.get(variable)
            if not folder or not os.path.isdir(folder):
                item.add_marker(pytest.mark.skip(
                    reason=f"data: {name} is absent ({variable} {'names no folder' if folder else 'is unset'})"))
