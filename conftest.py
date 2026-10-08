"""Test classes for the whole toolkit, so a slow suite can be told from a broken one.

Three classes, by marker:

* unmarked — fast and deterministic: synthetic or tiny inputs, no network, seconds per test.
  The default for every test.
* ``slow`` — runs on real observation data or a real model and takes minutes. Runs by default;
  ``pytest -m "not slow"`` gives the quick suite.
* ``live`` — reaches the network or an external service (an archive, a web API). Skipped unless
  asked for with ``--live`` or ``TOOLKIT_LIVE=1``; the authoring gate's test run has no network.

Mark a module with ``pytestmark = pytest.mark.slow`` (or ``live``), or one test with
``@pytest.mark.slow``.
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


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    if config.getoption("--live") or os.environ.get("TOOLKIT_LIVE") == "1":
        return
    skip = pytest.mark.skip(reason="live: run with --live or TOOLKIT_LIVE=1")
    for item in items:
        if item.get_closest_marker("live"):
            item.add_marker(skip)
