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

A test that reads or writes a data store takes the ``store`` fixture: the throwaway store a host gives
the test process. ``TOOLKIT_STORE`` names it (the name records carry); its producer's keys arrive in the
variables a tool process gets (``AWS_ACCESS_KEY_ID``, ``AWS_SECRET_ACCESS_KEY``, ``AWS_ENDPOINT_URL_S3``,
``AWS_REGION``, ``PGHOST``, ``PGPORT``, ``PGUSER``, ``PGPASSWORD``, ``PGDATABASE``); it holds the
buckets ``ixpe``, ``nicer`` and ``nustar``, versioned. Skipped, with the reason, when ``TOOLKIT_STORE``
is unset; the authoring gate starts a fresh one for each test run and fails a change whose store test
it skipped. A test's inputs reach that store only through the host's seed channel (the ``seed``
fixture, ``TOOLKIT_STORE_SEED``): a synthetic product made of given files, or a product or a file of a
named data set; a tool's own call is the only writer the test drives.
"""

from __future__ import annotations

import json
import os
import socket
import uuid

import pytest

#: The variable naming the throwaway store a test run has, and the variables its keys arrive in.
STORE_NAME = "TOOLKIT_STORE"
STORE_VARIABLES = ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_ENDPOINT_URL_S3", "PGHOST", "PGUSER",
                   "PGDATABASE")


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


@pytest.fixture
def store(monkeypatch):
    """The throwaway store the host names (see the module docstring): ``{"store": <its name>, "tag": <text
    unique to this test>}``; a product's subject or settings carry the tag, so a recipe is never one an
    earlier test stored. The provenance variables a product records are set to the test's own values
    when the host sets none."""
    if not os.environ.get(STORE_NAME):
        pytest.skip(f"store: no throwaway store ({STORE_NAME} is unset; the CI pipeline provides one)")
    missing = [n for n in STORE_VARIABLES if not os.environ.get(n)]
    if missing:
        pytest.fail(f"store: {STORE_NAME} names a store, but {', '.join(missing)} unset")
    monkeypatch.setenv("AR_TOOLKIT_COMMIT", os.environ.get("AR_TOOLKIT_COMMIT") or "0" * 40)
    monkeypatch.setenv("AR_ENV_FINGERPRINT", os.environ.get("AR_ENV_FINGERPRINT") or "test-environment")
    return {"store": os.environ[STORE_NAME], "tag": uuid.uuid4().hex[:12]}


#: The unix socket a test asks the host's seed channel through (one JSON line each way).
SEED_VARIABLE = "TOOLKIT_STORE_SEED"


@pytest.fixture
def seed(store):
    """``seed(**request) -> answer``: one request to the host's seed channel. A synthetic product:
    ``seed(kind=, subject=, bucket=, files={"<relative path>": <base64 text>})``; a data set's product:
    ``seed(data=, paths=[...], bucket=, kind=, subject=)``; one data-set file: ``seed(data=, path=,
    bucket=)``. The answer holds the ``record`` (a product's ``recipe.json`` record, or the file's) and,
    for a product, its ``objects``; a refusal raises ``RuntimeError`` with the host's words."""
    path = os.environ.get(SEED_VARIABLE)
    if not path:
        pytest.fail(f"store: {STORE_NAME} names a store, but {SEED_VARIABLE} is unset")

    def ask(**request):
        with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as conn:
            conn.connect(path)
            conn.sendall((json.dumps(request) + "\n").encode("utf-8"))
            with conn.makefile("rb") as answer:
                reply = json.loads(answer.readline() or b"{}")
        if "error" in reply or "record" not in reply:
            raise RuntimeError(f"the seed channel refused {sorted(request)}: {reply.get('error', reply)}")
        return reply
    return ask
