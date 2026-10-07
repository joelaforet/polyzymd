"""Settings for the whole test suite."""

from __future__ import annotations

import logging

import pytest


@pytest.fixture(autouse=True, scope="session")
def _isolated_hash_cache(tmp_path_factory: pytest.TempPathFactory) -> None:
    """Keep the file-hash cache of stored results out of the user's home during tests."""
    mp = pytest.MonkeyPatch()
    mp.setenv("POLYZYMD_CACHE_DIR", str(tmp_path_factory.mktemp("polyzymd_cache")))
    yield
    mp.undo()


@pytest.fixture(autouse=True)
def _restore_root_logging() -> None:
    """Give back the root logger's handlers and level after each test.

    The ``polyzymd`` command replaces the root handlers, which removes the
    file handler pytest keeps for the whole session. Without it, a later
    ``polyzymd analyze`` writes a log file into the current folder.
    """
    root = logging.getLogger()
    handlers, level = list(root.handlers), root.level
    yield
    root.handlers, root.level = handlers, level


@pytest.fixture()
def git_identity(monkeypatch: pytest.MonkeyPatch) -> None:
    """Give git a committer name and email, for tests that commit or tag."""
    for key, value in {
        "GIT_AUTHOR_NAME": "Test",
        "GIT_AUTHOR_EMAIL": "test@example.com",
        "GIT_COMMITTER_NAME": "Test",
        "GIT_COMMITTER_EMAIL": "test@example.com",
    }.items():
        monkeypatch.setenv(key, value)
