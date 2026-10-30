"""The builder environments install openff-pablo from its pinned tag."""

import importlib.metadata as metadata
import json

import pytest

PABLO_COMMIT = "fb364e9a9c0bdbde1262d1e822fe13d823033f40"


def test_openff_pablo_is_installed_from_the_pinned_tag() -> None:
    from openff.pablo import topology_from_pdb  # noqa: F401

    direct = metadata.distribution("openff-pablo").read_text("direct_url.json")
    if direct is None:
        pytest.skip("openff-pablo is not installed from git")
    vcs = json.loads(direct).get("vcs_info", {})
    assert vcs.get("commit_id") == PABLO_COMMIT
    assert vcs.get("requested_revision") == "v0.2.2"
