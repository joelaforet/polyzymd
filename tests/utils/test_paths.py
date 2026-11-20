"""Tests for printing and scripting paths without ``..`` parts."""

from __future__ import annotations

import os

from polyzymd.utils.paths import names_existing, shown_path, tidy_path


def test_dot_dot_parts_are_resolved(tmp_path):
    (tmp_path / "conditions" / "water").mkdir(parents=True)
    path = tmp_path / "conditions" / "water" / ".." / ".." / "runs" / "run_1"

    assert tidy_path(path) == (tmp_path / "runs" / "run_1").resolve()
    assert ".." not in shown_path(path)


def test_dot_dot_follows_a_symlinked_folder_as_the_system_does(tmp_path):
    """Behind a symlink, '..' is the parent of the target, not of the link."""
    (tmp_path / "real" / "deep").mkdir(parents=True)
    (tmp_path / "study").mkdir()
    os.symlink(tmp_path / "real" / "deep", tmp_path / "study" / "link")

    assert tidy_path(tmp_path / "study" / "link" / ".." / "runs") == (
        (tmp_path / "real" / "runs").resolve()
    )


def test_a_path_without_dot_dot_keeps_its_symlinks(tmp_path):
    (tmp_path / "real").mkdir()
    os.symlink(tmp_path / "real", tmp_path / "scratch")

    assert tidy_path(tmp_path / "scratch" / "run_1") == tmp_path / "scratch" / "run_1"


def test_names_existing_finds_files_and_folders(tmp_path):
    (tmp_path / "ref.pdb").write_text("END\n")

    assert names_existing(tmp_path / "ref.pdb", file=True)
    assert names_existing(tmp_path)
    assert not names_existing(tmp_path, file=True)
    assert not names_existing(tmp_path / "missing.pdb")


def test_text_that_cannot_be_a_path_names_nothing(tmp_path):
    """A long selection raised 'File name too long' when probed as a file."""
    selection = "protein and name CA and (resid " + " ".join(map(str, range(1, 112))) + ")"
    assert len(selection.encode()) > 255

    for text in (selection, tmp_path / selection, "ref\0.pdb"):
        assert not names_existing(text)
        assert not names_existing(text, file=True)
