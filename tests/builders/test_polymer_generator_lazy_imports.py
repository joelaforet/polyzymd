"""Subprocess checks for polymer generator lazy imports."""

from __future__ import annotations

import os
import subprocess
import sys
import textwrap
from pathlib import Path


def _run_import_guarded_script(script: str) -> subprocess.CompletedProcess[str]:
    """Run a Python snippet with the repository source path available.

    Parameters
    ----------
    script : str
        Python source code to run in a subprocess.

    Returns
    -------
    subprocess.CompletedProcess[str]
        Completed subprocess result for assertions.
    """
    repo_root = Path(__file__).resolve().parents[2]
    env = os.environ.copy()
    existing_pythonpath = env.get("PYTHONPATH")
    src_path = str(repo_root / "src")
    env["PYTHONPATH"] = (
        src_path if not existing_pythonpath else os.pathsep.join([src_path, existing_pythonpath])
    )

    return subprocess.run(
        [sys.executable, "-c", textwrap.dedent(script)],
        check=False,
        cwd=repo_root,
        env=env,
        capture_output=True,
        text=True,
    )


def test_polymer_generator_import_does_not_import_optional_chemistry_stack() -> None:
    """Importing polymer_generator should not import optional chemistry dependencies."""
    result = _run_import_guarded_script(
        r"""
        import importlib
        import sys

        blocked_prefixes = ("polymerist", "rdkit", "openff", "openmm")
        before = {name for name in sys.modules if name.startswith(blocked_prefixes)}
        module = importlib.import_module("polyzymd.builders.polymer_generator")
        if module.PolymerGenerator.__name__ != "PolymerGenerator":
            raise AssertionError("PolymerGenerator was not imported")
        imported = sorted(
            name for name in sys.modules if name.startswith(blocked_prefixes) and name not in before
        )
        if imported:
            raise AssertionError(f"Optional dependencies imported eagerly: {imported}")
        """,
    )

    assert result.returncode == 0, result.stderr


def test_unknown_label_validation_runs_before_optional_dependency_imports() -> None:
    """Unknown dynamic labels should fail before OpenFF or Polymerist imports are needed."""
    result = _run_import_guarded_script(
        r"""
        import builtins
        import importlib
        import tempfile
        from pathlib import Path

        blocked_prefixes = ("polymerist", "rdkit", "openff", "openmm")
        original_import = builtins.__import__

        def guarded_import(name, *args, **kwargs):
            if name.startswith(blocked_prefixes):
                raise RuntimeError(f"Unexpected optional dependency import: {name}")
            return original_import(name, *args, **kwargs)

        builtins.__import__ = guarded_import
        module = importlib.import_module("polyzymd.builders.polymer_generator")

        class FakeMonomerGroup:
            def __init__(self):
                self.monomers = {"EGMA_1-site": ["terminal"], "EGMA_2-site": ["middle"]}

        with tempfile.TemporaryDirectory() as tmpdir:
            generator = module.PolymerGenerator(
                monomer_group=FakeMonomerGroup(),
                cache_directory=Path(tmpdir),
            )
            try:
                generator.generate_polymer("AXA", {"A": "EGMA"})
            except ValueError as exc:
                if "No monomer name configured for sequence label: X" not in str(exc):
                    raise
            else:
                raise AssertionError("Unknown sequence label was accepted")
        """,
    )

    assert result.returncode == 0, result.stderr
