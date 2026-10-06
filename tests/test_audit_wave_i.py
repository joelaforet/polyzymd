"""Regressions for Joe's #163 review: users add other files to project and study folders.

PolyzyMD reads, hashes and deposits only files whose names or folders it
chooses in advance.
"""

from __future__ import annotations

from pathlib import Path


class TestStrayFiles:
    def test_stray_files_do_not_change_the_code_hash(self, tmp_path: Path, caplog) -> None:
        import logging

        from polyzymd.analyses.timeseries import folder_hash

        analyses = tmp_path / "analyses"
        (analyses / "data").mkdir(parents=True)
        (analyses / "f.py").write_text("def f(u):\n    return 1.0\n")
        (analyses / "data" / "table.csv").write_text("1\n")
        before = folder_hash(analyses / "f.py")
        (analyses / "notes.txt").write_text("my notes")
        (analyses / "copy.xtc").write_bytes(b"\0" * 1000)
        with caplog.at_level(logging.INFO):
            assert folder_hash(analyses / "f.py") == before
        assert "copy.xtc, notes.txt" in caplog.text and "data/" in caplog.text
        (analyses / "data" / "table.csv").write_text("2\n")
        assert folder_hash(analyses / "f.py") != before
