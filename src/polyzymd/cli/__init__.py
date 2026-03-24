"""Command-line interface for PolyzyMD."""

from __future__ import annotations

from typing import Any

__all__ = ["cli", "main"]


def __getattr__(name: str) -> Any:
    """Lazily expose CLI entry points."""
    if name in {"cli", "main"}:
        from polyzymd.cli.main import cli, main

        return {"cli": cli, "main": main}[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
