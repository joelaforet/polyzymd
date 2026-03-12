"""Per-module colored logging and CLI output for PolyzyMD.

This module provides:
- Automatic terminal color capability detection (truecolor / 256 / basic / none)
- A registry of per-module colors so each subsystem is visually distinct
- A ``ColoredFormatter`` that applies per-module tinting to INFO/DEBUG while
  keeping WARNING amber-yellow and ERROR red across all modules
- ``colored_echo()`` — a drop-in replacement for ``click.echo()`` that applies
  phase-aware coloring
- ``setup_colored_logging()`` — one-call initialization

Color output is automatically disabled when stdout is not a TTY, when the
``NO_COLOR`` environment variable is set (per https://no-color.org), or when
the ``--no-color`` CLI flag is passed.
"""

from __future__ import annotations

import logging
import os
import sys
from dataclasses import dataclass, field
from enum import Enum
from typing import IO

import click

# ---------------------------------------------------------------------------
# Terminal capability detection
# ---------------------------------------------------------------------------


class TerminalColorSupport(Enum):
    """Level of color support detected in the current terminal."""

    NONE = "none"  # No color (pipe, dumb terminal, NO_COLOR)
    BASIC = "basic"  # 16-color ANSI
    EXTENDED = "extended"  # 256-color xterm
    TRUECOLOR = "truecolor"  # 24-bit RGB


def detect_color_support(stream: IO[str] | None = None) -> TerminalColorSupport:
    """Detect the color capability of the given *stream* (default: stderr).

    The check respects ``NO_COLOR`` and ``TERM=dumb`` conventions.
    """
    if stream is None:
        stream = sys.stderr

    # NO_COLOR convention — any non-empty value disables color
    if os.environ.get("NO_COLOR", ""):
        return TerminalColorSupport.NONE

    # Must be a real TTY
    if not hasattr(stream, "isatty") or not stream.isatty():
        return TerminalColorSupport.NONE

    term = os.environ.get("TERM", "")
    if term == "dumb":
        return TerminalColorSupport.NONE

    colorterm = os.environ.get("COLORTERM", "").lower()
    if colorterm in ("truecolor", "24bit"):
        return TerminalColorSupport.TRUECOLOR

    # Many modern terminals set TERM=xterm-256color (or similar)
    if "256color" in term:
        return TerminalColorSupport.EXTENDED

    # Fallback — most terminals support at least basic 16 colors
    return TerminalColorSupport.BASIC


# Singleton — resolved once, then cached
_color_support: TerminalColorSupport | None = None


def get_color_support() -> TerminalColorSupport:
    """Return the cached color-support level (lazy-initialized)."""
    global _color_support  # noqa: PLW0603
    if _color_support is None:
        _color_support = detect_color_support()
    return _color_support


def set_color_support(level: TerminalColorSupport) -> None:
    """Override the detected color support (e.g. from ``--no-color``)."""
    global _color_support  # noqa: PLW0603
    _color_support = level


# ---------------------------------------------------------------------------
# Module color definitions
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ModuleColor:
    """Associates a logger-name prefix with an RGB color and fallback codes."""

    prefix: str
    rgb: tuple[int, int, int]
    xterm256: int  # Closest xterm-256 index
    basic_ansi: str  # 16-color ANSI escape for BASIC terminals (e.g. "\033[94m")


# Pre-computed nearest xterm-256 indices for each chosen RGB.
# (Verified with the 6x6x6 color cube formula.)
# Basic ANSI fallbacks use bright variants where possible for readability
# on dark terminal backgrounds (e.g. TERM=xterm-16color on HPC nodes).
_DEFAULT_MODULE_COLORS: list[ModuleColor] = [
    # CLI / Workflow — lavender → bright blue
    ModuleColor("polyzymd.cli", (180, 180, 215), 146, "\033[94m"),
    ModuleColor("polyzymd.workflow", (180, 180, 215), 146, "\033[94m"),
    # Building — sage green → bright green
    ModuleColor("polyzymd.builders", (175, 215, 175), 151, "\033[92m"),
    # Simulation (runner / continuation) — warm peach → bright magenta
    ModuleColor("polyzymd.simulation.runner", (215, 185, 175), 181, "\033[95m"),
    ModuleColor("polyzymd.simulation.continuation", (215, 185, 175), 181, "\033[95m"),
    # Progress / signals — steel blue → bright cyan
    ModuleColor("polyzymd.simulation.progress", (175, 200, 220), 152, "\033[96m"),
    ModuleColor("polyzymd.simulation.signals", (175, 200, 220), 152, "\033[96m"),
    # Core — parchment → dark yellow
    ModuleColor("polyzymd.core", (215, 210, 175), 187, "\033[33m"),
    # Data — aqua/mint → dark cyan
    ModuleColor("polyzymd.data", (175, 215, 205), 152, "\033[36m"),
    # Exporters — lilac → dark magenta
    ModuleColor("polyzymd.exporters", (205, 175, 215), 182, "\033[35m"),
    # Utils — neutral gray → white
    ModuleColor("polyzymd.utils", (195, 195, 195), 250, "\033[37m"),
]

# Level-override colors (apply regardless of module)
_LEVEL_COLORS: dict[int, tuple[int, int, int]] = {
    logging.WARNING: (255, 200, 50),  # Amber yellow
    logging.ERROR: (255, 80, 80),  # Red
    logging.CRITICAL: (255, 80, 80),  # Red (same as ERROR)
}

_LEVEL_XTERM256: dict[int, int] = {
    logging.WARNING: 220,  # Gold
    logging.ERROR: 196,  # Red
    logging.CRITICAL: 196,
}


# ---------------------------------------------------------------------------
# Color scheme registry
# ---------------------------------------------------------------------------


@dataclass
class ColorScheme:
    """Registry of module-to-color mappings.

    Colors are matched by *longest prefix* so that e.g.
    ``polyzymd.simulation.runner`` beats ``polyzymd.simulation``.
    """

    _module_colors: list[ModuleColor] = field(default_factory=list)
    _level_colors: dict[int, tuple[int, int, int]] = field(default_factory=dict)
    _level_xterm256: dict[int, int] = field(default_factory=dict)

    def register(
        self,
        prefix: str,
        rgb: tuple[int, int, int],
        xterm256: int,
        basic_ansi: str = "",
    ) -> None:
        """Add or replace a module color mapping."""
        # Remove any existing entry with the same prefix
        self._module_colors = [mc for mc in self._module_colors if mc.prefix != prefix]
        self._module_colors.append(ModuleColor(prefix, rgb, xterm256, basic_ansi))

    # -- colour resolution ---------------------------------------------------

    def _find_module_color(self, logger_name: str) -> ModuleColor | None:
        """Return the best-matching ``ModuleColor`` via longest-prefix match."""
        best: ModuleColor | None = None
        best_len = 0
        for mc in self._module_colors:
            if logger_name.startswith(mc.prefix) and len(mc.prefix) > best_len:
                best = mc
                best_len = len(mc.prefix)
        return best

    def get_ansi(self, logger_name: str, level: int) -> tuple[str, str]:
        """Return ``(open_seq, close_seq)`` ANSI codes for the given logger/level.

        Returns empty strings when color is disabled.
        """
        support = get_color_support()
        if support is TerminalColorSupport.NONE:
            return ("", "")

        reset = "\033[0m"

        # Level overrides first
        if level in self._level_colors:
            return (self._level_escape(level, support), reset)

        # Module color
        mc = self._find_module_color(logger_name)
        if mc is None:
            return ("", "")

        return (self._module_escape(mc, support), reset)

    # -- escape-sequence helpers ---------------------------------------------

    def _module_escape(self, mc: ModuleColor, support: TerminalColorSupport) -> str:
        if support is TerminalColorSupport.TRUECOLOR:
            r, g, b = mc.rgb
            return f"\033[38;2;{r};{g};{b}m"
        if support is TerminalColorSupport.EXTENDED:
            return f"\033[38;5;{mc.xterm256}m"
        # BASIC — use the 16-color ANSI fallback
        return mc.basic_ansi

    def _level_escape(self, level: int, support: TerminalColorSupport) -> str:
        if support is TerminalColorSupport.TRUECOLOR:
            r, g, b = self._level_colors[level]
            return f"\033[38;2;{r};{g};{b}m"
        if support is TerminalColorSupport.EXTENDED:
            idx = self._level_xterm256.get(level, 0)
            return f"\033[38;5;{idx}m"
        # BASIC fallback
        if level >= logging.ERROR:
            return "\033[31m"  # red
        if level >= logging.WARNING:
            return "\033[33m"  # yellow
        return ""


def default_scheme() -> ColorScheme:
    """Create a ``ColorScheme`` pre-loaded with the PolyzyMD palette."""
    scheme = ColorScheme(
        _module_colors=list(_DEFAULT_MODULE_COLORS),
        _level_colors=dict(_LEVEL_COLORS),
        _level_xterm256=dict(_LEVEL_XTERM256),
    )
    return scheme


# Singleton instance
_scheme: ColorScheme | None = None


def get_scheme() -> ColorScheme:
    """Return the global ``ColorScheme`` (lazy-created)."""
    global _scheme  # noqa: PLW0603
    if _scheme is None:
        _scheme = default_scheme()
    return _scheme


# ---------------------------------------------------------------------------
# Colored logging formatter
# ---------------------------------------------------------------------------


class ColoredFormatter(logging.Formatter):
    """Logging formatter that applies per-module ANSI colors.

    WARNING and ERROR messages override the module color with amber-yellow
    and red respectively.  DEBUG and INFO use the module's assigned tint.
    """

    def __init__(
        self,
        fmt: str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt: str | None = None,
        scheme: ColorScheme | None = None,
    ) -> None:
        super().__init__(fmt, datefmt)
        self._scheme = scheme or get_scheme()

    def format(self, record: logging.LogRecord) -> str:
        formatted = super().format(record)
        open_seq, close_seq = self._scheme.get_ansi(record.name, record.levelno)
        if open_seq:
            return f"{open_seq}{formatted}{close_seq}"
        return formatted


# ---------------------------------------------------------------------------
# colored_echo — click.echo replacement
# ---------------------------------------------------------------------------

# Phase-to-prefix mapping for colored_echo
_PHASE_PREFIXES: dict[str, str] = {
    "cli": "polyzymd.cli",
    "workflow": "polyzymd.workflow",
    "build": "polyzymd.builders",
    "simulation": "polyzymd.simulation.runner",
    "progress": "polyzymd.simulation.progress",
    "core": "polyzymd.core",
    "data": "polyzymd.data",
    "export": "polyzymd.exporters",
    "utils": "polyzymd.utils",
}


def colored_echo(
    message: object = "",
    *,
    phase: str = "cli",
    err: bool = False,
    nl: bool = True,
    level: int = logging.INFO,
) -> None:
    """Write *message* to the terminal with per-phase coloring.

    This is a thin wrapper around ``click.echo`` that applies the same
    per-module color scheme used by the logging formatter.

    Parameters
    ----------
    message : str
        Text to display.
    phase : str
        One of the keys in ``_PHASE_PREFIXES`` (e.g. ``"build"``,
        ``"simulation"``).  Determines the color.
    err : bool
        Write to stderr instead of stdout (passed to ``click.echo``).
    nl : bool
        Append a newline (passed to ``click.echo``).
    level : int
        Logging level — controls whether the level-override color is used
        (e.g. ``logging.WARNING`` → amber).
    """
    scheme = get_scheme()
    prefix = _PHASE_PREFIXES.get(phase, "polyzymd.cli")
    open_seq, close_seq = scheme.get_ansi(prefix, level)

    text = str(message) if message is not None else ""
    if open_seq:
        text = f"{open_seq}{text}{close_seq}"

    click.echo(text, err=err, nl=nl)


# ---------------------------------------------------------------------------
# One-call setup
# ---------------------------------------------------------------------------


def setup_colored_logging(
    *,
    verbose: bool = False,
    no_color: bool = False,
) -> None:
    """Configure the root logger with the colored formatter.

    Parameters
    ----------
    verbose : bool
        If True, set the root logger level to DEBUG.
    no_color : bool
        If True, forcibly disable all color output.
    """
    if no_color:
        set_color_support(TerminalColorSupport.NONE)

    level = logging.DEBUG if verbose else logging.INFO

    # Remove any existing handlers on the root logger so we don't double-print
    root = logging.getLogger()
    for handler in root.handlers[:]:
        root.removeHandler(handler)

    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(level)
    handler.setFormatter(ColoredFormatter())

    root.setLevel(level)
    root.addHandler(handler)
