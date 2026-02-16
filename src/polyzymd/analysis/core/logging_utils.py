"""Logging utilities for colorized terminal output.

This module provides a ColoredFormatter and setup function for consistent
logging with visual emphasis on warnings and errors in terminal output.
"""

from __future__ import annotations

import logging
import sys


class ColoredFormatter(logging.Formatter):
    """Formatter that adds ANSI color codes for WARNING and ERROR levels.

    Colors are only applied when output is to an interactive terminal (TTY).
    When redirecting to a file or pipe, plain text is used.

    Attributes
    ----------
    COLORS : dict
        Mapping of log levels to ANSI color codes.
    RESET : str
        ANSI code to reset text formatting.
    """

    COLORS = {
        logging.WARNING: "\033[93m",  # Yellow
        logging.ERROR: "\033[91m",  # Red
        logging.CRITICAL: "\033[91m",  # Red (bold could be added)
    }
    RESET = "\033[0m"

    def format(self, record: logging.LogRecord) -> str:
        """Format the log record with color if appropriate.

        Parameters
        ----------
        record : logging.LogRecord
            The log record to format.

        Returns
        -------
        str
            Formatted message, with ANSI color codes if outputting to TTY.
        """
        message = super().format(record)
        color = self.COLORS.get(record.levelno)
        if color and sys.stderr.isatty():
            return f"{color}{message}{self.RESET}"
        return message


def setup_logging(verbose: bool = False) -> None:
    """Set up logging with colored output for warnings and errors.

    This function configures the root logger with a ColoredFormatter
    that highlights WARNING and ERROR messages in yellow and red
    respectively when outputting to a terminal.

    Parameters
    ----------
    verbose : bool, optional
        If True, set log level to INFO. If False (default), set to WARNING.

    Examples
    --------
    >>> from polyzymd.analysis.core.logging_utils import setup_logging
    >>> setup_logging(verbose=True)  # Show INFO and above with colors
    >>> setup_logging(verbose=False)  # Show WARNING and above with colors
    """
    level = logging.INFO if verbose else logging.WARNING

    # Create handler with colored formatter
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(ColoredFormatter("%(message)s"))

    # Configure root logger
    logging.root.handlers = []
    logging.root.addHandler(handler)
    logging.root.setLevel(level)


def suppress_mdanalysis_info() -> None:
    """Suppress verbose MDAnalysis INFO-level log messages.

    MDAnalysis emits many INFO messages during Universe creation that
    clutter terminal output without providing actionable information:
    - "Setting segids from chainIDs..."
    - "The attribute(s) types have already been read..."
    - "attribute masses has been guessed successfully"

    This function sets the MDAnalysis logger to WARNING level to
    suppress these messages while preserving important warnings.
    """
    logging.getLogger("MDAnalysis").setLevel(logging.WARNING)
