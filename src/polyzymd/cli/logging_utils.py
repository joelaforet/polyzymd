"""Logging set-up for the analysis commands.

:func:`analysis_logging` writes every record to a log file and shows only
warnings and errors on the console unless ``verbose`` is set.
"""

from __future__ import annotations

import logging
import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

# Format of each line in the analysis log file
FULL_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

#: Loggers whose records go to the log file only, unless they are errors:
#: Python warnings captured from libraries (deprecations, missing unit cells)
#: and pymbar's notices, such as the banner it prints when JAX is absent.
FILE_ONLY_LOGGERS = ("py.warnings", "pymbar")


class _ConsoleFilter(logging.Filter):
    """Keep library chatter out of the console, and print each warning once.

    Library records below ERROR go to the log file only. A message already
    printed in this command is not printed again (loaders reached twice for
    one replicate would otherwise repeat theirs); the log file keeps every one.
    """

    def __init__(self) -> None:
        super().__init__()
        self._seen: set[tuple[str, int, str]] = set()

    def filter(self, record: logging.LogRecord) -> bool:
        name = record.name
        if any(name == prefix or name.startswith(prefix + ".") for prefix in FILE_ONLY_LOGGERS):
            return record.levelno >= logging.ERROR
        key = (name, record.levelno, record.getMessage())
        if key in self._seen:
            return False
        self._seen.add(key)
        return True


def analysis_logging(log_dir: "Path", command: str, verbose: bool = False) -> "Path":
    """Send the full log of an analysis command to a file and keep the console quiet.

    The console shows WARNING and above, without Python warnings captured from
    libraries or pymbar's notices; with ``verbose`` it keeps INFO. Everything,
    at INFO and above, with captured Python warnings, goes to
    ``log_dir/polyzymd-<command>-<time>.log``, which is returned so the command
    can print where it is. A quiet console keeps what an agent must read short.
    """
    from datetime import datetime
    from pathlib import Path

    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    # Array tasks can start in the same second; the task ID or process ID
    # keeps their logs apart.
    task = os.environ.get("SLURM_ARRAY_TASK_ID")
    tag = f"task{task}" if task else f"pid{os.getpid()}"
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    path = log_dir / f"polyzymd-{command}-{stamp}-{tag}.log"
    root = logging.getLogger()
    for handler in root.handlers:
        if isinstance(handler, logging.StreamHandler) and not isinstance(
            handler, logging.FileHandler
        ):
            handler.setLevel(logging.INFO if verbose else logging.WARNING)
            if not any(isinstance(f, _ConsoleFilter) for f in handler.filters):
                handler.addFilter(_ConsoleFilter())
    file_handler = logging.FileHandler(path)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter(FULL_FORMAT))
    root.addHandler(file_handler)
    if root.level > logging.INFO or root.level == logging.NOTSET:
        root.setLevel(logging.INFO)
    logging.captureWarnings(True)
    return path
