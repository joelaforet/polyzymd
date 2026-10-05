"""Read stored analysis results back without loading any trajectory.

``polyzymd analyze RUN --study study.yaml`` stores one run's results in
``<study>/results/<RUN>/``: each measured quantity under
``polyzymd_results/<name>/<condition>/replicate_<n>/`` (its ``record.json``
and ``series.npz`` or ``values.npz``), the report as ``report.json``, and the
figures under ``figures/``. :func:`read_results` reads them, so figure
scripts and notebooks can redraw a paper's figures from a study folder whose
trajectories are not on the machine.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from polyzymd.analyses.exceptions import ProtocolError

#: Name of the saved report in a run's results folder.
REPORT_FILE = "report.json"
#: Columns of :attr:`StoredResults.table`.
COLUMNS = ("name", "condition", "replicate", "part", "label", "frame", "time_ns", "value", "unit")


@dataclass
class StoredResults:
    """One analysis run's stored values and report.

    Attributes
    ----------
    folder : Path
        The run's results folder.
    report : ProtocolReport or None
        The report ``polyzymd analyze`` printed, or ``None`` when it saved none.
    table : pandas.DataFrame
        One row per stored value, with the columns :data:`COLUMNS`. ``name``
        is the measured quantity, as its folder under ``polyzymd_results/`` is
        named. A per-frame series has ``frame`` and ``time_ns`` and no
        ``label``; a per-replicate value has neither, and a labelled one (for
        example one value per residue) has its ``label`` and, when the
        function returned several parts, its ``part``.
    """

    folder: Path
    report: Any
    table: Any

    @property
    def names(self) -> list[str]:
        """The measured quantities stored for the run."""
        return sorted(set(self.table["name"]))


def _rows(record: dict[str, Any], folder: Path) -> list[tuple]:
    """Return the table rows of one replicate's stored values."""
    import numpy as np

    base = (record["name"], record["condition"], int(record["replicate"]))
    unit = record.get("unit")
    series = folder / "series.npz"
    if series.is_file():
        with np.load(series) as data:
            values, frames, times = data["values"], data["frames"], data["times"]
        if values.ndim == 2:
            # Several quantities per frame, one column per part.
            parts = json.loads((folder / "parts.json").read_text())
            return [
                (*base, part, None, int(f), float(t), float(v), unit)
                for i, part in enumerate(parts)
                for f, t, v in zip(frames, times, values[:, i], strict=True)
            ]
        return [
            (*base, None, None, int(f), float(t), float(v), unit)
            for f, t, v in zip(frames, times, values, strict=True)
        ]
    with np.load(folder / "values.npz") as data:
        values = np.asarray(data["values"], dtype=float)
    labels = record.get("labels")
    if labels == "returned":
        labels = json.loads((folder / "labels.json").read_text())
    if values.ndim == 0:
        return [(*base, None, None, None, None, float(values), unit)]
    parts_file = folder / "parts.json"
    if labels is None:
        # One value per part, with no labels.
        parts = json.loads(parts_file.read_text()) if parts_file.is_file() else range(len(values))
        return [(*base, p, None, None, None, float(v), unit) for p, v in zip(parts, values)]
    if values.ndim == 1:
        return [(*base, None, lab, None, None, float(v), unit) for lab, v in zip(labels, values)]
    parts = json.loads(parts_file.read_text()) if parts_file.is_file() else range(len(values))
    return [
        (*base, part, lab, None, None, float(v), unit)
        for part, row in zip(parts, values)
        for lab, v in zip(labels, row)
    ]


def read_results(folder: str | Path) -> StoredResults:
    """Read every stored value and the report of one analysis run's results folder.

    Raises
    ------
    ProtocolError
        If the folder holds no stored results.
    """
    import pandas as pd

    from polyzymd.analyses.protocols import ProtocolReport

    folder = Path(folder)
    records = sorted(folder.glob("polyzymd_results/*/*/replicate_*/record.json"))
    if not records:
        raise ProtocolError(
            f"No stored results under {folder}.",
            hint="Run polyzymd analyze RUN --study study.yaml first, on a machine with the "
            "trajectories.",
        )
    rows: list[tuple] = []
    for path in records:
        rows.extend(_rows(json.loads(path.read_text()), path.parent))
    report_path = folder / REPORT_FILE
    report = (
        ProtocolReport.model_validate_json(report_path.read_text())
        if report_path.is_file()
        else None
    )
    return StoredResults(folder, report, pd.DataFrame(rows, columns=list(COLUMNS)))
