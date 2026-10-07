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
import warnings
from dataclasses import dataclass, field
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
    warnings : list of str
        Where the report and the stored values disagree: no report, a
        partial report, or replicates stored but not in the report (or the
        reverse), as when a later run stored values and its report job
        failed. Empty when the report covers exactly the stored values.
        Each is also issued as a Python warning when the results are read.
    """

    folder: Path
    report: Any
    table: Any
    warnings: list[str] = field(default_factory=list)

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


def _missing_rows(filled: list[tuple[dict, list[tuple]]], rows: list[tuple]) -> list[tuple]:
    """Return the rows of labels a replicate lacks, with the run's ``missing`` value.

    A function with ``labels="returned"`` stores only the labels each
    replicate returned. The report gives a label that some replicate lacks
    the value ``missing`` (an unformed hydrogen bond is 0), so the stored
    table does the same: each replicate whose record names ``missing`` gets a
    row for every label of that quantity and part found in any replicate.
    """
    labels: dict[tuple, set] = {}
    for name, _c, _r, part, label, *_ in rows:
        if label is not None:
            labels.setdefault((name, part), set()).add(label)
    extra = []
    for record, own in filled:
        missing = float(record["missing"])
        base = (record["name"], record["condition"], int(record["replicate"]))
        have = {(part, label) for _n, _c, _r, part, label, *_ in own}
        for (name, part), every in labels.items():
            if name != record["name"]:
                continue
            for label in sorted(every - {lab for p, lab in have if p == part}, key=str):
                extra.append((*base, part, label, None, None, missing, record.get("unit")))
    return extra


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
    filled: list[tuple[dict, list[tuple]]] = []
    for path in records:
        record = json.loads(path.read_text())
        own = _rows(record, path.parent)
        rows.extend(own)
        if record.get("missing") is not None:
            filled.append((record, own))
    rows.extend(_missing_rows(filled, rows))
    report_path = folder / REPORT_FILE
    report = (
        ProtocolReport.model_validate_json(report_path.read_text())
        if report_path.is_file()
        else None
    )
    table = pd.DataFrame(rows, columns=list(COLUMNS))
    notes = _report_mismatch(report, table)
    for note in notes:
        warnings.warn(f"{folder}: {note}", UserWarning, stacklevel=2)
    return StoredResults(folder, report, table, notes)


def _report_mismatch(report: Any, table: Any) -> list[str]:
    """Return how a run's report and its stored values disagree, one note per difference."""
    if report is None:
        return ["no report.json: the values are stored but no report was written for them"]
    notes = []
    if report.status != "complete":
        notes.append(f"the report is {report.status}: " + "; ".join(report.problems))
    stored = {(str(c), int(r)) for c, r in zip(table["condition"], table["replicate"])}
    covered = {(c.label, int(r)) for c in report.conditions for r in c.replicates}
    if not covered:
        # A report without replicate numbers is compared by condition only.
        stored = {(c, 0) for c, _ in stored}
        covered = {(c.label, 0) for c in report.conditions}
    for what, pairs in (
        ("stored but not in the report", stored - covered),
        ("in the report but not stored", covered - stored),
    ):
        by_condition: dict[str, list[int]] = {}
        for label, index in sorted(pairs):
            by_condition.setdefault(label, []).append(index)
        for label, indices in by_condition.items():
            which = "" if indices == [0] else f" replicates {', '.join(map(str, indices))}"
            notes.append(
                f"condition {label}{which}: {what}; the report may be stale or partial "
                "(rerun polyzymd analyze to rewrite it)"
            )
    return notes
