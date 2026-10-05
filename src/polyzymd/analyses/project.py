"""``Project``: the studies of one paper, one per protein, read together.

:class:`Project` reads ``project.yaml`` (:mod:`polyzymd.analyses.project_file`)
and gives each study as a :class:`~polyzymd.analyses.study.Study`. Its
:meth:`Project.results` puts every study's stored results of one analysis in
one table with a ``study`` column, without loading any trajectory.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator

from polyzymd.analyses.exceptions import ProtocolError


@dataclass
class ProjectResults:
    """One analysis's stored results in every study that runs it.

    Attributes
    ----------
    table : pandas.DataFrame
        Every study's :attr:`StoredResults.table <polyzymd.analyses.results.StoredResults.table>`
        with a ``study`` column first, and a column for each factor any
        study's conditions declare (empty where a condition has none).
    reports : dict of str to ProtocolReport or None
        Each study's report, by study label. Comparisons stay within a study.
    folders : dict of str to Path
        Each study's results folder of the analysis.
    """

    table: Any
    reports: dict[str, Any]
    folders: dict[str, Path]


class Project:
    """The studies listed in a ``project.yaml``, one per protein.

    Parameters
    ----------
    path : str or Path
        The project folder or its ``project.yaml``.
    """

    def __init__(self, path: str | Path) -> None:
        from polyzymd.analyses.project_file import load_project_file

        self.protocol = load_project_file(path)
        self._studies: dict[str, Any] = {}

    @property
    def root(self) -> Path:
        """The project folder."""
        return self.protocol.root

    @property
    def labels(self) -> list[str]:
        """The study labels, in the order ``project.yaml`` lists them."""
        return list(self.protocol.studies)

    def __getitem__(self, label: str) -> Any:
        """Return the study ``label`` as a :class:`~polyzymd.analyses.study.Study`."""
        from polyzymd.analyses.study import Study

        if label not in self.protocol.studies:
            raise ProtocolError(
                f"{self.protocol.path} lists no study {label!r}.",
                hint=f"Use one of {', '.join(self.labels)}.",
            )
        if label not in self._studies:
            self._studies[label] = Study(self.protocol.studies[label])
        return self._studies[label]

    def __iter__(self) -> Iterator[Any]:
        return (self[label] for label in self.labels)

    def __len__(self) -> int:
        return len(self.protocol.studies)

    def __repr__(self) -> str:
        return f"Project({self.root}, studies={self.labels})"

    def runs_in(self, run: str) -> list[str]:
        """Return the labels of the studies that run the analysis ``run``."""
        return [label for label in self.labels if run in self[label].protocol.analyses]

    def results(self, run: str) -> ProjectResults:
        """Return the stored results of ``run`` in every study that runs it, in one table.

        Raises
        ------
        ProtocolError
            If no study runs ``run``, or a study that runs it has no stored
            results; the message names those studies, so no study is left out
            silently.
        """
        import pandas as pd

        labels = self.runs_in(run)
        if not labels:
            raise ProtocolError(
                f"No study of {self.protocol.path} runs {run!r}.",
                hint="Add it under analyses: in project.yaml or in a study.yaml.",
            )
        missing, stored = [], {}
        for label in labels:
            try:
                stored[label] = self[label].results(run)
            except ProtocolError:
                missing.append(label)
        if missing:
            raise ProtocolError(
                f"{run}: studies {', '.join(missing)} have no stored results.",
                hint=f"Run polyzymd analyze {run} --project {self.root} first.",
            )
        tables = [result.table.assign(study=label) for label, result in stored.items()]
        table = pd.concat(tables, ignore_index=True, sort=False)
        table = table[["study", *[c for c in table.columns if c != "study"]]]
        return ProjectResults(
            table=table,
            reports={label: result.report for label, result in stored.items()},
            folders={label: result.folder for label, result in stored.items()},
        )
