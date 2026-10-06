"""``Project``: the studies of one paper, one per protein, read together.

:class:`Project` reads ``project.yaml`` (:mod:`polyzymd.analyses.project_file`)
and gives each study as a :class:`~polyzymd.analyses.study.Study`.
:meth:`Project.results` puts every study's stored results of one analysis in
one table with a ``study`` column (a :class:`ProjectResults`), and
:meth:`Project.replicate_table` does the same for the per-replicate values.
Neither loads any trajectory.
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

    The file is read and checked when the object is made. Each study is
    opened as a :class:`~polyzymd.analyses.study.Study` the first time it is
    accessed and kept for later accesses. Iterating gives the studies in
    the order ``project.yaml`` lists them; ``len`` gives their number.

    Parameters
    ----------
    path : str or Path
        The project folder or its ``project.yaml``.

    Attributes
    ----------
    protocol : ProjectFile
        The checked contents of ``project.yaml``
        (:class:`~polyzymd.analyses.project_file.ProjectFile`).

    Raises
    ------
    ProtocolError
        If ``project.yaml`` cannot be read or fails its checks.
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
        """Return the study ``label`` as a :class:`~polyzymd.analyses.study.Study`.

        Parameters
        ----------
        label : str
            A study label from ``project.yaml``.

        Returns
        -------
        Study
            The study, opened on first access and cached.

        Raises
        ------
        ProtocolError
            If the project lists no study ``label``.
        """
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
        """Iterate over the studies in the order ``project.yaml`` lists them."""
        return (self[label] for label in self.labels)

    def __len__(self) -> int:
        """Return the number of studies."""
        return len(self.protocol.studies)

    def __repr__(self) -> str:
        """Return ``Project(<root>, studies=[...])``."""
        return f"Project({self.root}, studies={self.labels})"

    def runs_in(self, run: str) -> list[str]:
        """Return the labels of the studies that run the analysis ``run``.

        A study runs ``run`` when its analyses, those of its ``study.yaml``
        together with the project analyses that apply to it, include
        ``run``. Every study is opened to check.

        Parameters
        ----------
        run : str
            Name of the analysis run.

        Returns
        -------
        list of str
            The labels, in the order ``project.yaml`` lists the studies.
        """
        return [label for label in self.labels if run in self[label].protocol.analyses]

    def replicate_table(self, run: str) -> Any:
        """Return one row per replicate of ``run`` in every study that runs it.

        Concatenates each study's :meth:`Study.replicate_table
        <polyzymd.analyses.study.Study.replicate_table>`, each with its
        ``study`` column, in the order ``project.yaml`` lists the studies.
        No trajectory is loaded.

        Parameters
        ----------
        run : str
            Name of the analysis run.

        Returns
        -------
        pandas.DataFrame
            One row per replicate (and per label or part) of every study
            that runs ``run``; factor columns a study does not declare are
            empty for its rows.

        Raises
        ------
        ProtocolError
            If no study runs ``run``, or a study that runs it has no stored
            results; the message names those studies.
        """
        import pandas as pd

        self.results(run)  # names any study without stored results
        return pd.concat(
            [self[label].replicate_table(run) for label in self.runs_in(run)],
            ignore_index=True,
            sort=False,
        )

    def results(self, run: str) -> ProjectResults:
        """Return the stored results of ``run`` in every study that runs it, in one table.

        Reads each study's stored results with
        :meth:`Study.results <polyzymd.analyses.study.Study.results>`; no
        trajectory is loaded. The tables are concatenated with a ``study``
        column first.

        Parameters
        ----------
        run : str
            Name of the analysis run.

        Returns
        -------
        ProjectResults
            The combined table, and each study's report and results folder by
            study label.

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
