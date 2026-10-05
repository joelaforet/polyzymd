"""One-call, self-describing analysis protocol for scripts and coding agents.

:func:`analyze` takes the name of one of the analyses in
:data:`FUNCTION_ANALYSES` and one or more simulation config paths, builds a
:class:`~polyzymd.analyses.study.Study` of them, measures the analysis with
:meth:`~polyzymd.analyses.study.Study.timeseries` or
:meth:`~polyzymd.analyses.study.Study.per_replicate`, and returns a
:class:`ProtocolReport` in which every number states what it is. Every field
is described in ``docs/source/reference/analysis_protocol_report.md``.

The replicate is the sampling unit throughout: every mean, interval and test
uses replicates, never frames, as its sample.

References
----------
.. [1] Grossfield, A.; Patrone, P. N.; Roe, D. R.; Schultz, A. J.;
       Siderius, D. W.; Zuckerman, D. M. Best Practices for Quantification of
       Uncertainty and Sampling Quality in Molecular Simulations.
       Living J. Comput. Mol. Sci. 2018, 1 (1), 5067.
       https://doi.org/10.33011/livecoms.1.1.5067
.. [2] Welch, B. L. The Generalization of "Student's" Problem when Several
       Different Population Variances are Involved. Biometrika 1947, 34 (1-2),
       28-35. https://doi.org/10.1093/biomet/34.1-2.28
.. [3] Benjamini, Y.; Hochberg, Y. Controlling the False Discovery Rate: A
       Practical and Powerful Approach to Multiple Testing. J. R. Stat. Soc.
       Ser. B 1995, 57 (1), 289-300.
       https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Sequence

from pydantic import BaseModel, ConfigDict, Field

from polyzymd.analyses.exceptions import ProtocolError

#: Published page of the polyzymd analyze protocol.
ANALYZE_PROTOCOL_URL = (
    "https://polyzymd.readthedocs.io/en/latest/how_to/analysis_agent_protocol.html"
)

#: Agent skill that teaches the polyzymd analyze protocol.
ANALYZE_AGENT_SKILL = ".claude/skills/polyzymd-analyze/SKILL.md"

#: Published page on writing an analysis as a function for the study API.
ANALYSIS_API_URL = "https://polyzymd.readthedocs.io/en/latest/explanation/analysis_api.html"

#: Published page of the catalytic triad routine on the study API.
TRIAD_ROUTINE_URL = (
    "https://polyzymd.readthedocs.io/en/latest/how_to/analysis_triad_quickstart.html"
)

#: Where a retirement message sends a reader, and what to point an agent at.
RETIRED_DOCS_POINTER = (
    f"Read {ANALYZE_PROTOCOL_URL}, or point an agent at {ANALYZE_AGENT_SKILL} or that page "
    "to learn the protocol."
)

# Verdict vocabulary. Kept small so a caller can branch on it without parsing
# the rest of the sentence.
VERDICT_LARGER = "larger"
VERDICT_SMALLER = "smaller"
VERDICT_CHANGED = "changed"
VERDICT_NO_DIFFERENCE = "no significant difference"
VERDICT_NO_TEST = "no test recorded"
VERDICT_NOT_TESTABLE = "not testable"
VERDICT_VOCABULARY = (
    VERDICT_LARGER,
    VERDICT_SMALLER,
    VERDICT_CHANGED,
    VERDICT_NO_DIFFERENCE,
    VERDICT_NO_TEST,
    VERDICT_NOT_TESTABLE,
)

#: The analyses polyzymd analyze runs, each through Study.timeseries or
#: Study.per_replicate, with the settings each one takes and their defaults.
FUNCTION_ANALYSES = {
    "rg": {"selection": "protein"},
    "rmsd": {
        "selection": "protein and name CA",
        "alignment_selection": "protein and name CA",
        "reference_mode": "centroid",
        "reference_frame": 1,
        "reference_file": None,
    },
    "rmsf": {
        "selection": "protein and name CA",
        "alignment_selection": "protein and name CA",
        "reference_mode": "centroid",
        "reference_frame": 1,
        "reference_file": None,
        "highlight_residues": [],
        "core": None,
        "regions": {},
    },
    "rmsd_per_residue": {
        "selection": "protein and name CA",
        "alignment_selection": "protein and name CA",
        "reference_mode": None,
        "reference_frame": 1,
        "reference_file": None,
        "highlight_residues": [],
        "core": None,
        "regions": {},
    },
    "sasa": {
        "target": "protein",
        "contexts": {},
        "probe_radius_nm": 0.14,
        "n_sphere_points": 960,
    },
    "secondary_structure": {"selection": "protein", "scheme": "simplified"},
    "hydrogen_bonds": {
        "groups": {"protein": "chainid A", "polymer": "chainid C"},
        "summaries": {"protein_polymer": {"between": ["protein", "polymer"]}},
        "d_a_cutoff": 3.5,
        "d_h_a_angle_cutoff": 150.0,
        "donors": None,
        "hydrogens": None,
        "acceptors": None,
        "lifetime_key": "residue",
        "tolerance_ps": 0.0,
    },
    "native_contacts": {
        "selection": "protein and not element H",
        "reference_mode": None,
        "reference_frame": 1,
        "reference_file": None,
        "radius": 4.5,
        "min_separation": 3,
        "beta": 5.0,
        "lambda_constant": 1.8,
        "use_pbc": True,
        "regions": {},
    },
    "contacts": {
        "method": "occlusion",
        "polymer_selection": "chainid C",
        "protein_selection": "chainid A",
        "polymer_types": None,
        "use_pbc": True,
        "regions": {},
        "cutoff": 4.0,
        "heavy_atoms": True,
        "exposed_threshold": 0.2,
        "buried_threshold": 0.2,
        "max_asa": "theoretical",
        "probe_radius_nm": 0.14,
        "n_sphere_points": 960,
        "tolerance_ps": 0.0,
    },
    "distances": {"pairs": None, "threshold": 3.5, "use_pbc": True},
}


__all__ = [
    "FUNCTION_ANALYSES",
    "VERDICT_VOCABULARY",
    "ConditionReport",
    "PairwiseReport",
    "ProtocolProvenance",
    "ProtocolReport",
    "analyze",
]


# Report models. Field meanings live in reference/analysis_protocol_report.md,


class ConditionReport(BaseModel):
    """One condition's mean with the uncertainty and sample size behind it.

    ``ci_method`` is ``"student_t"`` from replicate values,
    ``"student_t_from_sem"`` when it was rebuilt from a stored standard error,
    and ``"not_estimable"`` when every replicate has the same value.
    ``replicates``, ``statistical_inefficiency`` and ``n_effective`` list, for
    each entry of ``replicate_values``, its replicate number and the pymbar
    statistical inefficiency and effective sample size of its time series.
    ``eq_detected_frame`` and ``eq_detected_ns`` give, for each, the start of
    the equilibrated region that pymbar ``detect_equilibration`` finds in the
    production series, as a production frame index from 0 and as simulation
    time. They are diagnostics and change no value. ``entry`` is the label of this row in a
    labelled result, such as a residue ID, and ``None`` otherwise.
    """

    model_config = ConfigDict(ser_json_inf_nan="strings")

    label: str
    entry: str | None = None
    n_replicates: int
    mean: float
    sem: float | None = None
    ci95: tuple[float, float] | None = None
    ci_method: str | None = None
    replicate_values: list[float] = Field(default_factory=list)
    replicates: list[int] = Field(default_factory=list)
    statistical_inefficiency: list[float] = Field(default_factory=list)
    n_effective: list[float] = Field(default_factory=list)
    eq_detected_frame: list[int] = Field(default_factory=list)
    eq_detected_ns: list[float] = Field(default_factory=list)


class PairwiseReport(BaseModel):
    """One comparison of the primary metric, control against one condition.

    ``delta`` is ``mean(b) - mean(a)`` and ``cohens_d`` is oriented to match it.
    ``p_adjusted`` of ``None`` means no corrected p value was computed, so
    the row describes a difference rather than deciding it; ``testable`` of
    ``False`` means a condition has fewer than two replicates. ``family_size``
    is the number of tests in the Benjamini-Hochberg family this row was
    corrected in, one family per outcome, and ``None`` when that is not known
    or the row was not tested. ``entry`` is the label compared in a labelled
    result, such as a residue ID, and ``None`` otherwise.
    """

    model_config = ConfigDict(ser_json_inf_nan="strings")

    a: str
    b: str
    entry: str | None = None
    delta: float
    delta_ci95: tuple[float, float] | None = None
    p: float | None = None
    p_adjusted: float | None = None
    test: str = "student_t"
    correction: str = "BH"
    family_size: int | None = None
    cohens_d: float | None = None
    hedges_g: float | None = None
    direction: str = "unchanged"
    significant: bool = False
    testable: bool = True


class TrendReport(BaseModel):
    """The slope of the replicate values against one numeric factor of the conditions.

    ``conditions`` are those that declare the factor; ``n_replicates`` the
    points fitted, one per replicate. ``slope`` is in the metric's unit per
    unit of the factor, with a 95 percent t interval; ``p`` tests zero slope
    and ``p_adjusted`` corrects it over the study's numeric factors
    (Benjamini-Hochberg). ``testable`` is ``False`` with fewer than two
    factor levels, three replicates, or values that all agree.
    """

    model_config = ConfigDict(ser_json_inf_nan="strings")

    factor: str
    conditions: list[str] = Field(default_factory=list)
    n_replicates: int = 0
    slope: float | None = None
    slope_ci95: tuple[float, float] | None = None
    p: float | None = None
    p_adjusted: float | None = None
    family_size: int | None = None
    r_squared: float | None = None
    significant: bool = False
    testable: bool = False


class ProtocolProvenance(BaseModel):
    """Versions, config hashes, output paths and settings of one protocol run.

    ``settings`` holds the analysis settings the analysis ran with and what
    they resolved to, such as the residues of an rmsf core. ``study`` is set
    for a run from a study file: its ``path``, ``sha256``, ``run`` and the
    ``settings`` the run was given (from the file and ``--set``), and
    ``git``, the study folder's commit and uncommitted files
    (:func:`~polyzymd.analyses.study_git.git_state`), or ``None`` outside a
    repository.
    """

    polyzymd_version: str
    mdanalysis_version: str | None = None
    config_hashes: dict[str, str] = Field(default_factory=dict)
    settings_fingerprint: str | None = None
    settings: dict[str, Any] = Field(default_factory=dict)
    output_paths: dict[str, str] = Field(default_factory=dict)
    study: dict[str, Any] | None = None


class ProtocolReport(BaseModel):
    """A validated answer to "what is this metric, and does it differ?".

    ``run`` names the selected value of an analysis that measures several,
    such as a pair label of ``distances``; ``all_metrics`` and ``all_runs``
    list the rest, the selected one first. ``status`` is ``"complete"``, or
    ``"partial"`` when a condition could not be measured or compared;
    ``problems`` then names each condition left out of the report, or the
    comparison that failed, with its error.
    """

    model_config = ConfigDict(ser_json_inf_nan="strings")

    analysis: str
    status: str = "complete"
    problems: list[str] = Field(default_factory=list)
    trends: list[TrendReport] = Field(default_factory=list)
    protocol_version: str
    metric: str
    unit: str | None = None
    run: str | None = None
    all_metrics: list[str] = Field(default_factory=list)
    all_runs: list[str] = Field(default_factory=list)
    equilibration: str
    stride: int = 1
    frames_per_replicate: dict[str, int | list[int] | None] = Field(default_factory=dict)
    conditions: list[ConditionReport] = Field(default_factory=list)
    pairwise: list[PairwiseReport] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    provenance: ProtocolProvenance
    verdict: list[str] = Field(default_factory=list)

    def to_agent_text(self) -> str:
        """Render the report as fixed-vocabulary text, one line per item.

        Every condition, comparison, warning and verdict gets its own line. A
        comparison of labelled values, such as residues, gets instead, per
        compared condition, one line with the number of labels tested, the
        family size and the number significantly lower and higher than the
        control, and one line each listing those labels with their difference
        and adjusted p value; every per-label row stays in the JSON form of
        the report. The output has no table borders, no colour and no blank
        lines.
        """
        run = f"  run {self.run}" if self.run else ""
        counts = {item.label: item.n_replicates for item in self.conditions}
        header = (
            f"# polyzymd analyze {self.analysis}  metric {self.metric}"
            f"  unit {self.unit or 'none'}{run}  eq {self.equilibration}"
            + (f"  stride {self.stride}" if self.stride != 1 else "")
            + f"  conditions {len(counts)}"
            f"  replicates {','.join(str(n) for n in counts.values()) or 'none'}"
            f"  protocol {self.analysis}/{self.protocol_version}"
            + (f"  status {self.status}" if self.status != "complete" else "")
        )
        if any(pair.entry is not None for pair in self.pairwise):
            body = [
                *_labelled_pairwise_lines(self.pairwise),
                f"note: the {len(self.conditions)} per-label condition rows and "
                f"{len(self.pairwise)} per-label comparison rows are in the JSON report",
            ]
        else:
            body = [
                *(_condition_line(item) for item in self.conditions),
                *(_pairwise_line(item) for item in self.pairwise),
            ]
        lines = [
            header,
            *(f"problem: {text}" for text in self.problems),
            *body,
            *(_trend_line(item) for item in self.trends),
            *(f"warning: {text}" for text in self.warnings),
            *(f"verdict: {text}" for text in self.verdict),
        ]
        return "\n".join(lines) + "\n"


# Public entry points


def analyze(
    name: str,
    configs: Sequence[Path | str],
    *,
    replicates: Sequence[int] | None = None,
    equilibration: str | None = None,
    settings: dict | None = None,
    labels: Sequence[str] | None = None,
    output_dir: Path | None = None,
    recompute: bool = False,
    run: str | None = None,
    eq_check: bool = True,
    plots: bool = True,
    stride: int = 1,
    data: dict[str, Path] | None = None,
    until: str | None = None,
) -> ProtocolReport:
    """Run one analysis over one or more simulation conditions.

    The first config is the control: every comparison is control against one
    other condition. With a single config no comparison is possible and
    ``pairwise`` is empty.

    Parameters
    ----------
    name : str
        Canonical analysis name, for example ``"rmsf"``.
    configs : sequence of Path or str
        Simulation ``config.yaml`` paths, control first.
    replicates : sequence of int, optional
        Replicates for every condition. Defaults to those found on disk.
    equilibration : str, optional
        Window discarded from every replicate, for example ``"10ns"``. Applied
        uniformly; defaults to the package default.
    settings : dict, optional
        Settings of the analysis; the keys it takes and their defaults are in
        :data:`FUNCTION_ANALYSES`.
    labels : sequence of str, optional
        One label per config. Defaults to each config's directory name.
    output_dir : Path, optional
        Where ``polyzymd_results/`` and ``figures/`` are written.
    recompute : bool, optional
        Recompute replicates instead of reusing cached results.
    run : str, optional
        Value to report, for an analysis that measures several (for example
        ``contact_fraction_residues`` of ``contacts``). Defaults to the first
        one; ``all_runs`` lists the others.
    eq_check : bool, optional
        Report the pymbar detected start of the equilibrated region of each
        replicate. ``False`` skips it. It changes no value either way.
    plots : bool, optional
        Draw the figures into ``<output_dir>/figures/<name>/`` and record that
        folder in ``provenance.output_paths["figures"]``. ``False`` draws none.
    stride : int, optional
        Measure every ``stride``-th production frame of every replicate, 1 by
        default; see :meth:`~polyzymd.analyses.study.Study.from_configs`.
    data : dict of str to Path, optional
        Condition label to the directory holding its run directories on this
        machine, in place of its config's ``scratch_directory``.
    until : str, optional
        End of a common analysis window, such as ``"38ns"``; see
        :class:`~polyzymd.analyses.study.Condition`.

    Returns
    -------
    ProtocolReport
        The validated report.

    Raises
    ------
    ProtocolError
        If the name is not in :data:`FUNCTION_ANALYSES` (``catalytic_triad``
        gets the triad routine instead), a config is missing, the labels do
        not match the configs, the settings are invalid, or no replicates are
        found.
    """
    _require_known(name)
    return _analyze_function(
        name,
        configs,
        replicates=replicates,
        equilibration=equilibration,
        settings=settings,
        labels=labels,
        output_dir=output_dir,
        recompute=recompute,
        run=run,
        eq_check=eq_check,
        plots=plots,
        stride=stride,
        data=data,
        until=until,
    )


def _require_known(name: str) -> None:
    """Raise ``ProtocolError`` unless ``name`` is in :data:`FUNCTION_ANALYSES`.

    ``catalytic_triad`` gets :func:`_refuse_retired`; any other unknown name
    gets the list of analyses and the page on writing a function instead.
    """
    _refuse_retired(name)
    if name not in FUNCTION_ANALYSES:
        raise ProtocolError(
            f"No analysis named {name!r}.",
            hint=(
                f"Use one of {', '.join(FUNCTION_ANALYSES)}. For another measurement, "
                f"write a function and run it with Study.timeseries or Study.per_replicate: "
                f"{ANALYSIS_API_URL}."
            ),
        )


def _refuse_retired(name: str) -> None:
    """Raise ``ProtocolError`` for ``catalytic_triad``, which is now a routine on the study API."""
    if name == "catalytic_triad":
        raise ProtocolError(
            "catalytic_triad is no longer a polyzymd analyze analysis: the triad is now a "
            "routine on the study API, which counts each triad hydrogen bond with "
            "functions.hbond_count and combines them with Timeseries.transform.",
            hint=(
                f"Follow {TRIAD_ROUTINE_URL} (docs/source/how_to/analysis_triad_quickstart.md), "
                "or point an agent at .claude/skills/polyzymd-analyze/SKILL.md or that page. "
                "For the triad distances run polyzymd analyze distances -c <config.yaml> "
                "--set pairs=<pairs.yaml>."
            ),
        )


def _analyze_function(
    name: str,
    configs: Sequence[Path | str],
    *,
    replicates: Sequence[int] | None,
    equilibration: str | None,
    settings: dict | None,
    labels: Sequence[str] | None,
    output_dir: Path | None,
    recompute: bool,
    run: str | None,
    eq_check: bool = True,
    plots: bool = True,
    stride: int = 1,
    data: dict[str, Path] | None = None,
    until: str | None = None,
) -> ProtocolReport:
    """Measure ``name`` on every production frame and report its per-replicate mean.

    ``rg`` measures :func:`~polyzymd.analyses.functions.radius_of_gyration`
    of ``selection``. ``rmsd`` measures :func:`~polyzymd.analyses.functions.rmsd`
    of ``selection`` from the reference that ``reference_mode``,
    ``reference_frame``, ``reference_file`` and ``alignment_selection`` give
    to :func:`~polyzymd.analyses.reference.reference`. Settings left out take
    the defaults in :data:`FUNCTION_ANALYSES`. With one config the report
    summarises it; with several it compares each one with the first by
    Welch's t test. With ``plots``, ``rg`` and ``rmsd`` draw
    ``<name>_timeseries`` and ``<name>_comparison``, and ``rg`` also
    ``rg_distribution``, into ``<output_dir>/figures/<name>/``. ``rmsf`` and ``rmsd_per_residue`` go to :func:`_analyze_rmsf`, and ``distances``
    to :func:`_analyze_pairs`. Raises
    ``ProtocolError`` for a setting the analysis does not take, or a ``run``
    for ``rg`` or ``rmsd``.
    """
    from polyzymd.analyses import functions
    from polyzymd.analyses.reference import reference
    from polyzymd.analyses.timeseries import select

    unknown = set(settings or {}) - set(FUNCTION_ANALYSES[name])
    pairs = name in (
        "distances",
        "rmsf",
        "rmsd_per_residue",
        "sasa",
        "secondary_structure",
        "contacts",
        "native_contacts",
        "hydrogen_bonds",
    )
    if unknown or (run is not None and not pairs):
        raise ProtocolError(
            f"{name} takes {'' if pairs else 'no run and '}no setting other than "
            f"{', '.join(FUNCTION_ANALYSES[name])}.",
            hint=f"Run polyzymd analyze {name} -c A/config.yaml "
            + (
                "--set pairs=pairs.yaml."
                if name == "distances"
                else f"--set {next(iter(FUNCTION_ANALYSES[name]))}=..., one of the settings above."
            ),
        )
    study = _study(configs, labels, equilibration, replicates, stride, data, until)
    if name in ("rmsf", "rmsd_per_residue"):
        return _analyze_rmsf(name, study, settings, run, recompute, output_dir, plots)
    if name == "sasa":
        return _analyze_sasa(study, settings, run, recompute, output_dir, eq_check, plots)
    if name == "secondary_structure":
        return _analyze_secondary_structure(study, settings, run, recompute, output_dir, plots)
    if name == "contacts":
        return _analyze_contacts(study, settings, run, recompute, output_dir, plots)
    if name == "hydrogen_bonds":
        return _analyze_hydrogen_bonds(study, settings, run, recompute, output_dir, plots)
    if name == "native_contacts":
        return _analyze_native_contacts(
            study, settings, run, recompute, output_dir, eq_check, plots
        )
    if pairs:
        return _analyze_pairs(
            name,
            study,
            settings,
            run,
            recompute=recompute,
            output_dir=output_dir,
            eq_check=eq_check,
            plots=plots,
        )
    settings = {**FUNCTION_ANALYSES[name], **(settings or {})}
    arguments = [select(str(settings["selection"]))]
    if name == "rmsd":
        arguments.append(
            reference(
                str(settings["reference_mode"]),
                str(settings["selection"]),
                frame=settings["reference_frame"],
                file=settings["reference_file"],
                alignment=str(settings["alignment_selection"]),
            )
        )
    function = functions.rmsd if name == "rmsd" else functions.radius_of_gyration
    series = study.timeseries(
        function,
        *arguments,
        unit="A",
        name=name,
        recompute=recompute,
        output_dir=output_dir,
        bounds=(0.0, None),
    )
    values = series.reduce("mean", detect_equilibration=eq_check)
    report = values.compare() if len(study) > 1 else values.summary()
    if plots:
        folder = _figures_dir(output_dir, name)
        series.plot(folder, f"{name}_timeseries")
        values.plot(folder, f"{name}_comparison")
        if name == "rg":
            series.plot_distribution(output_dir=folder, name="rg_distribution")
        report.provenance.output_paths["figures"] = str(folder)
    return report


def _analyze_rmsf(
    name: str,
    study: Any,
    settings: dict | None,
    run: str | None,
    recompute: bool,
    output_dir: Path | None,
    plots: bool,
) -> ProtocolReport:
    """Measure the per-residue RMS deviation, RMSF and offset of every replicate in one pass.

    :func:`~polyzymd.analyses.functions.rms_decomposition` superposes
    ``alignment_selection`` on the reference of ``reference_mode``,
    ``reference_frame`` and ``reference_file``, built by
    :func:`~polyzymd.analyses.reference.reference` for both selections
    together, and gives for each residue of ``selection`` its RMS deviation
    from the reference, its RMSF about the mean position and the offset of
    the mean position from the reference, labelled by residue ID, with
    their mean squares. For ``rmsd_per_residue`` a missing ``reference_mode``
    is ``"external"`` when a ``reference_file`` is given and ``"centroid"``
    otherwise; ``rmsf`` defaults to ``"centroid"``.

    Each replicate's headline value of a quantity is the root of its mean
    square over the core residues, ``core_rmsd_per_residue``, ``core_rmsf`` and
    ``core_offset``, so ``core_rmsd_per_residue`` squared equals ``core_rmsf``
    squared plus ``core_offset`` squared, because for every atom the mean
    square deviation is the variance about the mean position plus the squared
    distance of the mean from the reference. The core is the residues of
    ``selection`` that the MDAnalysis selection ``core`` also selects, by
    default all of them, and must be the same in every replicate. Each entry
    ``region: selection`` of ``regions`` gives ``<region>_rmsd_per_residue``,
    ``<region>_rmsf`` and ``<region>_offset`` the same way. The frames are
    superposed by ``alignment_selection`` whatever the core, so fit on the
    same core atoms to measure motion within that core. ``mean_rmsf`` and
    the other plain means over every residue are kept for reference, and
    ``rmsd_per_residue``, ``rmsf`` and ``offset`` are the profiles, compared
    residue by residue. ``run`` picks the reported result and defaults to
    ``core_<name>``. The settings, the resolved reference mode and the
    residues of the core and of each region are stored in
    ``provenance.settings``. With ``plots``, ``<part>_profile`` draws each
    profile with ``highlight_residues`` marked, ``rms_decomposition`` the
    three profiles of each condition together, ``rmsf_comparison`` the three
    core values, and with several conditions ``<part>_difference`` each
    condition's per-residue difference from the control with its interval
    and significant residues, into ``<output_dir>/figures/<name>/``.
    """
    import numpy as np

    from polyzymd.analyses import functions
    from polyzymd.analyses.figures import plot_decomposition, plot_differences, plot_values
    from polyzymd.analyses.reference import reference
    from polyzymd.analyses.timeseries import select

    settings = {**FUNCTION_ANALYSES[name], **(settings or {})}
    atoms, fit = str(settings["selection"]), str(settings["alignment_selection"])
    mode = settings["reference_mode"] or ("external" if settings["reference_file"] else "centroid")
    regions = settings["regions"] or {}
    if not isinstance(regions, dict) or {"core", "mean"} & set(regions):
        raise ProtocolError(
            f"{name}: regions must map names other than core and mean to selections, "
            f"got {regions!r}.",
            hint="Pass --set regions='{lid: resid 70-90}'.",
        )
    sets = {"core": settings["core"] or "all", **{str(k): str(v) for k, v in regions.items()}}
    runs = [f"{kind}_{part}" for kind in ("core", *regions, "mean") for part in functions.RMS_PARTS]
    runs += list(functions.RMS_PARTS)
    run = run or f"core_{name}"
    if run not in runs:
        raise ProtocolError(
            f"{name}: no result named {run!r}.", hint=f"Use --run with one of {runs}."
        )
    residues = {key: _residue_ids(study, f"({atoms}) and ({value})") for key, value in sets.items()}
    rows = study.per_replicate(
        functions.rms_decomposition,
        select(atoms),
        select(fit),
        reference(
            str(mode),
            f"({atoms}) or ({fit})",
            frame=settings["reference_frame"],
            file=settings["reference_file"],
            alignment=fit,
        ),
        unit="A",
        labels=lambda u: u.select_atoms(atoms).residues.resids,
        name="rms_decomposition",
        recompute=recompute,
        output_dir=output_dir,
        bounds=(0.0, None),
        parts=functions.RMS_PARTS + functions.MS_PARTS,
    )
    profiles = {part: rows[part] for part in functions.RMS_PARTS}
    results = {}
    for key, labels in residues.items():
        for part, square in zip(functions.RMS_PARTS, functions.MS_PARTS, strict=True):
            results[f"{key}_{part}"] = rows[square].over_labels(
                lambda values: float(np.sqrt(np.mean(values))), f"{key}_{part}", labels
            )
            results[f"{key}_{part}"].unit = "A"
            results[f"{key}_{part}"].bounds = (0.0, None)
    results.update(
        {f"mean_{part}": values.over_labels("mean") for part, values in profiles.items()}
    )
    results.update(profiles)
    values = results[run]
    report = values.compare() if len(study) > 1 else values.summary()
    report.provenance.settings = {
        **settings,
        "reference_mode": mode,
        "residues": {key: [int(r) for r in value] for key, value in residues.items()},
    }
    if plots:
        folder = _figures_dir(output_dir, name)
        highlight = settings["highlight_residues"] or []
        names = {"rmsd_per_residue": "RMS deviation", "rmsf": "RMSF", "offset": "offset"}
        for part, profile in profiles.items():
            title = f"Per-residue {names[part]}"
            profile.plot(folder, f"{part}_profile", title, None, highlight, "Residue")
            if len(study) > 1:
                compared = report if run == part else profile.compare()
                figure = f"{part}_difference"
                plot_differences(profile, compared, folder, figure, None, None, "Residue")
        shown = {names[part]: profile for part, profile in profiles.items()}
        plot_decomposition(shown, folder, "rms_decomposition", None, None, "Residue")
        core = [results[f"core_{part}"] for part in functions.RMS_PARTS]
        labels = [names[part] for part in functions.RMS_PARTS]
        plot_values(core, labels, folder, "rmsf_comparison", "Root mean square over the core")
        report.provenance.output_paths["figures"] = str(folder)
    return report.model_copy(update={"analysis": name, "run": run, "all_runs": runs})


def _empty_selections(study: Any, selections: dict[str, str]) -> dict[tuple[str, int], list[str]]:
    """Return, for each replicate where a named selection matches no atoms, those selections."""
    empty: dict[tuple[str, int], list[str]] = {}
    for condition in study:
        for replicate in condition.replicates:
            universe = replicate.universe()
            missing = [
                f"{name} {selection!r}"
                for name, selection in selections.items()
                if len(universe.select_atoms(selection)) == 0
            ]
            if missing:
                empty[(condition.label, replicate.index)] = missing
    return empty


def _first_universe(study: Any, empty: dict[tuple[str, int], list[str]], analysis: str) -> Any:
    """Return the universe of the first replicate whose selections all match atoms."""
    for condition in study:
        for replicate in condition.replicates:
            if (condition.label, replicate.index) not in empty:
                return replicate.universe()
    missing = sorted({name for names in empty.values() for name in names})
    raise ProtocolError(
        f"{analysis}: the selections {', '.join(missing)} match no atoms in any replicate.",
        hint="Choose selections that pick atoms, such as 'chainid A' for the protein and "
        "'chainid C' for the polymer.",
    )


def _report_skipping(
    values: Any, study: Any, empty: dict[tuple[str, int], list[str]], analysis: str
) -> ProtocolReport:
    """Summarise or compare ``values`` without the replicates where a selection matched no atoms.

    Those replicates are left out of every statistic, and a condition left
    without replicates is left out of the report, each with a warning. When
    the control is left out, the other conditions are summarised and not
    compared.
    """
    if empty:
        values.rows = {
            label: [row for row in rows if (label, row[0]) not in empty]
            for label, rows in values.rows.items()
        }
    labels = [condition.label for condition in study]
    kept = [label for label in labels if values.rows.get(label)]
    control = labels[0]
    if len(kept) > 1 and control in kept:
        report = values.compare(control=control, conditions=kept)
    else:
        report = values.summary(conditions=kept)
    by_condition: dict[str, list[str]] = {}
    for (label, index), missing in sorted(empty.items()):
        by_condition.setdefault(label, []).append(f"{index} ({'; '.join(missing)})")
    for label, entries in by_condition.items():
        left = "is left out" if label not in kept else "are left out"
        what = "the condition" if label not in kept else "those replicates"
        report.warnings.append(
            f"{analysis}: in condition {label}, replicate {', '.join(entries)} matched no "
            f"atoms, so {what} {left} of the statistics."
        )
    if len(labels) > 1 and control not in kept:
        report.warnings.append(
            f"{analysis}: the control {control} has no replicate where every selection matches "
            "atoms, so the other conditions are summarised and not compared. Give a condition "
            "with those atoms first to compare against it."
        )
    return report


def _hbond_summaries(settings: dict) -> dict[str, tuple[str, str | None]]:
    """Return each summary's name with the selections of its groups, the second ``None`` for within."""
    groups = settings["groups"] or {}
    summaries = settings["summaries"] or {}
    if isinstance(summaries, list):
        summaries = {str(item.get("name")): item for item in summaries if isinstance(item, dict)}
    if (
        not isinstance(groups, dict)
        or not groups
        or not isinstance(summaries, dict)
        or not summaries
    ):
        raise ProtocolError(
            "hydrogen_bonds: groups must map names to selections and summaries must map names "
            "to {between: [group, group]} or {within: group}.",
            hint="Pass --set groups='{protein: chainid A, polymer: chainid C}' --set "
            "summaries='{protein_polymer: {between: [protein, polymer]}}'.",
        )
    resolved: dict[str, tuple[str, str | None]] = {}
    for name, spec in summaries.items():
        spec = spec if isinstance(spec, dict) else {}
        between, within = spec.get("between"), spec.get("within")
        names = list(between) if between is not None else [within]
        if (between is None) == (within is None) or (between is not None and len(names) != 2):
            raise ProtocolError(
                f"hydrogen_bonds: summary {name!r} needs exactly one of between: [group, group] "
                "or within: group.",
                hint="Write summaries='{protein_polymer: {between: [protein, polymer]}}'.",
            )
        unknown = [group for group in names if group not in groups]
        if unknown:
            raise ProtocolError(
                f"hydrogen_bonds: summary {name!r} names groups {unknown} that groups does not "
                f"define; it defines {sorted(groups)}.",
                hint="Add the group to --set groups='{name: selection}'.",
            )
        first = str(groups[names[0]])
        resolved[str(name)] = (first, None if within is not None else str(groups[names[1]]))
    return resolved


def _hbond_residue_labels(universe: Any, selection: str) -> list:
    """Return labels for the residues of ``selection``: residue IDs, or ``chain:resid`` if IDs repeat."""
    residues = universe.select_atoms(selection).residues
    ids = [int(r.resid) for r in residues]
    if len(set(ids)) == len(ids):
        return ids
    labels = [f"{r.atoms[0].chainID}:{int(r.resid)}" for r in residues]
    if len(set(labels)) != len(labels):
        raise ProtocolError(
            f"hydrogen_bonds: the residues of {selection!r} repeat residue IDs even within a "
            "chain, so they cannot be told apart for the per-residue result.",
            hint="Make the summary's first group a selection of distinct residues, such as "
            "'chainid A'.",
        )
    return labels


def _analyze_hydrogen_bonds(
    study: Any,
    settings: dict | None,
    run: str | None,
    recompute: bool,
    output_dir: Path | None,
    plots: bool,
) -> ProtocolReport:
    """Count hydrogen bonds between or within named groups and report one result.

    ``groups`` maps names to MDAnalysis selections, and each entry of
    ``summaries`` is ``{between: [a, b]}``, hydrogen bonds with one partner in
    each group, or ``{within: a}``. :func:`~polyzymd.analyses.functions.hydrogen_bonds`
    runs MDAnalysis ``HydrogenBondAnalysis`` once per replicate for the chosen
    summary, with ``d_a_cutoff`` Å and ``d_h_a_angle_cutoff`` degrees, and
    donors, hydrogens and acceptors from
    :func:`~polyzymd.analyses.functions.hbond_atoms` unless the selections
    ``donors``, ``hydrogens`` or ``acceptors`` are given. Each summary ``s``
    gives ``s_mean_hbonds`` (the default for the first summary),
    ``s_mean_residue_pairs`` and ``s_any_fraction``. The atoms counted as
    hydrogens and acceptors are recorded under ``provenance.settings``, as
    counts per residue name and atom name. With ``plots``,
    ``hbonds_<run>_comparison`` goes to ``<output_dir>/figures/hydrogen_bonds/``.
    """
    from collections import Counter

    from polyzymd.analyses import functions
    from polyzymd.analyses.figures import plot_differences
    from polyzymd.analyses.timeseries import select

    settings = {**FUNCTION_ANALYSES["hydrogen_bonds"], **(settings or {})}
    summaries = _hbond_summaries(settings)
    parts = list(functions.HBOND_PARTS)
    life = {"mean_lifetime": 0, "lifetime_events": 1, "censored_fraction": 2}
    kinds = [*parts, *life, "residues", "pairs"]
    runs = [f"{name}_{kind}" for name in summaries for kind in kinds]
    run = run or runs[0]
    if run not in runs:
        raise ProtocolError(
            f"hydrogen_bonds: no result named {run!r}.", hint=f"Use --run with one of {runs}."
        )
    summary = next(
        name for name in summaries if run.startswith(f"{name}_") and run[len(name) + 1 :] in kinds
    )
    part = run[len(summary) + 1 :]
    if settings["lifetime_key"] not in ("residue", "atom"):
        raise ProtocolError(
            f"hydrogen_bonds: lifetime_key must be 'residue' or 'atom', got "
            f"{settings['lifetime_key']!r}.",
            hint="Pass --set lifetime_key=residue for residue pairs, or atom for atom pairs.",
        )
    if float(settings["tolerance_ps"]) < 0:
        raise ProtocolError(
            f"hydrogen_bonds: tolerance_ps must be at least 0, got {settings['tolerance_ps']}.",
            hint="Pass --set tolerance_ps=0 for bonds that end at the first absent frame.",
        )
    first, second = summaries[summary]
    both = first if second is None else f"({first}) or ({second})"
    group_selections = {
        "first group": first,
        **({} if second is None else {"second group": second}),
    }
    skipped = _empty_selections(study, group_selections)
    universe = _first_universe(study, skipped, "hydrogen_bonds")
    explicit = {
        key: select(f"({both}) and ({settings[key]})", allow_empty=True)
        for key in ("donors", "hydrogens", "acceptors")
        if settings[key]
    }
    atoms = universe.select_atoms(both)
    hydrogens, acceptors = (
        (None, None)
        if {"hydrogens", "acceptors"} <= set(explicit)
        else functions.hbond_atoms(atoms)
    )
    if "hydrogens" in explicit:
        hydrogens = universe.select_atoms(f"({both}) and ({settings['hydrogens']})")
    if "acceptors" in explicit:
        acceptors = universe.select_atoms(f"({both}) and ({settings['acceptors']})")

    def counts(group: Any) -> dict[str, int]:
        return dict(sorted(Counter(f"{a.resname} {a.name}" for a in group).items()))

    donor_atoms = (
        universe.select_atoms(f"({both}) and ({settings['donors']})")
        if "donors" in explicit
        else sum((h.bonded_atoms[:1] for h in hydrogens), universe.atoms[[]])
    )
    # Only non-default options, so a plain Python call reuses the stored records.
    options: dict[str, Any] = dict(explicit)
    if float(settings["d_a_cutoff"]) != functions.HBOND_DISTANCE:
        options["d_a_cutoff"] = float(settings["d_a_cutoff"])
    if float(settings["d_h_a_angle_cutoff"]) != functions.HBOND_ANGLE:
        options["d_h_a_angle_cutoff"] = float(settings["d_h_a_angle_cutoff"])
    arguments = [select(first, allow_empty=True)] + (
        [] if second is None else [select(second, allow_empty=True)]
    )
    if part in parts:
        rows = study.per_replicate(
            functions.hydrogen_bonds,
            *arguments,
            unit=None,
            name=f"hydrogen_bonds_{summary}",
            recompute=recompute,
            output_dir=output_dir,
            parts=parts,
            **options,
        )
        values = rows[part]
        values.bounds = (0.0, 1.0) if part == "any_fraction" else (0.0, None)
    elif part in life:
        lifetime_options = dict(options)
        if settings["lifetime_key"] != "residue":
            lifetime_options["key"] = settings["lifetime_key"]
        if float(settings["tolerance_ps"]) > 0:
            lifetime_options["tolerance_ps"] = float(settings["tolerance_ps"])
        rows = study.per_replicate(
            functions.hbond_lifetimes,
            *arguments,
            unit=None,
            name=f"hbond_lifetimes_{summary}",
            recompute=recompute,
            output_dir=output_dir,
            parts=list(functions.LIFETIME_PARTS),
            **lifetime_options,
        )
        values = rows[functions.LIFETIME_PARTS[life[part]]]
        values.unit, values.bounds = {
            "mean_lifetime": ("ns", (0.0, None)),
            "lifetime_events": (None, (0.0, None)),
            "censored_fraction": (None, (0.0, 1.0)),
        }[part]
    elif part == "pairs":
        values = study.per_replicate(
            functions.residue_pair_hbond_occupancy,
            *arguments,
            unit=None,
            labels="returned",
            missing=0.0,
            name=f"residue_pair_hbond_occupancy_{summary}",
            recompute=recompute,
            output_dir=output_dir,
            bounds=(0.0, 1.0),
            **options,
        )
    else:
        values = study.per_replicate(
            functions.residue_hbond_occupancy,
            *arguments,
            unit=None,
            labels=lambda u: _hbond_residue_labels(u, first),
            name=f"residue_hbond_occupancy_{summary}",
            recompute=recompute,
            output_dir=output_dir,
            bounds=(0.0, 1.0),
            **options,
        )
    values.metric = run
    report = _report_skipping(values, study, skipped, "hydrogen_bonds")
    if part in life:
        empty = [
            f"{label} replicate {row[0]}"
            for label, table in values.rows.items()
            for row in table
            if row[1] != row[1]
        ]
        if empty:
            report.warnings.append(
                f"hydrogen_bonds: {', '.join(empty)} have no hydrogen bond in summary "
                f"{summary!r}, so {run} is undefined (nan) there."
            )
    report.provenance.settings = {
        **settings,
        "summary": {"name": summary, "groups": [first] if second is None else [first, second]},
        "hbond_atoms": {
            "donors": counts(donor_atoms),
            "hydrogens": len(hydrogens),
            "acceptors": counts(acceptors),
        },
    }
    if plots:
        folder = _figures_dir(output_dir, "hydrogen_bonds")
        if part in ("residues", "pairs"):
            values.plot(
                folder,
                f"hbonds_{run}_profile",
                f"H-bond occupancy, {summary}",
                None,
                [],
                "Residue" if part == "residues" else "Residue pair",
            )
            if report.pairwise:
                plot_differences(
                    values,
                    report,
                    folder,
                    f"hbonds_{run}_difference",
                    None,
                    None,
                    "Residue" if part == "residues" else "Residue pair",
                )
        else:
            values.plot(folder, f"hbonds_{run}_comparison", title=run.replace("_", " "))
        report.provenance.output_paths["figures"] = str(folder)
    return report.model_copy(update={"analysis": "hydrogen_bonds", "run": run, "all_runs": runs})


def _analyze_native_contacts(
    study: Any,
    settings: dict | None,
    run: str | None,
    recompute: bool,
    output_dir: Path | None,
    eq_check: bool,
    plots: bool,
) -> ProtocolReport:
    """Measure the fraction of native contacts Q on every production frame and report its mean.

    :func:`~polyzymd.analyses.functions.native_contacts` measures ``selection``
    against the reference of ``reference_mode``, ``reference_frame`` and
    ``reference_file``, built by :func:`~polyzymd.analyses.reference.reference`.
    A missing ``reference_mode`` is ``"external"`` when a ``reference_file``
    is given, and otherwise ``"frame"``, production frame ``reference_frame``.
    ``radius``, ``min_separation``, ``beta`` and ``lambda_constant`` define
    the native pairs and the switching function; the defaults, on heavy
    atoms, are the definition of Best, Hummer and Eaton (2013). ``q``, the
    default result, counts every native pair; ``<region>_q`` for each entry
    ``region: selection`` of ``regions`` counts the pairs with at least one
    atom in the region. Only the chosen result is measured. Each replicate's
    value is its mean Q over production frames. With ``plots``,
    ``native_contacts_timeseries_<run>`` and ``native_contacts_comparison_<run>``
    go to ``<output_dir>/figures/native_contacts/``.
    """
    from polyzymd.analyses import functions
    from polyzymd.analyses.reference import reference
    from polyzymd.analyses.timeseries import select

    settings = {**FUNCTION_ANALYSES["native_contacts"], **(settings or {})}
    selection = str(settings["selection"])
    mode = settings["reference_mode"] or ("external" if settings["reference_file"] else "frame")
    regions = settings["regions"] or {}
    if not isinstance(regions, dict) or "q" in regions:
        raise ProtocolError(
            f"native_contacts: regions must map names other than 'q' to selections, got {regions!r}.",
            hint="Pass --set regions='{active_site: resid 70-90}'.",
        )
    runs = ["q", *(f"{name}_q" for name in regions)]
    run = run or "q"
    if run not in runs:
        raise ProtocolError(
            f"native_contacts: no result named {run!r}.", hint=f"Use --run with one of {runs}."
        )
    arguments = [
        select(selection),
        reference(
            mode,
            selection,
            frame=settings["reference_frame"],
            file=settings["reference_file"],
            alignment=selection,
        ),
    ]
    if run != "q":
        arguments.append(select(f"({selection}) and ({regions[run[: -len('_q')]]})"))
    # Only non-default options, so a plain Python call reuses the stored records.
    defaults = {
        "radius": functions.NATIVE_CONTACT_RADIUS,
        "min_separation": functions.NATIVE_CONTACT_SEPARATION,
        "beta": 5.0,
        "lambda_constant": 1.8,
    }
    options: dict[str, Any] = {
        key: settings[key] for key, value in defaults.items() if settings[key] != value
    }
    if not settings["use_pbc"]:
        options["pbc"] = False
    series = study.timeseries(
        functions.native_contacts,
        *arguments,
        unit=None,
        name=f"native_contacts_{run}",
        recompute=recompute,
        output_dir=output_dir,
        bounds=(0.0, 1.0),
        **options,
    )
    values = series.reduce("mean", detect_equilibration=eq_check)
    values.metric = f"mean_{run}"
    report = values.compare() if len(study) > 1 else values.summary()
    report.provenance.settings = {**settings, "reference_mode": mode}
    if plots:
        folder = _figures_dir(output_dir, "native_contacts")
        series.plot(folder, f"native_contacts_timeseries_{run}")
        values.plot(folder, f"native_contacts_comparison_{run}", title=f"Native contacts, {run}")
        report.provenance.output_paths["figures"] = str(folder)
    return report.model_copy(update={"analysis": "native_contacts", "run": run, "all_runs": runs})


def _analyze_sasa(
    study: Any,
    settings: dict | None,
    run: str | None,
    recompute: bool,
    output_dir: Path | None,
    eq_check: bool,
    plots: bool,
) -> ProtocolReport:
    """Measure the SASA of ``target`` in one context and report its total or its residues.

    ``contexts`` maps a name to the MDAnalysis selection of the atoms present
    in the calculation, for example ``{isolated: protein, with_polymer:
    protein or resname SBM EGM}``; each must contain every ``target`` atom.
    With no contexts, the target is measured alone under the name
    ``isolated``. Each context gives two results: ``<name>``, the mean over
    production frames of the target's total SASA from
    :func:`~polyzymd.analyses.functions.sasa`, and ``<name>_residues``, each
    target residue's mean SASA from
    :func:`~polyzymd.analyses.functions.residue_sasa`, compared residue by
    residue. Only the result that ``run`` picks is measured, by default the
    first context's total, because every context is a separate Shrake-Rupley
    pass over every frame. With ``plots``, a total draws
    ``sasa_timeseries_<name>``, ``sasa_comparison_<name>`` and
    ``sasa_distribution_<name>``, and a residue result ``sasa_profile_<name>``
    and, with several conditions, ``sasa_difference_<name>``, into
    ``<output_dir>/figures/sasa/``.
    """
    from polyzymd.analyses import functions
    from polyzymd.analyses.figures import plot_differences
    from polyzymd.analyses.timeseries import select

    settings = {**FUNCTION_ANALYSES["sasa"], **(settings or {})}
    target = str(settings["target"])
    contexts = settings["contexts"] or {"isolated": target}
    if not isinstance(contexts, dict) or not all(
        isinstance(key, str) and isinstance(value, str) and not key.endswith("_residues")
        for key, value in contexts.items()
    ):
        raise ProtocolError(
            f"sasa: contexts must map names, none ending in _residues, to selections, "
            f"got {contexts!r}.",
            hint="Pass --set contexts='{isolated: protein, with_polymer: protein or resname SBM EGM}'.",
        )
    runs = [key for name in contexts for key in (name, f"{name}_residues")]
    run = run or runs[0]
    if run not in runs:
        raise ProtocolError(
            f"sasa: no result named {run!r}.", hint=f"Use --run with one of {runs}."
        )
    residues = run.endswith("_residues") and run[: -len("_residues")] in contexts
    context = contexts[run[: -len("_residues")] if residues else run]
    # Pass only non-default options, so the stored record equals that of a
    # study.timeseries(functions.sasa, ...) call left at the defaults.
    defaults = {
        "probe_radius_nm": functions.SASA_PROBE_RADIUS_NM,
        "n_sphere_points": functions.SASA_SPHERE_POINTS,
    }
    options = {
        key: kind(settings[key])
        for key, kind in (("probe_radius_nm", float), ("n_sphere_points", int))
        if kind(settings[key]) != defaults[key]
    }
    folder = _figures_dir(output_dir, "sasa") if plots else None
    if residues:
        values = study.per_replicate(
            functions.residue_sasa,
            select(target),
            select(context),
            unit="A^2",
            labels=lambda u: u.select_atoms(target).residues.resids,
            name=f"sasa_{run}",
            recompute=recompute,
            output_dir=output_dir,
            bounds=(0.0, None),
            **options,
        )
        report = values.compare() if len(study) > 1 else values.summary()
        if plots:
            name = run[: -len("_residues")]
            values.plot(
                folder, f"sasa_profile_{name}", f"Per-residue SASA, {name}", None, [], "Residue"
            )
            if len(study) > 1:
                plot_differences(
                    values, report, folder, f"sasa_difference_{name}", None, None, "Residue"
                )
    else:
        series = study.timeseries(
            functions.sasa,
            select(target),
            select(context),
            unit="A^2",
            name=f"sasa_{run}",
            recompute=recompute,
            output_dir=output_dir,
            bounds=(0.0, None),
            **options,
        )
        values = series.reduce("mean", detect_equilibration=eq_check)
        values.metric = "mean_sasa"
        report = values.compare() if len(study) > 1 else values.summary()
        if plots:
            series.plot(folder, f"sasa_timeseries_{run}")
            values.plot(folder, f"sasa_comparison_{run}", title=f"SASA, {run}")
            series.plot_distribution(output_dir=folder, name=f"sasa_distribution_{run}")
    report.provenance.settings = {**settings, "contexts": dict(contexts)}
    if plots:
        report.provenance.output_paths["figures"] = str(folder)
    return report.model_copy(update={"analysis": "sasa", "run": run, "all_runs": runs})


def _analyze_secondary_structure(
    study: Any,
    settings: dict | None,
    run: str | None,
    recompute: bool,
    output_dir: Path | None,
    plots: bool,
) -> ProtocolReport:
    """Assign DSSP classes to every residue of ``selection`` and report one class.

    ``scheme`` picks MDTraj's DSSP: ``simplified`` (the default) gives helix,
    strand, coil and unassigned of
    :data:`~polyzymd.analyses.functions.DSSP_SIMPLIFIED`, and ``full`` the
    eight classes and unassigned of
    :data:`~polyzymd.analyses.functions.DSSP_CLASSES`.
    :func:`~polyzymd.analyses.functions.dssp_occupancy` gives, in one pass per
    replicate, each residue's fraction of production frames in each class.
    ``<name>`` is a class's mean over the residues, the fraction of
    residue-frames in it, and ``<name>_residues`` its per-residue profile,
    compared residue by residue. ``run`` defaults to the first class, helix or
    alpha_helix. A warning names the replicates with unassigned residues,
    which MDTraj gives ``"NA"`` when it cannot assign them. With ``plots``,
    ``ss_content_bars`` groups every class's fraction except unassigned, a
    total draws ``ss_<name>_comparison``, and a residue result
    ``ss_<name>_profile``, ``ss_classes_<name>`` (every class of each residue
    per condition) and, with several conditions, ``ss_<name>_difference``,
    into ``<output_dir>/figures/secondary_structure/``.
    """
    from polyzymd.analyses import functions
    from polyzymd.analyses.figures import plot_decomposition, plot_differences, plot_values
    from polyzymd.analyses.timeseries import select

    settings = {**FUNCTION_ANALYSES["secondary_structure"], **(settings or {})}
    atoms, scheme = str(settings["selection"]), settings["scheme"]
    if scheme not in ("simplified", "full"):
        raise ProtocolError(
            f"secondary_structure: scheme must be simplified or full, got {scheme!r}.",
            hint="Pass --set scheme=full for the eight DSSP classes.",
        )
    classes = list(functions.DSSP_SIMPLIFIED if scheme == "simplified" else functions.DSSP_CLASSES)
    runs = [key for name in classes for key in (name, f"{name}_residues")]
    run = run or classes[0]
    if run not in runs:
        raise ProtocolError(
            f"secondary_structure: no result named {run!r} in the {scheme} scheme.",
            hint=f"Use --run with one of {runs}, or --set scheme="
            + ("full" if scheme == "simplified" else "simplified")
            + " for the other classes.",
        )
    rows = study.per_replicate(
        functions.dssp_occupancy,
        select(atoms),
        unit=None,
        labels=lambda u: u.select_atoms(atoms).residues.resids,
        name=f"dssp_occupancy_{scheme}",
        recompute=recompute,
        output_dir=output_dir,
        bounds=(0.0, 1.0),
        parts=classes,
        simplified=scheme == "simplified",
    )
    totals = {name: rows[name].over_labels("mean", f"{name}_fraction") for name in classes}
    residues = run.endswith("_residues")
    values = rows[run[: -len("_residues")]] if residues else totals[run]
    report = values.compare() if len(study) > 1 else values.summary()
    unassigned = [
        f"{label} replicate {row[0]}"
        for label, table in totals["unassigned"].rows.items()
        for row in table
        if row[1] > 0
    ]
    if unassigned:
        report.warnings.append(
            "MDTraj could not assign a DSSP class to some residues (code NA) in "
            + ", ".join(unassigned)
            + "; they count in unassigned. Check for missing backbone atoms or "
            "non-standard residue names."
        )
    report.provenance.settings = dict(settings)
    if plots:
        folder = _figures_dir(output_dir, "secondary_structure")
        shown = [name for name in classes if name != "unassigned"]
        plot_values(
            [totals[name] for name in shown],
            shown,
            folder,
            "ss_content_bars",
            "Secondary structure",
        )
        if residues:
            name = run[: -len("_residues")]
            values.plot(folder, f"ss_{name}_profile", f"Per-residue {name}", None, [], "Residue")
            profiles = {part: rows[part] for part in shown}
            plot_decomposition(profiles, folder, f"ss_classes_{name}", None, None, "Residue")
            if len(study) > 1:
                plot_differences(
                    values, report, folder, f"ss_{name}_difference", None, None, "Residue"
                )
        else:
            values.plot(folder, f"ss_{run}_comparison", title=f"{run} fraction")
        report.provenance.output_paths["figures"] = str(folder)
    return report.model_copy(
        update={"analysis": "secondary_structure", "run": run, "all_runs": runs}
    )


#: Settings of polyzymd analyze contacts that only one method reads.
CONTACT_METHOD_SETTINGS = {
    "distance": ("cutoff", "heavy_atoms"),
    "occlusion": (
        "exposed_threshold",
        "buried_threshold",
        "max_asa",
        "probe_radius_nm",
        "n_sphere_points",
    ),
}


def _analyze_contacts(
    study: Any,
    settings: dict | None,
    run: str | None,
    recompute: bool,
    output_dir: Path | None,
    plots: bool,
) -> ProtocolReport:
    """Measure each protein residue's contact with the polymer and report one result.

    ``method`` picks what contact means on one frame:

    - ``occlusion`` (default):
      :func:`~polyzymd.analyses.functions.residue_occlusion` computes each
      residue's SASA with the protein alone and with the polymer. A residue is
      in contact when it is exposed without the polymer (relative SASA, over
      its maximum ASA from Tien et al. 2013, column ``max_asa``, at least
      ``exposed_threshold``) and buried by it (relative SASA below
      ``buried_threshold``, and lower than without the polymer). Residues without a maximum ASA,
      such as terminal caps, are not measured, and a warning names them.
    - ``distance``: :func:`~polyzymd.analyses.functions.residue_contacts`
      counts a contact when any polymer atom is within ``cutoff`` Å of the
      residue, comparing heavy atoms only when ``heavy_atoms`` is true.

    ``polymer_selection`` is narrowed to the residue names in
    ``polymer_types`` when given, and ``use_pbc`` uses the frame's box: the
    minimum image for ``distance``, and for ``occlusion`` each polymer
    molecule moved whole to its image nearest the protein. Results:

    - ``coverage``, the default: the fraction of residues in contact on at
      least one frame;
    - ``mean_contact_fraction``: the mean over residues of the fraction of
      frames in contact;
    - ``contact_fraction_residues``: each residue's contact fraction, compared
      residue by residue;
    - ``<type>_contact_fraction`` and ``<type>_contact_fraction_residues`` for
      each polymer residue name, such as one monomer type; for ``occlusion``,
      a contact with only that type's atoms present;
    - ``<class>_contact_fraction`` for each amino-acid class of
      :class:`~polyzymd.analyses.shared.groupings.base.ProteinAAClassification`
      present, and ``<region>_contact_fraction`` for each entry
      ``region: selection`` of ``regions``: the mean contact fraction of those
      residues;
    - for ``occlusion`` only, ``occluded_area``: the SASA in Å² the polymer
      removes from the measured residues per frame, ``occlusion_fraction``:
      that area over their SASA with the protein alone, and
      ``occluded_area_residues``: each residue's mean occluded area.

    With ``plots``, ``contacts_class_bars`` groups the classes, a one-value
    result draws ``contacts_<run>_comparison``, and a residue result
    ``contacts_<name>_profile`` and, with several conditions,
    ``contacts_<name>_difference``, into ``<output_dir>/figures/contacts/``.
    """
    import numpy as np

    from polyzymd.analyses import functions
    from polyzymd.analyses.figures import plot_differences, plot_values
    from polyzymd.analyses.shared.aa_classification import get_max_asa
    from polyzymd.analyses.shared.groupings.base import ProteinAAClassification
    from polyzymd.analyses.timeseries import ReplicateValues, select

    given = dict(settings or {})
    settings = {**FUNCTION_ANALYSES["contacts"], **given}
    method = settings["method"]
    if method not in CONTACT_METHOD_SETTINGS:
        raise ProtocolError(
            f"contacts: method must be 'occlusion' or 'distance', got {method!r}.",
            hint="Pass --set method=occlusion or --set method=distance.",
        )
    other = next(name for name in CONTACT_METHOD_SETTINGS if name != method)
    misplaced = sorted(set(given) & set(CONTACT_METHOD_SETTINGS[other]))
    if misplaced:
        raise ProtocolError(
            f"contacts: {', '.join(misplaced)} only apply to method={other}, and method is "
            f"{method}.",
            hint=f"Drop {', '.join(misplaced)}, or pass --set method={other}.",
        )
    if method == "occlusion" and settings["max_asa"] not in ("theoretical", "empirical"):
        raise ProtocolError(
            f"contacts: max_asa must be 'theoretical' or 'empirical', got {settings['max_asa']!r}.",
            hint="Pass --set max_asa=theoretical, the values Tien et al. 2013 recommend.",
        )
    tolerance = float(settings["tolerance_ps"])
    if tolerance < 0:
        raise ProtocolError(
            f"contacts: tolerance_ps must be at least 0, got {tolerance}.",
            hint="Pass --set tolerance_ps=0 for events that end at the first absent frame.",
        )
    protein = str(settings["protein_selection"])
    polymer = str(settings["polymer_selection"])
    types_filter = settings["polymer_types"]
    if types_filter:
        names = [types_filter] if isinstance(types_filter, str) else list(types_filter)
        polymer = f"({polymer}) and (resname {' '.join(str(name) for name in names)})"
    if method == "distance" and settings["heavy_atoms"]:
        protein = f"({protein}) and not element H"
        polymer = f"({polymer}) and not element H"
    regions = settings["regions"] or {}
    skipped = _empty_selections(study, {"protein_selection": protein, "polymer_selection": polymer})
    first = _first_universe(study, skipped, "contacts")
    protein_atoms = first.select_atoms(protein)
    polymer_atoms = first.select_atoms(polymer)
    types = sorted({str(name) for name in polymer_atoms.resnames})

    def measured(residues: Any) -> list:
        if method == "distance":
            return list(residues)
        return [r for r in residues if get_max_asa(str(r.resname), settings["max_asa"]) is not None]

    kept_residues = set(measured(protein_atoms.residues))
    unmeasured = [f"{r.resname}{r.resid}" for r in protein_atoms.residues if r not in kept_residues]
    grouping = ProteinAAClassification()
    by_class: dict[str, list[int]] = {}
    for residue in measured(protein_atoms.residues):
        by_class.setdefault(grouping.classify(str(residue.resname)), []).append(int(residue.resid))
    classes = [name for name in grouping.available_groups if name in by_class]
    reserved = {"coverage", "mean", "contact", "classes", "occluded", "occlusion", *types, *classes}
    if not isinstance(regions, dict) or reserved & set(regions):
        raise ProtocolError(
            f"contacts: regions must map names other than {sorted(reserved)} to selections, "
            f"got {regions!r}.",
            hint="Pass --set regions='{lid: resid 70-90}'.",
        )
    kept = {int(r.resid) for r in kept_residues}
    region_ids = {
        name: [i for i in _residue_ids(study, f"({protein}) and ({selection})") if i in kept]
        for name, selection in regions.items()
    }
    empty = [name for name, ids in region_ids.items() if not ids]
    if empty:
        raise ProtocolError(
            f"contacts: regions {empty} have no measured residue.",
            hint="Choose regions with standard amino acids; residues without a maximum ASA "
            "are not measured by method=occlusion.",
        )
    type_parts = [f"{name}_contact_fraction" for name in types]
    # Only non-default options, so a plain Python call reuses the stored records.
    options: dict[str, Any] = {"types": types}
    if not settings["use_pbc"]:
        options["pbc"] = False
    if method == "distance":
        function, name, parts = functions.residue_contacts, "residue_contacts", ["contact_fraction"]
        if float(settings["cutoff"]) != functions.CONTACT_CUTOFF:
            options["cutoff"] = float(settings["cutoff"])
    else:
        function, name, parts = (
            functions.residue_occlusion,
            "residue_occlusion",
            list(functions.OCCLUSION_PARTS),
        )
        defaults = {
            "exposed_threshold": functions.EXPOSED_THRESHOLD,
            "buried_threshold": functions.BURIED_THRESHOLD,
            "max_asa": "theoretical",
            "probe_radius_nm": functions.SASA_PROBE_RADIUS_NM,
            "n_sphere_points": functions.SASA_SPHERE_POINTS,
        }
        options |= {key: settings[key] for key, value in defaults.items() if settings[key] != value}
    occlusion = method == "occlusion"
    fraction_runs = [
        "coverage",
        "mean_contact_fraction",
        *(f"{name}_contact_fraction" for name in [*types, *classes, *region_ids]),
        *(["occluded_area", "occlusion_fraction"] if occlusion else []),
        "contact_fraction_residues",
        *(f"{name}_contact_fraction_residues" for name in types),
        *(["occluded_area_residues"] if occlusion else []),
    ]
    lifetime_runs = {
        "mean_lifetime": ("mean_lifetime", "polymer"),
        **{f"{name}_mean_lifetime": ("mean_lifetime", name) for name in types},
        "lifetime_events": ("n_events", "polymer"),
        "censored_fraction": ("censored_fraction", "polymer"),
    }
    all_runs = [*fraction_runs, *lifetime_runs]
    run = run or "coverage"
    if run not in all_runs:
        raise ProtocolError(
            f"contacts: no result named {run!r}.", hint=f"Use --run with one of {all_runs}."
        )
    if run in lifetime_runs:
        life = {key: value for key, value in options.items() if key != "types"}
        if not occlusion:
            life["method"] = method
        if tolerance > 0:
            life["tolerance_ps"] = tolerance
        table = study.per_replicate(
            functions.contact_lifetimes,
            select(protein, allow_empty=True),
            select(polymer, allow_empty=True),
            unit=None,
            labels=["polymer", *types],
            name="contact_lifetimes",
            recompute=recompute,
            output_dir=output_dir,
            parts=list(functions.LIFETIME_PARTS),
            types=types,
            **life,
        )
        part, group = lifetime_runs[run]
        values = table[part].over_labels(lambda v: float(v[0]), run, labels=[group])
        values.unit, values.bounds = {
            "mean_lifetime": ("ns", (0.0, None)),
            "n_events": (None, (0.0, None)),
            "censored_fraction": (None, (0.0, 1.0)),
        }[part]
        report = _report_skipping(values, study, skipped, "contacts")
        no_events = [
            f"{label} replicate {row[0]}"
            for label, table_rows in values.rows.items()
            for row in table_rows
            if not np.isfinite(row[1])
        ]
        if no_events:
            report.warnings.append(
                f"contacts: {', '.join(no_events)} have no contact event for {group}, so {run} "
                "is undefined (nan) there."
            )
        report.provenance.settings = {
            **{
                key: value
                for key, value in settings.items()
                if key not in CONTACT_METHOD_SETTINGS[other]
            },
            "protein_selection": protein,
            "polymer_selection": polymer,
            "polymer_types_found": types,
            "unmeasured_residues": unmeasured,
        }
        if plots:
            folder = _figures_dir(output_dir, "contacts")
            values.plot(folder, f"contacts_{run}_comparison", title=run.replace("_", " "))
            report.provenance.output_paths["figures"] = str(folder)
        return report.model_copy(update={"analysis": "contacts", "run": run, "all_runs": all_runs})
    rows = study.per_replicate(
        function,
        select(protein, allow_empty=True),
        select(polymer, allow_empty=True),
        unit=None,
        labels=lambda u: [int(r.resid) for r in measured(u.select_atoms(protein).residues)],
        name=name,
        recompute=recompute,
        output_dir=output_dir,
        bounds=(0.0, 1.0),
        parts=[*parts, *type_parts],
        **options,
    )
    profile = rows["contact_fraction"]
    totals = {
        "coverage": profile.over_labels(lambda v: float(np.mean(np.asarray(v) > 0)), "coverage"),
        "mean_contact_fraction": profile.over_labels("mean", "mean_contact_fraction"),
    }
    totals |= {
        f"{name}_contact_fraction": rows[f"{name}_contact_fraction"].over_labels(
            "mean", f"{name}_contact_fraction"
        )
        for name in types
    }
    totals |= {
        f"{name}_contact_fraction": profile.over_labels(
            "mean", f"{name}_contact_fraction", labels=by_class[name]
        )
        for name in classes
    }
    totals |= {
        f"{name}_contact_fraction": profile.over_labels(
            "mean", f"{name}_contact_fraction", labels=ids
        )
        for name, ids in region_ids.items()
    }
    for values in totals.values():
        values.bounds = (0.0, 1.0)
    residue_runs = {"contact_fraction_residues": profile} | {
        f"{name}_contact_fraction_residues": rows[f"{name}_contact_fraction"] for name in types
    }
    if method == "occlusion":
        area = rows["occluded_area"]
        exposed = rows["exposed_area"]
        area.unit, area.bounds = "A^2", (0.0, None)
        totals["occluded_area"] = area.over_labels(lambda v: float(np.sum(v)), "occluded_area")
        totals["occluded_area"].unit, totals["occluded_area"].bounds = "A^2", (0.0, None)
        ratio = {}
        for label, table in area.rows.items():
            alone = {row[0]: float(np.sum(row[1])) for row in exposed.rows[label]}
            ratio[label] = [
                (
                    row[0],
                    float(np.sum(row[1])) / alone[row[0]] if alone[row[0]] > 0 else 0.0,
                    *row[2:],
                )
                for row in table
            ]
        totals["occlusion_fraction"] = ReplicateValues(
            area.source, "occlusion_fraction", None, True, ratio
        )
        residue_runs["occluded_area_residues"] = area
    if [*totals, *residue_runs] != fraction_runs:
        raise RuntimeError("contacts: the planned results differ from the computed ones.")
    runs = all_runs
    values = residue_runs[run] if run in residue_runs else totals[run]
    report = _report_skipping(values, study, skipped, "contacts")
    if unmeasured:
        report.warnings.append(
            f"contacts: {len(unmeasured)} residues of the protein selection have no maximum "
            f"ASA and are not measured by method=occlusion: {', '.join(unmeasured)}. They "
            "still cover their neighbours."
        )
    report.provenance.settings = {
        **{
            key: value
            for key, value in settings.items()
            if key not in CONTACT_METHOD_SETTINGS[other]
        },
        "protein_selection": protein,
        "polymer_selection": polymer,
        "polymer_types_found": types,
        "unmeasured_residues": unmeasured,
        "residues": {"classes": by_class, **region_ids},
    }
    if plots:
        folder = _figures_dir(output_dir, "contacts")
        plot_values(
            [totals[f"{name}_contact_fraction"] for name in classes],
            classes,
            folder,
            "contacts_class_bars",
            "Contact fraction by amino-acid class",
        )
        if run in residue_runs:
            name = run[: -len("_residues")]
            values.plot(
                folder, f"contacts_{name}_profile", f"Per-residue {name}", None, [], "Residue"
            )
            if report.pairwise:
                plot_differences(
                    values, report, folder, f"contacts_{name}_difference", None, None, "Residue"
                )
        else:
            values.plot(folder, f"contacts_{run}_comparison", title=run.replace("_", " "))
        report.provenance.output_paths["figures"] = str(folder)
    return report.model_copy(update={"analysis": "contacts", "run": run, "all_runs": runs})


def _residue_ids(study: Any, selection: str) -> list[int]:
    """Return the residue IDs ``selection`` picks, the same in every replicate of ``study``."""
    found = {
        tuple(int(r) for r in replicate.universe().select_atoms(selection).residues.resids)
        for condition in study
        for replicate in condition.replicates
    }
    if len(found) != 1 or not next(iter(found)):
        raise ProtocolError(
            f"The selection {selection!r} picks {'no' if found == {()} else 'different'} "
            "residues in the replicates.",
            hint="Choose core and region selections that pick the same residues in every replicate.",
        )
    return list(found.pop())


def _figures_dir(output_dir: Path | None, name: str) -> Path:
    """Return ``<output_dir>/figures/<name>``, with the current directory by default."""
    return Path(output_dir or Path.cwd()).expanduser().resolve() / "figures" / name


def _study(
    configs: Sequence[Path | str],
    labels: Sequence[str] | None,
    equilibration: str | None,
    replicates: Sequence[int] | None,
    stride: int = 1,
    data: dict[str, Path] | None = None,
    until: str | None = None,
) -> Any:
    """Build the Study of ``configs``, with the package default equilibration window."""
    from polyzymd.analyses.study import Study
    from polyzymd.config.analysis_settings import AnalysisDefaults

    paths = [Path(item).expanduser().resolve() for item in configs]
    return Study.from_configs(
        dict(zip(_labels(paths, labels), paths, strict=True)),
        equilibration=equilibration or AnalysisDefaults().equilibration_time,
        replicates=replicates,
        stride=stride,
        data=data,
        until=until,
    )


def _analyze_pairs(
    name: str,
    study: Any,
    settings: dict | None,
    run: str | None,
    *,
    recompute: bool,
    output_dir: Path | None,
    eq_check: bool,
    plots: bool = True,
) -> ProtocolReport:
    """Measure every pair of ``distances`` and report one result.

    ``pairs`` is a list of mappings with ``label``, ``selection_a``,
    ``selection_b`` and optionally ``threshold``, ``below_label`` and
    ``above_label``, or the path of a YAML or JSON file holding that list. A
    selection may be wrapped in ``midpoint(...)`` or ``com(...)``. Each pair
    is measured once with :func:`~polyzymd.analyses.functions.pair_distance`,
    and its results are named ``<label>`` for the mean distance and
    ``<label> <below_label>`` for the fraction of frames strictly below the
    pair's threshold, which defaults to ``threshold``, computed from the
    stored distance with :func:`~polyzymd.analyses.functions.all_below`.
    ``run`` picks the result to report, by default the first, and
    ``all_runs`` lists them all. With ``plots``, every pair's distance
    distribution with its threshold is drawn as ``distance_kde_<label>`` and
    every fraction as ``distance_fraction_<result>`` into
    ``<output_dir>/figures/<name>/``. ``distance_kde_panel`` stacks every
    pair's distribution in one figure and ``distance_threshold_bars`` groups
    every fraction.
    """
    import yaml

    from polyzymd.analyses import functions
    from polyzymd.analyses.shared.selections import parse_selection_string
    from polyzymd.analyses.timeseries import select

    settings = {**FUNCTION_ANALYSES[name], **(settings or {})}
    pairs = settings["pairs"]
    if isinstance(pairs, (str, Path)):
        try:
            pairs = yaml.safe_load(Path(pairs).expanduser().read_text())
        except (OSError, yaml.YAMLError) as exc:
            raise ProtocolError(
                f"{name}: cannot read the pairs file {pairs}: {exc}",
                hint="Give a YAML or JSON list of pairs.",
            ) from exc
    keys = {"label", "selection_a", "selection_b", "threshold", "below_label", "above_label"}
    if (
        not isinstance(pairs, list)
        or not pairs
        or not all(
            isinstance(pair, dict) and {"label", "selection_a", "selection_b"} <= set(pair) <= keys
            for pair in pairs
        )
    ):
        raise ProtocolError(
            f"{name} needs pairs, a list of mappings with label, selection_a and selection_b, "
            f"and optionally threshold, below_label and above_label; got {pairs!r}.",
            hint="Write the list to pairs.yaml and pass --set pairs=pairs.yaml.",
        )
    results, distances, thresholds = {}, [], []
    for pair in pairs:
        a, b = (parse_selection_string(str(pair[key])) for key in ("selection_a", "selection_b"))
        threshold = float(
            settings["threshold"] if pair.get("threshold") is None else pair["threshold"]
        )
        distance = study.timeseries(
            functions.pair_distance,
            select(a.selection),
            select(b.selection),
            mode_a=a.mode.value,
            mode_b=b.mode.value,
            pbc=bool(settings["use_pbc"]),
            unit="A",
            name=f"{name}_{pair['label']}",
            recompute=recompute,
            output_dir=output_dir,
            bounds=(0.0, None),
        )
        below = pair.get("below_label") or f"below {threshold:g} A"
        fraction = distance.transform(
            functions.all_below,
            unit=None,
            bounds=(0.0, 1.0),
            name=f"{distance.name}_{below}",
            thresholds=[threshold],
        )
        results[str(pair["label"])] = (distance, "mean", "mean_distance")
        results[f"{pair['label']} {below}"] = (fraction, "fraction", "fraction_below_threshold")
        distances.append(distance)
        thresholds.append(threshold)
    run = next(iter(results)) if run is None else run
    if run not in results:
        raise ProtocolError(
            f"{name}: no result named {run!r}.", hint=f"Use --run with one of {list(results)}."
        )
    series, how, metric = results[run]
    values = series.reduce(how, detect_equilibration=eq_check)
    values.metric = metric
    report = values.compare() if len(study) > 1 else values.summary()
    if plots:
        folder, prefix = _figures_dir(output_dir, name), "distance"
        from polyzymd.analyses.figures import plot_distributions, plot_values

        titles = [f"{pair['label']} distance" for pair in pairs]
        for pair, distance, threshold, title in zip(pairs, distances, thresholds, titles):
            distance.plot_distribution(threshold, folder, f"{prefix}_kde_{pair['label']}", title)
        plot_distributions(distances, thresholds, titles, folder, f"{prefix}_kde_panel", "distance")
        fractions = {}
        for key, (series, how, metric) in results.items():
            if how == "fraction":
                fractions[key] = series.reduce(how, detect_equilibration=False)
                fractions[key].metric = metric
                fractions[key].plot(folder, f"{prefix}_fraction_{key}", title=key)
        plot_values(
            list(fractions.values()),
            list(fractions),
            folder,
            f"{prefix}_threshold_bars",
            "Distance contact fractions",
        )
        report.provenance.output_paths["figures"] = str(folder)
    return report.model_copy(update={"analysis": name, "run": run, "all_runs": list(results)})


# Condition labels


def _labels(configs: Sequence[Path], labels: Sequence[str] | None) -> list[str]:
    """Pick a unique label per condition, defaulting to the config directory name."""
    if labels is not None:
        chosen = [str(label) for label in labels]
        if len(chosen) != len(configs):
            raise ProtocolError(
                f"Got {len(chosen)} label(s) for {len(configs)} config(s).",
                hint="Pass one --label per -c, in the same order, or none at all.",
            )
        if len(set(chosen)) != len(chosen):
            raise ProtocolError(
                "Condition labels must be unique.", hint="Give each --label a different name."
            )
        return chosen

    seen: dict[str, int] = {}
    unique: list[str] = []
    for path in configs:
        base = path.parent.name if path.parent.name not in {"", ".", ".."} else path.stem
        seen[base] = seen.get(base, 0) + 1
        unique.append(base if seen[base] == 1 else f"{base}_{seen[base]}")
    return unique


# Report rows


def _condition(label: str, values: Sequence[float]) -> ConditionReport:
    """Summarise one condition's replicate values as n, mean, sem and 95 percent interval.

    The statistics come from :func:`~polyzymd.analyses.shared.statistics.mean_sem_ci`.
    A condition without values gets n 0, a NaN mean and no sem or interval.
    """
    from polyzymd.analyses.shared.statistics import mean_sem_ci

    if not values:
        return ConditionReport(label=label, n_replicates=0, mean=float("nan"))
    stats = mean_sem_ci(values)
    limits = None if stats.ci_low is None else (stats.ci_low, stats.ci_high)
    return ConditionReport(
        label=label,
        n_replicates=len(values),
        mean=stats.mean,
        sem=stats.sem,
        ci95=limits,
        ci_method=stats.ci_method,
        replicate_values=list(values),
    )


def _difference_ci(
    values_a: Sequence[float], values_b: Sequence[float], test: str
) -> tuple[float, float] | None:
    """Return the 95 percent interval on ``mean(b) - mean(a)``.

    The interval comes from ``scipy.stats.ttest_ind``, the same call that runs
    the test, so it matches the variance assumption of the reported test: a
    pooled variance for Student's t and separate variances with
    Welch-Satterthwaite degrees of freedom for Welch's t [2]_. Tukey HSD gets no
    interval, because its simultaneous intervals are not computed here. The
    interval covers this one difference and carries no multiplicity correction,
    so a comparison can be non-significant after the correction while its
    interval excludes zero. ``None`` when a condition has fewer than two values
    or both have zero variance, where no interval can be estimated.
    """
    if test == "tukey_hsd" or len(values_a) < 2 or len(values_b) < 2:
        return None
    if len(set(values_a)) == 1 and len(set(values_b)) == 1:
        return None

    from scipy import stats

    result = stats.ttest_ind(values_b, values_a, equal_var=test != "welch_t")
    low, high = (float(limit) for limit in result.confidence_interval(0.95))
    if not (math.isfinite(low) and math.isfinite(high)) or low == high:
        return None
    return (low, high)


def _trend_line(trend: TrendReport) -> str:
    """One report line per trend test."""
    if not trend.testable:
        return f"trend {trend.factor}  not testable  n {trend.n_replicates}"
    return (
        f"trend {trend.factor}  slope {_num(trend.slope)}  ci95 {_interval(trend.slope_ci95)}"
        f"  p {_num(trend.p)}  p_adj {_num(trend.p_adjusted)}  r2 {_num(trend.r_squared)}"
        f"  n {trend.n_replicates}  conditions {len(trend.conditions)}"
    )


def _verdict(
    metric: str,
    unit: str | None,
    conditions: Sequence[ConditionReport],
    pairwise: Sequence[PairwiseReport],
) -> list[str]:
    """Write one sentence per comparison, or one per condition when there are none."""
    unit_text = f" {unit}" if unit else ""
    if not pairwise:
        return [
            f"{item.label} {metric} {_num(item.mean)}{unit_text} "
            + (
                f"(95% CI {_interval(item.ci95)}, n {item.n_replicates})"
                if item.ci95
                else f"(no interval, n {item.n_replicates})"
            )
            for item in conditions
        ]

    counts = {item.label: item.n_replicates for item in conditions}
    sentences = []
    for pair in pairwise:
        n_text = f"n {counts.get(pair.a, 0)} vs {counts.get(pair.b, 0)}"
        evidence = (
            f"delta {_signed(pair.delta)}{unit_text}, 95% CI {_interval(pair.delta_ci95)}, "
            f"p_adj {_num(pair.p_adjusted)}, p {_num(pair.p)}, {n_text}"
        )
        if not pair.testable:
            sentences.append(
                f"{VERDICT_NOT_TESTABLE}: {metric} for {pair.a} vs {pair.b} needs at least two "
                f"replicates per condition and a value that varies ({n_text})"
            )
        elif pair.p_adjusted is None:
            sentences.append(
                f"{VERDICT_NO_TEST} for {metric} between {pair.a} and {pair.b}; no "
                f"multiplicity-corrected p value was computed ({evidence})"
            )
        elif not pair.significant:
            sentences.append(
                f"{VERDICT_NO_DIFFERENCE} in {metric} between {pair.a} and {pair.b} ({evidence})"
            )
        else:
            word = (
                VERDICT_CHANGED
                if pair.delta == 0
                else (VERDICT_LARGER if pair.delta > 0 else VERDICT_SMALLER)
            )
            sentences.append(f"{pair.b} {word} {metric} than {pair.a} ({evidence})")
    return sentences


# Formatting helpers


def _num(value: float | None) -> str:
    """Format one number with four significant digits, or ``na``."""
    if value is None:
        return "na"
    return "nan" if math.isnan(value) else f"{value:.4g}"


def _signed(value: float | None) -> str:
    """Format one number with an explicit sign, or ``na``."""
    if value is None:
        return "na"
    return "nan" if math.isnan(value) else f"{value:+.4g}"


def _interval(limits: Sequence[float] | None) -> str:
    """Format an interval as ``low to high``, or ``na``."""
    return "na" if limits is None else f"{_num(limits[0])} to {_num(limits[1])}"


def _condition_line(condition: ConditionReport) -> str:
    """Render one condition on a single line."""
    shown = ", ".join(_num(value) for value in condition.replicate_values)
    entry = "" if condition.entry is None else f"{condition.entry}  "
    line = (
        f"{entry}{condition.label}  n {condition.n_replicates}  mean {_num(condition.mean)}"
        f"  sem {_num(condition.sem)}  ci95 {_interval(condition.ci95)}"
        f"  values {shown or 'none'}"
    )
    if condition.statistical_inefficiency:
        line += (
            f"  replicates {', '.join(str(index) for index in condition.replicates)}"
            f"  g {', '.join(_num(value) for value in condition.statistical_inefficiency)}"
            f"  n_eff {', '.join(_num(value) for value in condition.n_effective)}"
        )
    if condition.eq_detected_ns:
        line += f"  eq_detected {_num(max(condition.eq_detected_ns))} ns"
    return line


def _pairwise_line(pair: PairwiseReport) -> str:
    """Render one comparison on a single line."""
    if not pair.testable:
        flag = "not_testable"
    elif pair.p_adjusted is None:
        flag = "no_test"
    else:
        flag = "significant" if pair.significant else "not_significant"
    family = "" if pair.family_size is None else f"  family {pair.family_size}"
    entry = "" if pair.entry is None else f"{pair.entry}  "
    return (
        f"{entry}{pair.a} vs {pair.b}  delta {_signed(pair.delta)}  ci95 {_interval(pair.delta_ci95)}"
        f"  p {_num(pair.p)}  p_adj {_num(pair.p_adjusted)}  test {pair.test}"
        f"  correction {pair.correction}{family}  d {_num(pair.cohens_d)}  {flag}"
    )


def _labelled_pairwise_lines(pairwise: Sequence[PairwiseReport]) -> list[str]:
    """Summarise per-label comparisons: counts per condition, then the significant labels."""
    lines = []
    for b in dict.fromkeys(pair.b for pair in pairwise):
        rows = [pair for pair in pairwise if pair.b == b]
        tested = [pair for pair in rows if pair.p_adjusted is not None]
        family = tested[0].family_size if tested else None
        found = {
            "lower": [pair for pair in tested if pair.significant and pair.delta < 0],
            "higher": [pair for pair in tested if pair.significant and pair.delta > 0],
        }
        lines.append(
            f"{rows[0].a} vs {b}  labels {len(rows)}  tested {len(tested)}  family "
            f"{family if family is not None else 'na'}  test {rows[0].test}  correction "
            f"{rows[0].correction}  lower {len(found['lower'])}  higher {len(found['higher'])}"
        )
        for side, items in found.items():
            if items:
                listed = ", ".join(
                    f"{pair.entry} delta {_signed(pair.delta)} p_adj {_num(pair.p_adjusted)}"
                    for pair in items
                )
                lines.append(f"{rows[0].a} vs {b}  {side}: {listed}")
    return lines
