"""Polymer bridging analysis plugin.

Quantifies oligomer-level multisite attachment to the protein surface directly
from trajectories. The core question is whether a single short oligomer
typically makes one-site or multi-site protein attachments.

**Status:** Experimental. Contributed as a proof-of-concept extensibility
exercise. Metric definitions, chemistry-aware profiling outputs, and
interpretation are subject to change.

Observation model
-----------------
One observation = one polymer fragment in one trajectory frame that contacts
at least one protein residue within the distance cutoff. All statistics are
computed over this observation population.

Primary metrics
---------------
- mean contacts per contacting oligomer
- fraction of contacting oligomer observations that are multisite
- fraction of contacting oligomer observations with high valency (3+ residues)

Chemistry-aware outputs (experimental)
--------------------------------------
- anchor protein class probabilities
- peripheral and all-multivalent protein class probabilities
- polymer contact/anchor monomer probabilities
- anchor-to-peripheral class matrices (row-normalized)
- polymer-anchor to protein-anchor matrices (row-normalized)
- ordered fragment signature probabilities (top 10)

Key definitions
---------------
- *Multisite*: effective eligible valency > 1. With positive
  ``min_ca_distance_angstrom``, eligibility is based on dynamic per-frame
  CA distances among contacted protein residues.
- *Anchor*: the polymer-protein residue pair with the minimum atom distance
  within a multivalent observation.
- *Peripheral*: all other eligible contacted residues in the observation.

Notes
-----
Interpretation is descriptive and experimental — observed frequencies over the
observation population, not normalized enrichment and not proof of mechanism.
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field, field_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ANOVAResult,
    BasePlotSettings,
    ComparisonResult,
    Condition,
    ConditionSummary,
    MetricValue,
    PairwiseResult,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.polymer_bridging._plot_settings import PolymerBridgingPlotSettings
from polyzymd.analyses.shared import apply_axis_style, get_colors, get_output_path, save_figure
from polyzymd.analyses.shared.config_hash import settings_fingerprint
from polyzymd.analyses.shared.plotting import (
    scatter_replicate_values,
    scatter_stacked_segment_replicates,
)
from polyzymd.analyses.stats import anova_test, pairwise_comparisons
from polyzymd.core.experimental import prefix_experimental_output

logger = logging.getLogger("polyzymd.analyses.polymer_bridging")

ResiduePairDistances = dict[tuple[int, int], float]

PROTEIN_GROUP_ORDER = [
    "aromatic",
    "charged_positive",
    "charged_negative",
    "polar",
    "nonpolar",
    "unknown",
]
EXPERIMENTAL_CHEMISTRY_FEATURES = ("polymer_bridging_chemistry",)


@dataclass(frozen=True)
class PolymerBridgingObservation:
    """One fragment-frame observation with chemistry metadata.

    Represents a single polymer fragment in a single trajectory frame that
    contacts at least one protein residue. This is the fundamental unit of
    data from which all bridging statistics are computed.

    Attributes
    ----------
    protein_residues : set[int]
        Residue IDs of all protein residues contacted by this fragment.
    protein_resnames : dict[int, str]
        Mapping from protein residue ID to residue name (3-letter code).
    protein_groups : dict[int, str]
        Mapping from protein residue ID to amino acid class
        (aromatic, charged_positive, etc.).
    contacting_polymer_resids : set[int]
        Residue IDs of polymer monomers that participate in contacts.
    polymer_resnames : dict[int, str]
        Mapping from polymer residue ID to monomer name.
    fragment_signature : tuple[str, ...]
        Ordered sequence of all monomer names in the fragment (full chain).
    ca_distances : ResiduePairDistances
        Frame-wise C-alpha distances for all pairs of contacted protein
        residues. Keys are ``(min_resid, max_resid)`` tuples.
    pair_min_distances : dict[tuple[int, int], float]
        Minimum atom-level distance for each (polymer_resid, protein_resid)
        pair. Used to determine the anchor.
    """

    protein_residues: set[int]
    protein_resnames: dict[int, str]
    protein_groups: dict[int, str]
    contacting_polymer_resids: set[int]
    polymer_resnames: dict[int, str]
    fragment_signature: tuple[str, ...]
    ca_distances: ResiduePairDistances
    pair_min_distances: dict[tuple[int, int], float]


FrameContactObservation = PolymerBridgingObservation | tuple[set[int], ResiduePairDistances]


class PolymerBridgingSettings(BaseModel):
    """Settings for oligomer bridging analysis.

    Parameters
    ----------
    protein_selection : str
        MDAnalysis atom selection string for the protein.
        Default ``"protein"``.
    polymer_selection : str
        MDAnalysis atom selection string for the polymer chain(s).
        Default ``"chainID C"`` (PolyzyMD chain convention).
    cutoff : float
        Contact distance cutoff in Angstroms. Atom pairs within this
        distance are considered in contact. Default ``4.5``.
    min_ca_distance_angstrom : float
        Minimum frame-wise C-alpha distance (Angstroms) between contacted
        protein residues for an observation to count as multisite. Set to
        ``0.0`` (default) to disable geometric filtering, or to a positive
        value (e.g., ``8.0``) to require spatially separated contacts.
        Must be >= 0.
    """

    protein_selection: str = Field(
        default="protein",
        description="Protein atom selection used for oligomer contact detection",
    )
    polymer_selection: str = Field(
        default="chainID C",
        description="Polymer selection used for compatibility with contacts filtering",
    )
    cutoff: float = Field(
        default=4.5,
        description="Contact cutoff in Angstroms for oligomer-protein contact detection",
    )
    min_ca_distance_angstrom: float = Field(
        default=0.0,
        description="Minimum frame-wise CA-CA distance between contacted residues to count as separated",
    )

    @field_validator("min_ca_distance_angstrom", mode="after")
    @classmethod
    def validate_distance(cls, v: float) -> float:
        if v < 0:
            raise ValueError("min_ca_distance_angstrom must be >= 0")
        return v


class PolymerBridgingReplicateResult(BaseModel):
    """Per-replicate oligomer bridging summary.

    Contains both the primary scalar metrics (multisite fraction, mean
    contacts, high-valency fraction) and the chemistry-aware probability
    distributions computed from one replicate trajectory.

    Attributes
    ----------
    replicate : int
        Replicate index.
    n_frames : int
        Number of trajectory frames analyzed (after equilibration).
    timestep_ps : float
        Trajectory timestep in picoseconds.
    min_ca_distance_angstrom : float
        The C-alpha distance threshold used for eligibility filtering.
    contacting_observations : int
        Total number of fragment-frame observations with at least one contact.
    multisite_observations : int
        Number of observations with eligible valency > 1.
    high_valency_observations : int
        Number of observations with eligible valency >= 3.
    mean_contacts_per_contacting_oligomer : float
        Average eligible valency across all contacting observations.
    multisite_fraction : float
        Fraction of contacting observations that are multisite.
    high_valency_fraction : float
        Fraction of contacting observations with valency >= 3.
    valency_probabilities : dict[str, float]
        Probability distribution over valency bins: ``"1"``, ``"2"``, ``"3+"``.
    anchor_protein_group_probabilities : dict[str, float]
        Amino acid class frequency of anchor residues in multivalent events.
    peripheral_protein_group_probabilities : dict[str, float]
        Amino acid class frequency of peripheral residues in multivalent events.
    multivalent_protein_group_probabilities : dict[str, float]
        Amino acid class frequency of all eligible residues in multivalent events.
    polymer_contact_type_probabilities : dict[str, float]
        Monomer-type frequency among all contacting polymer residues.
    polymer_anchor_type_probabilities : dict[str, float]
        Monomer-type frequency of the polymer anchor in multivalent events.
    anchor_to_peripheral_group_matrix : dict[str, dict[str, float]]
        Row-normalized matrix: anchor protein class -> peripheral protein class.
    polymer_anchor_to_protein_anchor_matrix : dict[str, dict[str, float]]
        Row-normalized matrix: polymer anchor monomer -> protein anchor class.
    fragment_signature_probabilities : dict[str, float]
        Top-10 ordered fragment signatures by frequency in multivalent events.
    """

    settings_fingerprint: str | None = None
    protein_selection: str | None = None
    polymer_selection: str | None = None
    cutoff: float | None = None
    config_hash: str = "unknown"
    equilibration_time: float = 0.0
    equilibration_unit: str = "ns"
    replicate: int
    n_frames: int
    timestep_ps: float
    min_ca_distance_angstrom: float
    contacting_observations: int
    multisite_observations: int
    high_valency_observations: int
    mean_contacts_per_contacting_oligomer: float
    multisite_fraction: float
    high_valency_fraction: float
    valency_probabilities: dict[str, float]
    anchor_protein_group_probabilities: dict[str, float] = Field(default_factory=dict)
    peripheral_protein_group_probabilities: dict[str, float] = Field(default_factory=dict)
    multivalent_protein_group_probabilities: dict[str, float] = Field(default_factory=dict)
    polymer_contact_type_probabilities: dict[str, float] = Field(default_factory=dict)
    polymer_anchor_type_probabilities: dict[str, float] = Field(default_factory=dict)
    anchor_to_peripheral_group_matrix: dict[str, dict[str, float]] = Field(default_factory=dict)
    polymer_anchor_to_protein_anchor_matrix: dict[str, dict[str, float]] = Field(
        default_factory=dict
    )
    fragment_signature_probabilities: dict[str, float] = Field(default_factory=dict)

    def save(self, path: Path | str) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.model_dump_json(indent=2))
        return path

    @classmethod
    def load(cls, path: Path | str) -> "PolymerBridgingReplicateResult":
        return cls.model_validate_json(Path(path).read_text())


class PolymerBridgingAggregatedResult(BaseModel):
    """Aggregated oligomer bridging result for one condition.

    Contains cross-replicate mean +/- SEM for all metrics and
    chemistry-aware probability distributions. Per-replicate values
    are preserved for downstream statistical testing.

    Attributes
    ----------
    n_replicates : int
        Number of replicates aggregated.
    replicates : list[int]
        Replicate indices.
    min_ca_distance_angstrom : float
        C-alpha distance threshold used.
    mean_contacts_per_contacting_oligomer : float
        Cross-replicate mean of per-replicate mean valency.
    mean_contacts_sem : float
        Standard error of the mean for valency.
    multisite_fraction : float
        Cross-replicate mean of per-replicate multisite fraction.
    multisite_fraction_sem : float
        Standard error of the mean for multisite fraction.
    high_valency_fraction : float
        Cross-replicate mean of per-replicate high-valency fraction.
    high_valency_fraction_sem : float
        Standard error of the mean for high-valency fraction.
    """

    n_replicates: int
    replicates: list[int]
    min_ca_distance_angstrom: float
    mean_contacts_per_contacting_oligomer: float
    mean_contacts_sem: float
    multisite_fraction: float
    multisite_fraction_sem: float
    high_valency_fraction: float
    high_valency_fraction_sem: float
    mean_contacts_per_contacting_oligomer_replicates: list[float]
    multisite_fraction_replicates: list[float]
    high_valency_fraction_replicates: list[float]
    valency_probabilities_mean: dict[str, float]
    valency_probabilities_sem: dict[str, float]
    valency_probabilities_per_replicate: dict[str, list[float]] = Field(default_factory=dict)
    anchor_protein_group_probabilities_mean: dict[str, float] = Field(default_factory=dict)
    anchor_protein_group_probabilities_sem: dict[str, float] = Field(default_factory=dict)
    anchor_protein_group_probabilities_per_replicate: dict[str, list[float]] = Field(
        default_factory=dict
    )
    peripheral_protein_group_probabilities_mean: dict[str, float] = Field(default_factory=dict)
    peripheral_protein_group_probabilities_sem: dict[str, float] = Field(default_factory=dict)
    peripheral_protein_group_probabilities_per_replicate: dict[str, list[float]] = Field(
        default_factory=dict
    )
    multivalent_protein_group_probabilities_mean: dict[str, float] = Field(default_factory=dict)
    multivalent_protein_group_probabilities_sem: dict[str, float] = Field(default_factory=dict)
    multivalent_protein_group_probabilities_per_replicate: dict[str, list[float]] = Field(
        default_factory=dict
    )
    polymer_contact_type_probabilities_mean: dict[str, float] = Field(default_factory=dict)
    polymer_contact_type_probabilities_sem: dict[str, float] = Field(default_factory=dict)
    polymer_contact_type_probabilities_per_replicate: dict[str, list[float]] = Field(
        default_factory=dict
    )
    polymer_anchor_type_probabilities_mean: dict[str, float] = Field(default_factory=dict)
    polymer_anchor_type_probabilities_sem: dict[str, float] = Field(default_factory=dict)
    polymer_anchor_type_probabilities_per_replicate: dict[str, list[float]] = Field(
        default_factory=dict
    )
    anchor_to_peripheral_group_matrix_mean: dict[str, dict[str, float]] = Field(
        default_factory=dict
    )
    anchor_to_peripheral_group_matrix_sem: dict[str, dict[str, float]] = Field(default_factory=dict)
    polymer_anchor_to_protein_anchor_matrix_mean: dict[str, dict[str, float]] = Field(
        default_factory=dict
    )
    polymer_anchor_to_protein_anchor_matrix_sem: dict[str, dict[str, float]] = Field(
        default_factory=dict
    )
    fragment_signature_probabilities_mean: dict[str, float] = Field(default_factory=dict)
    fragment_signature_probabilities_sem: dict[str, float] = Field(default_factory=dict)
    fragment_signature_probabilities_per_replicate: dict[str, list[float]] = Field(
        default_factory=dict
    )
    settings_fingerprint: str | None = None
    protein_selection: str | None = None
    polymer_selection: str | None = None
    cutoff: float | None = None
    config_hash: str = "unknown"
    equilibration_time: float = 0.0
    equilibration_unit: str = "ns"

    def save(self, path: Path | str) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.model_dump_json(indent=2))
        return path

    @classmethod
    def load(cls, path: Path | str) -> "PolymerBridgingAggregatedResult":
        return cls.model_validate_json(Path(path).read_text())


class PolymerBridgingAnalysis(Analysis):
    """Oligomer-level multisite protein attachment analysis.

    This plugin computes per-fragment, per-frame polymer bridging metrics
    directly from MD trajectories. It implements the full Analysis lifecycle:
    ``compute_replicate`` → ``aggregate`` → ``compare`` → ``plot`` → ``format``.

    **Status:** Experimental. Chemistry-aware outputs are labeled
    ``polymer_bridging_chemistry`` in the experimental feature system.

    Class Variables
    ---------------
    name : str
        Plugin name: ``"polymer_bridging"``.
    Settings : type
        ``PolymerBridgingSettings`` — Pydantic model for YAML configuration.
    AggregatedResultClass : type
        ``PolymerBridgingAggregatedResult`` — used by the default result loader.
    aliases : tuple[str, ...]
        ``("bridging",)`` — CLI shorthand.
    dependencies : tuple[str, ...]
        Empty — no prerequisite analyses.
    min_replicates : int
        ``2`` — minimum replicates for comparison statistics.

    Notes
    -----
    Conditions without polymer atoms are automatically filtered via
    ``filter_conditions()``. The comparison uses NaN-safe pairwise t-tests
    and ANOVA across three primary metrics.
    """

    name: ClassVar[str] = "polymer_bridging"
    Settings: ClassVar[type] = PolymerBridgingSettings
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = PolymerBridgingPlotSettings
    AggregatedResultClass: ClassVar[type] = PolymerBridgingAggregatedResult
    ReplicateResultClass: ClassVar[type | None] = PolymerBridgingReplicateResult
    aliases: ClassVar[tuple[str, ...]] = ("bridging",)
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2

    @staticmethod
    def _make_settings_cache_tag(settings: BaseModel) -> str:
        """Build a short cache tag from analysis settings.

        Parameters
        ----------
        settings : BaseModel
            Analysis settings model.

        Returns
        -------
        str
            First 8 hex characters from shared settings fingerprinting.
        """
        return settings_fingerprint(settings)

    @staticmethod
    def _make_window_cache_tag(equilibration: str) -> str:
        """Build a short cache tag from the analysis window.

        Parameters
        ----------
        equilibration : str
            Requested equilibration time string.

        Returns
        -------
        str
            Window tag safe for cache filenames.
        """

        from polyzymd.analyses.shared.loader import parse_time_string

        eq_value, eq_unit = parse_time_string(equilibration)
        return f"eq{eq_value:g}{eq_unit}"

    @staticmethod
    def _cache_matches_window(result: Any, equilibration: str) -> bool:
        """Return whether a cached result matches the requested analysis window.

        Parameters
        ----------
        result : Any
            Loaded cache result.
        equilibration : str
            Requested equilibration window from the comparison context.

        Returns
        -------
        bool
            ``True`` when stored equilibration metadata matches the requested
            window, otherwise ``False``.
        """

        import math

        from polyzymd.analyses.shared.loader import convert_time, parse_time_string

        expected_time, expected_unit = parse_time_string(equilibration)
        stored_time = getattr(result, "equilibration_time", None)
        stored_unit = getattr(result, "equilibration_unit", None)
        if stored_time is None or stored_unit is None:
            return False
        try:
            stored_time_ps = convert_time(float(stored_time), str(stored_unit), "ps")
            expected_time_ps = convert_time(float(expected_time), str(expected_unit), "ps")
        except (TypeError, ValueError):
            return False
        return math.isclose(stored_time_ps, expected_time_ps, rel_tol=0.0, abs_tol=1.0e-9)

    @staticmethod
    def _cache_has_settings_proof(result: Any, settings: BaseModel) -> bool:
        """Return whether cache content proves polymer bridging settings.

        Parameters
        ----------
        result : Any
            Loaded polymer bridging result.
        settings : BaseModel
            Active polymer bridging settings.

        Returns
        -------
        bool
            ``True`` when explicit content fields match settings without using
            filename-derived identity.
        """

        protein_selection = getattr(result, "protein_selection", None)
        polymer_selection = getattr(result, "polymer_selection", None)
        cutoff = getattr(result, "cutoff", None)
        min_ca_distance = getattr(result, "min_ca_distance_angstrom", None)
        try:
            cutoff_matches = abs(float(cutoff) - float(settings.cutoff)) <= 1e-6
            min_ca_matches = (
                abs(float(min_ca_distance) - float(settings.min_ca_distance_angstrom)) <= 1e-6
            )
        except (TypeError, ValueError):
            return False
        return (
            str(protein_selection).strip() == str(settings.protein_selection).strip()
            and str(polymer_selection).strip() == str(settings.polymer_selection).strip()
            and cutoff_matches
            and min_ca_matches
        )

    @staticmethod
    def _attach_settings_proof(result: Any, settings: BaseModel) -> None:
        """Attach current explicit settings proof fields to a result.

        Parameters
        ----------
        result : Any
            Polymer bridging result to update in place.
        settings : BaseModel
            Active polymer bridging settings.
        """

        result.settings_fingerprint = settings_fingerprint(settings)
        result.protein_selection = settings.protein_selection
        result.polymer_selection = settings.polymer_selection
        result.cutoff = float(settings.cutoff)

    def _cache_matches_context(
        self,
        result: Any,
        *,
        settings: BaseModel,
        equilibration: str,
        sim_config: Any,
        replicates: Sequence[int] | None = None,
        allow_replicate_subset: bool = False,
        source: Path | None = None,
    ) -> bool:
        """Return whether a cached aggregate matches the active context.

        Parameters
        ----------
        result : Any
            Loaded polymer-bridging result.
        settings : BaseModel
            Active analysis settings.
        equilibration : str
            Active equilibration/window request.
        sim_config : Any
            Active condition simulation configuration.
        replicates : Sequence[int] or None, optional
            Requested replicate tuple for aggregate cache validation.
        allow_replicate_subset : bool, optional
            Whether a stored subset of the requested replicates is accepted.
            This is used only when consuming finalized aggregates.
        source : Path or None, optional
            Cache source used for diagnostics.

        Returns
        -------
        bool
            ``True`` when settings, window, and known config identity match.
        """

        from polyzymd.analyses.shared.config_hash import validate_config_hash

        if not self._cache_matches_window(result, equilibration):
            logger.warning(
                "polymer_bridging: ignoring cache with mismatched equilibration window%s",
                f" at {source}" if source is not None else "",
            )
            return False

        stored_fingerprint = getattr(result, "settings_fingerprint", None)
        if stored_fingerprint is None:
            stored_fingerprint = getattr(result, "settings_fp", None)
        if stored_fingerprint is None:
            if self._cache_has_settings_proof(result, settings):
                self._attach_settings_proof(result, settings)
            else:
                logger.warning(
                    "polymer_bridging: ignoring cache without embedded settings fingerprint or "
                    "complete settings proof%s",
                    f" at {source}" if source is not None else "",
                )
                return False
        elif str(stored_fingerprint) != settings_fingerprint(settings):
            logger.warning(
                "polymer_bridging: ignoring cache with settings fingerprint mismatch%s: "
                "stored=%s, current=%s",
                f" at {source}" if source is not None else "",
                stored_fingerprint,
                settings_fingerprint(settings),
            )
            return False

        stored_hash = getattr(result, "config_hash", None)
        if stored_hash not in (None, "", "unknown") and not validate_config_hash(
            str(stored_hash),
            sim_config,
        ):
            return False

        if replicates is not None and not self._cache_matches_replicates(
            result,
            replicates,
            allow_subset=allow_replicate_subset,
            source=source,
        ):
            return False

        return True

    @staticmethod
    def _cache_matches_replicate_id(
        result: Any,
        replicate: int,
        *,
        source: Path | None = None,
    ) -> bool:
        """Return whether a cached per-replicate result matches the request.

        Parameters
        ----------
        result : Any
            Loaded per-replicate result.
        replicate : int
            Requested replicate ID.
        source : Path or None, optional
            Cache source path used for diagnostics.

        Returns
        -------
        bool
            ``True`` when the stored replicate ID equals ``replicate``.
        """

        stored_replicate = getattr(result, "replicate", None)
        try:
            stored = int(stored_replicate)
        except (TypeError, ValueError):
            logger.warning(
                "polymer_bridging: ignoring replicate cache without valid replicate identity%s",
                f" at {source}" if source is not None else "",
            )
            return False
        if stored != int(replicate):
            logger.warning(
                "polymer_bridging: ignoring replicate cache for replicate %d; requested %d%s",
                stored,
                int(replicate),
                f" at {source}" if source is not None else "",
            )
            return False
        return True

    @staticmethod
    def _align_replicate_results(
        results: Sequence[Any],
        replicates: Sequence[int],
    ) -> list[Any]:
        """Return replicate results ordered to match the requested IDs.

        Parameters
        ----------
        results : Sequence[Any]
            Per-replicate results returned by the orchestrator.
        replicates : Sequence[int]
            Requested replicate IDs from the aggregate context.

        Returns
        -------
        list[Any]
            Results ordered by ``replicates``.
        """

        requested = tuple(int(rep) for rep in replicates)
        by_replicate: dict[int, Any] = {}
        for result in results:
            stored_replicate = getattr(result, "replicate", None)
            try:
                replicate = int(stored_replicate)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    "polymer_bridging aggregate input lacks a valid replicate ID"
                ) from exc
            if replicate in by_replicate:
                raise ValueError(
                    f"polymer_bridging aggregate input duplicates replicate {replicate}"
                )
            by_replicate[replicate] = result

        requested_set = set(requested)
        stored_set = set(by_replicate)
        if stored_set != requested_set:
            raise ValueError(
                "polymer_bridging aggregate input replicate IDs do not match requested "
                f"replicates: stored={sorted(stored_set)}, requested={sorted(requested_set)}"
            )
        return [by_replicate[replicate] for replicate in requested]

    def _cache_matches_replicates(
        self,
        result: Any,
        replicates: Sequence[int],
        *,
        allow_subset: bool = False,
        source: Path | None = None,
    ) -> bool:
        """Return whether cached aggregate replicate identity matches the request.

        Parameters
        ----------
        result : Any
            Loaded polymer-bridging aggregate result.
        replicates : Sequence[int]
            Requested replicate numbers for the active aggregate.
        allow_subset : bool, optional
            Whether stored replicate identity may be a non-empty subset of the
            requested set.
        source : Path or None, optional
            Cache source path used for diagnostics.

        Returns
        -------
        bool
            ``True`` when the cached result declares the requested replicate
            set, or a subset when explicitly allowed. Missing or malformed
            identity is treated as incompatible because the aggregate cannot be
            proven current.
        """

        requested = tuple(sorted(int(rep) for rep in replicates))
        stored_replicates = getattr(result, "replicates", None)
        if not stored_replicates:
            logger.warning(
                "polymer_bridging: ignoring cache without replicate identity%s",
                f" at {source}" if source is not None else "",
            )
            return False

        try:
            stored = tuple(sorted(int(rep) for rep in stored_replicates))
        except (TypeError, ValueError):
            logger.warning(
                "polymer_bridging: ignoring cache with invalid replicate identity%s",
                f" at {source}" if source is not None else "",
            )
            return False

        if len(set(stored)) != len(stored):
            logger.warning(
                "polymer_bridging: ignoring cache with duplicate replicate identity%s",
                f" at {source}" if source is not None else "",
            )
            return False

        is_allowed_subset = allow_subset and set(stored).issubset(set(requested))
        if stored != requested and not is_allowed_subset:
            logger.warning(
                "polymer_bridging: ignoring cache for replicates %s; requested %s%s",
                stored,
                requested,
                f" at {source}" if source is not None else "",
            )
            return False

        stored_count = len(stored)
        if allow_subset and stored_count < self.min_replicates:
            logger.warning(
                "polymer_bridging: ignoring cache with %d replicates below minimum %d%s",
                stored_count,
                self.min_replicates,
                f" at {source}" if source is not None else "",
            )
            return False

        declared_count = getattr(result, "n_replicates", None)
        try:
            declared_count = int(declared_count)
        except (TypeError, ValueError):
            logger.warning(
                "polymer_bridging: ignoring cache with invalid n_replicates%s",
                f" at {source}" if source is not None else "",
            )
            return False
        if declared_count != stored_count:
            logger.warning(
                "polymer_bridging: ignoring cache with n_replicates=%d but %d stored "
                "replicate IDs%s",
                declared_count,
                stored_count,
                f" at {source}" if source is not None else "",
            )
            return False

        mismatched_vectors = self._replicate_vector_length_mismatches(result, stored_count)
        if mismatched_vectors:
            logger.warning(
                "polymer_bridging: ignoring cache with per-replicate vector length mismatch "
                "for %s%s",
                ", ".join(mismatched_vectors),
                f" at {source}" if source is not None else "",
            )
            return False

        return True

    @staticmethod
    def _replicate_vector_length_mismatches(result: Any, expected_count: int) -> list[str]:
        """Return per-replicate vector fields whose lengths are inconsistent.

        Parameters
        ----------
        result : Any
            Loaded polymer-bridging aggregate result.
        expected_count : int
            Number of stored replicate identities.

        Returns
        -------
        list[str]
            Field labels for vectors that do not match ``expected_count``.
        """

        mismatches: list[str] = []
        vector_fields = (
            "mean_contacts_per_contacting_oligomer_replicates",
            "multisite_fraction_replicates",
            "high_valency_fraction_replicates",
        )
        for field_name in vector_fields:
            values = getattr(result, field_name, None)
            if isinstance(values, (list, tuple)) and len(values) != expected_count:
                mismatches.append(field_name)

        mapping_fields = (
            "valency_probabilities_per_replicate",
            "anchor_protein_group_probabilities_per_replicate",
            "peripheral_protein_group_probabilities_per_replicate",
            "multivalent_protein_group_probabilities_per_replicate",
            "polymer_contact_type_probabilities_per_replicate",
            "polymer_anchor_type_probabilities_per_replicate",
            "fragment_signature_probabilities_per_replicate",
        )
        for field_name in mapping_fields:
            mapping = getattr(result, field_name, {}) or {}
            if not isinstance(mapping, dict):
                continue
            for key, values in mapping.items():
                if isinstance(values, (list, tuple)) and len(values) != expected_count:
                    mismatches.append(f"{field_name}[{key}]")

        return mismatches

    def _replicate_cache_path(
        self,
        output_dir: Path,
        settings: BaseModel,
        equilibration: str,
    ) -> Path:
        """Return the window-aware replicate sidecar path.

        Parameters
        ----------
        output_dir : Path
            Replicate output directory.
        settings : BaseModel
            Active analysis settings.
        equilibration : str
            Active equilibration/window request.

        Returns
        -------
        Path
            Fingerprinted sidecar path for the replicate cache.
        """

        settings_tag = self._make_settings_cache_tag(settings)
        window_tag = self._make_window_cache_tag(equilibration)
        return output_dir / f"polymer_bridging_{window_tag}_s{settings_tag}.json"

    def _aggregate_cache_path(
        self,
        output_dir: Path,
        settings: BaseModel,
        equilibration: str,
        replicates: Sequence[int],
    ) -> Path:
        """Return the window-aware aggregate sidecar path.

        Parameters
        ----------
        output_dir : Path
            Aggregated output directory.
        settings : BaseModel
            Active analysis settings.
        equilibration : str
            Active equilibration/window request.
        replicates : Sequence[int]
            Replicates represented by the aggregate.

        Returns
        -------
        Path
            Fingerprinted sidecar path for the aggregate cache.
        """

        settings_tag = self._make_settings_cache_tag(settings)
        window_tag = self._make_window_cache_tag(equilibration)
        rep_str = self._format_replicate_range(replicates)
        return output_dir / f"polymer_bridging_{window_tag}_{rep_str}_s{settings_tag}.json"

    def _load_validated_aggregated_result(
        self,
        aggregated_dir: Path,
        *,
        settings: BaseModel,
        equilibration: str,
        replicates: Sequence[int],
        sim_config: Any,
        recompute: bool,
        allow_replicate_subset: bool = False,
    ) -> Any | None:
        """Load an aggregate only when cache identity metadata matches.

        Parameters
        ----------
        aggregated_dir : Path
            Directory containing aggregated cache files.
        settings : BaseModel
            Active analysis settings.
        equilibration : str
            Active equilibration/window request.
        replicates : Sequence[int]
            Replicates expected in the aggregate.
        sim_config : Any
            Active condition simulation configuration.
        recompute : bool
            Whether cache loading should be skipped.
        allow_replicate_subset : bool, optional
            Whether finalized aggregates with successful replicate subsets may
            be loaded.

        Returns
        -------
        Any or None
            Valid aggregated result, or ``None`` when no compatible cache is
            available.
        """

        candidates = [
            self._aggregate_cache_path(aggregated_dir, settings, equilibration, replicates),
            self.aggregate_result_path(aggregated_dir),
        ]
        saw_json = aggregated_dir.exists() and any(aggregated_dir.glob("*.json"))

        for candidate in candidates:
            cached = self._check_cache(
                PolymerBridgingAggregatedResult,
                candidate,
                recompute=recompute,
                sim_config=None,
                settings=None,
            )
            if cached is None:
                continue
            if self._cache_matches_context(
                cached,
                settings=settings,
                equilibration=equilibration,
                sim_config=sim_config,
                replicates=replicates,
                allow_replicate_subset=allow_replicate_subset,
                source=candidate,
            ):
                return cached

        if saw_json:
            return None

        return self._load_aggregated_result(aggregated_dir)

    @staticmethod
    def _compute_frame_contacts(*args, **kwargs):
        """Compatibility facade for direct frame-contact computation."""

        return _compute_frame_contacts(*args, **kwargs)

    def compute_replicate(self, ctx: ReplicateContext, replicate: int) -> Any:
        """Compute oligomer bridging metrics through the runner seam."""

        cache_file = self._replicate_cache_path(ctx.output_dir, ctx.settings, ctx.equilibration)

        cached = self._check_cache(
            PolymerBridgingReplicateResult,
            cache_file,
            recompute=ctx.recompute,
            sim_config=ctx.sim_config,
            settings=None,
        )
        if (
            cached is not None
            and self._cache_matches_context(
                cached,
                settings=ctx.settings,
                equilibration=ctx.equilibration,
                sim_config=ctx.condition.sim_config,
                source=cache_file,
            )
            and self._cache_matches_replicate_id(cached, replicate, source=cache_file)
        ):
            if ctx.result_path is not None and not ctx.result_path.exists():
                self._attach_settings_proof(cached, ctx.settings)
                self.save_result(cached, ctx.result_path)
            return cached

        if ctx.result_path is not None:
            cached = self._check_cache(
                PolymerBridgingReplicateResult,
                ctx.result_path,
                recompute=ctx.recompute,
                sim_config=ctx.sim_config,
                settings=None,
            )
            if (
                cached is not None
                and self._cache_matches_context(
                    cached,
                    settings=ctx.settings,
                    equilibration=ctx.equilibration,
                    sim_config=ctx.condition.sim_config,
                    source=ctx.result_path,
                )
                and self._cache_matches_replicate_id(cached, replicate, source=ctx.result_path)
            ):
                if not cache_file.exists():
                    self._attach_settings_proof(cached, ctx.settings)
                    self.save_result(cached, cache_file)
                return cached

        result = super().compute_replicate(ctx, replicate)
        self.save_result(result, cache_file)
        if ctx.result_path is not None:
            self.save_result(result, ctx.result_path)
        return result

    def build_runner(
        self,
        ctx: ReplicateContext,
        replicate: int,
        universe: Any,
        window: Any,
    ) -> Any:
        """Build the runner-backed polymer bridging execution object."""

        from polyzymd.analyses.polymer_bridging._runner import PolymerBridgingReplicateRunner

        del replicate, window
        return PolymerBridgingReplicateRunner(
            universe=universe,
            protein_selection=ctx.settings.protein_selection,
            polymer_selection=ctx.settings.polymer_selection,
            cutoff=float(ctx.settings.cutoff),
            min_ca_distance_angstrom=float(ctx.settings.min_ca_distance_angstrom),
        )

    def summarize_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
        runner: Any,
        window: Any,
    ) -> Any:
        """Convert runner observations into the replicate result schema."""

        from polyzymd.analyses.shared.config_hash import compute_config_hash
        from polyzymd.analyses.shared.loader import parse_time_string

        runner_result = runner.results
        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        stats = _compute_bridging_statistics_from_frames(
            runner_result.observations,
            min_ca_distance_angstrom=ctx.settings.min_ca_distance_angstrom,
        )
        return PolymerBridgingReplicateResult(
            settings_fingerprint=self._make_settings_cache_tag(ctx.settings),
            protein_selection=ctx.settings.protein_selection,
            polymer_selection=ctx.settings.polymer_selection,
            cutoff=float(ctx.settings.cutoff),
            config_hash=compute_config_hash(ctx.sim_config),
            equilibration_time=eq_value,
            equilibration_unit=eq_unit,
            replicate=replicate,
            n_frames=int(getattr(runner_result, "n_frames", window.n_frames_selected)),
            timestep_ps=float(
                getattr(runner_result, "timestep_ps", window.timestep_ps * window.step)
            ),
            min_ca_distance_angstrom=ctx.settings.min_ca_distance_angstrom,
            contacting_observations=stats["contacting_observations"],
            multisite_observations=stats["multisite_observations"],
            high_valency_observations=stats["high_valency_observations"],
            mean_contacts_per_contacting_oligomer=stats["mean_contacts_per_contacting_oligomer"],
            multisite_fraction=stats["multisite_fraction"],
            high_valency_fraction=stats["high_valency_fraction"],
            valency_probabilities=stats["valency_probabilities"],
            anchor_protein_group_probabilities=stats["anchor_protein_group_probabilities"],
            peripheral_protein_group_probabilities=stats["peripheral_protein_group_probabilities"],
            multivalent_protein_group_probabilities=stats[
                "multivalent_protein_group_probabilities"
            ],
            polymer_contact_type_probabilities=stats["polymer_contact_type_probabilities"],
            polymer_anchor_type_probabilities=stats["polymer_anchor_type_probabilities"],
            anchor_to_peripheral_group_matrix=stats["anchor_to_peripheral_group_matrix"],
            polymer_anchor_to_protein_anchor_matrix=stats[
                "polymer_anchor_to_protein_anchor_matrix"
            ],
            fragment_signature_probabilities=stats["fragment_signature_probabilities"],
        )

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> Any:
        """Aggregate bridging metrics across replicates."""
        recompute = getattr(ctx, "recompute", False)
        cache_file = self._aggregate_cache_path(
            ctx.output_dir,
            ctx.settings,
            ctx.equilibration,
            ctx.replicates,
        )

        cached = self._check_cache(
            PolymerBridgingAggregatedResult,
            cache_file,
            recompute=recompute,
            sim_config=ctx.condition.sim_config,
            settings=None,
        )
        if cached is not None and self._cache_matches_context(
            cached,
            settings=ctx.settings,
            equilibration=ctx.equilibration,
            sim_config=ctx.condition.sim_config,
            replicates=ctx.replicates,
            source=cache_file,
        ):
            if ctx.result_path is not None and not ctx.result_path.exists():
                self._attach_settings_proof(cached, ctx.settings)
                self.save_result(cached, ctx.result_path)
            return cached

        if ctx.result_path is not None:
            cached = self._check_cache(
                PolymerBridgingAggregatedResult,
                ctx.result_path,
                recompute=recompute,
                sim_config=ctx.condition.sim_config,
                settings=None,
            )
            if cached is not None and self._cache_matches_context(
                cached,
                settings=ctx.settings,
                equilibration=ctx.equilibration,
                sim_config=ctx.condition.sim_config,
                replicates=ctx.replicates,
                source=ctx.result_path,
            ):
                if not cache_file.exists():
                    self._attach_settings_proof(cached, ctx.settings)
                    self.save_result(cached, cache_file)
                return cached

        aligned_results = self._align_replicate_results(results, ctx.replicates)

        mean_contacts = [float(r.mean_contacts_per_contacting_oligomer) for r in aligned_results]
        multisite = [float(r.multisite_fraction) for r in aligned_results]
        high_valency = [float(r.high_valency_fraction) for r in aligned_results]

        valency_keys = sorted({key for r in aligned_results for key in r.valency_probabilities})
        per_rep_valency = {
            key: [float(r.valency_probabilities.get(key, 0.0)) for r in aligned_results]
            for key in valency_keys
        }

        anchor_groups = _collect_probability_series(
            aligned_results,
            "anchor_protein_group_probabilities",
        )
        peripheral_groups = _collect_probability_series(
            aligned_results,
            "peripheral_protein_group_probabilities",
        )
        multivalent_groups = _collect_probability_series(
            aligned_results,
            "multivalent_protein_group_probabilities",
        )
        polymer_contact_types = _collect_probability_series(
            aligned_results,
            "polymer_contact_type_probabilities",
        )
        polymer_anchor_types = _collect_probability_series(
            aligned_results,
            "polymer_anchor_type_probabilities",
        )
        fragment_signatures = _collect_probability_series(
            aligned_results,
            "fragment_signature_probabilities",
        )
        anchor_peripheral_mean, anchor_peripheral_sem = _aggregate_nested_matrices(
            [r.anchor_to_peripheral_group_matrix for r in aligned_results]
        )
        polymer_anchor_protein_mean, polymer_anchor_protein_sem = _aggregate_nested_matrices(
            [r.polymer_anchor_to_protein_anchor_matrix for r in aligned_results]
        )

        from polyzymd.analyses.shared.config_hash import compute_config_hash
        from polyzymd.analyses.shared.loader import parse_time_string

        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        aggregated = PolymerBridgingAggregatedResult(
            settings_fingerprint=self._make_settings_cache_tag(ctx.settings),
            protein_selection=ctx.settings.protein_selection,
            polymer_selection=ctx.settings.polymer_selection,
            cutoff=float(ctx.settings.cutoff),
            config_hash=compute_config_hash(ctx.condition.sim_config),
            equilibration_time=eq_value,
            equilibration_unit=eq_unit,
            n_replicates=len(aligned_results),
            replicates=list(ctx.replicates),
            min_ca_distance_angstrom=float(ctx.settings.min_ca_distance_angstrom),
            mean_contacts_per_contacting_oligomer=float(np.mean(mean_contacts)),
            mean_contacts_sem=_sem(mean_contacts),
            multisite_fraction=float(np.mean(multisite)),
            multisite_fraction_sem=_sem(multisite),
            high_valency_fraction=float(np.mean(high_valency)),
            high_valency_fraction_sem=_sem(high_valency),
            mean_contacts_per_contacting_oligomer_replicates=mean_contacts,
            multisite_fraction_replicates=multisite,
            high_valency_fraction_replicates=high_valency,
            valency_probabilities_mean={
                key: float(np.mean(vals)) for key, vals in per_rep_valency.items()
            },
            valency_probabilities_sem={key: _sem(vals) for key, vals in per_rep_valency.items()},
            valency_probabilities_per_replicate=per_rep_valency,
            anchor_protein_group_probabilities_mean=_series_mean(anchor_groups),
            anchor_protein_group_probabilities_sem=_series_sem(anchor_groups),
            anchor_protein_group_probabilities_per_replicate=anchor_groups,
            peripheral_protein_group_probabilities_mean=_series_mean(peripheral_groups),
            peripheral_protein_group_probabilities_sem=_series_sem(peripheral_groups),
            peripheral_protein_group_probabilities_per_replicate=peripheral_groups,
            multivalent_protein_group_probabilities_mean=_series_mean(multivalent_groups),
            multivalent_protein_group_probabilities_sem=_series_sem(multivalent_groups),
            multivalent_protein_group_probabilities_per_replicate=multivalent_groups,
            polymer_contact_type_probabilities_mean=_series_mean(polymer_contact_types),
            polymer_contact_type_probabilities_sem=_series_sem(polymer_contact_types),
            polymer_contact_type_probabilities_per_replicate=polymer_contact_types,
            polymer_anchor_type_probabilities_mean=_series_mean(polymer_anchor_types),
            polymer_anchor_type_probabilities_sem=_series_sem(polymer_anchor_types),
            polymer_anchor_type_probabilities_per_replicate=polymer_anchor_types,
            anchor_to_peripheral_group_matrix_mean=anchor_peripheral_mean,
            anchor_to_peripheral_group_matrix_sem=anchor_peripheral_sem,
            polymer_anchor_to_protein_anchor_matrix_mean=polymer_anchor_protein_mean,
            polymer_anchor_to_protein_anchor_matrix_sem=polymer_anchor_protein_sem,
            fragment_signature_probabilities_mean=_series_mean(fragment_signatures),
            fragment_signature_probabilities_sem=_series_sem(fragment_signatures),
            fragment_signature_probabilities_per_replicate=fragment_signatures,
        )
        self.save_result(aggregated, cache_file)
        if ctx.result_path is not None:
            self.save_result(aggregated, ctx.result_path)
        return aggregated

    def filter_conditions(
        self,
        conditions: list[Condition],
        settings: BaseModel | None = None,
    ) -> list[Condition]:
        """Exclude conditions without polymer atoms."""
        polymer_selection = (
            settings.polymer_selection
            if settings is not None
            else self.Settings().polymer_selection
        )
        return [
            c for c in conditions if _condition_has_polymer(c, polymer_selection=polymer_selection)
        ]

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Expose scalar metrics for the default comparison pipeline."""
        if isinstance(summary, dict):
            data = summary
        else:
            data = summary.model_dump()
        return {
            "mean_contacts_per_contacting_oligomer": MetricValue(
                name="mean_contacts_per_contacting_oligomer",
                mean=float(data["mean_contacts_per_contacting_oligomer"]),
                sem=float(data["mean_contacts_sem"]),
                replicate_values=list(data["mean_contacts_per_contacting_oligomer_replicates"]),
                higher_is_better=True,
                direction_labels=("less bridging", "unchanged", "more bridging"),
            ),
            "multisite_fraction": MetricValue(
                name="multisite_fraction",
                mean=float(data["multisite_fraction"]),
                sem=float(data["multisite_fraction_sem"]),
                replicate_values=list(data["multisite_fraction_replicates"]),
                higher_is_better=True,
                direction_labels=("less multisite", "unchanged", "more multisite"),
            ),
            "high_valency_fraction": MetricValue(
                name="high_valency_fraction",
                mean=float(data["high_valency_fraction"]),
                sem=float(data["high_valency_fraction_sem"]),
                replicate_values=list(data["high_valency_fraction_replicates"]),
                higher_is_better=True,
                direction_labels=("less high-valency", "unchanged", "more high-valency"),
            ),
        }

    def compare(self, ctx) -> ComparisonResult | None:
        """Run a NaN-safe comparison across conditions."""
        metrics_by_condition: dict[str, dict[str, MetricValue]] = {}
        for cond in ctx.conditions:
            summary = ctx.aggregated_results.get(cond.label)
            if summary is not None and not self._cache_matches_context(
                summary,
                settings=ctx.settings,
                equilibration=ctx.equilibration,
                sim_config=cond.sim_config,
                replicates=cond.replicates,
                allow_replicate_subset=True,
            ):
                summary = None
            if summary is None:
                agg_dir_parent = ctx.analysis_dirs.get(cond.label)
                if agg_dir_parent is None:
                    continue
                summary = self._load_validated_aggregated_result(
                    agg_dir_parent / "aggregated",
                    settings=ctx.settings,
                    equilibration=ctx.equilibration,
                    replicates=cond.replicates,
                    sim_config=cond.sim_config,
                    recompute=ctx.recompute,
                    allow_replicate_subset=True,
                )
            if summary is None:
                continue
            extracted = self.extract_metrics(summary)
            if extracted:
                metrics_by_condition[cond.label] = extracted

        if not metrics_by_condition:
            return None

        metric_names: list[str] = []
        seen: set[str] = set()
        for cond_metrics in metrics_by_condition.values():
            for metric_name in cond_metrics:
                if metric_name not in seen:
                    seen.add(metric_name)
                    metric_names.append(metric_name)

        all_pairwise: list[PairwiseResult] = []
        all_anova: list[ANOVAResult] = []
        all_rankings: dict[str, list[str]] = {}

        for metric_name in metric_names:
            per_cond = {
                label: metrics[metric_name]
                for label, metrics in metrics_by_condition.items()
                if metric_name in metrics
            }
            if not per_cond:
                continue
            if len(per_cond) >= 2:
                all_pairwise.extend(
                    _safe_pairwise_comparisons(
                        per_cond,
                        ctx.effective_control,
                        ttest_method=ctx.ttest_method,
                        posthoc_method=ctx.posthoc_method,
                        fdr_alpha=ctx.fdr_alpha,
                    )
                )
            maybe_anova = _safe_anova(per_cond, metric_name, alpha=ctx.fdr_alpha)
            if maybe_anova is not None:
                all_anova.append(maybe_anova)
            all_rankings[metric_name] = sorted(
                per_cond.keys(), key=lambda lb: per_cond[lb].mean, reverse=True
            )

        condition_summaries: list[ConditionSummary] = []
        for label, metrics in metrics_by_condition.items():
            extra: dict[str, Any] = {}
            for metric_name, mv in metrics.items():
                extra[f"{metric_name}_mean"] = mv.mean
                extra[f"{metric_name}_sem"] = mv.sem
                extra[f"{metric_name}_replicate_values"] = mv.replicate_values
            n_reps = len(next(iter(metrics.values())).replicate_values) if metrics else 0
            condition_summaries.append(ConditionSummary(label=label, n_replicates=n_reps, **extra))

        # Apply Benjamini-Hochberg FDR correction across all pairwise tests
        # only when using t-tests
        fdr_alpha = getattr(ctx, "fdr_alpha", 0.05)
        posthoc_method = getattr(ctx, "posthoc_method", "ttest_bh")
        if all_pairwise and posthoc_method == "ttest_bh":
            from polyzymd.analyses.shared.inferential_statistics import benjamini_hochberg

            raw_p = [pw.p_value for pw in all_pairwise]
            bh_results = benjamini_hochberg(raw_p, alpha=fdr_alpha)
            for pw, bh in zip(all_pairwise, bh_results, strict=False):
                pw.p_value_adjusted = bh.adjusted_p_value
                p_for_significance = (
                    bh.adjusted_p_value if bh.adjusted_p_value is not None else pw.p_value
                )
                pw.significant = p_for_significance <= fdr_alpha

        from datetime import datetime

        from polyzymd import __version__

        primary_ranking = all_rankings.get(metric_names[0], []) if metric_names else []
        return ComparisonResult(
            analysis_type=self.name,
            name=ctx.name,
            control_label=ctx.effective_control,
            fdr_alpha=fdr_alpha,
            ttest_method=ctx.ttest_method,
            posthoc_method=ctx.posthoc_method,
            conditions=condition_summaries,
            pairwise_comparisons=all_pairwise,
            anova=all_anova if all_anova else None,
            ranking=primary_ranking,
            rankings_by_metric=all_rankings if len(metric_names) > 1 else None,
            equilibration_time=ctx.equilibration,
            created_at=datetime.now().isoformat(),
            polyzymd_version=__version__,
        )

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format bridging comparison for CLI display."""
        from polyzymd.analyses.stats import format_scalar_comparison

        if output_format == "json":
            return result.model_dump_json(indent=2)

        pieces = []
        metric_specs = [
            (
                "multisite_fraction",
                "Polymer Bridging Comparison",
                "Multisite Fraction",
                "",
            ),
            (
                "mean_contacts_per_contacting_oligomer",
                "Average Oligomer Valency",
                "Contacts / Oligomer",
                "",
            ),
            (
                "high_valency_fraction",
                "High-Valency Oligomers",
                "High-Valency Fraction",
                "",
            ),
        ]
        for metric_key, title, label, unit in metric_specs:
            filtered = _filter_comparison_result(result, metric_key)
            if filtered is None:
                continue
            pieces.append(
                format_scalar_comparison(
                    filtered,
                    title=title,
                    metric_label=label,
                    metric_unit=unit,
                    metric_key=metric_key,
                    output_format=output_format,
                    higher_is_better=True,
                )
            )
        formatted = "\n\n".join(pieces)
        return prefix_experimental_output(formatted, EXPERIMENTAL_CHEMISTRY_FEATURES, output_format)

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate comparison plots for oligomer bridging."""
        import matplotlib.pyplot as plt

        plot_settings = getattr(ctx.plot_settings, self.name)
        comparison = self._load_comparison_result(ctx.comparison_path)
        if comparison is None:
            return []

        labels = [cond.label for cond in comparison.conditions]
        condition_by_label = {cond.label: cond for cond in ctx.conditions}
        colors = get_colors(len(labels), ctx.plot_settings)
        output_paths: list[Path] = []

        def _ensure_replicate_bases(bases: list[float], n_values: int) -> list[float]:
            """Return stack bases aligned to the current component replicates."""
            if len(bases) < n_values:
                bases.extend([0.0] * (n_values - len(bases)))
            return bases[:n_values]

        def _advance_replicate_bases(bases: list[float], values: Sequence[Any]) -> None:
            """Accumulate unsigned replicate stack bases from finite segment values."""
            _ensure_replicate_bases(bases, len(values))
            for idx, raw_value in enumerate(values):
                try:
                    value = float(raw_value)
                except (TypeError, ValueError):
                    continue
                if np.isfinite(value):
                    bases[idx] += value

        if plot_settings.generate_multisite_bars:
            means = [
                getattr(cond, "multisite_fraction_mean", 0.0) for cond in comparison.conditions
            ]
            sems = [getattr(cond, "multisite_fraction_sem", 0.0) for cond in comparison.conditions]
            replicate_values = [
                getattr(cond, "multisite_fraction_replicate_values", [])
                for cond in comparison.conditions
            ]
            fig, ax = plt.subplots(figsize=plot_settings.figsize_bars)
            x = np.arange(len(labels))
            ax.bar(x, means, yerr=sems, color=colors, alpha=ctx.plot_settings.theme.bar_alpha)
            scatter_replicate_values(
                ax,
                x,
                replicate_values,
                ctx.plot_settings,
                orientation="vertical",
                bar_width=0.8,
            )
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=35, ha="right")
            ax.set_ylim(0, 1)
            apply_axis_style(
                ax,
                ctx.plot_settings,
                title="Multisite Oligomer Fraction",
                ylabel="Fraction of contacting oligomers",
            )
            out = get_output_path(
                ctx.output_dir, "polymer_bridging_multisite_fraction", ctx.plot_settings
            )
            output_paths.append(
                save_figure(
                    fig,
                    out,
                    ctx.plot_settings,
                    experimental_features=EXPERIMENTAL_CHEMISTRY_FEATURES,
                )
            )
            plt.close(fig)

        if plot_settings.generate_mean_contacts_bars:
            means = [
                getattr(cond, "mean_contacts_per_contacting_oligomer_mean", 0.0)
                for cond in comparison.conditions
            ]
            sems = [
                getattr(cond, "mean_contacts_per_contacting_oligomer_sem", 0.0)
                for cond in comparison.conditions
            ]
            replicate_values = [
                getattr(cond, "mean_contacts_per_contacting_oligomer_replicate_values", [])
                for cond in comparison.conditions
            ]
            fig, ax = plt.subplots(figsize=plot_settings.figsize_bars)
            x = np.arange(len(labels))
            ax.bar(x, means, yerr=sems, color=colors, alpha=ctx.plot_settings.theme.bar_alpha)
            scatter_replicate_values(
                ax,
                x,
                replicate_values,
                ctx.plot_settings,
                orientation="vertical",
                bar_width=0.8,
            )
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=35, ha="right")
            apply_axis_style(
                ax,
                ctx.plot_settings,
                title="Distinct Protein Residues Per Contacting Oligomer",
                ylabel="Mean contacted residues",
            )
            out = get_output_path(
                ctx.output_dir, "polymer_bridging_mean_contacts", ctx.plot_settings
            )
            output_paths.append(
                save_figure(
                    fig,
                    out,
                    ctx.plot_settings,
                    experimental_features=EXPERIMENTAL_CHEMISTRY_FEATURES,
                )
            )
            plt.close(fig)

        if plot_settings.generate_valency_stack:
            fig, ax = plt.subplots(figsize=plot_settings.figsize_stack)
            x = np.arange(len(labels))
            keys = ["1", "2", "3+"]
            bottoms = np.zeros(len(labels), dtype=float)
            replicate_bottoms: list[list[float]] = [[] for _ in labels]
            for idx, key in enumerate(keys):
                vals = []
                replicate_values = []
                for cond in comparison.conditions:
                    condition = condition_by_label.get(cond.label)
                    analysis_dir = ctx.analysis_dirs.get(cond.label)
                    summary = None
                    if analysis_dir is not None and condition is not None:
                        summary = self._load_validated_aggregated_result(
                            analysis_dir / "aggregated",
                            settings=ctx.settings,
                            equilibration=ctx.equilibration,
                            replicates=condition.replicates,
                            sim_config=condition.sim_config,
                            recompute=False,
                            allow_replicate_subset=True,
                        )
                    if summary is None:
                        vals.append(0.0)
                        replicate_values.append([])
                        continue
                    if isinstance(summary, dict):
                        probs = summary.get("valency_probabilities_mean", {})
                        per_replicate = summary.get("valency_probabilities_per_replicate", {})
                    else:
                        probs = summary.valency_probabilities_mean
                        per_replicate = summary.valency_probabilities_per_replicate
                    vals.append(float(probs.get(key, 0.0)))
                    key_replicates = per_replicate.get(key)
                    replicate_values.append(
                        list(key_replicates) if key_replicates is not None else []
                    )
                ax.bar(x, vals, bottom=bottoms, label=key, color=colors[idx % len(colors)])
                for label_idx, reps in enumerate(replicate_values):
                    if reps:
                        replicate_bases = _ensure_replicate_bases(
                            replicate_bottoms[label_idx], len(reps)
                        )
                        scatter_stacked_segment_replicates(
                            ax,
                            float(x[label_idx]),
                            float(bottoms[label_idx]),
                            reps,
                            ctx.plot_settings,
                            replicate_base_values=replicate_bases,
                        )
                bottoms += np.array(vals, dtype=float)
                for label_idx, reps in enumerate(replicate_values):
                    _advance_replicate_bases(replicate_bottoms[label_idx], reps)
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=35, ha="right")
            ax.set_ylim(0, 1)
            apply_axis_style(
                ax,
                ctx.plot_settings,
                title="Oligomer Attachment Valency Distribution",
                ylabel="Probability",
            )
            ax.legend(title="Residues contacted")
            out = get_output_path(
                ctx.output_dir, "polymer_bridging_valency_distribution", ctx.plot_settings
            )
            output_paths.append(
                save_figure(
                    fig,
                    out,
                    ctx.plot_settings,
                    experimental_features=EXPERIMENTAL_CHEMISTRY_FEATURES,
                )
            )
            plt.close(fig)

        summaries = {}
        for cond in comparison.conditions:
            condition = condition_by_label.get(cond.label)
            analysis_dir = ctx.analysis_dirs.get(cond.label)
            if analysis_dir is None or condition is None:
                continue
            summary = self._load_validated_aggregated_result(
                analysis_dir / "aggregated",
                settings=ctx.settings,
                equilibration=ctx.equilibration,
                replicates=condition.replicates,
                sim_config=condition.sim_config,
                recompute=False,
                allow_replicate_subset=True,
            )
            if summary is not None:
                summaries[cond.label] = summary

        if plot_settings.generate_anchor_group_bars:
            fig, ax = plt.subplots(figsize=plot_settings.figsize_bars)
            keys = [
                k
                for k in PROTEIN_GROUP_ORDER
                if any(
                    _get_summary_probs(
                        summaries.get(label), "anchor_protein_group_probabilities_mean"
                    ).get(k, 0.0)
                    > 0.0
                    for label in labels
                )
            ]
            x = np.arange(len(labels))
            width = 0.8 / max(len(keys), 1)
            for idx, key in enumerate(keys):
                means = [
                    _get_summary_probs(
                        summaries.get(label), "anchor_protein_group_probabilities_mean"
                    ).get(key, 0.0)
                    for label in labels
                ]
                replicate_values = [
                    _get_summary_probs(
                        summaries.get(label), "anchor_protein_group_probabilities_per_replicate"
                    ).get(key, [])
                    for label in labels
                ]
                bar_positions = x + idx * width - 0.4 + width / 2
                ax.bar(bar_positions, means, width=width, label=key)
                scatter_replicate_values(
                    ax,
                    bar_positions,
                    replicate_values,
                    ctx.plot_settings,
                    orientation="vertical",
                    bar_width=width,
                )
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=35, ha="right")
            ax.set_ylim(0, 1)
            apply_axis_style(
                ax, ctx.plot_settings, title="Anchor Protein Residue Class", ylabel="Probability"
            )
            if keys:
                ax.legend(title="Anchor class")
            out = get_output_path(
                ctx.output_dir, "polymer_bridging_anchor_groups", ctx.plot_settings
            )
            output_paths.append(
                save_figure(
                    fig,
                    out,
                    ctx.plot_settings,
                    experimental_features=EXPERIMENTAL_CHEMISTRY_FEATURES,
                )
            )
            plt.close(fig)

        if plot_settings.generate_protein_group_stack:
            fig, ax = plt.subplots(figsize=plot_settings.figsize_stack)
            x = np.arange(len(labels))
            bottoms = np.zeros(len(labels), dtype=float)
            replicate_bottoms = [[] for _ in labels]
            for idx, key in enumerate(PROTEIN_GROUP_ORDER):
                vals = [
                    _get_summary_probs(
                        summaries.get(label), "multivalent_protein_group_probabilities_mean"
                    ).get(key, 0.0)
                    for label in labels
                ]
                replicate_values = [
                    _get_summary_probs(
                        summaries.get(label),
                        "multivalent_protein_group_probabilities_per_replicate",
                    ).get(key, [])
                    for label in labels
                ]
                if not any(vals):
                    continue
                ax.bar(x, vals, bottom=bottoms, label=key, color=colors[idx % len(colors)])
                for label_idx, reps in enumerate(replicate_values):
                    if reps:
                        replicate_bases = _ensure_replicate_bases(
                            replicate_bottoms[label_idx], len(reps)
                        )
                        scatter_stacked_segment_replicates(
                            ax,
                            float(x[label_idx]),
                            float(bottoms[label_idx]),
                            reps,
                            ctx.plot_settings,
                            replicate_base_values=replicate_bases,
                        )
                bottoms += np.array(vals, dtype=float)
                for label_idx, reps in enumerate(replicate_values):
                    _advance_replicate_bases(replicate_bottoms[label_idx], reps)
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=35, ha="right")
            ax.set_ylim(0, 1)
            apply_axis_style(
                ax,
                ctx.plot_settings,
                title="Protein Residue Classes In Multivalent Events",
                ylabel="Probability",
            )
            if np.any(bottoms > 0):
                ax.legend(title="Protein class")
            out = get_output_path(
                ctx.output_dir, "polymer_bridging_protein_group_distribution", ctx.plot_settings
            )
            output_paths.append(
                save_figure(
                    fig,
                    out,
                    ctx.plot_settings,
                    experimental_features=EXPERIMENTAL_CHEMISTRY_FEATURES,
                )
            )
            plt.close(fig)

        if plot_settings.generate_anchor_peripheral_heatmap and labels:
            reference = summaries.get(labels[0])
            matrix = _get_nested_summary(reference, "anchor_to_peripheral_group_matrix_mean")
            if matrix:
                fig, ax = plt.subplots(figsize=plot_settings.figsize_heatmap)
                rows = sorted(matrix)
                cols = sorted({col for row in matrix.values() for col in row})
                data = np.array(
                    [[matrix.get(row, {}).get(col, 0.0) for col in cols] for row in rows]
                )
                im = ax.imshow(data, aspect="auto", cmap="viridis")
                ax.set_xticks(np.arange(len(cols)))
                ax.set_xticklabels(cols, rotation=35, ha="right")
                ax.set_yticks(np.arange(len(rows)))
                ax.set_yticklabels(rows)
                apply_axis_style(
                    ax,
                    ctx.plot_settings,
                    title=f"Anchor vs Peripheral Classes ({labels[0]})",
                    ylabel="Anchor class",
                )
                fig.colorbar(im, ax=ax, label="Row-normalized probability")
                out = get_output_path(
                    ctx.output_dir, "polymer_bridging_anchor_peripheral_heatmap", ctx.plot_settings
                )
                output_paths.append(
                    save_figure(
                        fig,
                        out,
                        ctx.plot_settings,
                        experimental_features=EXPERIMENTAL_CHEMISTRY_FEATURES,
                    )
                )
                plt.close(fig)

        if plot_settings.generate_polymer_anchor_heatmap and labels:
            reference = summaries.get(labels[0])
            matrix = _get_nested_summary(reference, "polymer_anchor_to_protein_anchor_matrix_mean")
            if matrix:
                fig, ax = plt.subplots(figsize=plot_settings.figsize_heatmap)
                rows = sorted(matrix)
                cols = sorted({col for row in matrix.values() for col in row})
                data = np.array(
                    [[matrix.get(row, {}).get(col, 0.0) for col in cols] for row in rows]
                )
                im = ax.imshow(data, aspect="auto", cmap="magma")
                ax.set_xticks(np.arange(len(cols)))
                ax.set_xticklabels(cols, rotation=35, ha="right")
                ax.set_yticks(np.arange(len(rows)))
                ax.set_yticklabels(rows)
                apply_axis_style(
                    ax,
                    ctx.plot_settings,
                    title=f"Polymer Anchor vs Protein Anchor ({labels[0]})",
                    ylabel="Polymer anchor type",
                )
                fig.colorbar(im, ax=ax, label="Row-normalized probability")
                out = get_output_path(
                    ctx.output_dir, "polymer_bridging_polymer_anchor_heatmap", ctx.plot_settings
                )
                output_paths.append(
                    save_figure(
                        fig,
                        out,
                        ctx.plot_settings,
                        experimental_features=EXPERIMENTAL_CHEMISTRY_FEATURES,
                    )
                )
                plt.close(fig)

        if plot_settings.generate_fragment_signature_bars and labels:
            reference = summaries.get(labels[0])
            signature_probs = _get_summary_probs(reference, "fragment_signature_probabilities_mean")
            if signature_probs:
                fig, ax = plt.subplots(figsize=plot_settings.figsize_bars)
                sig_labels = list(signature_probs.keys())[:10]
                values = [signature_probs[label] for label in sig_labels]
                x_sig = np.arange(len(sig_labels))
                ax.bar(x_sig, values, color=colors[: len(sig_labels)])
                signature_replicates = _get_summary_probs(
                    reference, "fragment_signature_probabilities_per_replicate"
                )
                scatter_replicate_values(
                    ax,
                    x_sig,
                    [signature_replicates.get(label, []) for label in sig_labels],
                    ctx.plot_settings,
                    orientation="vertical",
                    bar_width=0.8,
                )
                ax.set_xticks(x_sig)
                ax.set_xticklabels(sig_labels, rotation=35, ha="right")
                ax.set_ylim(0, max(values) * 1.15 if values else 1.0)
                apply_axis_style(
                    ax,
                    ctx.plot_settings,
                    title=f"Top Fragment Signatures ({labels[0]})",
                    ylabel="Probability",
                )
                out = get_output_path(
                    ctx.output_dir, "polymer_bridging_fragment_signatures", ctx.plot_settings
                )
                output_paths.append(
                    save_figure(
                        fig,
                        out,
                        ctx.plot_settings,
                        experimental_features=EXPERIMENTAL_CHEMISTRY_FEATURES,
                    )
                )
                plt.close(fig)

        return output_paths

    def _load_comparison_result(self, path: Path | None) -> ComparisonResult | None:
        if path is None or not path.exists():
            return None
        return ComparisonResult.load(path)


def _condition_has_polymer(cond: Condition, polymer_selection: str = "chainID C") -> bool:
    """Check whether a condition likely contains polymer atoms."""
    sim_config = cond.sim_config

    if hasattr(sim_config, "topology"):
        topo = sim_config.topology
        if hasattr(topo, "chains") and topo.chains:
            requested_chain = _simple_chain_id_selection(polymer_selection)
            chain_ids = [c.chain_id if hasattr(c, "chain_id") else c for c in topo.chains]
            if requested_chain is not None and requested_chain in chain_ids:
                return True

    try:
        import MDAnalysis as mda

        from polyzymd.analyses.shared.loader import TrajectoryLoader

        loader = TrajectoryLoader(sim_config)
        for rep in cond.replicates:
            run_dir = sim_config.get_working_directory(rep)
            topo_path = loader.find_topology(run_dir)
            universe = mda.Universe(str(topo_path))
            if len(universe.select_atoms(polymer_selection)) > 0:
                return True
    except ImportError:
        return False
    except (FileNotFoundError, OSError):
        return False
    except (ValueError, TypeError):
        return False

    return False


def _simple_chain_id_selection(selection: str) -> str | None:
    """Extract a simple chain ID from an MDAnalysis selection string.

    Parameters
    ----------
    selection : str
        MDAnalysis selection string.

    Returns
    -------
    str or None
        Requested chain ID for strictly simple ``chainID X`` selections,
        otherwise ``None``.
    """

    tokens = selection.strip().split()
    if len(tokens) == 2 and tokens[0].lower() == "chainid":
        return tokens[1]
    return None


def _compute_bridging_statistics_from_frames(
    frame_contacts: Sequence[set[int] | FrameContactObservation],
    *,
    min_ca_distance_angstrom: float,
    ca_distances: ResiduePairDistances | None = None,
) -> dict[str, Any]:
    """Compute per-replicate oligomer bridging statistics from frame contacts.

    Processes a list of observations and produces all primary metrics and
    chemistry-aware probability distributions.

    Parameters
    ----------
    frame_contacts : Sequence[set[int] | FrameContactObservation]
        List of observations. Each is either a
        ``PolymerBridgingObservation`` (with full chemistry metadata), a
        ``(residues, ca_distances)`` tuple, or a bare ``set[int]`` of
        contacted residue IDs.
    min_ca_distance_angstrom : float
        Minimum frame-wise C-alpha distance for eligibility.
    ca_distances : ResiduePairDistances or None, optional
        Fallback CA distances used when observations are bare sets.

    Returns
    -------
    dict[str, Any]
        Dictionary containing all computed metrics and probability
        distributions, keyed to match ``PolymerBridgingReplicateResult``
        field names.
    """
    eligible_counts: list[int] = []
    default_ca_distances = ca_distances or {}
    anchor_group_counter: Counter[str] = Counter()
    peripheral_group_counter: Counter[str] = Counter()
    multivalent_group_counter: Counter[str] = Counter()
    polymer_contact_counter: Counter[str] = Counter()
    polymer_anchor_counter: Counter[str] = Counter()
    anchor_to_peripheral_counter: Counter[tuple[str, str]] = Counter()
    polymer_anchor_to_protein_anchor_counter: Counter[tuple[str, str]] = Counter()
    fragment_signature_counter: Counter[str] = Counter()
    n_multivalent_events = 0

    for observation in frame_contacts:
        if isinstance(observation, PolymerBridgingObservation):
            residues = observation.protein_residues
            observation_distances = observation.ca_distances
        elif isinstance(observation, tuple):
            residues, observation_distances = observation
        else:
            residues = observation
            observation_distances = default_ca_distances
        if not residues:
            continue
        eligible = _count_eligible_residues(
            residues, min_ca_distance_angstrom, observation_distances
        )
        eligible_counts.append(eligible)

        if not isinstance(observation, PolymerBridgingObservation):
            continue

        eligible_residues = _eligible_residue_set(
            observation.protein_residues,
            min_ca_distance_angstrom,
            observation.ca_distances,
        )
        if len(eligible_residues) <= 1:
            continue

        n_multivalent_events += 1
        fragment_signature_counter[_format_signature(observation.fragment_signature)] += 1
        anchor = _find_anchor_pair(observation, eligible_residues)
        anchor_protein_group = (
            observation.protein_groups.get(anchor[1], "unknown") if anchor else "unknown"
        )
        if anchor:
            polymer_anchor_counter[observation.polymer_resnames.get(anchor[0], "unknown")] += 1
            polymer_anchor_to_protein_anchor_counter[
                (observation.polymer_resnames.get(anchor[0], "unknown"), anchor_protein_group)
            ] += 1
        anchor_group_counter[anchor_protein_group] += 1

        for polymer_resid in observation.contacting_polymer_resids:
            polymer_contact_counter[observation.polymer_resnames.get(polymer_resid, "unknown")] += 1

        for resid in eligible_residues:
            group = observation.protein_groups.get(resid, "unknown")
            multivalent_group_counter[group] += 1
            if anchor and resid == anchor[1]:
                continue
            peripheral_group_counter[group] += 1
            anchor_to_peripheral_counter[(anchor_protein_group, group)] += 1

    if not eligible_counts:
        probs = {"1": 0.0, "2": 0.0, "3+": 0.0}
        return {
            "contacting_observations": 0,
            "multisite_observations": 0,
            "high_valency_observations": 0,
            "mean_contacts_per_contacting_oligomer": 0.0,
            "multisite_fraction": 0.0,
            "high_valency_fraction": 0.0,
            "valency_probabilities": probs,
            "anchor_protein_group_probabilities": {},
            "peripheral_protein_group_probabilities": {},
            "multivalent_protein_group_probabilities": {},
            "polymer_contact_type_probabilities": {},
            "polymer_anchor_type_probabilities": {},
            "anchor_to_peripheral_group_matrix": {},
            "polymer_anchor_to_protein_anchor_matrix": {},
            "fragment_signature_probabilities": {},
        }

    arr = np.array(eligible_counts, dtype=float)
    probs = {
        "1": float(np.mean(arr == 1)),
        "2": float(np.mean(arr == 2)),
        "3+": float(np.mean(arr >= 3)),
    }
    return {
        "contacting_observations": len(eligible_counts),
        "multisite_observations": int(np.sum(arr > 1)),
        "high_valency_observations": int(np.sum(arr >= 3)),
        "mean_contacts_per_contacting_oligomer": float(np.mean(arr)),
        "multisite_fraction": float(np.mean(arr > 1)),
        "high_valency_fraction": float(np.mean(arr >= 3)),
        "valency_probabilities": probs,
        "anchor_protein_group_probabilities": _normalize_counter(anchor_group_counter),
        "peripheral_protein_group_probabilities": _normalize_counter(peripheral_group_counter),
        "multivalent_protein_group_probabilities": _normalize_counter(multivalent_group_counter),
        "polymer_contact_type_probabilities": _normalize_counter(polymer_contact_counter),
        "polymer_anchor_type_probabilities": _normalize_counter(polymer_anchor_counter),
        "anchor_to_peripheral_group_matrix": _normalize_nested_counter(
            anchor_to_peripheral_counter
        ),
        "polymer_anchor_to_protein_anchor_matrix": _normalize_nested_counter(
            polymer_anchor_to_protein_anchor_counter
        ),
        "fragment_signature_probabilities": _normalize_counter(
            fragment_signature_counter,
            top_n=10,
        ),
    }


def _count_eligible_residues(
    residues: set[int],
    min_ca_distance_angstrom: float,
    ca_distances: ResiduePairDistances,
) -> int:
    """Count residues that belong to at least one geometrically eligible pair.

    When ``min_ca_distance_angstrom <= 0``, returns ``len(residues)``
    (no filtering). Otherwise, returns the count of residues that
    participate in at least one pair with frame-wise CA distance >=
    the threshold.

    Parameters
    ----------
    residues : set[int]
        Contacted protein residue IDs.
    min_ca_distance_angstrom : float
        Minimum CA-CA distance for a pair to be eligible.
    ca_distances : ResiduePairDistances
        Frame-wise CA distances keyed by ``(min_resid, max_resid)``.

    Returns
    -------
    int
        Number of eligible residues (minimum 1 if any residues contacted).
    """
    if len(residues) <= 1 or min_ca_distance_angstrom <= 0:
        return len(residues)

    residues_sorted = sorted(residues)
    eligible: set[int] = set()
    for i, resid_i in enumerate(residues_sorted):
        for resid_j in residues_sorted[i + 1 :]:
            key = (min(resid_i, resid_j), max(resid_i, resid_j))
            if ca_distances.get(key, 0.0) >= min_ca_distance_angstrom:
                eligible.add(resid_i)
                eligible.add(resid_j)
    return len(eligible) if eligible else 1


def _eligible_residue_set(
    residues: set[int],
    min_ca_distance_angstrom: float,
    ca_distances: ResiduePairDistances,
) -> set[int]:
    """Return residues that belong to at least one eligible pair."""
    if len(residues) <= 1 or min_ca_distance_angstrom <= 0:
        return set(residues)

    residues_sorted = sorted(residues)
    eligible: set[int] = set()
    for i, resid_i in enumerate(residues_sorted):
        for resid_j in residues_sorted[i + 1 :]:
            key = (min(resid_i, resid_j), max(resid_i, resid_j))
            if ca_distances.get(key, 0.0) >= min_ca_distance_angstrom:
                eligible.add(resid_i)
                eligible.add(resid_j)
    return eligible or {min(residues_sorted)}


def _find_anchor_pair(
    observation: PolymerBridgingObservation,
    eligible_residues: set[int],
) -> tuple[int, int] | None:
    """Choose the closest polymer-protein residue pair within the eligible set.

    The anchor pair is the (polymer_resid, protein_resid) with the minimum
    atom-level distance, restricted to protein residues in the eligible set.

    Parameters
    ----------
    observation : PolymerBridgingObservation
        The full observation with pair minimum distances.
    eligible_residues : set[int]
        Protein residue IDs that passed the CA distance filter.

    Returns
    -------
    tuple[int, int] or None
        ``(polymer_resid, protein_resid)`` of the anchor pair, or None if
        no eligible pairs exist.
    """
    best_pair: tuple[int, int] | None = None
    best_distance = float("inf")
    for pair, distance in observation.pair_min_distances.items():
        polymer_resid, protein_resid = pair
        if protein_resid not in eligible_residues:
            continue
        if distance < best_distance:
            best_distance = distance
            best_pair = pair
    return best_pair


def _normalize_counter(counter: Counter[str], top_n: int | None = None) -> dict[str, float]:
    """Convert counts to probabilities."""
    if not counter:
        return {}
    items = counter.most_common(top_n) if top_n is not None else sorted(counter.items())
    total = sum(count for _, count in items)
    if total <= 0:
        return {}
    return {key: float(count / total) for key, count in items}


def _normalize_nested_counter(counter: Counter[tuple[str, str]]) -> dict[str, dict[str, float]]:
    """Convert pair counts into row-normalized probabilities."""
    rows: dict[str, Counter[str]] = {}
    for (left, right), count in counter.items():
        rows.setdefault(left, Counter())[right] += count
    result: dict[str, dict[str, float]] = {}
    for left, row_counter in rows.items():
        total = sum(row_counter.values())
        result[left] = {
            right: float(count / total) for right, count in sorted(row_counter.items()) if total > 0
        }
    return result


def _format_signature(signature: Sequence[str]) -> str:
    return "-".join(signature)


def _collect_probability_series(results: Sequence[Any], attr: str) -> dict[str, list[float]]:
    keys = sorted({key for r in results for key in getattr(r, attr, {}).keys()})
    return {key: [float(getattr(r, attr, {}).get(key, 0.0)) for r in results] for key in keys}


def _series_mean(series: dict[str, list[float]]) -> dict[str, float]:
    return {key: float(np.mean(vals)) for key, vals in series.items()}


def _series_sem(series: dict[str, list[float]]) -> dict[str, float]:
    return {key: _sem(vals) for key, vals in series.items()}


def _aggregate_nested_matrices(
    matrices: Sequence[dict[str, dict[str, float]]],
) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, float]]]:
    rows = sorted({row for matrix in matrices for row in matrix})
    cols = sorted({col for matrix in matrices for row in matrix.values() for col in row})
    mean_result: dict[str, dict[str, float]] = {}
    sem_result: dict[str, dict[str, float]] = {}
    for row in rows:
        mean_result[row] = {}
        sem_result[row] = {}
        for col in cols:
            values = [float(matrix.get(row, {}).get(col, 0.0)) for matrix in matrices]
            mean_result[row][col] = float(np.mean(values))
            sem_result[row][col] = _sem(values)
    return mean_result, sem_result


def _get_summary_probs(summary: Any, field_name: str) -> dict[str, float]:
    if summary is None:
        return {}
    if isinstance(summary, dict):
        return dict(summary.get(field_name, {}))
    return dict(getattr(summary, field_name, {}))


def _get_nested_summary(summary: Any, field_name: str) -> dict[str, dict[str, float]]:
    if summary is None:
        return {}
    if isinstance(summary, dict):
        return dict(summary.get(field_name, {}))
    return dict(getattr(summary, field_name, {}))


def _sem(values: Sequence[float]) -> float:
    arr = np.asarray(list(values), dtype=float)
    if arr.size <= 1:
        return 0.0
    return float(np.std(arr, ddof=1) / np.sqrt(arr.size))


def _filter_comparison_result(result: ComparisonResult, metric_key: str) -> ComparisonResult | None:
    conditions: list[ConditionSummary] = []
    for cond in result.conditions:
        payload = cond.model_dump()
        mean_key = f"{metric_key}_mean"
        sem_key = f"{metric_key}_sem"
        repl_key = f"{metric_key}_replicate_values"
        if mean_key not in payload:
            continue
        cond_summary = ConditionSummary(
            label=payload["label"],
            n_replicates=payload.get("n_replicates", 0),
            **{
                mean_key: payload.get(mean_key, 0.0),
                sem_key: payload.get(sem_key, 0.0),
                repl_key: payload.get(repl_key, []),
            },
        )
        conditions.append(cond_summary)

    if not conditions:
        return None

    pairwise = [p for p in result.pairwise_comparisons if p.metric == metric_key]
    anova = [a for a in (result.anova or []) if a.metric == metric_key]
    ranking = (result.rankings_by_metric or {}).get(metric_key, []) or result.ranking
    return ComparisonResult(
        analysis_type=result.analysis_type,
        name=result.name,
        control_label=result.control_label,
        conditions=conditions,
        pairwise_comparisons=pairwise,
        anova=anova or None,
        ranking=ranking,
        rankings_by_metric=None,
        equilibration_time=result.equilibration_time,
        created_at=result.created_at,
        polyzymd_version=result.polyzymd_version,
    )


def _safe_pairwise_comparisons(
    metrics_by_condition: dict[str, MetricValue],
    control_label: str | None,
    *,
    ttest_method: str = "student",
    posthoc_method: str = "ttest_bh",
    fdr_alpha: float = 0.05,
) -> list[PairwiseResult]:
    results = pairwise_comparisons(
        metrics_by_condition,
        control_label,
        ttest_method=ttest_method,
        posthoc_method=posthoc_method,
        fdr_alpha=fdr_alpha,
    )

    for result in results:
        has_nan = (
            not np.isfinite(result.t_statistic)
            or not np.isfinite(result.p_value)
            or not np.isfinite(result.cohens_d)
        )
        if has_nan:
            logger.debug(
                "Non-finite statistics for %s vs %s (%s); marking as not testable "
                "(likely zero-variance replicates)",
                result.condition_a,
                result.condition_b,
                result.metric,
            )
            result.t_statistic = 0.0
            result.p_value = 1.0
            result.cohens_d = 0.0
            result.significant = False
            result.direction = "not_testable"
            result.effect_size_interpretation = "not_testable"
        if result.p_value_adjusted is not None and not np.isfinite(result.p_value_adjusted):
            result.p_value_adjusted = result.p_value

    return results


def _safe_anova(
    metrics_by_condition: dict[str, MetricValue],
    metric_name: str,
    *,
    alpha: float,
) -> ANOVAResult | None:
    result = anova_test(metrics_by_condition, metric_name=metric_name, alpha=alpha)
    if result is None:
        return None

    has_nan = not np.isfinite(result.f_statistic) or not np.isfinite(result.p_value)
    if has_nan:
        logger.debug(
            "Non-finite ANOVA statistics for metric '%s'; marking as not testable "
            "(likely zero-variance replicates)",
            metric_name,
        )
        result.f_statistic = 0.0
        result.p_value = 1.0
        result.significant = False
    return result


def _compute_frame_contacts(
    condition: Condition,
    replicate: int,
    *,
    protein_selection: str,
    polymer_selection: str,
    cutoff: float,
    equilibration: str,
    min_ca_distance_angstrom: float = 0.0,
) -> tuple[list[FrameContactObservation], int, float]:
    """Compute observations through the runner compatibility facade.

    This helper is retained for tests and external callers that used the
    pre-runner API. The trajectory-native implementation lives in
    :mod:`polyzymd.analyses.polymer_bridging._runner`.
    """

    from polyzymd.analyses.polymer_bridging._runner import compute_frame_contacts

    return compute_frame_contacts(
        condition,
        replicate,
        protein_selection=protein_selection,
        polymer_selection=polymer_selection,
        cutoff=cutoff,
        equilibration=equilibration,
        min_ca_distance_angstrom=min_ca_distance_angstrom,
    )
