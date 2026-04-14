"""Shared helpers for computing binding preference from contacts data.

Used by the contacts, binding_free_energy, and polymer_affinity analysis
plugins to compute SASA-based binding preference enrichment.

Public functions
----------------
- find_enzyme_pdb()
- resolve_enzyme_pdb()
- try_load_cached_binding_preference()
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Protocol, Sequence

if TYPE_CHECKING:
    from polyzymd.analyses.shared.binding_preference import (
        AggregatedBindingPreferenceResult,
        BindingPreferenceResult,
    )

logger = logging.getLogger("polyzymd.analyses")


class ConditionLike(Protocol):
    """Minimal condition protocol required by BP helper functions.

    Attributes
    ----------
    label : str
        Condition label
    replicates : Sequence[int]
        Replicate IDs associated with this condition
    """

    label: str
    replicates: Sequence[int]


def find_enzyme_pdb(sim_config: Any) -> Path | None:
    """Find enzyme PDB file from simulation config.

    Searches common locations relative to the project directory.

    Parameters
    ----------
    sim_config : SimulationConfig
        Simulation configuration (must have ``output.projects_directory``).

    Returns
    -------
    Path or None
        Path to enzyme PDB, or None if not found.

    Raises
    ------
    ValueError
        If a glob pattern matches multiple candidate enzyme PDB files.
    """
    import glob as glob_module

    project_dir = sim_config.output.projects_directory
    possible_paths = [
        project_dir / "structures" / "enzyme.pdb",
        project_dir / "input" / "enzyme.pdb",
        project_dir.parent / "structures" / "enzyme.pdb",
        project_dir.parent / "enzyme.pdb",
    ]

    for path in possible_paths:
        if path.exists():
            return path

    # Try glob for any PDB with "enzyme" in name
    patterns = [
        str(project_dir / "**" / "*enzyme*.pdb"),
        str(project_dir.parent / "*enzyme*.pdb"),
    ]
    for pattern in patterns:
        matches = sorted(glob_module.glob(pattern, recursive=True))
        if matches:
            if len(matches) > 1:
                raise ValueError(
                    f"Ambiguous enzyme PDB auto-discovery for pattern '{pattern}': "
                    f"{len(matches)} matches found: " + ", ".join(matches)
                )
            return Path(matches[0])

    return None


def resolve_enzyme_pdb(
    enzyme_pdb_setting: str | None,
    source_path: Path | None,
    sim_config: Any,
) -> Path | None:
    """Resolve the enzyme PDB path from settings or auto-discovery.

    Parameters
    ----------
    enzyme_pdb_setting : str or None
        Explicit enzyme PDB path from analysis settings (e.g.,
        ``enzyme_pdb_for_sasa``). If relative, resolved against
        *source_path*'s parent directory.
    source_path : Path or None
        Path to the comparison.yaml file (used to resolve relative paths).
    sim_config : Any
        Simulation configuration for auto-discovery fallback.

    Returns
    -------
    Path or None
        Resolved enzyme PDB path, or None if not found.
    """
    if enzyme_pdb_setting:
        if source_path:
            enzyme_pdb = source_path.parent / enzyme_pdb_setting
        else:
            enzyme_pdb = Path(enzyme_pdb_setting)
        if enzyme_pdb.exists():
            return enzyme_pdb
        logger.warning(f"Explicit enzyme_pdb_for_sasa not found at {enzyme_pdb}")
        return None

    return find_enzyme_pdb(sim_config)


def try_load_cached_binding_preference(
    cond: ConditionLike,
    analysis_dir: Path,
    *,
    settings_fp: str | None = None,
) -> "AggregatedBindingPreferenceResult | BindingPreferenceResult | None":
    """Try to load cached binding preference results for a condition.

    Searches for binding preference files in order of preference:
    1. binding_preference_aggregated.json
    2. binding_preference_aggregated_reps*.json (glob pattern)
    3. binding_preference.json (single replicate)
    4. Per-replicate files (binding_preference_rep{N}.json)

    Parameters
    ----------
    cond : ConditionLike
        Condition to load.
    analysis_dir : Path
        Analysis directory for this condition.
    settings_fp : str or None, optional
        Settings fingerprint for cache lookup. When provided, fingerprinted
        cache files are searched first, then legacy filenames.

    Returns
    -------
    AggregatedBindingPreferenceResult | BindingPreferenceResult | None
        Loaded result, or None if not found.
    """
    import glob as glob_module

    from polyzymd.analyses.shared.binding_preference import (
        AggregatedBindingPreferenceResult,
        BindingPreferenceResult,
        aggregate_binding_preference,
    )

    if settings_fp is not None:
        fp_agg_path = analysis_dir / f"binding_preference_aggregated_s{settings_fp}.json"
        if fp_agg_path.exists():
            result = AggregatedBindingPreferenceResult.load(fp_agg_path)
            logger.debug(f"Loaded aggregated binding preference for {cond.label}")
            return result

        fp_agg_pattern = str(
            analysis_dir / f"binding_preference_aggregated_s{settings_fp}_reps*.json"
        )
        fp_agg_matches = sorted(glob_module.glob(fp_agg_pattern))
        if len(fp_agg_matches) == 1:
            result = AggregatedBindingPreferenceResult.load(fp_agg_matches[0])
            logger.debug(f"Loaded aggregated binding preference for {cond.label}")
            return result
        if len(fp_agg_matches) > 1:
            raise ValueError(
                f"Ambiguous binding preference cache for {cond.label}: "
                f"found {len(fp_agg_matches)} files matching '{fp_agg_pattern}': "
                + ", ".join(fp_agg_matches)
            )

        fp_single_path = analysis_dir / f"binding_preference_s{settings_fp}.json"
        if fp_single_path.exists():
            result = BindingPreferenceResult.load(fp_single_path)
            logger.debug(f"Loaded single binding preference for {cond.label}")
            return result

        fp_rep_results = []
        for rep in cond.replicates:
            fp_rep_path = analysis_dir / f"binding_preference_s{settings_fp}_rep{rep}.json"
            if fp_rep_path.exists():
                fp_rep_results.append(BindingPreferenceResult.load(fp_rep_path))

        if fp_rep_results:
            agg_result = aggregate_binding_preference(fp_rep_results)
            logger.debug(
                f"Aggregated {len(fp_rep_results)} replicate binding preference results "
                f"for {cond.label}"
            )
            return agg_result

    # Try aggregated result first (multi-replicate)
    agg_path = analysis_dir / "binding_preference_aggregated.json"
    if agg_path.exists():
        result = AggregatedBindingPreferenceResult.load(agg_path)
        logger.debug(f"Loaded aggregated binding preference for {cond.label}")
        return result

    # Try aggregated result with rep range in name (e.g., _reps1-3.json)
    agg_pattern = str(analysis_dir / "binding_preference_aggregated_reps*.json")
    agg_matches = sorted(glob_module.glob(agg_pattern))
    if len(agg_matches) == 1:
        result = AggregatedBindingPreferenceResult.load(agg_matches[0])
        logger.debug(f"Loaded aggregated binding preference for {cond.label}")
        return result
    if len(agg_matches) > 1:
        raise ValueError(
            f"Ambiguous binding preference cache for {cond.label}: "
            f"found {len(agg_matches)} files: " + ", ".join(agg_matches)
        )

    # Try single replicate result
    single_path = analysis_dir / "binding_preference.json"
    if single_path.exists():
        result = BindingPreferenceResult.load(single_path)
        logger.debug(f"Loaded single binding preference for {cond.label}")
        return result

    # Try per-replicate results and aggregate them
    rep_results = []
    for rep in cond.replicates:
        rep_path = analysis_dir / f"binding_preference_rep{rep}.json"
        if rep_path.exists():
            rep_results.append(BindingPreferenceResult.load(rep_path))

    if rep_results:
        agg_result = aggregate_binding_preference(rep_results)
        logger.debug(
            f"Aggregated {len(rep_results)} replicate binding preference results for {cond.label}"
        )
        return agg_result

    return None
