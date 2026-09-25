"""How much of the planned data a result was computed from.

Incomplete data is never refused. A comparison can run while some replicates
have not started and others are still being written, so users and agents can
look at what they have. Every result instead records what it was computed from,
so a preliminary number cannot be mistaken for a final one:

- each replicate records the frames it used, the simulated time they cover,
  the planned production length, and any production segment that was running,
  interrupted, failed or left out;
- each condition records the replicates its comparison lists, the ones used,
  why any were left out, and whether all of them are complete;
- each comparison is complete only when every condition is.

A replicate is complete when its trajectory covers the planned production
length and every production segment finished and was read. A replicate whose
config states no planned length is judged by its segments alone.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence

#: A trajectory covering this fraction of the planned length counts as finished,
#: which absorbs the one-frame difference between frame count and time span.
COMPLETE_FRACTION = 0.99


def replicate_completeness(
    frames: Any, provenance: Mapping[str, Any], sim_config: Any
) -> dict[str, Any]:
    """Describe how much of one replicate's planned production was analysed.

    Parameters
    ----------
    frames : FrameSelection
        Production window the plugin read.
    provenance : mapping
        Universe provider provenance, with ``segment_status`` and
        ``excluded_segments``.
    sim_config : SimulationConfig
        Condition config, for the planned production length.

    Returns
    -------
    dict
        ``frames_used``, ``analysed_ps`` (window start and end), ``covered_ns``
        (simulated time in the trajectory), ``planned_ns``,
        ``production_fraction``, ``unfinished_segments`` (segment index to its
        status), ``excluded_segments`` and ``complete``.
    """
    timestep = frames.timestep_ps
    total = frames.n_frames_total
    first = frames.first_frame_time_ps or 0.0
    covered_ns = total * timestep / 1000.0 if timestep and total else None
    start_ps = frames.selected_start_time_ps
    end_ps = first + (total - 1) * timestep if timestep and total else None
    planned_ns = _planned_ns(sim_config)
    fraction = covered_ns / planned_ns if covered_ns is not None and planned_ns else None
    status = provenance.get("segment_status") or {}
    unfinished = {
        str(index): str(state) for index, state in sorted(status.items()) if state != "completed"
    }
    excluded = sorted(int(index) for index in provenance.get("excluded_segments") or [])
    reached = fraction is None or fraction >= COMPLETE_FRACTION
    return {
        "frames_used": frames.n_frames_selected,
        "analysed_ps": [start_ps, end_ps],
        "covered_ns": covered_ns,
        "planned_ns": planned_ns,
        "production_fraction": fraction,
        "unfinished_segments": unfinished,
        "excluded_segments": excluded,
        "complete": bool(reached and not unfinished and not excluded),
    }


def condition_completeness(
    listed: Sequence[int],
    used: Mapping[int, Mapping[str, Any] | None],
    dropped: Mapping[int, str],
) -> dict[str, Any]:
    """Describe which of a condition's listed replicates a result used.

    Parameters
    ----------
    listed : sequence of int
        Replicates the comparison lists for the condition.
    used : mapping
        Replicate number to its :func:`replicate_completeness` record, or
        ``None`` for a result that predates the record.
    dropped : mapping
        Replicate number to the reason it was left out.

    Returns
    -------
    dict
        ``replicates_listed``, ``replicates_used``, ``dropped`` (replicate and
        reason), ``replicates`` (the per-replicate records),
        ``min_production_fraction`` and ``complete``.
    """
    records = {int(replicate): dict(record or {}) for replicate, record in used.items()}
    fractions = [
        record["production_fraction"]
        for record in records.values()
        if record.get("production_fraction") is not None
    ]
    all_used = set(records) >= {int(replicate) for replicate in listed}
    return {
        "replicates_listed": sorted(int(replicate) for replicate in listed),
        "replicates_used": sorted(records),
        "dropped": [
            {"replicate": int(replicate), "reason": str(reason)}
            for replicate, reason in sorted(dropped.items())
        ],
        "replicates": {str(replicate): record for replicate, record in sorted(records.items())},
        "min_production_fraction": min(fractions) if fractions else None,
        "complete": bool(
            all_used
            and records
            and all(record.get("complete", False) for record in records.values())
        ),
    }


def comparison_completeness(
    conditions: Mapping[str, Mapping[str, Any] | None],
    dropped_conditions: Sequence[str] = (),
) -> dict[str, Any]:
    """Describe a comparison from its conditions' completeness records."""
    records = {label: dict(record or {}) for label, record in conditions.items()}
    return {
        "complete": bool(
            not dropped_conditions and all(record.get("complete") for record in records.values())
        ),
        "conditions": records,
        "dropped_conditions": list(dropped_conditions),
    }


def summary(label: str, record: Mapping[str, Any] | None) -> str | None:
    """One line saying how a condition's result is partial, or ``None`` if it is not.

    For example ``"SBMA: replicates 1-3 of 1-5; replicate 2 at 64% of 100 ns"``.
    """
    if not record or record.get("complete"):
        return None
    listed = record.get("replicates_listed") or []
    used = record.get("replicates_used") or []
    parts = []
    if not record.get("replicates_used"):
        parts.append("no completeness record")
    elif set(used) != set(listed):
        parts.append(f"replicates {_ranges(used)} of {_ranges(listed)}")
    for replicate, item in (record.get("replicates") or {}).items():
        if item.get("complete", True):
            continue
        fraction, planned = item.get("production_fraction"), item.get("planned_ns")
        if fraction is not None and planned:
            detail = f"replicate {replicate} at {fraction:.0%} of {planned:g} ns"
        else:
            detail = f"replicate {replicate} incomplete"
        if "running" in (item.get("unfinished_segments") or {}).values():
            detail += " (still running)"
        parts.append(detail)
    return f"{label}: {'; '.join(parts) or 'incomplete'}"


def summaries(conditions: Mapping[str, Mapping[str, Any] | None]) -> list[str]:
    """The :func:`summary` line of every partial condition, in order."""
    return [line for label, record in conditions.items() if (line := summary(label, record))]


def _planned_ns(sim_config: Any) -> float | None:
    """Planned production length of a condition, or ``None`` if its config has none."""
    try:
        return float(sim_config.simulation_phases.production.duration)
    except (AttributeError, TypeError, ValueError):
        return None


def _ranges(values: Sequence[int]) -> str:
    """Compact replicate list, for example ``1-3, 5``."""
    ordered = sorted(int(value) for value in values)
    runs: list[str] = []
    start = previous = None
    for value in ordered:
        if start is None:
            start = previous = value
        elif value == previous + 1:
            previous = value
        else:
            runs.append(f"{start}" if start == previous else f"{start}-{previous}")
            start = previous = value
    if start is not None:
        runs.append(f"{start}" if start == previous else f"{start}-{previous}")
    return ", ".join(runs)
