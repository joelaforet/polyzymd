"""Aggregate result validation helpers for analysis plugins."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

from pydantic import BaseModel

from polyzymd.analyses._framework.cache_identity import settings_fingerprint

if TYPE_CHECKING:
    from polyzymd.analyses._framework.contexts import Condition


_UNKNOWN_VALUES = {None, "", "unknown"}


class AggregateValidationError(ValueError):
    """Raised when an aggregated result is stale for the active context."""


def aggregate_settings_fingerprint(settings: BaseModel | None) -> str | None:
    """Return the default settings fingerprint for aggregate validation.

    Parameters
    ----------
    settings : BaseModel or None
        Analysis settings for the active run.

    Returns
    -------
    str or None
        Short deterministic settings fingerprint, or ``None`` when no settings
        model is available.
    """

    if settings is None:
        return None
    return settings_fingerprint(settings)


def validate_aggregated_result(
    analysis: Any,
    result: Any,
    *,
    condition: Condition | None,
    settings: BaseModel | None,
    equilibration: str,
    source: str | Path | None = None,
    expected_replicates: Sequence[int] | None = None,
    allow_replicate_subset: bool = False,
) -> Any:
    """Validate and coerce an aggregated result for the active context.

    Parameters
    ----------
    analysis : Any
        Analysis plugin instance that owns result model configuration.
    result : Any
        Loaded or newly computed aggregate result.
    condition : Condition or None
        Condition providing simulation configuration and label context.
    settings : BaseModel or None
        Active analysis settings.
    equilibration : str
        Active equilibration window.
    source : str or Path or None, optional
        File path or description used in diagnostics.
    expected_replicates : sequence of int or None, optional
        Replicate IDs expected in the aggregate.
    allow_replicate_subset : bool, optional
        Accept a non-empty subset of ``expected_replicates`` when ``True``.

    Returns
    -------
    Any
        The validated result, coerced through ``AggregatedResultClass`` when
        applicable.

    Raises
    ------
    AggregateValidationError
        If the aggregate cannot prove compatibility with the active context.
    """

    coerced = _coerce_aggregate_result(analysis, result, source=source)
    _validate_config_hash(coerced, condition=condition, analysis_name=analysis.name, source=source)
    _validate_equilibration(
        coerced,
        equilibration=equilibration,
        analysis_name=analysis.name,
        source=source,
    )
    _validate_settings_fingerprint(
        analysis,
        coerced,
        settings=settings,
        source=source,
    )
    _validate_replicates(
        analysis,
        coerced,
        expected_replicates=expected_replicates,
        allow_replicate_subset=allow_replicate_subset,
        source=source,
    )
    return coerced


def _coerce_aggregate_result(analysis: Any, result: Any, *, source: str | Path | None) -> Any:
    """Coerce dict aggregate payloads through a plugin result model."""

    result_cls = type(analysis).AggregatedResultClass
    if result_cls is None or not isinstance(result, dict):
        return result
    try:
        if hasattr(result_cls, "model_validate"):
            return result_cls.model_validate(result)
        return result_cls(**result)
    except Exception as exc:
        raise _validation_error(
            analysis.name,
            source,
            f"could not parse aggregate as {result_cls.__name__}: {exc}",
        ) from exc


def _validate_config_hash(
    result: Any,
    *,
    condition: Condition | None,
    analysis_name: str,
    source: str | Path | None,
) -> None:
    """Validate stored config hash when the aggregate provides one."""

    stored_hash = _field_value(result, "config_hash")
    if stored_hash in _UNKNOWN_VALUES or condition is None:
        return

    from polyzymd.analyses._framework.cache_identity import validate_config_hash

    if not validate_config_hash(str(stored_hash), condition.sim_config):
        raise _validation_error(
            analysis_name,
            source,
            "config hash does not match the active condition config",
        )


def _validate_equilibration(
    result: Any,
    *,
    equilibration: str,
    analysis_name: str,
    source: str | Path | None,
) -> None:
    """Validate stored equilibration metadata when present."""

    stored_time = _field_value(result, "equilibration_time")
    stored_unit = _field_value(result, "equilibration_unit")
    if stored_time is None and stored_unit is None:
        return
    if stored_time is None or stored_unit is None:
        raise _validation_error(
            analysis_name,
            source,
            "aggregate has incomplete equilibration metadata",
        )

    from polyzymd.analyses.shared.loader import convert_time, parse_time_string

    try:
        expected_value, expected_unit = parse_time_string(equilibration)
        stored_ps = convert_time(float(stored_time), str(stored_unit), "ps")
        expected_ps = convert_time(float(expected_value), str(expected_unit), "ps")
    except (TypeError, ValueError) as exc:
        raise _validation_error(
            analysis_name,
            source,
            f"aggregate has invalid equilibration metadata: {exc}",
        ) from exc

    if abs(stored_ps - expected_ps) > 1.0e-9:
        raise _validation_error(
            analysis_name,
            source,
            f"equilibration mismatch: stored {float(stored_time):g}{stored_unit}, "
            f"requested {equilibration}",
        )


def _validate_settings_fingerprint(
    analysis: Any,
    result: Any,
    *,
    settings: BaseModel | None,
    source: str | Path | None,
) -> None:
    """Validate aggregate settings identity."""

    expected_fingerprint = analysis.aggregate_settings_fingerprint(settings)
    if expected_fingerprint is None:
        return

    stored_fingerprint = _settings_fingerprint_from_result(result)
    if stored_fingerprint is None:
        raise _validation_error(
            analysis.name,
            source,
            "aggregate is missing settings fingerprint metadata",
        )
    if str(stored_fingerprint) != str(expected_fingerprint):
        raise _validation_error(
            analysis.name,
            source,
            f"settings fingerprint mismatch: stored {stored_fingerprint}, "
            f"current {expected_fingerprint}",
        )


def _validate_replicates(
    analysis: Any,
    result: Any,
    *,
    expected_replicates: Sequence[int] | None,
    allow_replicate_subset: bool,
    source: str | Path | None,
) -> None:
    """Validate aggregate replicate identity and replicate count."""

    if expected_replicates is None:
        return

    expected = tuple(sorted(int(rep) for rep in expected_replicates))
    expected_set = set(expected)
    stored = _replicates_from_result(result)
    stored_count = _field_value(result, "n_replicates")
    if stored_count is None:
        stored_count = _metadata_value(result, "n_replicates")

    if stored is None and stored_count is None:
        raise _validation_error(
            analysis.name,
            source,
            "aggregate is missing replicate identity metadata",
        )

    if stored is not None:
        _validate_stored_replicate_ids(
            analysis,
            stored,
            expected=expected,
            expected_set=expected_set,
            allow_replicate_subset=allow_replicate_subset,
            source=source,
        )
        if stored_count is not None:
            _validate_replicate_count(
                analysis,
                stored_count,
                expected_count=len(stored),
                source=source,
                count_context="stored replicate IDs",
            )
        return

    if allow_replicate_subset:
        try:
            declared_count = int(stored_count)
        except (TypeError, ValueError) as exc:
            raise _validation_error(
                analysis.name,
                source,
                f"aggregate has invalid n_replicates={stored_count!r}",
            ) from exc
        if declared_count < analysis.min_replicates or declared_count > len(expected):
            raise _validation_error(
                analysis.name,
                source,
                f"n_replicates={declared_count} is not an allowed subset count for "
                f"expected replicates {list(expected)}",
            )
        return

    _validate_replicate_count(
        analysis,
        stored_count,
        expected_count=len(expected),
        source=source,
        count_context="expected replicates",
    )


def _validate_stored_replicate_ids(
    analysis: Any,
    stored_replicates: Sequence[int],
    *,
    expected: tuple[int, ...],
    expected_set: set[int],
    allow_replicate_subset: bool,
    source: str | Path | None,
) -> None:
    """Validate stored replicate IDs against expected IDs."""

    stored = tuple(sorted(int(rep) for rep in stored_replicates))
    if not stored:
        raise _validation_error(analysis.name, source, "aggregate has empty replicate identity")
    if len(set(stored)) != len(stored):
        raise _validation_error(analysis.name, source, "aggregate has duplicate replicate IDs")

    is_allowed_subset = allow_replicate_subset and set(stored).issubset(expected_set)
    if stored != expected and not is_allowed_subset:
        raise _validation_error(
            analysis.name,
            source,
            f"replicate mismatch: stored {list(stored)}, expected {list(expected)}",
        )
    if allow_replicate_subset and len(stored) < analysis.min_replicates:
        raise _validation_error(
            analysis.name,
            source,
            f"stored subset has {len(stored)} replicates below minimum {analysis.min_replicates}",
        )


def _validate_replicate_count(
    analysis: Any,
    stored_count: Any,
    *,
    expected_count: int,
    source: str | Path | None,
    count_context: str,
) -> None:
    """Validate a declared aggregate replicate count."""

    try:
        declared_count = int(stored_count)
    except (TypeError, ValueError) as exc:
        raise _validation_error(
            analysis.name,
            source,
            f"aggregate has invalid n_replicates={stored_count!r}",
        ) from exc
    if declared_count != expected_count:
        raise _validation_error(
            analysis.name,
            source,
            f"n_replicates={declared_count} does not match {count_context} count {expected_count}",
        )


def _settings_fingerprint_from_result(result: Any) -> str | None:
    """Return the settings fingerprint stored on an aggregate result."""

    for field_name in ("settings_fingerprint", "settings_fp"):
        value = _field_value(result, field_name)
        if value is not None:
            return str(value)
    for field_name in ("settings_fingerprint", "settings_fp"):
        value = _metadata_value(result, field_name)
        if value is not None:
            return str(value)
    return None


def _replicates_from_result(result: Any) -> tuple[int, ...] | None:
    """Return aggregate replicate IDs when the result declares them."""

    for field_name in ("replicates", "replicate_ids"):
        value = _field_value(result, field_name)
        if value is not None:
            try:
                return tuple(int(rep) for rep in value)
            except (TypeError, ValueError) as exc:
                raise AggregateValidationError(
                    "Aggregate result has invalid replicate identity. Recompute the analysis "
                    "or clear stale analysis cache files."
                ) from exc
    for field_name in ("replicates", "replicate_ids"):
        value = _metadata_value(result, field_name)
        if value is not None:
            try:
                return tuple(int(rep) for rep in value)
            except (TypeError, ValueError) as exc:
                raise AggregateValidationError(
                    "Aggregate result metadata has invalid replicate identity. Recompute the "
                    "analysis or clear stale analysis cache files."
                ) from exc
    return None


def _field_value(result: Any, field_name: str) -> Any:
    """Return a top-level result field from models or mappings."""

    if isinstance(result, dict):
        return result.get(field_name)
    return getattr(result, field_name, None)


def _metadata_value(result: Any, field_name: str) -> Any:
    """Return a metadata field from models or mappings."""

    metadata = _field_value(result, "metadata")
    if isinstance(metadata, dict):
        return metadata.get(field_name)
    return None


def _validation_error(
    analysis_name: str,
    source: str | Path | None,
    detail: str,
) -> AggregateValidationError:
    """Build a user-actionable aggregate validation error."""

    source_text = f" at {source}" if source is not None else ""
    return AggregateValidationError(
        f"{analysis_name}: stale or incompatible aggregated result{source_text}: {detail}. "
        "Recompute the analysis with --recompute or clear stale analysis cache files."
    )
