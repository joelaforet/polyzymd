"""Utilities for parsing and validating replicate range strings."""

from __future__ import annotations

import re


def parse_replicate_range(replicate_range: str) -> list[int]:
    """Parse a replicate range string into replicate numbers.

    Parameters
    ----------
    replicate_range : str
        Replicate range string (for example, ``"1-5"``, ``"1,3,5"``,
        ``"1-10:2"``).

    Returns
    -------
    list[int]
        Sorted and deduplicated replicate numbers.

    Raises
    ------
    ValueError
        If the range is empty, contains invalid entries, uses zero or negative
        replicate numbers, has reversed ranges, or has zero/negative step values.
    """
    if not replicate_range or not replicate_range.strip():
        raise ValueError("Replicate range cannot be empty")

    replicates: list[int] = []
    parts = replicate_range.split(",")

    for part in parts:
        part = part.strip()
        if not part:
            raise ValueError(f"Invalid replicate entry in range: {replicate_range}")

        if "-" in part:
            if ":" in part:
                range_part, step_text = part.split(":", maxsplit=1)
                step = int(step_text)
            else:
                range_part = part
                step = 1

            if step <= 0:
                raise ValueError(f"Replicate range step must be positive: {part}")

            bounds = range_part.split("-", maxsplit=1)
            if len(bounds) != 2 or not bounds[0] or not bounds[1]:
                raise ValueError(f"Invalid replicate range segment: {part}")

            start, end = map(int, bounds)
            if start <= 0 or end <= 0:
                raise ValueError(f"Replicate numbers must be positive: {part}")
            if start > end:
                raise ValueError(f"Replicate range start must be <= end: {part}")

            replicates.extend(range(start, end + 1, step))
        else:
            value = int(part)
            if value <= 0:
                raise ValueError(f"Replicate number must be positive: {part}")
            replicates.append(value)

    return sorted(list(set(replicates)))


def validate_replicate_range(replicate_range: str) -> bool:
    """Validate the format of a replicate range string.

    Parameters
    ----------
    replicate_range : str
        Range string to validate.

    Returns
    -------
    bool
        ``True`` when the format is valid.

    Raises
    ------
    ValueError
        If the range string format is invalid.
    """
    pattern = r"^(\d+(-\d+(:\d+)?)?)(,\d+(-\d+(:\d+)?)?)*$"
    if re.fullmatch(pattern, replicate_range) is None:
        raise ValueError(f"Invalid replicate range format: {replicate_range}")
    return True
