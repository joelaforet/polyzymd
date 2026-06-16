"""Advisory pixi environment warnings for CLI commands."""

from __future__ import annotations

import os
from collections.abc import Iterable

import click

DISABLE_ENV_WARNINGS_VALUES = {"1", "true", "yes", "on"}


def get_active_pixi_environment() -> str | None:
    """Return the active pixi environment name, if pixi set one.

    Returns
    -------
    str or None
        Value of ``PIXI_ENVIRONMENT_NAME`` when present and non-empty.
    """
    active_env = os.environ.get("PIXI_ENVIRONMENT_NAME")
    if active_env is None:
        return None
    active_env = active_env.strip()
    return active_env or None


def env_warnings_disabled() -> bool:
    """Return whether PolyzyMD CLI environment warnings are disabled.

    Returns
    -------
    bool
        ``True`` when ``POLYZYMD_DISABLE_ENV_WARNINGS`` is a recognized truthy value.
    """
    value = os.environ.get("POLYZYMD_DISABLE_ENV_WARNINGS", "")
    return value.strip().lower() in DISABLE_ENV_WARNINGS_VALUES


def warn_if_wrong_pixi_env(
    command: str,
    recommended: str | Iterable[str],
    accepted: Iterable[str] | None = None,
    detail: str | None = None,
    suggestion: str | None = None,
) -> None:
    """Warn when a CLI command is run from an unexpected pixi environment.

    This helper is intentionally advisory. It never raises, changes exit codes,
    or emits warnings when no pixi environment is active.

    Parameters
    ----------
    command : str
        User-facing command name to include in the warning.
    recommended : str or iterable of str
        Preferred environment or environments for this command.
    accepted : iterable of str, optional
        Environments that should not warn. Defaults to ``recommended``.
    detail : str, optional
        Extra context to append to the warning.
    suggestion : str, optional
        Explicit suggestion text. Defaults to ``pixi run`` guidance.
    """
    active_env = get_active_pixi_environment()
    if active_env is None or env_warnings_disabled():
        return

    recommended_envs = _as_tuple(recommended)
    accepted_envs = _as_tuple(accepted) if accepted is not None else recommended_envs
    if active_env in accepted_envs:
        return

    recommended_text = _format_env_list(recommended_envs)
    message = (
        f"Warning: polyzymd {command} is running in pixi environment "
        f"'{active_env}', but recommended environment is {recommended_text}."
    )
    if detail:
        message = f"{message} {detail}"
    if suggestion is None:
        first_env = recommended_envs[0] if recommended_envs else "build"
        suggestion = f"Use 'pixi run -e {first_env} polyzymd {command} ...' or switch shells."
    message = f"{message}\n{suggestion}\nSet POLYZYMD_DISABLE_ENV_WARNINGS=1 to hide this advisory warning."
    click.secho(message, fg="yellow", err=True)


def _as_tuple(values: str | Iterable[str]) -> tuple[str, ...]:
    """Normalize a string or iterable of strings to a tuple."""
    if isinstance(values, str):
        return (values,)
    return tuple(values)


def _format_env_list(values: tuple[str, ...]) -> str:
    """Format environment names for user-facing warning text."""
    quoted = [f"'{value}'" for value in values]
    if not quoted:
        return "a different pixi environment"
    if len(quoted) == 1:
        return quoted[0]
    if len(quoted) == 2:
        return f"{quoted[0]} or {quoted[1]}"
    return f"{', '.join(quoted[:-1])}, or {quoted[-1]}"
