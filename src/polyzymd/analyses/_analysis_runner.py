"""Runner-backed replicate execution helpers for analyses."""

from __future__ import annotations

import logging
from typing import Any

from polyzymd.analyses.exceptions import PluginContractError

logger = logging.getLogger("polyzymd.analyses")

_RUNNER_NOT_CONFIGURED = object()


def run_replicate_default(analysis: Any, ctx: Any, replicate: int, *, base_cls: type) -> Any:
    """Run the default per-replicate dispatch path.

    Parameters
    ----------
    analysis : Any
        Analysis instance.
    ctx : Any
        Replicate context.
    replicate : int
        Replicate number.
    base_cls : type
        Public facade class used for method identity checks.

    Returns
    -------
    Any
        Replicate result, or ``None`` when compute is disabled.
    """
    if not type(analysis).has_compute_stage:
        return None
    result = analysis._run_replicate_via_runner(ctx, replicate)
    if result is not _RUNNER_NOT_CONFIGURED:
        return result
    raise NotImplementedError(
        f"{type(analysis).__name__} public plugins must implement "
        "build_runner() + summarize_replicate() or set has_compute_stage = False; "
        "direct run_replicate() overrides are advanced/internal only."
    )


def run_replicate_via_runner(analysis: Any, ctx: Any, replicate: int, *, base_cls: type) -> Any:
    """Execute the opt-in runner-backed replicate path.

    Parameters
    ----------
    analysis : Any
        Analysis instance.
    ctx : Any
        Replicate context.
    replicate : int
        Replicate number.
    base_cls : type
        Public facade class used for method identity checks.

    Returns
    -------
    Any
        Replicate result, or a sentinel when the subclass did not opt in.
    """
    if type(analysis).build_runner is base_cls.build_runner:
        return _RUNNER_NOT_CONFIGURED

    loader = analysis._trajectory_loader_factory()(ctx.sim_config)
    universe = loader.load_universe(replicate)
    window = analysis.get_trajectory_window(ctx, replicate, loader, universe)
    if getattr(window, "warning_message", None):
        logger.warning(
            "%s: %s [condition=%s, replicate=%d]",
            analysis.name,
            window.warning_message,
            ctx.condition.label,
            replicate,
        )

    runner = analysis.build_runner(ctx, replicate, universe, window)
    if runner is None:
        raise PluginContractError(
            f"{type(analysis).__name__}.build_runner() returned None. "
            "Return a runner from build_runner(); direct run_replicate() overrides are "
            "advanced/internal only."
        )
    run_method = getattr(runner, "run", None)
    if not callable(run_method):
        raise PluginContractError(
            f"{type(analysis).__name__}.build_runner() must return an object with callable run(), "
            f"got {type(runner).__name__}"
        )

    executed_runner = run_method(**window.run_kwargs())
    if not hasattr(executed_runner, "results"):
        executed_runner = runner
    if not hasattr(executed_runner, "results"):
        raise PluginContractError(
            f"Runner returned by {type(analysis).__name__}.build_runner() must expose results "
            "after run()."
        )
    return analysis.summarize_replicate(ctx, replicate, executed_runner, window)
