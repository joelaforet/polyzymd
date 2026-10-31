"""The registered checks; importing this package registers every one of them."""

from polyzymd.audit.checks import (  # noqa: F401
    engine,
    geometry,
    modifications,
    parameters,
    restraints,
    system,
)
from polyzymd.audit.registry import CheckSpec, register


@register(
    CheckSpec(
        "modeling_exceptions.unused",
        "system",
        "WARN",
        "a modeling_exceptions entry that accepted nothing (stale in a published config)",
    )
)
def _unused_exceptions(context):  # reported by polyzymd.audit.report.audit
    return []
