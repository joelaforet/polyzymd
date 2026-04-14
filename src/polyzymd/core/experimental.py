"""Helpers for labeling experimental analyses and plots."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Final

EXPERIMENTAL_FEATURES: Final[dict[str, str]] = {
    "contacts_binding_preference": "Contacts binding preference",
    "exposure": "Exposure dynamics",
    "binding_free_energy": "Binding free energy",
    "polymer_affinity": "Polymer affinity score",
    "polymer_bridging_chemistry": "Polymer bridging chemistry profiling",
}

EXPERIMENTAL_COMPARISON_TYPES: Final[dict[str, tuple[str, ...]]] = {
    "exposure": ("exposure",),
    "binding_free_energy": ("binding_free_energy",),
    "polymer_affinity": ("polymer_affinity",),
}

EXPERIMENTAL_PLOT_TYPES: Final[dict[str, tuple[str, ...]]] = {
    "binding_preference_heatmap": ("contacts_binding_preference",),
    "binding_preference_bars": ("contacts_binding_preference",),
    "system_coverage_heatmap": ("contacts_binding_preference",),
    "system_coverage_bars": ("contacts_binding_preference",),
    "user_partition_bars": ("contacts_binding_preference",),
    "exposure_chaperone_fraction": ("exposure",),
    "exposure_enrichment_heatmap": ("exposure",),
    "bfe_heatmap": ("binding_free_energy",),
    "bfe_bars": ("binding_free_energy",),
    "affinity_stacked_bars": ("polymer_affinity",),
    "affinity_group_bars": ("polymer_affinity",),
}


def normalize_experimental_features(features: Iterable[str] | None) -> tuple[str, ...]:
    """Return de-duplicated experimental feature ids in stable order."""
    if features is None:
        return ()

    normalized: list[str] = []
    seen: set[str] = set()
    for feature in features:
        if feature not in EXPERIMENTAL_FEATURES or feature in seen:
            continue
        normalized.append(feature)
        seen.add(feature)
    return tuple(normalized)


def experimental_feature_labels(features: Iterable[str] | None) -> list[str]:
    """Return human-readable labels for experimental features."""
    return [EXPERIMENTAL_FEATURES[feature] for feature in normalize_experimental_features(features)]


def format_experimental_suffix(features: Iterable[str] | None) -> str:
    """Return a short suffix for experimental items in CLI listings."""
    return " [experimental]" if normalize_experimental_features(features) else ""


def experimental_features_for_comparison_type(
    comparison_type: str,
    analysis_settings: object | None = None,
) -> tuple[str, ...]:
    """Return experimental feature ids for a comparison type."""
    if comparison_type == "contacts" and getattr(
        analysis_settings, "compute_binding_preference", False
    ):
        return ("contacts_binding_preference",)
    return EXPERIMENTAL_COMPARISON_TYPES.get(comparison_type, ())


def experimental_features_for_plot_type(plot_type: str) -> tuple[str, ...]:
    """Return experimental feature ids for a registered plot type."""
    return EXPERIMENTAL_PLOT_TYPES.get(plot_type, ())


def experimental_warning_text(
    features: Iterable[str] | None,
    *,
    markdown: bool = False,
) -> str:
    """Build a standard warning banner for experimental analyses."""
    labels = experimental_feature_labels(features)
    if not labels:
        return ""

    affected = ", ".join(labels)
    note = "Definitions and interpretation may change after the presentation release."
    if markdown:
        return f"> WARNING: Experimental analysis\n> {note}\n> Affected: {affected}"
    return f"WARNING: Experimental analysis\n{note}\nAffected: {affected}"


def prefix_experimental_output(
    body: str,
    features: Iterable[str] | None,
    output_format: str,
) -> str:
    """Prefix formatted output with an experimental warning when appropriate."""
    normalized = normalize_experimental_features(features)
    if not normalized or output_format == "json":
        return body

    banner = experimental_warning_text(normalized, markdown=output_format == "markdown")
    if not body:
        return banner
    return f"{banner}\n\n{body}"


def echo_experimental_warning(features: Iterable[str] | None) -> None:
    """Print a standardized CLI warning for experimental analyses."""
    labels = experimental_feature_labels(features)
    if not labels:
        return

    import click

    click.secho("WARNING: Experimental analysis", fg="yellow", bold=True)
    click.echo("Definitions and interpretation may change after the presentation release.")
    click.echo(f"Affected: {', '.join(labels)}")
    click.echo()


def annotate_experimental_figure(fig: object, features: Sequence[str] | None) -> None:
    """Stamp a saved figure with a visible experimental tag."""
    if not normalize_experimental_features(features):
        return

    fig.text(
        0.01,
        0.99,
        "EXPERIMENTAL",
        fontsize=8,
        fontweight="bold",
        color="darkred",
        ha="left",
        va="top",
        bbox={"boxstyle": "round,pad=0.2", "facecolor": "mistyrose", "edgecolor": "darkred"},
    )
