"""PolyzyMD comparison module with lazy public exports."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_MODULE_EXPORTS: dict[str, list[str]] = {
    "polyzymd.compare.core.base": [
        "ANOVASummary",
        "BaseComparator",
        "BaseComparisonResult",
        "BaseConditionSummary",
        "PairwiseComparison",
    ],
    "polyzymd.compare.core.registry": ["ComparatorRegistry"],
    "polyzymd.compare.config": ["ComparisonConfig", "ConditionConfig"],
    "polyzymd.compare.settings": [
        "RMSFAnalysisSettings",
        "RMSFComparisonSettings",
        "CatalyticTriadAnalysisSettings",
        "CatalyticTriadComparisonSettings",
        "ContactsAnalysisSettings",
        "ContactsComparisonSettings",
        "BindingFreeEnergyAnalysisSettings",
        "BindingFreeEnergyComparisonSettings",
        "PolymerAffinityScoreSettings",
        "PolymerAffinityScoreComparisonSettings",
    ],
    "polyzymd.compare.comparators.rmsf": ["RMSFComparator"],
    "polyzymd.compare.comparators.triad": ["TriadComparator"],
    "polyzymd.compare.comparators.contacts": ["ContactsComparator"],
    "polyzymd.compare.comparators.binding_free_energy": ["BindingFreeEnergyComparator"],
    "polyzymd.compare.comparators.polymer_affinity": ["PolymerAffinityScoreComparator"],
    "polyzymd.compare.formatters": [
        "format_console_table",
        "format_markdown",
        "format_result",
        "to_json",
    ],
    "polyzymd.compare.triad_formatters": [
        "format_triad_console_table",
        "format_triad_markdown",
        "format_triad_result",
        "triad_to_json",
    ],
    "polyzymd.compare.contacts_formatters": [
        "contacts_to_json",
        "format_contacts_console_table",
        "format_contacts_markdown",
        "format_contacts_result",
    ],
    "polyzymd.compare.binding_free_energy_formatters": ["format_bfe_result"],
    "polyzymd.compare.polymer_affinity_formatters": ["format_affinity_result"],
    "polyzymd.compare.results": [
        "AffinityScoreConditionSummary",
        "AffinityScorePairwiseEntry",
        "AggregateComparisonResult",
        "BindingFreeEnergyResult",
        "ComparisonResult",
        "ContactsComparisonResult",
        "ContactsConditionSummary",
        "ContactsPairwiseComparison",
        "FreeEnergyConditionSummary",
        "FreeEnergyEntry",
        "FreeEnergyPairwiseEntry",
        "PolymerAffinityScoreResult",
        "TriadComparisonResult",
        "TriadConditionSummary",
    ],
}

_EXPORTS = {
    name: (module_name, name) for module_name, names in _MODULE_EXPORTS.items() for name in names
}

__all__ = [
    "BaseComparator",
    "BaseComparisonResult",
    "BaseConditionSummary",
    "PairwiseComparison",
    "ANOVASummary",
    "ComparatorRegistry",
    "ComparisonConfig",
    "ConditionConfig",
    "RMSFAnalysisSettings",
    "RMSFComparisonSettings",
    "CatalyticTriadAnalysisSettings",
    "CatalyticTriadComparisonSettings",
    "ContactsAnalysisSettings",
    "ContactsComparisonSettings",
    "BindingFreeEnergyAnalysisSettings",
    "BindingFreeEnergyComparisonSettings",
    "RMSFComparator",
    "ComparisonResult",
    "format_console_table",
    "format_markdown",
    "format_result",
    "to_json",
    "TriadComparator",
    "TriadComparisonResult",
    "TriadConditionSummary",
    "format_triad_console_table",
    "format_triad_markdown",
    "format_triad_result",
    "triad_to_json",
    "ContactsComparator",
    "ContactsComparisonResult",
    "ContactsConditionSummary",
    "ContactsPairwiseComparison",
    "AggregateComparisonResult",
    "format_contacts_console_table",
    "format_contacts_markdown",
    "format_contacts_result",
    "contacts_to_json",
    "BindingFreeEnergyComparator",
    "BindingFreeEnergyResult",
    "FreeEnergyConditionSummary",
    "FreeEnergyEntry",
    "FreeEnergyPairwiseEntry",
    "format_bfe_result",
    "PolymerAffinityScoreComparator",
    "PolymerAffinityScoreResult",
    "AffinityScoreConditionSummary",
    "AffinityScorePairwiseEntry",
    "PolymerAffinityScoreSettings",
    "PolymerAffinityScoreComparisonSettings",
    "format_affinity_result",
]


def __getattr__(name: str) -> Any:
    """Lazily import compare exports to avoid importing optional deps at install time."""
    try:
        module_name, attr_name = _EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc

    module = import_module(module_name)
    return getattr(module, attr_name)


def __dir__() -> list[str]:
    """Return module attributes for tab completion and introspection."""
    return sorted(set(globals()) | set(__all__))
