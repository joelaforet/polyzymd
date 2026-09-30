"""A comparison.yaml that still configures retired plugins keeps loading.

rg, rmsd, rmsf, sasa, secondary_structure, distances, catalytic_triad, contacts and
hydrogen_bonds left the plugin system for ``polyzymd analyze``. Each of their ``plugins``
blocks is ignored with one warning, and ``compare run <name>`` prints the equivalent
``polyzymd analyze`` command.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from polyzymd.cli.compare import compare
from polyzymd.config.comparison import ComparisonConfig
from tests._support.analysis_testkit import write_simulation_config

PAIR = {"label": "Ser-His", "selection_a": "resid 77 and name OG", "selection_b": "name NE2"}
RETIRED = {
    "rg": "radius_of_gyration",
    "rmsd": "rmsd",
    "rmsf": "rmsf",
    "sasa": "sasa",
    "secondary_structure": "dssp_occupancy",
    "distances": "pair_distance",
    "catalytic_triad": "pair_distance",
    "contacts": "residue_occlusion",
    "hydrogen_bonds": "hydrogen_bonds",
}


@pytest.fixture()
def comparison_file(tmp_path: Path) -> Path:
    """A comparison.yaml with a block for every retired plugin over two conditions."""
    conditions = []
    for label in ("No Polymer", "SBMA 50"):
        config = write_simulation_config(tmp_path / label.replace(" ", "_"), scratch=tmp_path)
        conditions.append({"label": label, "config": str(config), "replicates": [1, 2]})
    data = {
        "name": "retired_plugins",
        "control": "No Polymer",
        "conditions": conditions,
        "defaults": {"equilibration_time": "200ns"},
        "plugins": {
            "rg": {"runs": [{"label": "Protein", "selection": "protein"}]},
            "rmsd": {"runs": [{"label": "CA", "selection": "name CA"}]},
            "distances": {"pairs": [PAIR]},
            "catalytic_triad": {"threshold": 3.5, "pairs": [PAIR]},
            "rmsf": {"selection": "name CA"},
            "sasa": {"runs": [{"label": "Protein", "target_selection": "protein"}]},
            "secondary_structure": {"chain_id": "B"},
            "contacts": {"cutoff": 4.5, "polymer_selection": "resname SBM"},
            "hydrogen_bonds": {"distance_cutoff": 3.2},
        },
        "plot_settings": {
            "distances": {"use_kde": True},
            "catalytic_triad": {},
            "rmsf": {},
            "sasa": {},
            "secondary_structure": {},
            "contacts": {"enrichment_error_bar": "sem"},
            "hydrogen_bonds": {},
        },
    }
    path = tmp_path / "comparison.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    return path


def test_retired_blocks_are_ignored_with_one_warning_each(comparison_file: Path) -> None:
    """The file loads and each retired block warns once."""
    with pytest.warns(UserWarning, match="block, which is ignored") as record:
        config = ComparisonConfig.from_yaml(comparison_file)

    messages = [str(item.message) for item in record]
    for name, function in RETIRED.items():
        (message,) = [text for text in messages if f"plugins.{name} block" in text]
        assert f"polyzymd analyze {name} " in message
        ending = " (" if name in ("contacts", "hydrogen_bonds") else "."
        assert f"polyzymd.analyses.functions.{function}{ending}" in message
        assert ("--set pairs=" in message) == (function == "pair_distance")
        per_replicate = name in ("rmsf", "secondary_structure", "contacts", "hydrogen_bonds")
        method = "study.per_replicate" if per_replicate else "study.timeseries"
        assert f"in Python {method} with" in message
    (contacts,) = [text for text in messages if "plugins.contacts block" in text]
    assert (
        "functions.residue_occlusion (residue_contacts for --set method=distance, and "
        "contact_lifetimes for how long contacts last)." in contacts
    )
    (hbonds,) = [text for text in messages if "plugins.hydrogen_bonds block" in text]
    assert (
        "functions.hydrogen_bonds (hbond_lifetimes for how long hydrogen bonds last, "
        "residue_hbond_occupancy and residue_pair_hbond_occupancy" in hbonds
    )
    assert sum("plot_settings." in text for text in messages) == 7
    assert config.plugins.get_enabled_plugins() == []
    assert config.validate_config() == []


def test_compare_validate_passes(comparison_file: Path) -> None:
    """compare validate accepts a file whose plugin blocks are all retired."""
    validated = CliRunner().invoke(compare, ["validate", "-f", str(comparison_file)])
    assert validated.exit_code == 0, validated.output


@pytest.mark.parametrize("name", list(RETIRED))
def test_compare_run_prints_the_analyze_command(comparison_file: Path, name: str) -> None:
    """compare run of a retired name exits 1 with the polyzymd analyze command for this file."""
    result = CliRunner().invoke(compare, ["run", name, "-f", str(comparison_file)])

    assert result.exit_code == 1
    assert "no longer runs through 'polyzymd compare run'" in result.stderr
    assert "Unknown comparison type" not in result.stderr
    fix = next(line for line in result.stderr.splitlines() if line.startswith("fix: "))
    assert fix.startswith(f"fix: polyzymd analyze {name} -c ")
    assert "--label 'No Polymer'" in fix and "--label 'SBMA 50'" in fix
    pairs = " --set pairs=<pairs.yaml>" if RETIRED[name] == "pair_distance" else ""
    assert fix.endswith(f"--replicates 1,2 --eq 200ns{pairs}")


def test_compare_run_all_skips_retired(comparison_file, monkeypatch, tmp_path: Path) -> None:
    """compare run-all warns about the retired blocks and finds nothing left to run."""
    seen = {}

    def _run_all(config, **kwargs):
        seen["plugins"] = config.plugins.get_enabled_plugins()
        return {}

    monkeypatch.setattr("polyzymd.analyses.orchestrator.run_all_comparisons", _run_all)
    with pytest.warns(UserWarning, match="plugins.rmsd block, which is ignored"):
        result = CliRunner().invoke(compare, ["run-all", "-f", str(comparison_file)])

    assert result.exit_code == 1, result.output
    assert "No analyses are enabled in comparison.yaml." in result.output
    assert seen == {}
