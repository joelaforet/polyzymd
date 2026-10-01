"""The config hash that stored study results record must not change.

The expected values were computed with ``compute_config_hash`` before it moved
from ``polyzymd.analyses._framework.cache_identity`` to
``polyzymd.analyses.identity``. A different value would make every stored
result look as if it came from another config.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace as NS

from polyzymd.analyses.identity import compute_config_hash

#: The config.yaml of a LipA SBMA condition whose stored hydrogen-bond
#: results record the hash REAL_CONFIG_HASH. Only its paths are read, so the
#: files need not exist.
REAL_CONFIG = """\
name: sbma
engine: openmm
enzyme:
  name: LipA
  pdb_path: /home/joelaforet/Shirts-Lab-Linux/polyzymd-realdata/1ISP_clean_processed_moved_simulation_resids.pdb
thermodynamics:
  temperature: 363.0
simulation_phases:
  equilibration_stages:
    - name: heating
      duration: 0.3636
      temperature: 363.0
      ensemble: NVT
  production:
    ensemble: NPT
    duration: 1000.0
    samples: 2500
    checkpoint_interval: 60.0
output:
  projects_directory: /home/joelaforet/.claude/jobs/0e78795e/tmp/slice_rg_data/sbma
  scratch_directory: /home/joelaforet/Shirts-Lab-Linux/polyzymd-realdata
  naming_template: "LipA_ResorufinButyrate_SBMA-EGMA_A50_B50_1000ns_363K_run{replicate}"
"""
REAL_CONFIG_HASH = "4f53486e7bb4211f"


def _config(*, substrate: bool = True, polymers: bool = True) -> NS:
    """Return an object with every field compute_config_hash reads."""
    return NS(
        name="LipA_SBMA_363K",
        enzyme=NS(name="LipA", pdb_path=Path("/data/structures/LipA.pdb")),
        thermodynamics=NS(temperature=363.0, pressure=1.0),
        output=NS(
            projects_directory=Path("/projects/lipa"),
            effective_scratch_directory=Path("/scratch/lipa"),
            naming_template="{enzyme}_{temperature}K_run{replicate}",
        ),
        substrate=(
            NS(name="Resorufin-Butyrate", sdf_path=Path("/data/structures/rb.sdf"))
            if substrate
            else None
        ),
        polymers=NS(
            enabled=polymers,
            type_prefix="SBMA-EGPMA",
            length=5,
            count=2,
            monomers=[
                NS(label="A", probability=0.75, name="SBMA"),
                NS(label="B", probability=0.25, name="EGPMA"),
            ],
        ),
    )


def test_hash_matches_the_value_before_the_move() -> None:
    """The fields of a polymer condition with a substrate hash to the recorded value."""
    assert compute_config_hash(_config()) == "ebc33187070b0d6a"


def test_hash_without_substrate_or_polymers_matches_the_value_before_the_move() -> None:
    """Leaving out the substrate and disabling polymers hash to the recorded value."""
    assert compute_config_hash(_config(substrate=False, polymers=False)) == "bab071deeacfad10"


def test_real_config_hash_matches_stored_results(tmp_path: Path) -> None:
    """A real simulation config hashes to the value its stored results record."""
    from polyzymd.config.schema import SimulationConfig

    path = tmp_path / "config.yaml"
    path.write_text(REAL_CONFIG)
    assert compute_config_hash(SimulationConfig.from_yaml(path)) == REAL_CONFIG_HASH
