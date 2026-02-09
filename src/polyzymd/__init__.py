"""
PolyzyMD: Molecular dynamics simulation toolkit for enzyme-polymer systems.

A comprehensive package for building, running, and analyzing molecular dynamics
simulations of enzymes with co-polymers, designed for reproducible computational
research with HPC cluster support.

Example usage:
    >>> from polyzymd import SimulationConfig, SystemBuilder
    >>> config = SimulationConfig.from_yaml("my_simulation.yaml")
    >>> builder = SystemBuilder(config)
    >>> system = builder.build()
    >>> system.run()

Key modules:
    - config: Configuration management with YAML support
    - builders: System construction (enzyme, substrate, polymer, solvent)
    - simulation: MD execution and continuation
    - workflow: SLURM job submission and daisy-chaining
    - analysis: Trajectory analysis tools
"""

__version__ = "1.0.0"
__author__ = "Joe Laforet Jr."
__email__ = "jola3134@colorado.edu"

# Core imports for convenience
from polyzymd.config.schema import SimulationConfig
from polyzymd.builders.system_builder import SystemBuilder
from polyzymd.simulation.runner import SimulationRunner
from polyzymd.simulation.continuation import ContinuationManager
from polyzymd.workflow.slurm import SlurmConfig, SlurmScriptGenerator
from polyzymd.workflow.daisy_chain import DaisyChainSubmitter, submit_daisy_chain

__all__ = [
    # Configuration
    "SimulationConfig",
    # Building
    "SystemBuilder",
    # Simulation
    "SimulationRunner",
    "ContinuationManager",
    # Workflow
    "SlurmConfig",
    "SlurmScriptGenerator",
    "DaisyChainSubmitter",
    "submit_daisy_chain",
    # Version
    "__version__",
]
