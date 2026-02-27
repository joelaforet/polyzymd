"""
PolyzyMD: Molecular dynamics simulation toolkit for enzyme-polymer systems.

A comprehensive package for building, running, and analyzing molecular dynamics
simulations of enzymes with co-polymers, designed for reproducible computational
research with HPC cluster support.

Example usage:
    >>> from polyzymd.config import SimulationConfig
    >>> config = SimulationConfig.from_yaml("my_simulation.yaml")

    >>> from polyzymd import SystemBuilder
    >>> builder = SystemBuilder(config)
    >>> system = builder.build()

Key modules:
    - config: Configuration management with YAML support
    - builders: System construction (enzyme, substrate, polymer, solvent)
    - simulation: MD execution and continuation
    - workflow: SLURM job submission and daisy-chaining
    - exporters: GROMACS export functionality

Note:
    This package uses lazy imports for heavy dependencies (OpenMM, OpenFF, etc.).
    The config module can be imported without these dependencies, but simulation
    modules require a full conda environment with OpenMM and OpenFF installed.
"""

__version__ = "1.0.1"
__author__ = "Joseph R. Laforet Jr."
__email__ = "jola3134@colorado.edu"

# Define what's available for lazy import
__all__ = [
    # Version info
    "__version__",
    "__author__",
    "__email__",
    # Configuration (lightweight, always available)
    "SimulationConfig",
    # Building (requires OpenFF - lazy loaded)
    "SystemBuilder",
    # Simulation (requires OpenMM - lazy loaded)
    "SimulationRunner",
    "ContinuationManager",
    # Workflow (lightweight)
    "SlurmConfig",
    "SlurmScriptGenerator",
    "DaisyChainSubmitter",
    "submit_daisy_chain",
]


def __getattr__(name: str):
    """
    Lazy import heavy modules only when accessed.

    This allows the package to be imported without OpenMM/OpenFF installed,
    which is useful for:
    - Running basic CI tests (lint, config validation)
    - Building documentation
    - Inspecting package metadata

    The heavy imports only happen when you actually try to use classes
    like SystemBuilder or SimulationRunner.
    """
    # Configuration - lightweight, can always be imported
    if name == "SimulationConfig":
        from polyzymd.config.schema import SimulationConfig

        return SimulationConfig

    # Builders - require OpenFF
    if name == "SystemBuilder":
        from polyzymd.builders.system_builder import SystemBuilder

        return SystemBuilder

    # Simulation - require OpenMM
    if name == "SimulationRunner":
        from polyzymd.simulation.runner import SimulationRunner

        return SimulationRunner

    if name == "ContinuationManager":
        from polyzymd.simulation.continuation import ContinuationManager

        return ContinuationManager

    # Workflow - lightweight
    if name == "SlurmConfig":
        from polyzymd.workflow.slurm import SlurmConfig

        return SlurmConfig

    if name == "SlurmScriptGenerator":
        from polyzymd.workflow.slurm import SlurmScriptGenerator

        return SlurmScriptGenerator

    if name == "DaisyChainSubmitter":
        from polyzymd.workflow.daisy_chain import DaisyChainSubmitter

        return DaisyChainSubmitter

    if name == "submit_daisy_chain":
        from polyzymd.workflow.daisy_chain import submit_daisy_chain

        return submit_daisy_chain

    raise AttributeError(f"module 'polyzymd' has no attribute {name!r}")


def __dir__():
    """Return list of available attributes for tab completion."""
    return __all__
