API Overview
============

PolyzyMD is organized into the following modules:

Package Structure
-----------------

::

    polyzymd/
    ├── config/           # Configuration and YAML loading
    │   ├── schema.py     # Pydantic models for all config sections
    │   ├── comparison.py # ComparisonConfig, PlotSettings, condition models
    │   └── loader.py     # YAML loading utilities
    ├── builders/         # System building components
    │   ├── enzyme.py     # Enzyme/protein preparation
    │   ├── substrate.py  # Substrate/ligand handling
    │   ├── polymer.py    # Polymer chain generation
    │   ├── solvent.py    # Solvation and ion addition
    │   └── system_builder.py  # Main system assembly
    ├── simulation/       # MD simulation execution
    │   ├── runner.py     # Simulation runner
    │   ├── continuation.py    # Checkpoint continuation
    │   ├── progress.py   # Segment progress tracking
    │   └── signals.py    # SLURM signal handling
    ├── workflow/         # HPC workflow management
    │   ├── slurm.py      # SLURM script generation
    │   └── daisy_chain.py     # Job submission
    ├── core/             # Core utilities
    │   ├── parameters.py # Simulation parameters
    │   └── restraints.py # Restraint definitions
    ├── analyses/         # ★ Plugin system — unified analysis lifecycle
    │   ├── base.py       # Public facade, contexts, result models
    │   ├── discovery.py  # pkgutil-based auto-discovery
    │   ├── orchestrator.py  # Framework engine
    │   ├── stats.py      # Shared statistical utilities
    │   ├── _framework/   # Private/internal lifecycle and artifact internals
    │   ├── mda/          # Public MDAnalysis job and artifact layer
    │   ├── shared/       # Reusable utilities (TrajectoryLoader, alignment, etc.)
    │   ├── rmsd/         # RMSD plugin package
    │   ├── rg/           # Rg plugin package
    │   ├── rmsf/         # RMSF plugin package
    │   ├── contacts/     # Contacts plugin package
    │   ├── distances/    # Distance analysis plugin package
    │   ├── secondary_structure/  # Secondary structure plugin package
    │   └── ...           # Single-file or package plugins for each analysis type
    └── cli/              # Command-line interface
        ├── compare.py    # `polyzymd compare` subcommands
        └── main.py       # Click CLI

The ``analyses/_framework/`` package is private/internal implementation detail.
Contributor plugins should import public lifecycle symbols from
``polyzymd.analyses.base`` and MDAnalysis job/artifact utilities from
``polyzymd.analyses.mda``.


Key Classes
-----------

Configuration
~~~~~~~~~~~~~

- :py:class:`~polyzymd.config.schema.SimulationConfig` - Main configuration container
- :py:class:`~polyzymd.config.schema.EnzymeConfig` - Enzyme settings
- :py:class:`~polyzymd.config.schema.PolymerConfig` - Polymer settings
- :py:class:`~polyzymd.config.schema.OutputConfig` - Output directory settings

Building
~~~~~~~~

- :py:class:`~polyzymd.builders.system_builder.SystemBuilder` - Main system builder
- :py:class:`~polyzymd.builders.enzyme.EnzymeBuilder` - Enzyme preparation
- :py:class:`~polyzymd.builders.polymer.PolymerBuilder` - Polymer generation

Simulation
~~~~~~~~~~

- :py:class:`~polyzymd.simulation.runner.SimulationRunner` - Run simulations
- :py:class:`~polyzymd.simulation.continuation.ContinuationManager` - Continue from checkpoint
- :py:class:`~polyzymd.simulation.progress.SimulationProgress` - Track segment completion across jobs

Workflow
~~~~~~~~

- :py:class:`~polyzymd.workflow.daisy_chain.DaisyChainSubmitter` - Self-resubmitting SLURM job submission
- :py:class:`~polyzymd.workflow.slurm.SlurmConfig` - SLURM configuration

Restraints
~~~~~~~~~~

- :py:class:`~polyzymd.core.restraints.RestraintDefinition` - Restraint specification
- :py:class:`~polyzymd.core.restraints.AtomSelection` - Atom selection
- :py:class:`~polyzymd.core.restraints.RestraintFactory` - Create restraints from config

Analysis
~~~~~~~~

- :py:class:`~polyzymd.analyses.base.Analysis` - Plugin base class for all analyses
- :py:class:`~polyzymd.analyses.mda.job.MDAAnalysisJob` - Unit of trajectory-native MDAnalysis work
- :py:class:`~polyzymd.analyses.mda.frame_selection.FrameSelection` - Frame selection and equilibration-window descriptor
- :py:class:`~polyzymd.analyses.mda.artifacts.ReplicateArtifact` - Validated output for one replicate
- :py:class:`~polyzymd.analyses.mda.artifacts.ConditionArtifact` - Aggregated output for one condition
- :py:class:`~polyzymd.analyses.mda.artifacts.ComparisonArtifact` - Canonical comparison output for cross-condition results
- :py:class:`~polyzymd.analyses.mda.store.ArtifactStore` - Canonical artifact persistence and loading helper
- :py:class:`~polyzymd.analyses.base.MetricValue` - Scalar metric descriptor for default comparisons
- :py:class:`~polyzymd.analyses.contacts.ParallelContactAnalyzer` - Polymer-protein contacts

Comparison
~~~~~~~~~~

- :py:class:`~polyzymd.analyses.base.Analysis` - Plugin base class for all analyses
- :py:class:`~polyzymd.analyses.base.ComparisonResult` - Universal comparison result model
- :py:class:`~polyzymd.analyses.base.MetricValue` - Scalar metric descriptor
- :py:class:`~polyzymd.analyses.base.ReplicateContext` - Context for per-replicate computation
- :py:class:`~polyzymd.analyses.base.ComparisonContext` - Context for cross-condition comparison


Quick Reference
---------------

Load Configuration
~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from polyzymd.config.schema import SimulationConfig

    config = SimulationConfig.from_yaml("config.yaml")
    print(config.enzyme.name)

Build System
~~~~~~~~~~~~

.. code-block:: python

    from polyzymd.builders.system_builder import SystemBuilder

    builder = SystemBuilder(config)
    interchange = builder.build(replicate=1)

Run Simulation
~~~~~~~~~~~~~~

.. code-block:: python

    from polyzymd.simulation.runner import SimulationRunner

    runner = SimulationRunner(interchange, working_dir, config)
    runner.run_equilibration()
    runner.run_production(segment_index=0)

Submit to SLURM
~~~~~~~~~~~~~~~

.. code-block:: python

    from polyzymd.workflow.daisy_chain import submit_daisy_chain

    results = submit_daisy_chain(
        config_path="config.yaml",
        slurm_preset="aa100",
        replicates="1-5",
    )
