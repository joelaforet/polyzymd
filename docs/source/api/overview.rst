API Overview
============

PolyzyMD is organized into the following modules:

Package Structure
-----------------

::

    polyzymd/
    ├── config/           # Configuration and YAML loading
    │   ├── schema.py     # Pydantic models for all config sections
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
    │   ├── base.py       # Analysis ABC, context objects, result models
    │   ├── discovery.py  # pkgutil-based auto-discovery
    │   ├── orchestrator.py  # Framework engine
    │   ├── stats.py      # Shared statistical utilities
    │   ├── shared/       # Reusable utilities (TrajectoryLoader, alignment, etc.)
    │   ├── rmsd/         # RMSD plugin package
    │   ├── rg/           # Rg plugin package
    │   ├── rmsf/         # RMSF plugin package
    │   ├── contacts/     # Contacts plugin package
    │   ├── distances/    # Distance analysis plugin package
    │   ├── secondary_structure/  # Secondary structure plugin package
    │   └── ...           # One sub-package per analysis type
    ├── config/
    │   └── comparison.py # ComparisonConfig, PlotSettings, condition models
    ├── analyses/shared/
    │   ├── inferential_statistics.py # t-tests, ANOVA, effect sizes
    │   ├── result_io.py  # Comparison result discovery/loading
    │   └── paths.py      # sanitize_label() and path helpers
    ├── cli/
    │   └── compare.py    # `polyzymd compare` subcommands
    └── cli/              # Command-line interface
        └── main.py       # Click CLI


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
- :py:class:`~polyzymd.analyses.distances.DistanceCalculator` - Inter-group distances
- :py:class:`~polyzymd.analyses.contacts.ParallelContactAnalyzer` - Polymer-protein contacts
- :py:class:`~polyzymd.analyses._results_base.BaseAnalysisResult` - Serializable result base class

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
