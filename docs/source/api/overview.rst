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
    ├── analyses/         # ★ Study API, analysis functions and `polyzymd analyze`
    │   ├── study.py      # Study API: replicates as MDAnalysis universes
    │   ├── timeseries.py # study.timeseries / per_replicate, summaries and tests
    │   ├── functions.py  # Shipped analysis functions (rg, rmsd, rmsf, sasa, contacts,
    │   │                 # hydrogen bonds, distances, ...)
    │   ├── figures.py    # Figures drawn from stored values
    │   ├── protocols.py  # `polyzymd analyze` and the ProtocolReport
    │   ├── mda/          # MDAnalysis loading, file identity, job and artifact layer
    │   ├── shared/       # Reusable utilities (TrajectoryLoader, alignment, etc.)
    │   ├── base.py, discovery.py, orchestrator.py, stats.py, _framework/
    │   │                 # Analysis plugin framework: no shipped analysis is a
    │   │                 # plugin any more, and it is being removed
    └── cli/              # Command-line interface
        ├── analyze.py    # `polyzymd analyze`
        ├── compare.py    # `polyzymd compare` subcommands
        └── main.py       # Click CLI

No shipped analysis is a plugin: every analysis that ``polyzymd analyze``
offers is a function in ``polyzymd.analyses.functions`` run through
``polyzymd.analyses.study.Study``. The ``analyses/_framework/`` package is
private implementation detail of the plugin framework.


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

- :py:class:`~polyzymd.analyses.study.Study` - Every replicate of every condition as MDAnalysis universes
- :py:func:`~polyzymd.analyses.functions.hydrogen_bonds` - Hydrogen-bond counts between two groups
- :py:func:`~polyzymd.analyses.functions.residue_occlusion` - Polymer-protein contacts by occluded SASA
- :py:class:`~polyzymd.analyses.base.Analysis` - Base class of the plugin framework, which no shipped analysis uses
- :py:class:`~polyzymd.analyses.mda.job.MDAAnalysisJob` - Unit of trajectory-native MDAnalysis work
- :py:class:`~polyzymd.analyses.mda.frame_selection.FrameSelection` - Frame selection and equilibration-window descriptor
- :py:class:`~polyzymd.analyses.mda.artifacts.ReplicateArtifact` - Validated output for one replicate
- :py:class:`~polyzymd.analyses.mda.artifacts.ConditionArtifact` - Aggregated output for one condition
- :py:class:`~polyzymd.analyses.mda.artifacts.ComparisonArtifact` - Canonical comparison output for cross-condition results
- :py:class:`~polyzymd.analyses.mda.store.ArtifactStore` - Canonical artifact persistence and loading helper
- :py:class:`~polyzymd.analyses.base.MetricValue` - Scalar metric descriptor for default comparisons

Comparison
~~~~~~~~~~

- :py:class:`~polyzymd.analyses.protocols.ProtocolReport` - The report of ``polyzymd analyze`` and of ``summary()`` and ``compare()``
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
