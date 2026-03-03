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
    │   └── continuation.py    # Checkpoint continuation
    ├── workflow/         # HPC workflow management
    │   ├── slurm.py      # SLURM script generation
    │   └── daisy_chain.py     # Job chaining
    ├── core/             # Core utilities
    │   ├── parameters.py # Simulation parameters
    │   └── restraints.py # Restraint definitions
    ├── analysis/         # Post-simulation trajectory analysis
    │   ├── rmsf/         # Root Mean Square Fluctuation
    │   ├── distances/    # Inter-atomic distance analysis
    │   ├── triad/        # Catalytic triad integrity
    │   ├── contacts/     # Polymer-protein contacts, binding preference
    │   ├── sasa/         # Solvent Accessible Surface Area
    │   ├── exposure/     # Chaperone-like exposure dynamics
    │   ├── core/         # Statistics, autocorrelation, MetricType
    │   └── results/      # Serializable result models
    ├── compare/          # Multi-condition statistical comparison
    │   ├── comparators/  # RMSF, contacts, triad, BFE, affinity, exposure
    │   ├── results/      # Comparison result models
    │   └── plotters/     # Publication-quality figures
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

Workflow
~~~~~~~~

- :py:class:`~polyzymd.workflow.daisy_chain.DaisyChainSubmitter` - SLURM job submission
- :py:class:`~polyzymd.workflow.slurm.SlurmConfig` - SLURM configuration

Restraints
~~~~~~~~~~

- :py:class:`~polyzymd.core.restraints.RestraintDefinition` - Restraint specification
- :py:class:`~polyzymd.core.restraints.AtomSelection` - Atom selection
- :py:class:`~polyzymd.core.restraints.RestraintFactory` - Create restraints from config

Analysis
~~~~~~~~

- :py:class:`~polyzymd.analysis.rmsf.calculator.RMSFCalculator` - Per-residue RMSF
- :py:class:`~polyzymd.analysis.distances.calculator.DistanceCalculator` - Inter-group distances
- :py:class:`~polyzymd.analysis.triad.analyzer.CatalyticTriadAnalyzer` - Catalytic triad integrity
- :py:class:`~polyzymd.analysis.contacts.calculator_parallel.ParallelContactAnalyzer` - Polymer-protein contacts
- :py:class:`~polyzymd.analysis.results.base.BaseAnalysisResult` - Serializable result base class

Comparison
~~~~~~~~~~

- :py:class:`~polyzymd.compare.core.base.BaseComparator` - Template Method comparator base
- :py:class:`~polyzymd.compare.comparators.exposure.ExposureDynamicsComparator` - Chaperone analysis
- :py:class:`~polyzymd.compare.comparators.binding_free_energy.BindingFreeEnergyComparator` - Per-contact ΔG_sel
- :py:class:`~polyzymd.compare.comparators.polymer_affinity.PolymerAffinityScoreComparator` - Total interaction strength


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
