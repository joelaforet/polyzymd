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
