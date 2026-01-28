# Architecture Guide

This guide provides a detailed map of PolyzyMD's codebase, helping developers and power users understand where to find and modify specific functionality.

## Package Structure

```
src/polyzymd/
├── cli/                    # Command-line interface
│   └── main.py             # Click commands (build, run, submit, continue, validate)
├── config/                 # Configuration and validation
│   ├── schema.py           # Pydantic models for YAML config
│   └── loader.py           # YAML loading utilities
├── builders/               # System construction (Stage 1-4)
│   ├── system_builder.py   # Main orchestrator - coordinates all builders
│   ├── enzyme.py           # Protein/enzyme loading and preparation
│   ├── substrate.py        # Small molecule (ligand) handling
│   ├── polymer.py          # Random copolymer generation
│   └── solvent.py          # Solvation, ions, box setup
├── simulation/             # MD execution (Stage 5-6)
│   ├── runner.py           # Initial simulation (equilibration + production)
│   └── continuation.py     # Checkpoint-based continuation for daisy-chain
├── workflow/               # HPC job management
│   ├── slurm.py            # SLURM script templates and generation
│   └── daisy_chain.py      # Job submission with dependencies
├── core/                   # Shared utilities
│   ├── parameters.py       # Simulation parameter dataclasses
│   └── restraints.py       # Distance restraint definitions
└── __init__.py             # Package version and exports
```

## Data Flow Overview

```
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│   config.yaml   │────▶│  SimulationConfig│────▶│  SystemBuilder  │
└─────────────────┘     │  (schema.py)     │     │  (system_builder│
                        └─────────────────┘     │   .py)          │
                                                └────────┬────────┘
                                                         │
                        ┌────────────────────────────────┼────────────────────────────────┐
                        │                                │                                │
                        ▼                                ▼                                ▼
               ┌─────────────────┐             ┌─────────────────┐             ┌─────────────────┐
               │  EnzymeBuilder  │             │  SubstrateBuilder│            │  PolymerBuilder │
               │  (enzyme.py)    │             │  (substrate.py)  │            │  (polymer.py)   │
               └────────┬────────┘             └────────┬────────┘            └────────┬────────┘
                        │                                │                                │
                        └────────────────────────────────┼────────────────────────────────┘
                                                         │
                                                         ▼
                                                ┌─────────────────┐
                                                │  SolventBuilder │
                                                │  (solvent.py)   │
                                                └────────┬────────┘
                                                         │
                                                         ▼
                                                ┌─────────────────┐
                                                │   Interchange   │
                                                │   (OpenFF)      │
                                                └────────┬────────┘
                                                         │
                        ┌────────────────────────────────┼────────────────────────────────┐
                        │                                                                 │
                        ▼                                                                 ▼
               ┌─────────────────┐                                               ┌─────────────────┐
               │SimulationRunner │                                               │ContinuationMgr  │
               │  (runner.py)    │                                               │(continuation.py)│
               └─────────────────┘                                               └─────────────────┘
```

## Component Details

### CLI (`cli/main.py`)

The command-line interface is built with [Click](https://click.palletsprojects.com/). Each command maps to a specific workflow:

| Command | Function | Primary Module Used |
|---------|----------|---------------------|
| `polyzymd build` | Build system without running | `SystemBuilder` |
| `polyzymd run` | Build and run initial simulation | `SystemBuilder` + `SimulationRunner` |
| `polyzymd continue` | Continue from checkpoint | `ContinuationManager` |
| `polyzymd submit` | Submit daisy-chain to SLURM | `DaisyChainSubmitter` |
| `polyzymd validate` | Validate config file | `SimulationConfig` |

**Where to modify:**
- Add new CLI commands: `cli/main.py`
- Change CLI option defaults: Look for `@click.option` decorators

### Configuration (`config/schema.py`)

All configuration is defined as [Pydantic](https://docs.pydantic.dev/) models, providing automatic validation and type checking.

**Key classes:**
- `SimulationConfig` - Top-level container, has `from_yaml()` class method
- `EnzymeConfig` - Enzyme PDB path and settings
- `SubstrateConfig` - Ligand SDF path, charge method
- `PolymerConfig` - Polymer generation settings (monomers, length, count)
- `SolventConfig` - Water model, ion concentration, box shape
- `SimulationPhaseConfig` - Duration, timestep, ensemble settings
- `OutputConfig` - Directory paths for scratch/projects

**Where to modify:**
- Add new config options: Add fields to the appropriate Pydantic model
- Change validation rules: Add `@validator` decorators
- Change defaults: Modify field default values

### Builders (`builders/`)

The builders follow a pipeline pattern, each responsible for one component:

#### `system_builder.py` - Orchestrator

**Key methods:**
- `from_config(config)` - Create builder from SimulationConfig
- `build_from_config(config, working_dir, polymer_seed)` - Full pipeline
- `build_enzyme(pdb_path)` - Stage 1
- `build_substrate(sdf_path, ...)` - Stage 2
- `build_polymers(characters, probabilities, ...)` - Stage 3
- `pack_polymers(padding, working_dir)` - Uses PACKMOL
- `solvate(composition, padding, box_shape)` - Stage 4
- `create_interchange()` - Stage 5
- `get_openmm_components()` - Extract topology, system, positions

#### `enzyme.py` - Protein Loading

Uses OpenFF Toolkit to load PDB files with proper residue templates.

**Key method:** `build(pdb_path) -> Topology`

#### `substrate.py` - Ligand Handling

Loads docked conformers from SDF, assigns partial charges (NAGL by default).

**Key method:** `build(sdf_path, conformer_index, charge_method) -> Molecule`

#### `polymer.py` - Random Copolymer Generation

Generates random sequences based on monomer probabilities, loads pre-built SDFs.

**Key method:** `build(count, seed) -> Tuple[List[Molecule], List[int]]`

**Where to modify:**
- Change polymer loading: `_load_polymer_sdf()` method
- Change sequence generation: `_generate_sequences()` method

#### `solvent.py` - Solvation

Adds water (TIP3P), neutralizing ions, and optional co-solvents.

**Key method:** `solvate(topology, composition, padding, box_shape) -> Topology`

### Simulation (`simulation/`)

#### `runner.py` - Initial Runs

Handles energy minimization, equilibration (NVT), and production (NPT).

**Key methods:**
- `minimize(max_iterations, tolerance)` - Energy minimization
- `run_equilibration(temperature, duration_ns, ...)` - NVT equilibration
- `run_production(temperature, duration_ns, pressure, segment_index, ...)` - NPT production
- `run_from_config(config, segment_index)` - Run using config settings

**Reporters configured (lines 273-288, 398-413):**
```python
StateDataReporter(
    state_path,
    report_interval,
    step=True,
    time=True,
    potentialEnergy=True,
    kineticEnergy=True,
    totalEnergy=True,
    temperature=True,
    volume=True,
    density=True,
    speed=True,  # Reports performance in ns/day
)
```

**Where to modify:**
- Change reported quantities: Modify `StateDataReporter` arguments
- Change integrator settings: `_create_integrator()` method
- Change barostat settings: `_add_barostat()` method

#### `continuation.py` - Checkpoint Continuation

Loads state from previous segment checkpoint and continues simulation.

**Key methods:**
- `load_previous_state()` - Load checkpoint and system
- `run_segment(duration_ns, num_samples)` - Run next segment

### Workflow (`workflow/`)

#### `slurm.py` - SLURM Script Generation

Contains job script templates and preset configurations.

**Key components:**
- `SlurmConfig` - Dataclass with partition, QOS, time limit, etc.
- `SlurmConfig.from_preset(preset)` - Load preset (aa100, al40, testing, etc.)
- `INITIAL_JOB_TEMPLATE` - Template for first segment (builds system)
- `CONTINUATION_JOB_TEMPLATE` - Template for subsequent segments
- `SlurmJobGenerator` - Fills templates with runtime values

**Where to modify:**
- Add new SLURM presets: Add to `presets` dict in `from_preset()` method
- Change job script behavior: Modify `INITIAL_JOB_TEMPLATE` or `CONTINUATION_JOB_TEMPLATE`
- Change module loading: Edit the `module load` lines in templates

**Current template structure (simplified):**
```bash
#!/bin/bash
#SBATCH --partition={partition}
#SBATCH --qos={qos}
# ... other SBATCH directives

# Load conda environment (ignore module warnings)
module purge 2>/dev/null || true
module load miniforge 2>/dev/null || true
eval "$(conda shell.bash hook)"
mamba activate {conda_env}

# Enable strict error handling after environment setup
set -e

# Run simulation
polyzymd run -c "{config_path}" --replicate {replicate} ...
```

#### `daisy_chain.py` - Job Submission

Manages job dependencies for long simulations split across multiple SLURM jobs.

**Key components:**
- `DaisyChainConfig` - Settings for the chain (segments, time per segment)
- `DaisyChainSubmitter` - Generates scripts and submits with dependencies
- `submit_daisy_chain()` - Main entry point function

**Where to modify:**
- Change dependency type: Modify `--dependency=afterok:` in submission
- Change job naming: Modify `JobContext` creation

### Core Utilities (`core/`)

#### `parameters.py` - Simulation Parameters

Dataclasses for simulation settings, useful for programmatic use.

- `IntegratorParameters` - Timestep, friction, thermostat
- `EquilibrationParameters` - NVT settings
- `ProductionParameters` - NPT settings  
- `ReporterParameters` - Output frequency settings
- `SimulationParameters` - Combined container

#### `restraints.py` - Distance Restraints

Defines atom selections and harmonic distance restraints.

- `AtomSelection` - MDAnalysis-style selection with chain/resid/atom name
- `RestraintDefinition` - Two selections + distance + force constant
- `RestraintFactory` - Creates OpenMM forces from definitions

## Common Modifications

### Adding a New Reporter Output

**File:** `src/polyzymd/simulation/runner.py`

**Location:** Lines 275-288 (equilibration) and 399-413 (production)

Add parameters to `StateDataReporter`:
```python
StateDataReporter(
    ...,
    remainingTime=True,  # Add this for estimated time remaining
    elapsedTime=True,    # Add this for elapsed wall time
)
```

See [OpenMM StateDataReporter docs](http://docs.openmm.org/latest/api-python/generated/openmm.app.statedatareporter.StateDataReporter.html) for all options.

### Adding a New SLURM Preset

**File:** `src/polyzymd/workflow/slurm.py`

**Location:** `SlurmConfig.from_preset()` method (~line 50)

```python
presets = {
    "aa100": {...},
    "my_new_preset": {
        "partition": "my_partition",
        "qos": "my_qos",
        "time_limit": "12:00:00",
        "gpus": 2,
        ...
    },
}
```

### Changing the Force Field

**File:** `src/polyzymd/builders/system_builder.py`

**Location:** `__init__` method (line 54-58)

```python
def __init__(
    self,
    protein_forcefield: str = "ff14sb_off_impropers_0.0.4.offxml",  # Change this
    small_molecule_forcefield: str = "openff-2.0.0.offxml",         # Or this
) -> None:
```

Or in your config.yaml:
```yaml
force_field:
  protein: "ff19sb.offxml"
  small_molecule: "openff-2.1.0.offxml"
```

### Adding Custom Pre/Post-Processing Steps

**File:** `src/polyzymd/simulation/runner.py`

Add methods to `SimulationRunner` class, then call them from `cli/main.py`.

## Testing Your Changes

After modifying the code:

1. **Verify syntax:**
   ```bash
   python -m py_compile src/polyzymd/path/to/modified_file.py
   ```

2. **Reinstall package:**
   ```bash
   pip install -e .
   ```

3. **Test with dry-run:**
   ```bash
   polyzymd submit -c config.yaml --preset testing --dry-run
   ```

4. **Check generated scripts:**
   ```bash
   cat job_scripts/initial_seg0_rep1.sh
   ```

## Related Documentation

- [Configuration Guide](configuration.md) - YAML config options
- [HPC/SLURM Guide](hpc_slurm.md) - Job submission details
- [API Overview](../api/overview.rst) - Class and method reference
