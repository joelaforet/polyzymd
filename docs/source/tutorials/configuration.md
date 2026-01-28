# Configuration Reference

This document describes all configuration options for PolyzyMD YAML files.

## Configuration Structure

A complete configuration file has these sections:

```yaml
name: "simulation_name"
description: "optional description"

enzyme: { ... }           # Required
substrate: { ... }        # Optional (null for apo)
polymers: { ... }         # Optional (null to disable)
solvent: { ... }          # Required
restraints: [ ... ]       # Optional
thermodynamics: { ... }   # Required
simulation_phases: { ... } # Required
output: { ... }           # Required
force_field: { ... }      # Optional (has defaults)
```

---

## Enzyme Configuration

```yaml
enzyme:
  name: "LipA"                           # Identifier (required)
  pdb_path: "structures/enzyme.pdb"      # Path to PDB file (required)
  description: "Bacillus subtilis Lipase A"  # Optional description
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `name` | string | Yes | Short identifier for the enzyme |
| `pdb_path` | path | Yes | Path to prepared PDB file |
| `description` | string | No | Human-readable description |

---

## Substrate Configuration

```yaml
substrate:
  name: "Resorufin-Butyrate"             # Identifier (required)
  sdf_path: "structures/substrate.sdf"   # Path to SDF file (required)
  conformer_index: 0                     # Which conformer to use (default: 0)
  charge_method: "nagl"                  # Charge assignment method
  residue_name: "LIG"                    # 3-letter residue name
```

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `name` | string | Yes | - | Substrate identifier |
| `sdf_path` | path | Yes | - | Path to SDF with docked conformers |
| `conformer_index` | int | No | 0 | Index of conformer to use (0-indexed) |
| `charge_method` | string | No | "nagl" | Options: `nagl`, `espaloma`, `am1bcc` |
| `residue_name` | string | No | "LIG" | 3-letter code for topology |

### Charge Methods

| Method | Description | Speed |
|--------|-------------|-------|
| `nagl` | Graph neural network charges | Fast |
| `espaloma` | Machine learning charges | Medium |
| `am1bcc` | Semi-empirical QM charges | Slow |

### No Substrate (Apo Simulation)

```yaml
substrate: null
```

---

## Polymer Configuration

```yaml
polymers:
  enabled: true                          # Enable/disable polymers
  type_prefix: "SBMA-EGPMA"              # Polymer type identifier
  
  monomers:                              # Monomer definitions
    - label: "A"                         # Single character label
      probability: 0.98                  # Selection probability (0-1)
      name: "SBMA"                       # Full name (optional)
    - label: "B"
      probability: 0.02
      name: "EGPMA"
  
  length: 5                              # Monomers per chain
  count: 2                               # Number of polymer chains
  
  sdf_directory: null                    # Pre-built polymer SDFs (optional)
  cache_directory: ".polymer_cache"      # Cache for generated polymers
```

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `enabled` | bool | No | true | Enable polymer addition |
| `type_prefix` | string | Yes | - | Identifier for polymer type |
| `monomers` | list | Yes | - | Monomer specifications |
| `length` | int | Yes | - | Chain length (number of monomers) |
| `count` | int | Yes | - | Number of chains to add |
| `sdf_directory` | path | No | null | Directory with pre-built polymer SDFs |
| `cache_directory` | path | No | ".polymer_cache" | Cache directory |

### Monomer Specification

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `label` | string | Yes | Single character (A, B, C...) |
| `probability` | float | Yes | Selection probability (must sum to 1.0) |
| `name` | string | No | Full monomer name |

### No Polymers

```yaml
polymers: null
# or
polymers:
  enabled: false
```

---

## Solvent Configuration

```yaml
solvent:
  primary:
    type: "water"
    model: "tip3p"                       # Water model
  
  co_solvents: []                        # List of co-solvents (optional)
  
  ions:
    neutralize: true                     # Add counter-ions
    nacl_concentration: 0.15             # NaCl concentration (M)
  
  box:
    padding: 1.2                         # nm from solute to box edge
    shape: "rhombic_dodecahedron"        # Box shape
    target_density: 1.0                  # g/mL
    tolerance: 2.0                       # PACKMOL tolerance (Angstrom)
```

### Water Models

| Model | Description |
|-------|-------------|
| `tip3p` | TIP3P (default, fast) |
| `spce` | SPC/E |
| `tip4pew` | TIP4P-Ew |
| `opc` | OPC (accurate, slower) |

### Box Shapes

| Shape | Description |
|-------|-------------|
| `cube` | Cubic box |
| `rhombic_dodecahedron` | Space-efficient (default) |
| `truncated_octahedron` | Alternative space-efficient |

### Co-solvents

```yaml
co_solvents:
  - name: "dmso"
    smiles: "CS(=O)C"                    # SMILES string
    volume_fraction: 0.10                # 10% v/v
    residue_name: "DMS"                  # 3-letter code
```

---

## Restraints Configuration

```yaml
restraints:
  - type: "flat_bottom"                  # Restraint type
    name: "substrate_active_site"        # Identifier
    atom1:
      selection: "resid 77 and name OG"  # First atom selection
      description: "Catalytic serine"    # Optional description
    atom2:
      selection: "resname LIG and name C1"
      description: "Substrate carbon"
    distance: 3.3                        # Angstroms
    force_constant: 10000.0              # kJ/mol/nm²
    enabled: true                        # Enable/disable
```

See {doc}`restraints` for detailed selection syntax.

### Restraint Types

| Type | Description |
|------|-------------|
| `flat_bottom` | No force within threshold, harmonic beyond |
| `harmonic` | Harmonic potential at target distance |
| `upper_wall` | Prevent distance exceeding threshold |
| `lower_wall` | Prevent distance below threshold |

---

## Thermodynamics Configuration

```yaml
thermodynamics:
  temperature: 300.0                     # Kelvin
  pressure: 1.0                          # atmospheres
```

---

## Simulation Phases Configuration

```yaml
simulation_phases:
  equilibration:
    ensemble: "NVT"                      # NVT, NPT, or NVE
    duration: 1.0                        # nanoseconds
    samples: 100                         # frames to save
    time_step: 2.0                       # femtoseconds
    thermostat: "LangevinMiddle"
    thermostat_timescale: 1.0            # picoseconds
  
  production:
    ensemble: "NPT"
    duration: 100.0                      # nanoseconds total
    samples: 2500                        # total frames
    time_step: 2.0
    thermostat: "LangevinMiddle"
    thermostat_timescale: 1.0
    barostat: "MC"                       # Monte Carlo barostat
    barostat_frequency: 25               # steps between barostat moves
  
  segments: 10                           # Split production into segments
```

### Ensembles

| Ensemble | Description |
|----------|-------------|
| `NVT` | Constant volume, temperature |
| `NPT` | Constant pressure, temperature |
| `NVE` | Microcanonical (no thermostat) |

### Thermostats

| Thermostat | Description |
|------------|-------------|
| `LangevinMiddle` | Langevin integrator (recommended) |
| `Langevin` | Standard Langevin |
| `Andersen` | Andersen thermostat |
| `NoseHoover` | Nosé-Hoover chain |

### Barostats

| Barostat | Description |
|----------|-------------|
| `MC` | Monte Carlo barostat (recommended) |
| `MCA` | Monte Carlo anisotropic |

---

## Output Configuration

Environment variables (`$USER`, `$HOME`, `${VAR}`) and `~` are automatically expanded in path fields.

```yaml
output:
  # Directory structure - environment variables are expanded automatically
  projects_directory: "/projects/$USER/polyzymd"   # Scripts, logs
  scratch_directory: "/scratch/alpine/$USER/simulations"  # Trajectories
  
  # You can also use ~ for home directory
  # projects_directory: "~/polyzymd"
  
  # Subdirectories within projects_directory
  job_scripts_subdir: "job_scripts"
  slurm_logs_subdir: "slurm_logs"
  
  # Naming
  naming_template: "{enzyme}_{substrate}_{polymer_type}_{temperature}K_run{replicate}"
  
  # Output options
  save_checkpoint: true                  # Save restart files
  save_state_data: true                  # Save energy/temperature CSV
  trajectory_format: "dcd"               # dcd or xtc
```

### Naming Template Variables

| Variable | Description | Example |
|----------|-------------|---------|
| `{enzyme}` | Enzyme name | "LipA" |
| `{substrate}` | Substrate name | "ResorufinButyrate" |
| `{polymer_type}` | Polymer type | "SBMA-EGPMA" |
| `{temperature}` | Temperature in K | "300" |
| `{replicate}` | Replicate number | "1" |

---

## Force Field Configuration

```yaml
force_field:
  protein: "ff14sb_off_impropers_0.0.4.offxml"  # Protein force field
  small_molecule: "openff-2.0.0.offxml"          # Ligand/polymer force field
```

### Available Force Fields

**Protein:**
- `ff14sb_off_impropers_0.0.4.offxml` - Amber ff14SB (recommended)

**Small Molecule:**
- `openff-2.0.0.offxml` - OpenFF Sage 2.0 (recommended)
- `openff-2.1.0.offxml` - OpenFF Sage 2.1

---

## Complete Example

See the example configurations in `src/polyzymd/configs/examples/`:

- `enzyme_only.yaml` - Enzyme + substrate, no polymers
- `enzyme_polymer.yaml` - Full enzyme + polymer simulation
- `enzyme_cosolvent.yaml` - Enzyme with DMSO co-solvent
