# Dynamic Polymer Generation

This guide explains how to generate linear methacrylate polymer chains
on-the-fly from raw monomer SMILES strings.

:::{admonition} Environment Setup
:class: tip

All commands below assume you have activated the PolyzyMD pixi environment:

```bash
pixi shell -e build
```

Alternatively, prefix each command with `pixi run -e build`.
:::

## Overview

PolyzyMD supports generated polymers and explicit user-supplied molecule pools:

| Mode | Description | Use Case |
|------|-------------|----------|
| **Dynamic** | Generate linear methacrylate polymers from monomer SMILES | Flexibility, new monomers, rapid prototyping |
| **Fragments** | Assemble native linear explicit fragments from terminal/middle strings | Linear recipes where you provide fragment strings |
| **Provided molecules** | Pack explicit charged SDF molecules through `provided_molecules` | Branched, dendrimer, nonlinear, or externally authored molecules |
| **Cached** | Deprecated sequence-derived SDF filename libraries | Compatibility with old inventories only |

Dynamic mode with the bundled `"default"` reactions uses PolyzyMD's native
linear methacrylate backend to build, charge, and cache OpenFF molecules.

## When to Use Dynamic Mode

**Use Dynamic Mode when:**
- You want to test new monomer chemistries quickly
- You want native linear methacrylate generation from monomer SMILES
- You want PolyzyMD to charge and cache generated polymers automatically

**Use `provided_molecules` when:**
- You already have charged single-molecule SDFs
- The molecule is branched, dendrimeric, nonlinear, or externally authored
- You do not want PolyzyMD to regenerate or autocharge the molecule

---

## Quick Start

Configuration snippet for dynamic polymer generation; add the remaining required
simulation sections from your project template:

```yaml
name: "dynamic_polymer_test"

enzyme:
  name: "MyEnzyme"
  pdb_path: "structures/enzyme.pdb"

polymers:
  enabled: true
  generation_mode: "dynamic"      # Enable dynamic generation
  type_prefix: "SBMA-EGPMA"
  
  monomers:
    - label: "A"
      probability: 0.7
      name: "SBMA"
      smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])[N+](C([H])([H])[H])(C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])S(=O)(=O)[O-])C([H])([H])[H]"
    - label: "B"
      probability: 0.3
      name: "EGPMA"
      smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])Oc1c([H])c([H])c([H])c([H])c1[H])C([H])([H])[H]"
  
  length: 5
  count: 2
  charger: "nagl"
  reactions:
    initiation: "default"
    polymerization: "default"
    termination: "default"

# ... rest of configuration (solvent, thermodynamics, etc.)
```

Then build and run:

```bash
polyzymd build -c config.yaml
polyzymd run-segment -c config.yaml -r 1
```

---

## Configuration Reference

### Generation Mode

Snippet:

```yaml
polymers:
  generation_mode: "dynamic"
```

| Value | Description |
|-------|-------------|
| `dynamic` | Generate linear methacrylate polymers from monomer SMILES with bundled default reactions |
| `fragments` | Assemble native linear explicit fragments from `terminal` and `middle` strings |
| `cached` | Deprecated sequence-derived SDF library compatibility mode |

### Monomer SMILES

In dynamic mode, each monomer requires a `smiles` field:

Snippet:

```yaml
monomers:
  - label: "A"
    probability: 0.7
    name: "MMA"
    smiles: "C=C(C)C(=O)OC"    # Raw methacrylate SMILES
    residue_name: "MMA"        # Optional: 3-letter residue name
```

| Field | Required | Description |
|-------|----------|-------------|
| `label` | Yes | Single character identifier (A, B, C...) |
| `probability` | Yes | Selection probability (must sum to 1.0) |
| `name` | Yes | Monomer name for identification |
| `smiles` | Yes (dynamic) | Raw monomer SMILES string |
| `residue_name` | No | 3-letter residue name (auto-generated if not provided) |

```{important}
**SMILES Format**: Provide the raw, unactivated monomer SMILES. For methacrylates, this means the structure with the C=C double bond intact. The system will handle activation (chlorination) automatically.
```

### ATRP Reaction Configuration

By default, PolyzyMD uses bundled ATRP reaction templates. You can also specify custom reaction files:

Snippet:

```yaml
polymers:
  generation_mode: "dynamic"
  
  reactions:
    initiation: "default"        # Use bundled template
    polymerization: "default"    # Use bundled template
    termination: "default"       # Use bundled template
  
  # Or specify custom reaction files:
  # reactions:
  #   initiation: "/path/to/my_initiation.rxn"
  #   polymerization: "/path/to/my_polymerization.rxn"
  #   termination: "/path/to/my_termination.rxn"
```

### Charge Assignment

Snippet:

```yaml
polymers:
  charger: "nagl"    # Charge method for polymer atoms
```

| Method | Description | Speed |
|--------|-------------|-------|
| `nagl` | Graph neural network charges (recommended) | Fast |
| `espaloma` | Machine learning charges | Medium |
| `am1bcc` | Semi-empirical QM charges | Slow |

### Retry Configuration

Dynamic generation may occasionally produce structures with ring-piercing artifacts. The system automatically retries:

Snippet:

```yaml
polymers:
  max_retries: 10    # Maximum attempts before failing (default: 10)
```

---

## Complete Example Configuration

Here's a complete configuration for running a simulation with dynamically generated polymers:

```yaml
name: "LipA_dynamic_polymer_simulation"
description: "Lipase A with dynamically generated SBMA-EGPMA copolymers"

# Enzyme
enzyme:
  name: "LipA"
  pdb_path: "structures/enzyme.pdb"

# Substrate (optional)
substrate:
  name: "ResorufinButyrate"
  sdf_path: "structures/substrate.sdf"
  charge_method: "nagl"
  residue_name: "RBY"

# Dynamic Polymer Generation
polymers:
  enabled: true
  generation_mode: "dynamic"
  type_prefix: "SBMA-EGPMA"
  
  # ATRP reactions (use bundled defaults)
  reactions:
    initiation: "default"
    polymerization: "default"
    termination: "default"
  
  # Monomer definitions with SMILES
  monomers:
    - label: "A"
      probability: 0.7
      name: "SBMA"
      smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])[N+](C([H])([H])[H])(C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])S(=O)(=O)[O-])C([H])([H])[H]"
      residue_name: "SBM"
    - label: "B"
      probability: 0.3
      name: "EGPMA"
      smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])Oc1c([H])c([H])c([H])c([H])c1[H])C([H])([H])[H]"
      residue_name: "EGM"
  
  # Chain parameters
  length: 5
  count: 2
  
  # Charge assignment
  charger: "nagl"
  max_retries: 10
  
  # Caching
  cache_directory: ".polymer_cache"

# Solvent
solvent:
  primary:
    type: "water"
    model: "tip3p"
  co_solvents:
    - name: "dmso"
      mole_fraction: 0.10
      residue_name: "DMS"
  ions:
    neutralize: true
    nacl_concentration: 0.0
  box:
    padding: 1.2
    shape: "rhombic_dodecahedron"
    target_density: 1.05
    tolerance: 2.0

# Thermodynamics
thermodynamics:
  temperature: 300.0
  pressure: 1.0

# Simulation phases
simulation_phases:
  equilibration_stages:
    - name: "heating"
      duration: 0.2
      samples: 20
      ensemble: "NVT"
      temperature_start: 60.0
      temperature_end: 300.0
      temperature_increment: 1.0
      temperature_interval: 1200.0
      position_restraints:
        - group: "protein_heavy"
          force_constant: 4184.0
        - group: "polymer_heavy"
          force_constant: 4184.0
    - name: "free_equilibration"
      duration: 0.8
      samples: 80
      ensemble: "NPT"
      temperature: 300.0
  production:
    ensemble: "NPT"
    duration: 100.0
    samples: 2500
    time_step: 2.0
    thermostat: "LangevinMiddle"
    thermostat_timescale: 1.0
    barostat: "MC"
    barostat_frequency: 25
# Output
output:
  projects_directory: "."
  scratch_directory: null
  job_scripts_subdir: "job_scripts"
  slurm_logs_subdir: "slurm_logs"
  naming_template: "{enzyme}_{polymer_type}_dynamic_run{replicate}"
  save_checkpoint: true
  save_state_data: true
  trajectory_format: "dcd"
```

---

## How Dynamic Generation Works

### Step 1: Fragment Generation

When you run a simulation with `generation_mode: "dynamic"` and bundled
`"default"` reactions, the native backend:

1. **Loads raw methacrylate monomer SMILES** from your configuration
2. **Generates a linear native recipe** for the requested sequence
3. **Builds and validates OpenFF geometry**
4. **Assigns partial charges** using the configured charger
5. **Writes charged cache artifacts** under `cache_directory`

### Step 2: Polymer Building

For each polymer chain:

1. **Generate random sequence** based on monomer probabilities (e.g., "AABBA")
2. **Build a native linear methacrylate structure** for that exact sequence
3. **Validate geometry**
4. **Assign partial charges** using the configured charger

### Step 3: Caching

Generated fragments and polymers are cached for reuse:

```
.polymer_cache/
├── native-methacrylate-v*/             # Native dynamic artifacts
│   └── ..._charged.sdf                 # Charged polymer structures
└── native-fragments-v*/                # Native fragments artifacts, when used
    └── ..._charged.sdf
```

To regenerate from scratch, delete the cache:

```bash
rm -rf .polymer_cache
```

---

## Supported Chemistries

### ATRP (Atom-Transfer Radical Polymerization)

Currently, dynamic generation supports **methacrylate-based monomers** using ATRP chemistry:

- Sulfobetaine methacrylate (SBMA)
- Ethylene glycol phenyl ether methacrylate (EGPMA)
- Trimethylammonium ethyl methacrylate (TMAEMA)
- Sulfopropyl methacrylate (SPMA)
- Oligo(ethylene glycol) methacrylate (OEGMA)

### Adding New Monomers

To use a new methacrylate monomer:

1. Obtain the SMILES string (with C=C double bond intact)
2. Add it to your configuration:

```yaml
monomers:
  - label: "C"
    probability: 0.1
    name: "MyNewMonomer"
    smiles: "C=C(C)C(=O)OCCO"
```

```{note}
For non-methacrylate monomers or different polymerization chemistries, provide
explicit charged SDF molecules through `provided_molecules` or use the advanced
native fragments mode. Runtime CGSmiles authoring and a fragment-design notebook
are future work, not current user interfaces.
```

---

## Troubleshooting

### Cached artifact mismatch or geometry failure

**Cause**: Cached native artifacts were produced from an older recipe or the
current generated geometry failed validation.

**Solution**: Delete the cache and regenerate:
```bash
rm -rf .polymer_cache
```

### "Failed to build polymer after N attempts due to ring-piercing"

**Cause**: The polymer structure has atoms passing through rings.

**Solutions**:
1. Increase `max_retries` (default is 10)
2. Try shorter polymer chains
3. Check that monomer SMILES are correct

### "Symbol 'X' not present in monomer mapping"

**Cause**: Mismatch between generated sequence labels and configured monomers.

**Solution**: Ensure every sequence label (A, B, C...) appears in `monomers`.
Delete the cache and regenerate if the configuration changed.

### Slow first generation

**Cause**: First run builds, charges, validates, and writes cache artifacts.

**Solution**: This is normal for the first run. Subsequent runs reuse charged
cache artifacts and are much faster.

---

## Performance Tips

1. **Start small**: Test with short chains (length: 3-5) and few polymers (count: 1-2) first
2. **Use NAGL charger**: It's much faster than AM1BCC for polymer charging
3. **Reuse cache**: Keep `.polymer_cache` between runs for the same monomer set
4. **Parallelize later**: Dynamic generation currently runs sequentially; HPC parallelization is planned

---

## See Also

- {doc}`polymers` - General polymer configuration guide
- {doc}`../reference/configuration` - Complete configuration reference
- {doc}`gromacs_export` - Running with GROMACS
