# Dynamic Polymer Generation

This guide explains how to generate polymer chains on-the-fly from raw monomer SMILES strings, eliminating the need for pre-built polymer SDF files.

## Overview

PolyzyMD supports two modes for polymer generation:

| Mode | Description | Use Case |
|------|-------------|----------|
| **Cached** (default) | Load pre-built polymers from SDF files | Reproducibility, specific sequences |
| **Dynamic** | Generate polymers from monomer SMILES | Flexibility, new monomers, rapid prototyping |

Dynamic mode uses **ATRP (Atom-Transfer Radical Polymerization)** chemistry to:
1. Activate raw monomer SMILES via initiation reactions
2. Generate all possible polymer fragments
3. Build complete polymer chains with proper termination
4. Assign partial charges automatically

## When to Use Dynamic Mode

**Use Dynamic Mode when:**
- You want to test new monomer chemistries quickly
- You don't have pre-built polymer SDF files
- You want the system to handle fragment generation automatically

**Use Cached Mode when:**
- You need exact reproducibility of polymer structures
- You have validated, pre-built polymer SDFs
- You're running production simulations with known good structures

---

## Quick Start

Here's a minimal configuration for dynamic polymer generation:

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

# ... rest of configuration (solvent, thermodynamics, etc.)
```

Then run:

```bash
polyzymd run -c config.yaml -r 1
```

---

## Configuration Reference

### Generation Mode

```yaml
polymers:
  generation_mode: "dynamic"    # or "cached" (default)
```

| Value | Description |
|-------|-------------|
| `cached` | Load polymers from pre-built SDF files (default) |
| `dynamic` | Generate polymers from monomer SMILES using ATRP chemistry |

### Monomer SMILES

In dynamic mode, each monomer requires a `smiles` field:

```yaml
monomers:
  - label: "A"
    probability: 0.7
    name: "SBMA"
    smiles: "[H]C([H])=C(C(=O)..."    # Raw methacrylate SMILES
    residue_name: "SBM"               # Optional: 3-letter residue name
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
      volume_fraction: 0.30
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
  equilibration:
    ensemble: "NVT"
    duration: 1.0
    samples: 100
    time_step: 2.0
    thermostat: "LangevinMiddle"
    thermostat_timescale: 1.0
  production:
    ensemble: "NPT"
    duration: 100.0
    samples: 2500
    time_step: 2.0
    thermostat: "LangevinMiddle"
    thermostat_timescale: 1.0
    barostat: "MC"
    barostat_frequency: 25
  segments: 10

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

When you run a simulation with `generation_mode: "dynamic"`, the system:

1. **Loads raw monomer SMILES** from your configuration
2. **Runs initiation reactions** to activate the monomers (chlorination for ATRP)
3. **Runs polymerization reactions** to enumerate all possible fragments
4. **Runs termination reactions** on 1-site fragments to restore terminal alkenes
5. **Creates a MonomerGroup** with named fragments (e.g., `SBMA_2-site`, `SBMA_1-site`)

### Step 2: Polymer Building

For each polymer chain:

1. **Generate random sequence** based on monomer probabilities (e.g., "AABBA")
2. **Set terminal orientations** (head/tail use 1-site fragments)
3. **Build 3D structure** using Polymerist's `build_linear_polymer()`
4. **Validate no ring-piercing** (retry if necessary)
5. **Assign partial charges** using the configured charger

### Step 3: Caching

Generated fragments and polymers are cached for reuse:

```
.polymer_cache/
├── SBMA-EGPMA_monomer_group.json     # Fragment definitions
├── SBMA-EGPMA_AABBA_5-mer_charged.sdf  # Charged polymer structures
└── ...
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
    smiles: "[H]C([H])=C(C(=O)O...)C([H])([H])[H]"
```

```{note}
For non-methacrylate monomers or different polymerization chemistries (e.g., ring-opening, condensation), you would need to provide custom reaction files. This is an advanced use case.
```

---

## Troubleshooting

### "No 1-site terminal fragment found for monomer"

**Cause**: The cached MonomerGroup has old fragment names.

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

### "Symbol 'X' not present in monomer group mapping"

**Cause**: Mismatch between sequence labels and MonomerGroup fragments.

**Solution**: Ensure all monomer labels (A, B, C...) in your config have corresponding fragments generated. Delete the cache and regenerate.

### Slow fragment generation

**Cause**: First run generates all fragments, which can take time.

**Solution**: This is normal for the first run. Subsequent runs use the cached MonomerGroup and are much faster.

---

## Performance Tips

1. **Start small**: Test with short chains (length: 3-5) and few polymers (count: 1-2) first
2. **Use NAGL charger**: It's much faster than AM1BCC for polymer charging
3. **Reuse cache**: Keep `.polymer_cache` between runs for the same monomer set
4. **Parallelize later**: Dynamic generation currently runs sequentially; HPC parallelization is planned

---

## See Also

- {doc}`polymers` - General polymer configuration guide
- {doc}`configuration` - Complete configuration reference
- {doc}`gromacs_export` - Running with GROMACS
