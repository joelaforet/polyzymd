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

PolyzyMD supports two modes for polymer generation: **cached** (load from pre-built SDF files) and **dynamic** (generate on-the-fly from SMILES).

```{tip}
For a complete guide on dynamic polymer generation, see {doc}`dynamic_polymers`.
```

### Basic Configuration (Cached Mode)

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

### Dynamic Generation Configuration

To generate polymers on-the-fly from monomer SMILES (without pre-built SDF files):

```yaml
polymers:
  enabled: true
  generation_mode: "dynamic"             # Enable dynamic generation
  type_prefix: "SBMA-EGPMA"
  
  # ATRP reaction templates (use bundled defaults or custom paths)
  reactions:
    initiation: "default"                # or "/path/to/custom.rxn"
    polymerization: "default"
    termination: "default"
  
  monomers:
    - label: "A"
      probability: 0.7
      name: "SBMA"
      smiles: "[H]C([H])=C(C(=O)OC...)..."  # Required for dynamic mode
      residue_name: "SBM"                   # Optional 3-letter residue name
    - label: "B"
      probability: 0.3
      name: "EGPMA"
      smiles: "[H]C([H])=C(C(=O)OC...)..."
      residue_name: "EGM"
  
  length: 5
  count: 2
  charger: "nagl"                        # Charge method: nagl, espaloma, am1bcc
  max_retries: 10                        # Retries for ring-piercing detection
  cache_directory: ".polymer_cache"
```

### All Polymer Options

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `enabled` | bool | No | true | Enable polymer addition |
| `generation_mode` | string | No | "cached" | `cached` or `dynamic` |
| `type_prefix` | string | Yes | - | Identifier for polymer type |
| `monomers` | list | Yes | - | Monomer specifications |
| `length` | int | Yes | - | Chain length (number of monomers) |
| `count` | int | Yes | - | Number of chains to add |
| `sdf_directory` | path | No | null | Directory with pre-built polymer SDFs |
| `cache_directory` | path | No | ".polymer_cache" | Cache directory |
| `reactions` | object | No | all "default" | ATRP reaction templates (dynamic mode) |
| `charger` | string | No | "nagl" | Charge method for dynamic generation |
| `max_retries` | int | No | 10 | Max attempts for ring-piercing avoidance |

### Monomer Specification

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `label` | string | Yes | Single character (A, B, C...) |
| `probability` | float | Yes | Selection probability (must sum to 1.0) |
| `name` | string | No | Full monomer name |
| `smiles` | string | Dynamic only | Raw monomer SMILES (with C=C double bond) |
| `residue_name` | string | No | 3-letter residue code for topology |

### Charge Methods for Dynamic Generation

| Method | Description | Speed | Accuracy |
|--------|-------------|-------|----------|
| `nagl` | Graph neural network charges | Fast | Good |
| `espaloma` | Machine learning charges | Medium | Good |
| `am1bcc` | Semi-empirical QM charges | Slow | Best |

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

PolyzyMD supports adding co-solvents to your simulation system. You can specify co-solvents using either **volume fraction** (v/v) or **molar concentration**.

#### Specification Methods

| Method | Field | Description | Effect on Water |
|--------|-------|-------------|-----------------|
| Volume Fraction | `volume_fraction` | Fraction of box volume (0-1) | Reduces water proportionally |
| Concentration | `concentration` | Molar concentration (mol/L) | Additive (water unchanged) |

**Important:** Use exactly ONE method per co-solvent. Do not specify both `volume_fraction` and `concentration` for the same co-solvent.

#### Volume Fraction Method

Use this when you want a specific percentage of the solvent to be the co-solvent (e.g., "30% DMSO").

```yaml
co_solvents:
  - name: "dmso"
    volume_fraction: 0.30    # 30% v/v DMSO
```

**Formula:**

```
n = (V_box × phi × rho) / M

Where:
  n     = number of co-solvent molecules
  V_box = simulation box volume (L)
  phi   = volume fraction (e.g., 0.30 for 30%)
  rho   = co-solvent density (g/mL)
  M     = molar mass (g/mol)
```

**Source:** [`src/polyzymd/builders/solvent.py:267-287`](https://github.com/joelaforet/polyzymd/blob/main/src/polyzymd/builders/solvent.py#L267-L287)

The water count is reduced proportionally: if you specify 30% DMSO, water fills the remaining 70% of the box.

#### Concentration Method

Use this when you want a specific molar concentration (e.g., "2 M urea for protein denaturation studies").

```yaml
co_solvents:
  - name: "urea"
    concentration: 2.0       # 2 M urea
```

**Formula:**

```
n = C × V_box × N_A

Where:
  n     = number of co-solvent molecules
  C     = concentration (mol/L)
  V_box = simulation box volume (L)
  N_A   = Avogadro's number (implicit in OpenMM)
```

**Source:** [`src/polyzymd/builders/solvent.py:295-312`](https://github.com/joelaforet/polyzymd/blob/main/src/polyzymd/builders/solvent.py#L295-L312)

The water count is NOT reduced when using concentration. The co-solvent molecules are added to the existing water, which may slightly increase the effective density.

#### Built-in Co-solvent Library

PolyzyMD includes a library of common co-solvents with pre-defined SMILES and densities. Density values are sourced from [PubChem](https://pubchem.ncbi.nlm.nih.gov/), a public database of chemical compounds. Each compound has a unique Compound Identification Number (CID) that can be used to look up detailed information including density, structure, and safety data.

| Name | SMILES | Density (g/mL) | Reference |
|------|--------|----------------|-----------|
| `dmso` | `CS(=O)C` | 1.10 | [CID 679](https://pubchem.ncbi.nlm.nih.gov/compound/679) |
| `dmf` | `CN(C)C=O` | 0.95 | [CID 6228](https://pubchem.ncbi.nlm.nih.gov/compound/6228) |
| `acetonitrile` | `CC#N` | 0.786 | [CID 6342](https://pubchem.ncbi.nlm.nih.gov/compound/6342) |
| `urea` | `C(=O)(N)N` | 1.32 | [CID 1176](https://pubchem.ncbi.nlm.nih.gov/compound/1176) |
| `ethanol` | `CCO` | 0.789 | [CID 702](https://pubchem.ncbi.nlm.nih.gov/compound/702) |
| `methanol` | `CO` | 0.792 | [CID 887](https://pubchem.ncbi.nlm.nih.gov/compound/887) |
| `glycerol` | `C(C(CO)O)O` | 1.261 | [CID 753](https://pubchem.ncbi.nlm.nih.gov/compound/753) |
| `isopropanol` | `CC(C)O` | 0.786 | [CID 3776](https://pubchem.ncbi.nlm.nih.gov/compound/3776) |
| `acetone` | `CC(=O)C` | 0.784 | [CID 180](https://pubchem.ncbi.nlm.nih.gov/compound/180) |
| `thf` | `C1CCOC1` | 0.883 | [CID 8028](https://pubchem.ncbi.nlm.nih.gov/compound/8028) |
| `dioxane` | `C1COCCO1` | 1.033 | [CID 31275](https://pubchem.ncbi.nlm.nih.gov/compound/31275) |
| `ethylene_glycol` | `C(CO)O` | 1.114 | [CID 174](https://pubchem.ncbi.nlm.nih.gov/compound/174) |

For library co-solvents, you only need to specify the `name` and either `volume_fraction` or `concentration`:

```yaml
co_solvents:
  - name: "dmso"
    volume_fraction: 0.10    # 10% v/v DMSO - smiles and density auto-populated
```

#### Custom Co-solvents

For molecules not in the library, you must provide the SMILES string. Density is required only when using `volume_fraction`:

```yaml
co_solvents:
  # Custom co-solvent with volume fraction (density required)
  - name: "ethyl_acetate"
    smiles: "CCOC(=O)C"
    density: 0.902           # g/mL - required for volume_fraction
    volume_fraction: 0.15

  # Custom co-solvent with concentration (density not needed)
  - name: "my_additive"
    smiles: "CC(=O)NC"
    concentration: 0.5       # 0.5 M
```

#### Multiple Co-solvents

You can combine multiple co-solvents. Each can use either specification method independently:

```yaml
co_solvents:
  - name: "dmso"
    volume_fraction: 0.20    # 20% v/v DMSO
  - name: "urea"
    concentration: 1.0       # Plus 1 M urea
```

**Warning:** When using multiple co-solvents with `volume_fraction`, ensure the total does not exceed 1.0 (100%). The remaining fraction is filled with water.

```{warning}
**YAML List Syntax**

A common mistake is placing each field on a separate line with its own `-`, which creates multiple list items instead of one object with multiple fields.

**Incorrect** (creates 3 separate incomplete items):
~~~yaml
co_solvents:
  - name: "dmso"
  - volume_fraction: 0.30
  - residue_name: "DMS"
~~~

**Correct** (one item with 3 fields):
~~~yaml
co_solvents:
  - name: "dmso"
    volume_fraction: 0.30
    residue_name: "DMS"
~~~

The `-` character starts a **new list item**. All fields belonging to the same item must be indented to the same level *without* a leading `-`.
```

#### Assumptions and Limitations

- **Ideal mixing:** Volume fractions assume ideal mixing (volumes are additive). Real solutions may deviate.
- **Room temperature densities:** Library densities are approximate values at ~25C.
- **PACKMOL placement:** Co-solvent molecules are placed randomly by PACKMOL and may require equilibration to achieve uniform distribution.

### Solvent Parameterization

PolyzyMD uses **pre-computed partial charges** for all solvent molecules to ensure consistency and performance.

#### Why Pre-computed Charges?

When adding many copies of the same solvent molecule (e.g., 1000 DMSO molecules), each molecule should have **identical partial charges**. However, charge calculation methods like AM1BCC have numerical variability - running the calculation twice on the same molecule can produce slightly different charges.

If charges were computed independently for each solvent molecule:
1. **Inconsistency**: Identical molecules would have different parameters (physically incorrect)
2. **Performance**: AM1BCC is expensive; computing it 1000x is wasteful
3. **Force field issues**: Parameter variability can cause OpenFF Interchange errors

#### How It Works

PolyzyMD solves this by computing charges **once** and reusing them:

1. **Built-in solvents**: Pre-computed SDF files are bundled with the package (in `src/polyzymd/data/solvents/`)
2. **User cache**: Custom solvents are cached in `~/.polyzymd/solvent_cache/` after first use
3. **Lookup order**: Memory cache → Bundled SDFs → User cache → Generate and cache

```
# Lookup order for get_solvent_molecule("dmso")
1. Check in-memory cache (fastest)
2. Check bundled library: src/polyzymd/data/solvents/dmso.sdf
3. Check user cache: ~/.polyzymd/solvent_cache/dmso.sdf
4. Generate from SMILES + AM1BCC, save to user cache
```

#### Available Pre-computed Solvents

All 12 library co-solvents plus water models have pre-computed charges:

| Solvent | File | Charge Method |
|---------|------|---------------|
| TIP3P Water | `tip3p.sdf` | Literature values |
| DMSO | `dmso.sdf` | AM1BCC |
| DMF | `dmf.sdf` | AM1BCC |
| Acetonitrile | `acetonitrile.sdf` | AM1BCC |
| Urea | `urea.sdf` | AM1BCC |
| Ethanol | `ethanol.sdf` | AM1BCC |
| Methanol | `methanol.sdf` | AM1BCC |
| Glycerol | `glycerol.sdf` | AM1BCC |
| Isopropanol | `isopropanol.sdf` | AM1BCC |
| Acetone | `acetone.sdf` | AM1BCC |
| THF | `thf.sdf` | AM1BCC |
| Dioxane | `dioxane.sdf` | AM1BCC |
| Ethylene Glycol | `ethylene_glycol.sdf` | AM1BCC |

#### Custom Solvents

When you use a custom co-solvent (not in the library), PolyzyMD will:

1. Generate the molecule from your SMILES string
2. Compute AM1BCC partial charges (this may take a few seconds)
3. Cache the parameterized molecule to `~/.polyzymd/solvent_cache/`
4. Reuse the cached version for all future simulations

```yaml
co_solvents:
  - name: "my_custom_solvent"
    smiles: "CC(=O)OCC"        # First use: computes and caches charges
    concentration: 0.5         # Future uses: loads from cache instantly
```

#### Managing the Cache

You can inspect and manage the solvent cache programmatically:

```python
from polyzymd.data import list_available_solvents, clear_cache

# List all available solvents (bundled + cached)
solvents = list_available_solvents()
print(solvents)
# {'bundled': ['dmso', 'ethanol', ...], 'cached': ['my_custom_solvent']}

# Clear the user cache (does not affect bundled solvents)
clear_cache()
```

The user cache location is `~/.polyzymd/solvent_cache/`. You can safely delete this directory to force re-computation of custom solvents.

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

### Key Collision Warnings

When building systems with both proteins and small molecules, you may see warnings like:

```
Key collision with different parameters, fixing. Key is [#6X4:1]-[#1:2]
```

**This is expected behavior and does not indicate a problem.**

#### Why This Happens

PolyzyMD uses different force fields for different molecule types:
- **Proteins**: ff14SB (Amber force field ported to OpenFF format)
- **Small molecules**: OpenFF Sage 2.0 (general small molecule force field)

When these force fields are combined, the same SMIRKS pattern (e.g., `[#6X4:1]-[#1:2]` for sp³ carbon-hydrogen bonds) may appear in both, but with **different parameter values**. This is expected because:

1. ff14SB was optimized for protein behavior
2. OpenFF Sage was optimized for general organic molecules
3. Both are valid parameterizations for their respective domains

#### How OpenFF Handles This

OpenFF Interchange detects these collisions and resolves them by appending `_DUPLICATE` to the key, allowing both parameter sets to coexist:

```python
# Simplified OpenFF behavior
if key in existing_parameters:
    if parameters_are_identical:
        pass  # No action needed
    else:
        key.id += "_DUPLICATE"  # Keep both parameter sets
```

This ensures that:
- Protein atoms use ff14SB parameters
- Small molecule atoms use OpenFF Sage parameters
- The simulation runs correctly with appropriate parameters for each molecule type

#### What You'll See in Logs

With PolyzyMD's logging, you can identify which molecule combinations trigger collisions:

```
Combining 7 component Interchange(s)
  Components: LipA, ResorufinButyrate, EGPMA-SBMA_AAABA, ..., dmso, water/ions
[DEBUG] Combining 'LipA' with 'ResorufinButyrate'...
Key collision with different parameters, fixing. Key is [#6X4:1]-[#1:2]
...
```

Collisions typically occur when combining protein Interchanges (using ff14SB) with small molecule Interchanges (using OpenFF Sage).

#### Further Reading

For more details on this behavior, see the OpenFF Interchange documentation:
- [Sharp Edges: Combining Interchanges](https://docs.openforcefield.org/projects/interchange/en/stable/using/edges.html)

---

## Complete Example

See the example configurations in `src/polyzymd/configs/examples/`:

- `enzyme_only.yaml` - Enzyme + substrate, no polymers
- `enzyme_polymer.yaml` - Full enzyme + polymer simulation
- `enzyme_cosolvent.yaml` - Enzyme with DMSO co-solvent

---

## See Also

- {doc}`dynamic_polymers` - Dynamic polymer generation from SMILES
- {doc}`gromacs_export` - Running simulations with GROMACS
- {doc}`polymers` - Polymer setup guide
- {doc}`restraints` - Atom selection and restraints
- {doc}`cli_reference` - CLI documentation
