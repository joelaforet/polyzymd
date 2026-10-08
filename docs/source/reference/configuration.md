# Configuration Reference

This document describes all configuration options for PolyzyMD YAML files.

## Configuration Structure

A complete configuration file has these sections:

```yaml
name: "simulation_name"   # Required
engine: "openmm"          # Required: openmm or gromacs
description: "optional description"

enzyme: { ... }           # Required
substrate: { ... }        # Optional (null for apo)
polymers: { ... }         # Optional (null to disable)
solvent: { ... }          # Optional (default: TIP3P water, neutralizing ions)
restraints: [ ... ]       # Optional
thermodynamics: { ... }   # Required
simulation_phases: { ... } # Required
output: { ... }           # Optional (has defaults)
force_field: { ... }      # Optional (has defaults)
openmm: { ... }           # Optional (has defaults; read when engine is openmm)
gromacs: { ... }          # Optional (has defaults; read when engine is gromacs)
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
| `custom_substructures_path` | path | No | JSON file of residue templates for residues that OpenFF does not know, such as an N-terminal cystine. See {doc}`openff_pdb_ingestion` |
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
| `sdf_path` | path | Yes | - | Path to SDF with docked conformers. To make it from a crystal ligand, see {ref}`ligand-sdf-from-crystal` |
| `conformer_index` | int | No | 0 | Index of conformer to use (0-indexed) |
| `charge_method` | string | No | "nagl" | Options: `nagl`, `espaloma`, `am1bcc` |
| `residue_name` | string | No | "LIG" | 3-letter code for topology |

### Charge Methods

| Method | Description | Speed |
|--------|-------------|-------|
| `nagl` | Graph neural network charges | Fast |
| `espaloma` | Machine learning charges | Medium |
| `am1bcc` | Semi-empirical QM charges | Slow |

```{note}
`nagl`/OpenFF charging is the recommended default for routine PolyzyMD v1.3
workflows. `am1bcc` depends on AmberTools availability, which is
platform- and environment-specific and is not available in all default pixi
environments, especially macOS NumPy 2 environments. If you need AM1-BCC, use
pre-charged molecules or request/construct a dedicated AmberTools environment
for your platform.
```

### No Substrate (Apo Simulation)

```yaml
substrate: null
```

---

## Polymer Configuration

PolyzyMD supports two modes for polymer generation: **cached** (load from pre-built SDF files) and **dynamic** (generate on-the-fly from SMILES).

```{tip}
For a complete guide on dynamic polymer generation, see {doc}`../how_to/dynamic_polymers`.
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
  
  sdf_directory: "polymer_sdfs/SBMA-EGPMA" # Pre-built polymer SDFs (required in cached mode)
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
| `sdf_directory` | path | In cached mode | null | Directory with pre-built polymer SDFs |
| `cache_directory` | path | No | ".polymer_cache" | Cache directory |
| `reactions` | object | No | all "default" | ATRP reaction templates (dynamic mode) |
| `charger` | string | No | "nagl" | Charge method for dynamic generation |
| `max_retries` | int | No | 10 | Max attempts for ring-piercing avoidance |
| `random_seed` | int | No | null | Seed for polymer sequence draws (default: the replicate number) |
| `packing` | object | No | see below | PACKMOL placement settings |

### Packing Options (`polymers.packing`)

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `padding` | float (nm) | No | 2.0 | Room reserved for the polymer chains around the solute. Added to `solvent.box.padding` when the periodic cell is computed, and used as the padding of the confinement sphere |
| `tolerance` | float (Å) | No | 2.0 | PACKMOL minimum distance between any two atoms, including polymer to protein |
| `movebadrandom` | bool | No | false | Pass PACKMOL's `movebadrandom`; helps dense or heterogeneous systems converge |
| `confine_to_sphere` | bool | No | true | Confine chains to a sphere centred on the solute (radius = half the solute bounding-box diagonal + `padding`) while packing inside the final periodic brick. Set `false` to let chains fill the whole brick |
| `nloop` | int | No | 200 | Maximum PACKMOL GENCAN loops per molecule type |

PACKMOL is seeded with the replicate number for both polymer packing and
solvation, so replicates start from independent coordinates. The seeds are
recorded under `provenance` in `build_manifest.json`.

#### Where the polymers are packed

The periodic cell is computed **before** anything is packed, from the protein
and substrate alone:

```
edge = solute diameter + 2 * (polymers.packing.padding + solvent.box.padding)
box vectors = edge * shape_matrix
```

Chains are then packed inside the rectangular *brick* of that cell (shrunk by
`tolerance`, PACKMOL's own convention) plus the confinement sphere, and the
solvent fills the same brick afterwards. Because the cell depends only on the
enzyme and substrate, every replicate of a condition gets the same box volume,
the same water count and the same ion counts. The cell, the brick and the
sphere radius are recorded under `provenance` in `build_manifest.json`
(`box_vectors_nm`, `brick_nm`, `polymer_sphere_radius_nm`), so replicates can
be compared by diffing their manifests.

```{note}
A rhombic-dodecahedron brick is one edge long along `x` and `y`, but only
`sqrt(2)/2` (0.707) edges tall along `z`. The missing corners are supplied by
the periodic images. If the solute's bounding box would come closer than
`solvent.box.tolerance` to a brick face, the edge grows until it fits. The
build log and `build --dry-run` print the edge, the brick and the clearance to
each brick face.
```

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

```{note}
`nagl`/OpenFF charging is the recommended default for dynamic polymer generation.
`am1bcc` depends on AmberTools availability, which is platform- and
environment-specific and is not available in all default pixi environments,
especially macOS NumPy 2 environments. If you need AM1-BCC, use pre-charged
molecules or request/construct a dedicated AmberTools environment for your
platform.
```

### No Polymers

```yaml
polymers: null
```

You can also leave out the `polymers:` section. A section with
`enabled: false` must still contain `type_prefix`, `monomers`, `length` and
`count`, because the schema requires them.

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
    nacl_concentration: 0.15             # NaCl salt concentration (M)
  
  box:
    padding: 1.2                         # nm from solute to box edge (see below)
    shape: "rhombic_dodecahedron"        # Box shape
    target_density: 1.0                  # g/mL
    tolerance: 2.0                       # PACKMOL tolerance (Angstrom)
```

`nacl_concentration` sets the number of NaCl pairs. With `neutralize: true`,
the Na+ or Cl- ions that cancel the charge of the solute and co-solvents are
added on top of the salt, as OpenMM `Modeller` and `gmx genion -neutral` do.

`box.padding` is the distance from the **solute** to the box edge. The edge of
the cell is the solute diameter (its largest atom-to-atom distance) plus
`2 * padding`, as in `gmx editconf -d`. Every lattice vector of a cube or a
rhombic dodecahedron is at least one edge long. The solute
therefore starts at least `2 * padding` from each of its periodic copies, in
any orientation. The solute's bounding box is centred in the brick. When
polymers are configured, `polymers.packing.padding` is added to it and the
resulting cell is computed from the protein and substrate before any packing
happens, so it is identical across replicates of a condition; the packed
topology is not re-centred afterwards. Without polymers the box is computed
from the solute at solvation time. The number of waters and ions follows from the box
volume and `target_density`, so a deterministic box means deterministic
solvent counts.

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
| `cube` | Cube: three equal edges at right angles |
| `rhombic_dodecahedron` | Same edge, 71 % of the cube volume (default) |

Both shapes use the same edge, so they keep the solute equally far from its
periodic copies. The rhombic dodecahedron needs fewer waters.

### Co-solvents

PolyzyMD supports adding co-solvents to a water primary solvent. Give the amount of each as a **mole fraction**, a **molar concentration** or a **count** of molecules.

#### Specification Methods

| Method | Field | Description | Effect on Water |
|--------|-------|-------------|-----------------|
| Mole Fraction | `mole_fraction` | Fraction of neutral solvent molecules (0-1) | Replaces water in the neutral solvent mixture |
| Concentration | `concentration` | Molar concentration (mol/L) | Added on top of the water; only the ions of a charged co-solvent reduce the water (see below) |
| Count | `count` | Number of molecules in the box | Added on top of the water; only the ions of a charged co-solvent reduce the water (see below) |

**Important:** Use exactly ONE method per co-solvent. The previous `volume_fraction` key has been removed and is rejected instead of converted automatically. Existing configs that used `volume_fraction` must be updated explicitly to either `mole_fraction` or `concentration`; PolyzyMD does not infer mole fractions from volume fractions.

#### Mole Fraction Method

Use this when you want a specific molecule fraction of the neutral solvent mixture to be the co-solvent (e.g., "10 mol% DMSO").

```yaml
co_solvents:
  - name: "dmso"
    mole_fraction: 0.10      # 10 mol% DMSO
```

**Formula:**

```
x_water = 1 - sum(x_i)
M_avg   = x_water × M_water + sum(x_i × M_i)
N_total = neutral_solvent_mass / M_avg
n_i     = x_i × N_total
n_water = N_total - sum(n_i)

Where:
  x_i   = mole fraction for co-solvent i
  M_i   = molecular mass for co-solvent i
```

**Source:** [`src/polyzymd/builders/solvent.py`](https://github.com/joelaforet/polyzymd/blob/main/src/polyzymd/builders/solvent.py)

Mole-fraction co-solvents replace water in the neutral solvent mass budget. If you specify 10 mol% DMSO, DMSO molecules account for about 10% of neutral solvent molecules after rounding. Naming-template placeholders render mole fractions as mol-percent tokens with the `molpct` suffix; for example, `mole_fraction: 0.30` for DMSO renders as `dmso_30molpct`, and `{solvent_composition}` renders as `water_tip3p_dmso_30molpct`.

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
  N_A   = Avogadro's number
```

**Source:** [`src/polyzymd/builders/solvent.py`](https://github.com/joelaforet/polyzymd/blob/main/src/polyzymd/builders/solvent.py)

With `concentration` or `count`, the co-solvent molecules are added on top of the water, which slightly increases the density. Their mass does not reduce the water count. Ions do: the water fills the solvent mass that is left after all Na+ and Cl- ions. So a charged co-solvent, whose counter-ions or neutralizing ions are added as Na+ or Cl-, also reduces the water. For example, 8 dodecyl sulfate anions by `count` add 8 Na+ and give about 10 fewer waters than the same box without them.

#### Built-in Co-solvent Library

PolyzyMD includes a library of common co-solvents with pre-defined SMILES and densities. Density values are retained as metadata and are sourced from [PubChem](https://pubchem.ncbi.nlm.nih.gov/), a public database of chemical compounds. Each compound has a unique Compound Identification Number (CID) that can be used to look up detailed information including density, structure, and safety data.

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

For library co-solvents, you only need to specify the `name` and either `mole_fraction` or `concentration`:

```yaml
co_solvents:
  - name: "dmso"
    mole_fraction: 0.10      # 10 mol% DMSO - smiles auto-populated
```

#### Custom Co-solvents

For molecules not in the library, you must provide the SMILES string:

```yaml
co_solvents:
  # Custom co-solvent with mole fraction
  - name: "ethyl_acetate"
    smiles: "CCOC(=O)C"
    mole_fraction: 0.05

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
    mole_fraction: 0.10      # 10 mol% DMSO
  - name: "urea"
    concentration: 1.0       # Plus 1 M urea
```

**Warning:** When using multiple co-solvents with `mole_fraction`, ensure the total is less than 1.0 (100%). The remaining mole fraction is filled with water. Non-water primary solvents and 100% DMSO systems are out of scope for now; model solvent mixtures as water primary solvent plus `co_solvents` entries with `mole_fraction < 1.0`.

```{warning}
**YAML List Syntax**

A common mistake is placing each field on a separate line with its own `-`, which creates multiple list items instead of one object with multiple fields.

**Incorrect** (creates 3 separate incomplete items):
~~~yaml
co_solvents:
  - name: "dmso"
  - mole_fraction: 0.10
  - residue_name: "DMS"
~~~

**Correct** (one item with 3 fields):
~~~yaml
co_solvents:
  - name: "dmso"
    mole_fraction: 0.10
    residue_name: "DMS"
~~~

The `-` character starts a **new list item**. All fields belonging to the same item must be indented to the same level *without* a leading `-`.
```

#### Assumptions and Limitations

- **Ideal composition:** Mole fractions are converted to molecule counts by a weighted average molar mass. Real solutions may deviate from ideal density behavior.
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
2. Check bundled library: src/polyzymd/data/solvents/dmso.sdf (used when no SMILES or the library SMILES is given)
3. Check user cache: ~/.polyzymd/solvent_cache/<name>.<charge method>.<SMILES hash>.sdf
4. Generate from SMILES, charge with charge_method (default nagl), save to user cache
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

When you use a custom co-solvent (not in the library), PolyzyMD:

1. Generates the molecule from your SMILES string.
2. Assigns partial charges with `charge_method`: `nagl` by default, the
   OpenFF graph network trained on AM1-BCC. Set `charge_method: am1bcc` for
   AM1-BCC itself, which needs AmberTools. When you leave the default, the
   build writes a warning that asks you to check that NAGL charges suit the
   molecule.
3. Caches the charged molecule in `~/.polyzymd/solvent_cache/`, one file per
   SMILES and charge method, and reuses it. A changed SMILES gives a new
   molecule, even under the same `name`.

A charged molecule can carry Na+ or Cl- counter-ions in its SMILES
(`...[O-].[Na+]`). PolyzyMD splits them off: the co-solvent molecule has no
ion atoms, and each counter-ion is added as a Na+ or Cl- ion like the salt
ions. Without the counter-ion, `solvent.ions.neutralize: true` adds the ions
that cancel the net charge of the solute and co-solvents together. Both
spellings give a neutral system with the requested salt. Other counter-ions
are refused; leave them out of the SMILES.

```yaml
co_solvents:
  - name: "sds"
    smiles: "CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]"   # sodium dodecyl sulfate with Na+
    residue_name: "SDS"
    count: 8
    charge_method: nagl                          # the default; am1bcc needs AmberTools
```

#### Managing the Cache

You can inspect and manage the solvent cache programmatically:

```python
from polyzymd.data import list_available_solvents, clear_cache

# Map each available solvent name to its source
solvents = list_available_solvents()
print(solvents)
# {'tip3p': 'library', 'spce': 'built-in', 'dmso': 'library', ...,
#  'my_custom_solvent.nagl.9a7ebe92ac51': 'user_cache'}

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
      selection: "resid 76 and name OG"  # First atom selection
      description: "Catalytic serine"    # Optional description
    atom2:
      selection: "resname LIG and name C1"
      description: "Substrate carbon"
    distance: 3.3                        # Angstroms
    force_constant: 10000.0              # kJ/mol/nm²
    enabled: true                        # Enable/disable
```

See {doc}`../how_to/restraints` for detailed selection syntax.

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
  equilibration_stages:
    - name: "heating"
      samples: 20                        # frames to save
      ensemble: "NVT"
      temperature_start: 60.0            # starting temperature (K)
      temperature_end: 300.0             # final temperature (K)
      temperature_increment: 1.0         # increase per update (K)
      temperature_interval_steps: 600    # MD steps between updates
      position_restraints:
        - group: "protein_heavy"
          force_constant: 4184.0
    - name: "free_equilibration"
      duration: 0.8                      # nanoseconds
      samples: 80
      ensemble: "NPT"
      temperature: 300.0

  production:
    ensemble: "NPT"
    duration: 100.0                      # nanoseconds total
    samples: 2500                        # total frames
    time_step: 2.0
    thermostat: "LangevinMiddle"
    thermostat_timescale: 1.0
    barostat: "MC"                       # Monte Carlo barostat
    barostat_frequency: 25               # steps between barostat moves
  
```

Equilibration cannot be skipped: list at least one stage in
`equilibration_stages`, with a duration above 0 ns. Any duration above 0 is
accepted. The quickstart uses one stage of 0.002 ns.

### Minimization (`simulation_phases.minimization`)

Energy minimization runs before the first equilibration stage. It is a
runtime-only setting: changing it does not alter the built system or its
`build_manifest.json` hash.

```yaml
simulation_phases:
  minimization:
    freeze_solute: true      # hold protein + substrate heavy atoms fixed
    max_iterations: 1000     # 0 = run to convergence
    tolerance: 10.0          # kJ/mol/nm
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `freeze_solute` | bool | true | Freeze every protein and substrate **heavy** atom (the `solute_heavy` atom group) during minimization; only solvent, polymers, and the solute hydrogens relax. The prepared heavy-atom structure enters equilibration with its coordinates unchanged; PolyzyMD verifies a zero heavy-atom displacement and records it, together with `hydrogen_max_displacement_angstrom`, in `minimization/phase.json`. Solute hydrogens stay mobile on purpose so the minimizer places them on the force field's X–H constraint lengths rather than keeping the input PDB's |
| `max_iterations` | int | 1000 | Maximum minimizer iterations (0 = until convergence) |
| `tolerance` | float | 10.0 | Energy tolerance in kJ/mol/nm |

See {doc}`../explanation/simulation_safeguards` for why the solute is frozen.

```{note}
The trajectory frame interval is derived from `production.duration` and
`production.samples`. A config with `production.report_interval` is refused.
```

For a temperature ramp, omit `duration`. Set `temperature_increment` in K and
`temperature_interval_steps` in MD steps. PolyzyMD computes how many updates are
needed to reach the endpoint and derives the duration. Constant-temperature
stages continue to require an explicit `duration`.

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
| `Andersen` | Andersen thermostat. Only with `engine: gromacs`; the OpenMM engine runs `LangevinMiddle` instead, with a warning at run time. |
| `NoseHoover` | Nosé-Hoover chain. Only with `engine: gromacs`; the OpenMM engine runs `LangevinMiddle` instead, with a warning at run time. |

### Barostats

| Barostat | Description |
|----------|-------------|
| `MC` | Monte Carlo barostat (recommended) |
| `MCA` | anisotropic pressure coupling. Only with `engine: gromacs`; `validate` refuses it with `engine: openmm`. |

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
  
```

The schema also accepts `save_checkpoint`, `save_state_data` and
`trajectory_format`, but no code reads them. The OpenMM engine always writes
DCD trajectories and checkpoints; the GROMACS engine always writes XTC.

### Naming Template Variables

| Variable | Description | Example |
|----------|-------------|---------|
| `{enzyme}` | Enzyme name | "LipA" |
| `{substrate}` | Substrate name (hyphens removed), or `apo` without a substrate | "ResorufinButyrate" |
| `{polymer_type}` | Polymer type prefix and composition, or `none` without polymers | "SBMA-EGPMA_A70_B30" |
| `{duration}` | Production duration in ns: whole ns from 1 ns up, in full below 1 ns | "100", "0.005" |
| `{temperature}` | Temperature in K | "300" |
| `{replicate}` | Replicate number | "1" |
| `{primary_solvent}` | Primary solvent token | "water_tip3p" |
| `{cosolvent_composition}` | Co-solvents sorted by normalized name, or `none` | "dmso_30molpct_urea_2p5M" |
| `{solvent_composition}` | Primary solvent plus co-solvents when present | "water_tip3p_dmso_30molpct" |

Mole-fraction co-solvents use mol-percent naming tokens with the `molpct`
suffix. Concentration-based co-solvents use molarity tokens such as `2p5M`.
The removed `volume_fraction` field is not part of naming-template resolution.

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

## OpenMM Engine Configuration

The optional `openmm:` block selects the OpenMM platform. It is used when
`engine` is `openmm`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `platform` | `str` | `"CUDA"` | OpenMM platform: `CUDA`, `OpenCL`, `CPU` or `Reference`. PolyzyMD never falls back to another platform. |
| `device_index` | `str \| null` | `null` | GPU device index. |
| `precision` | `str` | `"mixed"` | Floating-point precision on CUDA. |

The replicate number fixes the starting structure and every random seed. On
the CPU, OpenCL and CUDA platforms, CPU threads, PME and GPUs add up forces in
a different order on each run, so two runs of a replicate differ from the
first minimization on (the slow Reference platform is deterministic). They agree
statistically, not frame by frame. Each production segment records the
platform and the property values it used under `openmm_platform` in
`progress.json`.

---

(config-gromacs)=
## GROMACS Engine Configuration

:::{versionadded} 1.3.0
:::

The optional `gromacs:` block sets how PolyzyMD calls GROMACS and which SLURM
resources a GROMACS job asks for. `analysis_topology` is read by the analyses
and by freeze. Every other field is used only by the SLURM scripts that
`polyzymd submit` and `polyzymd recover` write when the engine is GROMACS
(`engine: gromacs` in the config, or `--engine gromacs`). A local
`polyzymd run` and the `run_<prefix>_gromacs.sh` script that
`build --format gromacs` writes do not read those fields: they call `gmx` (or
`--gmx-path`) as `gmx mdrun -deffnm <stage> -v` with no extra flags. To use
other flags locally, run the GROMACS commands yourself.

### Minimal Example

```yaml
gromacs:
  module_load: "module load gcc/11.2.0 gromacs/2024.2"
  ntmpi: 1
  ntomp: 8
```

### GPU Example

```yaml
gromacs:
  gpu: true
  gpus: 1
  gmx_binary: "gmx"
  ntmpi: 1
  ntomp: 12
  module_load: "module load gcc/11.2.0 gromacs/2024.2"
  mdrun_flags: "-pin on"
```

`gpu: true` adds `-nb gpu`, `-pme gpu` and `-bonded gpu` unless
`mdrun_flags` sets those flags. It never adds `-update gpu`. See the notes.

### Full Field Reference

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `gmx_binary` | `str \| null` | `null` | GROMACS binary path or name. When null, resolved via `$GMX_BIN` environment variable or PATH discovery. |
| `analysis_topology` | `str \| null` | `null` | File name of the run's `.top` in the run folder. Analyses read it when MDAnalysis cannot read `prod.tpr`, and freeze deposits it with the files it includes. When null, PolyzyMD uses `<prefix>.top`, the file it writes, so another `.top` in the folder does not matter. |
| `mdrun_flags` | `str` | `""` | Extra flags passed to `gmx mdrun` for all stages of a SLURM job. |
| `mdrun_flags_equilibration` | `str \| null` | `null` | Override `mdrun_flags` for equilibration stages only. Falls back to `mdrun_flags` when null. |
| `mdrun_flags_production` | `str \| null` | `null` | Override `mdrun_flags` for production only. Falls back to `mdrun_flags` when null. |
| `grompp_flags` | `str` | `""` | Extra flags passed to `gmx grompp`, such as `-maxwarn 1` to accept a warning you have read. By default every warning stops the run. |
| `command_prefix` | `str \| null` | `null` | Prefix prepended to all GROMACS commands. Use for container wrappers (e.g., `singularity exec ...`). When set with a real-MPI binary, automatic `mpirun` wrapping is skipped. |
| `mpi_launcher_flags` | `str` | `""` | Extra flags for the MPI launcher (`mpirun`). Only used with real-MPI builds (`gmx_mpi`). |
| `module_load` | `str \| null` | `null` | Module load command inserted verbatim into SLURM scripts. The job runs it; the submitting host does not, so load the scheduler module in your shell before `submit`. List prerequisites before the GROMACS module. |
| `env_exports` | `dict[str, str]` | `{}` | Environment variables exported before GROMACS commands. Keys must be valid shell variable names. |
| `setup_commands` | `list[str]` | `[]` | Shell commands run after `module_load` and before GROMACS commands. |
| `ntmpi` | `int` | `1` | Number of MPI ranks for `gmx mdrun -ntmpi`. Also sets SLURM `--ntasks` unless `slurm_ntasks` overrides it. Must be >= 1. |
| `slurm_ntasks` | `int \| null` | `null` | Override SLURM `--ntasks` independently of GROMACS `-ntmpi`. For multi-node MPI+GPU workflows where scheduler tasks differ from thread-MPI ranks. Must be >= 1 when set. |
| `ntomp` | `int` | `8` | OpenMP threads per rank for `gmx mdrun -ntomp`. Sets SLURM `--cpus-per-task`. Must be >= 1. |
| `gpu` | `bool` | `false` | Request GPU via SLURM. When false, the `--gres=gpu` directive is omitted entirely. |
| `gpus` | `int` | `1` | Number of GPUs to request when `gpu` is true. Ignored when `gpu` is false. Must be >= 1. |
| `memory` | `str` | `"16G"` | SLURM `--mem` allocation for GROMACS jobs. |

### Notes

- Unsafe GPU flags (`-pme gpu`, `-bonded gpu`, `-update gpu`) are automatically
  stripped during energy minimization stages. Only `-nb gpu` is safe for EM.
- GROMACS updates on the GPU (`-update gpu`) only with `integrator = md`.
  The Langevin thermostats (the default) run as `integrator = sd`, and
  GROMACS then stops with *"Only the md integrator is supported"*. The GROMACS
  documentation lists the other conditions.
- When `gpu` is true and `ntmpi` > 1, a warning is emitted about GPU sharing.
- Set `slurm_ntasks` above `ntmpi` when the scheduler must reserve more tasks
  than GROMACS runs ranks, for example for a container or a multi-GPU allocation.
- If `mdrun_flags` contains `-ntmpi` or `-ntomp`, a warning is emitted when
  those values conflict with the explicit `ntmpi`/`ntomp` fields.

### Stage-specific mdrun flags

`mdrun_flags_equilibration` and `mdrun_flags_production` replace
`mdrun_flags` in their stages. When null (the default), the stage uses
`mdrun_flags`.

```yaml
gromacs:
  mdrun_flags: "-pin on"                      # all stages
  mdrun_flags_equilibration: "-pin on -dlb yes"
  mdrun_flags_production: "-pin on -dlb auto -nb gpu -pme gpu -bonded gpu"
```

### `command_prefix` and `mpi_launcher_flags`

Use one of the two, not both.

- `command_prefix` goes before every GROMACS command. Use it for a container
  or a site launcher:
  `command_prefix: "singularity exec --rocm /path/to/gromacs.sif"`. The value
  cannot hold shell variables such as `$PWD`; script generation refuses them.
- `mpi_launcher_flags` goes after the `mpirun` that PolyzyMD writes for a
  real-MPI binary: `mpi_launcher_flags: "-genv I_MPI_FABRICS shm:tcp"` gives
  `mpirun -genv I_MPI_FABRICS shm:tcp gmx_mpi mdrun ...`.

With `command_prefix` and a real-MPI binary, PolyzyMD writes no `mpirun`,
ignores `mpi_launcher_flags` and logs a warning.

### `env_exports` and `setup_commands`

The job script exports `env_exports`, then runs `setup_commands` in order,
after `module_load` and before the first GROMACS command:

```yaml
gromacs:
  env_exports:
    GMX_GPU_DD_COMMS: "true"
    OMP_PROC_BIND: "close"
  setup_commands:
    - "ulimit -s unlimited"
    - "source /opt/gromacs-2024/bin/GMXRC"
```

### mdrun flags

These `gmx mdrun` flags go in `mdrun_flags`, `mdrun_flags_equilibration` or
`mdrun_flags_production`.

| Flag | Description |
|------|-------------|
| `-ntmpi N` | Thread-MPI ranks (thread-MPI `gmx` builds). Set by `ntmpi`. |
| `-ntomp N` | OpenMP threads per rank. Set by `ntomp`. |
| `-npme N` | Dedicated PME ranks, for example `-npme 1` with 3 GPUs. |
| `-nb gpu` | Nonbonded forces on the GPU. Allowed in minimization. |
| `-pme gpu` | PME electrostatics on the GPU. Removed in minimization. |
| `-bonded gpu` | Bonded forces on the GPU. Removed in minimization. |
| `-update gpu` | Integration and constraints on the GPU. Needs `integrator = md`, so not with the Langevin thermostats. Removed in minimization. |
| `-pin on` | Pin threads to CPU cores. |
| `-pinstride N` | Stride between pinned threads. |
| `-dlb yes\|auto` | Dynamic load balancing. |
| `-gpu_id NNN` | GPU device IDs, for example `012` for 3 GPUs. |

For thread-MPI, real MPI and OpenMP threads, see
{doc}`../explanation/gromacs_parallelism`. For how to run and submit GROMACS
jobs, see {doc}`../how_to/run_gromacs`.

---

## Complete Example

Start from one of these complete configs:

- `polyzymd study add-condition "<label>" --new`, run in a study folder,
  writes `conditions/<label>/config.yaml`. This template has every section.
  The substrate, co-solvent, polymer and restraint sections are commented
  out, with a comment on each key.
- `examples/quickstart/config.yaml` (OpenMM) and
  `examples/quickstart/config_gromacs.yaml` (GROMACS), in the PolyzyMD
  repository, are a protein in water with NaCl. See
  {doc}`../get_started/quickstart`.

To add a substrate, co-solvents or polymers, copy the blocks of this page
into one of these configs.

---

## See Also

- {doc}`../how_to/dynamic_polymers` - Dynamic polymer generation from SMILES
- {doc}`../how_to/run_gromacs` - Run GROMACS simulations
- {doc}`../how_to/polymers` - Polymer setup guide
- {doc}`../how_to/restraints` - Atom selection and restraints
- {doc}`cli_reference` - CLI documentation
