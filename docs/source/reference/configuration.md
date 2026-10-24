# Configuration Reference

This document describes all configuration options for PolyzyMD YAML files.

## Configuration Structure

A complete configuration file has these sections:

```yaml
name: "simulation_name"     # Required
engine: "openmm"            # Required: openmm or gromacs
description: "optional description"

solutes: [ ... ]            # Required: what you add to the box (protein, ligand, polymers)
solvent: { ... }            # Required: force field, water model, co-solvents, ions, box
restraints: [ ... ]         # Optional
modeling_exceptions: [ ... ] # Optional: deliberate departures from a physics check
thermodynamics: { ... }     # Required
simulation_phases: { ... }  # Required
output: { ... }             # Optional (has defaults)
openmm: { ... }             # Optional (has defaults; read when engine is openmm)
gromacs: { ... }            # Optional (has defaults; read when engine is gromacs)
```

PolyzyMD assumes no force field, charge model or water model. Every solute
and the solvent name their own, so each choice is on record in the config
that built the system. The physics audit then checks that the built system
carries exactly these choices; see {doc}`physics_checks` and
{doc}`../explanation/physics_audit`.

A config written for 1.3.0rc4 or earlier (`enzyme:`, `substrate:`,
`polymers:`, `force_field:`) stops with an error that lists where each section
went. See {doc}`../how_to/rewrite_config_rc5`.

---

## Solutes

`solutes:` lists everything you add to the box. Each entry is one component.

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `type` | string | Yes | - | `protein`, `small_molecule` or `polymer`. It selects which keys the entry takes and how its physics is checked. It is never guessed from the other keys |
| `name` | string | No | the type | The handle that restraints, `modeling_exceptions`, the build record and analysis use. No spaces, `/` or `\`; not `solvent`, `ions`, `SOL` or `ION`. Names must be unique |
| `segid` | string | No | `PROT`, `LIG`, `POLY` | Segment ID of the component's atoms in the built PDB (columns 73-76), in analysis and in the GROMACS `components.ndx`. 1-4 characters, not `SOL` or `ION` |
| `force_field` | list of paths | Yes | - | The `.offxml` files that parameterize this component |
| `charges` | object | Yes | - | Where its partial charges come from; see below |

This release builds one solute of each type, and a protein is required. Every
solute but the protein uses the same force field file as `solvent`.
`validate` says so when a config goes beyond these limits; it never builds
something else.

### Charges

```yaml
charges: {method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}
```

| `method` | Meaning | Takes |
|----------|---------|-------|
| `library` | The force field's library charges (ff14SB's for a protein) | protein |
| `preset` | The charges stored in the component's SDF files, or shipped with PolyzyMD for a co-solvent such as DMSO | cached polymers, bundled co-solvents |
| `nagl` | OpenFF NAGL. `model` names the network file and is required | small molecules, dynamic polymers, co-solvents |
| `am1bcc` | AM1-BCC, computed for each molecule. Needs AmberTools | same as `nagl` |
| `espaloma` | Espaloma charges, computed for each molecule | same as `nagl` |

`model` is only for `nagl`. `openff-gnn-am1bcc-1.0.0.pt` is the model OpenFF
Sage 2.3.0 names in its own `NAGLCharges` section; PolyzyMD 1.3 release
candidates before rc.5 used `openff-gnn-am1bcc-0.1.0-rc.3.pt`. Every component
charged with `nagl` must name the same model in this release. Co-solvent
charges are cached per model, so a change never reuses charges from another
model.

`nagl` and OpenFF charging need no extra program. `am1bcc` depends on
AmberTools, which is not in every pixi environment (in particular macOS NumPy 2
environments). `validate` refuses `nagl` for elements other than H, C, N, O,
F, P, S, Cl, Br and I, and a method whose software is not installed.

### Protein (`type: protein`)

```yaml
solutes:
  - type: protein
    name: "LipA"                                  # optional; the run folder name uses it
    pdb_path: "structures/enzyme.pdb"
    force_field: ["ff14sb_off_impropers_0.0.4.offxml"]
    charges: {method: library}
    description: "Bacillus subtilis Lipase A"     # optional
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `pdb_path` | path | Yes | Prepared, protonated PDB file. `validate` refuses a missing file or one without ATOM/HETATM records |
| `charges` | object | Yes | `{method: library}`: a protein takes its force field's library charges |
| `custom_substructures_path` | path | No | JSON file of residue templates for residues that OpenFF does not know, such as an N-terminal cystine. See {doc}`openff_pdb_ingestion` |
| `description` | string | No | Human-readable description |

### Small molecule (`type: small_molecule`)

A ligand placed as given, such as a docked substrate. Leave the entry out for
an apo simulation.

```yaml
  - type: small_molecule
    name: "Resorufin-Butyrate"
    sdf_path: "structures/substrate.sdf"
    conformer_index: 0                            # which conformer to place
    residue_name: "LIG"
    force_field: ["openff-2.3.0.offxml"]
    charges: {method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}
```

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `sdf_path` | path | Yes | - | SDF with the docked, protonated conformers. To make it from a crystal ligand, see {ref}`ligand-sdf-from-crystal` |
| `conformer_index` | int | No | 0 | Index of the conformer to place (0-indexed). `validate` refuses an index the SDF does not hold |
| `residue_name` | string | No | "LIG" | 3-letter code for the topology |
| `charges` | object | Yes | - | `nagl` (with `model`), `am1bcc` or `espaloma` |

---

## Polymers (`type: polymer`)

PolyzyMD supports two modes for polymer generation: **cached** (load from pre-built SDF files) and **dynamic** (generate on-the-fly from SMILES).

```{tip}
For a complete guide on dynamic polymer generation, see {doc}`../how_to/dynamic_polymers`.
```

### Basic Configuration (Cached Mode)

```yaml
solutes:
  # ... the protein entry ...
  - type: polymer
    type_prefix: "SBMA-EGPMA"                # Polymer type identifier
    monomers:                                # Monomer definitions
      - label: "A"                           # Single character label
        probability: 0.98                    # Selection probability (0-1)
        name: "SBMA"                         # Full name (optional)
      - label: "B"
        probability: 0.02
        name: "EGPMA"
    length: 5                                # Monomers per chain
    count: 2                                 # Number of polymer chains
    sdf_directory: "polymer_sdfs/SBMA-EGPMA" # Pre-built polymer SDFs (required in cached mode)
    cache_directory: ".polymer_cache"        # Cache for generated polymers
    force_field: ["openff-2.3.0.offxml"]
    charges: {method: preset}                # the charges stored in the SDFs
```

The entry is one item of `solutes:`; all its keys are indented under the `-`.

### Dynamic Generation Configuration

To generate polymers on-the-fly from monomer SMILES (without pre-built SDF files):

```yaml
  - type: polymer
    generation_mode: "dynamic"               # Enable dynamic generation
    type_prefix: "SBMA-EGPMA"
    reactions:                               # ATRP reaction templates
      initiation: "default"                  # or "/path/to/custom.rxn"
      polymerization: "default"
      termination: "default"
    monomers:
      - label: "A"
        probability: 0.7
        name: "SBMA"
        smiles: "[H]C([H])=C(C(=O)OC...)..." # Required for dynamic mode
        residue_name: "SBM"                  # Optional 3-letter residue name
      - label: "B"
        probability: 0.3
        name: "EGPMA"
        smiles: "[H]C([H])=C(C(=O)OC...)..."
        residue_name: "EGM"
    length: 5
    count: 2
    force_field: ["openff-2.3.0.offxml"]
    charges: {method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}   # or am1bcc, espaloma
    max_retries: 10                          # Retries for ring-piercing detection
    cache_directory: ".polymer_cache"
```

### All Polymer Options

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `generation_mode` | string | No | "cached" | `cached` or `dynamic` |
| `type_prefix` | string | Yes | - | Identifier for polymer type |
| `monomers` | list | Yes | - | Monomer specifications |
| `length` | int | Yes | - | Chain length (number of monomers) |
| `count` | int | Yes | - | Number of chains to add |
| `sdf_directory` | path | In cached mode | null | Directory with pre-built polymer SDFs |
| `cache_directory` | path | No | ".polymer_cache" | Cache directory |
| `reactions` | object | No | all "default" | ATRP reaction templates (dynamic mode) |
| `force_field` | list | Yes | - | Force field files; the same file as `solvent.force_field` |
| `charges` | object | Yes | - | `preset` in cached mode (the charges in the SDFs); `nagl` (with `model`), `am1bcc` or `espaloma` in dynamic mode |
| `max_retries` | int | No | 10 | Max attempts for ring-piercing avoidance |
| `random_seed` | int | No | null | Seed for polymer sequence draws (default: the replicate number) |
| `packing` | object | No | see below | PACKMOL placement settings |

### Packing Options (`packing`)

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `padding` | float (nm) | No | 2.0 | Room reserved for the polymer chains around the solute. Added to `solvent.box.padding` when the periodic cell is computed, and used as the padding of the confinement sphere |
| `tolerance` | float (Å) | No | 2.0 | PACKMOL minimum distance between any two atoms, including polymer to protein |
| `movebadrandom` | bool | No | false | Pass PACKMOL's `movebadrandom`; helps dense or heterogeneous systems converge |
| `confine_to_sphere` | bool | No | true | Confine chains to a sphere centred on the solute (radius = half the solute bounding-box diagonal + `padding`) while packing inside the final periodic brick. Set `false` to let chains fill the whole brick |
| `nloop` | int | No | 200 | Maximum PACKMOL GENCAN loops per molecule type |

PACKMOL is seeded with the replicate number for both polymer packing and
solvation, so replicates start from independent coordinates. The seeds are
recorded under `provenance` in `build_provenance.json`.

#### Where the polymers are packed

The periodic cell is computed **before** anything is packed, from the protein
and substrate alone:

```
edge = solute diameter + 2 * (packing.padding + solvent.box.padding)
box vectors = edge * shape_matrix
```

Chains are then packed inside the rectangular *brick* of that cell (shrunk by
`tolerance`, PACKMOL's own convention) plus the confinement sphere, and the
solvent fills the same brick afterwards. Because the cell depends only on the
enzyme and substrate, every replicate of a condition gets the same box volume,
the same water count and the same ion counts. The cell, the brick and the
sphere radius are recorded under `provenance` in `build_provenance.json`
(`box_vectors_nm`, `brick_nm`, `polymer_sphere_radius_nm`), so replicates can
be compared by diffing their build records.

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

### Charges

Cached chains carry their charges in the SDF files: write
`charges: {method: preset}`. Dynamically generated chains are charged when
they are made: write `nagl` (with `model`), `am1bcc` or `espaloma`. `validate`
refuses any other combination, so the config never names a charge source the
build did not use. See [Charges](#charges).

### No Polymers

Leave the polymer entry out of `solutes:`.

---

## Solvent Configuration

```yaml
solvent:
  force_field: ["openff-2.3.0.offxml"]   # water, ions and co-solvents (required)
  water:
    model: "tip3p"                       # Water model (required)
  
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

`force_field` names the file that parameterizes water, ions and co-solvents;
Sage carries the TIP3P water and the ion parameters. `water.model` is
required: PolyzyMD assumes no water model. Water and co-solvents carry
segment ID `SOL` (co-solvents are told apart by residue name) and ions `ION`.

`nacl_concentration` sets the number of NaCl pairs. With `neutralize: true`,
the Na+ or Cl- ions that cancel the charge of the solute and co-solvents are
added on top of the salt, as OpenMM `Modeller` and `gmx genion -neutral` do.

`ionic_strength` (mol/L) replaces `nacl_concentration` when the target is the
ionic strength of the whole box, for example to match a buffer. The counter-ions
then count toward it: for a total of N monovalent ions and a solute charge Q,
the build adds (N - Q)/2 Na+ and (N + Q)/2 Cl-, the split of Machado and
Pantano (J. Chem. Theory Comput. 2020, 16, 1367; cite it when you use
`ionic_strength`, see {doc}`references`). N comes from
I = N / (2 V N_A), with V the solvent volume taken from the same
water-equivalent count as the NaCl pairs, and is never below |Q|, so the box
is always neutral. Charged co-solvents, such as a surfactant, count toward it
with their charge z (I = 1/2 sum c z^2), as do counter-ions written in their
SMILES; NaCl fills the rest and neutralizes the box. This is the formal ionic
strength of the molecules added: a surfactant above its critical micelle
concentration binds part of its counter-ions, which the formula does not
model. It needs `neutralize: true`. Example: 50 mM sodium phosphate at pH 7.5 has an
ionic strength of about 0.13 M, so `ionic_strength: 0.13` gives the box that
buffer's ionic strength with NaCl.

NaCl is the only salt the build adds: `validate`, `build`, `run` and `submit`
refuse a non-zero `kcl_concentration` or `mgcl2_concentration`.

`box.padding` is the distance from the **solute** to the box edge. The edge of
the cell is the solute diameter (its largest atom-to-atom distance) plus
`2 * padding`, as in `gmx editconf -d`. Every lattice vector of a cube or a
rhombic dodecahedron is at least one edge long. The solute
therefore starts at least `2 * padding` from each of its periodic copies, in
any orientation. The solute's bounding box is centred in the brick. When
a polymer solute is configured, its `packing.padding` is added to it and the
resulting cell is computed from the protein and substrate before any packing
happens, so it is identical across replicates of a condition; the packed
topology is not re-centred afterwards. Without polymers the box is computed
from the solute at solvation time. The number of waters and ions follows from the box
volume and `target_density`, so a deterministic box means deterministic
solvent counts.

### Water Models

`tip3p` (TIP3P) is the only water model PolyzyMD builds, and
`solvent.water.model` must name it. The force fields it
loads carry the TIP3P water parameters, so `validate`, `build`, `run` and
`submit` refuse any other `model`. A config of an existing run with another
model still loads for `status` and analysis.

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
    charges: {method: preset}
```

Every co-solvent names its charges. `preset` takes the charges PolyzyMD ships
for a library co-solvent (see the table below); `validate` refuses `preset`
for any other molecule, including a library name given another SMILES. A new
molecule takes `nagl` (with `model`), `am1bcc` or `espaloma`.

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
    charges: {method: preset}
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

| Name | SMILES | Density (g/mL) | Residue name | Reference |
|------|--------|----------------|--------------|-----------|
| `dmso` | `CS(=O)C` | 1.10 | `DMS` | [CID 679](https://pubchem.ncbi.nlm.nih.gov/compound/679) |
| `dmf` | `CN(C)C=O` | 0.95 | `DMF` | [CID 6228](https://pubchem.ncbi.nlm.nih.gov/compound/6228) |
| `acetonitrile` | `CC#N` | 0.786 | `CCN` | [CID 6342](https://pubchem.ncbi.nlm.nih.gov/compound/6342) |
| `urea` | `C(=O)(N)N` | 1.32 | `URE` | [CID 1176](https://pubchem.ncbi.nlm.nih.gov/compound/1176) |
| `ethanol` | `CCO` | 0.789 | `ETH` | [CID 702](https://pubchem.ncbi.nlm.nih.gov/compound/702) |
| `methanol` | `CO` | 0.792 | `MOH` | [CID 887](https://pubchem.ncbi.nlm.nih.gov/compound/887) |
| `glycerol` | `C(C(CO)O)O` | 1.261 | `GOL` | [CID 753](https://pubchem.ncbi.nlm.nih.gov/compound/753) |
| `isopropanol` | `CC(C)O` | 0.786 | `ISO` | [CID 3776](https://pubchem.ncbi.nlm.nih.gov/compound/3776) |
| `acetone` | `CC(=O)C` | 0.784 | `ACN` | [CID 180](https://pubchem.ncbi.nlm.nih.gov/compound/180) |
| `thf` | `C1CCOC1` | 0.883 | `THF` | [CID 8028](https://pubchem.ncbi.nlm.nih.gov/compound/8028) |
| `dioxane` | `C1COCCO1` | 1.033 | `DIO` | [CID 31275](https://pubchem.ncbi.nlm.nih.gov/compound/31275) |
| `ethylene_glycol` | `C(CO)O` | 1.114 | `EDO` | [CID 174](https://pubchem.ncbi.nlm.nih.gov/compound/174) |

A co-solvent's residue name must not be an amino-acid, nucleic-acid, water or ion name, because selections such as `protein` would then include it. Without `residue_name`, a custom co-solvent gets the first three letters of its name, or, when those clash, the first two letters and the next letter or digit that does not clash (`methylurea` becomes `MEH`, not `MET`); the dry run prints it. `validate` and `build` refuse a `residue_name` you set that clashes. Library co-solvents whose short name clashed use the PDB chemical-component code (glycerol is `GOL`, not `GLY`). Topologies built before PolyzyMD 1.3 still have the old names (`GLY` glycerol, `MET` methanol, `ACE` acetone and acetonitrile, `ETH` ethylene glycol); an older build has the clash when one of these residues has no `CA` atom, which MDAnalysis shows with `[r for r in u.select_atoms("resname GLY MET").residues if "CA" not in r.atoms.names]`.

For library co-solvents, you only need to specify the `name`, the amount and
`charges: {method: preset}`:

```yaml
co_solvents:
  - name: "dmso"
    mole_fraction: 0.10      # 10 mol% DMSO - smiles auto-populated
    charges: {method: preset}
```

#### Custom Co-solvents

For molecules not in the library, you must provide the SMILES string:

```yaml
co_solvents:
  # Custom co-solvent with mole fraction
  - name: "ethyl_acetate"
    smiles: "CCOC(=O)C"
    mole_fraction: 0.05
    charges: {method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}

  # Custom co-solvent with concentration (density not needed)
  - name: "my_additive"
    smiles: "CC(=O)NC"
    concentration: 0.5       # 0.5 M
    charges: {method: am1bcc}
```

#### Multiple Co-solvents

You can combine multiple co-solvents. Each can use either specification method independently:

```yaml
co_solvents:
  - name: "dmso"
    mole_fraction: 0.10      # 10 mol% DMSO
    charges: {method: preset}
  - name: "urea"
    concentration: 1.0       # Plus 1 M urea
    charges: {method: preset}
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

**Correct** (one item with 4 fields):
~~~yaml
co_solvents:
  - name: "dmso"
    mole_fraction: 0.10
    residue_name: "DMS"
    charges: {method: preset}
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
3. Check user cache: ~/.polyzymd/solvent_cache/<name>.<charge method>.<hash of SMILES, method and NAGL model>.sdf
4. Generate from SMILES, charge with the co-solvent's charges method, save to user cache
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
2. Assigns partial charges with the method its `charges` names: `nagl`
   (the OpenFF graph network trained on AM1-BCC, with the `model` you name),
   `am1bcc` (needs AmberTools) or `espaloma`. With `nagl`, the build writes a
   warning that asks you to check that NAGL charges suit the molecule.
3. Caches the charged molecule in `~/.polyzymd/solvent_cache/`, one file per
   SMILES, charge method and NAGL model, and reuses it. A changed SMILES gives a new
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
    charges: {method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}   # or am1bcc (needs AmberTools)
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
      component: "LipA"                  # The solute the atom belongs to
      selection: "resid 76 and name OG"  # Selection inside that component
      description: "Catalytic serine"    # Optional description
    atom2:
      component: "Resorufin-Butyrate"
      selection: "name C1"
      description: "Substrate carbon"
    distance: 3.3                        # Angstroms
    force_constant: 10000.0              # kJ/mol/nm²
    enabled: true                        # Enable/disable
```

Each atom names the solute it belongs to (`component`), and the selection is
matched only inside that component's atoms, so `resid 76` cannot pick a
polymer residue 76. Each selection must give exactly one atom. The physics
audit selects the atom again, with another library, and checks that the
restraint in the built system joins the atoms the config names
(`restraint.atom_identity`). See {doc}`../how_to/restraints` for detailed
selection syntax.

Only the OpenMM engine applies distance restraints. With `engine: gromacs`,
`validate` refuses an enabled restraint, and `build --format gromacs`,
`run --engine gromacs` and `submit --engine gromacs` refuse it too. Position
restraints in equilibration stages (`position_restraints:`) work on both
engines. On GROMACS they work for `protein_heavy`, `protein_backbone`,
`protein_calpha`, `ligand_heavy` and `polymer_heavy`; each group must have
the same `force_constant` in every stage, and all stages must use the same
protein group. `validate` refuses other position restraints on GROMACS, and
refuses `ligand_heavy` without a small-molecule solute or `polymer_heavy` without a polymer solute
on both engines.

### Restraint Types

| Type | Description |
|------|-------------|
| `flat_bottom` | No force within threshold, harmonic beyond |
| `harmonic` | Harmonic potential at target distance |
| `upper_wall` | Prevent distance exceeding threshold |
| `lower_wall` | Prevent distance below threshold |

---

## Modeling Exceptions

Every build runs the physics audit ({doc}`physics_checks`). A check that fails
stops the build; a check that warns is reported. When a failure or a warning
is deliberate, accept it here, with the reason. Crystal waters kept inside a
protein, for example, warn under `geometry.buried_water`:

```yaml
modeling_exceptions:
  - check: geometry.buried_water
    component: LipA          # optional; every component the check runs on when left out
    reason: "Crystal waters in the active-site pocket are kept on purpose."
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `check` | string | Yes | A check ID from {doc}`physics_checks`. `validate` refuses an unknown ID and names the closest ones |
| `component` | string | No | A solute name, `solvent`, `ions` or `system`. It must be a component the check runs on |
| `reason` | string | Yes | Why the departure is deliberate. No code reads it; it is there for the reader of the published config |

An accepted failure is recorded in `build_provenance.json` and
`build_provenance.md` with its reason and the atoms it covered. An entry that
accepted nothing is reported by the `modeling_exceptions.unused` check, so a
stale entry does not stay in a published config unnoticed.

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
`build_provenance.json` config hash.

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
`validate` refuses a phase with more `samples` than MD steps, or shorter than
one time step.
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
| `NVE` | Not run: both engines keep the thermostat on, so `validate` refuses it |

### Thermostats

| Thermostat | Description |
|------------|-------------|
| `LangevinMiddle` | Langevin integrator (recommended) |
| `Langevin` | Standard Langevin |
| `Andersen` | Andersen thermostat. Only with `engine: gromacs`; `validate` refuses it with `engine: openmm`. |
| `NoseHoover` | Nosé-Hoover chain. Only with `engine: gromacs`; `validate` refuses it with `engine: openmm`. |

### Barostats

| Barostat | Description |
|----------|-------------|
| `MC` | Monte Carlo barostat (recommended) |
| `MCA` | anisotropic pressure coupling. Only with `engine: gromacs`; `validate` refuses it with `engine: openmm`. |

---

## Output Configuration

Environment variables (`$USER`, `$HOME`, `${VAR}`) and `~` are automatically expanded in path fields. `validate` refuses a path that names an unset variable. Earlier versions did not expand `~` in input and output paths and wrote runs to a folder named `~` beside the config; while `~/...` does not exist, that folder is used, with a warning.

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

The schema still accepts `save_checkpoint`, `save_state_data` and
`trajectory_format` so the configs of existing runs load, but no code reads
them, and a config that gives one prints a warning. Both engines always write
checkpoints and energies; OpenMM writes DCD trajectories and GROMACS writes
XTC.

### Naming Template Variables

A config with another placeholder is refused, and `validate` refuses a
template without `{replicate}`, since the replicates would share one folder.

| Variable | Description | Example |
|----------|-------------|---------|
| `{enzyme}` | Name of the protein solute | "LipA" |
| `{substrate}` | Name of the small-molecule solute (hyphens removed), or `apo` without one | "ResorufinButyrate" |
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

## Force Fields

There is no `force_field:` section. Each solute names its own force field files
in `force_field:`, and `solvent.force_field` names the file for water, ions
and co-solvents:

```yaml
solutes:
  - type: protein
    force_field: ["ff14sb_off_impropers_0.0.4.offxml"]
    ...
solvent:
  force_field: ["openff-2.3.0.offxml"]
```

`validate` refuses a name that is not an installed force field or an
`.offxml` file. PolyzyMD recommends none: choose the force fields for your
system and cite them. Installed files include
`ff14sb_off_impropers_0.0.4.offxml` (Amber ff14SB, ported to OpenFF) and the
OpenFF Sage releases (`openff-2.0.0.offxml` to `openff-2.3.0.offxml`).
A config of 1.3.0rc4 or earlier that named no small-molecule force field was
built with `openff-2.0.0.offxml`; name that file to build the same physics.

Each component is parameterized by its own files alone. The
`parameters.bonded`, `parameters.nonbonded` and `charges.source` checks rebuild
every component from the files its entry names and compare the built system
with them term by term, so a term taken from another component's force field
fails the build.

### Key Collision Warnings

When building systems with both proteins and small molecules, you may see warnings like:

```
Key collision with different parameters, fixing. Key is [#6X4:1]-[#1:2]
```

**This is expected behavior and does not indicate a problem.**

#### Why This Happens

PolyzyMD uses different force fields for different molecule types:
- **Proteins**: ff14SB (Amber force field ported to OpenFF format)
- **Small molecules**: OpenFF Sage (general small molecule force field)

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
| `platform` | `str` | `"CUDA"` | OpenMM platform: `CUDA`, `OpenCL`, `CPU`, `Reference` or `HIP`; `validate` refuses other names. PolyzyMD never falls back to another platform. |
| `device_index` | `str \| null` | `null` | GPU device index. |
| `precision` | `str` | `"mixed"` | Floating-point precision on CUDA: `single`, `mixed` or `double`. |

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
| `memory` | `str` | `"16G"` | SLURM `--mem` allocation for GROMACS jobs: a number with an optional `K`, `M`, `G` or `T`. |

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
  cannot hold shell variables such as `$PWD`; `validate` refuses them.
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
  The small-molecule, polymer, co-solvent, restraint and
  `modeling_exceptions` entries are commented out, with a comment on each key.
- `examples/quickstart/config.yaml` (OpenMM) and
  `examples/quickstart/config_gromacs.yaml` (GROMACS), in the PolyzyMD
  repository, are a protein in water with NaCl. See
  {doc}`../get_started/quickstart`.

To add a small molecule, co-solvents or polymers, copy the entries of this
page into one of these configs.

---

## See Also

- {doc}`../how_to/dynamic_polymers` - Dynamic polymer generation from SMILES
- {doc}`../how_to/run_gromacs` - Run GROMACS simulations
- {doc}`../how_to/polymers` - Polymer setup guide
- {doc}`../how_to/restraints` - Atom selection and restraints
- {doc}`cli_reference` - CLI documentation
- {doc}`physics_checks` - The checks every build must pass
- {doc}`../how_to/rewrite_config_rc5` - Rewrite a 1.3.0rc4 config
