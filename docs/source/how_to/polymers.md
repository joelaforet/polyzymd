# Polymer Setup Guide

This guide covers configuring polymer chains in PolyzyMD simulations.

```{tip}
**Looking for dynamic polymer generation?** If you want native linear
methacrylate generation from SMILES strings, see the {doc}`dynamic_polymers`
guide. For the mBuild/OpenFF fragment adapter, see
{doc}`../reference/polymer_mbuild_openff`.
```

## Overview

PolyzyMD supports adding linear generated co-polymers and explicit charged SDF
molecules to your simulation box. Polymers are:

- Generated based on monomer probabilities
- Optionally assembled from explicit linear fragments
- Optionally combined with user-supplied charged SDF molecules through
  `provided_molecules`
- Placed around the enzyme using PACKMOL
- Parameterized with OpenFF force fields

## Basic Configuration

Polymer section snippet:

```yaml
polymers:
  enabled: true
  generation_mode: "dynamic"
  type_prefix: "SBMA-EGPMA"
  
  monomers:
    - label: "A"
      probability: 0.98
      name: "MMA"
      smiles: "C=C(C)C(=O)OC"
    - label: "B"
      probability: 0.02
      name: "EMA"
      smiles: "C=C(C)C(=O)OCC"
  
  length: 5      # 5-mer chains
  count: 2       # 2 polymer chains
  reactions:
    initiation: "default"
    polymerization: "default"
    termination: "default"
```

With the bundled `"default"` reactions, `generation_mode: "dynamic"` uses the
native linear methacrylate backend. Nonlinear, branched, and dendrimer molecules
are supplied as charged SDFs with `provided_molecules`; runtime CGSmiles and an
authoring notebook are future work, not current user interfaces.

---

## Monomer Definition

### Probability-Based Selection

Each chain is built by randomly selecting monomers based on their probabilities:

Snippet:

```yaml
monomers:
  - label: "A"
    probability: 0.98    # 98% chance
    name: "SBMA"
  - label: "B"
    probability: 0.02    # 2% chance
    name: "EGPMA"
```

```{important}
Probabilities must sum to 1.0 (100%).
```

### Multiple Monomers

You can define any number of monomer types:

Snippet:

```yaml
monomers:
  - label: "A"
    probability: 0.70
    name: "MonomerA"
  - label: "B"
    probability: 0.20
    name: "MonomerB"
  - label: "C"
    probability: 0.10
    name: "MonomerC"
```

### Homopolymers

For a homopolymer (single monomer type):

Snippet:

```yaml
monomers:
  - label: "A"
    probability: 1.0
    name: "PEG"
```

---

## Chain Configuration

### Chain Length

Number of monomers per chain:

Snippet:

```yaml
length: 5    # 5-mer (pentamer)
```

Typical values:
- **Short chains**: 3-5 monomers (faster simulations)
- **Medium chains**: 10-20 monomers
- **Long chains**: 50+ monomers (slower, more realistic)

### Number of Chains

Snippet:

```yaml
count: 2    # Add 2 polymer chains
```

More chains = larger system = slower simulation.

---

## Add Explicit Charged SDF Molecules

Use `provided_molecules` for explicit molecules that you have already authored,
charged, and validated. This is the preferred path for branched polymers,
dendrimers, nonlinear additives, and other molecules that PolyzyMD should pack
without regenerating or charging.

### Fixed Inventory

Use entry-level `count` values when the exact inventory is known:

Complete polymer section example; this assumes the referenced charged SDF files
already exist and satisfy the requirements in {doc}`../reference/data_requirements`:

```yaml
polymers:
  enabled: true
  generation_mode: "dynamic"
  type_prefix: "MMA"
  monomers:
    - label: "A"
      probability: 1.0
      name: "MethylMethacrylate"
      smiles: "C=C(C)C(=O)OC"
  length: 3
  count: 1
  reactions:
    initiation: "default"
    polymerization: "default"
    termination: "default"
  provided_molecules:
    - name: "charged_additives"
      entries:
        - sdf_path: "structures/branched_polymer_01.sdf"
          count: 1
        - sdf_path: "structures/branched_polymer_02.sdf"
          count: 2
```

Fixed pools must not set a pool-level `count` and must not use entry
probabilities.

### Probabilistic Pool

Use a pool-level `count`, per-entry probabilities, and optionally a pool-local
`seed` when PolyzyMD should sample from a validated molecule library:

Polymer section snippet:

```yaml
polymers:
  enabled: true
  generation_mode: "dynamic"
  type_prefix: "SBMA-EGPMA"
  # ... monomers, length, count, reactions ...
  random_seed: 2026
  provided_molecules:
    - name: "charged_oligomer_pool"
      count: 3
      seed: 17
      entries:
        - sdf_path: "structures/oligomer_a.sdf"
          probability: 0.7
        - sdf_path: "structures/oligomer_b.sdf"
          probability: 0.3
```

Probabilistic pools require probabilities on every entry, no entry-level counts,
and probabilities that sum to 1.0. Seed precedence is pool `seed`, then
`polymers.random_seed`, then the workflow or replicate seed.

### SDF Requirements

Every `provided_molecules` entry must be a charged single-molecule SDF. PolyzyMD
loads the molecule through OpenFF and validates that it contains exactly one
connected molecule, conformer coordinates, no dummy atoms, complete finite
partial charges, and finite coordinates. PolyzyMD does not autocharge these SDFs.

## Assemble Linear Explicit Fragments

Use `generation_mode: "fragments"` for native linear explicit-fragment assembly:

Complete polymer section example:

```yaml
polymers:
  enabled: true
  generation_mode: "fragments"
  type_prefix: "linear-carbon-fragment"
  monomers:
    - label: "A"
      probability: 1.0
      name: "CarbonFragment"
  fragments:
    A:
      terminal: "[*]C"
      middle: "[*:1]C[*:2]"
  length: 3
  count: 1
  charger: "nagl"
  cache_directory: ".polymer_cache"
```

Fragment mode is linear and sequence-directed:

- Schema validation checks that each monomer label has one `terminal` string and
  one `middle` string.
- Runtime parsing checks dummy atom counts, dummy atom degree, optional maps,
  RDKit parsing/embedding/sanitization, mBuild placement, and OpenFF
  conversion/charging.
- `terminal` strings must have exactly one `*`; `middle` strings must have
  exactly two `*` atoms.
- Optional `[*:1]` and `[*:2]` maps mark incoming and outgoing middle-fragment direction.
- PolyzyMD creates single inter-fragment bonds using mBuild `Port` objects and
  `force_overlap()`.
- The assembled molecule is converted to OpenFF, charged, and cached.

Do not add separate head or tail keys. PolyzyMD does not suggest chemistry for
fragment strings.

## Deprecated Sequence-Derived SDF Libraries

`sdf_directory` is deprecated compatibility behavior for old libraries whose
filenames are derived from generated sequences. It is not an alias for
`provided_molecules`.

Legacy compatibility snippet:

```yaml
polymers:
  enabled: true
  generation_mode: "cached"
  type_prefix: "SBMA-EGPMA"
  monomers:
    - label: "A"
      probability: 0.98
      name: "SBMA"
    - label: "B"
      probability: 0.02
      name: "EGPMA"
  length: 5
  count: 2
  sdf_directory: "legacy_polymer_sdfs"
```

Legacy files are looked up by generated sequence-derived names such as
`SBMA-EGPMA_seq=AABBA_5-mer_charged.sdf`. PolyzyMD emits a deprecation warning
when `sdf_directory` is used. For explicit charged SDF molecules, migrate to
`provided_molecules`.

---

## Polymer Cache

Generated polymers are cached for reuse:

Snippet:

```yaml
polymers:
  # ...
  cache_directory: ".polymer_cache"
```

This speeds up repeated runs with the same polymer sequences.

To clear the cache:

```bash
rm -rf .polymer_cache
```

---

## Example Configurations

### MMA-EMA Co-polymer

Two methacrylate monomers with different ester groups:

Polymer section snippet:

```yaml
polymers:
  enabled: true
  generation_mode: "dynamic"
  type_prefix: "MMA-EMA"
  monomers:
    - label: "A"
      probability: 0.98
      name: "MMA"           # Methyl methacrylate
      smiles: "C=C(C)C(=O)OC"
    - label: "B"
      probability: 0.02
      name: "EMA"           # Ethyl methacrylate
      smiles: "C=C(C)C(=O)OCC"
  length: 5
  count: 2
  reactions:
    initiation: "default"
    polymerization: "default"
    termination: "default"
```

### Methacrylate Homopolymer

Single-label native dynamic generation snippet:

```yaml
polymers:
  enabled: true
  generation_mode: "dynamic"
  type_prefix: "OEGMA"
  monomers:
    - label: "A"
      probability: 1.0
      name: "OEGMA"
      smiles: "C=C(C)C(=O)OCCO"
  length: 10
  count: 4
  reactions:
    initiation: "default"
    polymerization: "default"
    termination: "default"
```

### Block or Branched Additive

For a block-like, branched, or otherwise nonlinear structure that PolyzyMD should
not generate, provide a charged SDF through `provided_molecules`:

Polymer section snippet:

```yaml
polymers:
  enabled: true
  generation_mode: "dynamic"
  type_prefix: "SBMA-EGPMA"
  monomers:
    - label: "A"
      probability: 0.5
      name: "MMA"
      smiles: "C=C(C)C(=O)OC"
    - label: "B"
      probability: 0.5
      name: "EMA"
      smiles: "C=C(C)C(=O)OCC"
  length: 5
  count: 2
  reactions:
    initiation: "default"
    polymerization: "default"
    termination: "default"
  provided_molecules:
    - name: "block_additive"
      entries:
        - sdf_path: "structures/block_additive_charged.sdf"
          count: 1
```

---

## Disabling Polymers

### Control Simulations

For enzyme-only (control) simulations:

```yaml
polymers: null
```

Or explicitly:

```yaml
polymers:
  enabled: false
```

---

## Placement and Solvation

Polymers are placed in the simulation box using PACKMOL:

1. Enzyme (+ substrate) placed at center
2. Polymers placed around enzyme with minimum distance
3. Water molecules fill remaining space
4. Ions added to neutralize and reach target concentration

### Box Padding

The `solvent.box.padding` affects polymer placement:

```yaml
solvent:
  box:
    padding: 1.5    # nm - increase for more polymers
```

Larger padding = more space for polymers = larger system.

---

## Troubleshooting

### "PACKMOL failed"

Common causes:
- Box too small for all components
- Polymers too large
- Tolerance too tight

Solutions:
```yaml
solvent:
  box:
    padding: 2.0          # Increase padding
    tolerance: 2.5        # Increase tolerance (Angstrom)
```

### "Force field assignment failed"

The polymer structure may have issues:
- Check SDF files have correct bond orders
- For generated polymers, try a different polymer charge method:
  ```yaml
  polymers:
    charger: "am1bcc"    # More robust than NAGL, if AmberTools is available
  ```
- For `provided_molecules`, fix the input SDF. PolyzyMD requires complete finite
  partial charges and does not autocharge provided molecules.

### "Simulation unstable with polymers"

- Run longer equilibration:
  ```yaml
  simulation_phases:
    equilibration_stages:
      - name: "free_equilibration"
        duration: 2.0    # Increase from 1.0 ns
  ```
- Use softer restraints initially
- Check for clashes in initial structure

---

## Performance Considerations

| Configuration | System Size | Speed Impact |
|---------------|-------------|--------------|
| No polymers | Small | Fastest |
| 2 × 5-mer | Medium | ~10-20% slower |
| 4 × 10-mer | Large | ~30-50% slower |
| 10 × 20-mer | Very large | ~2-3× slower |

```{tip}
Start with small polymer systems (2 × 5-mer) to test your setup before scaling up.
```

---

## See Also

- {doc}`dynamic_polymers` - Generate native linear methacrylate polymers from SMILES
- {doc}`../reference/polymer_mbuild_openff` - mBuild/OpenFF adapter reference
- {doc}`gromacs_export` - Running simulations with GROMACS
- {doc}`../reference/configuration` - Complete configuration reference
