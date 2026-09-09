# Polymer Setup Guide

This guide covers configuring polymer chains in PolyzyMD simulations.

```{tip}
**Looking for dynamic polymer generation?** If you want to generate polymers on-the-fly from SMILES strings (without pre-built SDF files), see the {doc}`dynamic_polymers` tutorial.
```

## Overview

PolyzyMD supports adding random co-polymer chains to your simulation box. Polymers are:

- Generated based on monomer probabilities
- Placed around the enzyme using PACKMOL
- Parameterized with OpenFF force fields

## Basic Configuration

```yaml
polymers:
  enabled: true
  type_prefix: "SBMA-EGPMA"
  
  monomers:
    - label: "A"
      probability: 0.98
      name: "SBMA"
    - label: "B"
      probability: 0.02
      name: "EGPMA"
  
  length: 5      # 5-mer chains
  count: 2       # 2 polymer chains
```

---

## Monomer Definition

### Probability-Based Selection

Each chain is built by randomly selecting monomers based on their probabilities:

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

```yaml
length: 5    # 5-mer (pentamer)
```

Typical values:
- **Short chains**: 3-5 monomers (faster simulations)
- **Medium chains**: 10-20 monomers
- **Long chains**: 50+ monomers (slower, more realistic)

### Number of Chains

```yaml
count: 2    # Add 2 polymer chains
```

More chains = larger system = slower simulation.

---

## Pre-Built Polymer SDFs

For reproducibility, you can provide pre-built polymer structures instead of random generation.

### Directory Structure

```
polymer_sdfs/
└── SBMA-EGPMA/
    ├── AAAAA.sdf     # All A monomers
    ├── AAAAB.sdf     # 4 A's, 1 B
    ├── AAABA.sdf
    ├── AABAA.sdf
    └── ...
```

### Configuration

```yaml
polymers:
  enabled: true
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
  sdf_directory: "polymer_sdfs/SBMA-EGPMA"    # Path to pre-built SDFs
```

### Naming Convention

SDF files must be named with the monomer sequence:
- `AAAAA.sdf` - Sequence of 5 "A" monomers
- `AABBA.sdf` - Sequence A-A-B-B-A
- Labels must match those defined in `monomers`

---

## Polymer Cache

Generated polymers are cached for reuse:

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

### SBMA-EGPMA Co-polymer

Zwitterionic sulfobetaine with hydrophobic groups:

```yaml
polymers:
  enabled: true
  type_prefix: "SBMA-EGPMA"
  monomers:
    - label: "A"
      probability: 0.98
      name: "SBMA"          # Sulfobetaine methacrylate
    - label: "B"
      probability: 0.02
      name: "EGPMA"         # Ethylene glycol phenyl ether methacrylate
  length: 5
  count: 2
```

### PEG Homopolymer

Polyethylene glycol:

```yaml
polymers:
  enabled: true
  type_prefix: "PEG"
  monomers:
    - label: "A"
      probability: 1.0
      name: "EthyleneGlycol"
  length: 10
  count: 4
```

### Block Co-polymer (Approximate)

For a block-like structure, use pre-built SDFs:

```yaml
polymers:
  enabled: true
  type_prefix: "Block-AB"
  monomers:
    - label: "A"
      probability: 0.5
      name: "BlockA"
    - label: "B"
      probability: 0.5
      name: "BlockB"
  length: 10
  count: 2
  sdf_directory: "polymer_sdfs/block_copolymer"   # Pre-built block structures
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

1. Enzyme (+ substrate) placed at the center and held fixed
2. Polymers placed anywhere in the packing box; PACKMOL keeps every polymer
   atom at least `polymers.packing.tolerance` (default 2 Å) from the fixed
   solute and from other chains
3. Water molecules fill the remaining space (solvation box from
   `solvent.box.padding`)
4. Ions added to neutralize and reach target concentration

PACKMOL is seeded with the replicate number, so each replicate gets an
independent polymer arrangement and solvent configuration.

### Packing Box and Shell

`polymers.packing.padding` (default 2.0 nm) sets how far the polymer packing box
extends beyond the solute:

```yaml
polymers:
  packing:
    padding: 2.0          # nm around the solute for polymer placement
    tolerance: 2.0        # Å minimum atom-atom distance
    movebadrandom: true   # helps many unique chain types converge
```

Earlier versions also confined chains to a rectangular shell outside the
solute's bounding box. That shell was often thinner than the chains themselves
and made PACKMOL run to its loop limit without converging. It is now off by
default; set `packing.exclude_solute_bbox: true` only if you need the legacy
behaviour.

### After Packing

PACKMOL may exit with code 173 ("imperfect packing") for dense systems. PolyzyMD
accepts the best solution, and energy minimization resolves the residual
contacts. Two safeguards apply before any simulation starts:

- The build aborts with `SolvationClashError` if more than a handful of packed
  atoms overlap the solute, which indicates a coordinate-frame problem rather
  than imperfect packing.
- Minimization holds the protein and substrate fixed by default, so residual
  polymer or water contacts are resolved by moving the polymer or water, never
  the protein. See {doc}`../explanation/simulation_safeguards`.

---

## Troubleshooting

### "PACKMOL failed"

Common causes:
- Box too small for all components
- Polymers too large
- Tolerance too tight

Solutions:
```yaml
polymers:
  packing:
    padding: 2.5          # Larger polymer packing box
    movebadrandom: true   # Better convergence for many chain types
    nloop: 200            # Raise if PACKMOL stops at the loop limit
solvent:
  box:
    padding: 2.0          # Larger solvation box
```

If PACKMOL exits 173 (imperfect packing) the build still continues; the
residual contacts are removed during minimization. A `SolvationClashError`
means something else is wrong (see {doc}`troubleshooting`).

### "Force field assignment failed"

The polymer structure may have issues:
- Check SDF files have correct bond orders
- Try different charge method:
  ```yaml
  substrate:
    charge_method: "am1bcc"    # More robust than NAGL
  ```

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

- {doc}`dynamic_polymers` - Generate polymers from SMILES without pre-built SDF files
- {doc}`gromacs_export` - Running simulations with GROMACS
- {doc}`../reference/configuration` - Complete configuration reference
