# Restraints Guide

This guide covers how to define and apply distance restraints in PolyzyMD simulations.

## Overview

Restraints are commonly used to:

- Keep a substrate near the enzyme active site
- Maintain specific protein-ligand contacts
- Prevent dissociation during equilibration

PolyzyMD supports several restraint types with flexible atom selection.

## Restraint Types

### Flat-Bottom Restraint

The most common restraint type. No force is applied when atoms are within the threshold distance; a harmonic force kicks in beyond it.

```
U(r) = 0                      if r < r₀
U(r) = ½k(r - r₀)²           if r ≥ r₀
```

```yaml
restraints:
  - type: "flat_bottom"
    name: "substrate_active_site"
    atom1:
      selection: "resid 77 and name OG"
    atom2:
      selection: "resname LIG and name C1"
    distance: 3.3        # Angstroms - threshold
    force_constant: 10000.0  # kJ/mol/nm²
    enabled: true
```

### Harmonic Restraint

Standard harmonic potential - always applies force to maintain target distance.

```
U(r) = ½k(r - r₀)²
```

```yaml
restraints:
  - type: "harmonic"
    name: "distance_restraint"
    atom1:
      selection: "resid 50 and name CA"
    atom2:
      selection: "resid 100 and name CA"
    distance: 10.0       # Angstroms - target distance
    force_constant: 1000.0
    enabled: true
```

### Upper Wall

Prevents distance from exceeding threshold (same as flat-bottom).

```yaml
- type: "upper_wall"
```

### Lower Wall

Prevents distance from going below threshold.

```
U(r) = ½k(r₀ - r)²           if r < r₀
U(r) = 0                      if r ≥ r₀
```

```yaml
- type: "lower_wall"
```

---

## Atom Selection Syntax

PolyzyMD uses an MDAnalysis-inspired selection syntax.

### Available Keywords

| Keyword | Indexing | Description | Example |
|---------|----------|-------------|---------|
| `resid` | 1-indexed | Residue number (matches PDB/PyMOL) | `resid 77` |
| `resname` | - | Residue name | `resname SER` |
| `name` | - | Atom name | `name OG` |
| `index` | 0-indexed | OpenMM atom index | `index 2739` |
| `pdbindex` | 1-indexed | PDB atom serial (matches PyMOL) | `pdbindex 2740` |
| `chain` | - | Chain ID | `chain A` |

### Combining Selections

Use `and` / `or` to combine:

```yaml
# Both conditions must be true
selection: "resid 77 and name OG"

# Either condition
selection: "resname SER or resname THR"

# Complex combinations
selection: "(resid 77 and name OG) or (resid 110 and name NE2)"
```

### Index Conventions

```{important}
**Understanding indexing is critical for correct restraints!**
```

| What you see | Keyword to use | Value to use |
|--------------|----------------|--------------|
| PyMOL shows residue 77 | `resid` | 77 |
| PyMOL shows atom serial 2740 | `pdbindex` | 2740 |
| OpenMM internal index 2739 | `index` | 2739 |

**Rule of thumb:** Use `resid` and `pdbindex` - they match what you see in PyMOL.

---

## Workflow: Finding Atoms for Restraints

### Step 1: Open PDB in PyMOL

```bash
pymol structures/enzyme.pdb
```

### Step 2: Identify Target Atoms

In PyMOL, click on atoms to see their info:

```
You clicked /enzyme/A/SER`77/OG
 ID = 1192
```

This tells you:
- Residue: SER 77 (chain A)
- Atom name: OG
- Atom serial: would show in sequence window

### Step 3: Write Selection

```yaml
atom1:
  selection: "resid 77 and name OG"
  description: "Ser77 OG - catalytic nucleophile"
```

### Step 4: For Specific Atom Index

If you need a specific atom (e.g., when multiple atoms match):

1. In PyMOL, check the atom serial number
2. Use `pdbindex` with that number:

```yaml
atom1:
  selection: "resid 77 and pdbindex 1193"
  description: "Specific Ser77 atom"
```

---

## Common Restraint Patterns

### Substrate in Active Site

Keep substrate within bonding distance of catalytic residue:

```yaml
restraints:
  - type: "flat_bottom"
    name: "substrate_catalytic"
    atom1:
      selection: "resid 77 and name OG"
      description: "Catalytic serine"
    atom2:
      selection: "resname LIG and name C1"
      description: "Substrate carbonyl carbon"
    distance: 3.5
    force_constant: 10000.0
    enabled: true
```

### Multiple Active Site Contacts

```yaml
restraints:
  - type: "flat_bottom"
    name: "contact_1"
    atom1:
      selection: "resid 77 and name OG"
    atom2:
      selection: "resname LIG and name C1"
    distance: 3.5
    force_constant: 10000.0
    enabled: true
    
  - type: "flat_bottom"
    name: "contact_2"
    atom1:
      selection: "resid 110 and name NE2"
    atom2:
      selection: "resname LIG and name O2"
    distance: 3.5
    force_constant: 5000.0
    enabled: true
    
  - type: "flat_bottom"
    name: "contact_3"
    atom1:
      selection: "resid 12 and name N"
    atom2:
      selection: "resname LIG and name O1"
    distance: 4.0
    force_constant: 5000.0
    enabled: false  # Disabled for now
```

### Protein-Protein Distance

Monitor or restrain distance between protein regions:

```yaml
restraints:
  - type: "harmonic"
    name: "domain_distance"
    atom1:
      selection: "resid 50 and name CA"
    atom2:
      selection: "resid 150 and name CA"
    distance: 25.0
    force_constant: 100.0  # Weak restraint
    enabled: true
```

---

## Force Constant Guidelines

| Use Case | Force Constant (kJ/mol/nm²) |
|----------|---------------------------|
| Strong restraint (prevent dissociation) | 10000 - 50000 |
| Moderate restraint | 1000 - 5000 |
| Weak/guiding restraint | 100 - 500 |
| Very weak (biasing) | 10 - 50 |

```{note}
Start with stronger restraints during equilibration, then reduce or disable for production if desired.
```

---

## Troubleshooting

### "No atoms match selection"

```
ValueError: No atoms match selection: 'resid 77 and name OG'
```

**Causes:**
1. Residue numbering doesn't match (check PDB)
2. Atom name is different (check PDB atom names)
3. Typo in selection string

**Solution:** Open PDB in PyMOL and verify exact residue numbers and atom names.

### "Requires exactly one atom per selection"

```
ValueError: Restraint 'my_restraint' requires exactly one atom per selection.
Got 5 for atom1, 1 for atom2
```

**Cause:** Selection matches multiple atoms.

**Solution:** Make selection more specific:

```yaml
# Too broad - matches all atoms in residue
selection: "resid 77"

# Better - matches specific atom
selection: "resid 77 and name OG"

# Most specific - use atom index
selection: "resid 77 and pdbindex 1193"
```

### Restraint Not Taking Effect

1. Check `enabled: true` is set
2. Verify distance units (should be Angstroms)
3. Check force constant is reasonable (not too small)
4. Verify atom selections resolve correctly

---

## Disabling Restraints

### In Configuration

```yaml
restraints:
  - type: "flat_bottom"
    name: "my_restraint"
    # ... other settings ...
    enabled: false  # Disabled
```

### No Restraints at All

```yaml
restraints: []
```

Or simply omit the `restraints` section entirely.
