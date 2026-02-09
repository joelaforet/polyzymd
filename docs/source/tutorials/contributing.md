# Contributing Guide

This guide covers how to contribute to PolyzyMD, with a focus on adding new solvents to the library.

## Contributing New Solvents

The co-solvent library is designed to be easily extensible. Here's how to add a new solvent.

### Step 1: Add to the Co-solvent Library

Edit `src/polyzymd/data/cosolvent_library.py` and add your solvent to `COSOLVENT_LIBRARY`:

```python
COSOLVENT_LIBRARY: Dict[str, CoSolventData] = {
    # ... existing solvents ...
    
    "formamide": CoSolventData(
        name="formamide",
        smiles="C(=O)N",
        density=1.133,           # g/mL from PubChem
        residue_name="FOR",      # 3-letter code
        pubchem_cid=713,         # For reference
    ),
}
```

**Requirements:**
- **name**: Lowercase identifier (used in YAML configs)
- **smiles**: Valid SMILES string (verify with RDKit or OpenEye)
- **density**: Density in g/mL (required for volume fraction calculations)
- **residue_name**: 3-letter code for PDB/topology files
- **pubchem_cid**: PubChem Compound ID for traceability

### Step 2: Generate the Pre-computed SDF

Run the generator script to create the SDF file with embedded charges:

```bash
conda activate polymerist-env
cd /path/to/polyzymd
python src/polyzymd/data/solvents/_generator.py
```

This will:
1. Load the molecule from SMILES
2. Generate a 3D conformer
3. Compute AM1BCC partial charges
4. Save to `src/polyzymd/data/solvents/{name}.sdf`

### Step 3: Verify the SDF

Check that the SDF was created correctly:

```python
from openff.toolkit import Molecule

mol = Molecule.from_file("src/polyzymd/data/solvents/formamide.sdf")
print(f"Atoms: {mol.n_atoms}")
print(f"Charges: {mol.partial_charges}")
print(f"Residue: {mol.atoms[0].metadata.get('residue_name')}")
```

### Step 4: Test the New Solvent

Create a test configuration using your new solvent:

```yaml
solvent:
  primary:
    type: "water"
    model: "tip3p"
  co_solvents:
    - name: "formamide"
      volume_fraction: 0.10
```

Run the build to verify it works:

```bash
polyzymd build test_config.yaml --dry-run
```

### Step 5: Submit a Pull Request

1. Commit your changes:
   - `src/polyzymd/data/cosolvent_library.py` (library entry)
   - `src/polyzymd/data/solvents/{name}.sdf` (pre-computed charges)

2. Include in your PR description:
   - The solvent name and use case
   - PubChem CID reference for density value
   - Any special considerations (e.g., hydrogen bonding, coordination)

## SDF File Format

The SDF files in `src/polyzymd/data/solvents/` follow a specific format:

```
molecule_name
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    ...
M  CHG  ...
M  END
> <atom.dprop.PartialCharge>
-0.123456 0.234567 ...

> <residue_name>
DMS

$$$$
```

Key elements:
- **3D coordinates**: Generated conformer
- **PartialCharge property**: AM1BCC charges for each atom
- **residue_name property**: 3-letter residue code for topology

## Water Models

Water molecules use hardcoded literature charges (not AM1BCC) for accuracy:

| Model | O Charge | H Charge | Source |
|-------|----------|----------|--------|
| TIP3P | -0.834 | +0.417 | Jorgensen et al., 1983 |
| SPC/E | -0.8476 | +0.4238 | Berendsen et al., 1987 |

To add a new water model, edit `src/polyzymd/data/solvent_molecules.py`:

```python
def _create_tip4pew_water() -> Molecule:
    """Create TIP4P-Ew water with literature charges."""
    mol = Molecule.from_smiles("O")
    mol.generate_conformers(n_conformers=1)
    
    # TIP4P-Ew charges (note: virtual site not included)
    charges = [-0.84844, 0.42422, 0.42422] * unit.elementary_charge
    mol.partial_charges = charges
    
    for atom in mol.atoms:
        atom.metadata["residue_name"] = "HOH"
    
    return mol
```

## Code Style

- Follow PEP 8 and use `ruff` for linting
- Add type hints for all function signatures
- Include docstrings with Args/Returns sections
- Keep lines under 100 characters

## Testing

Before submitting:

```bash
# Run linting
ruff check src/

# Run type checking (if available)
mypy src/polyzymd/

# Test the build process
polyzymd build examples/enzyme_cosolvent.yaml --dry-run
```

## Questions?

Open an issue on GitHub: https://github.com/joelaforet/polyzymd/issues
