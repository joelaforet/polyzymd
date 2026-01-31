# Understanding Residue Assignment in PolyzyMD

*How PolyzyMD assigns unique residue identifiers to enable powerful trajectory analysis*

---

## The Problem: "One Molecule, One Residue"

When you build an MD system with most tools (including vanilla OpenFF/OpenMM workflows), you often encounter a frustrating limitation: **all water molecules end up with the same residue number**.

### What Happens by Default

In a typical OpenFF-based workflow:

```
# Default OpenFF behavior for water:
HOH  1   O     ...
HOH  1   H1    ...
HOH  1   H2    ...
HOH  1   O     ...   # Same residue number!
HOH  1   H1    ...
HOH  1   H2    ...
... (thousands more, all residue 1)
```

Every water molecule is assigned residue number 1. The same applies to ions (all Na+ are residue 1, all Cl- are residue 1) and often to co-solvents as well.

### Why This Is a Problem

This "one molecule type, one residue" approach loses critical information:

1. **You can't select individual molecules**
   ```python
   # This doesn't work - selects ALL water!
   water_near_protein = u.select_atoms("resname HOH and around 5 protein")
   
   # You want: "Give me the 10 water molecules closest to the active site"
   # But you can't distinguish water molecule #1 from #5000
   ```

2. **Distance calculations are ambiguous**
   ```python
   # Which water molecule is this distance for?
   distance = atom1.position - water.position  # Which water?!
   ```

3. **Visualization suffers**
   - PyMOL/VMD can't color individual water molecules
   - Can't track specific solvent molecules over time
   - Can't identify hydration shell residents

4. **Standard analysis tools break**
   - Radial distribution functions need molecule identifiers
   - Hydrogen bond analysis needs to track donor/acceptor molecules
   - Residence time calculations are impossible

---

## PolyzyMD's Solution: One Molecule = One Residue

PolyzyMD assigns **unique residue numbers to every molecule** in your system. This is a deliberate design choice that enables powerful downstream analysis.

### What PolyzyMD Produces

```
# PolyzyMD output for water:
HOH  1   O     ...   # Water molecule 1
HOH  1   H1    ...
HOH  1   H2    ...
HOH  2   O     ...   # Water molecule 2 (unique!)
HOH  2   H1    ...
HOH  2   H2    ...
HOH  3   O     ...   # Water molecule 3 (unique!)
...
```

Every molecule gets a unique (chain_id, residue_number) combination, enabling precise selection and analysis.

### Chain Assignment Strategy

PolyzyMD assigns chain IDs systematically:

| Component | Chain | Residue Numbers |
|-----------|-------|-----------------|
| Protein | A | Preserved from input PDB (e.g., 1-289) |
| Substrate | B | 1 |
| Polymers | C | Sequential per monomer (1, 2, 3, ...) |
| Solvent | D+ | Sequential per molecule (1-9999, then overflow to E, F, ...) |

**Overflow handling:** PDB format limits residue numbers to 9999. When exceeded, PolyzyMD automatically moves to the next chain letter:

```
Chain D: residues 1-9999  (waters 1-9999)
Chain E: residues 1-9999  (waters 10000-19999)
Chain F: residues 1-...   (and so on)
```

---

## Practical Analysis Examples

Now that every molecule has a unique identifier, here's what you can do:

### 1. Select Specific Solvent Molecules

```python
import MDAnalysis as mda

u = mda.Universe('solvated_system.pdb', 'trajectory.dcd')

# Select water molecules 1-10
first_ten_waters = u.select_atoms("resname HOH and resid 1:10")

# Select a specific ion
sodium_42 = u.select_atoms("resname NA and resid 42")

# Select water in chain D only (first 9999 waters)
chain_d_water = u.select_atoms("chainID D and resname HOH")
```

### 2. Find Waters Near the Active Site

```python
# Select waters within 5 Å of the catalytic serine
active_site = u.select_atoms("resid 77 and name OG")
nearby_waters = u.select_atoms("resname HOH and around 5.0 group active", active=active_site)

# Get unique water molecule IDs
water_resids = set(nearby_waters.resids)
print(f"Found {len(water_resids)} unique water molecules near active site")
print(f"Water residue IDs: {sorted(water_resids)[:10]}...")  # First 10
```

### 3. Track Hydration Shell Over Time

```python
# Track which waters are in the first hydration shell over the trajectory
hydration_residents = []

for ts in u.trajectory:
    shell = u.select_atoms("resname HOH and around 3.5 protein")
    # Get unique molecule IDs (not just atoms)
    molecule_ids = set(zip(shell.chainIDs, shell.resids))
    hydration_residents.append(molecule_ids)

# Find waters that stayed in the shell the entire time
persistent_waters = set.intersection(*hydration_residents)
print(f"{len(persistent_waters)} waters remained in hydration shell throughout")
```

### 4. Calculate Water Residence Times

```python
import numpy as np

# Define region of interest
active_site = u.select_atoms("resid 77 and name OG")

# Track presence of each water molecule
water_ids = u.select_atoms("resname HOH and name O").resids  # All water oxygens
residence_matrix = np.zeros((len(u.trajectory), max(water_ids) + 1), dtype=bool)

for i, ts in enumerate(u.trajectory):
    nearby = u.select_atoms("resname HOH and name O and around 5.0 group site", site=active_site)
    for resid in nearby.resids:
        residence_matrix[i, resid] = True

# Calculate residence time for each water
residence_times = residence_matrix.sum(axis=0)  # Frames present
top_residents = np.argsort(residence_times)[-10:][::-1]
print("Top 10 longest-residing waters near active site:")
for resid in top_residents:
    frames = residence_times[resid]
    print(f"  Water {resid}: {frames} frames ({frames * 0.1:.1f} ns)")  # assuming 100 ps/frame
```

### 5. Hydrogen Bond Analysis with Molecule Identity

```python
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis

# Analyze H-bonds between protein and specific water molecules
hbonds = HydrogenBondAnalysis(
    universe=u,
    donors_sel="protein and name N",
    hydrogens_sel="protein and name H", 
    acceptors_sel="resname HOH and name O",
)
hbonds.run()

# Get H-bond counts per water molecule
from collections import Counter
acceptor_resids = [hb[4] for hb in hbonds.results.hbonds]  # Column 4 is acceptor resid
water_hbond_counts = Counter(acceptor_resids)

print("Waters with most H-bonds to protein:")
for resid, count in water_hbond_counts.most_common(10):
    print(f"  Water {resid}: {count} H-bonds")
```

### 6. PyMOL Visualization

With unique residue numbers, you can do precise selections in PyMOL:

```python
# In PyMOL
# Color the 5 water molecules closest to the ligand
select nearby_water, (resn HOH within 4 of resn LIG)
# This now gives you individual molecules you can work with!

# Color specific waters by their behavior
# (after identifying "interesting" waters from your analysis)
select persistent_water, resn HOH and resi 1234+5678+9012
color cyan, persistent_water

# Show waters that bridge protein-ligand
select bridging, resn HOH and (within 3.5 of protein) and (within 3.5 of resn LIG)
show sticks, bridging
```

---

## Technical Details

### How It Works

PolyzyMD's `SystemBuilder._assign_pdb_identifiers()` method handles residue assignment:

1. **Protein:** Chain A, preserves original residue numbers from input PDB
2. **Substrate:** Chain B, residue 1
3. **Polymers:** Chain C, sequential numbering per monomer unit
4. **Solvent:** Starting from chain D, sequential numbering with overflow handling

```python
# Simplified logic for solvent assignment
residue_num = 1
chain_idx = 3  # Start at 'D'

for molecule in solvent_molecules:
    if residue_num > 9999:
        chain_idx += 1
        residue_num = 1
    
    for atom in molecule.atoms:
        atom.metadata["chain_id"] = CHAIN_LETTERS[chain_idx]
        atom.metadata["residue_number"] = str(residue_num)
    
    residue_num += 1
```

### Why Other Tools Don't Do This

The "one molecule type, one residue" approach is common because:

1. **Historical:** Traditional force field tools (AMBER, GROMACS) handled this differently
2. **Simplicity:** Easier to implement molecule type grouping
3. **File size:** Fewer unique residue IDs = smaller topology files
4. **Legacy compatibility:** Some analysis tools expect grouped residues

PolyzyMD prioritizes **analysis capability** over file size, recognizing that modern storage is cheap but lost information is irreplaceable.

### Verifying Your Topology

You can verify that residue assignment worked correctly:

```python
import MDAnalysis as mda

u = mda.Universe('solvated_system.pdb')

# Check water molecules have unique residues
waters = u.select_atoms("resname HOH")
water_molecules = waters.residues

print(f"Total water atoms: {len(waters)}")
print(f"Unique water residues: {len(water_molecules)}")
print(f"Expected (atoms/3): {len(waters) // 3}")

# Should match!
assert len(water_molecules) == len(waters) // 3, "Residue assignment issue!"

# Check residue number range
print(f"Water residue IDs range: {min(waters.resids)} to {max(waters.resids)}")

# Check chain distribution
from collections import Counter
chain_counts = Counter(waters.chainIDs)
print(f"Waters per chain: {dict(chain_counts)}")
```

---

## Comparison: Before and After

### Without Unique Residues (Default OpenFF)

```python
# Trying to find waters near active site
waters = u.select_atoms("resname HOH")
print(f"Found {len(waters)} water atoms")  # 30,000 atoms

# All have the same residue ID!
print(f"Unique resids: {set(waters.resids)}")  # {1}

# Can't distinguish molecules
# Can't track over time
# Can't calculate per-molecule properties
```

### With Unique Residues (PolyzyMD)

```python
# Find waters near active site
nearby = u.select_atoms("resname HOH and around 5 (resid 77 and name OG)")
molecule_ids = set(nearby.resids)
print(f"Found {len(molecule_ids)} unique water molecules")  # e.g., 23

# Track specific molecules
for ts in u.trajectory[::100]:  # Every 100th frame
    still_there = u.select_atoms(f"resid {' '.join(map(str, molecule_ids))} and around 5 (resid 77)")
    print(f"Frame {ts.frame}: {len(set(still_there.resids))} of {len(molecule_ids)} still nearby")
```

---

## Best Practices

### 1. Always Use the PolyzyMD-Generated PDB

The residue assignments are stored in `solvated_system.pdb`. Always use this file (not the input PDB) for analysis:

```python
# Correct
u = mda.Universe('output/solvated_system.pdb', 'output/production.dcd')

# Wrong - loses residue information!
u = mda.Universe('input/protein.pdb', 'output/production.dcd')
```

### 2. Mind the Chain Boundaries

If your system has >9999 solvent molecules, water spans multiple chains:

```python
# Select ALL water (across all chains)
all_water = u.select_atoms("resname HOH")

# Select water in first solvent chain only
chain_d_water = u.select_atoms("chainID D and resname HOH")
```

### 3. Use Molecule-Based Selections

When analyzing solvent, think in terms of molecules, not atoms:

```python
# Per-molecule analysis
for residue in u.select_atoms("resname HOH").residues:
    # residue.atoms gives you O, H1, H2 for this water
    com = residue.atoms.center_of_mass()
    # Now you can track this specific molecule
```

---

## Summary

PolyzyMD's "one molecule = one residue" approach enables:

- Selecting individual solvent molecules
- Tracking molecules over trajectory time
- Calculating per-molecule properties
- Proper hydrogen bond analysis
- Hydration shell characterization
- Residence time calculations
- Meaningful visualization in PyMOL/VMD

This is a deliberate design choice that prioritizes analysis capability. The topology file is slightly larger, but the analytical power gained is substantial.

---

*Last updated: January 2026*
