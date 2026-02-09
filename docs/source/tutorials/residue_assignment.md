# Understanding Residue Assignment in PolyzyMD

*How PolyzyMD assigns unique residue identifiers to enable powerful trajectory analysis*

---

## The Core Principle: One Repeat Unit = One Residue

In molecular dynamics, the concept of a "residue" comes from protein chemistry, where each amino acid is a distinct residue in the polymer chain. PolyzyMD extends this principle consistently across all molecule types:

> **One repeat unit = One residue**

This principle applies uniformly:

| Molecule Type | Repeat Unit | Residue Assignment |
|---------------|-------------|-------------------|
| **Protein** | Amino acid | Each amino acid is one residue (GLY, ALA, SER, ...) |
| **Polymer** | Monomer | Each monomer unit is one residue |
| **Solvent** (water, ions) | Entire molecule | Each molecule is one residue (the "repeat unit" is the molecule itself) |
| **Co-solvent** (DMSO, etc.) | Entire molecule | Each molecule is one residue |

For proteins and synthetic polymers, the repeat unit is the monomer. For small molecules like water or DMSO, the molecule itself is the fundamental unit—there's nothing smaller to repeat—so each molecule gets its own residue number.

---

## The Problem: Default Behavior Loses Information

When you build an MD system with most tools (including vanilla OpenFF/OpenMM workflows), you often encounter a frustrating limitation: **all molecules of the same type share a single residue number**.

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

This "one molecule **type**, one residue" approach conflates distinct molecules, losing critical information.

### Why This Is a Problem

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

## PolyzyMD's Solution: Consistent Repeat-Unit Residue Assignment

PolyzyMD assigns residue numbers based on the repeat-unit principle, ensuring every chemically meaningful unit has a unique identifier.

### What PolyzyMD Produces

**For proteins** (preserved from input):
```
# Each amino acid is its own residue (standard behavior)
ALA  1   N     ...
ALA  1   CA    ...
ALA  1   C     ...
GLY  2   N     ...   # Next residue
GLY  2   CA    ...
...
```

**For polymers** (one residue per monomer):
```
# Each monomer unit gets a unique residue number
SBMA 1   C1    ...   # First monomer
SBMA 1   C2    ...
SBMA 1   S     ...
SBMA 2   C1    ...   # Second monomer (unique!)
SBMA 2   C2    ...
...
```

**For solvent** (one residue per molecule):
```
# Each water molecule gets a unique residue number
HOH  1   O     ...   # Water molecule 1
HOH  1   H1    ...
HOH  1   H2    ...
HOH  2   O     ...   # Water molecule 2 (unique!)
HOH  2   H1    ...
HOH  2   H2    ...
HOH  3   O     ...   # Water molecule 3 (unique!)
...
```

Every repeat unit gets a unique (chain_id, residue_number) combination, enabling precise selection and analysis.

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
Chain D: residues 1-9999  (solvent molecules 1-9999)
Chain E: residues 1-9999  (solvent molecules 10000-19999)
Chain F: residues 1-...   (and so on)
```

---

## Practical Analysis Examples

Now that every repeat unit has a unique identifier, here's what you can do:

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

### 5. Analyze Polymer Monomers Individually

```python
# Select specific monomers in a polymer chain
# (Remember: each monomer is its own residue)
first_monomer = u.select_atoms("chainID C and resid 1")
last_monomer = u.select_atoms("chainID C and resid 10")  # For a 10-mer

# Calculate end-to-end distance over trajectory
e2e_distances = []
for ts in u.trajectory:
    com1 = first_monomer.center_of_mass()
    com2 = last_monomer.center_of_mass()
    e2e_distances.append(np.linalg.norm(com2 - com1))

print(f"End-to-end distance: {np.mean(e2e_distances):.2f} ± {np.std(e2e_distances):.2f} Å")
```

### 6. Hydrogen Bond Analysis with Molecule Identity

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

### 7. PyMOL Visualization

With unique residue numbers, you can do precise selections in PyMOL:

```
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

# Select specific monomers in polymer chain C
select monomer_5, chain C and resi 5
```

---

## Technical Details

### How It Works

PolyzyMD's `SystemBuilder._assign_pdb_identifiers()` method handles residue assignment:

1. **Protein:** Chain A, preserves original residue numbers from input PDB (amino acids already have correct residue assignments)
2. **Substrate:** Chain B, residue 1 (small molecule ligand is a single "unit")
3. **Polymers:** Chain C, sequential numbering per monomer unit (each monomer = one residue)
4. **Solvent:** Starting from chain D, sequential numbering per molecule (each molecule = one residue)

```python
# Simplified logic for solvent assignment
residue_num = 1
chain_idx = 3  # Start at 'D'

for molecule in solvent_molecules:
    if residue_num > 9999:
        chain_idx += 1
        residue_num = 1
    
    # Each molecule (the repeat unit for solvent) gets one residue number
    for atom in molecule.atoms:
        atom.metadata["chain_id"] = CHAIN_LETTERS[chain_idx]
        atom.metadata["residue_number"] = str(residue_num)
    
    residue_num += 1
```

### Why Other Tools Don't Do This

The "one molecule **type**, one residue" approach is common because:

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

### Without Proper Residue Assignment (Default OpenFF)

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

### With PolyzyMD's Repeat-Unit Residue Assignment

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

If your system has >9999 solvent molecules, they span multiple chains:

```python
# Select ALL water (across all chains)
all_water = u.select_atoms("resname HOH")

# Select water in first solvent chain only
chain_d_water = u.select_atoms("chainID D and resname HOH")
```

### 3. Think in Terms of Repeat Units

When analyzing your system, think about what the fundamental repeat unit is:

```python
# For proteins: iterate over amino acids (residues)
for residue in u.select_atoms("protein").residues:
    print(f"{residue.resname} {residue.resid}: {len(residue.atoms)} atoms")

# For polymers: iterate over monomers (residues)
for residue in u.select_atoms("chainID C").residues:
    print(f"Monomer {residue.resid}: {residue.resname}")

# For solvent: iterate over molecules (each is one residue)
for residue in u.select_atoms("resname HOH").residues:
    com = residue.atoms.center_of_mass()
    # Now you can track this specific water molecule
```

---

## Summary

PolyzyMD's **"one repeat unit = one residue"** principle provides a consistent, chemically meaningful approach to residue assignment:

- **Proteins:** Each amino acid is one residue (standard)
- **Polymers:** Each monomer is one residue (enables per-monomer analysis)
- **Solvent/Co-solvent:** Each molecule is one residue (enables molecule tracking)

This enables powerful analysis workflows:

- Selecting individual solvent molecules
- Tracking molecules over trajectory time
- Calculating per-molecule or per-monomer properties
- Proper hydrogen bond analysis
- Hydration shell characterization
- Residence time calculations
- Meaningful visualization in PyMOL/VMD

This is a deliberate design choice that prioritizes analysis capability. The topology file is slightly larger, but the analytical power gained is substantial.

---

*Last updated: January 2026*
