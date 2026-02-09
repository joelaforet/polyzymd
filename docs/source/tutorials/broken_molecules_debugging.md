# Help! My Simulation Shows a Ball of Bonds Instead of Reasonable Looking Molecules!

*A debugging case study from the PolyzyMD development team (January 2026)*

---

**You've waited hours (or days) for your MD simulation to finish. You excitedly load the trajectory in PyMOL or VMD, and... chaos. Bonds stretch across the entire simulation box. Molecules look like abstract art. Your beautiful enzyme-polymer system looks like someone threw spaghetti at the screen.**

If this sounds familiar, you're not alone. This document chronicles our journey debugging this exact issue in PolyzyMD, and we hope it saves future scientists the hours of head-scratching we went through.

---

## TL;DR (The Quick Fix)

**If you're in a hurry**, here's the essence:

- **Symptom:** Some molecules look broken (bonds span the box), but others look fine
- **Root Cause:** The atom order in your DCD trajectory doesn't match the atom order in your topology file (PDB)
- **Why it happens:** Any code that processes molecules by type (grouping waters together, ions together, etc.) can reorder atoms
- **Solution:** Ensure the code that creates your OpenMM system processes molecules in the SAME order as the topology file written for visualization

**Key diagnostic:** If only SOME molecules are broken (especially water/ions) while others (protein, ligand) look fine, you almost certainly have an atom order mismatch.

---

## The Problem We Encountered

### What We Saw

After running a production MD simulation with PolyzyMD, we loaded the trajectory in PyMOL:

```
PyMOL> load solvated_system.pdb
PyMOL> load_traj trajectory.dcd, solvated_system
```

Instead of a nice solvated enzyme system, we saw:

- **Bonds spanning the entire simulation box** - hundreds of angstroms long
- **Molecules that looked "exploded"** - atoms scattered everywhere
- **A tangled mess** that was completely unanalyzable

### The Crucial Observation

Here's what made this bug particularly insidious: **not everything was broken**.

- The **protein** looked perfect - proper secondary structure, reasonable motion
- The **substrate/ligand** looked fine - stayed near the active site
- The **polymers** looked correct - maintained their chain structure

But:

- **Water molecules** were completely scrambled
- **Ions** appeared in nonsensical positions
- The entire **solvent shell** was chaos

This selective breakage was the key clue that ultimately led us to the solution.

### Why This Matters

Broken trajectories aren't just ugly - they're scientifically dangerous:

1. **Analysis becomes meaningless** - RMSD, RMSF, distances are all wrong
2. **Subtle bugs might go unnoticed** - if you're only looking at the protein, you might miss the problem
3. **Wasted compute time** - hours/days of simulation that can't be used
4. **Potential for incorrect conclusions** - if you don't notice the problem

---

## Background: How DCD Trajectories Work

Before diving into our debugging journey, let's understand the fundamentals. This knowledge is essential for diagnosing trajectory problems.

### The DCD File Format

DCD (CHARMM/NAMD trajectory format) is one of the most common trajectory formats in MD. Here's what's inside:

```
DCD File Contents:
├── Header
│   ├── Number of atoms
│   ├── Number of frames
│   └── Timestep information
└── Frames (repeated for each saved frame)
    ├── Box vectors (if periodic)
    └── Coordinates
        ├── X coordinates for ALL atoms (atom 0, atom 1, ..., atom N-1)
        ├── Y coordinates for ALL atoms
        └── Z coordinates for ALL atoms
```

**What DCD files contain:**
- Atomic coordinates (X, Y, Z for each atom)
- Box dimensions (for periodic systems)
- Frame count and timing information

**What DCD files do NOT contain:**
- Atom names
- Residue names or numbers
- Chain identifiers
- Bond information
- Element types
- Charges
- **Any topology information whatsoever**

### The Critical Assumption

When you load a DCD file with a topology file (PDB, PSF, etc.), the visualization software makes a simple assumption:

> **Atom 0 in the DCD corresponds to Atom 0 in the topology.**
> **Atom 1 in the DCD corresponds to Atom 1 in the topology.**
> **And so on...**

There is no other way to map coordinates to atoms. The DCD format provides no atom identifiers - just coordinates indexed by position.

### A Simple Analogy

Imagine you have:
- A **guest list** (topology): "Seat 1: Alice, Seat 2: Bob, Seat 3: Carol..."
- A **seating chart update** (DCD frame): "Seat 1 is at position (x1,y1), Seat 2 is at (x2,y2)..."

If your seating chart was made with a different guest list order, you'd have Alice sitting where Bob should be, Bob where Carol should be, etc. The party would be chaos.

That's exactly what happens with mismatched atom order.

### Why DCD Works This Way

This design isn't a flaw - it's a feature:

1. **Efficiency:** No redundant topology data in every frame
2. **Speed:** Sequential coordinate arrays are fast to read/write
3. **Compatibility:** Works with any simulation package
4. **Size:** Trajectories stay small (coordinates only)

The tradeoff is that **you are responsible for ensuring atom order consistency**.

---

## Our Debugging Journey

Here's where we share our wrong turns, so you can avoid them.

### Attempt 1: "It Must Be PBC Wrapping!"

**Initial Hypothesis:** Periodic boundary conditions (PBC) are causing molecules to be split across box edges. The visualization software is drawing bonds across the box instead of wrapping them.

**Why this seemed reasonable:**
- PBC wrapping is a common source of visualization artifacts
- Molecules near box edges DO get split across boundaries
- Many trajectory post-processing discussions mention this

**What we tried:**

```python
# In the simulation runner, we tried various combinations:
state = simulation.context.getState(
    getPositions=True,
    enforcePeriodicBox=True   # Try to keep molecules whole
)

# Later tried:
enforcePeriodicBox=False  # Maybe this was the problem?
```

We also implemented custom unwrapping algorithms, trying to "make molecules whole" after the fact.

**Why it didn't work:**

PBC wrapping causes molecules to appear split at box boundaries, but the atoms themselves are still in the correct relative positions. If you zoom in on a "split" water molecule due to PBC, you'd see:
- One hydrogen on the left edge of the box
- The oxygen and other hydrogen on the right edge
- But the H-O distance would still be ~1 Å (correct)

Our problem was different: atoms were in completely wrong positions relative to each other. A water molecule might have one hydrogen near the protein and another hydrogen 50 Å away in the bulk solvent.

**The clue we missed:** If this were PBC wrapping, ALL molecules would be equally affected. But our protein looked perfect. Only the water/ions were scrambled.

### Attempt 2: "Maybe It's the DCD Reporter?"

**Second Hypothesis:** The OpenMM DCDReporter is writing coordinates incorrectly, perhaps with some buffering or ordering issue.

**What we tried:**
- Examining the reporter configuration
- Checking frame-by-frame coordinate outputs
- Comparing DCD output to PDB snapshots

**Why it didn't work:**

The DCDReporter was working correctly. When we extracted coordinates for a specific atom index from the DCD and from a PDB snapshot at the same timestep, they matched. The coordinates were being written correctly.

The problem wasn't that coordinates were wrong - it was that we were *interpreting* them with the wrong atom mapping.

### Attempt 3: "Is It the Topology File?"

**Third Hypothesis:** Maybe the PDB topology file we're using for visualization is corrupted or wrong.

**What we tried:**
- Regenerating the topology file
- Checking atom counts (they matched!)
- Visual inspection of the PDB

**What we found:**

The topology file was fine. Atom count matched the DCD perfectly. The PDB looked correct when loaded by itself.

**The key realization:** Both files were internally consistent. The problem was that they were consistent with DIFFERENT atom orderings.

### The Breakthrough: Pattern Recognition

We finally sat down and carefully observed which molecules were broken:

| Molecule Type | Appearance | Position in Old Processing Order |
|--------------|------------|----------------------------------|
| Protein | Perfect | 1st |
| Substrate | Perfect | 2nd |
| Polymers | Perfect | 3rd |
| Water | Broken | LAST (grouped) |
| Ions | Broken | LAST (grouped) |

**The pattern:** Everything that was processed early (in its "natural" order) looked fine. Everything that was grouped and processed at the end was broken.

This led us to examine our Interchange creation code...

---

## Root Cause: Atom Order Mismatch

### How OpenFF Interchange Works

[OpenFF Interchange](https://docs.openmm.org/latest/userguide/) is a powerful tool that converts OpenFF topology objects into OpenMM systems. The key operation for our story is `Interchange.combine()`:

```python
# Combining two interchanges appends atoms
combined = interchange_A.combine(interchange_B)

# Result atom order:
# [all atoms from A] + [all atoms from B]
```

The order of combination determines the final atom order.

### The Old Algorithm (Broken)

Our original `_create_interchange_batched()` method was designed for efficiency. It grouped molecules by type to minimize redundant parameterization:

```python
# OLD CODE (simplified) - THE PROBLEM

# Group molecules by SMILES (molecule type)
smiles_groups = defaultdict(list)
for molecule in solvated_topology.molecules:
    smiles = molecule.to_smiles()
    if smiles not in water_ion_smiles:
        smiles_groups[smiles].append(molecule)

# Create interchange for each molecule type
for smiles, molecules in smiles_groups.items():
    inc = ff.create_interchange(topology=molecules)
    all_interchanges.append(inc)

# Process water + ions LAST as a single batch
water_ion_mols = [m for m in topology if m.to_smiles() in water_ion_smiles]
water_ion_inc = ff.create_interchange(topology=water_ion_mols)
all_interchanges.append(water_ion_inc)  # <-- Always at the end!

# Combine all
combined = all_interchanges[0]
for inc in all_interchanges[1:]:
    combined = combined.combine(inc)
```

This produced an atom order like:
```
[protein atoms] [substrate atoms] [polymer atoms] [ALL water atoms] [ALL ion atoms]
```

### The Topology File Order

But `solvated_system.pdb` was written directly from the OpenFF topology BEFORE this reordering:

```python
# How the PDB was written
solvated_topology.to_file("solvated_system.pdb")
```

The OpenFF solvated topology had molecules in their original packing order:
```
[protein] [substrate] [polymer1] [water1] [Na+] [polymer2] [water2] [Cl-] [water3] ...
```

Water and ions were interspersed throughout, not grouped at the end!

### The Mismatch Visualized

```
PDB Topology (used for visualization):
Index:  0    1    2    3    4    5    6    7    8    9    ...
Atom:   N    CA   C    O    ...  OW   HW1  HW2  NA+  OW   ...
        └─protein─┘              └water┘     └ion┘ └water...

DCD Trajectory (after our processing):  
Index:  0    1    2    3    4    5    6    7    8    9    ...
Atom:   N    CA   C    O    ...  OW   OW   OW   OW   OW   ...
        └─protein─┘              └──ALL water grouped──┘

When PyMOL loads both:
- DCD atom 5 (a water oxygen) → mapped to PDB atom 5 (also water) ✓
- DCD atom 6 (a water oxygen) → mapped to PDB atom 6 (water hydrogen!) ✗
- DCD atom 7 (a water oxygen) → mapped to PDB atom 7 (water hydrogen!) ✗
- DCD atom 8 (a water oxygen) → mapped to PDB atom 8 (sodium ion!) ✗

Result: Complete chaos in the solvent
```

The protein matched because it was first in both orderings. The solvent was scrambled because it was reordered.

---

## The Solution

### New Algorithm: Preserve Order

The fix was conceptually simple: process molecules in their original order, only batching *consecutive* molecules of the same type:

```python
# NEW CODE (simplified) - THE FIX

current_smiles = None
current_batch = []

def flush_batch():
    """Create interchange for current batch."""
    if current_batch:
        inc = ff.create_interchange(topology=current_batch)
        all_interchanges.append(inc)
    current_batch.clear()

# Process molecules in ORIGINAL ORDER
for molecule in solvated_topology.molecules:  # Original iteration order!
    mol_smiles = molecule.to_smiles()
    
    if mol_smiles != current_smiles:
        # Different molecule type - flush previous batch
        flush_batch()
        current_smiles = mol_smiles
    
    current_batch.append(molecule)

flush_batch()  # Don't forget the last batch

# Combine all (now in correct order!)
combined = all_interchanges[0]
for inc in all_interchanges[1:]:
    combined = combined.combine(inc)
```

This preserves the exact molecule order while still batching consecutive same-type molecules for efficiency.

### Key Principle

> **The atom order in your trajectory MUST match the atom order in your topology file.**
> 
> Any code that manipulates molecular systems must preserve this order, or must regenerate the topology file to match the new order.

### Why Batching Consecutive Molecules Works

If the original order is:
```
[protein] [water1] [water2] [water3] [Na+] [water4] [water5] [Cl-]
```

The new algorithm produces batches:
```
Batch 1: [protein]           → Interchange 1
Batch 2: [water1,2,3]        → Interchange 2  (consecutive waters)
Batch 3: [Na+]               → Interchange 3
Batch 4: [water4,5]          → Interchange 4  (consecutive waters)
Batch 5: [Cl-]               → Interchange 5
```

Combined order: `[protein][water1,2,3][Na+][water4,5][Cl-]` - matches original!

---

## Lessons for Future Developers

### 1. DCD Files Are "Dumb"

Never forget that DCD files contain no topology information. They rely entirely on external files for atom identification. This means:

- The DCD doesn't know atom names, elements, or residues
- Coordinates are purely positional (atom 0, atom 1, ...)
- Any mismatch with the topology file causes silent corruption

### 2. Diagnosis Checklist

When you see broken molecules in a trajectory, ask yourself:

- [ ] **Are ALL molecules broken, or just some?**
  - All broken → probably PBC wrapping
  - Some broken → probably atom order mismatch

- [ ] **Which molecules are affected?**
  - Grouped at end of processing → likely reordering issue
  - Random distribution → look elsewhere

- [ ] **Do atom counts match?**
  ```python
  import MDAnalysis as mda
  u = mda.Universe('topology.pdb', 'trajectory.dcd')
  print(f"PDB atoms: {len(u.atoms)}")
  print(f"DCD atoms: {u.trajectory.n_atoms}")
  # These MUST match
  ```

- [ ] **Is there code that groups molecules by type?**
  - Dictionary grouping by SMILES/residue name
  - Separate processing of water/ions
  - Any sorting or reordering operations

### 3. Quick Verification Script

Use this script to sanity-check your trajectory:

```python
import MDAnalysis as mda

def verify_trajectory(topology_path, trajectory_path):
    """Verify trajectory is compatible with topology."""
    u = mda.Universe(topology_path, trajectory_path)
    
    # Check 1: Atom counts match
    topo_atoms = len(u.atoms)
    traj_atoms = u.trajectory.n_atoms
    print(f"Topology atoms: {topo_atoms}")
    print(f"Trajectory atoms: {traj_atoms}")
    
    if topo_atoms != traj_atoms:
        print("ERROR: Atom count mismatch!")
        return False
    
    # Check 2: Sample some atoms and their types
    print("\nFirst 10 atoms in topology:")
    for atom in u.atoms[:10]:
        print(f"  {atom.index}: {atom.name} ({atom.resname} {atom.resid})")
    
    # Check 3: Look at water molecules specifically
    waters = u.select_atoms("resname HOH or resname WAT or resname TIP3")
    if len(waters) > 0:
        print(f"\nFound {len(waters)} water atoms")
        print("First water molecule atoms:")
        first_water = u.select_atoms(f"resid {waters[0].resid}")
        for atom in first_water:
            print(f"  {atom.index}: {atom.name}")
    
    print("\nVisualize in PyMOL and check if molecules look whole!")
    return True

# Usage:
verify_trajectory('solvated_system.pdb', 'trajectory.dcd')
```

### 4. Prevention Strategies

1. **Test early:** Run a short simulation and check visualization BEFORE committing to long runs

2. **Document atom order assumptions:** If your code manipulates molecule order, document it clearly

3. **Write topology AFTER processing:** If you must reorder atoms, write the topology file from the same data structure that creates the simulation

4. **Use consistent iteration:** When iterating over molecules multiple times, ensure the order is deterministic

---

## How to Verify the Fix Worked

After applying a fix for atom order mismatch:

### Visual Verification

1. Load your files in PyMOL:
   ```
   load solvated_system.pdb
   load_traj trajectory.dcd, solvated_system
   ```

2. **Check the protein:** Should maintain secondary structure throughout

3. **Check water molecules:** Zoom in on individual waters
   - Should see proper H-O-H geometry
   - Bond lengths should be ~1 Å (not 50 Å!)
   
4. **Check ions:** Should appear as single atoms, not bonded to distant atoms

5. **Scrub through frames:** Problems often become more apparent with motion

### Programmatic Verification

```python
import MDAnalysis as mda

u = mda.Universe('solvated_system.pdb', 'trajectory.dcd')

# Check water geometry in first frame
waters = u.select_atoms("resname HOH or resname WAT")
for ts in u.trajectory[:1]:  # First frame
    # Get first water residue
    first_water_resid = waters[0].resid
    water_atoms = u.select_atoms(f"resid {first_water_resid}")
    
    if len(water_atoms) == 3:
        positions = water_atoms.positions
        # Calculate O-H distances
        oh1 = np.linalg.norm(positions[0] - positions[1])
        oh2 = np.linalg.norm(positions[0] - positions[2])
        
        print(f"Water O-H distances: {oh1:.2f} Å, {oh2:.2f} Å")
        
        if oh1 > 2.0 or oh2 > 2.0:
            print("WARNING: Water geometry looks wrong!")
        else:
            print("Water geometry looks correct!")
```

---

## Appendix: Technical Details

### Relevant Code Location

**File:** `src/polyzymd/builders/system_builder.py`

**Method:** `_create_interchange_batched()` (around line 500)

**Git commit:** `2047647` - "Fix atom order mismatch between solvated_system.pdb and Interchange"

### The Actual Fix (Abbreviated)

```python
def _create_interchange_batched(self, ff: ForceField, water_mol: Molecule) -> Interchange:
    """Create Interchange using batched molecule processing with preserved order.
    
    This implementation batches molecules by type for efficiency while preserving
    the exact molecule order from self._solvated_topology. This is critical for
    DCD trajectory compatibility - the atom order in the Interchange must match
    the atom order in solvated_system.pdb.
    """
    # ... setup code ...
    
    current_smiles = None
    current_batch: List[Molecule] = []

    def flush_batch():
        """Create interchange for current batch and add to list."""
        nonlocal current_batch, current_smiles
        if not current_batch:
            return
        # ... create interchange for batch ...
        current_batch = []

    # Process molecules in ORIGINAL ORDER
    for molecule in self._solvated_topology.molecules:
        mol_smiles = molecule.to_smiles()

        if mol_smiles != current_smiles:
            flush_batch()  # Different type - flush previous batch
            current_smiles = mol_smiles

        current_batch.append(molecule)

    flush_batch()  # Don't forget last batch
    
    # ... combine interchanges ...
```

### Related Resources

- [OpenFF Toolkit Documentation](https://docs.openforcefield.org/)
- [OpenMM User Guide](https://openmm.org/documentation)
- [MDAnalysis Documentation](https://docs.mdanalysis.org/)
- [DCD File Format Specification](https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html)

### When This Problem Might Appear Elsewhere

Watch out for atom order issues when:

1. **Using any molecule grouping/batching code** - efficiency optimizations often group by type
2. **Converting between file formats** - some formats sort atoms by element or residue
3. **Adding/removing atoms mid-simulation** - ensure all files are regenerated
4. **Combining systems from different sources** - verify consistent ordering
5. **Using topology manipulation tools** - some tools sort or reorder atoms

*Last updated: January 2026*

*Found this helpful? Found an error? Open an issue at https://github.com/joelaforet/polyzymd/issues*
