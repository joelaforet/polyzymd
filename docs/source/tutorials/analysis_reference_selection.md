# RMSF Analysis: Reference Structure Selection

This guide explains how to choose the appropriate reference structure for RMSF
(Root Mean Square Fluctuation) analysis in PolyzyMD. The choice of reference
affects the interpretation of your results and should match your scientific
question.

## Overview

RMSF measures how much each residue fluctuates around a reference position
during a simulation. Before computing RMSF, the trajectory must be **aligned**
to a reference structure to remove global translation and rotation.

The choice of reference structure determines *what* you're measuring:

| Reference Mode | What RMSF Measures |
|----------------|-------------------|
| `centroid` | Fluctuations around the equilibrium (most visited) conformation |
| `average` | Pure thermal fluctuations around the mathematical mean |
| `frame` | Fluctuations relative to a specific functional state |

## Reference Modes

### Centroid (Default)

```python
calc = RMSFCalculator(
    config,
    reference_mode="centroid",  # Default
    equilibration="100ns",
)
```

The **centroid mode** finds the most populated conformational state using
K-Means clustering on all protein atoms. It returns the actual frame from
the trajectory that best represents where the protein spends most of its time.

**When to use:**
- Standard protein flexibility analysis
- You want to know how much residues deviate from the equilibrium state
- You have a well-behaved single-state protein

**Implementation details:**
- Uses K-Means clustering with k=1 on all protein atoms (including side chains)
- Finds the frame closest to the cluster center
- This is a real sampled conformation, not a synthetic structure

### Average Structure

```python
calc = RMSFCalculator(
    config,
    reference_mode="average",
    equilibration="100ns",
)
```

The **average mode** computes the mathematical mean of all atomic positions
and aligns the trajectory to this synthetic structure.

**When to use:**
- You want a pure measure of thermal fluctuations
- Comparing to literature values that used average-based RMSF
- The protein samples a single conformational basin

**Caveats:**
- The average structure may have unphysical geometry (distorted bond lengths,
  angles) because it's a mathematical construct, not a real conformation
- If the protein transitions between distinct states, the average may not
  represent any physically meaningful conformation

### Specific Frame

```python
calc = RMSFCalculator(
    config,
    reference_mode="frame",
    reference_frame=500,  # 1-indexed frame number
    equilibration="100ns",
)
```

The **frame mode** aligns to a user-specified frame. This is powerful for
analyzing fluctuations relative to a known functional state.

**When to use:**
- Enzyme catalysis studies where you want to measure deviations from a
  "catalytically competent" conformation
- Comparing flexibility before/after a conformational change
- You have identified a specific frame of interest (e.g., substrate bound,
  active site properly configured)

**Example: Catalytic Competence Analysis**

For enzyme studies, you might want to find a frame where:
- The catalytic triad residues are properly hydrogen-bonded
- The substrate is positioned correctly
- The active site is in the "reactive" configuration

```python
# First, identify the catalytically competent frame
# (e.g., by analyzing hydrogen bond distances, active site geometry)
catalytic_frame = 1523  # Frame where catalytic triad is aligned

# Then compute RMSF relative to that state
calc = RMSFCalculator(
    config,
    reference_mode="frame",
    reference_frame=catalytic_frame,
    equilibration="100ns",
)
result = calc.compute(replicate=1)

# High RMSF values now indicate residues that frequently deviate
# from the catalytically competent geometry
```

## Alignment Selection

Independently of the reference mode, you can control which atoms are used
for the alignment superposition:

```python
calc = RMSFCalculator(
    config,
    reference_mode="centroid",
    alignment_selection="protein and name CA",  # Default: backbone alignment
    centroid_selection="protein",  # Default: all protein atoms for clustering
)
```

**Alignment selection** (`alignment_selection`):
- Controls which atoms are used to superimpose frames
- Default is `"protein and name CA"` (alpha carbons only)
- This removes global rotation/translation of the protein backbone

**Centroid selection** (`centroid_selection`):
- Controls which atoms are used for K-Means clustering (centroid mode only)
- Default is `"protein"` (all protein atoms)
- Using all atoms captures side chain conformations in the clustering

## Practical Recommendations

### For General Flexibility Analysis

Use the default centroid mode:

```bash
polyzymd analyze rmsf -c config.yaml -r 1 --eq-time 100ns
```

This gives you flexibility relative to the equilibrium state, which is
usually what you want for comparing flexibility across conditions or mutants.

### For Enzyme Mechanism Studies

Consider using frame mode with a catalytically relevant structure:

```python
# 1. Run initial analysis to find catalytically competent frames
# 2. Analyze hydrogen bonds, active site geometry, etc.
# 3. Select the best frame
# 4. Re-run RMSF with that frame as reference

calc = RMSFCalculator(
    config,
    reference_mode="frame",
    reference_frame=catalytic_frame_id,
    equilibration="100ns",
)
```

RMSF values will then tell you which residues deviate from the active
conformation, potentially disrupting catalysis.

### For Comparing to Published Data

Check what reference the published study used:
- Many older studies use average structures
- Some use the first frame or crystal structure
- Match the methodology for fair comparison

## Technical Notes

### Why Alignment Matters

Without alignment, RMSF measures the **total displacement** of atoms, which
includes:
1. Global translation (protein drifting through the box)
2. Global rotation (protein tumbling)
3. Internal flexibility (what we actually want)

Alignment removes (1) and (2), leaving only the internal flexibility signal.
This is why unaligned RMSF values can be ~50 Å while properly aligned values
are typically 0.5-5 Å for a stable protein.

### Memory Considerations

All alignment modes use in-memory trajectory alignment (`in_memory=True`).
For very long trajectories, this may require significant RAM. Consider:
- Using a larger equilibration time to reduce frames
- Processing replicates sequentially rather than in parallel

### Reproducibility

The result files store all alignment metadata:
- `reference_mode`: Which method was used
- `reference_frame`: The specific frame (1-indexed, or None for average)
- `alignment_selection`: Which atoms were used for superposition

This ensures your analysis is fully reproducible.
