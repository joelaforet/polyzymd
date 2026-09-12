# Secondary Structure Plugin Reference

For task-oriented guidance, start with {doc}`../how_to/analysis_chooser` and
then enable `secondary_structure` in `comparison.yaml`.

The plugin assigns a simplified DSSP class to every selected protein residue in
every production frame with `mdtraj.compute_dssp(simplified=True)`.

## Settings

Top-level plugin key: `plugins.secondary_structure`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `chain_id` | `str` | `"A"` | Protein chain letter passed to DSSP. Chain A is the PolyzyMD protein convention. |
| `selection` | `str \| null` | `null` | Explicit MDAnalysis selection for the protein residues. When set, this overrides `chain_id`. |

By default, PolyzyMD selects `protein and chainid A`. This matches the PDB and
PolyzyMD chain convention where chain A is the protein. GROMACS `.gro`
topologies may not preserve chain IDs, so `chainid` selections can fail in
MDAnalysis. For `.gro` inputs, set `selection` explicitly, for example:

```yaml
plugins:
  secondary_structure:
    selection: "protein"
```

Use MDAnalysis residue syntax to restrict the protein region:

```yaml
plugins:
  secondary_structure:
    selection: "protein and resid 1:269"
```

or zero-based residue indices:

```yaml
plugins:
  secondary_structure:
    selection: "protein and resindex 0:268"
```

DSSP requires complete, backbone-compatible protein residues. A CA-only
selection such as `protein and name CA` is refused with a `ReplicateError`, and
so is a selection that matches no atoms.

## Observables

| Observable | Kind | Unit | Values |
|---|---|---|---|
| `ss_helix` | `fraction` | fraction | Fraction of selected residues assigned H in each frame |
| `ss_strand` | `fraction` | fraction | Fraction assigned E in each frame |
| `ss_coil` | `fraction` | fraction | Fraction assigned C in each frame |
| `ss_unassigned` | `fraction` | fraction | Fraction mdtraj returns as NA in each frame |
| `helix_occupancy` | `profile` | fraction | Fraction of the window each residue spends in helix, indexed by residue ID |
| `strand_occupancy` | `profile` | fraction | Fraction of the window each residue spends in strand, indexed by residue ID |

The four fractions sum to 1.0 in every frame. Each replicate contributes the
mean over its frames, and the condition-level mean, SEM and 95 percent interval
are computed across replicates. Profiles are averaged element-wise across
replicates and are not tested pairwise.

`ss_unassigned` counts residues mdtraj cannot assign because they have no usable
backbone or carry a residue name mdtraj does not know. MDAnalysis accepts many
more residue names under the `protein` keyword than mdtraj does, so a non-zero
`ss_unassigned` means the selection needs attention rather than that the protein
is disordered. Releases before this one scored those residues as coil.

## Output files

Per-replicate results are written under
`analysis/<condition>/secondary_structure/run_<replicate>/`, with the per-frame
series in `observables.npz` beside `result.json`. Aggregated results are written
under `analysis/<condition>/secondary_structure/aggregated/`, and cross-condition
statistics to `comparison/secondary_structure/result.json`.

## References

- Kabsch, W. and Sander, C. (1983). Dictionary of protein secondary structure:
  pattern recognition of hydrogen-bonded and geometrical features.
  *Biopolymers*, 22(12), 2577-2637. doi:10.1002/bip.360221211
- McGibbon, R. T. et al. (2015). MDTraj: a modern open library for the analysis
  of molecular dynamics trajectories. *Biophysical Journal*, 109(8), 1528-1532.
  doi:10.1016/j.bpj.2015.08.015
