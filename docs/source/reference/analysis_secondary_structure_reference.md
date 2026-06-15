# Secondary Structure Plugin Reference

For task-oriented guidance, start with {doc}`../how_to/analysis_chooser` and
then enable `secondary_structure` in `comparison.yaml`.

## Settings

Top-level plugin key: `plugins.secondary_structure`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `chain_id` | `str` | `"A"` | Protein chain letter passed to DSSP. Chain A is the PolyzyMD protein convention. |
| `selection` | `str \| null` | `null` | Explicit MDAnalysis selection for the protein residues. When set, this overrides `chain_id`. |

The plugin computes simplified DSSP classes with MDTraj, aggregates helix,
strand, and coil persistence across replicates, and uses default scalar
comparison on `helix_fraction`.

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

DSSP requires complete, backbone-compatible protein residues. Do not use
CA-only selections such as `protein and name CA` for this plugin.

## Output files

Per-replicate results are written under
`analysis/<condition>/secondary_structure/run_<replicate>/`. Aggregated results
are written under `analysis/<condition>/secondary_structure/aggregated/`, and
cross-condition statistics are written to
`comparison/secondary_structure/result.json`.

## Plot outputs

| Plot output | Description |
|-------------|-------------|
| `ss_timeline_<condition>.png` | Residue-by-time secondary-structure heatmap for one condition |
| `ss_content_bars.png` | Condition-level helix, strand, and coil fractions |
| `ss_helix_bars.png`, `ss_strand_bars.png`, `ss_coil_bars.png` | One bar chart per secondary-structure class |
| `ss_persistence_diff_heatmap.png` | Difference in helix persistence relative to the control condition |

Time-axis plots assume uniformly saved frames. PolyzyMD maps frame index to time
as `frame_index * dt`; variable-timestep concatenated trajectories are not
supported.
