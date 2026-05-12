# Secondary Structure Plugin Reference

For task-oriented guidance, start with {doc}`../how_to/analysis_chooser` and
then enable `secondary_structure` in `comparison.yaml`.

## Settings

Top-level plugin key: `plugins.secondary_structure`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `chain_id` | `str` | `"A"` | Protein chain letter passed to DSSP. Chain A is the PolyzyMD protein convention. |

The plugin computes simplified DSSP classes with MDTraj, aggregates helix,
strand, and coil persistence across replicates, and uses default scalar
comparison on `helix_fraction`.

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
