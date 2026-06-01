# Hydrogen Bonds Plugin Reference

For task-oriented setup examples, see {doc}`../how_to/hydrogen_bonds`.

## Settings

Top-level plugin key and CLI name: `hydrogen_bonds`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `groups` | `dict[str, str]` | `{protein: "chainid A", polymer: "chainid C"}` | Named MDAnalysis selections used by summaries |
| `summaries` | list or mapping | `protein_polymer` between `protein` and `polymer` | Named H-bond summaries to compute |
| `distance_cutoff` | `float` | `3.0` | Donor-acceptor cutoff in Å |
| `angle_cutoff` | `float` | `150.0` | D-H...A angle cutoff in degrees |
| `update_selections` | `bool` | `true` | Re-evaluate selections on each frame |
| `top_n_pairs` | `int` | `15` | Top residue pairs shown in output and plots |
| `allow_empty_groups` | `bool` | `true` | Warn and skip summaries with empty groups instead of raising |
| `allow_overlapping_composition` | `bool` | `false` | Permit overlapping composition partitions |
| `composition` | mapping or `null` | `null` | Optional donor/acceptor partition analysis |
| `timestep_ps` | `float \| null` | `null` | Optional uniform frame spacing for time-axis plots |

Each summary defines exactly one of `between: [group_a, group_b]` or
`within: group_name`. Mapping-form summaries use the mapping key as the summary
name.

Hydrogen detection uses MDAnalysis `HydrogenBondAnalysis` with hydrogens
selected as `(<group union>) and element H`; topologies need explicit hydrogens
and reliable element metadata.

## Comparison mode

`hydrogen_bonds` overrides result loading to validate settings-sensitive cache
files, then uses the framework's default-style scalar comparison. It extracts
one `MetricValue` per configured summary named `mean_hbonds_<summary>` and runs
FDR-corrected pairwise tests and ANOVA per summary. `higher_is_better` is unset
because more or fewer H-bonds can be desirable depending on the system.

## Output files

Per-replicate cache files are named `hbonds_eq*.json` under
`analysis/<condition>/hydrogen_bonds/run_<replicate>/`. Aggregated results are
written under `analysis/<condition>/hydrogen_bonds/aggregated/`, and
cross-condition statistics are written to `comparison/hydrogen_bonds/result.json`.

## Plot outputs

| Plot output | Description |
|-------------|-------------|
| `hbond_summary_comparison.png` | Faceted mean H-bonds/frame bars for each summary |
| `hbond_timeseries_<summary>.png` | Mean H-bonds/frame over time for one summary |
| `hbond_top_pairs_<summary>.png` | Highest-occupancy residue-pair bars for one summary |
| `hbond_composition_absolute.png` | Stacked composition by donor/acceptor partition, when composition is enabled |
| `hbond_composition_fraction.png` | Fractional composition by donor/acceptor partition, when composition is enabled |

Time-axis plots assume uniformly saved frames. PolyzyMD maps frame index to time
as `frame_index * timestep_ps`; variable-timestep concatenated trajectories are
not supported.
