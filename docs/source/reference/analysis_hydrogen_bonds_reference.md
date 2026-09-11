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
| `allow_empty_groups` | `bool` | `false` | When `false`, an empty group selection raises `polyzymd.analyses.exceptions.SelectionError` naming the group and the selection. When `true`, summaries that use the empty group are warned about and skipped |
| `donor_acceptor_elements` | `tuple[str, ...]` | `["N", "O"]` | Elements allowed to act as donors and acceptors. Symbols are capitalized and de-duplicated; `H` is rejected |
| `allow_overlapping_composition` | `bool` | `false` | Permit overlapping composition partitions |
| `composition` | mapping or `null` | `null` | Optional donor/acceptor partition analysis |
| `hydrogens_selection` | `str \| null` | `null` | Advanced explicit-hydrogen selection override for unusual atom names |
| `timestep_ps` | `float \| null` | `null` | Optional uniform frame spacing for time-axis plots |

Each summary defines exactly one of `between: [group_a, group_b]` or
`within: group_name`. Mapping-form summaries use the mapping key as the summary
name.

Hydrogen detection uses MDAnalysis `HydrogenBondAnalysis` and requires explicit
hydrogen atoms in the topology. Both `donors_sel` and `acceptors_sel` are set to
`(<group union>) and element <donor_acceptor_elements>`, which is
`(<group union>) and element N O` with the defaults. PolyzyMD prefers
`(<group union>) and (element H)` for hydrogen selection. For GRO-like topologies
that do not provide MDAnalysis `elements`, PolyzyMD first tries to infer missing
elements safely from atom types or atom names. If element metadata remains
unavailable, `hydrogen_bonds` raises `SelectionError` instead of widening the
donor and acceptor selections to every atom. Set `hydrogens_selection` only for
unusual explicit-hydrogen naming schemes.

## Selection provenance

Replicate and condition artifacts record the selections that were used.

| Location | Key | Content |
|----------|-----|---------|
| `provenance` | `donor_acceptor_selection_policy.donors_selection` | Effective `donors_sel` string |
| `provenance` | `donor_acceptor_selection_policy.acceptors_selection` | Effective `acceptors_sel` string |
| `provenance` | `donor_acceptor_selection_policy.hydrogens_selection` | Effective `hydrogens_sel` string |
| `provenance` | `donor_acceptor_selection_policy.elements` | Elements the donor and acceptor selections were restricted to |
| `provenance` | `hydrogens_selection_policy.source` | `element`, `user`, or `name_fallback` |
| `metadata` | `donors_selection_string`, `acceptors_selection_string`, `hydrogens_selection_string` | The same strings, also written to the NPZ event sidecar |

`donor_acceptor_elements` is part of the settings fingerprint, so changing it
changes the cache identity and results computed with the previous default are
recomputed.

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
