# Hydrogen Bonds Plugin Reference

This page is lookup documentation for the `hydrogen_bonds` analysis plugin:
settings, selection behavior, the observables it reports, and the files it
writes.

For task-oriented setup examples, see {doc}`../how_to/hydrogen_bonds`.

## Plugin key

Top-level comparison YAML key: `plugins.hydrogen_bonds`.

```yaml
plugins:
  hydrogen_bonds:
    groups:
      protein: "chainid A"
      polymer: "chainid C"
    summaries:
      protein_polymer:
        between: [protein, polymer]
```

## Settings

### `plugins.hydrogen_bonds`

| Field | Type | Default | Constraints | Description |
|-------|------|---------|-------------|-------------|
| `groups` | mapping | `{protein: "chainid A", polymer: "chainid C"}` | | Named MDAnalysis selections the summaries read. |
| `summaries` | list or mapping | one `protein_polymer` summary | names unique; every named group defined | Partitions to report. A mapping uses its keys as the summary names. |
| `distance_cutoff` | float | `3.0` | `> 0` | Donor-acceptor distance cutoff in angstrom. |
| `angle_cutoff` | float | `150.0` | `> 0`, `<= 180` | D-H...A angle cutoff in degrees. |
| `donor_acceptor_elements` | list of string | `["N", "O"]` | known element symbols, not `H` | Elements allowed to donate and accept. Symbols are capitalized and de-duplicated, an unknown symbol is rejected when the config is read, and `H` is rejected because hydrogens are selected separately. Add `S` to include thiols and methionine sulfur. |
| `update_selections` | bool | `true` | | Re-evaluate MDAnalysis's donor, hydrogen and acceptor selections on each frame. It does not re-evaluate group membership, which is always fixed at the start of the window. |
| `allow_empty_groups` | bool | `false` | | When `false`, a group matching no atoms raises `SelectionError`. When `true`, the group is warned about and the summaries that read it report zero. |
| `top_n_pairs` | int | `15` | `>= 1` | Residue pairs kept in each occupancy profile. |
| `hydrogens_selection` | string or null | `null` | | Override for the hydrogen selection. `element H` is used when it is omitted. Set it only for topologies with unusual hydrogen naming. |
| `composition` | mapping or null | `null` | | Deprecated since v1.3 and ignored. The event sidecar carries every bond with its donor and acceptor atom indices, so a partition breakdown is a query over that table. Setting it raises a `DeprecationWarning`; it is rejected in v1.4. |
| `allow_overlapping_composition` | bool | `false` | | Deprecated since v1.3 and ignored, with `composition`. |
| `timestep_ps` | float or null | `null` | | Deprecated since v1.3 and ignored. The framework resolves the frame window and its time axis. |

### `summaries` entries

| Field | Type | Default | Constraints | Description |
|-------|------|---------|-------------|-------------|
| `name` | string | required | non-empty; unique | Name used in the observable names. A mapping-form summary takes its key. |
| `between` | pair of group names | `null` | exactly one of `between` and `within` | Bonds whose donor is in one group and whose acceptor is in the other, in either direction. |
| `within` | group name | `null` | exactly one of `between` and `within` | Bonds whose donor and acceptor are both in that group. |

## Selection behavior

MDAnalysis `HydrogenBondAnalysis` runs once over the union of the groups every
summary reads, and each summary is then a filter on the resulting event table.
Detection needs explicit hydrogen atoms and element metadata in the topology.

Both `donors_sel` and `acceptors_sel` are `(<group union>) and element <elements>`,
which is `(<group union>) and element N O` with the defaults. Without that
restriction MDAnalysis would treat any atom of `donors_sel` near a selected
hydrogen as a donor, which admits C-H donors and carbon acceptors that fall
outside the IUPAC definition. The configured symbols are matched against the
spellings the topology carries, so a topology that writes `CL` is still found by
an entry of `Cl`. The trajectory loader infers elements for GRO-like topologies
that carry only atom types or atom names, so those work too; a topology where
elements remain unavailable raises `SelectionError` rather than widening the
selection back to every atom.

Group membership is resolved once, from the first frame of the window, and then
held fixed. A coordinate-dependent group selection such as `around` or `sphzone`
therefore assigns a bond to the summary its atoms belonged to at that frame, not
the one they would belong to in the frame the bond was found in. The plugin logs
a warning naming any group whose selection uses one of those keywords.

Bonds whose donor and acceptor sit in the same residue are dropped from every
summary, because a residue hydrogen bonded to itself says nothing about the
interaction between the groups.

Groups may overlap. A bond in the overlap is counted by every summary whose
filter it passes, so overlapping groups double count by design.

## Observables

Each configured summary reports two observables.

| Name | Kind | Unit | Description |
|------|------|------|-------------|
| `hbonds_<summary>` | `mean_of_timeseries` | `count` | Number of hydrogen bonds matching that summary, one value per frame. |
| `pair_occupancy_<summary>` | `profile` | `fraction` | Occupancy of the `top_n_pairs` most persistent residue pairs, indexed by rank. |

Occupancy is the fraction of the window in which a residue pair held at least
one hydrogen bond, counted undirected, so a bond found in either direction
counts towards the same pair. The profile is indexed by rank rather than by
pair, because the pairs that reach the top differ between replicates: rank 0 is
the most occupied pair of that replicate. The pair itself is named in the
observable's `metadata["pair_labels"]`, in the form `TYR138(A)-SBM152(C)`,
ordered by chain and residue ID, and `index_label` is `occupancy rank`. Ranks
past the observed pairs are zero with an empty label, so every replicate of a
condition shares one index. `metadata["n_pairs_observed"]` records how many
pairs the summary saw in total.

Ranking inside a replicate biases the low ranks upward. Rank 0 is the largest of
many noisy occupancies, so its mean across replicates sits above the mean
occupancy of any one pair, and the bias grows with the number of candidate
pairs. Read the profile as the shape of a replicate's occupancy spectrum, and
use `metadata["pair_labels"]` to check whether the replicates agree on which
pairs are there. A specific pair worth testing belongs in its own summary, where
it becomes a `mean_of_timeseries` the framework compares properly.

## Event sidecar

Every replicate writes the raw MDAnalysis event table to
`sidecars/hydrogen_bond_events.npz` under the replicate directory, as one array
named `hydrogen_bond_events` with six columns: frame index, donor atom index,
hydrogen atom index, acceptor atom index, donor-acceptor distance in angstrom,
and D-H...A angle in degrees. Those names are written beside it as
`sidecars/hydrogen_bond_event_columns.npz`, so the table carries its own schema.
Atom indices are zero-based into the topology.
The table is the full detection output before any summary filter, so a question
the summaries do not answer can be asked of it directly.

## Canonical output paths

| Level | Path | Contents |
|-------|------|----------|
| Per replicate | `analysis/<condition>/hydrogen_bonds/run_<replicate>/result.json` | `ReplicateArtifact` whose `payload.observables` holds one reduced estimate per observable. |
| Per replicate sidecar | `analysis/<condition>/hydrogen_bonds/run_<replicate>/observables.npz` | Per-frame counts and per-rank occupancies, one array per observable name. |
| Per replicate sidecar | `analysis/<condition>/hydrogen_bonds/run_<replicate>/sidecars/hydrogen_bond_events.npz` | The raw event table. |
| Per condition | `analysis/<condition>/hydrogen_bonds/aggregated/result.json` | `ConditionArtifact` whose `payload.observables` holds one aggregate per observable. |
| Cross condition | `comparison/hydrogen_bonds/result.json` | `ComparisonArtifact` with the per-condition aggregates and the pairwise tests. |

## Artifact fields

A replicate estimate carries `name`, `kind`, `unit`, `n_frames`, the reduced
`value` for a time-series kind or the `profile`, `index` and `index_label` for a
profile, the plugin's `metadata`, and the correlation diagnostics `statistical_inefficiency` and `n_eff`.
The diagnostics are reported, never used to shrink an error bar.

A condition aggregate carries `replicate_values`, `mean`, `sem`, `ci95_low`,
`ci95_high`, `ci_method`, `coverage` and `n_replicates` for a time-series kind,
and `profile_mean`, `profile_sem`, `index` and `index_label` for a profile.
Every statistic is computed across replicates, never across frames. The
`metadata` on an aggregate is the first replicate's, so the pair labels there
name that replicate's ranks and not the condition's.

A comparison entry carries `control`, `condition`, `delta`, `percent_change`,
`test`, `p_value`, `p_adjusted`, `correction`, `cohens_d`, `significant`,
`testable` and `note`. Profiles are aggregated but not tested pairwise.

## Interpretation

A higher `hbonds_protein_polymer` in one condition than in another means the
polymer holds more hydrogen bonds to the protein per frame. The `delta` in a
comparison entry is in bonds per frame and `percent_change` is relative to the
control mean. `higher_is_better` is unset, because more or fewer hydrogen bonds
can be desirable depending on what the formulation is meant to do.

Counts are sensitive to the cutoffs. The defaults of 3.0 angstrom and 150
degrees follow Smith et al. 2019; a different pair of cutoffs gives a different
absolute count, so only compare conditions analyzed with the same settings.

## Figures

`compare run --plot` draws figures from the observable kind in the contract
runner. Until that work lands, a `hydrogen_bonds` run writes artifacts and the
text report but no figures.

## References

- Arunan, E. et al. (2011). Definition of the hydrogen bond (IUPAC
  Recommendations 2011). *Pure and Applied Chemistry*, 83(8), 1637-1641.
  doi:10.1351/PAC-REC-10-01-02
- Smith, P. et al. (2019). On the interaction of hyaluronic acid with synovial
  fluid lipid membranes. *Physical Chemistry Chemical Physics*, 21(19),
  9845-9857. doi:10.1039/C9CP01532A
- Michaud-Agrawal, N. et al. (2011). MDAnalysis: a toolkit for the analysis of
  molecular dynamics simulations. *Journal of Computational Chemistry*, 32(10),
  2319-2327. doi:10.1002/jcc.21787
- Gowers, R. J. et al. (2016). MDAnalysis: a Python package for the rapid
  analysis of molecular dynamics simulations. *Proceedings of the 15th Python in
  Science Conference*, 98-105. doi:10.25080/Majora-629e541a-00e

## See also

- {doc}`../how_to/hydrogen_bonds` (task recipes and commands)
- {doc}`comparison_yaml` (comparison file schema)
- {doc}`analysis_comparison_reference` (shared comparison behavior)
