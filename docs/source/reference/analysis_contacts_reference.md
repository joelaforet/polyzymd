# Contacts Plugin Reference

This page is lookup documentation for the `contacts` analysis plugin: settings,
the observables it reports, and the files it writes.

For a task-oriented setup and run workflow, see
{doc}`../how_to/analysis_contacts_quickstart`.

## Plugin key

Top-level comparison YAML key: `plugins.contacts`.

```yaml
plugins:
  contacts:
    protein_selection: "protein"
    polymer_selection: "resname SBM EGM"
    cutoff: 4.5
```

## Settings

### `plugins.contacts`

| Field | Type | Default | Constraints | Description |
|-------|------|---------|-------------|-------------|
| `protein_selection` | string | `"chainid A"` | must match at least one atom | MDAnalysis selection for protein atoms. |
| `polymer_selection` | string | `"chainid C"` | must match at least one atom | MDAnalysis selection for polymer atoms. |
| `cutoff` | float | `4.5` | `> 0` | Contact distance cutoff in angstrom. |
| `polymer_types` | list of string or null | `null` | | Restrict the polymer selection to these residue names. |
| `heavy_atoms_only` | bool | `false` | | Exclude hydrogens from both selections before the cutoff is applied. |
| `allow_single_fragment_fallback` | bool | `false` | | Put every polymer residue in chain 0 when the topology has no bonds, instead of raising `TopologyBondsMissingError`. |
| `residence_time_edges_ns` | list of float | `[0.0, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56, 5.12, 10.24, 20.48]` | at least two edges, increasing | Bin edges of the residence-time distribution in ns. |

### Retired settings

These were accepted before v1.3. They are ignored with a `DeprecationWarning`
so an existing comparison file still loads, and they are rejected in v1.4:
`grouping`, `compute_residence_times`, `protein_groups`, `protein_partitions`,
`fdr_alpha`, `min_effect_size`, `top_residues`, `compute_binding_preference`,
`surface_exposure_threshold`, `enzyme_pdb_for_sasa`, `include_default_aa_groups`,
`polymer_type_selections`, `polymer_chain`, `enrichment_normalization`.

The correction family and its alpha now come from the comparison file's own
`fdr_alpha`, which applies to every plugin in the run. Per-group and
per-partition summaries were a plotting concern; the per-residue profile carries
the same information at full resolution.

## What counts as a contact

A protein residue and a polymer residue are in contact in one frame when at
least one atom of the first is within `cutoff` of at least one atom of the
second. Distances use the minimum image convention with the box of that frame,
through the MDAnalysis `capped_distance` grid search.

Hydrogens count toward the cutoff by default, which is what this plugin has
always done and what the frozen parity reference reproduces. The literature
convention for a 4.5 A criterion is heavy atoms only; `heavy_atoms_only: true`
selects it and lowers every contact count. Hydrogens are excluded by element
where the topology has elements, and by a name test on a topology that does
not, which the run warns about.

A selection matching no atoms raises `SelectionError` rather than reporting
zero contacts. A window holding no frames raises `ReplicateError`.

## Polymer chain identity requires topology bonds

Chain identity comes from bonded fragments. A topology whose polymer atoms
carry no bonds raises `TopologyBondsMissingError`, because assigning every
polymer residue to chain 0 would distort the per-chain event stream without
failing. Set `allow_single_fragment_fallback: true` to opt in to the pre-1.3
behaviour; the run then emits a warning naming the topology.

## Observables

| Name | Kind | Unit | Description |
|------|------|------|-------------|
| `contact_count` | `mean_of_timeseries` | `count` | Number of distinct protein-polymer residue pairs in contact, one value per frame. |
| `coverage_per_frame` | `fraction` | `fraction` | Share of protein residues in contact with any polymer residue, one value per frame. |
| `coverage_any_frame` | `fraction` | `fraction` | Share of protein residues in contact at any point in the window, one value per replicate. It is a function of `contact_fraction`, so it is reported with its interval but declared `tested=False` and kept out of the pairwise tests and the correction family. |
| `contact_fraction` | `profile` | `fraction` | Share of frames each protein residue is in contact, indexed by residue ID. |
| `residence_time_distribution` | `profile` | `fraction` | Share of contact events whose duration falls in each bin, indexed by the lower bin edge in ns. |
| `mean_residence_time` | `profile` | `ns` | Mean duration of the contact events of each protein residue, indexed by residue ID. Zero for a residue with no events. |

Before v1.3 the plugin reported two replicate metrics, and neither name means
what the port reports under it. `mean_contact_fraction` is the mean of
`coverage_per_frame` over frames, the same number to floating-point noise. The
old `coverage` is `coverage_any_frame`.

An event runs from the first frame a residue pair is within the cutoff to the
last consecutive frame it stays within it. Durations are whole numbers of
frames, so the time axis of the window has to be evenly spaced; an uneven one
raises `ReplicateError`. An event already running when the window opens is
measured from the window's first frame and one still running when it closes is
measured to its last, so both are shortened, and neither is marked. An event
longer than the top bin edge is counted in the top bin; how many were is in
`metadata["residence_time_overflow_events"]`, so widen
`residence_time_edges_ns` when that count is not small.

Every observable carries the cutoff, the PBC policy, `heavy_atoms_only`, the
bond source of the chain identity, the number of polymer chains, the
residence-time overflow count and any retired settings the run ignored in
`metadata`, at replicate level and at condition level.

## Canonical output paths

| Level | Path | Contents |
|-------|------|----------|
| Per replicate | `analysis/<condition>/contacts/run_<replicate>/result.json` | `ReplicateArtifact` whose `payload.observables` holds one reduced estimate per observable. |
| Per replicate sidecar | `analysis/<condition>/contacts/run_<replicate>/observables.npz` | Full per-frame series and per-index vectors, one array per observable name. |
| Per replicate sidecar | `analysis/<condition>/contacts/run_<replicate>/sidecars/contact_events.npz` | The event table under the key `contact_events`. |
| Per condition | `analysis/<condition>/contacts/aggregated/result.json` | `ConditionArtifact` whose `payload.observables` holds one aggregate per observable. |
| Cross condition | `comparison/contacts/result.json` | `ComparisonArtifact` with the per-condition aggregates and the pairwise tests. |

### The event table

`contact_events.npz` holds one integer array of shape `(n_events, 4)`. The
columns are the protein residue ID, the polymer chain index, the first frame of
the event and its last frame, both as trajectory frame numbers.

## Artifact fields

A replicate estimate carries `name`, `kind`, `unit`, `n_frames`, the reduced
`value` for a time-series kind or the `profile` and `index` for a profile, and
the correlation diagnostics `statistical_inefficiency` and `n_eff`. The
diagnostics are reported, never used to shrink an error bar.

A condition aggregate carries `replicate_values`, `mean`, `sem`, `ci95_low`,
`ci95_high`, `ci_method`, the interval's own `coverage` field and `n_replicates`
for a time-series kind,
and `profile_mean`, `profile_sem` and `index` for a profile. Every statistic is
computed across replicates, never across frames.

A comparison entry carries `control`, `condition`, `delta`, `percent_change`,
`test`, `p_value`, `p_adjusted`, `correction`, `cohens_d`, `significant`,
`testable` and `note`. Profiles are aggregated but not tested pairwise.

## Units

| Quantity | Unit |
|----------|------|
| `cutoff` | angstrom |
| `contact_count` | count of residue pairs |
| `coverage_per_frame`, `coverage_any_frame`, `contact_fraction` | fraction |
| `residence_time_edges_ns`, `mean_residence_time` | ns |

## References

- Michaud-Agrawal, N., Denning, E. J., Woolf, T. B. & Beckstein, O. (2011).
  MDAnalysis: a toolkit for the analysis of molecular dynamics simulations.
  *Journal of Computational Chemistry*, 32(10), 2319-2327. doi:10.1002/jcc.21787
- Grossfield, A. et al. (2018). Best practices for quantifying the uncertainty
  in molecular simulations. *Living Journal of Computational Molecular Science*,
  1(1), 5067. doi:10.33011/livecoms.1.1.5067

## See also

- {doc}`../how_to/analysis_contacts_quickstart` (task recipes and commands)
- {doc}`comparison_yaml` (comparison file schema)
- {doc}`analysis_comparison_reference` (shared comparison behavior)
