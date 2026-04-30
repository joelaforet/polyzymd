# Contacts Plugin Reference

For a task-oriented setup and run workflow, see
{doc}`../how_to/analysis_contacts_quickstart`.

## Configuration Reference

Contacts plugin settings live under `plugins.contacts` in `comparison.yaml`.

### Core analysis fields (`ContactsSettings`)

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `polymer_selection` | `str` | `"chainID C"` | MDAnalysis selection for polymer atoms |
| `protein_selection` | `str` | `"protein"` | MDAnalysis selection for protein atoms |
| `cutoff` | `float` | `4.5` | Contact distance cutoff in Angstroms |
| `polymer_types` | `list[str] \| None` | `null` | Optional polymer residue-name filter |
| `grouping` | `str` | `"aa_class"` | Protein grouping mode: `aa_class`, `secondary_structure`, or `none` |
| `compute_residence_times` | `bool` | `true` | Compute aggregate residence-time summaries and plots |

Set `compute_residence_times: false` to skip aggregate residence-time summaries
and residence-time plotters. Per-replicate contact events remain stored because
they are the compressed representation used for contact fractions and
contacts-derived analyses. The setting participates in contacts cache identity;
older caches without this identity are accepted only for the default `true`
setting when their other metadata is valid.

### Binding preference and partition fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `compute_binding_preference` | `bool` | `false` | Enable binding-preference enrichment pipeline |
| `surface_exposure_threshold` | `float` | `0.2` | Relative SASA threshold for surface exposure |
| `enzyme_pdb_for_sasa` | `str \| None` | `null` | Optional enzyme PDB path for SASA computation |
| `include_default_aa_groups` | `bool` | `true` | Include default AA-class groups |
| `protein_groups` | `dict[str, list[int]] \| None` | `null` | Custom residue groups, e.g. `{active_site: [77, 133]}` |
| `protein_partitions` | `dict[str, list[str]] \| None` | `null` | Named partitions of `protein_groups` |
| `polymer_type_selections` | `dict[str, str] \| None` | `null` | Custom polymer type mappings by selection |
| `polymer_chain` | `str` | `"C"` | Polymer chain ID used for auto-detection |
| `enrichment_normalization` | `str` | `"residue"` | Deprecated backward-compatibility field (ignored) |

### Comparison output fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `fdr_alpha` | `float` | `0.05` | FDR alpha for Benjamini-Hochberg correction |
| `min_effect_size` | `float` | `0.5` | Minimum Cohen's d to flag/highlight |
| `top_residues` | `int` | `10` | Number of top residues shown in console output |

### Validation notes

- `grouping` must be one of `aa_class`, `secondary_structure`, or `none`
- `fdr_alpha` must be between 0 and 1
- If `protein_partitions` is provided, `protein_groups` must also be provided
- Partition group names must exist in `protein_groups`
- Residues cannot overlap across groups within the same partition

## Output Files

Contacts results are written under each condition's analysis directory.

```text
<projects_directory>/
└── analysis/
    └── <condition>/
        └── contacts/
            ├── run_1/
            │   ├── result.json
            │   └── contacts_eq10ns_cut4.5_s<fingerprint>_rep1.json
            ├── run_2/
            │   └── ...
            ├── run_3/
            │   └── ...
            ├── aggregated/
            │   └── result.json
            └── contacts_aggregated_eq10ns_cut4.5_s<fingerprint>_reps1-3.json
```

Depending on cache history and plugin version, legacy file names may also
exist (for example `contacts_rep1.json`).

### Per-replicate JSON structure (`ContactResult`)

Representative structure:

```python
{
    "analysis_type": "contacts",
    "residue_contacts": [
        {
            "protein_resid": 77,
            "protein_resname": "SER",
            "protein_group": "polar",
            "segment_contacts": [
                {
                    "polymer_resname": "SBM",
                    "polymer_resid": 403,
                    "polymer_chain_idx": 0,
                    "events": [{"start_frame": 120, "duration": 9}]
                }
            ],
            "statistical_inefficiency": 2.41,
            "n_effective": 3733.6
        }
    ],
    "n_frames": 9000,
    "timestep_ps": 10.0,
    "criteria_label": "any_atom_4.5A",
    "criteria_cutoff": 4.5,
    "start_frame": 1000,
    "schema_version": 2,
    "metadata": {
        "target_selector": "protein",
        "query_selector": "chainID C",
        "algorithm": "capped_distance"
    }
}
```

### Aggregated JSON structure (`AggregatedContactResult`)

Representative structure:

```python
{
    "analysis_type": "contacts_aggregated",
    "n_replicates": 3,
    "total_frames_per_replicate": [9000, 9000, 9000],
    "timestep_ps": 10.0,
    "criteria_label": "any_atom_4.5A",
    "criteria_cutoff": 4.5,
    "coverage_mean": 0.740,
    "coverage_sem": 0.011,
    "mean_contact_fraction": 0.180,
    "mean_contact_fraction_sem": 0.004,
    "group_stats": {
        "aromatic": [0.242, 0.013],
        "polar": [0.168, 0.009]
    },
    "residence_time_by_polymer_type": {
        "SBM": [9.60, 0.53],
        "EGM": [8.14, 0.56]
    },
    "residue_stats": [
        {
            "protein_resid": 77,
            "protein_group": "polar",
            "contact_fraction_mean": 0.211,
            "contact_fraction_sem": 0.016,
            "contact_fraction_per_replicate": [0.201, 0.232, 0.200],
            "by_polymer_type": {"SBM": [0.173, 0.012]},
            "residence_time_by_polymer_type": {"SBM": [7.2, 0.8]}
        }
    ],
    "metadata": {"aggregation_method": "mean_sem"}
}
```

## Plot Types

Contacts plots are generated through the comparison plotting workflow
(`polyzymd compare plot-all ...`) and controlled by `plot_settings.contacts`.

### Plot outputs

| Output stem | Description | Gate setting |
|-------------|-------------|--------------|
| `contact_fraction_profile` | Per-residue contact-fraction profile across conditions | `generate_contact_fraction_profile` |
| `contact_fraction_profile_<polymer_type>` | Per-residue profile split by polymer type (when multiple polymer types exist) | `generate_contact_fraction_profile` |
| `residence_time_profile` | Per-residue mean residence-time profile (ns) | `generate_residence_time_profile` |
| `residence_time_profile_<polymer_type>` | Per-residue residence-time profile by polymer type | `generate_residence_time_profile` |
| `cf_by_aa_class_bars` | Contact-fraction grouped bars by amino-acid class | `generate_cf_by_aa_class_bars` |
| `cf_by_partition_<partition>_bars` | Contact-fraction grouped bars by user-defined partition | `generate_cf_by_partition_bars` |
| `rt_by_aa_class_bars` | Residence-time grouped bars by amino-acid class | `generate_rt_by_aa_class_bars` |
| `rt_by_partition_<partition>_bars` | Residence-time grouped bars by user-defined partition | `generate_rt_by_partition_bars` |
| `system_coverage_bars` | Coverage-enrichment bars by AA class | `generate_system_coverage_bars` |
| `system_coverage_heatmap` | Coverage-enrichment heatmap | `generate_system_coverage_heatmap` |
| `user_partition_<partition>_bars` | Coverage-enrichment bars for user partition elements | `generate_user_partition_bars` |
| `binding_preference_bars` | Binding-preference enrichment bars | `generate_enrichment_bars` |
| `binding_preference_heatmap` | Binding-preference enrichment heatmap | `generate_enrichment_heatmap` |

### Contacts plot settings

| Field | Default | Description |
|-------|---------|-------------|
| `generate_enrichment_heatmap` | `true` | Enable binding-preference heatmap |
| `generate_enrichment_bars` | `true` | Enable binding-preference bars |
| `generate_system_coverage_heatmap` | `true` | Enable system-coverage heatmap |
| `generate_system_coverage_bars` | `true` | Enable system-coverage bars |
| `generate_user_partition_bars` | `true` | Enable user partition bar plots |
| `generate_contact_fraction_profile` | `true` | Enable per-residue contact-fraction profiles |
| `generate_residence_time_profile` | `true` | Enable per-residue residence-time profiles |
| `generate_cf_by_aa_class_bars` | `true` | Enable contact-fraction AA-class bars |
| `generate_cf_by_partition_bars` | `true` | Enable contact-fraction partition bars |
| `generate_rt_by_aa_class_bars` | `true` | Enable residence-time AA-class bars |
| `generate_rt_by_partition_bars` | `true` | Enable residence-time partition bars |
| `highlight_residues` | `[]` | Residues marked with vertical lines on profile plots |
| `contact_fraction_profile_threshold` | `null` | Optional threshold line on contact-fraction profile |

Figure-size and error-display fields are also available per plot type (for
example `figsize_contact_fraction_profile`,
`show_contact_fraction_profile_error`, `figsize_enrichment_bars`).

For global plotting keys (`style`, `dpi`, output format), see
{doc}`analysis_comparison_reference` and {doc}`comparison_yaml`.

## Common CLI Options

| Option | Default | Description |
|--------|---------|-------------|
| `-f, --file` | `comparison.yaml` | Path to comparison config |
| `--eq-time` | `0ns` | Equilibration time to skip |
| `--recompute` | off | Ignore cache and recompute |
| `--format` | `table` | Output format (`table` or `json`) |
| `-o, --output` | (none) | Write formatted output to file |
| `-q, --quiet` | off | Suppress INFO logs |
| `--debug` | off | Enable DEBUG logging |

Typical run command:

```bash
polyzymd compare run contacts -f comparison.yaml --eq-time 10ns
```

## Troubleshooting

### "No polymer atoms selected"

**Cause:** `polymer_selection` does not match any atoms.

**Fix:**

- Verify chain and residue naming in your topology
- Start with `polymer_selection: "chainID C"` and narrow incrementally
- Run with `--debug` to inspect selection behavior

### "Selection matched no atoms" (protein or polymer)

**Cause:** Selection syntax is valid but does not match this topology.

**Fix:**

- Check residue numbering and atom/residue naming
- Validate that your topology and trajectory belong together

### Missing replicate data / replicate skipped

**Message:** `Skipping replicate N: trajectory data not found`.

**Cause:** Missing files or incomplete simulation output for that replicate.

**Fix:**

- Confirm replicate output paths in the condition config
- Re-run after simulation completion
- Analysis continues with available replicates

### "protein_partitions requires protein_groups to be defined"

**Cause:** Partition references were configured without group definitions.

**Fix:** Add `protein_groups` and reference those names in
`protein_partitions`.

### Unexpected cache reuse after changing settings

**Cause:** Cached files from prior runs are still present.

**Fix:**

- Re-run with `--recompute`
- Or clear the relevant `analysis/<condition>/contacts/` directory

### Slow runtime

**Cause:** Large trajectories and large selections.

**Fix:**

- Increase `--eq-time` to skip equilibration frames
- Restrict `polymer_selection` and/or `protein_selection`
- Use cached results for repeated report generation
