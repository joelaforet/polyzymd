# Contacts Plugin Reference

For a task-oriented setup and run workflow, see
{doc}`../how_to/analysis_contacts_quickstart`.

## Configuration Reference

Contacts plugin settings live under `plugins.contacts` in `comparison.yaml`.

### Core analysis fields (`ContactsSettings`)

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `polymer_selection` | `str` | `"chainid C"` | MDAnalysis selection for polymer atoms |
| `protein_selection` | `str` | `"chainid A"` | MDAnalysis selection for protein atoms |
| `cutoff` | `float` | `4.5` | Contact distance cutoff in Angstroms |
| `polymer_types` | `list[str] \| None` | `null` | Optional polymer residue-name filter |
| `grouping` | `str` | `"aa_class"` | Protein grouping mode: `aa_class`, `secondary_structure`, or `none` |
| `compute_residence_times` | `bool` | `true` | Compute aggregate residence-time summaries and plots |

Set `compute_residence_times: false` to skip aggregate residence-time summaries
and residence-time plotters. Per-replicate contact events remain stored because
they are the compressed representation used for contact fractions and
contacts-derived analyses. The setting is validated through the canonical
contacts detection fingerprint recorded in replicate and condition artifacts.

### Partition fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `protein_groups` | `dict[str, list[int]] \| None` | `null` | Custom residue groups, e.g. `{active_site: [77, 133]}` |
| `protein_partitions` | `dict[str, list[str]] \| None` | `null` | Named partitions of `protein_groups` for contact-fraction and residence-time plots |

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
            │   └── sidecars/
            │       └── contact_events.npz
            ├── run_2/
            │   └── ...
            ├── run_3/
            │   └── ...
            ├── aggregated/
            │   ├── result.json
            │   └── sidecars/
            │       └── contact_profiles.npz
```

Legacy standalone JSON filenames from pre-artifact contacts runs are no longer
loaded by the v1.3 contacts workflow. Recompute contacts to produce canonical
artifact-store outputs.

### Per-replicate JSON structure (`ReplicateArtifact`)

Representative structure:

```python
{
    "analysis_name": "contacts",
    "replicate": 1,
    "payload": {
        "metrics": {"coverage": 0.74, "mean_contact_fraction": 0.18},
        "event_sidecar": "sidecars/contact_events.npz",
        "n_contact_events": 1240,
        "n_frames_used": 9000
    },
    "sidecars": [{"path": "sidecars/contact_events.npz", "metadata": {"kind": "contact_events"}}],
    "metadata": {
        "contacts_detection_fingerprint": "...",
        "equilibration": "10ns"
    }
}
```

### Aggregated JSON structure (`ConditionArtifact`)

Representative structure:

```python
{
    "analysis_name": "contacts",
    "condition_label": "PEGylated",
    "replicates": [1, 2, 3],
    "payload": {
        "metrics": {
            "coverage": {"values": [0.73, 0.75, 0.74], "mean": 0.74, "sem": 0.01},
            "mean_contact_fraction": {"values": [0.17, 0.19, 0.18], "mean": 0.18, "sem": 0.01}
        },
        "residue_stats": [
            {
                "protein_resid": 77,
                "protein_group": "polar",
                "contact_fraction_mean": 0.211,
                "contact_fraction_per_replicate": [0.201, 0.232, 0.200]
            }
        ],
        "profile_sidecar": "sidecars/contact_profiles.npz",
        "residence_time_by_polymer_type": {
            "SBM": {"mean_ns": 9.60, "sem_ns": 0.53}
        }
    },
    "metadata": {
        "contacts_detection_fingerprint": "...",
        "compute_residence_times": true,
        "equilibration": "10ns"
    }
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

### Contacts plot settings

| Field | Default | Description |
|-------|---------|-------------|
| `generate_contact_fraction_profile` | `true` | Enable per-residue contact-fraction profiles |
| `generate_residence_time_profile` | `true` | Enable per-residue residence-time profiles |
| `generate_cf_by_aa_class_bars` | `true` | Enable contact-fraction AA-class bars |
| `generate_cf_by_partition_bars` | `true` | Enable contact-fraction partition bars |
| `generate_rt_by_aa_class_bars` | `true` | Enable residence-time AA-class bars |
| `generate_rt_by_partition_bars` | `true` | Enable residence-time partition bars |
| `highlight_residues` | `[]` | Residues marked with vertical lines on profile plots |
| `contact_fraction_profile_threshold` | `null` | Optional threshold line on contact-fraction profile |

Figure-size and error-display fields are also available per plot type (for
example `figsize_contact_fraction_profile` and
`show_contact_fraction_profile_error`).

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
- Start with `polymer_selection: "chainid C"` and narrow incrementally
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
