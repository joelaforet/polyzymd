# SASA Plugin Reference

This page is lookup documentation for the `sasa` analysis plugin: settings,
selection behavior, output paths, artifact fields, comparison outputs, and plot
files.

For a guided workflow, see {doc}`../tutorials/sasa_analysis`. For practical
recipes and commands, see {doc}`../how_to/analysis_sasa_quickstart`.

## Plugin key

Top-level comparison YAML key: `plugins.sasa`.

```yaml
plugins:
  sasa:
    runs:
      - label: "protein_total"
        target_selection: "protein"
```

## Settings

### `plugins.sasa`

| Field | Type | Default | Constraints | Description |
|-------|------|---------|-------------|-------------|
| `runs` | list | required | at least one entry; labels must be unique | Named SASA computations to run. |
| `probe_radius_nm` | float | `0.14` | `> 0` | Shrake-Rupley probe radius in nanometers. |
| `n_sphere_points` | int | `960` | `>= 100` | Number of test points on each atom sphere. Higher is more accurate and slower. |
| `chunk_size` | int | `100` | `>= 1` | Frames processed per chunk for memory-managed computation. |

### `runs` entries

| Field | Type | Default | Constraints | Description |
|-------|------|---------|-------------|-------------|
| `label` | string | required | non-empty; must not contain `/` or `\` | Human-readable run label used in summaries and plot filenames. |
| `target_selection` | string | required | non-empty | MDAnalysis selection for atoms whose SASA is reported. |
| `context_selection` | string or null | `target_selection` | blank values become omitted | MDAnalysis selection for atoms included as surface blockers during SASA computation. |
| `stride` | int | `1` | `>= 1` | Analyze every Nth selected frame. |

## Target and context behavior

SASA runs separate what is reported from what can block the surface:

| Selection | Behavior |
|-----------|----------|
| `target_selection` | Defines the atoms/residues whose SASA values are summarized. |
| `context_selection` | Defines the atoms present in the Shrake-Rupley surface calculation. |

If `context_selection` is omitted, PolyzyMD sets it equal to
`target_selection`. This is useful for self-SASA measurements such as whole
protein SASA.

Examples:

| Goal | `target_selection` | `context_selection` |
|------|--------------------|---------------------|
| Whole-protein self-SASA | `protein` | `protein` or omitted |
| Protein SASA with polymer shielding | `protein` | `protein or chainid C` |
| Active-site SASA | `protein and (resid 77 or resid 156 or resid 262)` | `protein` |
| Active-site SASA with polymer shielding | `protein and (resid 77 or resid 156 or resid 262)` | `protein or chainid C` |
| Monomer-specific shielding | `protein` | `protein or resname SBMA` |

The project chain convention is A = protein, B = substrate, C = polymer, and
D+ = solvent/ions/other. Use lowercase `chainid` in MDAnalysis selections.

## Canonical output paths

SASA writes canonical artifact outputs for compute and aggregate stages, plus a
comparison output and plots.

| Level | Path | Contents |
|-------|------|----------|
| Per replicate | `analysis/<condition>/sasa/run_<replicate>/result.json` | `ReplicateArtifact` envelope with per-run payload summaries and sidecar references. |
| Per replicate sidecars | `analysis/<condition>/sasa/run_<replicate>/sidecars/*.npz` | Large arrays such as per-frame total SASA and per-residue SASA. |
| Per condition | `analysis/<condition>/sasa/aggregated/result.json` | `ConditionArtifact` envelope with aggregate per-run summaries across replicates. |
| Per condition sidecars | `analysis/<condition>/sasa/aggregated/sidecars/*.npz` | Aggregate arrays and supporting data, when written. |
| Cross condition | `comparison/sasa/result.json` | Comparison output with condition summaries, pairwise tests, ANOVA-by-run, rankings, and metadata. |
| Plots | `comparison/sasa/plots/` | SASA comparison, normalized-control, time-series, and profile plots. |

## Artifact envelope fields

Replicate and condition JSON files are artifact envelopes. The stable public
concepts are the artifact envelope and the canonical paths, not private helper
classes.

| Field | Meaning |
|-------|---------|
| `analysis_name` | Analysis plugin name, usually `sasa`. |
| `condition_label` | Comparison condition label. |
| `replicate` | Replicate number for replicate artifacts; absent or not meaningful for condition artifacts. |
| `payload` | JSON-compatible SASA summaries, metrics, run labels, and relative sidecar paths. |
| `metadata` | Settings fingerprints, software versions, equilibration labels, units, and related run metadata. |
| `provenance` | Input trajectory/topology identity and workflow details. |
| `sidecars` | Validated references to large sidecar files, including relative paths and integrity metadata. |

Common replicate `payload` keys include:

| Key | Meaning |
|-----|---------|
| `run_results` | List of per-run summaries for the replicate. |
| `n_runs` | Number of configured SASA runs. |
| `n_frames_total` | Total frames available after the workflow frame selection. |
| `n_frames_used` | Frames actually analyzed after per-run stride. |
| `metrics` / `replicate_metrics` | Scalar metrics extracted from run summaries. |
| `metric_metadata` | Units and labels for scalar metrics. |

Common per-run payload fields include:

| Key | Meaning |
|-----|---------|
| `label` | SASA run label. |
| `target_selection` | Selection whose SASA is reported. |
| `context_selection` | Selection used as the blocking context. |
| `mean_sasa` | Mean SASA for the run in A^2. |
| `sem_sasa` | Standard error estimate for the run in A^2. |
| `sidecar_path` | Relative path to the NPZ sidecar for arrays. |
| `probe_radius_nm` | Probe radius used for the calculation. |
| `n_sphere_points` | Sphere point count used for the calculation. |

## Loading artifacts with `ArtifactStore`

Use the public MDAnalysis artifact API to inspect canonical artifacts:

```python
from pathlib import Path

from polyzymd.analyses.mda import ArtifactStore

replicate_store = ArtifactStore(Path("analysis/With Polymer/sasa/run_1"))
replicate = replicate_store.read_replicate_result()
print(replicate.payload["run_results"])

condition_store = ArtifactStore(Path("analysis/With Polymer/sasa/aggregated"))
condition = condition_store.read_condition_result()
print(condition.payload)
```

Sidecar NPZ files are referenced from the artifact `sidecars` list and from
payload fields such as `sidecar_path`. Treat sidecars as large validated data
files linked by the artifact, not as independently discovered cache files.

## Comparison output

`comparison/sasa/result.json` contains cross-condition statistics organized by
configured run label.

| Field | Description |
|-------|-------------|
| `metric` | Comparison metric name, currently `mean_sasa`. |
| `name` | Comparison study name. |
| `n_runs` | Number of configured SASA runs. |
| `run_labels` | Ordered list of run labels. |
| `control_label` | Configured control condition, when present. |
| `conditions` | Per-condition summaries with per-run means and SEMs. |
| `pairwise_comparisons` | Per-run pairwise statistics between conditions. |
| `anova_by_run` | Per-run ANOVA results when testable. |
| `ranking_by_run` | Condition ranking for each run. |
| `equilibration_time` | Equilibration cutoff used for the comparison. |

Pairwise comparison entries include:

| Field | Description |
|-------|-------------|
| `run_label` | SASA run being compared. |
| `condition_a`, `condition_b` | Conditions in the comparison. |
| `p_value`, `p_value_adjusted` | Raw and adjusted p-values when available. |
| `cohens_d` | Effect size when available. |
| `direction` | `shielding`, `exposure`, or `unchanged` based on SASA change. |
| `significant` | Whether the comparison passed the configured significance rule. |
| `percent_change` | Percent change from condition A to condition B. |
| `testable` | Whether the comparison had enough data for a statistical test. |
| `note` | Explanation for non-testable or special cases. |

## Normalized-control formula

Normalized comparison plots use the configured control condition as the
denominator:

```text
percent change = (condition_mean - control_mean) / control_mean * 100
```

For a shielding run such as `protein_with_polymer`, negative values indicate
lower SASA than the control and are consistent with polymer shielding. Positive
values indicate increased exposure relative to the control.

## Plot outputs

For each configured run label, SASA may generate these plots under
`comparison/sasa/plots/`:

| Plot output | Description |
|-------------|-------------|
| `sasa_comparison_<run>.png` | Mean SASA bar chart with SEM and replicate points. |
| `sasa_normalized_comparison_<run>.png` | Percent change relative to the configured control. |
| `sasa_timeseries_<run>.png` | Per-frame SASA traces summarized across conditions. |
| `sasa_profile_<run>.png` | Per-residue mean SASA profile across conditions. |

Time-axis plots assume uniformly saved frames. PolyzyMD maps frame index to time
as `frame_index * dt`; variable-timestep concatenated trajectories are not
supported.

## Units

| Quantity | Unit |
|----------|------|
| `probe_radius_nm` | nm |
| SASA values in outputs | A^2 |
| Time sidecar arrays | ns |

## See also

- {doc}`../tutorials/sasa_analysis` — guided shielding tutorial
- {doc}`../how_to/analysis_sasa_quickstart` — task recipes and commands
- {doc}`comparison_yaml` — comparison file schema
- {doc}`analysis_comparison_reference` — shared comparison and plotting behavior
