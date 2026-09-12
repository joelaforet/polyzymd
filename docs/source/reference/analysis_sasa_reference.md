# SASA Plugin Reference

This page is lookup documentation for the `sasa` analysis plugin: settings,
selection behavior, the observables it reports, and the files it writes.

For a guided workflow, see {doc}`../tutorials/sasa_analysis`. For practical
recipes and commands, see {doc}`../how_to/analysis_sasa_quickstart`.

## Plugin key

Top-level comparison YAML key: `plugins.sasa`.

```yaml
plugins:
  sasa:
    runs:
      - label: "protein_isolated"
        target_selection: "protein"
```

## Settings

### `plugins.sasa`

| Field | Type | Default | Constraints | Description |
|-------|------|---------|-------------|-------------|
| `runs` | list | required | at least one entry; labels must be unique | Contexts to measure. |
| `probe_radius_nm` | float | `0.14` | `> 0` | Shrake-Rupley probe radius in nanometers. |
| `n_sphere_points` | int | `960` | `>= 100` | Test points on each atom sphere. Higher is more accurate and slower. |
| `chunk_size` | int | `100` | `>= 1` | Frames sent to MDTraj per call. It bounds memory and does not change the numbers. |

### `runs` entries

| Field | Type | Default | Constraints | Description |
|-------|------|---------|-------------|-------------|
| `label` | string | required | non-empty | Name used in the observable names. |
| `target_selection` | string | required | non-empty | MDAnalysis selection for the atoms whose area is reported. |
| `context_selection` | string or null | `target_selection` | must contain the target | MDAnalysis selection for the atoms allowed to block the surface. |
| `stride` | int | `1` | `>= 1` | Deprecated since v1.3 and ignored. The framework resolves one frame window for every observable from `--eq-time`. Setting it raises a `DeprecationWarning`; it is rejected in v1.4. |

## Target and context behavior

A SASA run separates what is reported from what can block the surface.

| Selection | Behavior |
|-----------|----------|
| `target_selection` | Atoms whose area is summed and grouped into residues. |
| `context_selection` | Atoms present in the Shrake-Rupley calculation. |

When `context_selection` is omitted it equals `target_selection`, which reports
the target's own surface. The target must be a subset of the context, otherwise
the run raises `ReplicateError`, because a target atom outside the context would
be measured without its own neighbors. A selection matching no atoms raises
`SelectionError` rather than reporting an area of zero.

Examples:

| Goal | `target_selection` | `context_selection` |
|------|--------------------|---------------------|
| Whole-protein self-SASA | `protein` | `protein` or omitted |
| Protein SASA with polymer shielding | `protein` | `protein or resname SBM EGM` |
| Protein SASA with substrate present | `protein` | `protein or resname RBY` |
| Active-site SASA | `protein and (resid 77 or resid 156 or resid 262)` | `protein` |
| Monomer-specific shielding | `protein` | `protein or resname SBM` |

The project chain convention is A = protein, B = substrate, C = polymer, and
D+ = solvent/ions/other. Solvent and ions are excluded from every calculation
because they are never named in a target or context selection.

## Observables

Each configured run reports two observables.

| Name | Kind | Unit | Description |
|------|------|------|-------------|
| `sasa_<label>` | `mean_of_timeseries` | `A^2` | Total area of the target in that context, one value per frame. |
| `relative_sasa_<label>` | `profile` | `fraction` | Mean over frames of each target residue's area divided by the maximum accessible area of its residue type, indexed by residue ID. |

The maximum accessible areas are the empirical tripeptide values of Tien et al.
2013, held in `polyzymd.analyses.shared.aa_classification.MAX_ASA_TABLE`. A
target residue whose name is not in that table raises `ReplicateError`, so the
per-residue profile applies to standard amino acids only.

Residues are grouped by topology residue index, so two residues that share a
chain, a residue ID and a residue name stay separate.

## Canonical output paths

| Level | Path | Contents |
|-------|------|----------|
| Per replicate | `analysis/<condition>/sasa/run_<replicate>/result.json` | `ReplicateArtifact` whose `payload.observables` holds one reduced estimate per observable. |
| Per replicate sidecar | `analysis/<condition>/sasa/run_<replicate>/observables.npz` | Full per-frame series and per-residue vectors, one array per observable name. |
| Per condition | `analysis/<condition>/sasa/aggregated/result.json` | `ConditionArtifact` whose `payload.observables` holds one aggregate per observable. |
| Cross condition | `comparison/sasa/result.json` | `ComparisonArtifact` with the per-condition aggregates and the pairwise tests. |

## Artifact fields

A replicate estimate carries `name`, `kind`, `unit`, `n_frames`, the reduced
`value` for a time-series kind or the `profile` and `index` for a profile, and
the correlation diagnostics `statistical_inefficiency` and `n_eff`. The
diagnostics are reported, never used to shrink an error bar.

A condition aggregate carries `replicate_values`, `mean`, `sem`, `ci95_low`,
`ci95_high`, `ci_method`, `coverage` and `n_replicates` for a time-series kind,
and `profile_mean`, `profile_sem` and `index` for a profile. Every statistic is
computed across replicates, never across frames.

A comparison entry carries `control`, `condition`, `delta`, `percent_change`,
`test`, `p_value`, `p_adjusted`, `correction`, `cohens_d`, `significant`,
`testable` and `note`. Profiles are aggregated but not tested pairwise.

## Interpretation

A lower `sasa_protein_with_polymer` than `sasa_protein_isolated` in the same
condition means polymer atoms cover protein surface. Comparing
`sasa_protein_with_polymer` across conditions against the no-polymer control
answers whether one formulation shields more than another. The `delta` in a
comparison entry is in square angstrom and `percent_change` is relative to the
control mean.

## Figures

`compare run --plot` draws figures from the observable kind in the contract
runner. Until that work lands, a `sasa` run writes artifacts and the text report
but no figures.

## Units

| Quantity | Unit |
|----------|------|
| `probe_radius_nm` | nm |
| `sasa_<label>` | A^2 |
| `relative_sasa_<label>` | fraction of the maximum accessible area |

## References

- Shrake, A. & Rupley, J. A. (1973). Environment and exposure to solvent of
  protein atoms. Lysozyme and insulin. *Journal of Molecular Biology*, 79(2),
  351-371. doi:10.1016/0022-2836(73)90011-9
- Tien, M. Z. et al. (2013). Maximum allowed solvent accessibilities of residues
  in proteins. *PLoS ONE*, 8(11), e80635. doi:10.1371/journal.pone.0080635
- McGibbon, R. T. et al. (2015). MDTraj: a modern open library for the analysis
  of molecular dynamics trajectories. *Biophysical Journal*, 109(8), 1528-1532.
  doi:10.1016/j.bpj.2015.08.015

## See also

- {doc}`../tutorials/sasa_analysis` (guided shielding tutorial)
- {doc}`../how_to/analysis_sasa_quickstart` (task recipes and commands)
- {doc}`comparison_yaml` (comparison file schema)
- {doc}`analysis_comparison_reference` (shared comparison behavior)
