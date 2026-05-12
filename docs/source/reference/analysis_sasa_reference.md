# SASA Plugin Reference

For a guided workflow, see {doc}`../tutorials/sasa_analysis`.

## Settings

Top-level plugin key: `plugins.sasa`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `runs` | `list[SASARunSettings]` | required | Named SASA computations to run |
| `probe_radius_nm` | `float` | `0.14` | Shrake-Rupley probe radius in nanometers |
| `n_sphere_points` | `int` | `960` | Sphere points for Shrake-Rupley SASA |
| `chunk_size` | `int` | `100` | Frames per chunk for memory-managed computation |

Each `runs` entry contains:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `label` | `str` | required | Human-readable run label; used in file names and plots |
| `target_selection` | `str` | required | Atoms whose SASA is reported |
| `context_selection` | `str \| null` | `target_selection` | Atoms included in the SASA context |
| `stride` | `int` | `1` | Analyze every Nth frame |

SASA uses a custom multi-run comparison model. Statistics are computed per run,
with pairwise tests, ANOVA by run, and optional normalized-to-control plots when
a control condition is configured.

## Output files

Per-replicate JSON and NPZ sidecars are written under
`analysis/<condition>/sasa/run_<replicate>/`. Aggregated results are written
under `analysis/<condition>/sasa/aggregated/`, and cross-condition statistics
are written to `comparison/sasa/result.json`.

## Plot outputs

For each configured run label, the plugin may generate:

| Plot output | Description |
|-------------|-------------|
| `sasa_comparison_<run>.png` | Mean SASA bar chart with replicate points |
| `sasa_normalized_comparison_<run>.png` | Percent change relative to control, when a control is configured |
| `sasa_timeseries_<run>.png` | Mean SASA over time with replicate shading |
| `sasa_profile_<run>.png` | Per-residue mean SASA profile |

Time-axis plots assume uniformly saved frames. PolyzyMD maps frame index to time
as `frame_index * dt`; variable-timestep concatenated trajectories are not
supported.
