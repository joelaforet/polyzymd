# Analysis Plugin Settings Reference

This page is a complete YAML key reference for plugin settings under
`comparison.yaml`:

```yaml
plugins:
  <plugin_name>:
    ...settings...
```

Use this as a lookup table for field names, types, defaults, and meanings.

FDR thresholds for comparison workflows are configured through top-level
`defaults.fdr_alpha` in `comparison.yaml` where supported, not through
plugin-local settings unless a plugin explicitly lists its own `fdr_alpha` field
below.

## Plot settings shared by every plugin

Each plugin's `plot_settings` block inherits this key.

| Key | Type | Default | Description |
|---|---|---|---|
| `error_bar` | `"ci95" \| "sem"` | `"ci95"` | Interval drawn on comparison bars and shaded bands. `ci95` draws the 95 percent Student t confidence interval across replicates. `sem` draws one standard error, which at `n = 3` is 4.3 times narrower |

```yaml
plugins:
  rmsd:
    plot_settings:
      error_bar: ci95
```

Whichever value is set, the figure carries a footnote naming the interval, the
number of replicates and the production window, and per-replicate points stay
overlaid on the bars. `hydrogen_bonds` has no plot settings model, so its
figures always use the default.

## Universe loading (`pbc_policy`)

Every plugin reads its coordinates through `TrajectoryLoader.load_universe()`
and `UniverseProvider.load_universe()`. Both accept a `pbc_policy` argument.

```{important}
`pbc_policy` is a Python argument, not yet a `comparison.yaml` key. Running
`polyzymd compare run` always loads with the default `"as_is"`. Code that calls
the loader or the universe provider directly can pass `"make_whole"` today, and
the policy in force is recorded in provenance either way.
```

| Value | Meaning |
|---|---|
| `"as_is"` (default) | Coordinates are used exactly as the trajectory stores them. No unwrap, centering, or make-whole step runs |
| `"make_whole"` | An MDAnalysis `unwrap` transformation is registered on the protein and polymer selection, so molecules split across a periodic boundary are rejoined before any measurement reads them |

`"make_whole"` walks the bond graph, so it raises
`polyzymd.analyses.exceptions.TopologyBondsMissingError` when the topology has
no bonds. The default stays `"as_is"`, so loading behaviour does not change
unless you ask for it.

The selection unwrapped by `"make_whole"` defaults to
`not (water or resname NA CL K MG ZN SOD CLA POT NA+ CL-)`, which is everything
that is not solvent or a monatomic ion.

### Universe provenance fields

`UniverseProvenance` (serialized into every MDAnalysis replicate artifact under
`universe_policy.provenance`) records how the coordinates were produced.

| Field | Type | Meaning |
|---|---|---|
| `pbc_policy` | `str` | The policy applied on load, `"as_is"` or `"make_whole"` |
| `topology_has_bonds` | `bool \| null` | Whether the loaded topology carries bonds. `null` before a universe has been loaded |
| `bond_source` | `str` | `"conect"` when bonds were read from the topology file, `"guessed"` when MDAnalysis inferred them, `"none"` when there are none |
| `trajectory_variant` | `str \| null` | Which trajectory the engine chose: `"centered"` for `prod_centered.xtc`, `"nojump"` for `prod_nojump.xtc`, `"raw"` otherwise. OpenMM segments are always `"raw"` |

`trajectory_variant` matters because the GROMACS job script post-processes the
production trajectory with `trjconv -pbc nojump` and then `-center -pbc mol -ur
compact`, and the engine prefers those files. OpenMM writes no such variant, so
the same plugin sees different coordinate semantics on the two engines. The
field records which one was read rather than leaving it implied by a filename.
## Plot settings shared by every plugin

Each plugin's `plot_settings` block inherits this key.

| Key | Type | Default | Description |
|---|---|---|---|
| `error_bar` | `"ci95" \| "sem"` | `"ci95"` | Interval drawn on comparison bars and shaded bands. `ci95` draws the 95 percent Student t confidence interval across replicates. `sem` draws one standard error, which at `n = 3` is 4.3 times narrower |

```yaml
plugins:
  rmsf:
    plot_settings:
      error_bar: ci95
```

Whichever value is set, the figure carries a footnote naming the interval, the
number of replicates and the production window, and per-replicate points stay
overlaid on the bars. `hydrogen_bonds` has no plot settings model, and neither
does a plugin written against the observable contract such as `rmsd`, so their
figures always use the default.

## `rmsf`

| Key | Type | Default | Description |
|---|---|---|---|
| `selection` | `str` | `"protein and name CA"` | MDAnalysis selection used for RMSF calculation |
| `reference_mode` | `str` | `"centroid"` | Reference mode: `centroid`, `average`, `frame`, or `external` |
| `reference_frame` | `int \| null` | `null` | Frame number used when `reference_mode: frame` (1-indexed) |
| `reference_file` | `str \| null` | `null` | External PDB path used when `reference_mode: external` |
| `alignment_selection` | `str` | `"protein and name CA"` | Selection used for trajectory alignment |
| `centroid_selection` | `str` | `"protein"` | Selection used to find centroid frame |

## `catalytic_triad`

| Key | Type | Default | Description |
|---|---|---|---|
| `name` | `str` | `"catalytic_triad"` | Name of the triad/active-site definition |
| `pairs` | `list[TriadPairSettings]` | required | Distance pairs to monitor |
| `threshold` | `float` | `3.5` | Contact threshold in Å |
| `description` | `str \| null` | `null` | Optional human-readable description |

`TriadPairSettings` entries in `pairs`:

| Key | Type | Default | Description |
|---|---|---|---|
| `label` | `str` | required | Human-readable pair name |
| `selection_a` | `str` | required | First atom/point selection |
| `selection_b` | `str` | required | Second atom/point selection |

## `distances`

| Key | Type | Default | Description |
|---|---|---|---|
| `threshold` | `float \| null` | `3.5` | Global distance threshold (Å) for contact-style state analysis |
| `pairs` | `list[DistancePairSettings]` | `[]` (must be non-empty) | Distance pairs to monitor |
| `use_pbc` | `bool` | `true` | Use minimum-image PBC-aware distances |
| `align_trajectory` | `bool` | `false` | Deprecated and ignored since 1.3.0; setting it to `true` raises a `DeprecationWarning` |
| `alignment_selection` | `str` | `"protein and name CA"` | Deprecated and ignored since 1.3.0 |
| `alignment_mode` | `str` | `"centroid"` | Deprecated and ignored since 1.3.0; still validated as `centroid`, `average`, or `frame` |
| `alignment_frame` | `int \| null` | `null` | Deprecated and ignored since 1.3.0 |

`DistancePairSettings` entries in `pairs`:

| Key | Type | Default | Description |
|---|---|---|---|
| `label` | `str` | required | Human-readable pair name |
| `selection_a` | `str` | required | First atom/point selection |
| `selection_b` | `str` | required | Second atom/point selection |
| `threshold` | `float \| null` | `null` | Per-pair threshold override (falls back to global `threshold`) |
| `below_label` | `str \| null` | `null` | Display label for below-threshold state |
| `above_label` | `str \| null` | `null` | Display label for above-threshold state |

## `contacts`

| Key | Type | Default | Description |
|---|---|---|---|
| `polymer_selection` | `str` | `"chainid C"` | MDAnalysis selection for polymer atoms |
| `protein_selection` | `str` | `"chainid A"` | MDAnalysis selection for protein atoms |
| `cutoff` | `float` | `4.5` | Contact cutoff distance in Å |
| `polymer_types` | `list[str] \| null` | `null` | Optional polymer residue-name filter |
| `grouping` | `str` | `"aa_class"` | Grouping mode: `aa_class`, `secondary_structure`, or `none` |
| `compute_residence_times` | `bool` | `true` | Compute aggregate residence-time summaries and plots; per-replicate contact events remain stored when disabled |
| `allow_single_fragment_fallback` | `bool` | `false` | Put every polymer residue in chain 0 when the topology has no bonds, instead of raising `TopologyBondsMissingError` |
| `protein_groups` | `dict[str, list[int]] \| null` | `null` | Custom residue groups, e.g. `{name: [resid, ...]}` |
| `protein_partitions` | `dict[str, list[str]] \| null` | `null` | Partition definitions over custom `protein_groups` |
| `fdr_alpha` | `float` | `0.05` | Benjamini-Hochberg false-discovery-rate alpha |
| `min_effect_size` | `float` | `0.5` | Minimum Cohen's d highlighted in output |
| `top_residues` | `int` | `10` | Number of top-contact residues shown in summaries |

## `secondary_structure`

| Key | Type | Default | Description |
|---|---|---|---|
| `chain_id` | `str` | `"A"` | Protein chain letter to analyze (PolyzyMD convention: chain A) |
| `selection` | `str \| null` | `null` | Explicit MDAnalysis protein-residue selection. Overrides `chain_id` when set |

The default is `protein and chainid A` for PDB/PolyzyMD chain-convention
compatibility. GROMACS `.gro` topologies may not preserve chain IDs; use
`selection: "protein"`, `selection: "protein and resid 1:269"`, or
`selection: "protein and resindex 0:268"` when chain IDs are unavailable.
DSSP requires complete residues; do not use CA-only selections such as
`protein and name CA`.

## `sasa`

| Key | Type | Default | Description |
|---|---|---|---|
| `runs` | `list[SASARunSettings]` | `[]` (must be non-empty) | SASA runs to compute |
| `probe_radius_nm` | `float` | `0.14` | MDTraj Shrake-Rupley probe radius (nm) |
| `n_sphere_points` | `int` | `960` | MDTraj Shrake-Rupley sphere point count |
| `chunk_size` | `int` | `100` | Frames per chunk for SASA computation |

`SASARunSettings` entries in `runs`:

| Key | Type | Default | Description |
|---|---|---|---|
| `label` | `str` | required | Human-readable run label |
| `target_selection` | `str` | required | Selection whose SASA is reported |
| `context_selection` | `str \| null` | `null` | Environment/context selection for SASA computation (`null` defaults to `target_selection`) |
| `stride` | `int` | `1` | Frame stride (1 means every frame) |

## `hydrogen_bonds`

| Key | Type | Default | Description |
|---|---|---|---|
| `groups` | `dict[str, str]` | `{protein: "chainid A", polymer: "chainid C"}` | Named MDAnalysis selections used by summaries |
| `summaries` | `list[HydrogenBondSummarySettings]` | `[{name: protein_polymer, between: [protein, polymer]}]` | Summary definitions to compute |
| `distance_cutoff` | `float` | `3.0` | Donor-acceptor cutoff (Å) |
| `angle_cutoff` | `float` | `150.0` | D-H...A angle cutoff (degrees) |
| `update_selections` | `bool` | `true` | Re-evaluate selections each frame |
| `top_n_pairs` | `int` | `15` | Number of top residue pairs reported |
| `allow_empty_groups` | `bool` | `false` | Raise `SelectionError` on an empty group; set `true` to warn and skip instead |
| `donor_acceptor_elements` | `tuple[str, ...]` | `["N", "O"]` | Elements allowed to act as donors and acceptors |
| `allow_overlapping_composition` | `bool` | `false` | Allow overlapping composition partitions (otherwise raise) |
| `composition` | `HydrogenBondCompositionSettings \| null` | `null` | Optional partitioning for composition analysis |
| `hydrogens_selection` | `str \| null` | `null` | Advanced explicit-hydrogen selection override for unusual atom names |
| `timestep_ps` | `float \| null` | `null` | Optional timestep override (ps) for time-axis plots |

Time-axis plots assume uniformly saved frames. PolyzyMD maps frame index to time
as `frame_index * timestep_ps`; variable-timestep concatenated trajectories are
not supported.

`HydrogenBondSummarySettings` entries in `summaries`:

| Key | Type | Default | Description |
|---|---|---|---|
| `name` | `str` | required | Unique summary name |
| `between` | `tuple[str, str] \| null` | `null` | Cross-group summary mode |
| `within` | `str \| null` | `null` | Intra-group summary mode |

Exactly one of `between` or `within` must be set for each summary.

Hydrogen detection uses MDAnalysis `HydrogenBondAnalysis` and requires explicit
hydrogens and element metadata. Donors and acceptors are
`(<group union>) and element <donor_acceptor_elements>` and hydrogens are
`(<group union>) and (element H)`. PolyzyMD infers missing elements for GRO-like
topologies when atom types or atom names are conservative enough. If elements
remain unavailable, or if the donor and acceptor selection matches no atoms, the
plugin raises `SelectionError` unless `allow_empty_groups` is true. Set
`hydrogens_selection` only for unusual explicit-hydrogen naming schemes.

`HydrogenBondCompositionSettings`:

| Key | Type | Default | Description |
|---|---|---|---|
| `partitions` | `dict[str, str]` | `{}` | Named composition partitions as MDAnalysis selections |

## `rg`

| Key | Type | Default | Description |
|---|---|---|---|
| `runs` | `list[RgRunSettings]` | `[]` (must be non-empty) | Named Rg runs to compute |

`RgRunSettings` entries in `runs`:

| Key | Type | Default | Description |
|---|---|---|---|
| `label` | `str` | required | Human-readable run label |
| `selection` | `str` | required | MDAnalysis selection for Rg calculation |
| `calculation_mode` | `"selection" \| "fragments"` | `"selection"` | Whole-selection vs fragment-reduced Rg mode |
| `fragment_weighting` | `"equal" \| "mass"` | `"equal"` | Fragment reduction weighting (fragment mode) |
| `save_fragment_distribution` | `bool` | `true` | Save per-fragment distribution sidecar outputs |
| `histogram_bins` | `int` | `50` | Histogram bins for fragment distribution summaries |
| `allow_single_fragment_fallback` | `bool` | `false` | Measure the whole selection as one fragment when the topology has no bonds, instead of raising `TopologyBondsMissingError` |

## `rmsd`

| Key | Type | Default | Description |
|---|---|---|---|
| `runs` | `list[RMSDRunSettings]` | required, at least one | Named RMSD runs to compute |

`RMSDRunSettings` entries in `runs`:

| Key | Type | Default | Description |
|---|---|---|---|
| `label` | `str` | required | Human-readable run label |
| `selection` | `str` | `"protein and name CA"` | Selection used for RMSD calculation |
| `alignment_selection` | `str` | `"protein and name CA"` | Selection used for alignment |
| `reference_mode` | `str` | `"centroid"` | Reference mode: `centroid`, `average`, `frame`, or `external` |
| `reference_frame` | `int` | `0` | 0-indexed frame index for `reference_mode: frame` |
| `reference_file` | `str \| null` | `null` | External PDB path for `reference_mode: external` |
| `centroid_selection` | `str \| null` | `null` | Optional centroid-mode selection (falls back to `alignment_selection`) |
