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
| `selection` | `str` | `"protein and name CA"` | MDAnalysis selection whose residues carry the profile |
| `alignment_selection` | `str` | `"protein and name CA"` | Selection superposed before the fluctuation is measured |
| `centroid_selection` | `str` | `"protein"` | Selection used to pick the representative frame in centroid mode |
| `reference_mode` | `str` | `"centroid"` | Reference mode: `centroid`, `average`, `frame`, or `external` |
| `reference_frame` | `int \| null` | `null` | Frame number used when `reference_mode: frame` (1-indexed) |
| `reference_file` | `str \| null` | `null` | External structure path used when `reference_mode: external` |

`rmsf` has no plot settings model.

## `catalytic_triad`

| Key | Type | Default | Description |
|---|---|---|---|
| `name` | `str` | `"catalytic_triad"` | Name of the triad/active-site definition |
| `pairs` | `list[PairSelection]` | required | Distance pairs to monitor, at least one |
| `threshold` | `float` | `3.5` | Contact threshold in Å |
| `description` | `str \| null` | `null` | Optional human-readable description |

`PairSelection` entries in `pairs`:

| Key | Type | Default | Description |
|---|---|---|---|
| `label` | `str` | required | Human-readable pair name |
| `selection_a` | `str` | required | First atom/point selection |
| `selection_b` | `str` | required | Second atom/point selection |

## `distances`

| Key | Type | Default | Description |
|---|---|---|---|
| `threshold` | `float \| null` | `3.5` | Global distance threshold (Å) for contact-style state analysis |
| `pairs` | `list[DistancePair]` | required | Distance pairs to monitor, at least one |
| `use_pbc` | `bool` | `true` | Use minimum-image PBC-aware distances |

The `align_trajectory` and `alignment_*` keys are accepted for one release,
ignored, and raise a `DeprecationWarning`.

`DistancePair` entries in `pairs`:

| Key | Type | Default | Description |
|---|---|---|---|
| `label` | `str` | required | Human-readable pair name |
| `selection_a` | `str` | required | First atom/point selection |
| `selection_b` | `str` | required | Second atom/point selection |
| `threshold` | `float \| null` | `null` | Per-pair threshold override (falls back to global `threshold`) |
| `below_label` | `str \| null` | `null` | Name of the below-threshold state; defaults to `"below <threshold> A"`. `above_label` is accepted for one release and ignored |

## `contacts`

| Key | Type | Default | Description |
|---|---|---|---|
| `protein_selection` | `str` | `"chainid A"` | MDAnalysis selection for protein atoms |
| `polymer_selection` | `str` | `"chainid C"` | MDAnalysis selection for polymer atoms |
| `cutoff` | `float` | `4.5` | Contact distance cutoff in Å |
| `polymer_types` | `list[str] \| null` | `null` | Restrict the polymer selection to these residue names |
| `heavy_atoms_only` | `bool` | `false` | Exclude hydrogens from both selections before the cutoff is applied |
| `allow_single_fragment_fallback` | `bool` | `false` | Put every polymer residue in chain 0 when the topology has no bonds, instead of raising `TopologyBondsMissingError` |
| `residence_time_edges_ns` | `list[float]` | `[0.0, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56, 5.12, 10.24, 20.48]` | Bin edges of the residence-time distribution in ns |

`grouping`, `compute_residence_times`, `protein_groups`, `protein_partitions`,
`fdr_alpha`, `min_effect_size` and `top_residues` are ignored with a
`DeprecationWarning` and rejected in v1.4. See
{doc}`analysis_contacts_reference`.

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

Outputs are four `fraction` observables (`ss_helix`, `ss_strand`, `ss_coil`,
`ss_unassigned`) and two `profile` observables (`helix_occupancy`,
`strand_occupancy`).

## `sasa`

| Key | Type | Default | Description |
|---|---|---|---|
| `runs` | `list[SASARun]` | required, non-empty, labels unique | Contexts to measure |
| `probe_radius_nm` | `float` | `0.14` | MDTraj Shrake-Rupley probe radius (nm) |
| `n_sphere_points` | `int` | `960` | MDTraj Shrake-Rupley sphere point count |
| `chunk_size` | `int` | `100` | Frames per MDTraj call; bounds memory and shifts areas by about 0.1 percent, so hold it fixed across a comparison |

`SASARun` entries in `runs`:

| Key | Type | Default | Description |
|---|---|---|---|
| `label` | `str` | required | Context name used in the observable names |
| `target_selection` | `str` | required | Selection whose area is reported |
| `context_selection` | `str \| null` | `null` | Selection allowed to block the surface (`null` defaults to `target_selection`) |
| `stride` | `int` | `1` | Deprecated since v1.3 and ignored; the framework resolves one frame window from `--eq-time` |

## `hydrogen_bonds`

| Key | Type | Default | Description |
|---|---|---|---|
| `groups` | `dict[str, str]` | `{protein: "chainid A", polymer: "chainid C"}` | Named MDAnalysis selections the summaries read |
| `summaries` | `list[HydrogenBondSummarySettings]` or mapping | `[{name: protein_polymer, between: [protein, polymer]}]` | Partitions to report; a mapping uses its keys as names |
| `distance_cutoff` | `float` | `3.0` | Donor-acceptor cutoff (A) |
| `angle_cutoff` | `float` | `150.0` | D-H...A angle cutoff (degrees) |
| `donor_acceptor_elements` | `tuple[str, ...]` | `["N", "O"]` | Elements allowed to donate and accept; capitalized and de-duplicated, `H` and unknown symbols rejected |
| `update_selections` | `bool` | `true` | Re-evaluate MDAnalysis's donor, hydrogen and acceptor selections each frame; group membership is always fixed at the start of the window |
| `allow_empty_groups` | `bool` | `false` | Raise `SelectionError` on an empty group; set `true` to report its summaries as zero |
| `top_n_pairs` | `int` | `15` | Residue pairs kept in each occupancy profile |
| `hydrogens_selection` | `str \| null` | `null` | Override for the hydrogen selection; `element H` by default |
| `composition` | mapping or `null` | `null` | Deprecated since v1.3 and ignored; query the event sidecar instead |
| `allow_overlapping_composition` | `bool` | `false` | Deprecated since v1.3 and ignored, with `composition` |
| `timestep_ps` | `float \| null` | `null` | Deprecated since v1.3 and ignored; the framework resolves the window |

`HydrogenBondSummarySettings` entries in `summaries`:

| Key | Type | Default | Description |
|---|---|---|---|
| `name` | `str` | required | Unique summary name; a mapping-form summary takes its key |
| `between` | `tuple[str, str] \| null` | `null` | Cross-group summary mode |
| `within` | `str \| null` | `null` | Intra-group summary mode |

Exactly one of `between` or `within` must be set for each summary, and every
group a summary names must be defined in `groups`.

Hydrogen detection uses MDAnalysis `HydrogenBondAnalysis` and requires explicit
hydrogens and element metadata. Donors and acceptors are
`(<group union>) and element <donor_acceptor_elements>` and hydrogens are
`(<group union>) and (element H)`. The trajectory loader infers elements for
GRO-like topologies that carry only atom types or atom names; a topology where
elements remain unavailable raises `SelectionError` rather than widening the
selection to every atom. Set
`hydrogens_selection` only for unusual explicit-hydrogen naming schemes.

## `rg`

| Key | Type | Default | Description |
|---|---|---|---|
| `runs` | `list[RgRunSettings]` | required, at least one | Named Rg runs to compute |

`RgRunSettings` entries in `runs`:

| Key | Type | Default | Description |
|---|---|---|---|
| `label` | `str` | required | Run label, slugged into the observable name |
| `selection` | `str` | required | MDAnalysis selection for Rg calculation |
| `calculation_mode` | `"selection" \| "fragments"` | `"selection"` | Whole-selection vs per-fragment Rg mode |
| `fragment_weighting` | `"equal" \| "mass"` | `"equal"` | Weighting of the per-frame mean over fragments (fragment mode) |
| `save_fragment_distribution` | `bool` | `true` | Report the fragment distribution as a profile over bins |
| `histogram_bins` | `int` | `50` | Bins in the fragment distribution |
| `histogram_range` | `[float, float]` | required when `save_fragment_distribution` is true | Range the fragment distribution covers, in A |
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
