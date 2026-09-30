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
  contacts:
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
