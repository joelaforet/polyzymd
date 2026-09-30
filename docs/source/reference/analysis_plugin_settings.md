# Analysis Settings Reference

Every shipped analysis takes its settings on the command line, as
`polyzymd analyze NAME -c <config.yaml> --set key=value`, the value read as
YAML, or as the `settings=` mapping of `polyzymd.analyses.analyze` in Python.
The settings each analysis takes, and their defaults, are listed in
`FUNCTION_ANALYSES` in `src/polyzymd/analyses/protocols.py`; a setting an
analysis does not take is refused with the list of the ones it does. The
`plugins:` section of `comparison.yaml` is retired: a block for a shipped
analysis is ignored with a warning, as {ref}`comparison-yaml-retired` lists.

| Analysis | Settings |
|---|---|
| `rg` | {doc}`../how_to/analysis_rg_quickstart` |
| `rmsd` | {doc}`../how_to/analysis_rmsd_quickstart` |
| `rmsf`, `residue_rmsd` | {doc}`../how_to/analysis_rmsf_quickstart` |
| `sasa` | {doc}`../how_to/analysis_sasa_quickstart` |
| `secondary_structure` | {doc}`../how_to/analysis_secondary_structure_quickstart` |
| `contacts` | {doc}`../how_to/analysis_contacts_quickstart` |
| `native_contacts` | {doc}`../how_to/analysis_native_contacts_quickstart` |
| `hydrogen_bonds` | {doc}`../how_to/hydrogen_bonds` |
| `distances` | {doc}`../how_to/analysis_distances_quickstart` |

What each setting changes in the measurement is described with the function
it is passed to in {doc}`analysis_functions`. The catalytic triad is a routine
on the study API rather than an analysis; see
{doc}`../how_to/analysis_triad_quickstart`.

## Error bars in figures

The figures of `polyzymd analyze` that draw each condition's mean draw its 95
percent Student t confidence interval across replicates, with the
per-replicate values overlaid, and a footnote naming the interval, the number
of replicates and the production window; with one replicate per condition the
footnote says no interval is drawn. The difference figures draw the 95 percent
interval of each difference from the control, named in their footnote with
the test and the Benjamini-Hochberg family size. No setting changes the
interval.

## Universe loading (`pbc_policy`)

Every analysis reads its coordinates through `UniverseProvider.load_universe()`,
which `Study` calls once per replicate, or `TrajectoryLoader.load_universe()`.
Both accept a `pbc_policy` argument.

```{important}
`pbc_policy` is a Python argument, not a `polyzymd analyze` setting. Running
`polyzymd analyze` always loads with the default `"as_is"`. Code that calls
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

`UniverseProvenance` records how the coordinates were produced; the plugin
framework serializes it into every MDAnalysis replicate artifact under
`universe_policy.provenance`.

| Field | Type | Meaning |
|---|---|---|
| `pbc_policy` | `str` | The policy applied on load, `"as_is"` or `"make_whole"` |
| `topology_has_bonds` | `bool \| null` | Whether the loaded topology carries bonds. `null` before a universe has been loaded |
| `bond_source` | `str` | `"conect"` when bonds were read from the topology file, `"guessed"` when MDAnalysis inferred them, `"none"` when there are none |
| `trajectory_variant` | `str \| null` | Which trajectory the engine chose: `"centered"` for `prod_centered.xtc`, `"nojump"` for `prod_nojump.xtc`, `"raw"` otherwise. OpenMM segments are always `"raw"` |

`trajectory_variant` matters because the GROMACS job script post-processes the
production trajectory with `trjconv -pbc nojump` and then `-center -pbc mol -ur
compact`, and the engine prefers those files. OpenMM writes no such variant, so
the same analysis sees different coordinate semantics on the two engines. The
field records which one was read rather than leaving it implied by a filename.
