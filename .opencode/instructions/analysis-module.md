# Analysis module rules

## The real tree

Verified against `src/polyzymd/analyses/` on 2026-09-30, after the
plugin framework was removed in the v1.3 analysis refactor. Every file named here exists. If you add or delete a
module, update this list in the same commit.

```
src/polyzymd/analyses/
├── study.py             # Study, Condition, Replicate: replicates as MDAnalysis universes
├── study_file.py        # study.yaml: conditions, equilibration, analysis runs; Study("study.yaml")
├── results.py           # read_results: stored values and report.json without trajectories
├── user_functions.py    # a study's own functions (function: file.py:name), hashed by whole file
├── study_scaffold.py    # polyzymd study init: the study folder layout
├── study_git.py         # the study folder's git state, recorded in reports
├── study_metadata.py    # metadata: block, CITATION.cff and .zenodo.json
├── study_freeze.py      # polyzymd study freeze: manifest, checklist, tag, deposit/, warnings
├── schemas/             # study-1 and manifest-1 JSON Schemas, shipped and deposited
├── study_upload_guide.py # deposit/upload/, trajectories.csv and UPLOAD.md; uploads nothing
├── timeseries.py        # Study.timeseries, Study.per_replicate, Timeseries, ReplicateValues
├── functions.py         # Shipped measurements: radius_of_gyration, rmsd, pair_distance,
│                        # all_below, native_contacts, rmsf, rmsd_per_residue,
│                        # rms_decomposition, sasa, residue_sasa, dssp_occupancy,
│                        # residue_contacts, residue_occlusion, contact_events,
│                        # restricted_mean_lifetime, contact_lifetimes, event_lifetimes,
│                        # hbond_atoms, hydrogen_bonds, hbond_lifetimes,
│                        # residue_hbond_occupancy, residue_pair_hbond_occupancy, hbond_count
├── reference.py         # pz.reference: external, frame, average and centroid references
├── figures.py           # Figures drawn from stored study results
├── protocols.py         # polyzymd analyze: FUNCTION_ANALYSES, _analyze_*, ProtocolReport
├── universe.py          # UniverseProvider, UniverseProvenance, FileIdentity
├── identity.py          # compute_config_hash: recorded by every stored result, never change it
├── exceptions.py        # Typed analysis errors
├── shared/              # aa_classification, autocorrelation, centroid, diagnostics,
│                        # gromacs, inferential_statistics, loader, plotting,
│                        # selections, statistics, topology, window, groupings/
```

There is no plugin class, registry, discovery or scaffold: an analysis is a
function. `polyzymd analyze NAME -c ... --submit` writes and submits the SLURM
jobs of one command through `workflow/analysis_submit.py`.

rg, rmsd, rmsf, rmsd_per_residue, distances, sasa, secondary_structure,
native_contacts, contacts and hydrogen_bonds are not plugins. They are functions in
`functions.py` run through the study API, listed with their settings in
`protocols.FUNCTION_ANALYSES`, and `polyzymd analyze <name>` runs them through
`protocols._analyze_function`. `polyzymd analyze contacts` is dispatched by
`protocols._analyze_contacts`: `functions.residue_occlusion` for
`method=occlusion`, the default, and `functions.residue_contacts` for
`method=distance` give coverage, contact fractions and occluded area, and
`functions.contact_lifetimes` gives `--run mean_lifetime`, `lifetime_events`
and `censored_fraction`. `polyzymd analyze hydrogen_bonds` is dispatched by
`protocols._analyze_hydrogen_bonds`: `functions.hydrogen_bonds` gives each
summary's counts, `functions.hbond_lifetimes` its lifetimes,
`functions.residue_hbond_occupancy` its `_residues` and
`functions.residue_pair_hbond_occupancy` (with `labels="returned"`) its `_pairs`.

Retired names and commands: `catalytic_triad` is a routine on the study API
(`docs/source/how_to/analysis_triad_quickstart.md`), not an analysis:
`polyzymd analyze catalytic_triad` raises `ProtocolError`
(`protocols._refuse_retired`) pointing to that page and to
`polyzymd analyze distances --set pairs=...`. Any other name outside
`FUNCTION_ANALYSES` raises `ProtocolError` listing them. `polyzymd analyze -f
comparison.yaml` raises `ProtocolError` with the equivalent `-c` command
(`cli.analyze._refuse_comparison_file`, which reads the file with
`yaml.safe_load`). The hidden `polyzymd compare` and `polyzymd new-analysis`
commands (`cli/retired.py`) accept any arguments and exit 2. Messages for
retired things name the replacement, the docs page to read and what to point
an agent at.

## Adding a measurement

Write a new measurement as a function, not a plugin, and run it through the
study API; see `docs/source/explanation/analysis_api.md`.

- A per-frame function takes MDAnalysis `AtomGroup`s at one frame and returns
  one number. Run it with `study.timeseries(fn, pz.select(...), ...)`.
- A per-replicate function also takes `frames` and returns a number or one value
  per label, such as per residue. Run it with `study.per_replicate(fn, ...,
  labels=..., parts=...)`.
- Call MDAnalysis or MDTraj for the measurement; do not reimplement them.
- To run it from a study folder without shipping it, list it in `study.yaml`
  as `function: analyses/file.py:name` with `kind`, `selections` and
  `settings`; see `docs/source/how_to/study_yaml.md`.
- To ship it in `polyzymd analyze`, add it to `FUNCTION_ANALYSES` with its
  settings and a dispatcher in `protocols.py`, as the existing analyses do, with
  its figures, a quick-start page and a real-data parity check.

## Loading trajectories

`polyzymd.analyses.shared.loader.TrajectoryLoader` is the canonical universe
loader. It resolves topology and trajectory files for a replicate, checks
segment lineage, and builds the MDAnalysis universe through
`loader.open_universe`.

Bonds and charges come from the run's force field: OpenMM runs from
`<segment>_system.xml` (`enrich_universe_force_field`), GROMACS runs from
`prod.tpr`. When MDAnalysis cannot read the TPR's version (GROMACS 2026 with
MDAnalysis 2.10), `open_universe` warns and builds the universe from the
run's `.top` with `shared/gromacs.universe_from_gromacs_top`, laid out
exactly as MDAnalysis lays out a TPR; `tests/analyses/test_gromacs_topology.py`
checks that against a GROMACS 2025 TPR. GROMACS universes take PolyzyMD's
chain IDs (A protein, B substrate, C polymer) from the build's
`solvated_system.pdb` (`apply_build_chain_ids`). `UniverseProvenance.bond_source`
names the bond source: `system_xml`, `tpr`, `top`, `conect`, `guessed` or
`none`.

`polyzymd.analyses.universe.UniverseProvider` wraps `TrajectoryLoader`. It
takes a `SimulationConfig`, instantiates the loader lazily, and adds input
provenance (`UniverseProvenance`) to each load. `Study` uses one provider per
condition. Nothing else should build a `Universe` directly, and no code should
construct file paths by hand.

## Statistical contract

The replicate is the sampling unit for every cross-condition test and every
condition-level uncertainty. Equilibration is one global value applied
uniformly, and no diagnostic is allowed to select data. Every metric carries a
unit and a stated uncertainty. Invoke the `livecoms-check` skill in
`.claude/skills/` before committing analysis code.
