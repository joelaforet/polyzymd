# Shipped analysis functions

These functions ship in `polyzymd.analyses.functions`. The per-frame functions
measure one frame and run through `Study.timeseries`, and the per-replicate
functions measure all production frames of one replicate and run through
`Study.per_replicate`, like functions you write yourself; see
{doc}`../explanation/analysis_api`.

## Per-frame functions

| Function | Arguments | Returns | Measurement | Command |
|---|---|---|---|---|
| `radius_of_gyration` | `atoms`, an `AtomGroup` | float, Å | MDAnalysis `AtomGroup.radius_of_gyration()`: mass-weighted, coordinates as loaded, molecules split across periodic boundaries not unwrapped | `polyzymd analyze rg` |
| `rmsd` | `atoms` and `reference`, `AtomGroup`s with the same atoms | float, Å | MDAnalysis `rms.rmsd(center=True, superposition=True)`: unweighted, `atoms` superposed on `reference` | `polyzymd analyze rmsd` |
| `pair_distance` | `atoms_a`, `atoms_b`, `mode_a` and `mode_b` (`single`, `midpoint`, `centroid` or `com`), `pbc` | float, Å | MDAnalysis `lib.distances.calc_bonds` between the two points, with the minimum image of the frame's box when `pbc` is true | `polyzymd analyze distances`, `polyzymd analyze catalytic_triad` |
| `all_below` | distance series, `thresholds` | 1 or 0 per frame | 1 for each frame in which every distance is strictly below its threshold; used with `Timeseries.transform` | `polyzymd analyze catalytic_triad` |
| `sasa` | `target` and `context`, `AtomGroup`s with every target atom in the context; `probe_radius_nm` (0.14) and `n_sphere_points` (960) | float, Å² | `mdtraj.shrake_rupley(mode="atom")` on the context atoms in a call of its own, with each atom's radius from MDTraj's element table plus the probe radius, summed over the target atoms; context atoms outside the target cover the target without being counted; periodic images not considered | `polyzymd analyze sasa` |

## Per-replicate functions

The three RMS functions each take `atoms`, the measured `AtomGroup`; `fit`, the `AtomGroup` superposed
on the reference; `reference`, the reference positions of `atoms | fit` from
`pz.reference` with the selection `"(<atoms>) or (<fit>)"`; and `frames`, the
production frame indices. Every frame is superposed on the reference by the
`fit` atoms, with the rotation of MDAnalysis `align.AlignTraj`, on a copy of
the coordinates. Each per-atom value is then averaged over the residue's atoms,
giving one value per residue in the order of `atoms.residues`.

| Function | Arguments | Returns | Measurement | Command |
|---|---|---|---|---|
| `rmsf` | `atoms`, `fit`, `reference`, `frames` | one value per residue, Å | MDAnalysis `rms.RMSF` of the superposed positions: the fluctuation of each atom about its mean position, as `gmx rmsf -o` gives | `polyzymd analyze rmsf` |
| `rms_deviation` | `atoms`, `fit`, `reference`, `frames` | one value per residue, Å | Root mean square deviation of each atom from its reference position, as `gmx rmsf -od` gives | `polyzymd analyze rms_deviation` |
| `rms_decomposition` | `atoms`, `fit`, `reference`, `frames` | six rows of one value per residue | Rows `RMS_PARTS = ("rms_deviation", "rmsf", "offset")`, the offset being the distance of each atom's mean position from its reference position, in Å; then rows `MS_PARTS = ("ms_deviation", "msf", "ms_offset")`, the means over each residue's atoms of the squared per-atom values, in Å². `ms_deviation = msf + ms_offset` for every residue | `polyzymd analyze rmsf`, `polyzymd analyze rms_deviation` |
| `residue_sasa` | `target`, `context` and `frames`, as `sasa` takes them, and the production frames | one value per target residue, Å² | Each frame's per-atom SASA as in `sasa`, one frame per MDTraj call, summed over each residue's target atoms and averaged over the frames | `polyzymd analyze sasa --run <context>_residues` |
| `dssp_occupancy` | `atoms`, whole residues, the production frames and `simplified` (default true) | one row per class, one value per residue, fraction of frames | `mdtraj.compute_dssp(simplified=simplified)` on every frame, 200 frames per call; the rows are helix, strand, coil and unassigned of `DSSP_SIMPLIFIED`, or with `simplified=False` the eight classes and unassigned of `DSSP_CLASSES`; one MDTraj chain per chain ID or segment | `polyzymd analyze secondary_structure` |
| `residue_contacts` | `protein` and `polymer`, `AtomGroup`s; the production frames; `cutoff` (4.0 Å), `types` (none) and `pbc` (true) | one row for the polymer, then one row per residue name in `types`, one value per protein residue, fraction of frames | On every frame, MDAnalysis `lib.distances.capped_distance` between the given polymer and protein atoms, with the minimum image of the frame's box when `pbc` is true; a residue is in contact when any of its atoms is within `cutoff` of a polymer atom, or for a row of `types`, of a polymer atom of that residue name | `polyzymd analyze contacts --set method=distance` |
| `residue_occlusion` | `protein` and `occluder`, `AtomGroup`s; the production frames; `threshold` (0.2), `types` (none), `max_asa` (`theoretical`), `pbc` (true), `probe_radius_nm` and `n_sphere_points` | rows `OCCLUSION_PARTS = ("contact_fraction", "exposed_fraction", "occluded_area", "exposed_area")`, then one contact-fraction row per residue name in `types`; one value per protein residue with a maximum ASA | Each residue's SASA with the protein alone and with the protein and the occluder, as `residue_sasa` computes it, one frame per MDTraj call, after each occluder molecule (bonded fragment) is moved whole to its periodic image nearest the protein when `pbc` is true; exposed when the SASA alone is at least `threshold` times the residue's maximum ASA of Tien et al. (2013), in contact when exposed and below that with the occluder; `occluded_area` is the mean of `max(0, alone - with)` in Å², `exposed_area` the mean SASA alone; a type row counts contact with only that residue name's occluder atoms present | `polyzymd analyze contacts` |
| `contact_lifetimes` | `protein` and `polymer`, `AtomGroup`s; the production frames; `method` (`occlusion`), `types` (none), `tolerance_ps` (0) and the method's options | rows `LIFETIME_PARTS = ("mean_lifetime", "n_events", "censored_fraction")`; one column for the polymer, then one per residue name in `types` | Each frame's contacts by `residue_occlusion`'s or `residue_contacts`'s rule; `contact_events` finds each residue's runs of frames in contact, after filling absences of at most `tolerance_ps` with MDAnalysis `lib.correlations.correct_intermittency`, a run including the first or last frame being censored; a run of `k` frames lasts `k` frame spacings; `mean_lifetime` is the Kaplan-Meier restricted mean in ns (`scipy.stats.ecdf` on `CensoredData`, area under the survival function up to the time the frames span), `nan` without events | `polyzymd analyze contacts --run mean_lifetime` |

The comparison with `gmx rmsf` on real trajectories, and the script that reruns
it, are in {doc}`../explanation/analysis_rmsf_verification`.

`polyzymd analyze rg` measures `radius_of_gyration` of `--set selection=...`
(default `protein`) on every production frame, reduces each replicate to its
mean, and compares every condition with the first by Welch's t test. See
{doc}`../how_to/analysis_rg_quickstart`.

`polyzymd analyze rmsd` measures `rmsd` of `--set selection=...` (default
`protein and name CA`) against
`pz.reference(reference_mode, selection, frame=reference_frame, file=reference_file, alignment=alignment_selection)`,
with `centroid` as the default mode. `reference_frame` counts production frames
from 1 after the equilibration window, and `alignment_selection` is the set of
atoms superposed to build the `average` and `centroid` references. Each frame
is superposed on the `selection` atoms themselves. See
{doc}`../how_to/analysis_rmsd_quickstart`.

`polyzymd analyze distances --set pairs=pairs.yaml` reads a YAML or JSON list of
pairs, each with `label`, `selection_a`, `selection_b` and optionally
`threshold` and `below_label`. It reports `<label>`, the mean distance, and
`<label> <below_label>`, the fraction of frames below the pair's threshold (the
analysis `threshold`, 3.5 Å by default, when the pair sets none).
`polyzymd analyze catalytic_triad` adds `simultaneous`, the fraction of frames in
which every pair is below its threshold, computed from the stored distances,
and reports it first. Pick a result with `--run`. See
{doc}`../how_to/analysis_distances_quickstart` and
{doc}`../how_to/analysis_triad_quickstart`.

`polyzymd analyze rmsf` and `polyzymd analyze rms_deviation` run
`rms_decomposition` once per replicate, with `selection` measured and
`alignment_selection` fitted (both `protein and name CA` by default) against
`pz.reference(reference_mode, "(selection) or (alignment_selection)", frame=reference_frame, file=reference_file, alignment=alignment_selection)`.
`reference_mode` defaults to `centroid` for `rmsf`; for `rms_deviation` it is
`external` when `reference_file` is set and `centroid` otherwise. Each
replicate's `core_<part>` value is the square root of the mean over the core
residues of the part's mean-square row, so
`core_rms_deviation² = core_rmsf² + core_offset²`. The core is the residues of
`selection` that `--set core=...` also selects, by default all of them. Each
entry of `--set regions='{name: selection}'` gives `<name>_<part>` the same way.
`mean_<part>` is the plain mean of the per-residue values, and `--run
rmsf`, `offset` or `rms_deviation` compares the profile residue by residue. The
default result is `core_rmsf` for `rmsf` and `core_rms_deviation` for
`rms_deviation`. The residues of the core and of each region are recorded
under `provenance.settings`. See {doc}`../how_to/analysis_rmsf_quickstart`.

`polyzymd analyze sasa` measures `sasa` of `--set target=...` (default
`protein`) in each context of `--set contexts='{name: selection}'` (default the
target alone, named `isolated`), with `probe_radius_nm` and `n_sphere_points`.
`--run <name>` reports the mean over production frames of the total, and
`--run <name>_residues` runs `residue_sasa` and compares every residue. Only the
chosen result is measured. The default is the first context's total. See
{doc}`../how_to/analysis_sasa_quickstart`. MDTraj is given one frame per call
because MDTraj 1.11.1 returns about 0.1 percent too much area for the later
frames each thread computes in a call; see {doc}`../explanation/analysis_sasa_verification`.

`polyzymd analyze secondary_structure` runs `dssp_occupancy` once per
replicate on `--set selection=...` (default `protein`, whole residues), in the
scheme of `--set scheme=...`: `simplified` (default) or `full`. `--run <class>`
reports a class's mean over the residues, the fraction of residue-frames in it,
and `--run <class>_residues` compares every residue. The default is the
scheme's first class, `helix` or `alpha_helix`. A warning names the replicates
with unassigned residues. See {doc}`../how_to/analysis_secondary_structure_quickstart`.

`polyzymd analyze contacts` runs, once per replicate between `--set
protein_selection=...` (default `chainid A`) and `polymer_selection` (default
`chainid C`, narrowed to the residue names of `polymer_types` when set),
`residue_occlusion` with `threshold`, `max_asa`, `probe_radius_nm` and
`n_sphere_points` for `method=occlusion` (default), or `residue_contacts` with
`cutoff` for `method=distance`, on heavy atoms only when `heavy_atoms` is true
(default). `use_pbc` sets `pbc`, and every polymer residue name of the
control's first replicate gets a row. `coverage`, the default, is the fraction
of measured residues with a contact fraction above zero;
`mean_contact_fraction`, `<type>_contact_fraction`, `<class>_contact_fraction`
for each amino-acid class of `ProteinAAClassification` and
`<region>_contact_fraction` for each entry of `--set regions='{name:
selection}'` are means of contact fractions over those residues;
`contact_fraction_residues` and `<type>_contact_fraction_residues` compare
every residue. For occlusion, `occluded_area` is the sum over residues of the
mean occluded area, `occlusion_fraction` is the summed occluded area over the
summed SASA alone, and `occluded_area_residues` compares every residue. The
selections, the types found, the residues without a maximum ASA and the
residues of each class and region are recorded under `provenance.settings`.
`--run mean_lifetime`, `<type>_mean_lifetime`, `lifetime_events` and `censored_fraction` run `contact_lifetimes` with the same method and `tolerance_ps`, and are measured only when chosen; see {doc}`../explanation/analysis_contact_lifetimes`. See
{doc}`../how_to/analysis_contacts_quickstart` and
{doc}`../explanation/analysis_contacts_verification`.

## Figures

`polyzymd analyze` writes these figures to `<output-dir>/figures/<analysis>/`
and records the folder under `output_paths.figures` in the JSON report;
`--no-plots`, or `plots=False` in `analyze`, skips them.

| Analysis | Figures |
|---|---|
| `rg` | `rg_timeseries` (every replicate against time), `rg_comparison` (condition means with every replicate value), `rg_distribution` (per-frame distribution) |
| `rmsd` | `rmsd_timeseries`, `rmsd_comparison` |
| `rmsf`, `rms_deviation` | `rmsf_profile`, `offset_profile` and `rms_deviation_profile` (each replicate and each condition's mean with its interval, `highlight_residues` marked), `rms_decomposition` (the three profiles of each condition together), `rmsf_comparison` (the three core values), and with several conditions `rmsf_difference`, `offset_difference` and `rms_deviation_difference` (each condition minus the control at every residue, with the interval of the difference and the significant residues marked) |
| `sasa` | For a total, `sasa_timeseries_<name>`, `sasa_comparison_<name>` and `sasa_distribution_<name>`; for `<name>_residues`, `sasa_profile_<name>` and with several conditions `sasa_difference_<name>` |
| `secondary_structure` | `ss_content_bars` (every class of the scheme but unassigned); for a total `ss_<name>_comparison`; for `<name>_residues`, `ss_<name>_profile`, `ss_classes_<name>` and with several conditions `ss_<name>_difference` |
| `contacts` | `contacts_class_bars` (each amino-acid class's mean contact fraction); for a one-value result `contacts_<run>_comparison`; for a residue result, `contacts_<name>_profile` and with several conditions `contacts_<name>_difference` |
| `distances` | `distance_kde_<pair>` for each pair with its threshold, `distance_fraction_<result>` for each fraction below threshold, and the grouped `distance_threshold_bars` (every pair's fraction) and `distance_kde_panel` (one panel per pair) |
| `catalytic_triad` | `triad_kde_<pair>` for each pair with its threshold, `triad_fraction_<result>` for `simultaneous` and each pair's fraction, and the grouped `triad_threshold_bars` (each pair's fraction, then all pairs) and `triad_kde_panel` (one panel per pair) |
