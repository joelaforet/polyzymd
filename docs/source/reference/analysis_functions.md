# Shipped analysis functions

These functions ship in `polyzymd.analyses.functions`. The per-frame functions
measure one frame and run through `Study.timeseries`, and the per-replicate
functions measure all production frames of one replicate and run through
`Study.per_replicate`, like functions you write yourself; see
{doc}`study_api`.

## Per-frame functions

| Function | Arguments | Returns | Measurement | Command |
|---|---|---|---|---|
| `radius_of_gyration` | `atoms`, an `AtomGroup` | float, Å | MDAnalysis `AtomGroup.radius_of_gyration()`: mass-weighted, coordinates as loaded, molecules split across periodic boundaries not unwrapped | `polyzymd analyze rg` |
| `rmsd` | `atoms` and `reference`, `AtomGroup`s with the same atoms | float, Å | MDAnalysis `rms.rmsd(center=True, superposition=True)`: unweighted, `atoms` superposed on `reference` | `polyzymd analyze rmsd` |
| `pair_distance` | `atoms_a`, `atoms_b`, `mode_a` and `mode_b` (`single`, `midpoint`, `centroid` or `com`), `pbc` | float, Å | MDAnalysis `lib.distances.calc_bonds` between the two points, with the minimum image of the frame's box when `pbc` is true | `polyzymd analyze distances` |
| `all_below` | distance series, `thresholds` | 1 or 0 per frame | 1 for each frame in which every distance is strictly below its threshold; used with `Timeseries.transform` | `polyzymd analyze distances`, for each pair's fraction below its threshold |
| `native_contacts` | `atoms` and `reference`, `AtomGroup`s with the same atoms; optional `region`; `radius` (4.5 Å), `min_separation` (3), `beta` (5 Å⁻¹), `lambda_constant` (1.8), `pbc` | float, from 0 to 1 | Native pairs: pairs of `atoms` more than `min_separation` residues apart whose `reference` positions are closer than `radius`, found once per reference without periodic images; on each frame, MDAnalysis `analysis.contacts.soft_cut_q` of their distances (`lib.distances.calc_bonds`, minimum image when `pbc` is true) and reference distances, the mean over pairs of `1/(1 + exp(beta (r - lambda_constant r0)))`; with `region`, only pairs with an atom in it. The defaults on heavy atoms are the definition of Best, Hummer and Eaton (2013) | `polyzymd analyze native_contacts` |
| `hbond_count` | `group_a`, optional `group_b`; `d_a_cutoff` (3.5 Å), `d_h_a_angle_cutoff` (150°), optional `donors`, `hydrogens` and `acceptors` | float, number of hydrogen bonds | MDAnalysis `HydrogenBondAnalysis` on the current frame alone, with the atoms and rules of `hydrogen_bonds` below; suits a few chosen atoms, such as the hydrogen bonds of a catalytic triad | none; the routine of {doc}`../how_to/analysis_triad_quickstart` |
| `sasa` | `target` and `context`, `AtomGroup`s with every target atom in the context; `probe_radius_nm` (0.14) and `n_sphere_points` (960) | float, Å² | `mdtraj.shrake_rupley(mode="atom")` on the context atoms in a call of its own, with each atom's radius from MDTraj's element table plus the probe radius, summed over the target atoms; context atoms outside the target cover the target without being counted; periodic images not considered | `polyzymd analyze sasa` |

## Per-replicate functions

The three RMS functions each take `atoms`, the measured `AtomGroup`; `fit`, the `AtomGroup` superposed
on the reference; `reference`, the reference positions of `atoms | fit` from
`pz.reference` with the selection `"(<atoms>) or (<fit>)"`; and `frames`, the
production frame indices. Every frame is superposed on the reference by the
`fit` atoms, with the rotation of MDAnalysis `align.AlignTraj`, on a copy of
the coordinates. Each residue's value is the square root of the mean over its
atoms of the squared per-atom values, as `gmx rmsf -res` combines atoms with
equal masses, giving one value per residue in the order of `atoms.residues`.

| Function | Arguments | Returns | Measurement | Command |
|---|---|---|---|---|
| `rmsf` | `atoms`, `fit`, `reference`, `frames` | one value per residue, Å | MDAnalysis `rms.RMSF` of the superposed positions: the fluctuation of each atom about its mean position, as `gmx rmsf -o` gives | `polyzymd analyze rmsf` |
| `rmsd_per_residue` | `atoms`, `fit`, `reference`, `frames` | one value per residue, Å | Root mean square over frames of each atom's distance from its reference position, as `gmx rmsf -od` gives; `rmsd_per_residue² = rmsf² + offset²` per atom. Not `rmsd`, which is one value per frame, the root mean square over atoms | `polyzymd analyze rmsd_per_residue` |
| `rms_decomposition` | `atoms`, `fit`, `reference`, `frames` | six rows of one value per residue | Rows `RMS_PARTS = ("rmsd_per_residue", "rmsf", "offset")`, the offset being the distance of each atom's mean position from its reference position, in Å; then rows `MS_PARTS = ("ms_deviation", "msf", "ms_offset")`, the means over each residue's atoms of the squared per-atom deviation, RMSF and offset, in Å². `ms_deviation = msf + ms_offset` for every residue | `polyzymd analyze rmsf`, `polyzymd analyze rmsd_per_residue` |
| `residue_sasa` | `target`, `context` and `frames`, as `sasa` takes them, and the production frames | one value per target residue, Å² | Each frame's per-atom SASA as in `sasa`, one frame per MDTraj call, summed over each residue's target atoms and averaged over the frames | `polyzymd analyze sasa --run <context>_residues` |
| `dssp_occupancy` | `atoms`, whole residues, the production frames and `simplified` (default true) | one row per class, one value per residue, fraction of frames | `mdtraj.compute_dssp(simplified=simplified)` on every frame, 200 frames per call; the rows are helix, strand, coil and unassigned of `DSSP_SIMPLIFIED`, or with `simplified=False` the eight classes and unassigned of `DSSP_CLASSES`; one MDTraj chain per chain ID or segment | `polyzymd analyze secondary_structure` |
| `residue_contacts` | `protein` and `polymer`, `AtomGroup`s; the production frames; `cutoff` (4.0 Å), `types` (none) and `pbc` (true) | one row for the polymer, then one row per residue name in `types`, one value per protein residue, fraction of frames | On every frame, MDAnalysis `lib.distances.capped_distance` between the given polymer and protein atoms, with the minimum image of the frame's box when `pbc` is true; a residue is in contact when any of its atoms is within `cutoff` of a polymer atom, or for a row of `types`, of a polymer atom of that residue name | `polyzymd analyze contacts --set method=distance` |
| `residue_occlusion` | `protein` and `occluder`, `AtomGroup`s; the production frames; `exposed_threshold` (0.2), `buried_threshold` (0.2), `types` (none), `max_asa` (`theoretical`), `pbc` (true), `probe_radius_nm` and `n_sphere_points` | rows `OCCLUSION_PARTS = ("contact_fraction", "exposed_fraction", "occluded_area", "exposed_area")`, then one contact-fraction row per residue name in `types`; one value per protein residue with a maximum ASA | Each residue's SASA with the protein alone and with the protein and the occluder, as `residue_sasa` computes it, one frame per MDTraj call, after each occluder molecule (bonded fragment) is moved whole to its periodic image nearest the protein when `pbc` is true; exposed when the SASA alone is at least `exposed_threshold` times the residue's maximum ASA of Tien et al. (2013), in contact when exposed and the SASA with the occluder is below `buried_threshold` times it and lower than alone; `occluded_area` is the mean of `max(0, alone - with)` in Å², `exposed_area` the mean SASA alone; a type row counts contact with only that residue name's occluder atoms present | `polyzymd analyze contacts` |
| `hbond_atoms` | `atoms`, an `AtomGroup` | the donatable hydrogens and the acceptors, two `AtomGroup`s | Reads the universe's bonds, which PolyzyMD takes from the run's OpenMM system: a donor is an N, O or S atom bonded to at least one hydrogen, and its hydrogens are those bonded to it; an acceptor is any O, or an N or S bonded to at most two atoms. Raises `ProtocolError` when a hydrogen of `atoms` has no bonded atom. Not a measurement: the atoms `hydrogen_bonds` and the functions below use unless given | `polyzymd analyze hydrogen_bonds`, recorded under `provenance.settings.hbond_atoms` |
| `hydrogen_bonds` | `group_a`, optional `group_b`, the production frames; `d_a_cutoff` (3.5 Å), `d_h_a_angle_cutoff` (150°), optional `donors`, `hydrogens` and `acceptors` | rows `HBOND_PARTS = ("mean_hbonds", "mean_residue_pairs", "any_fraction")` | MDAnalysis `HydrogenBondAnalysis` run once over the frames, with the hydrogens and acceptors of `hbond_atoms` unless given; a bond is a hydrogen of a donor within `d_a_cutoff` of an acceptor (minimum image of the frame's box) with a donor-hydrogen-acceptor angle of at least `d_h_a_angle_cutoff`. Each hydrogen is paired with its donor through the bonds, or by distance within 1.2 Å when `donors` is given. With `group_b`, only bonds with one partner in each group count, in either direction; without it, bonds within `group_a`. Bonds within one residue are left out. The rows are the mean hydrogen bonds per frame, the mean residue pairs joined by at least one per frame, and the fraction of frames with at least one | `polyzymd analyze hydrogen_bonds` |
| `hbond_lifetimes` | as `hydrogen_bonds`, and `key` (`residue`) and `tolerance_ps` (0) | rows `LIFETIME_PARTS = ("mean_lifetime", "n_events", "censored_fraction")` | The bonds of `hydrogen_bonds`; with `key="residue"` a pair is two residues, present on a frame when any hydrogen bond joins them, with `key="atom"` one donor atom and one acceptor atom; each pair's runs of frames are events, pooled as in `contact_lifetimes`: the Kaplan-Meier restricted mean lifetime in ns, the number of events and the fraction censored by the first or last frame | `polyzymd analyze hydrogen_bonds --run <s>_mean_lifetime` |
| `residue_hbond_occupancy` | as `hydrogen_bonds` | one value per residue of `group_a`, fraction of frames | The bonds of `hydrogen_bonds`; a residue counts on a frame when one of its atoms is the donor or acceptor of at least one | `polyzymd analyze hydrogen_bonds --run <s>_residues` |
| `residue_pair_hbond_occupancy` | as `hydrogen_bonds` | the pair labels, and one value per pair, fraction of frames | The bonds of `hydrogen_bonds`; a pair counts on a frame when at least one hydrogen bond joins it. A standard amino acid is named by its residue ID (`chain:resid` when IDs repeat) and any other residue by its residue name, so `149-SBM` is residue 149 with any SBM residue; with `group_b` the `group_a` residue comes first. Only pairs that form on some frame are returned, so run it with `labels="returned"` and `missing=0.0` | `polyzymd analyze hydrogen_bonds --run <s>_pairs` |
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
analysis `threshold`, 3.5 Å by default, when the pair sets none), computed from
the stored distance with `all_below`. Pick a result with `--run`. See
{doc}`../how_to/analysis_distances_quickstart`. The catalytic triad is a
routine on the study API rather than an analysis: it counts each triad
hydrogen bond with `hbond_count` and combines them with `Timeseries.transform`;
see {doc}`../how_to/analysis_triad_quickstart`.

`polyzymd analyze rmsf` and `polyzymd analyze rmsd_per_residue` run
`rms_decomposition` once per replicate, with `selection` measured and
`alignment_selection` fitted (both `protein and name CA` by default) against
`pz.reference(reference_mode, "(selection) or (alignment_selection)", frame=reference_frame, file=reference_file, alignment=alignment_selection)`.
Without `reference_mode`, the reference is `reference_file` when one is set
(`external`), otherwise the `centroid` frame. A `reference_file` with another
mode is refused. Each
replicate's `core_<part>` value is the square root of the mean over the core
residues of the part's mean-square row, so
`core_rmsd_per_residue² = core_rmsf² + core_offset²`. The core is the residues of
`selection` that `--set core=...` also selects, by default all of them. Each
entry of `--set regions='{name: selection}'` gives `<name>_<part>` the same way.
`mean_<part>` is the plain mean of the per-residue values, and `--run
rmsf`, `offset` or `rmsd_per_residue` compares the profile residue by residue. The
default result is `core_rmsf` for `rmsf` and `core_rmsd_per_residue` for
`rmsd_per_residue`. The residues of the core and of each region are recorded
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

`polyzymd analyze native_contacts` measures `native_contacts` of `--set
selection=...` (default `protein and not element H`) against
`pz.reference(reference_mode, selection, frame=reference_frame, file=reference_file)`,
on every production frame, and reduces each replicate to its mean. A missing
`reference_mode` is `external` when `reference_file` is set and `frame` with
production frame `reference_frame` (1) otherwise. `radius`, `min_separation`,
`beta`, `lambda_constant` and `use_pbc` are passed on. `--run q`, the default,
counts every native pair, and `--run <region>_q` for each entry of `--set
regions='{name: selection}'` only the pairs with an atom in the region. Only
the chosen result is measured. See
{doc}`../how_to/analysis_native_contacts_quickstart` and
{doc}`../explanation/analysis_native_contacts_verification`.

`polyzymd analyze secondary_structure` runs `dssp_occupancy` once per
replicate on `--set selection=...` (default `protein`, whole residues), in the
scheme of `--set scheme=...`: `simplified` (default) or `full`. `--run <class>`
reports a class's mean over the residues, the fraction of residue-frames in it,
and `--run <class>_residues` compares every residue. The default is the
scheme's first class, `helix` or `alpha_helix`. A warning names the replicates
with unassigned residues. See {doc}`../how_to/analysis_secondary_structure_quickstart`.

(dssp-classes)=
The two schemes use these classes. `simplified` calls
`mdtraj.compute_dssp(simplified=True)`, and `full` calls
`mdtraj.compute_dssp(simplified=False)`. `DSSP_GROUPS` gives the same mapping
in Python.

| `scheme=full` class | DSSP code | `scheme=simplified` class |
|---|---|---|
| `alpha_helix` | H | `helix` |
| `3_10_helix` | G | `helix` |
| `pi_helix` | I | `helix` |
| `extended_strand` | E | `strand` |
| `isolated_bridge` | B | `strand` |
| `turn` | T | `coil` |
| `bend` | S | `coil` |
| `loop` | blank | `coil` |
| `unassigned` | NA | `unassigned` |

`polyzymd analyze contacts` runs, once per replicate between `--set
protein_selection=...` (default null, the protein: `chainid A`) and
`polymer_selection` (default null, the polymer: `chainid C`; a replicate where it matches no atoms, such as a control without
polymer, has no contact: 0, and no comparison with its condition is tested),
`residue_occlusion` with `exposed_threshold`, `buried_threshold`, `max_asa`, `probe_radius_nm` and
`n_sphere_points` for `method=occlusion` (default), or `residue_contacts` with
`cutoff` for `method=distance`, on heavy atoms only when `heavy_atoms` is true
(default). `use_pbc` sets `pbc`, and every polymer residue name in any
condition gets a row (`polymer_types`), 0 in a replicate without it. `coverage`, the default, is the fraction
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

`polyzymd analyze hydrogen_bonds` counts hydrogen bonds for the summaries of
`--set summaries=...`, each `{between: [a, b]}` or `{within: a}` over the named
selections of `--set groups=...`; the defaults are `groups` `{protein: null,
polymer: null}`, where null selects the atoms of the role the group is named
after (`protein` chain A, `ligand` chain B, `polymer` chain C), and `summaries` `{protein_polymer: {between: [protein,
polymer]}}`. Only the chosen summary is measured. Each summary `<s>` gives
`<s>_mean_hbonds` (the default for the first summary), `<s>_mean_residue_pairs`
and `<s>_any_fraction` from `hydrogen_bonds`; `<s>_mean_lifetime`,
`<s>_lifetime_events` and `<s>_censored_fraction` from `hbond_lifetimes`, with
`lifetime_key` (`residue`) as `key` and `tolerance_ps` (0); `<s>_residues`
from `residue_hbond_occupancy`, one value per residue of the summary's first
group; and `<s>_pairs` from `residue_pair_hbond_occupancy`, with a pair that
forms in one replicate but not another counted 0 in the other. `d_a_cutoff`
(3.5 Å) and `d_h_a_angle_cutoff` (150°) set the geometry, and `donors`,
`hydrogens` and `acceptors` are selections, taken within the summary's
groups, that replace the atoms `hbond_atoms` chooses. The donors, hydrogens
and acceptors used are recorded under `provenance.settings.hbond_atoms`, as
counts per residue name and atom name. A replicate with no hydrogen bond in
the summary gets `nan` for a lifetime result and a warning. See
{doc}`../how_to/hydrogen_bonds` and
{doc}`../explanation/analysis_hydrogen_bonds_verification`.

## Figures

`polyzymd analyze` writes these figures to `<output-dir>/figures/<analysis>/`
and records the folder under `output_paths.figures` in the JSON report;
`--no-plots`, or `plots=False` in `analyze`, skips them.

| Analysis | Figures |
|---|---|
| `rg` | `rg_timeseries` (every replicate against time), `rg_comparison` (condition means with every replicate value), `rg_distribution` (per-frame distribution) |
| `rmsd` | `rmsd_timeseries`, `rmsd_comparison` |
| `rmsf`, `rmsd_per_residue` | `rmsf_profile`, `offset_profile` and `rmsd_per_residue_profile` (each replicate and each condition's mean with its interval, `highlight_residues` marked), `rms_decomposition` (the three profiles of each condition together), `rmsf_comparison` (the three core values), and with several conditions `rmsf_difference`, `offset_difference` and `rmsd_per_residue_difference` (each condition minus the control at every residue, with the interval of the difference and the significant residues marked) |
| `sasa` | For a total, `sasa_timeseries_<name>`, `sasa_comparison_<name>` and `sasa_distribution_<name>`; for `<name>_residues`, `sasa_profile_<name>` and with several conditions `sasa_difference_<name>` |
| `secondary_structure` | `ss_content_bars` (every class of the scheme but unassigned); for a total `ss_<name>_comparison`; for `<name>_residues`, `ss_<name>_profile`, `ss_classes_<name>` and with several conditions `ss_<name>_difference` |
| `contacts` | `contacts_class_bars` (each amino-acid class's mean contact fraction); for a one-value result `contacts_<run>_comparison`; for a residue result, `contacts_<name>_profile` and with several conditions `contacts_<name>_difference` |
| `native_contacts` | `native_contacts_timeseries_<run>` (every replicate against time) and `native_contacts_comparison_<run>` (condition means with every replicate value) |
| `distances` | `distance_kde_<pair>` for each pair with its threshold, `distance_fraction_<result>` for each fraction below threshold, and the grouped `distance_threshold_bars` (every pair's fraction) and `distance_kde_panel` (one panel per pair) |
| `hydrogen_bonds` | For a one-value result `hbonds_<run>_comparison` (condition means with every replicate value); for `<s>_residues` and `<s>_pairs`, `hbonds_<run>_profile` (each entry's value per replicate and each condition's mean with its interval) and with several conditions `hbonds_<run>_difference` (each condition minus the control at every residue or pair, with the interval of the difference and the significant entries marked) |
