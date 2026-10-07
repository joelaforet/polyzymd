# Data Requirements & Directory Layout

This page documents the directory structures, file formats, and naming
conventions that PolyzyMD uses for simulations and analysis. Use it as a
lookup reference when setting up new projects or troubleshooting missing-file
errors.

---

## Simulation Projects and Analysis Output

PolyzyMD keeps one simulation folder per condition, inside a {term}`study`.
Analysis reads those conditions and writes into a separate output directory:

| Directory | Created By | Holds |
|---|---|---|
| Simulation folder | `polyzymd study add-condition LABEL --new` | `config.yaml` and the inputs of one simulation condition; its replicates' trajectories live under the scratch or projects directory the config names, by default `runs/` of the project |
| Analysis output | `polyzymd analyze NAME -c A/config.yaml ...` or the study API | `polyzymd_results/` (every replicate's stored values and record), `figures/<analysis>/` and, with `--submit`, `slurm/`; the current directory unless `--output-dir` is given |

`polyzymd analyze` takes each condition's `config.yaml` with `-c`, control
first; no other file lists the conditions.

---

## Simulation Project Layout

`polyzymd project init my_paper --study lipa` followed by
`polyzymd study add-condition "No polymer" --new --study lipa` in `my_paper/`
creates:

```
my_paper/
├── project.yaml
├── lipa/
│   ├── study.yaml
│   └── conditions/no_polymer/
│       ├── config.yaml      # Simulation configuration (edit this)
│       └── structures/      # Input PDB/SDF files
└── runs/                    # Git-ignored; made by the first build or submit
    └── lipa/no_polymer/     # projects_directory of the config
        ├── job_scripts/     # Generated SLURM submission scripts
        └── ...              # One directory per replicate
```

`polyzymd submit` writes the SLURM logs into `slurm_logs/` of the folder you
run it from, and creates that folder when it submits.

The runs go into `runs/` unless you set `scratch_directory` in the config.
Trajectories can use a lot of disk space. On a cluster, set
`scratch_directory` to scratch storage.

After building and running a simulation, the **output directory** (on
scratch or in the projects directory) grows to:

```
{scratch_dir}/{naming_template}/       # One directory per replicate
├── solvated_system.pdb                # Viewer topology: names and coordinates (polyzymd build)
├── system.prmtop                      # Analysis topology: every atom and bond (polyzymd build)
├── system.xml                         # OpenMM System with restraints (polyzymd build)
├── build_manifest.json                # Hashes, config hash, versions, PACKMOL seeds
├── progress.json                      # Stage/segment records with runtime provenance
├── minimization/
│   ├── minimized_state.xml
│   └── phase.json                     # frozen_atoms, frozen_rmsd_angstrom, hydrogen_max_displacement_angstrom
├── equilibration_0_heating/           # Equilibration stage output
│   └── ...
├── production_0/                      # First production segment
│   ├── production_0_trajectory.dcd    # Trajectory
│   ├── production_0_topology.pdb      # Topology snapshot
│   └── production_0_parameters.json   # Parameters + "provenance" block
├── production_1/                      # Daisy-chain continuation segment
│   ├── production_1_trajectory.dcd
│   └── production_1_topology.pdb
├── production_2.hardkilled-20260909T181530Z/   # Retired unrecoverable segment
├── runtime_platform.json              # Pinned pixi env, OpenMM build, driver
├── STOP                               # Present only while the chain is stopped
└── ...                                # Additional segments if daisy-chained
```

Each replicate gets its own complete directory containing a topology file and
one or more trajectory segments.

### Control and marker files

| File | Written by | Meaning |
|------|------------|---------|
| `STOP` | `polyzymd cancel` | While present, the job wrapper submits no successor and a queued successor exits before starting a segment. Plain text: who stopped the chain, on which host, when, with which config and replicate, and how to undo it. Remove it (or run `polyzymd cancel --resume`) to allow resubmission |
| `runtime_platform.json` | first job of the chain | The runtime the chain is pinned to: `pixi_environment`, `openmm_version`, `platform`, `precision`, plus the observed driver and compute capability. A node that cannot reproduce the first four is rerouted, not run |
| `.successor-<job_id>` | job wrapper | Receipt proving a successor was already queued for that job, so a signal trap and normal exit cannot queue two |
| `.polyzymd.lock` | `run-segment` | Per-replicate `flock`; a second job on the same replicate exits with code 2 instead of running concurrently |
| `production_N.hardkilled-<ISO timestamp>` | `run-segment` | A segment that was hard-killed (no `INTERRUPTED` marker, no `restart_state.xml`, only a stale checkpoint) and could not be resumed. It is renamed out of the way, never deleted, so its frames survive; segment `N` is then re-run from the previous good state. Delete these once you no longer need the data |

A segment that *does* leave `restart_state.xml` behind is recoverable and is
never retired: the next job resumes from that portable state.

### Provenance files

| File | Key contents |
|------|--------------|
| `build_manifest.json` | SHA-256 of `solvated_system.pdb`, `system.xml` and `system.prmtop` when written, the config hash, `openmm_version`, `polyzymd_version`, and `provenance` (see below). `polyzymd submit --skip-build` refuses bundles whose config hash no longer matches `config.yaml` |
| `progress.json` | One record per equilibration stage and production segment with `polyzymd_version`, `openmm_version`, `pixi_environment` (null in files written by older versions). Rebuilt from a filesystem scan after a chain is cancelled and resubmitted; the segment provenance is then recovered from each `production_N_parameters.json` rather than lost. Published atomically (unique temporary file, fsync, rename) so a second writer cannot corrupt it |
| `production_N/production_N_parameters.json` | Simulation parameters plus a top-level `provenance` block (`polyzymd_version`, `openmm_version`, `pixi_environment`, `hostname`, `slurm_job_id`) |
| `minimization/phase.json` | Phase status, state path, `frozen_atoms` (the number of solute **heavy** atoms held fixed), `frozen_rmsd_angstrom` (0.0 when the solute was frozen), and `hydrogen_max_displacement_angstrom` (how far the furthest solute hydrogen moved onto its force-field constraint length; `null` for unfrozen minimization and for records written by older versions) |

#### `build_manifest.json` provenance keys

| Key | Meaning |
|-----|---------|
| `polymer_seed`, `packmol_seed` | Replicate seed used for polymer generation and PACKMOL |
| `polymer_packmol_seed` | Seed passed to the polymer-packing PACKMOL run |
| `solvent_packmol_seed` | Seed passed to the solvation PACKMOL run |
| `box_vectors_nm` | The periodic cell as a 3x3 row-major matrix in nm |
| `brick_nm` | Diagonal of the cell — the rectangular brick that PACKMOL fills, in nm |
| `deterministic_box` | `true` when the cell was computed from the protein + substrate before packing (the default for polymer builds); `false` for solute-only builds |
| `polymer_sphere_radius_nm` | Radius of the confinement sphere used for the chains (absent when `confine_to_sphere: false`) |

Replicates of one condition must agree on `box_vectors_nm`, `brick_nm` and
`polymer_sphere_radius_nm`, and therefore on their water and ion counts; they
differ only in the seeds. Diffing two manifests is the quickest way to confirm
that.

When loading multi-segment trajectories, the analysis loader verifies that the
segments form one contiguous, evenly spaced time line and raises
`TrajectoryLineageError` otherwise. A repeated last frame or one missing frame
at a segment boundary is repaired with a warning instead (see
{doc}`../how_to/troubleshooting`).

---

## Directory Naming Template

The `naming_template` field in the `output` section of `config.yaml` controls
how per-replicate directories are named.

**Default template:**

```
{enzyme}_{substrate}_{polymer_type}_{duration}ns_{temperature}K_run{replicate}
```

**Available placeholders:**

| Placeholder | Source | Example Value |
|---|---|---|
| `{enzyme}` | `enzyme.name` | `LipA` |
| `{substrate}` | `substrate.name` (hyphens removed), or `apo` if null | `ResorufinButyrate` |
| `{polymer_type}` | Derived from polymer config, or `none` if disabled | `SBMA-EGPMA_A70_B30` |
| `{temperature}` | `thermodynamics.temperature` (integer) | `300` |
| `{replicate}` | Replicate number (1-indexed) | `1` |
| `{duration}` | `simulation_phases.production.duration` in ns: whole ns from 1 ns up, in full below 1 ns | `100`, `0.005` |
| `{primary_solvent}` | Primary solvent token | `water_tip3p` |
| `{cosolvent_composition}` | Co-solvents sorted by normalized name, or `none` | `dmso_30molpct_urea_2p5M` |
| `{solvent_composition}` | Primary solvent plus co-solvents when present | `water_tip3p_dmso_30molpct` |

PolyzyMD uses the same resolved template for daisy-chain SLURM job names, so
job names match the per-replicate run directories after SLURM-safe
sanitization.

Solvent tokens are normalized for directory and SLURM job-name safety. Spaces,
slashes, parentheses, percent signs, and raw decimal points are replaced or
removed; concentration decimals use `p` (for example, `2.5 M` becomes
`urea_2p5M`). Mole-fraction co-solvents are rendered as mol-percent tokens
with the `molpct` suffix (for example, `mole_fraction: 0.30` becomes
`dmso_30molpct`). Ions are not included in solvent naming placeholders.

**Example resolved name:**

```
LipA_ResorufinButyrate_SBMA-EGPMA_A70_B30_100ns_300K_run1
```

---

## Scratch vs Projects Directories

PolyzyMD supports separating lightweight project files (scripts, logs) from
large simulation output (trajectories, checkpoints). This is common on HPC
systems where long-term storage and high-performance scratch are different
filesystems.

| Field | Purpose | Example |
|---|---|---|
| `projects_directory` | Scripts, configs, SLURM logs | `/projects/user/polyzymd` |
| `scratch_directory` | Trajectories, checkpoints, state data | `/scratch/alpine/user/simulations` |

If `scratch_directory` is `null` or omitted, all output goes to
`projects_directory`.

**Example `config.yaml` snippet:**

```yaml
output:
  projects_directory: "/projects/$USER/polyzymd"
  scratch_directory: "/scratch/alpine/$USER/simulations"
  naming_template: "{enzyme}_{substrate}_{polymer_type}_{duration}ns_{temperature}K_run{replicate}"
```

Environment variables (`$USER`, `$HOME`, `${VAR}`) and `~` are expanded
automatically in both path fields.

---

## What the Analysis Framework Expects

The `TrajectoryLoader` class resolves trajectory paths from a simulation
`config.yaml`. It uses the config's `scratch_directory` (or
`projects_directory` as fallback) combined with the `naming_template` to
locate each replicate's working directory.

### Topology and trajectory layout

Current OpenMM runs write two topologies in the replicate working directory
and production trajectories as indexed daisy-chain segments:
`production_N/production_N_trajectory.dcd`.

- `system.prmtop` is the analysis topology. It is built from the OpenMM
  topology and System together, so it carries every atom, residue, element,
  mass, charge and bond, including constrained bonds, with no column widths
  and no atom limit. Analyses load it when it is present. It has no chain
  IDs, so the loader takes them from `solvated_system.pdb` beside it, and
  `chainid A` selects the protein. For this it reads only the residue name and
  chain ID columns of the PDB's `ATOM` and `HETATM` lines, so it works at any
  system size.
- `solvated_system.pdb` is the viewer topology, for PyMOL or VMD with the DCD
  segments. Above 99,999 atoms OpenMM writes its serials in hex and MDAnalysis
  cannot read its CONECT records, and OpenMM writes CONECT records only for
  non-standard residues in any case, so as a topology it is a fallback for
  analysis, not the intended input.

Runs built before `system.prmtop` existed get one from their PDB and
`system.xml` with `polyzymd analysis-topology RUN_DIR...`. OpenMM's own PDB
reader accepts the hex serials it writes, so this works for any size.

GROMACS runs have the same split. `prod.tpr`, the compiled run input that
grompp writes before production, plays the role of `system.prmtop`: every
atom, bond, mass and charge, no atom limit, read natively by MDAnalysis.
Analyses prefer it over `solvated_system.pdb` and over any `.gro`, which
carries no bonds at all. No extra step is needed; the run directory keeps
`prod.tpr` after production.

MDAnalysis reads a TPR only up to the file version it knows: MDAnalysis 2.10
reads GROMACS 2025 files but not GROMACS 2026 ones. For a TPR it cannot read,
the loader warns and reads the run's `<prefix>.top`, with the `.itp` files it
includes, through MDAnalysis's `ITPParser`, laid out as MDAnalysis lays out a
TPR. Atoms, residue numbers, segments, charges, masses, elements (from the
atomic numbers in `[ atomtypes ]`), bonds, angles and dihedrals are then the
same as from the TPR, so results do not depend on the GROMACS version: on a
94,019-atom lysozyme-polymer run, `hydrogen_bonds`, `rmsd` and `contacts`
stored identical values from a GROMACS 2026 run read this way and from the
same run's TPR compiled by GROMACS 2025.

A TPR or `.top` names chains after molecule types (`MOL0`, `MOL1`, ...).
GROMACS universes take PolyzyMD's chain IDs (A protein, B substrate, C
polymer) from the build's `solvated_system.pdb` in the replicate directory,
when its atom count and residue names match, so `chainid A` and the other
default selections mean what they mean for an OpenMM run. Otherwise the
loader warns and the chains keep the molecule-type names.

| GROMACS file | Read for |
|---|---|
| `prod.tpr` | Atoms, bonds, charges, masses, elements |
| `<prefix>.top` and its `.itp` files | The same, when MDAnalysis cannot read `prod.tpr` |
| `solvated_system.pdb` (replicate directory) | Chain IDs |
| `prod_centered.xtc`, `prod_nojump.xtc` or `prod.xtc` | Coordinates |

Each completed production segment of an OpenMM run records, in
`progress.json`, the SHA-256 and size of its trajectory (`trajectory_sha256`,
`trajectory_bytes`) and the PolyzyMD and OpenMM versions that ran it.
Stored analysis results identify their input files by SHA-256 and size, taken
from there or computed once and cached. For GROMACS runs, downsampled
copies, and runs that finished before segments recorded their hashes,
`polyzymd hash-trajectories -c config.yaml` records them in
`trajectory_hashes.json` in the engine working directory, and never writes
`progress.json`; running it again changes nothing.

The provenance of each replicate names where the bonds came from, as
`bond_source`: `system_xml`, `tpr`, `top`, `conect`, `guessed` or `none`.

When multiple daisy-chain segments exist (e.g., `production_0/`,
`production_1/`, `production_2/`), they are automatically stitched together in
segment-index order using the MDAnalysis `ChainReader`. The resulting
`Universe` presents all segments as a single continuous trajectory.

---

## Input File Requirements

These are the input files placed in the simulation project's `structures/`
directory and referenced from `config.yaml`.

| File | Format | Config Field | Requirements |
|---|---|---|---|
| Protein structure | PDB (`.pdb`) | `enzyme.pdb_path` | Standard residue names, protonated at simulation pH, no missing heavy atoms in regions of interest |
| Substrate | SDF (`.sdf`) | `substrate.sdf_path` | 3D coordinates with docked pose, explicit hydrogens preferred |
| Polymer (if pre-built) | SDF (`.sdf`) | `polymers.sdf_directory` | One SDF per chain, or use dynamic generation from SMILES |
| Reaction templates | RXN or `"default"` | `polymers.reactions.*` | The string `"default"` loads bundled ATRP templates; a file path loads a custom template |

```{note}
The sentinel value `"default"` for reaction templates is **not** a file path.
Do not prepend a directory to it. PolyzyMD resolves `"default"` to bundled
reaction files at runtime.
```

---

## Analysis Output Layout

`polyzymd analyze hydrogen_bonds -c noPoly/config.yaml -c sbma/config.yaml
--output-dir results` writes:

```
results/
├── polyzymd_results/
│   └── <name>/                  # one folder per measurement name, e.g. hydrogen_bonds_protein_polymer
│       └── <condition>/
│           └── replicate_<n>/
│               ├── series.npz   # per-frame values (Study.timeseries), or
│               ├── values.npz   # one value or labelled array (Study.per_replicate)
│               └── record.json  # function, arguments, config hash, input files, frames, versions
├── figures/
│   └── hydrogen_bonds/          # the analysis's figures, unless --no-plots
└── slurm/                       # with --submit: one folder per submission
    └── hydrogen_bonds_<YYYYmmdd-HHMMSS>/
        ├── tasks.tsv
        ├── replicates.sbatch
        ├── report.sbatch
        ├── logs/
        └── report.txt           # report.json with --format json, or the -o path
```

A stored replicate result is read back instead of measured when every field
of its `record.json` matches the new call; see
{doc}`study_api`.

---

## Connecting It All Together

```
study add-condition -->  config.yaml  -->  polyzymd build  -->  polyzymd run  -->  trajectories/
                                                                                       |
polyzymd analyze NAME -c A/config.yaml -c B/config.yaml  -->  report + polyzymd_results/ + figures/
```

`polyzymd analyze` reads each condition's `config.yaml`, resolves the scratch
directory and naming template, then uses `TrajectoryLoader` to find the
topology and trajectory files of each replicate.

---

## Common Pitfalls

```{warning}
**Path resolution is config-relative, not CWD-relative.**
Relative paths in `config.yaml` (e.g., `enzyme.pdb_path: "structures/enzyme.pdb"`)
are resolved relative to the directory containing `config.yaml`, not your
shell's current working directory. Relative `-c` paths of `polyzymd analyze`
are resolved from the current directory.
```

- **Mismatched scratch directory.** If you built and ran simulations with one
  `scratch_directory` value but later changed it in `config.yaml`, analysis
  looks in the wrong location. The `scratch_directory` in
  `config.yaml` must match where the trajectory files actually reside.

- **The `"default"` sentinel for reactions.** Setting
  `polymers.reactions.initiation: "default"` tells PolyzyMD to use a bundled
  reaction template. Writing `"structures/default"` or any path containing
  `"default"` will fail because no such file exists.

- **Missing replicate directories.** Each replicate number given with
  `--replicates` must have a run directory on disk. If replicate 3 was never
  simulated, `polyzymd analyze` exits 2 and names the replicates it found.

- **Incomplete replicate directories.** Every replicate directory must contain
  at least a topology file (`system.prmtop` or `solvated_system.pdb`) and one or more production
  trajectory files. Partially completed simulations that crashed before writing
  a trajectory will cause load failures.

---

## See Also

- {doc}`configuration` -- Full configuration field reference
- {doc}`cli_reference` -- CLI command reference including `init` and `analyze`
- {doc}`../how_to/analysis_compare_conditions` -- How to set up and run a comparison
- {doc}`../get_started/quickstart` -- Run your first simulation end-to-end
