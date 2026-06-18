# Data Requirements & Directory Layout

This page documents the directory structures, file formats, and naming
conventions that PolyzyMD uses for simulations and analysis. Use it as a
lookup reference when setting up new projects or troubleshooting missing-file
errors.

---

## The Two-Project Pattern

PolyzyMD separates simulation execution from cross-condition analysis into two
distinct project types, each with its own directory scaffold and configuration
file:

| Project Type | Created By | Config File | Purpose |
|---|---|---|---|
| Simulation project | `polyzymd init -n <name>` | `config.yaml` | Build, run, and store one simulation condition |
| Comparison project | `polyzymd compare init -n <name>` | `comparison.yaml` | Analyze and compare results across conditions |

A comparison project does not contain trajectory data. Instead, its
`comparison.yaml` points to one or more simulation project `config.yaml` files,
which in turn resolve to the trajectory directories on disk.

---

## Simulation Project Layout

Running `polyzymd init -n my_simulation` creates:

```
my_simulation/
├── config.yaml              # Simulation configuration (edit this)
├── structures/              # Input PDB/SDF files
├── job_scripts/             # Generated SLURM submission scripts
└── slurm_logs/              # SLURM stdout/stderr logs
```

After building and running a simulation, the **output directory** (on
scratch or in the projects directory) grows to:

```
{scratch_dir}/{naming_template}/       # One directory per replicate
├── solvated_system.pdb                # Topology (created by polyzymd build)
├── equilibration_heating/             # Equilibration stage output
│   └── ...
├── production_0/                      # First production segment
│   ├── production_0_trajectory.dcd    # Trajectory
│   └── production_0_topology.pdb      # Topology snapshot
├── production_1/                      # Daisy-chain continuation segment
│   ├── production_1_trajectory.dcd
│   └── production_1_topology.pdb
└── ...                                # Additional segments if daisy-chained
```

Each replicate gets its own complete directory containing a topology file and
one or more trajectory segments.

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
| `{duration}` | `simulation_phases.production.duration` (integer ns) | `100` |
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

Current OpenMM runs write `solvated_system.pdb` in the replicate working
directory and production trajectories as indexed daisy-chain segments:
`production_N/production_N_trajectory.dcd`.

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

## Comparison Project Layout

Running `polyzymd compare init -n my_study` creates:

```
my_study/
├── comparison.yaml          # Analysis configuration (edit this)
├── comparison/              # Cross-condition comparison results
├── figures/                 # Generated plots
└── structures/              # (Optional) shared structure files (e.g., enzyme PDB for SASA)
```

Analysis runs also create and populate `analysis/` with canonical
`ReplicateArtifact` and `ConditionArtifact` outputs for per-replicate and
per-condition results. The `comparison/` directory is reserved for
cross-condition comparison results.

### comparison.yaml structure

The comparison config references simulation projects by pointing to their
`config.yaml` files:

```yaml
name: "polymer_study"
description: "Comparison of polymer conjugation effects"

control: "No Polymer"       # Label of the control condition, or null

conditions:
  - label: "No Polymer"
    config: "../no_polymer/config.yaml"
    replicates: [1, 2, 3]

  - label: "PEG 10k"
    config: "../peg_10k/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

plugins:
  rmsf:
    selection: "protein and name CA"
  # ... additional analysis plugins
```

Relative paths in `conditions[].config` are resolved relative to the directory
containing `comparison.yaml`, not the current working directory.

---

## Connecting It All Together

```
polyzymd init       -->  config.yaml  -->  polyzymd build  -->  polyzymd run  -->  trajectories/
                                                                                       |
polyzymd compare init  -->  comparison.yaml  -->  polyzymd compare run  -->  results + figures
                                  |
                    (points to config.yaml files)
```

The comparison framework reads each condition's `config.yaml`, resolves the
scratch directory and naming template, then uses `TrajectoryLoader` to find
topology and trajectory files for each replicate.

---

## Common Pitfalls

```{warning}
**Path resolution is config-relative, not CWD-relative.**
Relative paths in `config.yaml` (e.g., `enzyme.pdb_path: "structures/enzyme.pdb"`)
are resolved relative to the directory containing `config.yaml`, not your
shell's current working directory. The same applies to `conditions[].config`
paths in `comparison.yaml`.
```

- **Mismatched scratch directory.** If you built and ran simulations with one
  `scratch_directory` value but later changed it in `config.yaml`, the analysis
  framework will look in the wrong location. The `scratch_directory` in
  `config.yaml` must match where the trajectory files actually reside.

- **The `"default"` sentinel for reactions.** Setting
  `polymers.reactions.initiation: "default"` tells PolyzyMD to use a bundled
  reaction template. Writing `"structures/default"` or any path containing
  `"default"` will fail because no such file exists.

- **Missing replicate directories.** Each replicate number listed in
  `comparison.yaml` must have a corresponding directory on disk. If replicate 3
  was never simulated, the analysis will fail with a `FileNotFoundError`
  showing the expected path.

- **Incomplete replicate directories.** Every replicate directory must contain
  at least a topology file (`solvated_system.pdb`) and one or more production
  trajectory files. Partially completed simulations that crashed before writing
  a trajectory will cause load failures.

---

## See Also

- {doc}`configuration` -- Full configuration field reference
- {doc}`cli_reference` -- CLI command reference including `init` and `compare init`
- {doc}`../how_to/analysis_compare_conditions` -- How to set up and run a comparison
- {doc}`../get_started/quickstart` -- Run your first simulation end-to-end
