# Troubleshooting

Find the error message or the symptom below, and do the steps that follow it.

:::{admonition} Environment Setup
:class: tip

PolyzyMD has one pixi environment for each kind of work. Use the environment
of the command that fails:

| Task | Pixi environment |
|------|------------------|
| PDB preparation, `polyzymd validate`, `polyzymd build`, `polyzymd run` | `build` |
| OpenMM `polyzymd submit` or `recover --submit` | `build`. The SLURM job activates the site environment |
| OpenMM inside a SLURM job on an NVIDIA GPU | `sim-cuda-12-0`, `sim-cuda-12-4` or `sim-cuda-12-6` |
| Trajectory analysis and plots, `polyzymd analyze ...` | `analysis` |

Use `pixi shell -e <env>` to activate an environment, or put
`pixi run -e <env>` before a command.
:::

## Installation

### `ModuleNotFoundError: No module named 'openmm'`

The shell is not in a PolyzyMD environment. Activate the environment of the
task, and check OpenMM:

```bash
pixi shell -e build
python -c "import openmm; print(openmm.__version__)"
```

### `ImportError: cannot import name ... from 'polyzymd'`

The environment is out of date. From the repository root, update the clone
and install the environment again:

```bash
git pull
pixi install -e build
```

## Configuration

### `Field required`

```
enzyme.pdb_path
  Field required [type=missing, ...]
```

A required key is missing. Add it. For the required keys of each section, see
{doc}`../reference/configuration`.

### An unknown key

Every section of the config refuses a key that it does not define. A
misspelled key is therefore an error, not a silent default. Correct the
spelling. The error names the key and its section.

### `FileNotFoundError` for an input file

PolyzyMD reads a relative path in `config.yaml` relative to the folder of
`config.yaml`, not to the folder of the shell. Check the path from the folder
of the config. You can also give an absolute path:

```yaml
enzyme:
  pdb_path: "/full/path/to/enzyme.pdb"
```

### `yaml.scanner.ScannerError: mapping values are not allowed here`

The YAML syntax is wrong. Check these points:

- Use the same indentation for keys at one level. Two spaces is common.
- Put a space after each colon.
- Put quotes around values with special characters, such as `:` or `#`.

### Several `Field required` errors for the items of a list

```
error: solvent.co_solvents.0: Value error, Co-solvent 'dmso': give exactly one of mole_fraction, concentration and count, not none
error: solvent.co_solvents.1.name: Field required
error: solvent.co_solvents.2.name: Field required
fix: correct these keys in config.yaml and run polyzymd validate again.
```

Each key was written as a separate list item. A `-` starts a new item. Put
all keys of one item at the same indentation, with no `-` before the second
and later keys.

Incorrect (three items):

```yaml
co_solvents:
  - name: "dmso"
  - mole_fraction: 0.1
  - residue_name: "DMS"
```

Correct (one item with three keys):

```yaml
co_solvents:
  - name: "dmso"
    mole_fraction: 0.1
    residue_name: "DMS"
```

The same rule applies to every list, such as `co_solvents`, `monomers` and
`restraints`.

## System build

### Charge assignment fails for a substrate, co-solvent or polymer

1. Check that the SDF file has correct bond orders, formal charges and
   explicit hydrogens.
2. Use `charge_method: nagl` ({term}`NAGL`). It needs no extra program.
3. For a molecule that NAGL cannot charge, compute the charges outside
   PolyzyMD and store them in the SDF file.

The `am1bcc` method ({term}`AM1-BCC`) needs AmberTools, which the PolyzyMD
environments do not include. Use it only in an environment where you
installed AmberTools.

### `Packmol exited with return code N`

PACKMOL did not find a placement. The error names a log file with the
PACKMOL output. Do one or more of these steps:

1. Make the box larger:

   ```yaml
   solvent:
     box:
       padding: 2.0    # default 1.2 nm
   ```

2. Use fewer or shorter polymer chains (`polymers.count`, `polymers.length`).
3. For polymer packing, raise `polymers.packing.nloop`, or set
   `polymers.packing.movebadrandom: true`. See {doc}`polymers`.

Exit code 173 ("imperfect packing") is not an error. The build uses the best
placement that PACKMOL found, and energy minimization removes the remaining
close contacts.

### `SolvationClashError: N solvent atom(s) lie within ... A of the solute`

```
SolvationClashError: 1182 solvent atom(s) lie within 1.00 A of the solute ...
solute/solvent frame mismatch: the packed coordinates are shifted against the solute
```

The packed molecules sit inside the protein. The solute and the packed
coordinates are in different coordinate frames. This is a software or input
problem, not a packing problem. An imperfect PACKMOL run leaves at most a few
such contacts, and the check allows up to 20. Do not raise the limit, and do
not continue the build.

1. Update the `build` environment and PolyzyMD.
2. Check the coordinates of the input PDB: no atoms at the origin, and no
   duplicated models.
3. Report the problem on GitHub. Include the error message, and
   `packmol_input.txt` and `_PACKING_SOLUTE.pdb` from the build folder.

### `PeriodicImageClashError: N atom(s) ... lie within ... A of a periodic image`

```
PeriodicImageClashError: 34 atom(s) of the packed solute + polymers lie within
1.00 A of a periodic image (... minimum image separation 0.104 A between atoms
(14363, 6312) across lattice vector (0, 0, 1) ...)
```

Some molecules lie outside the rectangular brick of the periodic cell, so
they overlap their own periodic images. PACKMOL cannot see this, because it
runs without periodic boundaries. Minimization cannot remove an overlap this
close, so the build stops.

1. A system built by PolyzyMD 1.2 or by an earlier 1.3 release candidate must
   be built again. Do not edit it. Those versions could make a box whose brick
   was too short for the solute.
2. PolyzyMD sizes the box so that the solute fits inside the brick. The build
   log line `Box: edge ...` gives the clearance to each brick face. If all
   three clearances are at least the tolerance and the error still appears,
   report it with that line.

A contact between half the tolerance and the full tolerance gives only a
warning. Minimization removes it.

### `ValueError: No atoms match selection: 'resid 76 and name OG'`

A restraint selection matches no atom.

1. Use the residue numbers of the built system. PolyzyMD renumbers the
   protein residues from 1. See {doc}`../explanation/residue_assignment`.
2. Check the atom names in `solvated_system.pdb`.
3. Test the selection in MDAnalysis or PyMOL on `solvated_system.pdb`.

## Simulation

### `OpenMMException: Particle coordinate is NaN`

Possible causes are clashes in the starting structure, a time step that is
too large, or a system that is not stable.

1. Open `solvated_system.pdb` in a viewer, and look for clashes and for
   chains that pass through a ring.
2. Make the time step of the first equilibration stage smaller:

   ```yaml
   simulation_phases:
     equilibration_stages:
       - name: "heating"
         time_step: 1.0    # default 2.0 fs
   ```

3. Add a stage that holds the protein while the solvent and polymer relax.
   See {doc}`equilibration`.

### `CUDA out of memory`

The GPU has no more memory. Make the system smaller:

```yaml
solvent:
  box:
    padding: 1.0    # a smaller box
polymers:
  count: 1          # fewer polymer chains
```

### `slurmstepd: error: Detected 1 oom_kill event`

The SLURM job used more memory than it requested. This happens most often
during the energy minimization of a large system. The default request is 3G.
Request more with `--memory`:

```bash
polyzymd submit -c config.yaml --memory 8G
polyzymd recover -c config.yaml -r 1 --submit --memory 8G
```

Start points for the request:

| System | Memory |
|---|---|
| A small protein with one polymer chain | 3G (default) |
| Two to five polymer chains | 4G |
| More than five chains, or a large protein | 6G to 8G |

### The simulation is slow

1. Check that OpenMM loads its GPU platforms:

   ```python
   import openmm
   print(openmm.Platform.getPluginLoadFailures())
   ```

2. Set the platform in `config.yaml`, for example `openmm.platform: CUDA`.
   PolyzyMD never falls back to the CPU. An unavailable platform stops the
   simulation.
3. Save fewer frames. Lower `simulation_phases.production.samples`.

## SLURM

### The job stays pending with `REASON=Resources`

No GPU is free. Wait, or use another partition:

```bash
polyzymd submit -c config.yaml --preset al40
```

### A cancelled job comes back

```
scancel 1234567
# a few seconds later
squeue -u $USER   # the same replicate is queued again
```

`scancel` sends `SIGTERM`. The job reads it as an interruption with work
left, and submits a successor. Stop the chain, not the job:

```bash
polyzymd cancel -c config.yaml -r 1-3            # write STOP and cancel the jobs
polyzymd cancel -c config.yaml -r 1-3 --resume   # let the chain run again
```

OpenMM and GROMACS job scripts check for `STOP`. Job scripts written by
PolyzyMD 1.2 and earlier do not. Stop those chains with
`scancel --batch --signal=KILL <job_id>`. See
{ref}`Stop a chain <hpc-slurm-stop-a-chain>`.

### `ModuleNotFoundError` in the SLURM output

The job did not activate the pixi environment. The scripts of
`polyzymd submit` activate it with `pixi shell-hook`. If you edit a script by
hand, keep this line:

```bash
eval "$(pixi shell-hook -e sim-cuda-12-4 --manifest-path /path/to/polyzymd/pixi.toml)"
```

### `FATAL: CUDA routing failed after 3 retries`

The chain got three nodes in a row whose driver is too old for the pinned
`sim-cuda-*` environment, or that could not reproduce the recorded runtime.
The message names each node that it tried.

1. Submit again and exclude all of those nodes:

   ```bash
   polyzymd submit -c config.yaml -r 1 --preset blanca-shirts \
       --exclude <preset nodes>,<node 1>,<node 2>,<node 3>
   ```

   `--exclude` replaces the list of the preset, so include the nodes of the
   preset too. For the Blanca list, see {ref}`excluded-blanca-gpu-nodes`.

2. If the same nodes come back often, ask the maintainers to add them to the
   list of the preset.

A segment that runs resets the count. So the message means three failures in
a row, not three failures in the whole simulation.

### `Segment N was hard-killed ... Moved to production_N.hardkilled-...`

A segment stopped without a grace period, for example after `SIGKILL`, an
out-of-memory kill or a node failure. It left no `INTERRUPTED` marker and no
`restart_state.xml`. PolyzyMD cannot tell how much of the segment is good. It
renames the folder and runs segment `N` again from the last good state.

You need to do nothing. The renamed folder keeps the frames that were
written. Delete it when you no longer need it. If this happens often, lower
the time limit so that segments end before the scheduler kills them.

### `PermissionError: /scratch/...`

The scratch folder does not exist, or you cannot write to it:

```bash
mkdir -p /scratch/$USER/simulations
chmod 755 /scratch/$USER/simulations
```

### `submit --dry-run` wrote no scripts

`--dry-run` only prints the submission plan. It writes no file. To write the
scripts without submitting them, use `--generate-only`. See
{doc}`../reference/cli_reference`.

## Continuation

### The next segment cannot find the previous segment

1. Check that the previous segment finished. Its folder,
   `production_<N-1>/`, holds `production_<N-1>_state.xml` when it finished.
2. Check that the job uses the same config and the same scratch directory as
   the earlier segments.

### An error names two different particle counts

The files of the replicate folder come from different builds. Restore all
files from the same build. Do not copy single files until the counts agree.
See {doc}`hpc_slurm`.

### `Could not reload ...: OpenMM checkpoints are not portable`

```
RuntimeError: Could not reload .../production_3_checkpoint.chk: OpenMM
checkpoints are not portable and this process differs from the one that
wrote it (pixi_environment 'sim-cuda-12-4' -> 'sim-cuda-12-0')
```

A binary `.chk` checkpoint reloads only under the OpenMM build that wrote it.
The chain moved to another environment between segments, and the previous
segment left no portable state file.

1. Submit with the recorded environment. `runtime_platform.json` in the
   replicate folder names the `pixi_environment` and `openmm_version` of the
   chain. Give that environment with `--pixi-env`.
2. If the environment must change, rename the segment folder that cannot be
   recovered. Do not delete it. PolyzyMD then runs the segment again from the
   last portable state.

PolyzyMD always loads a portable state file (`production_N_state.xml`,
`interrupted_state.xml` or `restart_state.xml`) before a checkpoint, and logs
its choice. So this error appears only when a checkpoint was the only file
left.

### `TrajectoryLineageError: Trajectory segments do not form a single contiguous time line`

The analysis loader refuses production segments whose times overlap, run
backward or leave a gap, or whose frame intervals differ. It repairs two
boundary faults of the OpenMM restart chain, and prints a warning:

- A segment starts at the time of the last frame of the previous segment.
  The restart simulated that step again, so the loader leaves the last frame
  of the previous segment out of `Replicate.frames`.
- A segment starts two frame intervals after the previous segment ends. One
  frame was never written. The frames after the gap keep their recorded times
  in `Replicate.times`.

The replicate provenance records both repairs under `segment_join`. The
loader still refuses a larger overlap or gap, or a change of frame interval.
Causes are a missing segment, two chains that wrote into one replicate folder
(for example after a duplicated SLURM submission), or a report interval that
changed during the chain.

1. Read the error. It names the two segments, with their first and last
   times.
2. Read `progress.json`. Each segment records `polyzymd_version`,
   `openmm_version`, `pixi_environment` and `started_at`. Compare it with the
   `production_N/` folders, and find the second chain.
3. Move the segments of the second chain out of the replicate folder. Do not
   delete them. Then run the analysis again.

To inspect such a replicate anyway, call
`TrajectoryLoader.load_universe(replicate, verify_lineage=False)`.

## Visualization

### Bonds stretch across the whole box in PyMOL or VMD

- If **some** molecules are broken and others are not, the atom order of the
  trajectory does not match the topology. Load the trajectory with the
  topology PDB from the same `production_N/` folder.
- If **all** molecules are broken, the molecules are wrapped into the box.
  Make them whole: `pbc unwrap` or `pbc join` in VMD, or
  `transformations.unwrap` in MDAnalysis.

Before a long production campaign, run a short test and look at its
trajectory.

## Get help

Collect this information:

```bash
polyzymd info                                       # PolyzyMD and package versions
pixi list -e build | grep -E "openmm|openff|pydantic"
polyzymd validate -c config.yaml
polyzymd -v build -c config.yaml --dry-run          # verbose output
polyzymd --openff-logs build -c config.yaml         # OpenFF force-field logs
```

PolyzyMD hides the OpenFF Interchange and Toolkit logs by default. These
libraries write messages for each atom, which can be millions of lines for a
large system. Use `--openff-logs` when you look for a force-field or charge
problem.

Open an issue at <https://github.com/joelaforet/polyzymd/issues>. Include:

1. the full error message and traceback;
2. the output of `polyzymd info`;
3. your configuration file, with private paths removed;
4. the steps that reproduce the error.
