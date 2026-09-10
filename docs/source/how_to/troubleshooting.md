# Troubleshooting

Common issues and solutions for PolyzyMD.

:::{admonition} Environment Setup
:class: tip

PolyzyMD v1.3 uses split pixi environments. Pick the environment that matches
the command you are troubleshooting:

| Task | Pixi environment |
|------|------------------|
| Project setup, PDB preparation, `polyzymd validate`, `polyzymd build` | `build` |
| OpenMM `polyzymd submit` or `recover --submit` | `build`; the Slurm job activates the site runtime |
| Direct OpenMM `polyzymd run-segment` | `sim-cuda-12-4` or `sim-cuda-12-6` |
| Trajectory comparison, plotting, `polyzymd compare ...` | `analysis` |

Use `pixi shell -e <env>` to activate an environment, or prefix a command with
`pixi run -e <env>`.
:::

## Installation Issues

### "README.md not found"

```
OSError: Readme file does not exist: README.md
```

**Solution:** Pull latest changes or create the file:

```bash
git pull origin main
# or
echo "# PolyzyMD" > README.md
pip install -e .
```

### "Module not found: openmm"

```
ModuleNotFoundError: No module named 'openmm'
```

**Solution:** OpenMM is provided by pixi environments. Make sure you are inside
the environment for the task you are running:

```bash
pixi shell -e build            # setup, validation, system building
pixi shell -e sim-cuda-12-6    # CUDA simulation on a matching GPU cluster
python -c "import openmm; print(openmm.__version__)"
```

### "Import error after installation"

```
ImportError: cannot import name 'xxx' from 'polyzymd'
```

**Solution:** Reinstall in development mode:

```bash
pip uninstall polyzymd
pip install -e .
```

---

## Configuration Errors

### "Validation failed: field required"

```
ValidationError: enzyme -> pdb_path: field required
```

**Solution:** Check your YAML has all required fields. See {doc}`../reference/configuration` for required fields.

### "File not found"

```
FileNotFoundError: structures/enzyme.pdb
```

**Solution:** 
- Check the path is correct relative to where you run the command
- Use absolute paths if needed:
  ```yaml
  enzyme:
    pdb_path: "/full/path/to/enzyme.pdb"
  ```

### "Invalid YAML syntax"

```
yaml.scanner.ScannerError: mapping values are not allowed here
```

**Solution:** Check YAML formatting:
- Consistent indentation (2 spaces recommended)
- Colons followed by space
- Quotes around special characters

### "Multiple validation errors" / "Field required" for list items

```
Build failed: 3 validation errors for SimulationConfig
solvent.co_solvents.0
  Value error, Co-solvent 'dmso': Must specify either 'mole_fraction' or 'concentration'
solvent.co_solvents.1.name
  Field required [type=missing, input_value={'mole_fraction': 0.1}, input_type=dict]
solvent.co_solvents.2.name
  Field required [type=missing, input_value={'residue_name': 'DMS'}, input_type=dict]
```

**Cause:** Each field was written as a separate list item instead of grouping all fields under one item.

**Incorrect** (creates 3 separate items):
```yaml
co_solvents:
  - name: "dmso"
  - mole_fraction: 0.1
  - residue_name: "DMS"
```

**Correct** (one item with 3 fields):
```yaml
co_solvents:
  - name: "dmso"
    mole_fraction: 0.1
    residue_name: "DMS"
```

**Solution:** The `-` character starts a **new list item**. All fields belonging to the same item must be indented to the same level *without* a leading `-`. This applies to all list-based configurations:
- `co_solvents`
- `monomers`
- `restraints`

---

## System Building Errors

### "Charge assignment failed"

```
ChargeAssignmentError: Unable to assign charges
```

**Solutions:**

1. Try a different charge method:
   ```yaml
   substrate:
     charge_method: "am1bcc"    # Instead of "nagl"
   ```

2. Check your SDF file has correct bond orders and hydrogens

3. For complex molecules, pre-compute charges externally

### "PACKMOL failed to converge"

```
PackmolError: PACKMOL did not converge
```

**Solutions:**

1. Increase box size:
   ```yaml
   solvent:
     box:
       padding: 2.0    # Increase from 1.2
   ```

2. Increase tolerance:
   ```yaml
   solvent:
     box:
       tolerance: 3.0    # Increase from 2.0
   ```

3. Reduce number of polymers:
   ```yaml
   polymers:
     count: 1    # Reduce from higher number
   ```

4. For polymer packing specifically, raise `polymers.packing.nloop` or enable
   `polymers.packing.movebadrandom`. Exit code 173 ("imperfect packing") is
   not a failure: the build continues and minimization resolves the residual
   contacts.

### "SolvationClashError: N solvent atom(s) lie within 1.00 A of the solute"

```
SolvationClashError: 1182 solvent atom(s) lie within 1.00 A of the solute ...
solute/solvent frame mismatch, see d96b1fcd
```

The assembled solute and the packed solvent or polymer coordinates are not in
the same coordinate frame, so packed molecules sit inside the protein. This is
a software or input problem, not a packing-quality problem: PACKMOL never
places atoms this close to a fixed solute, and imperfect packing leaves at most
a handful of such contacts (the assertion allows up to 20). Do not lower the
limit or continue with the build.

**Solutions:**

1. Make sure you are running a PolyzyMD release that includes the
   BRICK-centered solute fix (1.2.1 or later) and that the `build` environment
   is up to date.
2. Check that the input PDB coordinates are sensible (no atoms at the origin,
   no duplicated models).
3. Report the `packmol_input.txt` and `_PACKING_SOLUTE.pdb` from the build
   directory together with the error message.

### "PeriodicImageClashError: N atom(s) ... lie within 1.00 A of a periodic image"

```
PeriodicImageClashError: 34 atom(s) of the packed solute + polymers lie within
1.00 A of a periodic image (... minimum image separation 0.104 A between atoms
(14363, 6312) across lattice vector (0, 0, 1) ...)
```

Molecules were packed outside the rectangular brick of the periodic cell, so
they overlap themselves across the cell boundary. PACKMOL cannot see this: it
runs without periodicity. Minimization cannot resolve a singular overlap
either, and the run would die with NaN, so the build stops.

**Solutions:**

1. Rebuild with PolyzyMD 1.3.0-rc5 or later, which computes the periodic cell
   before packing and packs the chains inside its brick. Bundles built by
   earlier versions must be rebuilt, not patched.
2. If you set `polymers.packing.box_vectors`, remove it: an explicit packing box
   opts out of the deterministic cell and can be larger than the brick.
3. If the message appears with `deterministic_box: true` in
   `build_manifest.json`, check the logged clearance to the brick faces. A
   solute that is long along `z` may not fit inside a rhombic-dodecahedron
   brick at the configured padding (the `z` clearance is
   `0.707 * padding - 0.146 * bbox_z`); raise `solvent.box.padding` or switch
   to `solvent.box.shape: cube`.
4. Contacts between half the tolerance and the full tolerance only produce a
   warning; those are resolved by minimization and need no action.

### "No atoms match selection"

```
ValueError: No atoms match selection: 'resid 77 and name OG'
```

**Solutions:**

1. Check residue numbering in your PDB
2. Verify atom names (case-sensitive in some contexts)
3. Open PDB in PyMOL to verify:
   ```
   PyMOL> select test, resid 77 and name OG
   ```

---

## Simulation Errors

### "NaN encountered"

```
OpenMMException: Particle coordinate is NaN
```

**Causes:**
- Bad initial structure (clashes)
- Time step too large
- Unstable system

**Solutions:**

1. Run energy minimization (should be automatic)

2. Reduce time step:
   ```yaml
   simulation_phases:
      equilibration_stages:
        - name: "heating"
          time_step: 1.0    # Reduce from 2.0 fs
   ```

3. Check initial structure for clashes in VMD/PyMOL

### "Out of memory"

There are two types of out-of-memory errors you may encounter:

#### GPU Memory (CUDA OOM)

```
CUDA out of memory
```

This means the GPU ran out of VRAM. **Solutions:**

1. Reduce system size:
   ```yaml
   solvent:
     box:
       padding: 1.0    # Smaller box
   polymers:
     count: 1          # Fewer polymers
   ```

2. Use single precision (default in OpenMM)

#### System Memory (SLURM OOM)

```
slurmstepd: error: Detected 1 oom_kill event in StepId=...
```

This means the job exceeded its allocated RAM. This often happens during energy minimization when loading large systems. **Solutions:**

1. Increase memory allocation using the `--memory` flag:
   ```bash
   # Default is 3G, increase for larger systems
   polyzymd submit -c config.yaml --memory 4G
   
   # For very large systems (many polymers, large proteins)
   polyzymd submit -c config.yaml --memory 8G

   # Also works with recover --submit
   polyzymd recover -c config.yaml -r 1 --submit --memory 8G
   ```

2. If using generated scripts directly, edit the `#SBATCH --mem` line:
   ```bash
   #SBATCH --mem=4G    # Increase from 3G
   ```

```{tip}
**Memory guidelines:**
- Small systems (1 polymer, small protein): 3G (default)
- Medium systems (2-5 polymers): 4G
- Large systems (5+ polymers, large proteins): 6-8G
```

### "Simulation too slow"

**Solutions:**

1. Verify GPU is being used:
   ```python
   import openmm
   print(openmm.Platform.getPluginLoadFailures())
   ```

2. Use CUDA platform explicitly (should be automatic)

3. Reduce output frequency:
   ```yaml
   simulation_phases:
     production:
       samples: 1000    # Fewer frames
   ```

---

## SLURM/HPC Errors

### "Job pending: Resources"

```
squeue shows REASON=Resources
```

**Solution:** Wait for GPUs to become available, or use different partition:

```bash
polyzymd submit -c config.yaml --preset al40    # Try different GPU type
```

### "Job failed: time limit"

```
TIMEOUT in job output
```

**Solution:** Use a shorter production duration for test runs or move to a
longer-walltime preset. PolyzyMD segments production automatically.

### "Module not found in job"

```
ModuleNotFoundError in SLURM output
```

**Solution:** Ensure your job script activates the pixi environment. PolyzyMD-generated
scripts handle this automatically via `pixi shell-hook`. If you are editing scripts
manually:

1. Edit generated script in `job_scripts/`
2. Add the pixi activation:
   ```bash
   eval "$(pixi shell-hook -e sim-cuda-12-4 --manifest-path /path/to/polyzymd/pixi.toml)"
   ```

### "A cancelled job comes straight back"

```
scancel 1234567
# ...seconds later
squeue -u $USER   # the same replicate is queued again
```

`scancel` sends `SIGTERM`. `run-segment` treats that as a graceful
interruption and exits 99, which the job wrapper reads as "interrupted, work
remains" — so it queues a successor, and each attempt creates another
`production_N` directory.

**Solution:** stop the chain rather than the job:

```bash
polyzymd cancel -c config.yaml -r 1-3          # marker + scancel
polyzymd cancel -c config.yaml -r 1-3 --resume # allow it to run again
```

Chains submitted before this feature existed keep their old script and do
not check the marker; stop those with
`scancel --batch --signal=KILL <job_id>`.

### "FATAL: CUDA routing failed after 3 retries"

The chain drew three consecutive nodes whose driver is too old for the pinned
`sim-cuda-*` environment, or that could not reproduce the recorded runtime.
The message names every node it tried.

**Solutions:**

1. Resubmit excluding all of them at once — the message gives you the list:
   ```bash
   polyzymd submit -c config.yaml -r 1 --preset blanca-shirts \
       --exclude bgpu-g4-u20,bgpu-g4-u24,bgpu-g4-u30
   ```
2. If the same nodes keep appearing, add them to the preset's exclusion list
   so no chain draws them.

A node that runs a segment successfully resets the budget, so this message
means three failures in a row, not three failures over the life of the run.

### "Segment N was hard-killed ... Moved to production_N.hardkilled-..."

A segment was killed without a grace period (SIGKILL, OOM, node failure), so
it left no `INTERRUPTED` marker and no `restart_state.xml` — only a stale
checkpoint at an unknown position. PolyzyMD cannot tell how much of that
segment is trustworthy, so it renames the directory out of the way and re-runs
segment `N` from the last good state.

**What to do:** nothing is required; the run continues. The retired directory
still holds the frames that were written, so inspect or salvage it if you
need to, and delete it when you no longer do. If you see many of these,
lower the wall-time or raise the preemption grace period so segments end
gracefully.

### "Permission denied on scratch"

```
PermissionError: /scratch/...
```

**Solution:** Check scratch directory exists and is writable:

```bash
mkdir -p /scratch/$USER/simulations
chmod 755 /scratch/$USER/simulations
```

---

## Continuation Errors

### "Checkpoint not found"

```
FileNotFoundError: checkpoint.chk not found
```

**Causes:**
- Previous segment didn't complete
- Wrong working directory

**Solutions:**

1. Check previous segment completed:
   ```bash
   ls -la /path/to/simulation/production_seg0/
   ```

2. Verify working directory path is correct

### "Could not reload ... OpenMM checkpoints are not portable"

```
RuntimeError: Could not reload .../production_3_checkpoint.chk: OpenMM
checkpoints are not portable and this process differs from the one that
wrote it (pixi_environment 'sim-cuda-12-4' -> 'sim-cuda-12-0')
```

Binary `.chk` checkpoints only reload under the OpenMM build that wrote them.
The chain was rerouted to a different environment between segments, and the
hard-killed previous segment left no portable state XML to recover from.

**Solutions:**

1. Prefer the recorded environment: the chain's `runtime_platform.json` names
   the `pixi_environment` and `openmm_version` it is pinned to. Submit with
   that `--pixi-env` (or `auto`, which reuses the recorded value).
2. If the environment genuinely has to change, retire the unrecoverable
   segment directory (rename it, do not delete it) so the segment is re-run
   from the last portable state.

PolyzyMD always prefers a portable serialized state
(`production_N_state.xml`, `interrupted_state.xml`, `restart_state.xml`) over
a checkpoint, and logs which one it chose, so this error only appears when a
checkpoint was the only surviving option.

### "State mismatch"

```
ValueError: System state doesn't match checkpoint
```

**Solution:** Don't modify the system between segments. If you need to change parameters, start fresh.

### "TrajectoryLineageError: Trajectory segments do not form a single contiguous time line"

The analysis loader refuses to concatenate production segments whose raw
timestamps overlap, run backwards, or leave gaps, or whose frame intervals
differ. This usually means two restart chains wrote into the same run
directory (for example a duplicated SLURM resubmission) or a segment is
missing.

**Solutions:**

1. Read the error: it names the offending segment pair with their first and
   last times.
2. Inspect `progress.json` (each segment records `polyzymd_version`,
   `openmm_version`, `pixi_environment`, `started_at`) and the
   `production_N/` directories to identify the branched chain.
3. Move the branched segments out of the run directory (do not delete them),
   then rerun the analysis. To inspect a bad chain anyway, call
   `TrajectoryLoader.load_universe(replicate, verify_lineage=False)`.

---

## Visualization Issues

### "Molecules appear broken/scrambled in PyMOL or VMD"

**Symptoms:**
- Bonds appear to span the entire simulation box (hundreds of angstroms)
- Molecules look like a "ball of bonds" or spaghetti
- Some molecules look correct while others are completely scrambled
- Water and ions particularly affected, but protein looks fine

**Quick Diagnosis:**
- If **ALL molecules** are broken → likely a PBC (periodic boundary conditions) wrapping issue
- If only **SOME molecules** are broken → likely an **atom order mismatch** between your trajectory and topology files

**Solutions:**

1. **For atom order mismatch:** This is a subtle but devastating bug where the atom order in your DCD trajectory doesn't match the atom order in your topology file. We encountered this exact issue during PolyzyMD development and wrote a detailed guide:
   
   See: {doc}`broken_molecules_debugging` - A complete debugging case study with diagnosis steps, root cause analysis, and solutions.

2. **For PBC wrapping:** Use your visualization software's unwrap/make-whole tools:
   - **PyMOL:** Use `intra_fit` command or external post-processing tools
   - **VMD:** Use `pbc unwrap` or `pbc join` commands
   - **MDAnalysis:** Use `transformations.unwrap()` 

**Prevention:** Always verify that new simulation pipelines produce correct trajectories by running a short test and checking visualization BEFORE committing to long production runs.

---

## Known Limitations

### Analysis supports OpenMM trajectories only

The `polyzymd compare` analysis workflow currently expects OpenMM-produced
trajectories (DCD format) in PolyzyMD's standard directory layout
(`production_N/production_N_trajectory.dcd`). GROMACS XTC trajectories
are not yet supported.

**Workarounds:**

- Use native GROMACS analysis tools (`gmx rms`, `gmx rmsf`, etc.)
- Use MDAnalysis directly with your GROMACS topology and XTC files

GROMACS trajectory support in `polyzymd compare` is planned
([#47](https://github.com/joelaforet/polyzymd/issues/47)).

---

## Submission Issues

### "I ran submit --dry-run but no scripts were created"

**Cause**: In v1.3.0, `submit --dry-run` is preview-only and writes nothing to disk.

**Fix**: Use `submit --generate-only` to create scripts without submitting them.
See the {doc}`/reference/cli_reference` for the full options reference.

---

## Getting Help

### Collect Debug Information

```bash
# Package version
polyzymd info

# Python environment
pixi list -e build | grep -E "openmm|openff|pydantic"

# Configuration validation
polyzymd validate -c config.yaml

# Debug mode for troubleshooting
polyzymd --debug build -c config.yaml --dry-run

# Enable OpenFF logs for debugging force field issues
polyzymd --openff-logs build -c config.yaml
```

```{note}
By default, verbose OpenFF Interchange/Toolkit logs are suppressed to keep log files readable. 
These libraries generate per-atom INFO messages that can produce millions of lines for large systems.
Use `--openff-logs` when debugging force field parameter or charge assignment issues.
```

### Report Issues

When reporting issues, include:

1. Full error message and traceback
2. Output of `polyzymd info`
3. Your configuration file (sanitized)
4. Steps to reproduce

Open issues at: https://github.com/joelaforet/polyzymd/issues
