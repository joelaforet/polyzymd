# Troubleshooting

Common issues and solutions for PolyzyMD.

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

**Solution:** Install via conda (not pip):

```bash
conda install -c conda-forge openmm
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

**Solution:** Check your YAML has all required fields. See {doc}`configuration` for required fields.

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
  Value error, Co-solvent 'dmso': Must specify either 'volume_fraction' or 'concentration'
solvent.co_solvents.1.name
  Field required [type=missing, input_value={'volume_fraction': 0.5}, input_type=dict]
solvent.co_solvents.2.name
  Field required [type=missing, input_value={'residue_name': 'DMS'}, input_type=dict]
```

**Cause:** Each field was written as a separate list item instead of grouping all fields under one item.

**Incorrect** (creates 3 separate items):
```yaml
co_solvents:
  - name: "dmso"
  - volume_fraction: 0.5
  - residue_name: "DMS"
```

**Correct** (one item with 3 fields):
```yaml
co_solvents:
  - name: "dmso"
    volume_fraction: 0.5
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
     equilibration:
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

**Solution:** Reduce segment duration:

```yaml
simulation_phases:
  production:
    duration: 100.0
  segments: 20    # More segments = shorter each
```

### "Module not found in job"

```
ModuleNotFoundError in SLURM output
```

**Solution:** Ensure your job script loads conda:

1. Edit generated script in `job_scripts/`
2. Add at the top:
   ```bash
   module load anaconda
   source activate polymerist-env
   ```

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

### "State mismatch"

```
ValueError: System state doesn't match checkpoint
```

**Solution:** Don't modify the system between segments. If you need to change parameters, start fresh.

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

## Getting Help

### Collect Debug Information

```bash
# Package version
polyzymd info

# Python environment
conda list | grep -E "openmm|openff|pydantic"

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
