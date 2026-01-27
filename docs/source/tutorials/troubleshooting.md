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

```
CUDA out of memory
```

**Solutions:**

1. Reduce system size:
   ```yaml
   solvent:
     box:
       padding: 1.0    # Smaller box
   polymers:
     count: 1          # Fewer polymers
   ```

2. Use single precision (default in OpenMM)

3. Request more GPU memory in SLURM:
   ```bash
   #SBATCH --mem=128G
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

## Getting Help

### Collect Debug Information

```bash
# Package version
polyzymd info

# Python environment
conda list | grep -E "openmm|openff|pydantic"

# Configuration validation
polyzymd validate -c config.yaml

# Verbose mode
polyzymd -v build -c config.yaml --dry-run
```

### Report Issues

When reporting issues, include:

1. Full error message and traceback
2. Output of `polyzymd info`
3. Your configuration file (sanitized)
4. Steps to reproduce

Open issues at: https://github.com/joelaforet/polyzymd/issues
