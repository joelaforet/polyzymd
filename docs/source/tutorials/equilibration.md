# Multi-Stage Equilibration

This tutorial explains how to use PolyzyMD's multi-stage equilibration feature with temperature ramping and position restraints. This is essential for properly equilibrating complex biomolecular systems like enzyme-polymer simulations.

## Overview

Proper equilibration of biomolecular systems typically requires a staged approach:

1. **Heating**: Gradually raise temperature while restraining heavy atoms
2. **Relaxation**: Release restraints on flexible components (e.g., polymers) while keeping protein restrained
3. **Free equilibration**: Remove all restraints and equilibrate at target conditions

This protocol prevents the system from exploding due to initial bad contacts and allows different components to relax at appropriate rates.

## Why Use Multi-Stage Equilibration?

Standard single-stage equilibration can fail for complex systems because:

- **Initial bad contacts**: Packed systems often have steric clashes that need gentle resolution
- **Temperature shock**: Suddenly assigning velocities at 300K can destabilize the system
- **Component flexibility**: Polymers and proteins have different relaxation timescales
- **Solvent penetration**: Water needs time to properly solvate all surfaces

Multi-stage equilibration addresses these issues by:

1. Starting at low temperature (60K) where kinetic energy is low
2. Gradually heating while heavy atoms are restrained to reference positions
3. Releasing restraints in stages, allowing flexible components to relax first
4. Final unrestrained equilibration ensures the system is properly relaxed

## Configuration

### Basic Structure

Multi-stage equilibration is configured using the `equilibration_stages` key instead of the simple `equilibration` key:

```yaml
simulation_phases:
  # Use equilibration_stages for multi-stage protocols
  equilibration_stages:
    - name: "heating"
      duration: 0.288
      ensemble: "NVT"
      # ... stage configuration
    
    - name: "relaxation"
      duration: 1.0
      ensemble: "NVT"
      # ... stage configuration
    
    - name: "free_equilibration"
      duration: 0.5
      ensemble: "NPT"
      # ... stage configuration
  
  production:
    # ... production configuration
```

```{important}
Use **either** `equilibration` (simple mode) **or** `equilibration_stages` (multi-stage mode), not both. If `equilibration_stages` is present, it takes precedence.
```

### Complete Example

Here's a full multi-stage equilibration configuration based on common MD literature protocols:

```yaml
simulation_phases:
  equilibration_stages:
    # Stage 1: Heat from 60K to 293K with position restraints
    - name: "heating"
      duration: 0.288              # nanoseconds
      samples: 50                  # trajectory frames to save
      ensemble: "NVT"
      temperature_start: 60.0      # starting temperature (K)
      temperature_end: 293.0       # final temperature (K)
      temperature_increment: 1.0   # temperature step (K)
      temperature_interval: 1200.0 # time between steps (fs)
      position_restraints:
        - group: "protein_heavy"
          force_constant: 4184.0   # kJ/mol/nm^2
        - group: "polymer_heavy"
          force_constant: 4184.0

    # Stage 2: Relax polymers while keeping protein restrained
    - name: "polymer_relaxation"
      duration: 1.0                # nanoseconds
      samples: 100
      ensemble: "NVT"
      temperature: 293.0           # constant temperature (K)
      position_restraints:
        - group: "protein_heavy"
          force_constant: 4184.0

    # Stage 3: Free equilibration with NPT ensemble
    - name: "free_equilibration"
      duration: 0.5                # nanoseconds
      samples: 100
      ensemble: "NPT"
      temperature: 293.0
      # No position_restraints = fully unrestrained

  production:
    ensemble: "NPT"
    duration: 100.0
    samples: 2500
    time_step: 2.0
    thermostat: "LangevinMiddle"
    barostat: "MC"
    barostat_frequency: 25

  segments: 10  # Split production into 10 segments for HPC
```

---

## Stage Configuration Options

### Basic Options

| Option | Type | Required | Description |
|--------|------|----------|-------------|
| `name` | string | Yes | Identifier for the stage (used in output filenames) |
| `duration` | float | Yes | Duration in nanoseconds |
| `samples` | int | No | Number of trajectory frames to save (default: 100) |
| `ensemble` | string | Yes | Thermodynamic ensemble: `NVT` or `NPT` |
| `time_step` | float | No | Integration timestep in femtoseconds (default: 2.0) |

### Temperature Options

You can specify either constant temperature or temperature ramping:

**Constant Temperature:**
```yaml
temperature: 293.0  # Kelvin
```

**Temperature Ramping:**
```yaml
temperature_start: 60.0       # Starting temperature (K)
temperature_end: 293.0        # Final temperature (K)
temperature_increment: 1.0    # Step size (K) - default: 1.0
temperature_interval: 1200.0  # Time between steps (fs) - default: 1200.0
```

### Position Restraints

Position restraints hold atoms near their reference positions (post-minimization coordinates) using a harmonic potential:

```yaml
position_restraints:
  - group: "protein_heavy"
    force_constant: 4184.0  # kJ/mol/nm^2
  - group: "polymer_heavy"
    force_constant: 4184.0
```

To have no restraints (free dynamics), simply omit the `position_restraints` key or set it to an empty list:

```yaml
position_restraints: []
```

---

## Atom Groups

Position restraints are applied to predefined atom groups. PolyzyMD automatically determines which atoms belong to each group based on the system topology.

### Available Groups

| Group Name | Description | Typical Use |
|------------|-------------|-------------|
| `protein_heavy` | All non-hydrogen protein atoms | Restrain entire protein structure |
| `protein_backbone` | Protein backbone atoms (N, CA, C, O) | Restrain protein fold, allow sidechain motion |
| `protein_calpha` | Protein alpha carbons only | Minimal protein restraint |
| `ligand_heavy` | All non-hydrogen substrate/ligand atoms | Restrain docked ligand position |
| `polymer_heavy` | All non-hydrogen polymer atoms | Restrain polymer chains |
| `solvent` | All solvent molecules (water + ions + cosolvents) | Rarely used |
| `water_only` | Only water molecules | Rarely used |
| `ions_only` | Only ions (Na+, Cl-, etc.) | Rarely used |
| `cosolvents_only` | Only cosolvent molecules | Rarely used |

### How Atom Groups Are Determined

PolyzyMD uses the PDB chain IDs to identify system components:

- **Chain A**: Protein (enzyme)
- **Chain B**: Substrate/ligand (if present)
- **Chain C**: Polymers (if present)
- **No chain ID**: Solvent, ions, cosolvents

Within each chain, atoms are classified as "heavy" (non-hydrogen) based on their element.

```{note}
The chain assignment assumes PolyzyMD built the system. If you're using an externally prepared system, ensure the chain IDs follow this convention.
```

---

## Force Constants

Position restraints use a harmonic potential:

```
U(r) = 0.5 * k * (r - r0)^2
```

where `k` is the force constant and `r0` is the reference position.

### Unit Conversion

Force constants are specified in **kJ/mol/nm^2**. A common literature value is 1.0 kcal/mol/A^2:

```
1.0 kcal/mol/A^2 = 4.184 kJ/mol/A^2 = 4184 kJ/mol/nm^2
```

### Recommended Values

| Use Case | Force Constant (kJ/mol/nm^2) |
|----------|------------------------------|
| Strong restraint (heating phase) | 4184.0 (= 1.0 kcal/mol/A^2) |
| Medium restraint | 2092.0 (= 0.5 kcal/mol/A^2) |
| Weak restraint | 418.4 (= 0.1 kcal/mol/A^2) |

---

## Temperature Ramping Details

### How It Works

When temperature ramping is enabled (both `temperature_start` and `temperature_end` are specified), PolyzyMD:

1. Sets initial velocities at `temperature_start`
2. Runs dynamics for `temperature_interval` femtoseconds
3. Increases temperature by `temperature_increment` K
4. Repeats until `temperature_end` is reached
5. Runs any remaining time at `temperature_end`

### Calculating Duration

For smooth heating, ensure the duration matches your temperature increment and interval:

```
duration = (temp_end - temp_start) / temp_increment * temp_interval / 1e6
```

**Example**: Heating from 60K to 293K with 1K increments every 1200 fs:

```
(293 - 60) / 1.0 * 1200 / 1e6 = 0.2796 ns ~ 0.288 ns
```

```{tip}
It's fine if the duration is slightly longer than the calculated heating time. The simulation will run at the final temperature for the remaining time.
```

---

## Periodic Boundary Conditions

### How Position Restraints Handle PBC

Position restraints use OpenMM's `periodicdistance()` function to correctly handle atoms that may be outside the primary periodic box (common after minimization or NPT equilibration):

```
U = 0.5 * k * periodicdistance(x, y, z, x0, y0, z0)^2
```

This ensures restraints work correctly regardless of where atoms are positioned relative to the periodic box boundaries.

```{important}
This is a critical implementation detail. Earlier versions used absolute distance calculations which could cause simulation crashes (NaN energies) in large systems. If you're using a custom restraint implementation, ensure it handles PBC correctly.
```

---

## Output Files

Each equilibration stage creates its own output directory with:

| File | Description |
|------|-------------|
| `equilibration_N_stagename_trajectory.dcd` | Trajectory file |
| `equilibration_N_stagename_state_data.csv` | Energy, temperature, density data |
| `equilibration_N_stagename_topology.pdb` | Reference PDB for trajectory |
| `equilibration_N_stagename_checkpoint.chk` | Checkpoint for restart |

Where `N` is the stage index (0, 1, 2...) and `stagename` is from your configuration.

---

## Typical Protocols

### Enzyme-Only System

For a simple enzyme simulation without polymers:

```yaml
equilibration_stages:
  - name: "heating"
    duration: 0.288
    ensemble: "NVT"
    temperature_start: 60.0
    temperature_end: 300.0
    position_restraints:
      - group: "protein_heavy"
        force_constant: 4184.0

  - name: "equilibration"
    duration: 1.0
    ensemble: "NPT"
    temperature: 300.0
    # No restraints
```

### Enzyme-Polymer System

For systems with polymers (the full protocol):

```yaml
equilibration_stages:
  - name: "heating"
    duration: 0.288
    ensemble: "NVT"
    temperature_start: 60.0
    temperature_end: 293.0
    position_restraints:
      - group: "protein_heavy"
        force_constant: 4184.0
      - group: "polymer_heavy"
        force_constant: 4184.0

  - name: "polymer_relaxation"
    duration: 1.0
    ensemble: "NVT"
    temperature: 293.0
    position_restraints:
      - group: "protein_heavy"
        force_constant: 4184.0

  - name: "free_equilibration"
    duration: 0.5
    ensemble: "NPT"
    temperature: 293.0
```

### Gradual Restraint Release

For very sensitive systems, you can gradually reduce restraint strength:

```yaml
equilibration_stages:
  - name: "heating"
    duration: 0.288
    ensemble: "NVT"
    temperature_start: 60.0
    temperature_end: 300.0
    position_restraints:
      - group: "protein_heavy"
        force_constant: 4184.0  # Full strength

  - name: "relax_1"
    duration: 0.5
    ensemble: "NVT"
    temperature: 300.0
    position_restraints:
      - group: "protein_heavy"
        force_constant: 2092.0  # Half strength

  - name: "relax_2"
    duration: 0.5
    ensemble: "NVT"
    temperature: 300.0
    position_restraints:
      - group: "protein_backbone"
        force_constant: 1046.0  # Backbone only, quarter strength

  - name: "free_equilibration"
    duration: 0.5
    ensemble: "NPT"
    temperature: 300.0
```

---

## Troubleshooting

### Simulation Crashes with NaN

**Symptoms**: Energy becomes NaN during equilibration, simulation terminates.

**Common causes**:
1. **Bad initial contacts**: Increase minimization iterations or add a gentle heating stage
2. **Timestep too large**: Reduce `time_step` to 1.0 fs for the heating stage
3. **Missing restraints**: Ensure restraints are applied during heating

### Temperature Spikes

**Symptoms**: Temperature overshoots target during heating.

**Solutions**:
1. Use smaller `temperature_increment` (e.g., 0.5 K instead of 1.0 K)
2. Increase `temperature_interval` (e.g., 2400 fs instead of 1200 fs)
3. Add more samples to monitor temperature closely

### Protein Unfolds During Polymer Relaxation

**Symptoms**: Protein structure changes when polymer restraints are released.

**Solutions**:
1. Keep protein restraints on during polymer relaxation (this is the default protocol)
2. Use a stronger force constant on the protein
3. Extend the heating phase to better equilibrate the system

### Slow Equilibration

**Symptoms**: System hasn't reached equilibrium after all stages.

**Solutions**:
1. Extend the `free_equilibration` stage
2. Add additional stages with gradually decreasing restraints
3. Monitor state_data.csv files to check for energy/temperature convergence

---

## References

The multi-stage equilibration protocol implemented in PolyzyMD is based on common practices in the MD literature:

> "The simulation protocol first involved 2000 steps of energy minimization, after which the system was heated from 60 to 300 K over the course of 288 ps MD simulation by incrementing the temperature by 1 K every 600 time steps (1200 fs). The system was then equilibrated at 300 K for 1 ns before production MD simulation. The heating of the system and equilibration were performed using 1.0 kcal/mol/A^2 harmonic positional restraints on all protein backbone atoms."

This approach ensures stable equilibration of complex biomolecular systems while preventing artifacts from temperature shock or bad initial contacts.
