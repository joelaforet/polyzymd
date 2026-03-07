#!/usr/bin/env python
"""Create a synthetic test environment for smart-restart testing.

This script builds a tiny waterbox (~64 waters) using OpenMM directly,
runs a short production segment 0, and writes the exact directory structure
that ContinuationManager and the recover CLI command expect.

Output structure:
    test_restart_workdir/
        solvated_system.pdb
        production_0/
            production_0_state.xml
            production_0_system.xml
            production_0_parameters.json
            production_0_trajectory.dcd
            production_0_state_data.csv
            production_0_topology.pdb
            production_0_checkpoint.chk

Usage:
    # On the cluster, inside the polyzymd-env conda environment:
    mamba run -n polyzymd-env python setup_test_env.py

    # Or if already activated:
    python setup_test_env.py

Requirements:
    - OpenMM (available in polyzymd-env)
    - No PolyzyMD imports needed (pure OpenMM)

Takes ~10s on GPU, ~30s on CPU.
"""

import json
import time
from pathlib import Path

import openmm
from openmm import XmlSerializer, unit
from openmm.app import (
    DCDReporter,
    ForceField,
    Modeller,
    PDBFile,
    Simulation,
    StateDataReporter,
    Topology,
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
WORKDIR = Path("test_restart_workdir")
TEMPERATURE = 300.0  # Kelvin
PRESSURE = 1.0  # atm
TIMESTEP_FS = 2.0
FRICTION = 1.0  # 1/ps
NUM_PRODUCTION_STEPS = 1000  # Very short segment 0
NUM_SAMPLES = 10  # 10 frames for segment 0
BOX_SIZE_NM = 2.0  # Tiny box


def main():
    t0 = time.time()
    print("=" * 70)
    print("Smart-Restart Test Environment Setup")
    print("=" * 70)

    # Clean up previous test if exists
    if WORKDIR.exists():
        import shutil

        print(f"Removing existing {WORKDIR}/ ...")
        shutil.rmtree(WORKDIR)

    WORKDIR.mkdir(parents=True)
    seg_dir = WORKDIR / "production_0"
    seg_dir.mkdir()

    # -----------------------------------------------------------------
    # Step 1: Build a tiny waterbox
    # -----------------------------------------------------------------
    print("\n[1/5] Building tiny waterbox...")

    # Use TIP3P water model from amber14 force field
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

    # Create an empty topology and add solvent
    topology = Topology()
    topology.setUnitCellDimensions(openmm.Vec3(BOX_SIZE_NM, BOX_SIZE_NM, BOX_SIZE_NM))

    modeller = Modeller(topology, [])
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        boxSize=openmm.Vec3(BOX_SIZE_NM, BOX_SIZE_NM, BOX_SIZE_NM) * unit.nanometers,
    )

    n_atoms = modeller.topology.getNumAtoms()
    n_residues = sum(1 for _ in modeller.topology.residues())
    print(f"  Created waterbox: {n_residues} residues, {n_atoms} atoms")

    # -----------------------------------------------------------------
    # Step 2: Save solvated PDB (topology reference for ContinuationManager)
    # -----------------------------------------------------------------
    print("[2/5] Saving solvated_system.pdb...")
    pdb_path = WORKDIR / "solvated_system.pdb"
    with open(pdb_path, "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    print(f"  Saved to {pdb_path}")

    # -----------------------------------------------------------------
    # Step 3: Create system and run minimization + short equilibration
    # -----------------------------------------------------------------
    print("[3/5] Creating system and minimizing...")

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=openmm.app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=openmm.app.HBonds,
    )

    # Add barostat
    barostat = openmm.MonteCarloBarostat(
        PRESSURE * unit.atmospheres,
        TEMPERATURE * unit.kelvin,
        25,
    )
    system.addForce(barostat)

    # Create integrator
    integrator = openmm.LangevinMiddleIntegrator(
        TEMPERATURE * unit.kelvin,
        FRICTION / unit.picoseconds,
        TIMESTEP_FS * unit.femtoseconds,
    )

    # Try GPU first, fall back to CPU
    try:
        platform = openmm.Platform.getPlatformByName("CUDA")
        print("  Using CUDA platform")
    except Exception:
        try:
            platform = openmm.Platform.getPlatformByName("OpenCL")
            print("  Using OpenCL platform")
        except Exception:
            platform = openmm.Platform.getPlatformByName("CPU")
            print("  Using CPU platform (slower but works)")

    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    # Minimize
    simulation.minimizeEnergy()
    print("  Energy minimized")

    # Short equilibration (500 steps)
    simulation.context.setVelocitiesToTemperature(TEMPERATURE * unit.kelvin)
    simulation.step(500)
    print("  Equilibrated (500 steps)")

    # -----------------------------------------------------------------
    # Step 4: Run short production (segment 0)
    # -----------------------------------------------------------------
    print("[4/5] Running production segment 0...")

    report_interval = max(1, NUM_PRODUCTION_STEPS // NUM_SAMPLES)

    # Add reporters
    traj_path = seg_dir / "production_0_trajectory.dcd"
    simulation.reporters.append(DCDReporter(str(traj_path), report_interval))

    state_data_path = seg_dir / "production_0_state_data.csv"
    simulation.reporters.append(
        StateDataReporter(
            str(state_data_path),
            report_interval,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            speed=True,
        )
    )

    # Run production
    simulation.step(NUM_PRODUCTION_STEPS)
    print(f"  Completed {NUM_PRODUCTION_STEPS} steps ({NUM_SAMPLES} frames)")

    # -----------------------------------------------------------------
    # Step 5: Save end-of-segment files (what ContinuationManager needs)
    # -----------------------------------------------------------------
    print("[5/5] Saving end-of-segment files...")

    # Save topology PDB
    topo_pdb_path = seg_dir / "production_0_topology.pdb"
    state_for_pdb = simulation.context.getState(getPositions=True)
    with open(topo_pdb_path, "w") as f:
        PDBFile.writeFile(
            simulation.topology,
            state_for_pdb.getPositions(),
            f,
        )

    # Save state XML (portable)
    state = simulation.context.getState(
        getPositions=True,
        getVelocities=True,
        getForces=True,
        getEnergy=True,
        getParameters=True,
    )
    state_xml_path = seg_dir / "production_0_state.xml"
    with open(state_xml_path, "w") as f:
        f.write(XmlSerializer.serialize(state))
    print(f"  Saved {state_xml_path}")

    # Save system XML
    system_xml_path = seg_dir / "production_0_system.xml"
    with open(system_xml_path, "w") as f:
        f.write(XmlSerializer.serialize(system))
    print(f"  Saved {system_xml_path}")

    # Save checkpoint
    chk_path = seg_dir / "production_0_checkpoint.chk"
    simulation.saveCheckpoint(str(chk_path))
    print(f"  Saved {chk_path}")

    # Save parameters JSON (matches runner.py format exactly)
    duration_ns = (NUM_PRODUCTION_STEPS * TIMESTEP_FS) / 1e6  # fs -> ns
    params_dict = {
        "__class__": "SimulationParameters",
        "__values__": {
            "thermo_params": {
                "__class__": "ThermoParameters",
                "__values__": {
                    "temperature": {
                        "__class__": "Quantity",
                        "__values__": {"value": TEMPERATURE, "unit": "kelvin"},
                    },
                    "thermostat_params": {
                        "__class__": "ThermostatParameters",
                        "__values__": {
                            "temperature": {
                                "__class__": "Quantity",
                                "__values__": {
                                    "value": TEMPERATURE,
                                    "unit": "kelvin",
                                },
                            },
                            "timescale": {
                                "__class__": "Quantity",
                                "__values__": {
                                    "value": FRICTION,
                                    "unit": "/picosecond",
                                },
                            },
                        },
                    },
                    "barostat_params": {
                        "__class__": "BarostatParameters",
                        "__values__": {
                            "pressure": {
                                "__class__": "Quantity",
                                "__values__": {
                                    "value": PRESSURE,
                                    "unit": "atmosphere",
                                },
                            },
                            "temperature": {
                                "__class__": "Quantity",
                                "__values__": {
                                    "value": TEMPERATURE,
                                    "unit": "kelvin",
                                },
                            },
                            "frequency": 25,
                        },
                    },
                },
            },
            "integ_params": {
                "__class__": "IntegratorParameters",
                "__values__": {
                    "time_step": {
                        "__class__": "Quantity",
                        "__values__": {"value": TIMESTEP_FS, "unit": "femtosecond"},
                    },
                    "total_time": {
                        "__class__": "Quantity",
                        "__values__": {"value": duration_ns, "unit": "nanosecond"},
                    },
                    "num_samples": NUM_SAMPLES,
                },
            },
            "reporter_params": {
                "__class__": "ReporterParameters",
                "__values__": {
                    "report_interval": report_interval,
                    "report_trajectory": True,
                    "report_state_data": True,
                },
            },
        },
    }
    param_path = seg_dir / "production_0_parameters.json"
    with open(param_path, "w") as f:
        json.dump(params_dict, f, indent=2)
    print(f"  Saved {param_path}")

    # -----------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------
    elapsed = time.time() - t0
    print()
    print("=" * 70)
    print(f"Setup complete in {elapsed:.1f}s")
    print("=" * 70)
    print()
    print("Directory structure:")
    print(f"  {WORKDIR}/")
    print(f"    solvated_system.pdb          ({pdb_path.stat().st_size // 1024} KB)")
    print(f"    production_0/")
    for fp in sorted(seg_dir.iterdir()):
        print(f"      {fp.name:<45} ({fp.stat().st_size // 1024} KB)")
    print()
    print("Next steps:")
    print("  1. Run interactive tests:  bash test_restart_interactive.sh")
    print("  2. Or submit SLURM chain:  bash submit_test_chain.sh")
    print()


if __name__ == "__main__":
    main()
