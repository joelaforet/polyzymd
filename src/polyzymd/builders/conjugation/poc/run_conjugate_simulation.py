#!/usr/bin/env python
"""
Full simulation workflow for polymer-protein conjugate.

Steps:
  1. Run conjugation_poc.py to build the conjugate and create the Interchange
  2. Solvate with TIP3P water + 137 mM NaCl
  3. Create solvated Interchange (ff14SB + Sage + TIP3P)
  4. Energy minimization (full system)
  5. Multi-stage equilibration:
     - Stage 1: NVT heating 60K -> 310K (0.3 ns) with protein+polymer heavy atoms restrained
     - Stage 2: NVT polymer relaxation at 310K (0.5 ns) with protein heavy atoms restrained
     - Stage 3: NPT free equilibration at 310K (0.5 ns) no restraints
  6. Production NPT (1 ns) at 310K, 1 atm
  7. Validate: print energy, temperature, density from state data

Usage:
    pixi run -e conjugation python run_conjugate_simulation.py
"""

from __future__ import annotations

import importlib.util
import logging
import sys
import time
from pathlib import Path

import numpy as np

# ── Configure logging ──────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("conjugate_sim")

# ── Paths ──────────────────────────────────────────────────────────────────
_POC_DIR = Path(__file__).resolve().parent
NOTEBOOK_PATH = _POC_DIR / "conjugation_poc.py"
OUTPUT_DIR = _POC_DIR / "output" / "simulation_output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Simulation parameters (matching polyzymd sample_config.yaml) ───────────
TARGET_TEMP = 310.0  # K
PRESSURE = 1.0  # atm
TIMESTEP = 2.0  # fs
FRICTION = 1.0  # 1/ps (= 1 ps timescale)
BAROSTAT_FREQ = 25  # steps
RESTRAINT_K = 4184.0  # kJ/mol/nm^2 (= 1 kcal/mol/A^2)

# Box / solvation
PADDING_NM = 1.5  # nm (smaller than 2.5 to keep system manageable for CPU)
NACL_CONC = 0.137  # mol/L
TARGET_DENSITY = 1.0  # g/mL

# Equilibration stages
HEATING_DURATION_NS = 0.030  # 30 ps (reduced from 300 ps for POC)
HEATING_TEMP_START = 60.0  # K
HEATING_TEMP_END = 310.0  # K
HEATING_TEMP_INCREMENT = 10.0  # K per update (faster ramp for POC)
HEATING_TEMP_INTERVAL_FS = 1200.0  # fs between temperature updates

POLYMER_RELAX_DURATION_NS = 0.050  # 50 ps (reduced from 500 ps for POC)

FREE_EQUIL_DURATION_NS = 0.050  # 50 ps (reduced from 500 ps for POC)

# Production
PRODUCTION_DURATION_NS = 0.100  # 100 ps for POC validation (500 ns for real)
PRODUCTION_SAMPLES = 200  # trajectory frames

# Minimization
MINIMIZATION_MAX_ITER = 10000
MINIMIZATION_TOLERANCE = 50.0  # kJ/mol/nm


# ═══════════════════════════════════════════════════════════════════════════
# Step 1: Run conjugation_poc.py to get the Interchange
# ═══════════════════════════════════════════════════════════════════════════
def run_conjugation_notebook():
    """Execute conjugation_poc.py and return key objects."""
    logger.info("=" * 70)
    logger.info("STEP 1: Running conjugation_poc.py to build conjugate")
    logger.info("=" * 70)

    t0 = time.perf_counter()

    # Suppress noisy logs from the notebook
    for name in [
        "openff.toolkit",
        "openff.interchange",
        "openff.units",
        "conjugation_poc",
    ]:
        logging.getLogger(name).setLevel(logging.WARNING)

    if not NOTEBOOK_PATH.exists():
        raise FileNotFoundError(f"Conjugation notebook not found: {NOTEBOOK_PATH}")

    spec = importlib.util.spec_from_file_location("conjugation_poc", str(NOTEBOOK_PATH))
    if spec is None:
        raise ImportError(f"Could not create module spec for notebook: {NOTEBOOK_PATH}")

    if spec.loader is None:
        raise ImportError(f"Module spec has no loader for notebook: {NOTEBOOK_PATH}")

    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    outputs, defs = mod.app.run()

    # Restore log levels
    for name in [
        "openff.toolkit",
        "openff.interchange",
        "openff.units",
        "conjugation_poc",
    ]:
        logging.getLogger(name).setLevel(logging.INFO)

    interchange = defs["interchange"]
    conjugate_off = defs["conjugate_off"]
    energy_before = defs["energy_before"]
    energy_after = defs["energy_after"]
    component_metadata = defs["component_metadata"]
    free_polymer_offs = defs["free_polymer_offs"]
    conjugated_sequences = defs["conjugated_sequences"]
    free_sequences = defs["free_sequences"]

    assert component_metadata is not None, "component_metadata missing from notebook"
    # Mirror notebook configuration: N_FREE_POLYMERS = 5
    assert len(free_polymer_offs) == 5, f"Expected 5 free polymers, got {len(free_polymer_offs)}"
    # Mirror notebook configuration: N_CONJUGATED = 2
    assert len(component_metadata["conjugated_polymers"]) == 2
    # Mirror notebook configuration: N_CONJUGATED = 2
    assert len(conjugated_sequences) == 2
    # Mirror notebook configuration: CONJUGATED_LENGTH = 10
    assert all(len(s) == 10 for s in conjugated_sequences)
    # Mirror notebook configuration: N_FREE_POLYMERS = 5
    assert len(free_sequences) == 5
    # Mirror notebook configuration: FREE_POLYMER_LENGTH = 5
    assert all(len(s) == 5 for s in free_sequences)

    dt = time.perf_counter() - t0
    logger.info(
        "Conjugation complete in %.1fs: %d atoms, E_min=%.1f kJ/mol",
        dt,
        conjugate_off.n_atoms,
        energy_after,
    )

    return interchange, conjugate_off, free_polymer_offs, component_metadata, defs


# ═══════════════════════════════════════════════════════════════════════════
# Step 2: Solvate the conjugate
# ═══════════════════════════════════════════════════════════════════════════
def solvate_system(interchange, conjugate_off, free_polymer_offs):
    """Solvate the conjugate with TIP3P water + NaCl ions.

    Returns a solvated OpenFF Topology with positions set as conformers.
    """
    logger.info("=" * 70)
    logger.info("STEP 2: Solvating conjugate")
    logger.info("=" * 70)

    from openff.interchange.components._packmol import UNIT_CUBE, solvate_topology
    from openff.toolkit import Topology
    from openff.units import Quantity, unit

    t0 = time.perf_counter()

    # solvate_topology needs a plain OpenFF Topology with positions stored
    # as molecule conformers (not the InterchangeTopology wrapper).
    # Build a fresh Topology from the conjugate molecule which already has
    # correct conformers.
    topology = Topology.from_molecules([conjugate_off, *free_polymer_offs])

    # Verify positions are accessible
    positions = topology.get_positions()
    assert positions is not None, "Topology must have positions for solvation"

    logger.info("Solute: %d atoms, %d molecules", topology.n_atoms, topology.n_molecules)

    # Solvate using OpenFF's built-in solvate_topology
    # This adds TIP3P water + NaCl ions in a cubic box
    logger.info(
        "Solvating with TIP3P water, %.3f M NaCl, padding=%.1f nm, cube box",
        NACL_CONC,
        PADDING_NM,
    )

    solvated_topology = solvate_topology(
        topology=topology,
        nacl_conc=Quantity(NACL_CONC, unit.mole / unit.liter),
        padding=Quantity(PADDING_NM, unit.nanometer),
        box_shape=UNIT_CUBE,
        target_density=Quantity(TARGET_DENSITY, unit.gram / unit.milliliter),
        tolerance=Quantity(0.2, unit.nanometer),  # 2 Angstrom default
    )

    # Count components
    n_total = solvated_topology.n_atoms
    n_solute = topology.n_atoms
    n_solvent = n_total - n_solute

    # Count water and ions
    n_water = 0
    n_na = 0
    n_cl = 0
    for mol in solvated_topology.molecules:
        smiles = mol.to_smiles()
        if "O" in smiles and mol.n_atoms == 3:
            n_water += 1
        elif "[Na+]" in smiles:
            n_na += 1
        elif "[Cl-]" in smiles:
            n_cl += 1

    logger.info(
        "Solvation complete: %d total atoms (%d solute + %d solvent)",
        n_total,
        n_solute,
        n_solvent,
    )
    logger.info("  Water: %d, Na+: %d, Cl-: %d", n_water, n_na, n_cl)

    box_vectors = solvated_topology.box_vectors
    if box_vectors is not None:
        bv = box_vectors.m_as(unit.nanometer)
        logger.info(
            "  Box: %.2f x %.2f x %.2f nm",
            bv[0][0],
            bv[1][1],
            bv[2][2],
        )

    dt = time.perf_counter() - t0
    logger.info("Solvation took %.1fs", dt)

    return solvated_topology


# ═══════════════════════════════════════════════════════════════════════════
# Step 3: Create solvated Interchange
# ═══════════════════════════════════════════════════════════════════════════
def create_solvated_interchange(solvated_topology, conjugate_off, free_polymer_offs):
    """Create Interchange for the solvated system.

    Uses ff14SB for protein, Sage for polymer, and pre-computed charges
    from the conjugation notebook.
    """
    logger.info("=" * 70)
    logger.info("STEP 3: Creating solvated Interchange")
    logger.info("=" * 70)

    from openff.toolkit import ForceField

    t0 = time.perf_counter()

    ff = ForceField(
        "ff14sb_off_impropers_0.0.4.offxml",
        "openff-2.0.0.offxml",
    )

    logger.info("Creating Interchange with ff14SB + Sage...")
    logger.info(
        "  Topology: %d atoms, %d molecules",
        solvated_topology.n_atoms,
        solvated_topology.n_molecules,
    )

    # The conjugate_off molecule has pre-computed charges (ff14SB for protein,
    # NAGL for polymer). Pass it via charge_from_molecules so the FF doesn't
    # try to re-compute charges.
    # Suppress per-atom logging flood from OpenFF nonbonded assignment
    _nonbonded_logger = logging.getLogger("openff.interchange.smirnoff._nonbonded")
    _prev_level = _nonbonded_logger.level
    _nonbonded_logger.setLevel(logging.WARNING)
    try:
        charge_molecules = [conjugate_off, *free_polymer_offs]
        solvated_interchange = ff.create_interchange(
            solvated_topology,
            charge_from_molecules=charge_molecules,
        )
    finally:
        _nonbonded_logger.setLevel(_prev_level)

    dt = time.perf_counter() - t0
    logger.info("Interchange created in %.1fs", dt)

    return solvated_interchange


# ═══════════════════════════════════════════════════════════════════════════
# Step 4-7: OpenMM simulation
# ═══════════════════════════════════════════════════════════════════════════
def run_openmm_simulation(
    solvated_interchange, conjugate_off, free_polymer_offs, component_metadata
):
    """Run the full OpenMM simulation: minimize -> equilibrate -> produce."""
    import openmm
    from openff.interchange.interop.openmm._positions import to_openmm_positions
    from openmm import XmlSerializer
    from openmm import app as omm_app
    from openmm import unit as omm_unit

    # ── Convert Interchange to OpenMM ──────────────────────────────────────
    logger.info("=" * 70)
    logger.info("STEP 4: Converting to OpenMM and minimizing")
    logger.info("=" * 70)

    t0 = time.perf_counter()

    logger.info("Converting Interchange to OpenMM topology...")
    omm_topology = solvated_interchange.to_openmm_topology()
    logger.info("  Topology: %d atoms", sum(1 for _ in omm_topology.atoms()))

    logger.info("Converting Interchange to OpenMM system...")
    omm_system = solvated_interchange.to_openmm(combine_nonbonded_forces=True)
    logger.info("  System forces: %d", omm_system.getNumForces())

    logger.info("Converting Interchange to OpenMM positions...")
    omm_positions = to_openmm_positions(solvated_interchange, include_virtual_sites=True)

    dt = time.perf_counter() - t0
    logger.info("OpenMM conversion took %.1fs", dt)

    # Save system XML for reproducibility
    system_xml_path = OUTPUT_DIR / "system.xml"
    with open(system_xml_path, "w") as f:
        f.write(XmlSerializer.serialize(omm_system))
    logger.info("Saved system XML: %s", system_xml_path)

    # ── Atom groups from metadata ──────────────────────────────────────────
    groups = component_metadata["restraint_groups"]
    protein_heavy_indices = groups["protein_heavy"]
    backbone_indices = groups["protein_backbone_heavy"]
    conjugate_heavy_indices = groups["conjugate_heavy"]
    free_polymer_heavy_indices = groups["free_polymer_heavy"]

    assert not (set(free_polymer_heavy_indices) & set(conjugate_heavy_indices)), (
        "Free polymer indices overlap with conjugate indices!"
    )

    logger.info(
        "Atom groups (from metadata): %d protein heavy, %d backbone heavy, "
        "%d conjugate heavy, %d free polymer heavy, %d total atoms",
        len(protein_heavy_indices),
        len(backbone_indices),
        len(conjugate_heavy_indices),
        len(free_polymer_heavy_indices),
        sum(1 for _ in omm_topology.atoms()),
    )

    # ── Helper: add/remove positional restraints ───────────────────────────
    def add_position_restraints(system, atom_indices, positions, k_kj_mol_nm2):
        """Add harmonic positional restraints to specified atoms."""
        force = openmm.CustomExternalForce("0.5*k*periodicdistance(x, y, z, x0, y0, z0)^2")
        force.addPerParticleParameter("k")
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")

        for idx in atom_indices:
            pos = positions[idx]
            x0 = pos[0].value_in_unit(omm_unit.nanometer)
            y0 = pos[1].value_in_unit(omm_unit.nanometer)
            z0 = pos[2].value_in_unit(omm_unit.nanometer)
            force.addParticle(idx, [k_kj_mol_nm2, x0, y0, z0])

        force_idx = system.addForce(force)
        logger.info(
            "  Added restraints to %d atoms (k=%.0f kJ/mol/nm^2, force_idx=%d)",
            len(atom_indices),
            k_kj_mol_nm2,
            force_idx,
        )
        return force_idx

    def remove_forces(system, force_indices):
        """Remove forces by index (in reverse order to keep indices valid)."""
        for idx in sorted(force_indices, reverse=True):
            system.removeForce(idx)
        logger.info("  Removed %d force(s)", len(force_indices))

    def log_energy_decomposition(context, system, label):
        """Log potential energy contributions for each force group.

        Parameters
        ----------
        context : openmm.Context
            OpenMM context used to evaluate energies.
        system : openmm.System
            OpenMM system containing force objects.
        label : str
            Label for the current decomposition snapshot.
        """
        original_groups = []
        for i in range(system.getNumForces()):
            force = system.getForce(i)
            original_groups.append(force.getForceGroup())
            force.setForceGroup(i)

        context.reinitialize(preserveState=True)

        logger.info("  Energy decomposition (%s):", label)
        for i in range(system.getNumForces()):
            state_i = context.getState(getEnergy=True, groups={i})
            energy_i = state_i.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)
            force_name = system.getForce(i).__class__.__name__
            logger.info("    Force %d (%s): %.2f kJ/mol", i, force_name, energy_i)

        for i in range(system.getNumForces()):
            system.getForce(i).setForceGroup(original_groups[i])

        context.reinitialize(preserveState=True)

    def get_platform():
        """Get the best available platform, validated by creating a test Context."""
        # Build a tiny test system to validate the platform actually works
        test_system = openmm.System()
        test_system.addParticle(1.0)
        test_force = openmm.CustomExternalForce("0")
        test_force.addParticle(0, [])
        test_system.addForce(test_force)

        def test_integrator_factory():
            """Create a test integrator for platform validation.

            Returns
            -------
            openmm.VerletIntegrator
                Minimal integrator used to probe platform availability.
            """

            return openmm.VerletIntegrator(0.001)

        for name in ["CUDA", "OpenCL", "CPU"]:
            try:
                plat = openmm.Platform.getPlatformByName(name)
                # Actually create a Context to verify the platform works
                test_ctx = openmm.Context(test_system, test_integrator_factory(), plat)
                del test_ctx
                logger.info("Using %s platform (validated)", name)
                return plat
            except Exception as e:
                logger.warning("Platform %s failed validation: %s", name, e)
                continue
        raise RuntimeError("No OpenMM platform available")

    platform = get_platform()

    # ── Step 4: Multi-stage energy minimization ────────────────────────────
    logger.info("=" * 70)
    logger.info("STEP 4: Multi-stage energy minimization")
    logger.info("=" * 70)
    t0 = time.perf_counter()

    logger.info("Running multi-stage minimization (tol=%.1f kJ/mol/nm)...", MINIMIZATION_TOLERANCE)

    # Initial unrestrained state
    integrator = openmm.VerletIntegrator(1.0 * omm_unit.femtosecond)
    sim = omm_app.Simulation(omm_topology, omm_system, integrator, platform)
    sim.context.setPositions(omm_positions)

    state_initial = sim.context.getState(getEnergy=True)
    e_initial = state_initial.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)
    logger.info("  E_initial = %.2f kJ/mol", e_initial)

    log_energy_decomposition(sim.context, omm_system, "before minimization")
    del sim

    # Stage A: restrain all conjugate heavy atoms
    logger.info("  Stage A: Restrained minimization (conjugate heavy atoms)")
    restrained_a = conjugate_heavy_indices
    restraint_idx_a = add_position_restraints(omm_system, restrained_a, omm_positions, RESTRAINT_K)

    integrator = openmm.VerletIntegrator(1.0 * omm_unit.femtosecond)
    sim = omm_app.Simulation(omm_topology, omm_system, integrator, platform)
    sim.context.setPositions(omm_positions)

    state_a_before = sim.context.getState(getEnergy=True)
    e_a_before = state_a_before.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)

    sim.minimizeEnergy(
        tolerance=MINIMIZATION_TOLERANCE * omm_unit.kilojoule_per_mole / omm_unit.nanometer,
        maxIterations=5000,
    )

    state_a_after = sim.context.getState(getEnergy=True, getPositions=True)
    e_a_after = state_a_after.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)
    positions_a = state_a_after.getPositions()
    logger.info("    E_stage_A: %.2f -> %.2f kJ/mol", e_a_before, e_a_after)

    omm_system.removeForce(restraint_idx_a)
    del sim

    # Stage B: restrain protein backbone atoms
    logger.info("  Stage B: Backbone-restrained minimization")
    restraint_idx_b = add_position_restraints(
        omm_system, backbone_indices, positions_a, RESTRAINT_K
    )

    integrator = openmm.VerletIntegrator(1.0 * omm_unit.femtosecond)
    sim = omm_app.Simulation(omm_topology, omm_system, integrator, platform)
    sim.context.setPositions(positions_a)

    state_b_before = sim.context.getState(getEnergy=True)
    e_b_before = state_b_before.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)

    sim.minimizeEnergy(
        tolerance=MINIMIZATION_TOLERANCE * omm_unit.kilojoule_per_mole / omm_unit.nanometer,
        maxIterations=5000,
    )

    state_b_after = sim.context.getState(getEnergy=True, getPositions=True)
    e_b_after = state_b_after.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)
    positions_b = state_b_after.getPositions()
    logger.info("    E_stage_B: %.2f -> %.2f kJ/mol", e_b_before, e_b_after)

    omm_system.removeForce(restraint_idx_b)
    del sim

    # Stage C: fully unrestrained minimization
    logger.info("  Stage C: Unrestrained minimization (%d iter)", MINIMIZATION_MAX_ITER)
    integrator = openmm.VerletIntegrator(1.0 * omm_unit.femtosecond)
    sim = omm_app.Simulation(omm_topology, omm_system, integrator, platform)
    sim.context.setPositions(positions_b)

    state_c_before = sim.context.getState(getEnergy=True)
    e_c_before = state_c_before.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)

    sim.minimizeEnergy(
        tolerance=MINIMIZATION_TOLERANCE * omm_unit.kilojoule_per_mole / omm_unit.nanometer,
        maxIterations=MINIMIZATION_MAX_ITER,
    )

    state_c_after = sim.context.getState(getEnergy=True, getPositions=True)
    e_c_after = state_c_after.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)
    positions_min = state_c_after.getPositions()
    box_vectors = state_c_after.getPeriodicBoxVectors()
    logger.info("    E_stage_C: %.2f -> %.2f kJ/mol", e_c_before, e_c_after)

    log_energy_decomposition(sim.context, omm_system, "after minimization")

    dt = time.perf_counter() - t0
    logger.info(
        "  Total minimization: %.2f -> %.2f kJ/mol (took %.1fs)",
        e_initial,
        e_c_after,
        dt,
    )
    e_after = e_c_after

    # Save minimized PDB
    min_pdb_path = OUTPUT_DIR / "solvated_minimized.pdb"
    with open(min_pdb_path, "w") as f:
        omm_app.PDBFile.writeFile(omm_topology, positions_min, f)
    logger.info("  Saved: %s", min_pdb_path)

    # Delete minimization simulation (can't reuse with a different integrator)
    del sim

    # ── Reference positions for restraints ─────────────────────────────────
    reference_positions = positions_min

    # Current state tracking
    current_positions = positions_min
    current_velocities = None
    current_box_vectors = box_vectors

    # ═══════════════════════════════════════════════════════════════════════
    # Step 5: Multi-stage equilibration
    # ═══════════════════════════════════════════════════════════════════════
    logger.info("=" * 70)
    logger.info("STEP 5: Multi-stage equilibration")
    logger.info("=" * 70)

    # ── Stage 1: NVT Heating ──────────────────────────────────────────────
    logger.info("-" * 50)
    logger.info(
        "Stage 1: NVT Heating %.0f -> %.0f K (%.3f ns)",
        HEATING_TEMP_START,
        HEATING_TEMP_END,
        HEATING_DURATION_NS,
    )
    logger.info("  Restraints: protein_heavy + polymer_heavy @ %.0f kJ/mol/nm^2", RESTRAINT_K)
    logger.info("-" * 50)

    t0 = time.perf_counter()

    # Remove any existing barostat (NVT)
    forces_to_remove = []
    for i in range(omm_system.getNumForces()):
        if isinstance(omm_system.getForce(i), openmm.MonteCarloBarostat):
            forces_to_remove.append(i)
    remove_forces(omm_system, forces_to_remove) if forces_to_remove else None

    # Add restraints on protein + polymer heavy atoms
    all_restrained = conjugate_heavy_indices
    restraint_idx_s1 = add_position_restraints(
        omm_system,
        all_restrained,
        reference_positions,
        RESTRAINT_K,
    )

    # Create integrator at starting temperature
    integrator = openmm.LangevinMiddleIntegrator(
        HEATING_TEMP_START * omm_unit.kelvin,
        FRICTION / omm_unit.picosecond,
        TIMESTEP * omm_unit.femtosecond,
    )

    sim = omm_app.Simulation(omm_topology, omm_system, integrator, platform)
    if current_box_vectors is not None:
        sim.context.setPeriodicBoxVectors(*current_box_vectors)
    sim.context.setPositions(current_positions)
    sim.context.setVelocitiesToTemperature(HEATING_TEMP_START * omm_unit.kelvin)

    # Set up reporters
    stage1_dir = OUTPUT_DIR / "equilibration_0_heating"
    stage1_dir.mkdir(exist_ok=True)

    total_steps_s1 = int(HEATING_DURATION_NS * 1e6 / TIMESTEP)
    report_interval_s1 = max(1, total_steps_s1 // 10)

    sim.reporters.append(
        omm_app.StateDataReporter(
            str(stage1_dir / "state_data.csv"),
            report_interval_s1,
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
    sim.reporters.append(
        omm_app.DCDReporter(str(stage1_dir / "trajectory.dcd"), report_interval_s1)
    )

    # Temperature ramping
    steps_per_update = int(HEATING_TEMP_INTERVAL_FS / TIMESTEP)
    current_temp = HEATING_TEMP_START
    steps_done = 0

    while current_temp < HEATING_TEMP_END:
        integrator.setTemperature(current_temp * omm_unit.kelvin)
        steps_to_run = min(steps_per_update, total_steps_s1 - steps_done)
        if steps_to_run <= 0:
            break
        sim.step(steps_to_run)
        steps_done += steps_to_run
        current_temp += HEATING_TEMP_INCREMENT

    # Run remaining steps at final temperature
    integrator.setTemperature(HEATING_TEMP_END * omm_unit.kelvin)
    remaining = total_steps_s1 - steps_done
    if remaining > 0:
        sim.step(remaining)
        steps_done += remaining

    state = sim.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    current_positions = state.getPositions()
    current_velocities = state.getVelocities()
    current_box_vectors = state.getPeriodicBoxVectors()
    e_s1 = state.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)

    dt = time.perf_counter() - t0
    logger.info("  Stage 1 complete: %d steps, E=%.2f kJ/mol, took %.1fs", steps_done, e_s1, dt)

    # Save checkpoint
    sim.saveCheckpoint(str(stage1_dir / "checkpoint.chk"))

    # Remove restraints and simulation
    remove_forces(omm_system, [restraint_idx_s1])
    del sim

    # ── Stage 2: NVT Polymer Relaxation ───────────────────────────────────
    logger.info("-" * 50)
    logger.info(
        "Stage 2: NVT Polymer Relaxation at %.0f K (%.3f ns)",
        TARGET_TEMP,
        POLYMER_RELAX_DURATION_NS,
    )
    logger.info("  Restraints: protein_heavy only @ %.0f kJ/mol/nm^2", RESTRAINT_K)
    logger.info("-" * 50)

    t0 = time.perf_counter()

    # Ensure no barostat (NVT)
    forces_to_remove = []
    for i in range(omm_system.getNumForces()):
        if isinstance(omm_system.getForce(i), openmm.MonteCarloBarostat):
            forces_to_remove.append(i)
    remove_forces(omm_system, forces_to_remove) if forces_to_remove else None

    # Add restraints on protein heavy atoms only (polymer free)
    restraint_idx_s2 = add_position_restraints(
        omm_system,
        protein_heavy_indices,
        reference_positions,
        RESTRAINT_K,
    )

    integrator = openmm.LangevinMiddleIntegrator(
        TARGET_TEMP * omm_unit.kelvin,
        FRICTION / omm_unit.picosecond,
        TIMESTEP * omm_unit.femtosecond,
    )

    sim = omm_app.Simulation(omm_topology, omm_system, integrator, platform)
    if current_box_vectors is not None:
        sim.context.setPeriodicBoxVectors(*current_box_vectors)
    sim.context.setPositions(current_positions)
    sim.context.setVelocities(current_velocities)

    stage2_dir = OUTPUT_DIR / "equilibration_1_polymer_relaxation"
    stage2_dir.mkdir(exist_ok=True)

    total_steps_s2 = int(POLYMER_RELAX_DURATION_NS * 1e6 / TIMESTEP)
    report_interval_s2 = max(1, total_steps_s2 // 10)

    sim.reporters.append(
        omm_app.StateDataReporter(
            str(stage2_dir / "state_data.csv"),
            report_interval_s2,
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
    sim.reporters.append(
        omm_app.DCDReporter(str(stage2_dir / "trajectory.dcd"), report_interval_s2)
    )

    logger.info("  Running %d steps...", total_steps_s2)
    sim.step(total_steps_s2)

    state = sim.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    current_positions = state.getPositions()
    current_velocities = state.getVelocities()
    current_box_vectors = state.getPeriodicBoxVectors()
    e_s2 = state.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)

    dt = time.perf_counter() - t0
    logger.info("  Stage 2 complete: %d steps, E=%.2f kJ/mol, took %.1fs", total_steps_s2, e_s2, dt)

    sim.saveCheckpoint(str(stage2_dir / "checkpoint.chk"))
    remove_forces(omm_system, [restraint_idx_s2])
    del sim

    # ── Stage 3: NPT Free Equilibration ──────────────────────────────────
    logger.info("-" * 50)
    logger.info(
        "Stage 3: NPT Free Equilibration at %.0f K, %.1f atm (%.3f ns)",
        TARGET_TEMP,
        PRESSURE,
        FREE_EQUIL_DURATION_NS,
    )
    logger.info("  Restraints: none")
    logger.info("-" * 50)

    t0 = time.perf_counter()

    # Add barostat for NPT
    barostat = openmm.MonteCarloBarostat(
        PRESSURE * omm_unit.atmosphere,
        TARGET_TEMP * omm_unit.kelvin,
        BAROSTAT_FREQ,
    )
    omm_system.addForce(barostat)
    logger.info(
        "  Added MC barostat: %.1f atm, %.0f K, freq=%d", PRESSURE, TARGET_TEMP, BAROSTAT_FREQ
    )

    integrator = openmm.LangevinMiddleIntegrator(
        TARGET_TEMP * omm_unit.kelvin,
        FRICTION / omm_unit.picosecond,
        TIMESTEP * omm_unit.femtosecond,
    )

    sim = omm_app.Simulation(omm_topology, omm_system, integrator, platform)
    if current_box_vectors is not None:
        sim.context.setPeriodicBoxVectors(*current_box_vectors)
    sim.context.setPositions(current_positions)
    sim.context.setVelocities(current_velocities)

    stage3_dir = OUTPUT_DIR / "equilibration_2_free"
    stage3_dir.mkdir(exist_ok=True)

    total_steps_s3 = int(FREE_EQUIL_DURATION_NS * 1e6 / TIMESTEP)
    report_interval_s3 = max(1, total_steps_s3 // 10)

    sim.reporters.append(
        omm_app.StateDataReporter(
            str(stage3_dir / "state_data.csv"),
            report_interval_s3,
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
    sim.reporters.append(
        omm_app.DCDReporter(str(stage3_dir / "trajectory.dcd"), report_interval_s3)
    )

    logger.info("  Running %d steps...", total_steps_s3)
    sim.step(total_steps_s3)

    state = sim.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    current_positions = state.getPositions()
    current_velocities = state.getVelocities()
    current_box_vectors = state.getPeriodicBoxVectors()
    e_s3 = state.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)

    dt = time.perf_counter() - t0
    logger.info("  Stage 3 complete: %d steps, E=%.2f kJ/mol, took %.1fs", total_steps_s3, e_s3, dt)

    sim.saveCheckpoint(str(stage3_dir / "checkpoint.chk"))

    # Save equilibrated PDB
    eq_pdb_path = OUTPUT_DIR / "equilibrated.pdb"
    with open(eq_pdb_path, "w") as f:
        omm_app.PDBFile.writeFile(omm_topology, current_positions, f)
    logger.info("  Saved equilibrated PDB: %s", eq_pdb_path)

    del sim
    # Note: keep the barostat for production (also NPT)

    # ═══════════════════════════════════════════════════════════════════════
    # Step 6: Production NPT
    # ═══════════════════════════════════════════════════════════════════════
    logger.info("=" * 70)
    logger.info(
        "STEP 6: Production NPT (%.3f ns, %.0f K, %.1f atm)",
        PRODUCTION_DURATION_NS,
        TARGET_TEMP,
        PRESSURE,
    )
    logger.info("=" * 70)

    t0 = time.perf_counter()

    # Update barostat temperature if needed
    for i in range(omm_system.getNumForces()):
        force = omm_system.getForce(i)
        if isinstance(force, openmm.MonteCarloBarostat):
            force.setDefaultTemperature(TARGET_TEMP * omm_unit.kelvin)
            force.setDefaultPressure(PRESSURE * omm_unit.atmosphere)
            force.setFrequency(BAROSTAT_FREQ)

    integrator = openmm.LangevinMiddleIntegrator(
        TARGET_TEMP * omm_unit.kelvin,
        FRICTION / omm_unit.picosecond,
        TIMESTEP * omm_unit.femtosecond,
    )

    sim = omm_app.Simulation(omm_topology, omm_system, integrator, platform)
    if current_box_vectors is not None:
        sim.context.setPeriodicBoxVectors(*current_box_vectors)
    sim.context.setPositions(current_positions)
    sim.context.setVelocities(current_velocities)

    prod_dir = OUTPUT_DIR / "production_0"
    prod_dir.mkdir(exist_ok=True)

    total_steps_prod = int(PRODUCTION_DURATION_NS * 1e6 / TIMESTEP)
    report_interval_prod = max(1, total_steps_prod // PRODUCTION_SAMPLES)

    sim.reporters.append(
        omm_app.StateDataReporter(
            str(prod_dir / "state_data.csv"),
            report_interval_prod,
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
    sim.reporters.append(
        omm_app.DCDReporter(str(prod_dir / "trajectory.dcd"), report_interval_prod)
    )
    sim.reporters.append(
        omm_app.CheckpointReporter(str(prod_dir / "checkpoint.chk"), report_interval_prod)
    )

    # Log initial state
    _state = sim.context.getState(getEnergy=True)
    _e = _state.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)
    logger.info("  Production initial PE: %.2f kJ/mol", _e)

    # Save topology PDB
    with open(prod_dir / "topology.pdb", "w") as f:
        omm_app.PDBFile.writeFile(omm_topology, current_positions, f)

    logger.info("  Running %d steps (%d frames)...", total_steps_prod, PRODUCTION_SAMPLES)
    sim.step(total_steps_prod)

    state = sim.context.getState(getPositions=True, getEnergy=True)
    e_prod = state.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)
    final_positions = state.getPositions()

    dt = time.perf_counter() - t0
    logger.info(
        "  Production complete: %d steps, E=%.2f kJ/mol, took %.1fs", total_steps_prod, e_prod, dt
    )

    # Save final state
    sim.saveCheckpoint(str(prod_dir / "final_checkpoint.chk"))

    # Save final PDB
    final_pdb_path = OUTPUT_DIR / "production_final.pdb"
    with open(final_pdb_path, "w") as f:
        omm_app.PDBFile.writeFile(omm_topology, final_positions, f)
    logger.info("  Saved final PDB: %s", final_pdb_path)

    del sim

    return {
        "e_minimized": e_after,
        "e_stage1": e_s1,
        "e_stage2": e_s2,
        "e_stage3": e_s3,
        "e_production_final": e_prod,
    }


# ═══════════════════════════════════════════════════════════════════════════
# Step 7: Validation
# ═══════════════════════════════════════════════════════════════════════════
def validate_results(energies):
    """Parse state data CSVs and validate simulation stability."""
    logger.info("=" * 70)
    logger.info("STEP 7: Validation")
    logger.info("=" * 70)

    import csv

    all_ok = True

    for stage_name in [
        "equilibration_0_heating",
        "equilibration_1_polymer_relaxation",
        "equilibration_2_free",
        "production_0",
    ]:
        csv_path = OUTPUT_DIR / stage_name / "state_data.csv"
        if not csv_path.exists():
            logger.warning("  Missing: %s", csv_path)
            continue

        with open(csv_path, "r") as f:
            # OpenMM StateDataReporter prefixes the header with '#'
            # Strip the '#' from the header but keep it; skip other comment lines
            lines = []
            for line in f:
                if line.startswith("#"):
                    # Could be the header row (e.g. #"Step","Time (ps)",...)
                    stripped = line.lstrip("#")
                    if stripped.strip():
                        lines.append(stripped)
                else:
                    lines.append(line)

        if not lines:
            logger.warning("  Empty state data: %s", stage_name)
            continue

        reader = csv.DictReader(lines)
        temps = []
        pe_vals = []
        densities = []
        volumes = []

        for row in reader:
            try:
                temps.append(float(row["Temperature (K)"]))
                pe_vals.append(float(row["Potential Energy (kJ/mole)"]))
                if "Density (g/mL)" in row:
                    densities.append(float(row["Density (g/mL)"]))
                if "Box Volume (nm^3)" in row:
                    volumes.append(float(row["Box Volume (nm^3)"]))
            except (KeyError, ValueError):
                continue

        if not temps:
            logger.warning("  No data in: %s", stage_name)
            continue

        temp_mean = np.mean(temps)
        temp_std = np.std(temps)
        pe_mean = np.mean(pe_vals)
        pe_std = np.std(pe_vals)

        logger.info("  %s:", stage_name)
        logger.info("    Temperature: %.1f +/- %.1f K (n=%d)", temp_mean, temp_std, len(temps))
        logger.info("    Potential E: %.0f +/- %.0f kJ/mol", pe_mean, pe_std)

        if densities:
            d_mean = np.mean(densities)
            d_std = np.std(densities)
            logger.info("    Density: %.4f +/- %.4f g/mL", d_mean, d_std)

        if volumes:
            v_mean = np.mean(volumes)
            v_std = np.std(volumes)
            logger.info("    Volume: %.1f +/- %.1f nm^3", v_mean, v_std)

        # Check for NaN or Inf in energies
        if any(np.isnan(pe_vals)) or any(np.isinf(pe_vals)):
            logger.error("    FAIL: NaN/Inf detected in potential energy!")
            all_ok = False

        # Check temperature stability (production only)
        if "production" in stage_name:
            if abs(temp_mean - TARGET_TEMP) > 20:
                logger.warning("    WARNING: Temperature drift > 20K from target")
            if temp_std > 30:
                logger.warning("    WARNING: Temperature fluctuations > 30K")

    logger.info("-" * 50)
    logger.info("Energy summary:")
    for key, val in energies.items():
        logger.info("  %s: %.2f kJ/mol", key, val)

    if all_ok:
        logger.info("VALIDATION PASSED: All energies finite, simulation stable")
    else:
        logger.error("VALIDATION FAILED: Check errors above")

    return all_ok


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════
def main():
    """Run the full workflow."""
    total_t0 = time.perf_counter()

    logger.info("=" * 70)
    logger.info("CONJUGATE SIMULATION WORKFLOW")
    logger.info("=" * 70)
    logger.info("Output directory: %s", OUTPUT_DIR)

    # Step 1: Build conjugate
    interchange, conjugate_off, free_polymer_offs, component_metadata, defs = (
        run_conjugation_notebook()
    )

    # Step 2: Solvate
    solvated_topology = solvate_system(interchange, conjugate_off, free_polymer_offs)

    # Step 3: Create solvated Interchange
    solvated_interchange = create_solvated_interchange(
        solvated_topology,
        conjugate_off,
        free_polymer_offs,
    )

    # Steps 4-6: Minimize, equilibrate, produce
    energies = run_openmm_simulation(
        solvated_interchange,
        conjugate_off,
        free_polymer_offs,
        component_metadata,
    )

    # Step 7: Validate
    success = validate_results(energies)

    total_dt = time.perf_counter() - total_t0
    logger.info("=" * 70)
    logger.info("TOTAL TIME: %.1f seconds (%.1f minutes)", total_dt, total_dt / 60)
    logger.info("=" * 70)

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
