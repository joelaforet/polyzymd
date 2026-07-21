"""
GROMACS export functionality for PolyzyMD.

This module provides comprehensive GROMACS export capabilities including:
- MDP file generation from PolyzyMD config parameters
- Position restraint file generation for equilibration
- Topology modification with #ifdef POSRES blocks
- Automated run script generation

The generated files are designed to match OpenFF Interchange force field
parameters (rcoulomb=0.9, rvdw=0.9, PME, etc.) for 1:1 parity with OpenMM.

Made by PolyzyMD, by Joseph R. Laforet Jr.
"""

from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

from polyzymd.core.branding import FULL_CREDIT_LINE, SHORT_CREDIT_LINE
from polyzymd.exporters.exact_openmm import (
    AtomIdentity,
    ExactTopologyBondRecord,
    gromacs_atom_order_hash,
    gromacs_topology_hash,
)
from polyzymd.utils.templates import render_package_template
from polyzymd.workflow.slurm import _validate_script_value

if TYPE_CHECKING:
    from openff.interchange import Interchange
    from openmm.app import Topology as OpenMMTopology

    from polyzymd.config.schema import (
        EquilibrationStageConfig,
        SimulationConfig,
        SimulationPhaseConfig,
    )
    from polyzymd.core.atom_groups import AtomGroupResolver, SystemComponentInfo

logger = logging.getLogger(__name__)

_EXACT_PAIR_TOLERANCE = 1e-12

# =============================================================================
# Constants and Mappings
# =============================================================================

# Thermostat mapping: PolyzyMD -> GROMACS
THERMOSTAT_MAP = {
    "LangevinMiddle": "v-rescale",
    "Langevin": "v-rescale",
    "NoseHoover": "nose-hoover",
    "Andersen": "andersen",
}

# Barostat mapping: PolyzyMD -> GROMACS (tcoupl, pcoupltype)
BAROSTAT_MAP = {
    "MC": ("c-rescale", "isotropic"),
    "MCA": ("c-rescale", "anisotropic"),
}

# OpenFF-compatible default parameters for 1:1 parity with OpenMM
OPENFF_DEFAULTS = {
    "cutoff-scheme": "verlet",
    "pbc": "xyz",
    "verlet-buffer-tolerance": 0.005,
    "coulombtype": "PME",
    "rcoulomb": 0.9,
    "fourier-spacing": 0.12,
    "ewald-rtol": 1e-5,
    "vdwtype": "cutoff",
    "rvdw": 0.9,
    "vdw-modifier": "Potential-switch",
    "rvdw-switch": 0.8,
    "DispCorr": "EnerPres",
    "constraints": "h-bonds",
    "constraint_algorithm": "lincs",
    "nstlist": 20,
}

# Water compressibility at 300K (bar^-1)
WATER_COMPRESSIBILITY = 4.5e-5

# Branding
POLYZYMD_BRANDING = FULL_CREDIT_LINE
_GROMACS_TEMPLATE_PACKAGE = "polyzymd.engines.gromacs"
_LOCAL_RUN_TEMPLATE = "run_gromacs.sh.jinja"


def _validate_local_shell_token(value: str, field_name: str) -> str:
    """Reject unsafe local run-script values expanded as shell tokens.

    Parameters
    ----------
    value : str
        The configured value to validate.
    field_name : str
        Name of the field for error messages.

    Returns
    -------
    str
        The validated value unchanged.

    Raises
    ------
    ValueError
        If the value contains unsafe shell characters or whitespace.
    """

    _validate_script_value(value, field_name)
    if any(char.isspace() for char in value):
        raise ValueError(
            f"Local GROMACS script field '{field_name}' contains whitespace: {value!r}. "
            "Values expanded as shell tokens must not contain whitespace."
        )
    return value


# =============================================================================
# Data Classes
# =============================================================================


@dataclass
class MDPParameters:
    """Container for MDP file parameters.

    This dataclass holds all parameters needed to generate a GROMACS MDP file.
    Parameters are organized by category for clarity.
    """

    # Identification
    title: str = ""
    stage_type: str = ""  # em, eq, prod

    # Integration
    integrator: str = "md"
    dt: float = 0.002  # ps
    nsteps: int = 0

    # Output control
    nstxout: int = 0
    nstvout: int = 0
    nstlog: int = 5000
    nstenergy: int = 5000
    nstxout_compressed: int = 5000
    compressed_x_grps: str = "System"

    # Neighbor searching (OpenFF defaults)
    cutoff_scheme: str = "verlet"
    nstlist: int = 20
    pbc: str = "xyz"
    verlet_buffer_tolerance: float = 0.005

    # Electrostatics (OpenFF defaults)
    coulombtype: str = "PME"
    rcoulomb: float = 0.9
    fourier_spacing: float = 0.12
    ewald_rtol: float = 1e-5

    # Van der Waals (OpenFF defaults)
    vdwtype: str = "cutoff"
    rvdw: float = 0.9
    vdw_modifier: str = "Potential-switch"
    rvdw_switch: float = 0.8
    dispcorr: str = "EnerPres"

    # Constraints
    continuation: bool = False
    constraints: str = "h-bonds"
    constraint_algorithm: str = "lincs"

    # Temperature coupling
    tcoupl: str = "no"
    tc_grps: str = "System"
    tau_t: float = 0.5
    ref_t: float = 300.0

    # Pressure coupling
    pcoupl: str = "no"
    pcoupltype: str = "isotropic"
    tau_p: float = 5.0
    ref_p: float = 1.0
    compressibility: float = WATER_COMPRESSIBILITY

    # Velocity generation
    gen_vel: bool = False
    gen_temp: float = 300.0
    gen_seed: int = -1

    # Energy minimization specific
    emtol: float = 500.0
    emstep: float = 0.01

    # Annealing (for temperature ramping)
    annealing: str = "no"
    annealing_npoints: int = 0
    annealing_time: List[float] = field(default_factory=list)
    annealing_temp: List[float] = field(default_factory=list)

    # Position restraints
    define: str = ""  # e.g., "-DPOSRES_PROTEIN -DPOSRES_LIGAND"

    def to_mdp_string(self) -> str:
        """Convert parameters to MDP file format string.

        Returns:
            String containing the MDP file contents.
        """
        lines = []

        # Header comment
        lines.append(f"; {self.title}")
        lines.append(f"; Generated by PolyzyMD on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append(f"; {POLYZYMD_BRANDING}")
        lines.append("")

        if self.stage_type == "em":
            lines.extend(self._generate_em_section())
        else:
            lines.extend(self._generate_md_section())

        return "\n".join(lines)

    def _generate_em_section(self) -> List[str]:
        """Generate MDP content for energy minimization."""
        lines = []

        # Integration
        lines.append("; Run control")
        lines.append(f"integrator      = {self.integrator}")
        lines.append(f"emtol           = {self.emtol}")
        lines.append(f"emstep          = {self.emstep}")
        lines.append(f"nsteps          = {self.nsteps}")
        lines.append("")

        # Output
        lines.append("; Output control")
        lines.append("nstxout         = 0")
        lines.append(f"nstlog          = {self.nstlog}")
        lines.append(f"nstenergy       = {self.nstenergy}")
        lines.append("")

        # Neighbor searching
        lines.extend(self._generate_neighbor_section())

        # Electrostatics
        lines.extend(self._generate_electrostatics_section())

        # Van der Waals
        lines.extend(self._generate_vdw_section())

        # No temperature/pressure coupling
        lines.append("; No temperature/pressure coupling during minimization")
        lines.append("tcoupl          = no")
        lines.append("pcoupl          = no")

        return lines

    def _generate_md_section(self) -> List[str]:
        """Generate MDP content for MD simulations."""
        lines = []

        # Run control
        lines.append("; Run control")
        lines.append(f"integrator      = {self.integrator}")
        lines.append(f"dt              = {self.dt}")
        lines.append(f"nsteps          = {self.nsteps}")
        lines.append("")

        # Output control
        lines.append("; Output control")
        lines.append(f"nstxout         = {self.nstxout}")
        lines.append(f"nstvout         = {self.nstvout}")
        lines.append(f"nstlog          = {self.nstlog}")
        lines.append(f"nstenergy       = {self.nstenergy}")
        lines.append(f"nstxout-compressed = {self.nstxout_compressed}")
        if self.compressed_x_grps:
            lines.append(f"compressed-x-grps = {self.compressed_x_grps}")
        lines.append("")

        # Constraints
        lines.append("; Bond constraints")
        lines.append(f"continuation    = {'yes' if self.continuation else 'no'}")
        lines.append(f"constraints     = {self.constraints}")
        lines.append(f"constraint_algorithm = {self.constraint_algorithm}")
        lines.append("")

        # Neighbor searching
        lines.extend(self._generate_neighbor_section())

        # Electrostatics
        lines.extend(self._generate_electrostatics_section())

        # Van der Waals
        lines.extend(self._generate_vdw_section())

        # Temperature coupling
        lines.append("; Temperature coupling")
        lines.append(f"tcoupl          = {self.tcoupl}")
        lines.append(f"tc-grps         = {self.tc_grps}")
        lines.append(f"tau_t           = {self.tau_t}")
        lines.append(f"ref_t           = {self.ref_t}")
        lines.append("")

        # Annealing (if enabled)
        if self.annealing != "no" and self.annealing_npoints > 0:
            lines.append("; Temperature annealing")
            lines.append(f"annealing       = {self.annealing}")
            lines.append(f"annealing-npoints = {self.annealing_npoints}")
            lines.append(f"annealing-time  = {' '.join(str(t) for t in self.annealing_time)}")
            lines.append(f"annealing-temp  = {' '.join(str(t) for t in self.annealing_temp)}")
            lines.append("")

        # Pressure coupling
        lines.append("; Pressure coupling")
        if self.pcoupl != "no":
            lines.append(f"pcoupl          = {self.pcoupl}")
            lines.append(f"pcoupltype      = {self.pcoupltype}")
            lines.append(f"tau_p           = {self.tau_p}")
            lines.append(f"ref_p           = {self.ref_p}")
            lines.append(f"compressibility = {self.compressibility}")
        else:
            lines.append("pcoupl          = no")
        lines.append("")

        # Velocity generation
        lines.append("; Velocity generation")
        if self.gen_vel:
            lines.append("gen_vel         = yes")
            lines.append(f"gen_temp        = {self.gen_temp}")
            lines.append(f"gen_seed        = {self.gen_seed}")
        else:
            lines.append("gen_vel         = no")

        # Position restraints define
        if self.define:
            lines.append("")
            lines.append("; Position restraints")
            lines.append(f"define          = {self.define}")

        return lines

    def _generate_neighbor_section(self) -> List[str]:
        """Generate neighbor searching parameters."""
        lines = [
            "; Neighbor searching",
            f"cutoff-scheme   = {self.cutoff_scheme}",
            f"nstlist         = {self.nstlist}",
            f"pbc             = {self.pbc}",
            f"verlet-buffer-tolerance = {self.verlet_buffer_tolerance}",
            "",
        ]
        return lines

    def _generate_electrostatics_section(self) -> List[str]:
        """Generate electrostatics parameters."""
        lines = [
            "; Electrostatics",
            f"coulombtype     = {self.coulombtype}",
            f"rcoulomb        = {self.rcoulomb}",
            f"fourier-spacing = {self.fourier_spacing}",
            f"ewald-rtol      = {self.ewald_rtol}",
            "",
        ]
        return lines

    def _generate_vdw_section(self) -> List[str]:
        """Generate van der Waals parameters."""
        lines = [
            "; Van der Waals",
            f"vdwtype         = {self.vdwtype}",
            f"rvdw            = {self.rvdw}",
            f"vdw-modifier    = {self.vdw_modifier}",
            f"rvdw-switch     = {self.rvdw_switch}",
            f"DispCorr        = {self.dispcorr}",
            "",
        ]
        return lines


# =============================================================================
# MDP Generator
# =============================================================================


class MDPGenerator:
    """Generates GROMACS MDP files from PolyzyMD configuration.

    This class handles the conversion of PolyzyMD simulation parameters to
    GROMACS MDP format, supporting:
    - Energy minimization
    - Multi-stage equilibration with position restraints
    - Temperature ramping via GROMACS annealing
    - Production MD

    The generated MDP files use OpenFF-compatible parameters for force field
    parity with OpenMM simulations.

    Example:
        >>> generator = MDPGenerator(config)
        >>> em_mdp = generator.generate_energy_minimization()
        >>> eq_mdps = generator.generate_equilibration_stages()
        >>> prod_mdp = generator.generate_production()
    """

    def __init__(self, config: "SimulationConfig"):
        """Initialize the MDP generator.

        Args:
            config: PolyzyMD SimulationConfig object containing all simulation
                parameters including thermodynamics, simulation phases, etc.
        """
        self._config = config
        self._temperature = config.thermodynamics.temperature
        self._pressure = config.thermodynamics.pressure

    def generate_energy_minimization(
        self,
        nsteps: int = 50000,
        emtol: float = 500.0,
    ) -> MDPParameters:
        """Generate MDP parameters for energy minimization.

        Args:
            nsteps: Maximum number of minimization steps.
            emtol: Energy tolerance (kJ/mol/nm) for convergence.

        Returns:
            MDPParameters configured for steepest descent minimization.
        """
        return MDPParameters(
            title="Energy Minimization - Run FIRST before dynamics",
            stage_type="em",
            integrator="steep",
            nsteps=nsteps,
            emtol=emtol,
            emstep=0.01,
            nstlog=500,
            nstenergy=500,
        )

    def generate_equilibration_stages(self) -> List[Tuple[str, MDPParameters]]:
        """Generate MDP parameters for all equilibration stages.

        Generates one MDP file per configured equilibration stage.

        Returns:
            List of (filename, MDPParameters) tuples. Filenames follow the
            pattern "eq_01_{name}.mdp", "eq_02_{name}.mdp", etc.
        """
        return self._generate_staged_equilibration(
            self._config.simulation_phases.equilibration_stages
        )

    def _generate_staged_equilibration(
        self,
        stages: List["EquilibrationStageConfig"],
    ) -> List[Tuple[str, MDPParameters]]:
        """Generate MDPs for multi-stage equilibration protocol.

        Args:
            stages: List of EquilibrationStageConfig objects.

        Returns:
            List of (filename, MDPParameters) tuples, one per stage.
        """
        result = []

        for i, stage in enumerate(stages):
            stage_num = i + 1
            is_first = i == 0

            # Build position restraint defines
            posres_defines = self._build_posres_defines(stage)

            # Handle temperature (constant or ramping)
            if stage.is_temperature_ramping:
                params = self._create_annealing_params(
                    stage=stage,
                    stage_num=stage_num,
                    is_first_stage=is_first,
                    posres_defines=posres_defines,
                )
            else:
                # Use stage temperature or fall back to global
                temperature = stage.temperature or self._temperature

                params = self._create_equilibration_params(
                    name=stage.name,
                    duration_ns=stage.duration,
                    samples=stage.samples,
                    time_step_fs=stage.time_step or 2.0,
                    temperature=temperature,
                    ensemble=stage.ensemble.value,
                    thermostat=(stage.thermostat.value if stage.thermostat else "LangevinMiddle"),
                    thermostat_timescale=stage.thermostat_timescale or 1.0,
                    barostat=stage.barostat.value if stage.barostat else None,
                    is_first_stage=is_first,
                    continuation=not is_first,
                    posres_defines=posres_defines,
                )

            filename = f"eq_{stage_num:02d}_{stage.name}.mdp"
            result.append((filename, params))

        return result

    def _create_equilibration_params(
        self,
        name: str,
        duration_ns: float,
        samples: int,
        time_step_fs: float,
        temperature: float,
        ensemble: str,
        thermostat: str,
        thermostat_timescale: float,
        barostat: Optional[str],
        is_first_stage: bool,
        continuation: bool,
        posres_defines: str = "",
    ) -> MDPParameters:
        """Create MDPParameters for an equilibration stage.

        Args:
            name: Stage name for title.
            duration_ns: Duration in nanoseconds.
            samples: Number of trajectory frames to save.
            time_step_fs: Time step in femtoseconds.
            temperature: Target temperature in Kelvin.
            ensemble: "NVT" or "NPT".
            thermostat: Thermostat type name.
            thermostat_timescale: Thermostat coupling timescale in ps.
            barostat: Barostat type name or None.
            is_first_stage: Whether this is the first equilibration stage.
            continuation: Whether continuing from previous simulation.
            posres_defines: Position restraint defines string.

        Returns:
            Configured MDPParameters object.
        """
        dt_ps = time_step_fs / 1000.0
        nsteps = int(duration_ns * 1e6 / time_step_fs)
        output_interval = max(1, nsteps // samples) if samples > 0 else 5000

        # Map thermostat
        tcoupl = THERMOSTAT_MAP.get(thermostat, "v-rescale")

        # Map barostat
        is_npt = ensemble == "NPT"
        if is_npt and barostat:
            pcoupl, pcoupltype = BAROSTAT_MAP.get(barostat, ("c-rescale", "isotropic"))
        else:
            pcoupl, pcoupltype = "no", "isotropic"

        # Pressure in bar (config is in atm)
        ref_p = self._pressure * 1.01325

        return MDPParameters(
            title=f"{ensemble} Equilibration: {name}",
            stage_type="eq",
            integrator="md",
            dt=dt_ps,
            nsteps=nsteps,
            nstxout=output_interval,
            nstvout=output_interval,
            nstlog=output_interval,
            nstenergy=output_interval,
            nstxout_compressed=output_interval,
            continuation=continuation,
            tcoupl=tcoupl,
            tau_t=thermostat_timescale,
            ref_t=temperature,
            pcoupl=pcoupl,
            pcoupltype=pcoupltype,
            ref_p=ref_p,
            gen_vel=is_first_stage,
            gen_temp=temperature,
            define=posres_defines,
        )

    def _create_annealing_params(
        self,
        stage: "EquilibrationStageConfig",
        stage_num: int,
        is_first_stage: bool,
        posres_defines: str = "",
    ) -> MDPParameters:
        """Create MDPParameters with GROMACS annealing for temperature ramping.

        Args:
            stage: EquilibrationStageConfig with temperature ramping.
            stage_num: Stage number (1-indexed).
            is_first_stage: Whether this is the first equilibration stage.
            posres_defines: Position restraint defines string.

        Returns:
            MDPParameters configured with annealing.
        """
        time_step_fs = stage.time_step or 2.0
        dt_ps = time_step_fs / 1000.0
        nsteps = int(stage.duration * 1e6 / time_step_fs)
        samples = stage.samples
        output_interval = max(1, nsteps // samples) if samples > 0 else 5000

        # Temperature ramping parameters
        t_start = stage.temperature_start
        t_end = stage.temperature_end
        duration_ps = stage.duration * 1000  # ns to ps

        # Map thermostat
        thermostat = stage.thermostat.value if stage.thermostat else "LangevinMiddle"
        tcoupl = THERMOSTAT_MAP.get(thermostat, "v-rescale")

        # Map barostat
        is_npt = stage.ensemble.value == "NPT"
        barostat = stage.barostat.value if stage.barostat else None
        if is_npt and barostat:
            pcoupl, pcoupltype = BAROSTAT_MAP.get(barostat, ("c-rescale", "isotropic"))
        else:
            pcoupl, pcoupltype = "no", "isotropic"

        ref_p = self._pressure * 1.01325

        # GROMACS annealing: simple linear ramp
        # annealing-time: times in ps
        # annealing-temp: temperatures at those times
        annealing_time = [0.0, duration_ps]
        annealing_temp = [t_start, t_end]

        return MDPParameters(
            title=f"Temperature Ramping: {stage.name} ({t_start}K -> {t_end}K)",
            stage_type="eq",
            integrator="md",
            dt=dt_ps,
            nsteps=nsteps,
            nstxout=output_interval,
            nstvout=output_interval,
            nstlog=output_interval,
            nstenergy=output_interval,
            nstxout_compressed=output_interval,
            continuation=not is_first_stage,
            tcoupl=tcoupl,
            tau_t=stage.thermostat_timescale or 0.5,
            ref_t=t_end,  # Reference temp is final temp
            pcoupl=pcoupl,
            pcoupltype=pcoupltype,
            ref_p=ref_p,
            gen_vel=is_first_stage,
            gen_temp=t_start,  # Generate velocities at start temp
            annealing="single",
            annealing_npoints=2,
            annealing_time=annealing_time,
            annealing_temp=annealing_temp,
            define=posres_defines,
        )

    def _build_posres_defines(self, stage: "EquilibrationStageConfig") -> str:
        """Build position restraint defines string from stage config.

        Args:
            stage: EquilibrationStageConfig with position_restraints.

        Returns:
            Define string like "-DPOSRES_PROTEIN -DPOSRES_LIGAND".
        """
        if not stage.position_restraints:
            return ""

        defines = []
        for posres in stage.position_restraints:
            # Convert group name to define macro
            # protein_heavy -> POSRES_PROTEIN
            # ligand_heavy -> POSRES_LIGAND
            # polymer_heavy -> POSRES_POLYMER
            group = posres.group.upper()
            if "PROTEIN" in group:
                defines.append("-DPOSRES_PROTEIN")
            elif "LIGAND" in group:
                defines.append("-DPOSRES_LIGAND")
            elif "POLYMER" in group:
                defines.append("-DPOSRES_POLYMER")

        # Remove duplicates while preserving order
        seen = set()
        unique_defines = []
        for d in defines:
            if d not in seen:
                seen.add(d)
                unique_defines.append(d)

        return " ".join(unique_defines)

    def generate_production(self) -> MDPParameters:
        """Generate MDP parameters for production MD.

        Returns:
            MDPParameters configured for production simulation.
        """
        prod = self._config.simulation_phases.production

        dt_ps = prod.time_step / 1000.0
        nsteps = int(prod.duration * 1e6 / prod.time_step)
        output_interval = max(1, nsteps // prod.samples) if prod.samples > 0 else 5000

        # Map thermostat
        thermostat = prod.thermostat.value if prod.thermostat else "LangevinMiddle"
        tcoupl = THERMOSTAT_MAP.get(thermostat, "v-rescale")

        # Map barostat
        is_npt = prod.ensemble.value == "NPT"
        barostat = prod.barostat.value if prod.barostat else None
        if is_npt and barostat:
            pcoupl, pcoupltype = BAROSTAT_MAP.get(barostat, ("c-rescale", "isotropic"))
        else:
            pcoupl, pcoupltype = "no", "isotropic"

        ref_p = self._pressure * 1.01325  # GROMACS uses bar internally, convert from atm

        return MDPParameters(
            title=f"Production MD ({prod.duration} ns)",
            stage_type="prod",
            integrator="md",
            dt=dt_ps,
            nsteps=nsteps,
            nstxout=0,  # Don't write full precision coords
            nstvout=0,  # Don't write velocities
            nstlog=output_interval,
            nstenergy=output_interval,
            nstxout_compressed=output_interval,
            compressed_x_grps="System",
            continuation=True,
            tcoupl=tcoupl,
            tau_t=prod.thermostat_timescale,
            ref_t=self._temperature,
            pcoupl=pcoupl,
            pcoupltype=pcoupltype,
            ref_p=ref_p,
            gen_vel=False,
        )


# =============================================================================
# Position Restraint Generator
# =============================================================================


class PositionRestraintGenerator:
    """Generates GROMACS position restraint sections for molecule ITP files.

    This class adds position restraint sections directly to molecule .itp files
    using local atom indices (1 to N within each molecule), which is required
    by GROMACS.

    For protein and ligand components (single molecule type each), it uses the
    AtomGroupResolver to identify restrained atoms by global index and maps
    them to local indices within that molecule's ITP.

    For polymers, OpenFF Interchange may create multiple unique molecule types
    (one per unique monomer sequence). This class discovers ALL polymer ITP
    files and adds position restraints to each by parsing the ITP's [ atoms ]
    section to identify heavy (non-hydrogen) atoms directly, rather than
    relying on global-to-local index mapping.

    GROMACS requires position restraints to be:
    1. Inside the [ moleculetype ] section of each molecule
    2. Using local atom indices (1 to N) within that molecule

    Example:
        >>> generator = PositionRestraintGenerator(topology, component_info)
        >>> generator.add_posres_to_itp_files(config, output_dir, "MySystem")
    """

    # Maps atom group names to (component_type, posres_define)
    GROUP_MAPPING: Dict[str, Tuple[str, str]] = {
        "protein_heavy": ("protein", "POSRES_PROTEIN"),
        "protein_backbone": ("protein", "POSRES_PROTEIN"),
        "protein_calpha": ("protein", "POSRES_PROTEIN"),
        "ligand_heavy": ("ligand", "POSRES_LIGAND"),
        "polymer_heavy": ("polymer", "POSRES_POLYMER"),
    }

    def __init__(
        self,
        topology: "OpenMMTopology",
        component_info: "SystemComponentInfo",
    ):
        """Initialize the position restraint generator.

        Args:
            topology: OpenMM Topology object from Interchange.
            component_info: SystemComponentInfo with atom counts/chain assignments.
        """
        from polyzymd.core.atom_groups import AtomGroupResolver

        self._topology = topology
        self._component_info = component_info
        self._resolver = AtomGroupResolver(topology, component_info)

    def add_posres_to_itp_files(
        self,
        config: "SimulationConfig",
        output_dir: Path,
        prefix: str,
    ) -> Dict[str, str]:
        """Add position restraint sections to molecule .itp files.

        For each restraint group in the config's equilibration stages:
        1. Determines the component type (protein, ligand, or polymer)
        2. Finds the relevant ITP file(s)
        3. Identifies restrained atoms (heavy atoms for *_heavy groups)
        4. Appends #ifdef POSRES blocks to each ITP

        For protein/ligand, there is exactly one ITP per component. For
        polymers, there may be many unique molecule types (one per unique
        monomer sequence), and ALL polymer ITPs receive position restraints.

        Args:
            config: SimulationConfig with equilibration stages.
            output_dir: Directory containing .itp files.
            prefix: Filename prefix used for .itp files.

        Returns:
            Dictionary mapping component types to POSRES define names used.
            E.g., {"protein": "POSRES_PROTEIN", "polymer": "POSRES_POLYMER"}
        """
        phases = config.simulation_phases
        posres_defines: Dict[str, str] = {}

        # Collect all unique (group, force_constant) pairs, keeping max fc
        restraints_needed: Dict[str, float] = {}
        for stage in phases.equilibration_stages:
            if not stage.position_restraints:
                continue
            for posres in stage.position_restraints:
                group = posres.group
                fc = posres.force_constant
                if group not in restraints_needed or fc > restraints_needed[group]:
                    restraints_needed[group] = fc

        if not restraints_needed:
            logger.info("No position restraints configured")
            return posres_defines

        # Process each restraint group
        for group_name, force_constant in restraints_needed.items():
            component_type, posres_define = self.GROUP_MAPPING.get(group_name, (None, None))
            if component_type is None or posres_define is None:
                logger.warning(f"Cannot determine component for group '{group_name}' - skipping")
                continue

            if component_type == "polymer":
                count = self._add_posres_to_polymer_itps(
                    output_dir, prefix, force_constant, posres_define, group_name
                )
                if count > 0:
                    posres_defines[component_type] = posres_define
            else:
                added = self._add_posres_to_single_itp(
                    output_dir,
                    prefix,
                    component_type,
                    group_name,
                    force_constant,
                    posres_define,
                )
                if added:
                    posres_defines[component_type] = posres_define

        return posres_defines

    # ------------------------------------------------------------------
    # Single-ITP path (protein / ligand)
    # ------------------------------------------------------------------

    def _add_posres_to_single_itp(
        self,
        output_dir: Path,
        prefix: str,
        component_type: str,
        group_name: str,
        force_constant: float,
        posres_define: str,
    ) -> bool:
        """Add position restraints to a single-molecule-type ITP (protein or ligand).

        Uses global-to-local index mapping via the AtomGroupResolver.

        Returns:
            True if restraints were added, False otherwise.
        """
        global_indices = self._resolver.resolve(group_name)
        if not global_indices:
            logger.warning(f"No atoms found for group '{group_name}' - skipping")
            return False

        itp_path = self._find_single_itp(output_dir, prefix, component_type)
        if itp_path is None:
            logger.warning(
                f"Cannot find ITP file for {component_type} - skipping position restraints"
            )
            return False

        local_indices = self._global_to_local_indices(global_indices, component_type)
        if not local_indices:
            logger.warning(f"No local indices for group '{group_name}' - skipping")
            return False

        self._append_posres_to_itp(
            itp_path, local_indices, force_constant, posres_define, group_name
        )
        logger.info(
            f"Added {len(local_indices)} position restraints to "
            f"{itp_path.name} (#ifdef {posres_define})"
        )
        return True

    def _find_single_itp(
        self,
        output_dir: Path,
        prefix: str,
        component_type: str,
    ) -> Optional[Path]:
        """Find the ITP file for a single-instance component (protein or ligand).

        OpenFF Interchange names molecule types MOL0, MOL1, ... based on
        PolyzyMD's building order:
        - MOL0 = protein
        - MOL1 = substrate/ligand (if protein present, else MOL0)

        Args:
            output_dir: Directory containing ITP files.
            prefix: Filename prefix.
            component_type: "protein" or "ligand"

        Returns:
            Path to the ITP file, or None if not found.
        """
        if component_type == "protein":
            mol_index = 0
        elif component_type == "ligand":
            mol_index = 1 if self._component_info.has_protein else 0
        else:
            return None

        itp_path = output_dir / f"{prefix}_MOL{mol_index}.itp"
        if itp_path.exists():
            return itp_path

        logger.warning(f"ITP file not found: {itp_path}")
        return None

    def _global_to_local_indices(
        self,
        global_indices: List[int],
        component_type: str,
    ) -> List[int]:
        """Convert global atom indices to local 1-indexed indices.

        Args:
            global_indices: 0-indexed global atom indices.
            component_type: "protein" or "ligand"

        Returns:
            Sorted list of 1-indexed local atom indices for GROMACS.
        """
        if component_type == "protein":
            start_index = 0
            n_atoms = self._component_info.n_protein_atoms
        elif component_type == "ligand":
            start_index = self._component_info.n_protein_atoms
            n_atoms = self._component_info.n_substrate_atoms
        else:
            return []

        end_index = start_index + n_atoms
        local_indices = [
            global_idx - start_index + 1
            for global_idx in global_indices
            if start_index <= global_idx < end_index
        ]
        return sorted(local_indices)

    # ------------------------------------------------------------------
    # Multi-ITP path (polymers)
    # ------------------------------------------------------------------

    def _add_posres_to_polymer_itps(
        self,
        output_dir: Path,
        prefix: str,
        force_constant: float,
        posres_define: str,
        group_name: str,
    ) -> int:
        """Add position restraints to ALL polymer ITP files.

        OpenFF Interchange creates one molecule type (and ITP file) per unique
        polymer sequence. For random copolymers, 25 chains might produce 12+
        unique molecule types. Every polymer ITP needs its own position
        restraints.

        Instead of mapping global indices (which requires knowing the exact
        molecule ordering), this method parses each ITP's [ atoms ] section
        directly and identifies heavy atoms by atom name (names NOT starting
        with 'H' are heavy atoms).

        Args:
            output_dir: Directory containing ITP files.
            prefix: Filename prefix.
            force_constant: Force constant in kJ/mol/nm^2.
            posres_define: POSRES define name (e.g., "POSRES_POLYMER").
            group_name: Atom group name for comment headers.

        Returns:
            Number of ITP files that received position restraints.
        """
        polymer_itps = self._find_all_polymer_itps(output_dir, prefix)
        if not polymer_itps:
            logger.warning("No polymer ITP files found - skipping polymer position restraints")
            return 0

        count = 0
        for itp_path in polymer_itps:
            heavy_indices = self._get_heavy_atom_indices_from_itp(itp_path)
            if not heavy_indices:
                logger.warning(f"No heavy atoms found in {itp_path.name} - skipping")
                continue

            self._append_posres_to_itp(
                itp_path, heavy_indices, force_constant, posres_define, group_name
            )
            count += 1
            logger.info(
                f"Added {len(heavy_indices)} position restraints to "
                f"{itp_path.name} (#ifdef {posres_define})"
            )

        logger.info(f"Position restraints added to {count}/{len(polymer_itps)} polymer ITP files")
        return count

    def _find_all_polymer_itps(
        self,
        output_dir: Path,
        prefix: str,
    ) -> List[Path]:
        """Discover all polymer molecule ITP files in the output directory.

        Polymer ITPs start at the first MOL index after protein and ligand.
        We glob for all ``{prefix}_MOL*.itp`` files, exclude known non-polymer
        molecules (water, ions), and return the rest sorted by MOL index.

        Args:
            output_dir: Directory containing ITP files.
            prefix: Filename prefix (e.g., "FnIII_OEGMA-SBMA-1-1").

        Returns:
            Sorted list of polymer ITP file paths.
        """
        # Determine the first polymer MOL index
        first_polymer_mol = 0
        if self._component_info.has_protein:
            first_polymer_mol += 1
        if self._component_info.has_substrate:
            first_polymer_mol += 1

        # Glob all MOL ITP files
        all_mol_itps = sorted(
            output_dir.glob(f"{prefix}_MOL*.itp"),
            key=lambda p: self._extract_mol_index(p) or 0,
        )

        # Keep only those at or above the first polymer index, excluding
        # water/ion molecule types (which have very few atoms, typically 3-4)
        polymer_itps = []
        for itp_path in all_mol_itps:
            mol_idx = self._extract_mol_index(itp_path)
            if mol_idx is None or mol_idx < first_polymer_mol:
                continue

            # Skip water/ion ITPs by checking molecule name and atom count
            mol_name = self._get_molecule_name_from_itp(itp_path)
            if mol_name and mol_name.lower() in ("water", "sol", "tip3p", "tip4p", "spce"):
                continue

            atom_count = self._get_atom_count_from_itp(itp_path)
            if atom_count <= 5:
                # Likely water (3 atoms) or ion (1 atom), not a polymer
                logger.debug(
                    f"Skipping {itp_path.name} ({atom_count} atoms) - too few atoms for a polymer"
                )
                continue

            polymer_itps.append(itp_path)

        return polymer_itps

    @staticmethod
    def _extract_mol_index(itp_path: Path) -> Optional[int]:
        """Extract the MOL index from a filename like ``prefix_MOL3.itp``.

        Returns:
            Integer MOL index, or None if the filename doesn't match.
        """
        match = re.search(r"_MOL(\d+)\.itp$", itp_path.name)
        return int(match.group(1)) if match else None

    @staticmethod
    def _get_molecule_name_from_itp(itp_path: Path) -> Optional[str]:
        """Parse the molecule name from the [ moleculetype ] section of an ITP.

        Returns:
            The molecule type name, or None if parsing fails.
        """
        try:
            in_moleculetype = False
            for line in itp_path.read_text().splitlines():
                stripped = line.strip()
                if not stripped or stripped.startswith(";"):
                    continue
                if stripped.startswith("["):
                    section = stripped.strip("[] \t").lower()
                    in_moleculetype = section == "moleculetype"
                    continue
                if in_moleculetype:
                    parts = stripped.split()
                    if parts:
                        return parts[0]
            return None
        except (OSError, UnicodeDecodeError) as exc:
            logger.warning(
                "Could not read molecule name from %s (%s). "
                "The exporter will fall back to filename-based molecule ordering.",
                itp_path,
                exc,
            )
            return None

    @staticmethod
    def _get_heavy_atom_indices_from_itp(itp_path: Path) -> List[int]:
        """Parse an ITP file and return 1-indexed local indices of heavy atoms.

        Heavy atoms are identified by atom name: any atom whose name does NOT
        start with 'H' is considered heavy. This convention is universal across
        AMBER, CHARMM, OPLS, and SMIRNOFF force fields. It is also invariant
        to hydrogen mass repartitioning (HMR), which changes atom masses but
        never changes atom names.

        The ITP [ atoms ] section format (from OpenFF Interchange) is::

            ;index, atom type, resnum, resname, name, cgnr, charge, mass
                 1 MOL2_0       1 RES1     C1      1  -0.123456789012  12.010780000000

        Column indices (0-based, whitespace-split):
        - 0: atom index (1-indexed)
        - 4: atom name

        Args:
            itp_path: Path to the ITP file.

        Returns:
            Sorted list of 1-indexed atom indices for non-hydrogen atoms.
        """
        heavy_indices: List[int] = []
        try:
            lines = itp_path.read_text().splitlines()
        except (OSError, UnicodeDecodeError) as exc:
            logger.warning(f"Failed to read ITP file {itp_path}: {exc}")
            return []

        in_atoms_section = False
        for line in lines:
            stripped = line.strip()
            if not stripped or stripped.startswith(";"):
                continue
            if stripped.startswith("["):
                section = stripped.strip("[] \t").lower()
                in_atoms_section = section == "atoms"
                continue
            if in_atoms_section:
                parts = stripped.split()
                if len(parts) >= 5 and parts[0].isdigit():
                    atom_index = int(parts[0])
                    atom_name = parts[4]
                    if not atom_name.startswith("H"):
                        heavy_indices.append(atom_index)
        return sorted(heavy_indices)

    # ------------------------------------------------------------------
    # Shared helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _get_atom_count_from_itp(itp_path: Path) -> int:
        """Parse an ITP file to count the number of atoms in the molecule.

        Reads the [ atoms ] section and counts data lines.

        Args:
            itp_path: Path to the ITP file.

        Returns:
            Number of atoms in the molecule type, or 0 on error.
        """
        try:
            lines = itp_path.read_text().splitlines()
        except (OSError, UnicodeDecodeError) as exc:
            logger.warning(f"Failed to read ITP file {itp_path}: {exc}")
            return 0

        in_atoms_section = False
        atom_count = 0
        for line in lines:
            stripped = line.strip()
            if not stripped or stripped.startswith(";"):
                continue
            if stripped.startswith("["):
                section = stripped.strip("[] \t").lower()
                in_atoms_section = section == "atoms"
                continue
            if in_atoms_section:
                parts = stripped.split()
                if parts and parts[0].isdigit():
                    atom_count += 1
        return atom_count

    @staticmethod
    def _append_posres_to_itp(
        itp_path: Path,
        local_indices: List[int],
        force_constant: float,
        posres_define: str,
        group_name: str,
    ) -> None:
        """Append position restraint section to an ITP file.

        Args:
            itp_path: Path to the ITP file.
            local_indices: List of 1-indexed local atom indices.
            force_constant: Force constant in kJ/mol/nm^2.
            posres_define: Name of the POSRES define (e.g., "POSRES_PROTEIN").
            group_name: Name of the atom group for comments.
        """
        content = itp_path.read_text()

        lines = [
            "",
            f"; Position restraints for {group_name}",
            f"; Generated by PolyzyMD - {POLYZYMD_BRANDING}",
            f"#ifdef {posres_define}",
            "[ position_restraints ]",
            "; ai   funct   fcx         fcy         fcz",
        ]

        for idx in local_indices:
            lines.append(
                f"{idx:6d}  1     {force_constant:.1f}    {force_constant:.1f}    {force_constant:.1f}"
            )

        lines.append("#endif")
        lines.append("")

        new_content = content.rstrip() + "\n" + "\n".join(lines)
        itp_path.write_text(new_content)


# =============================================================================
# Run Script Generator
# =============================================================================


class RunScriptGenerator:
    """Generates GROMACS run scripts for the full simulation workflow.

    Creates shell scripts that execute the complete GROMACS workflow:
    1. Energy minimization
    2. Equilibration stages
    3. Production MD
    4. Trajectory post-processing (PBC unwrapping)
    """

    def __init__(
        self,
        prefix: str,
        equilibration_mdps: List[str],
        gmx_command: str = "gmx",
    ):
        """Initialize the run script generator.

        Args:
            prefix: System prefix for file naming.
            equilibration_mdps: List of equilibration MDP filenames.
            gmx_command: GROMACS command (default "gmx", can be path).
        """
        self._prefix = prefix
        self._eq_mdps = equilibration_mdps
        self._gmx = gmx_command

    def generate(self, output_path: Path) -> Path:
        """Generate the GROMACS run script.

        Parameters
        ----------
        output_path : Path
            Path to write the script.

        Returns
        -------
        Path
            Path to the generated script.
        """
        self._validate_script_values()
        content = render_package_template(
            _GROMACS_TEMPLATE_PACKAGE,
            _LOCAL_RUN_TEMPLATE,
            {
                "branding": POLYZYMD_BRANDING,
                "short_credit_line": SHORT_CREDIT_LINE,
                "generated_on": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "gmx_command": self._gmx,
                "prefix": self._prefix,
                "energy_minimization_block": "\n".join(
                    self._generate_energy_minimization()
                ).rstrip(),
                "equilibration_block": "\n".join(self._generate_equilibration()).rstrip(),
                "production_block": "\n".join(self._generate_production()).rstrip(),
                "postprocessing_block": "\n".join(self._generate_postprocessing()).rstrip(),
                "footer_block": "\n".join(self._generate_footer()).rstrip(),
            },
        ).rstrip("\n")
        output_path.write_text(content)
        output_path.chmod(0o755)  # Make executable

        logger.info(f"Generated GROMACS run script: {output_path}")
        return output_path

    def _validate_script_values(self) -> None:
        """Validate values interpolated into the local GROMACS shell script.

        Raises
        ------
        ValueError
            If a configured value contains unsafe shell metacharacters.
        """
        _validate_local_shell_token(self._gmx, "gmx_command")
        _validate_local_shell_token(self._prefix, "prefix")
        for mdp_name in self._eq_mdps:
            _validate_local_shell_token(mdp_name, "equilibration_mdp")

    def _generate_header(self) -> List[str]:
        """Generate script header with configuration."""
        return [
            "#!/bin/bash",
            "# GROMACS Workflow Script",
            f"# {POLYZYMD_BRANDING}",
            f"# Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "#",
            "# Usage: ./run_gromacs.sh",
            "#",
            "# This script runs the complete MD workflow:",
            "#   1. Energy minimization",
            "#   2. Equilibration (multiple stages if configured)",
            "#   3. Production MD",
            "#   4. Trajectory post-processing",
            "",
            "set -e  # Exit on error",
            "",
            "# Configuration",
            f'GMX="{self._gmx}"',
            f'PREFIX="{self._prefix}"',
            "",
            'echo "========================================"',
            'echo "GROMACS Workflow for ${PREFIX}"',
            'echo "========================================"',
            f'echo "{SHORT_CREDIT_LINE}"',
            'echo "Using GROMACS: $GMX"',
            'echo ""',
            "",
        ]

    def _generate_energy_minimization(self) -> List[str]:
        """Generate energy minimization commands."""
        return [
            "# ========================================",
            "# Step 1: Energy Minimization",
            "# ========================================",
            'echo "=== Step 1: Energy Minimization ==="',
            "",
            "if [ ! -f em.gro ]; then",
            '    echo "=== Step 1: Energy Minimization ==="',
            "",
            "    $GMX grompp -f em.mdp -c ${PREFIX}.gro -r ${PREFIX}.gro -p ${PREFIX}.top -o em.tpr -maxwarn 1",
            "    $GMX mdrun -deffnm em -v",
            "",
            "    # Check if standard minimization succeeded",
            "    if [ ! -f em.gro ]; then",
            '        echo "ERROR: Energy minimization failed!"',
            "        exit 1",
            "    fi",
            "",
            "    # Verify energy minimization health",
            '    if grep -qi "force.*not finite\\|inf.*atom" em.log 2>/dev/null; then',
            '        echo "FATAL: Energy minimization failed — infinite forces detected in em.log"',
            '        echo "This usually indicates severe atomic overlaps in the initial structure."',
            '        echo "Try increasing packing padding or box size."',
            "        exit 1",
            "    fi",
            '    echo "Minimization complete: em.gro"',
            '    echo ""',
            "else",
            '    echo "Skipping energy minimization (em.gro exists)."',
            "fi",
            "",
        ]

    def _generate_equilibration(self) -> List[str]:
        """Generate equilibration stage commands."""
        lines = [
            "# ========================================",
            "# Step 2: Equilibration",
            "# ========================================",
        ]

        prev_output = "em"

        for i, mdp_file in enumerate(self._eq_mdps):
            stage_num = i + 1
            mdp_name = Path(mdp_file).stem  # Remove .mdp extension
            output_name = f"eq_{stage_num:02d}"

            lines.extend(
                [
                    f'echo "=== Step 2.{stage_num}: Equilibration - {mdp_name} ==="',
                    "",
                ]
            )

            # Build grompp command
            # Always use -r em.gro for reference coordinates (needed for position restraints)
            grompp_cmd = f"$GMX grompp -f {mdp_file} -c {prev_output}.gro -r em.gro"
            if prev_output != "em":
                grompp_cmd += f" -t {prev_output}.cpt"
            grompp_cmd += f" -p ${{PREFIX}}.top -o {output_name}.tpr -maxwarn 1"

            lines.extend(
                [
                    grompp_cmd,
                    f"$GMX mdrun -deffnm {output_name} -v",
                    "",
                    f"if [ ! -f {output_name}.gro ]; then",
                    f'    echo "ERROR: Equilibration stage {stage_num} failed!"',
                    "    exit 1",
                    "fi",
                    f'echo "Equilibration stage {stage_num} complete: {output_name}.gro"',
                    'echo ""',
                    "",
                ]
            )

            prev_output = output_name

        # Store last equilibration output for production
        lines.append(f'LAST_EQ="{prev_output}"')
        lines.append("")

        return lines

    def _generate_production(self) -> List[str]:
        """Generate production MD commands."""
        return [
            "# ========================================",
            "# Step 3: Production MD",
            "# ========================================",
            'echo "=== Step 3: Production MD ==="',
            "",
            "$GMX grompp -f prod.mdp -c ${LAST_EQ}.gro -t ${LAST_EQ}.cpt -p ${PREFIX}.top -o prod.tpr -maxwarn 1",
            "$GMX mdrun -deffnm prod -v",
            "",
            'echo "Production complete: prod.xtc"',
            'echo ""',
            "",
        ]

    def _generate_postprocessing(self) -> List[str]:
        """Generate trajectory post-processing commands."""
        return [
            "# ========================================",
            "# Step 4: Trajectory Post-processing",
            "# ========================================",
            'echo "=== Step 4: Trajectory Post-processing ==="',
            'echo "Removing PBC jumps and centering..."',
            "",
            "# Remove periodic boundary jumps",
            'echo "System" | $GMX trjconv -s prod.tpr -f prod.xtc -o prod_nojump.xtc -pbc nojump',
            "",
            "# Center and make molecules whole",
            'echo -e "Protein\\nSystem" | $GMX trjconv -s prod.tpr -f prod_nojump.xtc -o prod_centered.xtc -center -pbc mol -ur compact',
            "",
            'echo "Post-processing complete."',
            'echo "  - prod_nojump.xtc: Trajectory without PBC jumps"',
            'echo "  - prod_centered.xtc: Centered trajectory for visualization"',
            'echo ""',
            "",
        ]

    def _generate_footer(self) -> List[str]:
        """Generate script footer with summary."""
        return [
            "# ========================================",
            "# Workflow Complete",
            "# ========================================",
            'echo "========================================"',
            'echo "Workflow Complete!"',
            'echo "========================================"',
            'echo ""',
            'echo "Output files:"',
            'echo "  em.gro             - Minimized structure"',
            'echo "  eq_*.gro           - Equilibrated structures"',
            'echo "  prod.xtc           - Production trajectory"',
            'echo "  prod.edr           - Energy file"',
            'echo "  prod.log           - Log file"',
            'echo "  prod_centered.xtc  - Centered trajectory for visualization"',
            'echo ""',
            f'echo "{POLYZYMD_BRANDING}"',
        ]


# =============================================================================
# Main Exporter Class
# =============================================================================


def patch_gromacs_topology_with_exact_exceptions(top_path: Path, sidecar: Any) -> dict[str, Any]:
    """Patch a GROMACS topology with exact local exception semantics."""
    lines = top_path.read_text(encoding="utf-8").splitlines()
    topology = _parse_gromacs_topology_layout(lines)
    _preflight_exact_gromacs_identity(topology, sidecar)
    grouped = _group_exact_exceptions_by_molecule_type(topology, sidecar)
    raw_pair_set = _expand_local_pair_set(topology)
    exact_nonzero_set = _exact_pair_set(sidecar.nonzero_exceptions)
    if not raw_pair_set:
        raise RuntimeError("Exact GROMACS export failed closed: raw baseline has no [ pairs ] rows")
    if _exact_pair_set(sidecar.zero_exceptions) & raw_pair_set:
        raise RuntimeError("Exact GROMACS export failed closed: raw pairs include zero exclusions")
    if raw_pair_set != exact_nonzero_set:
        raise RuntimeError(
            "Exact GROMACS export failed closed: raw pair set does not match exact nonzero "
            f"exception set (raw={len(raw_pair_set)}, exact={len(exact_nonzero_set)})"
        )

    patched_lines = _replace_defaults(lines)
    patched_lines = _patch_local_exception_sections(patched_lines, topology, grouped, sidecar)
    top_path.write_text("\n".join(patched_lines) + "\n", encoding="utf-8")

    patched_topology = _parse_gromacs_topology_layout(patched_lines)
    validation = _validate_expanded_local_semantics(patched_topology, sidecar)
    return {
        "schema_version": 1,
        "route": "native_openmm_glycam_exact_gromacs",
        "topology": str(top_path),
        "raw_pair_count": len(raw_pair_set),
        "raw_local_pair_row_count": sum(
            len(block.sections.get("pairs", _TopologySection(-1, -1, ())).entries)
            for block in topology.molecule_types.values()
        ),
        "exact_nonzero_exception_count": len(sidecar.nonzero_exceptions),
        "exact_zero_exception_count": len(sidecar.zero_exceptions),
        "patched_pair_count": validation["expanded_pair_count"],
        "patched_local_pair_row_count": validation["local_pair_row_count"],
        "patched_exclusion_count": validation["expanded_exclusion_count"],
        "patched_local_exclusion_row_count": validation["local_exclusion_row_count"],
        "pair_mismatch_count": validation["pair_mismatch_count"],
        "exclusion_mismatch_count": validation["exclusion_mismatch_count"],
        "zero_pairs_in_patched_pairs": validation["zero_pairs_in_patched_pairs"],
        "constraint_count": len(sidecar.constraints),
        "exception_hash": sidecar.exception_hash,
        "atom_order_hash": sidecar.atom_order_hash,
        "gromacs_atom_order_hash": getattr(sidecar, "gromacs_atom_order_hash", None),
        "gromacs_topology_hash": getattr(sidecar, "gromacs_topology_hash", None),
        "route_invariants": dict(getattr(sidecar, "route_invariants", {}) or {}),
        "nonbonded_metadata": sidecar.nonbonded_metadata.model_dump(mode="json"),
    }


def _section_entries(lines: list[str], section: str) -> list[list[str]]:
    """Return non-comment token rows from a GROMACS topology section."""
    entries: list[list[str]] = []
    active = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("["):
            active = stripped.lower() == f"[ {section.lower()} ]"
            continue
        if active and stripped and not stripped.startswith(";"):
            entries.append(stripped.split(";", 1)[0].split())
    return entries


@dataclass(frozen=True)
class _TopologySection:
    """A local section inside one GROMACS molecule type."""

    header_index: int
    end_index: int
    entries: tuple[tuple[str, ...], ...]


@dataclass(frozen=True)
class _MoleculeTypeBlock:
    """A parsed GROMACS molecule type block."""

    name: str
    start_index: int
    end_index: int
    atom_count: int
    sections: dict[str, _TopologySection]


@dataclass(frozen=True)
class _AtomInstance:
    """A global atom mapped into GROMACS molecule-copy coordinates."""

    molecule_type: str
    copy_index: int
    local_index: int


@dataclass(frozen=True)
class _TopologyLayout:
    """Parsed molecule layout and global atom map for a GROMACS topology."""

    molecule_types: dict[str, _MoleculeTypeBlock]
    molecule_counts: dict[str, int]
    molecule_order: tuple[tuple[str, int], ...]
    atom_instances: dict[int, _AtomInstance]


@dataclass(frozen=True)
class _GroupedLocalExceptions:
    """Exact exception rows grouped by molecule type and local atom pair."""

    pairs: dict[str, dict[tuple[int, int], Any]]
    exclusions: dict[str, dict[tuple[int, int], Any]]


def _preflight_exact_gromacs_identity(topology: _TopologyLayout, sidecar: Any) -> None:
    """Validate parsed GROMACS atom order and topology before mutation."""

    atoms, bonds = _expanded_gromacs_atom_and_bond_identity(topology)
    _validate_parsed_gromacs_against_authoritative(atoms, bonds, sidecar)
    expected_atoms = tuple(getattr(sidecar, "gromacs_atoms", ()) or ())
    expected_bonds = tuple(getattr(sidecar, "gromacs_topology_bonds", ()) or ())
    if not expected_atoms and not expected_bonds:
        return
    if len(atoms) != sidecar.particle_count:
        raise RuntimeError(
            "Exact GROMACS export failed closed before patching: particle count mismatch "
            f"GROMACS={len(atoms)} sidecar={sidecar.particle_count}"
        )
    if tuple((a.index, a.name, a.residue_name, a.residue_id) for a in atoms) != tuple(
        (a.index, a.name, a.residue_name, a.residue_id) for a in expected_atoms
    ):
        raise RuntimeError(
            "Exact GROMACS export failed closed before patching: atom order/identity mismatch"
        )
    if gromacs_atom_order_hash(tuple(atoms)) != sidecar.gromacs_atom_order_hash:
        raise RuntimeError(
            "Exact GROMACS export failed closed before patching: atom-order hash mismatch"
        )
    if _topology_bond_pairs(bonds) != _topology_bond_pairs(expected_bonds):
        raise RuntimeError(
            "Exact GROMACS export failed closed before patching: topology bond identity mismatch"
        )
    if gromacs_topology_hash(tuple(atoms), tuple(bonds)) != sidecar.gromacs_topology_hash:
        raise RuntimeError(
            "Exact GROMACS export failed closed before patching: topology bond hash mismatch"
        )


def _validate_parsed_gromacs_against_authoritative(
    atoms: list[AtomIdentity], bonds: list[ExactTopologyBondRecord], sidecar: Any
) -> None:
    """Validate a raw GROMACS topology against authoritative OpenMM identity.

    GROMACS topology rows do not preserve chain identifiers. Interchange may also
    renumber repeated water residues and write common water residue-name aliases.
    The only accepted normalization is therefore: ignore chain IDs, require exact
    atom order and atom names, require exact non-water residue names and residue
    IDs, normalize known water residue names to ``HOH``, and ignore only water
    residue IDs.
    """

    if len(atoms) != sidecar.particle_count:
        raise RuntimeError(
            "Exact GROMACS export failed closed before patching: particle count mismatch "
            f"GROMACS={len(atoms)} sidecar={sidecar.particle_count}"
        )
    expected_atoms = _normalized_authoritative_atoms_for_gromacs(sidecar.atoms, atoms)
    parsed_atoms = tuple(_normalized_gromacs_atom(atom) for atom in atoms)
    if parsed_atoms != expected_atoms:
        raise RuntimeError(
            "Exact GROMACS export failed closed before patching: atom order/identity mismatch"
        )
    if _gromacs_explicit_bond_pairs(bonds, sidecar) != _authoritative_gromacs_explicit_bond_pairs(
        sidecar
    ):
        raise RuntimeError(
            "Exact GROMACS export failed closed before patching: topology bond identity mismatch"
        )


def _normalized_authoritative_atoms_for_gromacs(
    authoritative_atoms: tuple[AtomIdentity, ...], parsed_atoms: list[AtomIdentity]
) -> tuple[tuple[int, str, str, str], ...]:
    """Return authoritative atom identities under documented GROMACS normalization."""

    residue_id_offset = _gromacs_residue_id_offset(authoritative_atoms, parsed_atoms)
    return tuple(
        _normalized_authoritative_atom(authoritative, parsed, residue_id_offset=residue_id_offset)
        for authoritative, parsed in zip(authoritative_atoms, parsed_atoms, strict=True)
    )


def _normalized_authoritative_atom(
    authoritative: AtomIdentity, parsed: AtomIdentity, *, residue_id_offset: int | None = None
) -> tuple[int, str, str, str]:
    """Normalize one authoritative atom for comparison with parsed GROMACS rows."""

    if _is_water_residue(authoritative.residue_name) and _is_water_residue(parsed.residue_name):
        return (authoritative.index, authoritative.name, "HOH", "")
    if _is_standard_ion_residue(authoritative.residue_name) and _is_standard_ion_residue(
        parsed.residue_name
    ):
        return (authoritative.index, authoritative.name.upper(), authoritative.residue_name, "")
    residue_id = authoritative.residue_id
    if residue_id_offset is not None and _integer_string(residue_id):
        residue_id = str(int(residue_id) + residue_id_offset)
    return (
        authoritative.index,
        authoritative.name,
        authoritative.residue_name,
        residue_id,
    )


def _normalized_gromacs_atom(atom: AtomIdentity) -> tuple[int, str, str, str]:
    """Normalize one parsed GROMACS atom under the accepted export policy."""

    if _is_water_residue(atom.residue_name):
        return (atom.index, atom.name, "HOH", "")
    if _is_standard_ion_residue(atom.residue_name):
        return (atom.index, atom.name.upper(), atom.residue_name, "")
    return (atom.index, atom.name, atom.residue_name, atom.residue_id)


def _gromacs_residue_id_offset(
    authoritative_atoms: tuple[AtomIdentity, ...], parsed_atoms: list[AtomIdentity]
) -> int | None:
    """Return a uniform GROMACS residue-id offset when Interchange renumbers residues."""

    offsets = set()
    for authoritative, parsed in zip(authoritative_atoms, parsed_atoms, strict=True):
        if _is_water_residue(authoritative.residue_name) or _is_standard_ion_residue(
            authoritative.residue_name
        ):
            continue
        if authoritative.name != parsed.name or authoritative.residue_name != parsed.residue_name:
            return None
        if not (_integer_string(authoritative.residue_id) and _integer_string(parsed.residue_id)):
            return None
        offsets.add(int(parsed.residue_id) - int(authoritative.residue_id))
        if len(offsets) > 1:
            return None
    return offsets.pop() if offsets else None


def _integer_string(value: str) -> bool:
    """Return whether a string contains an integer value."""

    try:
        int(value)
    except ValueError:
        return False
    return True


def _is_water_residue(residue_name: str) -> bool:
    """Return whether a residue name is an accepted water alias."""

    return residue_name.upper() in {"HOH", "WAT", "SOL", "TIP3", "TIP3P"}


def _is_standard_ion_residue(residue_name: str) -> bool:
    """Return whether a residue name is a standard monoatomic ion alias."""

    return residue_name.upper() in {"NA", "NA+", "CL", "CL-", "K", "K+"}


def _topology_bond_pairs(
    bonds: list[ExactTopologyBondRecord] | tuple[ExactTopologyBondRecord, ...],
) -> set[tuple[int, int]]:
    """Return exact global topology bond pairs independent of section origin."""

    return {tuple(sorted((bond.i, bond.j))) for bond in bonds}


def _authoritative_gromacs_explicit_bond_pairs(sidecar: Any) -> set[tuple[int, int]]:
    """Return authoritative bond pairs expected as explicit GROMACS topology rows."""

    topology_pairs = _topology_bond_pairs(sidecar.topology_bonds)
    constrained_pairs = {tuple(sorted((record.i, record.j))) for record in sidecar.constraints}
    return topology_pairs - constrained_pairs


def _gromacs_explicit_bond_pairs(
    bonds: list[ExactTopologyBondRecord] | tuple[ExactTopologyBondRecord, ...], sidecar: Any
) -> set[tuple[int, int]]:
    """Return parsed GROMACS pairs excluding constraints already proven by OpenMM."""

    constrained_pairs = {tuple(sorted((record.i, record.j))) for record in sidecar.constraints}
    return _topology_bond_pairs(bonds) - constrained_pairs


def _expanded_gromacs_atom_and_bond_identity(
    topology: _TopologyLayout,
) -> tuple[list[AtomIdentity], list[ExactTopologyBondRecord]]:
    """Expand parsed GROMACS atom rows and topology bonds across molecule copies."""

    atoms: list[AtomIdentity] = []
    bonds: list[ExactTopologyBondRecord] = []
    for molecule_type, copies, copy_start in _expanded_molecule_rows_with_ordinals(
        topology.molecule_order
    ):
        block = topology.molecule_types[molecule_type]
        atom_rows = block.sections.get("atoms", _TopologySection(-1, -1, ())).entries
        bond_rows = block.sections.get("bonds", _TopologySection(-1, -1, ())).entries
        constraint_rows = block.sections.get("constraints", _TopologySection(-1, -1, ())).entries
        settle_rows = block.sections.get("settles", _TopologySection(-1, -1, ())).entries
        for _ in range(copies):
            offset = len(atoms)
            local_atoms: dict[int, AtomIdentity] = {}
            for row in atom_rows:
                if len(row) < 5:
                    raise RuntimeError("Exact GROMACS export found malformed [ atoms ] row")
                local_index = int(row[0])
                atom = AtomIdentity(
                    index=offset + local_index,
                    name=row[4],
                    residue_name=row[3],
                    residue_id=row[2],
                    chain_id="",
                )
                local_atoms[local_index] = atom
                atoms.append(atom)
            for row in (*bond_rows, *constraint_rows):
                if len(row) < 2:
                    continue
                left, right = sorted((int(row[0]), int(row[1])))
                bonds.append(
                    ExactTopologyBondRecord(
                        i=offset + left,
                        j=offset + right,
                        atom_i=local_atoms[left],
                        atom_j=local_atoms[right],
                    )
                )
            for row in settle_rows:
                if len(row) < 1:
                    continue
                first = int(row[0])
                for other in (first + 1, first + 2):
                    if other not in local_atoms:
                        continue
                    bonds.append(
                        ExactTopologyBondRecord(
                            i=offset + first,
                            j=offset + other,
                            atom_i=local_atoms[first],
                            atom_j=local_atoms[other],
                        )
                    )
    return atoms, bonds


def _sidecar_with_parsed_gromacs_identity(top_path: Path, sidecar: Any) -> Any:
    """Return a sidecar augmented with validated GROMACS-visible identity."""

    topology = _parse_gromacs_topology_layout(top_path.read_text(encoding="utf-8").splitlines())
    atoms, bonds = _expanded_gromacs_atom_and_bond_identity(topology)
    _validate_parsed_gromacs_against_authoritative(atoms, bonds, sidecar)
    return sidecar.model_copy(
        update={
            "gromacs_atoms": tuple(atoms),
            "gromacs_topology_bonds": tuple(bonds),
            "gromacs_atom_order_hash": gromacs_atom_order_hash(tuple(atoms)),
            "gromacs_topology_hash": gromacs_topology_hash(tuple(atoms), tuple(bonds)),
        }
    )


def _parse_gromacs_topology_layout(lines: list[str]) -> _TopologyLayout:
    """Parse molecule types, atom counts, molecule copies, and global atom map."""
    blocks = _parse_molecule_type_blocks(lines)
    molecule_order = tuple(_parse_molecules_from_lines(lines))
    if not blocks or not molecule_order:
        raise RuntimeError("Exact GROMACS export requires molecule types and [ molecules ] layout")
    counts: dict[str, int] = {}
    atom_instances: dict[int, _AtomInstance] = {}
    global_index = 1
    for molecule_type, copies, copy_start in _expanded_molecule_rows_with_ordinals(molecule_order):
        block = blocks.get(molecule_type)
        if block is None:
            raise RuntimeError(f"[ molecules ] references unknown molecule type {molecule_type!r}")
        counts[molecule_type] = counts.get(molecule_type, 0) + copies
        for copy_offset in range(copies):
            copy_index = copy_start + copy_offset
            for local_index in range(1, block.atom_count + 1):
                atom_instances[global_index] = _AtomInstance(
                    molecule_type=molecule_type,
                    copy_index=copy_index,
                    local_index=local_index,
                )
                global_index += 1
    return _TopologyLayout(
        molecule_types=blocks,
        molecule_counts=counts,
        molecule_order=molecule_order,
        atom_instances=atom_instances,
    )


def _parse_molecule_type_blocks(lines: list[str]) -> dict[str, _MoleculeTypeBlock]:
    """Return molecule type blocks keyed by name."""
    starts = [idx for idx, line in enumerate(lines) if line.strip().lower() == "[ moleculetype ]"]
    blocks: dict[str, _MoleculeTypeBlock] = {}
    for ordinal, start in enumerate(starts):
        end = (
            starts[ordinal + 1]
            if ordinal + 1 < len(starts)
            else _system_section_index(lines, start)
        )
        name = _molecule_type_name(lines, start, end)
        sections = _sections_in_block(lines, start, end)
        atom_count = len(sections.get("atoms", _TopologySection(-1, -1, ())).entries)
        if atom_count <= 0:
            raise RuntimeError(f"Molecule type {name!r} has no parseable [ atoms ] rows")
        blocks[name] = _MoleculeTypeBlock(
            name=name,
            start_index=start,
            end_index=end,
            atom_count=atom_count,
            sections=sections,
        )
    return blocks


def _system_section_index(lines: list[str], start: int) -> int:
    """Return the section index ending the last molecule type block."""
    for idx in range(start + 1, len(lines)):
        if lines[idx].strip().lower() == "[ system ]":
            return idx
    return len(lines)


def _molecule_type_name(lines: list[str], start: int, end: int) -> str:
    """Return the first data token in a molecule type section."""
    for line in lines[start + 1 : end]:
        stripped = line.strip()
        if not stripped or stripped.startswith(";"):
            continue
        if stripped.startswith("["):
            break
        return stripped.split()[0]
    raise RuntimeError("Could not parse [ moleculetype ] name")


def _sections_in_block(lines: list[str], start: int, end: int) -> dict[str, _TopologySection]:
    """Return local sections in one molecule type block."""
    section_headers = [idx for idx in range(start + 1, end) if lines[idx].strip().startswith("[")]
    sections: dict[str, _TopologySection] = {}
    for ordinal, header in enumerate(section_headers):
        name = lines[header].strip().strip("[] \t").lower()
        section_end = section_headers[ordinal + 1] if ordinal + 1 < len(section_headers) else end
        entries = tuple(
            tuple(stripped.split(";", 1)[0].split())
            for stripped in (line.strip() for line in lines[header + 1 : section_end])
            if stripped and not stripped.startswith(";")
        )
        sections[name] = _TopologySection(header, section_end, entries)
    return sections


def _parse_molecules_from_lines(lines: list[str]) -> list[tuple[str, int]]:
    """Parse the global GROMACS [ molecules ] layout."""
    entries = _section_entries(lines, "molecules")
    layout = []
    for entry in entries:
        if len(entry) >= 2:
            layout.append((entry[0], int(entry[1])))
    return layout


def _expanded_molecule_rows_with_ordinals(
    molecule_order: tuple[tuple[str, int], ...] | list[tuple[str, int]],
) -> list[tuple[str, int, int]]:
    """Return molecule rows with cumulative per-type copy ordinal starts."""

    copy_ordinals: dict[str, int] = {}
    expanded: list[tuple[str, int, int]] = []
    for molecule_type, copies in molecule_order:
        start = copy_ordinals.get(molecule_type, 0)
        expanded.append((molecule_type, copies, start))
        copy_ordinals[molecule_type] = start + copies
    return expanded


def _replace_defaults(lines: list[str]) -> list[str]:
    """Return topology lines with ``gen-pairs`` disabled."""
    new_lines: list[str] = []
    index = 0
    while index < len(lines):
        line = lines[index]
        stripped = line.strip().lower()
        if stripped == "[ defaults ]":
            new_lines.append(line)
            index += 1
            while index < len(lines) and (
                not lines[index].strip() or lines[index].strip().startswith(";")
            ):
                new_lines.append(lines[index])
                index += 1
            if index >= len(lines):
                raise RuntimeError("Exact GROMACS export could not find [ defaults ] data row")
            parts = lines[index].split()
            if len(parts) < 5:
                raise RuntimeError("Exact GROMACS export found malformed [ defaults ] row")
            new_lines.append(f"{parts[0]:>6} {parts[1]:>5} no {parts[3]:>12} {parts[4]:>12}")
            index += 1
            continue
        new_lines.append(line)
        index += 1
    return new_lines


def _group_exact_exceptions_by_molecule_type(
    topology: _TopologyLayout, sidecar: Any
) -> _GroupedLocalExceptions:
    """Group exact global exceptions by molecule type and local pair."""
    pairs = _group_records_by_local_pair(topology, sidecar.nonzero_exceptions, sidecar)
    exclusions = _group_records_by_local_pair(
        topology, sidecar.zero_exceptions + sidecar.nonzero_exceptions, sidecar
    )
    return _GroupedLocalExceptions(pairs=pairs, exclusions=exclusions)


def _group_records_by_local_pair(
    topology: _TopologyLayout, records: tuple[Any, ...], sidecar: Any
) -> dict[str, dict[tuple[int, int], Any]]:
    """Return representative records grouped by molecule type and local pair."""
    grouped: dict[str, dict[tuple[int, int], Any]] = {}
    counts: dict[tuple[str, tuple[int, int]], int] = {}
    particles = {particle.index: particle for particle in sidecar.particles}
    signatures: dict[tuple[str, tuple[int, int]], tuple[Any, ...]] = {}
    for record in records:
        left = topology.atom_instances.get(record.i)
        right = topology.atom_instances.get(record.j)
        if left is None or right is None:
            raise RuntimeError("Exact sidecar atom index is outside the GROMACS topology")
        if left.molecule_type != right.molecule_type or left.copy_index != right.copy_index:
            raise RuntimeError(
                "Exact GROMACS export does not support inter-copy or inter-type exceptions: "
                f"{record.i}-{record.j}"
            )
        local_pair = tuple(sorted((left.local_index, right.local_index)))
        key = (left.molecule_type, local_pair)
        signature = _record_semantic_signature(record, particles)
        previous = signatures.setdefault(key, signature)
        if previous != signature:
            raise RuntimeError(
                "Repeated molecule type has inconsistent local exception semantics for "
                f"{left.molecule_type} {local_pair}"
            )
        grouped.setdefault(left.molecule_type, {})[local_pair] = record
        counts[key] = counts.get(key, 0) + 1
    for molecule_type, local_pair in counts:
        expected = topology.molecule_counts[molecule_type]
        observed = counts[(molecule_type, local_pair)]
        if observed != expected:
            raise RuntimeError(
                "Repeated molecule type local exception count mismatch for "
                f"{molecule_type} {local_pair}: observed={observed}, expected={expected}"
            )
    return grouped


def _record_semantic_signature(record: Any, particles: dict[int, Any]) -> tuple[Any, ...]:
    """Return a comparable local exception and particle signature."""
    left = particles[record.i]
    right = particles[record.j]
    return (
        f"{record.charge_product_e2:.16g}",
        f"{record.sigma_nm:.16g}",
        f"{record.epsilon_kj_mol:.16g}",
        f"{left.charge_e:.16g}",
        f"{right.charge_e:.16g}",
        f"{left.sigma_nm:.16g}",
        f"{right.sigma_nm:.16g}",
        f"{left.epsilon_kj_mol:.16g}",
        f"{right.epsilon_kj_mol:.16g}",
    )


def _patch_local_exception_sections(
    lines: list[str],
    topology: _TopologyLayout,
    grouped: _GroupedLocalExceptions,
    sidecar: Any,
) -> list[str]:
    """Replace each molecule type's local pairs and exclusions."""
    new_lines: list[str] = []
    index = 0
    blocks_by_start = {block.start_index: block for block in topology.molecule_types.values()}
    while index < len(lines):
        block = blocks_by_start.get(index)
        if block is None:
            new_lines.append(lines[index])
            index += 1
            continue
        new_lines.extend(_patched_molecule_type_block(lines, block, grouped, sidecar))
        index = block.end_index
    return new_lines


def _patched_molecule_type_block(
    lines: list[str],
    block: _MoleculeTypeBlock,
    grouped: _GroupedLocalExceptions,
    sidecar: Any,
) -> list[str]:
    """Return one molecule type block with exact local pair/exclusion sections."""
    pair_lines = _exact_pairs_section(grouped.pairs.get(block.name, {}), sidecar)
    exclusion_lines = _exact_exclusions_section(grouped.exclusions.get(block.name, {}))
    out: list[str] = []
    index = block.start_index
    first_section_index = min(
        (section.header_index for section in block.sections.values()), default=block.end_index
    )
    pair_inserted = False
    exclusion_inserted = False
    while index < block.end_index:
        stripped = lines[index].strip().lower()
        if block.start_index < index < first_section_index and _is_molecule_type_data_line(
            lines[index]
        ):
            out.append(_molecule_type_line_with_nrexcl_zero(lines[index]))
            index += 1
            continue
        if stripped == "[ pairs ]":
            out.extend(pair_lines)
            pair_inserted = True
            index += 1
            while index < block.end_index and not lines[index].strip().startswith("["):
                index += 1
            continue
        if stripped == "[ exclusions ]":
            out.extend(exclusion_lines)
            exclusion_inserted = True
            index += 1
            while index < block.end_index and not lines[index].strip().startswith("["):
                index += 1
            continue
        if pair_lines and not pair_inserted and _should_insert_pair_section_before(stripped):
            out.extend(pair_lines)
            pair_inserted = True
        out.append(lines[index])
        index += 1
    if pair_lines and not pair_inserted:
        out.extend(pair_lines)
    if exclusion_lines and not exclusion_inserted:
        out.extend(exclusion_lines)
    return out


def _should_insert_pair_section_before(section: str) -> bool:
    """Return whether exact pairs should precede this local section."""
    return section in {"[ bonds ]", "[ settles ]", "[ angles ]", "[ dihedrals ]"}


def _is_molecule_type_data_line(line: str) -> bool:
    """Return whether a line contains molecule type data tokens."""
    stripped = line.strip()
    return bool(stripped and not stripped.startswith(";") and not stripped.startswith("["))


def _molecule_type_line_with_nrexcl_zero(line: str) -> str:
    """Return a molecule type row with automatic bonded exclusions disabled."""
    data, separator, comment = line.partition(";")
    parts = data.split()
    if len(parts) < 2:
        raise RuntimeError("Exact GROMACS export found malformed [ moleculetype ] row")
    rendered = f"{parts[0]:<12} 0"
    if separator:
        return f"{rendered} ;{comment}"
    return rendered


def _exact_pairs_section(records: dict[tuple[int, int], Any], sidecar: Any) -> list[str]:
    """Return exact local GROMACS function-2 pair rows."""
    if not records:
        return []
    charges = {particle.index: particle.charge_e for particle in sidecar.particles}
    lines = [
        "[ pairs ]",
        "; ai aj funct fudgeQQ qi qj sigma epsilon ; exact OpenMM NonbondedForce exceptions",
    ]
    for local_pair, record in sorted(records.items()):
        qi = charges[record.i]
        qj = charges[record.j]
        if abs(qi * qj) < 1e-14:
            qi_out = record.charge_product_e2
            qj_out = 1.0
            fudge = 1.0
        else:
            qi_out = qi
            qj_out = qj
            fudge = record.charge_product_e2 / (qi * qj)
        lines.append(
            f"{local_pair[0]:6d} {local_pair[1]:6d} 2 {fudge: .16g} {qi_out: .16g} "
            f"{qj_out: .16g} {record.sigma_nm: .16g} {record.epsilon_kj_mol: .16g}"
        )
    return lines


def _exact_exclusions_section(records: dict[tuple[int, int], Any]) -> list[str]:
    """Return exact local GROMACS regular nonbonded exclusion rows."""
    if not records:
        return []
    lines = [
        "[ exclusions ]",
        "; ai aj ; exact OpenMM NonbondedForce exceptions",
    ]
    lines.extend(f"{pair[0]:6d} {pair[1]:6d}" for pair in sorted(records))
    return lines


def _expand_local_pairs(
    topology: _TopologyLayout, sidecar: Any
) -> dict[tuple[int, int], tuple[float, float, float]]:
    """Expand local GROMACS pair rows to global pair parameters."""
    expanded: dict[tuple[int, int], tuple[float, float, float]] = {}
    for molecule_type, copies, copy_start in _expanded_molecule_rows_with_ordinals(
        topology.molecule_order
    ):
        block = topology.molecule_types[molecule_type]
        entries = block.sections.get("pairs", _TopologySection(-1, -1, ())).entries
        for copy_offset in range(copies):
            copy_index = copy_start + copy_offset
            offset = _copy_global_offset(topology, molecule_type, copy_index)
            for entry in entries:
                if len(entry) < 8:
                    continue
                local_i, local_j = int(entry[0]), int(entry[1])
                pair = tuple(sorted((offset + local_i, offset + local_j)))
                fudge, qi, qj, sigma, epsilon = (float(value) for value in entry[3:8])
                expanded[pair] = (fudge * qi * qj, sigma, epsilon)
    return expanded


def _expand_local_pair_set(topology: _TopologyLayout) -> set[tuple[int, int]]:
    """Expand local GROMACS pair rows to global atom pairs without parameters."""
    expanded: set[tuple[int, int]] = set()
    for molecule_type, copies, copy_start in _expanded_molecule_rows_with_ordinals(
        topology.molecule_order
    ):
        block = topology.molecule_types[molecule_type]
        entries = block.sections.get("pairs", _TopologySection(-1, -1, ())).entries
        for copy_offset in range(copies):
            copy_index = copy_start + copy_offset
            offset = _copy_global_offset(topology, molecule_type, copy_index)
            for entry in entries:
                if len(entry) < 2:
                    continue
                expanded.add(tuple(sorted((offset + int(entry[0]), offset + int(entry[1])))))
    return expanded


def _expand_local_exclusions(topology: _TopologyLayout) -> set[tuple[int, int]]:
    """Expand local GROMACS exclusion rows to global atom pairs."""
    expanded: set[tuple[int, int]] = set()
    for molecule_type, copies, copy_start in _expanded_molecule_rows_with_ordinals(
        topology.molecule_order
    ):
        block = topology.molecule_types[molecule_type]
        entries = block.sections.get("exclusions", _TopologySection(-1, -1, ())).entries
        for copy_offset in range(copies):
            copy_index = copy_start + copy_offset
            offset = _copy_global_offset(topology, molecule_type, copy_index)
            for entry in entries:
                if len(entry) < 2:
                    continue
                first = int(entry[0])
                for other in entry[1:]:
                    expanded.add(tuple(sorted((offset + first, offset + int(other)))))
    return expanded


def _copy_global_offset(topology: _TopologyLayout, molecule_type: str, copy_index: int) -> int:
    """Return the zero-based global atom offset for a molecule copy."""
    offset = 0
    remaining = copy_index
    for current_type, copies in topology.molecule_order:
        atom_count = topology.molecule_types[current_type].atom_count
        if current_type == molecule_type:
            if remaining < copies:
                return offset + remaining * atom_count
            remaining -= copies
        offset += copies * atom_count
    raise RuntimeError(f"Could not resolve molecule copy {molecule_type}[{copy_index}]")


def _validate_expanded_local_semantics(topology: _TopologyLayout, sidecar: Any) -> dict[str, int]:
    """Validate patched local pairs and exclusions by global re-expansion."""
    expanded_pairs = _expand_local_pairs(topology, sidecar)
    expanded_exclusions = _expand_local_exclusions(topology)
    exact_pairs = _exact_pair_parameter_table(sidecar.nonzero_exceptions)
    exact_zero = _exact_pair_set(sidecar.zero_exceptions)
    exact_exclusions = exact_zero | set(exact_pairs)
    pair_mismatches = 0
    if set(expanded_pairs) != set(exact_pairs):
        pair_mismatches += len(set(expanded_pairs) ^ set(exact_pairs))
    for pair, expected in exact_pairs.items():
        got = expanded_pairs.get(pair)
        if got is None:
            continue
        if any(
            abs(got_value - expected_value) > _EXACT_PAIR_TOLERANCE
            for got_value, expected_value in zip(got, expected, strict=True)
        ):
            pair_mismatches += 1
    exclusion_mismatches = len(expanded_exclusions ^ exact_exclusions)
    zero_pairs_in_patched = len(exact_zero & set(expanded_pairs))
    if pair_mismatches or exclusion_mismatches or zero_pairs_in_patched:
        raise RuntimeError(
            "Exact GROMACS export failed closed after local semantic validation: "
            f"pair_mismatches={pair_mismatches}, "
            f"exclusion_mismatches={exclusion_mismatches}, "
            f"zero_pairs_in_patched={zero_pairs_in_patched}"
        )
    return {
        "expanded_pair_count": len(expanded_pairs),
        "expanded_exclusion_count": len(expanded_exclusions),
        "local_pair_row_count": sum(
            len(block.sections.get("pairs", _TopologySection(-1, -1, ())).entries)
            for block in topology.molecule_types.values()
        ),
        "local_exclusion_row_count": sum(
            len(block.sections.get("exclusions", _TopologySection(-1, -1, ())).entries)
            for block in topology.molecule_types.values()
        ),
        "pair_mismatch_count": pair_mismatches,
        "exclusion_mismatch_count": exclusion_mismatches,
        "zero_pairs_in_patched_pairs": zero_pairs_in_patched,
    }


def _exact_pair_set(records: tuple[Any, ...]) -> set[tuple[int, int]]:
    """Return global sorted atom pairs from exact exception records."""
    return {tuple(sorted((record.i, record.j))) for record in records}


def _exact_pair_parameter_table(
    records: tuple[Any, ...],
) -> dict[tuple[int, int], tuple[float, float, float]]:
    """Return exact global pair parameters keyed by sorted atom pair."""
    return {
        tuple(sorted((record.i, record.j))): (
            record.charge_product_e2,
            record.sigma_nm,
            record.epsilon_kj_mol,
        )
        for record in records
    }


def _exact_mdp_string(params: MDPParameters, sidecar: Any) -> str:
    """Return MDP text aligned narrowly to exact OpenMM nonbonded metadata."""
    metadata = sidecar.nonbonded_metadata
    params.rcoulomb = metadata.cutoff_nm
    params.rvdw = metadata.cutoff_nm
    params.ewald_rtol = metadata.ewald_error_tolerance
    params.coulombtype = "PME" if metadata.method.upper() == "PME" else metadata.method
    params.dispcorr = "EnerPres" if metadata.use_dispersion_correction else "no"
    if metadata.use_switching_function and metadata.switching_distance_nm is not None:
        params.vdw_modifier = "Potential-switch"
        params.rvdw_switch = metadata.switching_distance_nm
    else:
        params.vdw_modifier = "None"
    return params.to_mdp_string()


class GromacsExporter:
    """Coordinates full GROMACS export from PolyzyMD systems.

    This is the main entry point for GROMACS export. It coordinates:
    - Coordinate/topology export via Interchange
    - MDP file generation from config
    - Position restraint file generation
    - Topology modification for posres includes
    - Run script generation

    Example:
        >>> exporter = GromacsExporter(interchange, config, component_info)
        >>> exporter.export(output_dir, prefix="my_system")
    """

    def __init__(
        self,
        interchange: "Interchange",
        config: "SimulationConfig",
        component_info: Optional["SystemComponentInfo"] = None,
    ):
        """Initialize the GROMACS exporter.

        Args:
            interchange: OpenFF Interchange object with parameterized system, or
                a PolyzyMD ``ExactExportBundle`` for the native GLYCAM
                exact-exception compatibility bridge.
            config: PolyzyMD SimulationConfig with simulation parameters.
            component_info: Optional SystemComponentInfo for position restraints.
                If not provided, position restraints will be skipped.
        """
        self._interchange = interchange
        self._config = config
        self._component_info = component_info

    def export(
        self,
        output_dir: Union[str, Path],
        prefix: Optional[str] = None,
        gmx_command: str = "gmx",
    ) -> Dict[str, Path]:
        """Export the complete GROMACS simulation setup.

        Args:
            output_dir: Directory to write all files.
            prefix: Filename prefix. If None, generates from config.
            gmx_command: GROMACS command for run script.

        Returns:
            Dictionary mapping file types to paths:
            - "gro": Coordinate file
            - "top": Topology file
            - "em_mdp": Energy minimization MDP
            - "eq_mdps": List of equilibration MDP paths
            - "prod_mdp": Production MDP
            - "posres": Dictionary of position restraint files
            - "run_script": Run script path
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if prefix is None:
            prefix = self._generate_prefix()

        result: Dict[str, any] = {}

        if getattr(self._interchange, "is_exact_export_bundle", False):
            return self._export_exact_bundle(output_dir, prefix, gmx_command)

        # Step 1: Export coordinates and topology via Interchange
        logger.info("Exporting coordinates and topology...")
        gro_path, top_path = self._export_interchange(output_dir, prefix)
        result["gro"] = gro_path
        result["top"] = top_path

        # Step 2: Generate MDP files
        logger.info("Generating MDP files...")
        mdp_generator = MDPGenerator(self._config)

        # Energy minimization
        em_params = mdp_generator.generate_energy_minimization()
        em_path = output_dir / "em.mdp"
        em_path.write_text(em_params.to_mdp_string())
        result["em_mdp"] = em_path

        # Equilibration stages
        eq_mdps = mdp_generator.generate_equilibration_stages()
        eq_paths = []
        for filename, params in eq_mdps:
            mdp_path = output_dir / filename
            mdp_path.write_text(params.to_mdp_string())
            eq_paths.append(mdp_path)
        result["eq_mdps"] = eq_paths

        # Production
        prod_params = mdp_generator.generate_production()
        prod_path = output_dir / "prod.mdp"
        prod_path.write_text(prod_params.to_mdp_string())
        result["prod_mdp"] = prod_path

        # Step 3: Add position restraints to molecule ITP files (if component_info available)
        result["posres"] = {}
        result["posres_defines"] = {}
        if self._component_info is not None:
            logger.info("Adding position restraints to molecule ITP files...")
            topology = self._interchange.to_openmm_topology()
            posres_generator = PositionRestraintGenerator(topology, self._component_info)
            posres_defines = posres_generator.add_posres_to_itp_files(
                self._config, output_dir, prefix
            )
            result["posres_defines"] = posres_defines

            # Track which components have posres for user info
            if posres_defines:
                logger.info(f"Position restraints added for: {', '.join(posres_defines.keys())}")
                logger.info(f"POSRES defines: {', '.join(posres_defines.values())}")

        # Step 4: Generate run script
        logger.info("Generating run script...")
        eq_mdp_names = [p.name for p in eq_paths]
        script_generator = RunScriptGenerator(prefix, eq_mdp_names, gmx_command)
        script_path = output_dir / f"run_{prefix}_gromacs.sh"
        script_generator.generate(script_path)
        result["run_script"] = script_path

        # Log summary
        self._log_summary(result, output_dir)

        return result

    def _export_exact_bundle(
        self,
        output_dir: Path,
        prefix: str,
        gmx_command: str,
    ) -> Dict[str, Any]:
        """Export a native GLYCAM exact bundle with explicit OpenMM exceptions.

        This branch is a PolyzyMD compatibility bridge and not a vanilla
        Interchange export. A private baseline Interchange writes raw GROMACS
        files, then the topology is patched to ``gen-pairs = no`` with explicit
        function-2 rows for every nonzero OpenMM exception in the sidecar.
        """
        if self._component_info is not None:
            raise ValueError(
                "GROMACS exact native OpenMM GLYCAM export does not support component_info "
                "position-restraint postprocessing because it cannot be proven without changing "
                "exact topology semantics. Disable component_info or use the OpenFF/Sage route."
            )
        bundle = self._interchange
        baseline = bundle.require_private_baseline_interchange()
        original_interchange = self._interchange
        self._interchange = baseline
        try:
            self._fix_zero_indexed_residues()
        finally:
            self._interchange = original_interchange
        output_prefix = str(output_dir / prefix)
        baseline.to_gromacs(prefix=output_prefix, monolithic=True, _merge_atom_types=True)
        gro_path = output_dir / f"{prefix}.gro"
        top_path = output_dir / f"{prefix}.top"
        stub_mdp = output_dir / f"{prefix}.mdp"
        if stub_mdp.exists():
            stub_mdp.unlink()

        self._fix_gro_residue_numbering(gro_path, top_path, output_dir, prefix)

        sidecar = _sidecar_with_parsed_gromacs_identity(top_path, bundle.sidecar)
        if bundle.sidecar_path is not None:
            sidecar.save(bundle.sidecar_path)

        audit = patch_gromacs_topology_with_exact_exceptions(
            top_path=top_path,
            sidecar=sidecar,
        )
        audit_path = output_dir / f"{prefix}_exact_gromacs_audit.json"
        audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")

        mdp_generator = MDPGenerator(self._config)
        em_path = output_dir / "em.mdp"
        em_path.write_text(
            _exact_mdp_string(mdp_generator.generate_energy_minimization(), bundle.sidecar)
        )
        eq_paths = []
        for filename, params in mdp_generator.generate_equilibration_stages():
            mdp_path = output_dir / filename
            mdp_path.write_text(_exact_mdp_string(params, bundle.sidecar))
            eq_paths.append(mdp_path)
        prod_path = output_dir / "prod.mdp"
        prod_path.write_text(_exact_mdp_string(mdp_generator.generate_production(), bundle.sidecar))

        script_generator = RunScriptGenerator(prefix, [p.name for p in eq_paths], gmx_command)
        script_path = output_dir / f"run_{prefix}_gromacs.sh"
        script_generator.generate(script_path)

        result: Dict[str, Any] = {
            "gro": gro_path,
            "top": top_path,
            "em_mdp": em_path,
            "eq_mdps": eq_paths,
            "prod_mdp": prod_path,
            "posres": {},
            "posres_defines": {},
            "run_script": script_path,
            "exact_exception_sidecar": bundle.sidecar_path,
            "exact_gromacs_audit": audit_path,
        }
        self._log_summary(result, output_dir)
        return result

    def _export_interchange(
        self,
        output_dir: Path,
        prefix: str,
    ) -> Tuple[Path, Path]:
        """Export Interchange to GROMACS format.

        Args:
            output_dir: Output directory.
            prefix: Filename prefix.

        Returns:
            Tuple of (gro_path, top_path).
        """
        # Fix 0-indexed residues for GROMACS
        self._fix_zero_indexed_residues()

        # Export using Interchange
        output_prefix = str(output_dir / prefix)
        self._interchange.to_gromacs(
            prefix=output_prefix,
            monolithic=False,
            _merge_atom_types=True,
        )

        gro_path = output_dir / f"{prefix}.gro"
        top_path = output_dir / f"{prefix}.top"

        # Remove the stub MDP generated by Interchange
        stub_mdp = output_dir / f"{prefix}.mdp"
        if stub_mdp.exists():
            stub_mdp.unlink()

        # Fix include paths in topology file
        # Interchange writes absolute-looking paths relative to CWD, but GROMACS
        # needs paths relative to the .top file location
        self._fix_topology_include_paths(top_path, output_dir, prefix)

        # Fix residue numbering in .gro for multi-residue molecule copies.
        # OpenFF Interchange adds copy_index to each atom's residue_index,
        # which is correct for single-residue molecules (water) but creates
        # a sliding offset for multi-residue molecules (polymers).
        self._fix_gro_residue_numbering(gro_path, top_path, output_dir, prefix)

        return gro_path, top_path

    def _fix_topology_include_paths(
        self,
        top_path: Path,
        output_dir: Path,
        prefix: str,
    ) -> None:
        """Fix #include paths in topology file to be relative to .top location.

        OpenFF Interchange writes include paths relative to the current working
        directory at export time (e.g., "output_files/replicate_1/gromacs/file.itp").
        GROMACS expects paths relative to the .top file location, so we need to
        convert them to just the filename (e.g., "file.itp").

        Args:
            top_path: Path to the topology file.
            output_dir: Directory containing the topology file.
            prefix: Filename prefix used for the export.
        """
        import re

        content = top_path.read_text()

        # Pattern to match #include directives with paths containing our prefix
        # This matches: #include "any/path/prefix_something.itp"
        # And converts to: #include "prefix_something.itp"
        include_pattern = re.compile(r'#include\s+"([^"]*/' + re.escape(prefix) + r'_[^"]+\.itp)"')

        def fix_path(match: re.Match) -> str:
            full_path = match.group(1)
            # Extract just the filename
            filename = Path(full_path).name
            return f'#include "{filename}"'

        new_content = include_pattern.sub(fix_path, content)

        # Only write if we made changes
        if new_content != content:
            top_path.write_text(new_content)
            logger.info(f"Fixed include paths in topology file: {top_path.name}")

    def _fix_gro_residue_numbering(
        self,
        gro_path: Path,
        top_path: Path,
        output_dir: Path,
        prefix: str,
    ) -> None:
        """Fix residue numbering in .gro file for multi-residue molecule copies.

        OpenFF Interchange's GRO writer computes residue numbers as::

            (atom.residue_index + copy_index) % 100_000

        This is correct for single-residue molecules (water, ions) where each
        copy should have a unique residue number. For multi-residue molecules
        (polymers with N monomers), it creates a sliding offset:
        - Copy 0: resid 1, 2, 3, 4, 5
        - Copy 1: resid 2, 3, 4, 5, 6  (wrong — should be 6, 7, 8, 9, 10)
        - Copy 2: resid 3, 4, 5, 6, 7  (wrong — should be 11, 12, 13, 14, 15)

        This method post-processes the .gro to assign globally sequential
        residue numbers: each copy's residues continue from the previous copy.
        This enables unique residue-based selection in analysis tools like
        MDAnalysis (e.g., ``u.select_atoms("resid 11:15")`` for chain 3).

        For single-residue molecules, the original numbering is preserved since
        ``resid + copy_index`` is already correct (each copy gets one unique id).

        Args:
            gro_path: Path to the .gro file.
            top_path: Path to the .top file (for molecule layout).
            output_dir: Directory containing ITP files.
            prefix: Filename prefix for ITP files.
        """
        # Parse molecule layout from .top: [(mol_name, n_copies), ...]
        mol_layout = self._parse_molecules_from_top(top_path)
        if not mol_layout:
            logger.debug("No molecule layout found in .top - skipping residue fix")
            return

        # Build per-molecule-type info: atoms per molecule and residues per
        # molecule (number of unique residue indices in the template).
        mol_info: Dict[str, Tuple[int, int]] = {}  # name -> (n_atoms, n_residues)
        for mol_name, _ in mol_layout:
            if mol_name in mol_info:
                continue
            itp_path = output_dir / f"{prefix}_{mol_name}.itp"
            if not itp_path.exists():
                # Fallback: water/ions may not have separate ITP files
                mol_info[mol_name] = (0, 1)
                continue
            n_atoms, n_residues = self._parse_itp_atom_and_residue_counts(itp_path)
            mol_info[mol_name] = (n_atoms, n_residues)

        # Check if any molecule type has >1 residue (otherwise no fix needed)
        needs_fix = any(n_res > 1 for _, n_res in mol_info.values())
        if not needs_fix:
            logger.debug("No multi-residue molecule types - skipping residue fix")
            return

        # Read and fix the .gro file
        lines = gro_path.read_text().splitlines()
        if len(lines) < 3:
            return

        title = lines[0]
        n_atoms_line = lines[1]
        atom_lines = lines[2:-1]  # Everything between header and box vector
        box_line = lines[-1]

        n_atoms_total = int(n_atoms_line.strip())
        if len(atom_lines) != n_atoms_total:
            logger.warning(
                f"GRO atom count mismatch: header says {n_atoms_total}, "
                f"found {len(atom_lines)} lines - skipping residue fix"
            )
            return

        # Walk through atom lines, assigning residue numbers per molecule layout
        fixed_lines = []
        atom_idx = 0
        next_resid = 1  # Global residue counter (1-indexed)

        for mol_name, n_copies in mol_layout:
            n_atoms_per_mol, n_residues_per_mol = mol_info.get(mol_name, (0, 1))

            if n_atoms_per_mol == 0:
                # Unknown molecule type — count atoms from .gro lines by
                # looking for the next molecule type change. This shouldn't
                # happen with well-formed exports, but handle gracefully.
                logger.warning(
                    f"Unknown atom count for molecule type '{mol_name}' - "
                    "preserving original residue numbers"
                )
                # Skip n_copies worth of atoms with original numbering
                # We can't determine atom count, so just leave remaining lines
                break

            if n_residues_per_mol <= 1:
                # Single-residue molecule (water, ion, small molecule).
                # Each copy gets one unique residue number — the original
                # copy_index arithmetic is already correct for these.
                for copy_idx in range(n_copies):
                    for _ in range(n_atoms_per_mol):
                        if atom_idx < len(atom_lines):
                            fixed_lines.append(atom_lines[atom_idx])
                            atom_idx += 1
                    next_resid += 1
            else:
                # Multi-residue molecule (polymer). Read the template residue
                # indices from copy 0, then apply globally sequential numbering
                # to all copies.
                if atom_idx + n_atoms_per_mol > len(atom_lines):
                    logger.warning(
                        f"Not enough atoms for template of '{mol_name}' - stopping residue fix"
                    )
                    break

                # Parse template: extract per-atom residue index from copy 0
                template_resids: List[int] = []
                for i in range(n_atoms_per_mol):
                    line = atom_lines[atom_idx + i]
                    resid_str = line[:5].strip()
                    try:
                        template_resids.append(int(resid_str))
                    except ValueError:
                        template_resids.append(1)

                # Determine unique residue indices in template order
                seen_resids: List[int] = []
                for r in template_resids:
                    if r not in seen_resids:
                        seen_resids.append(r)
                template_min = seen_resids[0] if seen_resids else 1

                # Process all copies
                for copy_idx in range(n_copies):
                    # Offset for this copy: template residues are rebased to
                    # start at next_resid
                    resid_offset = next_resid - template_min

                    for atom_in_mol in range(n_atoms_per_mol):
                        if atom_idx >= len(atom_lines):
                            break
                        line = atom_lines[atom_idx]
                        new_resid = (template_resids[atom_in_mol] + resid_offset) % 100_000
                        # GRO format: first 5 chars are residue number
                        fixed_line = f"{new_resid:5d}{line[5:]}"
                        fixed_lines.append(fixed_line)
                        atom_idx += 1

                    next_resid += n_residues_per_mol

        # Append any remaining atoms (shouldn't happen, but be safe)
        while atom_idx < len(atom_lines):
            fixed_lines.append(atom_lines[atom_idx])
            atom_idx += 1

        # Write fixed .gro
        output_lines = [title, n_atoms_line] + fixed_lines + [box_line]
        gro_path.write_text("\n".join(output_lines) + "\n")

        logger.info(
            f"Fixed residue numbering in {gro_path.name} "
            f"(globally sequential for multi-residue molecules)"
        )

    @staticmethod
    def _parse_molecules_from_top(top_path: Path) -> List[Tuple[str, int]]:
        """Parse the [ molecules ] section from a GROMACS .top file.

        Returns:
            List of (molecule_name, n_copies) tuples in topology order.
        """
        molecules: List[Tuple[str, int]] = []
        try:
            lines = top_path.read_text().splitlines()
        except (OSError, UnicodeDecodeError) as exc:
            logger.warning(f"Failed to read molecules from {top_path}: {exc}")
            return molecules

        in_molecules = False
        for line in lines:
            stripped = line.strip()
            if not stripped or stripped.startswith(";"):
                continue
            if stripped.startswith("["):
                section = stripped.strip("[] \t").lower()
                in_molecules = section == "molecules"
                continue
            if in_molecules:
                parts = stripped.split()
                if len(parts) >= 2:
                    mol_name = parts[0]
                    try:
                        n_copies = int(parts[1])
                    except ValueError:
                        continue
                    molecules.append((mol_name, n_copies))
        return molecules

    @staticmethod
    def _parse_itp_atom_and_residue_counts(
        itp_path: Path,
    ) -> Tuple[int, int]:
        """Parse an ITP file to get atom count and unique residue count.

        Args:
            itp_path: Path to the ITP file.

        Returns:
            Tuple of (n_atoms, n_unique_residues).
        """
        n_atoms = 0
        residue_indices: set = set()
        try:
            lines = itp_path.read_text().splitlines()
        except (OSError, UnicodeDecodeError) as exc:
            logger.warning(f"Failed to read ITP file {itp_path}: {exc}")
            return 0, 1

        in_atoms = False
        for line in lines:
            stripped = line.strip()
            if not stripped or stripped.startswith(";"):
                continue
            if stripped.startswith("["):
                section = stripped.strip("[] \t").lower()
                in_atoms = section == "atoms"
                continue
            if in_atoms:
                parts = stripped.split()
                if len(parts) >= 3 and parts[0].isdigit():
                    n_atoms += 1
                    try:
                        residue_indices.add(int(parts[2]))
                    except ValueError:
                        pass
        return n_atoms, max(len(residue_indices), 1)

    def _fix_zero_indexed_residues(self) -> None:
        """Fix 0-indexed residues for GROMACS compatibility."""
        found_zero_indexed = False

        for molecule in self._interchange.topology.molecules:
            for atom in molecule.atoms:
                residue_num = atom.metadata.get("residue_number")
                if residue_num is not None:
                    if isinstance(residue_num, str):
                        try:
                            if int(residue_num) == 0:
                                found_zero_indexed = True
                                break
                        except ValueError:
                            pass
                    elif residue_num == 0:
                        found_zero_indexed = True
                        break
            if found_zero_indexed:
                break

        if not found_zero_indexed:
            return

        logger.info("Fixing 0-indexed residues for GROMACS compatibility")
        for molecule in self._interchange.topology.molecules:
            for atom in molecule.atoms:
                residue_num = atom.metadata.get("residue_number")
                if residue_num is not None:
                    if isinstance(residue_num, str):
                        try:
                            atom.metadata["residue_number"] = str(int(residue_num) + 1)
                        except ValueError:
                            pass
                    else:
                        atom.metadata["residue_number"] = residue_num + 1

    def _generate_prefix(self) -> str:
        """Generate filename prefix from config."""
        parts = []

        if self._config.enzyme:
            parts.append(self._config.enzyme.name)

        if self._config.polymers and self._config.polymers.enabled:
            parts.append(self._config.polymers.type_prefix)

        return "_".join(parts) if parts else "system"

    def _log_summary(self, result: Dict[str, any], output_dir: Path) -> None:
        """Log summary of generated files."""
        logger.info("")
        logger.info("=" * 60)
        logger.info("GROMACS Export Complete")
        logger.info("=" * 60)
        logger.info(f"Output directory: {output_dir}")
        logger.info("")
        logger.info("Generated files:")
        logger.info(f"  Coordinates:      {result['gro'].name}")
        logger.info(f"  Topology:         {result['top'].name}")
        logger.info(f"  Energy min MDP:   {result['em_mdp'].name}")
        for eq_path in result["eq_mdps"]:
            logger.info(f"  Equilibration:    {eq_path.name}")
        logger.info(f"  Production MDP:   {result['prod_mdp'].name}")
        if result.get("posres_defines"):
            logger.info("  Position restraints added to molecule ITP files:")
            for component, define in result["posres_defines"].items():
                logger.info(f"    - {component}: #ifdef {define}")
        logger.info(f"  Run script:       {result['run_script'].name}")
        logger.info("")
        logger.info(f"To run: cd {output_dir} && ./{result['run_script'].name}")
        logger.info("")
        logger.info(POLYZYMD_BRANDING)


# =============================================================================
# GROMACS Runner (Local Execution)
# =============================================================================


class GromacsError(Exception):
    """Exception raised when a GROMACS command fails."""

    def __init__(self, command: str, returncode: int, message: str = ""):
        self.command = command
        self.returncode = returncode
        self.message = message
        super().__init__(f"GROMACS command failed: {command} (exit code {returncode})")


class GromacsRunner:
    """Execute GROMACS workflow from exported files with streaming output.

    This class runs the complete GROMACS MD workflow locally, streaming
    all output to stdout in real-time for familiar GROMACS user experience.

    The workflow consists of:
        1. Energy minimization (grompp + mdrun)
        2. Equilibration stages (loop of grompp + mdrun)
        3. Production MD (grompp + mdrun)
        4. Trajectory post-processing (trjconv)

    On any failure, execution stops immediately and all intermediate files
    are preserved for debugging.

    Example:
        >>> runner = GromacsRunner(
        ...     working_dir=Path("./gromacs"),
        ...     prefix="lysozyme_PEG",
        ...     equilibration_mdps=["eq_01_nvt.mdp", "eq_02_npt.mdp"],
        ... )
        >>> runner.run_full_workflow()
    """

    def __init__(
        self,
        working_dir: Path,
        prefix: str,
        equilibration_mdps: List[str],
        gmx_command: str = "gmx",
    ):
        """Initialize the GROMACS runner.

        Args:
            working_dir: Directory containing GROMACS input files.
            prefix: System prefix for filenames (e.g., "lysozyme_PEG").
            equilibration_mdps: List of equilibration MDP filenames in order.
            gmx_command: Path to GROMACS executable (default: "gmx").
        """
        self._working_dir = Path(working_dir)
        self._prefix = prefix
        self._eq_mdps = equilibration_mdps
        self._gmx = gmx_command

        # Track state for workflow
        self._last_eq_output: Optional[str] = None

    def run_full_workflow(self) -> None:
        """Run the complete GROMACS workflow.

        Executes energy minimization, equilibration, production, and
        post-processing in sequence. Stops immediately on any failure.

        Raises:
            GromacsError: If any GROMACS command fails.
            FileNotFoundError: If required input files are missing.
        """
        self._print_banner("GROMACS Workflow", f"System: {self._prefix}")

        # Verify input files exist
        self._verify_inputs()

        # Run workflow steps
        self._run_energy_minimization()
        self._run_equilibration()
        self._run_production()
        self._run_postprocessing()

        self._print_summary()

    def _verify_inputs(self) -> None:
        """Verify that required input files exist."""
        required_files = [
            f"{self._prefix}.gro",
            f"{self._prefix}.top",
            "em.mdp",
            "prod.mdp",
        ]
        required_files.extend(self._eq_mdps)

        missing = []
        for filename in required_files:
            filepath = self._working_dir / filename
            if not filepath.exists():
                missing.append(filename)

        if missing:
            raise FileNotFoundError(
                f"Missing required files in {self._working_dir}: {', '.join(missing)}"
            )

    def _run_command(
        self,
        cmd: List[str],
        description: str,
        stdin_input: Optional[str] = None,
    ) -> None:
        """Run a command with streaming stdout/stderr.

        Args:
            cmd: Command and arguments as list.
            description: Human-readable description for logging.
            stdin_input: Optional string to pipe to stdin.

        Raises:
            GromacsError: If command returns non-zero exit code.
        """
        import subprocess
        import sys

        logger.info(f"\n>>> {description}")
        logger.info(f">>> Command: {' '.join(cmd)}")
        logger.info("-" * 60)

        try:
            process = subprocess.Popen(
                cmd,
                cwd=self._working_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                stdin=subprocess.PIPE if stdin_input else None,
                text=True,
                bufsize=1,  # Line buffered
            )

            # If we have stdin input, write it and close
            if stdin_input:
                process.stdin.write(stdin_input)
                process.stdin.close()

            # Stream output line by line
            for line in process.stdout:
                sys.stdout.write(line)
                sys.stdout.flush()

            process.wait()

            if process.returncode != 0:
                raise GromacsError(
                    command=" ".join(cmd),
                    returncode=process.returncode,
                    message=f"Step failed: {description}",
                )

            logger.info("-" * 60)
            logger.info(f">>> {description} completed successfully")

        except FileNotFoundError:
            raise GromacsError(
                command=cmd[0],
                returncode=-1,
                message=f"GROMACS executable not found: {self._gmx}. "
                "Ensure GROMACS is installed and in your PATH, or use --gmx-path.",
            )

    def _run_energy_minimization(self) -> None:
        """Run energy minimization step."""
        self._print_banner("Step 1: Energy Minimization")

        # grompp
        grompp_cmd = [
            self._gmx,
            "grompp",
            "-f",
            "em.mdp",
            "-c",
            f"{self._prefix}.gro",
            "-r",
            f"{self._prefix}.gro",
            "-p",
            f"{self._prefix}.top",
            "-o",
            "em.tpr",
            "-maxwarn",
            "1",
        ]

        self._run_command(
            grompp_cmd,
            "Preparing energy minimization (grompp)",
        )

        # mdrun
        self._run_command(
            [self._gmx, "mdrun", "-deffnm", "em", "-v"],
            "Running energy minimization (mdrun)",
        )

        # Verify output
        if not (self._working_dir / "em.gro").exists():
            raise GromacsError(
                command="mdrun -deffnm em",
                returncode=1,
                message="Energy minimization failed: em.gro not produced",
            )

        self._check_energy_minimization_health()

    def _check_energy_minimization_health(self) -> None:
        """Check EM log for non-finite force failures.

        Raises
        ------
        RuntimeError
            If the EM log reports non-finite forces.
        """
        em_log_path = self._working_dir / "em.log"
        if not em_log_path.exists():
            return

        em_log_text = em_log_path.read_text(errors="ignore")
        if re.search(r"force.*not finite|inf.*atom", em_log_text, flags=re.IGNORECASE):
            raise RuntimeError(
                "Energy minimization failed — infinite forces detected in em.log. "
                "This usually indicates severe atomic overlaps in the initial structure. "
                "Try increasing packing padding or box size."
            )

    def _run_equilibration(self) -> None:
        """Run all equilibration stages."""
        self._print_banner("Step 2: Equilibration", f"{len(self._eq_mdps)} stage(s)")

        prev_output = "em"

        for i, mdp_file in enumerate(self._eq_mdps):
            stage_num = i + 1
            output_name = f"eq_{stage_num:02d}"

            logger.info(
                f"\n=== Equilibration Stage {stage_num}/{len(self._eq_mdps)}: {mdp_file} ==="
            )

            # Build grompp command
            # Always use -r em.gro for reference coordinates (needed for position restraints)
            grompp_cmd = [
                self._gmx,
                "grompp",
                "-f",
                mdp_file,
                "-c",
                f"{prev_output}.gro",
                "-r",
                "em.gro",
                "-p",
                f"{self._prefix}.top",
                "-o",
                f"{output_name}.tpr",
                "-maxwarn",
                "1",
            ]

            # Add checkpoint from previous stage (except after EM)
            if prev_output != "em":
                grompp_cmd.extend(["-t", f"{prev_output}.cpt"])

            self._run_command(
                grompp_cmd,
                f"Preparing equilibration stage {stage_num} (grompp)",
            )

            # mdrun
            self._run_command(
                [self._gmx, "mdrun", "-deffnm", output_name, "-v"],
                f"Running equilibration stage {stage_num} (mdrun)",
            )

            # Verify output
            if not (self._working_dir / f"{output_name}.gro").exists():
                raise GromacsError(
                    command=f"mdrun -deffnm {output_name}",
                    returncode=1,
                    message=f"Equilibration stage {stage_num} failed: {output_name}.gro not produced",
                )

            prev_output = output_name

        self._last_eq_output = prev_output

    def _run_production(self) -> None:
        """Run production MD."""
        self._print_banner("Step 3: Production MD")

        if self._last_eq_output is None:
            raise RuntimeError("Equilibration must complete before production")

        # grompp
        self._run_command(
            [
                self._gmx,
                "grompp",
                "-f",
                "prod.mdp",
                "-c",
                f"{self._last_eq_output}.gro",
                "-t",
                f"{self._last_eq_output}.cpt",
                "-p",
                f"{self._prefix}.top",
                "-o",
                "prod.tpr",
                "-maxwarn",
                "1",
            ],
            "Preparing production MD (grompp)",
        )

        # mdrun
        self._run_command(
            [self._gmx, "mdrun", "-deffnm", "prod", "-v"],
            "Running production MD (mdrun)",
        )

        # Verify output
        if not (self._working_dir / "prod.xtc").exists():
            raise GromacsError(
                command="mdrun -deffnm prod",
                returncode=1,
                message="Production MD failed: prod.xtc not produced",
            )

    def _run_postprocessing(self) -> None:
        """Run trajectory post-processing with trjconv."""
        self._print_banner("Step 4: Trajectory Post-processing")

        logger.info("Removing periodic boundary jumps and centering trajectory...")

        # Step 1: Remove PBC jumps
        # echo "System" | gmx trjconv -s prod.tpr -f prod.xtc -o prod_nojump.xtc -pbc nojump
        self._run_command(
            [
                self._gmx,
                "trjconv",
                "-s",
                "prod.tpr",
                "-f",
                "prod.xtc",
                "-o",
                "prod_nojump.xtc",
                "-pbc",
                "nojump",
            ],
            "Removing PBC jumps (trjconv)",
            stdin_input="System\n",
        )

        # Step 2: Center on protein and make molecules whole
        # echo -e "Protein\nSystem" | gmx trjconv -s prod.tpr -f prod_nojump.xtc -o prod_centered.xtc -center -pbc mol -ur compact
        self._run_command(
            [
                self._gmx,
                "trjconv",
                "-s",
                "prod.tpr",
                "-f",
                "prod_nojump.xtc",
                "-o",
                "prod_centered.xtc",
                "-center",
                "-pbc",
                "mol",
                "-ur",
                "compact",
            ],
            "Centering trajectory (trjconv)",
            stdin_input="Protein\nSystem\n",
        )

    def _print_banner(self, title: str, subtitle: str = "") -> None:
        """Print a section banner."""
        logger.info("\n" + "=" * 60)
        logger.info(title)
        if subtitle:
            logger.info(subtitle)
        logger.info("=" * 60)

    def _print_summary(self) -> None:
        """Print workflow completion summary."""
        self._print_banner("Workflow Complete!", POLYZYMD_BRANDING)
        logger.info("\nOutput files:")
        logger.info("  em.gro              - Minimized structure")
        for i in range(len(self._eq_mdps)):
            logger.info(f"  eq_{i + 1:02d}.gro           - Equilibrated structure (stage {i + 1})")
        logger.info("  prod.xtc            - Production trajectory")
        logger.info("  prod.edr            - Energy file")
        logger.info("  prod.log            - Log file")
        logger.info("  prod_nojump.xtc     - Trajectory without PBC jumps")
        logger.info("  prod_centered.xtc   - Centered trajectory for visualization")
        logger.info(f"\nAll files in: {self._working_dir}")
