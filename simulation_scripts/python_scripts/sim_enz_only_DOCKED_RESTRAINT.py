import openmm
from openff.interchange.components import _packmol as packmol
from openmm.unit import Quantity, nanometer, angstrom
from polymerist.mdtools.openfftools import boxvectors

from openmm.unit import gram, centimeter
from polymerist.mdtools.openfftools.solvation import packing, solvents

from openff.toolkit import Molecule, Topology, ForceField

from pathlib import Path
from polymerist.genutils.fileutils.pathutils import assemble_path
from polymerist.mdtools.openfftools.partition import partition
from polymerist.mdtools.openfftools.topology import topology_to_sdf, topology_from_sdf, get_largest_offmol

import argparse

parser = argparse.ArgumentParser(add_help=True)

parser.add_argument("-p", "--protein_name", type=str, required=True, help="Name of the protein (str)")
parser.add_argument("-t", "--temperature", type=float, required=True, help="Temperature in Kelvin (float)")
parser.add_argument("-e", "--eq_time", type=float, required=True, help="Equilibration time in nanoseconds (float)")
parser.add_argument("-r", "--prod_time", type=float, required=True, help="Production time in nanoseconds (float)")
parser.add_argument("-pp", "--path_to_pro_pdb", type=str, required=True, help="Path to the protein PDB file (str)")
parser.add_argument("-rep", "--replicate", type=int, required=True, help="Replicate number (int)")

# comment below arguments out if not doing daisy chain
parser.add_argument("--frames", type=int, default=2500,
                    help="""Number of frames to save for this simulation. NOTE: since this \\
                          workflow is configured to send this off once, the number in \\
                          this function should be 1/N what you want for the entire simulation. \\
                          I.E, you want 2500 frames and 10 segments, so this is 250.""")
parser.add_argument("--num_segs", type=int, default=1, help="If breaking simulation into daisy chains, number of breaks.")

args = parser.parse_args()

protein_name = args.protein_name #"LipA"
ligand = "Resorufin-Butyrate"

# FIX THIS WHEN NOT LAZY!!!!
if protein_name == "LipA":
    file_containing_docked_structures = "docked.sdf"
elif protein_name == "LipA-dBur-100":
    file_containing_docked_structures = "docked_RBY_for_LipA_dBur.sdf"
elif protein_name == "LipA-dSA-100":
    file_containing_docked_structures = "docked_RBY_for_LipA_dSA-100.sdf"
elif protein_name == "LipA-dSA-50":
    file_containing_docked_structures = "docked_RBY_for_LipA_dSA-50.sdf"

temperature =  args.temperature #293 # Kelvin
eq_time = args.eq_time # 0.5 # nanoseconds
prod_time = args.prod_time #200 # nanoseconds
replicate_num = args.replicate

#comment out if no daisy chain
num_prod_frames = args.frames

# for iterating over docked substrate idxs
best_substrate_idx = 6

# Setup working directory
working_dir_tag = f"3,3A_RESTRAINT_oxyRestraint_{protein_name}_{ligand}_conf{best_substrate_idx+1}_{temperature}K_{eq_time}ns-NVT_{prod_time*args.num_segs}ns-NPT_run{replicate_num}"
WORKING_DIR: Path = Path(working_dir_tag)
WORKING_DIR.mkdir(exist_ok=True)


# working_dir_tag = f"10A_RESTRAINT_{protein_name}_{ligand}_{temperature}K_{eq_time}ns-NVT_{prod_time}ns-NPT_run{replicate_num}"
# WORKING_DIR : Path = Path(working_dir_tag)
# WORKING_DIR.mkdir(exist_ok=True)

path_to_pro_pdb = args.path_to_pro_pdb #"NH3_terminal_updated.pdb"

off_pro_top = Topology.from_pdb(path_to_pro_pdb)
was_partitioned = partition(off_pro_top)
assert was_partitioned

# define periodic box vectors, based on a padded version of the tight bounding box of all your molecules (the topology)
box_padding : Quantity = 1.2 * nanometer # THIS IS PADDING BETWEEN PROTEIN AND POLYMERS
box_vecs_tight = boxvectors.get_topology_bbox(off_pro_top)
box_vecs = boxvectors.pad_box_vectors_uniform(box_vecs_tight, box_padding)
target_density : Quantity = 1.0 * gram * centimeter**-3

from polymerist.mdtools.openfftools.solvation import physprops
from typing import Optional, Union
import logging
import numpy as np

LOGGER = logging.getLogger(__name__)

from polymerist.mdtools.openfftools import topology
from openff.toolkit import ForceField
from polymerist.mdtools.openmmtools.forcegroups import impose_unique_force_groups

#################################################################

# Prepare substrate structure, for generating charges

from polymerist.mdtools.openfftools.partialcharge.molchargers import NAGLCharger, EspalomaCharger

charger = NAGLCharger()

off_lig_mol = Molecule.from_file(file_containing_docked_structures)
off_lig_mol_best = charger.charge_molecule(off_lig_mol[best_substrate_idx])  # I looked at the structures in PyMol and
                                                            # selected the one I thought most amenable to Sn2 chemistry
                                                            # did 8 for first run, now considering 1
                                                            # conf 1 also looked pretty good with the 3.3 A restraint
                                                            # conf 1 looks the best IMO

##############################################################################
# Combine both ligand and protein into one Topology
combined_top = Topology.from_molecules([off_pro_top.molecule(0), off_lig_mol_best])

##################################################################
# Re-number the Chains

total_mols = combined_top.n_molecules
chain_idxs = range(1, total_mols+1)

i = 0
for mol in combined_top.molecules:
    for atom in mol.atoms:
        atom.metadata['chain_id'] = f"{i+1}"
    i+=1

##################################################################
# Solvate the topology and add ions
from openff.interchange.components._packmol import *
import numpy as np
import MDAnalysis as mda
from polymerist.mdtools.openfftools.unitsys import openmm_to_openff

def solvate_topology_with_ions(
    topology: Topology,
    nacl_conc: Quantity = Quantity(0.1, "mole / liter"),
    padding: Quantity = Quantity(1.2, "nanometer"),
    box_shape: NDArray = RHOMBIC_DODECAHEDRON,
    target_density: Quantity = Quantity(0.9, "gram / milliliter"),
    tolerance: Quantity = Quantity(2.0, "angstrom"),

) -> Topology:

    """
    Add water and ions to neutralise and solvate a topology.
    Parameters
    ----------
    topology
        The OpenFF Topology to solvate.
    nacl_conc
        The bulk concentration of NaCl in the solvent, in units compatible with
        molarity. This is used to calculate a mass fraction for the bulk
        solvent and does not represent the actual concentration in the final
        box.

    padding : Scalar with dimensions of length
        The desired distance between the solute and the edge of the box. Ignored
        if the topology already has box vectors. The usual recommendation is
        that this equals or exceeds the VdW cut-off distance, so that the
        solute is isolated by its periodic images by twice the cut-off.
    box_shape : Array with shape (3, 3)
        An array defining the box vectors of a box with the desired shape and
        unit periodic image distance. This shape will be scaled to satisfy the
        padding given above. Some typical shapes are provided in this module.
    target_density : Scalar with dimensions of mass density

        The target mass density for the packed box.
    tolerance: Scalar with dimensions of distance

        The minimum spacing between molecules during packing in units of
        distance. The default is large so that added waters do not disrupt the
        structure of proteins; when constructing a mixture of small molecules,
        values as small as 0.5 Å will converge faster and can still produce
        stable simulations after energy minimisation.

    Returns
    -------
    Topology
        An OpenFF ``Topology`` with the solvated system.
    Raises
    ------
    PACKMOLRuntimeError
        When packmol fails to execute / converge.
    Notes
    -----
    Returned topologies may have larger box vectors than what would be defined
    by the target density.
    """

    #_check_box_shape_shape(box_shape)
    #box_vecs = boxvectors.box_vectors_flexible(box_vecs) # up-convert to full box vectors if only dimensions are given
    min_box_vecs = boxvectors.get_topology_bbox(topology)
    box_vecs = boxvectors.pad_box_vectors_uniform(min_box_vecs, padding)

    if not np.all(box_vecs >= min_box_vecs): # TODO : change this to only check XYZ dimensions
        raise boxvectors.BoxVectorError(f'Desired box dimensions ({box_vecs}) are smaller than minimum excluded Topology dimensions ({min_box_vecs})')

    box_vecs = box_shape @ box_vecs
    print(box_vecs)

    box_center = np.sum(box_vecs, axis=0) / 2

    # Translate all molecules in the topology so they are centered at the box_center
    # NOTE SUPER HACKY, MUST CHANGE FOR A DIFFERENT ENZYME, BUT FOR NOW
    # THIS WORKS AND WILL RECENTER LipA AND THE LIGAND AT CENTER OF BOX
    com_protein = openmm_to_openff(mda.Universe(path_to_pro_pdb).select_atoms("protein").center_of_mass()*angstrom)
    box_center = (np.sum(box_vecs, axis=0) / 2)
    trans_vec = box_center - com_protein
    old_positions = topology.get_positions()
    new_positions = old_positions + trans_vec
    topology.set_positions(new_positions)

    # assign properties from desired box vectors if size checks pass
    box_vol = boxvectors.get_box_volume(box_vecs, units_as_openm=False)

    print(f"Computed new solvent box vectors of {box_vecs}")

    # Compute target masses of solvent

    #box_volume = numpy.linalg.det(box_vectors.m) * box_vectors.u**3
    box_volume = box_vol
    target_mass = box_volume * target_density
    solute_mass = sum(sum([atom.mass for atom in molecule.atoms]) for molecule in topology.molecules)
    solvent_mass = target_mass - solute_mass

    #_check_add_positive_mass(solvent_mass)

    # Get reference data and prepare solvent molecules
    water = Molecule.from_smiles("O")
    na = Molecule.from_smiles("[Na+]")
    cl = Molecule.from_smiles("[Cl-]")
    nacl_mass = sum([atom.mass for atom in na.atoms]) + sum(
        [atom.mass for atom in cl.atoms],

    )

    water_mass = sum([atom.mass for atom in water.atoms])
    molarity_pure_water = Quantity(55.5, "mole / liter")

    # Compute the number of salt "molecules" to add from the mass and concentration
    nacl_mass_fraction = (nacl_conc * nacl_mass) / (molarity_pure_water * water_mass)
    nacl_mass_to_add = solvent_mass * nacl_mass_fraction
    nacl_to_add = (nacl_mass_to_add / nacl_mass).m_as("dimensionless").round()

    # Compute the number of water molecules to add to make up the remaining mass
    water_mass_to_add = solvent_mass - nacl_mass
    water_to_add = (water_mass_to_add / water_mass).m_as("dimensionless").round()

    # Neutralise the system by adding and removing salt
    solute_charge = sum([molecule.total_charge for molecule in topology.molecules])
    na_to_add = np.ceil(nacl_to_add - solute_charge.m / 2.0)
    cl_to_add = np.floor(nacl_to_add + solute_charge.m / 2.0)

    print(f"Computed new solvent box vectors of {box_vecs}")

    # Pack the box
    return pack_box(
        [water, na, cl],
        [int(water_to_add), int(na_to_add), int(cl_to_add)],
        solute=topology,
        tolerance=tolerance,
        box_vectors=box_vecs,
    )


#################################################################################################
import numpy as np
from polymerist.mdtools.openfftools.unitsys import openmm_to_openff
from openff.interchange.components import _packmol as packmol

solvated_combined_top = solvate_topology_with_ions(topology=combined_top,
                                                 padding = openmm_to_openff(0.9 * nanometer),
                                                 nacl_conc = Quantity(0.1, "mole / liter"),
                                                 box_shape=packmol.RHOMBIC_DODECAHEDRON)#,
                                                 #box_shape=packmol.UNIT_CUBE)

for mol in solvated_combined_top.molecules:
    if mol.to_smiles() == "[H][O][H]":
        for atom in mol.atoms:
            atom.metadata['residue_name'] = "HOH"

    elif mol.to_smiles() == "[Na+]":
        for atom in mol.atoms:
            atom.metadata['residue_name'] = "Na+"

    elif mol.to_smiles() == "[Cl-]":
        for atom in mol.atoms:
            atom.metadata['residue_name'] = "Cl-"

    elif mol.to_smiles() == "[H][C]1=[C]([H])[C]2=[N][c]3[c]([H])[c]([H])[c]([O][C](=[O])[C]([H])([H])[C]([H])([H])[C]([H])([H])[H])[c]([H])[c]3[O][C]2=[C]([H])[C]1=[O]":
        for atom in mol.atoms:
            atom.metadata['residue_name'] = "RBY"

solvated_combined_top.to_file(assemble_path(WORKING_DIR, working_dir_tag, postfix='solvated', extension='pdb')) # saving extra PDB purely for visualization (SDF only lets you see molecules one-at-a-time)


##################################################################
# create interchange object
ff = ForceField("ff14sb_off_impropers_0.0.4.offxml", 'openff-2.0.0.offxml') #, )
inc = ff.create_interchange(solvated_combined_top, charge_from_molecules=[off_lig_mol_best, solvents.water_TIP3P])

# define how you want to run your simulation

from openmm.unit import kelvin, nanoseconds, picoseconds, femtoseconds
from openmm.unit import femtosecond, picosecond, nanosecond, kelvin, atmosphere
from polymerist.mdtools.openmmtools.parameters import (
    ThermoParameters,
    ThermostatParameters,
    BarostatParameters,
    IntegratorParameters,
    ReporterParameters,
    SimulationParameters
)
from polymerist.mdtools.openmmtools.reporters import DEFAULT_STATE_DATA_PROPS

# Legacy way of describing simulation parameters
# sim_params_equil = SimulationParameters(
#     thermo_params=ThermoParameters(
#         ensemble='NVT',
#         temperature=temperature*kelvin,
#         friction_coeff=1*picoseconds**-1,
#     ),
#     integ_params=IntegratorParameters(
#         time_step=2*femtoseconds,
#         total_time=eq_time*nanoseconds,
#         num_samples=10,
#     ),
#     reporter_params=ReporterParameters(
#         traj_ext='dcd',
#         state_data=DEFAULT_STATE_DATA_PROPS, # NOTE: you can tweak which of these you want to write out during a simulation
#     )
# )

# New Spec compliant with Polymerist 1.0
sim_params_equil = SimulationParameters(
    thermo_params=ThermoParameters(
        thermostat_params=ThermostatParameters(
            temperature=temperature*kelvin,
            timescale=1*picoseconds**-1,
            thermostat='LangevinMiddle'
        ),
    ),
    integ_params=IntegratorParameters(
        time_step=2*femtosecond,
        total_time=eq_time*nanoseconds,
        num_samples=10,
    ),
    reporter_params=ReporterParameters(
        state_data=DEFAULT_STATE_DATA_PROPS,
        traj_ext='dcd',
    ),
)

# Legacy
# sim_params_prod = SimulationParameters(
#     thermo_params=ThermoParameters(
#         ensemble='NPT',
#         temperature=temperature*kelvin,
#         friction_coeff=1*picoseconds**-1,
#     ),
#     integ_params=IntegratorParameters(
#         time_step=2*femtoseconds,
#         total_time=prod_time*nanoseconds,
#         num_samples=num_prod_frames, #2500,
#     ),
#     reporter_params=ReporterParameters(
#         traj_ext='dcd',
#         state_data=DEFAULT_STATE_DATA_PROPS, # NOTE: you can tweak which of these you want to write out during a simulation
#     )
# )
# sim_params.to_file(WORKING_DIR / 'equil.json')

sim_params_prod = SimulationParameters(
    thermo_params=ThermoParameters(
        thermostat_params=ThermostatParameters(
            temperature=temperature*kelvin,
            timescale=1*picoseconds**-1,
            thermostat='LangevinMiddle'
        ),
        barostat_params=BarostatParameters(
            pressure=1*atmosphere,
            update_frequency = 25,
            barostat='MC', # Monte-Carlo barostat
        )
    ),
    integ_params=IntegratorParameters(
        time_step=2*femtosecond,
        total_time=prod_time*nanoseconds,
        num_samples=num_prod_frames,
    ),
    reporter_params=ReporterParameters(
        traj_ext='dcd',
        state_data=DEFAULT_STATE_DATA_PROPS
    ),
)

# extract OpenMM stuff from the Interchange, ensuring you can decompose your FF energy contributions later
from openff.interchange.interop.openmm._positions import to_openmm_positions

omm_topology = inc.to_openmm_topology()
omm_system  = inc.to_openmm(combine_nonbonded_forces=False)
omm_positions = to_openmm_positions(inc, include_virtual_sites=True)

##########################################################################################
from openmm import CustomBondForce
from openmm.unit import nanometer, angstrom, kilojoule_per_mole

ser_O_group = []
RBY_grp = []

# Old selection way, works for LipA but no guarantees
# for other proteins

# for atom in omm_topology.atoms():
#     if atom.name == "OG" and atom.residue.index == 76:
#         print(atom)
#         ser_O_group.append(atom.index)

#     if atom.residue.name == "RBY" and atom.index == 2739:
#         print(atom)
#         RBY_grp.append(atom.index)


MET_78_group = []
ILE_12_group = []

for atom in omm_topology.atoms():
    if atom.name == "OG" and atom.residue.index == 76:
        print(atom)
        ser_O_group.append(atom.index)

    if atom.residue.name == "RBY" and atom.index == 2739:
        print(atom)
        RBY_grp.append(atom.index)

    if  atom.residue.index == 77:
        if atom.index == 1192: # in the pdb the atom is 1193 - 1 == 1192
            MET_78_group.append(atom.index)
            print(atom)

    if atom.residue.index == 11:
        if atom.index == 171: # in the pdb atom is 172 - 1 == 171
            print(atom)
            ILE_12_group.append(atom.index)


"""
Adds a flat-bottom potential between two selected atoms:
U(r) = 0 if r < r0
        0.5 * k * (r - r0)^2 if r >= r0
"""

restraint_dist = 3.3 # angstroms

k = 10000 * kilojoule_per_mole / nanometer**2
r0 = restraint_dist * angstrom

# OpenMM expression: r is the distance between atoms
expression = "step(r - r0) * 0.5 * k * (r - r0)^2"
fb_force = CustomBondForce(expression)
fb_force.addGlobalParameter("k", k)
fb_force.addGlobalParameter("r0", r0)

atom1_index = ser_O_group[0]
atom2_index = RBY_grp[0]

# Step 4: Add the bond between selected atoms
fb_force.addBond(atom1_index, atom2_index, [])

# Step 5: Add the force to the system
omm_system.addForce(fb_force)

print(f"Flat-bottom restraint added between atom {atom1_index} and {atom2_index} with r0={r0}Å, k={k}")

#########################################################################################

"""
Adds a flat-bottom potential between two selected atoms:
U(r) = 0 if r < r1
        0.5 * k * (r - r1)^2 if r >= r1
"""

restraint_dist = 6.0 # angstroms

k = 10000 * kilojoule_per_mole / nanometer**2
r1 = restraint_dist * angstrom

# OpenMM expression: r is the distance between atoms
expression = "step(r - r1) * 0.5 * k * (r - r1)^2"
fb_force_oxy = CustomBondForce(expression)
fb_force_oxy.addGlobalParameter("k", k)
fb_force_oxy.addGlobalParameter("r1", r1)

atom1_index = MET_78_group[0]
atom2_index = RBY_grp[0]

# Step 4: Add the bond between selected atoms
fb_force_oxy.addBond(atom1_index, atom2_index, [])

# Step 5: Add the force to the system
omm_system.addForce(fb_force_oxy)

print(f"Flat-bottom restraint added between atom {atom1_index} and {atom2_index} with r1={r1}Å, k={k}")

atom1_index = ILE_12_group[0]
atom2_index = RBY_grp[0]

# Step 4: Add the bond between selected atoms
fb_force_oxy.addBond(atom1_index, atom2_index, [])

# Step 5: Add the force to the system
#omm_system.addForce(fb_force_oxy)

print(f"Flat-bottom restraint added between atom {atom1_index} and {atom2_index} with r1={r1}Å, k={k}")



#########################################################################################

impose_unique_force_groups(omm_system) # ensure each Force is separate to enable mapping of energy contributions

from polymerist.mdtools.openmmtools.execution import run_simulation_schedule
from polymerist.genutils.logutils.IOHandlers import MSFHandlerFlex, get_active_loggers

logpath = assemble_path(WORKING_DIR, 'simulation_status', extension='log')
with MSFHandlerFlex(filename=logpath, loggers=get_active_loggers(), proc_name='AF3_LipA_sims') as logger:
    history = run_simulation_schedule(
        working_dir=WORKING_DIR,
        schedule={
            'equilibration' : sim_params_equil,
            'production_0' : sim_params_prod
        },
        init_top=omm_topology,
        init_sys=omm_system,
        init_pos=omm_positions,
        return_history=True
    )


