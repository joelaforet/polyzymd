import argparse
from pathlib import Path
from openff.toolkit import Topology
from polymerist.mdtools.openfftools.topology import get_largest_offmol
from polymerist.mdtools.openfftools.topology import topology_from_sdf
from polymerist.mdtools.openfftools.solvation import physprops
from polymerist.mdtools.openfftools import boxvectors
from polymerist.mdtools.openfftools.partition import partition
from typing import Optional, Union
import numpy as np
import logging
from openff.units import Quantity
import openmm
import openmm.unit
from openff.toolkit import Molecule

from openmm.unit import nanometer, gram, centimeter

from openff.interchange.components import _packmol as packmol

from polymerist.genutils.fileutils.pathutils import assemble_path
from polymerist.mdtools.openfftools.solvation import packing, solvents

LOGGER = logging.getLogger(__name__)

parser = argparse.ArgumentParser(add_help=True)

# Original arguments
parser.add_argument("-p", "--protein_name", type=str, required=True, help="Name of the protein (str)")
parser.add_argument("-o", "--oligomer_name", type=str, required=True, help="Name of the oligomer (str)")
parser.add_argument("-n", "--num_oligomers", type=int, required=True, help="Number of oligomers (int)")
parser.add_argument("-t", "--temperature", type=float, required=True, help="Temperature in Kelvin (float)")
parser.add_argument("-e", "--eq_time", type=float, required=True, help="Equilibration time in nanoseconds (float)")
parser.add_argument("-r", "--prod_time", type=float, required=True, help="Production time in nanoseconds (float)")
parser.add_argument("-pp", "--path_to_pro_pdb", type=str, required=True, help="Path to the protein PDB file (str)")
parser.add_argument("-ps", "--path_to_poly_sdf", type=str, required=True, help="Path to the *directory* containing polymer SDFs (str)")
parser.add_argument("-rep", "--replicate", type=int, required=True, help="Replicate number (int)")
parser.add_argument("--num_segs", type=int, default=1, help="If breaking simulation into daisy chains, number of breaks.")

# New arguments for co-polymer randomization
parser.add_argument("-chars", "--characters", type=str, nargs="+", required=True,
                    help="List of monomer unit labels (e.g., A B C)")
parser.add_argument("-probs", "--probabilities", type=float, nargs="+", required=True,
                    help="List of monomer probabilities (must sum to 1.0)")
parser.add_argument("-plen", "--polymer_length", type=int, required=True,
                    help="Length of each polymer (number of monomer units)")
parser.add_argument("-prefix", "--polymer_type_prefix", type=str, required=True,
                    help="Prefix of polymer type in filenames (e.g., SBMA-EGPMA)")
parser.add_argument("--frames", type=int, default=2500,
                    help="""Number of frames to save for this simulation. NOTE: since this \\
                          workflow is configured to send this off once, the number in \\
                          this function should be 1/N what you want for the entire simulation. \\
                          I.E, you want 2500 frames and 10 segments, so this is 250.""")
args = parser.parse_args()

# Assign input arguments
protein_name = args.protein_name
oligomer_name = args.oligomer_name
num_oligomers = args.num_oligomers
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

temperature = args.temperature
eq_time = args.eq_time
prod_time = args.prod_time
replicate_num = args.replicate
path_to_pro_pdb = args.path_to_pro_pdb
path_to_poly_sdf_dir = args.path_to_poly_sdf  # NOTE: now a directory, not a file

# New variables
characters = args.characters
probabilities = args.probabilities
polymer_length = args.polymer_length
polymer_type_prefix = args.polymer_type_prefix

num_prod_frames = args.frames

# Setup working directory
working_dir_tag = f"3,3A_RESTRAINT_{protein_name}_{ligand}_{oligomer_name}-{int(probabilities[1]*100)}%_{num_oligomers}x_{temperature}K_{eq_time}ns-NVT_{prod_time*args.num_segs}ns-NPT_run{replicate_num}"
WORKING_DIR: Path = Path(working_dir_tag)
WORKING_DIR.mkdir(exist_ok=True)

# Load protein topology
off_pro_top = Topology.from_pdb(path_to_pro_pdb)
was_partitioned = partition(off_pro_top)
assert was_partitioned

# Define box
# Changed padding from 1.2 to 1.5nm 5/16/25
box_padding: Quantity = 1.5 * nanometer  # Padding between protein and polymers
box_vecs_tight = boxvectors.get_topology_bbox(off_pro_top)
box_vecs = boxvectors.pad_box_vectors_uniform(box_vecs_tight, box_padding)

# Target density (if you need it later)
target_density: Quantity = 1.0 * gram * centimeter**-3

import random
from collections import Counter
import os
from typing import List, Union, Optional
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
# off_lig_mol_best = charger.charge_molecule(off_lig_mol[8]) # I looked at the structures in PyMol and 
#                                                            # selected the one I thought most amenable to Sn2 chemistry

# FIX THIS WHEN NOT LAZY!!!!
if protein_name == "LipA":
    #best_mol_idx = 8 #previously used 8, switching to 1
    best_mol_idx = 1
elif protein_name == "LipA-dBur-100":
    best_mol_idx = 0
elif protein_name == "LipA-dSA-100":
    best_mol_idx = 0
elif protein_name == "LipA-dSA-50":
    best_mol_idx = 0

off_lig_mol_best = charger.charge_molecule(off_lig_mol[best_mol_idx]) # I looked at the structures in PyMol and 
                                                           # selected the one I thought most amenable to Sn2 chemistry


##############################################################################

##############################################################################

# Combine both ligand and protein into one Topology

combined_top = Topology.from_molecules([off_pro_top.molecule(0), off_lig_mol_best])

for mol in combined_top.molecules:
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

##################################################################


# TODO ADD IN CODE THAT REJECTS DUPLICATE SEQUENCES
# I.E. SEQUENCE 'ABAAA' VS 'AAABA' SHOULD MAP TO THE SAME MOLECULE!
# EDIT fixed this :)

def pack_topology_with_random_polymers_from_sequences(offtop,
                                                       characters: List[str],
                                                       probabilities: List[float],
                                                       box_vecs,
                                                       N: int,
                                                       molecule_directory: str,
                                                       exclusion: Optional = None,
                                                       working_directory=None,
                                                       retain_working_files=True,
                                                       polymer_length: int = 5,
                                                       polymer_type_prefix: str = "SBMA-EGPMA"):
    '''
    Pack a Topology with N randomly generated polymer sequences, selecting matching .sdf files.

    Parameters:
    - offtop: Topology to solvate.
    - characters: Possible monomer units ("A", "B", etc).
    - probabilities: Probabilities associated with each monomer unit.
    - box_vecs: Desired box vectors.
    - N: Total number of polymer chains to insert.
    - molecule_directory: Directory where .sdf files are stored.
    - polymer_length: Number of monomers per chain (default 5).
    - polymer_type_prefix: Prefix for molecule filenames (e.g., "SBMA-EGPMA").
    '''
    

    assert abs(sum(probabilities) - 1.0) < 1e-8, "Probabilities must sum to 1."
    assert len(characters) == len(probabilities), "Characters and probabilities must match."

    box_vecs = boxvectors.box_vectors_flexible(box_vecs)
    min_box_vecs = boxvectors.get_topology_bbox(offtop) * 1.1
    if exclusion is not None:
        min_box_vecs = boxvectors.pad_box_vectors_uniform(min_box_vecs, exclusion)

    if not np.all(box_vecs >= min_box_vecs):
        raise boxvectors.BoxVectorError(f'Desired box dimensions ({box_vecs}) are smaller than minimum excluded Topology dimensions ({min_box_vecs})')

    box_vol = boxvectors.get_box_volume(box_vecs, units_as_openm=True)

    # Step 1: Randomly generate N polymer sequences
    def generate_random_polymer(N, characters, probabilities):
        return ''.join(random.choices(characters, weights=probabilities, k=N))

    polymer_sequences = [generate_random_polymer(polymer_length, characters, probabilities) for _ in range(N)]

    # Define a canonical form: min of the string and its reverse
    def canonical(s):
        return min(s, s[::-1])

    # Convert all strings to their canonical form
    canonical_polys = [canonical(p) for p in polymer_sequences]

    # Step 2: Count how many times each unique polymer was generated
    sequence_counts = Counter(canonical_polys)
    print(sequence_counts)

    # Step 3: Load molecule objects matching each unique sequence
    loaded_mols = {}
    for sequence in sequence_counts:
        filename = f"{polymer_type_prefix}_{sequence}_{polymer_length}-mer_charged.sdf"
        filepath = os.path.join(molecule_directory, filename)
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Could not find file for sequence {sequence}: {filepath}")
        off_poly_top = topology_from_sdf(filepath)
        off_poly_mol = get_largest_offmol(off_poly_top)
        loaded_mols[sequence] = off_poly_mol

    # Step 4: Prepare lists for packmol
    solvent_list = []
    num_list = []
    seq_list = []

    for sequence, count in sequence_counts.items():
        solvent_list.append(loaded_mols[sequence])
        num_list.append(count)
        seq_list.append(sequence)

    # Step 5: Pack the box
    LOGGER.info(f"Packing {sum(num_list)} total polymers: {sequence_counts}")
    packed_top = packmol.pack_box(
        solvent_list,
        num_list,
        offtop,
        box_vectors=box_vecs,
        box_shape=packmol.UNIT_CUBE,
        center_solute='BRICK',
        working_directory=working_directory,
        retain_working_files=retain_working_files
    )
    LOGGER.info('Packmol packing converged')

    packed_top.box_vectors = box_vecs
    print("Sequences generated:")
    print(seq_list)
    return packed_top, solvent_list


from polymerist.mdtools.openfftools import topology
from openff.toolkit import ForceField
from polymerist.mdtools.openmmtools.forcegroups import impose_unique_force_groups

print("Adding oligomers to box!")

poly_pro_top, polymer_mol_list = pack_topology_with_random_polymers_from_sequences(
    offtop=combined_top,
    characters=characters,
    probabilities=probabilities,
    box_vecs=box_vecs,
    N=num_oligomers,
    molecule_directory=path_to_poly_sdf_dir,
    polymer_length=polymer_length,
    polymer_type_prefix=oligomer_name,  # or whatever your base prefix is
    working_directory=str(WORKING_DIR),
    retain_working_files=True
)

print("Finished adding oligomers!")

##################################################################

# Re-number the Chains

total_mols = poly_pro_top.n_molecules
chain_idxs = range(1, total_mols+1)

i = 0
for mol in poly_pro_top.molecules:
    for atom in mol.atoms:
        atom.metadata['chain_id'] = f"{i+1}"
    i+=1

##################################################################

##################################################################

# Utility function for computing center of mass

from openff.toolkit import Quantity, unit

def calculate_center_of_mass(topology):
    # Dictionary of atomic masses for common atoms used in organic chemistry
    common_atoms = {
        1: 1.008,   # Hydrogen
        6: 12.011,  # Carbon
        7: 14.007,  # Nitrogen
        8: 15.999,  # Oxygen
        9: 18.998,  # Fluorine
        15: 30.974, # Phosphorus
        16: 32.06,  # Sulfur
        17: 35.45,  # Chlorine
    }

    all_conf_coords = []
    all_masses = []

    for molecule in topology.molecules:
        conf_coords = molecule.conformers[0].m_as(unit.angstrom)
        all_conf_coords.append(conf_coords)

        for atom in molecule.atoms:
            atomic_number = atom.atomic_number
            atomic_mass = common_atoms.get(atomic_number, 0) # Get atomic mass from dictionary, default to 0 if not found
            all_masses.append(atomic_mass)

    all_positions = np.concatenate(all_conf_coords, axis=0)
    all_masses = np.array(all_masses)

    # Calculate the center of mass
    center_of_mass = np.average(all_positions, axis=0, weights=all_masses)*unit.angstrom

    return center_of_mass


##################################################################


# Solvate the topology and add ions

import numpy as np
from polymerist.mdtools.openfftools.unitsys import openmm_to_openff
from openff.interchange.components import _packmol as packmol

from openff.interchange.components._packmol import *

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
    
    print(min_box_vecs)

    box_vecs = boxvectors.pad_box_vectors_uniform(min_box_vecs, padding)

    print("Pre-transform")
    print(box_vecs)

    if not np.all(box_vecs >= min_box_vecs): # TODO : change this to only check XYZ dimensions
        raise boxvectors.BoxVectorError(f'Desired box dimensions ({box_vecs}) are smaller than minimum excluded Topology dimensions ({min_box_vecs})')

    box_vecs = box_shape @ box_vecs

    print(box_vecs)
    
    # Translate all molecules in the topology so they are centered at the box_center
    # NOTE SUPER HACKY, MUST CHANGE FOR A DIFFERENT ENZYME, BUT FOR NOW
    # THIS WORKS AND WILL RECENTER LipA AND THE LIGAND AT CENTER OF BOX

    #com_protein = openmm_to_openff(mda.Universe("NH3_terminal_updated.pdb").select_atoms("all").center_of_mass()*angstrom)
    box_center = (np.sum(box_vecs, axis=0) / 2)
    
    com_traj = calculate_center_of_mass(topology)

    print("COM")
    print(com_traj)
    
    print("Box center")
    print(box_center)

    trans_vec = box_center - com_traj

    print(trans_vec)

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


print("Adding water and ions to the box!")
solvated_poly_pro_top = solvate_topology_with_ions(topology=poly_pro_top,
                                                 padding = openmm_to_openff(0.9 * nanometer),
                                                 nacl_conc = Quantity(0.1, "mole / liter"),
                                                 box_shape=packmol.RHOMBIC_DODECAHEDRON)#,
                                                 #box_shape=packmol.UNIT_CUBE)

print("Done adding solvent+ions!")

for mol in solvated_poly_pro_top.molecules:
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
        #print("meep")
        for atom in mol.atoms:
            atom.metadata['residue_name'] = "RBY"
        

solvated_poly_pro_top.to_file(assemble_path(WORKING_DIR, working_dir_tag, postfix='solvated', extension='pdb')) # saving extra PDB purely for visualization (SDF only lets you see molecules one-at-a-time)

##################################################################

# EDIT CIRCA 5/16/25
# NOTE making an interchange with many unique molecules is SLOW
# It is faster to make one interchange per unique molecule class
# and then run the Interchange.combine() function to make one final
# interchange that contains all of the chemical structures.
# The following hotfix does exactly that.


# create interchange object
ff = ForceField("ff14sb_off_impropers_0.0.4.offxml", 'openff-2.0.0.offxml') #, )

poly_mapper = {}

for poly in polymer_mol_list:
    poly_mapper[poly.to_smiles()] = poly

poly_mapper[off_lig_mol_best.to_smiles()] = off_lig_mol_best


all_incs = []
all_mols = []

for molecule in solvated_poly_pro_top.molecules:
    if molecule.to_smiles() == "[H][O][H]":
        continue

    if molecule.to_smiles() == "[Na+]":
        continue

    if molecule.to_smiles() == "[Cl-]":
        continue

    if molecule.to_smiles() in poly_mapper.keys():
        all_mols = [poly_mapper[molecule.to_smiles()]]
        print("found oligomer")
    else:
        all_mols = []

    #all_mols.append(solvents.water_TIP3P)

    all_incs.append(ff.create_interchange(molecule.to_topology(), charge_from_molecules=all_mols))

## Process the water + ions all at once, this saves computational overhead

print("Processing solvent + ions into an Interchange...")

water_ions_mols = []
prot_mol = []

for mol in solvated_poly_pro_top.molecules:
    if mol.to_smiles() == "[H][O][H]":
        water_ions_mols.append(mol)
    elif mol.to_smiles() == "[Na+]":
        water_ions_mols.append(mol)
    elif mol.to_smiles() == "[Cl-]":
        water_ions_mols.append(mol)

water_ion_inc = ff.create_interchange(topology=water_ions_mols, charge_from_molecules=[solvents.water_TIP3P])

print("Made solvent/ion Interchange")

# Combine all of the non-solvent components into one interchange
# Add in the solvent at the end to save compute

combined_inc = all_incs[0]

print("Combining components of Interchange...")
for inc in all_incs[1:]:
    combined_inc = combined_inc.combine(inc)

print("Adding in solvent...")

combined_inc = combined_inc.combine(water_ion_inc)

print("Re-setting the box vectors")
combined_inc.box = solvated_poly_pro_top.box_vectors

##################################################################

# TODO figure out how to iterate over all the unique molecule types in the topology
# I.E. grab all water molecules, export this to a new, empty topology, then make an
# interchange of just that component. Do this for all of the sub-components, and then
# combine the interchanges into one final, big interchange of the whole system

inc = combined_inc # should be good now

print("Done making Interchange!")
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

# Code to add in a flat bottom potential restraint between the Carbonyl carbon of the substrate
# And the O of the Catalytic Serine
# NOTE Must change this based on the enzyme/ substrate pair you are investigating

##########################################################################################
from openmm import CustomBondForce
from openmm.unit import nanometer, angstrom, kilojoule_per_mole

ser_O_group = []
RBY_grp = []

for atom in omm_topology.atoms():
    if atom.name == "OG" and atom.residue.index == 76:
        print(atom)
        ser_O_group.append(atom.index)

    if atom.residue.name == "RBY" and atom.index == 2739:
        print(atom)
        RBY_grp.append(atom.index)

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

# #################################################

# import MDAnalysis as mda
# from biopandas.pdb import PandasPdb

# def fix_prot_atoms_names(universe, prot_path = "NH3_terminal_His_proton_updated.pdb"):
    
#     # Takes in a MDAnalysis universe object and applies the correct
#     # naming convention to atoms in a protein

#     # Load the PDB file to match atom names
#     ppdb = PandasPdb().read_pdb(path=prot_path)
#     df = ppdb.df['ATOM']

#     prot_atom_grp = universe.select_atoms("protein")
#     assert len(df) == len(prot_atom_grp)

#     for k in range(len(prot_atom_grp)):
#         atom = prot_atom_grp[k]
#         df_atom = df.iloc[k]
#         assert df_atom['residue_number'] == atom.resid
#         assert df_atom['residue_name'] == atom.resname
#         atom.type = df_atom['atom_name']

#     return universe

# # Write the updated topology with the correct names for protein atoms:
# u = mda.Universe(
#                 f"{WORKING_DIR}/production/production_topology.pdb",
#                 f"{WORKING_DIR}/production/production_trajectory.dcd"
#             )

# u = fix_prot_atoms_names(u, prot_path=path_to_pro_pdb)

# all_atoms = u.select_atoms("all")

# for atom in u.select_atoms("protein"):
#     atom.name = atom.type

# for atom in u.select_atoms("not protein"):
#     if atom.chainID == "A":
#         atom.chainID = "B"

# with mda.Writer(f"{WORKING_DIR}/production/updated_topology.pdb") as pdb:
#     pdb.write(all_atoms)


