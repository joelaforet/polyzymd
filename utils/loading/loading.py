'''Functions used to load polymer-enzyme MD trajectories into MDAnalysis Universes'''

__author__ = 'Joe Laforet Jr.'
__email__ = 'jola3134@colorado.edu'

from sklearn.cluster import KMeans
from MDAnalysis.analysis import align
import numpy as np

from biopandas.pdb import PandasPdb

import MDAnalysis as mda
from MDAnalysis import transformations
from tqdm import tqdm
import os
import warnings

########################################################################################################################################################

# Hack to re-name all of the atoms in a universe that correspond to a protein
# to the names found in the file prot_path
# this SHOULD be a .pdb file that is Chemical Component Dictionary (CCD) compliant.
# NOTE not tested for systems with more than one protein...

def fix_prot_atoms_names(universe: mda.Universe,
                         prot_path: str="NH3_terminal_His_proton_updated.pdb") -> mda.Universe:
    
    # Takes in a MDAnalysis universe object and applies the correct
    # naming convention to atoms in a protein

    # Load the PDB file to match atom names
    ppdb = PandasPdb().read_pdb(path=prot_path)
    df = ppdb.df['ATOM']

    prot_atom_grp = universe.select_atoms("protein")
    assert len(df) == len(prot_atom_grp)

    for k in range(len(prot_atom_grp)):
        atom = prot_atom_grp[k]
        df_atom = df.iloc[k]
        assert df_atom['residue_number'] == atom.resid
        assert df_atom['residue_name'] == atom.resname
        atom.type = df_atom['atom_name']
        atom.name = atom.type

    return universe


##################################################################################################################################################

# Used to take in a MDAnalysis universe, and a PDB file containing the
# protein atoms names correctly, and writes a new PDB file for the 
# universe that contains the corrected names for protein atoms.
# NOTE: if you are visualizing a PDB in PyMol and the SS is not displaying
# correctly, try running this function and visualizing the new topology.

def write_corrected_protein_topology(universe: mda.Universe,
                                     prot_path: str="NH3_terminal_His_proton_updated.pdb",
                                     out_name: str= "Protein_Corrected_Topology.pdb"):

    u = fix_prot_atoms_names(universe)
    all_atoms = u.select_atoms("all")

    for atom in u.select_atoms("protein"):
        atom.name = atom.type

    for atom in u.select_atoms("not protein"):
        if atom.chainID == "A":
            atom.chainID = "B"

    with mda.Writer(out_name) as pdb:
        pdb.write(all_atoms)

    print(f"Successfully wrote new topology to: {out_name}")

    return None

###########################################################################################################################################


# This function identifies a centroid structure by applying KMeans clustering on the coordinates of the protein
# atoms in the trajectory.
# You can set equil_offset to only look at the traj[equil_offset:] frames.
# NOTE universe passed in must be a MDAnalysis universe and have valid protein atoms
# Consider using the fix_prot_atom_names() function above if your protein atoms are not in the convetion outlined by
# the Chemical Component Dictionary (CCD) https://www.wwpdb.org/data/ccd
# This site is an interactive way of looking up common naming conventions https://www.ebi.ac.uk/pdbe-srv/pdbechem/
# ALSO NOTE, PyMol won't display proteins correctly + you won't be able to perform protein specific functionality
# Like align if the atoms aren't named properly.

def get_representative_frame(universe: mda.Universe,
                             equil_offset: int=0) -> int:

    """Function computes the most representative configuration of a protein in a universe
       if equil_offset is specified, will only look at the frames after offset steps"""

    # Select alpha carbons of the protein
    alpha_carbons = universe.select_atoms("protein")

    print("grabbing coords")
    # Get the coordinates of alpha carbons for each frame
    coordinates = np.array([alpha_carbons.positions for ts in universe.trajectory[equil_offset:]])

    # Reshape coordinates to 2D array for K-means clustering
    n_frames, n_atoms, n_dim = coordinates.shape
    coordinates_reshaped = coordinates.reshape(n_frames, n_atoms * n_dim)

    print("doing Kmeans")
    # Perform K-means clustering with 1 cluster
    kmeans = KMeans(n_clusters=1)
    kmeans.fit(coordinates_reshaped)

    # Get the cluster center
    cluster_center = kmeans.cluster_centers_[0]

    # Reshape cluster center back to original shape
    cluster_center_reshaped = cluster_center.reshape(n_atoms, n_dim)

    print("Finding most representative frame")
    # Find the frame that is most representative of the cluster (closest to the cluster center)
    distances = np.linalg.norm(coordinates - cluster_center_reshaped, axis=(1, 2))
    representative_frame_index = np.argmin(distances) + equil_offset

    print(f"The most representative frame is frame number {representative_frame_index}.")

    return representative_frame_index

#######################################################################################################################################################

# Suppress some MDAnalysis warnings about writing PDB files
warnings.filterwarnings('ignore')

# Takes in a dictionary structured like:
# temp_reps = {temp_1: num_reps,
#              temp_2: num_reps...}
# Will process all replicates of that given temperature and append the universes to list

# Each entry in this list will be a tuple of universes structured like [(full_traj, ref_frame)]
# The full_traj will be all frames of the simulation with the following transforms applied:
# 1. Center protein and ligand in the Periodic Cell NOTE the ligand is hard-coded as 'resname RBY', change this for
# a different ligand.
# 2. Performs KMeans clustering on the trajectory to identify a centroid structure. 
# 3. Given the centroid structure, align all other frames to this one.
# 4. Write newly centered+aligned trajectory to 'prot_centered_trajectory.dcd'

# Depending on your applications, it may be desired to use the unaligned version of the universe.

# The entries of the list will be [temp_1[0], temp_1[1]...temp_1[num_reps-1], temp_2[0], temp_2[1]...temp_2[num_reps-1]]
# for all temp/reps present in the temp_reps dictionary

# "simulation_id": "LipA_Resorufin-Butyrate" (stuff before thermo_params)
# "thermo_id": "0.5ns-NVT_1000.0ns-NPT", (formatted like timeEQ-EQ_ensemble_timePROD-PROD_ensemble)
# "dir_id": "Docked_Substrate/Resorufin-Butyrate/Simulating_just_LipA", (directory path where sim files live)
# "temp_reps": {293: 5, 323: 5, 363: 5}, {temp_1: num_reps_1, temp_2: num_reps_2, etc}
# "equilibrium_offset": 500, 
# "prot_path" : "NH3_terminal_His_proton_updated.pdb" relative path to correctly named pdb file

# NOTE Universes have the FULL amount of frames loaded in them,
# in this case, equilibrium offset is used to section off what frames
# get included in the clustering calculation
# For later analysis, you must filter out the first equilibrium_offset frames
# I may change this later, but for now this is what we have.

def load_or_transform_universes(simulation_id: str,
                                thermo_id: str,
                                dir_id: str,
                                temperatures_replicates: dict,
                                equilibrium_offset: int=0,
                                prot_path: str="NH3_terminal_His_proton_updated.pdb",
                                ligand_name: str="RBY")->list:
    
    u_arr = []
    for temp, num_replicates in temperatures_replicates.items():
        for run in range(1, num_replicates + 1):
            transformed_trajectory = f"{dir_id}/{simulation_id}_{float(temp)}K_{thermo_id}_run{run}/production/prot_centered_trajectory.dcd"
            ref_output_topology = f"{dir_id}/{simulation_id}_{float(temp)}K_{thermo_id}_run{run}/production/production_topology.pdb"
            ref_output_trajectory = f"{dir_id}/{simulation_id}_{float(temp)}K_{thermo_id}_run{run}/production/ref_prot_centered_trajectory.dcd"

            if os.path.exists(transformed_trajectory) and os.path.exists(ref_output_topology) and os.path.exists(ref_output_trajectory):
                # Load the transformed universe and reference structure from file
                print("Cached data found!")
                universe = mda.Universe(f"{dir_id}/{simulation_id}_{float(temp)}K_{thermo_id}_run{run}/production/production_topology.pdb", transformed_trajectory)
                ref = mda.Universe(ref_output_topology, ref_output_trajectory)

                universe = fix_prot_atoms_names(universe=universe, prot_path=prot_path)
                ref = fix_prot_atoms_names(universe=ref, prot_path=prot_path)


            else:
                # Original loading and transformation process
                universe = mda.Universe(
                    f"{dir_id}/{simulation_id}_{float(temp)}K_{thermo_id}_run{run}/production/production_topology.pdb",
                    f"{dir_id}/{simulation_id}_{float(temp)}K_{thermo_id}_run{run}/production/production_trajectory.dcd"
                )

                # Fix atom names so MDAnalysis recognizes keywords like protein
                universe = fix_prot_atoms_names(universe=universe, prot_path=prot_path)

                # Workflow for transformations to center protein in box
                center_prot = transformations.center_in_box(ag=universe.select_atoms(f'protein or resname {ligand_name}'), center='mass')
                workflow = [center_prot]
                universe.trajectory.add_transformations(*workflow)

                # Get representative frame index
                representative_frame_index = get_representative_frame(universe=universe, equil_offset=equilibrium_offset)

                # Calculate average structure and align trajectory
                print("Running AverageStructure")
                average = align.AverageStructure(universe, universe, select='protein and type CA', ref_frame=representative_frame_index).run(verbose=True)
                ref = average.results.universe
                print("Starting AlignTraj")
                aligner = align.AlignTraj(universe, ref, select='protein and type CA', in_memory=False).run(verbose=True)
                del aligner   
                del average


                all_atoms = universe.select_atoms("all")
                all_atoms.write(transformed_trajectory, frames='all')
                all_atoms.write(ref_output_trajectory, frames=universe.trajectory[[representative_frame_index]])

            u_arr.append((universe, ref))

    return u_arr


# DEPRECATED: Loading universe function, but without the load_from_cache feature
# may be useful for some other purpose so I won't delete it.

'''def load_universes(simulation_id, thermo_id, dir_id, temperatures_replicates, equilibrium_offset=0, prot_path = "NH3_terminal_updated.pdb"):
    u_arr = []
    for temp, num_replicates in temperatures_replicates.items():
        for run in range(1, num_replicates + 1):
            universe = mda.Universe(
                f"{dir_id}/{simulation_id}_{float(temp)}K_{thermo_id}_run{run}/production/production_topology.pdb",
                f"{dir_id}/{simulation_id}_{float(temp)}K_{thermo_id}_run{run}/production/production_trajectory.dcd"
            )

            # need to fix atom names so MDAnalysis recognizes keywords like protein
            universe = fix_prot_atoms_names(universe=universe, prot_path=prot_path)

            # workflow for transformations to center protein in box
            center_prot = transformations.center_in_box(ag = universe.select_atoms('protein'), center='mass')

            #workflow = [noJump_trans, center_prot_poly]

            workflow=[center_prot]

            universe.trajectory.add_transformations(*workflow)

            #ag = universe.select_atoms("protein")
            
            representative_frame_index = get_representative_frame(universe=universe, equil_offset=equilibrium_offset)

            print("running AverageStructure")
            average = align.AverageStructure(universe, universe, select='protein and type CA', ref_frame=representative_frame_index).run(verbose=True)
            ref = average.results.universe
            print("starting AlignTraj")
            aligner = align.AlignTraj(universe, ref, select='protein and type CA', in_memory=False).run(verbose=True)
            del aligner   
            del average
             
            u_arr.append((universe, ref))

    return u_arr
    '''

##################################################################################################################################################################

# Below function is useful for examining one universe at a time
# Use this when you have a specfic temp/rep/simulation combination
# you want to analyze. This fxn is also called by the parallel 
# submit_distance_jobs.py script.

# NOTE this function does NOT align the trajectory, only loads it in place
def load_universe(simulation_id: str,
                  thermo_id: str,
                  dir_id: str,
                  temp: int,
                  run: int,
                  prot_path: str="NH3_terminal_updated.pdb",
                  ligand_name: str="RBY")->mda.Universe:
    """Load and prepare a single universe with transformations applied"""
    
    # Try to load pre-processed trajectory first
    transformed_trajectory = f"{dir_id}/{simulation_id}_{float(temp)}K_{thermo_id}_run{run}/production/prot_centered_trajectory.dcd"
    topology_file = f"{dir_id}/{simulation_id}_{float(temp)}K_{thermo_id}_run{run}/production/production_topology.pdb"
    
    if os.path.exists(transformed_trajectory) and os.path.exists(topology_file):
        print(f"Loading pre-processed trajectory: {transformed_trajectory}")
        universe = mda.Universe(topology_file, transformed_trajectory)
    else:
        # Load original trajectory and apply transformations
        print(f"Loading original trajectory and applying transformations")
        universe = mda.Universe(
            f"{dir_id}/{simulation_id}_{float(temp)}K_{thermo_id}_run{run}/production/production_topology.pdb",
            f"{dir_id}/{simulation_id}_{float(temp)}K_{thermo_id}_run{run}/production/production_trajectory.dcd"
        )
        
        # Apply transformations
        center_prot = transformations.center_in_box(ag=universe.select_atoms(f'protein or resname {ligand_name}'), center='mass')
        workflow = [center_prot]
        universe.trajectory.add_transformations(*workflow)
    
    # Fix atom names
    universe = fix_prot_atoms_names(universe=universe, prot_path=prot_path)
    
    return universe

###########################################################################################################################################

def load_or_transform_daisychain_universes(simulation_id: str,
                                thermo_id: str,
                                dir_id: str,
                                temperatures_replicates: dict,
                                equilibrium_offset: int=0,
                                prot_path: str="NH3_terminal_His_proton_updated.pdb",
                                ligand_name: str="RBY")->list:
    
    u_arr = []
    for temp, num_replicates in temperatures_replicates.items():
        for run in range(1, num_replicates + 1):
            base_path = f"{dir_id}/{simulation_id}_{float(temp)}K_{thermo_id}_run{run}"
            
            # Define cached file paths (assuming we cache in the first production directory)
            transformed_trajectory = f"{base_path}/production_0/prot_centered_trajectory.dcd"
            ref_output_topology = f"{base_path}/production_0/production_0_topology.pdb"
            ref_output_trajectory = f"{base_path}/production_0/ref_prot_centered_trajectory.dcd"

            legacy_transformed_trajectory = f"{base_path}/production/prot_centered_trajectory.dcd"
            legacy_ref_output_topology = f"{base_path}/production/production_topology.pdb"
            legacy_ref_output_trajectory = f"{base_path}/production/ref_prot_centered_trajectory.dcd"

            if os.path.exists(transformed_trajectory) and os.path.exists(ref_output_topology) and os.path.exists(ref_output_trajectory):
                # Load the transformed universe and reference structure from file
                print("Cached data found!")
                universe = mda.Universe(ref_output_topology, transformed_trajectory)
                ref = mda.Universe(ref_output_topology, ref_output_trajectory)

                universe = fix_prot_atoms_names(universe=universe, prot_path=prot_path)
                ref = fix_prot_atoms_names(universe=ref, prot_path=prot_path)

            elif os.path.exists(legacy_transformed_trajectory) and os.path.exists(legacy_ref_output_topology) and os.path.exists(legacy_ref_output_trajectory):
                # Load the transformed universe and reference structure from file
                print("Cached data found!")
                print("Loading from Legacy format (no Daisy Chain)")
                universe = mda.Universe(legacy_ref_output_topology, legacy_transformed_trajectory)
                ref = mda.Universe(legacy_ref_output_topology, legacy_ref_output_trajectory)

                universe = fix_prot_atoms_names(universe=universe, prot_path=prot_path)
                ref = fix_prot_atoms_names(universe=ref, prot_path=prot_path)

            else:
                # Find all production_N directories and collect all DCD files
                import glob
                production_dirs = glob.glob(f"{base_path}/production_*")
                
                # Sort directories numerically
                def extract_dir_number(dirpath):
                    import re
                    match = re.search(r'production_(\d+)$', dirpath)
                    return int(match.group(1)) if match else 0
                
                production_dirs.sort(key=extract_dir_number)
                
                if not production_dirs:
                    print(f"Warning: No production directories found in {base_path}")
                    continue
                
                # Collect all DCD files from all production directories
                all_dcd_files = []
                topology_file = None
                
                for prod_dir in production_dirs:
                    # Get all DCD files in this production directory
                    dcd_pattern = f"{prod_dir}/production_*.dcd"
                    dcd_files = glob.glob(dcd_pattern)
                    
                    # Sort DCD files within this directory numerically
                    def extract_dcd_number(filepath):
                        import re
                        match = re.search(r'production_(\d+)\.dcd$', filepath)
                        return int(match.group(1)) if match else 0
                    
                    dcd_files.sort(key=extract_dcd_number)
                    all_dcd_files.extend(dcd_files)
                    
                    # Get topology file from the first directory (assuming it's the same across all)
                    if topology_file is None:
                        potential_topology = f"{prod_dir}/production_0_topology.pdb"
                        if os.path.exists(potential_topology):
                            topology_file = potential_topology
                
                if not all_dcd_files:
                    print(f"Warning: No production DCD files found in any production directories")
                    continue
                    
                if topology_file is None:
                    print(f"Warning: No topology file found in production directories")
                    continue
                
                print(f"Loading {len(all_dcd_files)} DCD files from {len(production_dirs)} production directories")
                
                # Load universe with topology and all DCD files
                universe = mda.Universe(topology_file, *all_dcd_files)

                # Fix atom names so MDAnalysis recognizes keywords like protein
                universe = fix_prot_atoms_names(universe=universe, prot_path=prot_path)

                # Workflow for transformations to center protein in box
                center_prot = transformations.center_in_box(ag=universe.select_atoms(f'protein or resname {ligand_name}'), center='mass')
                workflow = [center_prot]
                universe.trajectory.add_transformations(*workflow)

                # Get representative frame index
                #representative_frame_index = get_representative_frame(universe=universe, equil_offset=equilibrium_offset)

                # Calculate average structure and align trajectory
                print("Running AverageStructure")
                representative_frame_index = 0 # reference is just the first frame, want to refer to equilibrated structure
                average = align.AverageStructure(universe, universe, select='protein and type CA', ref_frame=representative_frame_index).run(verbose=True)
                ref = average.results.universe
                print("Starting AlignTraj")
                aligner = align.AlignTraj(universe, ref, select='protein and type CA', in_memory=False).run(verbose=True)
                del aligner   
                del average

                # Save transformed trajectories
                all_atoms = universe.select_atoms("all")
                all_atoms.write(transformed_trajectory, frames='all')
                all_atoms.write(ref_output_trajectory, frames=universe.trajectory[[representative_frame_index]])

            u_arr.append((universe, ref))

    return u_arr
