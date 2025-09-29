# Populate the workspace for CoPolymer simulations with Enzymes (including controls)
import signac
import itertools

project = signac.init_project("project_Proteins_and_CoPolymers")

# https://signac.readthedocs.io/en/latest/recipes.html#defining-a-grid-of-state-point-values

def grid(gridspec):
    """Yields the Cartesian product of a `dict` of iterables.

    The input ``gridspec`` is a dictionary whose keys correspond to
    parameter names. Each key is associated with an iterable of the
    values that parameter could take on. The result is a sequence of
    dictionaries where each dictionary has one of the unique combinations
    of the parameter values.
    """
    for values in itertools.product(*gridspec.values()):
        yield dict(zip(gridspec.keys(), values))


def generate_all_statepoints():
    """Generate both polymer and control simulation statepoints."""
    
    # Common parameters for both polymer and control simulations
    temperatures = [293.0, 323.0, 363.0]
    proteins = ["LipA"]
    ligands = ["Resorufin-Butyrate"]
    ligand_restraint_dist = [10.0]  # angstroms
    
    # Simulation parameters
    oligomer_lengths = [5]  # just 5-mers for now
    number_of_oligomers = [38]  # only testing 38 in the box
    
    equilibration_times = [0.5]  # time in ns to do equilibration
    production_times = [1000.0]  # time in ns to do production
    
    # Ensembles to use for MD
    equilibration_ensemble = ['NVT']
    production_ensemble = ['NPT']
    
    # Number of replicates
    N = 5
    replicate_nums = [x+1 for x in range(N)]
    
    # Trajectory logging frequency (in ns)
    log_freq = [0.4]
    
    statepoints = []
    
    # Generate polymer simulations
    polymer_1_SMILES = ['CC(Cl)(CCl)C(=O)OCC[N+](C)(C)CCC[S+]([O-])([O-])=O']  # SBMA
    polymer_2_SMILES = ['CC(Cl)(CCl)C(=O)OCCOC1=CC=CC=C1']  # EGPMA
    polymer_ratios = [0.0, 0.01, 0.02, 0.05, 0.10]  # By convention, this is fraction of co-polymer that is polymer_2
    
    polymer_statepoint_grid = {
        "SMILES_poly_1": polymer_1_SMILES,
        "SMILES_poly_2": polymer_2_SMILES,
        "protein": proteins,
        "ligand": ligands,
        "ligand_restraint_dist": ligand_restraint_dist,
        "polymer_ratio": polymer_ratios,
        "olig_length": oligomer_lengths,
        "num_oligomers": number_of_oligomers,
        "temperature": temperatures,
        "eq_time": equilibration_times,
        "eq_ensemble": equilibration_ensemble,
        "prod_time": production_times,
        "prod_ensemble": production_ensemble,
        "replicate_num": replicate_nums,
        "traj_log_freq": log_freq
    }
    
    # Add polymer statepoints
    for sp in grid(polymer_statepoint_grid):
        statepoints.append(sp)
    
    # Generate control simulations (no polymer)
    control_statepoint_grid = {
        "SMILES_poly_1": [None],  # No polymer 1
        "SMILES_poly_2": [None],  # No polymer 2
        "protein": proteins,
        "ligand": ligands,
        "ligand_restraint_dist": ligand_restraint_dist,
        "polymer_ratio": [None],  # Not applicable for controls
        "olig_length": [None],    # Not applicable for controls
        "num_oligomers": [0],     # No oligomers in control
        "temperature": temperatures,
        "eq_time": equilibration_times,
        "eq_ensemble": equilibration_ensemble,
        "prod_time": production_times,
        "prod_ensemble": production_ensemble,
        "replicate_num": replicate_nums,
        "traj_log_freq": log_freq
    }
    
    # Add control statepoints
    for sp in grid(control_statepoint_grid):
        statepoints.append(sp)
    
    return statepoints


def main():
    """Initialize all statepoints in the Signac project."""
    
    count = 0
    
    # Generate all statepoints (polymer + control)
    all_statepoints = generate_all_statepoints()
    
    for sp in all_statepoints:
        print("Initializing job", sp)
        project.open_job(sp).init()
        count += 1
    
    # Count polymer vs control jobs for summary
    polymer_jobs = len([sp for sp in all_statepoints if sp["SMILES_poly_1"] is not None])
    control_jobs = len([sp for sp in all_statepoints if sp["SMILES_poly_1"] is None])
    
    print(f"\nMade {count} total statepoints!")
    print(f"  - Polymer simulations: {polymer_jobs}")
    print(f"  - Control simulations: {control_jobs}")
    
    # Print some example queries you can use
    print(f"\nExample queries:")
    print(f"  Find all control jobs: project.find_jobs({{'SMILES_poly_1': None}})")
    print(f"  Find all polymer jobs: project.find_jobs({{'SMILES_poly_1': {{'$ne': None}}}})")
    print(f"  Find jobs at 293K: project.find_jobs({{'temperature': 293.0}})")
    print(f"  Find 1% polymer ratio jobs: project.find_jobs({{'polymer_ratio': 0.01}})")


if __name__ == "__main__":
    main()