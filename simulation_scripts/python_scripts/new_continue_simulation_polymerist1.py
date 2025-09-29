"""
Continuation script for MD simulations using native OpenMM.
Loads state from previous segment and continues simulation.
"""

import sys
import json
import argparse
from pathlib import Path
import logging
from openmm.app import PDBFile, DCDReporter, StateDataReporter, CheckpointReporter, Simulation
from openmm.app.topology import Topology
from openmm import XmlSerializer, System, State, LangevinMiddleIntegrator, MonteCarloBarostat
from openmm.unit import Quantity, nanoseconds, kelvin, bar, femtoseconds, picoseconds
import openmm.unit as u

def quantity_from_dict(qdict):
    """Convert serialized quantity dictionary back to OpenMM Quantity."""
    value = qdict["__values__"]["value"]
    unit_str = qdict["__values__"]["unit"]

    if unit_str.startswith("/"):
        base_unit = getattr(u, unit_str[1:])
        return value / base_unit
    else:
        return value * getattr(u, unit_str)

def setup_logging(log_file):
    """Setup logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

def create_integrator_from_params(param_dict):
    """Create OpenMM integrator from parameter dictionary."""
    # Extract integrator parameters
    integ_raw = param_dict["__values__"]["integ_params"]["__values__"]
    time_step = quantity_from_dict(integ_raw["time_step"])
    
    # Extract thermodynamic parameters
    thermo_raw = param_dict["__values__"]["thermo_params"]["__values__"]
    temperature = quantity_from_dict(thermo_raw["temperature"])
    friction_coeff = quantity_from_dict(thermo_raw["friction_coeff"])
    
    # Create Langevin integrator
    integrator = LangevinMiddleIntegrator(temperature, friction_coeff, time_step)
    
    return integrator

def add_barostat_if_needed(system, param_dict):
    """Add barostat to system if pressure is specified."""
    thermo_raw = param_dict["__values__"]["thermo_params"]["__values__"]
    ensemble = thermo_raw["ensemble"]
    
    if ensemble.upper() == "NPT":
        temperature = quantity_from_dict(thermo_raw["temperature"])
        pressure = quantity_from_dict(thermo_raw["pressure"])
        barostat_freq = thermo_raw.get("barostat_freq", 25)
        
        # Check if barostat already exists
        has_barostat = any(isinstance(force, MonteCarloBarostat) 
                          for force in system.getForces())
        
        if not has_barostat:
            barostat = MonteCarloBarostat(pressure, temperature, barostat_freq)
            system.addForce(barostat)
            print(f"Added barostat: {pressure} at {temperature}, frequency {barostat_freq}")

def setup_reporters(simulation, working_dir, segment_index, param_dict):
    """Setup reporters for the simulation."""
    # Extract reporter parameters
    reporter_raw = param_dict["__values__"]["reporter_params"]["__values__"]
    
    report_checkpoint = reporter_raw.get("report_checkpoint", True)
    report_state = reporter_raw.get("report_state", True)
    report_trajectory = reporter_raw.get("report_trajectory", True)
    report_state_data = reporter_raw.get("report_state_data", True)
    traj_ext = reporter_raw.get("traj_ext", "dcd")
    num_steps = reporter_raw.get("num_steps", 100000)
    
    # Extract integrator parameters for step calculation
    integ_raw = param_dict["__values__"]["integ_params"]["__values__"]
    total_time = quantity_from_dict(integ_raw["total_time"])
    time_step = quantity_from_dict(integ_raw["time_step"])
    num_samples = integ_raw["num_samples"]
    
    # Calculate total steps and reporting interval
    total_steps = int(total_time / time_step)
    report_interval = max(1, total_steps // num_samples)
    
    output_dir = working_dir / f"production_{segment_index}"
    output_dir.mkdir(exist_ok=True)
    
    # Trajectory reporter
    if report_trajectory:
        traj_file = output_dir / f"production_{segment_index}_trajectory.{traj_ext}"
        if traj_ext.lower() == "dcd":
            simulation.reporters.append(DCDReporter(str(traj_file), report_interval))
        else:
            print(f"Warning: Trajectory format {traj_ext} not supported, using DCD")
            simulation.reporters.append(DCDReporter(str(traj_file), report_interval))
    
    # State data reporter
    if report_state_data:
        state_file = output_dir / f"production_{segment_index}_state_data.csv"
        simulation.reporters.append(StateDataReporter(
            str(state_file), 
            report_interval,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            speed=True
        ))
    
    # Checkpoint reporter
    if report_checkpoint:
        checkpoint_file = output_dir / f"production_{segment_index}_checkpoint.chk"
        simulation.reporters.append(CheckpointReporter(str(checkpoint_file), report_interval))
    
    return total_steps, output_dir

def save_final_state(simulation, output_dir, segment_index):
    """Save the final state and system after simulation."""
    # Save final state
    state_file = output_dir / f"production_{segment_index}_state.xml"
    state = simulation.context.getState(
        getPositions=True, 
        getVelocities=True, 
        getForces=True, 
        getEnergy=True,
        getParameters=True,
        enforcePeriodicBox=True
    )
    
    with open(state_file, 'w') as f:
        f.write(XmlSerializer.serialize(state))
    
    # Save system
    system_file = output_dir / f"production_{segment_index}_system.xml"
    with open(system_file, 'w') as f:
        f.write(XmlSerializer.serialize(simulation.system))
    
    print(f"Saved final state to: {state_file}")
    print(f"Saved system to: {system_file}")

def find_solvated_pdb(working_dir):
    """Find the solvated PDB file in the working directory."""
    working_dir = Path(working_dir)
    
    # Look for common solvated PDB file patterns
    patterns = [
        "*solvated.pdb",
        "*_solvated.pdb", 
        "solvated_*.pdb",
        "*.pdb"  # fallback to any PDB file
    ]
    
    for pattern in patterns:
        pdb_files = list(working_dir.glob(pattern))
        if pdb_files:
            return pdb_files[0]
    
    raise FileNotFoundError(f"Could not find solvated PDB file in {working_dir}")

def main():
    parser = argparse.ArgumentParser(description="Continue MD simulation from previous segment")
    parser.add_argument("-s", "--segment_index", type=int, required=True,
                        help="Current segment index (1-based)")
    parser.add_argument("-w", "--working_dir", type=str, required=True,
                        help="Working directory path")
    parser.add_argument("-t", "--segment_time", type=float, required=True,
                        help="Time for this segment in nanoseconds")
    parser.add_argument("-n", "--num_samples", type=int, default=50,
                        help="Number of frames to save for this segment (default: 50)")
    
    args = parser.parse_args()
    
    # Setup paths
    segment_index = args.segment_index
    working_dir = Path(args.working_dir)
    prev_segment = segment_index - 1
    
    print(f"Starting continuation of segment {segment_index}")
    print(f"Working directory: {working_dir}")
    print(f"Previous segment: {prev_segment}")
    
    # Paths to prior outputs
    state_path = working_dir / f"production_{prev_segment}/production_{prev_segment}_state.xml"
    sys_path = working_dir / f"production_{prev_segment}/production_{prev_segment}_system.xml"
    param_path = working_dir / f"production_{prev_segment}/production_{prev_segment}_parameters.json"
    
    # Check if required files exist
    for path in [state_path, sys_path, param_path]:
        print(path)
        if not path.exists():
            raise FileNotFoundError(f"Required file not found: {path}")

    sim_params_path = assemble_path(working_dir / f"production_{prev_segment}", f'production_{prev_segment}_parameters', extension='json')
    sim_params_from_file = SimulationParameters.from_file(sim_params_path)

    sim_paths = SimulationPaths.from_dir_and_parameters(working_dir/f'production_{prev_segment}', sim_params=sim_params_from_file)
    # print(sim_paths.parameters_path)
    # print(sim_paths.topology_path)
    # print(sim_paths.system_path)
    # print(sim_paths.trajectory_path)

    # Find solvated PDB file
    solvated_pdb_path = find_solvated_pdb(working_dir)
    print(f"Using solvated PDB: {solvated_pdb_path}")
    
    # Setup logging
    log_file = working_dir / f'simulation_status_segment_{segment_index}.log'
    logger = setup_logging(log_file)
    
    try:
        # Load OpenMM components
        logger.info("Loading OpenMM system...")
        # with open(sys_path, 'r') as f:
        #     omm_system: System = XmlSerializer.deserialize(f.read())

        with open(sim_paths.system_path, 'r') as f:
            omm_system: System = XmlSerializer.deserialize(f.read())
        
        logger.info("Loading topology...")
        omm_topology: Topology = PDBFile(str(solvated_pdb_path)).topology
        
        # # Load simulation parameters
        # logger.info("Loading simulation parameters...")
        # with open(param_path, 'r') as f:
        #     param_dict = json.load(f)
        
        # # Update parameters for this segment
        # param_dict["__values__"]["integ_params"]["__values__"]["total_time"] = {
        #     "__values__": {"value": args.segment_time, "unit": "nanoseconds"}
        # }
        # param_dict["__values__"]["integ_params"]["__values__"]["num_samples"] = args.num_samples
        
        # logger.info(f"Segment time: {args.segment_time} ns")
        # logger.info(f"Number of samples: {args.num_samples}")
        
        # # Add barostat if needed
        # add_barostat_if_needed(omm_system, param_dict)
        
        # # Create integrator
        # logger.info("Creating integrator...")
        # integrator = create_integrator_from_params(param_dict)
        
        # # Create simulation
        # logger.info("Creating simulation...")
        # simulation = Simulation(omm_topology, omm_system, integrator)
        
        # # Load state from previous segment
        # logger.info("Loading state from previous segment...")
        # simulation.loadState(str(state_path))
        
        # # Setup reporters
        # logger.info("Setting up reporters...")
        # total_steps, output_dir = setup_reporters(simulation, working_dir, segment_index, param_dict)
        
        # # Save parameters for this segment
        # param_file = output_dir / f"production_{segment_index}_parameters.json"
        # with open(param_file, 'w') as f:
        #     json.dump(param_dict, f, indent=2)
        
        # # Run simulation
        # logger.info(f"Starting simulation for segment {segment_index}...")
        # logger.info(f"Total steps: {total_steps}")
        
        # simulation.step(total_steps)
        
        # # Save final state
        # save_final_state(simulation, output_dir, segment_index)
        

        from polymerist.genutils.logutils.IOHandlers import MSFHandlerFlex 
        from polymerist.mdtools.openmmtools.execution import run_simulation_schedule

        # MSFHandler compiles log and error events into single, "master" log file
        logpath = assemble_path(POLYMER_SIM_DIR, 'simulation_status', 'log')
        with MSFHandlerFlex(filename=logpath, proc_name=f'{polymer_name}_sims') as logger:
            history = run_simulation_schedule(
                POLYMER_SIM_DIR, 
                schedule={ # simulations will be run in the order they appear here
                    f'production_{segment_index}'  : sim_params_from_file,
                }, 
                init_top=omm_topology,
                init_sys=omm_system,
                init_pos=interchange.get_positions(include_virtual_sites=True).to_openmm(),
                return_history=True
            )



        logger.info(f"Segment {segment_index} completed successfully!")
        return 0
        
    except Exception as e:
        logger.error(f"Error during simulation: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())