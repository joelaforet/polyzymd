#!/usr/bin/env python3
"""
Daisy-chain MD simulation submission script for HPC SLURM scheduler.
Breaks long simulations into smaller dependent jobs.
"""

"""
python daisy_chain_submission_NOPOLY.py -p LipA -o SBMA-EGPMA -n 2 -t 293 -e 0.5 -r 1 -pp NH3_terminal_His_proton_updated.pdb -ps ATRP_SBMA_EGPMA_5-mer -rep 1-2 -chars A B -probs 0.98 0.02 -plen 5 -prefix SBMA-EGPMA -seg 5 --preset blanca-shirts --frames 100
"""

import argparse
import subprocess
import sys
import os
from pathlib import Path
import re

def parse_replicate_range(replicate_range):
    """Parse SLURM array range into a list of individual replicate numbers."""
    replicates = []
    
    # Handle comma-separated ranges
    parts = replicate_range.split(',')
    
    for part in parts:
        part = part.strip()
        if '-' in part:
            # Handle range like "1-5" or "1-5:2"
            if ':' in part:
                range_part, step = part.split(':')
                step = int(step)
            else:
                range_part = part
                step = 1
            
            start, end = map(int, range_part.split('-'))
            replicates.extend(range(start, end + 1, step))
        else:
            # Single number
            replicates.append(int(part))
    
    return sorted(list(set(replicates)))  # Remove duplicates and sort

INITIAL_NOPOLY_JOB_SCRIPT_TEMPLATE ="""#!/bin/bash \n
#SBATCH --partition={partition} \n
#SBATCH --job-name=i_{job_name} \n
#SBATCH --output={job_output} \n
#SBATCH --qos={qos} \n
#SBATCH --nodes=1 \n
#SBATCH --ntasks=1 \n
#SBATCH --mem=3G \n
#SBATCH --time={time_limit} \n
#SBATCH --gres=gpu:1 \n
#SBATCH --mail-type=FAIL \n
#SBATCH --mail-user={email} \n
#SBATCH --account={account} \n

module purge
module load miniforge
mamba activate polymerist-env

# Run the original simulation script with modified production time
python sim_enz_only_DOCKED_RESTRAINT.py \\
    -p {protein_name} \\
    -t {temperature} \\
    -e {eq_time} \\
    -r {segment_prod_time} \\
    -pp {path_to_pro_pdb} \\
    -rep {replicate_num} \\
    --frames {segment_frames_to_save} \\
    --num_segs {total_segments}

echo "Segment {segment_index} completed successfully"
"""


def create_initial_job_script(args, replicate_num, segment_index=0, total_segments=10):
    """Fill in the initial simulation job script using a predefined template."""

    segment_prod_time = args.prod_time / total_segments
    segment_frames_to_save = int(args.frames / total_segments)

    job_script = INITIAL_NOPOLY_JOB_SCRIPT_TEMPLATE.format(
        partition=args.partition,
        job_name=f"s{segment_index}_r{replicate_num}_{args.temperature}K_{args.protein_name}",
        job_output=f"s{segment_index}_r{replicate_num}_{args.temperature}K_{args.protein_name}%.%A_%a.out",
        qos=args.qos,
        time_limit=args.time_limit,
        email=args.email,
        account=args.account,
        protein_name=args.protein_name,
        temperature=args.temperature,
        eq_time=args.eq_time,
        segment_prod_time=segment_prod_time,
        path_to_pro_pdb=args.path_to_pro_pdb,
        replicate_num=replicate_num,
        segment_frames_to_save=segment_frames_to_save,
        total_segments=total_segments,
        segment_index=segment_index
    )

    return job_script

# CONTINUATION_JOB_SCRIPT_TEMPLATE = """#!/bin/bash \n
# #SBATCH --partition={partition} \n
# #SBATCH --job-name=c_{job_name} \n
# #SBATCH --output={job_output} \n
# #SBATCH --qos={qos} \n
# #SBATCH --nodes=1 \n
# #SBATCH --ntasks=1 \n
# #SBATCH --mem=3G \n
# #SBATCH --time={time_limit} \n
# #SBATCH --gres=gpu:1 \n
# #SBATCH --mail-type=FAIL \n
# #SBATCH --mail-user={email} \n
# #SBATCH --account={account} \n

# module purge
# module load miniforge
# mamba activate polymerist-env

# # Run the continuation script
# python new_continue_simulation.py \\
#     -s {segment_index} \\
#     -w {working_dir} \\
#     -t {segment_prod_time} \\
#     -n {segment_frames_to_save}

# echo "Segment {segment_index} completed successfully"
# """

# UPDATED for Polymerist 1.0 API Changes
CONTINUATION_JOB_SCRIPT_TEMPLATE = """#!/bin/bash \n
#SBATCH --partition={partition} \n
#SBATCH --job-name=c_{job_name} \n
#SBATCH --output={job_output} \n
#SBATCH --qos={qos} \n
#SBATCH --nodes=1 \n
#SBATCH --ntasks=1 \n
#SBATCH --mem=3G \n
#SBATCH --time={time_limit} \n
#SBATCH --gres=gpu:1 \n
#SBATCH --mail-type=FAIL \n
#SBATCH --mail-user={email} \n
#SBATCH --account={account} \n

module purge
module load miniforge
mamba activate polymerist-env

# Run the continuation script
python new_continue_simulation_polymerist2.py \\
    -s {segment_index} \\
    -w {working_dir} \\
    -t {segment_prod_time} \\
    -n {segment_frames_to_save}

echo "Segment {segment_index} completed successfully"
"""


def create_continuation_job_script(args, replicate_num, segment_index, total_segments=10, ligand="Resorufin-Butyrate"):
    """Fill in the continuation simulation job script using a predefined template."""

    # Calculate the production time for this segment
    segment_prod_time = args.prod_time / total_segments
    segment_frames_to_save = int(args.frames / total_segments)

    working_dir = (
        f"3,3A_RESTRAINT_{args.protein_name}_{ligand}_{args.temperature}K_"
        f"{args.eq_time}ns-NVT_{args.prod_time}ns-NPT_run{replicate_num}"
    )

    job_script = CONTINUATION_JOB_SCRIPT_TEMPLATE.format(
        partition=args.partition,
        job_name=f"s{segment_index}_r{replicate_num}_{args.temperature}_{args.protein_name}",
        job_output=f"s{segment_index}_r{replicate_num}_{args.temperature}_{args.protein_name}%.%A_%a.out",
        qos=args.qos,
        time_limit=args.time_limit,
        email=args.email,
        account=args.account,
        segment_index=segment_index,
        working_dir=working_dir,
        segment_prod_time=segment_prod_time,
        segment_frames_to_save=segment_frames_to_save
    )

    return job_script






def validate_replicate_range(replicate_range):
    """Validate that replicate_range is in proper SLURM array format."""
    # Check for common SLURM array formats: 1-5, 1,3,5, 1-5:2, etc.
    import re
    pattern = r'^(\d+(-\d+(:\d+)?)?)(,\d+(-\d+(:\d+)?)?)*$'
    if not re.match(pattern, replicate_range):
        raise ValueError(f"Invalid replicate range format: {replicate_range}")
    return True

def main():
    parser = argparse.ArgumentParser(add_help=True, description="Submit daisy-chained MD simulation jobs")

    # Original arguments
    parser.add_argument("-p", "--protein_name", type=str, required=True, help="Name of the protein (str)")
    parser.add_argument("-t", "--temperature", type=float, required=True, help="Temperature in Kelvin (float)")
    parser.add_argument("-e", "--eq_time", type=float, required=True, help="Equilibration time in nanoseconds (float)")
    parser.add_argument("-r", "--prod_time", type=float, required=True, help="Production time in nanoseconds (float)")
    parser.add_argument("-pp", "--path_to_pro_pdb", type=str, required=True, help="Path to the protein PDB file (str)")
    parser.add_argument("-rep", "--replicate_range", type=str, required=True, 
                        help="SLURM array range for replicates (e.g., '1-5', '1,3,5', '1-10:2')")
    
    # New arguments for daisy-chaining
    parser.add_argument("-seg", "--total_segments", type=int, default=10,
                        help="Total number of segments to split the simulation into (default: 10)")
    parser.add_argument("--dry-run", action="store_true", 
                        help="Generate job scripts but don't submit them")

    # HPC partition configuration arguments
    parser.add_argument("--partition", type=str, default="aa100",
                        help="SLURM partition to use (default: aa100)")
    parser.add_argument("--qos", type=str, default="normal",
                        help="SLURM QoS to use (default: normal)")
    parser.add_argument("--account", type=str, default="ucb625_asc1",
                        help="SLURM account to use (default: ucb625_asc1)")
    parser.add_argument("--email", type=str, default="jola3134@colorado.edu",
                        help="Email for job notifications (default: jola3134@colorado.edu)")
    parser.add_argument("--time-limit", type=str, default="23:59:59",
                        help="Time limit for each job (default: 23:59:59)")
    
    # Preset partition configurations
    parser.add_argument("--preset", type=str, choices=["aa100", "blanca-shirts", "testing", "al40"], 
                        help="Use preset partition configuration (overrides individual partition settings)")

    parser.add_argument("--frames", type=int, default=2500,
                        help="Total number of frames to save for the production part of the simulation.")

    args = parser.parse_args()

    # Apply preset configurations if specified
    if args.preset == "aa100":
        args.partition = "aa100"
        args.qos = "normal"
        args.account = "ucb625_asc1"
        args.time_limit = "23:59:59"

    elif args.preset == "al40":
        args.partition = "al40"
        args.qos = "normal"
        args.account = "ucb625_asc1"
        args.time_limit = "23:59:59"

    elif args.preset == "blanca-shirts":
        args.partition = "blanca,blanca-shirts"
        args.qos = "preemptable"
        args.account = "blanca-shirts"
        args.time_limit = "23:59:59"

    elif args.preset == 'testing':
        args.partition= "atesting_a100"
        args.qos = "testing"
        args.account = "ucb625_asc1"
        args.time_limit = "0:59:59"


    def submit_job(script_content, script_name, dependency_job_id=None, output_dir="daisy_chain_scripts"):
        """Submit a job to SLURM and return the job ID."""
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Full path for the script
        script_name = os.path.join(output_dir, script_name)

        # Write the script to a file
        with open(script_name, 'w') as f:
            f.write(script_content)
        
        # Make the script executable
        os.chmod(script_name, 0o755)
        
        # Prepare the sbatch command
        cmd = ['sbatch']
        if dependency_job_id:
            cmd.extend(['--dependency', f'afterok:{dependency_job_id}'])
        if args.preset == "blanca-shirts":
            cmd.extend(['--exclude=bgpu-bortz1'])
        cmd.append(script_name)
        
        # Submit the job
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            # Extract job ID from output (format: "Submitted batch job 12345")
            job_id = result.stdout.strip().split()[-1]
            print(f"Submitted job {job_id} from script {script_name}")
            return job_id
        except subprocess.CalledProcessError as e:
            print(f"Error submitting job: {e}")
            print(f"STDOUT: {e.stdout}")
            print(f"STDERR: {e.stderr}")
            sys.exit(1)


    # Validate replicate range format
    try:
        validate_replicate_range(args.replicate_range)
        replicate_list = parse_replicate_range(args.replicate_range)
    except ValueError as e:
        print(f"Error: {e}")
        print("Valid formats: '1-5', '1,3,5', '1-10:2', etc.")
        sys.exit(1)

    total_jobs = len(replicate_list) * args.total_segments
    
    print(f"Preparing {args.total_segments} segment simulation jobs for {args.protein_name}")
    print(f"Total production time: {args.prod_time} ns")
    print(f"Time per segment: {args.prod_time / args.total_segments} ns")
    print(f"Replicates: {replicate_list} ({len(replicate_list)} total)")
    print(f"Total jobs to submit: {total_jobs}")
    print(f"Dependency chains: {len(replicate_list)} independent chains of {args.total_segments} jobs each")
    print(f"Partition: {args.partition}")
    print(f"QoS: {args.qos}")
    print(f"Account: {args.account}")
    print(f"Time limit: {args.time_limit}")
    print(f"Email: {args.email}")
    print()

    # # Check that continue_simulation.py exists
    # if not Path("new_continue_simulation.py").exists():
    #     print("Error: new_continue_simulation.py not found in current directory")
    #     print("Please ensure new_continue_simulation.py is in your working directory")
    #     sys.exit(1)

    # Check that continue_simulation.py exists
    if not Path("new_continue_simulation_polymerist2.py").exists():
        print("Error: new_continue_simulation_polymerist2.py not found in current directory")
        print("Please ensure new_continue_simulation_polymerist2.py is in your working directory")
        sys.exit(1)


    replicate_job_chains = {}
    
    # Submit jobs for each replicate
    for replicate_num in replicate_list:
        print(f"Submitting job chain for replicate {replicate_num}...")
        replicate_job_chains[replicate_num] = []
        
        # Submit jobs for each segment of this replicate
        for segment_index in range(args.total_segments):

            if segment_index == 0:
                # Initial job
                script_content = create_initial_job_script(args, replicate_num, segment_index, args.total_segments)
                script_name = f"initial_seg{segment_index}_rep{replicate_num}_{args.protein_name}_RBY.sh"
                dependency_job_id = None
            else:
                # Continuation job
                script_content = create_continuation_job_script(args, replicate_num, segment_index, args.total_segments)
                script_name = f"continue_seg{segment_index}_rep{replicate_num}_{args.protein_name}_RBY.sh"
                # Depend on the previous segment of the SAME replicate
                dependency_job_id = replicate_job_chains[replicate_num][segment_index - 1]
            
            if args.dry_run:
                print(f"  Would create script: {script_name}")
                with open(script_name, 'w') as f:
                    f.write(script_content)
                # For dry run, use dummy job ID
                job_id = f"DUMMY_{replicate_num}_{segment_index}"
            else:
                job_id = submit_job(script_content,
                                           script_name,
                                           dependency_job_id)
                print(f"  Submitted segment {segment_index}: job {job_id}")
            
            replicate_job_chains[replicate_num].append(job_id)
    
    if args.dry_run:
        print(f"\nDry run completed. {total_jobs} job scripts created.")
        print("Review the scripts and run without --dry-run to submit them.")
    else:
        print(f"\nAll {total_jobs} jobs submitted successfully!")
        print("\nDependency chains:")
        for replicate_num, job_ids in replicate_job_chains.items():
            print(f"  Replicate {replicate_num}: {' -> '.join(job_ids)}")
        
        print(f"\nMonitor progress with: squeue -u $USER")
        print(f"Check job details with: scontrol show job <job_id>")


if __name__ == "__main__":
    main()
