#!/usr/bin/env python3
"""
Script to migrate data from old directory naming convention to Signac structure.
Uses Signac project structure to iterate over jobs and copy corresponding old simulation data.

python src/migrate_signac.py /pl/active/shirts_archive/LaforetJoe/Enzyme_Immobilization_2025/Substrate_Restrained_Sims --project-path project_Proteins_and_CoPolymers/ --dry-run

"""

import os
import shutil
from pathlib import Path
import argparse
import logging
import signac

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def construct_old_directory_name(job):
    """
    Construct old directory name from Signac job statepoint data.
    Handles both polymer simulations and control simulations (no polymer).
    
    Args:
        job: Signac job object
        
    Returns:
        str: Old directory name following the convention
    """
    # Extract parameters from job statepoint
    restraint_dist = job.sp.ligand_restraint_dist
    protein = job.sp.protein
    ligand = job.sp.ligand
    temperature = job.sp.temperature
    eq_time = job.sp.eq_time
    eq_ensemble = job.sp.eq_ensemble
    prod_time = job.sp.prod_time
    prod_ensemble = job.sp.prod_ensemble
    replicate_num = job.sp.replicate_num
    
    # Check if this is a control simulation (no polymer)
    is_control = job.sp.SMILES_poly_1 is None
    
    # Convert restraint distance to string with A_RESTRAINT
    restraint_str = f"{restraint_dist:g}A_RESTRAINT"
    
    # Construct temperature string
    temp_str = f"{temperature}K"
    
    # Construct equilibration string
    eq_str = f"{eq_time}ns-{eq_ensemble}"
    
    # Construct production string
    prod_str = f"{prod_time}ns-{prod_ensemble}"
    
    # Construct replicate string
    replicate_str = f"run{replicate_num}"
    
    if is_control:
        # Control simulation directory pattern:
        # "10A_RESTRAINT_LipA_Resorufin-Butyrate_293.0K_0.5ns-NVT_1000.0ns-NPT_run1"
        old_dir_name = f"{restraint_str}_{protein}_{ligand}_{temp_str}_{eq_str}_{prod_str}_{replicate_str}"
    else:
        # Polymer simulation directory pattern:
        # "10A_RESTRAINT_LipA_Resorufin-Butyrate_SBMA-EGPMA-0%_38x_293.0K_0.5ns-NVT_1000.0ns-NPT_run1"
        polymer_ratio = job.sp.polymer_ratio
        num_oligomers = job.sp.num_oligomers
        
        # Convert polymer ratio to percentage
        polymer_percent = f"{polymer_ratio*100:g}%"
        
        # Construct polymer string (assuming SBMA-EGPMA pattern)
        polymer_str = f"SBMA-EGPMA-{polymer_percent}"
        
        # Construct oligomer string
        oligomer_str = f"{num_oligomers}x"
        
        old_dir_name = f"{restraint_str}_{protein}_{ligand}_{polymer_str}_{oligomer_str}_{temp_str}_{eq_str}_{prod_str}_{replicate_str}"
    
    return old_dir_name

def copy_job_data(job, old_data_root, dry_run=False):
    """
    Copy data from old directory structure to Signac job directory.
    
    Args:
        job: Signac job object
        old_data_root (Path): Root directory containing old simulation data
        dry_run (bool): If True, only print what would be copied
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Construct old directory name and path
        old_dir_name = construct_old_directory_name(job)
        old_dir_path = old_data_root / old_dir_name
        
        logger.info(f"\nProcessing job: {job.id}")
        logger.info(f"Looking for old directory: {old_dir_name}")
        
        if not old_dir_path.exists():
            logger.warning(f"Old directory does not exist: {old_dir_path}")
            return False
        
        # Get job directory path
        job_dir = Path(job.path)
        
        # Get list of files/directories to copy from inside the old directory
        items_to_copy = []
        for item in os.listdir(old_dir_path):
            src_item = old_dir_path / item
            dst_item = job_dir / item
            
            # Skip if destination already exists
            if dst_item.exists():
                logger.info(f"Skipping {item} (already exists in job directory)")
                continue
                
            items_to_copy.append((src_item, dst_item, item))
        
        if dry_run:
            logger.info(f"DRY RUN: Would copy {len(items_to_copy)} items to job {job.id}")
            for src_item, dst_item, item in items_to_copy:
                item_type = "directory" if src_item.is_dir() else "file"
                logger.info(f"  Would copy {item_type}: {item}")
            return True
        
        # Actually copy the items from inside the old directory
        for src_item, dst_item, item in items_to_copy:
            if src_item.is_dir():
                shutil.copytree(src_item, dst_item)
                logger.info(f"Copied directory: {item}")
            else:
                shutil.copy2(src_item, dst_item)
                logger.info(f"Copied file: {item}")
        
        logger.info(f"Successfully copied {len(items_to_copy)} items to job {job.id}")
        return True
        
    except Exception as e:
        logger.error(f"Error processing job {job.id}: {e}")
        return False

def migrate_data(*jobs, old_data_root, dry_run=False):
    """
    Migrate data for all specified jobs.
    
    Args:
        *jobs: Variable number of Signac job objects
        old_data_root (str): Root directory containing old simulation data
        dry_run (bool): If True, only show what would be copied
    """
    old_data_path = Path(old_data_root).resolve()
    
    if not old_data_path.exists():
        logger.error(f"Old data root directory does not exist: {old_data_path}")
        return
    
    logger.info(f"Old data root: {old_data_path}")
    logger.info(f"Dry run mode: {dry_run}")
    logger.info(f"Processing {len(jobs)} jobs")
    
    successful_copies = 0
    failed_copies = 0
    
    for job in jobs:
        if copy_job_data(job, old_data_path, dry_run):
            successful_copies += 1
        else:
            failed_copies += 1
    
    logger.info(f"\nMigration complete!")
    logger.info(f"Successful copies: {successful_copies}")
    logger.info(f"Failed copies: {failed_copies}")

def migrate_all_jobs_with_project(project, old_data_root, dry_run=False):
    """
    Migrate data for all jobs in the specified Signac project.
    
    Args:
        project: Signac project object
        old_data_root (str): Root directory containing old simulation data
        dry_run (bool): If True, only show what would be copied
    """
    jobs = list(project)
    
    logger.info(f"Found {len(jobs)} jobs in project")
    migrate_data(*jobs, old_data_root=old_data_root, dry_run=dry_run)

def migrate_all_jobs(old_data_root, dry_run=False):
    """
    Migrate data for all jobs in the current Signac project.
    
    Args:
        old_data_root (str): Root directory containing old simulation data
        dry_run (bool): If True, only show what would be copied
    """
    project = signac.get_project()
    jobs = list(project)
    
    logger.info(f"Found {len(jobs)} jobs in project")
    migrate_data(*jobs, old_data_root=old_data_root, dry_run=dry_run)

def main():
    parser = argparse.ArgumentParser(description="Migrate data from old naming convention to Signac jobs")
    parser.add_argument("old_data_root", help="Root directory containing old simulation data")
    parser.add_argument("--project-path", help="Path to Signac project directory (default: current directory)")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be copied without actually copying")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    parser.add_argument("--jobs", nargs="*", help="Specific job IDs to process (if not specified, processes all jobs)")
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # Get project, optionally from specified path
        if args.project_path:
            project_path = Path(args.project_path).resolve()
            logger.info(f"Looking for Signac project at: {project_path}")
            
            # Check if the project directory exists
            if not project_path.exists():
                logger.error(f"Project directory does not exist: {project_path}")
                return 1
            
            # Check if it contains signac.rc or other project indicators
            if not (project_path / "signac.rc").exists() and not (project_path / "workspace").exists():
                logger.warning(f"Directory may not be a Signac project (no signac.rc or workspace found): {project_path}")
            
            # Change to the project directory temporarily to open the project
            original_cwd = os.getcwd()
            try:
                os.chdir(project_path)
                project = signac.get_project()
                logger.info(f"Successfully opened Signac project at: {project_path}")
            finally:
                os.chdir(original_cwd)
        else:
            project = signac.get_project()
            logger.info(f"Using Signac project at: {os.getcwd()}")
        
        if args.jobs:
            # Process specific jobs
            jobs = [project.open_job(id=job_id) for job_id in args.jobs]
            migrate_data(*jobs, old_data_root=args.old_data_root, dry_run=args.dry_run)
        else:
            # Process all jobs
            migrate_all_jobs_with_project(project, args.old_data_root, args.dry_run)
            
    except Exception as e:
        logger.error(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())