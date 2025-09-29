#!/bin/bash
#SBATCH --partition=blanca-shirts
#SBATCH --job-name=submit_things
#SBATCH --output=meep.%A_%a.out
#SBATCH --qos=blanca-shirts
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0:3:59
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jola3134@colorado.edu
#SBATCH --account=blanca-shirts

module purge
module load miniforge
mamba activate polymerist-env

python daisy_chain_submission_NOPOLY.py -p LipA -t 363 -e 0.5 -r 1000 -pp NH3_terminal_His_proton_updated.pdb -rep 1-3 -seg 10 --preset blanca-shirts --frames 2500

