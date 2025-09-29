#!/bin/bash
#SBATCH --partition=blanca-shirts
#SBATCH --job-name=submit_things
#SBATCH --output=meep.%A_%a.out
#SBATCH --qos=blanca-shirts
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=23:59:59
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jola3134@colorado.edu
#SBATCH --account=blanca-shirts

module purge
module load miniforge
mamba activate polymerist-env

python daisy_chain_submission.py -p LipA -o EGPMA-SBMA -n 38 -t 363 -e 0.5 -r 1000 -pp NH3_terminal_His_proton_updated.pdb -ps ATRP_EGPMA_SBMA_5-mer -rep 1-3 -chars A B -probs 0.10 0.90 -plen 5 -prefix EGPMA-SBMA -seg 10 --preset blanca-shirts --frames 2500

