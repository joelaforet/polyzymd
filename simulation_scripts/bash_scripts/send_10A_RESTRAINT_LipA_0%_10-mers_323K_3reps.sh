#!/bin/bash
#SBATCH --partition=blanca,blanca-shirts
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

python daisy_chain_submission.py -p LipA -o SBMA-EGPMA -n 19 -t 323 -e 0.5 -r 1000 -pp NH3_terminal_His_proton_updated.pdb -ps ATRP_SBMA_EGPMA_10-mer -rep 1-3 -chars A B -probs 1.00 0.00 -plen 10 -prefix SBMA-EGPMA -seg 10 --preset blanca-shirts --frames 2500

