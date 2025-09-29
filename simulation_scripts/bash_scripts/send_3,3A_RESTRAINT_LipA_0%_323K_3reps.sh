#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --job-name=submit_things
#SBATCH --output=meep.%A_%a.out
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=23:59:59
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jola3134@colorado.edu
#SBATCH --account=ucb625_asc1

module purge
module load miniforge
mamba activate polymerist-env

python daisy_chain_submission.py -p LipA -o EGPMA-SBMA -n 38 -t 323 -e 0.5 -r 1000 -pp NH3_terminal_His_proton_updated.pdb -ps ATRP_EGPMA_SBMA_5-mer -rep 1-3 -chars A B -probs 0.00 1.00 -plen 5 -prefix EGPMA-SBMA -seg 10 --preset al40 --frames 2500

