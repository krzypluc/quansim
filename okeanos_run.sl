#!/bin/bash -l 
#SBATCH -N 1
#SBATCH --ntasks-per-node 32       # Numer of tasks per node
#SBATCH --time=00:05:00            # Required time
#SBATCH -A g86-1055                  # Account
#SBATCH --output=output.out

srun hostname > nodes.txt

srun -N 1 -n 1 -c 48 java -jar accquansim-0.1-beta-shaded.jar config/config.yml
