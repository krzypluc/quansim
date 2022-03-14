#!/bin/bash -l 
#SBATCH -N 1
#SBATCH --ntasks-per-node 32       # Numer of tasks per node
#SBATCH --time=00:05:00            # Required time
#SBATCH -A g86-1055                  # Account
#SBATCH --output=output.out

srun hostname > nodes.txt

javac -d classes -sourcepath ./src/main/java -classpath ./packages/pcj.jar ./src/main/java/compute/E1aClass.java

srun -N 1 -n 1 -c 48 java -cp "classes:./packages/pcj.jar" compute.E1aClass
