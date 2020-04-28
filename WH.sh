#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --account=blanca-appm-student
#SBATCH --job-name=WH
#SBATCH --nodes=1
#SBATCH --ntasks=36
#SBATCH --output=WH.out

module purge
module load gcc/8.2.0
module load openmpi_ucx/4.0.0

#mpicxx RMA_test.cpp
mpirun ./main
