#!/bin/bash

#SBATCH --job-name=test_1node
#SBATCH --partition=express
##SBATCH --exclude=""
#SBATCH -N 1
#SBATCH --tasks-per-node=48
#SBATCH --output=equilNpT.out
#SBATCH --error=equilNpT.err
#SBATCH --time=0-01:00:00

module load GCC/11.2.0 OpenMPI/4.1.1 NAMD/2.14-mpi


cd $PWD
charmrun +p $SLURM_NTASKS namd2 equilibrationNpT.conf > equilNpT.log

