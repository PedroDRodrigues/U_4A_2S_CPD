#!/bin/bash

#SBATCH --job-name=life3d-mpi
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBTACH --ntasks=16
#SBTACH --ntasks-per-node=4
#SBATCH --nodes=4

ARGS=(10, 512, 0.4, 0) 

srun ./life3d-mpi ${ARGS[@]}