#!/usr/bin/env bash

#SBATCH --job-name=life3d-mpi-1000-g20-1node
#SBATCH --output=mpi_%j-mpi-1000-g20-1node.out
#SBATCH --error=mpi_%j-mpi-1000-g20-1node.err
#SBTACH --ntasks=1
#SBTACH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4

srun ./life3d-mpi 1000 64 0.4 0
