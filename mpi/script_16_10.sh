#!/usr/bin/env bash

#SBATCH --job-name=life3d-mpi-10-g20-16node
#SBATCH --output=mpi_%j-mpi-10-g20-16node.out
#SBATCH --error=mpi_%j-mpi-10-g20-16node.err
#SBTACH --ntasks=16
#SBTACH --ntasks-per-node=1
#SBATCH --nodes=16
#SBATCH --cpus-per-task=4

srun ./life3d-mpi 10 512 0.4 0
