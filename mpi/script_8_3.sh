#!/usr/bin/env bash

#SBATCH --job-name=life3d-mpi-3-g20-8node
#SBATCH --output=mpi_%j-mpi-3-g20-8node.out
#SBATCH --error=mpi_%j-mpi-3-g20-8node.err
#SBTACH --ntasks=8
#SBTACH --ntasks-per-node=1
#SBATCH --nodes=8
#SBATCH --cpus-per-task=4

srun ./life3d-mpi 3 1024 0.4 100
