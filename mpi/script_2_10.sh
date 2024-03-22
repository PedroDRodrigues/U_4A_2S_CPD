#!/usr/bin/env bash

#SBATCH --job-name=life3d-mpi-10-g20-2node
#SBATCH --output=mpi_%j-mpi-10-g20-2node.out
#SBATCH --error=mpi_%j-mpi-10-g20-2node.err
#SBTACH --ntasks=2
#SBTACH --ntasks-per-node=1
#SBATCH --nodes=2
#SBATCH --cpus-per-task=4

srun ./life3d-mpi 10 512 0.4 0
