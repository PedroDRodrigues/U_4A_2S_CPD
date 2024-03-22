#!/usr/bin/env bash

#SBATCH --job-name=life3d-mpi-10-g20-8node
#SBATCH --output=mpi_%j-mpi-10-g20-8node.out
#SBATCH --error=mpi_%j-mpi-10-g20-8node.err
#SBTACH --ntasks=8
#SBTACH --ntasks-per-node=1
#SBATCH --nodes=8
#SBATCH --cpus-per-task=4

srun ./life3d-mpi 10 512 0.4 0
