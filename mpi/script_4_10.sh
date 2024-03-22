#!/usr/bin/env bash

#SBATCH --job-name=life3d-mpi-10-g20-4node
#SBATCH --output=mpi_%j-mpi-10-g20-4node.out
#SBATCH --error=mpi_%j-mpi-10-g20-4node.err
#SBTACH --ntasks=4
#SBTACH --ntasks-per-node=1
#SBATCH --nodes=4
#SBATCH --cpus-per-task=4

srun ./life3d-mpi 10 512 0.4 0
