#!/usr/bin/env bash

#SBATCH --job-name=life3d-mpi-200-g20-serial
#SBATCH --output=mpi_%j-mpi-200-g20-serial.out
#SBATCH --error=mpi_%j-mpi-200-g20-serial.err
#SBTACH --ntasks=1
#SBTACH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

srun ./life3d-mpi 200 128 0.5 1000
