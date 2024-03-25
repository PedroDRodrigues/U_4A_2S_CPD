#!/usr/bin/env bash

#SBATCH --job-name=life3d-mpi-1000-g20-serial
#SBATCH --output=mpi_%j-mpi-1000-g20-serial.out
#SBATCH --error=mpi_%j-mpi-1000-g20-serial.err
#SBTACH --ntasks=1
#SBTACH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -x lab2p[1-20]
#SBATCH -C lab4

srun ./life3d-mpi 1000 64 0.4 0
