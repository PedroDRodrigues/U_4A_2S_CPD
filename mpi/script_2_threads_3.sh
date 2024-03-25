#!/usr/bin/env bash

#SBATCH --job-name=life3d-mpi-3-g20-2threads
#SBATCH --output=mpi_%j-mpi-3-g20-2threads.out
#SBATCH --error=mpi_%j-mpi-3-g20-2threads.err
#SBTACH --ntasks=1
#SBTACH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH -x lab2p[1-20]
#SBATCH -C lab4

srun ./life3d-mpi 3 1024 0.4 100
