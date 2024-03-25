#!/usr/bin/env bash

#SBATCH --job-name=life3d-mpi-200-g20-16node
#SBATCH --output=mpi_%j-mpi-200-g20-16node.out
#SBATCH --error=mpi_%j-mpi-200-g20-16node.err
#SBTACH --ntasks=16
#SBTACH --ntasks-per-node=1
#SBATCH --nodes=16
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH -x lab2p[1-20]

srun ./life3d-mpi 200 128 0.5 1000
