#!/usr/bin/env bash

#SBATCH --job-name=life3d-mpi-200-g20-2node
#SBATCH --output=mpi_%j-mpi-200-g20-2node.out
#SBATCH --error=mpi_%j-mpi-200-g20-2node.err
#SBTACH --ntasks=2
#SBTACH --ntasks-per-node=1
#SBATCH --nodes=2
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH -x lab2p[1-20]

srun ./life3d-mpi 200 128 0.5 1000
