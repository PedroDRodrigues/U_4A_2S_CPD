#!/usr/bin/env bash

#SBATCH --job-name=life3d-mpi-1000-g20-2node
#SBATCH --output=mpi_%j-mpi-1000-g20-2node.out
#SBATCH --error=mpi_%j-mpi-1000-g20-2node.err
#SBTACH --ntasks=1
#SBTACH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -x lab2p[1-20]
ARGS=(3, 1024, 0.4, 1000)

srun ./life3d-mpi ${ARGS[@]}
