#!/bin/bash

# Name job 'poisson'
#PBS -N poisson

# Allocate two nodes with 12 processors from the default resources
#PBS -lnodes=2:ppn=12:default

# Expect to run up to 5 minutes
#PBS -lwalltime=00:45:00

# Memory per process
#PBS -lpmem=2000MB

# Run on the freecycle account
#PBS -A freecycle

# Run in the optimist queue by default
#PBS -q optimist

#Join stdout and stderr output to one file
#PBS -j oe

# Change directory to dir with the job script
cd ${PBS_O_WORKDIR}

# Load needed modules
module load intelcomp
module load openmpi/1.4.3-intel

# Set thread affinity
KMP_AFFINITY="granularity=fine,compact"

# Run with 8 MPI processes, each with 3 threads

N=1024
while [ $N -lt 16385 ]; do
    echo .
    echo Running sequential solver with n = $N
    time ./Release/poisson-single $N
    let N=N*2
done

N=1024
while [ $N -lt 16385 ]; do
    echo .
    echo Running MPI solver with n = $N, nodes = 2 and procs per node = 12
    OMP_NUM_THREADS=1 mpirun -npernode 12 ./Release/poisson $N
    let N=N*2
done

N=1024
while [ $N -lt 16385 ]; do
    echo .
    echo Running MPI solver with n = $N, nodes = 2 and procs per node = 2
    OMP_NUM_THREADS=6 mpirun -npernode 2 ./Release/poisson $N
    let N=N*2
done
