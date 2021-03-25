#!/bin/bash

module purge
module load hdf5/1.10.6
module load openmpi/4.0.1-cuda
module load cuda10.2/10.2.89
module load gsl/2.6
module list

export MPI_HOME='/cm/shared/apps/openmpi/openmpi-4.0.1.cuda/'
export GRAKLE_HOME='/home/brvillas/code/grackle'
export POISSON_SOLVER='-DSOR'
export SUFFIX='.sor'
make clean
make -j 40

source ~/.bashrc
