#!/usr/bin/bash

module load mpi

# Build the crystal using conventional or primitive axes
#  Create the map file for fix-phonon

mpirun -np 4 lmp_g++_openmpi \
-in in.main \
-log log.main \
-var TSIM 350 \
-var TEQ 200 \
-var NUC 20