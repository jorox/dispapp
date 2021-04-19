#!/usr/bin/bash

module load mpi

WORKDIR=workspace
SRCDIR=src

export PYTHONPATH=$SRCDIR/python:$PYTHONPATH # add mdlib to python modules

# Build the crystal and create map file
python $SRCDIR/python/pylmp pre $WORKDIR/argon_prim_ucell.lmp 5 10 20 $WORKDIR/argon_prim.lmp

# Run sim
for P in 0 -500 -200 -100 -50 -20 -10 10 20 50 100 200 500
do
mpirun -np 4 lmp_g++_openmpi \
-in $SRCDIR/lammps/in.main \
-log log.main \
-var PRESS $P \
-var TEQ 500 \
-var TCOR 100 \
-var NUC 20 \
-var WORKDIR $WORKDIR \
-log log.p$P
done

# Post-process
python $SRCDIR/python/pylmp post $WORKDIR "ArPrim.*.log"
