#!/bin/bash
#
mpirun -use-hwthread-cpus -np 10 ./cms cots_1999
mpirun -use-hwthread-cpus -np 10 ./cms cots_1998
mpirun -use-hwthread-cpus -np 10 ./cms cots_1997
mpirun -use-hwthread-cpus -np 10 ./cms cots_1996
mpirun -use-hwthread-cpus -np 10 ./cms cots_1995
mpirun -use-hwthread-cpus -np 10 ./cms cots_1994
mpirun -use-hwthread-cpus -np 10 ./cms cots_1993

