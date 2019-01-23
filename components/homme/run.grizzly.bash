#!/bin/bash
#
#   Jobscript for launching dcmip2016 test 1 on the NERSC Edison machine
#
# usage: sbatch jobscript-...

#SBATCH -J h30v200X30ne        # job name
#SBATCH -o bw.o%j  # output and error file name (%j expands to jobID)
#SBATCH -N 30                # total number of mpi tasks requested
#SBATCH --ntasks-per-node 36
##SBATCH -p debug              # queue (partition) -- normal, development, etc.
#SBATCH -t 00:00:30           # run time (hh:mm:ss)
#SBATCH -A w19_teleconnections               # charge hours to account 1
#SBATCH --qos=interactive

source $HOME/homme.bash

date

srun -N 30 ../../../test_execs/theta-l-nlev30/theta-l-nlev30 < input.nl

date
