#!/bin/tcsh -f
#PBS -l nodes=4:ppn=8
#PBS -l walltime=8:00:00
#PBS -N swtc2
#PBS -j oe
#PBS -A FY081407
#PBX -a 1758
#XXX -W depend=afterany:jobid

#
#  Shallow water test case 6 Runge-Kutta test
#  used to check stability of m-stage RK schemes
#
set wdir = /lustre/scratch3/turquoise/jrub/acme/lorenz2
set HOMME = /turquoise/usr/projects/climate/jrub/FRONTS/acme/components/homme
set MACH = /turquoise/usr/projects/climate/jrub/FRONTS/acme-balu/grizzly.cmake
#set MACH = $HOMME/cmake/machineFiles/darwin.cmake
#set MACH = $HOMME/cmake/machineFiles/rhel5.cmake
set input = $HOMME/test/sw_conservative

set builddir = $wdir/build
set rundir = $wdir/swtc2-peixoto/NE80-t144-smooth.00
mkdir -p $rundir
mkdir -p $builddir

cd $builddir


set NCPU = 36
if ( ${?PBS_NODEFILE} ) then
   set NCPU = `wc $PBS_NODEFILE | awk '{print $1}' - `
endif
if ( ${?SLURM_NNODES} ) then
   # SLURM_NNODES  = number of nodes
   # hard to tell how many cores per nodes
   # set NCPU to zero, and mpirun will use the max allowed
   set NCPU = 0
endif
echo NCPU = $NCPU

set test_case = swtc2

#configure the model
#balu build=1 leads to cmake failing. So, instead configure and build in usual fashion rather than here.
#balu use this script only to run the SW test case instead
set build = 0
set make = 0
if ( $#argv >= 1) then
  if ( $1 == 'build' ) set build = 1
endif
#cmake:
cd $builddir
if ( $build == 1 ) then
   rm -rf CMakeFiles CMakeCache.txt
   cmake -C $MACH -DSWEQX_PLEV=1  -DSWEQX_NP=4 $HOMME
   exit
endif
if ( $make == 1 ) then
   make -j4 sweqx
    if ($status) exit
endif
set exe = $builddir/src/sweqx/sweqx



cd $rundir
mkdir -p movies


# defaults:
set smooth=0
set rk_stage=0
set nu = 0
set nu_s = -1  # defaults to nu
set LFTfreq = 0

set hypervis_subcycle =  1
set integration = explicit

set NE = 90
#set nu = 5.1e15


### leapfrog
<<<<<<< HEAD
set smooth = 0.05 ; set LFTfreq = 0
set ndays = 36 #12
set output_end_time  = 36
set tstep = 144
=======
# balu set smooth = 0.05 ; set LFTfreq = 0
# balu set tstep = 144
>>>>>>> nonhydro/homme/lorenz2-sw-peixoto


### leapfrog-trapazoidal 
set smooth = 0 ; set LFTfreq = 1
set tstep = 144

### RK2-M (what's used in 3D code)
#set smooth = 0 ; set LFTfreq = 4
#set tstep = 144

### RK2-SSP 
#set integration = runge_kutta
#set rk_stage = 2
#set tstep = 75

set filter_freq = 0
set name = ${test_case}-NE${NE}-t${tstep}


set sfreq = 6
@ sfreq *= 3600
set sfreq = `echo "$sfreq / $tstep" | bc`


sed s/ne=.\*/"ne = $NE"/  $input/swtc2.nl |\
sed s/ndays.\*/"ndays = $ndays"/  |\
sed s/output_end_time.\*/"output_end_time = $output_end_time"/  |\
sed s/tstep.\*/"tstep = $tstep"/  |\
sed s/smooth.\*/"smooth = $smooth"/  |\
sed s/test_case.\*/"test_case = \'$test_case\'"/  |\
sed s/integration.\*/"integration = '$integration'"/  |\
sed s/rk_stage_user.\*/"rk_stage_user = $rk_stage"/  |\
sed s/LFTfreq.\*/"LFTfreq = $LFTfreq"/  |\
sed s/nu=.\*/"nu= $nu"/  |\
sed s/nu_s=.\*/"nu_s= $nu_s"/  |\
sed s/filter_freq.\*/"filter_freq = $filter_freq"/  |\
sed s/hypervis_subcycle.\*/"hypervis_subcycle = $hypervis_subcycle"/  |\
sed s/statefreq.\*/"statefreq = $sfreq"/  \
> input.nl

date
mpirun -np $NCPU $exe < input.nl | tee  sweq.out
date

mv -f sweq.mass $name.mass
mv -f sweq.out $name.out

mv -f swtc2.l1.errors  $name.l1.errors
mv -f swtc2.l2.errors  $name.l2.errors
mv -f swtc2.linf.errors  $name.linf.errors

echo "day,l1 error: "  `tail -1 $name.l1.errors`
echo "day,l2 error: "  `tail -1 $name.l2.errors`
echo "day,l8 error: "  `tail -1 $name.linf.errors`


