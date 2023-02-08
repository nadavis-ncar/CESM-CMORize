#!/bin/bash
#PBS -N CMOR
#PBS -j oe
#PBS -M nadavis@ucar.edu
#PBS -l select=1:ncpus=1
#PBS -A P93300607
#PBS -l walltime=12:00:00
#PBS -q premium

export TMPDIR=/glade/scratch/$USER/tmp
mkdir -p $TMPDIR

module rm ncarenv
module load matlab/R2021b

mkdir -p output

export MPSTASKS=36
export MPSACCOUNT=P93300607
export MPSQUEUE=regular
export MPSWALLTIME="12:00:00"
SECONDS=0

matlab  -nosplash -nodesktop -r "cmor_main"

echo "Time elapsed = $SECONDS s"
