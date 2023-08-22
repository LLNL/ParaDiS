#!/bin/sh

rm -rf Tantalum.log Tantalum_results slurm*.out

exe=../bin/paradis
ctl=Tantalum.ctrl
dat=Tantalum.data
log=Tantalum.log

# serial...
# $exe -d $dat $ctl | tee -a $log

# parallel...
n=8
mpirun -n $n $exe -d $dat $ctl | tee -a $log

