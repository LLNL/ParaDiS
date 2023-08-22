#!/bin/sh

rm -rf Be.log Be_results slurm*.out

exe=../bin/paradis
ctl=Be.ctrl
dat=Be.data
log=Be.log

# serial...
# $exe -d $dat $ctl | tee -a $log

# parallel...
n=8
mpirun -n $n $exe -d $dat $ctl | tee -a $log

