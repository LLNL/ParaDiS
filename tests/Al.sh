#!/bin/sh

rm -rf Al.log Al_results slurm*.out

exe=../bin/paradis
ctl=Al.ctrl
dat=Al.data
log=Al.log

# serial...
# $exe -d $dat $ctl | tee -a $log

# parallel...
n=8
mpirun -n $n $exe -d $dat $ctl | tee -a $log

