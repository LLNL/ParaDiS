#!/bin/sh

rm -rf *.log *_results .DS_Store *.core

tst=ThreeSpheres

exe=../bin/paradis
ctl=./$tst.ctrl
dat=./$tst.data
log=./$tst.log

# serial...

#$exe -d $dat $ctl | tee -a $log 

# parallel...

n=8
mpirun -n $n $exe -d $dat $ctl | tee -a $log 
