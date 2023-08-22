#!/bin/sh 

nm=ta_bcc_lin_fmm

rm -rf ./$nm.log ./results  slurm*.out

# serial...
#$exe -d $nm.data $nm.ctrl | tee -a $nm.log 

# parallel...
mpirun -n 8 ../../../bin/paradis -d $nm.data $nm.ctrl | tee -a $nm.log 
