#!/bin/sh 

nm=ta_pos_tbl_2

rm -rf ./$nm.log ./results  slurm*.out

# serial...
#$exe -d $nm.data $nm.ctrl | tee -a $nm.log 

# parallel...
mpirun -n 8 ../../../bin/paradis -d $nm.data $nm.ctrl | tee -a $nm.log 

