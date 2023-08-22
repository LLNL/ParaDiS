#!/bin/bash 

#------------------------------------------------------------------------------
# This script drives an anisotropic test on tantalum.
#
# Note - you need to compile ParaDiS with anisotropic forces enabled, qmax=10, 
#        Taylor FMM disabled, and uniform FMM enabled...
#
# DEFS += -DANISOTROPIC    (enabled aniso forces)
# DEFS += -DANISO_QMAX=10  (set qmax=10)
# TAYLOR_FMM_MODE=OFF      (disabled)
# DEFS += -DUNIFORMFMM     (use uniform FMM solution)
#
# The FMM table for these runs were generated as follows...
# ../bin/StanfordFMMsTableGen -cubesize 232 -c11 258.85e9 -c12 169.44e9 -c44 93.27e9 -corder 4 -mattype bcc -outfile ta_aniso.fmm 
#------------------------------------------------------------------------------

if [ -z $1 ] ; then mins=10 ; else mins=$1 ; fi

n=32

exe=../bin/paradis
dat=ta_aniso.data
ctl=ta_aniso_neg.ctrl
log=ta_aniso_neg.log
dir=ta_aniso_neg.rslts

rm -rf $log $dir .DS_Store *.core slurm*.out

sys=`uname -n`

if [ "${sys:0:7}" == "rzmanta"  ] ; then sys='lc-sierra' ; fi
if [ "${sys:0:7}" == "rzansel"  ] ; then sys='lc-sierra' ; fi
if [ "${sys:0:6}" == "lassen"   ] ; then sys='lc-sierra' ; fi
if [ "${sys:0:8}" == "rzhasgpu" ] ; then sys='lc-linux'  ; fi
if [ "${sys:0:7}" == "rztrona"  ] ; then sys='lc-linux'  ; fi
if [ "${sys:0:7}" == "rztopaz"  ] ; then sys='lc-linux'  ; fi
if [ "${sys:0:6}" == "pascal"   ] ; then sys='lc-linux'  ; fi
if [ "${sys:0:6}" == "quartz"   ] ; then sys='lc-linux'  ; fi
if [ "${sys:0:5}" == "syrah"    ] ; then sys='lc-linux'  ; fi
if [ "${sys:0:5}" == "borax"    ] ; then sys='lc-linux'  ; fi
if [ "${sys:0:7}" == "surface"  ] ; then sys='lc-linux'  ; fi
if [ "${sys:0:5}" == "wagon"    ] ; then sys='osx'       ; fi
if [ "${sys:0:8}" == "spurlake" ] ; then sys='osx'       ; fi

echo start job `date`

# serial...
# $exe -d $dat $ctl | tee -a $log 

# parallel on OSX...

if [ $sys == "osx" ] ; then
   mpirun -n $n $exe -d $dat $ctl | tee -a $log 
fi

# parallel on lc linux systems...

if [ $sys == "lc-linux" ] ; then
   salloc -n $n -c 1 -t $mins -p pdebug srun -n $n $exe -d $dat $ctl | tee -a $log 
fi

# parallel on lc sierra systems...

if [ $sys == "lc-sierra" ] ; then
   lalloc 1 -W $mins -q pdebug jsrun -n 1 -r 1 -a $n -c $n $exe -d $dat $ctl | tee -a $log 
fi

echo job ended: `date`
