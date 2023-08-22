#!/bin/bash 

n=16

if [ -z $1 ] ; then mins=4 ; else mins=$1 ; fi

exe=../bin/paradis
dat=ghost2_fail.data
ctl=ghost2_fail.ctrl
log=ghost2_fail.log
dir=ghost2_fail.rslts

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
   n=32
   salloc -n $n -c 1 -t $mins -p pdebug srun -n $n $exe -maxstep 100 -doms 4 4 2 -cells 8 8 8 -d $dat $ctl | tee -a $log 
fi

# parallel on lc sierra systems...

if [ $sys == "lc-sierra" ] ; then
   n=32
   lalloc 1 -W $mins -q pdebug jsrun -n$n -r$n -a1 -c1 $exe -maxstep 100 -doms 4 4 2 -cells 8 8 8 -d $dat $ctl | tee -a $log 
fi

#----------------------------------------------------------------------
# generate a restart file animation...

(rm -rf $dir/gnu ; mkdir -p $dir/gnu)
(cd $dir/gnu; ../../../bin/paradis_gnu -src ../restart -gnu rs )

echo job ended: `date`
