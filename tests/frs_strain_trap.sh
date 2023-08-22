#!/bin/bash

if [ -z $1 ] ; then tmin=4 ; else tmin=$1 ; fi

exe=../bin/paradis
dat=frs_strain.data
ctl=frs_strain_trap.ctrl
log=frs_strain_trap.log
dir=frs_strain_trap.rslts

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

# set the domain and cell configurations...

ndx=2;                  # domains(x)
ndy=2;                  # domains(y)
ndz=2;                  # domains(z)
n=$((ndx*ndy*ndz))      # total domains(x,y,z)

ncx=8;                  # cells(x)
ncy=8;                  # cells(y)
ncz=8;                  # cells(z)
maxstep=1000            # maximum number of steps to run

# serial on OSX...
# $exe -d $dat $ctl | tee -a $log

# parallel on OSX...

if [ $sys == "osx" ] ; then
   mpirun -n $n $exe -maxstep $maxstep -doms $ndx $ndy $ndz -cells $ncx $ncy $ncz -d $dat $ctl | tee -a $log
fi

# on larger machines, use 4x4x2 (32) domains...

ndx=4;
ndy=4;
ndz=2; 
n=$((ndx*ndy*ndz))

# parallel on lc linux systems...

if [ $sys == "lc-linux" ] ; then
   salloc -n $n -c 1 -t $tmin -p pdebug srun -n $n $exe -maxstep $maxstep -doms $ndx $ndy $ndz -cells $ncx $ncy $ncz -d $dat $ctl | tee -a $log
fi

# parallel on lc sierra systems...

if [ $sys == "lc-sierra" ] ; then
   lalloc 1 -W $tmin -q pdebug jsrun -n$n -r$n -a1 -c1 $exe -maxstep $maxstep -doms $ndx $ndy $ndz -cells $ncx $ncy $ncz -d $dat $ctl | tee -a $log
fi

#----------------------------------------------------------------------
# generate a restart file animation...

#(rm -rf $dir/gnu ; mkdir -p $dir/gnu)
#(cd $dir/gnu; ../../../bin/paradis_gnu -src ../restart -gnu rs )

echo job ended: `date`
