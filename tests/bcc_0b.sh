#!/bin/bash 

#----------------------------------------------------------------------
# setup...

n=8

if [ -z $1 ] ; then mins=10 ; else mins=$1 ; fi

exe=../bin/paradis
ctl=bcc_0b.ctrl
log=bcc_0b.log
dbg=bcc_0b.lldb
dir=bcc_0b.rslts

# there are three datasets that can be used with this script...

dat=bcc_0b.data
dat=bcc_0b_screws.data
dat=bcc_0b_loops.data

#rm -rf $log $dir .DS_Store *.core slurm*.out

#----------------------------------------------------------------------
# check for restart...

if [ -f $dir/latest_restart ] ; then
   restart=`tail -1 ./$dir/latest_restart`
   ctl=$dir/$restart
   dat=$dir/$restart.data

   if [ $ctl == "$dir/restart/restart.cn" ] ; then
      ctl=$dir/restart/restart.cn
      dat=$dir/restart/restart.data
   fi
fi

#----------------------------------------------------------------------
# identify the system we are on...

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

#----------------------------------------------------------------------
# run paradis...

echo start job `date`

# serial...
# $exe -d $dat $ctl | tee -a $log 

# parallel debug...
# mpirun -n 8 xterm -e lldb $exe -s $dbg

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

#----------------------------------------------------------------------
# generate a restart file animation...

#(rm -rf $dir/gnu ; mkdir -p $dir/gnu)
#(cd $dir/gnu; ../../../bin/paradis_gnu -src ../restart -gnu rs )

echo job ended: `date`

