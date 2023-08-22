#!/bin/bash 

n=8

if [ -z $1 ] ; then mins=4 ; else mins=$1 ; fi

paradis=../bin/paradis
dat=bcc_screws.data
ctl=bcc_screws.ctrl
log=bcc_screws.log
dir=bcc_screws.rslts

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

# check for restart...

if [ -f $dir/latest_restart ] ; then
   set restart=`tail -1 ./$dir/latest_restart`
   set ctl=$dir/$restart
   set dat=$dir/$restart.data

   if [ $ctl == "$dir/restart/restart.cn" ] ; then
      set ctl=$dir/restart/restart.cn
      set dat=$dir/restart/restart.data
   fi

fi

echo start job `date`

# serial...
# $paradis -d $dat $ctl | tee -a $log 

# parallel on OSX...

if [ $sys == "osx" ] ; then
   mpirun -n $n $paradis -d $dat $ctl | tee -a $log 
fi

# parallel on lc linux systems...

if [ $sys == "lc-linux" ] ; then
   salloc -n $n -c 1 -t $mins -p pdebug srun -n $n $paradis -d $dat $ctl | tee -a $log 
fi

# parallel on lc sierra systems...

if [ $sys == "lc-sierra" ] ; then
   lalloc 1 -W $mins -q pdebug jsrun -n 1 -r 1 -a $n -c $n $paradis -d $dat $ctl | tee -a $log 
fi

echo job ended: `date`
