#!/bin/csh -x

set paradis=../bin/paradis

set tst=Copper

set dat=$tst.data
set ctl=$tst.ctrl
set dir=$tst.rslts

# check for restart...

if ( -f $dir/latest_restart ) then
   set restart=`tail -1 ./$dir/latest_restart`
   set ctl=$dir/$restart
   set dat=$dir/$restart.data

   if ( $ctl == "$dir/restart/restart.cn" ) then
      set ctl=$dir/restart/restart.cn
      set dat=$dir/restart/restart.data
   endif

endif

# serial...
#$paradis -d $dat $ctl | tee -a $tst.log 

# parallel...
mpirun -n 8 $paradis -d $dat $ctl | tee -a $tst.log 

