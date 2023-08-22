#!/bin/bash

# The new Sierra clusters at LLNL do not have GPUs installed on the 
# interactive nodes that used to build and submit jobs.  To get access 
# to a compute node with a GPU, you much submit a batch job or request 
# interactive access to a compute node.  This script was composed to 
# assist with spawning interactive sessions.
#
# Note that The LSF system is deployed on the early access systems for 
# Sierra. LSF uses the bsub command to send jobs to the queues. 
# Note - the current early access systems do not have user accounts
# tied to banks. We use the 'guests' account for submitting jobs.

mins=30          # set default to 30 minutes

if [ $1 ]; then  # override default with cmd line arg 1
   mins=$1       #
fi               #

# bsub argument notes...
# -x            : exclusive access to node
# -n 20         : allocate 20 tasks (should be coupled with an mpirun command)
# -Is           : interactive session
# -G guests     : use guest account resources
# -W $mins      : time limit in minutes (or hh:mm)
# -XF           : uses SSH X11 forwarding
# xterm -sb     : spawns an xterm
# /usr/bin/bash : spawns a bash session in the current terminal

# some common invocations...

#cmd="bsub -x -n 20     -XF -W $mins -G guests xterm -sb"
#cmd="bsub -x -n 20 -Is -XF -W $mins -G guests /usr/bin/bash"
#cmd="bsub          -Is -XF -W $mins -G guests /usr/bin/bash"
 cmd="bsub          -Is -XF -W $mins -G guests /usr/bin/bash"

# launch the session...

echo $cmd; $cmd
