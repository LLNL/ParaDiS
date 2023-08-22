#!/bin/bash
#MSUB -q pdebug
#MSUB -A bdivp
#MSUB -l nodes=2
#MSUB -l resfailpolicy=ignore # don't kill job allocation upon a node failure
#MSUB -l qos=normal
#MSUB -l gres=ignore
#MSUB -l walltime=00:10:00

#
#   This is a sample script to demonstrate how to run a simulation using
#   the Scalable Checkpoint/Restart (SCR) library.
#
#   Much of the functionality of SCR is controlled by environment
#   variables.  Many of the common ones and their uses are described
#   below.
#

# 
#   Use the dotkit utility to modify the environment with
#   correct paths to the SCR libraries, executables, etc.
#
. /usr/local/tools/dotkit/init.sh; use scr-1.1

#
#   Define the number of nodes and tasks needed for the simualtion.
#
#   Note: If this number of nodes is less than the number of nodes
#   allocated to the job and a node fails, the simulation can be
#   restarted on the remaining nodes available in the job allocation.
#
numNodes=1
numTasks=8


#
#   SCR Environment Variables:
#
#   SCR_CACHE_BASE
#       Specifies the node local path for checkpoints.
#       These checkpoints files will be lost when the parallel job ends.
#       Use the SCR_FLUSH environment variable to periodically flush a
#       local checkpoint to permanent storage.
#
#   SCR_CHECKPOINT_INTERVAL
#       Defines the frequency with which SCR checkpoints are taken.
#       For example, assuming an application checks with the SCR library
#       once per timestep, setting SCR_CHECKPOINT_INTERVAL=10 will
#       cause SCR to write a checkpoint file (to local storage)
#
#   SCR_COPY_TYPE 
#       Indicates node local checkpoint data is stored.
#       Valid options are:
#
#           SCR_COPY_TYPE=LOCAL
#           SCR_COPY_TYPE=PARTNER
#           SCR_COPY_TYPE=XOR
#           SCR_COPY_TYPE=FILE
#
#       See the SCR documentation for details on the behavior of these
#       various settings.
#
#   SCR_DISTRIBUTE
#       Is a toggle enabling/disabling file distribution
#       during SCR_Init().  File distribution enables scr_srun to replace
#       failed nodes with spare nodes during a restart from the checkpoint
#       cache.
#
#   SCR_ENABLE
#       Is a toggle enabling/disabling SCR at run time.
#
#   SCR_FETCH
#       Enables/disables SCR from fetching checkpoint files
#       from the file system during startup.
#
#   SCR_FLUSH
#       Specifies the number of checkpoints written between
#       periodic writes of SCR checkpoints to the file system.  Set
#       to zero disables flushes to the file system. 
#
#   SCR_HALT_SECONDS
#       Set to a positive integer to instruct SCR to halt the job
#       after a successful checkpoint if the remaining time in the
#       job allocation is less than the specified number of seconds.
#
#   SCR_HOP_DISTANCE
#       When SCR_COPY_TYPE is set to PARTNER, each tasks's checkpoint
#       data is replicated on another tasks's local space.  This allows
#       the application to withstand some degree of node failures.  When
#       this environment variable is set, each task's checkpoint data
#       will be replicated on the node owing the task with taskID
#       SCR_HOP_DISTANCE higher (when this would exceed the number of
#       tasks, the numbers wrap back to the first task). This value
#       should be set such that the duplicated data is on a different
#       node.  (i.e. If a simulation is using 8 tasks per node, the
#       hop distance should be set to 8).
#
#   SCR_RUNS
#       Sets the number of times scr_srun will attempt to run the
#       simulation within the job allocation.  Setting this to -1
#       specifies an unlimited number of retries (that is, until
#       the job allocation ends due to a job time limit, the job
#       being explicitly cancelled via scancel, etc).
#
#   SCR_PREFIX
#       Specifes the global path (i.e. permanent storage on the
#       parallel file system) where checkpoint directories should
#       be read from or written to when SCR_FETCH and/or SCR_FLUSH
#       are enabled
#
export SCR_ENABLE=1
export SCR_LOG_ENABLE=0
export SCR_JOB_NAME=paradis
export SCR_DISTRIBUTE=0
export SCR_PREFIX=/p/lscratchrza/hommes/SCR
export SCR_FETCH=1
export SCR_FLUSH=10
export SCR_CACHE_BASE=$TMP
export SCR_COPY_TYPE=LOCAL
#export SCR_HOP_DISTANCE=8
export SCR_RUNS=3
export SCR_CHECKPOINT_INTERVAL=10
export SCR_HALT_SECONDS=300

#
#   Use scr_srun rather than execute 'srun' directly.  The scr_srun
#   is just a wapper for srun that does some more SCR specific setup,
#   handles rerunning the application (if applicable) in the event
#   of node failures, etc.
#
scr_srun -l -N ${numNodes} -n ${numTasks} ./bin/paradis ./tests/Copper.ctrl

