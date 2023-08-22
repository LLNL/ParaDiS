/*-------------------------------------------------------------------------
 *
 *      Module:      GenerateOutput.c
 *      Description: 
 *
 *      Includes public functions:
 *
 *          DoParallelIO()
 *          GenerateOutput()
 *          GetOutputTypes()
 *          GetParIOGroup()
 *
 *      Includes private functions:
 *          FreeBinFileArrays()
 *          SetBinFileArrays()
 *
 *------------------------------------------------------------------------*/
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Node.h"
#include "Comm.h"
#include "Util.h"
#include "WriteProp.h"
#include "DisplayC.h"
#include "Restart.h"
#include "UTC_Time.h"
#include "DeltaTime.h"
#include "Narms_Diag.h"
#include "Malloc.h"

/*-------------------------------------------------------------------------
 *
 *      Function:    SendWriteToken
 *      Description: Send the write token to the specified task
 *
 *      Arguments:
 *          taskRank Rank of the task (in MPI_COMM_WORLD) to which
 *                   the write token is to be sent.
 *
 *------------------------------------------------------------------------*/
static void SendWriteToken(int taskRank)
{
#ifdef PARALLEL
        int writeToken = 0;

        MPI_Send((int *) &writeToken, 1, MPI_INT, taskRank, MSG_TOKEN_RING, MPI_COMM_WORLD);
#endif
        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    RecvWriteToken
 *      Description: Wait for the write token from the specified task
 *
 *      Arguments:
 *          taskRank Rank of the task (in MPI_COMM_WORLD) from which
 *                   the write token is expected.
 *
 *------------------------------------------------------------------------*/
static void RecvWriteToken(int taskRank)
{
#ifdef PARALLEL
        int        writeToken=0;
        MPI_Status reqStatus;

        MPI_Recv((int *) &writeToken, 1, MPI_INT, taskRank, MSG_TOKEN_RING, MPI_COMM_WORLD, &reqStatus);
#endif
        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GetOutputTypes
 *      Description: Given the stage of program execution, determine
 *                   which types of output need to be generated and
 *                   increment any applicable file sequence counters for 
 *                   those output types. -- Note the counters are not to
 *                   be updated while dumping data during program termination.
 *      Args:
 *          stage        indicates the program execution stage which
 *                       determines the types of output this function
 *                       might produce.  Valid values for this field are:
 *
 *                               STAGE_INIT
 *                               STAGE_CYCLE
 *                               STAGE_TERM
 *
 *          outputTypes  pointer to an integer which on exit from this
 *                       subroutine will contain a bitmask with 1 bit
 *                       set for each type of output to be generated.
 *
 *------------------------------------------------------------------------*/
static void GetOutputTypes(Home_t *home, int stage, int *outputTypes)
{

#ifndef NO_XWINDOW
        static int xWinAlive = 1;
#endif

        Param_t *param    = home->param;
        real8    timeNow  = param->timeNow;
        real8    dumpTime = 0.0;

/*
 *      If X-window plotting is enabled (and active) and the code
 *      is not in the process of terminating, set the flag to do
 *      the plotting.
 */
#ifndef NO_XWINDOW
        if ((stage == STAGE_INIT) || (stage == STAGE_CYCLE)) 
        {
            if (xWinAlive) 
            {
                if (home->myDomain == 0) 
                {
                    xWinAlive = WinAlive();
                }
#ifdef PARALLEL
                MPI_Bcast((int *) &xWinAlive, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
                if (xWinAlive)
                {
                    *outputTypes |= GEN_XWIN_DATA;
                }
            }
        }
#endif

/*
 *      If restart file dumps are enabled, do so if the code is
 *      in the termination stage, or it is at the end of a cycle
 *      and either the elapsed simulation time since the last dump
 *      has exceeded the allowable delta time between restart dumps
 *      OR the current cycle is a multiple of the specified restart
 *      dump frequency.
 *
 *      NOTE: The last two conditions are mutually exclusive. If
 *            a delta time has been specified to determine the 
 *            time to dump files, no check of the current
 *            cycle against the dump frequency will be done.
 */
        if (param->savecn) 
        {
            if (stage == STAGE_TERM) {
                *outputTypes |= GEN_RESTART_DATA;
            } else if (stage == STAGE_CYCLE) {
                if (param->savecndt > 0.0) {
                    dumpTime = param->savecntime + param->savecndt;
                    if (timeNow >= dumpTime) {
                        param->savecntime = timeNow;
                        param->savecncounter++;
                        *outputTypes |= GEN_RESTART_DATA;
                    }
                } else if ((home->cycle % param->savecnfreq) == 0) {
                    param->savecncounter++;
                    *outputTypes |= GEN_RESTART_DATA;
                }
            }
        }

#ifdef USE_SCR
/*
 *     Check if it is time to dump a checkpoint/restart file
 *     via the SCR library.  The SCR_Need_checkpoint() function
 *     returns SCR_FAILURE if SCR is not enabled or initialized
 *     and returns SCR_SUCCESS if a checkpoint file shuold be
 *     dumped.  Environment variables SCR_CHECKPOINT_INTERVAL,
 *     SCR_CHECKPOINT_SECONDS, and SCR_CHECKPOINT_OVERHEAD
 *     are used in SCR_Need_checkpoint() to make the determination.
 */
        if (stage == STAGE_TERM) 
        {
            int flag=0;
            if (SCR_Need_checkpoint(&flag) == SCR_SUCCESS) 
            {
                *outputTypes |= GEN_SCR_DATA;
            }
        } 
        else if (stage == STAGE_CYCLE) 
        {
            int flag=0;
            if (SCR_Need_checkpoint(&flag) == SCR_SUCCESS) 
            {
                if (flag) *outputTypes |= GEN_SCR_DATA;
            }
        }
#endif /* USE_SCR */
/*
 *      If GNUPLOT file dumps are enabled, do so if the code is
 *      in the termination stage, or at the end of a cycle and
 *      either the elapsed simulation time since the last dump
 *      has exceeded the allowable delta time between file dumps
 *      OR the current cycle is a multiple of the specified dump
 *      frequency.
 *
 *      NOTE: The last two conditions are mutually exclusive. If
 *            a delta time has been specified to determine the 
 *            time to dump files, no check of the current
 *            cycle against the dump frequency will be done.
 */
        if (param->gnuplot) {
            if (stage == STAGE_TERM) {
                *outputTypes |= GEN_GNUPLOT_DATA;
            } else if (stage == STAGE_CYCLE) {
                if (param->gnuplotdt > 0.0) {
                    dumpTime = param->gnuplottime + param->gnuplotdt;
                    if (timeNow >= dumpTime) {
                        param->gnuplottime = timeNow;
                        param->gnuplotcounter++;
                        *outputTypes |= GEN_GNUPLOT_DATA;
                    }
                } else if ((home->cycle % param->gnuplotfreq) == 0) {
                    param->gnuplotcounter++;
                    *outputTypes |= GEN_GNUPLOT_DATA;
                }
            }
        }

/*
 *      If TECPLOT file dumps are enabled, do so if the code is
 *      in the termination stage, or it is at the end of a cycle
 *      and either the elapsed simulation time since the last dump
 *      has exceeded the allowable delta time between file dumps
 *      OR the current cycle is a multiple of the specified 
 *      dump frequency.
 *
 *      NOTE: The last two conditions are mutually exclusive. If
 *            a delta time has been specified to determine the 
 *            time to dump files, no check of the current
 *            cycle against the dump frequency will be done.
 */
        if (param->tecplot) {
            if (stage == STAGE_TERM) {
                *outputTypes |= GEN_TECPLOT_DATA;
            } else if (stage == STAGE_CYCLE) {
                if (param->tecplotdt > 0.0) {
                    dumpTime = param->tecplottime + param->tecplotdt;
                    if (timeNow >= dumpTime) {
                        param->tecplottime = timeNow;
                        param->tecplotcounter++;
                        *outputTypes |= GEN_TECPLOT_DATA;
                    }
                } else if ((home->cycle % param->tecplotfreq) == 0) {
                    param->tecplotcounter++;
                    *outputTypes |= GEN_TECPLOT_DATA;
                }
            }
        }

/*
 *      If velocity file dumps are enabled, do so if the code is
 *      in the termination stage, or it is at the end of a cycle
 *      and either the elapsed simulation time since the last dump
 *      has exceeded the allowable delta time between file dumps
 *      OR the current cycle is a multiple of the specified 
 *      dump frequency.
 *
 *      NOTE: The last two conditions are mutually exclusive. If
 *            a delta time has been specified to determine the 
 *            time to dump files, no check of the current
 *            cycle against the dump frequency will be done.
 */
        if (param->velfile) {
            if (stage == STAGE_TERM) {
                *outputTypes |= GEN_VELOCITY_DATA;
            } else if (stage == STAGE_CYCLE) {
                if (param->velfiledt > 0.0) {
                    dumpTime = param->velfiletime + param->velfiledt;
                    if (timeNow >= dumpTime) {
                        param->velfiletime = timeNow;
                        param->velfilecounter++;
                        *outputTypes |= GEN_VELOCITY_DATA;
                    }
                } else if ((home->cycle % param->velfilefreq) == 0) {
                    param->velfilecounter++;
                    *outputTypes |= GEN_VELOCITY_DATA;
                }
            }
        }

/*
 *      If force file dumps are enabled, do so if the code is
 *      in the termination stage, or it is at the end of a cycle
 *      and either the elapsed simulation time since the last dump
 *      has exceeded the allowable delta time between file dumps
 *      OR the current cycle is a multiple of the specified 
 *      dump frequency.
 *
 *      NOTE: The last two conditions are mutually exclusive. If
 *            a delta time has been specified to determine the 
 *            time to dump files, no check of the current
 *            cycle against the dump frequency will be done.
 */
        if (param->writeForce) {
            if (stage == STAGE_TERM) {
                *outputTypes |= GEN_FORCE_DATA;
            } else if (stage == STAGE_CYCLE) {
                if (param->writeForceDT > 0.0) {
                    dumpTime = param->writeForceTime + param->writeForceDT;
                    if (timeNow >= dumpTime) {
                        param->writeForceTime = timeNow;
                        param->writeForceCounter++;
                        *outputTypes |= GEN_FORCE_DATA;
                    }
                } else if ((home->cycle % param->writeForceFreq) == 0) {
                    param->writeForceCounter++;
                    *outputTypes |= GEN_FORCE_DATA;
                }
            }
        }

/*
 *      If arm file dumps are enabled, do so if the code is
 *      in the termination stage, or it is at the end of a cycle
 *      and either the elapsed simulation time since the last dump
 *      has exceeded the allowable delta time between file dumps
 *      OR the current cycle is a multiple of the specified 
 *      dump frequency.
 *
 *      NOTE: The last two conditions are mutually exclusive. If
 *            a delta time has been specified to determine the 
 *            time to dump files, no check of the current
 *            cycle against the dump frequency will be done.
 */
        if (param->armfile) {
            if (stage == STAGE_TERM) {
                *outputTypes |= GEN_ARM_DATA;
            } else if (stage == STAGE_CYCLE) {
                if (param->armfiledt > 0.0) {
                    dumpTime = param->armfiletime + param->armfiledt;
                    if (timeNow >= dumpTime) {
                        param->armfiletime = timeNow;
                        param->armfilecounter++;
                        *outputTypes |= GEN_ARM_DATA;
                    }
                } else if ((home->cycle % param->armfilefreq) == 0) {
                    param->armfilecounter++;
                    *outputTypes |= GEN_ARM_DATA;
                }
            }
        }

/*
 *      If POLEFIG file dumps are enabled, do so if the code is
 *      in the termination stage, or it is at the end of a cycle
 *      and either the elapsed simulation time since the last dump
 *      has exceeded the allowable delta time between file dumps
 *      OR the current cycle is a multiple of the specified 
 *      dump frequency.
 *
 *      NOTE: The last two conditions are mutually exclusive. If
 *            a delta time has been specified to determine the 
 *            time to dump files, no check of the current
 *            cycle against the dump frequency will be done.
 */
        if (param->polefigfile) {
            if (stage == STAGE_TERM) {
                *outputTypes |= GEN_POLEFIG_DATA;
            } else if (stage == STAGE_CYCLE) {
                if (param->polefigdt > 0.0) {
                    dumpTime = param->polefigtime + param->polefigdt;
                    if (timeNow >= dumpTime) {
                        param->polefigtime = timeNow;
                        param->polefigcounter++;
                        *outputTypes |= GEN_POLEFIG_DATA;
                    }
                } else if ((home->cycle % param->polefigfreq) == 0) {
                    param->polefigcounter++;
                    *outputTypes |= GEN_POLEFIG_DATA;
                }
            }
        }
    
/*
 *      If POVRAY file dumps are enabled, do so if the code is
 *      in the termination stage, or it is at the end of a cycle
 *      and either the elapsed simulation time since the last dump
 *      has exceeded the allowable delta time between file dumps
 *      OR the current cycle is a multiple of the specified 
 *      dump frequency.
 *
 *      NOTE: The last two conditions are mutually exclusive. If
 *            a delta time has been specified to determine the 
 *            time to dump files, no check of the current
 *            cycle against the dump frequency will be done.
 */
        if (param->povray) {
            if (stage == STAGE_TERM) {
                *outputTypes |= GEN_POVRAY_DATA;
            } else if (stage == STAGE_CYCLE) {
                if (param->povraydt > 0.0) {
                    dumpTime = param->povraytime + param->povraydt;
                    if (timeNow >= dumpTime) {
                        param->povraytime = timeNow;
                        param->povraycounter++;
                        *outputTypes |= GEN_POVRAY_DATA;
                    }
                } else if ((home->cycle % param->povrayfreq) == 0) {
                    param->povraycounter++;
                    *outputTypes |= GEN_POVRAY_DATA;
                }
            }
        }

/*
 *      If VisIt file dumps are enabled, do so if the code is
 *      in the termination stage, or it is at the end of a cycle
 *      and either the elapsed simulation time since the last dump
 *      has exceeded the allowable delta time between file dumps
 *      OR the current cycle is a multiple of the specified 
 *      dump frequency.
 *
 *      NOTE: The last two conditions are mutually exclusive. If
 *            a delta time has been specified to determine the 
 *            time to dump files, no check of the current
 *            cycle against the dump frequency will be done.
 */
        if (param->writeVisit) {
            if (stage == STAGE_TERM) {
                *outputTypes |= GEN_VISIT_DATA;
            } else if (stage == STAGE_CYCLE) {
                if (param->writeVisitDT > 0.0) {
                    dumpTime = param->writeVisitTime + param->writeVisitDT;
                    if (timeNow >= dumpTime) {
                        param->writeVisitTime = timeNow;
                        param->writeVisitCounter++;
                        *outputTypes |= GEN_VISIT_DATA;
                    }
                } else if ((home->cycle % param->writeVisitFreq) == 0) {
                    param->writeVisitCounter++;
                    *outputTypes |= GEN_VISIT_DATA;
                }
            }
        }

/*
 *      Check if creation of the density field file is enabled.  Only
 *      done at code termination time...
 */
        if ((stage == STAGE_TERM) &&
            ((param->savedensityspec[0] > 0) &&
             (param->savedensityspec[1] > 0) &&
             (param->savedensityspec[2] > 0))) {
            *outputTypes |= GEN_DENSITY_DATA;
        }

/*
 *      If timer file dumps are enabled, do so if the code
 *      is in the termination stage, or it is at the end
 *      of a cycle and either the elapsed simulation time
 *      since the last dump has exceeded the allowable delta time
 *      between file dumps OR the current cycle is a multiple of
 *      the specified dump frequency.
 *
 *      NOTE: The last two conditions are mutually exclusive. If
 *            a delta time has been specified to determine the 
 *            time to dump files, no check of the current
 *            cycle against the dump frequency will be done.
 */
        if (param->savetimers) {
            if (stage == STAGE_TERM) {
                *outputTypes |= GEN_TIMER_DATA;
            } else if (stage == STAGE_CYCLE) {
                if (param->savetimersdt > 0.0) {
                    dumpTime = param->savetimerstime + param->savetimersdt;
                    if (timeNow >= dumpTime) {
                        param->savetimerstime = timeNow;
                        param->savetimerscounter++;
                        *outputTypes |= GEN_TIMER_DATA;
                    }
                } else if ((home->cycle % param->savetimersfreq) == 0) {
                    param->savetimerscounter++;
                    *outputTypes |= GEN_TIMER_DATA;
                }
            }
        }

/*
 *      If FLUX file dumps are enabled, do so if the code is
 *      at the end of a cycle and either the elapsed simulation
 *      time since the last dump has exceeded the allowable delta
 *      time between file dumps OR the current cycle is a multiple
 *      of the specified dump frequency.
 *
 *      NOTE: The last two conditions are mutually exclusive. If
 *            a delta time has been specified to determine the 
 *            time to dump files, no check of the current
 *            cycle against the dump frequency will be done.
 */
        if (param->writeFlux && (stage == STAGE_CYCLE)) {
            if (param->writeFluxDT > 0.0) {
                dumpTime = param->writeFluxTime + param->writeFluxDT;
                if (timeNow >= dumpTime) {
                    param->writeFluxTime = timeNow;
                    param->writeFluxCounter++;
                    *outputTypes |= GEN_FLUX_DATA;
                }
            } else if ((home->cycle % param->writeFluxFreq) == 0) {
                param->writeFluxCounter++;
                *outputTypes |= GEN_FLUX_DATA;
            }
        }

/*
 *      If properties file dumps are enabled, do so if the code is
 *      at the end of a cycle and either the elapsed simulation
 *      time since the last dump has exceeded the allowable delta
 *      time between file dumps OR the current cycle is a multiple
 *      of the specified dump frequency.
 *
 *      NOTE: The last two conditions are mutually exclusive. If
 *            a delta time has been specified to determine the 
 *            time to dump files, no check of the current
 *            cycle against the dump frequency will be done.
 */
        if (param->saveprop && (stage == STAGE_CYCLE)) {
            if (param->savepropdt > 0.0) {
                dumpTime = param->saveproptime + param->savepropdt;
                if (timeNow >= dumpTime) {
                    param->saveproptime = timeNow;
                    *outputTypes |= GEN_PROPERTIES_DATA;
                }
            } else if ((home->cycle % param->savepropfreq) == 0) {
                *outputTypes |= GEN_PROPERTIES_DATA;
            }
        }

/*
 *      If POSTSCRIPT file dumps are enabled, do so if the code is
 *      at the end of a cycle and either the elapsed simulation
 *      time since the last dump has exceeded the allowable delta
 *      time between file dumps OR the current cycle is a multiple
 *      of the specified dump frequency.
 *
 *      NOTE: The last two conditions are mutually exclusive. If
 *            a delta time has been specified to determine the 
 *            time to dump files, no check of the current
 *            cycle against the dump frequency will be done.
 */
        if ((param->psfile) && (stage == STAGE_CYCLE)) {
            if (param->psfiledt > 0.0) {
                dumpTime = param->psfiletime + param->psfiledt;
                if (timeNow > dumpTime) {
                    *outputTypes |= GEN_POSTSCRIPT_DATA;
                }
            } else {
                if ((home->cycle % param->psfilefreq) == 0) {
                    *outputTypes |= GEN_POSTSCRIPT_DATA;
                }
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GetParIOGroup
 *      Description: Determine the parallel I/O group of which the current
 *                   task is a member, as well as the previous and
 *                   next tasks in the same group, etc.
 *
 *------------------------------------------------------------------------*/
void GetParallelIOGroup(Home_t *home)
{
        int     taskID, group;
        int     thisGroupSize;
        int     firstInGroup, lastInGroup, thisGroup;

        Param_t *param = home->param;

        int  numIOGroups    = param->numIOGroups;
        int  numTasks       = home->numDomains;
        int  thisTask       = home->myDomain;

        int  smallGroupSize = numTasks / numIOGroups;
        int  numLargeGroups = numTasks % numIOGroups;
        int  largeGroupSize = smallGroupSize + (numLargeGroups > 0);

        int *taskList       = (int *)malloc(numIOGroups * sizeof(int));

/*
 *      First make a list including only the processes that are
 *      the last in their respective IO groups
 */
        for (taskID = 0, group = 0; group < numIOGroups; group++) 
        {
            thisGroupSize = (group < numLargeGroups) ?  largeGroupSize : smallGroupSize;

            lastInGroup = taskID + (thisGroupSize-1);

            taskList[group] = lastInGroup;
            taskID += thisGroupSize;
        }

/*
 *      Now calculate the task specific stuff like group ID, first
 *      and last tasks in the group, and so on.
 */
        if ((numLargeGroups * largeGroupSize) > thisTask) {
            thisGroup = thisTask / largeGroupSize;
            firstInGroup = thisGroup * largeGroupSize;
            lastInGroup = firstInGroup + largeGroupSize - 1;
        } else {
            thisGroup = numLargeGroups +
                        ((thisTask - (numLargeGroups * largeGroupSize)) /
                         smallGroupSize);
            firstInGroup = numLargeGroups * largeGroupSize +
                           (thisGroup - numLargeGroups) * smallGroupSize;
            lastInGroup = firstInGroup + smallGroupSize -1;
        }

        home->ioGroupNum       = thisGroup;
        home->firstInIOGroup   = firstInGroup;
        home->lastInIOGroup    = lastInGroup;
        home->prevInIOGroup    = MAX(thisTask-1, firstInGroup);
        home->nextInIOGroup    = MIN(thisTask+1, lastInGroup);
        home->isLastInIOGroup  = (thisTask == lastInGroup);
        home->isFirstInIOGroup = (thisTask == firstInGroup);

#ifdef PARALLEL
/*
 *      Create a new MPI group containing the processes identified above,
 *      and a new communicator encompasing those processes.
 */
        MPI_Group groupCommWorld, groupLastInIOGroup;

        MPI_Comm_group (MPI_COMM_WORLD, &groupCommWorld);
        MPI_Group_incl (groupCommWorld, numIOGroups, taskList, &groupLastInIOGroup);
        MPI_Comm_create(MPI_COMM_WORLD, groupLastInIOGroup, &home->commLastInIOGroup);

/*
 *      Free the MPI group definitions which are no longer needed
 */
        MPI_Group_free(&groupCommWorld);
        MPI_Group_free(&groupLastInIOGroup);
#endif

        free(taskList);

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GetCounts
 *
 *      Description: Obtain both local and global node/segment counts.
 *
 *      Args:
 *          numNodes      location in which to return to caller the total
 *                        number of nodes in the simulation.
 *          numSegs       location in which to return to caller the total
 *                        number of unique segments in the simulation.
 *          numArms       location in which to return the sum arm count
 *                        of all nodes in the simulation
 *          numLocalNodes location in which to return the number of nodes
 *                        in the local domain.
 *          numLocalSeg   location in which to return the number of unique
 *                        segments in the local domain.
 *
 *      Last Modified: 09/10/2008 - gh  Modified to add the <numLocalNodes>
 *                                  and <numLocalSegs> parameters.
 *
 *------------------------------------------------------------------------*/
static void GetCounts(Home_t *home, int *numNodes, int *numSegs, int *numArms,
                      int *numLocalNodes, int *numLocalSegs)
{
        int   locVals   [3] = { 0,0,0 };
        int   globalVals[3] = { 0,0,0 };

        for (int i=0; i<home->newNodeKeyPtr; i++) 
        {
            Node_t *node = home->nodeKeys[i];
             
            if (node)
            {
                locVals[0] += 1;
                locVals[2] += node->numNbrs;
    
                for (int j=0; j < node->numNbrs; j++) 
                {
                    if ((home->myDomain == node->nbrTag[j].domainID) &&
                        (node->nbrTag[j].index < i)) {
                        continue;
                    }
                    locVals[1] += 1;
                }
            }
        }

#ifdef PARALLEL
        MPI_Reduce((int *) locVals, (int *) globalVals, 3, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        *numNodes = globalVals[0];
        *numSegs  = globalVals[1];
        *numArms  = globalVals[2];
#else
        *numNodes = locVals[0];
        *numSegs  = locVals[1];
        *numArms  = locVals[2];
#endif

        *numLocalNodes = locVals[0];
        *numLocalSegs  = locVals[1];

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    SetBinFileArrays
 *
 *      Description: When writing HDF5 files, the nodal data is aggregated
 *                   into arrays of like elements.  This function allocates
 *                   the needed arrays, pulls node and segment data from the
 *                   Node_t structures, and fills in the new arrays so they
 *                   are ready to write.
 *
 *      Args:
 *          binData  pointer to structure with arrays for binary writes.
 *                   This structure will be populated for the caller.
 *
 *      Last Modified: 09/10/2008 - original version
 *
 *------------------------------------------------------------------------*/
static void SetBinFileArrays(Home_t *home, BinFileData_t *binData)
{
        int    size;

        if (binData->nodeCount == 0) { return; }

/*
 *      Allocate 1 element per node for these arrays
 */
        size = binData->nodeCount * sizeof(int);
        binData->nodeIndex      = (int *)malloc(size);
        binData->nodeConstraint = (int *)malloc(size);
        binData->nodeNumSegs    = (int *)malloc(size);

/*
 *      Need 3 elements per node for position data
 */
        size = binData->nodeCount * sizeof(double) * 3;
        binData->nodePos = (real8 *)malloc(size);

/*
 *      Need 3 ints to define endpoints for each segment.  First
 *      int is index of first node (always the local domain, so
 *      no need to store domain id), last two ints for tag of
 *      2nd endpoint which may or may not be in the local domain.
 */
        size = binData->segCount * sizeof(int) * 3;
        binData->segTags = (int *)malloc(size);

/*
 *      Need 3 elements per segment for these arrays
 */
        size = binData->segCount * sizeof(real8) * 3;
        binData->burgersVec = (real8 *)malloc(size);
        binData->glidePlane = (real8 *)malloc(size);

/*
 *      Go through each node and store it's data in the appropriate
 *      arrays.
 */
        int nOffset = 0;
        int sOffset = 0;

        for (int i=0; i < home->newNodeKeyPtr; i++) 
        {
            Node_t *node = home->nodeKeys[i];

            if (node == (Node_t *)NULL) { continue; }

            binData->nodeIndex[nOffset] = node->myTag.index;
            binData->nodeConstraint[nOffset] = node->constraint;
            binData->nodeNumSegs[nOffset] = node->numNbrs;
            binData->nodePos[nOffset*3  ] = node->x;
            binData->nodePos[nOffset*3+1] = node->y;
            binData->nodePos[nOffset*3+2] = node->z;

            nOffset++;
/*
 *          Add any segment with a remote neighbor to the segment
 *          list.  Completely local segments only get added if the
 *          current node 'owns' the segment.
 */
            for (int j=0; j < node->numNbrs; j++) 
            {
                if ((home->myDomain == node->nbrTag[j].domainID) &&
                    (node->nbrTag[j].index < i)) {
                    continue;
                }

                binData->segTags[sOffset*3  ] = node->myTag.index;
                binData->segTags[sOffset*3+1] = node->nbrTag[j].domainID;
                binData->segTags[sOffset*3+2] = node->nbrTag[j].index;

                binData->burgersVec[sOffset*3  ] = node->burgX[j];
                binData->burgersVec[sOffset*3+1] = node->burgY[j];
                binData->burgersVec[sOffset*3+2] = node->burgZ[j];

                binData->glidePlane[sOffset*3  ] = node->nx[j];
                binData->glidePlane[sOffset*3+1] = node->ny[j];
                binData->glidePlane[sOffset*3+2] = node->nz[j];

                sOffset++;
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    FreeBinFileArrays
 *
 *      Description: This function frees the arrays that had been written
 *                   to the HDF5 file.
 *
 *      Args:
 *          binData  pointer to structure with array pointers.  On
 *                   completion, the arrays will have been freed and the
 *                   array pointers in <binData> set to NULL.
 *
 *      Last Modified: 09/10/2008 - original version
 *
 *------------------------------------------------------------------------*/
static void FreeBinFileArrays(BinFileData_t *binData)
{
        if (binData->nodeIndex      != (int   *)NULL) { free(binData->nodeIndex);      }
        if (binData->nodeConstraint != (int   *)NULL) { free(binData->nodeConstraint); }
        if (binData->nodeNumSegs    != (int   *)NULL) { free(binData->nodeNumSegs);    }
        if (binData->segTags        != (int   *)NULL) { free(binData->segTags);        }
        if (binData->nodePos        != (real8 *)NULL) { free(binData->nodePos);        }
        if (binData->burgersVec     != (real8 *)NULL) { free(binData->burgersVec);     }
        if (binData->glidePlane     != (real8 *)NULL) { free(binData->glidePlane);     }

        binData->nodeIndex      = (int   *)NULL;
        binData->nodeConstraint = (int   *)NULL;
        binData->nodeNumSegs    = (int   *)NULL;
        binData->segTags        = (int   *)NULL;
        binData->nodePos        = (real8 *)NULL;
        binData->burgersVec     = (real8 *)NULL;
        binData->glidePlane     = (real8 *)NULL;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    DoParallelIO
 *      Description: 
 *
 *      Args:
 *          outputTypes  integer bitfield.  All non-zero bits correspond
 *                       to specific types of output that are to be
 *                       generated at this stage of program execution.
 *                       See Util.h for the mapping of output types
 *                       to bit positions.
 *          stage        indicates the program execution stage which
 *                       determines the types of output this function
 *                       might produce.
 *
 *------------------------------------------------------------------------*/
static void DoParallelIO(Home_t *home, int outputTypes, int stage)
{
        int     ioGroup, prevInGroup, nextInGroup, numIOGroups;
        int     thisDomain, isFirstInGroup, isLastInGroup;
        int     writePrologue, writeEpilogue;
        int     sendToken, recvToken;
        int     numSegs = 0, numArms = 0;
        int     nodesWritten = 0, segsWritten = 0;
        int     countInGroup[2] = {0, 0};
        char    baseName[128];
        Param_t *param;
        BinFileData_t binData;
#ifdef PARALLEL
        MPI_Request req;
        MPI_Status  reqStatus;
#endif

        thisDomain  = home->myDomain;
        param       = home->param;

        numIOGroups = param->numIOGroups;

        memset(&binData, 0, sizeof(BinFileData_t));

        binData.firstInGroup = home->firstInIOGroup;
        binData.lastInGroup = home->lastInIOGroup;

        ioGroup = home->ioGroupNum;

        isFirstInGroup = home->isFirstInIOGroup;
        isLastInGroup = home->isLastInIOGroup;

        prevInGroup = home->prevInIOGroup;
        nextInGroup = home->nextInIOGroup;

        recvToken = (thisDomain != prevInGroup);
        sendToken = (thisDomain != nextInGroup);

/*
 *      Domain zero (first member of first I/O group) will always
 *      do output file creation and initialization, and the last
 *      member of the last I/O group always adds any necessary
 *      trailers/whatever to the file
 */
        writePrologue = (thisDomain == 0);
        writeEpilogue = ((ioGroup == (numIOGroups-1)) && isLastInGroup);

/*
 *      Certain output routines require total counts of nodes,
 *      segments, etc.  If we're doing any of these types of output
 *      then do a global operation to get these values before
 *      beginning the I/O.
 */
        if ((outputTypes & (GEN_RESTART_DATA|GEN_SCR_DATA)) ||
            (outputTypes & GEN_TECPLOT_DATA)) {
            GetCounts(home, &param->nodeCount, &numSegs, &numArms,
                      &binData.nodeCount, &binData.segCount);
        }

/*
 *      If we're about to write a binary restart file, set up the
 *      arrays that will be dumped directly to the binary file.
 *      Do this now so that when it is the task's turn to write, the
 *      data is prepped and ready to be written.
 */
        if (((outputTypes & (GEN_RESTART_DATA|GEN_SCR_DATA)) != 0) &&
             (param->writeBinRestart != 0)) {
            SetBinFileArrays(home, &binData);
        }

/*
 *      Call all I/O functions...
 *
 *      Before writing any specific output file, wait for the write
 *      token from the predecesor in the I/O group.  When done
 *      write the specific output type, if this task is not the
 *      last in its I/O group, then send the write token to the
 *      next process in the I/O group so the other task can start
 *      on that output type while the current task continues on
 *      with other output types.
 */
        if ((outputTypes & GEN_RESTART_DATA) != 0) {
            if (stage == STAGE_CYCLE) {
                snprintf(baseName, sizeof(baseName), "rs%04d", param->savecncounter);
            } else {
                snprintf(baseName, sizeof(baseName), "restart.cn");
            }
            if (recvToken) RecvWriteToken(prevInGroup);
            if (param->writeBinRestart) {
                WriteBinaryRestart(home, baseName, ioGroup, isFirstInGroup,
                                   writePrologue, writeEpilogue, &binData);
            } else {
                WriteRestart(home, baseName, ioGroup, isFirstInGroup,
                             writePrologue, writeEpilogue);
            }
            if (sendToken) SendWriteToken(nextInGroup);
        }

#ifdef USE_SCR
        if ((outputTypes & GEN_SCR_DATA) != 0) {
            int statusC = writePrologue, statusD = -1;
            int fileSegCount, maxstep;
            char name[32];
            char file[SCR_MAX_FILENAME];

            fileSegCount = param->numFileSegments; /* save */
            param->numFileSegments = home->numDomains; /* modify */

            SCR_Start_checkpoint();

            if (writePrologue) {
                if (SCR_Route_file("ctrl.0", file) == SCR_SUCCESS) {
                    statusC = WriteCtrl(home, file);
                }
            }

            sprintf(name, "data.%d", thisDomain);
            if (SCR_Route_file(name, file) == SCR_SUCCESS) {
                if (param->writeBinRestart) {
                    binData.firstInGroup = thisDomain;
                    binData.lastInGroup = thisDomain;
                    statusD = WriteBinData(home, file, 1,
                                           writePrologue, writeEpilogue,
                                           &binData);
                } else {
                    statusD = WriteData(home, file, 1,
                                        writePrologue, writeEpilogue);
                }
            }

            SCR_Complete_checkpoint(statusC == 0 && statusD == 0);

            param->numFileSegments = fileSegCount; /* restore */
        }
#endif /* USE_SCR */

/*
 *      Free the checkpoint data.
 */
        if (((outputTypes & (GEN_RESTART_DATA|GEN_SCR_DATA)) != 0) &&
             (param->writeBinRestart != 0)) {
            FreeBinFileArrays(&binData);
        }


        if ((outputTypes & GEN_GNUPLOT_DATA) != 0) {
            if (stage == STAGE_CYCLE) {
                snprintf(baseName, sizeof(baseName), "0t%04d", param->gnuplotcounter);
            } else {
                snprintf(baseName, sizeof(baseName), "gnuplot.final");
            }
            if (recvToken) RecvWriteToken(prevInGroup);
            Gnuplot(home, baseName, ioGroup, isFirstInGroup,
                    writePrologue, writeEpilogue);
            if (sendToken) SendWriteToken(nextInGroup);
        }


        if ((outputTypes & GEN_TECPLOT_DATA) != 0) {
            if (stage == STAGE_CYCLE) {
                snprintf(baseName, sizeof(baseName), "tecdata%04d", param->tecplotcounter);
            } else {
                snprintf(baseName, sizeof(baseName), "tecdata.final");
            }
            if (recvToken) RecvWriteToken(prevInGroup);
            Tecplot(home, baseName, ioGroup, isFirstInGroup,
                    writePrologue, writeEpilogue, numSegs);
            if (sendToken) SendWriteToken(nextInGroup);
        }


        if ((outputTypes & GEN_VELOCITY_DATA) != 0) {
            if (stage == STAGE_CYCLE) {
                snprintf(baseName, sizeof(baseName), "vel%04d", param->velfilecounter);
            } else {
                snprintf(baseName, sizeof(baseName), "vel.final");
            }
            if (recvToken) RecvWriteToken(prevInGroup);
            WriteVelocity(home, baseName, ioGroup, isFirstInGroup,
                          writePrologue, writeEpilogue);
            if (sendToken) SendWriteToken(nextInGroup);
        }


        if ((outputTypes & GEN_FORCE_DATA) != 0) {
            if (stage == STAGE_CYCLE) {
                snprintf(baseName, sizeof(baseName), "force%04d", param->writeForceCounter);
            } else {
                snprintf(baseName, sizeof(baseName), "force.final");
            }
            if (recvToken) RecvWriteToken(prevInGroup);
            WriteForce(home, baseName, ioGroup, isFirstInGroup,
                       writePrologue, writeEpilogue);
            if (sendToken) SendWriteToken(nextInGroup);
        }


        if ((outputTypes & GEN_ARM_DATA) != 0) {
            if (stage == STAGE_CYCLE) {
                snprintf(baseName, sizeof(baseName), "arm%04d", param->armfilecounter);
            } else {
                snprintf(baseName, sizeof(baseName), "arm.final");
            }
            if (recvToken) RecvWriteToken(prevInGroup);
            WriteArms(home, baseName, ioGroup, isFirstInGroup,
                      writePrologue, writeEpilogue);
            if (sendToken) SendWriteToken(nextInGroup);
        }


        if ((outputTypes & GEN_POLEFIG_DATA) != 0) {
            if (stage == STAGE_CYCLE) {
                snprintf(baseName, sizeof(baseName), "polefig%04d", param->polefigcounter);
            } else {
                snprintf(baseName, sizeof(baseName), "polefig.final");
            }
            if (recvToken) RecvWriteToken(prevInGroup);
            WritePoleFig(home, baseName, ioGroup, isFirstInGroup,
                         writePrologue, writeEpilogue);
            if (sendToken) SendWriteToken(nextInGroup);
        }


        if ((outputTypes & GEN_POVRAY_DATA) != 0) {
            if (stage == STAGE_CYCLE) {
                snprintf(baseName, sizeof(baseName), "povframe%04d", param->povraycounter);
            } else {
                snprintf(baseName, sizeof(baseName), "povframe.final");
            }
            if (recvToken) RecvWriteToken(prevInGroup);
            WritePovray(home, baseName, ioGroup, isFirstInGroup,
                        writePrologue, writeEpilogue);
            if (sendToken) SendWriteToken(nextInGroup);
        }


        if ((outputTypes & GEN_VISIT_DATA) != 0) {
            if (stage == STAGE_CYCLE) {
                snprintf(baseName, sizeof(baseName), "visit%04d", param->writeVisitCounter);
            } else {
                snprintf(baseName, sizeof(baseName), "visit.final");
            }
            if (recvToken) RecvWriteToken(prevInGroup);
            WriteVisit(home, baseName, writePrologue, writeEpilogue,
                       &nodesWritten, &segsWritten);
            if (sendToken) SendWriteToken(nextInGroup);
/*
 *          For binary format files, the last process in each IO group
 *          can calculate the total number of segments written by the
 *          group, but for text files, each process must receive the
 *          number of segments written by the predecessor processes in
 *          the group, increment the value by the number it wrote, and
 *          pass it along until the total reaches the last processor in
 *          the IO group.
 */
#ifdef PARALLEL
            if (param->writeVisitSegmentsAsText ||
                param->writeVisitNodesAsText) {

                if (recvToken) {
                    MPI_Recv((int *) countInGroup, 2, MPI_INT, prevInGroup, MSG_VISIT_COUNTS, MPI_COMM_WORLD, &reqStatus);
                }

                countInGroup[0] += nodesWritten;
                countInGroup[1] += segsWritten;

                if (sendToken) {
                    MPI_Isend((int *) countInGroup, 2, MPI_INT, nextInGroup, MSG_VISIT_COUNTS, MPI_COMM_WORLD, &req);
                }
            } else {
                countInGroup[0] += nodesWritten;
                countInGroup[1] += segsWritten;
            }
#else
            countInGroup[0] = nodesWritten;
            countInGroup[1] = segsWritten;
#endif

/*
 *          Last thing is that we have to create a metadata file to go
 *          along with the visit data files.  However the needed per-file
 *          segment/node counts for a given file portion are known
 *          only to the final task in the I/O group responsible
 *          for writing the file.  
 *
 *          Do a reduce among all tasks in new group, gathering
 *          data to the last task in the last IO group (i.e. this is
 *          the same task as the highest numbered task in COMM_WORLD)
 *
 *          Once the data is collected, write the metadata file.
 */
            if (isLastInGroup) {
                int *visitDataCounts;

                visitDataCounts = (int *)malloc(numIOGroups * 2 * sizeof(int));

#ifdef PARALLEL
                MPI_Gather((int *) countInGroup   , 2, MPI_INT, 
                           (int *) visitDataCounts, 2, MPI_INT, numIOGroups-1, home->commLastInIOGroup);
#else
                visitDataCounts[0] = countInGroup[0];
                visitDataCounts[1] = countInGroup[1];
#endif

                if (thisDomain == (home->numDomains-1)) {

                    if (stage == STAGE_CYCLE) {
                        snprintf(baseName, sizeof(baseName), "visit%04d", param->writeVisitCounter);
                    } else {
                        snprintf(baseName, sizeof(baseName), "visit.final");
                    }

                    WriteVisitMetaDataFile(home, baseName, visitDataCounts);
                }

                free(visitDataCounts);
            }
        }

/*
 *      Just to be on the safe side, try flushing all output
 *      to disk before writing the name of the latest restart file
 */
        sync();

/*
 *      If we wrote the text restart file we still need to store the
 *      name of the recently written restart file to disk.  This
 *      involves an explicit syncronization point in the code since
 *      we don't want to do this until all processes have completed
 *      writing their restart data.
 */
        if ((outputTypes & GEN_RESTART_DATA) != 0) {
            if (stage == STAGE_CYCLE) {
                snprintf(baseName, sizeof(baseName), "rs%04d", param->savecncounter);
            } else {
                snprintf(baseName, sizeof(baseName), "restart.cn");
            }
#ifdef PARALLEL
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            if (thisDomain == 0) {
                SetLatestRestart(baseName);
            }
        }

/*
 *      Last thing to do is wait for the completion of the 
 *      asynchronous send (if any) of data counts associated with
 *      the visit output to the next processor in the IO group.
 */
#ifdef PARALLEL
        if (((outputTypes & GEN_VISIT_DATA) != 0) && sendToken &&
             (param->writeVisitSegmentsAsText||param->writeVisitNodesAsText)) {
            MPI_Wait(&req, &reqStatus);
        }
#endif


        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GenerateOutput
 *      Description: This subroutine controls generation of all types
 *                   of output appropriate to the current stage of program
 *                   execution.
 *      Args:
 *          stage    indicates the program execution stage which determines
 *                   the types of output this function might produce.  Valid
 *                   values for this field are:
 *
 *                           STAGE_INIT
 *                           STAGE_CYCLE
 *                           STAGE_TERM
 *
 *------------------------------------------------------------------------*/
void GenerateOutput(Home_t *home, int stage)
{
        int     outputTypes = 0;
        char    fileName[128];
        time_t  tp;
        Param_t *param;


        TimerStart(home, GENERATE_IO);

        param = ( home ? home->param : 0 );

        static int  init=1;

        if (init && (home->myDomain==0))
        {
           home->utc_cycle_start = Current_UTC_Time();
           home->utc_cycle_prev  = home->utc_cycle_start;
           home->utc_cycle_curr  = home->utc_cycle_start;
           init=0;
        }

/*
 *      If we're at the end of a cycle, update the cycle number and
 *      print a status message to stdout.
 */
        if (stage==STAGE_CYCLE) 
        {
            param->timeNow += param->realdt;
            home->cycle++;

            if ((param->loadType != 2) && (home->myDomain == 0)) 
            {
                time(&tp);

                home->utc_cycle_curr = Current_UTC_Time();

                long long dt0 = home->utc_cycle_curr - home->utc_cycle_prev;
                long long dt1 = home->utc_cycle_curr - home->utc_cycle_start;

                int hh=0, mm=0; double ss=0; 
                UTC_To_HHMMSS(dt1,hh,mm,ss);

                char timestr[64];  
                asctime_r(localtime(&tp),timestr);
                timestr[strlen(timestr)-1] = 0;

                int mem_curr = (int)(PDS_Memory_Diagnostics::Get_RSS_Current() >> 20);  // mem_curr = current memory allocation
              //int mem_peak = (int)(PDS_Memory_Diagnostics::Get_RSS_Peak   () >> 20);  // mem_peak = peak    memory allocation

                printf("cycle=%-8d  realdt=%e timeNow=%e date=%s elapsed=%02d:%02d:%04.1lf walldt=%4.1lfs mem=%dMB\n",
                        home->cycle, param->timeNow-param->timeStart, param->timeNow, timestr, hh, mm, ss, (double) dt0/1000000.0, mem_curr );

                home->utc_cycle_prev = home->utc_cycle_curr;
            }
        }

/*
 *      Determine what types of output need to be generated at
 *      this stage.
 */
        GetOutputTypes(home, stage, &outputTypes);

/*
 *      If we are writing either properties data or a restart file,
 *      accumulate the total dislocation density for the entire system
 *      to be written further on.
 */
#ifdef PARALLEL
        if (((outputTypes & GEN_PROPERTIES_DATA) != 0) ||
            ((outputTypes & GEN_RESTART_DATA   ) != 0)) 
        {
            real8 localDensity  = param->disloDensity;
            real8 globalDensity = 0.0;
            MPI_Allreduce((real8 *) &localDensity, (real8 *) &globalDensity, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            param->disloDensity = globalDensity;
        }
#endif

#ifndef NO_XWINDOW
/*
 *      The X-Window plotting is a different beast since all
 *      the data must get filtered through processor zero...
 *      so for now we handle it the old way.
 */
        if ((outputTypes & GEN_XWIN_DATA) != 0) {
            CommSendMirrorNodes(home, stage);
        }

/*
 *      And now the postscript files which are also created from
 *      the X-Window info on processor zero.
 */
        if ((outputTypes & GEN_POSTSCRIPT_DATA) != 0) 
        {
            if (home->myDomain == 0) 
            { printf(" +++ Writing PostScript file\n"); }

            sleep(2);  /* to insure all segments are drawn */
            WinWritePS();
        }
#endif

/*
 *      Dump slip amount and strain decomposition?
 */
        if ((outputTypes & GEN_FLUX_DATA) != 0) 
        {
            sprintf(fileName,  "densflux%04d",home->param->writeFluxCounter);
            WriteDensFlux(fileName,home);
        }

        if (outputTypes & GEN_DENSITY_DATA) 
        {
            WriteDensityField(home, (char *) "densityfield.out");
        }
/*
 *      Preserve various properites i.e. density, etc.
 *               (property vs. time files)
 */
        if ((outputTypes & GEN_PROPERTIES_DATA) != 0)
        {
            WriteProp(home, DENSITY);  /* Actually, length of disloctions   */
                                       /* both total and broken down by     */
                                       /* groups of burgers vector types    */

            WriteProp(home, DENSITY_DELTA);  /* Per-burgers vector density */
                                             /* gain/loss                  */

            WriteProp(home, EPSDOT);   /* Resultant strain rate (scalar)    */

            WriteProp(home, ALL_EPS);  /* Plastic strain tensor             */

            WriteProp(home, EPS);      /* Resultant plastic strain (scalar) */
        }

/*
 *      A number of the types of output will be done (potentially)
 *      in parallel.  Handle those in a separate routine but only
 *      as long as the user-provided flag <skipIO> is not set.
 */
        if (((outputTypes & GEN_RESTART_DATA)  != 0) ||
            ((outputTypes & GEN_SCR_DATA)      != 0) ||
            ((outputTypes & GEN_GNUPLOT_DATA)  != 0) ||
            ((outputTypes & GEN_VELOCITY_DATA) != 0) ||
            ((outputTypes & GEN_FORCE_DATA)    != 0) ||
            ((outputTypes & GEN_ARM_DATA)      != 0) ||
            ((outputTypes & GEN_POLEFIG_DATA)  != 0) ||
            ((outputTypes & GEN_POVRAY_DATA)   != 0) ||
            ((outputTypes & GEN_VISIT_DATA)    != 0) ||
            ((outputTypes & GEN_TECPLOT_DATA)  != 0)) 
        {
            if (!param->skipIO) {
                DoParallelIO(home, outputTypes, stage);
            }
        }

/*
 *      Accumulate/save the deltaTime output.  Note that this is
 *      somewhat different that the other output in that the deltaTime
 *      data is accumulated at each time step and appended in aggregate
 *      whenever the savedtfreq is reached.
 */
        if ( (stage==STAGE_CYCLE) && (param->savedt) )
        {
           DeltaTime_Save(home);
        }

/*      Accumulate/save the node arm histogram diagnostic (if active)
 */
        if ( param->save_narms )
        {
           if (stage==STAGE_CYCLE)  Narms_Diag_Save(home);
           if (stage==STAGE_TERM )  Narms_Diag_Save(home,1);
        }

        TimerStop(home, GENERATE_IO);

/*
 *      Measure dead time waiting for all processes to complete I/O
 */
#ifdef PARALLEL
#ifdef SYNC_TIMERS
        TimerStart (home, IO_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop  (home, IO_BARRIER);
#endif
#endif

/*
 *      If necessary dump latest timing data into a file.  If we're at
 *      process termination, don't restart the timers.
 */
        if ((outputTypes & GEN_TIMER_DATA) != 0) 
        {
            TimerStop    (home, TOTAL_TIME);
            TimeAtRestart(home, stage);

            if (stage == STAGE_CYCLE) 
            {
                TimerStart(home, TOTAL_TIME);
            }
        }

        return;
}
