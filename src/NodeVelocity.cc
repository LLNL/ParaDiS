//-------------------------------------------------------------------------
//
//      Module:       NodeVelocity.c
//      Description:  Contains functions to control setting nodal
//                    velocities, generating velocity statistics,
//                    and setting the timestep.
//
//-------------------------------------------------------------------------

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Mobility.h"
#include "Util.h"

// Define the types of velocity statistics to be accumulated.
// NOTE: V_MAXSTATS should always be the final item in the
// enumerated list below.

typedef enum
{
    V_NODE_COUNT = 0,
    V_AVERAGE_X     ,
    V_AVERAGE_Y     ,
    V_AVERAGE_Z     ,
    V_VAR           ,
    V_MAXSTATS
} VStatTypes_t;

#ifdef VEL_STATISTICS
//-------------------------------------------------------------------------
//
//      Function:    GetVelocityStatistics
//      Description: If gathering of velocity statistics is enabled,
//                   gather the statistics (defined by the VStatTypes_t
//                   above)
//
//-------------------------------------------------------------------------

void GetVelocityStatistics(Home_t *home)
{
        real8    velStatsLocal [V_MAXSTATS] = { 0.0 };
        real8    velStatsGlobal[V_MAXSTATS] = { 0.0 };

        // Loop over all the nodes, accumulating all the necessary data

        for (int i=0; i<home->newNodeKeyPtr; i++)
        {
            Node_t *node = home->nodeKeys[i];

            if (node)
            {
                real8 vx = node->vX;
                real8 vy = node->vY;
                real8 vz = node->vZ;

                real8 v2 = vx*vx + vy*vy + vz*vz;

                // If we're gathering statistics, accumulate the necessary data.
                // Otherwise just find the highest nodal velocity.

                velStatsLocal[V_AVERAGE_X] += vx;
                velStatsLocal[V_AVERAGE_Y] += vy;
                velStatsLocal[V_AVERAGE_Z] += vz;
                velStatsLocal[V_VAR]       += v2;
                velStatsLocal[V_NODE_COUNT]++;
            }
        }

        if (velStatsLocal[V_NODE_COUNT] > 0)
        {
            MPI_Allreduce(velStatsLocal, velStatsGlobal, V_MAXSTATS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            int   nodeCount   = velStatsGlobal[V_NODE_COUNT];

            real8 vAveragex   = velStatsGlobal[V_AVERAGE_X] / nodeCount;
            real8 vAveragey   = velStatsGlobal[V_AVERAGE_Y] / nodeCount;
            real8 vAveragez   = velStatsGlobal[V_AVERAGE_Z] / nodeCount;

            real8 vAveragesq  =   vAveragex*vAveragex 
                                + vAveragey*vAveragey 
                                + vAveragez*vAveragez;

            real8 vStDev      = sqrt(velStatsGlobal[V_VAR] / nodeCount) - vAveragesq;

            home->param->vStDev     = vStDev;
            home->param->vAverage   = sqrt(vAveragesq);
        }

}
#endif  // VEL_STATISTICS

//-------------------------------------------------------------------------
//
//      Function:    ApplyNodeConstraintsToVelocity
//      Description: Determines if the node has any constraints that
//                   affect the node's velocity and applies the appropriate
//                   velocity adjustments.
//
//      Arguments:
//          IN:  node  Pointer to node to be adjusted
//
//-------------------------------------------------------------------------

void ApplyNodeConstraintsToVelocity(Home_t *home, Node_t *node)
{
   if (HAS_ANY_OF_CONSTRAINTS(node->constraint, PINNED_X)) { node->vX = 0.0; }
   if (HAS_ANY_OF_CONSTRAINTS(node->constraint, PINNED_Y)) { node->vY = 0.0; }
   if (HAS_ANY_OF_CONSTRAINTS(node->constraint, PINNED_Z)) { node->vZ = 0.0; }
}


//-------------------------------------------------------------------------
//      Function:    CalcNodeVelocities
//      Description: Driver function that will invoke the appropriate
//                   mobility function to update the velocity of every
//                   native node, then apply the velocity cutoff if
//                   applicable.
//
//      Arguments:
//          zeroOnErr  Flag indicating if nodal velocity should be
//                     zeroed for any node for which the mobility
//                     function was unable to calculate a velocity.
//          doAll      Flag indicating if ALL nodes are to have
//                     veolcity recalculated, or just those nodes
//                     that have had their forces updated.
//                     NOTE: Unless we're using a mobility function that
//                     handles segments intersecting inclusions in a
//                     special way, the caller should likely always set
//                     this flag to 1.
//
//      Returns:  0 if velocity was successfully calculated for all nodes
//                1 if the mobility functions were unable to converge
//                  on a velocity for one or more nodes.
//-------------------------------------------------------------------------

int CalcNodeVelocities(Home_t *home, int zeroOnErr, int doAll)
{
        TimerStart(home, CALC_VELOCITY);

        int domainMobError = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+ : domainMobError) shared(doAll, zeroOnErr)
#endif
        {
            int  threadMobError=0;
            int  threadID=0, threadIterStart=0, threadIterEnd=0;

            GetThreadIterationIndices(home->newNodeKeyPtr, &threadID, &threadIterStart, &threadIterEnd);

            for (int j=threadIterStart; j < threadIterEnd; j++)
            {
                // If we encountered a mobility error on a previous node,
                // we'll probably be cutting the timestep, so don't
                // bother restting the velocity of any subsequent nodes.
                //
                // Note: continue the loop rather than breaking out, though,
                // because a 'break' from the loop would prevent the compiler
                // from threading the loop.

                if (threadMobError != 0) 
                { continue; }

                Node_t  *node = home->nodeKeys[j];

                if ( node == (Node_t *) NULL ) 
                { continue; }

                // If we do not need to recalculate velocity on ALL nodes,
                // skip this node unless the forces have just been updated.

                if ((doAll == 0) && ((node->flags & NODE_RESET_FORCES) == 0)) 
                { continue; }

                // We set a pointer to the appropriate mobility function
                // during initialization, so just invoke the function now by
                // calling EvaluvateMobility(), which calls the mobility
                // function pointer.

                MobArgs_t   mobArgs;
                int         nodeMobError = EvaluateMobility(home, node, &mobArgs);

                threadMobError |= nodeMobError;

                // If we had problems calculating the mobility for
                // this node, do any special handling.

                if (nodeMobError) 
                {
#ifdef DEBUG_TIMESTEP
                    printf("Mobility error on node (%d,%d)\n", node->myTag.domainID, node->myTag.index);
                    PrintNode(node);
#endif
                    if (zeroOnErr) 
                    {
                        node->vX = 0.0;
                        node->vY = 0.0;
                        node->vZ = 0.0;
                    }
                }

                // We also used the 'reset forces' flag to determine if we needed
                // to recalculate velocity, but now it can be reset.

                node->flags &= ~NODE_RESET_FORCES;
            }

            domainMobError += threadMobError;

        }  // end omp parallel section


        // Zero out the 'reset forces' flag for all ghost nodes 

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i=0; i < home->ghostNodeCount; i++)
        {
            Node_t *ghost = home->ghostNodeList[i];

            if (ghost)
            {  ghost->flags &= ~NODE_RESET_FORCES; }
        }

#ifdef VEL_STATISTICS
        GetVelocityStatistics(home);
#endif

#ifdef ESHELBY
        // If we're simulating eshelby inclusiuons and were making use of the
        // segment/particle intersection list, we're now done with it and need
        // to clear it out.

        SegPartListClear(home);
#endif

        TimerStop(home, CALC_VELOCITY);

#ifdef PARALLEL
#ifdef SYNC_TIMERS
        TimerStart (home, CALC_VELOCITY_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop  (home, CALC_VELOCITY_BARRIER);
#endif
#endif

        return(domainMobError);
}
