
#include "mpi_portability.h"

#include "Home.h"
#include "Comm.h"


/*---------------------------------------------------------------------------
 *
 *      Function:     ResetGlidePlanes
 *      Description:  If needed, this function will recalculate the glide
 *                    plane info for all locally owned segments, then 
 *                    distribute the plane information to the remote
 *                    domains which have these segments as ghosts.
 *
 *-------------------------------------------------------------------------*/
void ResetGlidePlanes(Home_t *home)
{
        Param_t *param;


/*
 *      Currently, the only situation in which we need to go through this
 *      process is when glide planes are being enforced but allowed some
 *      'fuzziness'.
 */
        param = home->param;

        if (!(param->enforceGlidePlanes && param->allowFuzzyGlidePlanes)) {
            return;
        }
 
/*
 *      Loop through all local nodes and recalculate the glide plane
 *      for any segment owned by one of those nodes.
 *
 *      NOTE: Although the loop over nodes is threaded and the
 *            glide-plane update modifies the neighbor node as
 *            well as the local node, we do not actually need
 *            to do any sort of locking.  This is because the
 *            glide plane info for a given segment will be 
 *            accessed and updated by only one thread, so there
 *            is no chance of contention between threads.
 */
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            int nodeID, threadID, threadIterStart, threadIterEnd;

            GetThreadIterationIndices(home->newNodeKeyPtr, &threadID,
                                      &threadIterStart, &threadIterEnd);

/*
 *          Loop through all the nodes...
 */
            for (nodeID = threadIterStart; nodeID < threadIterEnd; nodeID++) {
                int     numNbrs, segIndex;
                Node_t  *node;

                if ((node = home->nodeKeys[nodeID]) == (Node_t *)NULL) {
                    continue;
                }

                numNbrs = node->numNbrs;

                for (segIndex = 0; segIndex < numNbrs; segIndex++) {
                    int     nbrSegID;
                    real8   lDir[3], burg[3], newPlane[3];
                    Node_t  *nbrNode;

                    nbrNode = GetNeighborNode(home, node, segIndex);

                    if (nbrNode == (Node_t *)NULL) {
                        continue;
                    }

/*
 *                  Only check the segment from one end so there's no
 *                  duplication of work done.
 */
                    if (OrderNodes(node, nbrNode) < 0) {
                        continue;
                    }

                    burg[X] = node->burgX[segIndex];
                    burg[Y] = node->burgY[segIndex];
                    burg[Z] = node->burgZ[segIndex];

                    lDir[X] = nbrNode->x - node->x;
                    lDir[Y] = nbrNode->y - node->y;
                    lDir[Z] = nbrNode->z - node->z;

                    ZImage(param, &lDir[X], &lDir[Y], &lDir[Z]);
                    NormalizeVec(lDir);

/*
 *                  Calculate the new glide plane.  If it is undefined,
 *                  the segment is screw and should simply maintain the
 *                  current glide plane.
 */
                    FindPreciseGlidePlane(home, burg, lDir, newPlane,
                                          param->allowFuzzyGlidePlanes);

                    if (fabs(DotProduct(newPlane, newPlane)) < 1.0e-03) {
                        continue;
                    }

/*
 *                  Segment is not screw, so reset the segment' plane
 *                  to the newly calculated one.
 */
                    node->nx[segIndex] = newPlane[X];
                    node->ny[segIndex] = newPlane[Y];
                    node->nz[segIndex] = newPlane[Z];

                    nbrSegID = GetArmID(nbrNode, node);

                    nbrNode->nx[nbrSegID] = newPlane[X];
                    nbrNode->ny[nbrSegID] = newPlane[Y];
                    nbrNode->nz[nbrSegID] = newPlane[Z];
                }
            }

        }  /* end "omp parallel" section */

/*
 *      Now send the glide plane info to all the domains that use it in the
 *      ghost data
 */
        CommSendGhostPlanes(home);
                     
        return;
}
