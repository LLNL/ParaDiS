/**************************************************************************
 *
 *      Module:      ForwardEulerIntegrator.c
 *      Description: Contains functions necessary to implement an
 *                   explicit forwardeuler timestep integration method
 *                   
 *
 ***************************************************************************/

#include "mpi_portability.h"

#include "Home.h"
#include "Comm.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif  /* ifdef _ARLFEM */

#define MIN_DELTA_T 1.0e-11


/*------------------------------------------------------------------------
 *
 *      Function:    PreserveNodalData
 *      Description: Copies the current nodal position/velocity
 *                   data to <old*> values for all local and ghost nodes.
 *                   These old values are needed during the timestep
 *                   integration routines and for calculating plastic
 *                   strain.
 *
 *-----------------------------------------------------------------------*/
static void PreserveNodalData(Home_t *home)
{
        int    i;

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < home->newNodeKeyPtr; i++) {
            Node_t *node;

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            node->oldx = node->x;
            node->oldy = node->y;
            node->oldz = node->z;

            node->oldvX = node->vX;
            node->oldvY = node->vY;
            node->oldvZ = node->vZ;
        }

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < home->ghostNodeCount; i++) {
            Node_t *node;

            node = home->ghostNodeList[i];

            node->oldx = node->x;
            node->oldy = node->y;
            node->oldz = node->z;

            node->oldvX = node->vX;
            node->oldvY = node->vY;
            node->oldvZ = node->vZ;
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       AdvanceAllNodes
 *      Description:    Advance all nodes (local and ghost) in time
 *                      and space by a specified amount.  Note: the
 *                      function assumes current nodal positions and
 *                      velocities have already been preserved in
 *                      the old* variables by a call to PreserveNodalData().
 *
 *-------------------------------------------------------------------------*/
static void AdvanceAllNodes(Home_t *home)
{
        int     i;
        real8   dt;
        Param_t *param;

        param = home->param;
        param->realdt = param->deltaTT;

        dt = param->realdt;

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < home->newNodeKeyPtr; i++){
                real8   x, y, z;
                Node_t  *node;

                if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;

                x = node->oldx + node->oldvX * dt;
                y = node->oldy + node->oldvY * dt;
                z = node->oldz + node->oldvZ * dt;

                FoldBox(param,&x,&y,&z);

                node->x = x;
                node->y = y;
                node->z = z;
        }

/*
 *      Also need to move Ghost nodes
 */
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < home->ghostNodeCount; i++){
                real8   x, y, z;
                Node_t  *node;

                node = home->ghostNodeList[i];

                x = node->oldx + node->oldvX * dt;
                y = node->oldy + node->oldvY * dt;
                z = node->oldz + node->oldvZ * dt;

                FoldBox(param,&x,&y,&z);

                node->x = x;
                node->y = y;
                node->z = z;
        }

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:    ForwardEulerIntegrator()
 *      Description: Use the current nodal velocities, the maximum flight
 *                   distance (param->rmax), and the change in nodal
 *                   velocities since the previous step to determine
 *                   the duration of the next timestep.
 *
 *                   Before returning to the caller, the subroutine
 *                   will advance all nodes to their new positions
 *                   and recalculate the force/velocity of all nodes
 *                   at their new positions.
 *
 *------------------------------------------------------------------------*/
void ForwardEulerIntegrator(Home_t *home)
{
        int      i, doAll = 1;
        real8    vx, vy, vz, v2, vmax, vmax2;
        real8    oldDT, newDT, deltaTT;
        real8    temp1, temp2, rmaxFact;
        real8    dirChange, magV, magOldV, magDeltaV;
        real8    xDelta, yDelta, zDelta, deltaVLimit, deltaAdjFact;
        Param_t  *param;
        Node_t   *node;
        static int firstTime = 1;

/*
 *      If this is the first time into this function, make sure
 *      we initialize the *old* force/velocities before we attempt to 
 *      use them.
 */
        if (firstTime) {
            PreserveNodalData(home);
            firstTime = 0;
        }

        param=home->param;

        vmax  = 0.0;
        vmax2 = 0.0;

        oldDT = param->deltaTT;
        newDT = oldDT * 1.25;
        rmaxFact = param->rmax * 0.20;

        deltaVLimit = 0.15;
        deltaAdjFact = 0.25 / deltaVLimit;

/*
 *      Loop over all the nodes, obtain the maximum velocity from among
 *      all nodes, and calculate what the new timestep would be based
 *      solely on the velocty changes since the previous timestep.
 */
        for (i=0; i<home->newNodeKeyPtr;i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;

            vx = node->vX;
            vy = node->vY;
            vz = node->vZ;

            v2 = vx*vx+vy*vy+vz*vz;

            if (vmax2<v2) {
                vmax2 = v2;
            }

/*
 *          Calculate the node's change in velocity, and (possibly)
 *          adjust the timestep accordingly.
 *
 *          NOTE: Use a minimum (MIN_DELTA_T) beyond which we won't drop
 *          the timestep for velocity changes.  (tests further on may
 *          drop the timestep below this limit to keep nodal movement
 *          below the rmax value)
 */

            magV = sqrt(node->vX*node->vX +
                        node->vY*node->vY +
                        node->vZ*node->vZ);

            magOldV = sqrt(node->oldvX*node->oldvX +
                           node->oldvY*node->oldvY +
                           node->oldvZ*node->oldvZ);

/*
 *          If new velocity is non-zero and the node has not
 *          reversed its direction of movement, check if there
 *          has been a significant enough change in velocity to
 *          warrant altering the timestep duration.
 */
            if ((magV > 0.1) && (magOldV > 0.1)) {

                dirChange = ((node->oldvX*node->vX < 0) ||
                             (node->oldvY*node->vY < 0) ||
                             (node->oldvZ*node->vZ < 0));

                if (!dirChange) {

                    if (magV > magOldV) temp1 = magV/magOldV - 1;
                    else temp1 = magOldV/magV - 1;

                    if (temp1 < deltaVLimit) {
                        temp2 = oldDT * (1.25 - (temp1 * deltaAdjFact));
                        if ((temp2 < newDT) && (magV*newDT > rmaxFact)) {
                            newDT = temp2;
                        }
                    } else {
                        temp2 = oldDT / (1 + temp1);
                        if ((temp2 < newDT) && (magV*newDT > rmaxFact)) {
                            newDT = temp2;
                        }
                    }
                } else {
/*
 *                  FIX ME!  Node has reversed direction... what do we want
 *                  to do here???
 */
                    xDelta = node->vX-node->oldvX;
                    yDelta = node->vY-node->oldvY;
                    zDelta = node->vZ-node->oldvZ;

                    magDeltaV = sqrt(xDelta*xDelta + yDelta*yDelta + zDelta*zDelta);

                    if ((magDeltaV * newDT) > (rmaxFact)) {
                        temp2 = 0.40 * oldDT;
                        if (temp2 < newDT) {
                            newDT = temp2;
                        }
                    }
                }
            }

            if (newDT < MIN_DELTA_T) newDT = MIN_DELTA_T;

        }

        if (newDT < 0.1*oldDT) newDT = 0.1*oldDT;

        vmax = sqrt(vmax2);

/*
 *      No node is permitted to move more than a distance of
 *      param->rmax in a single timestep.  Use the highest nodal
 *      velocity to calculate the maximum time step delta that
 *      can be used while still limiting the flight distance of
 *      the fastest node to rmax.
 */
        if (vmax > 0.0)
            deltaTT = param->rmax / vmax;
        else
            deltaTT = param->maxDT;

        if (newDT < deltaTT) deltaTT = newDT;

#ifdef PARALLEL
        MPI_Allreduce(&deltaTT, &param->deltaTT, 1, MPI_DOUBLE, MPI_MIN,
                      MPI_COMM_WORLD);
#else
        param->deltaTT = deltaTT;
#endif

	param->timeStart = param->timeNow;

/*
 *      Update the force and velocity for ghost nodes, and then reposition
 *      all nodes based on their velocities and the new timestep.
 */
        PreserveNodalData(home);
        AdvanceAllNodes(home);

#ifdef _ARLFEM
/*
 *      Now that the nodes have been repositioned, we need to deal
 *      with any nodes that may have moved outside the free surfaces,
 *      and also rebuild the list of surface-intersecting segments.
 */
        FixSurfaceTopology(home);
#endif

/*
 *      Recalculate forces and velocities at the new positions for all nodes
 */
        NodeForce(home, FULL);
        CalcNodeVelocities(home, 0, doAll);
        CommSendVelocity(home);

        return;
}
