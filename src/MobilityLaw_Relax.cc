/**************************************************************************
 *
 *      Module:       MobilityLaw_Relax
 *      Description:  Contains functions for calculating mobility of nodes
 *                    using a simple steepest descent.
 *
 *      Includes functions:
 *
 *            MobilityLaw_Relax()
 *                
 ***************************************************************************/
#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "Mobility.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif

/**************************************************************************
 *
 *      Function:     MobilityLaw_Relax
 *      Description:  This function calculates the velocity for a single
 *                    specified node.
 *
 *      Arguments:
 *          IN:  node      Pointer to the node for which to
 *                         calculate velocity
 *      IN/OUT: mobArgs    Structure containing additional
 *                         parameters for conveying information
 *                         to/from the mobility function.
 *
 *      Returns:  0 on success
 *                1 if velocity could not be determined
 *
 *
 *************************************************************************/
int Mobility_Relax(Home_t *home, Node_t *node, MobArgs_t *mobArgs)
{
        int     i, numNbrs;
        real8   Mx, My, Mz;
        real8   dx, dy, dz;
        real8   totLen, invLen;
        Node_t  *nbr;
        Param_t *param;
  
        real8 (*invDragMatrix)[3]    =  mobArgs->invDragMatrix;
     // real8 (*glideConstraints)[3] =  mobArgs->glideConstraints;
        int    *numGlideConstraints  = &mobArgs->numGlideConstraints;

        Matrix33_Zero(invDragMatrix);

/*
 *      This mobility module imposes no glide constraints on the node.
 */
        *numGlideConstraints = 0;

/*
 *      If the node is pinned in ALL dimensions it cannot be moved
 *      so just zero the velocity and return
 */
        if (HAS_ALL_OF_CONSTRAINTS(node->constraint, PINNED_NODE)) {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }
  
/*
 *      If any of the arms of the node has a burgers vector that
 *      has explicitly been set to be sessile (via control file
 *      inputs), the node may not be moved.
 */
        if (NodeHasSessileBurg(home, node)) {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

        param = home->param;
        numNbrs = node->numNbrs;
  
        Mx = param->MobRelaxX;
        My = param->MobRelaxY;
        Mz = param->MobRelaxZ;
  
/*
 *      Since velocity is proportional to total force per unit length
 *      we need to sum 1/2 the lengths of all the node's segments.
 */
        totLen = 0.0;

        for (i = 0; i < numNbrs; i++) {

            if ((nbr = GetNeighborNode(home, node, i)) == (Node_t *)NULL) {
                continue;
            }

            dx = nbr->x - node->x;
            dy = nbr->y - node->y;
            dz = nbr->z - node->z;

            ZImage(param, &dx, &dy, &dz) ;

            totLen += 0.5 * sqrt(dx*dx + dy*dy + dz*dz);
        }

/*
 *      It's possible this function was called for a node which had only zero-
 *      length segments (during SplitSurfaceNodes() for example).  If that is
 *      the case, just set the velocity to zero and return.
 */
        if (totLen == 0.0) {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

/*
 *      By default, velocity is simply proportional to total force per
 *      unit length.  However, user has the option to remove the length
 *      factor from the calculations...
 */
        if (param->MobRelaxScaleByLength) {
            invLen = 1.0 / totLen;
        } else {
            invLen = 1.0;
        }

        node->vX = node->fX * Mx * invLen;
        node->vY = node->fY * My * invLen;
        node->vZ = node->fZ * Mz * invLen;


	invDragMatrix[0][0] = Mx * invLen;
	invDragMatrix[1][1] = My * invLen;
	invDragMatrix[2][2] = Mz * invLen;

        return(0);
}
