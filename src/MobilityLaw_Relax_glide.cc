/**************************************************************************
 *
 *      Module:       MobilityLaw_Relax_glide
 *      Description:  Contains functions for calculating mobility of nodes
 *                    using a simple steepest descent, but restricting
 *                    motion to the glide planes defined by the node's
 *                    segments.
 *                    
 *
 *      Includes functions:
 *
 *            MobilityLaw_Relax_glide()
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
 *      Function:     MobilityLaw_Relax_glide
 *      Description:  This function calculates the velocity for a single
 *                    specified node.
 *
 *      Arguments:
 *          IN:  node          Pointer to the node for which to
 *                             calculate velocity
 *          IN/OUT: mobArgs    Structure containing additional
 *                             parameters for conveying information
 *                             to/from the mobility function.
 *      Returns:  0 on success
 *                1 if velocity could not be determined
 *
 *
 *************************************************************************/
int Mobility_Relax_glide(Home_t *home, Node_t *node, MobArgs_t *mobArgs)
{
        int     i;
        int     numNbrs,numNonZeroLenSegs;
        real8   Mx, My, Mz, eps = 1.0e-12;
        real8   totLen, invLen;
        real8   nVel[3], segLen[3];
        Node_t  *nbrNode;
        Param_t *param;

        real8 (*invDragMatrix)[3]    =  mobArgs->invDragMatrix;
     // real8 (*glideConstraints)[3] =  mobArgs->glideConstraints;
        int    *numGlideConstraints  = &mobArgs->numGlideConstraints;

        Matrix33_Zero(invDragMatrix);
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
 *      Copy and normalize glide plane constraints for any non-zero length
 *      segments.  Compute dislocation lengths
 */
        totLen = 0.0;

        numNonZeroLenSegs = 0;
        for (i = 0; i < numNbrs; i++) {
            real8 oneSegLen;

            nbrNode = GetNeighborNode(home, node, i);

            if (nbrNode == (Node_t *)NULL) {
                continue;
            }

/*
 *          Skip any zero length segments.
 */
            segLen[X] = nbrNode->x - node->x;
            segLen[Y] = nbrNode->y - node->y;
            segLen[Z] = nbrNode->z - node->z;

            ZImage(param, &segLen[X], &segLen[Y], &segLen[Z]);

            oneSegLen = sqrt(DotProduct(segLen, segLen));

            if (oneSegLen < eps) {
                continue;
            }

            totLen += 0.5 * oneSegLen;

            numNonZeroLenSegs++;
        }

/*
 *      It's possible this function was called for a node which had only zero-
 *      length segments (during SplitSurfaceNodes() for example).  If that is
 *      the case, just set the velocity to zero and return.
 */
        if (numNonZeroLenSegs == 0) {
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

        nVel[X] = Mx * node->fX * invLen;
        nVel[Y] = My * node->fY * invLen;
        nVel[Z] = Mz * node->fZ * invLen;

        invDragMatrix[0][0] = Mx * invLen;
        invDragMatrix[1][1] = My * invLen;
        invDragMatrix[2][2] = Mz * invLen;

        node->vX = nVel[X];
        node->vY = nVel[Y];
        node->vZ = nVel[Z];

        return(0);
}
