/**************************************************************************
 *
 *  Function:    Mobility_FCC_climb
 *  Author:      Tom Arsenlis, Sylvie Aubry
 *  Description: This is alternate version of a generic mobility function
 *               for FCC materials.  It is very similar in function and
 *               structure to the BCC_glide mobility although glide planes
 *               rather than burgers vectors are used to identify junctions.
 *               Additionally, it does permit dislocations a small amount
 *               of climb which means that the motion is not completely glide
 *               restricted.  As a result, the segment glide planes must
 *               be recalculated during each timestep after the nodes are
 *               repositioned.
 *
 *               Note: Use of this mobility function explicitly enables
 *               the 'enforceGlidePlanes' flag in the code, but also permits
 *               some 'fuzziness' in those glide planes by setting the
 *               'allowFuzzyGlidePlanes' flag.
 *
 *  Arguments:
 *      IN:  node         Pointer to the node for which to
 *                        calculate velocity
 *      IN/OUT: mobArgs   Structure containing additional
 *                        parameters for conveying information
 *                        to/from the mobility function.
 *
 *  Returns:  0 on success
 *            1 if velocity could not be determined
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

//#define M33_INVERSE(a,b) Matrix33_Inverse(a,b)
//#define M33_INVERSE(a,b) Matrix33_SVD_Inverse(a,b)
  #define M33_INVERSE(a,b) Matrix33_PseudoInverse(a,b)

#define MAX_SEG_PER_NODE 10


int Mobility_FCC_climb(Home_t *home, Node_t *node, MobArgs_t *mobArgs)
{
        int     i, j, numNbrs; 
        int     numNonZeroLenSegs = 0;
        real8   bDotb;
        real8   halfLen;
        real8   angle1, angle2;
        real8   costheta, costheta2, cosCritical1, cosCritical2;
        real8   Bline, Bscrew, Bedge, Bglide, Bclimb, Beclimb;
        real8   Bscrew2, Bedge2;
        real8   invBscrew2, invBedge2;
        real8   BlmBecl;
        real8   eps = 1.0e-12;
        real8   segLen[MAX_SEG_PER_NODE];
        real8   burgList[MAX_SEG_PER_NODE][3];
        real8   segVector[MAX_SEG_PER_NODE][3];
        real8   nDir[MAX_SEG_PER_NODE][3];
        real8   lineDir[MAX_SEG_PER_NODE][3];
        real8   nForce[3], nVel[3];
        real8   Btotal[3][3] = {{0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0}};

        Param_t *param  = home->param;

        real8 (*invDragMatrix)[3]    =  mobArgs->invDragMatrix;
     // real8 (*glideConstraints)[3] =  mobArgs->glideConstraints;
        int    *numGlideConstraints  = &mobArgs->numGlideConstraints;

        Matrix33_Zero(invDragMatrix);
        *numGlideConstraints = 0;

        Bscrew = 1.0 / param->MobScrew;
        Bedge  = 1.0 / param->MobEdge;

/*
 *      Climb is set this way to help with convergence and the fact that
 *      this is a glide restricted mobility
 */
        Beclimb    = 10000.0 * Bedge;
       
        Bedge2     = Bedge * Bedge;
        Bscrew2    = Bscrew * Bscrew;

        Bline      = 1.0 * MIN(Bscrew, Bedge);
        BlmBecl    = Bline - Beclimb; 

        invBscrew2 = 1.0 / (Bscrew2);
        invBedge2  = 1.0 / (Bedge2);

        numNbrs = node->numNbrs;

        angle1 = 1.0;
        angle2 = 5.0;

/*
 *      If node is 'pinned' in ALL dimensions, or the node has any arms 
 *      with a burgers vector that has explicitly been set to be sessile (via
 *      control file inputs), the node may not be moved so just zero the
 *      velocity and return
 */
        if (HAS_ALL_OF_CONSTRAINTS(node->constraint, PINNED_NODE) ||
            NodeHasSessileBurg(home, node)) {

            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;

            return(0);
        }

/*
 *      have to go through all node's segments and set up the
 *      appropriate array of glide planes.  While we're at it, might
 *      as well save some additional info we'll need later.
 */
        cosCritical1 = (M_PI * angle1 / 180.0);
        cosCritical2 = (M_PI * angle2 / 180.0);

        for (i = 0; i < numNbrs; i++) {
            real8   segLen2;
            Node_t  *nbrNode;

            if ((nbrNode = GetNeighborNode(home, node, i)) == (Node_t *)NULL) {
                continue;
            }

            segVector[i][X] = node->x - nbrNode->x;
            segVector[i][Y] = node->y - nbrNode->y;
            segVector[i][Z] = node->z - nbrNode->z;

            ZImage(param, &segVector[i][X], &segVector[i][Y],
                   &segVector[i][Z]);

/*
 *          Skip the zero length segments
 */
            segLen2 = DotProduct(segVector[i], segVector[i]);

            if ((segLen2) < eps) {
                segLen[i] = 0.0;
                continue;
            }

            numNonZeroLenSegs++;

            segLen[i] = sqrt(segLen2);

/*
 *          Save some segment specific stuff for later rather
 */
            burgList[i][X] = node->burgX[i];
            burgList[i][Y] = node->burgY[i];
            burgList[i][Z] = node->burgZ[i];

            VECTOR_COPY(lineDir[i], segVector[i]);
            NormalizeVec(lineDir[i]);

            nDir[i][X] = node->nx[i];
            nDir[i][Y] = node->ny[i];
            nDir[i][Z] = node->nz[i];

/*
 *          If this is one arm of a 2-node, see if it is a screw segment
 *
 *          The segment glide planes are calculated based on whether
 *          the segment is screw or not.  Two angles (angle1 and angle2)
 *          are defined which determine the 'screwness' of the segment.
 *          Any segment within <angle1> degrees of being screw is
 *          considered screw and its glide plane will simply be
 *          maintained.  A segment within <angle2> degrees of screw (but
 *          more than <angle1>) is considered to be in a 'mixed' region
 *          where a glide plane is calculated from a linear combination
 *          of the original screw glide plane and l cross b.  Glide planes
 *          for all other segments are simply calculated via l cross b.
 */
            if (numNbrs == 2) {

                bDotb = DotProduct(burgList[i], burgList[i]);
                costheta = DotProduct(lineDir[i], burgList[i]);
                costheta2 = (costheta*costheta) / bDotb;

                if (costheta2 <= (cosCritical1 * cosCritical1)) {
                    /* not screw */
                    NormalizedCrossVector(burgList[i], lineDir[i], nDir[i]);
                } else if (costheta2 <= (cosCritical2 * cosCritical2)) {
                    /* mixed region */
                    real8  acostheta2, temp;
                    real8  nScrew[3], nNoScrew[3];
                    VECTOR_COPY(nScrew, nDir[i]);
                    NormalizedCrossVector(burgList[i], lineDir[i], nNoScrew);
                    acostheta2 = acos(sqrt(costheta2));
                    temp = cos(0.5*(acostheta2-angle1)*M_PI/(angle2-angle1));
                    temp = 1.0 - (temp * temp);
                    for (j = 0; j < 3; j++) {
                        nDir[i][j] = nScrew[j] + (nNoScrew[j]-nScrew[j])*temp;
                    }
                }
            }
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
 *      Now initialize the velocity
 */
        node->vX = 0.0e0;
        node->vY = 0.0e0;
        node->vZ = 0.0e0;

/*
 *      Begin construction of the node drag matrix
 *
 *      Loop over all arms of the node, adding each arm's contribution
 *      to the drag matrix.
 */
        for (i = 0; i < numNbrs; i++) {
            real8 b[3], d[3], m[3], n[3], nCryst[3];

/*
 *          If the segment is zero length (which can happen when
 *          the mobility function is being called from SplitMultiNodes())
 *          just skip the segment.
 */
            if (segLen[i] < eps) {
                continue;
            }

            VECTOR_COPY(b, burgList[i]);
            VECTOR_COPY(d, lineDir[i]);
            VECTOR_COPY(n, nDir[i]);

            halfLen = segLen[i] * 0.5;

            bDotb = DotProduct(b, b);
            costheta = DotProduct(d, b);
            costheta2 = costheta * costheta / bDotb;

/*
 *          If needed, rotate a copy of the glide plane vector from the
 *          laboratory frame to the crystal frame.
 */
            if (param->useLabFrame) {
                Matrix33Vector3Multiply(param->rotMatrixInverse, n, nCryst);
            } else {
                VECTOR_COPY(nCryst, n);
            }

/*
 *          arms not on [1 1 1] planes don't move as readily as
 *          other arms, so must be handled specially.
 */
            if (fabs(nCryst[X] * nCryst[Y] * nCryst[Z]) < 0.13608) {
                if (numNbrs == 2) {
                    Btotal[0][0] += halfLen * Beclimb;
                    Btotal[1][1] += halfLen * Beclimb;
                    Btotal[2][2] += halfLen * Beclimb;
                } else {
                    Btotal[0][0] += halfLen * (d[X]*d[X]*BlmBecl + Beclimb);
                    Btotal[0][1] += halfLen * (d[X]*d[Y]*BlmBecl);
                    Btotal[0][2] += halfLen * (d[X]*d[Z]*BlmBecl);
                    Btotal[1][1] += halfLen * (d[Y]*d[Y]*BlmBecl + Beclimb);
                    Btotal[1][2] += halfLen * (d[Y]*d[Z]*BlmBecl);
                    Btotal[2][2] += halfLen * (d[Z]*d[Z]*BlmBecl + Beclimb);
                }
            } else {
                NormalizedCrossVector(n, d, m);
                Bglide = sqrt(invBedge2 + (invBscrew2-invBedge2) *
                              costheta2);
                Bglide = 1.0 / Bglide;
                Bclimb = Beclimb;
                Btotal[0][0] += (m[X]*m[X] * Bglide +
                                 n[X]*n[X] * Bclimb +
                                 d[X]*d[X] * Bline) * halfLen;
                Btotal[0][1] += (m[X]*m[Y] * Bglide +
                                 n[X]*n[Y] * Bclimb +
                                 d[X]*d[Y] * Bline) * halfLen;
                Btotal[0][2] += (m[X]*m[Z] * Bglide +
                                 n[X]*n[Z] * Bclimb +
                                 d[X]*d[Z] * Bline) * halfLen;
                Btotal[1][1] += (m[Y]*m[Y] * Bglide +
                                 n[Y]*n[Y] * Bclimb +
                                 d[Y]*d[Y] * Bline) * halfLen;
                Btotal[1][2] += (m[Y]*m[Z] * Bglide +
                                 n[Y]*n[Z] * Bclimb +
                                 d[Y]*d[Z] * Bline) * halfLen;
                Btotal[2][2] += (m[Z]*m[Z] * Bglide +
                                 n[Z]*n[Z] * Bclimb +
                                 d[Z]*d[Z] * Bline) * halfLen;
            }

        }  /* End loop over arms */

        Btotal[1][0] = Btotal[0][1];
        Btotal[2][0] = Btotal[0][2];
        Btotal[2][1] = Btotal[1][2];

/*
 *      At this point we should check if the matrix is invertable and
 *      if it isn't, find the eigen values and eigen vectors of the drag
 *      matrix, then invert the drag matrix keeping zero eigen values
 *      as zero.
 *
 *      FIX ME!  For now, we're assuming the matrix is invertible.
 */
        nForce[0] = node->fX;
        nForce[1] = node->fY;
        nForce[2] = node->fZ;

        if ( M33_INVERSE(invDragMatrix, Btotal) < 0 ) 
        { Fatal("%s::%s(%d) : Cannot invert drag matrix!", __FILE__, __func__, __LINE__ ); }

        Matrix33Vector3Multiply(invDragMatrix, nForce, nVel);

        node->vX = nVel[0];
        node->vY = nVel[1];
        node->vZ = nVel[2];

        return(0);
}
