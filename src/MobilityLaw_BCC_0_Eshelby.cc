/**************************************************************************
 *
 *      Module:       MobilityLaw_BCC_0_Eshelby
 *      Description:  Contains functions for calculating mobility of nodes
 *                    in BCC metals with adjustments for dislocation
 *                    segments that intersect Eshelby inclusions.
 *                    Based on the Arsenlis matlab code mobbcc0_Eshelbyy.m
 *
 *      Includes functions:
 *
 *            Mobility_BCC_0_Eshelby()
 *                
 ***************************************************************************/
#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "Mobility.h"
#include "M33.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif

#define NUMLIMITS_MAX (2*MAX_SEGPART_INTERSECTIONS+1)

//#define M33_INVERSE(a,b) Matrix33_Inverse(a,b)
//#define M33_INVERSE(a,b) Matrix33_SVD_Inverse(a,b)
  #define M33_INVERSE(a,b) Matrix33_PseudoInverse(a,b)


/**************************************************************************
 *
 *      Function:     Mobility_BCC_0_Eshelby
 *      Description:  This function calculates the velocity for a single
 *                    specified node.
 *
 *      Arguments:
 *          IN:  node                Pointer to the node for which to
 *                                   calculate velocity
 *          IN/OUT: mobArgs          Structure containing additional
 *                                   parameters for conveying information
 *                                   to/from the mobility function.
 *
 *      Returns:  0 on success
 *                1 if velocity could not be determined
 *
 *************************************************************************/
int  Mobility_BCC_0_Eshelby(Home_t *home, Node_t *node, MobArgs_t *mobArgs)
{
        int     i, j, m, n, nbrs;
        int     numNonZeroLenSegs = 0;
        real8   tmp, tmp3[3], vOld[3], totLength;
        real8   massMult=0, massMatrix[3][3];
        real8   bx, by, bz;
        real8   dx, dy, dz;
        real8   mx, my, mz;
        real8   nx, ny, nz;
        real8   mag, mag2, halfMag, invMag;
        real8   invbMag2, bMag2, costheta, costheta2, invsqrt1mcostheta2;
        real8   Bline, Bscrew, Bedge, Bglide, Bclimb, Beclimb;
        real8   Bscrew2, Beclimb2;
        real8   invBscrew2, invBedge2;
        real8   BlmBsc, BclmBsc, BglmBsc, BlmBecl;
        real8   eps = 1.0e-12;
        real8   nForce[3], nVel[3];
        real8   burgCryst[3];
        real8   Bseg[3][3];
        real8   Btotal[3][3];
#ifdef ESHELBY
        real8   BsegAlt[3][3];
#endif
        Node_t  *nbrNode;
        Param_t *param;

        real8 (*invDragMatrix)[3]   =  mobArgs->invDragMatrix;
        int    *numGlideConstraints = &mobArgs->numGlideConstraints;

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

        param = home->param;

        Bscrew     = 1.0 / param->MobScrew;
        Bedge      = 1.0 / param->MobEdge;
        Beclimb    = 1.0 / param->MobClimb;

        Bscrew2    = Bscrew * Bscrew;
        Beclimb2   = Beclimb * Beclimb;

        Bline      = 1.0e-2 * MIN(Bscrew, Bedge);
        BlmBsc     = Bline - Bscrew;
        BlmBecl    = Bline - Beclimb;

        invBscrew2 = 1.0 / (Bscrew*Bscrew);
        invBedge2  = 1.0 / (Bedge*Bedge);

#ifdef ESHELBY
        real8 alpha           = param->MobEshelbyResist;
        real8 Bscrew_Eshelby  = alpha * Bscrew;
        real8 Bedge_Eshelby   = alpha * Bedge;
        real8 Beclimb_Eshelby = alpha * Beclimb;

        real8 BlmBscE         = Bline - Bscrew_Eshelby;
        real8 BlmBeclE        = Bline - Beclimb_Eshelby;

        real8 BscrewE2        = Bscrew_Eshelby  * Bscrew_Eshelby ;
        real8 BeclimbE2       = Beclimb_Eshelby * Beclimb_Eshelby;

        real8 invBscrewE2     = 1.0 / (Bscrew_Eshelby*Bscrew_Eshelby);
        real8 invBedgeE2      = 1.0 / (Bedge_Eshelby *Bedge_Eshelby );
#endif

        nbrs = node->numNbrs;

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                Btotal[i][j]  = 0.0;
            }
        }
/*
 *      If we need to include inertial terms...
 */
        if (param->includeInertia) {
            massMult   = 0.25 * param->massDensity * 
                         (param->burgMag * param->burgMag);

            vOld[0] = node->oldvX;
            vOld[1] = node->oldvY;
            vOld[2] = node->oldvZ;
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

/*
 *      Loop over all arms of the node, adding each arm's contribution
 *      to the drag matrix.
 */
        totLength = 0.0;

#ifdef ESHELBY
        real8  pos1[3];
        pos1[X] = node->x;
        pos1[Y] = node->y;
        pos1[Z] = node->z;
#endif

        for (i = 0; i < nbrs; i++) {  
            for (m = 0; m < 3; m++) {
                for (n = 0; n < 3; n++) {
                    Bseg[m][n]    = 0.0;
#ifdef ESHELBY
                    BsegAlt[m][n] = 0.0;
#endif
                }
            }

            bx = node->burgX[i];
            by = node->burgY[i];
            bz = node->burgZ[i];

            bMag2 = (bx*bx + by*by + bz*bz);
            invbMag2 = 1.0 / bMag2;

/*
 *          Calculate the length of the arm and its tangent line direction
 */
            nbrNode = GetNeighborNode(home, node, i);

            if (nbrNode == (Node_t *)NULL) {
                printf("WARNING: Neighbor not found at %s line %d\n",
                       __FILE__, __LINE__);
                continue;
            }

            dx = nbrNode->x - node->x;
            dy = nbrNode->y - node->y;
            dz = nbrNode->z - node->z;

            ZImage(param, &dx, &dy, &dz);


#ifdef ESHELBY
            real8  pos2[3];
            pos2[X] = node->x + dx;
            pos2[Y] = node->y + dy;
            pos2[Z] = node->z + dz;
#endif

            mag2    = dx*dx + dy*dy + dz*dz;
/*
 *          If the segment is zero length (which can happen when
 *          the mobility function is being called from SplitMultiNodes())
 *          just skip the segment.
 */
            if (mag2 < eps) {
                continue;
            }

            numNonZeroLenSegs++;

            mag     = sqrt(mag2);
            halfMag = 0.5 * mag;
            invMag  = 1.0 / mag;

            dx *= invMag;
            dy *= invMag;
            dz *= invMag;

/*
 *          Calculate how close to screw the arm is
 */
            costheta = (dx*bx + dy*by + dz*bz);
            costheta2 = (costheta*costheta) * invbMag2;

/*
 *          [0 0 1] arms don't move as readily as other arms, so must be
 *          handled specially.
 *
 *          If needed, rotate a copy of the burgers vector from the
 *          laboratory frame to the crystal frame.
 */
            if (param->useLabFrame) {
                real8 bTmp[3] = {bx, by, bz};
                Matrix33Vector3Multiply(param->rotMatrixInverse,bTmp,burgCryst);
            } else {
                burgCryst[X] = bx;
                burgCryst[Y] = by;
                burgCryst[Z] = bz;
            }

            if (fabs(burgCryst[X]*burgCryst[Y]*burgCryst[Z]) < eps) {
                Bseg[0][0] += (dx*dx * BlmBecl + Beclimb);
                Bseg[0][1] += (dx*dy * BlmBecl);
                Bseg[0][2] += (dx*dz * BlmBecl);
                Bseg[1][1] += (dy*dy * BlmBecl + Beclimb);
                Bseg[1][2] += (dy*dz * BlmBecl);
                Bseg[2][2] += (dz*dz * BlmBecl + Beclimb);
            } else  {
/*
 *              Arm is not [0 0 1], so build the drag matrix assuming the
 *              dislocation is screw type
 */
                Bseg[0][0] += (dx*dx * BlmBsc + Bscrew);
                Bseg[0][1] += (dx*dy * BlmBsc);
                Bseg[0][2] += (dx*dz * BlmBsc);
                Bseg[1][1] += (dy*dy * BlmBsc + Bscrew);
                Bseg[1][2] += (dy*dz * BlmBsc);
                Bseg[2][2] += (dz*dz * BlmBsc + Bscrew);

/*
 *              Now correct the drag matrix for dislocations that are
 *              not screw
 */
                if ((1.0 - costheta2) > eps) {

                    invsqrt1mcostheta2 = 1.0 / sqrt((1.0 - costheta2) * bMag2);
#ifdef NAN_CHECK
                    if (isnan(invsqrt1mcostheta2) != 0) {
                        Fatal("Mobility_BCC_0: invsqrt1mcostheta2 = "
                              "NaN\n  invsqrt1mcostheta2 = 1.0 / "
                              "sqrt((1.0 - costheta2) * bMag2)\n  where"
                              "costheta2 = %lf and bMag2 = %lf", costheta2,
                              bMag2);
                    }
#endif
                    xvector(bx, by, bz, dx, dy, dz, &nx, &ny, &nz);
                    nx *= invsqrt1mcostheta2;
                    ny *= invsqrt1mcostheta2;
                    nz *= invsqrt1mcostheta2;

                    xvector(nx, ny, nz, dx, dy, dz, &mx, &my, &mz);


                    Bglide = sqrt(invBedge2+(invBscrew2-invBedge2)*costheta2);
                    Bglide = 1.0 / Bglide;
                    Bclimb = sqrt(Beclimb2 + (Bscrew2 - Beclimb2) * costheta2);

#ifdef NAN_CHECK
                    if (isnan(Bglide) != 0) {
                        Fatal("Mobility_BCC_0: Bglide = NaN\n"
                              "  Bglide = sqrt(invBedge2 + "
                              "(invBscrew2-invBedge2)*costheta2)\n"
                              "  where invBedge2 = %lf, invBscrew2 = %lf, "
                              "costheta2 = %lf", invBedge2, invBscrew2,
                              costheta2);
                    }

                    if (isnan(Bclimb) != 0) {
                        Fatal("Mobility_BCC_0: Bclimb = NaN\n"
                              "  Bclimb = sqrt(Beclimb2 + "
                              "(Bscrew2-Beclimb2)*costheta2)\n"
                              "  where Beclimb2 = %lf, Bscrew2 = %lf, "
                              "costheta2 = %lf", Beclimb2, Bscrew2,
                              costheta2);
                    }
#endif
                    BclmBsc = Bclimb - Bscrew;
                    BglmBsc = Bglide - Bscrew;

                    Bseg[0][0] += (nx*nx * BclmBsc + mx*mx * BglmBsc);
                    Bseg[0][1] += (nx*ny * BclmBsc + mx*my * BglmBsc);
                    Bseg[0][2] += (nx*nz * BclmBsc + mx*mz * BglmBsc);
                    Bseg[1][1] += (ny*ny * BclmBsc + my*my * BglmBsc);
                    Bseg[1][2] += (ny*nz * BclmBsc + my*mz * BglmBsc);
                    Bseg[2][2] += (nz*nz * BclmBsc + mz*mz * BglmBsc);
                }
            }  /* End non-[0 0 1] arm */

#ifdef ESHELBY
/*
 *          if the segment intersects any eshelby inclusions, add
 *          in the mobility corrections from the inclusion.
 */
            SegPartIntersect_t *intersection = SegPartListLookup(home, &node->myTag, &nbrNode->myTag);

            if (intersection != (SegPartIntersect_t *)NULL) {
/*
 *              Have to modify the drag matrix since the segment
 *              intersects some particles.
 */
                if (fabs(burgCryst[X]*burgCryst[Y]*burgCryst[Z]) < eps) {
                    BsegAlt[0][0] = (dx*dx*BlmBeclE + Beclimb_Eshelby);
                    BsegAlt[0][1] = (dx*dy*BlmBeclE);
                    BsegAlt[0][2] = (dx*dz*BlmBeclE);
                    BsegAlt[1][1] = (dy*dy*BlmBeclE + Beclimb_Eshelby);
                    BsegAlt[1][2] = (dy*dz*BlmBeclE);
                    BsegAlt[2][2] = (dz*dz*BlmBeclE + Beclimb_Eshelby);
                } else {
                    if ((1.0 - costheta2) > eps) {
                        Bglide = sqrt(invBedgeE2 + (invBscrewE2 - invBedgeE2) * costheta2);
                        Bglide = 1.0 / Bglide;
                        Bclimb = sqrt(BeclimbE2 + (BscrewE2 - BeclimbE2) * costheta2);
                        BsegAlt[0][0] = (nx*nx * Bclimb + mx*mx * Bglide + dx*dx * Bline);
                        BsegAlt[0][1] = (nx*ny * Bclimb + mx*my * Bglide + dx*dy * Bline);
                        BsegAlt[0][2] = (nx*nz * Bclimb + mx*mz * Bglide + dx*dz * Bline);
                        BsegAlt[1][1] = (ny*ny * Bclimb + my*my * Bglide + dy*dy * Bline);
                        BsegAlt[1][2] = (ny*nz * Bclimb + my*mz * Bglide + dy*dz * Bline);
                        BsegAlt[2][2] = (nz*nz * Bclimb + mz*mz * Bglide + dz*dz * Bline);
                    } else {
                        BsegAlt[0][0] = (dx*dx*BlmBscE+Bscrew_Eshelby);
                        BsegAlt[0][1] = (dx*dy*BlmBscE);
                        BsegAlt[0][2] = (dx*dz*BlmBscE);
                        BsegAlt[1][1] = (dy*dy*BlmBscE+Bscrew_Eshelby);
                        BsegAlt[1][2] = (dy*dz*BlmBscE);
                        BsegAlt[2][2] = (dz*dz*BlmBscE+Bscrew_Eshelby);
                    }
                }

/*
 *              Loop over every inclusion intersected by this segment
 */
                int numLimits = 0;
                real8 ratioList[2][NUMLIMITS_MAX];

                for (int iCnt=0; iCnt < intersection->numIntersections; iCnt++) 
                {
                    int ninc;
                    real8 dx, dy, dz;
                    real8 ratio[2];

                    EInclusion_t *inclusion = &home->eshelbyInclusions[intersection->inclusionIndex[iCnt]];

                    dx = inclusion->position[X] - pos1[X];
                    dy = inclusion->position[Y] - pos1[Y];
                    dz = inclusion->position[Z] - pos1[Z];
                    
                    ZImage(param, &dx, &dy, &dz);
                    
                    real8 partPos[3];
                    partPos[X] = pos1[X] + dx;
                    partPos[Y] = pos1[Y] + dy;
                    partPos[Z] = pos1[Z] + dz;

                    int IsIntercept;
                    IsIntercept = IntersectionEllipseSegment(partPos, inclusion->radius,
                                                             inclusion->rotation, pos1, pos2, 
                                                             1.0, &ninc, ratio);
 
                    ratio[0] = MIN(1.0, MAX(ratio[0], 0.0));
                    ratio[1] = MIN(1.0, MAX(ratio[1], 0.0));

/*
 *                  This basically breaks the segment into regions that
 *                  are overlapped by inclusions and regions that do not.
 *                  Regions overlapped by multiple inclusions are
 *                  aggregated together.
 */
                    if (numLimits == 0) {
                        ratioList[0][0] = ratio[0];
                        ratioList[1][0] = ratio[1];
                        numLimits = 1;
                    } else {
                        int newRegion = 0;
                        for (int k=0; k < numLimits; k++) {
                            if ((ratio[1] < ratioList[0][k]) ||
                                (ratio[0] > ratioList[1][k])) {
                                newRegion += 1;
                            } else {
                                if (ratio[0] < ratioList[0][k]) {
                                    ratioList[0][k] = ratio[0];
                                }
                                if (ratio[1] > ratioList[1][k]) {
                                    ratioList[1][k] = ratio[1];
                                }
                            }
                        }

                        if (newRegion == numLimits) {
                            ratioList[0][numLimits] = ratio[0];
                            ratioList[1][numLimits] = ratio[1];
                            numLimits += 1;
                            if (numLimits >= NUMLIMITS_MAX) {
                                Fatal("numLimits (%d) exceeds max of %d", numLimits, NUMLIMITS_MAX);
                            }
                        }
                    }
                }  /* loop over intersecting particles */

                real8 factor = 0.0;

                for (j = 0; j < numLimits; j++) {
                    factor += (1.0 - ratioList[0][j]) * (1.0 - ratioList[0][j]) -
                              (1.0 - ratioList[1][j]) * (1.0 - ratioList[1][j]);
                }

                Bseg[0][0] = Bseg[0][0] * (1.0 - factor) + BsegAlt[0][0] * factor;
                Bseg[0][1] = Bseg[0][1] * (1.0 - factor) + BsegAlt[0][1] * factor;
                Bseg[0][2] = Bseg[0][2] * (1.0 - factor) + BsegAlt[0][2] * factor;
                Bseg[1][1] = Bseg[1][1] * (1.0 - factor) + BsegAlt[1][1] * factor;
                Bseg[1][2] = Bseg[1][2] * (1.0 - factor) + BsegAlt[1][2] * factor;
                Bseg[2][2] = Bseg[2][2] * (1.0 - factor) + BsegAlt[2][2] * factor;

            }  /* if (intersection != NULL) */

#endif

            totLength += mag;

            Btotal[0][0] += halfMag * Bseg[0][0];
            Btotal[0][1] += halfMag * Bseg[0][1];
            Btotal[0][2] += halfMag * Bseg[0][2];
            Btotal[1][1] += halfMag * Bseg[1][1];
            Btotal[1][2] += halfMag * Bseg[1][2];
            Btotal[2][2] += halfMag * Bseg[2][2];

        }  /* End loop over arms */


        Btotal[1][0] = Btotal[0][1];
        Btotal[2][0] = Btotal[0][2];
        Btotal[2][1] = Btotal[1][2];

/*
 *      It's possible this function was called for a node which only
 *      had zero length segments (during SplitSurfaceNodes() for example).
 *      If that is the case, just set the velocity to zero and return;
 */
        if (numNonZeroLenSegs == 0) {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

/*
 *      If we need to include inertial terms...
 */
        if (param->includeInertia) {

            tmp = totLength * massMult / param->deltaTT;

            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    massMatrix[i][j] = tmp * (i == j);
                    Btotal[i][j] += tmp * (i == j);
                }
            }

            Matrix33Vector3Multiply(massMatrix, vOld, tmp3);
        } else {
            tmp3[0] = 0.0;
            tmp3[1] = 0.0;
            tmp3[2] = 0.0;
        }

/*
 *      At this point we should check if the matrix is invertible and
 *      if it isn't, find the eigen values and eigen vectors of the drag
 *      matrix, then invert the drag matrix keeping zero eigen values as zero.
 */
        
        nForce[0] = node->fX + tmp3[0];
        nForce[1] = node->fY + tmp3[1];
        nForce[2] = node->fZ + tmp3[2];

        if ( M33_INVERSE(invDragMatrix, Btotal) < 0 )
        { Fatal("%s::%s(%d) : Cannot invert drag matrix!", __FILE__, __func__, __LINE__ ); }

        Matrix33Vector3Multiply(invDragMatrix, nForce, nVel);

        node->vX = nVel[0];
        node->vY = nVel[1];
        node->vZ = nVel[2];

        return(0);

}

