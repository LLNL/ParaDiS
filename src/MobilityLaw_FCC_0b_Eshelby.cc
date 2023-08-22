/**************************************************************************
 *
 *  Function    : Mobility_FCC_0b_Eshelby
 *  Author      : Sylvie Aubry from Tom Arsenlis' FCC_0b
 *  Description : This is alternate version of a generic mobility function
 *                for FCC materials.  It is very similar in function and
 *                structure to the BCC_glide mobility although glide planes
 *                rather than burgers vectors are used to identify junctions.
 *                Additionally, this function attempts to identify nodes
 *                that both have short segments and have reversed the 
 *                their direction of velocity.  When found, these nodes
 *                are artifically slowed for a single timestep in an effort
 *                to dampen motion of flickering/oscillating nodes.
 *
 *  Arguments:
 *          IN:  node        Pointer to the node for which to
 *                           calculate velocity
 *          IN/OUT: mobArgs  Structure containing additional
 *                           parameters for conveying information
 *                           to/from the mobility function.
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

#define NUMLIMITS_MAX (2*MAX_SEGPART_INTERSECTIONS+1)


int Mobility_FCC_0b_Eshelby(Home_t *home, Node_t *node, MobArgs_t *mobArgs)
{
    int     i, numNbrs, numNorms; 
    int     numNonZeroLenSegs=0;
    real8   dx, dy, dz;
    real8   nx, ny, nz, nx2, ny2, nz2; 
    real8   temp;
    real8   mag, mag2, halfMag, invMag;
    real8   invbMag2, bMag2, costheta, costheta2;
    real8   Bline, Bscrew, Bedge, Bglide, Bclimb, Beclimb;
    real8   Bscrew2, Beclimb2;
    real8   invBscrew2, invBedge2;
    real8   BlmBsc, BclmBsc, BglmBsc, BlmBecl;
    real8   shortSegCutoff;
    real8   eps = 1.0e-12;
    real8   tor=1.e-5;
    real8   shortSegRatio = 1.0;
    real8   nForce[3], nVel[3];
    real8   Bseg[3][3]= {{0.0, 0.0, 0.0},
			 {0.0, 0.0, 0.0},
			 {0.0, 0.0, 0.0}};
    real8   Btotal[3][3] = {{0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0}};
#ifdef ESHELBY
    real8   BsegAlt[3][3]= {{0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0}};
#endif
    real8   normCrystal[3][3];
    Node_t  *nbrNode;
    Param_t *param;

    param  = home->param;
    Bscrew = 1.0 / param->MobScrew;
    Bedge  = 1.0 / param->MobEdge;

    real8 (*invDragMatrix)[3]    =  mobArgs->invDragMatrix;
    real8 (*glideConstraints)[3] =  mobArgs->glideConstraints;
    int   *numGlideConstraints   = &mobArgs->numGlideConstraints;

    Matrix33_Zero(invDragMatrix);
    *numGlideConstraints= 0;

/*
 *  Climb is set this way to help with convergence and the fact that this
 *  is a glide restricted mobility
 */
    Beclimb    = 1.0 / param->MobClimb;
    
    /* Climb mobility for junction nodes */
    if (param->MobClimbJunc > 0.0) {
      if (node->numNbrs > 2)
	Beclimb = 1.0 / param->MobClimbJunc;
    }
       
    Bscrew2    = Bscrew * Bscrew;
    Beclimb2   = Beclimb * Beclimb;

    Bline      = 1.0 * MIN(Bscrew, Bedge);
    BlmBsc     = Bline - Bscrew;
    BlmBecl    = Bline - Beclimb; 

    invBscrew2 = 1.0 / (Bscrew*Bscrew);
    invBedge2  = 1.0 / (Bedge*Bedge);

#ifdef ESHELBY
    real8 pos1[3], pos2[3];
    real8 alpha      = param->MobEshelbyResist;
    
    real8 BscrewE    = alpha * Bscrew;
    real8 BedgeE     = alpha * Bedge;
    real8 BeclimbE   = alpha * Beclimb;
    
    real8 BscrewE2   = BscrewE  * BscrewE ;
    real8 BeclimbE2  = BeclimbE * BeclimbE;
    
    real8 BlineE     = 1.0e-2 * MIN(BscrewE, BedgeE);
    real8 BlmBscE    = BlineE - BscrewE;
    real8 BlmBeclE   = BlineE - BeclimbE;
    
    real8 invBscrewE2   = 1.0 / (BscrewE*BscrewE);
    real8 invBedgeE2    = 1.0 / (BedgeE *BedgeE );
#endif


    numNbrs = node->numNbrs;

    shortSegCutoff = 0.5 * param->minSeg;

/*
 *  If node is 'pinned' in ALL dimensions, or the node has any arms 
 *  with a burgers vector that has explicitly been set to be sessile (via
 *  control file inputs), the node may not be moved so just zero the velocity
 *  and return
 */
    if (HAS_ALL_OF_CONSTRAINTS(node->constraint, PINNED_NODE) ||
        NodeHasSessileBurg(home, node)) {

        node->vX = 0.0;
        node->vY = 0.0;
        node->vZ = 0.0;

        return(0);
    }

/*
 *  It's possible this function was called for a node which had only zero-
 *  length segments (during SplitSurfaceNodes() for example).  If that is
 *  the case, just set the velocity to zero and return.
 */
    for (i = 0; i < numNbrs; i++) {

        if ((nbrNode = GetNeighborNode(home, node, i)) == (Node_t *)NULL) {
            continue;
        }

        dx = node->x - nbrNode->x;
        dy = node->y - nbrNode->y;
        dz = node->z - nbrNode->z;

        ZImage(param, &dx, &dy, &dz);

/*
 *      This assumes that lengths are given in terms of burgers vector units
 */
        if ((dx*dx + dy*dy + dz*dz) > eps) {
            numNonZeroLenSegs++;
        }
    }

    if (numNonZeroLenSegs == 0) {
        node->vX = 0.0;
        node->vY = 0.0;
        node->vZ = 0.0;
        return(0);
    }

/*
 *  Get the normal for the first segment in both the lab and crystal frames
 */
    numNorms = 0; 

    glideConstraints[0][X] = node->nx[0];
    glideConstraints[0][Y] = node->ny[0];
    glideConstraints[0][Z] = node->nz[0]; 

    NormalizeVec(glideConstraints[0]);

    if (param->useLabFrame) {
        Matrix33Vector3Multiply(param->rotMatrixInverse, glideConstraints[0],
                                normCrystal[0]);
    } else {
        VECTOR_COPY(normCrystal[0], glideConstraints[0]);
    }

/*
 *  if this test passes then the dislocation segment is a junction
 *  segment which will be constrained to grow and shrink but not glide.
 *  Here glide planes are assumed to be of 111 type and any zero in the
 *  plane data is assumed to be associated with a non-glide plane
 */
    temp = normCrystal[0][X] *
           normCrystal[0][Y] *
           normCrystal[0][Z];

    if (fabs(temp) < eps) {

        numNorms = 1;
        nbrNode = GetNeighborNode(home, node, 0);

        glideConstraints[2][X] = node->x - nbrNode->x;
        glideConstraints[2][Y] = node->y - nbrNode->y;
        glideConstraints[2][Z] = node->z - nbrNode->z;

        ZImage(param, &glideConstraints[2][X], &glideConstraints[2][Y],
               &glideConstraints[2][Z]);

        NormalizeVec(glideConstraints[2]);
        cross(glideConstraints[2], glideConstraints[0], glideConstraints[1]);
    }
    

    if (numNbrs > 1) { 

        i = 1;

        while ((numNorms < 2) && (i < numNbrs)){

            nx = node->nx[i];
            ny = node->ny[i];
            nz = node->nz[i];

            Normalize(&nx, &ny, &nz);
    
            if (numNorms == 0) {

                temp = fabs(glideConstraints[0][X] * nx +
                            glideConstraints[0][Y] * ny +
                            glideConstraints[0][Z] * nz);

                if (fabs(1.0e0-temp) > tor) {

                    numNorms = 1;

                    nx -= temp * glideConstraints[0][X];
                    ny -= temp * glideConstraints[0][Y];
                    nz -= temp * glideConstraints[0][Z];

                    Normalize(&nx, &ny, &nz);

                    glideConstraints[1][X] = nx;
                    glideConstraints[1][Y] = ny; 
                    glideConstraints[1][Z] = nz;

                    cross(glideConstraints[0],
                          glideConstraints[1],
                          glideConstraints[2]);
                }

            } else {/* numNorms==1*/

                temp = glideConstraints[2][X] * nx +
                       glideConstraints[2][Y] * ny +
                       glideConstraints[2][Z] * nz;

                if (fabs(temp) > tor ) {
                    numNorms = 2;
                }
            }

/*
 *          check to see if the normal is a non-glide plane and then add to
 *          the checks so that it constrains the junction dislocation to
 *          only move along its line
 */
            temp = normCrystal[0][X] *
                   normCrystal[0][Y] *
                   normCrystal[0][Z];

            if (fabs(temp) < eps) {

                nbrNode = GetNeighborNode(home, node, i);

                dx = node->x - nbrNode->x;
                dy = node->y - nbrNode->y;
                dx = node->z - nbrNode->z;

                ZImage(param, &dx, &dy, &dz);
                Normalize(&dx, &dy, &dz);
                xvector(nx, ny, nz, dx, dy, dz, &nx2, &ny2, &nz2);

                if (numNorms == 0) {

                    temp = fabs(glideConstraints[0][X] * nx2 +
                                glideConstraints[0][Y] * ny2 +
                                glideConstraints[0][Z] * nz2);

                    if (fabs(1.0e0 - temp) > tor ) {

                        numNorms = 1;

                        nx2 -= temp * glideConstraints[0][X];
                        ny2 -= temp * glideConstraints[0][Y];
                        nz2 -= temp * glideConstraints[0][Z];

                        Normalize(&nx2, &ny2, &nz2);

                        glideConstraints[1][X] = nx2;
                        glideConstraints[1][Y] = ny2; 
                        glideConstraints[1][Z] = nz2;

                        cross(glideConstraints[0],
                              glideConstraints[1],
                              glideConstraints[2]);
                    }

                } else {/* numNorms==1*/

                    temp = glideConstraints[2][X] * nx2 +
                           glideConstraints[2][Y] * ny2 +
                           glideConstraints[2][Z] * nz2;

                    if (fabs(temp) > tor) {
                        numNorms = 2;
                    }
                }
            }

            i++;

        } /* end while */

    } /* end if (numNbrs > 1) */   

/*
 *      Glide constraints have been set, now initialize the velocity
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
#ifdef ESHELBY
        pos1[X] = node->x;
        pos1[Y] = node->y;
        pos1[Z] = node->z;
#endif

	if (numNorms < 2) {
            for (i = 0; i < numNbrs; i++) {

#ifdef ESHELBY
 
	        Matrix33_Zero(Bseg);	      
	        Matrix33_Zero(BsegAlt);
#endif	      
                real8 b[3], d[3], m[3], n[3], nCryst[3];

                b[X] = node->burgX[i];
                b[Y] = node->burgY[i];
                b[Z] = node->burgZ[i];

                bMag2 = (b[X]*b[X] + b[Y]*b[Y] + b[Z]*b[Z]);
                invbMag2 = 1.0 / bMag2;

/*
 *              Calculate the length of the arm and its tangent line direction
 */
                nbrNode = GetNeighborNode(home, node, i);

                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                d[X] = nbrNode->x - node->x;
                d[Y] = nbrNode->y - node->y;
                d[Z] = nbrNode->z - node->z;

                ZImage(param, &d[X], &d[Y], &d[Z]);

#ifdef ESHELBY
		pos2[X] = node->x + d[X];
		pos2[Y] = node->y + d[Y];
		pos2[Z] = node->z + d[Z];
#endif

		mag2 = d[X]*d[X] + d[Y]*d[Y] + d[Z]*d[Z];

/*
 *              If the segment is zero length (which can happen when
 *              the mobility function is being called from SplitMultiNodes())
 *              just skip the segment.
 */
                if (mag2 < eps) {
                    continue;
                }

                mag     = sqrt(mag2);
                halfMag = mag/2.0;
                invMag  = 1.0 / mag;

                d[X] *= invMag;
                d[Y] *= invMag;
                d[Z] *= invMag;

/*
 *              If the node has a very short segment, we need to flag it
 *              for later
 */
                if (mag < shortSegCutoff) {
                    shortSegRatio = MIN(shortSegRatio, mag/param->minSeg);
                }

/*
 *              Calculate how close to screw the arm is
 */
                costheta = (d[X]*b[X] + d[Y]*b[Y] + d[Z]*b[Z]);
                costheta2 = (costheta*costheta) * invbMag2;
    
/*
 *              arms not on [1 1 1] planes don't move as readily as
 *              other arms, so must be handled specially.
 *
 *              If needed, rotate a copy of the glide plane vector from the
 *              laboratory frame to the crystal frame.
 */
                n[X] = node->nx[i];
                n[Y] = node->ny[i];
                n[Z] = node->nz[i];

                if (param->useLabFrame) {
                    Matrix33Vector3Multiply(param->rotMatrixInverse, n, nCryst);
                } else {
                    VECTOR_COPY(nCryst, n);
                }

                if (fabs(nCryst[X] * nCryst[Y] * nCryst[Z]) < eps) {
                    if (numNbrs == 2) {
                        Bseg[0][0] = Beclimb;
                        Bseg[1][1] = Beclimb;
                        Bseg[2][2] = Beclimb;
                    } else {
                        Bseg[0][0] = (d[X]*d[X]*BlmBecl + Beclimb);
                        Bseg[0][1] = (d[X]*d[Y]*BlmBecl);
                        Bseg[0][2] = (d[X]*d[Z]*BlmBecl);
                        Bseg[1][1] = (d[Y]*d[Y]*BlmBecl + Beclimb);
                        Bseg[1][2] = (d[Y]*d[Z]*BlmBecl);
                        Bseg[2][2] = (d[Z]*d[Z]*BlmBecl + Beclimb);
                    }
                } else  {
/*
 *                  Arm is a regular glide arm, so build the drag matrix
 *                  assuming the dislocation is screw type
 */
                    Bseg[0][0] = (d[X]*d[X] * BlmBsc + Bscrew);
                    Bseg[0][1] = (d[X]*d[Y] * BlmBsc);
                    Bseg[0][2] = (d[X]*d[Z] * BlmBsc);
                    Bseg[1][1] = (d[Y]*d[Y] * BlmBsc + Bscrew);
                    Bseg[1][2] = (d[Y]*d[Z] * BlmBsc);
                    Bseg[2][2] = (d[Z]*d[Z] * BlmBsc + Bscrew);

/*
 *                  Now correct the drag matrix for dislocations that are
 *                  not screw
 */
                    if ((1.0 - costheta2) > eps) {

                        cross(n, d, m);

                        Bglide = sqrt(invBedge2 + (invBscrew2-invBedge2) *
                                      costheta2);
                        Bglide = 1.0 / Bglide;
                        Bclimb = sqrt(Beclimb2 + (Bscrew2 - Beclimb2) *
                                      costheta2);

#ifdef NAN_CHECK
                        if (isnan(Bglide) != 0) {
                            Fatal("Mobility_FCC_0b_Eshelby: Bglide = NaN\n"
                                "  Bglide = sqrt(invBedge2 + "
                                "(invBscrew2-invBedge2)*costheta2)\n"
                                "  where invBedge2 = %lf, invBscrew2 = %lf, "
                                "costheta2 = %lf", invBedge2, invBscrew2,
                                costheta2);
                        }

                        if (isnan(Bclimb) != 0) {
                            Fatal("Mobility_FCC_0b_Eshelby: Bclimb = NaN\n"
                                "  Bclimb = sqrt(Beclimb2 + "
                                "(Bscrew2-Beclimb2)*costheta2)\n"
                                "  where Beclimb2 = %lf, Bscrew2 = %lf, "
                                "costheta2 = %lf", Beclimb2, Bscrew2,
                                costheta2);
                        }
#endif
                        BclmBsc = Bclimb - Bscrew;
                        BglmBsc = Bglide - Bscrew;


                        Bseg[0][0] += (n[X]*n[X] * BclmBsc +
                                        m[X]*m[X] * BglmBsc);
                        Bseg[0][1] += (n[X]*n[Y] * BclmBsc +
                                        m[X]*m[Y] * BglmBsc);
                        Bseg[0][2] += (n[X]*n[Z] * BclmBsc +
                                        m[X]*m[Z] * BglmBsc);
                        Bseg[1][1] += (n[Y]*n[Y] * BclmBsc +
                                        m[Y]*m[Y] * BglmBsc);
                        Bseg[1][2] += (n[Y]*n[Z] * BclmBsc +
                                        m[Y]*m[Z] * BglmBsc);
                        Bseg[2][2] += (n[Z]*n[Z] * BclmBsc +
                                        m[Z]*m[Z] * BglmBsc);
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
                if (fabs(nCryst[X] * nCryst[Y] * nCryst[Z]) < eps) {
                    if (numNbrs == 2) {
                        BsegAlt[0][0] = BeclimbE;
                        BsegAlt[1][1] = BeclimbE;
                        BsegAlt[2][2] = BeclimbE;
                    } else {
                        BsegAlt[0][0] = (d[X]*d[X]*BlmBeclE + BeclimbE);
                        BsegAlt[0][1] = (d[X]*d[Y]*BlmBeclE);
                        BsegAlt[0][2] = (d[X]*d[Z]*BlmBeclE);
                        BsegAlt[1][1] = (d[Y]*d[Y]*BlmBeclE + BeclimbE);
                        BsegAlt[1][2] = (d[Y]*d[Z]*BlmBeclE);
                        BsegAlt[2][2] = (d[Z]*d[Z]*BlmBeclE + BeclimbE);
                    }
                } else  {
/*
 *                  Arm is a regular glide arm, so build the drag matrix
 *                  assuming the dislocation is screw type
 */
                    BsegAlt[0][0] = (d[X]*d[X] * BlmBscE + BscrewE);
                    BsegAlt[0][1] = (d[X]*d[Y] * BlmBscE);
                    BsegAlt[0][2] = (d[X]*d[Z] * BlmBscE);
                    BsegAlt[1][1] = (d[Y]*d[Y] * BlmBscE + BscrewE);
                    BsegAlt[1][2] = (d[Y]*d[Z] * BlmBscE);
                    BsegAlt[2][2] = (d[Z]*d[Z] * BlmBscE + BscrewE);

/*
 *                  Now correct the drag matrix for dislocations that are
 *                  not screw
 */
                    if ((1.0 - costheta2) > eps) {

                        cross(n, d, m);

                        Bglide = sqrt(invBedgeE2 + (invBscrewE2-invBedgeE2) *
                                      costheta2);
                        Bglide = 1.0 / Bglide;
                        Bclimb = sqrt(BeclimbE2 + (BscrewE2 - BeclimbE2) *
                                      costheta2);

#ifdef NAN_CHECK
                        if (isnan(Bglide) != 0) {
                            Fatal("Mobility_FCC_0b_Eshelby: BglideE = NaN\n"
                                "  BglideE = sqrt(invBedgeE2 + "
                                "(invBscrewE2-invBedgeE2)*costheta2)\n"
                                "  where invBedgeE2 = %lf, invBscrewE2 = %lf, "
                                "costheta2 = %lf", invBedgeE2, invBscrewE2,
                                costheta2);
                        }

                        if (isnan(Bclimb) != 0) {
                            Fatal("Mobility_FCC_0b_Eshelby: BclimbE = NaN\n"
                                "  BclimbE = sqrt(BeclimbE2 + "
                                "(BscrewE2-BeclimbE2)*costheta2)\n"
                                "  where BeclimbE2 = %lf, BscrewE2 = %lf, "
                                "costheta2 = %lf", BeclimbE2, BscrewE2,
                                costheta2);
                        }
#endif
                        BclmBsc = Bclimb - BscrewE;
                        BglmBsc = Bglide - BscrewE;


                        BsegAlt[0][0] += (n[X]*n[X] * BclmBsc + m[X]*m[X] * BglmBsc);
                        BsegAlt[0][1] += (n[X]*n[Y] * BclmBsc + m[X]*m[Y] * BglmBsc);
                        BsegAlt[0][2] += (n[X]*n[Z] * BclmBsc + m[X]*m[Z] * BglmBsc);
                        BsegAlt[1][1] += (n[Y]*n[Y] * BclmBsc + m[Y]*m[Y] * BglmBsc);
                        BsegAlt[1][2] += (n[Y]*n[Z] * BclmBsc + m[Y]*m[Z] * BglmBsc);
                        BsegAlt[2][2] += (n[Z]*n[Z] * BclmBsc + m[Z]*m[Z] * BglmBsc);
                    }
                }  /* End non-[0 0 1] arm */

		
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
                                                             inclusion->rotation, pos1, pos2, 1.0,
                                                             &ninc, ratio);
 
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

                for (int j = 0; j < numLimits; j++) {
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
 *          At this point we should check if the matrix is invertable and
 *          if it isn't, find the eigen values and eigen vectors of the drag
 *          matrix, then invert the drag matrix keeping zero eigen values
 *          as zero.
 *
 *          FIX ME!  For now, we're assuming the matrix is invertible.
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

	}

        return(0);
}
