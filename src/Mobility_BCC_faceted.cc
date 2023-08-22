/**************************************************************************
 *
 *      Module:       Mobility_BCC_faceted.c
 *      Description:  This mobility is intended for bcc crystals with
 *                    reduced mobility in certain line directions.
 *                    The mobility is linear in terms of the velocity
 *                    of a dislocation being proportional to the force
 *                    for every orientation; however the value of the
 *                    linear constant, Drag or Mobility, is different
 *                    for different line directions.  The mobility
 *                    is a planar mobility in that it restricts dislocations
 *                    to move on their glide planes identified by the normal
 *                    direction in the segment information.  
 *
 *                    Based on the Arsenlis matlab code mobbcc_facet.m
 *
 *                    The standard faceting within this mobility is
 *                    controlled by adding an additional drag to dislocation
 *                    oriented within a specific angular window of a
 *                    specified angle.  One period of the (1/2*cos(theta/)+1)
 *                    from -pi <theta< +pi is used to specify the additional
 *                    drag where theta/pi=(dislocation character angle -
 *                    window center) / window width.  Up to 10 windows each
 *                    with a unique center and width and additional drag
 *                    may be specified.
 *
 *                    A "cusping" variant of the faceting can be selected
 *                    which is a different way to describe the mobility
 *                    change across segment characters.  A 'cusp' shape will
 *                    have a flat background that changes smoothly to a
 *                    sharp tip at a specified angle.  It's similar to
 *                    the trench mobility in that both describe mobility
 *                    change across segment character say from mixed to
 *                    edge to screw.  However, it's different in the
 *                    sense that for trench mobility, one adds additional
 *                    terms at the trench to slow the mobility inside the
 *                    trench,  while for cusps the mobility changes smoothly
 *                    from fast to VERY slow values at the tip of the cusp.
 *                    
 *      Includes functions:
 *
 *            InitMob_BCC_faceted()
 *            Mobility_BCC_faceted()
 *                
 ***************************************************************************/
#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Mobility.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif

//#define M33_INVERSE(a,b) Matrix33_Inverse(a,b)
//#define M33_INVERSE(a,b) Matrix33_SVD_Inverse(a,b)
  #define M33_INVERSE(a,b) Matrix33_PseudoInverse(a,b)

static real8  Bglide, Bglide_p, Bclimb, Bline;
static real8  BlmBc, BlmBg, BcmBg;
static real8  Btrench[MAX_MOB_TRENCHES];
static real8  trenchAngle[MAX_MOB_TRENCHES];
static real8  trenchWidth[MAX_MOB_TRENCHES];


/*
 *      Some of the values we need never change throughout the
 *      execution, so set some static values once during initialization
 *      so they'll be preserved for subsequent calls.
 */
void InitMob_BCC_faceted(Home_t *home)
{
        int i;

        Bglide = 1.0 / home->param->MobGlide; /* floor on drag that applies */
                                              /* to all glide dislocations  */

        Bclimb = 1.0 / home->param->MobClimb; /* we actually use the climb  */
                                              /* mobility here to determine */
                                              /* the drag on junctions      */

        Bglide_p = Bglide * sqrt(2.0) / 4.0;  /* Only needed when 'cusps' */
                                              /* are enabled              */
        Bline  = Bglide;

        BlmBc = Bline - Bclimb;
        BlmBg = Bline - Bglide;
        BcmBg = Bclimb - Bglide;

        for (i = 0; i < home->param->MobNumTrenches; i++) {
            Btrench[i] = 1.0 / home->param->MobTrench[i];
        }

/*
 *      Convert trench angles and widths from degrees to radians
 */
        for (i = 0; i < home->param->MobNumTrenches; i++) {
            trenchAngle[i] = home->param->MobTrenchAngle[i] * M_PI / 180.0;
            trenchWidth[i] = home->param->MobTrenchWidth[i] * M_PI / 180.0;
        }

        return;
}


/**************************************************************************
 *
 *      Function:     Mobility_BCC_faceted
 *      Description:  This function calculates the velocity for a single
 *                    specified node.
 *
 *      Arguments:
 *          IN:  node         Pointer to the node for which to
 *                            calculate velocity
 *          IN/OUT: mobArgs   Structure containing additional
 *                            parameters for conveying information
 *                            to/from the mobility function.
 *
 *      Returns:  0 on success
 *                1 if velocity could not be determined
 *
 *************************************************************************/
int  Mobility_BCC_faceted(Home_t *home, Node_t *node, MobArgs_t *mobArgs)
{
        int     i, j;
        int     totLength;
        int     numNbrs, numNonZeroLenSegs;
        real8   trenchMin, trenchMax;
        real8   mag, halfMag, invMag, mag2;
        real8   theta, cth, factor;
        real8   tmp, tmp3[3];
        real8   burgNorm, burgVec[3];
        real8   vOld[3];
        real8   eps = 1.0e-06;
        real8   massMult=0.0, massMatrix[3][3];
        real8   gpNormal, glidePlane[3];
        real8   lDir[3], lDirTM[3][3], lDirNorm;
        real8   nDir[3], nDirTM[3][3];
        real8   mDir[3], mDirTM[3][3];
        real8   nForce[3], nVel[3];
        real8   Btotal[3][3] = {{0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0}};
        Param_t *param;
        Node_t  *nbrNode;

        real8 (*invBtotal)[3]        =  mobArgs->invDragMatrix;
     // real8 (*glideConstraints)[3] =  mobArgs->glideConstraints;
        int    *numGlideConstraints  = &mobArgs->numGlideConstraints;

        Matrix33_Zero(invBtotal);
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

        numNonZeroLenSegs = 0;

        numNbrs = node->numNbrs;

/*
 *      If we need to include inertial terms...
 */
        if (param->includeInertia) {
            massMult   = 0.25 * param->massDensity * 
                         (param->burgMag * param->burgMag);

            vOld[0] = node->oldvX;
            vOld[1] = node->oldvY;
            vOld[2] = node->oldvZ;
/*
 *          If necessary, rotate the old velocity vector from the
 *          laboratory frame to the crystal frame
 */
            if (param->useLabFrame) {
                real8 rotatedVec[3];
                Matrix33Vector3Multiply(param->rotMatrixInverse, vOld,
                                        rotatedVec);
                VECTOR_COPY(vOld, rotatedVec);
            }
        }


/*
 *      Loop over all arms of the node, adding each arm's contribution
 *      to the drag matrix.
 */
        totLength = 0.0;

        for (i = 0; i < numNbrs; i++) {

/*
 *          Calculate the length of the arm and its tangent line direction
 */
            nbrNode = GetNeighborNode(home, node, i);

            if (nbrNode == (Node_t *)NULL) {
                printf("WARNING: Neighbor not found at %s line %d\n",
                       __FILE__, __LINE__);
                continue;
            }

            lDir[X] = nbrNode->x - node->x;
            lDir[Y] = nbrNode->y - node->y;
            lDir[Z] = nbrNode->z - node->z;

            ZImage(param, &lDir[X], &lDir[Y], &lDir[Z]);

            mag2    = lDir[X]*lDir[X] +
                      lDir[Y]*lDir[Y] +
                      lDir[Z]*lDir[Z];
/*
 *          If the segment is zero length (which can happen when
 *          the mobility function is being called from SplitMultiNodes())
 *          just skip the segment.
 */
            if (mag2 < eps) {
                continue;
            }

            numNonZeroLenSegs++;

            burgVec[X] = node->burgX[i];
            burgVec[Y] = node->burgY[i];
            burgVec[Z] = node->burgZ[i];

/*
 *          If needed, rotate the burgers vector and line sense from the
 *          laboratory frame to the crystal frame.
 */
            if (param->useLabFrame) {
                real8 bRot[3], lDirRot[3];

                Matrix33Vector3Multiply(param->rotMatrixInverse, burgVec, bRot);
                Matrix33Vector3Multiply(param->rotMatrixInverse, lDir, lDirRot);

                VECTOR_COPY(burgVec, bRot);
                VECTOR_COPY(lDir, lDirRot);
            }

            burgNorm = Normal(burgVec);

            mag = sqrt(mag2);
            halfMag = 0.5 * mag;
            invMag  = 1.0 / mag;

            totLength += mag;

            lDir[X] *= invMag;
            lDir[Y] *= invMag;
            lDir[Z] *= invMag;

            lDirNorm = Normal(lDir);

            glidePlane[X] = node->nx[i];
            glidePlane[Y] = node->ny[i];
            glidePlane[Z] = node->nz[i];

/*
 *          If needed, rotate the glide plane from the laboratory frame
 *          to the crystal frame.
 */
            if (param->useLabFrame) {
                real8 glidePlaneRot[3];

                Matrix33Vector3Multiply(param->rotMatrixInverse, glidePlane,
                                        glidePlaneRot);

                VECTOR_COPY(glidePlane, glidePlaneRot);
            }

            gpNormal = Normal(glidePlane);

            Vec3TransposeAndMult(lDir, lDirTM);

            if ((fabs(burgVec[X] * burgVec[Y] * burgVec[Z]) /
                 (burgNorm*burgNorm*burgNorm)) < eps) {
/*
 *              Junctions are treated such that they have climb mobility
 *              perpendicular to their line directions, and glide/line
 *              mobility in their line direction.
 */
                Btotal[0][0] += 2 * mag * (BlmBc * lDirTM[0][0] + Bclimb);
                Btotal[0][1] += 2 * mag * (BlmBc * lDirTM[0][1]);
                Btotal[0][2] += 2 * mag * (BlmBc * lDirTM[0][2]);
                Btotal[1][1] += 2 * mag * (BlmBc * lDirTM[1][1] + Bclimb);
                Btotal[1][2] += 2 * mag * (BlmBc * lDirTM[1][2]);
                Btotal[2][2] += 2 * mag * (BlmBc * lDirTM[2][2] + Bclimb);

            } else {
                nDir[X] = glidePlane[X] / gpNormal;
                nDir[Y] = glidePlane[Y] / gpNormal;
                nDir[Z] = glidePlane[Z] / gpNormal;
                Vec3TransposeAndMult(nDir, nDirTM);
                cross(lDir, nDir, mDir);

                if (!param->MobFacetedUseCusps) {
/*
 *                  Not using cusps, so treat trenches in the standard
 *                  fashion.  (Bglide is separate from Btrench in the trench
 *                  mobility)
 */
                    Btotal[0][0] += halfMag * (BlmBg * lDirTM[0][0] +
                                               BcmBg * nDirTM[0][0] + Bglide);
                    Btotal[0][1] += halfMag * (BlmBg * lDirTM[0][1] +
                                               BcmBg * nDirTM[0][1]);
                    Btotal[0][2] += halfMag * (BlmBg * lDirTM[0][2] +
                                               BcmBg * nDirTM[0][2]);
                    Btotal[1][1] += halfMag * (BlmBg * lDirTM[1][1] +
                                               BcmBg * nDirTM[1][1] + Bglide);
                    Btotal[1][2] += halfMag * (BlmBg * lDirTM[1][2] +
                                               BcmBg * nDirTM[1][2]);
                    Btotal[2][2] += halfMag * (BlmBg * lDirTM[2][2] +
                                               BcmBg * nDirTM[2][2] + Bglide);
                } else {
/*
 *                  We're using cusps rather than the standard trenches.
 *                  Bglide will be included in tmp below together with
 *                  Btrench in cusp mobility.  For cusps, we don't
 *                  subtract Bglide from Bline (to get BlmBg) and Bclimb
 *                  (to get BcmBg). Also, make sure climb mobility is smaller
 *                  than the slowest trench mobility at screw/edge
 *                  orientation
 */
                    Btotal[0][0] += halfMag * (Bclimb * nDirTM[0][0]);
                    Btotal[0][1] += halfMag * (Bclimb * nDirTM[0][1]);
                    Btotal[0][2] += halfMag * (Bclimb * nDirTM[0][2]);
                    Btotal[1][1] += halfMag * (Bclimb * nDirTM[1][1]);
                    Btotal[1][2] += halfMag * (Bclimb * nDirTM[1][2]);
                    Btotal[2][2] += halfMag * (Bclimb * nDirTM[2][2]);
    
                }    

/*
 *              <cth> *should* be in the range  -1.0 <= cth <= +1.0
 *              however, due to machine precision issues, the calculation
 *              used can sometimes result in values minutely outside that
 *              range, so explicitly force the result back inside the range.
 */
                cth = DotProduct(burgVec, lDir) / (burgNorm * lDirNorm);
                cth = MIN(MAX(cth, -1.0), 1.0);

                cross(burgVec, lDir, tmp3);
                tmp = DotProduct(tmp3, nDir);
                theta = acos(cth) * Sign(tmp);

                Vec3TransposeAndMult(mDir, mDirTM);

                if (theta > (M_PI/2.0)) {
                    theta -= M_PI;
                } else if (theta <= -(M_PI/2.0)) {
                    theta += M_PI;
                }

                if (!param->MobFacetedUseCusps) {
/*
 *                  Use standard trench mobility
 */ 
                    for (j = 0; j < param->MobNumTrenches; j++) {

                        trenchMin = trenchAngle[j] - trenchWidth[j];
                        trenchMax = trenchAngle[j] + trenchWidth[j];

                        if (trenchMin <= -(M_PI/2.0)) {
                            trenchMin += M_PI;
                        }

                        if (trenchMax > (M_PI/2.0)) {
                            trenchMax -= M_PI;
                        }

                        if (     ((theta > trenchMin) && (theta < trenchMax)) 
                             || (((theta > trenchMin) || (theta < trenchMax)) && (trenchMin > trenchMax)) ) 
                        {
                            tmp = sin((theta - trenchAngle[j]) *
                                      (0.5e0*M_PI / trenchWidth[j]));
                            tmp *= tmp;
                            factor = 1.0e0 - pow(tmp, 10.0);
                
                            tmp = halfMag * Btrench[j] * factor;

                            Btotal[0][0] += tmp * (lDirTM[0][0] * 1.0e-03 +
                                                   mDirTM[0][0]);
                            Btotal[0][1] += tmp * (lDirTM[0][1] * 1.0e-03 +
                                                   mDirTM[0][1]);
                            Btotal[0][2] += tmp * (lDirTM[0][2] * 1.0e-03 +
                                                   mDirTM[0][2]);
                            Btotal[1][1] += tmp * (lDirTM[1][1] * 1.0e-03 +
                                                   mDirTM[1][1]);
                            Btotal[1][2] += tmp * (lDirTM[1][2] * 1.0e-03 +
                                                   mDirTM[1][2]);
                            Btotal[2][2] += tmp * (lDirTM[2][2] * 1.0e-03 +
                                                   mDirTM[2][2]);
                        }

                    }  /* end for (j = 0; j < param->MobNumTrenches; ...) */

                } else {
/*
 *                  Use the cusp mobility
 */ 
                    if (theta < 0.0 ) {
                        theta += M_PI;
                    }

                    if (theta > (M_PI/2.0)) {
                        theta = M_PI - theta; 
                    }

                    for (j = 0; j < param->MobNumTrenches; j++) {


                        tmp = sin(theta-trenchAngle[j]);  

                        factor = 1.0 / sqrt(1.0/Btrench[j]/Btrench[j] +
                                            1.0/Bglide_p/Bglide_p*tmp*tmp); 

                        tmp = halfMag * factor;

                        Btotal[0][0] += tmp * (lDirTM[0][0] * 1.0e-03 +
                                               mDirTM[0][0]);
                        Btotal[0][1] += tmp * (lDirTM[0][1] * 1.0e-03 +
                                               mDirTM[0][1]);
                        Btotal[0][2] += tmp * (lDirTM[0][2] * 1.0e-03 +
                                               mDirTM[0][2]);
                        Btotal[1][1] += tmp * (lDirTM[1][1] * 1.0e-03 +
                                               mDirTM[1][1]);
                        Btotal[1][2] += tmp * (lDirTM[1][2] * 1.0e-03 +
                                               mDirTM[1][2]);
                        Btotal[2][2] += tmp * (lDirTM[2][2] * 1.0e-03 +
                                               mDirTM[2][2]);

                    }  /* end for (j = 0; j < param->MobNumTrenches; ...) */

                }  /* end code for cusps */

            }  /* end else */

        }  /* end for (i = 0; i < numNbrs; ...) */

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

        nForce[0] = node->fX;
        nForce[1] = node->fY;
        nForce[2] = node->fZ;

/*
 *      If needed, rotate the force vector from the laboratory frame to the
 *      crystal frame
 */
        if (param->useLabFrame) {
            real8 rotForce[3];
            Matrix33Vector3Multiply(param->rotMatrixInverse, nForce, rotForce);
            VECTOR_COPY(nForce, rotForce);
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

            nForce[0] += tmp3[0];
            nForce[1] += tmp3[1];
            nForce[2] += tmp3[2];
        }

/*
 *      At this point we should check if the matrix is invertable and
 *      if it isn't, find the eigen values and eigen vectors of the drag
 *      matrix, then invert the drag matrix keeping zero eigen values as zero.
 *
 *      FIX ME!  For now, we're assuming the matrix is invertable.
 */
        if ( M33_INVERSE(invBtotal, Btotal) < 0 ) 
        { Fatal("%s::%s(%d) : Cannot invert matrix for node (%d,%d)!", __FILE__, __func__, __LINE__, node->myTag.domainID, node->myTag.index); }

        Matrix33Vector3Multiply(invBtotal, nForce, nVel);

/*
 *      If needed, rotate the velocity vector back to the laboratory frame
 *      from the crystal frame
 */
        if (param->useLabFrame) {
            real8 vTemp[3];
            Matrix33Vector3Multiply(param->rotMatrix, nVel, vTemp);
            VECTOR_COPY(nVel, vTemp);
        }


#ifdef NAN_CHECK
        if ((isnan(nVel[0]) || isinf(nVel[0])) ||
            (isnan(nVel[1]) || isinf(nVel[1])) ||
            (isnan(nVel[2]) || isinf(nVel[2]))) {
            Fatal("(%d,%d) v = (%e %e %e) in BCC_faceted",
                  node->myTag.domainID, node->myTag.index,
                  nVel[0], nVel[1], nVel[2]);
        }
#endif
        node->vX = nVel[0];
        node->vY = nVel[1];
        node->vZ = nVel[2];

        return(0);
}
