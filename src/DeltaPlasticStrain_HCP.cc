/***************************************************************************
 *   
 *      Author:       Sylvie Aubry
 *
 *      Function:     DeltaPlasticStrain_HCP
 *
 *      Description:  Calculate the plastic strain increment
 *                    by dislocation motion in HCP materials.
 *
 *      Public functions:
 *          DeltaPlasticStrain_HCP()
 *
 *      Private functions:
 *          InitDeltaPlasticStrain_HCP()
 *
 ***************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Node.h"
#include "Util.h"
#include "Mobility.h"

/***************************************************************************
 *
 *      Function:    InitDeltaPlasticStrain_HCP
 *
 *      Description: Several of the arrays used during the delta plastic
 *                   strain calculations are dynamically calculated, but
 *                   once set, never need to be changed again for the
 *                   duration of the simulation.  This function should be
 *                   invoked during the first call to DeltaPlasticStrain_HCP()
 *                   to initialize those arrays, and never used again
 *                   after that.
 *
 *                   NOTE: upon return to the caller, the <nanv>, <tanv>
 *                         and <bnorm> arrays are in the laboratory frame
 *                         not the crystalographic frame.
 *
 ***************************************************************************/
void InitDeltaPlasticStrain_HCP(Home_t *home,
                real8 nanv[HCP_NUM_GLIDE_PLANES][3][HCP_NUM_GLIDE_BURG],
                real8 tanv[HCP_NUM_GLIDE_PLANES][3][HCP_NUM_GLIDE_BURG],
                real8 bnorm[HCP_NUM_GLIDE_BURG][3])
{
        int       i, j, k, m;
        int      *numGlissilePlanesPerBurg, *burgFirstPlaneIndex;
        real8   (*burgList)[3], (*planeList)[3];
        Param_t  *param = home->param;

        burgList = home->burgData.burgList;
        planeList = home->burgData.planeList;
        burgFirstPlaneIndex = home->burgData.burgFirstPlaneIndex;
        numGlissilePlanesPerBurg = home->burgData.numGlissilePlanesPerBurg;

        for (k = 0; k < HCP_NUM_GLIDE_BURG; k++) {
            int start, end;

            for (j = 0; j < 3; j++)  {
                bnorm[k][j] = burgList[k][j];
            }
            NormalizeVec(bnorm[k]);

            m = 0;

            start = burgFirstPlaneIndex[k];
            end = burgFirstPlaneIndex[k] + numGlissilePlanesPerBurg[k];

            for (i = start; i < end; i++) {
                for (j = 0; j < 3; j++) {
                    nanv[m][j][k] = planeList[i][j];
                }
                m += 1;
            }
        }


/*
 *      Get edge directions
 */
        for (i = 0; i < HCP_NUM_GLIDE_BURG; i++) {
            real8 bb[3], pp[3], ee[3];

            VECTOR_COPY(bb, bnorm[i]);

            for (j = 0; j < HCP_NUM_GLIDE_PLANES; j++) {

                pp[0] = nanv[j][0][i];
                pp[1] = nanv[j][1][i];
                pp[2] = nanv[j][2][i];

                cross(bb, pp, ee);
                NormalizeVec(ee);

                tanv[j][0][i] = ee[0];
                tanv[j][1][i] = ee[1];
                tanv[j][2][i] = ee[2];
            }
        }

/*
 *      Burgers vector 9 has 3 planes only. Initialize normal
 */
        nanv[3][0][9]=0;
        nanv[3][1][9]=0;
        nanv[3][2][9]=0;

/*
 *      Burgers vector 9 has 3 planes only. Initialize edge
 */
        tanv[3][0][9]=0;
        tanv[3][1][9]=0;
        tanv[3][2][9]=0;

#if 0
/*
 *      Some stuff for debugging
 */
        printf("cOVERa=%f\n", home->param->cOVERa);

        for (i = 0; i < HCP_NUM_GLIDE_BURG; i++) {
            printf("bnorm[%d]=%f %f %f\n", i,
                   bnorm[i][0], bnorm[i][1], bnorm[i][2]);
        }
        printf("\n");

        for (i = 0; i < HCP_NUM_GLIDE_BURG; i++) {
            for (j = 0; j < HCP_NUM_GLIDE_PLANES; j++) {
                printf("nanv[%d][%d]=%f %f %f\n", i, j,
                       nanv[j][0][i], nanv[j][1][i], nanv[j][2][i]);
            }
            printf("\n");
        }

        for (i = 0; i < HCP_NUM_GLIDE_BURG; i++) {
            for (j = 0; j < HCP_NUM_GLIDE_PLANES; j++) {
                printf("tanv[%d][%d]=%f %f %f\n", i, j,
                       tanv[j][0][i], tanv[j][1][i], tanv[j][2][i]);
            }
            printf("\n");
        }
#endif

/*
 *      If necessary, rotate the nanv, tanv and bnorm arrays from the
 *      crystal frame to the lab frame
 */
        if (param->useLabFrame) {
            int plane, bIndex;

            for (plane = 0; plane < HCP_NUM_GLIDE_PLANES; plane++) {

                for (bIndex = 0; bIndex < HCP_NUM_GLIDE_BURG; bIndex++) {
                    real8 nanvCryst[3], nanvLab[3], tanvCryst[3], tanvLab[3];

                    nanvCryst[X] = nanv[plane][X][bIndex];
                    nanvCryst[Y] = nanv[plane][Y][bIndex];
                    nanvCryst[Z] = nanv[plane][Z][bIndex];

                    tanvCryst[X] = tanv[plane][X][bIndex];
                    tanvCryst[Y] = tanv[plane][Y][bIndex];
                    tanvCryst[Z] = tanv[plane][Z][bIndex];

                    Matrix33Vector3Multiply(param->rotMatrix,nanvCryst,nanvLab);
                    Matrix33Vector3Multiply(param->rotMatrix,tanvCryst,tanvLab);

                    nanv[plane][X][bIndex] = nanvLab[X];
                    nanv[plane][Y][bIndex] = nanvLab[Y];
                    nanv[plane][Z][bIndex] = nanvLab[Z];

                    tanv[plane][X][bIndex] = tanvLab[X];
                    tanv[plane][Y][bIndex] = tanvLab[Y];
                    tanv[plane][Z][bIndex] = tanvLab[Z];
                }
            }

            for (bIndex = 0; bIndex < HCP_NUM_GLIDE_BURG; bIndex++) {
                real8 bnormLab[3];

                Matrix33Vector3Multiply(param->rotMatrix,bnorm[bIndex],bnormLab);
                VECTOR_COPY(bnorm[bIndex], bnormLab);
            }
        }

        return;
}


void DeltaPlasticStrain_HCP(Home_t *home)
{
        int          i, j;
        real8        factor0, bmagsq;
        real8        gdspn[6], gdstn[6];
        real8        dstn[3][3], dspn[3][3];
        real8        Ltot[HCP_NUM_GLIDE_BURG][17];
        real8        areaSwept[HCP_NUM_GLIDE_BURG][9];
        real8        (*burgList)[3];
        Param_t      *param;
        static int   initDone = 0;
        static real8 bnorm[HCP_NUM_GLIDE_BURG][3];
//      static real8 burgMagLocal[HCP_NUM_GLIDE_BURG];
        static real8 nanv[HCP_NUM_GLIDE_PLANES][3][HCP_NUM_GLIDE_BURG];
        static real8 tanv[HCP_NUM_GLIDE_PLANES][3][HCP_NUM_GLIDE_BURG];
     
//   Ltot components:
//   Ltot[*][0]  total edge_density (based on total length - screw length)
//   Ltot[*][1]  edge_density plane 1
//   Ltot[*][2]  edge_density plane 2
//   Ltot[*][3]  edge_density plane 3
//   Ltot[*][4]  edge_density plane 4
//   Ltot[*][5]  edge_density all other planes
//   Ltot[*][6]  total screw_density (based on screw length)
//   Ltot[*][7]  screw_density plane 1
//   Ltot[*][8]  screw_density plane 2
//   Ltot[*][9]  screw_density plane 3
//   Ltot[*][10] screw_density plane 4
//   Ltot[*][11] screw_density all other planes
//   Ltot[*][12] total_density plane 1 = sqrt(edgeDen1^2-screwDen1^2)
//   Ltot[*][13] total_density plane 2 = sqrt(edgeDen2^2-screwDen2^2)
//   Ltot[*][14] total_density plane 3 = sqrt(edgeDen3^2-screwDen3^2)
//   Ltot[*][15] total_density plane 4 = sqrt(edgeDen4^2-screwDen4^2)
//   Ltot[*][16] total_density all other planes = sqrt(edgeOther^2-screwOther^2)
//

        param = home->param;
        bmagsq = param->burgMag * param->burgMag;

        Matrix33_Zero(dstn);
        Matrix33_Zero(dspn);

        memset(Ltot, 0, sizeof(Ltot));
        memset(areaSwept, 0, sizeof(areaSwept));

        param->disloDensity = 0.0;

        for (i = 0; i < param->numBurgGroups; i++) {
            param->partialDisloDensity[i] = 0.0;
        }

        for (i=0; i < 6; i++) {
            param->delpStrain[i] = 0.0;                                         
        }                                                                

/*
 *      Set some pointers to the burgers vector and glide plane lists
 *      that were created during initialization.
 */
        burgList = home->burgData.burgList;

/*
 *      The first time we enter this function, initialize a few
 *      arrays that will remain static for the duration of the
 *      simulation once they are set.
 */
        if (!initDone++) {
            InitDeltaPlasticStrain_HCP(home, nanv, tanv, bnorm);
        }

/*
 *      Keep a list of the magnitude of the glide Burgers vectors
 *      of the simulation. They are used for comparison purposes only between
 *      gdstn and pstn at the end of this file.
 */
    //  for (i=0; i< HCP_NUM_GLIDE_BURG; i++) {
    //      burgMagLocal[i] = 0;
    //  }


/*
 *      Loop through all the nodes calculating the flux contribution
 *      from each of the node's segments.  We want the loop threaded,
 *      but we also need to sum thread-specific totals when we're
 *      done, so we can't just automatically distribute the loop
 *      iterations.  Instead, we'll set up a threaded execution block
 *      where each thread calculates it's own loop indices and then
 *      safely adds its contributions to the summed totals in a 'critical'
 *      section after its loop iterations have been completed.
 */
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        int    k, l;
        int    nodeIndex;
        int    threadID, threadIterStart, threadIterEnd;
        real8  thread_disloDensity;
        real8  *thread_partialDensity;
        real8  thread_dstn[3][3], thread_dspn[3][3];
        real8  thread_Ltot[HCP_NUM_GLIDE_BURG][17];
        real8  thread_areaSwept[HCP_NUM_GLIDE_BURG][9];

        thread_disloDensity = 0.0;
        thread_partialDensity = (real8 *)calloc(1, param->numBurgGroups *
                                                sizeof(real8));

        Matrix33_Zero(thread_dstn);
        Matrix33_Zero(thread_dspn);

        memset(thread_Ltot, 0, sizeof(thread_Ltot));
        memset(thread_areaSwept, 0, sizeof(thread_areaSwept));

        GetThreadIterationIndices(home->newNodeKeyPtr, &threadID,
                                  &threadIterStart, &threadIterEnd);

        for (nodeIndex=threadIterStart; nodeIndex<threadIterEnd; nodeIndex++) {
            int    m;
            int    numSeg, segIndex;
            real8  localstrain[3][3];
            Node_t *node;

            if ((node = home->nodeKeys[nodeIndex]) == (Node_t *)NULL) {
                continue;
            }

            numSeg = node->numNbrs;

            Matrix33_Zero(localstrain);
        
/*
 *          Loop over every segment attached to the node
 */        
            for (segIndex = 0; segIndex < numSeg; segIndex++) {
               int    index2, index3;
               int    bIndex, burgGroup;
               int    offset1, offset2;
               int    aOffset1, aOffset2;
               int    planeIndex;
               real8  bx, by, bz;
               real8  ex, ey, ez;
               real8  nx, ny, nz;
               real8  hhx, hhy, hhz;
               real8  delxnode, delynode, delznode;
               real8  delxnbr, delynbr, delznbr;
               real8  deltax1x, deltax1y, deltax1z;
               real8  deltax2x, deltax2y, deltax2z;
               real8  tmpx, tmpy, tmpz;
               real8  deltaDensity;
               real8  edgeDensity;
               real8  maxVal;
               real8  size;
               real8  sbsign, sb;
               real8  segd, segd2;
               real8  seg[3], seg2[3];
               real8  segleft[3], segleft2[3], qs[3];
               real8  Ltemp2[2];
               real8  Ltemp[HCP_NUM_GLIDE_PLANES];
               real8  stemp[HCP_NUM_GLIDE_PLANES];
               real8  bCryst[3];
               real8  zeta[3];
               real8  vec1[3], vec2[3], outVec[2];
               real8  thread_dyad[3][3];
               real8  segPlane[3];
               Node_t *nbr;

               nbr = GetNeighborNode(home, node, segIndex);

               if (nbr == (Node_t *)NULL) {
                   printf("WARNING: Neighbor not found at %s line %d\n",
                          __FILE__, __LINE__);
                   continue;
               }

               /* Avoid double counting */

               if (OrderNodes(node, nbr) != -1) {
                   continue;
               }

/*
 *             This check is in support of some local site-specific
 *             post-processing only and all *normal* ParaDiS simulations
 *             will simply continue on through this check.
 */
               if (param->ignoreNonLoopFlux) {
                   if (!HAS_ANY_OF_CONSTRAINTS(node->constraint, LOOP_NODE) ||
                       !HAS_ANY_OF_CONSTRAINTS(nbr->constraint,  LOOP_NODE)) {
                       continue;
                   }
               }
/*
 *             For later, we need a copy of the burgers vector guaranteed
 *             to be in the crystal frame
 */
               bx = node->burgX[segIndex];
               by = node->burgY[segIndex];
               bz = node->burgZ[segIndex];

               if (param->useLabFrame) {
                   real8 bLab[3] = {bx, by, bz};
                   Matrix33Vector3Multiply(param->rotMatrixInverse,bLab,bCryst);
               } else {
                   bCryst[X] = bx;
                   bCryst[Y] = by;
                   bCryst[Z] = bz;
               }

               ex = nbr->x - node->oldx; 
               ey = nbr->y - node->oldy; 
               ez = nbr->z - node->oldz; 

               hhx = node->x - nbr->oldx;
               hhy = node->y - nbr->oldy;
               hhz = node->z - nbr->oldz;

               ZImage(param, &ex, &ey, &ez);
               ZImage(param, &hhx, &hhy, &hhz);

               nx = (ey*hhz - ez*hhy) * 0.5;
               ny = (ez*hhx - ex*hhz) * 0.5;
               nz = (ex*hhy - ey*hhx) * 0.5;

               thread_dyad[0][0] = bx*nx;
               thread_dyad[0][1] = bx*ny;
               thread_dyad[0][2] = bx*nz;
               thread_dyad[1][0] = by*nx;
               thread_dyad[1][1] = by*ny;
               thread_dyad[1][2] = by*nz;
               thread_dyad[2][0] = bz*nx;
               thread_dyad[2][1] = bz*ny;
               thread_dyad[2][2] = bz*nz;
           
               for (k = 0; k < 3; k++){
                  for (l = 0; l < 3; l++){
                      real8 tmpDstn;
                      tmpDstn = (thread_dyad[l][k]+thread_dyad[k][l]) *
                                0.5/param->simVol;
                      thread_dstn[l][k] += tmpDstn;
                      thread_dspn[l][k] += (thread_dyad[l][k] -
                                            thread_dyad[k][l]) *
                                           0.5/param->simVol;
                      localstrain[l][k] += tmpDstn;
                  }
               }

               delxnode = node->x - node->oldx;
               delynode = node->y - node->oldy;
               delznode = node->z - node->oldz;

               delxnbr = nbr->x - nbr->oldx;
               delynbr = nbr->y - nbr->oldy;
               delznbr = nbr->z - nbr->oldz;

               ZImage(param, &delxnode, &delynode, &delznode);
               ZImage(param, &delxnbr, &delynbr, &delznbr);

               tmpx = nbr->x - node->x;
               tmpy = nbr->y - node->y;
               tmpz = nbr->z - node->z;

               ZImage(param, &tmpx, &tmpy, &tmpz);

               zeta[0] = tmpx;
               zeta[1] = tmpy;
               zeta[2] = tmpz;

               size = sqrt(DotProduct(zeta, zeta));

               deltaDensity = size * param->burgVolFactor;
               thread_disloDensity += deltaDensity;
               Normalize(&nx, &ny, &nz);
           
               deltax1x= node->x - node->oldx;
               deltax1y= node->y - node->oldy;
               deltax1z= node->z - node->oldz;

               ZImage(param, &deltax1x, &deltax1y, &deltax1z);

               deltax2x = nbr->x - nbr->oldx;
               deltax2y = nbr->y - nbr->oldy;
               deltax2z = nbr->z - nbr->oldz;

               ZImage(param, &deltax2x, &deltax2y, &deltax2z);

               seg[0] = nbr->x - node->x;                                  
               seg[1] = nbr->y - node->y;    
               seg[2] = nbr->z - node->z;

               ZImage(param, &seg[0], &seg[1], &seg[2]);

               seg2[0]= nbr->oldx - node->oldx;
               seg2[1]= nbr->oldy - node->oldy;
               seg2[2]= nbr->oldz - node->oldz;

               ZImage(param, &seg2[0], &seg2[1], &seg2[2]);

/* 
 *             Get the index in the reference Burgers vectors list
 *             of current Burgers vector (bx, by, bz).  We only
 *             Check the 1st 10 glide burgers vectors here.
 */
               GetBurgIndexHCP(bCryst, 0, HCP_NUM_GLIDE_BURG,
                               burgList, &bIndex);

/*
 *             For HCP mobility types we track density of the following
 *             groupings of burgers vectors.  Any burgers vectors not
 *             mentioned are ignored.
 *
 *                 group #     burgers vector types
 *                     0       a/2[-1 sqrt(3) 0]    a/2[1 -sqrt(3) 0]
 *                     1          [a 0 0]              [-a 0 0]
 *                     2       a/2[1 sqrt(3) 0 ]    a/2[-1 -sqrt(3) 0 ]
 *                     3       1/2[-a sqrt(3)a 2c]  1/2[a -sqrt(3)a -2c]
 *                     4          [a 0 c]              [-a 0 -c]
 *                     5       1/2[a sqrt(3)a -2c]  1/2[-a -sqrt(3)a 2c]
 *                     6       1/2[-a sqrt(3)a -2c] 1/2[a -sqrt(3)a 2c]
 *                     7          [-a 0 c]             [a 0 -c]
 *                     8       1/2[a sqrt(3) 2c]    1/2[-a -sqrt(3) -2c]
 *                     9          [0 0 c]              [0 0 -c]
 *                    10       all others
 */

	       if (bIndex < 0) 
		 burgGroup = 10;
	       else
		 burgGroup = bIndex;


/*
 *             And, if the burger's vector fell into one of groups of
 *             burgers vectors whose density we're tracking, increment
 *             that group's density value.
 */
               if (burgGroup >= 0) {
                   thread_partialDensity[burgGroup] += deltaDensity;
               }

/*
 *             We're really only interested in the glide Burgers vectors (the
 *             first 10 Burgers vectors), so skip the rest.
 */
	       if (bIndex < 0) {
                   continue;
               }
       
/* 
 *             Keep the Burgers vector magnitude for comparisons done later.
 */
	   //  burgMagLocal[bIndex] = sqrt(bx*bx + by*by + bz*bz);

/*
 *             Compute the plastic strain rate directly
 */

               sbsign = bnorm[bIndex][0]*bx +
                        bnorm[bIndex][1]*by +
                        bnorm[bIndex][2]*bz;

               sb = ((sbsign < 0) ? -1.0 : 1.0);

/*
 *             segx vector and delta2 vector defined above
 *
 *             segd is the dotted value of seg length with
 *             burgers vector => Screw density
 *
 *             segleft is the total edge density for the segment
 */
               segd = bnorm[bIndex][0]*seg[0] +
                      bnorm[bIndex][1]*seg[1] +
                      bnorm[bIndex][2]*seg[2];

               // Update total screw length for the burgers vector
               thread_Ltot[bIndex][6] += fabs(segd);

               segleft[0] = seg[0] - segd * bnorm[bIndex][0];
               segleft[1] = seg[1] - segd * bnorm[bIndex][1];
               segleft[2] = seg[2] - segd * bnorm[bIndex][2];

/*
 *             For the first 10 burgers vectors find the screw density
 *             associated with the glissile glide planes.  Screw density
 *             for all other planes is lumped together into a single
 *             value.
 *
 *             Also track the total density per plane.
 */

               if (bIndex < 10) {
                   int    numGlissileGlidePlanes;
                   int   *numGlissilePlanesPerBurg;

                   segPlane[0] = node->nx[segIndex];
                   segPlane[1] = node->ny[segIndex];
                   segPlane[2] = node->nz[segIndex];

                   FindGlidePlaneHCP(home, bIndex, segPlane, &planeIndex);

                   numGlissilePlanesPerBurg =
                           home->burgData.numGlissilePlanesPerBurg;
                   numGlissileGlidePlanes  =
                           numGlissilePlanesPerBurg[bIndex];

/*
 *                 If <planeIndex> corresponds to one of the glissile planes, 
 *                 update the screw density associated with the plane.
 *                 In all other cases, just update the sum total for all the
 *                 sessile planes.
 */
                   if (planeIndex < numGlissileGlidePlanes) {
                       thread_Ltot[bIndex][7+planeIndex] += fabs(segd);
                   } else {
                       thread_Ltot[bIndex][11] += fabs(segd);
                   }

/*
 *                 Increment the total edge density for this burgers vector
 */
                   edgeDensity = sqrt(DotProduct(segleft, segleft));
                   thread_Ltot[bIndex][0] += edgeDensity; 
               }

                for (m = 0; m < HCP_NUM_GLIDE_PLANES; m++) {
                    Ltemp[m] = fabs(tanv[m][0][bIndex]*segleft[0] +
                                    tanv[m][1][bIndex]*segleft[1] +
                                    tanv[m][2][bIndex]*segleft[2]);
                }

/*
 *              Find the indices of the two highest values and select
 *              offsets in the nanv, tanv and areaSwept arrays.  Note that
 *              the 10th burgers vector (bIndex == 9) has 1 fewer planes
 *              than the other
 */
                if (bIndex == 9) {
                    FindMax(Ltemp, HCP_NUM_GLIDE_PLANES-1, &maxVal, &index2);
                    index3 = Get2ndMaxIndex(Ltemp, HCP_NUM_GLIDE_PLANES-1, index2);
                } else {
                    FindMax(Ltemp, HCP_NUM_GLIDE_PLANES, &maxVal, &index2);
                    index3 = Get2ndMaxIndex(Ltemp, HCP_NUM_GLIDE_PLANES, index2);
                }

                if ((index2 < 0) || (index3 < 0)) {
                    continue;
                }

                offset1 = index2;
                offset2 = index3;

                aOffset1 = offset1 + 1;
                aOffset2 = offset2 + 1;

/*
 *              Lower Triangle:
 *
 *              For climb
 */
                xvector(seg[0], seg[1], seg[2], deltax2x, deltax2y, deltax2z,
                        &tmpx, &tmpy, &tmpz);

                thread_areaSwept[bIndex][0] += sb*0.5*(tmpx*bnorm[bIndex][0] +
                                                       tmpy*bnorm[bIndex][1] +
                                                       tmpz*bnorm[bIndex][2]);


                vec1[0] = tanv[offset1][0][bIndex];
                vec1[1] = tanv[offset1][1][bIndex];
                vec1[2] = tanv[offset1][2][bIndex];

                vec2[0] = tanv[offset2][0][bIndex];
                vec2[1] = tanv[offset2][1][bIndex];
                vec2[2] = tanv[offset2][2][bIndex];

                DecompVec(segleft, vec1, vec2, Ltemp2);
                     
                thread_Ltot[bIndex][aOffset1] += fabs(Ltemp2[0]);
                thread_Ltot[bIndex][aOffset2] += fabs(Ltemp2[1]);
                     
                xvector(tanv[offset1][0][bIndex], tanv[offset1][1][bIndex],
                        tanv[offset1][2][bIndex], deltax2x, deltax2y,
                        deltax2z, &tmpx, &tmpy, &tmpz);
                     
                thread_areaSwept[bIndex][aOffset1] += sb * 0.5 * Ltemp2[0] *
                                    (tmpx*nanv[offset1][0][bIndex] +
                                     tmpy*nanv[offset1][1][bIndex] +
                                     tmpz*nanv[offset1][2][bIndex]);
                     
                xvector(tanv[offset2][0][bIndex], tanv[offset2][1][bIndex],
                        tanv[offset2][2][bIndex], deltax2x, deltax2y,
                        deltax2z, &tmpx, &tmpy, &tmpz);
                     
                thread_areaSwept[bIndex][aOffset2] += sb * 0.5 * Ltemp2[1] *
                                    (tmpx*nanv[offset2][0][bIndex] +
                                     tmpy*nanv[offset2][1][bIndex] +
                                     tmpz*nanv[offset2][2][bIndex]);

/*
 *             Screw portion, first...
 */
               xvector(bnorm[bIndex][0], bnorm[bIndex][1], bnorm[bIndex][2],
                       deltax2x, deltax2y, deltax2z, &tmpx, &tmpy, &tmpz);

               qs[0] = sb * 0.5 * segd * tmpx;
               qs[1] = sb * 0.5 * segd * tmpy;
               qs[2] = sb * 0.5 * segd * tmpz;

               stemp[0] = fabs(nanv[0][0][bIndex]*qs[0] +
                               nanv[0][1][bIndex]*qs[1] +
                               nanv[0][2][bIndex]*qs[2]);

               stemp[1] = fabs(nanv[1][0][bIndex]*qs[0] +
                               nanv[1][1][bIndex]*qs[1] +
                               nanv[1][2][bIndex]*qs[2]);
 
               stemp[2] = fabs(nanv[2][0][bIndex]*qs[0] +
                               nanv[2][1][bIndex]*qs[1] +
                               nanv[2][2][bIndex]*qs[2]);

               stemp[3] = fabs(nanv[3][0][bIndex]*qs[0] +
                               nanv[3][1][bIndex]*qs[1] +
                               nanv[3][2][bIndex]*qs[2]);

/*
 *             Find the indices of the two highest values and select
 *             offsets in the nanv, tanv and areaSwept arrays.  Note that
 *             the 10th burgers vector (bIndex == 9) has 1 fewer planes
 *             than the other
 */
               if (bIndex == 9) {
                   FindMax(stemp, HCP_NUM_GLIDE_PLANES-1, &maxVal, &index2);
                   index3 = Get2ndMaxIndex(stemp, HCP_NUM_GLIDE_PLANES-1, index2);
                   offset1 = index2;
                   offset2 = index3;
                   aOffset1 = offset1 + 4;
                   aOffset2 = offset2 + 4;
               } else {
                   FindMax(stemp, HCP_NUM_GLIDE_PLANES, &maxVal, &index2);
                   index3 = Get2ndMaxIndex(stemp, HCP_NUM_GLIDE_PLANES, index2);
                   offset1 = index2;
                   offset2 = index3;
                   aOffset1 = offset1 + 5;
                   aOffset2 = offset2 + 5;
               }

               if ((index2 < 0) || (index3 < 0)) {
                   continue;
               }

               vec1[0] = nanv[offset1][0][bIndex];
               vec1[1] = nanv[offset1][1][bIndex];
               vec1[2] = nanv[offset1][2][bIndex];

               vec2[0] = nanv[offset2][0][bIndex];
               vec2[1] = nanv[offset2][1][bIndex];
               vec2[2] = nanv[offset2][2][bIndex];

               DecompVec(qs, vec1, vec2, outVec);
  
               thread_areaSwept[bIndex][aOffset1] += outVec[0];
               thread_areaSwept[bIndex][aOffset2] += outVec[1];

 
/********************* the other half ****************************/
 
/*
 *              segx vector and delta2 vector defined above
 */
                segd2 = bnorm[bIndex][0]*seg2[0] +
                        bnorm[bIndex][1]*seg2[1] +
                        bnorm[bIndex][2]*seg2[2];

                segleft2[0] = seg2[0] - segd2 * bnorm[bIndex][0];
                segleft2[1] = seg2[1] - segd2 * bnorm[bIndex][1];
                segleft2[2] = seg2[2] - segd2 * bnorm[bIndex][2];
 
                for (m = 0; m < HCP_NUM_GLIDE_PLANES; m++) {
                    Ltemp[m] = fabs(tanv[m][0][bIndex]*segleft2[0] +
                                    tanv[m][1][bIndex]*segleft2[1] +
                                    tanv[m][2][bIndex]*segleft2[2]);
                }

/*
 *              Find the indices of the two highest values and select
 *              offsets in the nanv, tanv and areaSwept arrays.  Note that
 *              the 10th burgers vector (bIndex == 9) has 1 fewer planes
 *              than the other
 */
                if (bIndex == 9) {
                    FindMax(Ltemp, HCP_NUM_GLIDE_PLANES-1, &maxVal, &index2);
                    index3 = Get2ndMaxIndex(Ltemp, HCP_NUM_GLIDE_PLANES-1, index2);
                } else {
                    FindMax(Ltemp, HCP_NUM_GLIDE_PLANES, &maxVal, &index2);
                    index3 = Get2ndMaxIndex(Ltemp, HCP_NUM_GLIDE_PLANES, index2);
                }

                if ((index2 < 0) || (index3 < 0)) {
                    continue;
                }

                offset1 = index2;
                offset2 = index3;

                aOffset1 = offset1 + 1;
                aOffset2 = offset2 + 1;

/*
 *              For climb first
 */
                xvector(seg2[0], seg2[1], seg2[2],
                        deltax1x, deltax1y, deltax1z,
                        &tmpx, &tmpy, &tmpz);

                thread_areaSwept[bIndex][0] += sb * 0.5 * (tmpx*bnorm[bIndex][0] +
                                                          tmpy*bnorm[bIndex][1] +
                                                          tmpz*bnorm[bIndex][2]);
 

                vec1[0] = tanv[offset1][0][bIndex];
                vec1[1] = tanv[offset1][1][bIndex];
                vec1[2] = tanv[offset1][2][bIndex];

                vec2[0] = tanv[offset2][0][bIndex];
                vec2[1] = tanv[offset2][1][bIndex];
                vec2[2] = tanv[offset2][2][bIndex];

                DecompVec(segleft2, vec1, vec2, Ltemp2);

                xvector(tanv[offset1][0][bIndex], tanv[offset1][1][bIndex],
                        tanv[offset1][2][bIndex], deltax1x, deltax1y,
                        deltax1z, &tmpx, &tmpy, &tmpz);
                     
                thread_areaSwept[bIndex][aOffset1] += sb * 0.5 * Ltemp2[0] *
                                    (tmpx*nanv[offset1][0][bIndex] +
                                     tmpy*nanv[offset1][1][bIndex] +
                                     tmpz*nanv[offset1][2][bIndex]);
                     
                xvector(tanv[offset2][0][bIndex], tanv[offset2][1][bIndex],
                        tanv[offset2][2][bIndex], deltax1x, deltax1y,
                        deltax1z, &tmpx, &tmpy, &tmpz);
                     
                thread_areaSwept[bIndex][aOffset2] += sb * 0.5 * Ltemp2[1] *
                                    (tmpx*nanv[offset2][0][bIndex] +
                                     tmpy*nanv[offset2][1][bIndex] +
                                     tmpz*nanv[offset2][2][bIndex]);

/*
 *              Screw portion, second...
 */
                xvector(bnorm[bIndex][0], bnorm[bIndex][1], bnorm[bIndex][2],
                        deltax1x, deltax1y, deltax1z, &tmpx, &tmpy, &tmpz);

                qs[0] = sb * 0.5 * segd2 * tmpx;
                qs[1] = sb * 0.5 * segd2 * tmpy;
                qs[2] = sb * 0.5 * segd2 * tmpz;

                stemp[0] = fabs(nanv[0][0][bIndex]*qs[0] +
                                nanv[0][1][bIndex]*qs[1] +
                                nanv[0][2][bIndex]*qs[2]);

                stemp[1] = fabs(nanv[1][0][bIndex]*qs[0] +
                                nanv[1][1][bIndex]*qs[1] +
                                nanv[1][2][bIndex]*qs[2]);
 
                stemp[2] = fabs(nanv[2][0][bIndex]*qs[0] +
                                nanv[2][1][bIndex]*qs[1] +
                                nanv[2][2][bIndex]*qs[2]);

                stemp[3] = fabs(nanv[3][0][bIndex]*qs[0] +
                                nanv[3][1][bIndex]*qs[1] +
                                nanv[3][2][bIndex]*qs[2]);

/*
 *              Find the indices of the two highest values and select
 *              offsets in the nanv, tanv and areaSwept arrays.  Note that
 *              the 10th burgers vector (bIndex == 9) has 1 fewer planes
 *              than the other
 */
                if (bIndex == 9) {
                    FindMax(stemp, HCP_NUM_GLIDE_PLANES-1, &maxVal, &index2);
                    index3 = Get2ndMaxIndex(stemp, HCP_NUM_GLIDE_PLANES-1, index2);
                    offset1 = index2;
                    offset2 = index3;
                    aOffset1 = offset1 + 4;
                    aOffset2 = offset2 + 4;
                } else {
                    FindMax(stemp, HCP_NUM_GLIDE_PLANES, &maxVal, &index2);
                    index3 = Get2ndMaxIndex(stemp, HCP_NUM_GLIDE_PLANES, index2);
                    offset1 = index2;
                    offset2 = index3;
                    aOffset1 = offset1 + 5;
                    aOffset2 = offset2 + 5;
                }

                if ((index2 < 0) || (index3 < 0)) {
                    continue;
                }

                vec1[0] = nanv[offset1][0][bIndex];
                vec1[1] = nanv[offset1][1][bIndex];
                vec1[2] = nanv[offset1][2][bIndex];

                vec2[0] = nanv[offset2][0][bIndex];
                vec2[1] = nanv[offset2][1][bIndex];
                vec2[2] = nanv[offset2][2][bIndex];

                DecompVec(qs, vec1, vec2, outVec);

                thread_areaSwept[bIndex][aOffset1] += outVec[0];
                thread_areaSwept[bIndex][aOffset2] += outVec[1];

            }  /* for (segIndex = 0; segIndex < numSeg; ...) */
           
#if 0
/*
 *          Some stuff for debugging
 */
            printf("thread_Ltot ---------------- \n");
            for (int ti = 0; ti < HCP_NUM_GLIDE_BURG; ti++) {
                printf("%e %e %e %e %e\n", thread_Ltot[ti][0],
                       thread_Ltot[ti][1], thread_Ltot[ti][2],
                       thread_Ltot[ti][3], thread_Ltot[ti][4]);
            }

            printf("thread_areaSwept printing ---- \n");
            for (int ti = 0; ti < HCP_NUM_GLIDE_BURG; ti++) {
                printf(" %e %e %e  %e %e %e  %e %e %e\n",
                       thread_areaSwept[ti][0], thread_areaSwept[ti][1],
                       thread_areaSwept[ti][2], thread_areaSwept[ti][3],
                       thread_areaSwept[ti][4], thread_areaSwept[ti][5],
                       thread_areaSwept[ti][6], thread_areaSwept[ti][7],
                       thread_areaSwept[ti][8]);
            }
#endif
        }  /* for (nodeIndex = 0; nodeIndex < home->newNodeKeyPtr; ...) */

/*
 *      At this point, each thread has a portion of several arrays.
 *      The threads now update the global arrays with their contributions
 *      in a 'critical' block so there's no conflicts among the threads.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_DPS_HCP_A)
#endif
        {
            for (k = 0; k < 3; k++) {
                for (l = 0; l < 3; l++) {
                    dstn[k][l] += thread_dstn[k][l];
                    dspn[k][l] += thread_dspn[k][l];
                }
            }

            for (k = 0; k < HCP_NUM_GLIDE_BURG; k++) {

                for (l = 0; l < 17; l++) {
                    Ltot[k][l] += thread_Ltot[k][l];
                }

                for (l = 0; l < 9; l++) {
                    areaSwept[k][l] += thread_areaSwept[k][l];
                }
            }

            for (k = 0; k < param->numBurgGroups; k++) {
                param->partialDisloDensity[k] += thread_partialDensity[k];
            }

            param->disloDensity += thread_disloDensity;

        }  /* end "omp critical" section */

        if (thread_partialDensity != NULL) {
            free(thread_partialDensity);
        }

    }  /* end "omp parallel" section */

/*
 *      We have per-plane density for edges and screws, but we want 
 *      need to sum the edge and screw values for per-plane totals
 */
        for (i = 0; i < HCP_NUM_GLIDE_BURG; i++) {
            for (j = 0; j < 5; j++) {
                Ltot[i][12+j] = sqrt((Ltot[i][1+j] * Ltot[i][1+j]) +
                                     (Ltot[i][7+j] * Ltot[i][7+j]));
            }
        }

/*
 *      Density decomposition
 */
        factor0 = param->simVol*bmagsq;

        for (i = 0; i < HCP_NUM_GLIDE_BURG; i++) {

            for (j = 0; j < 9; j++) {
                param->HCP_dfluxtot[i][j] = areaSwept[i][j] / param->simVol / 
                                            param->realdt;
            }

            for (j = 0; j < 17; j++) {
                param->HCP_dLtotDecomp[i][j] = Ltot[i][j] / factor0; 
            }
        }
 
#if 0
/*
 *      Some stuff for debugging
 */
        printf("Ltot ----------------, Factor0=%f\n", factor0);
        for (i = 0; i < HCP_NUM_GLIDE_BURG; i++) {
            printf("%e %e %e %e %e\n", param->HCP_dLtot[i][0],
                   param->HCP_dLtot[i][1], param->HCP_dLtot[i][2],
                   param->HCP_dLtot[i][3], param->HCP_dLtot[i][4]);
        }

        printf("fluxtot printing ---- \n");
        for (i = 0; i < HCP_NUM_GLIDE_BURG; i++) {
            printf(" %e %e %e  %e %e %e  %e %e %e\n",
                   param->HCP_dfluxtot[i][0], param->HCP_dfluxtot[i][1],
                   param->HCP_dfluxtot[i][2], param->HCP_dfluxtot[i][3],
                   param->HCP_dfluxtot[i][4], param->HCP_dfluxtot[i][5],
                   param->HCP_dfluxtot[i][6], param->HCP_dfluxtot[i][7],
                   param->HCP_dfluxtot[i][8]);
        }
#endif

/*
 *      Accumulate delta strain this subcycle into delta strain this cycle
 */
        param->delpStrain[0] = dstn[0][0];
        param->delpStrain[1] = dstn[1][1];
        param->delpStrain[2] = dstn[2][2];
        param->delpStrain[3] = dstn[1][2];
        param->delpStrain[4] = dstn[0][2];
        param->delpStrain[5] = dstn[0][1];
        
        param->delpSpin[0] = dspn[0][0];
        param->delpSpin[1] = dspn[1][1];
        param->delpSpin[2] = dspn[2][2];
        param->delpSpin[3] = dspn[1][2];
        param->delpSpin[4] = dspn[0][2];
        param->delpSpin[5] = dspn[0][1];
                                                                                
#ifdef PARALLEL
/*
 *      We've calculated processor specific values, now sum the
 *      delta strain from all processors, and accumulate into net strain
 */
        MPI_Allreduce(param->delpStrain, gdstn, 6, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(param->delpSpin, gdspn, 6, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);

        for (i = 0; i < 6; i++) {
            param->delpStrain[i] = gdstn[i];
            param->delpSpin[i] = gdspn[i];
        }

/*
 *      Flux decomposition
 */
        MPI_Allreduce((real8 *) param->HCP_dLtotDecomp, (real8 *) param->HCP_LtotDecomp, HCP_NUM_GLIDE_BURG*17, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce((real8 *) param->HCP_dfluxtot   , (real8 *) param->HCP_fluxtot   , HCP_NUM_GLIDE_BURG*9 , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
/*
 *      For serial compilation, no need to accumulate values from remote
 *      processors, just copy the data into the param struct.
 */
        for (i = 0; i < HCP_NUM_GLIDE_BURG; i++) {

            for (j = 0; j < 9; j++) {
                param->HCP_fluxtot[i][j] = param->HCP_dfluxtot[i][j];
            }

            for (j = 0; j < 17; j++) {
                param->HCP_LtotDecomp[i][j] = param->HCP_dLtotDecomp[i][j];
            }
        } 

        for (i = 0; i < 6; i++) {
            gdstn[i] = param->delpStrain[i];
            gdspn[i] = param->delpSpin[i];
        }
#endif /* ifdef PARALLEL */

#if 0
/*
 *      For debugging only
 */
        if (home->myDomain == 0) {
            int   k, l;
            real8 pstn[3][3], dyad[3][3];
            real8 dyad0[3][3], dyad1[3][3], dyad2[3][3], dyad3[3][3];

            Matrix33_Zero(pstn);
/*
 *          Strain increment using decomposition here
 */
            for (i = 0; i < HCP_NUM_GLIDE_BURG; i++) {

                /* climb */
                dyad[0][0] = bnorm[i][0]*bnorm[i][0];
                dyad[0][1] = bnorm[i][0]*bnorm[i][1];
                dyad[0][2] = bnorm[i][0]*bnorm[i][2];

                dyad[1][0] = bnorm[i][1]*bnorm[i][0];
                dyad[1][1] = bnorm[i][1]*bnorm[i][1];
                dyad[1][2] = bnorm[i][1]*bnorm[i][2];

                dyad[2][0] = bnorm[i][2]*bnorm[i][0];
                dyad[2][1] = bnorm[i][2]*bnorm[i][1];
                dyad[2][2] = bnorm[i][2]*bnorm[i][2];

                for (k = 0; k < 3; k++){
                    for (l = 0; l < 3; l++){
                        pstn[k][l] += param->HCP_fluxtot[i][0] *
                                      (dyad[l][k]+dyad[k][l]) *
                                      burgMagLocal[i] * param->realdt * 0.5;
                    }
                }
            }
              
            for (i = 0; i < HCP_NUM_GLIDE_BURG-1; i++) {
                /* Four planes */
                dyad0[0][0] = bnorm[i][0]*nanv[0][0][i];
                dyad0[0][1] = bnorm[i][0]*nanv[0][1][i];
                dyad0[0][2] = bnorm[i][0]*nanv[0][2][i]; 
                
                dyad0[1][0] = bnorm[i][1]*nanv[0][0][i];
                dyad0[1][1] = bnorm[i][1]*nanv[0][1][i];
                dyad0[1][2] = bnorm[i][1]*nanv[0][2][i]; 

                dyad0[2][0] = bnorm[i][2]*nanv[0][0][i];
                dyad0[2][1] = bnorm[i][2]*nanv[0][1][i];
                dyad0[2][2] = bnorm[i][2]*nanv[0][2][i]; 


                dyad1[0][0] = bnorm[i][0]*nanv[1][0][i];
                dyad1[0][1] = bnorm[i][0]*nanv[1][1][i];
                dyad1[0][2] = bnorm[i][0]*nanv[1][2][i]; 

                dyad1[1][0] = bnorm[i][1]*nanv[1][0][i];
                dyad1[1][1] = bnorm[i][1]*nanv[1][1][i];
                dyad1[1][2] = bnorm[i][1]*nanv[1][2][i]; 

                dyad1[2][0] = bnorm[i][2]*nanv[1][0][i];
                dyad1[2][1] = bnorm[i][2]*nanv[1][1][i];
                dyad1[2][2] = bnorm[i][2]*nanv[1][2][i]; 


                dyad2[0][0] = bnorm[i][0]*nanv[2][0][i];
                dyad2[0][1] = bnorm[i][0]*nanv[2][1][i];
                dyad2[0][2] = bnorm[i][0]*nanv[2][2][i]; 

                dyad2[1][0] = bnorm[i][1]*nanv[2][0][i];
                dyad2[1][1] = bnorm[i][1]*nanv[2][1][i];
                dyad2[1][2] = bnorm[i][1]*nanv[2][2][i]; 

                dyad2[2][0] = bnorm[i][2]*nanv[2][0][i];
                dyad2[2][1] = bnorm[i][2]*nanv[2][1][i];
                dyad2[2][2] = bnorm[i][2]*nanv[2][2][i]; 


                dyad3[0][0] = bnorm[i][0]*nanv[3][0][i];
                dyad3[0][1] = bnorm[i][0]*nanv[3][1][i];
                dyad3[0][2] = bnorm[i][0]*nanv[3][2][i]; 

                dyad3[1][0] = bnorm[i][1]*nanv[3][0][i];
                dyad3[1][1] = bnorm[i][1]*nanv[3][1][i];
                dyad3[1][2] = bnorm[i][1]*nanv[3][2][i]; 

                dyad3[2][0] = bnorm[i][2]*nanv[3][0][i];
                dyad3[2][1] = bnorm[i][2]*nanv[3][1][i];
                dyad3[2][2] = bnorm[i][2]*nanv[3][2][i]; 

                for (k = 0; k < 3; k++) {
                    for(l = 0; l < 3; l++) {
                    
                        pstn[k][l] += ((param->HCP_fluxtot[i][1]+param->HCP_fluxtot[i][5]) * (dyad0[l][k]+dyad0[k][l]) +
                                       (param->HCP_fluxtot[i][2]+param->HCP_fluxtot[i][6]) * (dyad1[l][k]+dyad1[k][l]) +
                                       (param->HCP_fluxtot[i][3]+param->HCP_fluxtot[i][7]) * (dyad2[l][k]+dyad2[k][l]) +
                                       (param->HCP_fluxtot[i][4]+param->HCP_fluxtot[i][8]) * (dyad3[l][k]+dyad3[k][l]) ) * burgMagLocal[i] * param->realdt * 0.5;
                    }
                }
            }  /* for (i=0;i<HCP_NUM_GLIDE_BURG-1; etc. */


            /* Treat here the special case of Burgers vector 9 */
            
            i = 9;
            
            /* Three planes */
            dyad0[0][0] = bnorm[i][0]*nanv[0][0][i];
            dyad0[0][1] = bnorm[i][0]*nanv[0][1][i];
            dyad0[0][2] = bnorm[i][0]*nanv[0][2][i]; 
            
            dyad0[1][0] = bnorm[i][1]*nanv[0][0][i];
            dyad0[1][1] = bnorm[i][1]*nanv[0][1][i];
            dyad0[1][2] = bnorm[i][1]*nanv[0][2][i]; 
            
            dyad0[2][0] = bnorm[i][2]*nanv[0][0][i];
            dyad0[2][1] = bnorm[i][2]*nanv[0][1][i];
            dyad0[2][2] = bnorm[i][2]*nanv[0][2][i]; 

            
            dyad1[0][0] = bnorm[i][0]*nanv[1][0][i];
            dyad1[0][1] = bnorm[i][0]*nanv[1][1][i];
            dyad1[0][2] = bnorm[i][0]*nanv[1][2][i]; 
            
            dyad1[1][0] = bnorm[i][1]*nanv[1][0][i];
            dyad1[1][1] = bnorm[i][1]*nanv[1][1][i];
            dyad1[1][2] = bnorm[i][1]*nanv[1][2][i]; 
            
            dyad1[2][0] = bnorm[i][2]*nanv[1][0][i];
            dyad1[2][1] = bnorm[i][2]*nanv[1][1][i];
            dyad1[2][2] = bnorm[i][2]*nanv[1][2][i]; 
            

            dyad2[0][0] = bnorm[i][0]*nanv[2][0][i];
            dyad2[0][1] = bnorm[i][0]*nanv[2][1][i];
            dyad2[0][2] = bnorm[i][0]*nanv[2][2][i]; 
            
            dyad2[1][0] = bnorm[i][1]*nanv[2][0][i];
            dyad2[1][1] = bnorm[i][1]*nanv[2][1][i];
            dyad2[1][2] = bnorm[i][1]*nanv[2][2][i]; 
            
            dyad2[2][0] = bnorm[i][2]*nanv[2][0][i];
            dyad2[2][1] = bnorm[i][2]*nanv[2][1][i];
            dyad2[2][2] = bnorm[i][2]*nanv[2][2][i]; 
            
            for (k = 0; k < 3; k++) {
                for(l = 0; l < 3; l++) {
                
                    pstn[k][l] += ((param->HCP_fluxtot[i][1]+param->HCP_fluxtot[i][4]) * (dyad0[l][k]+dyad0[k][l]) +
                                   (param->HCP_fluxtot[i][2]+param->HCP_fluxtot[i][5]) * (dyad1[l][k]+dyad1[k][l]) +
                                   (param->HCP_fluxtot[i][3]+param->HCP_fluxtot[i][6]) * (dyad2[l][k]+dyad2[k][l])) * burgMagLocal[i] * param->realdt * 0.5;
                }
            }

/* 
 *          gdstn tensor the total sum of all decomposed strain flux components
 *          The sum should be identical to delpStrain if the decomposition is
 *          done correctly .. To check if both methods produce the same
 *          result ..  Do not delete the print statements.
 */
#if 0
            {
                real8 ratio[6];

                printf("--simVoll = %e -----\n", param->simVol);
                printf("--realdt = %e -----\n", param->realdt);

                printf("gdstn= %.8e %.8e %.8e %.8e %.8e %.8e \n",
                       gdstn[0], gdstn[1], gdstn[2], gdstn[3],
                       gdstn[4], gdstn[5]);

                printf("pstn = %.8e %.8e %.8e %.8e %.8e %.8e \n",
                       pstn[0][0], pstn[1][1], pstn[2][2], pstn[1][2],
                       pstn[0][2], pstn[0][1]);
           
                ratio[0] = (pstn[0][0] == 0.0) ? 0.0 : gdstn[0]/pstn[0][0];
                ratio[1] = (pstn[1][1] == 0.0) ? 0.0 : gdstn[1]/pstn[1][1];
                ratio[2] = (pstn[2][2] == 0.0) ? 0.0 : gdstn[2]/pstn[2][2];
                ratio[3] = (pstn[1][2] == 0.0) ? 0.0 : gdstn[3]/pstn[1][2];
                ratio[4] = (pstn[0][2] == 0.0) ? 0.0 : gdstn[4]/pstn[0][2];
                ratio[5] = (pstn[0][1] == 0.0) ? 0.0 : gdstn[5]/pstn[0][1];

                printf("gdstn/pstn = %.8e %.8e %.8e %.8e %.8e %.8e \n",
                       ratio[0], ratio[1], ratio[2],
                       ratio[3], ratio[4], ratio[5]);
            }
#endif
        } /*if masternode */
#endif

#if 0
        if (home->myDomain == 0) {
            printf("param->delpStrain[0]= %e\n",param->delpStrain[0]);
            printf("param->delpStrain[1]= %e\n",param->delpStrain[1]);
            printf("param->delpStrain[2]= %e\n",param->delpStrain[2]);
            printf("param->delpStrain[3]= %e\n",param->delpStrain[3]);
            printf("param->delpStrain[4]= %e\n",param->delpStrain[4]);
            printf("param->delpStrain[5]= %e\n",param->delpStrain[5]);
        }
#endif                                                                  
#if 0
        if (home->myDomain == 0) {
            printf("param->delpSpin[0]= %e\n",param->delpSpin[0]);
            printf("param->delpSpin[1]= %e\n",param->delpSpin[1]);
            printf("param->delpSpin[2]= %e\n",param->delpSpin[2]);
            printf("param->delpSpin[3]= %e\n",param->delpSpin[3]);
            printf("param->delpSpin[4]= %e\n",param->delpSpin[4]);
            printf("param->delpSpin[5]= %e\n",param->delpSpin[5]);
        }

#endif                                                                  

        return;
}
