/***************************************************************************
 *   
 *      Author:       Moono Rhee
 *
 *      Function:     DeltaPlasticStrain_BCC
 *
 *      Description:  Calculate the plastic strain increment
 *                    by dislocation motion in BCC materials.
 *
 *      Last Modified:  
 *        01/08/2002 - original version
 *        05/20/2002 - M. Rhee: revisited to include dislocation density
 *        03/19/2004 - M. Hiratani: Double counting of nodes
 *                     is removed.  Accordingly, dislocationdensity was
 *                     recorrected.
 *        03/25/2011 - MRhee, re-wrote to treat the lab-frame dislocation
 *                     geometry.
 *                   - reason:2x2 matrices can't be inverted for some frames
 *                     when L,flux decomp is done
 *                   - burg mags in 'pstn' are scaled for diff b mags in
 *                     input file
 *
 *      Public functions:
 *          DeltaPlasticStrain_BCC()
 *
 *      Private functions:
 *          InitDeltaPlasticStrain_BCC()
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
 *      Function:    InitDeltaPlasticStrain_BCC
 *
 *      Description: Several of the arrays used during the delta plastic
 *                   strain calculations are dynamically calculated, but
 *                   once set, never need to be changed again for the
 *                   duration of the simulation.  This function should be
 *                   invoked during the first call to DeltaPlasticStrain_BCC()
 *                   to initialize those arrays, and never used again
 *                   after that.
 *
 ***************************************************************************/
static void InitDeltaPlasticStrain_BCC(Home_t *home,
                                       real8 nanv[3][3][4],
                                       real8 tanv[3][3][4],
                                       real8 burgv[4][3],
                                       real8 bnorm[4][3])
{
        int      i, j, k, bIndex;
        real8    sqt2, sqt6;
        Param_t *param = home->param;

        nanv[0][0][0]=  0.; nanv[0][1][0]=  1.; nanv[0][2][0]= -1.;
        nanv[1][0][0]= -1.; nanv[1][1][0]=  0.; nanv[1][2][0]=  1.;
        nanv[2][0][0]=  1.; nanv[2][1][0]= -1.; nanv[2][2][0]=  0.;
 
        nanv[0][0][1]=  0.; nanv[0][1][1]= -1.; nanv[0][2][1]=  1.;
        nanv[1][0][1]=  1.; nanv[1][1][1]=  0.; nanv[1][2][1]=  1.;
        nanv[2][0][1]= -1.; nanv[2][1][1]= -1.; nanv[2][2][1]=  0.;
 
        nanv[0][0][2]=  0.; nanv[0][1][2]= -1.; nanv[0][2][2]= -1.;
        nanv[1][0][2]=  1.; nanv[1][1][2]=  0.; nanv[1][2][2]= -1.;
        nanv[2][0][2]=  1.; nanv[2][1][2]=  1.; nanv[2][2][2]=  0.;
 
        nanv[0][0][3]=  0.; nanv[0][1][3]=  1.; nanv[0][2][3]=  1.;
        nanv[1][0][3]= -1.; nanv[1][1][3]=  0.; nanv[1][2][3]= -1.;
        nanv[2][0][3]= -1.; nanv[2][1][3]=  1.; nanv[2][2][3]=  0.;

 
        tanv[0][0][0] = -2.0; tanv[0][1][0] =  1.0; tanv[0][2][0] =  1.0;
        tanv[1][0][0] =  1.0; tanv[1][1][0] = -2.0; tanv[1][2][0] =  1.0;
        tanv[2][0][0] =  1.0; tanv[2][1][0] =  1.0; tanv[2][2][0] = -2.0;
 
        tanv[0][0][1] =  2.0; tanv[0][1][1] =  1.0; tanv[0][2][1] =  1.0;
        tanv[1][0][1] =  1.0; tanv[1][1][1] =  2.0; tanv[1][2][1] = -1.0;
        tanv[2][0][1] =  1.0; tanv[2][1][1] = -1.0; tanv[2][2][1] =  2.0;
 
        tanv[0][0][2] =  2.0; tanv[0][1][2] =  1.0; tanv[0][2][2] = -1.0;
        tanv[1][0][2] =  1.0; tanv[1][1][2] =  2.0; tanv[1][2][2] =  1.0;
        tanv[2][0][2] = -1.0; tanv[2][1][2] =  1.0; tanv[2][2][2] =  2.0;
 
        tanv[0][0][3] =  2.0; tanv[0][1][3] = -1.0; tanv[0][2][3] =  1.0;
        tanv[1][0][3] = -1.0; tanv[1][1][3] =  2.0; tanv[1][2][3] =  1.0;
        tanv[2][0][3] =  1.0; tanv[2][1][3] =  1.0; tanv[2][2][3] =  2.0;
 

        sqt2 = M_SQRT1_2;
        sqt6 = 1.0 / sqrt(6.0);
 
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                 for (k = 0; k < 4; k++) {
                     nanv[i][j][k] = nanv[i][j][k] * sqt2;
                     tanv[i][j][k] = tanv[i][j][k] * sqt6;
                 }
            }
        }

#if 0
/*
 *      Just some stuff for debugging
 */
        for (j = 0; j < 4; j++) {
            for (i = 0; i < 3; i++) {
                printf(" nanv = %e %e %e \n",
                       nanv[i][0][j], nanv[i][1][j], nanv[i][2][j]);
            }
            printf("\n");
        }

        printf("\n");

        for (j = 0; j < 4; j++) {
            for (i = 0; i < 3; i++) {
                printf(" tanv = %e %e %e \n",
                       tanv[i][0][j], tanv[i][1][j], tanv[i][2][j]);
            }
            printf("\n");
        }
#endif

/*
 *      If necessary, rotate the nanv and tanv arrays from the crystal frame
 *      to the lab frame
 */
        if (param->useLabFrame) {
            int plane;

            for (plane = 0; plane < 3; plane++) {

                for (bIndex = 0; bIndex < 4; bIndex++) {
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
        }

#if 0
/*
 *      Debug only...
 */
        printf("---------------------------------------------------------\n");
        printf("After rotation \n");
        printf("---------------------------------------------------------\n");

        for (j = 0; j < 4; j++){
            for (i = 0; i < 3; i++) {
                printf(" nanv = %e %e %e \n",
                       nanv[i][0][j], nanv[i][1][j], nanv[i][2][j]);
            }
            printf("\n");
        }

        printf("\n");

        for (j = 0; j < 4; j++){
            for (i = 0; i < 3; i++) {
                printf(" tanv = %e %e %e \n",
                       tanv[i][0][j], tanv[i][1][j], tanv[i][2][j]);
            }
            printf("\n");
        }
        printf("---------------------------------------------------------\n");
#endif

        burgv[0][0] =  0.5; burgv[0][1] =  0.5; burgv[0][2] =  0.5;
        burgv[1][0] = -0.5; burgv[1][1] =  0.5; burgv[1][2] =  0.5;
        burgv[2][0] =  0.5; burgv[2][1] = -0.5; burgv[2][2] =  0.5;
        burgv[3][0] =  0.5; burgv[3][1] =  0.5; burgv[3][2] = -0.5;
    
        for (i = 0; i < 4; i++) { /* 111 type only; otherwise, would be 7 */
            for (j = 0; j < 3; j++) {
                 bnorm[i][j]=burgv[i][j] / sqrt(0.75);
            }
        }

/*
 *      If necessary, rotate the bnorm array from the crystal frame to
 *      the lab frame
 */
        if (param->useLabFrame) {

            for (bIndex = 0; bIndex < 4; bIndex++) {
                real8 bnormLab[3];

                Matrix33Vector3Multiply(param->rotMatrix,bnorm[bIndex],bnormLab);
                VECTOR_COPY(bnorm[bIndex], bnormLab);
            }
        }
 
#if 0
/*
 *      Debug only. Only look at the first 4 (the 111 types)
 */
        for (i = 0; i < 4; i++) {
            printf(" bnorm = %e %e %e \n",
                   bnorm[i][0], bnorm[i][1], bnorm[i][2]);
        }
#endif

        return;
}


void DeltaPlasticStrain_BCC(Home_t *home)
{
        int          i, j;
        real8        factor0, bmagsq;
        real8        eps=1.e-12;
        real8        gdspn[6], gdstn[6];
        real8        areaSwept[4][7];
        real8        Ltot[4][11];
        real8        dstn[3][3], dspn[3][3];
        Param_t      *param;
        static int   initDone = 0;
        static real8 nanv[3][3][4], tanv[3][3][4];
        static real8 burgv[4][3], bnorm[4][3];

// Ltot contents:
//   [burg][0]  total edge_density (based on total length - screw length)
//   [burg][1]  edge_density plane 1
//   [burg][2]  edge_density plane 2
//   [burg][3]  edge_density plane 3
//   [burg][4]  total screw_density (based on screw length)
//   [burg][5]  screw_density plane 1
//   [burg][6]  screw_density plane 2
//   [burg][7]  screw_density plane 3
//   [burg][8]  total_density plane 1 = sqrt(edge_density1^2 + screw_density1^2)
//   [burg][9]  total_density plane 2 = sqrt(edge_density2^2 + screw_density2^2)
//   [burg][10] total_density plane 3 = sqrt(edge_density3^2 + screw_density3^2)

        param = home->param;

        bmagsq = param->burgMag * param->burgMag;
        param->disloDensity = 0.0;

        memset(Ltot, 0, sizeof(Ltot));
        memset(areaSwept, 0, sizeof(areaSwept));
 
        Matrix33_Zero(dstn);
        Matrix33_Zero(dspn);

        for (i = 0; i < param->numBurgGroups; i++) {
            param->partialDisloDensity[i] = 0.0;
        }

        for (i=0; i < 6; i++) {
            param->delpStrain[i] = 0.0;                                         
        }                                                                

/*
 *      The first time we enter this function, initialize a few
 *      arrays that will remain static for the duration of the
 *      simulation once they are set.
 */
        if (!initDone++) {
            InitDeltaPlasticStrain_BCC(home, nanv, tanv, burgv, bnorm);
        }
             
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
        int    tk, tl;
        int    nodeIndex;
        int    threadID, threadIterStart, threadIterEnd;
        real8  thread_disloDensity;
        real8  *thread_partialDensity;
        real8  thread_dstn[3][3], thread_dspn[3][3];
        real8  thread_dyad[3][3];
        real8  thread_areaSwept[4][7];
        real8  thread_Ltot[4][11];

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
            int    burgGroup, zeroCount, index, index2;
            int    planeIndex;
            real8  bx, by, bz;
            real8  nx, ny, nz;
            real8  tx, ty, tz;
            real8  tmpx, tmpy, tmpz;
            real8  deltax1x, deltax1y, deltax1z;
            real8  deltax2x, deltax2y, deltax2z;
            real8  size, tmpmax;
            real8  deltaDensity;
            real8  sb, sbsign;
            real8  segd, segd2;
            real8  qs[3], stemp[3];
            real8  vec1[3], vec2[3], outVec[2];
            real8  seg[3], seg2[3];
            real8  Ltemp[3], Ltemp2[2];
            real8  segleft[3], segleft2[3];
            real8  localstrain[3][3];
            Node_t *node, *nbr;

            if ((node = home->nodeKeys[nodeIndex]) == (Node_t *)NULL) {
                continue;
            }

            numSeg = node->numNbrs;

            Matrix33_Zero(localstrain);
 
/*
 *          Loop over every segment attached to the node
 */        
            for (segIndex = 0; segIndex < numSeg; segIndex++) {
               int   offset1, offset2;
               real8 minVal;
               real8 edgeDensity;
               real8 ex, ey, ez;
               real8 hhx, hhy, hhz;
               real8 delxnode, delynode, delznode;
               real8 delxnbr, delynbr, delznbr;
               real8 bCryst[3], zeta[3];

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

/* 
 *             0.5 for cross(e, hh) => area
 */
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
           
               for (tk = 0; tk < 3; tk++){
                  for (tl = 0; tl < 3; tl++){
                      real8 tmpDstn;
                      tmpDstn = (thread_dyad[tl][tk] + thread_dyad[tk][tl]) *
                                0.5 / param->simVol;
                      thread_dstn[tl][tk] += tmpDstn;
                      thread_dspn[tl][tk] += (thread_dyad[tl][tk] -
                                            thread_dyad[tk][tl]) *
                                           0.5 / param->simVol;
                      localstrain[tl][tk] += tmpDstn;
                  }
               }

               /* Modified for flux decomposition calculation 06/22/04 M.Rhee */
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

               size=sqrt(zeta[0]*zeta[0] + zeta[1]*zeta[1] + zeta[2]*zeta[2]);

/*
 *             deltaDensity = size/param->simVol/bmagsq;
 */
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
 *             In addition to tracking total dislocation density, we need
 *             track dislocation density for specific sets of burgers
 *             vectors.  The mobility type (BCC, FCC, etc) determine
 *             the number of such groupings and which burgers vectors
 *             are of interest.
 *
 *             For BCC mobility types we track density of the following
 *             groupings of burgers vectors.  Any burgers vectors not
 *             mentioned are ignored.
 *
 *                     group #     burgers vector types
 *                       0         [ 1 1 1] [-1-1-1]
 *                       1         [-1 1 1] [ 1-1-1]
 *                       2         [ 1-1 1] [-1 1-1]
 *                       3         [ 1 1-1] [-1-1 1]
 *                       4         [ 1 0 0] [-1 0 0]
 *                                 [ 0 1 0] [ 0-1 0]
 *                                 [ 0 0 1] [ 0 0-1]
 *
 *             *************************************************
 *             ***                                           ***
 *             ***                  IMPORTANT!               ***
 *             ***   If you change any groupings of burgers  ***
 *             ***   vectors below, you must change the      ***
 *             ***   DENSITY_FILE_VERSION number defined     ***
 *             ***   in WriteProp.c!                         ***
 *             ***                                           ***
 *             *************************************************
 */
               burgGroup = -1;

               zeroCount = (bCryst[X] == 0.0) +
                           (bCryst[Y] == 0.0) +
                           (bCryst[Z] == 0.0);

               switch(zeroCount) {

                  case 0:
                    /* all types  [-]1 [-]1 [-]1 */

                    tx = bCryst[X];  ty = bCryst[Y];  tz = bCryst[Z];

                    if ((tx * ty * tz) < 0.0) {
                         tx = -bCryst[X];
                         ty = -bCryst[Y];
                         tz = -bCryst[Z];
                    }

                    /* types -1  1  1 */
                    if ((ty < 0.0) && (tz < 0.0)) burgGroup = 1;
 
                    /* types  1 -1  1 */
                    if ((tx < 0.0) && (tz < 0.0)) burgGroup = 2;

                    /* types  1  1 -1 */
                    if ((tx < 0.0) && (ty < 0.0)) burgGroup = 3;

                    /* types  1  1  1 */
                    if ((tx > 0.0) && (ty > 0.0) && (tz > 0.0))
                        burgGroup = 0;

                    break;

                 case 2:
/*
 *                  one of 100 types; group all these into a single category
 */
                    burgGroup = 4;
                    break;

                  default:
/*
 *                  Ignore all other types of burgers vectors
 */
                    break;
               } /* end switch (zeroCount) */

/*
 *             And, if the burger's vector fell into one of groups of
 *             burgers vectors whose density we're tracking, increment
 *             that group's density value.
 */
               if (burgGroup >= 0) {
                   thread_partialDensity[burgGroup] += deltaDensity;
               }

/*
 *             If this is not a [1 1 1] type burgers vector, skip it.
 */
               if ((fabs(bCryst[X])*fabs(bCryst[Y])*fabs(bCryst[Z])) < eps) {
                   continue;
               }

/*
 *             max index for 4 burgers vector for now
 */
               index = -1;
               tmpmax=0.0;

               for (m = 0; m < 4; m++) {
                   real8 tmpVal;
                   tmpVal = fabs(bnorm[m][0]*bx +
                                 bnorm[m][1]*by +
                                 bnorm[m][2]*bz);
                   if (tmpVal > tmpmax) {
                       tmpmax = tmpVal;
                       index = m;
                   }
               }

               if (index < 0) {
                   continue;
               }

               sbsign = bnorm[index][0]* bx +
                        bnorm[index][1]* by +
                        bnorm[index][2]* bz;

               sb = ((sbsign < 0) ? -1.0 : 1.0);

/*
 *             Find the index of the appropriate glide plane for
 *             this segment.
 */
	       planeIndex = -1;
	       tmpmax = 0.0;
	       
	       for (m = 0; m < 3; m++) {
                    real8 tmpVal;
                    tmpVal = fabs(nanv[m][0][index]*nx +
                                  nanv[m][1][index]*ny +
                                  nanv[m][2][index]*nz);
                    if (tmpVal > tmpmax) {
                        tmpmax = tmpVal;
                        planeIndex = m;
                    }
	       }
	       
	       if (planeIndex < 0) {
		 continue;
	       }

/*
 *             segx vector and delta2 vector defined above
 *
 *             segd is the dotted value of seg length with
 *             burgers vector => Screw density
 *
 *             segleft is the total edge density for the segment
 */
               segd = bnorm[index][0]*seg[0] +
                      bnorm[index][1]*seg[1] +
                      bnorm[index][2]*seg[2];

               segleft[0] = seg[0] - segd * bnorm[index][0];
               segleft[1] = seg[1] - segd * bnorm[index][1];
               segleft[2] = seg[2] - segd * bnorm[index][2];

               edgeDensity = sqrt(DotProduct(segleft, segleft));

/*
 *             Update the total edge density for this burgers vector,
 *             the total screw density for this burgers vector and the
 *             screw density for this plane.
 */
               thread_Ltot[index][0] += edgeDensity; // total edge density
               thread_Ltot[index][4] += fabs(segd);  // total screw density
               thread_Ltot[index][5+planeIndex] += fabs(segd);

/*
 *             min index2
 */
               for (m = 0; m < 3; m++) {
                   Ltemp[m] = fabs(tanv[m][0][index]*segleft[0] +
                                   tanv[m][1][index]*segleft[1] +
                                   tanv[m][2][index]*segleft[2]);
               }

               FindMin(Ltemp, 3, &minVal, &index2);
 
/*
 *             Lower Triangle:
 *
 *             For climb  
 */
               xvector(seg[0], seg[1], seg[2], deltax2x, deltax2y, deltax2z,
                       &tmpx, &tmpy, &tmpz);

               thread_areaSwept[index][0] += sb*0.5*(tmpx*bnorm[index][0] +
                                                     tmpy*bnorm[index][1] +
                                                     tmpz*bnorm[index][2]);
 
               switch (index2) {
                   case 0:
                       offset1 = 1;
                       offset2 = 2;
                       break;
                   case 1:
                       offset1 = 0;
                       offset2 = 2;
                       break;
                   default:
                       offset1 = 0;
                       offset2 = 1;
                       break;
               }

               vec1[0] = tanv[offset1][0][index];
               vec1[1] = tanv[offset1][1][index];
               vec1[2] = tanv[offset1][2][index];

               vec2[0] = tanv[offset2][0][index];
               vec2[1] = tanv[offset2][1][index];
               vec2[2] = tanv[offset2][2][index];

               DecompVec(segleft, vec1, vec2, Ltemp2);

               thread_Ltot[index][offset1+1] += fabs(Ltemp2[0]);
               thread_Ltot[index][offset2+1] += fabs(Ltemp2[1]);

               xvector(tanv[offset1][0][index], tanv[offset1][1][index],
                       tanv[offset1][2][index], deltax2x, deltax2y, deltax2z, 
                       &tmpx, &tmpy, &tmpz);

               thread_areaSwept[index][offset1+1] += sb * 0.5 * Ltemp2[0] *
                                      (tmpx*nanv[offset1][0][index] +
                                       tmpy*nanv[offset1][1][index] +
                                       tmpz*nanv[offset1][2][index]);
       
               xvector(tanv[offset2][0][index], tanv[offset2][1][index],
                       tanv[offset2][2][index], deltax2x, deltax2y,
                       deltax2z, &tmpx, &tmpy, &tmpz);

               thread_areaSwept[index][offset2+1] += sb * 0.5 * Ltemp2[1] *
                                      (tmpx*nanv[offset2][0][index] +
                                       tmpy*nanv[offset2][1][index] +
                                       tmpz*nanv[offset2][2][index]);

/*
 *             For screw (first part, lower triangle), decompose 'qs' vector 
 *             for two highest planes.
 */
               xvector(bnorm[index][0], bnorm[index][1], bnorm[index][2],
                       deltax2x, deltax2y, deltax2z, &tmpx, &tmpy, &tmpz);

               qs[0] = sb * 0.5 * segd * tmpx;
               qs[1] = sb * 0.5 * segd * tmpy;
               qs[2] = sb * 0.5 * segd * tmpz;

               stemp[0] = fabs(nanv[0][0][index]*qs[0] +
                               nanv[0][1][index]*qs[1] +
                               nanv[0][2][index]*qs[2]);

               stemp[1] = fabs(nanv[1][0][index]*qs[0] +
                               nanv[1][1][index]*qs[1] +
                               nanv[1][2][index]*qs[2]);

               stemp[2] = fabs(nanv[2][0][index]*qs[0] +
                               nanv[2][1][index]*qs[1] +
                               nanv[2][2][index]*qs[2]);
 
/*
 *             Find min index2
 */
               FindMin(stemp, 3, &minVal, &index2);

               switch (index2) {
                   case 0:
                       offset1 = 1;
                       offset2 = 2;
                       break;
                   case 1:
                       offset1 = 0;
                       offset2 = 2;
                       break;
                   default:
                       offset1 = 0;
                       offset2 = 1;
                       break;
               }

               vec1[0] = nanv[offset1][0][index];
               vec1[1] = nanv[offset1][1][index];
               vec1[2] = nanv[offset1][2][index];

               vec2[0] = nanv[offset2][0][index];
               vec2[1] = nanv[offset2][1][index];
               vec2[2] = nanv[offset2][2][index];

               DecompVec(qs, vec1, vec2, outVec); 

               thread_areaSwept[index][offset1+4] += outVec[0];
               thread_areaSwept[index][offset2+4] += outVec[1];

/*
 *             The other half, upper triangle:
 *
 *             segx vector and delta2 vector defined above
 */
               segd2 = bnorm[index][0]*seg2[0] +
                       bnorm[index][1]*seg2[1] +
                       bnorm[index][2]*seg2[2];

               segleft2[0] = seg2[0] - segd2 * bnorm[index][0];
               segleft2[1] = seg2[1] - segd2 * bnorm[index][1];
               segleft2[2] = seg2[2] - segd2 * bnorm[index][2];
 
/*
 *             min index2
 */
               for (m = 0; m < 3; m++) {
                   Ltemp[m] = fabs(tanv[m][0][index]*segleft2[0] +
                                   tanv[m][1][index]*segleft2[1] +
                                   tanv[m][2][index]*segleft2[2]);
               }

               FindMin(Ltemp, 3, &minVal, &index2);
 
/*
 *             for climb first
 */
               xvector(seg2[0], seg2[1], seg2[2],
                       deltax1x, deltax1y, deltax1z,
                       &tmpx, &tmpy, &tmpz);

               thread_areaSwept[index][0] += sb * 0.5 * (tmpx*bnorm[index][0] +
                                                         tmpy*bnorm[index][1] +
                                                         tmpz*bnorm[index][2]);
 
               switch (index2) {
                   case 0:
                       offset1 = 1;
                       offset2 = 2;
                       break;
                   case 1:
                       offset1 = 0;
                       offset2 = 2;
                       break;
                   default:
                       offset1 = 0;
                       offset2 = 1;
                       break;
               }

               vec1[0] = tanv[offset1][0][index];
               vec1[1] = tanv[offset1][1][index];
               vec1[2] = tanv[offset1][2][index];

               vec2[0] = tanv[offset2][0][index];
               vec2[1] = tanv[offset2][1][index];
               vec2[2] = tanv[offset2][2][index];

               DecompVec(segleft2, vec1, vec2, Ltemp2);

               xvector(tanv[offset1][0][index], tanv[offset1][1][index],
                       tanv[offset1][2][index], deltax1x, deltax1y,
                       deltax1z, &tmpx, &tmpy, &tmpz);

               thread_areaSwept[index][offset1+1] += sb * 0.5 * Ltemp2[0] *
                                      (tmpx*nanv[offset1][0][index] +
                                       tmpy*nanv[offset1][1][index] +
                                       tmpz*nanv[offset1][2][index]);
 
               xvector(tanv[offset2][0][index], tanv[offset2][1][index],
                       tanv[offset2][2][index], deltax1x, deltax1y,
                       deltax1z, &tmpx, &tmpy, &tmpz);

               thread_areaSwept[index][offset2+1] += sb * 0.5 * Ltemp2[1] *
                                      (tmpx*nanv[offset2][0][index] +
                                       tmpy*nanv[offset2][1][index] +
                                       tmpz*nanv[offset2][2][index]);
 
/*
 *             For screw (second part, upper triangle), decompose 'qs' vector 
 *             for two highest planes.
 */
               xvector(bnorm[index][0], bnorm[index][1], bnorm[index][2],
                       deltax1x, deltax1y, deltax1z, &tmpx, &tmpy, &tmpz);

               qs[0] = sb * 0.5 * segd2 * tmpx;
               qs[1] = sb * 0.5 * segd2 * tmpy;
               qs[2] = sb * 0.5 * segd2 * tmpz;

               stemp[0] = fabs(nanv[0][0][index]*qs[0] +
                               nanv[0][1][index]*qs[1] +
                               nanv[0][2][index]*qs[2]);

               stemp[1] = fabs(nanv[1][0][index]*qs[0] +
                               nanv[1][1][index]*qs[1] +
                               nanv[1][2][index]*qs[2]);

               stemp[2] = fabs(nanv[2][0][index]*qs[0] +
                               nanv[2][1][index]*qs[1] +
                               nanv[2][2][index]*qs[2]);
 
/*
 *             Find min index2
 */
               FindAbsMin(stemp, 3, &minVal, &index2);
 
               switch (index2) {
                   case 0:
                       offset1 = 1;
                       offset2 = 2;
                       break;
                   case 1:
                       offset1 = 0;
                       offset2 = 2;
                       break;
                   default:
                       offset1 = 0;
                       offset2 = 1;
                       break;
               }

               vec1[0] = nanv[offset1][0][index];
               vec1[1] = nanv[offset1][1][index];
               vec1[2] = nanv[offset1][2][index];

               vec2[0] = nanv[offset2][0][index];
               vec2[1] = nanv[offset2][1][index];
               vec2[2] = nanv[offset2][2][index];

               DecompVec(qs, vec1, vec2, outVec);

               thread_areaSwept[index][offset1+4] += outVec[0];
               thread_areaSwept[index][offset2+4] += outVec[1];

#if 0
/*
 *             Some stuff for debugging
 */
               printf("thread_Ltot ---------------- \n");
               for (m = 0; m < 4; m++) {
                   printf("%e %e %e %e \n",
                          thread_Ltot[m][0], thread_Ltot[m][1],
                          thread_Ltot[m][2], thread_Ltot[m][3]);
               }
   
               printf("thread_areaSwept ----------- \n");
               for (m = 0; m < 4; m++) {
                   printf(" %e %e %e %e %e %e %e\n",
                          thread_areaSwept[m][0], thread_areaSwept[m][1],
                          thread_areaSwept[m][2], thread_areaSwept[m][3],
                          thread_areaSwept[m][4], thread_areaSwept[m][5],
                          thread_areaSwept[m][6]);
               }
#endif
           }  /* for (segIndex = 0; segIndex < numSeg; ...) */
        }  /* for (nodeIndex = 0; nodeIndex < home->newNodeKeyPtr; ...) */

/*
 *      At this point, each thread has a portion of several arrays.
 *      The threads now update the global arrays with their contributions
 *      in a 'critical' block so there's no conflicts among the threads.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_DPS_BCC_A)
#endif
        {
            for (tk = 0; tk < 3; tk++) {
                for (tl = 0; tl < 3; tl++) {
                    dstn[tk][tl] += thread_dstn[tk][tl];
                    dspn[tk][tl] += thread_dspn[tk][tl];
                }
            }

            for (tk = 0; tk < 4; tk++) {

                for (tl = 0; tl < 7; tl++) {
                    areaSwept[tk][tl] += thread_areaSwept[tk][tl];
                }

                for (tl = 0; tl < 11; tl++) {
                    Ltot[tk][tl] += thread_Ltot[tk][tl];
                }
            }

            for (tk = 0; tk < param->numBurgGroups; tk++) {
                param->partialDisloDensity[tk] += thread_partialDensity[tk];
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
        for (i = 0; i < 4; i++) {
            // plane 1 edge + screw
            Ltot[i][8]  = sqrt((Ltot[i][1] * Ltot[i][1])+
                               (Ltot[i][5] * Ltot[i][5]));

            // plane 2 edge + screw
            Ltot[i][9]  = sqrt((Ltot[i][2] * Ltot[i][2])+
                               (Ltot[i][6] * Ltot[i][6]));

            // plane 3 edge + screw
            Ltot[i][10] = sqrt((Ltot[i][3] * Ltot[i][3])+
                               (Ltot[i][7] * Ltot[i][7]));
        }

/*
 *      Density decomposition
 */
        factor0 = param->simVol*bmagsq;

        for (i = 0; i < 4; i++) {

            for (j = 0; j < 11; j++) {
                param->BCC_dLtotDecomp[i][j] = Ltot[i][j] / factor0; 
            }

            for (j = 0; j < 7; j++) {
                param->BCC_dfluxtot[i][j] = areaSwept[i][j] / param->simVol /
                                            param->realdt;
            }
        }
 
#if 0
        printf("Ltot ---------------- \n");
        for (i = 0; i < 4; i++) {
            printf("%e %e %e %e \n", Ltot[i][0], Ltot[i][1],
                   Ltot[i][2], Ltot[i][3]);
        }
#endif


/*
 *      accumulate delta strain this subcycle into delta strain this cycle
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
        MPI_Allreduce(param->delpStrain, gdstn, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(param->delpSpin  , gdspn, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        for (i = 0; i < 6; i++) {
            param->delpStrain[i] = gdstn[i];
            param->delpSpin  [i] = gdspn[i];
        }

/*
 *      Flux decomposition
 *
 *      NOTE: dfluxtot is the 'areaSwept'. Divided by realdt, gives
 *            real "flux" unit.
 *      28 = 4 Burgers vectors  x 7 components of flux
 *      44 = 4 Burgers vectors  x 11 components of density
 */
        MPI_Allreduce((real8 *) param->BCC_dLtotDecomp, (real8 *) param->BCC_LtotDecomp, 44, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce((real8 *) param->BCC_dfluxtot   , (real8 *) param->BCC_fluxtot   , 28, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
/*
 *      For serial compilation, no need to accumulate values from remote
 *      processors, just copy the data into the param struct.
 */
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 11; j++) {
                param->BCC_LtotDecomp[i][j] = param->BCC_dLtotDecomp[i][j];
            }
            for (j = 0; j < 7; j++) {
                param->BCC_fluxtot[i][j] = param->BCC_dfluxtot[i][j];
            }
        }

        for (i = 0; i < 6; i++) {
            gdstn[i] = param->delpStrain[i];
            gdspn[i] = param->delpSpin[i];
        }
#endif /* ifdef PARALLEL */


/*
 *      For debugging only.  See comments further on...
 */
#if 0
        if (home->myDomain == 0) {
            int   k, l;
            real8 pstn[3][3];
            real8 dyad[3][3], dyad0[3][3], dyad1[3][3], dyad2[3][3];

            Matrix33_Zero(pstn);

/*
 *          Strain increment using decomposition here
 */
            for (i = 0; i < 4; i++) {

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
                         pstn[k][l] += param->BCC_fluxtot[i][0] *
                                      (dyad[l][k]+dyad[k][l]) *
                                      0.5 * param->realdt;
                    }
                }
              
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
                  
/*
 *                      Volume is already taken care of when dfluxtot is
 *                      calculated above
 */
                        pstn[k][l] += ((param->BCC_fluxtot[i][1]+param->BCC_fluxtot[i][4]) * (dyad0[l][k]+dyad0[k][l]) +
                                       (param->BCC_fluxtot[i][2]+param->BCC_fluxtot[i][5]) * (dyad1[l][k]+dyad1[k][l]) +
                                       (param->BCC_fluxtot[i][3]+param->BCC_fluxtot[i][6]) * (dyad2[l][k]+dyad2[k][l])) * param->realdt * 0.5;
                    }
                }
            }  /* for (i=0;i<4; etc. */

/* 
 *          The dstn tensor is the total sum of all decomposed strain flux
 *          components.  The sum should be identical to delpStrain if the
 *          decomposition is done correctly .. To check if both methods
 *          produce the same result ..  Do not delete the print statements.
 *
 *          NOTE: If b's are not normalized in the input deck, care must be
 *                taken when comparing these two values
 */
#if 0
            printf("--simVoll = %e -----\n", param->simVol);
            printf("dstn[3][3]----------\n");
            for (i = 0; i < 3; i++) {
                printf("  %.12e %.12e %.12e \n",
                       dstn[i][0], dstn[i][1], dstn[i][2]);                        
            }

            printf("pstn[3][3]-simVol=%e--------\n", param->simVol);
            for (i = 0; i < 3; i++) {
                printf("  %.12e %.12e %.12e \n",
                       pstn[i][0], pstn[i][1], pstn[i][2]);
            }
#endif

#ifdef PARALLEL
#if 0
           printf("gdstn= %.8e %.8e %.8e %.8e %.8e %.8e \n",gdstn[0],
                  gdstn[1],gdstn[2],gdstn[3],gdstn[4],gdstn[5]);
           printf("pstn = %.8e %.8e %.8e %.8e %.8e %.8e \n",pstn[0][0],
                  pstn[1][1],pstn[2][2],pstn[1][2],pstn[0][2],pstn[0][1]);
#endif
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
