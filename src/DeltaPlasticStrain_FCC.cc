/***************************************************************************
 *   
 *      Author:       Sylvie Aubry (based on DeltaPlasticStrain_BCC.c)
 *
 *      Function:     DeltaPlasticStrain_FCC
 *
 *      Description:  Calculate the plastic strain increment
 *                    by dislocation motion in FCC materials.
 *
 *      Public functions:
 *          DeltaPlasticStrain_FCC()
 *
 *      Private functions:
 *          InitDeltaPlasticStrain_FCC()
 *
 ***************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Node.h"
#include "Util.h"
#include "Mobility.h"

/*
 *      We only deal with specific burgers vectors in this function.
 *      Define the number of main burgers vectors we handle and the
 *      number of planes for each burgers vector.
 */
#define NUM_BURG    6
#define NUM_PLANES  3


/***************************************************************************
 *
 *      Function:    InitDeltaPlasticStrain_FCC
 *
 *      Description: Several of the arrays used during the delta plastic
 *                   strain calculations are dynamically calculated, but
 *                   once set, never need to be changed again for the
 *                   duration of the simulation.  This function should be
 *                   invoked during the first call to DeltaPlasticStrain_FCC()
 *                   to initialize those arrays, and never used again
 *                   after that.
 *
 ***************************************************************************/
static void InitDeltaPlasticStrain_FCC(Home_t *home,
                                       real8 nanv[NUM_PLANES][3][NUM_BURG],
                                       real8 tanv[NUM_PLANES][3][NUM_BURG],
                                       real8 burgv[NUM_BURG][3],
                                       real8 bnorm[NUM_BURG][3])
{
        int      bIndex;
        Param_t *param = home->param;


        nanv[0][0][0]=0.5773502692; nanv[0][1][0]=-0.5773502692; nanv[0][2][0]=0.5773502692;
        nanv[1][0][0]=-0.5773502692; nanv[1][1][0]=0.5773502692; nanv[1][2][0]=0.5773502692;
        nanv[2][0][0]=0.0000000000; nanv[2][1][0]=0.0000000000; nanv[2][2][0]=1.0000000000;

        nanv[0][0][1]=0.5773502692; nanv[0][1][1]=0.5773502692; nanv[0][2][1]=0.5773502692;
        nanv[1][0][1]=0.5773502692; nanv[1][1][1]=0.5773502692; nanv[1][2][1]=-0.5773502692;
        nanv[2][0][1]=0.0000000000; nanv[2][1][1]=0.0000000000; nanv[2][2][1]=1.0000000000;

        nanv[0][0][2]=0.5773502692; nanv[0][1][2]=0.5773502692; nanv[0][2][2]=-0.5773502692;
        nanv[1][0][2]=-0.5773502692; nanv[1][1][2]=0.5773502692; nanv[1][2][2]=0.5773502692;
        nanv[2][0][2]=0.0000000000; nanv[2][1][2]=1.0000000000; nanv[2][2][2]=0.0000000000;

        nanv[0][0][3]=0.5773502692; nanv[0][1][3]=0.5773502692; nanv[0][2][3]=0.5773502692;
        nanv[1][0][3]=0.5773502692; nanv[1][1][3]=-0.5773502692; nanv[1][2][3]=0.5773502692;
        nanv[2][0][3]=0.0000000000; nanv[2][1][3]=1.0000000000; nanv[2][2][3]=0.0000000000;

        nanv[0][0][4]=0.5773502692; nanv[0][1][4]=0.5773502692; nanv[0][2][4]=-0.5773502692;
        nanv[1][0][4]=0.5773502692; nanv[1][1][4]=-0.5773502692; nanv[1][2][4]=0.5773502692;
        nanv[2][0][4]=1.0000000000; nanv[2][1][4]=0.0000000000; nanv[2][2][4]=0.0000000000;

        nanv[0][0][5]=0.5773502692; nanv[0][1][5]=0.5773502692; nanv[0][2][5]=0.5773502692;
        nanv[1][0][5]=-0.5773502692; nanv[1][1][5]=0.5773502692; nanv[1][2][5]=0.5773502692;
        nanv[2][0][5]=1.0000000000; nanv[2][1][5]=0.0000000000; nanv[2][2][5]=0.0000000000;

        tanv[0][0][0]=-0.4082482905; tanv[0][1][0]=0.4082482905; tanv[0][2][0]=0.8164965809;
        tanv[1][0][0]=-0.4082482905; tanv[1][1][0]=0.4082482905; tanv[1][2][0]=-0.8164965809;
        tanv[2][0][0]=-0.7071067812; tanv[2][1][0]=0.7071067812; tanv[2][2][0]=0.0000000000;

        tanv[0][0][1]=0.4082482905; tanv[0][1][1]=0.4082482905; tanv[0][2][1]=-0.8164965809;
        tanv[1][0][1]=-0.4082482905; tanv[1][1][1]=-0.4082482905; tanv[1][2][1]=-0.8164965809;
        tanv[2][0][1]=0.7071067812; tanv[2][1][1]=0.7071067812; tanv[2][2][1]=-0.0000000000;

        tanv[0][0][2]=0.4082482905; tanv[0][1][2]=-0.8164965809; tanv[0][2][2]=-0.4082482905;
        tanv[1][0][2]=0.4082482905; tanv[1][1][2]=0.8164965809; tanv[1][2][2]=-0.4082482905;
        tanv[2][0][2]=0.7071067812; tanv[2][1][2]=0.0000000000; tanv[2][2][2]=-0.7071067812;

        tanv[0][0][3]=-0.4082482905; tanv[0][1][3]=0.8164965809; tanv[0][2][3]=-0.4082482905;
        tanv[1][0][3]=0.4082482905; tanv[1][1][3]=0.8164965809; tanv[1][2][3]=0.4082482905;
        tanv[2][0][3]=-0.7071067812; tanv[2][1][3]=0.0000000000; tanv[2][2][3]=-0.7071067812;

        tanv[0][0][4]=0.8164965809; tanv[0][1][4]=-0.4082482905; tanv[0][2][4]=0.4082482905;
        tanv[1][0][4]=-0.8164965809; tanv[1][1][4]=-0.4082482905; tanv[1][2][4]=0.4082482905;
        tanv[2][0][4]=0.0000000000; tanv[2][1][4]=-0.7071067812; tanv[2][2][4]=0.7071067812;

        tanv[0][0][5]=-0.8164965809; tanv[0][1][5]=0.4082482905; tanv[0][2][5]=0.4082482905;
        tanv[1][0][5]=-0.8164965809; tanv[1][1][5]=-0.4082482905; tanv[1][2][5]=-0.4082482905;
        tanv[2][0][5]=-0.0000000000; tanv[2][1][5]=0.7071067812; tanv[2][2][5]=0.7071067812;

        bnorm[0][0]=0.7071067812; bnorm[0][1]=0.7071067812; bnorm[0][2]=0.0000000000;
        bnorm[1][0]=0.7071067812; bnorm[1][1]=-0.7071067812; bnorm[1][2]=0.0000000000;
        bnorm[2][0]=0.7071067812; bnorm[2][1]=0.0000000000; bnorm[2][2]=0.7071067812;
        bnorm[3][0]=0.7071067812; bnorm[3][1]=0.0000000000; bnorm[3][2]=-0.7071067812;
        bnorm[4][0]=0.0000000000; bnorm[4][1]=0.7071067812; bnorm[4][2]=0.7071067812;
        bnorm[5][0]=0.0000000000; bnorm[5][1]=0.7071067812; bnorm[5][2]=-0.7071067812;

        burgv[0][0]=0.7071067812; burgv[0][1]=0.7071067812; burgv[0][2]=0.0000000000;
        burgv[1][0]=0.7071067812; burgv[1][1]=-0.7071067812; burgv[1][2]=0.0000000000;
        burgv[2][0]=0.7071067812; burgv[2][1]=0.0000000000; burgv[2][2]=0.7071067812;
        burgv[3][0]=0.7071067812; burgv[3][1]=0.0000000000; burgv[3][2]=-0.7071067812;
        burgv[4][0]=0.0000000000; burgv[4][1]=0.7071067812; burgv[4][2]=0.7071067812;
        burgv[5][0]=0.0000000000; burgv[5][1]=0.7071067812; burgv[5][2]=-0.7071067812;

/*
 *      If necessary, rotate the nanv and tanv arrays from the crystal frame
 *      to the lab frame
 */
        if (param->useLabFrame) {
            int plane;

            for (plane = 0; plane < NUM_PLANES; plane++) {

                for (bIndex = 0; bIndex < NUM_BURG; bIndex++) {
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

/*
 *      If necessary, rotate the bnorm array from the crystal frame to
 *      the lab frame
 */
        if (param->useLabFrame) {

            for (bIndex = 0; bIndex < NUM_BURG; bIndex++) {
                real8 bnormLab[3];

                Matrix33Vector3Multiply(param->rotMatrix,bnorm[bIndex],bnormLab);
                VECTOR_COPY(bnorm[bIndex], bnormLab);
            }
        }

#if 0
/*
 *      Some stuff for debugging
 */
        for (int i = 0; i < NUM_BURG; i++) {
            printf("bnorm[%d]=%f %f %f\n", i,
                   bnorm[i][0], bnorm[i][1], bnorm[i][2]);
        }
        printf("\n");

        for (int i = 0; i < NUM_BURG; i++) {
            for (int j = 0; j < NUM_PLANES; j++) {
                printf("nanv[%d][%d]=%f %f %f\n", i, j,
                       nanv[j][0][i], nanv[j][1][i], nanv[j][2][i]);
            }
            printf("\n");
        }

        for (int i = 0; i < NUM_BURG; i++) {
            for (int j = 0; j < NUM_PLANES; j++) {
                printf("tanv[%d][%d]=%f %f %f\n", i, j,
                       tanv[j][0][i], tanv[j][1][i], tanv[j][2][i]);
            }
            printf("\n");
        }

#endif
        return;
}


void DeltaPlasticStrain_FCC(Home_t *home)
{
        int          i, j;
        real8        factor0, bmagsq;
        real8        eps = 1.e-12;
        real8        gdspn[6], gdstn[6];
        real8        dstn[3][3], dspn[3][3];
        real8        Ltot[NUM_BURG][11], areaSwept[NUM_BURG][7];
        Param_t      *param;
        static int   initDone = 0;
        static real8 burgv[NUM_BURG][3], bnorm[NUM_BURG][3];
        static real8 nanv[NUM_PLANES][3][NUM_BURG];
        static real8 tanv[NUM_PLANES][3][NUM_BURG];

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
//

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
            InitDeltaPlasticStrain_FCC(home, nanv, tanv, burgv, bnorm);
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
        int    k, l;
        int    nodeIndex;
        int    threadID, threadIterStart, threadIterEnd;
        real8  thread_disloDensity;
        real8  *thread_partialDensity;
        real8  thread_dstn[3][3], thread_dspn[3][3], thread_dyad[3][3];
        real8  thread_Ltot[NUM_BURG][11], thread_areaSwept[NUM_BURG][7];

        thread_disloDensity = 0.0;
        thread_partialDensity = (real8 *)calloc(1, param->numBurgGroups *
                                               sizeof(real8));

        Matrix33_Zero(thread_dstn);
        Matrix33_Zero(thread_dspn);

        memset(thread_Ltot, 0, sizeof(thread_Ltot));
        memset(thread_areaSwept, 0, sizeof(thread_areaSwept));

        GetThreadIterationIndices(home->newNodeKeyPtr, &threadID,
                                  &threadIterStart, &threadIterEnd);

/*
 *      Each thread now loops over its portion of the node list
 */
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
               int    offset1, offset2;
               int    index, index2;
               int    burgGroup, negCount, zeroCount;
               int    planeIndex;
               real8  bx, by, bz;
               real8  ex, ey, ez;
               real8  nx, ny, nz;
               real8  tx, ty, tz;
               real8  hhx, hhy, hhz;
               real8  tmpx, tmpy, tmpz;
               real8  deltax1x, deltax1y, deltax1z;
               real8  deltax2x, deltax2y, deltax2z;
               real8  delxnode, delynode, delznode;
               real8  delxnbr, delynbr, delznbr;
               real8  minVal;
               real8  deltaDensity;
               real8  edgeDensity;
               real8  sbsign, sb;
               real8  segd, segd2;
               real8  size, tmpmax;
               real8  bCryst[3];
               real8  zeta[3];
               real8  seg[3], seg2[3];
               real8  segleft[3], segleft2[3], qs[3];
               real8  vec1[3], vec2[3], outVec[2];
               real8  Ltemp[NUM_PLANES], Ltemp2[2];
               real8  stemp[NUM_PLANES];
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
 *             For FCC mobility types we track density of the following
 *             groupings of burgers vectors.  Any burgers vectors not
 *             mentioned are ignored.
 *
 *                 group #     burgers vector types
 *                   0         [ 1 1 0] [-1-1 0]
 *                   1         [-1 1 0] [ 1-1 0]
 *                   2         [ 1 0 1] [-1 0-1]
 *                   3         [-1 0 1] [ 1 0-1]
 *                   4         [ 0 1 1] [ 0-1-1]
 *                   5         [ 0-1 1] [ 0 1-1]
 *                   6         all others
 */
               burgGroup = -1;

               zeroCount = (bCryst[X] == 0.0) +
                           (bCryst[Y] == 0.0) +
                           (bCryst[Z] == 0.0);

               negCount  = (bCryst[X] < 0.0) +
                           (bCryst[Y] < 0.0) +
                           (bCryst[Z] < 0.0);

               switch (zeroCount) {
               case 1:

                   tx = bCryst[X];
                   ty = bCryst[Y];
                   tz = bCryst[Z];

                   if (negCount > 1) {
                       tx = -bCryst[X];
                       ty = -bCryst[Y];
                       tz = -bCryst[Z];
                   }

                   /* types  [ 1  1  0] [-1 -1  0] */
                   if ((tx > 0.0) && (ty > 0.0)) burgGroup = 0;

                   /* types  [ 1 -1  0] [-1  1  0] */
                   else if ((tx != 0.0) && (ty != 0.0) &&
                            ((tx * ty) < 0.0)) burgGroup = 1;

                   /* types  [ 1  0  1] [-1  0 -1] */
                   else if ((tx > 0.0) && (tz > 0.0)) burgGroup = 2;

                   /* types  [-1  0  1] [ 1  0 -1] */
                   else if ((tx != 0.0) && (tz != 0.0) &&
                            ((tx * tz) < 0.0)) burgGroup = 3;

                   /* types  [ 0  1  1] [ 0 -1 -1] */
                   else if ((ty > 0.0) && (tz > 0.0)) burgGroup = 4;

                   /* types  [ 0  1 -1] [ 0 -1  1] */
                   else if ((ty != 0.0) && (tz != 0.0) &&
                            ((ty * tz) < 0.0)) burgGroup = 5;

                   break;

               default:
/*
 *                 All other types of burgers vectors get grouped together
 */
                   burgGroup = 6;
                   break;

               }  /* end switch (zeroCount) */


/*
 *             And, if the burger's vector fell into one of groups of
 *             burgers vectors whose density we're tracking, increment
 *             that group's density value.
 */
               if (burgGroup >= 0) {
                   thread_partialDensity[burgGroup] += deltaDensity;
               }

/*
 *             If this is not a [0 1 1] type burgers vector, skip it.
 */
               if ((fabs(bCryst[X])*fabs(bCryst[Y])*fabs(bCryst[Z])) > eps) {
                   continue;
               }

               if (fabs(bCryst[X]) < eps && fabs(bCryst[Y]) < eps){
                   continue;
               }

               if (fabs(bCryst[Y]) < eps && fabs(bCryst[Z]) < eps){
                   continue;
               }

               if (fabs(bCryst[X]) < eps && fabs(bCryst[Z]) < eps){
                   continue;
               }

/*
 *             max index for NUM_BURG burgers vector for now
 */
               index = -1;
               tmpmax = 0.0;

               for (m = 0; m < NUM_BURG; m++) {
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

               sbsign = bnorm[index][0]*bx +
                        bnorm[index][1]*by +
                        bnorm[index][2]*bz;

               sb = ((sbsign < 0) ? -1.0 : 1.0);

/*
 *             Find the index of the appropriate glide plane for
 *             this segment.
 */
               planeIndex = -1;
               tmpmax = 0.0;

               for (m = 0; m < NUM_PLANES; m++) {
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
               thread_Ltot[index][0] += edgeDensity;  // total edge density
               thread_Ltot[index][4] += fabs(segd);   // total screw density
               thread_Ltot[index][5+planeIndex] += fabs(segd);

/*
 *             min index2
 */
               for (m = 0; m < NUM_PLANES; m++) {
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

               areaSwept[index][0] += sb*0.5*(tmpx*bnorm[index][0] +
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
                       tanv[offset1][2][index], deltax2x, deltax2y,
                       deltax2z, &tmpx, &tmpy, &tmpz);

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
 *             Screw portion, first...
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
               FindAbsMin(stemp, NUM_PLANES, &minVal, &index2);

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


/********************* the other half ****************************/
 
/*
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
               for (m = 0; m < NUM_PLANES; m++) {
                   Ltemp[m] = fabs(tanv[m][0][index]*segleft2[0] +
                                   tanv[m][1][index]*segleft2[1] +
                                   tanv[m][2][index]*segleft2[2]);
               }

               FindMin(Ltemp, NUM_PLANES, &minVal, &index2);
 
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
 *             Screw portion, second...
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
               FindMin(stemp, NUM_PLANES, &minVal, &index2);
 
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

	    }  /* for (segIndex = 0; segIndex < nc; ...) */
        }  /* for (nodeIndex = 0; nodeIndex < home->newNodeKeyPtr; ...) */

/*
 *      At this point, each thread has a portion of several arrays.
 *      The threads now update the global arrays with their contributions
 *      in a 'critical' block so there's no conflicts among the threads.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_DPS_FCC_A)
#endif
        {
            for (k = 0; k < 3; k++) {
                for (l = 0; l < 3; l++) {
                    dstn[k][l] += thread_dstn[k][l];
                    dspn[k][l] += thread_dspn[k][l];
                }
            }

            for (k = 0; k < NUM_BURG; k++) {
                for (l = 0; l < 7; l++) {
                    areaSwept[k][l] += thread_areaSwept[k][l];
                }
                for (l = 0; l < 11; l++) {
                    Ltot[k][l] += thread_Ltot[k][l];
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

#if 0
/*
 *      Some stuff for debugging
 */
        printf("Ltot ---------------- \n");
        for (i = 0; i < NUM_BURG; i++) {
            printf("%e %e %e %e \n", Ltot[i][0], Ltot[i][1],
                   Ltot[i][2], Ltot[i][3]);
        }
        printf("areaSwept printing ---- \n");
        for (i = 0; i < NUM_BURG; i++) {
            printf(" %e %e %e %e %e %e %e\n",
                   areaSwept[i][0], areaSwept[i][1],
                   areaSwept[i][2], areaSwept[i][3],
                   areaSwept[i][4], areaSwept[i][5],
                   areaSwept[i][6]);
        }
#endif

/*
 *      We have per-plane density for edges and screws, but we want 
 *      need to sum the edge and screw values for per-pane totals
 */
        for (i = 0; i < NUM_BURG; i++) {
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

        for (i = 0; i < NUM_BURG; i++) {

            for (j = 0; j < 11; j++) {
                param->FCC_dLtotDecomp[i][j] = Ltot[i][j] / factor0; 
            }

            for (j = 0; j < 7; j++) {
                param->FCC_dfluxtot[i][j] = areaSwept[i][j] / param->simVol / 
                                            param->realdt;
            }
        }
 
#if 0
        printf("Ltot ---------------- \n");
        for (i = 0; i < NUM_BURG; i++) {
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
        MPI_Allreduce((real8 *) param->FCC_dLtotDecomp, (real8 *) param->FCC_LtotDecomp, NUM_BURG*11, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce((real8 *) param->FCC_dfluxtot   , (real8 *) param->FCC_fluxtot   , NUM_BURG*7 , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
/*
 *      For serial compilation, no need to accumulate values from remote
 *      processors, just copy the data into the param struct.
 */
        for (i = 0; i < NUM_BURG; i++) {
            for (j = 0; j < 11; j++) {
                param->FCC_LtotDecomp[i][j] = param->FCC_dLtotDecomp[i][j];
            }
            for (j = 0; j < 7; j++) {
                param->FCC_fluxtot[i][j] = param->FCC_dfluxtot[i][j];
            }
        } 

        for (i = 0; i < 6; i++) {
            gdstn[i] = param->delpStrain[i];
            gdspn[i] = param->delpSpin[i];
        }
#endif /* ifdef PARALLEL */


#if 0
/*
 *      For debugging only.  See comments further on
 */
        if (home->myDomain == 0) {
            int   k, l;
            real8 pstn[3][3];
            real8 val[3][3];
            real8 dyad[3][3], dyad0[3][3], dyad1[3][3], dyad2[3][3];

/*
 *          Strain increment using decomposition here
 */
            for (i = 0; i < NUM_BURG; i++) {

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
                        pstn[k][l] += param->FCC_fluxtot[i][0] *
                                      (dyad[l][k]+dyad[k][l]) * 0.5 /
                                      param->simVol;
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
                  
                        pstn[k][l] += ((param->FCC_fluxtot[i][1]+param->FCC_fluxtot[i][4]) * (dyad0[l][k]+dyad0[k][l]) +
                                       (param->FCC_fluxtot[i][2]+param->FCC_fluxtot[i][5]) * (dyad1[l][k]+dyad1[k][l]) +
                                       (param->FCC_fluxtot[i][3]+param->FCC_fluxtot[i][6]) * (dyad2[l][k]+dyad2[k][l])) * 0.5 / param->simVol;
                    }
                }
            }  /* for (i=0;i<NUM_BURG; etc. */

/* 
 *          dstn2 tensor the total sum of all decomposed strain flux components
 *          The sum should be identical to delpStrain if the decomposition is
 *          done correctly .. To check if both methods produce the same
 *          result ..  Do not delete the print statements.
 */
#if 0
            printf("--simVoll = %e -----\n", param->simVol);
            printf("--realdt = %e -----\n", param->realdt);

            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    val[i][j] = pstn[i][j]*(param->simVol * param->realdt);
                }
            }

            printf("gdstn= %.8e %.8e %.8e %.8e %.8e %.8e \n",
                   gdstn[0], gdstn[1], gdstn[2],
                   gdstn[3], gdstn[4], gdstn[5]);
            printf("pstn = %.8e %.8e %.8e %.8e %.8e %.8e \n",
                   val[0][0], val[1][1], val[2][2],
                   val[1][2], val[0][2], val[0][1]);

            printf("gdstn/pstn = %.8e %.8e %.8e %.8e %.8e %.8e \n",
                   gdstn[0]/val[0][0], gdstn[1]/val[1][1], gdstn[2]/val[2][2],
                   gdstn[3]/val[1][2], gdstn[4]/val[0][2], gdstn[5]/val[0][1]);

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
