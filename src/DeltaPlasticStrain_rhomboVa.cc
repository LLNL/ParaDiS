/***************************************************************************
 *   
 *      Author:       Moono Rhee
 *
 *      Function:     DeltaPlasticStrain_rhomboVa.c
 *
 *      Description:  Calculate the plastic strain increment
 *                    by dislocation motion in rhombohedral
 *                    vanadium.
 *
 *      Public functions:
 *          DeltaPlasticStrain_rhomboVa()
 *
 *      Private functions:
 *          InitDeltaPlasticStrain_rhomboVa()
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
 *      Function:    InitDeltaPlasticStrain_rhomboVa
 *
 *      Description: Several of the arrays used during the delta plastic
 *                   strain calculations are dynamically calculated, but
 *                   once set, never need to be changed again for the
 *                   duration of the simulation.  This function should be
 *                   invoked during the first call to
 *                   DeltaPlasticStrain_rhomboVa() to initialize those
 *                   arrays, and never used again after that.
 *
 ***************************************************************************/
static void InitDeltaPlasticStrain_rhomboVa(Home_t *home,
                                            real8 nanv[3][3][4],
                                            real8 tanv[3][3][4],
                                            real8 burgv[13][3],
                                            real8 bnorm[13][3])
{
        int      i, j, k;
        real8    x110, psmall, ppp;
        real8    rbmag, brmag;
        Param_t *param = home->param;


        burgv[0][0] = 1.450000000000000e+000;
        burgv[0][1] = 1.450000000000000e+000;
        burgv[0][2] = 1.450000000000000e+000;

        rbmag = sqrt(burgv[0][0]*burgv[0][0] * 3.00000000000e+000);

        for (i = 0; i < 3; i++) {
            burgv[0][i] = burgv[0][i] / rbmag;
        }
      
        burgv[1][0] =  -1.187000000000000e+000 / rbmag;
        burgv[1][1] =   1.318500000000000e+000 / rbmag;
        burgv[1][2] =   1.318500000000000e+000 / rbmag;

        burgv[2][0] =   1.318500000000000e+000 / rbmag;
        burgv[2][1] =  -1.187000000000000e+000 / rbmag;
        burgv[2][2] =   1.318500000000000e+000 / rbmag;

        burgv[3][0] =   1.318500000000000e+000 / rbmag;
        burgv[3][1] =   1.318500000000000e+000 / rbmag;
        burgv[3][2] =  -1.187000000000000e+000 / rbmag;



#if 0
/*
 *      The following junction burgers vectors are from Tom
 */
burgv(5:13,:)=[ 1.31500000000e-01,  1.31500000000e-01,  2.63700000000e+00;
                1.31500000000e-01,  2.63700000000e+00,  1.31500000000e-01;
                2.63700000000e+00,  1.31500000000e-01,  1.31500000000e-01;
                                0,  2.50550000000e+00, -2.50550000000e+00; 
                2.50550000000e+00,                  0, -2.50550000000e+00;
                2.50550000000e+00, -2.50550000000e+00,                  0;
                2.76850000000e+00,  2.76850000000e+00,  2.63000000000e-01;
                2.76850000000e+00,  2.63000000000e-01,  2.76850000000e+00;
                2.63000000000e-01,  2.76850000000e+00,  2.76850000000e+00]
#endif

        burgv[4][0]  =  1.31500000000e-01;
        burgv[4][1]  =  1.31500000000e-01;
        burgv[4][2]  =  2.63700000000e+00;

        burgv[5][0]  =  1.31500000000e-01;
        burgv[5][1]  =  2.63700000000e+00;
        burgv[5][2]  =  1.31500000000e-01;

        burgv[6][0]  =  2.63700000000e+00;
        burgv[6][1]  =  1.31500000000e-01;
        burgv[6][2]  =  1.31500000000e-01;

        burgv[7][0]  =  0.0;
        burgv[7][1]  =  2.50550000000e+00;
        burgv[7][2]  = -2.50550000000e+00; 

        burgv[8][0]  =  2.50550000000e+00;
        burgv[8][1]  =  0.0;
        burgv[8][2]  = -2.50550000000e+00;

        burgv[9][0]  =  2.50550000000e+00;
        burgv[9][1]  = -2.50550000000e+00;
        burgv[9][2]  =  0.0;

        burgv[10][0] =  2.76850000000e+00;
        burgv[10][1] =  2.76850000000e+00;
        burgv[10][2] =  2.63000000000e-01;

        burgv[11][0] =  2.76850000000e+00;
        burgv[11][1] =  2.63000000000e-01;
        burgv[11][2] =  2.76850000000e+00;

        burgv[12][0] =  2.63000000000e-01;
        burgv[12][1] =  2.76850000000e+00;
        burgv[12][2] =  2.76850000000e+00;

/*
 *      bnorm is normalized burgers vector, for numerical purpose
 */
        for (i = 0; i < 13; i++) { 

            brmag = sqrt(burgv[i][0]*burgv[i][0] +
                         burgv[i][1]*burgv[i][1] +
                         burgv[i][2]*burgv[i][2]);

            for (j = 0; j < 3; j++) {
                bnorm[i][j] = burgv[i][j] / brmag ; 
            }
        }

/*
 *      Normalized plane directions
 */
        x110   = 7.071067811865475e-001;
        psmall = 7.034825319468779e-002;
        ppp    = 7.053549189140372e-001;

/*
 *           |--- 3 planes
 *           |  |-- x,y,z
 *           |  |
 *           |  |  |--- 4 burgers vectors
 *           |  |  |
 *      nanv[x][y][z]    
 *
 *     b1=[111]/sqrt(3)
 */
        nanv[0][0][0] = -x110;  /*b4,p1 */
        nanv[0][1][0] =  x110;
        nanv[0][2][0] =   0.0;

        nanv[1][0][0] =  x110;  /*b3,p1*/
        nanv[1][1][0] =   0.0;
        nanv[1][2][0] = -x110;

        nanv[2][0][0] =   0.0;  /*b2,p1*/
        nanv[2][1][0] = -x110;
        nanv[2][2][0] =  x110;

/* 
 *      normalized: b2= -0.537007682227198 0.596499266231306 0.596499266231306
 *      actual:     b2= -0.4726266298      0.5249858563      0.5249858563
 */
        nanv[0][0][1] =   0.0; /*b1,p3*/
        nanv[0][1][1] = -x110;
        nanv[0][2][1] =  x110;

        nanv[1][0][1] =  ppp;  /*b4,p3*/
        nanv[1][1][1] = -psmall;
        nanv[1][2][1] =  ppp;

        nanv[2][0][1] = -ppp;  /*b3,p3 */ 
        nanv[2][1][1] = -ppp;
        nanv[2][2][1] =  psmall;

/*
 *      normalized: b3=  0.596499266231306 -0.537007682227198 0.596499266231306
 *      actual:     b3 = 0.5249858563      -0.4726266298      0.5249858563 
 */
        nanv[0][0][2] =  x110;   /*b1,p2*/
        nanv[0][1][2] =   0.0;
        nanv[0][2][2] = -x110;

        nanv[1][0][2] =  psmall; /*b2,p3*/
        nanv[1][1][2] = -ppp;
        nanv[1][2][2] = -ppp;

        nanv[2][0][2] = -ppp;    /*b2,p3*/
        nanv[2][1][2] = -ppp;
        nanv[2][2][2] =  psmall;

/*
 *      normalized: b4=  0.596499266231306 0.596499266231306 -0.537007682227198
 *      actual:     b4=  0.5249858563      0.5249858563      -0.4726266298
 */
        nanv[0][0][3] = -x110;   /*b1,p1*/
        nanv[0][1][3] =  x110;
        nanv[0][2][3] =   0.0;

        nanv[1][0][3] =  psmall; /*b3,p2*/
        nanv[1][1][3] = -ppp;
        nanv[1][2][3] = -ppp;

        nanv[2][0][3] =  ppp;    /*b2,p2 */
        nanv[2][1][3] = -psmall;
        nanv[2][2][3] =  ppp;

/*
 *      Tangent vectors: along the other screw direction, analogous to
 *      edge directions for bcc0
 *
 *      Actual angle between short b's = 73.45 degrees
 *      Angle between long b & short b = 67.7459 deg
 *
 *           |--- 3 planes for each burgers vector
 *           |  |-- x,y,z
 *           |  |  4 burgers vectors
 *           |  |  |              
 *      tanv[x][y][z]    
 *
 *      for b1: three tangent directions;
 *      for the other b's: b4,b3,b2, should be done
 */
        tanv[0][0][0] = bnorm[3][0]; /*b4*/
        tanv[0][1][0] = bnorm[3][1];
        tanv[0][2][0] = bnorm[3][2];

        tanv[1][0][0] = bnorm[2][0]; /*b3*/
        tanv[1][1][0] = bnorm[2][1];
        tanv[1][2][0] = bnorm[2][2];

        tanv[2][0][0] = bnorm[1][0]; /*b2*/
        tanv[2][1][0] = bnorm[1][1];
        tanv[2][2][0] = bnorm[1][2];

/*
 *     b2: b1,b4,b3
 */
        tanv[0][0][1] =  bnorm[0][0]; /*b1*/
        tanv[0][1][1] =  bnorm[0][1];
        tanv[0][2][1] =  bnorm[0][2];

        tanv[1][0][1] = -bnorm[3][0]; /*b4*/
        tanv[1][1][1] = -bnorm[3][1];
        tanv[1][2][1] = -bnorm[3][2];

        tanv[2][0][1] = -bnorm[2][0]; /*b3*/
        tanv[2][1][1] = -bnorm[2][1];
        tanv[2][2][1] = -bnorm[2][2];

/*
 *      b3: b1,b4,b2
 */
        tanv[0][0][2] =  bnorm[0][0]; /*b1*/
        tanv[0][1][2] =  bnorm[0][1];
        tanv[0][2][2] =  bnorm[0][2];

        tanv[1][0][2] = -bnorm[3][0]; /*b4*/
        tanv[1][1][2] = -bnorm[3][1];
        tanv[1][2][2] = -bnorm[3][2];

        tanv[2][0][2] = -bnorm[1][0]; /*b2*/
        tanv[2][1][2] = -bnorm[1][1];
        tanv[2][2][2] = -bnorm[1][2];

/*
 *      b4: b1,b3,b2
 */
        tanv[0][0][3] =  bnorm[0][0]; /*b1*/
        tanv[0][1][3] =  bnorm[0][1];
        tanv[0][2][3] =  bnorm[0][2];

        tanv[1][0][3] = -bnorm[2][0]; /*b3*/
        tanv[1][1][3] = -bnorm[2][1];
        tanv[1][2][3] = -bnorm[2][2];

        tanv[2][0][3] = -bnorm[1][0]; /*b2*/
        tanv[2][1][3] = -bnorm[1][1];
        tanv[2][2][3] = -bnorm[1][2];

/*
 *      Check if they are making 75 degrees instead of 105
 */
        for (i = 0; i < 4; i++){

            for (j = 0; j < 3; j++){
                int checkSgn;

                checkSgn = tanv[j][0][i]*bnorm[i][0] +
                           tanv[j][1][i]*bnorm[i][1] +
                           tanv[j][2][i]*bnorm[i][2];

                if (checkSgn < 0.0 ) {

                    for (k = 0; k < 3; k++) {
                        tanv[j][k][i] = -tanv[j][k][i];
                    }
                }
            }
        }

/*
 *      If necessary, rotate the nanv, tanv, and bnorm arrays from the
 *      crystal frame to the lab frame
 */
        if (param->useLabFrame) {
            int plane, bIndex;

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

            for (bIndex = 0; bIndex < 13; bIndex++) {
                real8 bnormLab[3];

                Matrix33Vector3Multiply(param->rotMatrix,bnorm[bIndex],bnormLab);
                VECTOR_COPY(bnorm[bIndex], bnormLab);
            }
        }

        return;
}


void DeltaPlasticStrain_rhomboVa(Home_t *home)
{
        int          i, j, k;
        real8        eps = 1.0e-12;
        real8        dstn[3][3],dspn[3][3];
        real8        Ltot[4][3][2], areaSwept[4][3][3];
        Param_t      *param;
        static int   initDone = 0;
        static real8 burgv[13][3], bnorm[13][3];
        static real8 nanv[3][3][4], tanv[3][3][4];


        param = home->param;

/*
 *      For rhombohedral, we want to keep only two for each slip system,
 *      one for the screw and the other for 75 degree dislocation
 *      The index should change here.
 */ 
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 3; j++) {
                areaSwept[i][j][2] = 0.0;
                for(k = 0; k < 2; k++) {
                    Ltot[i][j][k]    = 0.0;
                    areaSwept[i][j][k] = 0.0;
                }
            }
        }
 
        Matrix33_Zero(dstn);
        Matrix33_Zero(dspn);

        param->disloDensity = 0.0;

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
            InitDeltaPlasticStrain_rhomboVa(home, nanv, tanv, burgv, bnorm);
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
        real8  thread_dstn[3][3], thread_dspn[3][3];
        real8  thread_Ltot[4][3][2], thread_areaSwept[4][3][3];

        thread_disloDensity = 0.0;

        Matrix33_Zero(thread_dstn);
        Matrix33_Zero(thread_dspn);

        memset(thread_Ltot, 0, sizeof(thread_Ltot));
        memset(thread_areaSwept, 0, sizeof(thread_areaSwept));

        GetThreadIterationIndices(home->newNodeKeyPtr, &threadID,
                                  &threadIterStart, &threadIterEnd);

        for (nodeIndex=threadIterStart; nodeIndex<threadIterEnd; nodeIndex++) {
            int    numSeg, segIndex;
            real8  localstrain[3][3];
            Node_t *node;


            if ((node = home->nodeKeys[nodeIndex]) == (Node_t *)NULL) {
                continue;
            }

            Matrix33_Zero(localstrain);

/*
 *          Loop over every segment attached to the node
 */        
            numSeg = node->numNbrs;

            for (segIndex = 0; segIndex < numSeg; segIndex++) {
                int    m, bIndex, normIndex;
                real8  bx, by, bz;
                real8  ex, ey, ez;
                real8  nx, ny, nz;
                real8  hhx, hhy, hhz;
                real8  tmpx, tmpy, tmpz;
                real8  deltax1x, deltax1y, deltax1z;
                real8  deltax2x, deltax2y, deltax2z;
                real8  delxnode, delynode, delznode;
                real8  delxnbr, delynbr, delznbr;
                real8  bxnorm, bynorm, bznorm;
                real8  size, sbsign, sb;
                real8  cc;
                real8  l_ts, l_t75, ldd_ts, ldd_t75;
                real8  l2_ts, l2_t75, ldd2_ts, ldd2_t75;
                real8  tmpmax;
                real8  burgmagLocal, deltaDensity;
                real8  bCryst[3];
                real8  vect[3];
                real8  zeta[3];
                real8  seg[3], seg2[3];
                real8  thread_dyad[3][3];
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
 *              This check is in support of some local site-specific
 *              post-processing only and all *normal* ParaDiS simulations
 *              will simply continue on through this check.
 */
                if (param->ignoreNonLoopFlux) {
                    if (!HAS_ANY_OF_CONSTRAINTS(node->constraint, LOOP_NODE) ||
                        !HAS_ANY_OF_CONSTRAINTS(nbr->constraint,  LOOP_NODE)) {
                        continue;
                    }
                }

/*
 *              For later, we need a copy of the burgers vector guaranteed
 *              to be in the crystal frame
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

                burgmagLocal = sqrt(bx*bx + by *by + bz* bz);

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

#if 0
                for (tk=0;tk<3;tk++) {
                    printf("thread_dyad = %g %g %g \n", thread_dyad[tk][0],
                           thread_dyad[tk][1], thread_dyad[tk][2]);
                }
                printf(" \n");
#endif
           
                for (tk = 0; tk < 3; tk++){
                   for (tl = 0; tl < 3; tl++){
                       real8 tmpDstn;
                       tmpDstn = (thread_dyad[tl][tk] + thread_dyad[tk][tl]) *
                                 0.5 /param->simVol;
                       thread_dstn[tl][tk] += tmpDstn;
                       thread_dspn[tl][tk] += (thread_dyad[tl][tk] -
                                               thread_dyad[tk][tl]) *
                                              0.5/param->simVol;
                       localstrain[tl][tk] += tmpDstn;
                   }
                }
 
                /* Modified for flux decomp calculation 06/22/04 M.Rhee */
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

                deltaDensity = size /
                               (param->simVol*param->burgMag*param->burgMag);

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
 *              Safety check: check for  burg close to zero.
 *
 *              This may not be needed for the rhombohedral, but leave 
 *              it it just in case.
 */

                if ((fabs(bCryst[X]) *
                     fabs(bCryst[Y]) *
                     fabs(bCryst[Z])) < eps) {
                    continue;
                }

/*
 *              flux decomposition
 */
                bxnorm = bx / burgmagLocal;
                bynorm = by / burgmagLocal;
                bznorm = bz / burgmagLocal;

                bIndex = -1;
                tmpmax=0.0;

                for (m = 0; m < 13; m++) {
                    real8 burgTmp;
                    burgTmp = fabs(bnorm[m][0]*bxnorm +
                                   bnorm[m][1]*bynorm +
                                   bnorm[m][2]*bznorm);
                    if (burgTmp > tmpmax) {
                        tmpmax = burgTmp;
                        bIndex = m;
                    }
                }

/*
 *              Skip if bIndex is negative or corresponds to a junction burg
 */
                if ((bIndex > 3) || (bIndex < 0))  {
                    continue;
                }

                sbsign = bnorm[bIndex][0]*bx +
                         bnorm[bIndex][1]*by +
                         bnorm[bIndex][2]*bz;

/*
 *              Determine sign for sign consistency for flux calculation
 */
                sb = ((sbsign < 0) ? -1.0 : 1.0);

/*
 *              Plane normal: nx,ny,nz *should* already normalized but
 *              make sure...
 */
                nx = node->nx[segIndex];
                ny = node->ny[segIndex];
                nz = node->nz[segIndex];

                Normalize(&nx,&ny,&nz);

/*
 *              Corresponding plane normal index
 */
                normIndex = -1;
                tmpmax=0.0;

                for (m = 0; m < 3; m++) {
                    real8 normTmp;
/*
 *                  nanv are all normalized above
 */
                    normTmp = fabs(nanv[m][0][bIndex]*nx +
                                   nanv[m][1][bIndex]*ny +
                                   nanv[m][2][bIndex]*nz);
                    if (normTmp > tmpmax) {
                        tmpmax = normTmp;
                        normIndex = m;
                    }
                }

/*
 *              Sanity check...
 */
                if (normIndex < 0) {
                    continue;
                }

/*
                deltax2 is the vector from old nbr end to new nbr position
                deltax2x = nbr->x - nbr->oldx; 
                deltax2y = nbr->y - nbr->oldy;
                deltax2z = nbr->z - nbr->oldz;
 */

                cc =  bnorm[bIndex][0]*tanv[normIndex][0][bIndex] 
                    + bnorm[bIndex][1]*tanv[normIndex][1][bIndex] 
                    + bnorm[bIndex][2]*tanv[normIndex][2][bIndex];

                ldd_ts = seg[0]*bnorm[bIndex][0] 
                       + seg[1]*bnorm[bIndex][1] 
                       + seg[2]*bnorm[bIndex][2];

                ldd_t75 = seg[0]*tanv[normIndex][0][bIndex] 
                        + seg[1]*tanv[normIndex][1][bIndex] 
                        + seg[2]*tanv[normIndex][2][bIndex];

                l_ts = ((ldd_ts  - ldd_t75 * cc )/ (1.-cc * cc));
                l_t75= ((ldd_t75 - ldd_ts * cc )/ (1.-cc * cc));

                /* printf("l_ts= %g,l_t75=%g\n",l_ts,l_t75); */

                thread_Ltot[bIndex][normIndex][0] += fabs(l_ts);
                thread_Ltot[bIndex][normIndex][1] += fabs(l_t75);

/*
 *              climb flux (nonglide term) first, non glide fluxes Eq. 14
 *
 *                 (l_DD(t+dt) cross delx2) burg/burgmag / 2.0 / delta time
 *
 *              cross of (seg, delatx2), 
 */
/*
 *              l_DD x dx2
 */
                xvector(seg[0], seg[1], seg[2],
                        deltax2x, deltax2y, deltax2z,
                        &tmpx, &tmpy, &tmpz);

                thread_areaSwept[bIndex][normIndex][0] +=
                        sb * 0.5 *
                        (tmpx * bnorm[bIndex][0] +
                         tmpy * bnorm[bIndex][1] +
                         tmpz * bnorm[bIndex][2]);

/*
 *              l_ts part
 */
                vect[0] = bnorm[bIndex][0] * l_ts;
                vect[1] = bnorm[bIndex][1] * l_ts;
                vect[2] = bnorm[bIndex][2] * l_ts;

                xvector(vect[0], vect[1], vect[2],
                        deltax2x, deltax2y, deltax2z,
                        &tmpx, &tmpy, &tmpz);

                thread_areaSwept[bIndex][normIndex][1] +=
                        sb * 0.5 *
                        (tmpx * nanv[normIndex][0][bIndex] +
                         tmpy * nanv[normIndex][1][bIndex] +
                         tmpz * nanv[normIndex][2][bIndex]);

/*
 *              l_t75, cross vector is to get area
 */
                vect[0] = tanv[normIndex][0][bIndex] * l_t75;
                vect[1] = tanv[normIndex][1][bIndex] * l_t75;
                vect[2] = tanv[normIndex][2][bIndex] * l_t75;

                xvector(vect[0], vect[1], vect[2],
                        deltax2x, deltax2y, deltax2z,
                        &tmpx, &tmpy, &tmpz);

                thread_areaSwept[bIndex][normIndex][2] +=
                        sb * 0.5 *
                        (tmpx * nanv[normIndex][0][bIndex] +
                         tmpy * nanv[normIndex][1][bIndex] +
                         tmpz * nanv[normIndex][2][bIndex]);
 
/*
 *              note: flux second  bracket 0: climb
 *                         1: t75, flux caused by 75 degree dislocation
 *                         2: ts, flux caused by screw dislocation
 */

/********************* the other half ****************************/
 

                ldd2_ts  = seg2[0]*bnorm[bIndex][0] 
                         + seg2[1]*bnorm[bIndex][1] 
                         + seg2[2]*bnorm[bIndex][2];

                ldd2_t75 = seg2[0]*tanv[normIndex][0][bIndex] 
                         + seg2[1]*tanv[normIndex][1][bIndex] 
                         + seg2[2]*tanv[normIndex][2][bIndex] ;

                l2_ts = ((ldd2_ts - ldd2_t75 * cc ) / (1.-cc * cc));
                l2_t75= ((ldd2_t75 - ldd2_ts * cc ) / (1.-cc * cc));

                xvector(seg2[0], seg2[1], seg2[2],
                        deltax1x, deltax1y, deltax1z,
                        &tmpx, &tmpy, &tmpz);
/*
 *              climb
 */
                thread_areaSwept[bIndex][normIndex][0] += sb * 0.5 *
                        (tmpx * bnorm[bIndex][0] +
                         tmpy * bnorm[bIndex][1] +
                         tmpz * bnorm[bIndex][2]);
 
/*
 *              l2_ts,dx1
 */
                vect[0] = bnorm[bIndex][0] * l2_ts ;
                vect[1] = bnorm[bIndex][1] * l2_ts ;
                vect[2] = bnorm[bIndex][2] * l2_ts ;

                xvector(vect[0], vect[1], vect[2],
                        deltax1x, deltax1y, deltax1z,
                        &tmpx, &tmpy, &tmpz);
 
                thread_areaSwept[bIndex][normIndex][1] += sb * 0.5 *
                        (tmpx * nanv[normIndex][0][bIndex] +
                         tmpy * nanv[normIndex][1][bIndex] +
                         tmpz * nanv[normIndex][2][bIndex]);

/*
 *              l2_t75,dx1
 */
                vect[0] = tanv[normIndex][0][bIndex] * l2_t75 ;
                vect[1] = tanv[normIndex][1][bIndex] * l2_t75 ;
                vect[2] = tanv[normIndex][2][bIndex] * l2_t75 ;

                xvector(vect[0], vect[1], vect[2],
                        deltax1x, deltax1y,deltax1z,
                        &tmpx, &tmpy, &tmpz);
 
                thread_areaSwept[bIndex][normIndex][2] += sb * 0.5 *
                        (tmpx * nanv[normIndex][0][bIndex] +
                         tmpy * nanv[normIndex][1][bIndex] +
                         tmpz * nanv[normIndex][2][bIndex]);
 
/*
 *              note: flux second index 0: climb
 *                    1: ts, flux caused by screw dislocation
 *                    2: t75, flux caused by 75 degree dislocation
 */

#if 0
                printf("Ltot ---------------- \n");
                for (tk = 0; tk < 4; tk++) {
                    for (tl = 0; tl < 3; tl++) {
                        printf("%e %e\n", thread_Ltot[tk][tl][0],
                               thread_Ltot[tk][tl][1]);
                    }
                }
   
                printf("fluxtot printing ---- \n");
                for (tk = 0; tk < 4; tk++) {
                    for (tl = 0; tl < 3; tl++) {
                        printf(" %e %e %e \n",
                               thread_areaSwept[tk][tl][0],
                               thread_areaSwept[tk][tl][1],
                               thread_areaSwept[tk][tl][2]);
                    }
                    printf("\n");
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
#pragma omp critical (CRIT_DPS_RHOMBO_A)
#endif
        {
            for (tk = 0; tk < 3; tk++) {
                for (tl = 0; tl < 3; tl++) {
                    dstn[tk][tl] += thread_dstn[tk][tl];
                    dspn[tk][tl] += thread_dspn[tk][tl];
                }
            }

            for (tk = 0; tk < 4; tk++) {
                for (tl = 0; tl < 3; tl++) {
                    Ltot[tk][tl][0] += thread_Ltot[tk][tl][0];
                    Ltot[tk][tl][1] += thread_Ltot[tk][tl][1];
                    areaSwept[tk][tl][0] += thread_areaSwept[tk][tl][0];
                    areaSwept[tk][tl][1] += thread_areaSwept[tk][tl][1];
                    areaSwept[tk][tl][2] += thread_areaSwept[tk][tl][2];
                }
            }

            param->disloDensity += thread_disloDensity;

        }  /* end "omp critical" section */

    }  /* end "omp parallel section" */


/*
 *      Density decomposition
 *
 *      Note: The units of param->rhombo_dLtot[i][j][k] are in units of
 *            reference length.
 *
 *            The units of param->rhombo_dfluxtot[i][j][k] are in units of
 *            inverse reference length inverse seconds
 */


        for (i = 0; i < 4; i++) {

           for (j = 0; j < 3; j++) {

              param->rhombo_dLtot[i][j][0] = Ltot[i][j][0];
              param->rhombo_dLtot[i][j][1] = Ltot[i][j][1];

              param->rhombo_dfluxtot[i][j][0] = areaSwept[i][j][0] /
                                                param->simVol/param->realdt;
              param->rhombo_dfluxtot[i][j][1] = areaSwept[i][j][1] /
                                                param->simVol/param->realdt;
              param->rhombo_dfluxtot[i][j][2] = areaSwept[i][j][2] /
                                                param->simVol/param->realdt;
           }
        }
 
#if 0
        printf("Ltot ---------------- \n");
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 3; j++) {
                printf("%e %e\n", Ltot[i][j][0], Ltot[i][j][1]);
            }
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
                                                                                

/*
 *      We've calculated processor specific values for several items,
 *      now we need to sum the values from all processors.
 */
#ifdef PARALLEL
/*
 *      Sum the delta strain from all processors, and accumulate
 *      into net strain
 */
        real8 gdspn[6] = { 0.0 };
        real8 gdstn[6] = { 0.0 };

        MPI_Allreduce(param->delpStrain, gdstn, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(param->delpSpin  , gdspn, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        for (i = 0; i < 6; i++) {
            param->delpStrain[i] = gdstn[i];
            param->delpSpin  [i] = gdspn[i];
        }

/*
 *      Flux decomposition, 4x3x2=24, 4x3x3=36
 */
        MPI_Allreduce((real8 *) param->rhombo_dLtot   , (real8 *) param->rhombo_Ltot   , 24, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce((real8 *) param->rhombo_dfluxtot, (real8 *) param->rhombo_fluxtot, 36, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
/*
 *      For serial compilation, no need to accumulate values from remote
 *      processors, just copy the data into the param struct.
 */
        for (i = 0; i < 4; i++) {

           for (j = 0; j < 3; j++) {

              param->rhombo_Ltot[i][j][0] = param->rhombo_dLtot[i][j][0];
              param->rhombo_Ltot[i][j][1] = param->rhombo_dLtot[i][j][1];

              param->rhombo_fluxtot[i][j][0] = param->rhombo_dfluxtot[i][j][0];
              param->rhombo_fluxtot[i][j][1] = param->rhombo_dfluxtot[i][j][1];
              param->rhombo_fluxtot[i][j][2] = param->rhombo_dfluxtot[i][j][2];
           }
        }

#endif /* ifdef PARALLEL */


#if 0
/*
 *      For debug only
 *
 *      WARNING!  This code will not function properly!  It depends
 *                on <burgmagLocal> for the pstn calculation, but
 *                <burgmagLocal> is only set in the above loop over
 *                segments of each node, and is not calculated within
 *                this loop.
 */
        if (home->myDomain == 0) {
            int   l;
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
                       pstn[k][l] += ((areaSwept[i][0][0]) * (dyad[l][k]+dyad[k][l]) 
                                     +(areaSwept[i][1][0]) * (dyad[l][k]+dyad[k][l])
                                     +(areaSwept[i][2][0]) * (dyad[l][k]+dyad[k][l]) ) * 0.5 * burgmagLocal / param->simVol;
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


#if 0
                for (i=0;i<3;i++) {
                    printf("b=%d,dyad0= %g %g %g \n", i,
                           dyad0[i][0],dyad0[i][1],dyad0[i][2]);
                }
                printf(" \n");
#endif


                dyad1[0][0] = bnorm[i][0]*nanv[1][0][i];
                dyad1[0][1] = bnorm[i][0]*nanv[1][1][i];
                dyad1[0][2] = bnorm[i][0]*nanv[1][2][i]; 
                dyad1[1][0] = bnorm[i][1]*nanv[1][0][i];
                dyad1[1][1] = bnorm[i][1]*nanv[1][1][i];
                dyad1[1][2] = bnorm[i][1]*nanv[1][2][i]; 
                dyad1[2][0] = bnorm[i][2]*nanv[1][0][i];
                dyad1[2][1] = bnorm[i][2]*nanv[1][1][i];
                dyad1[2][2] = bnorm[i][2]*nanv[1][2][i]; 

#if 0
                for (i=0;i<3;i++) {
                    printf("b=%d,dyad1= %g %g %g \n", i,
                           dyad1[i][0],dyad1[i][1],dyad1[i][2]);
                }
                printf(" \n");
#endif

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
                  
                     pstn[k][l] +=((areaSwept[i][0][1]+areaSwept[i][0][2]) * (dyad0[l][k]+dyad0[k][l]) +
                                   (areaSwept[i][1][1]+areaSwept[i][1][2]) * (dyad1[l][k]+dyad1[k][l]) +
                                   (areaSwept[i][2][1]+areaSwept[i][2][2]) * (dyad2[l][k]+dyad2[k][l]))*0.5 * burgmagLocal / param->simVol;
                    }
                }
            }  /* for (i=0;i<4; etc. */

/* 
 *          dstn tensor the total sum of all decomposed strain flux components
 *          The sum should be identical to delpStrain if the decomposition is
 *          done correctly .. To check if both methods produce the same
 *          result ..  Do not delete the print statements.
 */
#if 0
            printf("--simVoll = %e -----\n", param->simVol);
            printf("dstn[3][3]----------\n");
            for (i = 0; i < 3; i++) {
                printf("  %.12e %.12e %.12e \n", dstn[i][0], dstn[i][1],
                       dstn[i][2]);                        
            }

            printf("pstn[3][3]-simVol=%e--------\n", param->simVol);
            for (i = 0; i < 3; i++) {
                printf("  %.12e %.12e %.12e \n", pstn[i][0], pstn[i][1],
                       pstn[i][2]);
            }

           printf("gdstn= %.8e %.8e %.8e %.8e %.8e %.8e \n",gdstn[0],
                  gdstn[1],gdstn[2],gdstn[3],gdstn[4],gdstn[5]);
           printf("pstn = %.8e %.8e %.8e %.8e %.8e %.8e \n",pstn[0][0],
                  pstn[1][1],pstn[2][2],pstn[1][2],pstn[0][2],pstn[0][1]);

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
