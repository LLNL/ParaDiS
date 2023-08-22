/**************************************************************************
 *
 *                               ALERT!
 *
 *    THE FUNCTIONS IN THIS MODULE ARE NOT CURRENTLY INVOKED BY ANY
 *    OF THE PARADIS CODE.  THE FUNCTIONS HAVE BEEN DEFINED FOR FUTURE
 *    USE.
 *
 *************************************************************************/


/***************************************************************************
 *
 *      Module:      StresssDueToSemiInfSegAnisotropic.cc
 *
 *      Description: Contains functions needed to calculate anisotropic
 *                   stress at a given point from a specified dislocation
 *                   segment.
 *
 *
 *                   This routines uses only the single integral
 *
 *                   Jz_ijk = \int (R . c12)^i (R . rotMatrixRow3)^j dz
 *                                ----------------------
 *                                        Ra^k
 *
 **************************************************************************/

#include "mpi_portability.h"

#include "Home.h"

#if defined ANISOTROPIC && defined _ARLFEM
/***************************************************************************
 *
 *      Function:    StressDueToSemiInfSegAnisotropic()
 *
 *      Description: Uses anisotropic elasticity to calculate the stress
 *                   at a given point from a specified dislocation
 *                   segment.
 *
 *      Parameters:
 *          IN: point{xyz} Position at which to calculate stress from the
 *                         <p1>,<p2> segment.
 *          IN: p1{xyz}    First endpoint of segment
 *          IN: p2{xyz}    Position of a point along the semi-infinite segment
 *                         starting at p1
 *          IN: b{xyz}     Segment beurgers vector from <p1> to <p2>
 *          IN: a          Core radius
 *          IN: qMax       Base factor for determining the number of terms
 *                         in the expansion for spherical harmonics.  The
 *                         total number of terms is 2 * qMax + 1.
 *          OUT: sigma     Resulting stress. Components ordered as follows:
 *                           [0] = stressxx
 *                           [1] = stressyy
 *                           [2] = stresszz
 *                           [3] = stressyz
 *                           [4] = stressxz
 *                           [5] = stressxy
 *
 **************************************************************************/
void StresssDueToSemiInfSegAnisotropic(Home_t *home,
                                       real8 pointx, real8 pointy, real8 pointz,
                                       real8 p1x, real8 p1y, real8 p1z,
                                       real8 p2x, real8 p2y, real8 p2z,
                                       real8 bx,  real8 by,  real8 bz,
                                       real8 a, int qMax,
                                       real8 sigma[6])
{
        int           qMaxp1 = qMax+1;
        int           qMaxp2 = qMax+2;
        int           twoqMaxp3 = 2*qMax+3;
        int           i, j, k, m, n, ii, jj, kk, iii, jjj, twoi;
        int           local, local2, local3, local4;
        int           offsetG;
        int           ind, ind1, indm1;
        int           FqRealNumElem, FqImagNumElem;
        short         *RemArray, *ImmArray;
        short         *RelArray, *ImlArray;
        short         *Re18Array, *Im18Array;
        short         *ReiiArray, *ImiiArray;
        real8         *FqReal, *FqImag, *rotMatrixRow3;
        real8         dX, dY, dZ;
        real8         invL;
        real8         d2;
        real8         a2, a2_d2, a2d2inv;
        real8         delta;
        real8         thetad;
        real8         Rasum;
        real8         onemfactor;
        real8         Ap12, Ap23, Ap31, Ap41, Ap43, Ap51, Ap52, Ap62, Ap63;
        real8         bt11, bt12, bt13;
        real8         bt21, bt22, bt23;
        real8         bt31, bt32, bt33;
        real8         tmp1, tmp2;
        real8         t[3];
        real8         y0;
        real8         R[3];
        real8         R2;
        real8         Ra;
        real8         factor[qMaxp2];
        real8         factorinv[qMaxp2];
        real8         dfactor[qMaxp2];
        real8         Rdsin[twoqMaxp3];
        real8         Rainv[qMaxp1];
        real8         Ra2inv;
        real8         indd[twoqMaxp3];
        real8         indtd[twoqMaxp3];
        real8         Ap[6][3];
        real8         G[18];
        real8         ApGaux[3][3];
        real8         ApG[6];
        real8         (*ecMatrix)[6];
        real8          reJyzt, imJyzt;
        complex8 *rotMatrix12;
        complex8 alpha;
        complex8 gammad;
        complex8 ctmp;
        complex8 afactor[qMaxp2];
        complex8 inda[twoqMaxp3];
        complex8 indgd[twoqMaxp3];
        complex8 Rdpin[twoqMaxp3];
        complex8 RdpRds[twoqMaxp3][twoqMaxp3];
        complex8 Jyzt0[twoqMaxp3][2];
        complex8 Jyzt1[twoqMaxp3][2];
        complex8 Jyzt2[twoqMaxp3][2];

/*
 *      Need some stuff from 'home'
 */
        FqReal = home->anisoVars.FqReal;
        FqImag = home->anisoVars.FqImag;

        RemArray = home->anisoVars.sphericalHarmonicsRemKept;
        ImmArray = home->anisoVars.sphericalHarmonicsImmKept;
        RelArray = home->anisoVars.sphericalHarmonicsRelKept;
        ImlArray = home->anisoVars.sphericalHarmonicsImlKept;
        Re18Array= home->anisoVars.sphericalHarmonicsRe18Kept;
        Im18Array= home->anisoVars.sphericalHarmonicsIm18Kept;
        ReiiArray= home->anisoVars.sphericalHarmonicsReiKept;
        ImiiArray= home->anisoVars.sphericalHarmonicsImiKept;

        rotMatrix12 = home->anisoVars.anisoRotMatrix12;
        rotMatrixRow3  = &home->anisoVars.anisoRotMatrix[2][0];

        ecMatrix = home->param->elasticConstantMatrix;

/*
 *      Initialization
 */
        for (i = 0; i < 18; i++) {
            G[i] = 0.0;
        }

/*
 *      Not all elements are set, so initialize all to zero at start
 */

        for (i = 0; i < twoqMaxp3; i++) {
            Jyzt0[i][0] = 0.0; Jyzt0[i][1] = 0.0;
            Jyzt1[i][0] = 0.0; Jyzt1[i][1] = 0.0;
            Jyzt2[i][0] = 0.0; Jyzt2[i][1] = 0.0;
        }

        for (i = 0; i < twoqMaxp3; i++) {
            for (j = 0; j < twoqMaxp3; j++) {
                RdpRds[i][j] = 0.0;
            }
        }

/*
 *      Define line direction
 */
        t[0] = p2x - p1x;
        t[1] = p2y - p1y;
        t[2] = p2z - p1z;

        NormalizeVec(t);

/*
 *      Compute integrals via recurrence
 */
        alpha = DotProduct(t, rotMatrix12);
        delta = DotProduct(t, rotMatrixRow3);

        factorinv[0]=-1;
        for (i = 1; i < qMaxp2; i++) {
            factorinv[i] = factorinv[i-1]+2;
        }

        for (i = 0; i < qMaxp2; i++) {
            factor[i] = 1.0 / factorinv[i];
        }

        for (i = 0; i < qMaxp2; i++) {
            afactor[i] = factor[i] * alpha;
            dfactor[i] = factor[i] * delta;
        }

        R[0] = pointx - p1x;
        R[1] = pointy - p1y;
        R[2] = pointz - p1z;

        y0 = DotProduct(R, t);

        Rdpin[0] = 1.0;
        Rdpin[1] = DotProduct(R, rotMatrix12);

        Rdsin[0] = 1.0;
        Rdsin[1] = DotProduct(R, rotMatrixRow3);

        gammad = Rdpin[1] - y0 * alpha;
        thetad = Rdsin[1] - y0 * delta;

        for (i = 2; i < twoqMaxp3; i++) {
            Rdpin[i] = Rdpin[i-1] * Rdpin[1];
            Rdsin[i] = Rdsin[i-1] * Rdsin[1];
        }


        inda[0] = 0.0;
        indd[0] = 0.0;
        indgd[0] = 0.0;
        indtd[0] = 0.0;

        for (i = 1; i < twoqMaxp3; i++) {
            inda[i] = inda[i-1] + alpha;
            indd[i] = indd[i-1] + delta;
            indgd[i] = indgd[i-1] + gammad;
            indtd[i] = indtd[i-1] + thetad;
        }


        for (i = 0; i < twoqMaxp3; i++) {
            for (j = 0; j < 2*qMaxp2-(i+1); j++) {
                RdpRds[i][j] = Rdpin[i] * Rdsin[j];
            }
        }

        R2 = DotProduct(R, R);

        d2 = R2 - y0*y0;

        a2 = a * a;
        a2_d2 = a2 + d2;

        Ra = sqrt(R2 + a2);

        a2d2inv = 1.0 / a2_d2;
        Rainv[0] = 1.0 / Ra;

        Ra2inv = Rainv[0] * Rainv[0];

        for (i = 1; i < qMaxp1; i++) {
            Rainv[i] = Rainv[i-1] * Ra2inv;
        }

 
/*
 *      Beginning terms of the reccurrence
 */ 

        Jyzt1[1][0] = -(Rainv[0]*y0+1)*a2d2inv;

        Jyzt2[1][0] = Rainv[0]*delta + thetad*Jyzt1[1][0];
        Jyzt2[2][0] = Rainv[0]*alpha + gammad*Jyzt1[1][0];


        local3 = 0;
        local4 = 0;

        reJyzt = creal(Jyzt2[1][0]);

        for (ii = 0; ii < 18; ii++) {
            G[ii] += reJyzt*FqReal[local3];
            local3++;
        }

        imJyzt = cimag(Jyzt2[1][0]);

#if 0
//
//      TEMP FIX BEFORE FULL REMOVAL OF COMPLEX VARIABLES
//
//      The new FqImag array does not contain the first 18 elements
//      because they were all zero, so this code is obsolete?
//
        for (ii = 0; ii < 18; ii++) {
            G[ii] -= imJyzt*FqImag[local4];
            local4++;
        }
#endif


        reJyzt = creal(Jyzt2[2][0]);

        for (ii = 0; ii < 18; ii++) {
            G[ii] += reJyzt*FqReal[local3];
            local3++;
        }

        imJyzt = cimag(Jyzt2[2][0]);

        for (ii = 0; ii < 18; ii++) 
        {
            G[ii] -= imJyzt*FqImag[local4];
            local4++;
        }

/*
 *      Recursive integrals
 */
        local  = 2;
        local2 = 2;

        ind = 1;
        indm1 = 0; 

        for (i = 2; i < qMaxp2; i++) {
            twoi = 2*i;

            onemfactor=(1.0 - (2*i-2)*factor[i]);

            kk = twoi - 2;
            k = kk;

            for (j = 0; j < kk; j++) {

                k--;

                ctmp =-y0 * RdpRds[j][k] * Rainv[i-1];

                Jyzt0[j+1][ind] = a2d2inv * (factor[i] * (ctmp + 
                                  indgd[j] * Jyzt1[j][indm1] +
                                  indtd[k] * Jyzt1[j+1][indm1]) +
                                  onemfactor * Jyzt2[j+1][indm1]);
            }

            jj = twoi-4;

            /* n = 0 */
            j = jj+1;

            tmp2 =-Rdsin[j] * Rainv[i-1];

            Jyzt1[1][ind] = dfactor[i] * (indd[j] * Jyzt1[1][indm1] - tmp2) +
                            thetad * Jyzt0[1][ind];

            kk = twoi-2;
            k = kk;
            for (j = 0; j < kk; j++) {

                k--;

                ctmp =-RdpRds[j][k] * Rainv[i-1];

                Jyzt1[j+2][ind] = afactor[i] * (indd[k] * Jyzt1[j+1][indm1] +
                                                inda[j] * Jyzt1[j][indm1] -ctmp) +
                                  gammad * Jyzt0[j+1][ind];
            }

            /* n = 1 */
            j = jj+2;

            tmp2 =-Rdsin[j] * Rainv[i-1];

            Jyzt2[1][ind] = dfactor[i] * (indd[j] * Jyzt2[1][indm1] - tmp2) +
                            thetad * Jyzt1[1][ind];

            kk = twoi-1;
            k = kk;

            for (j = 0; j < kk; j++) {

                k--;

                ctmp =-RdpRds[j][k] * Rainv[i-1];

                Jyzt2[j+2][ind] = afactor[i] * (indd[k] * Jyzt2[j+1][indm1] +
                                                inda[j] * Jyzt2[j][indm1] -ctmp) +
                                  gammad * Jyzt1[j+1][ind];
            }


            for (jjj = 1; jjj < RemArray[i]; jjj++) {

                j = RelArray[local];
                reJyzt = creal(Jyzt2[j][ind]);

                for (iii = 0; iii < Re18Array[local]; iii++) {
                    ii = ReiiArray[local3];
                    G[ii] += reJyzt*FqReal[local3];
                    local3++;
                }

                local++;  
            }

            for (jjj = 1; jjj < ImmArray[i]; jjj++) {

                j = ImlArray[local2];
                imJyzt = cimag(Jyzt2[j][ind]);

                for (iii = 0; iii < Im18Array[local2]; iii++) {
                    ii = ImiiArray[local4];
                    G[ii] -= imJyzt*FqImag[local4];
                    local4++;
                }

                local2++;
            }

            indm1 = ind;
            ind = 1 - ind;


        }  /* end for (i = 2; i < q; i++) */


/*
 *      Compute pre-factor eps C C b t bp tp in the form of
 *      matrices Ap, B, Bp, A
 */
        bt11 = bx*t[0];
        bt22 = by*t[1];
        bt33 = bz*t[2];

        bt12 = bx*t[1];
        bt23 = by*t[2];
        bt31 = bz*t[0];

        bt21 = by*t[0];
        bt32 = bz*t[1];
        bt13 = bx*t[2];

        Ap12 = bt13;
        Ap62 = bt23;
        Ap51 = bt12;
        Ap23 = bt21;
        Ap43 = bt31;
        Ap31 = bt32;

        Ap63 = bt11 - bt22;
        Ap41 = bt22 - bt33;
        Ap52 = bt33 - bt11;



/*
 *      Compute pre-factor eps C C b t bp tp in the form of
 *      matrices Ap, B, Bp, A
 */
        bt11 = bx*t[0];
        bt22 = by*t[1];
        bt33 = bz*t[2];

        bt12 = bx*t[1];
        bt23 = by*t[2];
        bt31 = bz*t[0];

        bt21 = by*t[0];
        bt32 = bz*t[1];
        bt13 = bx*t[2];

        Ap12 = bt13;
        Ap62 = bt23;
        Ap51 = bt12;
        Ap23 = bt21;
        Ap43 = bt31;
        Ap31 = bt32;

        Ap63 = bt11 - bt22;
        Ap41 = bt22 - bt33;
        Ap52 = bt33 - bt11;

        switch(home->param->materialType) {
            case MAT_TYPE_BCC:
            case MAT_TYPE_FCC:
            case MAT_TYPE_HCP:
/*
 *              For hexagonal or cubic elastic constant matrix
 */
                for (i = 0; i < 6; i++) {
                    for (j = 0; j < 3; j++) {
                        Ap[i][j] = 0.0;
                    }
                }

                for (j = 0; j < 3; j++) {
                    Ap[j][0] = -ecMatrix[j][1]*Ap62 + ecMatrix[j][2]*Ap31;
                    Ap[j][1] =  ecMatrix[j][0]*Ap12 - ecMatrix[j][2]*Ap43;
                    Ap[j][2] = -ecMatrix[j][0]*Ap51 + ecMatrix[j][1]*Ap23;
                }

                Ap[3][0] =  ecMatrix[3][3]*Ap41;
                Ap[3][1] = -ecMatrix[3][3]*Ap23;
                Ap[3][2] =  ecMatrix[3][3]*Ap43;

                Ap[4][0] =  ecMatrix[4][4]*Ap51;
                Ap[4][1] =  ecMatrix[4][4]*Ap52;
                Ap[4][2] = -ecMatrix[4][4]*Ap31;

                Ap[5][0] = -ecMatrix[5][5]*Ap12;
                Ap[5][1] =  ecMatrix[5][5]*Ap62;
                Ap[5][2] =  ecMatrix[5][5]*Ap63;

                break;

            case MAT_TYPE_RHOMBOHEDRAL_VA:
/*
 *              For both general or symmetric elastic constant matrix
 */
                for (j = 0; j < 6; j++) {
                    Ap[j][0] = -ecMatrix[j][1]*Ap62 + ecMatrix[j][2]*Ap31 +
                                ecMatrix[j][3]*Ap41 + ecMatrix[j][4]*Ap51 -
                                ecMatrix[j][5]*Ap12;
                    Ap[j][1] =  ecMatrix[j][0]*Ap12 - ecMatrix[j][2]*Ap43 -
                                ecMatrix[j][3]*Ap23 + ecMatrix[j][4]*Ap52 +
                                ecMatrix[j][5]*Ap62;
                    Ap[j][2] = -ecMatrix[j][0]*Ap51 + ecMatrix[j][1]*Ap23 +
                                ecMatrix[j][3]*Ap43 - ecMatrix[j][4]*Ap31 +
                                ecMatrix[j][5]*Ap63;
                }

                break;

        }  /* end switch(home->param->materialType) */


        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {

                ApGaux[i][j] = 0.0;

                for (k = 0; k < 6; k++) {
                    offsetG = 6*j+k;
                    ApGaux[i][j] += Ap[k][i]*G[offsetG];
                }
            }
        }

        ApG[0] = ApGaux[0][0];
        ApG[1] = ApGaux[1][1];
        ApG[2] = ApGaux[2][2];

        ApG[3] = (ApGaux[1][2] + ApGaux[2][1]);
        ApG[4] = (ApGaux[2][0] + ApGaux[0][2]);
        ApG[5] = (ApGaux[0][1] + ApGaux[1][0]);

        for (i = 0; i < 6; i++) {

            sigma[i] = 0.0;

            for (j = 0; j < 6; j++) {

                sigma[i] += ecMatrix[i][j] * ApG[j];
            }

            sigma[i] *= INV_4_PI_SQUARED;
        }

        return;
}
#endif  // end if defined ANISOTROPIC && _ARLFEM
