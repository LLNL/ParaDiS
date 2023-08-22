/***************************************************************************
 *
 *      Module:      SemiInfiniteSegSegForceAnisotropicParallel.cc
 *
 *      Description: Contains function to calcluate the anisotropic force
 *                   on a dislocation from another parallel dislocation
 *
 ***************************************************************************/

#include "mpi_portability.h"

#include "Home.h"

#if defined ANISOTROPIC && defined _ARLFEM
#if 0  // Anisotropic forces need to be updated using table for version 3.
/***************************************************************************
 *
 *      Function:    SemiInfiniteSegSegForceAnisotropicParallel()
 *
 *      Description: Calculates the anisotropic force on segment 2 from
 *                   the near-parallel segment 1.
 *
 *      Parameters:
 *          IN:   p1{xyz}       Position of first segment's first endpoint
 *          IN:   p2{xyz}       Point on the semi-infinite dislocation starting
 *                              at p1.
 *          IN:   p3{xyz}       Position of second segment's first endpoint
 *          IN:   p4{xyz}       Position of second segment's second endpoint
 *          IN:   bp{xyz}       Burgers vector of segment 1 from p1 to p2
 *          IN:   b{xyz}        Burgers vector of segment 2 from p3 to p4
 *          IN:   a             Core radius
 *          IN:   qMax          Base factor for determining the number of terms
 *                              in the expansion for spherical harmonics.  The
 *                              total number of terms is 2 * qMax + 1.
 *          OUT:  f3{xyz}       Force on segment 2 at position <p3>
 *          OUT:  f4{xyz}       Force on segment 2 at position <p4>
 *
 ***************************************************************************/

void SemiInfiniteSegSegForceAnisotropicParallel(
                                        Home_t *home,
                                        real8 p1x, real8 p1y, real8 p1z,
                                        real8 p2x, real8 p2y, real8 p2z,
                                        real8 p3x, real8 p3y, real8 p3z,
                                        real8 p4x, real8 p4y, real8 p4z,
                                        real8 bpx, real8 bpy, real8 bpz,
                                        real8 bx,  real8 by,  real8 bz,
                                        real8 a, int qMax, 
                                        real8 *fp3x, real8 *fp3y, real8 *fp3z,
                                        real8 *fp4x, real8 *fp4y, real8 *fp4z)
{
        int           qMaxp1 = qMax+1;
        int           qMaxp2 = qMax+2;
        int           twoqMaxp2 = 2*qMax+2;
        int           twoqMaxp3 = 2*qMax+3;
        int           twoqMaxp4 = 2*qMax+4;

        int           local, local2, local3, local4, offsetG;
        int           i, j, k, m, n, ii, jj, kk, iii;
        int           ind, ind1, indm1;

        short        *RemArray, *ImmArray;
        short        *RelArray, *ImlArray;
        short        *Re18Array, *Im18Array;
        short        *ReiiArray, *ImiiArray;

        real8         (*ecMatrix)[6];
        real8         *FqReal, *FqImag, *rotMatrixRow3;
        real8         dX, dY, dZ, invL;
        real8         coeff34;
        real8         onemfactor, onem2factor;
        real8         a2, d2, a2_d2, a2d2inv;
        real8         temp[4], lntemp[4];

        real8         delta;
        real8         thetad;
        real8         din[twoqMaxp4];
        real8         dfactorTmp, commonsum, tempsuminf;

        real8         bptp11, bptp12, bptp13;
        real8         bptp21, bptp22, bptp23;
        real8         bptp31, bptp32, bptp33;
        real8         Ap12, Ap23, Ap31, Ap41, Ap43, Ap51, Ap52, Ap62, Ap63;
        real8         B12, B14, B16, B23, B24, B25, B31, B35, B36;
        real8         bt11, bt12, bt13;
        real8         bt21, bt22, bt23;
        real8         bt31, bt32, bt33;
        real8         center[3];
        real8         t[3];
        real8         y0, y1, z0;
        real8         ypz[2], ymz[2];
        real8         factor[qMaxp2];
        real8         factorinv[qMaxp2];
        real8         dfactor[qMaxp2];
        real8         Rdsin0[twoqMaxp3], Rdsin1[twoqMaxp3];
        real8         indd[twoqMaxp3];
        real8         indtd[twoqMaxp3];
        real8         R2[2];
        real8         Ra[2];
        real8         Ra2inv[2];
        real8         Rainv0[qMaxp1], Rainv1[qMaxp1];
        real8         Rasum;
        real8         R[2][3];
        real8         A[6][3], Ap[6][3];
        real8         B[3][6], Bp[3][6];
        real8         G3[18], G4[18];
        real8         ApG3[6], ApG4[6];
        real8         ApGaux3[3][3], ApGaux4[3][3];
        real8         reHy, reHt, reHy1, reHy2;
        real8         imHy, imHt, imHy1, imHy2;
        real8         ReFqval, ImFqval;

        complex8 *rotMatrix12;
        complex8 tmp1, tmp2, tmp3, tmpsum, tmpinf1, tmpinf2;
        complex8 alpha, gammad;
        complex8 afactorTmp;

        complex8 Rasuminf, zJyzsuminf;
        complex8 ain[twoqMaxp4], adin[twoqMaxp4][twoqMaxp4];


        complex8 Ht1_11,  Ht1_21;
        complex8 Hyt2_11, Hyt2_21;

        complex8 Jyzsuminf;
        complex8 afactor[qMaxp2];
        complex8 Rdpin0[twoqMaxp3], Rdpin1[twoqMaxp3];
        complex8 inda[twoqMaxp3];
        complex8 indgd[twoqMaxp3];
        complex8 RdpRds0[twoqMaxp3][twoqMaxp3];
        complex8 RdpRds1[twoqMaxp3][twoqMaxp3];

        complex8 Htinf0[twoqMaxp4][2], Htinf1[twoqMaxp4][2];
        complex8 Htinf2[twoqMaxp4][2];
        complex8 Hytinf0[twoqMaxp3][2], Hytinf1[twoqMaxp3][2];
        complex8 Hytinf2[twoqMaxp3][2];
        complex8 Hztinf0[twoqMaxp3][2], Hztinf1[twoqMaxp3][2];
        complex8 Hztinf2[twoqMaxp3][2];

        complex8 Kyztsum0[twoqMaxp3][2], Kyztsum1[twoqMaxp3][2];
        complex8 Kyztsum2[twoqMaxp3][2];
        complex8 Lyztsum0[twoqMaxp3][2], Lyztsum1[twoqMaxp3][2];
        complex8 Lyztsum2[twoqMaxp3][2];
        complex8 Kyztinfsum0[twoqMaxp3][2], Kyztinfsum1[twoqMaxp3][2];

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
            G3[i] = 0.0;
            G4[i] = 0.0;
        }


        for (i = 0; i < twoqMaxp4; i++) {
            Htinf0[i][0] = 0.0;   Htinf0[i][1] = 0.0;
            Htinf1[i][0] = 0.0;   Htinf1[i][1] = 0.0;
            Htinf2[i][0] = 0.0;   Htinf2[i][1] = 0.0;
        }

        for (i = 0; i < twoqMaxp3; i++) {
            Hytinf0[i][0] = 0.0;
            Hytinf0[i][1] = 0.0;
            Hytinf1[i][0] = 0.0;
            Hytinf1[i][1] = 0.0;
            Hytinf2[i][0] = 0.0;
            Hytinf2[i][1] = 0.0;

            Hztinf0[i][0] = 0.0;
            Hztinf0[i][1] = 0.0;
            Hztinf1[i][0] = 0.0;
            Hztinf1[i][1] = 0.0;
            Hztinf2[i][0] = 0.0;
            Hztinf2[i][1] = 0.0;


            Kyztsum0[i][0] = 0.0;
            Kyztsum0[i][1] = 0.0;
            Kyztsum1[i][0] = 0.0;
            Kyztsum1[i][1] = 0.0;
            Kyztsum2[i][0] = 0.0;
            Kyztsum2[i][1] = 0.0;

            Lyztsum0[i][0] = 0.0;
            Lyztsum0[i][1] = 0.0;
            Lyztsum1[i][0] = 0.0;
            Lyztsum1[i][1] = 0.0;
            Lyztsum2[i][0] = 0.0;
            Lyztsum2[i][1] = 0.0;

            Kyztinfsum0[i][0] = 0.0;
            Kyztinfsum0[i][1] = 0.0;
            Kyztinfsum1[i][0] = 0.0;
            Kyztinfsum1[i][1] = 0.0;
        }

/*
 *      Not some of the  elements, the other will be set
 */
        for (i = 0; i < twoqMaxp3; i++) {
            for (j = twoqMaxp4-(i+1); j < twoqMaxp3; j++) {
                RdpRds0[i][j] = 0.0;
                RdpRds1[i][j] = 0.0;
            }
        }

/*
 *      Get line directions and set some line-specific coefficients.
 *      Note: segments are parallel, so we only need line direction
 *      for 1 of the segments.
 */
        dX = p4x - p3x;
        dY = p4y - p3y;
        dZ = p4z - p3z;

        invL = 1.0 / sqrt(dX*dX + dY*dY + dZ*dZ);
        coeff34 = invL * INV_4_PI_SQUARED;

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

        R[0][0] = p3x - p1x;
        R[0][1] = p3y - p1y;
        R[0][2] = p3z - p1z;

        R[1][0] = p4x - p1x;
        R[1][1] = p4y - p1y;
        R[1][2] = p4z - p1z;

        for (i = 0; i < 2; i++) {
            R2[i] = DotProduct(R[i], R[i]);
        }

        center[0] = 0.5*(p1x + 0.5 * (p3x + p4x));
        center[1] = 0.5*(p1y + 0.5 * (p3y + p4y));
        center[2] = 0.5*(p1z + 0.5 * (p3z + p4z));

        y0 =   (p3x - center[0]) * t[0] +
               (p3y - center[1]) * t[1] +
               (p3z - center[2]) * t[2];

        y1 =   (p4x - center[0]) * t[0] +
               (p4y - center[1]) * t[1] +
               (p4z - center[2]) * t[2];

        z0 = -((p1x - center[0]) * t[0] +
               (p1y - center[1]) * t[1] +
               (p1z - center[2]) * t[2]);

        ypz[0] = y0 + z0;
        ypz[1] = y1 + z0;

        ymz[0] = y0 - z0;
        ymz[1] = y1 - z0;

        Rdpin0[0] = 1.0;
        Rdpin0[1] = DotProduct(R[0], rotMatrix12);

        Rdsin0[0] = 1.0;
        Rdsin0[1] = DotProduct(R[0], rotMatrixRow3);

        Rdpin1[0] = 1.0;
        Rdpin1[1] = DotProduct(R[1], rotMatrix12);

        Rdsin1[0] = 1.0;
        Rdsin1[1] = DotProduct(R[1], rotMatrixRow3);

        for (i = 2; i < twoqMaxp3; i++) {
            Rdpin0[i] = Rdpin0[i-1] * Rdpin0[1];
            Rdpin1[i] = Rdpin1[i-1] * Rdpin1[1];

            Rdsin0[i] = Rdsin0[i-1] * Rdsin0[1];
            Rdsin1[i] = Rdsin1[i-1] * Rdsin1[1];
        }

        gammad = Rdpin0[1] - ypz[0] * alpha;
        thetad = Rdsin0[1] - ypz[0] * delta;

        inda[0] = 0.0;
        indd[0] = 0.0;
        indgd[0] = 0.0;
        indtd[0] = 0.0;

        ain[0] = 0.0;
        din[0] = 0.0;
        ain[1] = 1.0;
        din[1] = 1.0;
        ain[2] = alpha;
        din[2] = delta;

        for (i = 1; i < twoqMaxp3; i++) {
            inda[i] = inda[i-1] + alpha;
            indd[i] = indd[i-1] + delta;
            indgd[i] = indgd[i-1] + gammad;
            indtd[i] = indtd[i-1] + thetad;
        }

        for (i = 2; i < twoqMaxp3; i++) {
            ain[i+1] = ain[i] * ain[2];
            din[i+1] = din[i] * din[2];
        }

        adin[0][0] = 0.0;

        for (i = 0; i < twoqMaxp3; i++) {
            int jMax = twoqMaxp4-(i+1);

            for (j = 0; j < jMax; j++) {
                RdpRds0[i][j] = Rdpin0[i] * Rdsin0[j];
                RdpRds1[i][j] = Rdpin1[i] * Rdsin1[j];

                adin[i+1][j+1] = ain[i+1] * din[j+1];
            }
        }

        d2 = 0.5 * (R2[0] - ypz[0]*ypz[0] + 
                    R2[1] - ypz[1]*ypz[1]);

        a2 = a * a;
        a2_d2 = a2 + d2;
        a2d2inv = 1.0 / a2_d2;

        Ra[0] = sqrt(R2[0] + a2);
        Ra[1] = sqrt(R2[1] + a2);

        Rainv0[0] = 1.0 / Ra[0];
        Rainv1[0] = 1.0 / Ra[1];

        Ra2inv[0] = Rainv0[0] * Rainv0[0];
        Ra2inv[1] = Rainv1[0] * Rainv1[0];

        for (i = 1; i < qMaxp1; i++) {
            Rainv0[i] = Rainv0[i-1] * Ra2inv[0];
            Rainv1[i] = Rainv1[i-1] * Ra2inv[1];
        }

        Kyztinfsum1[1][0] = y1 - y0;
        Jyzsuminf=log((Ra[0]+ypz[0])/((Ra[1]+ypz[1])));

        Rasum=Ra[0]-Ra[1];
        Rasuminf = Rasum - Kyztinfsum1[1][0];

        Kyztsum1[1][0] = z0*Jyzsuminf;
        Kyztsum2[2][0] = alpha*z0*Rasum+gammad*Kyztsum1[1][0];
        Kyztsum2[1][0] = delta*z0*Rasum+thetad*Kyztsum1[1][0];

        Lyztsum1[1][0] = y0*log(Ra[0]+ypz[0]) -
                         y1*log(Ra[1]+ypz[1]) +
                         log(a2_d2*0.5)*Kyztinfsum1[1][0];

        commonsum=(y0*(Ra[0]+y0)-y1*(Ra[1]+y1));
        Lyztsum2[1][0] = delta*commonsum+thetad*Lyztsum1[1][0];
        Lyztsum2[2][0] = alpha*commonsum+gammad*Lyztsum1[1][0];

        tempsuminf=(ymz[0]*Ra[0]-ymz[1]*Ra[1]+y0*y0-y1*y1)*a2d2inv;

        zJyzsuminf=Kyztsum1[1][0]-Kyztinfsum1[1][0];

        Hytinf1[1][1] = -0.5 * (Jyzsuminf-tempsuminf);
        Hztinf1[1][1] = -0.5 * (Jyzsuminf+tempsuminf);

        tmp1 = zJyzsuminf-Rasuminf;
        tmp2 = Lyztsum1[1][0]-Rasuminf;

        Htinf0[1][1] = Rasuminf*a2d2inv;

        Hytinf2[1][1] = delta*tmp1 + thetad*Hytinf1[1][1];
        Hytinf2[2][1] = alpha*tmp1 + gammad*Hytinf1[1][1];

        Hztinf2[1][1] = delta*tmp2 + thetad*Hztinf1[1][1];
        Hztinf2[2][1] = alpha*tmp2 + gammad*Hztinf1[1][1];

        Htinf1[1][1] = delta*(Hytinf1[1][1] + Hztinf1[1][1]) + thetad*Htinf0[1][1];
        Htinf2[1][1] = delta*(Hytinf2[1][1] + Hztinf2[1][1]) + thetad*Htinf1[1][1];
        Htinf1[2][1] = alpha*(Hytinf1[1][1] + Hztinf1[1][1]) + gammad*Htinf0[1][1];

        Htinf2[2][1] = alpha*(Hytinf2[1][1] + Hztinf2[1][1]) + gammad*Htinf1[1][1];
        Htinf2[3][1] = alpha*(Hytinf2[2][1] + Hztinf2[2][1]) + gammad*Htinf1[2][1];

        Hyt2_11 = Hytinf2[1][1];
        Hyt2_21 = Hytinf2[2][1];
        Ht1_11  = Htinf1[1][1];
        Ht1_21  = Htinf1[2][1];

        local = 0;
        local3 = 0;
        local4 = 0;


        /* H023 */
        reHy = creal(Hyt2_11);
        reHt = creal(Ht1_11);
        reHy1 = reHy - y0 * reHt;
        reHy2 = reHy - y1 * reHt;

        imHy = cimag(Hyt2_11);
        imHt = cimag(Ht1_11);
        imHy1 = imHy - y0 * imHt;
        imHy2 = imHy - y1 * imHt;


        for (i = 0; i < 18; i++) {
            ReFqval = FqReal[local3];
            G3[i] -= reHy2*ReFqval;
            G4[i] += reHy1*ReFqval;
            local3++;

#if 0
//
//          TEMP FIX BEFORE FULL REMOVAL OF COMPLEX VARIABLES
//
//          The new FqImag array does not contain the first 18 elements
//          because they were all zero, so this code is obsolete?
//
            ImFqval = FqImag[local4];
            G3[i] += imHy2*ImFqval;
            G4[i] -= imHy1*ImFqval;
            local4++;
#endif
        }

        /* H203 */
        local = 1;
        reHy = creal(Hyt2_21); 
        reHt = creal(Ht1_21);
        reHy1 = reHy - y0 * reHt;
        reHy2 = reHy - y1 * reHt;

        imHy = cimag(Hyt2_21); 
        imHt = cimag(Ht1_21);
        imHy1 = imHy - y0 * imHt;
        imHy2 = imHy - y1 * imHt;

        for (i = 0; i < 18; i++) {
            ReFqval = FqReal[local3];
            G3[i] -= reHy2*ReFqval;
            G4[i] += reHy1*ReFqval;

            ImFqval = FqImag[local4];
            G3[i] += imHy2*ImFqval;
            G4[i] -= imHy1*ImFqval;
            local3++;
            local4++;
        }

/*
 *      Recursive integrals
 */
        ind = 1;
        indm1 = 0;
        local = 2;
        local2= 2;

        for (i = 2; i < qMaxp2; i++) {
            real8 factor_i   = factor[i];
            real8 factor_im1 = factor[i-1];

            kk = 2*i-2;
            k = kk;
            onemfactor = 1.0 - kk * factor_im1;

            dfactorTmp = dfactor[i-1];
            afactorTmp = afactor[i-1];


            for (j = 0; j < kk; j++) {
                k--;

                tmp1 = ypz[0]*RdpRds0[j][k]*Rainv0[i-2];
                tmp2 = ypz[1]*RdpRds1[j][k]*Rainv1[i-2];
                tmp3 = adin[j][k+1]*indgd[j]+adin[j+1][k]*indtd[k];

                tmpsum = z0*(tmp1-tmp2);
        
                Kyztsum0[j+1][ind] = a2d2inv * (factor_im1 * (tmpsum +
                                     indgd[j]*Kyztsum1[j][indm1] + 
                                     indtd[k]*Kyztsum1[j+1][indm1]) +
                                     onemfactor*Kyztsum2[j+1][indm1]);


                tmpsum = y0*(tmp1+y0*adin[j+1][k+1]+tmp3) - 
                         y1*(tmp2+y1*adin[j+1][k+1]+tmp3);

                Lyztsum0[j+1][ind] = a2d2inv * (factor_im1 * (tmpsum +
                                     indgd[j]*Lyztsum1[j][indm1] + 
                                     indtd[k]*Lyztsum1[j+1][indm1]) +
                                     onemfactor*Lyztsum2[j+1][indm1]);
            }

            /* n = 0 */
            jj = 2*i-4;
            j = jj+1;

            Kyztinfsum1[1][ind]=dfactorTmp*(indd[j]*Kyztinfsum1[1][indm1]);


            tmp1=Rdsin0[j]*Rainv0[i-2];
            tmp2=Rdsin1[j]*Rainv1[i-2];

            tmpsum = z0*(tmp1-tmp2);

            Kyztsum1[1][ind] = dfactorTmp *
                               (indd[j] * Kyztsum1[1][indm1] - tmpsum) +
                               thetad * Kyztsum0[1][ind];

            tmpsum = y0*(tmp1+din[j+1])-y1*(tmp2+din[j+1]);

            Lyztsum1[1][ind] = dfactorTmp *
                               (indd[j] * Lyztsum1[1][indm1] - tmpsum) +
                               thetad * Lyztsum0[1][ind];

            kk = 2*i-2;
            k = kk;

            for (j = 0; j < kk; j++) {
                k--;

                Kyztinfsum1[j+2][ind] = afactor[i-1] *
                                        (indd[k]*Kyztinfsum1[j+1][indm1] + 
                                         inda[j]*Kyztinfsum1[j][indm1]);
         
                tmp1 = RdpRds0[j][k]*Rainv0[i-2];
                tmp2 = RdpRds1[j][k]*Rainv1[i-2];

                tmpsum = z0*(tmp1-tmp2);

                Kyztsum1[j+2][ind] = afactorTmp *
                                     (indd[k] * Kyztsum1[j+1][indm1] +
                                      inda[j] * Kyztsum1[j][indm1] - tmpsum) +
                                     gammad*Kyztsum0[j+1][ind];

                tmpsum = y0*(tmp1+adin[j+1][k+1])-y1*(tmp2+adin[j+1][k+1]);

                Lyztsum1[j+2][ind] = afactorTmp *
                                     (indd[k] * Lyztsum1[j+1][indm1] +
                                      inda[j] * Lyztsum1[j][indm1] - tmpsum) +
                                     gammad*Lyztsum0[j+1][ind];
            }

            /* n = 1 */
            j = jj+2;

            tmp1=Rdsin0[j]*Rainv0[i-2];
            tmp2=Rdsin1[j]*Rainv1[i-2];
            tmp3=din[j]*indtd[j];

            tmpsum=z0*(tmp1-Rdsin1[j]*Rainv1[i-2]);


            Kyztsum2[1][ind] = dfactorTmp *
                               (indd[j] * Kyztsum2[1][indm1] - tmpsum) +
                               thetad * Kyztsum1[1][ind];

            tmpsum = y0*(tmp1+din[j+1]*y0+tmp3) - 
                     y1*(tmp2+din[j+1]*y1+tmp3);


            Lyztsum2[1][ind] = dfactorTmp *
                               (indd[j] * Lyztsum2[1][indm1] - tmpsum) +
                               thetad * Lyztsum1[1][ind];

            kk = 2*i-1;
            k = kk;

            for (j = 0; j < kk; j++) {
                k--;

                tmp1=RdpRds0[j][k]*Rainv0[i-2];
                tmp2=RdpRds1[j][k]*Rainv1[i-2];
                tmp3=adin[j][k+1]*indgd[j]+adin[j+1][k]*indtd[k];

                tmpsum = z0*(tmp1-tmp2);

                Kyztsum2[j+2][ind] = afactorTmp *
                                     (indd[k] * Kyztsum2[j+1][indm1] +
                                      inda[j] * Kyztsum2[j][indm1] - tmpsum) +
                                     gammad*Kyztsum1[j+1][ind];

                tmpsum = y0*(tmp1+adin[j+1][k+1]*y0+tmp3) - 
                         y1*(tmp2+adin[j+1][k+1]*y1+tmp3);


                Lyztsum2[j+2][ind] = afactorTmp *
                                     (indd[k] * Lyztsum2[j+1][indm1] +
                                      inda[j] * Lyztsum2[j][indm1] - tmpsum) +
                                      gammad*Lyztsum1[j+1][ind];
            }  /* end for (j = 0; j < kk; j++) */


            indm1 = ind;
            ind = 1 - ind;

            onem2factor = 1.0 - (2*i-1) * factor_i;

            kk = 2*i-2;
            k = kk;

            for (j = 0; j < kk; j++) {
                complex8 tmpa, tmpb;

                k--;

                tmp1 = RdpRds0[j][k]*Rainv0[i-2];
                tmp2 = RdpRds1[j][k]*Rainv1[i-2];

                tmp3 = y0*(tmp1+adin[j+1][k+1]) -
                       y1*(tmp2+adin[j+1][k+1]);

                tmpinf1 = factor[i-1]*(inda[j]*Lyztsum1[j][ind] + 
                                       indd[k]*Lyztsum1[j+1][ind] - tmp3);

                tmp3 = z0*(tmp1-tmp2);
                tmpinf2 = factor[i-1] *
                          (inda[j]*(Kyztsum1[j][ind]   - Kyztinfsum1[j][ind]) + 
                           indd[k]*(Kyztsum1[j+1][ind] - Kyztinfsum1[j+1][ind]) -
                           tmp3);

                Hytinf0[j+1][ind] = a2d2inv *
                                    (onem2factor * Hytinf2[j+1][indm1] + factor_i *
                                     (tmpinf1 +
                                      indgd[j] * Hytinf1[j][indm1] +
                                      indtd[k] * Hytinf1[j+1][indm1] -
                                      Hztinf2[j+1][indm1]));

                Hztinf0[j+1][ind] = a2d2inv *
                                    (onem2factor * Hztinf2[j+1][indm1] + factor_i *
                                     (tmpinf2 + 
                                      indgd[j] * Hztinf1[j][indm1] +
                                      indtd[k] * Hztinf1[j+1][indm1] -
                                      Hytinf2[j+1][indm1]));

            }  /* end for (j = 0; j < kk; j++) */

            /* n = 0 */
            jj = 2*i-4;
            j = jj+1;

            Hytinf1[1][ind] = dfactor[i] *
                              (Htinf1[1][indm1] + indd[j]*Hytinf1[1][indm1] - 
                               Lyztsum0[1][indm1]) +  thetad*Hytinf0[1][ind];

            Hztinf1[1][ind] = dfactor[i] *
                              (Htinf1[1][indm1] + indd[j]*Hztinf1[1][indm1] - 
                               Kyztsum0[1][indm1] + Kyztinfsum0[1][indm1]) + 
                               thetad*Hztinf0[1][ind];

            kk = 2*i-1;

            for (j = 1; j < kk; j++) {
                k = kk - (j+1);

                Hytinf1[j+1][ind] = afactor[i] *
                                    (Htinf1[j][indm1] +
                                     inda[j-1] * Hytinf1[j-1][indm1]   + 
                                     indd[k]   * Hytinf1[j][indm1]     - 
                                     Lyztsum0[j][indm1]) +
                                    gammad*Hytinf0[j][ind];

                Hztinf1[j+1][ind] = afactor[i] *
                                    (Htinf1[j][indm1]  +
                                     inda[j-1] * Hztinf1[j-1][indm1] + 
                                     indd[k] * Hztinf1[j][indm1] -
                                     Kyztsum0[j][indm1] + Kyztinfsum0[j][indm1]) +
                                    gammad*Hztinf0[j][ind];            
            }

            /* n = 1 */
            j = jj+2;

            Hytinf2[1][ind] = dfactor[i] *
                              (Htinf2[1][indm1] + indd[j]*Hytinf2[1][indm1] -
                               Lyztsum1[1][indm1]) + thetad*Hytinf1[1][ind];

            Hztinf2[1][ind] = dfactor[i] *
                              (Htinf2[1][indm1] + indd[j]*Hztinf2[1][indm1] -
                               Kyztsum1[1][indm1] + Kyztinfsum1[1][indm1]) + 
                              thetad*Hztinf1[1][ind];

            kk = 2*i;

            for (j = 1; j < kk; j++) {
                k = kk - (j+1);

                Hytinf2[j+1][ind] = afactor[i] *
                                    (Htinf2[j][indm1] +
                                     inda[j-1]*Hytinf2[j-1][indm1] +
                                     indd[k]*Hytinf2[j][indm1] -
                                     Lyztsum1[j][indm1]) +
                                    gammad*Hytinf1[j][ind];

                Hztinf2[j+1][ind] = afactor[i] *
                                    (Htinf2[j][indm1] + 
                                     inda[j-1]*Hztinf2[j-1][indm1] +
                                     indd[k]*Hztinf2[j][indm1] - 
                                     Kyztsum1[j][indm1] + Kyztinfsum1[j][indm1]) +
                                    gammad*Hztinf1[j][ind];
            }
 
/*
 *          Construct H
 */
            onem2factor = 1.0 - factor_i * 2 * i;

            kk = 2*i;

            for (j = 1; j < kk; j++) {

                k = kk - (j+1);

                tmpsum = Kyztsum1[j][indm1] - Kyztinfsum1[j][indm1] + 
                         Lyztsum1[j][indm1];

                Htinf0[j][ind] = a2d2inv * (factor_i *  (tmpsum + 
                                 indgd[j-1] * Htinf1[j-1][indm1] + 
                                 indtd[k]   * Htinf1[j][indm1])  +
                                 onem2factor*Htinf2[j][indm1]); 
            }

            jj = 2*i-3;

            /* m=1 */
            Htinf1[1][ind] = delta  * (Hytinf1[1][ind] + Hztinf1[1][ind]) +
                             thetad *  Htinf0[1][ind];

            for (j = 1; j < jj+3; j++) {

                Htinf1[j+1][ind] = alpha  * (Hytinf1[j][ind] + Hztinf1[j][ind]) +
                                   gammad *  Htinf0[j][ind];
            }

            /* m=2 */
            Htinf2[1][ind] = delta  * (Hytinf2[1][ind] + Hztinf2[1][ind]) +
                             thetad *   Htinf1[1][ind];

            for (j = 1; j < jj+4; j++) {

                Htinf2[j+1][ind] = alpha  * (Hytinf2[j][ind] + Hztinf2[j][ind]) +
                                   gammad * Htinf1[j][ind];
            }

            for (jj = 1; jj < RemArray[i]; jj++) {

                j = RelArray[local];

                reHy = creal(Hytinf2[j][ind]);
                reHt = creal(Htinf1[j][ind]);

                reHy1 = reHy - y0 * reHt;
                reHy2 = reHy - y1 * reHt;

                for (iii = 0; iii < Re18Array[local]; iii++) {
                    ii = ReiiArray[local3];
                    ReFqval = FqReal[local3];
                    G3[ii] -= reHy2*ReFqval;
                    G4[ii] += reHy1*ReFqval;
                    local3++;
                }

                local++;  
            }

            for (jj = 1; jj < ImmArray[i]; jj++) {

                j = ImlArray[local2];

                imHy = cimag(Hytinf2[j][ind]);
                imHt = cimag(Htinf1[j][ind]);
                imHy1 = imHy - y0 * imHt;
                imHy2 = imHy - y1 * imHt;

                for (iii = 0; iii < Im18Array[local2]; iii++) {
                    ii = ImiiArray[local4];
                    ImFqval = FqImag[local4];
                    G3[ii] += imHy2*ImFqval;
                    G4[ii] -= imHy1*ImFqval;
                    local4++;
                }

                local2++;
            }

        }  /* end for (i = 2; i < qMaxp2; i++) */
 
/*
 *      Compute pre-factor eps C C b t bp in the form of
 *      matrices Ap, B, Bp, A
 *
 *      Note: segments are parallel, so <t> is used in calculations
 *      for <bptp*> variables as well as <bt*> variables.
 */
        bptp11 = bpx*t[0];
        bptp22 = bpy*t[1];
        bptp33 = bpz*t[2];

        bptp12 = bpx*t[1];
        bptp23 = bpy*t[2];
        bptp31 = bpz*t[0];

        bptp21 = bpy*t[0];
        bptp32 = bpz*t[1];
        bptp13 = bpx*t[2];

        Ap12 = bptp13;
        Ap62 = bptp23;
        Ap51 = bptp12;
        Ap23 = bptp21;
        Ap43 = bptp31;
        Ap31 = bptp32;

        Ap63 = bptp11 - bptp22;
        Ap41 = bptp22 - bptp33;
        Ap52 = bptp33 - bptp11;

        bt12 = bx*t[1];
        bt23 = by*t[2];
        bt31 = bz*t[0];

        bt21 = by*t[0];
        bt32 = bz*t[1];
        bt13 = bx*t[2];

        bt11 = bx*t[0];
        bt22 = by*t[1];
        bt33 = bz*t[2];

        B12 = bt23;
        B35 = bt32;
        B31 = bt12;
        B16 = bt13;
        B23 = bt31;
        B24 = bt21;

        B36 = bt22 - bt11;
        B14 = bt33 - bt22;
        B25 = bt11 - bt33;

        switch(home->param->materialType) {
            case MAT_TYPE_BCC:
            case MAT_TYPE_FCC:
            case MAT_TYPE_HCP:
/*
 *              For hexagonal or cubic elastic constant matrix
 */
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

                for (j = 0; j < 3; j++) {
                    A[j][0] = -ecMatrix[j][1]*B12 + ecMatrix[j][2]*B35;
                    A[j][1] =  ecMatrix[j][0]*B16 - ecMatrix[j][2]*B23;
                    A[j][2] = -ecMatrix[j][0]*B31 + ecMatrix[j][1]*B24;
                }

                A[3][0] = -ecMatrix[3][3]*B14;
                A[3][1] = -ecMatrix[3][3]*B24;
                A[3][2] =  ecMatrix[3][3]*B23;

                A[4][0] =  ecMatrix[4][4]*B31;
                A[4][1] = -ecMatrix[4][4]*B25;
                A[4][2] = -ecMatrix[4][4]*B35;

                A[5][0] = -ecMatrix[5][5]*B16;
                A[5][1] =  ecMatrix[5][5]*B12;
                A[5][2] = -ecMatrix[5][5]*B36;


                for (j = 0; j < 6; j++) {
                    B[0][j]  = -A[j][0];
                    B[1][j]  = -A[j][1];
                    B[2][j]  = -A[j][2];
                }

                break;

            case MAT_TYPE_RHOMBOHEDRAL_VA:
/*
 *              For general elastic constant matrix
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

                    B[0][j] =  ecMatrix[1][j]*B12 - ecMatrix[2][j]*B35 +
                               ecMatrix[3][j]*B14 - ecMatrix[4][j]*B31 +
                               ecMatrix[5][j]*B16;
                    B[1][j] = -ecMatrix[0][j]*B16 + ecMatrix[2][j]*B23 +
                               ecMatrix[3][j]*B24 + ecMatrix[4][j]*B25 -
                               ecMatrix[5][j]*B12;
                    B[2][j] =  ecMatrix[0][j]*B31 - ecMatrix[1][j]*B24 -
                               ecMatrix[3][j]*B23 + ecMatrix[4][j]*B35 +
                               ecMatrix[5][j]*B36;
                }

                break;
#if 0 
            case XXX:
/*
 *              For symetric elastic constant matrix
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

                    A[j][0] = -ecMatrix[j][1]*B12 + ecMatrix[j][2]*B35 -
                               ecMatrix[j][3]*B14 + ecMatrix[j][4]*B31 -
                               ecMatrix[j][5]*B16;

                    A[j][1] =  ecMatrix[j][0]*B16 - ecMatrix[j][2]*B23 -
                               ecMatrix[j][3]*B24 - ecMatrix[j][4]*B25 +
                               ecMatrix[j][5]*B12;

                    A[j][2] = -ecMatrix[j][0]*B31 + ecMatrix[j][1]*B24 +
                               ecMatrix[j][3]*B23 - ecMatrix[j][4]*B35 -
                               ecMatrix[j][5]*B36;
                }

                for (j = 0; j < 6; j++) {
                    B[0][j]  = -A[j][0];
                    B[1][j]  = -A[j][1];
                    B[2][j]  = -A[j][2];
                }

                break;
#endif  /* end for symetric C */

        }  /* end switch(home->param->materialType) */

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                ApGaux3[i][j] = 0.0;
                ApGaux4[i][j] = 0.0;

                for (k = 0; k < 6; k++) {

                    offsetG = 6*j+k;
                    ApGaux3[i][j] = ApGaux3[i][j] + Ap[k][i]*G3[offsetG];
                    ApGaux4[i][j] = ApGaux4[i][j] + Ap[k][i]*G4[offsetG];
                }
            }
        }

        ApG3[0] = ApGaux3[0][0];
        ApG3[1] = ApGaux3[1][1];
        ApG3[2] = ApGaux3[2][2];

        ApG3[3] = (ApGaux3[1][2] + ApGaux3[2][1]);
        ApG3[4] = (ApGaux3[2][0] + ApGaux3[0][2]);
        ApG3[5] = (ApGaux3[0][1] + ApGaux3[1][0]);

        ApG4[0] = ApGaux4[0][0];
        ApG4[1] = ApGaux4[1][1];
        ApG4[2] = ApGaux4[2][2];

        ApG4[3] = (ApGaux4[1][2] + ApGaux4[2][1]);
        ApG4[4] = (ApGaux4[2][0] + ApGaux4[0][2]);
        ApG4[5] = (ApGaux4[0][1] + ApGaux4[1][0]);

        *fp3x = 0.0;
        *fp3y = 0.0;
        *fp3z = 0.0;

        *fp4x = 0.0;
        *fp4y = 0.0;
        *fp4z = 0.0;

        for (j = 0; j < 6; j++) {

            *fp3x += B[0][j] * ApG3[j];
            *fp3y += B[1][j] * ApG3[j];
            *fp3z += B[2][j] * ApG3[j];

            *fp4x += B[0][j] * ApG4[j];
            *fp4y += B[1][j] * ApG4[j];
            *fp4z += B[2][j] * ApG4[j];
        }

        *fp3x *= coeff34;
        *fp3y *= coeff34;
        *fp3z *= coeff34;

        *fp4x *= coeff34;
        *fp4y *= coeff34;
        *fp4z *= coeff34;

        return;
}
#endif
#endif  // end if defined ANISOTROPIC && _ARLFEM
