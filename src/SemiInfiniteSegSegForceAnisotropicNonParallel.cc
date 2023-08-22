
#include "mpi_portability.h"

#include "Home.h"

#ifdef ANISO_QMAX
#define QMAX ANISO_QMAX
#else
#define QMAX 20
#endif

#if defined ANISOTROPIC && defined _ARLFEM
/***************************************************************************
 *
 *      Function:    SemiInfiniteSegSegForceAnisotropicNonParallel()
 *
 *      Description: Calculates the anisotropic force on segment 2 from
 *                   the semi-inf non-parallel segment 1.
 *
 *      Parameters:
 *          IN:   p1{xyz}       Position of first segment's first endpoint
 *          IN:   p2{xyz}       Point on the semi-infinite dislocation
 *                              with an origin at point p1
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
#if 0  // Anisotropic forces need to be updated using table for version 3.


void SemiInfiniteSegSegForceAnisotropicNonParallel(Home_t *home,
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
        int qMaxp1 = qMax+1;
        int qMaxp2 = qMax+2;
        int twoqMaxp2 = 2*qMax+2; 
        int twoqMaxp3 = 2*qMax+3;
        int twoqMaxp4 = 2*qMax+4;

        int   i, j, k, m, n, ii, jj, kk, iii;
        int   im1,twoi;
        int   local, local2, local3, local4;
        int   ilist, ind, ind1, indm1;

        short   *RemArray, *ImmArray;
        short   *RelArray, *ImlArray;
        short   *Re18Array, *Im18Array;
        short   *ReiiArray, *ImiiArray;

        real8 (*ecMatrix)[6];
        real8 *FqReal, *FqImag, *rotMatrixRow3;
        real8 a2, a2_d2, a2d2inv;
        real8 adW_003inf;
        real8 Ap12, Ap23, Ap31, Ap41, Ap43, Ap51, Ap52, Ap62, Ap63;
        real8 B12, B14, B16, B23, B24, B25, B31, B35, B36;
        real8 bt11,bt22,bt33;
        real8 bptp11,bptp22,bptp33;
        real8 c, c2;
        real8 cdenom;
        real8 coeff12, coeff34;
        real8 commonW223, commonW223inf;
        real8 d;
        real8 denom;
        real8 dX, dY, dZ;
        real8 delta, epsilon;
        real8 dy[2];
        real8 dmce, emcd;
        real8 invL;
        real8 Rasum, Rasuminf;
        real8 thetad;
        real8 onemc2, onemc2inv;
        real8 onem2factor;
        real8 tp[3], t[3];
        real8 A[6][3], Ap[6][3];
        real8 B[3][6];
        real8 G3[18], G4[18];
        real8 ApGaux3[3][3], ApGaux4[3][3];
        real8 ApG3[6], ApG4[6];
        real8 indd[twoqMaxp3];
        real8 inde[twoqMaxp3];
        real8 indtd[twoqMaxp3];

        real8 indemcdzptd[twoqMaxp3];
        real8 inddmceyptd[2][twoqMaxp3];

        real8 inddy[2][twoqMaxp3];
        real8 inddyptd[2][twoqMaxp3];
        real8 dyptd[2];

        real8 kdmce[twoqMaxp2];
        real8 kemcd[twoqMaxp2];
        real8 Rdt[2], Rdtp[2];
        real8 R2[2], Ra[2];
        real8 Rainv0[qMaxp1];
        real8 Rainv1[qMaxp1];

        real8 Ra2inv[2];
        real8 RdtRainv[qMaxp1][2];
        real8 RdtpRainv[qMaxp1][2];
        real8 temp1, temp2, tmpVeca[2], tmpVecb[2], tmpVecinf[2];
        real8 factorinv[qMaxp2];
        real8 factor[qMaxp2], factor1[2], factor2[2];
        real8 newfactor[qMaxp2];
        real8 R[2][3];
        real8 temp,tempbis;
        real8 temp3a, temp3b, temp3c, temp3d;
        real8 y0,y1,z0,z1;

        real8 yin[2], zin[2];

        real8 dmceyptd[2], emcdzptd;
        real8 txtp[3];
        real8 Rdsin0[twoqMaxp3];
        real8 Rdsin1[twoqMaxp3];

        real8 ein[twoqMaxp4];

        real8 r1, r2;
        real8 i1, i2;
        real8 tmpre0, tmpre1, tmpreinf0, tmpreinf1;
        real8 onemfactor;
        real8 dfactor, efactor;

        real8 reHy1, reHy2, reHz1, reHz2;
        real8 imHy1, imHy2, imHz1, imHz2;
        real8 ReFqval, ImFqval;
        real8 reHzt, reHyt, reHt;
        real8 imHzt, imHyt, imHt;

        complex8 *rotMatrix12;
        complex8 tempt;
        complex8 afactor, bfactor;
        complex8 tmpcpx0, tmpcpx1;
        complex8 tmpcpxinf0, tmpcpxinf1;
        complex8 alpha, beta;
        complex8 amcb, bmca;
        complex8 gammad;
        complex8 aypgd[2];
        complex8 I013inf;
        complex8 I023inf;
        complex8 I103inf;
        complex8 I113inf;
        complex8 I203inf;
        complex8 Jytsuminf, Jztsuminf;

        complex8 amcbypgd[2], bmcazpgd;

        complex8 inda[twoqMaxp3];
        complex8 indb[twoqMaxp3];
        complex8 indgd[twoqMaxp3];

        complex8 ay[2], inday[2][twoqMaxp3];

        complex8 indbmcazpgd[twoqMaxp3];
        complex8 indamcbypgd[2][twoqMaxp3];
        complex8 indaypgd[2][twoqMaxp3];

        complex8 jamcb[twoqMaxp2];
        complex8 jbmca[twoqMaxp2];
        complex8 Rdpin0[twoqMaxp3];
        complex8 Rdpin1[twoqMaxp3];

        complex8 bin[twoqMaxp4];
        complex8 bein[twoqMaxp4][twoqMaxp4];

        complex8 RdpRds0[twoqMaxp3][twoqMaxp3];
        complex8 RdpRds1[twoqMaxp3][twoqMaxp3];

        complex8 Htinf0[2][twoqMaxp4], Htinf1[2][twoqMaxp4];
        complex8 Htinf2[2][twoqMaxp4];

        complex8 Hytinf0[twoqMaxp3],Hytinf1[twoqMaxp3];
        complex8 Hztinf0[twoqMaxp3],Hztinf1[twoqMaxp3];

        complex8 Jyt00[2][twoqMaxp3],Jyt10[2][twoqMaxp3],Jyt20[2][twoqMaxp3];
        complex8 Jyt01[2][twoqMaxp3],Jyt11[2][twoqMaxp3],Jyt21[2][twoqMaxp3];
        complex8 Jzt00[2][twoqMaxp3],Jzt10[2][twoqMaxp3],Jzt20[2][twoqMaxp3];
        complex8 Jzt01[2][twoqMaxp3],Jzt11[2][twoqMaxp3],Jzt21[2][twoqMaxp3];

        complex8 Jytinf00[2][twoqMaxp3], Jytinf10[2][twoqMaxp3];
        complex8 Jytinf20[2][twoqMaxp3];
        complex8 Jytinf01[2][twoqMaxp3], Jytinf11[2][twoqMaxp3];
        complex8 Jytinf21[2][twoqMaxp3];
        complex8 Jztinf00[2][twoqMaxp3], Jztinf10[2][twoqMaxp3];
        complex8 Jztinf20[2][twoqMaxp3];
        complex8 Jztinf01[2][twoqMaxp3], Jztinf11[2][twoqMaxp3];
        complex8 Jztinf21[2][twoqMaxp3];

        complex8 Lytinf0[2][twoqMaxp3],Lytinf1[2][twoqMaxp3];
        complex8 Kytinf0[2][twoqMaxp3],Kytinf1[2][twoqMaxp3];

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

        AnisotropicVars_t *avars = &home->anisoVars;

        real8    *fqr    = 0;                                 
        real8    *fqi    = 0;                                 
        fqr    = (real8 *) avars->FqReal_v3;
        fqi    = (real8 *) avars->FqImag_v3;
        real8   (*FqReal)[(qMax+2)*(qMax+1)] = (real8 (*)[(qMax+2)*(qMax+1)]) fqr;
        real8   (*FqImag)[(qMax+2)*(qMax+1)] = (real8 (*)[(qMax+2)*(qMax+1)]) fqi;

        rotMatrix12 = avars->anisoRotMatrix12;
        rotMatrixRow3  = &(avars->anisoRotMatrix[2][0]);

        //c66    = (real8 *) avars->elasticConstantMatrix2D;
        ecMatrix = home->param->elasticConstantMatrix;


/*
 *      Initialization
 */
        for (i = 0; i < 18; i++) {
            G3[i] = 0.0;
            G4[i] = 0.0;
        }

/*
 *      Initialize all large arrays
 */
        for (i = 0; i < 2; i++) {
            for (j = 0; j < twoqMaxp4; j++) {
                Htinf0[i][j] = 0.0;
                Htinf1[i][j] = 0.0;
                Htinf2[i][j] = 0.0;
            }
        }

        for (i = 0; i < twoqMaxp3; i++) {
            Hytinf0[i] = 0.0; Hytinf1[i] = 0.0;
            Hztinf0[i] = 0.0; Hztinf1[i] = 0.0;
        }

        for (i = 0; i < 2; i++) {
            for (j = 0; j < twoqMaxp3; j++) {
                Jyt00[i][j] = 0.0; Jyt10[i][j] = 0.0; Jyt20[i][j] = 0.0;
                Jyt01[i][j] = 0.0; Jyt11[i][j] = 0.0; Jyt21[i][j] = 0.0;

                Jzt00[i][j] = 0.0; Jzt10[i][j] = 0.0; Jzt20[i][j] = 0.0;
                Jzt01[i][j] = 0.0; Jzt11[i][j] = 0.0; Jzt21[i][j] = 0.0;

                Jytinf00[i][j] = 0.0; Jytinf10[i][j] = 0.0; Jytinf20[i][j] = 0.0;
                Jytinf01[i][j] = 0.0; Jytinf11[i][j] = 0.0; Jytinf21[i][j] = 0.0;

                Jztinf00[i][j] = 0.0; Jztinf10[i][j] = 0.0; Jztinf20[i][j] = 0.0;
                Jztinf01[i][j] = 0.0; Jztinf11[i][j] = 0.0; Jztinf21[i][j] = 0.0;

                Lytinf0[i][j] = 0.0; Lytinf1[i][j] = 0.0;
                Kytinf0[i][j] = 0.0; Kytinf1[i][j] = 0.0;
            }
        }

        for (i = 0; i < twoqMaxp3; i++) {
            for (j = 0; j < twoqMaxp3; j++) {
                RdpRds0[i][j] = 0.0;
                RdpRds1[i][j] = 0.0;

                bein[i+1][j+1] = 0.0;
            }
        }

/* 
 *      Define line directions
 */
        tp[X] = p2x - p1x;
        tp[Y] = p2y - p1y;
        tp[Z] = p2z - p1z;

        NormalizeVec(tp);

        dX = p4x - p3x;
        dY = p4y - p3y;
        dZ = p4z - p3z;

        temp = dX * dX + dY * dY + dZ * dZ;
        invL = 1.0 / sqrt(temp);

        coeff34 = invL * INV_4_PI_SQUARED;

        t[X] = dX * invL;
        t[Y] = dY * invL;
        t[Z] = dZ * invL;

/*
 *      Compute double and single integrals via recurrence
 */
        c = DotProduct(tp, t);
        c2 = c * c;
        onemc2 = 1.0 - c2;
        onemc2inv = 1.0/onemc2;

        txtp[0] = t[1] * tp[2] - t[2] * tp[1];
        txtp[1] = t[2] * tp[0] - t[0] * tp[2];
        txtp[2] = t[0] * tp[1] - t[1] * tp[0];

        alpha = DotProduct(t,  rotMatrix12);
        beta  = DotProduct(tp, rotMatrix12);

        delta   = DotProduct(t,  rotMatrixRow3);
        epsilon = DotProduct(tp, rotMatrixRow3);

        amcb = alpha - (c * beta);
        bmca = beta - (c * alpha);

        dmce = delta - (c * epsilon);
        emcd = epsilon - (c * delta);

        R[0][0] = p3x - p1x;
        R[0][1] = p3y - p1y;
        R[0][2] = p3z - p1z;

        R[1][0] = p4x - p1x;
        R[1][1] = p4y - p1y;
        R[1][2] = p4z - p1z;


        d = DotProduct(R[0], txtp) * onemc2inv;

        tmpVeca[0] = DotProduct(R[0],t);
        tmpVeca[1] = DotProduct(R[1],t);

        tmpVecb[0] = DotProduct(R[0],tp);
        tmpVecb[1] = DotProduct(R[1],tp);

        y0 = (tmpVeca[0] - c * tmpVecb[0]) * onemc2inv;
        y1 = (tmpVeca[1] - c * tmpVecb[1]) * onemc2inv;

        z0 = (tmpVecb[0] - c * tmpVeca[0]) * onemc2inv;

        yin[0] = y0;
        yin[1] = y1;

        Rdpin0[0] = 1.0;
        Rdsin0[0] = 1.0;
        Rdpin0[1] = DotProduct(R[0], rotMatrix12);
        Rdsin0[1] = DotProduct(R[0], rotMatrixRow3);

        Rdpin1[0] = 1.0;
        Rdsin1[0] = 1.0;
        Rdpin1[1] = DotProduct(R[1], rotMatrix12);
        Rdsin1[1] = DotProduct(R[1], rotMatrixRow3);

        gammad = 0.5 * (Rdpin0[1] + Rdpin1[1] - (y0+y1) * alpha) -
                        beta * z0;
        thetad = 0.5 * (Rdsin0[1] + Rdsin1[1] - (y0+y1) * delta) -
                        epsilon * z0;

        bmcazpgd = (bmca * z0 + gammad);
        emcdzptd = (emcd * z0 + thetad);

        amcbypgd[0] = (amcb * y0 + gammad);
        dmceyptd[0] = (dmce * y0 + thetad);

        amcbypgd[1] = (amcb * y1 + gammad);
        dmceyptd[1] = (dmce * y1 + thetad);

        dy[0]=delta * y0;
        dy[1]=delta * y1;

        ay[0]=alpha * y0;
        ay[1]=alpha * y1;

        dyptd[0]=(dy[0]+thetad);
        dyptd[1]=(dy[1]+thetad);

        aypgd[0]=(ay[0]+gammad);
        aypgd[1]=(ay[1]+gammad);

        bin[0]=0.0e0;
        ein[0]=0.0e0;
        bin[1]=1;
        ein[1]=1;
        bin[2]=beta;
        ein[2]=epsilon;


        for (i = 2; i < twoqMaxp3; i++) {

             Rdpin0[i] = Rdpin0[i-1] * Rdpin0[1];
             Rdpin1[i] = Rdpin1[i-1] * Rdpin1[1];

             Rdsin0[i] = Rdsin0[i-1] * Rdsin0[1];
             Rdsin1[i] = Rdsin1[i-1] * Rdsin1[1];

             bin[i+1] = bin[i] * bin[2];
             ein[i+1] = ein[i] * ein[2];
        }

        for (i = 0; i < twoqMaxp3; i++) {
            for (j = 0; j < twoqMaxp3-i; j++) {
                RdpRds0[i][j] = Rdpin0[i] * Rdsin0[j];
                RdpRds1[i][j] = Rdpin1[i] * Rdsin1[j];

                bein[i+1][j+1] = bin[i+1] * ein[j+1];
            }
        }

        a2 = a * a;
        a2_d2 = a2 + d * d * onemc2;

        Rdt[0] = y0 + z0 * c;
        Rdt[1] = y1 + z0 * c;

        Rdtp[0] = z0 + y0 * c;
        Rdtp[1] = z0 + y1 * c;

        cdenom = onemc2  *  a2_d2;

        for (i = 0; i < 2; i++) {
            R2[i] = DotProduct(R[i], R[i]);

            temp1=R2[i] + a2;

            tmpVeca[i] = temp1 - Rdt[i] * Rdt[i];
            tmpVecb[i] = temp1 - Rdtp[i] * Rdtp[i];

            Ra[i] = sqrt(temp1);
        }

        temp = sqrt(cdenom);     

        factorinv[0]=-1;
        for (i = 1; i < qMaxp2; i++) {
            factorinv[i] = factorinv[i-1]+2;
        }

        Rainv0[0] = 1.0 / Ra[0];
        Rainv1[0] = 1.0 / Ra[1];

        a2d2inv = 1.0 / a2_d2;
        denom = 1.0 /temp;

        for (i = 0; i < 2; i++) {
            factor1[i] = 1.0 / tmpVeca[i];
            factor2[i] = 1.0 / tmpVecb[i];
        }

        for (i = 0; i < qMaxp2; i++) {
            factor[i]  = 1.0 / factorinv[i];
        }

        for (i = 0; i < qMaxp2; i++) {
            newfactor[i] = factor[i] * onemc2inv;
        }

        Ra2inv[0] = Rainv0[0] * Rainv0[0];
        Ra2inv[1] = Rainv1[0] * Rainv1[0];

        RdtRainv[0][0]  = Rainv0[0] * Rdt[0];
        RdtpRainv[0][0] = Rainv0[0] * Rdtp[0];

        RdtRainv[0][1]  = Rainv1[0] * Rdt[1];
        RdtpRainv[0][1] = Rainv1[0] * Rdtp[1];

        cdenom = (1.0+c) * denom;

        for (i = 1; i < qMaxp1; i++) {
            ii = i-1;
            Rainv0[i] = Rainv0[ii] * Ra2inv[0];
            RdtRainv[i][0] = Rainv0[i] * Rdt[0];
            RdtpRainv[i][0] = Rainv0[i] * Rdtp[0];

            Rainv1[i] = Rainv1[ii] * Ra2inv[1];
            RdtRainv[i][1] = Rainv1[i] * Rdt[1];
            RdtpRainv[i][1] = Rainv1[i] * Rdtp[1];
        }

        inda[0]  =  0.0e0;
        indb[0]  =  0.0e0;
        indd[0]  =  0.0e0;
        inde[0]  =  0.0e0;
        indgd[0] =  0.0e0;
        indtd[0] =  0.0e0;

        for (m = 1; m < twoqMaxp3; m++) {
            inda[m]  = inda[m-1] + alpha;
            indb[m]  = indb[m-1] + beta;
            indd[m]  = indd[m-1] + delta;
            inde[m]  = inde[m-1] + epsilon;
            indgd[m] = indgd[m-1] + gammad;
            indtd[m] = indtd[m-1] + thetad;
        }

        indemcdzptd[0] = 0.0e0;
        indbmcazpgd[0] = 0.0e0;

        for (i = 0; i < 2; i++) {
            indamcbypgd[i][0] = 0.0e0;
            inddmceyptd[i][0] = 0.0e0;

            inday[i][0] = 0.0e0;
            inddy[i][0] = 0.0e0;
            inddyptd[i][0] = 0.0e0;
            indaypgd[i][0] = 0.0e0;
        }

        for (m = 1; m < twoqMaxp3; m++) {

            indemcdzptd[m] = indemcdzptd[m-1]+emcdzptd;
            indbmcazpgd[m] = indbmcazpgd[m-1]+bmcazpgd;

            for (i = 0; i < 2; i++) {
                indamcbypgd[i][m] = indamcbypgd[i][m-1]+amcbypgd[i];
                inddmceyptd[i][m] = inddmceyptd[i][m-1]+dmceyptd[i];

                inday[i][m]    = inday[i][m-1]    + ay[i];
                inddy[i][m]    = inddy[i][m-1]    + dy[i];
                inddyptd[i][m] = inddyptd[i][m-1] + dyptd[i];
                indaypgd[i][m] = indaypgd[i][m-1] + aypgd[i];
            }
        }

        jamcb[0] =  0.0e0;
        jbmca[0] =  0.0e0;
        kdmce[0] =  0.0e0;
        kemcd[0] =  0.0e0;

        for (i = 1; i < twoqMaxp2; i++) {
            ii = i-1;
            jamcb[i] = jamcb[ii] + amcb;
            jbmca[i] = jbmca[ii] + bmca;
            kdmce[i] = kdmce[ii] + dmce;
            kemcd[i] = kemcd[ii] + emcd;
        }



        Jyt10[0][1] = log(Ra[0]+Rdt[0]);
        Jzt10[0][1] = log(Ra[0]+Rdtp[0]);

        Jyt11[0][1] = log(Ra[1]+Rdt[1]);
        Jzt11[0][1] = log(Ra[1]+Rdtp[1]);

        Jztinf10[0][1]=log((1-c2) * y0 * y0+a2_d2);
        Kytinf0[0][1]=-y0;

        Jztinf11[0][1]=log((1-c2) * y1 * y1+a2_d2);
        Kytinf1[0][1]=-y1;

        /* m = 0 */
        Jyt20[0][2] =   alpha * Ra[0] + bmcazpgd    * Jyt10[0][1];
        Jyt20[0][1] =   delta * Ra[0] + emcdzptd    * Jyt10[0][1];
        Jzt20[0][2] =    beta * Ra[0] + amcbypgd[0] * Jzt10[0][1];
        Jzt20[0][1] = epsilon * Ra[0] + dmceyptd[0] * Jzt10[0][1];

        Jytinf20[0][2] = beta    * Kytinf0[0][1];
        Jytinf20[0][1] = epsilon * Kytinf0[0][1];
        Jztinf20[0][1] =-epsilon * (c * y0) + dmceyptd[0] * Jztinf10[0][1];
        Jztinf20[0][2] =-beta    * (c * y0) + amcbypgd[0] * Jztinf10[0][1];

        /* m = 1 */
        Jyt21[0][2] =   alpha * Ra[1] + bmcazpgd * Jyt11[0][1];
        Jyt21[0][1] =   delta * Ra[1] + emcdzptd * Jyt11[0][1];
        Jzt21[0][2] =    beta * Ra[1] + amcbypgd[1] * Jzt11[0][1];
        Jzt21[0][1] = epsilon * Ra[1] + dmceyptd[1] * Jzt11[0][1];

        Jytinf21[0][2] = beta    * Kytinf1[0][1];
        Jytinf21[0][1] = epsilon * Kytinf1[0][1];
        Jztinf21[0][1] =-epsilon * (c * y1) + dmceyptd[1] * Jztinf11[0][1];
        Jztinf21[0][2] =-beta    * (c * y1) + amcbypgd[1] * Jztinf11[0][1];

        Jytsuminf =   Jyt10[0][1] -    Jyt11[0][1] - 
                   Jytinf10[0][1] + Jytinf11[0][1];


        Jztsuminf =   Jzt10[0][1] -    Jzt11[0][1] - 
                   Jztinf10[0][1] + Jztinf11[0][1];

/*
 *      Initialize Hysumint, Hzsumint, Hsumint
 *
 *      Preliminary double integrals involving unique terms
 */
        tmpVeca[0]=(Ra[0]+y0+z0) * cdenom;
        tmpVeca[1]=(Ra[1]+y1+z0) * cdenom;

        tmpVecinf[0]=(1-c) * y0 * cdenom;
        tmpVecinf[1]=(1-c) * y1 * cdenom;

        Htinf0[1][1] = atan(tmpVeca[1]) - atan(tmpVeca[0]) -
                       atan(tmpVecinf[1]) + atan(tmpVecinf[0]);

        Htinf0[1][1] = Htinf0[1][1] * 2.0 * denom;

        Rasuminf = Ra[0] - Ra[1] + c * (y0-y1);
        adW_003inf=a2_d2 * creal(Htinf0[1][1]);

        commonW223inf=(c * Rasuminf - adW_003inf ) * onemc2inv;

        I113inf = (c * adW_003inf - Rasuminf) * onemc2inv;
        I203inf = z0 * (Jyt10[0][1]-Jyt11[0][1]) -
                  Kytinf0[0][1] + Kytinf1[0][1] + commonW223inf;
        I023inf = y0 * (Jzt10[0][1]-Jztinf10[0][1]) -
                  y1 * (Jzt11[0][1]-Jztinf11[0][1]) + commonW223inf;

        I103inf = ( c * Jytsuminf  - Jztsuminf ) * onemc2inv;
        I013inf = ( c * Jztsuminf  - Jytsuminf ) * onemc2inv;

        Hytinf0[1] = I103inf;
        Hytinf1[1] = delta * I203inf + epsilon * I113inf+thetad * Hytinf0[1];
        Hytinf1[2] = alpha * I203inf    + beta * I113inf+gammad * Hytinf0[1];


        Hztinf0[1] = I013inf;
        Hztinf1[1] = delta * I113inf + epsilon * I023inf+thetad * Hztinf0[1];
        Hztinf1[2] = alpha * I113inf    + beta * I023inf+gammad * Hztinf0[1];

        Htinf1[1][1] = delta   * I103inf +
                       epsilon * I013inf +
                       thetad  * Htinf0[1][1];

        Htinf1[1][2] = alpha   * I103inf +
                       beta    * I013inf +
                       gammad  * Htinf0[1][1];

        Htinf2[1][1] = delta   * Hytinf1[1] +
                       epsilon * Hztinf1[1] +
                       thetad  * Htinf1[1][1];

        Htinf2[1][2] = alpha   * Hytinf1[1] +
                       beta    * Hztinf1[1] +
                       gammad  * Htinf1[1][1];

        Htinf2[1][3] = alpha   * Hytinf1[2] +
                       beta    * Hztinf1[2] +
                       gammad  * Htinf1[1][2];

        local3 = 0;
        local4 = 0;

        reHyt = creal(Hytinf1[1]);
        reHt  = creal(Htinf1[1][1]); 

        imHyt = cimag(Hytinf1[1]);
        imHt  = cimag(Htinf1[1][1]); 

        reHy1 = reHyt - y0 * reHt; 
        reHy2 = reHyt - y1 * reHt; 

        imHy1 = imHyt - y0 * imHt; 
        imHy2 = imHyt - y1 * imHt; 

        local = 0;

        for (i = 0; i < 18; i++) {
            ReFqval = FqReal[local3];

            G3[i] -= reHy2 * ReFqval;
            G4[i] += reHy1 * ReFqval;
            local3++;
#if 0
//
//          TEMP FIX BEFORE FULL REMOVAL OF COMPLEX VARIABLES
//
//          The new FqImag array does not contain the first 18 elements
//          because they were all zero, so this code is obsolete?
//
            ImFqval = FqImag[local4];

            G3[i] += imHy2 * ImFqval;
            G4[i] -= imHy1 * ImFqval;
            local4++;
#endif
        }

        reHyt = creal(Hytinf1[2]);
        reHt  = creal(Htinf1[1][2]); 

        imHyt = cimag(Hytinf1[2]);
        imHt  = cimag(Htinf1[1][2]); 

        reHy1 = reHyt - y0 * reHt; 
        reHy2 = reHyt - y1 * reHt; 

        imHy1 = imHyt - y0 * imHt; 
        imHy2 = imHyt - y1 * imHt; 

        local = 1;

        for (i = 0; i < 18; i++) {
            ReFqval = FqReal[local3];
            G3[i] -= reHy2 * ReFqval;
            G4[i] += reHy1 * ReFqval;

            ImFqval = FqImag[local4];
            G3[i] += imHy2 * ImFqval;
            G4[i] -= imHy1 * ImFqval;
            local3++;
            local4++;
        }

/*
 *      Start recursion
 */
        ind = 1;
        indm1 = 0;
        local = 2;
        local2= 2;

        for (i = 2; i < qMaxp2; i++) {
            complex8 tmp1, tmp2;

/*
 *          Calculate J00k
 */
            im1 = i-1;
            twoi = 2 * i;
    
            kk = twoi-2;
            onemfactor = 1.0 - kk * factor[im1];

            afactor = factor[im1] * alpha;
            bfactor = factor[im1] * beta;
            dfactor = factor[im1] * delta;
            efactor = factor[im1] * epsilon;

/*
 *          Calculate Ji00k J0ik
 */
            k = kk;

            tmpVeca[0]=(2-kk) * c * y0; 
            tmpVeca[1]=(2-kk) * c * y1; 

            tmpVecb[0]=((kk-1) * c2-1) * y0;
            tmpVecb[1]=((kk-1) * c2-1) * y1;

            for (j = 0; j < kk; j++ ) {

                k--;

                /* m = 0 */
                Jyt00[ind][j+1] = factor1[0] * 
                        (factor[i-1] * (RdpRds0[j][k] * RdtRainv[i-2][0]    + 
                                        indbmcazpgd[j] * Jyt10[indm1][j]    +
                                        indemcdzptd[k] * Jyt10[indm1][j+1]) + 
                         onemfactor * Jyt20[indm1][j+1]);

                Jzt00[ind][j+1] = factor2[0]  * 
                        (factor[i-1] * (RdpRds0[j][k] * RdtpRainv[i-2][0]      +
                                        indamcbypgd[0][j] * Jzt10[indm1][j]    +
                                        inddmceyptd[0][k] * Jzt10[indm1][j+1]) + 
                         onemfactor * Jzt20[indm1][j+1]);

                tempt = -bein[j+1][k+1] * tmpVeca[0] - 
                        (indaypgd[0][j] * bein[j][k+1] +
                         inddyptd[0][k] * bein[j+1][k]);


                Jztinf00[ind][j+1] = factor2[0] * 
                                     (factor[i-1] *
                                      (tempt +
                                       indamcbypgd[0][j] * Jztinf10[indm1][j]+
                                       inddmceyptd[0][k] * Jztinf10[indm1][j+1]) +
                                      onemfactor * Jztinf20[indm1][j+1]);


                tempt = bein[j+1][k+1] * tmpVecb[0] - c *
                        (inday[0][j] * bein[j][k+1] +
                         inddy[0][k] * bein[j+1][k]);


                Lytinf0[ind][j+1] = onemc2inv *
                                    (onemfactor * Jytinf20[indm1][j+1] +
                                     factor[i-1] *
                                     (jbmca[j] * Kytinf0[indm1][j]+
                                      kemcd[k] * Kytinf0[indm1][j+1]+tempt));

                /* m = 1 */
                Jyt01[ind][j+1] = factor1[1]  * 
                        (factor[i-1] * (RdpRds1[j][k] * RdtRainv[i-2][1]    + 
                                        indbmcazpgd[j] * Jyt11[indm1][j]    +
                                        indemcdzptd[k] * Jyt11[indm1][j+1]) + 
                         onemfactor * Jyt21[indm1][j+1]);


                Jzt01[ind][j+1] = factor2[1]  * 
                        (factor[i-1] * (RdpRds1[j][k] * RdtpRainv[i-2][1]      + 
                                        indamcbypgd[1][j] * Jzt11[indm1][j]    +
                                        inddmceyptd[1][k] * Jzt11[indm1][j+1]) + 
                         onemfactor * Jzt21[indm1][j+1]);

                tempt = -bein[j+1][k+1] * tmpVeca[1] -
                        (indaypgd[1][j] * bein[j][k+1]+
                         inddyptd[1][k] * bein[j+1][k]);

                Jztinf01[ind][j+1] = factor2[1] * 
                                     (factor[i-1] *
                                      (tempt +
                                       indamcbypgd[1][j] * Jztinf11[indm1][j]+
                                       inddmceyptd[1][k] * Jztinf11[indm1][j+1]) +
                                      onemfactor * Jztinf21[indm1][j+1]);

                tempt = bein[j+1][k+1] * tmpVecb[1] - c * 
                        (inday[1][j] * bein[j][k+1]+
                         inddy[1][k] * bein[j+1][k]);


                Lytinf1[ind][j+1] = onemc2inv *
                                    (onemfactor*Jytinf21[indm1][j+1]+factor[i-1] *
                                     (jbmca[j] * Kytinf1[indm1][j]+
                                      kemcd[k] * Kytinf1[indm1][j+1]+tempt));

            }  /* end for [j = 0; j < kk; ...) */

            jj = twoi-5;

            tmpVeca[0]= (jj+2) * c * y0;
            tmpVeca[1]= (jj+2) * c * y1;

            j = jj + 2;


            tmpre0 = Rdsin0[j] * Rainv0[i-2];
            tmpre1 = Rdsin1[j] * Rainv1[i-2];

            /* n = 1 */
            Jyt10[ind][1] = dfactor  *  (indd[j] * Jyt10[indm1][1]-tmpre0) +
                            emcdzptd * Jyt00[ind][1];
            Jzt10[ind][1] = efactor  *  (inde[j] * Jzt10[indm1][1]-tmpre0) +
                            dmceyptd[0] * Jzt00[ind][1];

            Jyt11[ind][1] = dfactor  *  (indd[j] * Jyt11[indm1][1]-tmpre1) +
                            emcdzptd * Jyt01[ind][1];
            Jzt11[ind][1] = efactor  *  (inde[j] * Jzt11[indm1][1]-tmpre1) +
                            dmceyptd[1] * Jzt01[ind][1];


            tmpre0 = ein[j+1] * tmpVeca[0]-inddy[0][j] * ein[j];


            Kytinf0[ind][1] = dfactor * (indd[j] * Kytinf0[indm1][1]-tmpre0) +
                              emcd * Lytinf0[ind][1]; 
      
            Jztinf10[ind][1] = efactor * (inde[j] * Jztinf10[indm1][1]+ein[j+1]) +
                               dmceyptd[0] * Jztinf00[ind][1];

            tmpre1 = ein[j+1] * tmpVeca[1]-inddy[1][j] * ein[j];


            Kytinf1[ind][1] = dfactor * (indd[j] * Kytinf1[indm1][1]-tmpre1) +
                              emcd * Lytinf1[ind][1]; 
      
            Jztinf11[ind][1] = efactor * (inde[j] * Jztinf11[indm1][1]+ein[j+1]) + 
                               dmceyptd[1] * Jztinf01[ind][1];


            kk = twoi-2;
            k = kk;

            for (j = 0; j < kk; j++) {

                k--;

                tmpcpx0 = RdpRds0[j][k] * Rainv0[i-2];
                tmpcpx1 = RdpRds1[j][k] * Rainv1[i-2];

                Jyt10[ind][j+2] = afactor  * 
                                  (indd[k] * Jyt10[indm1][j+1] +
                                   inda[j] * Jyt10[indm1][j] - tmpcpx0) +
                                  bmcazpgd * Jyt00[ind][j+1];
                Jzt10[ind][j+2] = bfactor  * 
                                  (inde[k] * Jzt10[indm1][j+1] +
                                   indb[j] * Jzt10[indm1][j] - tmpcpx0) +
                                  amcbypgd[0] * Jzt00[ind][j+1];

                Jyt11[ind][j+2] = afactor  * 
                                  (indd[k] * Jyt11[indm1][j+1] +
                                   inda[j] * Jyt11[indm1][j] - tmpcpx1) +
                                  bmcazpgd * Jyt01[ind][j+1];
                Jzt11[ind][j+2] = bfactor  * 
                                  (inde[k] * Jzt11[indm1][j+1] +
                                   indb[j] * Jzt11[indm1][j] - tmpcpx1) +
                                  amcbypgd[1] * Jzt01[ind][j+1];

                tmpcpxinf0 = bein[j+1][k+1] * tmpVeca[0] -
                             (inday[0][j] * bein[j][k+1]  +
                              inddy[0][k] * bein[j+1][k]);
                tmpcpxinf1 = bein[j+1][k+1] * tmpVeca[1] -
                             (inday[1][j] * bein[j][k+1] + 
                              inddy[1][k] * bein[j+1][k]);

                Kytinf0[ind][j+2] = afactor *
                                    (indd[k] * Kytinf0[indm1][j+1] +
                                     inda[j] * Kytinf0[indm1][j]   -
                                     tmpcpxinf0) + bmca * Lytinf0[ind][j+1];

                Jztinf10[ind][j+2] = bfactor *
                                     (inde[k] * Jztinf10[indm1][j+1] +
                                      indb[j] * Jztinf10[indm1][j]   +
                                      bein[j+1][k+1]) +
                                     amcbypgd[0] * Jztinf00[ind][j+1];


                Kytinf1[ind][j+2] = afactor *
                                    (indd[k] * Kytinf1[indm1][j+1] +
                                     inda[j] * Kytinf1[indm1][j]   - 
                                     tmpcpxinf1) +
                                    bmca * Lytinf1[ind][j+1];

                Jztinf11[ind][j+2] = bfactor *
                                     (inde[k] * Jztinf11[indm1][j+1] +
                                      indb[j] * Jztinf11[indm1][j]   +
                                      bein[j+1][k+1]) +
                                     amcbypgd[1] * Jztinf01[ind][j+1];

            }  /* END n = 1 */

            /* n = 2 */
            j = jj + 3;

            tmpre0 = Rdsin0[j] * Rainv0[i-2];
            tmpre1 = Rdsin1[j] * Rainv1[i-2];

            Jyt20[ind][1] = dfactor  *  (indd[j] * Jyt20[indm1][1]-tmpre0) +
                            emcdzptd * Jyt10[ind][1];
            Jzt20[ind][1] = efactor  *  (inde[j] * Jzt20[indm1][1]-tmpre0) +
                            dmceyptd[0] * Jzt10[ind][1];

            Jyt21[ind][1] = dfactor  *  (indd[j] * Jyt21[indm1][1]-tmpre1) +
                            emcdzptd * Jyt11[ind][1];
            Jzt21[ind][1] = efactor  *  (inde[j] * Jzt21[indm1][1]-tmpre1) +
                            dmceyptd[1] * Jzt11[ind][1];


            tmpre0 = ein[j+1] * tmpVeca[0]-inddy[0][j] * ein[j];
            tmpre1 = ein[j+1] * tmpVeca[1]-inddy[1][j] * ein[j];

            tmpreinf0 = ein[j+1] * tmpVeca[0]-inddyptd[0][j] * ein[j];
            tmpreinf1 = ein[j+1] * tmpVeca[1]-inddyptd[1][j] * ein[j];

            Jytinf20[ind][1] = dfactor * (indd[j] * Jytinf20[indm1][1]-tmpre0) +
                               emcd * Kytinf0[ind][1]; 
      
            Jztinf20[ind][1] = efactor * (inde[j] * Jztinf20[indm1][1]-tmpreinf0) +
                               dmceyptd[0] * Jztinf10[ind][1];


            Jytinf21[ind][1] = dfactor * (indd[j] * Jytinf21[indm1][1]-tmpre1) +
                               emcd * Kytinf1[ind][1]; 
      
            Jztinf21[ind][1] = efactor * (inde[j] * Jztinf21[indm1][1]-tmpreinf1) + 
                               dmceyptd[1] * Jztinf11[ind][1];

            kk = twoi-1;
            k=kk;

            for (j = 0; j < kk; j++) {

                k--;

                tmpcpx0 = RdpRds0[j][k] * Rainv0[i-2];
                tmpcpx1 = RdpRds1[j][k] * Rainv1[i-2];

                Jyt20[ind][j+2] = afactor  * 
                                  (indd[k] * Jyt20[indm1][j+1] +
                                   inda[j] * Jyt20[indm1][j] - tmpcpx0) +
                                  bmcazpgd * Jyt10[ind][j+1];
                Jzt20[ind][j+2] = bfactor  * 
                                  (inde[k] * Jzt20[indm1][j+1] +
                                   indb[j] * Jzt20[indm1][j] - tmpcpx0) +
                                  amcbypgd[0] * Jzt10[ind][j+1];

                Jyt21[ind][j+2] = afactor  * 
                                  (indd[k] * Jyt21[indm1][j+1] +
                                   inda[j] * Jyt21[indm1][j] - tmpcpx1) +
                                  bmcazpgd * Jyt11[ind][j+1];
                Jzt21[ind][j+2] = bfactor  * 
                                  (inde[k] * Jzt21[indm1][j+1] +
                                   indb[j] * Jzt21[indm1][j] - tmpcpx1) +
                                  amcbypgd[1] * Jzt11[ind][j+1];


                tmpcpx0    = bein[j+1][k+1]  *  tmpVeca[0] -
                             (inday[0][j] * bein[j][k+1] +
                              inddy[0][k] * bein[j+1][k]);

                tmpcpx1    = bein[j+1][k+1]  *  tmpVeca[1] -
                             (inday[1][j] * bein[j][k+1] +
                              inddy[1][k] * bein[j+1][k]);

                tmpcpxinf0 = bein[j+1][k+1]  *  tmpVeca[0] - 
                             (indaypgd[0][j] * bein[j][k+1] +
                              inddyptd[0][k] * bein[j+1][k]);

                tmpcpxinf1 = bein[j+1][k+1]  *  tmpVeca[1] - 
                             (indaypgd[1][j] * bein[j][k+1] +
                              inddyptd[1][k] * bein[j+1][k]);

                Jytinf20[ind][j+2] = afactor  * 
                                     (indd[k] * Jytinf20[indm1][j+1]+
                                      inda[j] * Jytinf20[indm1][j] - tmpcpx0)+ 
                                     bmca * Kytinf0[ind][j+1];

                Jztinf20[ind][j+2] = bfactor  * 
                                     (inde[k] * Jztinf20[indm1][j+1]+
                                      indb[j] * Jztinf20[indm1][j] - tmpcpxinf0)+ 
                                     amcbypgd[0] * Jztinf10[ind][j+1];


                Jytinf21[ind][j+2] = afactor  * 
                                     (indd[k] * Jytinf21[indm1][j+1]+
                                      inda[j] * Jytinf21[indm1][j] - tmpcpx1)+ 
                                     bmca * Kytinf1[ind][j+1];

                Jztinf21[ind][j+2] = bfactor  * 
                                     (inde[k] * Jztinf21[indm1][j+1]+
                                      indb[j] * Jztinf21[indm1][j] - tmpcpxinf1)+ 
                                     amcbypgd[1] * Jztinf11[ind][j+1];

           }  /*  END n = 2 */

/*
 *          Recursive double integrals
 *
 *          First construct Hyt and Hzt
 */
            indm1 = ind;
            ind   = 1-ind;

/*
 *          M = 1
 */
            kk = twoi;

            for (j = 1; j < kk; j++) {

                k = kk - (j+1);

                tmp1 =    Jzt10[indm1][j] -    Jzt11[indm1][j] -
                       Jztinf10[indm1][j] + Jztinf11[indm1][j];
                tmp2 =    Jyt10[indm1][j] -    Jyt11[indm1][j] -
                       Jytinf10[indm1][j] + Jytinf11[indm1][j];

                Hytinf0[j] = newfactor[i]  * 
                       (jamcb[j-1] * Htinf1[indm1][j-1] +
                        kdmce[k]   * Htinf1[indm1][j]   -
                        tmp1 + c * tmp2);

                Hztinf0[j] = newfactor[i]  * 
                       (jbmca[j-1] * Htinf1[indm1][j-1] +
                        kemcd[k]   * Htinf1[indm1][j]   -
                        tmp2 + c * tmp1);

            }  /* END M = 1 */

/*
 *          M = 2
 */
            kk = twoi+1;

            for (j = 1; j < kk; j++) {

                k = kk - (j+1);

                tmp1 =    Jzt20[indm1][j] -    Jzt21[indm1][j] -
                       Jztinf20[indm1][j] + Jztinf21[indm1][j];
                tmp2 =    Jyt20[indm1][j] -    Jyt21[indm1][j] -
                       Jytinf20[indm1][j] + Jytinf21[indm1][j];

                Hytinf1[j] = newfactor[i]  * 
                          (jamcb[j-1] * Htinf2[indm1][j-1] +
                           kdmce[k]   * Htinf2[indm1][j] -
                           tmp1 + c * tmp2);

                Hztinf1[j] = newfactor[i]  * 
                          (jbmca[j-1] * Htinf2[indm1][j-1] +
                           kemcd[k]   * Htinf2[indm1][j] - 
                           tmp2 + c * tmp1);

            }  /* END M = 2 */

/*
 *          Construct Ht
 */
            onem2factor = 1.0 - factor[i]  *  twoi;
            kk = twoi;
 
            for (j = 1; j < kk; j++) {

                k = kk - (j+1);

                tempt = y0 * (Jzt10[indm1][j] - Jztinf10[indm1][j]) + 
                        z0 * (Jyt10[indm1][j] -    Jyt11[indm1][j]) -
                        y1 * (Jzt11[indm1][j] - Jztinf11[indm1][j]) -
                            Kytinf0[indm1][j] +  Kytinf1[indm1][j];

                Htinf0[ind][j] = a2d2inv  * 
                                 (factor[i]  *
                                  (tempt +
                                   indgd[j-1]*Htinf1[indm1][j-1] +
                                   indtd[k]    * Htinf1[indm1][j]) +
                                  onem2factor * Htinf2[indm1][j]);
            }

            /* m = 0 */
            Htinf1[ind][1] = delta   * Hytinf0[1] +
                             epsilon * Hztinf0[1] +
                             thetad  * Htinf0[ind][1];

            for (j = 1; j < twoi; j++) {

                Htinf1[ind][j+1] = alpha  * Hytinf0[j] +
                                   beta   * Hztinf0[j] +
                                   gammad * Htinf0[ind][j];
            }           

            /* m = 1 */
            Htinf2[ind][1] = delta   * Hytinf1[1] +
                             epsilon * Hztinf1[1] +
                             thetad  * Htinf1[ind][1];

            for (j = 1; j < twoi+1; j++) {

                Htinf2[ind][j+1] = alpha  * Hytinf1[j] +
                                   beta   * Hztinf1[j] +
                                   gammad * Htinf1[ind][j];
            }           

            for (jj = 1; jj < RemArray[i]; jj++) {

                j = RelArray[local];

                reHyt = creal(Hytinf1[j]);
                reHt = creal(Htinf1[ind][j]);

                reHy1 = reHyt- y0 * reHt;
                reHy2 = reHyt- y1 * reHt;

                for (iii = 0; iii < Re18Array[local]; iii++) {
                    ii = ReiiArray[local3];
                    ReFqval = FqReal[local3];
                    G3[ii] -= reHy2 * ReFqval;
                    G4[ii] += reHy1 * ReFqval;
                    local3++;
                }

                local++;  
            }

            for (jj = 1; jj < ImmArray[i]; jj++) {

                j = ImlArray[local2];

                imHyt = cimag(Hytinf1[j]);
                imHt  = cimag(Htinf1[ind][j]);

                imHy1 = imHyt- y0 * imHt;
                imHy2 = imHyt- y1 * imHt;

                for (iii = 0; iii < Im18Array[local2]; iii++) {
                    ii = ImiiArray[local4];
                    ImFqval = FqImag[local4];
                    G3[ii] += imHy2 * ImFqval;
                    G4[ii] -= imHy1 * ImFqval;
                    local4++;
                }

                local2++;
            }

        }  /* end for (i = 2; i < qMaxp2; ...) */

/*
 *      Compute pre-factor eps C C b t bp tp in the form of
 *      matrices Ap, B, Bp, A
 */
        bptp11 = bpx * tp[0];
        bptp22 = bpy * tp[1];
        bptp33 = bpz * tp[2];

        Ap51 = bpx * tp[1]; 
        Ap62 = bpy * tp[2]; 
        Ap43 = bpz * tp[0]; 

        Ap23  = bpy * tp[0];
        Ap31  = bpz * tp[1];
        Ap12  = bpx * tp[2];

        Ap63 = bptp11 - bptp22;
        Ap41 = bptp22 - bptp33;
        Ap52 = bptp33 - bptp11;

        B31 = bx * t[1]; 
        B12 = by * t[2]; 
        B23 = bz * t[0]; 

        B24 = by * t[0]; 
        B35 = bz * t[1]; 
        B16 = bx * t[2]; 

        bt11 = bx * t[0];
        bt22 = by * t[1];
        bt33 = bz * t[2];

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

                    local = 6*j+k;

                    ApGaux3[i][j] += Ap[k][i]*G3[local];
                    ApGaux4[i][j] += Ap[k][i]*G4[local];
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

/*
 *      Deduce forces
 */
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
