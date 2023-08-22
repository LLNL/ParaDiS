/***************************************************************************
 *
 *      Module:      StressDueToSegAnisotropic.c
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

#ifdef ANISOTROPIC
#ifdef ANISO_QMAX
#define QMAX ANISO_QMAX
#else
#define QMAX 20
#endif

//---------------------------------------------------------------------------
// The IBM compiler continues to manifest errors when casting variable
// length and multidimensional matrices.  The ptrcpy() macro copies pointers
// via memcpy() to avoid the casting problems.
//---------------------------------------------------------------------------

#define ptrcpy(a,b)  memcpy(&(a),&(b),sizeof(void *))

/***************************************************************************
 *
 *      Function:    StressDueToSegAnisotropic()
 *
 *      Description: Uses anisotropic elasticity to calculate the stress
 *                   at a given point from a specified dislocation
 *                   segment.
 *
 *      Parameters:
 *          IN: point{xyz} Position at which to calculate stress from the
 *                         <p1>,<p2> segment.
 *          IN: p1{xyz}    First endpoint of segment
 *          IN: p2{xyz}    Second endpoint of segment
 *          IN: b{xyz}     Segment beurgers vector from <p1> to <p2>
 *          IN: a          Core radius
 *          IN: qmax       Base factor for determining the number of terms
 *                         in the expansion for spherical harmonics.  The
 *                         total number of terms is 2 * qmax + 1.
 *          OUT: sigma     Resulting stress
 *                         CAREFUL : The order of the stress is different than the
 *                         isotropic one
 *                         [0] = stressxx
 *                         [1] = stressyy
 *                         [2] = stresszz
 *                         [3] = stressxy
 *                         [4] = stressyz
 *                         [5] = stressxz
 *
 **************************************************************************/

void StressDueToSegAnisotropic
(
    Home_t *home,
    real8 pointx, real8 pointy, real8 pointz,
    real8 p1x, real8 p1y, real8 p1z,
    real8 p2x, real8 p2y, real8 p2z,
    real8 bx,  real8 by,  real8 bz,
    real8 a, int qmax,
    real8 sigma[6]
)
{
    int        i, j, k, ii, jj, kk, twoi;
    int        local;
    int        offsetG;
    int        ind, indm1;
    real8      d2;
    real8      a2, a2_d2, a2d2inv;
    real8      delta;
    real8      thetad;
    real8      Rasum;
    real8      onemfactor;
    real8      Ap12, Ap23, Ap31, Ap41, Ap43, Ap51, Ap52, Ap62, Ap63;
    real8      bt11, bt12, bt13;
    real8      bt21, bt22, bt23;
    real8      bt31, bt32, bt33;
    real8      tmp2;
    real8      t[3];
    real8      y[2];
    real8      R[2][3];
    real8      R2[2];
    real8      Ra[2];
    real8      Ap[6][3];
    real8      G[18];
    real8      ApGaux[3][3];
    real8      ApG[6];
    real8      factor   [  QMAX+2];
    real8      factorinv[  QMAX+2];
    real8      dfactor  [  QMAX+2];
    real8      Rdsin0   [2*QMAX+3];
    real8      Rdsin1   [2*QMAX+3];
    real8      Rainv0   [  QMAX+1];
    real8      Rainv1   [  QMAX+1];
    real8      Ra2inv   [2];
    real8      indd     [2*QMAX+3];
    real8      indtd    [2*QMAX+3];
    real8      reJyzt, imJyzt;

    complex8   alpha;
    complex8   gammad;
    complex8   ctmp1, ctmp2;
    complex8   afactor  [  QMAX+2];
    complex8   inda     [2*QMAX+3];
    complex8   indgd    [2*QMAX+3];
    complex8   Rdpin0   [2*QMAX+3];
    complex8   Rdpin1   [2*QMAX+3];
    complex8   RdpRds0  [2*QMAX+3][2*QMAX+3];
    complex8   RdpRds1  [2*QMAX+3][2*QMAX+3];
    complex8   Jyzt0    [2*QMAX+3][2];
    complex8   Jyzt1    [2*QMAX+3][2];
    complex8   Jyzt2    [2*QMAX+3][2];

    // Deal with case where qmax exceeds allocated space...

    if ( qmax>QMAX )
    {
       printf("%s::%s() ln=%d error - qmax(%d) larger than QMAX(%d) - recompile with larger QMAX!\n", __FILE__, __func__, __LINE__, qmax, QMAX );
       exit(0);
    }

    // Need some stuff from 'home'

    AnisotropicVars_t *avars = &home->anisoVars;

    real8    *fqr = (real8 *) avars->FqReal_v3;
    real8    *fqi = (real8 *) avars->FqImag_v3;
    real8   (*FqReal)[(qmax+2)*(qmax+1)]=0;  ptrcpy(FqReal,fqr);
    real8   (*FqImag)[(qmax+2)*(qmax+1)]=0;  ptrcpy(FqImag,fqi);

    complex8 *rotMatrix12   =   avars->anisoRotMatrix12;
    real8    *rotMatrixRow3 = &(avars->anisoRotMatrix[2][0]);

    real8   (*ecMatrix)[6]  = home->param->elasticConstantMatrix;

    // Initialization

    for (i = 0; i < 18; i++)
    {
        G[i] = 0.0;
    }

    // Define line direction

    t[0] = p2x - p1x;
    t[1] = p2y - p1y;
    t[2] = p2z - p1z;

    NormalizeVec(t);

    // Compute integrals via recurrence

    alpha = DotProduct(t, rotMatrix12  );
    delta = DotProduct(t, rotMatrixRow3);

    for (i = 0; i < (qmax+2); i++)
    {
        factorinv[i] = 2*i-1;
        factor[i] = 1.0 / factorinv[i];
        afactor[i] = factor[i] * alpha;
        dfactor[i] = factor[i] * delta;
    }

    R[0][0] = pointx - p1x;
    R[0][1] = pointy - p1y;
    R[0][2] = pointz - p1z;

    R[1][0] = pointx - p2x;
    R[1][1] = pointy - p2y;
    R[1][2] = pointz - p2z;

    y[0] = DotProduct(R[0], t);
    y[1] = DotProduct(R[1], t);


    Rdpin0[0] = 1.0;
    Rdpin1[0] = 1.0;
    Rdpin0[1] = DotProduct(R[0], rotMatrix12);
    Rdpin1[1] = DotProduct(R[1], rotMatrix12);

    Rdsin0[0] = 1.0;
    Rdsin1[0] = 1.0;
    Rdsin0[1] = DotProduct(R[0], rotMatrixRow3);
    Rdsin1[1] = DotProduct(R[1], rotMatrixRow3);

    for (i = 2; i < (2*qmax+3); i++)
    {
        int im1 = i-1;

        Rdpin0[i] = Rdpin0[im1] * Rdpin0[1];
        Rdpin1[i] = Rdpin1[im1] * Rdpin1[1];

        Rdsin0[i] = Rdsin0[im1] * Rdsin0[1];
        Rdsin1[i] = Rdsin1[im1] * Rdsin1[1];
    }

    gammad = Rdpin0[1] - y[0] * alpha;
    thetad = Rdsin0[1] - y[0] * delta;

    for (i = 0; i < (2*qmax+3); i++)
    {
        inda[i] = alpha * i;
        indd[i] = delta * i;
        indgd[i] = gammad * i;
        indtd[i] = thetad * i;
    }

    // Not all elements are set, so initialize all to zero at start

    for (i = 0; i < (2*qmax+3); i++)
    {
        Jyzt0[i][0] = 0.0; Jyzt0[i][1] = 0.0;
        Jyzt1[i][0] = 0.0; Jyzt1[i][1] = 0.0;
        Jyzt2[i][0] = 0.0; Jyzt2[i][1] = 0.0;
    }

    for (i = 0; i < (2*qmax+3); i++)
    {
        for (j = 0; j < (2*qmax+3); j++)
        {
            RdpRds0[i][j] = 0.0;
            RdpRds1[i][j] = 0.0;
        }
    }

    for (i = 0; i < (2*qmax+3); i++)
    {
        for (j = 0; j < 2*(qmax+2)-(i+1); j++)
        {
            RdpRds0[i][j] = Rdpin0[i] * Rdsin0[j];
            RdpRds1[i][j] = Rdpin1[i] * Rdsin1[j];
        }
    }

    R2[0] = DotProduct(R[0], R[0]);
    R2[1] = DotProduct(R[1], R[1]);

    d2 = 0.5 * (R2[0] - y[0]*y[0] +
                R2[1] - y[1]*y[1]);

    a2 = a * a;
    a2_d2 = a2 + d2;
    a2d2inv = 1.0 / a2_d2;

    Ra[0] = sqrt(R2[0] + a2);
    Ra[1] = sqrt(R2[1] + a2);
    Rasum = Ra[1] - Ra[0];

    Rainv0[0] = 1.0 / Ra[0];
    Rainv1[0] = 1.0 / Ra[1];

    Ra2inv[0] = Rainv0[0] * Rainv0[0];
    Ra2inv[1] = Rainv1[0] * Rainv1[0];

    for (i = 1; i < (qmax+1); i++)
    {
        int im1 = i-1;

        Rainv0[i] = Rainv0[im1] * Ra2inv[0];
        Rainv1[i] = Rainv1[im1] * Ra2inv[1];
    }

    Jyzt1[1][0] = log((Ra[1]+y[1]) / (Ra[0]+y[0]));
    Jyzt2[2][0] = alpha*Rasum + gammad*Jyzt1[1][0];
    Jyzt2[1][0] = delta*Rasum + thetad*Jyzt1[1][0];

    ind = 1;
    indm1 = 0;
    onemfactor=(1.0 - 2.0*factor[1]);

    ctmp1 = y[1] * RdpRds1[0][1] * Rainv1[0] -
            y[0] * RdpRds0[0][1] * Rainv0[0];

    ctmp2 = y[1] * RdpRds1[1][0] * Rainv1[0] -
            y[0] * RdpRds0[1][0] * Rainv0[0];

    Jyzt0[1][ind] = a2d2inv *
                    (factor[1] *
                     (ctmp1 + indgd[0] * Jyzt1[0][indm1] +
                      indtd[1] * Jyzt1[1][indm1]) +
                     onemfactor * Jyzt2[1][indm1]);

    Jyzt0[2][ind] = a2d2inv *
                    (factor[1] *
                     (ctmp2 + indgd[1] * Jyzt1[1][indm1] +
                      indtd[0] * Jyzt1[2][indm1]) +
                     onemfactor * Jyzt2[2][indm1]);


    local = 0;
    reJyzt = creal(Jyzt0[1][ind]);
    imJyzt = cimag(Jyzt0[1][ind]);
    for (ii = 0; ii < 18; ii++)
    {
        G[ii] += reJyzt*FqReal[ii][local];
        G[ii] -= imJyzt*FqImag[ii][local];
    }

    local = 1;
    reJyzt = creal(Jyzt0[2][ind]);
    imJyzt = cimag(Jyzt0[2][ind]);
    for (ii = 0; ii < 18; ii++)
    {
        G[ii] += reJyzt*FqReal[ii][local];
        G[ii] -= imJyzt*FqImag[ii][local];
    }

    // Recursive integrals

    local  = 2;

    for (i = 2; i < (qmax+2); i++)
    {
        twoi = 2*i;
        jj = twoi-4;

        /* n = 0 */
        j = jj+1;

        tmp2 = Rdsin1[j] * Rainv1[i-2] -
               Rdsin0[j] * Rainv0[i-2];

        Jyzt1[1][ind] = dfactor[i-1] *
                        (indd[j] * Jyzt1[1][indm1] - tmp2) +
                        thetad * Jyzt0[1][ind];

        kk = twoi-2;
        k = kk;
        for (j = 0; j < kk; j++)
        {

            k--;

            ctmp1 = RdpRds1[j][k] * Rainv1[i-2] -
                    RdpRds0[j][k] * Rainv0[i-2];

            Jyzt1[j+2][ind] = afactor[i-1] *
                              (indd[k] * Jyzt1[j+1][indm1] +
                               inda[j] * Jyzt1[j][indm1] -ctmp1) +
                              gammad * Jyzt0[j+1][ind];
        }

        /* n = 1 */
        j = jj+2;

        tmp2 = Rdsin1[j] * Rainv1[i-2] -
               Rdsin0[j] * Rainv0[i-2];

        Jyzt2[1][ind] = dfactor[i-1] *
                        (indd[j] * Jyzt2[1][indm1] - tmp2) +
                        thetad * Jyzt1[1][ind];

        kk = twoi-1;
        k = kk;

        for (j = 0; j < kk; j++)
        {

            k--;

            ctmp1 = RdpRds1[j][k] * Rainv1[i-2] -
                    RdpRds0[j][k] * Rainv0[i-2];

            Jyzt2[j+2][ind] = afactor[i-1] *
                              (indd[k] * Jyzt2[j+1][indm1] +
                               inda[j] * Jyzt2[j][indm1] -ctmp1) +
                              gammad * Jyzt1[j+1][ind];
        }


        indm1 = ind;
        ind = 1 - ind;

        kk = twoi;
        onemfactor=(1.0e0-kk*factor[i]);

        k = kk;

        for (j = 0; j < kk; j++)
        {

            k--;

            ctmp1 = y[1] * RdpRds1[j][k] * Rainv1[i-1] -
                    y[0] * RdpRds0[j][k] * Rainv0[i-1];

            Jyzt0[j+1][ind] = a2d2inv *
                              (factor[i] * (ctmp1 + indgd[j] * Jyzt1[j][indm1] +
                                            indtd[k] * Jyzt1[j+1][indm1]) +
                               onemfactor * Jyzt2[j+1][indm1]);
        }


        for (j = 1; j < 2*i+1; j++)
        {
            reJyzt = creal(Jyzt0[j][ind]);
            imJyzt = cimag(Jyzt0[j][ind]);

            for (ii = 0; ii < 18; ii++)
            {
                G[ii] += reJyzt*FqReal[ii][local];
                G[ii] -= imJyzt*FqImag[ii][local];
            }

            local++;
        }

    }  /* end for (i = 2; i < q; i++) */


    // Compute pre-factor eps C C b t bp tp in the form of
    // matrices Ap, B, Bp, A

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

    // For hexagonal or cubic elastic constant matrix

    for (i = 0; i < 6; i++)
    {
        for (j = 0; j < 3; j++)
        {
            Ap[i][j] = 0.0;
        }
    }

    for (j = 0; j < 3; j++)
    {
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

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            ApGaux[i][j] = 0.0;

            for (k = 0; k < 6; k++)
            {
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

    for (i = 0; i < 6; i++)
    {
        sigma[i] = 0.0;

        for (j = 0; j < 6; j++)
        {
            sigma[i] += ecMatrix[i][j] * ApG[j];
        }

        sigma[i] *= INV_4_PI_SQUARED;

        // printf("sigma[%d]=%g\n",i,sigma[i]);
    }

    return;
}
#endif // ifdef ANISOTROPIC
