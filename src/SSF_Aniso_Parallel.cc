#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mpi_portability.h"
#include "cuda_portability.h"

#include "Typedefs.h"
#include "Complex_t.h"
#include "V3.h"

#include "SSF_Aniso.h"

//------------------------------------------------------------------------------------------------------------

#ifdef ANISO_QMAX
#define QMAX ANISO_QMAX
#else
#define QMAX 20
#endif

//------------------------------------------------------------------------------------------------------------

#define VZERO(a,n)   for (int i=0; (i<(n)); ++i) { (a)[i]=0.0; }

#define MZERO(a,n,m) for (int i=0; (i<(n)); ++i) \
                     for (int j=0; (j<(m)); ++j) { (a)[i][j]=0.0; }

#define M33_ZERO(a)  for (int i=0; (i<3); ++i)   { (a)[i][0]=0.0;  (a)[i][1]=0.0;  (a)[i][2]=0.0;  }

//------------------------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Aniso_Parallel
(
          real8    *f3    ,   ///< force on node 3 <xyz>   (returned)
          real8    *f4    ,   ///< force on node 4 <xyz>   (returned)
    const real8    *p1    ,   ///< position of node 1 <xyz>
    const real8    *p2    ,   ///< position of node 2 <xyz>
    const real8    *p3    ,   ///< position of node 3 <xyz>
    const real8    *p4    ,   ///< position of node 4 <xyz>
    const real8    *b12   ,   ///< burgers vector (p1->p2) <xyz>
    const real8    *b34   ,   ///< burgers vector (p3->p4) <xyz>
    const real8     a     ,   ///< core radius (b)
    const int       qmax  ,   ///< spherical harmonics expansion factor
    const int       mtype ,   ///< material type index
    const real8    *c66   ,   ///< elastic constants matrix (6x6)
    const complex8 *c12   ,   ///< complex elements of aniso rotation matrix (1x3,complex)
    const real8    *e3    ,   ///< real    elements of aniso rotation matrix (1x3,real   )
    const real8    *fqr   ,   ///< Fq table (real)      size=18*(qmax+2)*(qmax+1)
    const real8    *fqi       ///< Fq table (imaginary) size=18*(qmax+2)*(qmax+1)
)
{
    real8   (*ecm   )[6]                 = (real8 (*)[6])                 c66;
    real8   (*FqReal)[(QMAX+2)*(QMAX+1)] = (real8 (*)[(QMAX+2)*(QMAX+1)]) fqr;
    real8   (*FqImag)[(QMAX+2)*(QMAX+1)] = (real8 (*)[(QMAX+2)*(QMAX+1)]) fqi;

    int        k, jj, kk;
    int        ind, indm1;
    real8      y  [2];
    real8      z  [2];
    real8      yin[4];
    real8      zin[4];
    real8      ypz[4];
    real8      ymz[4];
    real8      R2 [4];
    real8      Ra [4];
    real8      Ra2inv[4];

    real8      G3       [18];
    real8      G4       [18];
    real8      factor   [  QMAX+2];
    real8      dfactor  [  QMAX+2];
    real8      Rdsin    [2*QMAX+3][4];
    real8      indd     [2*QMAX+3];
    real8      indtd    [2*QMAX+3];
    real8      Rainv    [  QMAX+1][4];

    complex8   afactor  [  QMAX+2];
    complex8   Rdpin    [2*QMAX+3][4];
    complex8   inda     [2*QMAX+3];
    complex8   indgd    [2*QMAX+3];
    complex8   Ht       [2*QMAX+4][3][2];
    complex8   Hyt      [2*QMAX+3][3][2];
    complex8   Hzt      [2*QMAX+3][3][2];
    complex8   Jyzt     [2*QMAX+3][3][2][4];

    VZERO(G3,18);
    VZERO(G4,18);

    real8 p12[3]  = { p2[0]-p1[0],
                      p2[1]-p1[1],
                      p2[2]-p1[2] };

    real8 p34[3]  = { p4[0]-p3[0],
                      p4[1]-p3[1],
                      p4[2]-p3[2] };

    real8 p12_mag = V3_MAG(p12);
    real8 p34_mag = V3_MAG(p34);

    real8 invL    = ( (p12_mag>0) ? (1.0/p12_mag) : 0.0 );

    real8 tp[3]   = { p12[0]*invL,
                      p12[1]*invL,
                      p12[2]*invL };

          invL    = ( (p34_mag>0) ? (1.0/p34_mag) : 0.0 );
    real8 coeff34 = invL / (4.0*M_PI*M_PI);

    real8 t[3]    = { p34[0]*invL,
                      p34[1]*invL,
                      p34[2]*invL };

    // Compute integrals via recurrence

    complex8 alpha = V3_DOT(t,c12);
    real8    delta = V3_DOT(t,e3 );

    for (int i=0; i<(qmax+2); ++i) { factor [i] = 1.0/(2.0*(i+1)-3.0); }
    for (int i=0; i<(qmax+2); ++i) { afactor[i] = alpha * factor[i]; }
    for (int i=0; i<(qmax+2); ++i) { dfactor[i] = delta * factor[i]; }

    real8 R[4][3]   = { { p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2] },
                        { p3[0]-p2[0], p3[1]-p2[1], p3[2]-p2[2] },
                        { p4[0]-p1[0], p4[1]-p1[1], p4[2]-p1[2] },
                        { p4[0]-p2[0], p4[1]-p2[1], p4[2]-p2[2] } };

    real8 center[3] = { 0.25*(p1[0]+p2[0]+p3[0]+p4[0]),
                        0.25*(p1[1]+p2[1]+p3[1]+p4[1]),
                        0.25*(p1[2]+p2[2]+p3[2]+p4[2]) };

    y[0] =  (((p3[0] - center[0]) * t[0]) +
             ((p3[1] - center[1]) * t[1]) +
             ((p3[2] - center[2]) * t[2]));

    y[1] =  (((p4[0] - center[0]) * t[0]) +
             ((p4[1] - center[1]) * t[1]) +
             ((p4[2] - center[2]) * t[2]));

    z[0] = -(((p1[0] - center[0]) * t[0]) +
             ((p1[1] - center[1]) * t[1]) +
             ((p1[2] - center[2]) * t[2]));

    z[1] = -(((p2[0] - center[0]) * t[0]) +
             ((p2[1] - center[1]) * t[1]) +
             ((p2[2] - center[2]) * t[2]));

    yin[0] = y[0];
    yin[1] = y[0];
    yin[2] = y[1];
    yin[3] = y[1];

    zin[0] = z[0];
    zin[1] = z[1];
    zin[2] = z[0];
    zin[3] = z[1];

    for (int i=0; (i<4); ++i) { ypz[i] = yin[i]+zin[i]; }
    for (int i=0; (i<4); ++i) { ymz[i] = yin[i]-zin[i]; }

    for (int i=0; (i<4); ++i)
    {
        Rdpin[0][i] = 1.0;
        Rdsin[0][i] = 1.0;

        Rdpin[1][i] = V3_DOT(R[i],c12);
        Rdsin[1][i] = V3_DOT(R[i],e3 );
    }

    for (int i=2; i<(2*qmax+3); ++i)
    {
        Rdpin[i][0] = Rdpin[i-1][0] * Rdpin[1][0];
        Rdpin[i][1] = Rdpin[i-1][1] * Rdpin[1][1];
        Rdpin[i][2] = Rdpin[i-1][2] * Rdpin[1][2];
        Rdpin[i][3] = Rdpin[i-1][3] * Rdpin[1][3];

        Rdsin[i][0] = Rdsin[i-1][0] * Rdsin[1][0];
        Rdsin[i][1] = Rdsin[i-1][1] * Rdsin[1][1];
        Rdsin[i][2] = Rdsin[i-1][2] * Rdsin[1][2];
        Rdsin[i][3] = Rdsin[i-1][3] * Rdsin[1][3];
    }

    complex8  gammad = Rdpin[1][0] - ypz[0] * alpha;
    real8     thetad = Rdsin[1][0] - ypz[0] * delta;

    for (int i=0; i<(2*qmax+3); ++i) { inda [i] = i*alpha ; }
    for (int i=0; i<(2*qmax+3); ++i) { indd [i] = i*delta ; }
    for (int i=0; i<(2*qmax+3); ++i) { indgd[i] = i*gammad; }
    for (int i=0; i<(2*qmax+3); ++i) { indtd[i] = i*thetad; }

    for (int i=0; (i<4); ++i) { R2[i] = V3_DOT(R[i],R[i]); }

    real8 d2      = 0.5 * (R2[0] - ypz[0]*ypz[0] +
                           R2[3] - ypz[3]*ypz[3]);

    real8 a2      = a*a;
    real8 a2_d2   = a2+d2;
    real8 a2d2inv = 1.0/a2_d2;

    for (int i=0; (i<4); ++i) { Ra      [i] = sqrt(R2[i] + a2); }
    for (int i=0; (i<4); ++i) { Rainv[0][i] = 1.0 / Ra[i]; }
    for (int i=0; (i<4); ++i) { Ra2inv  [i] = Rainv[0][i] * Rainv[0][i]; }

    for (int i=1; i<(qmax+1); ++i)
    {
        Rainv[i][0] = Rainv[i-1][0] * Ra2inv[0];
        Rainv[i][1] = Rainv[i-1][1] * Ra2inv[1];
        Rainv[i][2] = Rainv[i-1][2] * Ra2inv[2];
        Rainv[i][3] = Rainv[i-1][3] * Ra2inv[3];
    }

    for (int i=0; (i<4); ++i) { Jyzt[0][0][0][i] = 0.0; }

    for (int i=0; (i<4); ++i)
    {
        Jyzt[1][1][0][i] = log(Ra[i]+ypz[i]);
        Jyzt[2][2][0][i] = alpha*Ra[i] + gammad*Jyzt[1][1][0][i];
        Jyzt[1][2][0][i] = delta*Ra[i] + thetad*Jyzt[1][1][0][i];
    }

    real8    Rasum   = Ra[0] - Ra[1] - Ra[2] + Ra[3];
    real8    tmpsum  = a2d2inv * (ymz[0]*Ra[0] - ymz[1]*Ra[1] -
                                  ymz[2]*Ra[2] + ymz[3]*Ra[3]);

    complex8 Jyzsum  = Jyzt[1][1][0][0] - Jyzt[1][1][0][1] -
                       Jyzt[1][1][0][2] + Jyzt[1][1][0][3];

    complex8 yJyzsum = yin[0]*Jyzt[1][1][0][0] - yin[1]*Jyzt[1][1][0][1] -
                       yin[2]*Jyzt[1][1][0][2] + yin[3]*Jyzt[1][1][0][3];

    complex8 zJyzsum = zin[0]*Jyzt[1][1][0][0] - zin[1]*Jyzt[1][1][0][1] -
                       zin[2]*Jyzt[1][1][0][2] + zin[3]*Jyzt[1][1][0][3];

    for (int i=0; i<(2*qmax+3); ++i)
    {
        for (int j=0; j<3; ++j)
        {
            Hyt[i][j][0] = 0.0;
            Hyt[i][j][1] = 0.0;
            Hzt[i][j][0] = 0.0;
            Hzt[i][j][1] = 0.0;
        }
    }

    Hyt[1][1][1] = -0.5 * (Jyzsum + tmpsum);
    Hzt[1][1][1] = -0.5 * (Jyzsum - tmpsum);

    complex8 tmp1 = zJyzsum - Rasum;
    complex8 tmp2 = yJyzsum - Rasum;

    Ht [1][0][1] = Rasum * a2d2inv;

    Hyt[2][2][1] = alpha*tmp1 + gammad*Hyt[1][1][1];
    Hyt[1][2][1] = delta*tmp1 + thetad*Hyt[1][1][1];
    Hzt[2][2][1] = alpha*tmp2 + gammad*Hzt[1][1][1];
    Hzt[1][2][1] = delta*tmp2 + thetad*Hzt[1][1][1];

    Ht [1][1][1] = delta*(Hyt[1][1][1] + Hzt[1][1][1]) + thetad*Ht[1][0][1];
    Ht [1][2][1] = delta*(Hyt[1][2][1] + Hzt[1][2][1]) + thetad*Ht[1][1][1];
    Ht [2][1][1] = alpha*(Hyt[1][1][1] + Hzt[1][1][1]) + gammad*Ht[1][0][1];
    Ht [2][2][1] = alpha*(Hyt[1][2][1] + Hzt[1][2][1]) + gammad*Ht[1][1][1];
    Ht [3][2][1] = alpha*(Hyt[2][2][1] + Hzt[2][2][1]) + gammad*Ht[2][1][1];

    int      local  = 0;
    complex8 alpha0 = Hyt[1][2][1]-y[0]*Ht[1][1][1];
    complex8 beta0  = Hyt[1][2][1]-y[1]*Ht[1][1][1];

    real8    reHy1  = creal(alpha0);
    real8    reHy2  = creal(beta0 );
    real8    imHy1  = cimag(alpha0);
    real8    imHy2  = cimag(beta0 );

    for (int i=0; (i<18); ++i)
    {
        real8 fqr = FqReal[i][local];
        real8 fqi = FqImag[i][local];

        G3[i] += (imHy2*fqi - reHy2*fqr);
        G4[i] += (reHy1*fqr - imHy1*fqi);
    }

    local  = 1;
    alpha0 = Hyt[2][2][1]-y[0]*Ht[2][1][1];
    beta0  = Hyt[2][2][1]-y[1]*Ht[2][1][1];

    reHy1  = creal(alpha0);
    reHy2  = creal(beta0 );
    imHy1  = cimag(alpha0);
    imHy2  = cimag(beta0 );

    for (int i=0; (i<18); ++i)
    {
        real8 fqr = FqReal[i][local];
        real8 fqi = FqImag[i][local];

        G3[i] += (imHy2*fqi - reHy2*fqr);
        G4[i] += (reHy1*fqr - imHy1*fqi);
    }

    ind   = 1;
    indm1 = 0;
    local = 2;

    // Recursive integrals

    for (int i=2; i<(qmax+2); i++)
    {
        real8  onemfactor = 1.0 - (2*(i+1)-4) * factor[i-1];

        kk = 2*i-2;
        for (int j=0; j<kk; ++j)
        {
            k = kk-(j+1);
            for (int m=0; m<4; ++m)
            {
                Jyzt[j+1][0][ind][m] = a2d2inv *
                                       (factor[i-1] *
                                        (ypz[m]*Rdpin[j][m]*Rdsin[k][m]*Rainv[i-2][m] +
                                         indgd[j]  *Jyzt[j  ][1][indm1][m] +
                                         indtd[k]  *Jyzt[j+1][1][indm1][m]) +
                                         onemfactor*Jyzt[j+1][2][indm1][m]);
            }
        }

        jj = 2*i-4;
        for (int n=0; n<2; ++n)
        {
            int j = jj+n+1;
            for (int m=0; m<4; ++m)
            {
                tmp2 = Rdsin[j][m]*Rainv[i-2][m];
                Jyzt[1][n+1][ind][m] = dfactor[i-1] * (indd[j] * Jyzt[1][n+1][indm1][m] - tmp2) + thetad * Jyzt[1][n][ind][m];
            }
            kk = 2*i-3+(n+1);
            for (int j=0; j<kk; ++j)
            {
                k = kk - (j+1);
                for (int m=0; m<4; ++m)
                {
                    tmp1 = Rdpin[j][m]*Rdsin[k][m]*Rainv[i-2][m];
                    Jyzt[j+2][n+1][ind][m] = afactor[i-1] * (indd[k] * Jyzt[j+1][n+1][indm1][m] + inda[j] * Jyzt[j][n+1][indm1][m] - tmp1) + gammad*Jyzt[j+1][n][ind][m];
                }
            }
        }

        indm1 =   ind;
        ind   = 1-ind;

        real8 onem2factor = 1.0 - (2*i-1) * factor[i];

        kk = 2*i-2;
        for (int j=0; (j<kk); ++j)
        {
            complex8 tmpv4[4];

            k = kk - (j+1);
            for (int m=0; (m<4); ++m)
            {
                tmpv4[m] = inda[j] * Jyzt[j][1][ind][m] + indd[k] * Jyzt[j+1][1][ind][m] - Rdpin[j][m]*Rdsin[k][m] * Rainv[i-2][m];
            }

            tmp1 = factor[i-1] * (yin[0] * tmpv4[0] - yin[1] * tmpv4[1] - yin[2] * tmpv4[2] + yin[3] * tmpv4[3]);
            tmp2 = factor[i-1] * (zin[0] * tmpv4[0] - zin[1] * tmpv4[1] - zin[2] * tmpv4[2] + zin[3] * tmpv4[3]);

            Hyt[j+1][0][ind] = a2d2inv * (onem2factor * Hyt[j+1][2][indm1] + factor[i] * (tmp1 + indgd[j] * Hyt[j][1][indm1] + indtd[k] * Hyt[j+1][1][indm1] - Hzt[j+1][2][indm1]));
            Hzt[j+1][0][ind] = a2d2inv * (onem2factor * Hzt[j+1][2][indm1] + factor[i] * (tmp2 + indgd[j] * Hzt[j][1][indm1] + indtd[k] * Hzt[j+1][1][indm1] - Hyt[j+1][2][indm1]));
        }

        jj = 2*i-4;
        for (int n=0; n<2; n++)
        {
            int j = jj+n+1;

            tmp1 = yin[0]*Jyzt[1][n][indm1][0] -
                   yin[1]*Jyzt[1][n][indm1][1] -
                   yin[2]*Jyzt[1][n][indm1][2] +
                   yin[3]*Jyzt[1][n][indm1][3];

            tmp2 = zin[0]*Jyzt[1][n][indm1][0] -
                   zin[1]*Jyzt[1][n][indm1][1] -
                   zin[2]*Jyzt[1][n][indm1][2] +
                   zin[3]*Jyzt[1][n][indm1][3];

            Hyt[1][n+1][ind] = dfactor[i] * (Ht[1][n+1][indm1] + indd[j]*Hyt[1][n+1][indm1] - tmp1) + thetad*Hyt[1][n][ind];
            Hzt[1][n+1][ind] = dfactor[i] * (Ht[1][n+1][indm1] + indd[j]*Hzt[1][n+1][indm1] - tmp2) + thetad*Hzt[1][n][ind];

            kk = 2*i-2+(n+1);

            for (int j=1; j<kk; ++j)
            {
                k = kk - (j+1);

                tmp1 = yin[0]*Jyzt[j][n][indm1][0] -
                       yin[1]*Jyzt[j][n][indm1][1] -
                       yin[2]*Jyzt[j][n][indm1][2] +
                       yin[3]*Jyzt[j][n][indm1][3];

                Hyt[j+1][n+1][ind] = afactor[i] * (Ht[j][n+1][indm1] + inda[j-1]*Hyt[j-1][n+1][indm1] + indd[k]*Hyt[j][n+1][indm1] - tmp1) + gammad*Hyt[j][n][ind];

                tmp2 = zin[0]*Jyzt[j][n][indm1][0] -
                       zin[1]*Jyzt[j][n][indm1][1] -
                       zin[2]*Jyzt[j][n][indm1][2] +
                       zin[3]*Jyzt[j][n][indm1][3];

                Hzt[j+1][n+1][ind] = afactor[i] * (Ht[j][n+1][indm1] + inda[j-1]*Hzt[j-1][n+1][indm1] + indd[k]*Hzt[j][n+1][indm1] - tmp2) + gammad*Hzt[j][n][ind];
            }
        }

        // Construct H

        onem2factor = 1.0 - factor[i] * (2*(i+1)-2);

        kk = 2*i;

        for (int j=1; (j<kk); ++j)
        {
            k = kk - (j+1);

            tmp1 = ypz[0]*Jyzt[j][1][indm1][0] -
                   ypz[1]*Jyzt[j][1][indm1][1] -
                   ypz[2]*Jyzt[j][1][indm1][2] +
                   ypz[3]*Jyzt[j][1][indm1][3];

            Ht[j][0][ind] = a2d2inv * (factor[i] *
                                       (tmp1 + indgd[j-1]*Ht[j-1][1][indm1] +
                                        indtd[k]*Ht[j][1][indm1]) +
                                       onem2factor*Ht[j][2][indm1]);
        }

        jj = 2*i-3;

        for (int m=1; (m<3); ++m)
        {
            Ht[1][m][ind] = delta  * (Hyt[1][m][ind] + Hzt[1][m][ind]) +
                            thetad * Ht[1][m-1][ind];

            for (int j=1; j<(jj+m+2); ++j)
            {
                Ht[j+1][m][ind] =   alpha  * (Hyt[j][m  ][ind] + Hzt[j][m][ind])
                                  + gammad * (Ht [j][m-1][ind]                 );
            }
        }

        for (int j=1; j<(2*i+1); ++j)
        {
            alpha0 = Hyt[j][2][ind] - y[0]*Ht[j][1][ind];
            beta0  = Hyt[j][2][ind] - y[1]*Ht[j][1][ind];

            real8 reHy1 = creal(alpha0);
            real8 reHy2 = creal(beta0 );
            real8 imHy1 = cimag(alpha0);
            real8 imHy2 = cimag(beta0 );

            for (int m=0; (m<18); ++m)
            {
                real8 fqr = FqReal[m][local];
                real8 fqi = FqImag[m][local];

                G3[m] += imHy2*fqi - reHy2*fqr;
                G4[m] += reHy1*fqr - imHy1*fqi;
            }

            local++;
        }
    }

    // Compute pre-factor eps C C b t bp tp in the form of
    // matrices Ap, B, Bp, A

    real8 bptp11 = b12[0]*tp[0];
    real8 bptp22 = b12[1]*tp[1];
    real8 bptp33 = b12[2]*tp[2];

    real8 bptp12 = b12[0]*tp[1];
    real8 bptp23 = b12[1]*tp[2];
    real8 bptp31 = b12[2]*tp[0];

    real8 bptp21 = b12[1]*tp[0];
    real8 bptp32 = b12[2]*tp[1];
    real8 bptp13 = b12[0]*tp[2];

    real8 Ap12   = bptp13;
    real8 Ap62   = bptp23;
    real8 Ap51   = bptp12;
    real8 Ap23   = bptp21;
    real8 Ap43   = bptp31;
    real8 Ap31   = bptp32;

    real8 Ap63   = bptp11 - bptp22;
    real8 Ap41   = bptp22 - bptp33;
    real8 Ap52   = bptp33 - bptp11;

    real8 bt12   = b34[0]*t[1];
    real8 bt23   = b34[1]*t[2];
    real8 bt31   = b34[2]*t[0];

    real8 bt21   = b34[1]*t[0];
    real8 bt32   = b34[2]*t[1];
    real8 bt13   = b34[0]*t[2];

    real8 bt11   = b34[0]*t[0];
    real8 bt22   = b34[1]*t[1];
    real8 bt33   = b34[2]*t[2];

    real8 B12    = bt23;
    real8 B35    = bt32;
    real8 B31    = bt12;
    real8 B16    = bt13;
    real8 B23    = bt31;
    real8 B24    = bt21;

    real8 B36    = bt22 - bt11;
    real8 B14    = bt33 - bt22;
    real8 B25    = bt11 - bt33;

    real8 A [6][3];
    real8 Ap[6][3];
    real8 B [3][6];
    real8 ApG3[6];
    real8 ApG4[6];

#if 0
    // For general C
    {
        for (int j=0; j<6; ++j)
        {
            Ap[j][0] = -ecm[j][1]*Ap62 + ecm[j][2]*Ap31 + ecm[j][3]*Ap41 + ecm[j][4]*Ap51 - ecm[j][5]*Ap12;
            Ap[j][1] =  ecm[j][0]*Ap12 - ecm[j][2]*Ap43 - ecm[j][3]*Ap23 + ecm[j][4]*Ap52 + ecm[j][5]*Ap62;
            Ap[j][2] = -ecm[j][0]*Ap51 + ecm[j][1]*Ap23 + ecm[j][3]*Ap43 - ecm[j][4]*Ap31 + ecm[j][5]*Ap63;

            B[0][j]  =  ecm[1][j]*B12  - ecm[2][j]*B35  + ecm[3][j]*B14  - ecm[4][j]*B31  + ecm[5][j]*B16;
            B[1][j]  = -ecm[0][j]*B16  + ecm[2][j]*B23  + ecm[3][j]*B24  + ecm[4][j]*B25  - ecm[5][j]*B12;
            B[2][j]  =  ecm[0][j]*B31  - ecm[1][j]*B24  - ecm[3][j]*B23  + ecm[4][j]*B35  + ecm[5][j]*B36;
        }
    }
#endif

#if 0
    // For symetric C
    {
        for (int j=0; j<6; ++j)
        {
            Ap[j][0] = -ecm[j][1]*Ap62 + ecm[j][2]*Ap31 + ecm[j][3]*Ap41 +
                        ecm[j][4]*Ap51 - ecm[j][5]*Ap12;
            Ap[j][1] =  ecm[j][0]*Ap12 - ecm[j][2]*Ap43 - ecm[j][3]*Ap23 +
                        ecm[j][4]*Ap52 + ecm[j][5]*Ap62;
            Ap[j][2] = -ecm[j][0]*Ap51 + ecm[j][1]*Ap23 + ecm[j][3]*Ap43 -
                        ecm[j][4]*Ap31 + ecm[j][5]*Ap63;

            A[j][0]  = -ecm[j][1]*B12 + ecm[j][2]*B35 - ecm[j][3]*B14 +
                        ecm[j][4]*B31 - ecm[j][5]*B16;
            A[j][1]  =  ecm[j][0]*B16 - ecm[j][2]*B23 - ecm[j][3]*B24 -
                        ecm[j][4]*B25 + ecm[j][5]*B12;
            A[j][2]  = -ecm[j][0]*B31 + ecm[j][1]*B24 + ecm[j][3]*B23 -
                        ecm[j][4]*B35 - ecm[j][5]*B36;
        }

        for (int i=0; i<3; ++i)
        for (int j=0; j<6; ++j)
           B[i][j]  = -A[j][i];
    }
#endif


#if 1
    // For hexa or cubic C
    {
        for (int i=0; i<6; ++i)
        {
           V3_ZERO(A [i]);
           V3_ZERO(Ap[i]);
        }

        for (int j=0; j<3; ++j)
        {
            Ap[j][0] = -ecm[j][1]*Ap62 + ecm[j][2]*Ap31;
            Ap[j][1] =  ecm[j][0]*Ap12 - ecm[j][2]*Ap43;
            Ap[j][2] = -ecm[j][0]*Ap51 + ecm[j][1]*Ap23;

            A [j][0] = -ecm[j][1]*B12  + ecm[j][2]*B35;
            A [j][1] =  ecm[j][0]*B16  - ecm[j][2]*B23;
            A [j][2] = -ecm[j][0]*B31  + ecm[j][1]*B24;
        }

        Ap[3][0] =  ecm[3][3]*Ap41;
        Ap[3][1] = -ecm[3][3]*Ap23;
        Ap[3][2] =  ecm[3][3]*Ap43;

        A [3][0] = -ecm[3][3]*B14;
        A [3][1] = -ecm[3][3]*B24;
        A [3][2] =  ecm[3][3]*B23;

        Ap[4][0] =  ecm[4][4]*Ap51;
        Ap[4][1] =  ecm[4][4]*Ap52;
        Ap[4][2] = -ecm[4][4]*Ap31;

        A [4][0] =  ecm[4][4]*B31;
        A [4][1] = -ecm[4][4]*B25;
        A [4][2] = -ecm[4][4]*B35;

        Ap[5][0] = -ecm[5][5]*Ap12;
        Ap[5][1] =  ecm[5][5]*Ap62;
        Ap[5][2] =  ecm[5][5]*Ap63;

        A [5][0] = -ecm[5][5]*B16;
        A [5][1] =  ecm[5][5]*B12;
        A [5][2] = -ecm[5][5]*B36;

        for (int i=0; i<3; ++i)
        for (int j=0; j<6; ++j)
           B[i][j] = -A[j][i];
    }
#endif

    real8  ApGaux3[3][3];   M33_ZERO(ApGaux3);
    real8  ApGaux4[3][3];   M33_ZERO(ApGaux4);

    for (int i=0; (i<3); ++i)
    {
        for (int j=0; (j<3); ++j)
        {
            ApGaux3[i][j] = 0.0;
            ApGaux4[i][j] = 0.0;

            for (int k=0; (k<6); ++k)
            {
                local = 6*j+k;
                ApGaux3[i][j] = ApGaux3[i][j] + Ap[k][i]*G3[local];
                ApGaux4[i][j] = ApGaux4[i][j] + Ap[k][i]*G4[local];
            }
        }
    }

    ApG3[0] = ApGaux3[0][0];
    ApG3[1] = ApGaux3[1][1];
    ApG3[2] = ApGaux3[2][2];

    ApG3[3] = ApGaux3[1][2] + ApGaux3[2][1];
    ApG3[4] = ApGaux3[2][0] + ApGaux3[0][2];
    ApG3[5] = ApGaux3[0][1] + ApGaux3[1][0];

    ApG4[0] = ApGaux4[0][0];
    ApG4[1] = ApGaux4[1][1];
    ApG4[2] = ApGaux4[2][2];

    ApG4[3] = ApGaux4[1][2] + ApGaux4[2][1];
    ApG4[4] = ApGaux4[2][0] + ApGaux4[0][2];
    ApG4[5] = ApGaux4[0][1] + ApGaux4[1][0];

    if (f3)
    {
        V3_ZERO(f3);

        for (int j=0; (j<6); ++j)
        {
            f3[0] += B[0][j] * ApG3[j];
            f3[1] += B[1][j] * ApG3[j];
            f3[2] += B[2][j] * ApG3[j];
        }

        f3[0] *= coeff34;
        f3[1] *= coeff34;
        f3[2] *= coeff34;
    }

    if (f4)
    {
        V3_ZERO(f4);

        for (int j=0; (j<6); ++j)
        {
            f4[0] += B[0][j] * ApG4[j];
            f4[1] += B[1][j] * ApG4[j];
            f4[2] += B[2][j] * ApG4[j];
        }

        f4[0] *= coeff34;
        f4[1] *= coeff34;
        f4[2] *= coeff34;
    }
}

// SSF_Aniso_Para_Approx()
//
// Approximate the force on segment 2 from the near-parallel segment 1 using a
// quadratic shape function.
//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Aniso_Parallel_Approx
(
          real8    *f3    ,   ///< approximate force on node 3 <xyz> (returned)
          real8    *f4    ,   ///< approximate force on node 4 <xyz> (returned)
    const real8    *p1    ,   ///< postion of node 1 <xyz>
    const real8    *p2    ,   ///< postion of node 2 <xyz>
    const real8    *p3    ,   ///< postion of node 3 <xyz>
    const real8    *p4    ,   ///< postion of node 4 <xyz>
    const real8    *b12   ,   ///< burgers vector (p1 to p2) <xyz>
    const real8    *b34   ,   ///< burgers vector (p3 to p4) <xyz>
    const real8     a     ,   ///< core radius
    const int       qmax  ,   ///< base factor for spherical harmonics expansion
    const int       mtype ,   ///< material type index
    const real8     ecrit ,   ///< critical angle for parallelism test
    const real8    *c66   ,
    const complex8 *c12   ,
    const real8    *e3    ,
    const real8    *fqr   ,
    const real8    *fqi
)
{
   real8 t[3] = { p4[0]-p3[0],
                  p4[1]-p3[1],
                  p4[2]-p3[2] };

   V3_NORMALIZE(t,t);

   real8 p12[3]     = { p2[0]-p1[0],
                        p2[1]-p1[1],
                        p2[2]-p1[2] };

   real8 temp       = V3_DOT(p12,p12);
   real8 oneoverLp  = 1.0/sqrt(temp);

   real8 proj       = V3_DOT(p12,t);

   real8 diff[3]    = { p12[0] - proj*t[0],
                        p12[1] - proj*t[1],
                        p12[2] - proj*t[2] };

   real8 magdiff    = V3_DOT(diff,diff);
         magdiff    = sqrt(magdiff);

   real8 diffdir[3] = { 0,0,0 };

   real8 eps = 1e-12;

   if (magdiff*oneoverLp > eps)
   {
       diffdir[0] = diff[0]/magdiff;
       diffdir[1] = diff[1]/magdiff;
       diffdir[2] = diff[2]/magdiff;
   }
   else
   {
       // Find another direction

       real8 M[3][3] = { { 1.0-t[0]*t[0],    -t[0]*t[1],    -t[0]*t[2] },
                         {    -t[1]*t[0], 1.0-t[1]*t[1],    -t[1]*t[2] },
                         {    -t[2]*t[0],    -t[2]*t[1], 1.0-t[2]*t[2] } };

       real8 s =   (M[0][0]*M[0][0])
                 + (M[1][0]*M[1][0])
                 + (M[2][0]*M[2][0]);

       if (s > eps)
       {
           s = 1.0/sqrt(s);

           diffdir[0] = s*M[0][0];
           diffdir[1] = s*M[1][0];
           diffdir[2] = s*M[2][0];
       }
       else
       {
           s =   (M[0][1]*M[0][1])
               + (M[1][1]*M[1][1])
               + (M[2][1]*M[2][1]);

           s = 1.0/sqrt(s);

           diffdir[0] = s*M[0][1];
           diffdir[1] = s*M[1][1];
           diffdir[2] = s*M[2][1];
       }
   }

   real8 lmag     = proj*0.5;
   real8 dmag     = lmag*(1.0+1e-8)*sqrt(ecrit/(1.0-ecrit));

   real8 lpmid[3] = { 0.5*(p1[0]+p2[0]),
                      0.5*(p1[1]+p2[1]),
                      0.5*(p1[2]+p2[2]) };

   real8 x1a[3]   = { lpmid[0]-lmag*t[0] - diffdir[0]*dmag,
                      lpmid[1]-lmag*t[1] - diffdir[1]*dmag,
                      lpmid[2]-lmag*t[2] - diffdir[2]*dmag };

   real8 x2a[3]   = { lpmid[0]+lmag*t[0] + diffdir[0]*dmag,
                      lpmid[1]+lmag*t[1] + diffdir[1]*dmag,
                      lpmid[2]+lmag*t[2] + diffdir[2]*dmag };

   real8 x1b[3]   = { lpmid[0]-lmag*t[0] + diffdir[0]*dmag,
                      lpmid[1]-lmag*t[1] + diffdir[1]*dmag,
                      lpmid[2]-lmag*t[2] + diffdir[2]*dmag };

   real8 x2b[3]   = { lpmid[0]+lmag*t[0] - diffdir[0]*dmag,
                      lpmid[1]+lmag*t[1] - diffdir[1]*dmag,
                      lpmid[2]+lmag*t[2] - diffdir[2]*dmag };

   real8 x1p[3]   = { lpmid[0]-lmag*t[0],
                      lpmid[1]-lmag*t[1],
                      lpmid[2]-lmag*t[2] };

   real8 x2p[3]   = { lpmid[0]+lmag*t[0],
                      lpmid[1]+lmag*t[1],
                      lpmid[2]+lmag*t[2] };

   real8 f3a[3]   = { 0,0,0 };
   real8 f4a[3]   = { 0,0,0 };
   real8 f3b[3]   = { 0,0,0 };
   real8 f4b[3]   = { 0,0,0 };
   real8 f3p[3]   = { 0,0,0 };
   real8 f4p[3]   = { 0,0,0 };

   SSF_Aniso_NonParallel ( 0,0,f3a,f4a, x1a,x2a,p3,p4,b12,b34, a,qmax,mtype, c66,c12,e3,fqr,fqi );
   SSF_Aniso_NonParallel ( 0,0,f3b,f4b, x1b,x2b,p3,p4,b12,b34, a,qmax,mtype, c66,c12,e3,fqr,fqi );
   SSF_Aniso_Parallel    (     f3p,f4p, x1p,x2p,p3,p4,b12,b34, a,qmax,mtype, c66,c12,e3,fqr,fqi );

   real8 shape    =  0.25*magdiff/dmag;
   real8 shapea   = (2.0*shape+1.0)*shape;
   real8 shapeb   = (2.0*shape-1.0)*shape;
   real8 shapep   = (1.0+2.0*shape)*(1.0-2.0*shape);

   if (f3)
   {
      f3[0] = (shapea*f3a[0]) + (shapep*f3p[0]) + (shapeb*f3b[0]);
      f3[1] = (shapea*f3a[1]) + (shapep*f3p[1]) + (shapeb*f3b[1]);
      f3[2] = (shapea*f3a[2]) + (shapep*f3p[2]) + (shapeb*f3b[2]);
   }

   if (f4)
   {
      f4[0] = (shapea*f4a[0]) + (shapep*f4p[0]) + (shapeb*f4b[0]);
      f4[1] = (shapea*f4a[1]) + (shapep*f4p[1]) + (shapeb*f4b[1]);
      f4[2] = (shapea*f4a[2]) + (shapep*f4p[2]) + (shapeb*f4b[2]);
   }
}

// SSF_Aniso_Parallel()
//
// This version is the main entry point for the parallel anisotropic forces.
// This routine assumes two parallel segments have been provided and employs
// the shape-function based approximation based on the line directions of
// the two segments.
//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Aniso_Parallel
(
          real8    *f1    ,   ///< force    on node 1 <xyz>   (returned)
          real8    *f2    ,   ///< force    on node 2 <xyz>   (returned)
          real8    *f3    ,   ///< force    on node 3 <xyz>   (returned)
          real8    *f4    ,   ///< force    on node 4 <xyz>   (returned)
    const real8    *p1    ,   ///< position of node 1 <xyz>
    const real8    *p2    ,   ///< position of node 2 <xyz>
    const real8    *p3    ,   ///< position of node 3 <xyz>
    const real8    *p4    ,   ///< position of node 4 <xyz>
    const real8    *b12   ,   ///< burgers vector (p1->p2) <xyz>
    const real8    *b34   ,   ///< burgers vector (p3->p4) <xyz>
    const real8     a     ,   ///< core radius (b)
    const int       qmax  ,   ///< spherical harmonics expansion factor
    const int       mtype ,   ///< material type index
    const real8     ecrit ,   ///< critical angle for parallelism test
    const real8    *c66   ,   ///<
    const complex8 *c12   ,   ///<
    const real8    *e3    ,   ///<
    const real8    *fqr   ,   ///<
    const real8    *fqi       ///<
)
{
    // Define line directions...

    real8 t12[3] = { p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2] };
    real8 t34[3] = { p4[0]-p3[0], p4[1]-p3[1], p4[2]-p3[2] };

    V3_NORMALIZE(t12,t12);
    V3_NORMALIZE(t34,t34);

    real8 c = V3_DOT(t34,t12);

    if (c>0.0)
    {
        // (t12 == t34)

        SSF_Aniso_Parallel_Approx( f1,f2, p3,p4,p1,p2, b34,b12, a,qmax,mtype,ecrit, c66,c12,e3,fqr,fqi );
        SSF_Aniso_Parallel_Approx( f3,f4, p1,p2,p3,p4, b12,b34, a,qmax,mtype,ecrit, c66,c12,e3,fqr,fqi );
    }
    else
    {
        // (t12 == -t34)

        real8 b12n[3] = { -b12[0], -b12[1], -b12[2] };    // b12n = -b12

        SSF_Aniso_Parallel_Approx( f2,f1, p3,p4,p2,p1, b34 ,b12n, a,qmax,mtype,ecrit, c66,c12,e3,fqr,fqi );
        SSF_Aniso_Parallel_Approx( f3,f4, p2,p1,p3,p4, b12n,b34 , a,qmax,mtype,ecrit, c66,c12,e3,fqr,fqi );
    }
}

