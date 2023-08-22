//-----------------------------------------------------------------------------------------------
//
// Module:      SSF_Aniso_NonParallel.cc
//
// Description: Contains functions for calculating anisotropic elastic
//              forces between two non-parallel dislocation sgements
//
//
// H, Hy and Hz are double integrals computed using recurrence.
//
// H_ijk = \int \int (R . c12)^i (R . rotMatrix3)^j dy dz
//                   ----------------------
//                     Ra^k
// Hy_ijk = \int \int y (R . c12)^i (R . rotMatrix3)^j dy dz
//                   ----------------------
//                     Ra^k
// H_ijk = \int \int  z (R . c12)^i (R . rotMatrix3)^j dy dz
//                   ----------------------
//                     Ra^k
// The J integrals are the corresponding single integrals
//
// Jy_ijk = \int (R . c12)^i (R . rotMatrix3)^j dy
//              ----------------------
//                     Ra^k
//
// Jz_ijk = \int (R . c12)^i (R . rotMatrix3)^j dz
//              ----------------------
//                     Ra^k
//
// R = yt + zt' + du
//-----------------------------------------------------------------------------------------------

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "mpi_portability.h"
#include "cuda_portability.h"

#include "Typedefs.h"
#include "Complex_t.h"
#include "V3.h"
#include "Anisotropic.h"

#include "SSF_Aniso.h"

//-----------------------------------------------------------------------------------------------

#ifdef ANISO_QMAX
#define QMAX ANISO_QMAX
#else
#define QMAX 20
#endif

//-----------------------------------------------------------------------------------------------

#define VZERO(a,n)   for (int i=0; (i<(n)); ++i) { (a)[i]=0.0; }

#define MZERO(a,n,m) for (int i=0; (i<(n)); ++i) \
                     for (int j=0; (j<(m)); ++j) { (a)[i][j]=0.0; }

#define M33_ZERO(a)  for (int i=0; (i<3); ++i)   { (a)[i][0]=0.0;  (a)[i][1]=0.0;  (a)[i][2]=0.0;  }
#define M63_ZERO(a)  for (int i=0; (i<6); ++i)   { (a)[i][0]=0.0;  (a)[i][1]=0.0;  (a)[i][2]=0.0;  }

// SSF_Aniso_Stack_Bytes()
//
// Returns an estimate to the size of the stack in bytes for the anisotropic routines.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
size_t SSF_Aniso_Stack_Bytes (const int qmax)
{
    size_t bytes=0;

    bytes += sizeof(real8    [4]              );
    bytes += sizeof(real8    [4]              );
    bytes += sizeof(real8    [4]              );
    bytes += sizeof(real8    [4]              );
    bytes += sizeof(real8    [4]              );
    bytes += sizeof(real8    [4]              );
    bytes += sizeof(real8    [4]              );
    bytes += sizeof(real8    [4]              );
    bytes += sizeof(real8    [4]              );
    bytes += sizeof(real8    [4]              );
    bytes += sizeof(real8    [4]              );

    bytes += sizeof(complex8 [4]              );
    bytes += sizeof(complex8 [4]              );

    bytes += sizeof(real8    [   2*qmax+3]    );
    bytes += sizeof(real8    [   2*qmax+3]    );
    bytes += sizeof(real8    [   2*qmax+3]    );
    bytes += sizeof(real8    [   2*qmax+3][4] );
    bytes += sizeof(real8    [   2*qmax+3][4] );
    bytes += sizeof(real8    [   2*qmax+2]    );
    bytes += sizeof(real8    [   2*qmax+2]    );
    bytes += sizeof(real8    [     qmax+1][4] );
    bytes += sizeof(real8    [     qmax+1][4] );
    bytes += sizeof(real8    [     qmax+2]    );
    bytes += sizeof(real8    [     qmax+2]    );
    bytes += sizeof(real8    [     qmax+2]    );

    bytes += sizeof(real8    [     qmax+1]    );
    bytes += sizeof(real8    [     qmax+1]    );
    bytes += sizeof(real8    [     qmax+1]    );
    bytes += sizeof(real8    [     qmax+1]    );

    bytes += sizeof(real8    [   2*qmax+3]    );
    bytes += sizeof(real8    [   2*qmax+3]    );
    bytes += sizeof(real8    [   2*qmax+3]    );
    bytes += sizeof(real8    [   2*qmax+3]    );

    bytes += sizeof(complex8 [   2*qmax+3]    );
    bytes += sizeof(complex8 [   2*qmax+3]    );
    bytes += sizeof(complex8 [   2*qmax+3]    );
    bytes += sizeof(complex8 [   2*qmax+3]    );

    bytes += sizeof(complex8 [   2*qmax+3]    );
    bytes += sizeof(complex8 [   2*qmax+3]    );
    bytes += sizeof(complex8 [   2*qmax+3]    );
    bytes += sizeof(complex8 [   2*qmax+3][4] );
    bytes += sizeof(complex8 [   2*qmax+3][4] );
    bytes += sizeof(complex8 [   2*qmax+2]    );
    bytes += sizeof(complex8 [   2*qmax+2]    );
    bytes += sizeof(complex8 [2][2*qmax+4]    );
    bytes += sizeof(complex8 [2][2*qmax+4]    );
    bytes += sizeof(complex8 [2][2*qmax+4]    );
    bytes += sizeof(complex8 [   2*qmax+3]    );
    bytes += sizeof(complex8 [   2*qmax+3]    );
    bytes += sizeof(complex8 [   2*qmax+3]    );
    bytes += sizeof(complex8 [   2*qmax+3]    );

    bytes += sizeof(complex8 [2][2*qmax+3]    );
    bytes += sizeof(complex8 [2][2*qmax+3]    );
    bytes += sizeof(complex8 [2][2*qmax+3]    );
    bytes += sizeof(complex8 [2][2*qmax+3]    );
    bytes += sizeof(complex8 [2][2*qmax+3]    );
    bytes += sizeof(complex8 [2][2*qmax+3]    );
    bytes += sizeof(complex8 [2][2*qmax+3]    );
    bytes += sizeof(complex8 [2][2*qmax+3]    );

    bytes += sizeof(complex8 [2][2*qmax+3]    );
    bytes += sizeof(complex8 [2][2*qmax+3]    );
    bytes += sizeof(complex8 [2][2*qmax+3]    );
    bytes += sizeof(complex8 [2][2*qmax+3]    );
    bytes += sizeof(complex8 [2][2*qmax+3]    );
    bytes += sizeof(complex8 [2][2*qmax+3]    );
    bytes += sizeof(complex8 [2][2*qmax+3]    );
    bytes += sizeof(complex8 [2][2*qmax+3]    );

    bytes += sizeof(real8    [18]             );
    bytes += sizeof(real8    [18]             );
    bytes += sizeof(real8    [18]             );
    bytes += sizeof(real8    [18]             );

    return(bytes);
}

__cuda_hdev__
void SSF_Aniso_NonParallel
(
          real8    *f1   ,       ///< force on node 1 <xyz>   (returned)
          real8    *f2   ,       ///< force on node 2 <xyz>   (returned)
          real8    *f3   ,       ///< force on node 3 <xyz>   (returned)
          real8    *f4   ,       ///< force on node 4 <xyz>   (returned)
    const real8    *p1   ,       ///< position of node 1 <xyz>
    const real8    *p2   ,       ///< position of node 2 <xyz>
    const real8    *p3   ,       ///< position of node 3 <xyz>
    const real8    *p4   ,       ///< position of node 4 <xyz>
    const real8    *b12  ,       ///< burgers vector (p1->p2) <xyz>
    const real8    *b34  ,       ///< burgers vector (p3->p4) <xyz>
    const real8     a    ,       ///< core radius (b)
    const int       qmax ,       ///< spherical harmonics expansion factor
    const int       mtype,       ///< material type index
    const real8    *c66  ,       ///< elastic constants matrix (6x6)
    const complex8 *c12  ,       ///< complex elements of aniso rotation matrix (1x3,complex)
    const real8    *e3   ,       ///< real    elements of aniso rotation matrix (1x3,real   )
    const real8    *fqr  ,       ///< Fq table (real)      size=18*(qmax+2)*(qmax+1)
    const real8    *fqi          ///< Fq table (imaginary) size=18*(qmax+2)*(qmax+1)
)
{
    real8   (*ecm)[6]                    = (real8 (*)[6])                 c66;
    real8   (*FqReal)[(QMAX+2)*(QMAX+1)] = (real8 (*)[(QMAX+2)*(QMAX+1)]) fqr;
    real8   (*FqImag)[(QMAX+2)*(QMAX+1)] = (real8 (*)[(QMAX+2)*(QMAX+1)]) fqi;

    real8    Rdt         [4];
    real8    Rdtp        [4];
    real8    R2          [4];
    real8    Ra          [4];
    real8    Ra2inv      [4];
    real8    tmp_v4a     [4];
    real8    tmp_v4b     [4];
    real8    factor1     [4];
    real8    factor2     [4];
    real8    dmceyptd    [4];
    real8    emcdzptd    [4];

    complex8 amcbypgd    [4];
    complex8 bmcazpgd    [4];

    real8    indd        [    2*QMAX+3];
    real8    inde        [    2*QMAX+3];
    real8    indtd       [    2*QMAX+3];
    real8    indemcdzptd [    2*QMAX+3][4];
    real8    inddmceyptd [    2*QMAX+3][4];
    real8    kdmce       [    2*QMAX+2];
    real8    kemcd       [    2*QMAX+2];
    real8    RdtRainv    [      QMAX+1][4];
    real8    RdtpRainv   [      QMAX+1][4];
    real8    factorinv   [      QMAX+2];
    real8    factor      [      QMAX+2];
    real8    newfactor   [      QMAX+2];

    real8    Rainv0      [      QMAX+1];
    real8    Rainv1      [      QMAX+1];
    real8    Rainv2      [      QMAX+1];
    real8    Rainv3      [      QMAX+1];

    real8    Rdsin0      [    2*QMAX+3];
    real8    Rdsin1      [    2*QMAX+3];
    real8    Rdsin2      [    2*QMAX+3];
    real8    Rdsin3      [    2*QMAX+3];

    complex8 Rdpin0      [    2*QMAX+3];
    complex8 Rdpin1      [    2*QMAX+3];
    complex8 Rdpin2      [    2*QMAX+3];
    complex8 Rdpin3      [    2*QMAX+3];

    complex8 inda        [    2*QMAX+3];
    complex8 indb        [    2*QMAX+3];
    complex8 indgd       [    2*QMAX+3];
    complex8 indbmcazpgd [    2*QMAX+3][4];
    complex8 indamcbypgd [    2*QMAX+3][4];
    complex8 jamcb       [    2*QMAX+2];
    complex8 jbmca       [    2*QMAX+2];
    complex8 Ht0         [ 2][2*QMAX+4];
    complex8 Ht1         [ 2][2*QMAX+4];
    complex8 Ht2         [ 2][2*QMAX+4];
    complex8 Hyt0        [    2*QMAX+3];
    complex8 Hyt1        [    2*QMAX+3];
    complex8 Hzt0        [    2*QMAX+3];
    complex8 Hzt1        [    2*QMAX+3];

    complex8 Jyt10       [ 2][2*QMAX+3];   MZERO(Jyt10,2,2*QMAX+3);
    complex8 Jyt20       [ 2][2*QMAX+3];   MZERO(Jyt20,2,2*QMAX+3);
    complex8 Jyt11       [ 2][2*QMAX+3];   MZERO(Jyt11,2,2*QMAX+3);
    complex8 Jyt21       [ 2][2*QMAX+3];   MZERO(Jyt21,2,2*QMAX+3);
    complex8 Jyt12       [ 2][2*QMAX+3];   MZERO(Jyt12,2,2*QMAX+3);
    complex8 Jyt22       [ 2][2*QMAX+3];   MZERO(Jyt22,2,2*QMAX+3);
    complex8 Jyt13       [ 2][2*QMAX+3];   MZERO(Jyt13,2,2*QMAX+3);
    complex8 Jyt23       [ 2][2*QMAX+3];   MZERO(Jyt23,2,2*QMAX+3);

    complex8 Jzt10       [ 2][2*QMAX+3];   MZERO(Jzt10,2,2*QMAX+3);
    complex8 Jzt20       [ 2][2*QMAX+3];   MZERO(Jzt20,2,2*QMAX+3);
    complex8 Jzt11       [ 2][2*QMAX+3];   MZERO(Jzt11,2,2*QMAX+3);
    complex8 Jzt21       [ 2][2*QMAX+3];   MZERO(Jzt21,2,2*QMAX+3);
    complex8 Jzt12       [ 2][2*QMAX+3];   MZERO(Jzt12,2,2*QMAX+3);
    complex8 Jzt22       [ 2][2*QMAX+3];   MZERO(Jzt22,2,2*QMAX+3);
    complex8 Jzt13       [ 2][2*QMAX+3];   MZERO(Jzt13,2,2*QMAX+3);
    complex8 Jzt23       [ 2][2*QMAX+3];   MZERO(Jzt23,2,2*QMAX+3);

    real8    G1[18];                       VZERO(G1,18);
    real8    G2[18];                       VZERO(G2,18);
    real8    G3[18];                       VZERO(G3,18);
    real8    G4[18];                       VZERO(G4,18);

    // Define line directions

    real8 dX = p2[0] - p1[0];
    real8 dY = p2[1] - p1[1];
    real8 dZ = p2[2] - p1[2];

    real8 invL    = 1.0 / sqrt(dX*dX + dY*dY + dZ*dZ);
    real8 coeff12 = invL / (4*M_PI*M_PI);

    real8 tp[3] = { dX * invL,
                    dY * invL,
                    dZ * invL };

          dX = p4[0] - p3[0];
          dY = p4[1] - p3[1];
          dZ = p4[2] - p3[2];

          invL = 1.0 / sqrt(dX*dX + dY*dY + dZ*dZ);

    real8 coeff34 = invL / (4*M_PI*M_PI);

    real8 t[3] = { dX * invL,
                   dY * invL,
                   dZ * invL };

    // Compute double and single integrals via recurrence

    real8 c         = V3_DOT(tp, t);
    real8 c2        = c * c;
    real8 onemc2    = 1.0 - c2;
    real8 onemc2inv = 1.0/onemc2;

    real8 txtp[3]   = { t[1]*tp[2] - t[2]*tp[1],
                        t[2]*tp[0] - t[0]*tp[2],
                        t[0]*tp[1] - t[1]*tp[0] };

    complex8 alpha   = V3_DOT(t, c12);
    complex8 beta    = V3_DOT(tp,c12);
    real8    delta   = V3_DOT(t, e3);
    real8    epsilon = V3_DOT(tp,e3);

    complex8 amcb    = alpha   - (c * beta   );
    complex8 bmca    = beta    - (c * alpha  );

    real8    dmce    = delta   - (c * epsilon);
    real8    emcd    = epsilon - (c * delta  );

    real8    R[4][3] = { { (p3[0]-p1[0]), (p3[1]-p1[1]), (p3[2]-p1[2]) },
                         { (p3[0]-p2[0]), (p3[1]-p2[1]), (p3[2]-p2[2]) },
                         { (p4[0]-p1[0]), (p4[1]-p1[1]), (p4[2]-p1[2]) },
                         { (p4[0]-p2[0]), (p4[1]-p2[1]), (p4[2]-p2[2]) } };

    real8  d = V3_DOT(R[0], txtp) * onemc2inv;

    real8  temp1[2] = { V3_DOT(R[0],t ),
                        V3_DOT(R[3],t ) };

    real8  temp2[2] = { V3_DOT(R[0],tp),
                        V3_DOT(R[3],tp) };

    real8  y0 = (temp1[0] - c * temp2[0]) * onemc2inv;
    real8  y1 = (temp1[1] - c * temp2[1]) * onemc2inv;

    real8  z0 = (temp2[0] - c * temp1[0]) * onemc2inv;
    real8  z1 = (temp2[1] - c * temp1[1]) * onemc2inv;

    real8  yin[4] = { y0, y0, y1, y1 };
    real8  zin[4] = { z0, z1, z0, z1 };

    Rdpin0[0] = 1.0;
    Rdsin0[0] = 1.0;
    Rdpin0[1] = V3_DOT(R[0], c12);
    Rdsin0[1] = V3_DOT(R[0], e3 );

    Rdpin1[0] = 1.0;
    Rdsin1[0] = 1.0;
    Rdpin1[1] = V3_DOT(R[1], c12);
    Rdsin1[1] = V3_DOT(R[1], e3 );

    Rdpin2[0] = 1.0;
    Rdsin2[0] = 1.0;
    Rdpin2[1] = V3_DOT(R[2], c12);
    Rdsin2[1] = V3_DOT(R[2], e3 );

    Rdpin3[0] = 1.0;
    Rdsin3[0] = 1.0;
    Rdpin3[1] = V3_DOT(R[3], c12);
    Rdsin3[1] = V3_DOT(R[3], e3 );

    complex8 gammad = 0.5 * (Rdpin0[1] + Rdpin3[1] - (y0+y1) * alpha - beta    * (z0+z1));
    real8    thetad = 0.5 * (Rdsin0[1] + Rdsin3[1] - (y0+y1) * delta - epsilon * (z0+z1));

    for (int i=0; (i<4); ++i)
    {
        bmcazpgd[i] = (bmca * zin[i] + gammad);
        emcdzptd[i] = (emcd * zin[i] + thetad);
        amcbypgd[i] = (amcb * yin[i] + gammad);
        dmceyptd[i] = (dmce * yin[i] + thetad);

        R2[i] = V3_DOT(R[i], R[i]);
    }

    real8 a2      = a * a;
    real8 a2_d2   = a2 + d * d * onemc2;

    real8 cdenom  = onemc2 * a2_d2;
    real8 a2d2inv = 1.0 / a2_d2;
    real8 denom   = 1.0 / sqrt(cdenom);
          cdenom  = (1.0+c) * denom;

    for (int i=0; (i<4); ++i) { Rdt    [i] = yin[i] + zin[i]*c;  }
    for (int i=0; (i<4); ++i) { Rdtp   [i] = zin[i] + yin[i]*c;  }
    for (int i=0; (i<4); ++i) { tmp_v4a[i] = R2 [i] + a2;        }
    for (int i=0; (i<4); ++i) { Ra     [i] = sqrt(tmp_v4a[i]);   }
    for (int i=0; (i<4); ++i) { tmp_v4a[i] = Ra[i]+Rdt [i]; }
    for (int i=0; (i<4); ++i) { tmp_v4b[i] = Ra[i]+Rdtp[i]; }

    Rainv0[0] = 1.0 / Ra[0];
    Rainv1[0] = 1.0 / Ra[1];
    Rainv2[0] = 1.0 / Ra[2];
    Rainv3[0] = 1.0 / Ra[3];

    Ra2inv[0] = Rainv0[0] * Rainv0[0];
    Ra2inv[1] = Rainv1[0] * Rainv1[0];
    Ra2inv[2] = Rainv2[0] * Rainv2[0];
    Ra2inv[3] = Rainv3[0] * Rainv3[0];

    RdtRainv [0][0] = Rainv0[0] * Rdt [0];
    RdtpRainv[0][0] = Rainv0[0] * Rdtp[0];

    RdtRainv [0][1] = Rainv1[0] * Rdt [1];
    RdtpRainv[0][1] = Rainv1[0] * Rdtp[1];

    RdtRainv [0][2] = Rainv2[0] * Rdt [2];
    RdtpRainv[0][2] = Rainv2[0] * Rdtp[2];

    RdtRainv [0][3] = Rainv3[0] * Rdt [3];
    RdtpRainv[0][3] = Rainv3[0] * Rdtp[3];

    // Initialize Jy, Jz, Jysum, Jzsum

    Jyt10[0][0] = 0.0;  Jzt10[0][0] = 0.0;
    Jyt11[0][0] = 0.0;  Jzt11[0][0] = 0.0;
    Jyt12[0][0] = 0.0;  Jzt12[0][0] = 0.0;
    Jyt13[0][0] = 0.0;  Jzt13[0][0] = 0.0;

    Jyt20[0][0] = 0.0;  Jzt20[0][0] = 0.0;
    Jyt21[0][0] = 0.0;  Jzt21[0][0] = 0.0;
    Jyt22[0][0] = 0.0;  Jzt22[0][0] = 0.0;
    Jyt23[0][0] = 0.0;  Jzt23[0][0] = 0.0;

    Jyt10[0][1] = log(tmp_v4a[0]);
    Jzt10[0][1] = log(tmp_v4b[0]);

    Jyt11[0][1] = log(tmp_v4a[1]);
    Jzt11[0][1] = log(tmp_v4b[1]);

    Jyt12[0][1] = log(tmp_v4a[2]);
    Jzt12[0][1] = log(tmp_v4b[2]);

    Jyt13[0][1] = log(tmp_v4a[3]);
    Jzt13[0][1] = log(tmp_v4b[3]);

    Jyt20[0][2] =   alpha * Ra[0] + bmcazpgd[0] * Jyt10[0][1];
    Jyt20[0][1] =   delta * Ra[0] + emcdzptd[0] * Jyt10[0][1];
    Jzt20[0][2] =    beta * Ra[0] + amcbypgd[0] * Jzt10[0][1];
    Jzt20[0][1] = epsilon * Ra[0] + dmceyptd[0] * Jzt10[0][1];

    Jyt21[0][2] =   alpha * Ra[1] + bmcazpgd[1] * Jyt11[0][1];
    Jyt21[0][1] =   delta * Ra[1] + emcdzptd[1] * Jyt11[0][1];
    Jzt21[0][2] =    beta * Ra[1] + amcbypgd[1] * Jzt11[0][1];
    Jzt21[0][1] = epsilon * Ra[1] + dmceyptd[1] * Jzt11[0][1];

    Jyt22[0][2] =   alpha * Ra[2] + bmcazpgd[2] * Jyt12[0][1];
    Jyt22[0][1] =   delta * Ra[2] + emcdzptd[2] * Jyt12[0][1];
    Jzt22[0][2] =    beta * Ra[2] + amcbypgd[2] * Jzt12[0][1];
    Jzt22[0][1] = epsilon * Ra[2] + dmceyptd[2] * Jzt12[0][1];

    Jyt23[0][2] =   alpha * Ra[3] + bmcazpgd[3] * Jyt13[0][1];
    Jyt23[0][1] =   delta * Ra[3] + emcdzptd[3] * Jyt13[0][1];
    Jzt23[0][2] =    beta * Ra[3] + amcbypgd[3] * Jzt13[0][1];
    Jzt23[0][1] = epsilon * Ra[3] + dmceyptd[3] * Jzt13[0][1];

    for (int i=0; (i<4); ++i) { tmp_v4a[i] = R2[i] + a2 - Rdt [i]*Rdt [i]; }
    for (int i=0; (i<4); ++i) { tmp_v4b[i] = R2[i] + a2 - Rdtp[i]*Rdtp[i]; }

    complex8 Jytsum = Jyt10[0][1] - Jyt11[0][1] - Jyt12[0][1] + Jyt13[0][1];
    complex8 Jztsum = Jzt10[0][1] - Jzt11[0][1] - Jzt12[0][1] + Jzt13[0][1];

    // Initialize Hysumint, Hzsumint, Hsumint
    //
    // Preliminary double integrals involving unique terms

    real8 temp3a = atan((Ra[1] + y0 + z1) * cdenom);
    real8 temp3b = atan((Ra[0] + y0 + z0) * cdenom);
    real8 temp3c = atan((Ra[3] + y1 + z1) * cdenom);
    real8 temp3d = atan((Ra[2] + y1 + z0) * cdenom);

             Ht0[1][1]    = temp3a - temp3b - temp3c + temp3d;

             Ht0[1][1]    = Ht0[1][1] * 2.0 * denom;

    real8    Rasum        = Ra[0] - Ra[1] - Ra[2] + Ra[3];

    real8    adW_003      = a2_d2 * creal(Ht0[1][1]);

    real8    commonW223   = (c * Rasum - adW_003) * onemc2inv;

    complex8 I113 = (c * adW_003 - Rasum) * onemc2inv;
    complex8 I203 = z0*Jyt10[0][1] - z1*Jyt11[0][1] - z0*Jyt12[0][1] + z1*Jyt13[0][1] + commonW223;
    complex8 I023 = y0*Jzt10[0][1] - y0*Jzt11[0][1] - y1*Jzt12[0][1] + y1*Jzt13[0][1] + commonW223;
    complex8 I103 = (c * Jytsum - Jztsum) * onemc2inv;
    complex8 I013 = (c * Jztsum - Jytsum) * onemc2inv;

    Hyt0[1]   = I103;
    Hzt0[1]   = I013;

    Hyt1[1]   = delta*I203    + epsilon*I113    + thetad*Hyt0[1];
    Hyt1[2]   = alpha*I203    +    beta*I113    + gammad*Hyt0[1];

    Hzt1[1]   = delta*I113    + epsilon*I023    + thetad*Hzt0[1];
    Hzt1[2]   = alpha*I113    +    beta*I023    + gammad*Hzt0[1];

    Ht1[1][1] = delta*I103    + epsilon*I013    + thetad*Ht0[1][1];
    Ht2[1][1] = delta*Hyt1[1] + epsilon*Hzt1[1] + thetad*Ht1[1][1];
    Ht1[1][2] = alpha*I103    + beta   *I013    + gammad*Ht0[1][1];
    Ht2[1][2] = alpha*Hyt1[1] + beta   *Hzt1[1] + gammad*Ht1[1][1];
    Ht2[1][3] = alpha*Hyt1[2] + beta   *Hzt1[2] + gammad*Ht1[1][2];

    real8 reHzt = creal(Hzt1[1]);
    real8 reHyt = creal(Hyt1[1]);
    real8 reHt  = creal(Ht1[1][1]);

    real8 imHzt = cimag(Hzt1[1]);
    real8 imHyt = cimag(Hyt1[1]);
    real8 imHt  = cimag(Ht1[1][1]);

    real8 reHz1 = reHzt - z0*reHt;
    real8 reHz2 = reHzt - z1*reHt;
    real8 reHy1 = reHyt - y0*reHt;
    real8 reHy2 = reHyt - y1*reHt;

    real8 imHz1 = imHzt - z0*imHt;
    real8 imHz2 = imHzt - z1*imHt;
    real8 imHy1 = imHyt - y0*imHt;
    real8 imHy2 = imHyt - y1*imHt;

    for (int i=0; (i<18); ++i)
    {
        real8 a = FqReal[i][0];
        real8 b = FqImag[i][0];

        G1[i] += b*imHz2 - a*reHz2;
        G2[i] += a*reHz1 - b*imHz1;
        G3[i] += b*imHy2 - a*reHy2;
        G4[i] += a*reHy1 - b*imHy1;
    }

    reHzt = creal(Hzt1[2]   );
    reHyt = creal(Hyt1[2]   );
    reHt  = creal(Ht1 [1][2]);

    imHzt = cimag(Hzt1[2]   );
    imHyt = cimag(Hyt1[2]   );
    imHt  = cimag(Ht1 [1][2]);

    reHz1 = reHzt - z0*reHt;
    reHz2 = reHzt - z1*reHt;
    reHy1 = reHyt - y0*reHt;
    reHy2 = reHyt - y1*reHt;

    imHz1 = imHzt - z0*imHt;
    imHz2 = imHzt - z1*imHt;
    imHy1 = imHyt - y0*imHt;
    imHy2 = imHyt - y1*imHt;

    for (int i=0; (i<18); ++i)
    {
        real8 a = FqReal[i][1];
        real8 b = FqImag[i][1];

        G1[i] += b*imHz2 - a*reHz2;
        G2[i] += a*reHz1 - b*imHz1;
        G3[i] += b*imHy2 - a*reHy2;
        G4[i] += a*reHy1 - b*imHy1;
    }

    for (int i=2; i<(2*qmax+3); ++i)
    {
        Rdpin0[i] = Rdpin0[i-1] * Rdpin0[1];
        Rdpin1[i] = Rdpin1[i-1] * Rdpin1[1];
        Rdpin2[i] = Rdpin2[i-1] * Rdpin2[1];
        Rdpin3[i] = Rdpin3[i-1] * Rdpin3[1];

        Rdsin0[i] = Rdsin0[i-1] * Rdsin0[1];
        Rdsin1[i] = Rdsin1[i-1] * Rdsin1[1];
        Rdsin2[i] = Rdsin2[i-1] * Rdsin2[1];
        Rdsin3[i] = Rdsin3[i-1] * Rdsin3[1];
    }

    for (int i=0; (i<4); i++)
    {
        factor1[i] = 1.0 / tmp_v4a[i];
        factor2[i] = 1.0 / tmp_v4b[i];
    }

    factorinv[0]=-1;
    for (int i=1; i<(qmax+2); i++) { factorinv[i] = factorinv[i-1]+2;      }
    for (int i=0; i<(qmax+2); i++) { factor   [i] = 1.0 / factorinv[i];    }
    for (int i=0; i<(qmax+2); i++) { newfactor[i] = factor[i] * onemc2inv; }

    for (int i=1,j=0; i<(qmax+1); i++,j++)
    {
        Rainv0   [i]    = Rainv0[j] * Ra2inv[0];
        RdtRainv [i][0] = Rainv0[i] * Rdt   [0];
        RdtpRainv[i][0] = Rainv0[i] * Rdtp  [0];
    }

    for (int i=1,j=0; i<(qmax+1); i++,j++)
    {
        Rainv1   [i]    = Rainv1[j] * Ra2inv[1];
        RdtRainv [i][1] = Rainv1[i] * Rdt   [1];
        RdtpRainv[i][1] = Rainv1[i] * Rdtp  [1];
    }

    for (int i=1,j=0; i<(qmax+1); i++,j++)
    {
        Rainv2   [i]    = Rainv2[j] * Ra2inv[2];
        RdtRainv [i][2] = Rainv2[i] * Rdt   [2];
        RdtpRainv[i][2] = Rainv2[i] * Rdtp  [2];
    }

    for (int i=1,j=0; i<(qmax+1); i++,j++)
    {
        Rainv3   [i]    = Rainv3[j] * Ra2inv[3];
        RdtRainv [i][3] = Rainv3[i] * Rdt   [3];
        RdtpRainv[i][3] = Rainv3[i] * Rdtp  [3];
    }

    cuda_syncthreads();

    inda [0] = 0.0;
    indb [0] = 0.0;
    indd [0] = 0.0;
    inde [0] = 0.0;
    indgd[0] = 0.0;
    indtd[0] = 0.0;

    for (int i=1; i<(2*qmax+3); i++)
    {
        inda [i] = inda [i-1] + alpha  ;
        indb [i] = indb [i-1] + beta   ;
        indd [i] = indd [i-1] + delta  ;
        inde [i] = inde [i-1] + epsilon;
        indgd[i] = indgd[i-1] + gammad ;
        indtd[i] = indtd[i-1] + thetad ;
    }

    for (int i=0; (i<4); i++)
    {
        indemcdzptd[0][i] = 0.0;
        indbmcazpgd[0][i] = 0.0;
        indamcbypgd[0][i] = 0.0;
        inddmceyptd[0][i] = 0.0;
    }

    for (int m=1; m<(2*qmax+3); m++)
    {
        for (int i=0; (i<4); i++)
        {
            indemcdzptd[m][i] = indemcdzptd[m-1][i]+emcdzptd[i];
            indbmcazpgd[m][i] = indbmcazpgd[m-1][i]+bmcazpgd[i];
            indamcbypgd[m][i] = indamcbypgd[m-1][i]+amcbypgd[i];
            inddmceyptd[m][i] = inddmceyptd[m-1][i]+dmceyptd[i];
        }
    }

    jamcb[0] = 0.0;
    jbmca[0] = 0.0;
    kdmce[0] = 0.0;
    kemcd[0] = 0.0;

    for (int i=1,j=0; i<(2*qmax+2); i++,j++)
    {
        jamcb[i] = jamcb[j] + amcb;
        jbmca[i] = jbmca[j] + bmca;
        kdmce[i] = kdmce[j] + dmce;
        kemcd[i] = kemcd[j] + emcd;
    }

    cuda_syncthreads();

    // Start recursion

    int ind   = 1;
    int indm1 = 0;
    int local = 2;
    for (int i=2; i<(qmax+2); i++)
    {
        real8 onemfactor = 1.0 - (2*i-2) * factor[i-1];

        // Calculate J00k

        complex8 afactor = factor[i-1] * alpha  ;
        complex8 bfactor = factor[i-1] * beta   ;
        real8    dfactor = factor[i-1] * delta  ;
        real8    efactor = factor[i-1] * epsilon;

        cuda_syncthreads();

        // m=0
        {
            complex8 Jyt00[2*QMAX+3]; Jyt00[0] = 0.0;
            complex8 Jzt00[2*QMAX+3]; Jzt00[0] = 0.0;

            for (int j=0,k=2*i-3; j<(2*i-2); j++,k-- ) { Jyt00[j+1] = factor1[0] * (factor[i-1]*(Rdpin0[j]*Rdsin0[k]*RdtRainv [i-2][0] + indbmcazpgd[j][0]*Jyt10[indm1][j] + indemcdzptd[k][0]*Jyt10[indm1][j+1]) + onemfactor*Jyt20[indm1][j+1]); }
            for (int j=0,k=2*i-3; j<(2*i-2); j++,k-- ) { Jzt00[j+1] = factor2[0] * (factor[i-1]*(Rdpin0[j]*Rdsin0[k]*RdtpRainv[i-2][0] + indamcbypgd[j][0]*Jzt10[indm1][j] + inddmceyptd[k][0]*Jzt10[indm1][j+1]) + onemfactor*Jzt20[indm1][j+1]); }

            int   j     = 2*i-3;
            real8 tmpre = Rdsin0[j]*Rainv0[i-2];

            Jyt10[ind][1] = dfactor * (indd[j]*Jyt10[indm1][1]-tmpre) + emcdzptd[0]*Jyt00[1];
            Jzt10[ind][1] = efactor * (inde[j]*Jzt10[indm1][1]-tmpre) + dmceyptd[0]*Jzt00[1];

            for (int j=0,k=2*i-3; j<(2*i-2); j++,k--)
            {
                complex8 cpx = Rdpin0[j]*Rdsin0[k]*Rainv0[i-2];

                Jyt10[ind][j+2] = afactor * (indd[k]*Jyt10[indm1][j+1] + inda[j]*Jyt10[indm1][j] - cpx) + bmcazpgd[0]*Jyt00[j+1];
                Jzt10[ind][j+2] = bfactor * (inde[k]*Jzt10[indm1][j+1] + indb[j]*Jzt10[indm1][j] - cpx) + amcbypgd[0]*Jzt00[j+1];
            }
        }

        // m=1
        {
            complex8 Jyt01[2*QMAX+3]; Jyt01[0]=0.0;
            complex8 Jzt01[2*QMAX+3]; Jzt01[0]=0.0;

            for (int j=0,k=2*i-3; j<(2*i-2); j++,k-- ) { Jyt01[j+1] = factor1[1] * (factor[i-1]*(Rdpin1[j]*Rdsin1[k]*RdtRainv [i-2][1] + indbmcazpgd[j][1]*Jyt11[indm1][j] + indemcdzptd[k][1]*Jyt11[indm1][j+1]) + onemfactor*Jyt21[indm1][j+1]); }
            for (int j=0,k=2*i-3; j<(2*i-2); j++,k-- ) { Jzt01[j+1] = factor2[1] * (factor[i-1]*(Rdpin1[j]*Rdsin1[k]*RdtpRainv[i-2][1] + indamcbypgd[j][1]*Jzt11[indm1][j] + inddmceyptd[k][1]*Jzt11[indm1][j+1]) + onemfactor*Jzt21[indm1][j+1]); }

            int   j     = 2*i-3;
            real8 tmpre = Rdsin1[j]*Rainv1[i-2];

            Jyt11[ind][1] = dfactor * (indd[j]*Jyt11[indm1][1]-tmpre) + emcdzptd[1]*Jyt01[1];
            Jzt11[ind][1] = efactor * (inde[j]*Jzt11[indm1][1]-tmpre) + dmceyptd[1]*Jzt01[1];

            for (int j=0,k=2*i-3; j<(2*i-2); j++,k--)
            {
                complex8 cpx = Rdpin1[j]*Rdsin1[k]*Rainv1[i-2];

                Jyt11[ind][j+2] = afactor * (indd[k]*Jyt11[indm1][j+1] + inda[j]*Jyt11[indm1][j] - cpx) + bmcazpgd[1]*Jyt01[j+1];
                Jzt11[ind][j+2] = bfactor * (inde[k]*Jzt11[indm1][j+1] + indb[j]*Jzt11[indm1][j] - cpx) + amcbypgd[1]*Jzt01[j+1];
            }
        }

        // m=2
        {
            complex8 Jyt02[2*QMAX+3]; Jyt02[0]=0.0;
            complex8 Jzt02[2*QMAX+3]; Jzt02[0]=0.0;

            for (int j=0,k=2*i-3; j<(2*i-2); j++,k-- ) { Jyt02[j+1] = factor1[2] * (factor[i-1]*(Rdpin2[j]*Rdsin2[k]*RdtRainv [i-2][2] + indbmcazpgd[j][2]*Jyt12[indm1][j] + indemcdzptd[k][2]*Jyt12[indm1][j+1]) + onemfactor*Jyt22[indm1][j+1]); }
            for (int j=0,k=2*i-3; j<(2*i-2); j++,k-- ) { Jzt02[j+1] = factor2[2] * (factor[i-1]*(Rdpin2[j]*Rdsin2[k]*RdtpRainv[i-2][2] + indamcbypgd[j][2]*Jzt12[indm1][j] + inddmceyptd[k][2]*Jzt12[indm1][j+1]) + onemfactor*Jzt22[indm1][j+1]); }

            int   j     = 2*i-3;
            real8 tmpre = Rdsin2[j]*Rainv2[i-2];

            Jyt12[ind][1] = dfactor * (indd[j]*Jyt12[indm1][1]-tmpre) + emcdzptd[2]*Jyt02[1];
            Jzt12[ind][1] = efactor * (inde[j]*Jzt12[indm1][1]-tmpre) + dmceyptd[2]*Jzt02[1];

            for (int j=0,k=2*i-3; j<(2*i-2); j++,k--)
            {
                complex8 cpx = Rdpin2[j]*Rdsin2[k]*Rainv2[i-2];

                Jyt12[ind][j+2] = afactor * (indd[k]*Jyt12[indm1][j+1] + inda[j]*Jyt12[indm1][j] - cpx) + bmcazpgd[2]*Jyt02[j+1];
                Jzt12[ind][j+2] = bfactor * (inde[k]*Jzt12[indm1][j+1] + indb[j]*Jzt12[indm1][j] - cpx) + amcbypgd[2]*Jzt02[j+1];
            }
        }

        // m=3
        {
            complex8 Jyt03[2*QMAX+3]; Jyt03[0]=0.0;
            complex8 Jzt03[2*QMAX+3]; Jzt03[0]=0.0;

            for (int j=0,k=2*i-3; j<(2*i-2); j++,k-- ) { Jyt03[j+1] = factor1[3] * (factor[i-1]*(Rdpin3[j]*Rdsin3[k]*RdtRainv [i-2][3] + indbmcazpgd[j][3]*Jyt13[indm1][j] + indemcdzptd[k][3]*Jyt13[indm1][j+1]) + onemfactor*Jyt23[indm1][j+1]); }
            for (int j=0,k=2*i-3; j<(2*i-2); j++,k-- ) { Jzt03[j+1] = factor2[3] * (factor[i-1]*(Rdpin3[j]*Rdsin3[k]*RdtpRainv[i-2][3] + indamcbypgd[j][3]*Jzt13[indm1][j] + inddmceyptd[k][3]*Jzt13[indm1][j+1]) + onemfactor*Jzt23[indm1][j+1]); }

            int   j     = 2*i-3;
            real8 tmpre = Rdsin3[j]*Rainv3[i-2];

            Jyt13[ind][1] = dfactor * (indd[j]*Jyt13[indm1][1]-tmpre) + emcdzptd[3]*Jyt03[1];
            Jzt13[ind][1] = efactor * (inde[j]*Jzt13[indm1][1]-tmpre) + dmceyptd[3]*Jzt03[1];

            for (int j=0,k=2*i-3; j<(2*i-2); j++,k--)
            {
                complex8 cpx = Rdpin3[j]*Rdsin3[k]*Rainv3[i-2];

                Jyt13[ind][j+2] = afactor * (indd[k]*Jyt13[indm1][j+1] + inda[j]*Jyt13[indm1][j] - cpx) + bmcazpgd[3]*Jyt03[j+1];
                Jzt13[ind][j+2] = bfactor * (inde[k]*Jzt13[indm1][j+1] + indb[j]*Jzt13[indm1][j] - cpx) + amcbypgd[3]*Jzt03[j+1];
            }
        }

        cuda_syncthreads();

        // m=0
        {
            int   j     = 2*i-2;
            real8 tmpre = Rdsin0[j]*Rainv0[i-2];

            Jyt20[ind][1] = dfactor * (indd[j]*Jyt20[indm1][1]-tmpre) + emcdzptd[0]*Jyt10[ind][1];
            Jzt20[ind][1] = efactor * (inde[j]*Jzt20[indm1][1]-tmpre) + dmceyptd[0]*Jzt10[ind][1];

            for (int j=0,k=2*i-2; j<(2*i-1); j++,k--)
            {
                complex8 cpx = Rdpin0[j]*Rdsin0[k]*Rainv0[i-2];

                Jyt20[ind][j+2] = afactor * (indd[k]*Jyt20[indm1][j+1] + inda[j]*Jyt20[indm1][j] - cpx) + bmcazpgd[0]*Jyt10[ind][j+1];
                Jzt20[ind][j+2] = bfactor * (inde[k]*Jzt20[indm1][j+1] + indb[j]*Jzt20[indm1][j] - cpx) + amcbypgd[0]*Jzt10[ind][j+1];
            }
        }

        // m=1
        {
            int   j     = 2*i-2;
            real8 tmpre = Rdsin1[j]*Rainv1[i-2];

            Jyt21[ind][1] = dfactor * (indd[j]*Jyt21[indm1][1]-tmpre) + emcdzptd[1]*Jyt11[ind][1];
            Jzt21[ind][1] = efactor * (inde[j]*Jzt21[indm1][1]-tmpre) + dmceyptd[1]*Jzt11[ind][1];

            for (int j=0,k=2*i-2; j<(2*i-1); j++,k--)
            {
                complex8 cpx = Rdpin1[j]*Rdsin1[k]*Rainv1[i-2];

                Jyt21[ind][j+2] = afactor * (indd[k]*Jyt21[indm1][j+1] + inda[j]*Jyt21[indm1][j] - cpx) + bmcazpgd[1]*Jyt11[ind][j+1];
                Jzt21[ind][j+2] = bfactor * (inde[k]*Jzt21[indm1][j+1] + indb[j]*Jzt21[indm1][j] - cpx) + amcbypgd[1]*Jzt11[ind][j+1];
            }
        }

        // m=2
        {
            int   j     = 2*i-2;
            real8 tmpre = Rdsin2[j]*Rainv2[i-2];

            Jyt22[ind][1] = dfactor * (indd[j]*Jyt22[indm1][1]-tmpre) + emcdzptd[2]*Jyt12[ind][1];
            Jzt22[ind][1] = efactor * (inde[j]*Jzt22[indm1][1]-tmpre) + dmceyptd[2]*Jzt12[ind][1];

            for (int j=0,k=2*i-2; j<(2*i-1); j++,k--)
            {
                complex8 cpx = Rdpin2[j]*Rdsin2[k]*Rainv2[i-2];

                Jyt22[ind][j+2] = afactor * (indd[k]*Jyt22[indm1][j+1] + inda[j]*Jyt22[indm1][j] - cpx) + bmcazpgd[2]*Jyt12[ind][j+1];
                Jzt22[ind][j+2] = bfactor * (inde[k]*Jzt22[indm1][j+1] + indb[j]*Jzt22[indm1][j] - cpx) + amcbypgd[2]*Jzt12[ind][j+1];
            }
        }

        // m=3
        {
            int   j     = 2*i-2;
            real8 tmpre = Rdsin3[j]*Rainv3[i-2];

            Jyt23[ind][1] = dfactor * (indd[j]*Jyt23[indm1][1]-tmpre) + emcdzptd[3]*Jyt13[ind][1];
            Jzt23[ind][1] = efactor * (inde[j]*Jzt23[indm1][1]-tmpre) + dmceyptd[3]*Jzt13[ind][1];

            for (int j=0,k=2*i-2; j<(2*i-1); j++,k--)
            {
                complex8 cpx = Rdpin3[j]*Rdsin3[k]*Rainv3[i-2];

                Jyt23[ind][j+2] = afactor * (indd[k]*Jyt23[indm1][j+1] + inda[j]*Jyt23[indm1][j] - cpx) + bmcazpgd[3]*Jyt13[ind][j+1];
                Jzt23[ind][j+2] = bfactor * (inde[k]*Jzt23[indm1][j+1] + indb[j]*Jzt23[indm1][j] - cpx) + amcbypgd[3]*Jzt13[ind][j+1];
            }
        }

        cuda_syncthreads();

        // Recursive double integrals
        //
        // First construct Hyt and Hzt

        indm1 =   ind;
        ind   = 1-ind;

        // m=1
        for (int j=1,k=2*i-2; j<(2*i); j++,k--)
        {
            complex8 tmp1 = Jzt10[indm1][j] - Jzt11[indm1][j] - Jzt12[indm1][j] + Jzt13[indm1][j];
            complex8 tmp2 = Jyt10[indm1][j] - Jyt11[indm1][j] - Jyt12[indm1][j] + Jyt13[indm1][j];

            Hyt0[j] = newfactor[i] * (jamcb[j-1]*Ht1[indm1][j-1] + kdmce[k]*Ht1[indm1][j] - tmp1+c*tmp2);
            Hzt0[j] = newfactor[i] * (jbmca[j-1]*Ht1[indm1][j-1] + kemcd[k]*Ht1[indm1][j] - tmp2+c*tmp1);
        }

        // m=2
        for (int j=1,k=2*i-1; j<(2*i+1); j++,k--)
        {
            complex8 tmp1 = Jzt20[indm1][j] - Jzt21[indm1][j] - Jzt22[indm1][j] + Jzt23[indm1][j];
            complex8 tmp2 = Jyt20[indm1][j] - Jyt21[indm1][j] - Jyt22[indm1][j] + Jyt23[indm1][j];

            Hyt1[j] = newfactor[i] * (jamcb[j-1]*Ht2[indm1][j-1] + kdmce[k]*Ht2[indm1][j] - tmp1+c*tmp2);
            Hzt1[j] = newfactor[i] * (jbmca[j-1]*Ht2[indm1][j-1] + kemcd[k]*Ht2[indm1][j] - tmp2+c*tmp1);
        }

        // Construct Ht

        real8 onem2factor = 1.0 - factor[i] * 2*i;

        for (int j=1,k=2*i-2; j<(2*i); j++,k--)
        {
            complex8 cpx;

            cpx    = (y0 * Jzt10[indm1][j]);
            cpx   += (z0 * Jyt10[indm1][j]);
            cpx   -= (y0 * Jzt11[indm1][j]);
            cpx   -= (z1 * Jyt11[indm1][j]);
            cpx   -= (y1 * Jzt12[indm1][j]);
            cpx   -= (z0 * Jyt12[indm1][j]);
            cpx   += (y1 * Jzt13[indm1][j]);
            cpx   += (z1 * Jyt13[indm1][j]);

            Ht0[ind][j] = a2d2inv * (factor[i] * (cpx + indgd[j-1]*Ht1[indm1][j-1] + indtd[k]*Ht1[indm1][j]) + onem2factor*Ht2[indm1][j]);
        }

        Ht1[ind][1] = delta   * Hyt0[1] +
                      epsilon * Hzt0[1] +
                      thetad  * Ht0[ind][1];

        for (int j=1; j<(2*i); j++)
        {
            Ht1[ind][j+1] = alpha  * Hyt0[j] + beta   * Hzt0[j] + gammad * Ht0[ind][j];
        }

        Ht2[ind][1] = delta   * Hyt1[1] +
                      epsilon * Hzt1[1] +
                      thetad  * Ht1[ind][1];

        for (int j=1; j<(2*i+1); j++)
        {
            Ht2[ind][j+1] = alpha  * Hyt1[j] + beta   * Hzt1[j] + gammad * Ht1[ind][j];
        }

        for (int j=1; j<(2*i+1); j++)
        {
            reHzt = creal(Hzt1[j]);
            reHyt = creal(Hyt1[j]);

            imHzt = cimag(Hzt1[j]);
            imHyt = cimag(Hyt1[j]);

            reHt  = creal(Ht1[ind][j]);
            imHt  = cimag(Ht1[ind][j]);

            reHy1 = reHyt- y0*reHt;
            reHy2 = reHyt- y1*reHt;
            reHz1 = reHzt- z0*reHt;
            reHz2 = reHzt- z1*reHt;

            imHy1 = imHyt- y0*imHt;
            imHy2 = imHyt- y1*imHt;
            imHz1 = imHzt- z0*imHt;
            imHz2 = imHzt- z1*imHt;

            for (int k=0; (k<18); ++k)
            {
                real8 a = FqReal[k][local];
                real8 b = FqImag[k][local];

                G1[k] += b*imHz2 - a*reHz2;
                G2[k] += a*reHz1 - b*imHz1;
                G3[k] += b*imHy2 - a*reHy2;
                G4[k] += a*reHy1 - b*imHy1;
            }

            local++;
        }
    }  /* end for [i = 2; i<q; ...] */

    cuda_syncthreads();

    // End of recurrence : H integrals computed.
    // Multiply H integral with precomputed vector F

#if 0
    complex8 HzmH1[(QMAX+1)*(QMAX+2)];
    complex8 HzmH2[(QMAX+1)*(QMAX+2)];
    complex8 HymH1[(QMAX+1)*(QMAX+2)];
    complex8 HymH2[(QMAX+1)*(QMAX+2)];

    VZERO(G1,18);
    VZERO(G2,18);
    VZERO(G3,18);
    VZERO(G4,18);

    local = 0;
    for (int i=0; i<((qmax+1)*(qmax+2)); ++i)
    {
        reHz1 = creal(HzmH1[local]);
        reHz2 = creal(HzmH2[local]);
        reHy1 = creal(HymH1[local]);
        reHy2 = creal(HymH2[local]);

        imHz1 = cimag(HzmH1[local]);
        imHz2 = cimag(HzmH2[local]);
        imHy1 = cimag(HymH1[local]);
        imHy2 = cimag(HymH2[local]);

        for (int j=0; (j<18); j++)
        {
            real8 ReFqval = FqReal[j][local];
            real8 ImFqval = FqImag[j][local];

            G1[j] -= reHz2*ReFqval;
            G2[j] += reHz1*ReFqval;
            G3[j] -= reHy2*ReFqval;
            G4[j] += reHy1*ReFqval;

            G1[j] += imHz2*ImFqval;
            G2[j] -= imHz1*ImFqval;
            G3[j] += imHy2*ImFqval;
            G4[j] -= imHy1*ImFqval;
        }
        local++;
    }
#endif

    // Compute pre-factor eps C C b t bp tp in the form of
    // matrices Ap, B, Bp, A

    real8 bptp11 = b12[0]*tp[0];
    real8 bptp22 = b12[1]*tp[1];
    real8 bptp33 = b12[2]*tp[2];

    real8 Ap51   = b12[0]*tp[1];
    real8 Ap62   = b12[1]*tp[2];
    real8 Ap43   = b12[2]*tp[0];

    real8 Ap23   = b12[1]*tp[0];
    real8 Ap31   = b12[2]*tp[1];
    real8 Ap12   = b12[0]*tp[2];

    real8 Ap63   = bptp11 - bptp22;
    real8 Ap41   = bptp22 - bptp33;
    real8 Ap52   = bptp33 - bptp11;

    real8 B31    = b34[0]*t[1];
    real8 B12    = b34[1]*t[2];
    real8 B23    = b34[2]*t[0];

    real8 B24    = b34[1]*t[0];
    real8 B35    = b34[2]*t[1];
    real8 B16    = b34[0]*t[2];

    real8 bt11   = b34[0]*t[0];
    real8 bt22   = b34[1]*t[1];
    real8 bt33   = b34[2]*t[2];

    real8 B36    = bt22 - bt11;
    real8 B14    = bt33 - bt22;
    real8 B25    = bt11 - bt33;

    real8 AG1 [6]   ; VZERO(AG1 ,6);
    real8 AG2 [6]   ; VZERO(AG2 ,6);
    real8 ApG3[6]   ; VZERO(ApG3,6);
    real8 ApG4[6]   ; VZERO(ApG4,6);
    real8 A   [6][3]; MZERO(A   ,6,3);
    real8 Ap  [6][3]; MZERO(Ap  ,6,3);
    real8 b0  [6]   ; VZERO(b0  ,6);
    real8 bp0 [6]   ; VZERO(bp0 ,6);
    real8 b1  [6]   ; VZERO(b1  ,6);
    real8 bp1 [6]   ; VZERO(bp1 ,6);
    real8 b2  [6]   ; VZERO(b2  ,6);
    real8 bp2 [6]   ; VZERO(bp2 ,6);

    // For general elastic constant matrix

    for (int i=0; (i<6); ++i)
    {
        Ap[i][0] = -ecm[i][1]*Ap62 + ecm[i][2]*Ap31 + ecm[i][3]*Ap41 + ecm[i][4]*Ap51 - ecm[i][5]*Ap12;
        Ap[i][1] =  ecm[i][0]*Ap12 - ecm[i][2]*Ap43 - ecm[i][3]*Ap23 + ecm[i][4]*Ap52 + ecm[i][5]*Ap62;
        Ap[i][2] = -ecm[i][0]*Ap51 + ecm[i][1]*Ap23 + ecm[i][3]*Ap43 - ecm[i][4]*Ap31 + ecm[i][5]*Ap63;
    }

    for (int i=0; (i<6); ++i)
    {
        bp0[i]   =  ecm[1][i]*Ap62 - ecm[2][i]*Ap31 - ecm[3][i]*Ap41 - ecm[4][i]*Ap51 + ecm[5][i]*Ap12;
        bp1[i]   = -ecm[0][i]*Ap12 + ecm[2][i]*Ap43 + ecm[3][i]*Ap23 - ecm[4][i]*Ap52 - ecm[5][i]*Ap62;
        bp2[i]   =  ecm[0][i]*Ap51 - ecm[1][i]*Ap23 - ecm[3][i]*Ap43 + ecm[4][i]*Ap31 - ecm[5][i]*Ap63;
    }

    for (int i=0; (i<6); ++i)
    {
        b0[i]    =  ecm[1][i]*B12  - ecm[2][i]*B35  + ecm[3][i]*B14  - ecm[4][i]*B31  + ecm[5][i]*B16;
        b1[i]    = -ecm[0][i]*B16  + ecm[2][i]*B23  + ecm[3][i]*B24  + ecm[4][i]*B25  - ecm[5][i]*B12;
        b2[i]    =  ecm[0][i]*B31  - ecm[1][i]*B24  - ecm[3][i]*B23  + ecm[4][i]*B35  + ecm[5][i]*B36;
    }

    for (int i=0; (i<6); ++i)
    {
        A[i][0]  = -ecm[i][1]*B12  + ecm[i][2]*B35  - ecm[i][3]*B14  + ecm[i][4]*B31  - ecm[i][5]*B16;
        A[i][1]  =  ecm[i][0]*B16  - ecm[i][2]*B23  - ecm[i][3]*B24  - ecm[i][4]*B25  + ecm[i][5]*B12;
        A[i][2]  = -ecm[i][0]*B31  + ecm[i][1]*B24  + ecm[i][3]*B23  - ecm[i][4]*B35  - ecm[i][5]*B36;
    }

    cuda_syncthreads();

    real8 AGaux1 [3][3]; M33_ZERO(AGaux1 );
    real8 AGaux2 [3][3]; M33_ZERO(AGaux2 );
    real8 ApGaux3[3][3]; M33_ZERO(ApGaux3);
    real8 ApGaux4[3][3]; M33_ZERO(ApGaux4);

    for (int i=0; (i<3); ++i)
    {
        for (int j=0; (j<3); ++j)
        {
            for (int k=0,m=6*j; (k<6); ++k,++m)
            {
                AGaux1 [i][j] += A [k][i]*G1[m];
                AGaux2 [i][j] += A [k][i]*G2[m];
            }

            for (int k=0,m=6*j; (k<6); ++k,++m)
            {
                ApGaux3[i][j] += Ap[k][i]*G3[m];
                ApGaux4[i][j] += Ap[k][i]*G4[m];
            }
        }
    }

    AG1 [0] = AGaux1[0][0];
    AG1 [1] = AGaux1[1][1];
    AG1 [2] = AGaux1[2][2];

    AG1 [3] = AGaux1[1][2] + AGaux1[2][1];
    AG1 [4] = AGaux1[2][0] + AGaux1[0][2];
    AG1 [5] = AGaux1[0][1] + AGaux1[1][0];

    AG2 [0] = AGaux2[0][0];
    AG2 [1] = AGaux2[1][1];
    AG2 [2] = AGaux2[2][2];

    AG2 [3] = AGaux2[1][2] + AGaux2[2][1];
    AG2 [4] = AGaux2[2][0] + AGaux2[0][2];
    AG2 [5] = AGaux2[0][1] + AGaux2[1][0];

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

    // Deduce forces

    if (f1)
    {
        V3_ZERO(f1);

        for (int i=0; (i<6); ++i)
        {
            f1[0] += bp0[i] * AG1 [i];
            f1[1] += bp1[i] * AG1 [i];
            f1[2] += bp2[i] * AG1 [i];
        }

        f1[0] *= coeff12;
        f1[1] *= coeff12;
        f1[2] *= coeff12;
    }

    if (f2)
    {
        V3_ZERO(f2);

        for (int i=0; (i<6); ++i)
        {
            f2[0] += bp0[i] * AG2 [i];
            f2[1] += bp1[i] * AG2 [i];
            f2[2] += bp2[i] * AG2 [i];
        }

        f2[0] *= coeff12;
        f2[1] *= coeff12;
        f2[2] *= coeff12;
    }

    if (f3)
    {
        V3_ZERO(f3);

        for (int i=0; (i<6); ++i)
        {
            f3[0] += b0 [i] * ApG3[i];
            f3[1] += b1 [i] * ApG3[i];
            f3[2] += b2 [i] * ApG3[i];
        }

        f3[0] *= coeff34;
        f3[1] *= coeff34;
        f3[2] *= coeff34;
    }

    if (f4)
    {
        V3_ZERO(f4);

        for (int i=0; (i<6); ++i)
        {
            f4[0] += b0 [i] * ApG4[i];
            f4[1] += b1 [i] * ApG4[i];
            f4[2] += b2 [i] * ApG4[i];
        }

        f4[0] *= coeff34;
        f4[1] *= coeff34;
        f4[2] *= coeff34;
    }
}

