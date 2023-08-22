#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi_portability.h"
#include "cuda_portability.h"

#include "Typedefs.h"
#include "V3.h"
#include "SSF_Iso_Rational.h"

// SSF_Iso_Rational_Correction_Half()
//
// Applies the isotropic rational force correction for a single segment.
//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Iso_Rational_Correction_Half
(
         real8  *f1,    ///< corrected force on node 1 <xyz> (updated)
         real8  *f2,    ///< corrected force on node 2 <xyz> (updated)
   const real8  *p1,    ///< position of node 1 <xyz>
   const real8  *p2,    ///< position of node 2 <xyz>
   const real8  *p3,    ///< position of node 3 <xyz>
   const real8  *p4,    ///< position of node 4 <xyz>
   const real8  *b1,    ///< burgers vector for segment p1->p2
   const real8  *b3,    ///< burgers vector for segment p3->p4
   const real8   a ,    ///< core radius
   const real8   mu,    ///< shear modulus
   const real8   nu,    ///< poisson ratio
   const real8   ecrit  ///< critical angle
)
{
   const real8 p12[3]     = { (p2[0]-p1[0]),
                              (p2[1]-p1[1]),
                              (p2[2]-p1[2]) };

   const real8 invVecLen  = 1.0 / sqrt(V3_DOT(p12,p12));

   const real8 unitVec[3] = { p12[0] * invVecLen,
                              p12[1] * invVecLen,
                              p12[2] * invVecLen };

   const real8 R[4][3]    = { { (p1[0]-p3[0]), (p1[1]-p3[1]), (p1[2]-p3[2]) },
                              { (p2[0]-p3[0]), (p2[1]-p3[1]), (p2[2]-p3[2]) },
                              { (p1[0]-p4[0]), (p1[1]-p4[1]), (p1[2]-p4[2]) },
                              { (p2[0]-p4[0]), (p2[1]-p4[1]), (p2[2]-p4[2]) } };

   const real8 y[4]       = { V3_DOT(R[0], unitVec),
                              V3_DOT(R[1], unitVec),
                              V3_DOT(R[2], unitVec),
                              V3_DOT(R[3], unitVec) };

   const real8 nd[2][3]   = { { R[0][0] - (y[0] * unitVec[0]),
                                R[0][1] - (y[0] * unitVec[1]),
                                R[0][2] - (y[0] * unitVec[2]) },
                              { R[2][0] - (y[2] * unitVec[0]),
                                R[2][1] - (y[2] * unitVec[1]),
                                R[2][2] - (y[2] * unitVec[2]) } };

   const real8 d2[2]      = { V3_DOT(nd[0], nd[0]),
                              V3_DOT(nd[1], nd[1]) };

   const real8 tdb        = V3_DOT(unitVec, b1);
   const real8 tdbp       = V3_DOT(unitVec, b3);

   const real8 txbpxt[3]  = { b3[0] - (tdbp * unitVec[0]),
                              b3[1] - (tdbp * unitVec[1]),
                              b3[2] - (tdbp * unitVec[2]) };

   real8 bxt[3] = { 0.0 };
   V3_CROSS(bxt, b1, unitVec);

   const real8 bxtdbp  = V3_DOT(bxt,b3);
   const real8 txbpdb  = bxtdbp;

   const real8 nddb[2] = { V3_DOT(nd[0],b1),
                           V3_DOT(nd[1],b1) };

   real8 ndxt   [2][3] = {{ 0.0 }};
   real8 ndxbp  [2][3] = {{ 0.0 }};
   real8 ndxbpxt[2][3] = {{ 0.0 }};

   V3_CROSS(ndxt[0], nd[0], unitVec);
   V3_CROSS(ndxt[1], nd[1], unitVec);

   V3_CROSS(ndxbp[0], nd[0], b3);
   V3_CROSS(ndxbp[1], nd[1], b3);

   V3_CROSS(ndxbpxt[0], ndxbp[0], unitVec);
   V3_CROSS(ndxbpxt[1], ndxbp[1], unitVec);

   const real8 ndxbpdb[2] = { V3_DOT(ndxbp[0], b1),
                              V3_DOT(ndxbp[1], b1) };

   real8 iCorr003a[3], iCorr103a[3], iCorr203a[3];
   real8 iCorr003b[3], iCorr103b[3], iCorr203b[3];

   for (int i=0; (i<3); i++)
   {
       iCorr203a[i] = txbpxt     [i] * tdb;

       iCorr103a[i] = (ndxt   [0][i] * txbpdb    ) +
                      (txbpxt    [i] * nddb   [0]) +
                      (ndxbpxt[0][i] * tdb);

       iCorr003a[i] = (ndxt   [0][i] * ndxbpdb[0]) +
                      (ndxbpxt[0][i] * nddb   [0]);

       iCorr203b[i] = iCorr203a[i];

       iCorr103b[i] = (ndxt   [1][i] * txbpdb    ) +
                      (txbpxt    [i] * nddb   [1]) +
                      (ndxbpxt[1][i] * tdb);

       iCorr003b[i] = (ndxt   [1][i] * ndxbpdb[1]) +
                      (ndxbpxt[1][i] * nddb   [1]);
   }

   real8 f003[4], f103[4], f203[4], f303[4];

   for (int i=0; (i<4); i++)
   {
       real8 a2_d2    = a*a + d2[(i>1)];
       real8 a2_d2Inv = 1.0/a2_d2;
       real8 ysqr     = y[i]*y[i];

       real8 ra       = sqrt(ysqr + a2_d2);
       real8 ra_inv   = 1.0/ra;

       f003[i]  = y[i] * ra_inv * a2_d2Inv;
       f103[i]  = -ra_inv;
       f203[i]  = log(ra+y[i]) - y[i] * ra_inv;
       f303[i]  = ra + ra_inv * a2_d2;
   }

   real8 fCorr203a = (f303[1] - f303[0]) - (y[0] * (f203[1] - f203[0]));
   real8 fCorr103a = (f203[1] - f203[0]) - (y[0] * (f103[1] - f103[0]));
   real8 fCorr003a = (f103[1] - f103[0]) - (y[0] * (f003[1] - f003[0]));

   real8 fCorr203b = (f303[3] - f303[2]) - (y[2] * (f203[3] - f203[2]));
   real8 fCorr103b = (f203[3] - f203[2]) - (y[2] * (f103[3] - f103[2]));
   real8 fCorr003b = (f103[3] - f103[2]) - (y[2] * (f003[3] - f003[2]));

   real8 f1Corr[3] = { 0.0, 0.0, 0.0 };
   real8 f2Corr[3] = { 0.0, 0.0, 0.0 };

   for (int i=0; (i<3); i++)
   {
       f2Corr[i] = iCorr003a[i] * fCorr003a +
                   iCorr103a[i] * fCorr103a +
                   iCorr203a[i] * fCorr203a -
                   iCorr003b[i] * fCorr003b -
                   iCorr103b[i] * fCorr103b -
                   iCorr203b[i] * fCorr203b;
   }

   fCorr203a = (f303[1] - f303[0]) - (y[1] * (f203[1] - f203[0]));
   fCorr103a = (f203[1] - f203[0]) - (y[1] * (f103[1] - f103[0]));
   fCorr003a = (f103[1] - f103[0]) - (y[1] * (f003[1] - f003[0]));

   fCorr203b = (f303[3] - f303[2]) - (y[3] * (f203[3] - f203[2]));
   fCorr103b = (f203[3] - f203[2]) - (y[3] * (f103[3] - f103[2]));
   fCorr003b = (f103[3] - f103[2]) - (y[3] * (f003[3] - f003[2]));

   for (int i=0; (i<3); i++)
   {
       f1Corr[i] = iCorr003b[i] * fCorr003b +
                   iCorr103b[i] * fCorr103b +
                   iCorr203b[i] * fCorr203b -
                   iCorr003a[i] * fCorr003a -
                   iCorr103a[i] * fCorr103a -
                   iCorr203a[i] * fCorr203a;
   }

   const real8 factor = 0.125 * invVecLen * mu / M_PI;

   f1[0] += f1Corr[0] * factor;
   f1[1] += f1Corr[1] * factor;
   f1[2] += f1Corr[2] * factor;

   f2[0] += f2Corr[0] * factor;
   f2[1] += f2Corr[1] * factor;
   f2[2] += f2Corr[2] * factor;
}

// SSF_Iso_Rational_Correction()
//
// Applies the isotropic rational force correction for pair of segments.
//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Iso_Rational_Correction
(
         real8  *f1,    ///< corrected force on node 1 <xyz> (updated)
         real8  *f2,    ///< corrected force on node 2 <xyz> (updated)
         real8  *f3,    ///< corrected force on node 3 <xyz> (updated)
         real8  *f4,    ///< corrected force on node 4 <xyz> (updated)
   const real8  *p1,    ///< position of node 1 <xyz>
   const real8  *p2,    ///< position of node 2 <xyz>
   const real8  *p3,    ///< position of node 3 <xyz>
   const real8  *p4,    ///< position of node 4 <xyz>
   const real8  *b1,    ///< burgers vector for segment p1->p2
   const real8  *b3,    ///< burgers vector for segment p3->p4
   const real8   a ,    ///< core radius
   const real8   mu,    ///< shear modulus
   const real8   nu,    ///< poisson ratio
   const real8   ecrit  ///< critical angle
)
{
   SSF_Iso_Rational_Correction_Half(f1,f2, p1,p2,p3,p4, b1,b3, a,mu,nu,ecrit);
   SSF_Iso_Rational_Correction_Half(f3,f4, p3,p4,p1,p2, b3,b1, a,mu,nu,ecrit);
}

