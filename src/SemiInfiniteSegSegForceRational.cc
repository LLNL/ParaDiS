/**************************************************************************
 *
 *      Module:  SemiInfiniteSegSegForceRational.c
 *
 *      Authors:  Tom Arsenlis, Meijie Tang
 *
 *      Description: This module contains the functions needed for
 *               calculating forces from interactions between a
 *               finite length dislocation and a semi-infinite dislocation.
 *               This is only applicable when the FEM code is hooked into
 *               ParaDiS for doing simulations with free surfaces conditions.
 *
 *      Includes public functions:
 *              SemiInfiniteSegSegForceRationalCorrection()
 *
 *************************************************************************/

#include "mpi_portability.h"

#include "Home.h"


#ifdef _ARLFEM
/*-------------------------------------------------------------------------
 *
 *      Function:       SemiInfiniteSegSegForceRationalCorrection
 *      Description:    Used to calculate the 'rational' interaction force 
 *                      correction between a finite-length dislocation
 *                      segment and a semi-infinite-length dislocation segment.
 *
 *      Arguments:
 *              x1,y1,z1     Coordinates of the surface node which is the
 *                           endpoint of the semi-infinite dislocation segment.
 *              x2,y2,z2     Coordinates of a point which is outside the
 *                           free surface and on the semi-infinite dislocation
 *                           segment.
 *              x3,y3,z3,
 *              x4,y4,z4     endpoints for a finite-length dislocation segment
 *                           starting at p3x,p3y,p3z and ending at p4x,p4y,p4z
 *              bp{xyz}      burgers vector for semi-infinite segment p1 to p2
 *              b{xyz}       burgers vector for finite segment p3 to p4
 *              a            core parameter
 *              MU           shear modulus
 *              NU           poisson ratio
 *              fp3*,fp4*    Locations in which to return forces on nodes
 *                           located at p3 and p4 respectively
 *                      
 *-----------------------------------------------------------------------*/
void SemiInfiniteSegSegForceRationalCorrection(Home_t *home,
                 real8 x1, real8 y1, real8 z1,
                 real8 x2, real8 y2, real8 z2,
                 real8 x3, real8 y3, real8 z3,
                 real8 x4, real8 y4, real8 z4,
                 real8 bpx, real8 bpy, real8 bpz,
                 real8 bx, real8 by, real8 bz,
                 real8 a, real8 MU, real8 NU,
                 real8 *fp3x, real8 *fp3y, real8 *fp3z,
                 real8 *fp4x, real8 *fp4y, real8 *fp4z)
{
        int   i;
        real8 a2_d2, a2_d2Inv;
        real8 d2, invVecLen, factor;
        real8 fCorr003a, fCorr103a, fCorr203a;
        real8 tdb, tdbp, txbpdb, txbpxt[3];
        real8 nddb, ndxbpdb, nd[3], ndxt[3], ndxbp[3], ndxbpxt[3];
        real8 b[3], bp[3], bxt[3];
        real8 R[2][3], Ra[2], RaInv[2];
        real8 vec[3], unitVec[3];
        real8 y[2], ySquared[2];
        real8 f003[2], f103[2], f203[2], f303[2];
        real8 f3Corr[3], f4Corr[3];
        real8 iCorr003a[3], iCorr103a[3], iCorr203a[3];


        b[0] = bx;
        b[1] = by;
        b[2] = bz;

        bp[0] = bpx;
        bp[1] = bpy;
        bp[2] = bpz;

        vec[X] = x4 - x3;
        vec[Y] = y4 - y3;
        vec[Z] = z4 - z3;

        invVecLen = 1.0 / sqrt(DotProduct(vec, vec));

        unitVec[X] = vec[X] * invVecLen;
        unitVec[Y] = vec[Y] * invVecLen;
        unitVec[Z] = vec[Z] * invVecLen;

        R[0][0] = x3 - x1;
        R[0][1] = y3 - y1;
        R[0][2] = z3 - z1;

        R[1][0] = x4 - x1;
        R[1][1] = y4 - y1;
        R[1][2] = z4 - z1;

        y[0] = DotProduct(R[0], unitVec);
        y[1] = DotProduct(R[1], unitVec);

        for (i = 0; i < 3; i++) {
            nd[i] = R[0][i] - (y[0] * unitVec[i]);
        }

        cross(b, unitVec, bxt);
        cross(nd, unitVec, ndxt);
        cross(nd, bp, ndxbp);
        cross(ndxbp, unitVec, ndxbpxt);

        d2      = DotProduct(nd, nd);
        tdb     = DotProduct(unitVec, b);
        tdbp    = DotProduct(unitVec, bp);
        nddb    = DotProduct(nd, b);
        txbpdb  = DotProduct(bxt, bp);
        ndxbpdb = DotProduct(ndxbp, b);

        for (i = 0; i < 3; i++) {
            txbpxt[i] = bp[i] - (tdbp * unitVec[i]);
        }

        for (i = 0; i < 3; i++) {
            iCorr203a[i] = txbpxt[i] * tdb;
            iCorr103a[i] = (ndxt[i] * txbpdb)  + (txbpxt[i] * nddb) +
                           (ndxbpxt[i] * tdb);
            iCorr003a[i] = (ndxt[i] * ndxbpdb) + (ndxbpxt[i] * nddb);
        }

        a2_d2 = a * a + d2;
        a2_d2Inv = 1.0 / a2_d2;

        for (i = 0; i < 2; i++) {
            ySquared[i] = y[i] * y[i];
            Ra[i] = sqrt(ySquared[i] + a2_d2);
            RaInv[i] = 1.0 / Ra[i];
            f003[i] = y[i] * RaInv[i] * a2_d2Inv;
            f103[i] = -RaInv[i];
            f203[i] = log(Ra[i]+y[i]) - y[i] * RaInv[i];
            f303[i] = Ra[i] + a2_d2 * RaInv[i];
        }

        factor = 0.125 * invVecLen * MU / M_PI;

        fCorr203a = (f303[1] - f303[0]) - (y[0] * (f203[1] - f203[0]));
        fCorr103a = (f203[1] - f203[0]) - (y[0] * (f103[1] - f103[0]));
        fCorr003a = (f103[1] - f103[0]) - (y[0] * (f003[1] - f003[0]));

        for (i = 0; i < 3; i++) {
            f4Corr[i] = iCorr003a[i] * fCorr003a +
                        iCorr103a[i] * fCorr103a +
                        iCorr203a[i] * fCorr203a;
        }

        fCorr203a = (f303[1] - f303[0]) - (y[1] * (f203[1] - f203[0]));
        fCorr103a = (f203[1] - f203[0]) - (y[1] * (f103[1] - f103[0]));
        fCorr003a = (f103[1] - f103[0]) - (y[1] * (f003[1] - f003[0]));

        for (i = 0; i < 3; i++) {
            f3Corr[i] = -iCorr003a[i] * fCorr003a -
                         iCorr103a[i] * fCorr103a -
                         iCorr203a[i] * fCorr203a;
        }

        *fp3x = f3Corr[0] * factor;
        *fp3y = f3Corr[1] * factor;
        *fp3x = f3Corr[2] * factor;

        *fp4x = f4Corr[0] * factor;
        *fp4y = f4Corr[1] * factor;
        *fp4z = f4Corr[2] * factor;

        return;
}
#endif  /* ifdef _ARLFEM */
