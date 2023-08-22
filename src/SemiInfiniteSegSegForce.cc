/**************************************************************************
 *
 *      Module:  SemiInfiniteSegSegForce.c
 *
 *      Description: This module contains a simple wrapper function
 *               to invoke the semi-infinite seg/seg force function
 *               appropriate to the manner in which the code has 
 *               been configured.
 *
 *               This module is only applicable when the ARL FEM code
 *               has been coupled to ParaDiS.
 *
 *      Includes public functions:
 *              SemiInfiniteSegSegForce()
 *
 *************************************************************************/

#include "mpi_portability.h"

#include "Home.h"


#ifdef _ARLFEM
/*-------------------------------------------------------------------------
 *
 *      Function:       SemiInfiniteSegSegForce
 *      Description:    Wrapper function which invokes an appropriate
 *                      force function based on whether the code has
 *                      been compiled to use rational or non-rational
 *                      dislocation segment forces.
 *
 *      Arguments:
 *              p1*          Coordinates of the surface node which is the
 *                           endpoint of the semi-infinite dislocation segment.
 *              p2*          Coordinates of a point which is outside the
 *                           free surface and on the semi-infinite dislocation
 *                           segment.
 *              p3*,p4*      endpoints for a finite-length dislocation segment
 *                           starting at p3x,p3y,p3z and ending at p4x,p4y,p4z
 *              bxp,byp,bzp  burgers vector for segment p1 to p2
 *              bx,by,bz     burgers vector for segment p3 to p4
 *              a            core parameter
 *              MU           shear modulus
 *              NU           poisson ratio
 *              fp1*         Location at which to return forces at the endpoint
 *                           of the semi-infinite segment.  WARNING! Currently
 *                           force at this point is always zeroed out!
 *              fp3*,fp4*    Locations in which to return forces on nodes
 *                           located at p3 and p4 respectively
 *                      
 *-----------------------------------------------------------------------*/
void SemiInfiniteSegSegForce(Home_t *home,
                 real8 p1x, real8 p1y, real8 p1z,
                 real8 p2x, real8 p2y, real8 p2z,
                 real8 p3x, real8 p3y, real8 p3z,
                 real8 p4x, real8 p4y, real8 p4z,
                 real8 bpx, real8 bpy, real8 bpz,
                 real8 bx, real8 by, real8 bz,
                 real8 a, real8 MU, real8 NU,
                 real8 *fp1x, real8 *fp1y, real8 *fp1z,
                 real8 *fp3x, real8 *fp3y, real8 *fp3z,
                 real8 *fp4x, real8 *fp4y, real8 *fp4z)
{

   Fatal("Semi-infinite forces for parallel and non-parallel segments need to be updated\n");

#if 0 //def ANISOTROPIC
        real8 c, c2, onemc2;
        real8 lineDir1[3], lineDir2[3];

/*
 *      For now, force at p1 is explicitly zeroed...
 */
        *fp1x = 0.0;
        *fp1y = 0.0;
        *fp1z = 0.0;

/*
 *      Define the line directions and determine if the segments are
 *      parallel or not.
 */
        lineDir1[X] = p2x - p1x;
        lineDir1[Y] = p2y - p1y;
        lineDir1[Z] = p2z - p1z;

        NormalizeVec(lineDir1);

        lineDir2[X] = p4x - p3x;
        lineDir2[Y] = p4y - p3y;
        lineDir2[Z] = p4z - p3z;

        NormalizeVec(lineDir2);

        c = DotProduct(lineDir1, lineDir2);
        c2 = c * c;
        onemc2 = 1.0 - c2;

        if (onemc2 > 1.0e-04) {

            SemiInfiniteSegSegForceAnisotropicNonParallel(home,
                        p1x, p1y, p1z, p2x, p2y, p2z,
                        p3x, p3y, p3z, p4x, p4y, p4z,
                        bpx, bpy, bpz, bx, by, bz,
                        a, home->param->anisoHarmonicsNumTermsBase,
                        fp3x, fp3y, fp3z,
                        fp4x, fp4y, fp4z);
        } else {

            SemiInfiniteSegSegForceAnisotropicParallel(home,
                        p1x, p1y, p1z, p2x, p2y, p2z,
                        p3x, p3y, p3z, p4x, p4y, p4z,
                        bpx, bpy, bpz, bx, by, bz,
                        a, home->param->anisoHarmonicsNumTermsBase,
                        fp3x, fp3y, fp3z,
                        fp4x, fp4y, fp4z);
        }

#else  // ANISOTROPIC not defined
        real8 fp3xCorr, fp3yCorr, fp3zCorr;
        real8 fp4xCorr, fp4yCorr, fp4zCorr;

        SemiInfiniteSegSegForceNonRational(home,
                    p1x, p1y, p1z, p2x, p2y, p2z,
                    p3x, p3y, p3z, p4x, p4y, p4z,
                    bpx, bpy, bpz, bx, by, bz,
                    a, MU, NU,
                    fp1x, fp1y, fp1z,
                    fp3x, fp3y, fp3z,
                    fp4x, fp4y, fp4z);

#ifdef USE_RATIONAL_SEG_FORCES
/*
 *      When using 'rational' forces, the standard isotropic force calulation
 *      is done above, and corrections are calculated below and added
 *      to the force values for the finite segment.
 */
        SemiInfiniteSegSegForceRationalCorrection(home,
                    p1x, p1y, p1z, p2x, p2y, p2z,
                    p3x, p3y, p3z, p4x, p4y, p4z,
                    bpx, bpy, bpz, bx, by, bz,
                    a, MU, NU,
                    &fp3xCorr, &fp3yCorr, &fp3zCorr,
                    &fp4xCorr, &fp4yCorr, &fp4zCorr);

        *fp3x += fp3xCorr;
        *fp3y += fp3yCorr;
        *fp3z += fp3zCorr;

        *fp4x += fp4xCorr;
        *fp4y += fp4yCorr;
        *fp4z += fp4zCorr;

#endif  // end #ifdef USE_RATIONAL_SEG_FORCES
#endif  // end #else use isotropic functions

        return;
}
#endif  /* ifdef _ARLFEM */
