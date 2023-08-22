/**************************************************************************
 *
 *      Module:  This module contains a simple wrapper function
 *               to invoke the seg/seg force function appropriate
 *               to the manner in which the code has been configured.
 *
 *      Includes public functions:
 *              SegSegForce()
 *
 *************************************************************************/

#include "mpi_portability.h"
#include "cuda_portability.h"

#include "Home.h"
#include "V3.h"
#include "SSF_Iso.h"
#include "SSF_Aniso.h"

// SegSegForce()
//
// Wrapper function which invokes an appropriate force function based on the manner
// in which the code has been compiled.
//-------------------------------------------------------------------------------------------------------

void SegSegForce
(
   Home_t *home,                             ///< points to home data structure
   real8 p1x, real8 p1y, real8 p1z,          ///< endpoint for first  dislocation segment (node 1)
   real8 p2x, real8 p2y, real8 p2z,          ///< endpoint for first  dislocation segment (node 2)
   real8 p3x, real8 p3y, real8 p3z,          ///< endpoint for second dislocation segment (node 3)
   real8 p4x, real8 p4y, real8 p4z,          ///< endpoint for second dislocation segment (node 4)
   real8 b1x, real8 b1y, real8 b1z,          ///< burgers vector for segment p1 to p2
   real8 b3x, real8 b3y, real8 b3z,          ///< burgers vector for segment p3 to p4
   real8 a  ,                                ///< core radius
   real8 mu ,                                ///< shear modulus
   real8 nu ,                                ///< poisson ratio
   int seg12Local,                           ///< 1 if either node of segment p1->p2 is local to the current domain, zero otherwise
   int seg34Local,                           ///< 1 if either node of segment p3->p4 is local to the current domain, zero otherwise
   real8 *f1x, real8 *f1y, real8 *f1z,       ///< pointers to receive resulting force on node 1
   real8 *f2x, real8 *f2y, real8 *f2z,       ///< pointers to receive resulting force on node 2
   real8 *f3x, real8 *f3y, real8 *f3z,       ///< pointers to receive resulting force on node 3
   real8 *f4x, real8 *f4y, real8 *f4z        ///< pointers to receive resulting force on node 4
)
{
              real8 f1[3] = { 0.0, 0.0, 0.0 };
              real8 f2[3] = { 0.0, 0.0, 0.0 };
              real8 f3[3] = { 0.0, 0.0, 0.0 };
              real8 f4[3] = { 0.0, 0.0, 0.0 };

        const real8 p1[3] = { p1x, p1y, p1z };
        const real8 p2[3] = { p2x, p2y, p2z };
        const real8 p3[3] = { p3x, p3y, p3z };
        const real8 p4[3] = { p4x, p4y, p4z };

        const real8 b1[3] = { b1x, b1y, b1z };
        const real8 b3[3] = { b3x, b3y, b3z };

        const real8 ecrit = 1.0e-4;

#ifdef ANISOTROPIC
        SSF_Aniso(f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, home);
#else
        SSF_Iso  (f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, a,mu,nu,ecrit);
#endif

        *f1x=f1[0];   *f1y=f1[1];   *f1z=f1[2];
        *f2x=f2[0];   *f2y=f2[1];   *f2z=f2[2];
        *f3x=f3[0];   *f3y=f3[1];   *f3z=f3[2];
        *f4x=f4[0];   *f4y=f4[1];   *f4z=f4[2];

}

//-------------------------------------------------------------------------------------------------------

void SegSegForce
(
         Home_t   *home ,  ///< points to home data structure
         real8    *f1   ,  ///< resulting force on node 1 <xyz>
         real8    *f2   ,  ///< resulting force on node 2 <xyz>
         real8    *f3   ,  ///< resulting force on node 3 <xyz>
         real8    *f4   ,  ///< resulting force on node 4 <xyz>
   const real8    *p1   ,  ///< position of node 1 <xyz>
   const real8    *p2   ,  ///< position of node 2 <xyz>
   const real8    *p3   ,  ///< position of node 3 <xyz>
   const real8    *p4   ,  ///< position of node 4 <xyz>
   const real8    *b1   ,  ///< burgers vector for segment x1->x2
   const real8    *b3   ,  ///< burgers vector for segment x3->x4
   const real8     a    ,  ///< core radius
   const real8     mu   ,  ///< shear modulus
   const real8     nu   ,  ///< poisson ratio
   const real8     ecrit   ///< critical angle for parallelism test
)
{
        V3_ZERO(f1);
        V3_ZERO(f2);
        V3_ZERO(f3);
        V3_ZERO(f4);

#ifdef ANISOTROPIC
        SSF_Aniso(f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, home);
#else
        SSF_Iso  (f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, a,mu,nu,ecrit);
#endif
}

