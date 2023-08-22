#pragma once

#ifndef _PDS_ANISOTROPIC_H
#define _PDS_ANISOTROPIC_H

//-----------------------------------------------------------------------------------------------
//  Anisotropic.h  Contains prototypes for functions specific to the
//                 code for calculating aforces using anisotropic
//                 elasticity and other related structures, definitions etc.
//-----------------------------------------------------------------------------------------------

//--------------------------------------------------
// Prototypes for anisotropic elasticity functions
//--------------------------------------------------

void  AnisotropicInit(Home_t *home);

#ifdef _ARLFEM
void SemiInfiniteSegSegForceAnisotropicParallel
(
        Home_t *home,                                    ///< points to home data structure
        real8 p1x, real8 p1y, real8 p1z,                 ///< position of node 1 <xyz>
        real8 p2x, real8 p2y, real8 p2z,                 ///< position of node 2 <xyz>
        real8 p3x, real8 p3y, real8 p3z,                 ///< position of node 3 <xyz>
        real8 p4x, real8 p4y, real8 p4z,                 ///< position of node 4 <xyz>
        real8 bpx, real8 bpy, real8 bpz,                 ///< burgers vector (p1->p2) <xyz>
        real8 bx,  real8 by,  real8 bz,                  ///< burgers vector (p3->p4) <xyz>
        real8 a,                                         ///< core radius
        int   qmax,                                      ///< spherical harmonic expansion factor
        real8 *fp3x, real8 *fp3y, real8 *fp3z,           ///< resulting force on node 3
        real8 *fp4x, real8 *fp4y, real8 *fp4z            ///< resulting force on node 4
);

void SemiInfiniteSegSegForceAnisotropicNonParallel
(
        Home_t *home,                                    ///< points to home data structure
        real8 p1x, real8 p1y, real8 p1z,                 ///< position of node 1 <xyz>
        real8 p2x, real8 p2y, real8 p2z,                 ///< position of node 2 <xyz>
        real8 p3x, real8 p3y, real8 p3z,                 ///< position of node 3 <xyz>
        real8 p4x, real8 p4y, real8 p4z,                 ///< position of node 4 <xyz>
        real8 bpx, real8 bpy, real8 bpz,                 ///< burgers vector (p1->p2) <xyz>
        real8 bx,  real8 by,  real8 bz,                  ///< burgers vector (p3->p4) <xyz>
        real8 a,                                         ///< core radius
        int   qmax,                                      ///< spherical harmonic expansion factor
        real8 *fp3x, real8 *fp3y, real8 *fp3z,           ///< resulting force on node 3
        real8 *fp4x, real8 *fp4y, real8 *fp4z            ///< resulting force on node 4
);
#endif // _ARLFEM

void CoreForceandDerivativeAnisotropic
(
        Home_t *home,                                    ///< points to home data structure
        real8 x1, real8 y1, real8 z1,                    ///< position of node 1 <xyz>
        real8 x2, real8 y2, real8 z2,                    ///< position of node 1 <xyz>
        real8 b[3],                                      ///< burgers vector (p1->p2) <xyz>
        real8 Ecore,                                     ///<
        int   qmax,                                      ///< spherical harmonic expansion factor
        real8 F[3], real8 dFdx[3][3]                     ///<
);

void CoreForceAnisotropic
(
        Home_t *home,                                    ///< points to home data structure
        real8 x1, real8 y1, real8 z1,                    ///< position of node 1 <xyz>
        real8 x2, real8 y2, real8 z2,                    ///< position of node 1 <xyz>
        real8 b[3],                                      ///< burgers vector (p1->p2) <xyz>
        real8 Ecore,                                     ///<
        int   qmax,                                      ///< spherical harmonic expansion factor
        real8 fCore[3]                                   ///<
);

void  SelfForceAnisotropic
(
        Home_t *home, int coreOnly,                      ///< points to home data structure
        real8 bx, real8 by, real8 bz,                    ///< burgers vector (p1->p2) <xyz>
        real8 x1, real8 y1, real8 z1,                    ///< position of node 1 <xyz>
        real8 x2, real8 y2, real8 z2,                    ///< position of node 1 <xyz>
        real8 a,                                         ///< core radius
        real8 Ecore,                                     ///< 
        real8 f1[3],                                     ///< resulting force on node 1
        real8 f2[3]                                      ///< resulting force on node 2
);

void  StressDueToSegAnisotropic
(
        Home_t *home,                                    ///< points to home data structure
        real8 pointx, real8 pointy, real8 pointz,                                ///<
        real8 p1x, real8 p1y, real8 p1z,                 ///< position of node 1 <xyz>
        real8 p2x, real8 p2y, real8 p2z,                 ///< position of node 2 <xyz>
        real8 bx,  real8 by,  real8 bz,                  ///< burgers vector (p1->p2) <xyz>
        real8 a,                                         ///< core radius
        int   qmax,                                      ///< spherical harmonic expansion factor
        real8 sigma[6]                                   ///<
);

void SetElasticConstantMatrix(Home_t *home);

void SetElasticConstantMatrix4D
(
        real8 ecMatrix  [6][6],
        real8 ecMatrix4D[3][3][3][3]
);

void ComputeLegendre(int P, real8 x, real8 *legendre);

#endif // _PDS_ANISOTROPIC_H
