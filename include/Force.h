#pragma once

#ifndef _PDS_FORCE_H
#define _PDS_FORCE_H

/***************************************************************************
 *
 *      Module:         Force.h
 *      Description:    This header is primarily for prototypes of various
 *                      functions used in calculating or updating
 *                      forces or stresses.
 *
 ***************************************************************************/

#include "Home.h"

void AddtoArmForce(Node_t *node, int arm, real8 f[3]);
void AddtoNodeForce(Node_t *node, real8 f[3]);
void ComputeForces(Home_t *home, Node_t *node1, Node_t *node2,
        Node_t *node3, Node_t *node4, real8 *f1, real8 *f2,
        real8 *f3, real8 *f4);
void ComputeSegSigbRem(Home_t *home, int reqType);
void deWitStress(real8 mu, real8 nu,
        real8 Sigma[][3],
        real8 burgX, real8 burgY, real8 burgZ,
        real8 xA, real8 yA, real8 zA,
        real8 xB, real8 yB, real8 zB,
        real8 x0, real8 y0, real8 z0);
void deWitInteraction(real8 mu, real8 nu, real8  Sigma[][3],
        real8  px, real8  py, real8  pz,
        real8  tp1, real8  tp2, real8  tp3,
        real8  burgX, real8  burgY, real8  burgZ,
        real8  *cntx, real8  *cnty, real8  *cntz);
void dSegImgStress(Home_t *home, real8 Sigma[][3],
        real8 px, real8 py, real8 pz,
        real8 dlx, real8 dly, real8 dlz,
        real8 burgX, real8 burgY, real8 burgZ,
        real8 rx, real8 ry, real8 rz, int pbc);

#ifdef ESHELBY

#ifdef ESHELBYFORCE
void EshelbyForce(Home_t *home, real8 EPos[3], real8 radius[3], real8 EStrain[6],
                  real8 newPos1[3], real8 newPos2[3], real8 burg[3],
                  real8 f1[3], real8 f2[3]);
#endif

#ifdef ESHELBYCORE
void SegPartCoreForce(Home_t *home, real8 EPos[3], real8 radius[3],
                      real8 rotation[2][3], real8 pos1[3], real8 pos2[3],
                      real8 b[3], real8 f1[3], real8 f2[3]);
#endif

#ifdef ESHELBYSTEP
void InclusionStepForce(real8 gamma, real8  C[3], real8 Radius[3], real8 Rotation[2][3],
                        real8 newPos1[3], real8 newPos2[3], real8 vb[3],
                        real8 f1[3], real8 f2[3]);
#endif

#endif   // ESHELBY

void ExtPKForce(real8 str[3][3],
        real8 bx, real8 by, real8 bz, real8 X1, real8 Y1, real8 Z1,
        real8 X2, real8 Y2, real8 Z2, real8 f1[3], real8 f2[3]);
void FindFSegComb(Home_t *home, real8 p0[3], real8 p1[3], real8 p2[3],
        real8 burg1[3], real8 burg2[3], real8 fp0seg1[3],
        real8 fp1seg1[3], real8 fp1seg2[3], real8 fp2seg2[3],
        real8 f0new[3], real8 f1new[3]);
void FindSubFSeg(Home_t *home, real8 p1[3], real8 p2[3], real8 burg[3],
        real8 oldfp1[3], real8 oldfp2[3], real8 newpos[3],
        real8 f0seg1[3], real8 f1seg1[3], real8 f0seg2[3],
        real8 f1seg2[3]);
void GetFieldPointStress(Home_t *home, real8 x, real8 y, real8 z,
        real8 totStress[3][3]);
void GetFieldPointStressRem(Home_t *home, real8 x, real8 y, real8 z,
        int cellX, int cellY, int cellZ, real8 totStress[3][3]);
real8 HCPEcoreFactor(Home_t *home, real8 burg[3]);
void LineTensionForce(Home_t *home, real8 x1, real8 y1, real8 z1,
        real8 x2, real8 y2, real8 z2, real8 bx, real8 by, real8 bz,
        real8 f1[3], real8 f2[3]);
void LocalSegForces(Home_t *home, int reqType);
void SubcyclingSegForces(Home_t *home, int reqType);
void NodeForce(Home_t *home, int reqType);
void OsmoticForce(Home_t *home, real8 x1, real8 y1, real8 z1,
        real8 x2, real8 y2, real8 z2, real8 bx, real8 by, real8 bz,
        real8 f1[3], real8 f2[3]);
void PKForce(real8 sigb[3],
        real8 X1, real8 Y1, real8 Z1, real8 X2, real8 Y2, real8 Z2,
        real8 f1[3], real8 f2[3]);
void ReevaluateForces(Home_t *home);

#ifndef ANISOTROPIC
void CoreForceIsotropic
(
    Home_t *home,
    real8 x1, real8 y1, real8 z1,
    real8 x2, real8 y2, real8 z2,
    real8 b[3], real8 Ecore, real8 nu,
    real8 f2Core[3]
);

void CoreForceandDerivativeIsotropic
(
    Home_t *home,
    real8 x1, real8 y1, real8 z1,
    real8 x2, real8 y2, real8 z2,
    real8 b[3], real8 Ecore, real8 nu,
    real8 fCore[3], real8 dfCoredx[3][3]
);

void SelfForceIsotropic
(
    Home_t *home,
    int coreOnly, real8 mu, real8 nu,
    real8 bx, real8 by, real8 bz,
    real8 x1, real8 y1, real8 z1,
    real8 x2, real8 y2, real8 z2,
    real8 a,  real8 Ecore,
    real8 f1[3], real8 f2[3]
);

#endif  // !ANISOTROPIC

void SegSegForce
(
    Home_t *home,
    real8 p1x, real8 p1y, real8 p1z,
    real8 p2x, real8 p2y, real8 p2z,
    real8 p3x, real8 p3y, real8 p3z,
    real8 p4x, real8 p4y, real8 p4z,
    real8 bpx, real8 bpy, real8 bpz,
    real8 bx , real8 by , real8 bz ,
    real8 a, real8 mu, real8 nu,
    int seg12Local, int seg34Local,
    real8 *fp1x, real8 *fp1y, real8 *fp1z,
    real8 *fp2x, real8 *fp2y, real8 *fp2z,
    real8 *fp3x, real8 *fp3y, real8 *fp3z,
    real8 *fp4x, real8 *fp4y, real8 *fp4z
);

void SegSegForce
(
         Home_t   *home ,         ///< points to home data structure
         real8    *f1   ,         ///< resulting force on node 1 <xyz>
         real8    *f2   ,         ///< resulting force on node 2 <xyz>
         real8    *f3   ,         ///< resulting force on node 3 <xyz>
         real8    *f4   ,         ///< resulting force on node 4 <xyz>
   const real8    *p1   ,         ///< position of node 1 <xyz>
   const real8    *p2   ,         ///< position of node 2 <xyz>
   const real8    *p3   ,         ///< position of node 3 <xyz>
   const real8    *p4   ,         ///< position of node 4 <xyz>
   const real8    *b1   ,         ///< burgers vector for segment x1->x2
   const real8    *b3   ,         ///< burgers vector for segment x3->x4
   const real8     a    ,         ///< core radius
   const real8     mu   ,         ///< shear modulus
   const real8     nu   ,         ///< poisson ratio
   const real8     ecrit=1.0e-4   ///< critical angle for parallelism test
);

void SelfForce
(
         Home_t *home,
         int coreOnly, real8 mu, real8 nu,
         real8 bx, real8 by, real8 bz,
         real8 x1, real8 y1, real8 z1,
         real8 x2, real8 y2, real8 z2,
         real8 a,  real8 Ecore,
         real8 f1[3], real8 f2[3]
);

#ifdef _ARLFEM
void SemiInfiniteSegSegForce(Home_t *home,
        real8 p1x, real8 p1y, real8 p1z,
        real8 p2x, real8 p2y, real8 p2z,
        real8 p3x, real8 p3y, real8 p3z,
        real8 p4x, real8 p4y, real8 p4z,
        real8 bpx, real8 bpy, real8 bpz,
        real8 bx, real8 by, real8 bz,
        real8 a, real8 mu, real8 nu,
        real8 *fp1x, real8 *fp1y, real8 *fp1z,
        real8 *fp3x, real8 *fp3y, real8 *fp3z,
        real8 *fp4x, real8 *fp4y, real8 *fp4z);

void SemiInfiniteSegSegForceNonRational(Home_t *home,
        real8 p1x, real8 p1y, real8 p1z,
        real8 p2x, real8 p2y, real8 p2z,
        real8 p3x, real8 p3y, real8 p3z,
        real8 p4x, real8 p4y, real8 p4z,
        real8 bpx, real8 bpy, real8 bpz,
        real8 bx, real8 by, real8 bz,
        real8 a, real8 mu, real8 nu,
        real8 *fp1x, real8 *fp1y, real8 *fp1z,
        real8 *fp3x, real8 *fp3y, real8 *fp3z,
        real8 *fp4x, real8 *fp4y, real8 *fp4z);

void SemiInfiniteSegSegForceRational(Home_t *home,
        real8 p1x, real8 p1y, real8 p1z,
        real8 p2x, real8 p2y, real8 p2z,
        real8 p3x, real8 p3y, real8 p3z,
        real8 p4x, real8 p4y, real8 p4z,
        real8 bpx, real8 bpy, real8 bpz,
        real8 bx, real8 by, real8 bz,
        real8 a, real8 mu, real8 nu,
        real8 *fp1x, real8 *fp1y, real8 *fp1z,
        real8 *fp3x, real8 *fp3y, real8 *fp3z,
        real8 *fp4x, real8 *fp4y, real8 *fp4z);

void VirtualSegForce(Home_t *home, real8 p1[3], real8 p2[3],
        real8 burg[3], real8 f1[3], real8 f2[3]);

#endif  /* if defined _ARLFEM */

void SetOneNodeForce(Home_t *home, Node_t *node1);
void StressDueToSeg(Home_t *home, real8 px, real8 py, real8 pz,
        real8 p1x, real8 p1y, real8 p1z, real8 p2x, real8 p2y, real8 p2z,
        real8 bx, real8 by, real8 bz, real8 a, real8 mu, real8 nu,
        real8 *stress);

void ZeroNodeForces(Home_t *home, Node_t **nodes, int ncnt, int reqType);
void ZeroNodeForces(Home_t *home,                           int reqType);

#endif  // _PDS_FORCE_H
