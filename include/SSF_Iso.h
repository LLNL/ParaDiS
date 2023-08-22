#pragma once

#ifndef _PDS_SSF_ISO_H
#define _PDS_SSF_ISO_H

#include "cuda_portability.h"
#include "Typedefs.h"
#include "SSF_FV_t.h"
#include "SSF_PV_t.h"

__cuda_hdev__ 
int SSF_Segs_Parallel
(
   const real8  *p1         ,  ///< position of node 1 <xyz> (segment 1)
   const real8  *p2         ,  ///< position of node 2 <xyz> (segment 1)
   const real8  *p3         ,  ///< position of node 3 <xyz> (segment 2)
   const real8  *p4         ,  ///< position of node 4 <xyz> (segment 2)
   const real8   ecrit         ///< critical angle for parallelism test 
);

__cuda_hdev__ 
void SSF_Iso_Force_Integrals
(
         real8  *Fyintegrals,  ///< points to an array to receive the Y-integrals (y0...y15)
         real8  *Fzintegrals,  ///< points to an array to receive the Z-integrals (z0...z15)
         real8  *y          ,  ///< points to an array to receive the mean yin values
         real8  *z          ,  ///< points to an array to receive the mean zin values
   const real8   a2         ,  ///<     (a*a)
   const real8   d          ,  ///< 
   const real8   c          ,  ///< 
   const real8   c2         ,  ///<     (c*c)
   const real8   onemc2     ,  ///<   (1-c*c)
   const real8   onemc2inv  ,  ///< 1/(1-c*c)
   const real8  *R          ,  ///< 
   const real8  *t          ,  ///< 
   const real8  *tp            ///< 
);

__cuda_hdev__ 
void SSF_Iso_Forces 
(
         real8  *f3         ,  ///< resulting force on node 3 <xyz>
         real8  *f4         ,  ///< resulting force on node 4 <xyz>
   const real8  *Fintegrals ,  ///< integrals array
   const real8  *t          ,  ///< 
   const real8  *tp         ,  ///< 
   const real8  *txtp       ,  ///< Cross(t,tp)
   const real8  *b          ,  ///< burgers vector <xyz>
   const real8  *bp         ,  ///< burgers vector <xyz>
   const real8   c          ,  ///<
   const real8   d          ,  ///<
   const real8   m4p        ,  ///< constants of integration
   const real8   m4pn       ,  ///<
   const real8   a2m4pn     ,  ///<
   const real8   a2m8p      ,  ///<
   const real8  *y          ,  ///<
   const real8   oneoverL      ///<
);

__cuda_hdev__ 
void SSF_Iso_Special_Force_Integrals 
(
         real8  *Fintegrals ,  ///<
   const real8   a2_d2      ,  ///<
   const real8  *y          ,  ///<
   const real8  *z             ///<
);

__cuda_hdev__ 
void SSF_Iso_Special_Forces
(
         real8  *f1         ,  ///< force on node 1 <xyz>
         real8  *f2         ,  ///< force on node 2 <xyz>
   const real8  *Fintegrals ,  ///< integrals array
   const real8  *t          ,  ///< 
   const real8  *b          ,  ///< burgers vector <xyz>
   const real8  *bp         ,  ///< burgers vector <xyz>
   const real8  *nd         ,  ///<
   const real8   m4p        ,  ///< constants of integration
   const real8   m4pn       ,  ///<
   const real8   a2m4pn     ,  ///<
   const real8   a2m8p      ,  ///<
   const real8  *y          ,  ///<
   const real8   oneoverL      ///<
);

__cuda_hdev__ 
void SSF_Iso_Special_Remote_Node_Force_Half 
(
         real8  *f3         ,  ///< force on node 3 <xyz>  
         real8  *f4         ,  ///< force on node 4 <xyz>  
   const real8  *p1         ,  ///< position of node 1 <xyz>
   const real8  *p2         ,  ///< position of node 2 <xyz>
   const real8  *p3         ,  ///< position of node 3 <xyz>
   const real8  *p4         ,  ///< position of node 4 <xyz>
   const real8  *b1         ,  ///< burgers vector (node 1->2)
   const real8  *b3         ,  ///< burgers vector (node 3->4)
   const real8   a          ,  ///< core radius
   const real8   mu         ,  ///< shear modulus
   const real8   nu         ,  ///< poisson ratio
   const real8   ecrit  
);

__cuda_hdev__ 
void SSF_Iso_Non_Parallel
(
         real8  *f1         ,  ///< force on node 1 <xyz>  (returned if not null)
         real8  *f2         ,  ///< force on node 2 <xyz>  (returned if not null)
         real8  *f3         ,  ///< force on node 3 <xyz>  (returned if not null)
         real8  *f4         ,  ///< force on node 4 <xyz>  (returned if not null)
   const real8  *p1         ,  ///< position of node 1 <xyz>
   const real8  *p2         ,  ///< position of node 2 <xyz>
   const real8  *p3         ,  ///< position of node 3 <xyz>
   const real8  *p4         ,  ///< position of node 4 <xyz>
   const real8  *b1         ,  ///< burgers vector (node 1->2)
   const real8  *b3         ,  ///< burgers vector (node 3->4)
   const real8   a          ,  ///< core radius
   const real8   mu         ,  ///< shear modulus
   const real8   nu         ,  ///< poisson ratio
   const real8   ecrit         ///< critical angle for parallelism test 
);

__cuda_hdev__ 
void SSF_Iso_Parallel 
(
         real8  *f1         ,  ///< resulting force on node 1 <xyz>
         real8  *f2         ,  ///< resulting force on node 2 <xyz>
         real8  *f3         ,  ///< resulting force on node 3 <xyz>
         real8  *f4         ,  ///< resulting force on node 4 <xyz>
   const real8  *p1         ,  ///< position of node 1 <xyz>
   const real8  *p2         ,  ///< position of node 2 <xyz>
   const real8  *p3         ,  ///< position of node 3 <xyz>
   const real8  *p4         ,  ///< position of node 4 <xyz>
   const real8  *b1         ,  ///< burgers vector for segment x1->x2
   const real8  *b3         ,  ///< burgers vector for segment x3->x4
   const real8   a          ,  ///< core radius
   const real8   mu         ,  ///< shear modulus
   const real8   nu         ,  ///< poisson ratio
   const real8   ecrit         ///< critical angle for parallelism test
);

__cuda_hdev__
void SSF_Iso
(
         real8    *f1   ,  ///< resulting force on node 1 (xyz)
         real8    *f2   ,  ///< resulting force on node 2 (xyz)
         real8    *f3   ,  ///< resulting force on node 3 (xyz)
         real8    *f4   ,  ///< resulting force on node 4 (xyz)
   const real8    *p1   ,  ///< node 1 position (xyz)
   const real8    *p2   ,  ///< node 2 position (xyz)
   const real8    *p3   ,  ///< node 3 position (xyz)
   const real8    *p4   ,  ///< node 4 position (xyz)
   const real8    *b1   ,  ///< burgers vector (p1->p2) (xyz)
   const real8    *b3   ,  ///< burgers vector (p3->p4) (xyz)
   const real8     a    ,  ///< core radius
   const real8     mu   ,  ///< shear modulus
   const real8     nu   ,  ///< poisson ratio
   const real8     ecrit,  ///< critical angle for parallelism test
   const real8     lx   ,  ///< simulation box size (x)
   const real8     ly   ,  ///< simulation box size (y)
   const real8     lz   ,  ///< simulation box size (z)
   const real8     sx   ,  ///< reciprocal of simulation box size (1.0/lx) (zero if PBC not active)
   const real8     sy   ,  ///< reciprocal of simulation box size (1.0/ly) (zero if PBC not active)
   const real8     sz      ///< reciprocal of simulation box size (1.0/lz) (zero if PBC not active)
);

__cuda_hdev__
void SSF_Iso
(
         real8    *fv   ,  ///< resulting forces                     (f1,f2,f3,f4,...)
   const real8    *pv   ,  ///< source positions and burgers vectors (p1,p2,p3,p4,b12,b34,...)
   const real8     a    ,  ///< core radius
   const real8     mu   ,  ///< shear modulus
   const real8     nu   ,  ///< poisson ratio
   const real8     ecrit,  ///< critical angle for parallelism test
   const real8     lx   ,  ///< simulation box size (x)
   const real8     ly   ,  ///< simulation box size (y)
   const real8     lz   ,  ///< simulation box size (z)
   const real8     sx   ,  ///< reciprocal of simulation box size (1.0/lx) (zero if PBC not active)
   const real8     sy   ,  ///< reciprocal of simulation box size (1.0/ly) (zero if PBC not active)
   const real8     sz      ///< reciprocal of simulation box size (1.0/lz) (zero if PBC not active)
);

__cuda_hdev__
void SSF_Iso
(
         SSF_FV_t *fv   ,  ///< returned forces        (f1,f2,f3,f4,...)
   const real8    *nv   ,  ///< source node positions  (n1,n2,n3,n4,n5,n6,n7,...)
   const SSF_PV_t *pv   ,  ///< source pvec structures (pv1,pv2,pv3,pv4,pv5,...)
   const real8     a    ,  ///< core radius
   const real8     mu   ,  ///< shear modulus
   const real8     nu   ,  ///< poisson ratio
   const real8     ecrit,  ///< critical angle for parallelism test
   const real8     lx   ,  ///< simulation box size (x)
   const real8     ly   ,  ///< simulation box size (y)
   const real8     lz   ,  ///< simulation box size (z)
   const real8     sx   ,  ///< reciprocal of simulation box size (1.0/lx) (zero if PBC not active)
   const real8     sy   ,  ///< reciprocal of simulation box size (1.0/ly) (zero if PBC not active)
   const real8     sz      ///< reciprocal of simulation box size (1.0/lz) (zero if PBC not active)
);

__cuda_hdev__ 
void SSF_Iso
(
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
);

__cuda_hdev__
void SSF_Iso
(
         real8  &f1x,       real8  &f1y,       real8  &f1z,    ///< resulting force on node 1 <xyz>
         real8  &f2x,       real8  &f2y,       real8  &f2z,    ///< resulting force on node 2 <xyz>
         real8  &f3x,       real8  &f3y,       real8  &f3z,    ///< resulting force on node 3 <xyz>
         real8  &f4x,       real8  &f4y,       real8  &f4z,    ///< resulting force on node 4 <xyz>

   const real8   p1x, const real8   p1y, const real8   p1z,    ///< position of node 1 <xyz>
   const real8   p2x, const real8   p2y, const real8   p2z,    ///< position of node 2 <xyz>
   const real8   p3x, const real8   p3y, const real8   p3z,    ///< position of node 3 <xyz>
   const real8   p4x, const real8   p4y, const real8   p4z,    ///< position of node 4 <xyz>

   const real8   b1x, const real8   b1y, const real8   b1z,    ///< burgers vector for segment p1->p2
   const real8   b3x, const real8   b3y, const real8   b3z,    ///< burgers vector for segment p3->p4

   const real8   a     ,                                       ///< core radius
   const real8   mu    ,                                       ///< shear modulus
   const real8   nu    ,                                       ///< poisson ratio
   const real8   ecrit                                         ///< critical angle for parallelism test
);

__cuda_hdev__
void SSF_Iso
(
         real8  *f1x,       real8  *f1y,       real8  *f1z,    ///< resulting force on node 1 <xyz>
         real8  *f2x,       real8  *f2y,       real8  *f2z,    ///< resulting force on node 2 <xyz>
         real8  *f3x,       real8  *f3y,       real8  *f3z,    ///< resulting force on node 3 <xyz>
         real8  *f4x,       real8  *f4y,       real8  *f4z,    ///< resulting force on node 4 <xyz>

   const real8   p1x, const real8   p1y, const real8   p1z,    ///< position of node 1 <xyz>
   const real8   p2x, const real8   p2y, const real8   p2z,    ///< position of node 2 <xyz>
   const real8   p3x, const real8   p3y, const real8   p3z,    ///< position of node 3 <xyz>
   const real8   p4x, const real8   p4y, const real8   p4z,    ///< position of node 4 <xyz>

   const real8   b1x, const real8   b1y, const real8   b1z,    ///< burgers vector for segment p1->p2
   const real8   b3x, const real8   b3y, const real8   b3z,    ///< burgers vector for segment p3->p4

   const real8   a     ,                                       ///< core radius
   const real8   mu    ,                                       ///< shear modulus
   const real8   nu    ,                                       ///< poisson ratio
   const real8   ecrit                                         ///< critical angle for parallelism test
);

#endif // _PDS_SSF_ISO_H
