#pragma once

#ifndef _PDS_SSF_ANISO_H
#define _PDS_SSF_ANISO_H

#include <stdlib.h>

#include "cuda_portability.h"

#include "Typedefs.h"
#include "SSF_FV_t.h"
#include "SSF_PV_t.h"

//------------------------------------------------------------------------------------------------------------

__cuda_host__
size_t SSF_Aniso_Stack_Bytes (const int qmax);

//------------------------------------------------------------------------------------------------------------

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
);

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
);

//------------------------------------------------------------------------------------------------------------

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
);

//------------------------------------------------------------------------------------------------------------

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
);

//------------------------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Aniso
(
         real8  *f1        ,    ///< resulting force on node 1 <xyz>
         real8  *f2        ,    ///< resulting force on node 2 <xyz>
         real8  *f3        ,    ///< resulting force on node 3 <xyz>
         real8  *f4        ,    ///< resulting force on node 4 <xyz>
   const real8  *p1        ,    ///< position of node 1 <xyz>
   const real8  *p2        ,    ///< position of node 2 <xyz>
   const real8  *p3        ,    ///< position of node 3 <xyz>
   const real8  *p4        ,    ///< position of node 4 <xyz>
   const real8  *b1        ,    ///< burgers vector for segment p1->p2
   const real8  *b3        ,    ///< burgers vector for segment p3->p4
         Home_t *home           ///< points to home data structure
);

//------------------------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Aniso
(
         real8         *f1   ,  ///< resulting force on node 1 (xyz)
         real8         *f2   ,  ///< resulting force on node 2 (xyz)
         real8         *f3   ,  ///< resulting force on node 3 (xyz)
         real8         *f4   ,  ///< resulting force on node 4 (xyz)
   const real8         *p1   ,  ///< node 1 position (xyz)
   const real8         *p2   ,  ///< node 2 position (xyz)
   const real8         *p3   ,  ///< node 3 position (xyz)
   const real8         *p4   ,  ///< node 4 position (xyz)
   const real8         *b1   ,  ///< burgers vector (p1->p2) (xyz)
   const real8         *b3   ,  ///< burgers vector (p3->p4) (xyz)
   const real8          a    ,  ///< core radius
   const int            qmax ,  ///< base factor for spherical harmonics expansion
   const int            mtype,  ///< material type index
   const real8          ecrit,  ///< critical angle for parallelism test
   const real8          lx   ,  ///< simulation box size (x)
   const real8          ly   ,  ///< simulation box size (y)
   const real8          lz   ,  ///< simulation box size (z)
   const real8          sx   ,  ///< reciprocal of simulation box size (1.0/lx) (zero if PBC not active)
   const real8          sy   ,  ///< reciprocal of simulation box size (1.0/ly) (zero if PBC not active)
   const real8          sz   ,  ///< reciprocal of simulation box size (1.0/lz) (zero if PBC not active)
   const real8         *c66  ,  ///< elastic constants matrix (6x6)
   const complex8      *c12  ,  ///< complex elements of aniso rotation matrix (1x3,complex)
   const real8         *e3   ,  ///< real    elements of aniso rotation matrix (1x3,real   )
   const real8         *fqr  ,  ///< Fq table (real)
   const real8         *fqi     ///< Fq table (imaginary)
);

//------------------------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Aniso
(
         real8         *fv   ,  ///< resulting forces                     (f1,f2,f3,f4,...)
   const real8         *pv   ,  ///< source positions and burgers vectors (p1,p2,p3,p4,b12,b34,...)
   const real8          a    ,  ///< core radius
   const int            qmax ,  ///< base factor for spherical harmonics expansion
   const int            mtype,  ///< material type index
   const real8          ecrit,  ///< critical angle for parallelism test
   const real8          lx   ,  ///< simulation box size (x)
   const real8          ly   ,  ///< simulation box size (y)
   const real8          lz   ,  ///< simulation box size (z)
   const real8          sx   ,  ///< reciprocal of simulation box size (1.0/lx) (zero if PBC not active)
   const real8          sy   ,  ///< reciprocal of simulation box size (1.0/ly) (zero if PBC not active)
   const real8          sz   ,  ///< reciprocal of simulation box size (1.0/lz) (zero if PBC not active)
   const real8         *c66  ,  ///< elastic constants matrix (6x6)
   const complex8      *c12  ,  ///< complex elements of aniso rotation matrix (1x3,complex)
   const real8         *e3   ,  ///< real    elements of aniso rotation matrix (1x3,real   )
   const real8         *fqr  ,  ///< Fq table (real)
   const real8         *fqi     ///< Fq table (imaginary)
);

//------------------------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Aniso
(
         SSF_FV_t      *fv   ,  ///< returned forces        (f1,f2,f3,f4,...)
   const real8         *nv   ,  ///< source node positions  (n1,n2,n3,n4,n5,n6,n7,...)
   const SSF_PV_t      *pv   ,  ///< source pvec structures (pv1,pv2,pv3,pv4,pv5,...)
   const real8          a    ,  ///< core radius
   const int            qmax ,  ///< base factor for spherical harmonics expansion
   const int            mtype,  ///< material type index
   const real8          ecrit,  ///< critical angle for parallelism test
   const real8          lx   ,  ///< simulation box size (x)
   const real8          ly   ,  ///< simulation box size (y)
   const real8          lz   ,  ///< simulation box size (z)
   const real8          sx   ,  ///< reciprocal of simulation box size (1.0/lx) (zero if PBC not active)
   const real8          sy   ,  ///< reciprocal of simulation box size (1.0/ly) (zero if PBC not active)
   const real8          sz   ,  ///< reciprocal of simulation box size (1.0/lz) (zero if PBC not active)
   const real8         *c66  ,  ///< elastic constants matrix (6x6)
   const complex8      *c12  ,  ///< complex elements of aniso rotation matrix (1x3,complex)
   const real8         *e3   ,  ///< real    elements of aniso rotation matrix (1x3,real   )
   const real8         *fqr  ,  ///< Fq table (real)
   const real8         *fqi     ///< Fq table (imaginary)
);

#endif   //  _PDS_SSF_ANISO_H
