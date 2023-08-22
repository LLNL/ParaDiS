#pragma once

#ifndef _PDS_SSF_DRIVER_H
#define _PDS_SSF_DRIVER_H

#include <stdlib.h>

#include "cuda_portability.h"

#include "Typedefs.h"
#include "Home.h"
#include "Node.h"
#include "Segment.h"
#include "AnisotropicVars_t.h"
#include "SSF_FV_t.h"
#include "SSF_PV_t.h"
#include "SSF_Node_t.h"

//------------------------------------------------------------------------------------------------------------
// expose module statics...
//------------------------------------------------------------------------------------------------------------

__cuda_host__ int            SSF_Initialized     (void);

__cuda_host__ size_t         SSF_Mbuf_CPU_Bytes  (void);
__cuda_host__ size_t         SSF_Mbuf_GPU_Bytes  (void);
__cuda_host__ unsigned char *SSF_Mbuf_CPU        (void);
__cuda_host__ unsigned char *SSF_Mbuf_GPU        (void);

__cuda_host__ real8          SSF_A               (void);
__cuda_host__ real8          SSF_MU              (void);
__cuda_host__ real8          SSF_NU              (void);
__cuda_host__ int            SSF_Qmax            (void);
__cuda_host__ int            SSF_Mtype           (void);
__cuda_host__ real8          SSF_Ecrit           (void);

__cuda_host__ int            SSF_PBX             (void);
__cuda_host__ int            SSF_PBY             (void);
__cuda_host__ int            SSF_PBZ             (void);

__cuda_host__ real8          SSF_LX              (void);
__cuda_host__ real8          SSF_LY              (void);
__cuda_host__ real8          SSF_LZ              (void);

__cuda_host__ real8          SSF_SX              (void);
__cuda_host__ real8          SSF_SY              (void);
__cuda_host__ real8          SSF_SZ              (void);

__cuda_host__ real8         *SSF_Aniso_C66C      (void);
__cuda_host__ complex8      *SSF_Aniso_C12C      (void);
__cuda_host__ real8         *SSF_Aniso_E3C       (void);
__cuda_host__ real8         *SSF_Aniso_FQRC      (void);
__cuda_host__ real8         *SSF_Aniso_FQIC      (void);

__cuda_host__ real8         *SSF_Aniso_C66G      (void);
__cuda_host__ complex8      *SSF_Aniso_C12G      (void);
__cuda_host__ real8         *SSF_Aniso_E3G       (void);
__cuda_host__ real8         *SSF_Aniso_FQRG      (void);
__cuda_host__ real8         *SSF_Aniso_FQIG      (void);

__cuda_host__ size_t         SSF_Nodes_CPU_Bytes (void);
__cuda_host__ size_t         SSF_Nodes_GPU_Bytes (void);
__cuda_host__ unsigned char *SSF_Nodes_CPU       (void);
__cuda_host__ unsigned char *SSF_Nodes_GPU       (void);

//------------------------------------------------------------------------------------------------------------
// memory allocation...
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Allocate
(
   const size_t         bytes   ///< host and device buffer sizes
);

__cuda_host__
void SSF_Allocate
(
   const int            nn,   ///< number of nodes
   const int            np    ///< number of segment pairs
);

__cuda_host__
void SSF_Allocate_Nodes
(
         Node_t       **narr     ,    ///< array  of node pointers (may be sparse)
   const int            narr_cnt      ///< number of node pointers
);

__cuda_host__
void SSF_Free (void);

//------------------------------------------------------------------------------------------------------------
// initialization...
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Initialize_CUDA
(
   const int            nthrd ,  ///< number of cuda threads per stream block
   const int            npblk ,  ///< number of pairs per stream block
   const int            nstrm    ///< number of cuda streams
);

__cuda_host__
void SSF_Initialize_Iso
(
   const real8          a     ,  ///< core radius (b)
   const real8          mu    ,  ///< shear modulus
   const real8          nu    ,  ///< poisson ratio
   const real8          ecrit ,  ///< critical angle for parallelism test
   const int            pbx   ,  ///< pbc active (x) (1=yes,0=no)
   const int            pby   ,  ///< pbc active (y) (1=yes,0=no)
   const int            pbz   ,  ///< pbc active (z) (1=yes,0=no)
   const real8          lx    ,  ///< simulation box size (x)
   const real8          ly    ,  ///< simulation box size (y)
   const real8          lz    ,  ///< simulation box size (z)
   const int            nthrd ,  ///< number of cuda threads per stream block
   const int            npblk ,  ///< number of pairs per stream block
   const int            nstrm    ///< number of cuda streams
);

__cuda_host__
void SSF_Initialize_Aniso
(
   const real8          a    ,   ///< core radius (b)
   const real8          mu   ,   ///< shear modulus
   const real8          nu   ,   ///< poisson ratio
   const real8          ecrit,   ///< critical angle for parallelism test
   const int            pbx  ,   ///< pbc active (x) (1=yes,0=no)
   const int            pby  ,   ///< pbc active (y) (1=yes,0=no)
   const int            pbz  ,   ///< pbc active (z) (1=yes,0=no)
   const real8          lx   ,   ///< simulation box size (x)
   const real8          ly   ,   ///< simulation box size (y)
   const real8          lz   ,   ///< simulation box size (z)
   const int            qmax ,   ///< spherical harmonics expansion factor
   const int            mtype,   ///< material type index
   const real8         *c66  ,   ///< elastic constants matrix (6x6)
   const complex8      *c12  ,   ///< complex elements of aniso rotation matrix (1x3,complex)
   const real8         *e3   ,   ///< real    elements of aniso rotation matrix (1x3,real   )
   const real8         *fqr  ,   ///< Fq table (real)      size=18*(qmax+2)*(qmax+1)
   const real8         *fqi  ,   ///< Fq table (imaginary) size=18*(qmax+2)*(qmax+1)
   const int            nthrd,   ///< number of cuda threads per stream block
   const int            npblk,   ///< number of pairs per stream block
   const int            nstrm    ///< number of cuda streams
);

__cuda_host__
void SSF_Initialize
(
   Home_t              *home  ///< points to Home data structure
);

//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Update_Constants
(
   Home_t *home  ///< points to Home data structure
);

__cuda_host__
void SSF_Update_Constants
(
   const real8          a      ,  ///< core radius (b)
   const real8          mu     ,  ///< shear modulus
   const real8          nu     ,  ///< poisson ratio
   const real8          ecrit  ,  ///< critical angle for parallelism test
   const int            pbx    ,  ///< pbc active (x) (1=yes,0=no)
   const int            pby    ,  ///< pbc active (y) (1=yes,0=no)
   const int            pbz    ,  ///< pbc active (z) (1=yes,0=no)
   const real8          lx     ,  ///< simulation box size (x)
   const real8          ly     ,  ///< simulation box size (y)
   const real8          lz     ,  ///< simulation box size (z)
   const int            qmax =0,  ///< spherical harmonics expansion factor
   const int            mtype=0   ///< material type index
);

//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Nodes_Serialize
(
   const SSF_Node_t    *nodes,    ///< array of nodes
   const int            ncnt      ///< number of nodes
);


__cuda_host__
void SSF_Nodes_Serialize
(
   const Node_t       **narr     ,    ///< array of node pointers (may be sparse)
   const int            narr_cnt      ///< length of the node pointer array
);

//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Aniso_CPU
(
         real8         *fv,  ///< returned forces                      (f1,f2,f3,f4,....)
         real8         *pv,  ///< source positions and burgers vectors (p1,p2,p3,p4,b12,b34,...)
   const int            np   ///< number of segment pairs
);

__cuda_host__
void SSF_Aniso_CPU
(
         SSF_FV_t      *fv,  ///< returned forces        (f1,f2,f3,f4,...)
   const real8         *nv,  ///< source node positions  (n1,n2,n3,n4,n5,n6,n7,...)
   const SSF_PV_t      *pv,  ///< source pvec structures (pv1,pv2,pv3,pv4,pv5,pv6,...)
   const int            np   ///< number of segment pairs
);

__cuda_host__
void SSF_Iso_CPU
(
         real8         *fv,  ///< returned forces                      (f1,f2,f3,f4,....)
         real8         *pv,  ///< source positions and burgers vectors (p1,p2,p3,p4,b12,b34,...)
   const int            np   ///< number of segment pairs
);

__cuda_host__
void SSF_Iso_CPU
(
          SSF_FV_t     *fv,  ///< returned forces        (f1,f2,f3,f4,...)
   const real8         *nv,  ///< source node positions  (n1,n2,n3,n4,n5,n6,n7,...)
   const SSF_PV_t      *pv,  ///< source pvec structures (pv1,pv2,pv3,pv4,pv5,pv6,...)
   const int            np   ///< number of segment pairs
);

__cuda_host__
void SSF_Aniso_Force_Reset_CPU
(
          SSF_Node_t   *nv   ,  ///< array  of node structures
   const int            ncnt    ///< number of nodes
);

__cuda_host__
void SSF_Aniso_Force_Sum_CPU
(
          SSF_Node_t   *nv   ,  ///< node structures
   const SSF_FV_t      *fv   ,  ///< source computed forces     (f1,f2,f3,f4,...)
   const SSF_PV_t      *pv   ,  ///< source pair vector structs
   const int            ncnt ,  ///< number of nodes in the node vector
   const int            pcnt    ///< number of segment pairs
);

__cuda_host__
void SSF_GPU
(
         real8         *fv   ,  ///< resulting forces            (f1,f2,f3,f4,....)
   const real8         *nv   ,  ///< source node positions       (n1,n2,n3,n4,n5,n6,...)
   const real8         *pv   ,  ///< source pair vector structs  (p1,p2,p3,p4,b1,b3,...)
   const int            nn   ,  ///< number of source nodes
   const int            np   ,  ///< number of segment pairs
   const int            mode    ///< isotropy mode (0=isotropic,1=anisotropic)
);

__cuda_host__
void SSF_GPU
(
          SSF_FV_t     *fv  ,  ///< resulting forces            (f1,f2,f3,f4,....)
   const real8         *nv  ,  ///< source node positions       (n1,n2,n3,n4,n5,n6,...)
   const SSF_PV_t      *pv  ,  ///< source pair vector structs  (p1,p2,p3,p4,p5,p6,...)
   const int            nn  ,  ///< number of source nodes
   const int            np  ,  ///< number of segment pairs
   const int            mode   ///< isotropy mode (0=isotropic,1=anisotropic)
);

__cuda_host__
void SSF_GPU_Strm
(
         SSF_FV_t      *fv  ,   ///< resulting forces        (f1,f2,f3,f4,....)
   const real8         *nv  ,   ///< source node positions   (n1,n2,n3,n4,n5,n6,...)
   const SSF_PV_t      *pv  ,   ///< source positions        (p1,p2,p3,p4,b12,b34,...)
   const int            nn  ,   ///< number of source nodes
   const int            np  ,   ///< number of segment pairs
   const int            mode    ///< isotropy mode (0=isotropic,1=anisotropic)
);

__cuda_host__
void SSF_Index_Nodes
(
         Node_t       **nodes ,  ///< sparse array of node pointers (can include null pointers)
   const int            ncnt     ///< number of node pointers
);

__cuda_host__
Node_t **SSF_Gather_Nodes
(
         int          & nn     ,  ///< returns the actual node count (native+ghost)
         Node_t       **native ,  ///< sparse array of native node pointers (can include null pointers)
   const int            ncnt   ,  ///< number of native node pointers
         Node_t       **ghosts ,  ///< sparse array of ghost  node pointers (can include null pointers)
   const int            gcnt      ///< number of ghost node pointers
);

__cuda_host__
void SSF_Gather_Nodes
(
         real8         *nv   ,  ///< source node positions   (n1,n2,n3,n4,n5,n6,...)
         Node_t       **nodes,  ///< dense array of node pointers (assembled above)
   const int            nn      ///< number of nodes
);

__cuda_host__
void SSF_Send_Nodes
(
         Node_t       **narr   ,  ///< array of node pointers
   const int            ncnt      ///< node count
);

__cuda_host__
void SSF_Recv_Nodes
(
         Node_t       **narr   ,  ///< array of node pointers
   const int            ncnt      ///< node count
);

__cuda_host__
void SSF_Gather_Pairs
(
         real8         *pv   ,  ///< source positions        (p1,p2,p3,p4,b12,b34,...)
         SegmentPair_t *sp   ,  ///< array  of source segment pairs
   const int            np      ///< number of source segment pairs
);

__cuda_host__
void SSF_Gather_Pairs
(
         SSF_PV_t      *pv   ,  ///< source positions        (p1,p2,p3,p4,b12,b34,...)
         SegmentPair_t *sp   ,  ///< array  of source segment pairs
   const int            np      ///< number of source segment pairs
);

__cuda_host__
void SSF_Scatter_Forces
(
         SegmentPair_t *sp  , ///< array  of source segment pairs
         SSF_FV_t      *fv  , ///< computed segment pair forces
   const int            np    ///< number of source segment pairs
);

//------------------------------------------------------------------------------------------------------------
// main driver...
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Compute_Forces
(
         Home_t        *home,  ///< points to Home data structure
         SegmentPair_t *sp  ,  ///< array  of source segment pairs
   const int            np  ,  ///< number of source segment pairs
   const int            mode   ///< isotropy mode (0=isotropic,1=anisotropic)
);

#endif  //  _PDS_SSF_DRIVER_H
