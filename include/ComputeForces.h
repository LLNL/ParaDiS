#pragma once

#ifndef _PDS_COMPUTE_FORCES_H
#define _PDS_COMPUTE_FORCES_H

#include "Home.h"

extern void Update_Forces
(
   Segment_t *seg,       ///< points to segment to update
   Node_t    *n1 ,       ///< points to node 1
   Node_t    *n2 ,       ///< points to node 2
   real8     *f1 ,       ///< accumulated force on node 1 <xyz>
   real8     *f2         ///< accumulated force on node 2 <xyz>
);

extern void ComputeForces
(
   Home_t *home,         ///< points to the home data structure (used to extract a,MU,NU parameters).
   Node_t *n1  ,         ///< points to node 1 (segment 1, node 1)
   Node_t *n2  ,         ///< points to node 2 (segment 1, node 2)
   Node_t *n3  ,         ///< points to node 3 (segment 2, node 1)
   Node_t *n4  ,         ///< points to node 4 (segment 2, node 2)
   real8  *f1  ,         ///< resulting force on node 1 <xyz> (returned)
   real8  *f2  ,         ///< resulting force on node 2 <xyz> (returned)
   real8  *f3  ,         ///< resulting force on node 3 <xyz> (returned)
   real8  *f4            ///< resulting force on node 4 <xyz> (returned)
);

extern void ComputeForces
(
         Home_t *home,   ///< points to the home data structure (used to extract a,MU,NU parameters).
         Node_t *n1  ,   ///< points to node 1 (segment 1, node 1)
         Node_t *n2  ,   ///< points to node 2 (segment 1, node 2)
         Node_t *n3  ,   ///< points to node 3 (segment 2, node 1)
         Node_t *n4  ,   ///< points to node 4 (segment 2, node 2)
         real8  *f1  ,   ///< resulting force on node 1 <xyz> (returned)
         real8  *f2  ,   ///< resulting force on node 2 <xyz> (returned)
         real8  *f3  ,   ///< resulting force on node 3 <xyz> (returned)
         real8  *f4  ,   ///< resulting force on node 4 <xyz> (returned)
   const real8   a   ,   ///< core radius (b)
   const real8   mu  ,   ///< shear modulus
   const real8   nu  ,   ///< poisson ratio
   const real8   lx  ,   ///< simulation box size (x)
   const real8   ly  ,   ///< simulation box size (y)
   const real8   lz  ,   ///< simulation box size (z)
   const real8   sx  ,   ///< reciprocal of simulation box size (1.0/lx) (or zero if PBC not active)
   const real8   sy  ,   ///< reciprocal of simulation box size (1.0/ly) (or zero if PBC not active)
   const real8   sz      ///< reciprocal of simulation box size (1.0/lz) (or zero if PBC not active)
);

extern void ComputeForces
(
   Home_t        *home,  ///< points to home data structure
   SegmentPair_t *sp  ,  ///< points to an array of segment pairs
   int            np     ///< number of segment pairs
);

#endif  // ifndef _COMPUTEFORCES_H
