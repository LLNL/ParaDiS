#pragma once

#ifndef _PDS_PBC_H
#define _PDS_PBC_H

#include <math.h>

#include "cuda_portability.h"
#include "Typedefs.h"
#include "Param.h"

//-----------------------------------------------------------------------------------------------

class PBC_t
{
   public :
      BoundType_t  bndx;   ///< identifies the X-axis periodic boundary type (periodic, free, or reflecting)
      BoundType_t  bndy;   ///< identifies the Y-axis periodic boundary type (periodic, free, or reflecting)
      BoundType_t  bndz;   ///< identifies the Z-axis periodic boundary type (periodic, free, or reflecting)

      real8        x0  ;   ///< simulation box minimum extant (x)
      real8        y0  ;   ///< simulation box minimum extant (y)
      real8        z0  ;   ///< simulation box minimum extant (z)

      real8        x1  ;   ///< simulation box maximum extant (x)
      real8        y1  ;   ///< simulation box maximum extant (y)
      real8        z1  ;   ///< simulation box maximum extant (z)

      real8        cx  ;   ///< simulation center (x)
      real8        cy  ;   ///< simulation center (y)
      real8        cz  ;   ///< simulation center (z)

      real8        lx  ;   ///< simulation box size (x)
      real8        ly  ;   ///< simulation box size (y)
      real8        lz  ;   ///< simulation box size (z)

      real8        sx  ;   ///< reciprocal of simulation box size (1.0/lx) (zero if PBC not periodic)
      real8        sy  ;   ///< reciprocal of simulation box size (1.0/ly) (zero if PBC not periodic)
      real8        sz  ;   ///< reciprocal of simulation box size (1.0/lz) (zero if PBC not periodic)

   public :
     __cuda_hdev__  PBC_t(void);
     __cuda_hdev__  PBC_t(const PBC_t   &pbc  );
     __cuda_hdev__  PBC_t(const Param_t *param);
     __cuda_hdev__  PBC_t(const BoundType_t _bndx, const BoundType_t _bndy, const BoundType_t _bndz,
                          const real8       _x0  , const real8       _y0,   const real8       _z0  ,
                          const real8       _x1  , const real8       _y1,   const real8       _z1   );
   
     __cuda_hdev__ ~PBC_t() {}  // (no local allocations)

     __cuda_hdev__ const PBC_t & operator =  (const PBC_t   &pbc  );
     __cuda_hdev__ const PBC_t & operator =  (const Param_t *param); 

     __cuda_hdev__ void  Reposition     (const real8 *p1, real8 *p2)                                        const;
     __cuda_hdev__ void  Reposition     (const real8 *p1, real8 *p2, real8 *p3, real8 *p4)                  const;
     __cuda_hdev__ void  Reposition     (const real8 *p1, real8 *p2, real8 *p3, real8 *p4, const real8 *pc) const;

     __cuda_hdev__ void  FoldBox        (real8 *x, real8 *y, real8 *z) const;
     __cuda_hdev__ void  FoldBox        (real8 &x, real8 &y, real8 &z) const;
     __cuda_hdev__ void  FoldBox        (real8 *v)                     const;

     __cuda_hdev__ void  Reposition     (const real8 x0, const real8 y0, const real8 z0, real8 *x, real8 *y, real8 *z) const;
     __cuda_hdev__ void  Reposition     (const real8 x0, const real8 y0, const real8 z0, real8 &x, real8 &y, real8 &z) const;

     __cuda_hdev__ void  ZImage         (real8 *x, real8 *y, real8 *z) const;
     __cuda_hdev__ void  ZImage         (real8 &x, real8 &y, real8 &z) const;
     __cuda_hdev__ void  ZImage         (real8 *v)                     const;

};

// PBC_Position()
// 
// All the routines below are various overloaded version of the PBC positioning function.
// When two endpoints are provided (p1,p2) the second endpoint is moved to the closest
// periodic image of the position with respect to the first endpoint (p1).
//
// When called with 4 endpoints (p1,p2,p3,p4) the endpoints to (p2,p3,p4) are moved to
// the closest periodic image with respect to the first endpoint (p1).  Note that the 
// resulting endpoints may be located outside the primary simulation box and in one of 
// the nearby periodic images.  Further note that the images of the second segment (p3,p4)
// are moved with respect to the distance between the midpoints of the two segments.
//
// When called with 4 endpoints (p1,p2,p3,p4) and a reference (pc) the endpoints are 
// moved with respect to the center of the cell containing the first endpoint (p1).
//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void PBC_Position
(
   const PBC_t *pbc, ///< points to PBC structure
   const real8 *p1 , ///< position of node 1 <xyz> (unchanged)
         real8 *p2   ///< position of node 2 <xyz> (moved to nearest image)
);

__cuda_hdev__
void PBC_Position
(
   const PBC_t *pbc, ///< points to PBC structure
   const real8 *p1 , ///< position of node 1 <xyz> (unchanged)
         real8 *p2 , ///< position of node 2 <xyz> (moved to nearest image)
         real8 *p3 , ///< position of node 3 <xyz> (moved to nearest image)
         real8 *p4   ///< position of node 4 <xyz> (moved to nearest image)
);

__cuda_hdev__
void PBC_Position
(
   const PBC_t *pbc, ///< points to PBC structure
   const real8 *p1 , ///< position of node 1 <xyz> (unchanged)
         real8 *p2 , ///< position of node 2 <xyz> (moved to nearest image)
         real8 *p3 , ///< position of node 3 <xyz> (moved to nearest image)
         real8 *p4 , ///< position of node 4 <xyz> (moved to nearest image)
   const real8 *pc   ///< center of cell containing segment 1 (p1,p2)
);

__cuda_hdev__
void PBC_Position
(
   const real8 *p1,  ///< position of node 1 <xyz> (unchanged)
         real8 *p2,  ///< position of node 2 <xyz> (moved to nearest image)
   const real8  lx,  ///< simulation box size (x)
   const real8  ly,  ///< simulation box size (y)
   const real8  lz   ///< simulation box size (z)
                     ///< NB - this version assumes PBC always active!
);

__cuda_hdev__
void PBC_Position
(
   const real8 *p1,  ///< position of node 1 <xyz> (unchanged)
         real8 *p2,  ///< position of node 2 <xyz> (moved to nearest image)
   const real8  lx,  ///< simulation box size (x)
   const real8  ly,  ///< simulation box size (y)
   const real8  lz,  ///< simulation box size (z)
   const real8  sx,  ///< reciprocal of simulation box size (1.0/lx) (zero if PBC not active)
   const real8  sy,  ///< reciprocal of simulation box size (1.0/ly) (zero if PBC not active)
   const real8  sz   ///< reciprocal of simulation box size (1.0/lz) (zero if PBC not active)
);

__cuda_hdev__
void PBC_Position
(
   const real8 *p1,  ///< position of node 1 <xyz> (unchanged)
         real8 *p2,  ///< position of node 2 <xyz> (moved to nearest image)
         real8 *p3,  ///< position of node 3 <xyz> (moved to nearest image)
         real8 *p4,  ///< position of node 4 <xyz> (moved to nearest image)
   const real8  lx,  ///< simulation box size (x)
   const real8  ly,  ///< simulation box size (y)
   const real8  lz,  ///< simulation box size (z)
   const real8  sx,  ///< reciprocal of simulation box size (1.0/lx) (zero if PBC not active)
   const real8  sy,  ///< reciprocal of simulation box size (1.0/ly) (zero if PBC not active)
   const real8  sz   ///< reciprocal of simulation box size (1.0/lz) (zero if PBC not active)
);

__cuda_hdev__
void PBC_Position
(
   const real8 *p1,  ///< position of node 1 <xyz> (unchanged)
         real8 *p2,  ///< position of node 2 <xyz> (moved to nearest image)
         real8 *p3,  ///< position of node 3 <xyz> (moved to nearest image)
         real8 *p4,  ///< position of node 4 <xyz> (moved to nearest image)
   const real8 *pc,  ///< center of cell containing segment 1 (p1,p2)
   const real8  lx,  ///< simulation box size (x)
   const real8  ly,  ///< simulation box size (y)
   const real8  lz,  ///< simulation box size (z)
   const real8  sx,  ///< reciprocal of simulation box size (1.0/lx) (zero if PBC not active)
   const real8  sy,  ///< reciprocal of simulation box size (1.0/ly) (zero if PBC not active)
   const real8  sz   ///< reciprocal of simulation box size (1.0/lz) (zero if PBC not active)
);

// PBC_FoldBox()
//
// The FoldBox routines will move a given endpoint location from outside the simulation box 
// to its physical position inside the simulation box. Endpoints that are already within the
// simulation box are unchanged. Multiple overloaded versions are available.
//
// Note that the fold box routines use the simulation center as the anchor point for the
// adjustments.
//-----------------------------------------------------------------------------------------------

__cuda_hdev__ void PBC_FoldBox  (const PBC_t *pbc, real8 *x, real8 *y, real8 *z);   ///< repositions endpoint <xyz>
__cuda_hdev__ void PBC_FoldBox  (const PBC_t *pbc, real8 &x, real8 &y, real8 &z);   ///< repositions endpoint <xyz>
__cuda_hdev__ void PBC_FoldBox  (const PBC_t *pbc, real8 *v);                       ///< repositions endpoint vector <xyz>

// PBC_ZImage()
//
//
//-----------------------------------------------------------------------------------------------

__cuda_hdev__ void PBC_ZImage   (const PBC_t *pbc, real8 *x, real8 *y, real8 *z);
__cuda_hdev__ void PBC_ZImage   (const PBC_t *pbc, real8 &x, real8 &y, real8 &z);
__cuda_hdev__ void PBC_ZImage   (const PBC_t *pbc, real8 *v);

// PBC_Position()
//
// Given the position of two endpoints, will update the position of the second endpoint
// to the location of the nearest image with respect to the first.
//
//    <x0,y0,z0>  : reference position
//    <x,y,z>     : position to be updated
//    <p0>        : reference position     <xyz>
//    <p1>        : position to be updated <xyz>
//-----------------------------------------------------------------------------------------------

__cuda_hdev__ void PBC_Position (const PBC_t *pbc, const real8 x0, const real8 y0, const real8 z0, real8 *x, real8 *y, real8 *z);
__cuda_hdev__ void PBC_Position (const PBC_t *pbc, const real8 x0, const real8 y0, const real8 z0, real8 &x, real8 &y, real8 &z);
__cuda_hdev__ void PBC_Position (const PBC_t *pbc, const real8 *p0, real8 *p1);

#endif  // _PDS_PBC_H
