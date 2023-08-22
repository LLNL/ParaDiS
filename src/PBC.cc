#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"
#include "cuda_portability.h"

#include "Typedefs.h"
#include "Param.h"
#include "PBC.h"

//------------------------------------------------------------------------------------------------------------

#define __inline__ inline __attribute__((__always_inline__))

//------------------------------------------------------------------------------------------------------------
// This module implements several Periodic Boundary Correction (PBC) modules used by the GPU implemention
// of forces. There are several versions depending on whether Cell information is available or not.
//------------------------------------------------------------------------------------------------------------

// constructors...
//------------------------------------------------------------------------------------------------------------

__cuda_hdev__  
PBC_t::PBC_t(void)
{
   bndx=bndy=bndz=Periodic;

   x0=y0=z0=0.0;
   x1=y1=z1=0.0;
   cx=cy=cz=0.0;
   lx=ly=lz=0.0;
   sx=sy=sz=0.0;
}

__cuda_hdev__ PBC_t::PBC_t
(
   const BoundType_t _bndx, const BoundType_t _bndy, const BoundType_t _bndz,
   const real8       _x0  , const real8       _y0,   const real8       _z0  ,
   const real8       _x1  , const real8       _y1,   const real8       _z1
)
{
      bndx = _bndx;
      bndy = _bndy;
      bndz = _bndz;

      x0   = _x0;
      y0   = _y0;
      z0   = _z0;

      x1   = _x1;
      y1   = _y1;
      z1   = _z1;

      cx   = (x1+x0)/2.0;
      cy   = (y1+y0)/2.0;
      cz   = (z1+z0)/2.0;

      lx   = fabs(x1-x0);
      ly   = fabs(y1-y0);
      lz   = fabs(z1-z0);

      sx   = ( ((bndx == Periodic) && (lx>0.0)) ? (1.0/lx) : 0.0 );
      sy   = ( ((bndy == Periodic) && (ly>0.0)) ? (1.0/ly) : 0.0 );
      sz   = ( ((bndz == Periodic) && (lz>0.0)) ? (1.0/lz) : 0.0 );
}

__cuda_hdev__ PBC_t::PBC_t (const PBC_t   &pbc  ) { *this = pbc;   }
__cuda_hdev__ PBC_t::PBC_t (const Param_t *param) { *this = param; }

__cuda_hdev__ 
const PBC_t & PBC_t::operator = (const PBC_t & pbc)
{
   bndx = pbc.bndx;
   bndy = pbc.bndy;
   bndz = pbc.bndz;

   x0   = pbc.x0;
   y0   = pbc.y0;
   z0   = pbc.z0;

   x1   = pbc.x1;
   y1   = pbc.y1;
   z1   = pbc.z1;

   cx   = pbc.cx;
   cy   = pbc.cy;
   cz   = pbc.cz;

   lx   = pbc.lx;
   ly   = pbc.ly;
   lz   = pbc.lz;

   sx   = pbc.sx;
   sy   = pbc.sy;
   sz   = pbc.sz;

   return(*this);
}

__cuda_hdev__ 
const PBC_t & PBC_t::operator = (const Param_t *param)
{
   if (param)
   {
      bndx = param->xBoundType;
      bndy = param->yBoundType;
      bndz = param->zBoundType;

      x0   = param->minSideX;
      y0   = param->minSideY;
      z0   = param->minSideZ;

      x1   = param->maxSideX;
      y1   = param->maxSideY;
      z1   = param->maxSideZ;

      cx   = (x1+x0)/2.0;
      cy   = (y1+y0)/2.0;
      cz   = (z1+z0)/2.0;

      lx   = fabs(x1-x0);
      ly   = fabs(y1-y0);
      lz   = fabs(z1-z0);

      sx   = ( ((bndx == Periodic) && (lx>0.0)) ? (1.0/lx) : 0.0 );
      sy   = ( ((bndy == Periodic) && (ly>0.0)) ? (1.0/ly) : 0.0 );
      sz   = ( ((bndz == Periodic) && (lz>0.0)) ? (1.0/lz) : 0.0 );
   }

   return(*this);
}

// Reposition()
//
// Updates position 2 to the position of the nearest image of p2 to p1.
//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void PBC_t::Reposition
(
   const real8 *p1,   ///< position of node 1 <xyz> (unchanged)
         real8 *p2    ///< position of node 2 <xyz> (moved to nearest image)
) const
{
   p2[0] -= rint((p2[0]-p1[0])*sx)*lx;   // adjust p2 relative to p1
   p2[1] -= rint((p2[1]-p1[1])*sy)*ly;   //
   p2[2] -= rint((p2[2]-p1[2])*sz)*lz;   //
}

// Reposition()
//
// Updates the endpoints both segments based on the midpoints of the two segments.
// Note that the anchor for the system is p1. Endpoints p2 and p4 are adjusted to
// the nearest images to p1 and p3 (respectively), then both endpoints of segment 2
// (p3,p4)are adjusted to the nearest image of segment 1 (p1,p2) using the
// midpoints of the two segments.
//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void PBC_t::Reposition
(
   const real8 *p1,   ///< position of node 1 <xyz> (anchor, unchanged)
         real8 *p2,   ///< position of node 2 <xyz> (moved to nearest image)
         real8 *p3,   ///< position of node 3 <xyz> (moved to nearest image)
         real8 *p4    ///< position of node 4 <xyz> (moved to nearest image)
) const
{
   p2[0] -= rint((p2[0]-p1[0])*sx)*lx;           // adjust p2 relative to p1
   p2[1] -= rint((p2[1]-p1[1])*sy)*ly;           //
   p2[2] -= rint((p2[2]-p1[2])*sz)*lz;           //

   p4[0] -= rint((p4[0]-p3[0])*sx)*lx;           // adjust p4 relative to p3
   p4[1] -= rint((p4[1]-p3[1])*sy)*ly;           //
   p4[2] -= rint((p4[2]-p3[2])*sz)*lz;           //

   real8 pa[3] = { (p1[0]+p2[0])/2.0,            // pa = midpoint of segment 1 (p1->p2)
                   (p1[1]+p2[1])/2.0,            //
                   (p1[2]+p2[2])/2.0 };          //

   real8 pb[3] = { (p3[0]+p4[0])/2.0,            // pb = midpoint of segment 2 (p3->p4)
                   (p3[1]+p4[1])/2.0,            //
                   (p3[2]+p4[2])/2.0 };          //

   real8 pc[3] = { -rint((pb[0]-pa[0])*sx)*lx,   // pc = pbc correction for segment 2 (p3->p4)
                   -rint((pb[1]-pa[1])*sy)*ly,   //
                   -rint((pb[2]-pa[2])*sz)*lz }; //

   p3[0] += pc[0];                               // adjust p3
   p3[1] += pc[1];                               //
   p3[2] += pc[2];                               //

   p4[0] += pc[0];                               // adjust p4
   p4[1] += pc[1];                               //
   p4[2] += pc[2];                               //
}

//------------------------------------------------------------------------------------------------------------

__cuda_hdev__
void PBC_t::Reposition
(
   const real8 *p1,   ///< position of node 1 <xyz> (unchanged)
         real8 *p2,   ///< position of node 2 <xyz> (moved to nearest image)
         real8 *p3,   ///< position of node 3 <xyz> (moved to nearest image)
         real8 *p4,   ///< position of node 4 <xyz> (moved to nearest image)
   const real8 *pc    ///< center of cell containing segment 1 (p1,p2)
) const
{
   p2[0] -= rint((p2[0]-p1[0])*sx)*lx;   // adjust p2 relative to p1
   p2[1] -= rint((p2[1]-p1[1])*sy)*ly;   //
   p2[2] -= rint((p2[2]-p1[2])*sz)*lz;   //

   p3[0] -= rint((p3[0]-pc[0])*sx)*lx;   // adjust p3 relative to cell center
   p3[1] -= rint((p3[1]-pc[1])*sy)*ly;   //
   p3[2] -= rint((p3[2]-pc[2])*sz)*lz;   //

   p4[0] -= rint((p4[0]-p3[0])*sx)*lx;   // adjust p4 relative to p3
   p4[1] -= rint((p4[1]-p3[1])*sy)*ly;   //
   p4[2] -= rint((p4[2]-p3[2])*sz)*lz;   //
}

// FoldBox()
//
// Will reposition a given position <xyz> to a point within the primary simulation box
// Note that the adjustment is made with respect to the simulation center (cx,cy,cz).
//------------------------------------------------------------------------------------------------------------

__cuda_hdev__
void PBC_t::FoldBox (real8 *x, real8 *y, real8 *z) const
{
   *x -= rint(((*x)-cx)*sx)*lx;
   *y -= rint(((*y)-cy)*sy)*ly;
   *z -= rint(((*z)-cz)*sz)*lz;
}

__cuda_hdev__
void PBC_t::FoldBox (real8 &x, real8 &y, real8 &z) const
{
   x -= rint((x-cx)*sx)*lx;
   y -= rint((y-cy)*sy)*ly;
   z -= rint((z-cz)*sz)*lz;
}

__cuda_hdev__
void PBC_t::FoldBox (real8 *v) const
{
   v[0] -= rint((v[0]-cx)*sx)*lx;
   v[1] -= rint((v[1]-cy)*sy)*ly;
   v[2] -= rint((v[2]-cz)*sz)*lz;
}

//------------------------------------------------------------------------------------------------------------

__cuda_hdev__
void  PBC_t::Reposition
(
   const real8  xref, const real8  yref, const real8  zref,    ///< reference position <xyz>
         real8 *x   ,       real8 *y   ,       real8 *z        ///< position to be updated <xyz>
) const
{
   *x -= rint(((*x)-xref)*sx)*lx;
   *y -= rint(((*y)-yref)*sy)*ly;
   *z -= rint(((*z)-zref)*sz)*lz;
}

__cuda_hdev__
void  PBC_t::Reposition
(
   const real8  xref, const real8  yref, const real8  zref,    ///< reference position <xyz>
         real8 &x   ,       real8 &y   ,       real8 &z        ///< position to be updated <xyz>
) const
{
   x -= rint((x-xref)*sx)*lx;
   y -= rint((y-yref)*sy)*ly;
   z -= rint((z-zref)*sz)*lz;
}

// ZImage()
//------------------------------------------------------------------------------------------------------------

__cuda_hdev__
void PBC_t::ZImage (real8 *x, real8 *y, real8 *z) const
{
   *x -= rint((*x)*sx)*lx;
   *y -= rint((*y)*sy)*ly;
   *z -= rint((*z)*sz)*lz;
}

__cuda_hdev__
void PBC_t::ZImage(real8 &x, real8 &y, real8 &z) const
{
   x -= rint(x*sx)*lx;
   y -= rint(y*sy)*ly;
   z -= rint(z*sz)*lz;
}

__cuda_hdev__
void PBC_t::ZImage (real8 *v) const
{
   v[0] -= rint(v[0]*sx)*lx;
   v[1] -= rint(v[1]*sy)*ly;
   v[2] -= rint(v[2]*sz)*lz;
}

// PBC_Position()
//
// Will update the position of the second node to the position of the nearest image with respect
// to the first node.
//------------------------------------------------------------------------------------------------------------

__inline__ __cuda_hdev__ real8 PBC_Position(const real8 x, const real8 xr, const real8 lx                ) { return(x - rint((x-xr)/lx)*lx ); }
__inline__ __cuda_hdev__ real8 PBC_Position(const real8 x, const real8 xr, const real8 lx, const real8 sx) { return(x - rint((x-xr)*sx)*lx ); }

__cuda_hdev__
void PBC_Position
(
   const real8 *p1,  ///< position of node 1 <xyz> (unchanged)
         real8 *p2,  ///< position of node 2 <xyz> (moved to nearest image)
   const real8  lx,  ///< simulation box size (x)
   const real8  ly,  ///< simulation box size (y)
   const real8  lz   ///< simulation box size (z)
                     ///< NB - this version assumes PBC always active!
)
{
                                          // adjust p2 relative to p1
   p2[0] = PBC_Position(p2[0],p1[0],lx);  // (assumes PBC active)
   p2[1] = PBC_Position(p2[1],p1[1],ly);  // (assumes PBC active)
   p2[2] = PBC_Position(p2[2],p1[2],lz);  // (assumes PBC active)
}

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
)
{
   p2[0] = PBC_Position(p2[0],p1[0],lx,sx);  // adjust p2 relative to p1
   p2[1] = PBC_Position(p2[1],p1[1],ly,sy);  //
   p2[2] = PBC_Position(p2[2],p1[2],lz,sz);  //
}

// PBC_Position()
//
// Will update the positions of nodes 2,3,4 as necessary for periodic boundary conditions.
// The position of node 1 is unchanged.
//
// This version adjusts based on the midpoint of the two segments.
//------------------------------------------------------------------------------------------------------------

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
)
{
   p2[0] = PBC_Position(p2[0],p1[0],lx,sx);      // adjust p2 relative to p1
   p2[1] = PBC_Position(p2[1],p1[1],ly,sy);      //
   p2[2] = PBC_Position(p2[2],p1[2],lz,sz);      //

   p4[0] = PBC_Position(p4[0],p3[0],lx,sx);      // adjust p4 relative to p3
   p4[1] = PBC_Position(p4[1],p3[1],ly,sy);      //
   p4[2] = PBC_Position(p4[2],p3[2],lz,sz);      //

   real8 pa[3] = { (p1[0]+p2[0])/2.0,            // pa = midpoint of segment 1
                   (p1[1]+p2[1])/2.0,            //
                   (p1[2]+p2[2])/2.0 };          //

   real8 pb[3] = { (p3[0]+p4[0])/2.0,            // pb = midpoint of segment 2
                   (p3[1]+p4[1])/2.0,            //
                   (p3[2]+p4[2])/2.0 };          //

   real8 pc[3] = { -rint((pb[0]-pa[0])*sx)*lx,   // pc = pbc correction for segment 2
                   -rint((pb[1]-pa[1])*sy)*ly,   //
                   -rint((pb[2]-pa[2])*sz)*lz }; //

   p3[0] += pc[0];                               // apply pbc correction to p3 (segment 2)
   p3[1] += pc[1];                               //
   p3[2] += pc[2];                               //

   p4[0] += pc[0];                               // apply pbc correction to p4 (segment 2)
   p4[1] += pc[1];                               //
   p4[2] += pc[2];                               //
}

// PBC_Position()
//
// Will update the positions of nodes 2,3,4 as necessary for periodic boundary conditions.
// The position of node 1 is unchanged.
//
// This version adjusts based on the center of the cell assigned to segment 1 (p1,p2).
//------------------------------------------------------------------------------------------------------------

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
)
{
   p2[0] = PBC_Position(p2[0],p1[0],lx,sx);   // adjust p2 relative to p1
   p2[1] = PBC_Position(p2[1],p1[1],ly,sy);   //
   p2[2] = PBC_Position(p2[2],p1[2],lz,sz);   //

   p3[0] = PBC_Position(p3[0],pc[0],lx,sx);   // adjust p3 relative to cell center
   p3[1] = PBC_Position(p3[1],pc[1],ly,sy);   //
   p3[2] = PBC_Position(p3[2],pc[2],lz,sz);   //

   p4[0] = PBC_Position(p4[0],p3[0],lx,sx);   // adjust p4 relative to p3
   p4[1] = PBC_Position(p4[1],p3[1],ly,sy);   //
   p4[2] = PBC_Position(p4[2],p3[2],lz,sz);   //
}

//------------------------------------------------------------------------------------------------------------
// Function:     FoldBox
// Description:  Given a set of coordinates, adjusts the
//               coordinates to the corresponding point within
//               the primary (non-periodic) image of the problem
//               space.  If the provided coordinates are already
//               within the primary image, no adjustments are made.
//
// Arguments:
//    param       Pointer to structure containing control parameters
//    x, y, z     Pointers to components of coordinate.  On return
//                to the caller, these coordinates will have been
//                adjusted (if necessary) to the coordinates of the
//                corresponding point within the primary image of
//                the problem space.
//
//                Note that if periodic boundaries are not in use, the position will not be updated.
//------------------------------------------------------------------------------------------------------------

__inline__
__cuda_hdev__
real8 FoldBox(const real8 x, const real8 xmin, const real8 xmax, const BoundType_t bndx)
{
   const real8 cx = (xmax+xmin)/2.0;                                       // cx = center of simulation box
   const real8 lx = fabs(xmax-xmin);                                       // lx = simulation box size
   const real8 sx = ( ((bndx==Periodic) && (lx>0.0)) ? (1.0/lx) : 0.0 );   // sx = reciprocol of simulation box size, or zero (if PBC not enabled)

   return( x - rint((x-cx)*sx)*lx );
}

void FoldBox
(
   const Param_t *param,                   ///< points to current control parameters
         real8   *x, real8 *y, real8 *z    ///< position to be updated <xyx>
)
{
   *x = FoldBox( *x, param->minSideX, param->maxSideX, param->xBoundType );
   *y = FoldBox( *y, param->minSideY, param->maxSideY, param->yBoundType );
   *z = FoldBox( *z, param->minSideZ, param->maxSideZ, param->zBoundType );
}

void FoldBox
(
   const Param_t *param,                   ///< points to current control parameters
         real8   &x, real8 &y, real8 &z    ///< position to be updated <xyx>
)
{
   x = FoldBox( x, param->minSideX, param->maxSideX, param->xBoundType );
   y = FoldBox( y, param->minSideY, param->maxSideY, param->yBoundType );
   z = FoldBox( z, param->minSideZ, param->maxSideZ, param->zBoundType );
}

void FoldBox
(
   const Param_t *param,    ///< points to current control parameters
         real8   *p         ///< position to be updated <xyx>
)
{
   p[0] = FoldBox( p[0], param->minSideX, param->maxSideX, param->xBoundType );
   p[1] = FoldBox( p[1], param->minSideY, param->maxSideY, param->yBoundType );
   p[2] = FoldBox( p[2], param->minSideZ, param->maxSideZ, param->zBoundType );
}

__inline__
__cuda_hdev__
real8 PBC_FoldBox(const real8 x, const real8 cx, const real8 lx, const real8 sx)
{
   return( x - rint((x-cx)*sx)*lx );
}

__cuda_hdev__
void PBC_FoldBox (const PBC_t *pbc, real8 *x, real8 *y, real8 *z)
{
   *x = PBC_FoldBox(*x, pbc->cx, pbc->lx, pbc->sx );
   *y = PBC_FoldBox(*y, pbc->cy, pbc->ly, pbc->sy );
   *z = PBC_FoldBox(*z, pbc->cz, pbc->lz, pbc->sz );
}

__cuda_hdev__
void PBC_FoldBox (const PBC_t *pbc, real8 &x, real8 &y, real8 &z)
{
   x = PBC_FoldBox(x, pbc->cx, pbc->lx, pbc->sx );
   y = PBC_FoldBox(y, pbc->cy, pbc->ly, pbc->sy );
   z = PBC_FoldBox(z, pbc->cz, pbc->lz, pbc->sz );
}

__cuda_hdev__
void PBC_FoldBox (const PBC_t *pbc, real8 *v)
{
   v[0] = PBC_FoldBox(v[0], pbc->cx, pbc->lx, pbc->sx );
   v[1] = PBC_FoldBox(v[1], pbc->cy, pbc->ly, pbc->sy );
   v[2] = PBC_FoldBox(v[2], pbc->cz, pbc->lz, pbc->sz );
}


//------------------------------------------------------------------------------------------------------------
// Function:      PBCPOSITION
// Description:   Finds the position of the image of point (x,y,z)
//                that is closest to point (x0,y0,z0).
//
//                NOTE: The returned position is not required
//                to be in the primary image but may actually be in
//                one of the periodic images
//
// Arguments:
//    param       Pointer to structure containing control parameters
//    x0, y0, z0  Components of position used as a base.
//    x, y, z     Pointers to components of secondary point.  These
//                values will be overwritten with coordinates of the
//                image of this point closest to (x0,y0,z0).
//
//                Note that if periodic boundaries are not in use, the position will not be updated.
//------------------------------------------------------------------------------------------------------------

__inline__
__cuda_hdev__
real8 PBC_Position(const real8 x, const real8 xr, const real8 xmin, const real8 xmax, const BoundType_t bndx)
{
   const real8 lx = fabs(xmax-xmin);                                        // lx = simulation box size
   const real8 sx = ( ((bndx==Periodic) && (lx>0.0)) ? (1.0/lx) : 0.0 );    // sx = reciprocol of simulation box size, or zero (if PBC not enabled)

   return( x - rint((x-xr)*sx)*lx );
}

void PBCPOSITION
(
   const Param_t *param,                                       ///< points to control parameters
   const real8  xref, const real8  yref, const real8  zref,    ///< reference position     <xyz>
         real8 *x   ,       real8 *y   ,       real8 *z        ///< position to be updated <xyz>
)
{
   *x = PBC_Position( *x, xref, param->minSideX, param->maxSideX, param->xBoundType );
   *y = PBC_Position( *y, yref, param->minSideY, param->maxSideY, param->yBoundType );
   *z = PBC_Position( *z, zref, param->minSideZ, param->maxSideZ, param->zBoundType );
}

void PBCPOSITION
(
   const Param_t *param,                                       ///< points to control parameters
   const real8  xref, const real8  yref, const real8  zref,    ///< reference position     <xyz>
         real8 &x   ,       real8 &y   ,       real8 &z        ///< position to be updated <xyz>
)
{
   x = PBC_Position( x, xref, param->minSideX, param->maxSideX, param->xBoundType );
   y = PBC_Position( y, yref, param->minSideY, param->maxSideY, param->yBoundType );
   z = PBC_Position( z, zref, param->minSideZ, param->maxSideZ, param->zBoundType );
}

void PBCPOSITION
(
   const Param_t *param,   ///< points to control parameters
   const real8   *vref ,   ///< reference position     <xyz>
         real8   *v        ///< position to be updated <xyz>
)
{
   v[0] = PBC_Position( v[0], vref[0], param->minSideX, param->maxSideX, param->xBoundType );
   v[1] = PBC_Position( v[1], vref[1], param->minSideY, param->maxSideY, param->yBoundType );
   v[2] = PBC_Position( v[2], vref[2], param->minSideZ, param->maxSideZ, param->zBoundType );
}

__cuda_hdev__
void PBC_Position (const PBC_t *pbc, const real8 xref, const real8 yref, const real8 zref, real8 *x, real8 *y, real8 *z)
{
   *x = PBC_Position( *x, xref, pbc->lx, pbc->sx );
   *y = PBC_Position( *y, yref, pbc->ly, pbc->sy );
   *z = PBC_Position( *z, zref, pbc->lz, pbc->sz );
}

__cuda_hdev__
void PBC_Position (const PBC_t *pbc, const real8 xref, const real8 yref, const real8 zref, real8 &x, real8 &y, real8 &z)
{
   x = PBC_Position( x, xref, pbc->lx, pbc->sx );
   y = PBC_Position( y, yref, pbc->ly, pbc->sy );
   z = PBC_Position( z, zref, pbc->lz, pbc->sz );
}

__cuda_hdev__
void PBC_Position (const PBC_t *pbc, const real8 *vref, real8 *v)
{
   v[0] = PBC_Position( v[0], vref[0], pbc->lx, pbc->sx );
   v[1] = PBC_Position( v[1], vref[1], pbc->ly, pbc->sy );
   v[2] = PBC_Position( v[2], vref[2], pbc->lz, pbc->sz );
}

//------------------------------------------------------------------------------------------------------------
// Function:     ZImage
// Description:  Finds the minimum image of (x,y,z) in the problem box.
//
//               Typical use of this function specifies (x,y,z)
//               as vector from a source point to a secondary
//               point.  Upon return to the caller, the vector
//               will have been adjusted to be avector from the
//               source point to the closest image of the secondary
//               point.
//
// Arguments:
//    param      Pointer to structure containing control parameters
//    x, y, z    Pointers to components of point (or vector).
//
//               Note that if periodic boundaries are not in use, the position will not be updated.
//------------------------------------------------------------------------------------------------------------

__inline__
__cuda_hdev__
real8 PBC_ZImage(const real8 x, const real8 xmin, const real8 xmax, const BoundType_t bndx )
{
   const real8 lx = fabs(xmax-xmin);                                        // lx = simulation box size
   const real8 sx = ( ((bndx==Periodic) && (lx>0.0)) ? (1.0/lx) : 0.0 );    // sx = reciprocol of simulation box size, or zero (if PBC not enabled)

   return( x - rint(x*sx)*lx );
}

void ZImage
(
   const Param_t *param,                  ///< points to control parameters
         real8 *x, real8 *y, real8 *z     ///< position to be updated <xyz>
)
{
   *x = PBC_ZImage( *x, param->minSideX, param->maxSideX, param->xBoundType );
   *y = PBC_ZImage( *y, param->minSideY, param->maxSideY, param->yBoundType );
   *z = PBC_ZImage( *z, param->minSideZ, param->maxSideZ, param->zBoundType );
}

void ZImage
(
   const Param_t *param,                  ///< points to control parameters
         real8 &x, real8 &y, real8 &z     ///< position to be updated <xyz>
)
{
   x = PBC_ZImage( x, param->minSideX, param->maxSideX, param->xBoundType );
   y = PBC_ZImage( y, param->minSideY, param->maxSideY, param->yBoundType );
   z = PBC_ZImage( z, param->minSideZ, param->maxSideZ, param->zBoundType );
}

void ZImage
(
   const Param_t *param,    ///< points to control parameters
         real8   *p         ///< position to be updated <xyz>
)
{
   p[0] = PBC_ZImage( p[0], param->minSideX, param->maxSideX, param->xBoundType );
   p[1] = PBC_ZImage( p[1], param->minSideY, param->maxSideY, param->yBoundType );
   p[2] = PBC_ZImage( p[2], param->minSideZ, param->maxSideZ, param->zBoundType );
}

//------------------------------------------------------------------------------------------------------------

__inline__
__cuda_hdev__
real8 PBC_ZImage(const real8 x, const real8 lx, const real8 sx)
{
   return( x - rint(x*sx)*lx );
}

__cuda_hdev__
void PBC_ZImage (const PBC_t *pbc, real8 *x, real8 *y, real8 *z)
{
   *x = PBC_ZImage( *x, pbc->lx, pbc->sx );
   *y = PBC_ZImage( *y, pbc->ly, pbc->sy );
   *z = PBC_ZImage( *z, pbc->lz, pbc->sz );
}

__cuda_hdev__
void PBC_ZImage (const PBC_t *pbc, real8 &x, real8 &y, real8 &z)
{
   x = PBC_ZImage( x, pbc->lx, pbc->sx );
   y = PBC_ZImage( y, pbc->ly, pbc->sy );
   z = PBC_ZImage( z, pbc->lz, pbc->sz );
}

__cuda_hdev__
void PBC_ZImage (const PBC_t *pbc, real8 *v)
{
   v[0] = PBC_ZImage( v[0], pbc->lx, pbc->sx );
   v[1] = PBC_ZImage( v[1], pbc->ly, pbc->sy );
   v[2] = PBC_ZImage( v[2], pbc->lz, pbc->sz );
}

