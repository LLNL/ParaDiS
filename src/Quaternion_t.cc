#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Typedefs.h"
#include "Quaternion_t.h"
#include "V3.h"


Quaternion_t Quaternion_t::operator *  (const Quaternion_t & q)
{
   real8 va[3] = {   qx,   qy,   qz };
   real8 vb[3] = { q.qx, q.qy, q.qz };
   real8 vx[3] = {  0.0,  0.0,  0.0 };

   V3_CROSS(vx,va,vb);

   vx[0] += va[0]*q.qs;
   vx[1] += va[1]*q.qs;
   vx[2] += va[2]*q.qs;

   vx[0] += vb[0]*  qs;
   vx[1] += vb[1]*  qs;
   vx[2] += vb[2]*  qs;

   return ( Quaternion_t(vx[0], vx[1], vx[2], (q.qs*qs)-V3_DOT(va,vb) ) );
}

int Quaternion_t::operator == (const Quaternion_t  & q) const
{
   return (     ( fabs(qx-q.qx) > 1.0e-8 )
             || ( fabs(qy-q.qy) > 1.0e-8 )
             || ( fabs(qz-q.qz) > 1.0e-8 )
             || ( fabs(qs-q.qs) > 1.0e-8 ) ? 0 : 1 );
}

int Quaternion_t::operator == (const real8           s) const
{
   return (     ( fabs(qx-s) > 1.0e-8 )
             || ( fabs(qy-s) > 1.0e-8 )
             || ( fabs(qz-s) > 1.0e-8 )
             || ( fabs(qs-s) > 1.0e-8 ) ? 0 : 1 );
}

int Quaternion_t::operator != (const Quaternion_t  & q) const
{
   return (     ( fabs(qx-q.qx) > 1.0e-8 )
             || ( fabs(qy-q.qy) > 1.0e-8 )
             || ( fabs(qz-q.qz) > 1.0e-8 )
             || ( fabs(qs-q.qs) > 1.0e-8 ) ? 1 : 0 );
}

int Quaternion_t::operator != (const real8           s) const
{
   return (     ( fabs(qx-s) > 1.0e-8 )
             || ( fabs(qy-s) > 1.0e-8 )
             || ( fabs(qz-s) > 1.0e-8 )
             || ( fabs(qs-s) > 1.0e-8 ) ? 1 : 0 );
}

