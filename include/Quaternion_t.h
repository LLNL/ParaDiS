#pragma once

#ifndef __PDS_QUATERNION_H
#define __PDS_QUATERNION_H

#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Typedefs.h"

class Quaternion_t
{
   public : 
      union
      {
         real8     qv[4];  ///< quaternion as 4-element vector (x,y,z,s)

         struct 
         {
            real8  qx;   ///< basis vector (x)
            real8  qy;   ///< basis vector (y)
            real8  qz;   ///< basis vector (z)
            real8  qs;   ///< scale (s)
         };
      };

   public : 
       Quaternion_t(void)                                                       { qx=0.0; qy=0.0; qz=0.0; qs=1.0; }
       Quaternion_t(const int   x, const int   y, const int   z, const int   s) { qx=(real8) x; qy=(real8) y; qz=(real8) z; qs=(real8) s; }
       Quaternion_t(const real8 x, const real8 y, const real8 z, const real8 s) { qx=x; qy=y; qz=z; qs=s; }
      ~Quaternion_t() {} // (nothing allocated)

      const Quaternion_t & operator = (const Quaternion_t & q) { qx=q.qx; qy=q.qy; qz=q.qz; qs=q.qs; return(*this); } 

      void operator += (const Quaternion_t & q) { qx += q.qx; qy += q.qy; qz += q.qz; qs += q.qs; } 
      void operator -= (const Quaternion_t & q) { qx -= q.qx; qy -= q.qy; qz -= q.qz; qs -= q.qs; } 
      void operator += (const real8          s) { qx += s   ; qy += s   ; qz += s   ; qs += s   ; } 
      void operator -= (const real8          s) { qx -= s   ; qy -= s   ; qz -= s   ; qs -= s   ; } 
      void operator *= (const real8          s) { qx *= s   ; qy *= s   ; qz *= s   ; qs *= s   ; } 
      void operator /= (const real8          s) { real8 t = ( (fabs(s)>0.0) ? (1.0/s) : 0.0 ); qx *= t; qy *= t; qz *= t; qs *= t; } 

      Quaternion_t operator +  (const Quaternion_t & q) { return( Quaternion_t( (qx+q.qx), (qy+q.qy), (qz+q.qz), (qs+q.qs) ) ); }
      Quaternion_t operator -  (const Quaternion_t & q) { return( Quaternion_t( (qx-q.qx), (qy-q.qy), (qz-q.qz), (qs-q.qs) ) ); }

      Quaternion_t operator *  (const Quaternion_t & q);
      Quaternion_t operator *  (const real8          s) { return( Quaternion_t( (s*qx*s), (s*qy), (s*qz), (s*qs) ) ); }

      Quaternion_t operator /  (const real8          s) { real8 t = ( (fabs(s)>0.0) ? (1.0/s) : 0.0 ); 
                                                          return( Quaternion_t( (t*qx), (t*qy), (t*qz), (t*qs) ) ); }

      Quaternion_t operator - () const { return( Quaternion_t(-qx,-qy,-qz,-qs) ); }

      int          operator == (const Quaternion_t  & q) const;
      int          operator == (const real8           s) const;
      int          operator != (const Quaternion_t  & q) const;
      int          operator != (const real8           s) const;

      void         Abs       (void)                         { qx=fabs(qx); qy=fabs(qy); qz=fabs(qz); qs=fabs(qs); } 
      real8        Dot       (const Quaternion_t & q) const { return( (qx*q.qx) + (qy*q.qy) + (qz*q.qz) + (qs*q.qs) ); }
      real8        Norm      (void)                   const { return( sqrt((qx*qx) + (qy*qy) + (qz*qz) + (qs*qs)) ); }
      real8        Mag       (void)                   const { return( sqrt((qx*qx) + (qy*qy) + (qz*qz) + (qs*qs)) ); }
      void         Normalize (void)                         { real8 s = sqrt((qx*qx) + (qy*qy) + (qz*qz) + (qs*qs)); 
                                                                    s = ( (fabs(s)>0.0) ? (1.0/s) : 0.0 ); qx*=s; qy*=s; qz*=s; qs*=s; }
};

inline Quaternion_t Abs       (const Quaternion_t & q) { return( Quaternion_t( fabs(q.qx), fabs(q.qy), fabs(q.qz), fabs(q.qs)) ); } 
inline real8        Norm      (const Quaternion_t & q) { return( sqrt((q.qx*q.qx) + (q.qy*q.qy) + (q.qz*q.qz) + (q.qs*q.qs)) ); }
inline real8        Mag       (const Quaternion_t & q) { return( sqrt((q.qx*q.qx) + (q.qy*q.qy) + (q.qz*q.qz) + (q.qs*q.qs)) ); }

inline Quaternion_t Normalize (const Quaternion_t & q) { real8 s = sqrt((q.qx*q.qx) + (q.qy*q.qy) + (q.qz*q.qz) + (q.qs*q.qs)); 
                                                               s = ( (fabs(s)>0.0) ? (1.0/s) : 0.0 ); 
                                                         return( Quaternion_t( (s*q.qx), (s*q.qy), (s*q.qz), (s*q.qs) ) ); }

inline real8        Dot       (const Quaternion_t & qa, const Quaternion_t & qb) { return( (qa.qx*qb.qx) + (qa.qy*qb.qy) + (qa.qz*qb.qz) + (qa.qs*qb.qs) ); }

#endif
