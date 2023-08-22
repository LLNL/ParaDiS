/**************************************************************************
 *  Module:  V3 - implements a simple 1x3 vector operations
 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "V3.h"
#include "M33.h"

/*===========================================================================================*/

int V3_Near (const V3 va, const V3 vb)
{ 
   if (va && vb) 
   {
      if (    (fabs(va[0]-vb[0])>1e-12) 
           || (fabs(va[1]-vb[1])>1e-12) 
           || (fabs(va[2]-vb[2])>1e-12) ) return(0);
      else
         return(1);
   }

   return(0); 
}

/*===========================================================================================*/

int V3_Near (const V3 va, const V3 vb, const real8 eps)
{ 
   if (va && vb) 
   {
      if (    (fabs(va[0]-vb[0])>eps)
           || (fabs(va[1]-vb[1])>eps)
           || (fabs(va[2]-vb[2])>eps) ) return(0);
      else
         return(1);
   }

   return(0); 
}

/*===========================================================================================*/

int V3_Near (const V3 v, const real8 x, const real8 y, const real8 z, const real8 eps) 
{ 
   if (v) 
   {
      if (    (fabs(v[0]-x)>eps) 
           || (fabs(v[1]-y)>eps) 
           || (fabs(v[2]-z)>eps) ) return(0);
      else
         return(1);
   }

   return(0); 
}

/*===========================================================================================*/

void V3_Normalize (real8 *v, const int vn)
{
   if (v && (vn>0))
   {
      for (int i=0; (i<vn); i++, v+=3)
         V3_NORMALIZE(v,v);
   } 
}

/*===========================================================================================*/

void V3_Scale (real8 *v, const int vn, const real8 s)
{
   if (v && (vn>0))
   {
      for (int i=0; (i<vn); i++, v+=3)
         V3_SCALE(v,v,s);
   } 
}

void V3_Scale (real8 *v, const int vn, const real8 sx, const real8 sy, const real8 sz)
{
   if (v && (vn>0))
   {
      for (int i=0; (i<vn); i++, v+=3)
         V3_SCALE_S3(v,v,sx,sy,sz);
   } 
}

// V3_Aligned()
//
// Returns whether two normalized vectors are aligned.
//-----------------------------------------------------------------------------------------------

int V3_Aligned (const real8 *v, const real8 x, const real8 y, const real8 z, const real8 eps)
{
   real8 va[3] = { v[0], v[1], v[2] };
   real8 vb[3] = { x, y, z };

   real8 sa = sqrt( va[0]*va[0] + va[1]*va[1] + va[2]*va[2] );   
         sa = ( (fabs(sa)>0.0) ? (1.0/sa) : 0.0 );

   real8 sb = sqrt( vb[0]*vb[0] + vb[1]*vb[1] + vb[2]*vb[2] );   
         sb = ( (fabs(sb)>0.0) ? (1.0/sb) : 0.0 );

         va[0]*=sa;  va[1]*=sa;  va[2]*=sa;
         vb[0]*=sb;  vb[1]*=sb;  vb[2]*=sb;

   return (    (fabs(va[0]-vb[0])>eps)
            || (fabs(va[1]-vb[1])>eps)
            || (fabs(va[2]-vb[2])>eps) ? 0 : 1 );
}

int V3_Aligned (const real8 *v0, const real8 *v1, const real8 eps)
{
   real8 va[3] = { v0[0], v0[1], v0[2] };
   real8 vb[3] = { v1[0], v1[1], v1[2] };

   real8 sa = sqrt( va[0]*va[0] + va[1]*va[1] + va[2]*va[2] );   
         sa = ( (fabs(sa)>0.0) ? (1.0/sa) : 0.0 );

   real8 sb = sqrt( vb[0]*vb[0] + vb[1]*vb[1] + vb[2]*vb[2] );   
         sb = ( (fabs(sb)>0.0) ? (1.0/sb) : 0.0 );

         va[0]*=sa;  va[1]*=sa;  va[2]*=sa;
         vb[0]*=sb;  vb[1]*=sb;  vb[2]*=sb;

   return (    (fabs(va[0]-vb[0])>eps)
            || (fabs(va[1]-vb[1])>eps)
            || (fabs(va[2]-vb[2])>eps) ? 0 : 1 );
}

/*-----------------------------------------------------------------------------------------------
 * Ray_Plane_Intersect()
 *
 * Calculates the intersection of the ray defined by PA-PB and the plane that includes three 
 * non-colinear points (P0,P1,and P2).
 *
 * Inputs,
 *   pa  - starting point of ray
 *   pb  - point on ray other than pa
 *   p0  - vertex of triangle on plane
 *   p1  - vertex of triangle on plane
 *   p2  - vertex of triangle on plane (note - p0,p1,and p2 cannot be collinear)
 *
 * Returns
 *   pt  - the intercept point
 *   tuv - the solution to,
 *
 *      |t|   | (xa-xb)  (x1-x0)  (x2-x0) |-1   |(xa-x0)|
 *      |u| = | (ya-yb)  (y1-y0)  (y2-y0) |   * |(ya-y0)|
 *      |v|   | (za-zb)  (z1-z0)  (z2-z0) |     |(za-z0)|
 *
 *  where,
 *    t - identifies the parametric that solves for the intercept. 
 *          t=[0-1] : intercept is between PA and PB
 *          t<0     : intercept is behind PA
 *          t>1     : intercept is beyond PB
 *  
 *  additionally, if...
 *    u=[0,1] and v=[0,1] and (u+v)<=1 : the intercept lies within the triangle specified by P0,P1, and P2.
 *
 * Note that the solution uses an 3x3 matrix for the inversion (by Cramer's rule).  I haven't 
 * figured out what happens in the degenerate case where there is no intercept because the 
 * ray is parallel to the plane or the caller has provided bad input (PA==PB or P0-P2 collinear).
 *
 * reference : http://en.wikipedia.org/wiki/Line-plane_intersection
 *-----------------------------------------------------------------------------------------------*/

void V3_Ray_Plane_Intersect (
   const V3 pa, const V3 pb, 
   const V3 p0, const V3 p1, const V3 p2, 
   V3 pt,  
   V3 tuv )
{
   V3    pr; 
   M33   m33, m33i;

   M33_SET(m33, (pa[0]-pb[0]), (p1[0]-p0[0]), (p2[0]-p0[0]),
                (pa[1]-pb[1]), (p1[1]-p0[1]), (p2[1]-p0[1]),
                (pa[2]-pb[2]), (p1[2]-p0[2]), (p2[2]-p0[2]) );

   M33_INV(m33i,m33);

   V3_SUB(pr,pa,p0);

   V3_M33_V3_MUL(tuv,m33i,pr);

   V3_SUB   (pt,pb,pa);
   V3_MUL_VS(pt,pt,tuv[0]);
   V3_ADD   (pt,pt,pa);
}

