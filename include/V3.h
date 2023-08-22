#pragma once

#ifndef _PDS_V3_H
#define _PDS_V3_H

/*--------------------------------------------------------------------------
 *	V3.h Declarations for 1x3 vector class
 *------------------------------------------------------------------------*/

#include <math.h>
#include "Typedefs.h"

typedef real8 V3[3];

#define V3_SET(a,x,y,z)         { (a)[0]=(x);  (a)[1]=(y);  (a)[2]=(z);  }
#define V3_ZERO(a)              { (a)[0]=0.0;  (a)[1]=0.0;  (a)[2]=0.0;  }
#define V3_ONE(a)               { (a)[0]=1.0;  (a)[1]=1.0;  (a)[2]=1.0;  }
#define V3_PI(a)                { (a)[0]=M_PI; (a)[1]=M_PI; (a)[2]=M_PI; }

#define V3_FROM_SCALAR(v,x,y,z) { (v)[0]=(x); \
                                  (v)[1]=(y); \
                                  (v)[2]=(z); }

#define V3_TO_SCALAR(x,y,z,v)   { (x)=(v)[0]; \
                                  (y)=(v)[1]; \
                                  (z)=(v)[2]; }

#define V3_COPY(a,b)            { (a)[0]=(b)[0]; \
                                  (a)[1]=(b)[1]; \
                                  (a)[2]=(b)[2]; }

#define V3_ADD(a,b,c)           { (a)[0]=((b)[0]+(c)[0]); \
                                  (a)[1]=((b)[1]+(c)[1]); \
                                  (a)[2]=((b)[2]+(c)[2]); }

#define V3_ADD_VS(a,b,s)        { (a)[0]=((b)[0]+(s)); \
                                  (a)[1]=((b)[1]+(s)); \
                                  (a)[2]=((b)[2]+(s)); }

#define V3_ADD_SV(a,s,b)        { (a)[0]=((s)+(b)[0]); \
                                  (a)[1]=((s)+(b)[1]); \
                                  (a)[2]=((s)+(b)[2]); }

#define V3_SUB(a,b,c)           { (a)[0]=((b)[0]-(c)[0]); \
                                  (a)[1]=((b)[1]-(c)[1]); \
                                  (a)[2]=((b)[2]-(c)[2]); }

#define V3_SUB_VS(a,b,s)        { (a)[0]=((b)[0]-(s)); \
                                  (a)[1]=((b)[1]-(s)); \
                                  (a)[2]=((b)[2]-(s)); }

#define V3_SUB_SV(a,s,b)        { (a)[0]=((s)-(b)[0]); \
                                  (a)[1]=((s)-(b)[1]); \
                                  (a)[2]=((s)-(b)[2]); }

#define V3_MUL(a,b,c)           { (a)[0]=((b)[0]*(c)[0]); \
                                  (a)[1]=((b)[1]*(c)[1]); \
                                  (a)[2]=((b)[2]*(c)[2]); }

#define V3_MUL_VS(a,b,s)        { (a)[0]=((b)[0]*(s)); \
                                  (a)[1]=((b)[1]*(s)); \
                                  (a)[2]=((b)[2]*(s)); }

#define V3_MUL_SV(a,s,b)        { (a)[0]=((s)*(b)[0]); \
                                  (a)[1]=((s)*(b)[1]); \
                                  (a)[2]=((s)*(b)[2]); }

#define V3_SCALE(a,b,s)         { (a)[0]=((b)[0]*(s)); \
                                  (a)[1]=((b)[1]*(s)); \
                                  (a)[2]=((b)[2]*(s)); }

#define V3_SCALE_S3(a,b,x,y,z)  { (a)[0]=((b)[0]*(x)); \
                                  (a)[1]=((b)[1]*(y)); \
                                  (a)[2]=((b)[2]*(z)); }

#define V3_ACCUM(a,b)           { (a)[0]+=((b)[0]); \
                                  (a)[1]+=((b)[1]); \
                                  (a)[2]+=((b)[2]); }

#define V3_REMOVE(a,b)          { (a)[0]-=((b)[0]); \
                                  (a)[1]-=((b)[1]); \
                                  (a)[2]-=((b)[2]); }

#define V3_HSUM(a)              (  (a)[0]+(a)[1]+(a)[2] )

#define V3_DOT(a,b)             ( ((a)[0]*(b)[0]) + ((a)[1]*(b)[1]) + ((a)[2]*(b)[2]) )

#define V3_NORMAL(a)            ( sqrt( ((a)[0]*(a)[0]) + ((a)[1]*(a)[1]) + ((a)[2]*(a)[2]) ) )
#define V3_MAG(a)               ( sqrt( ((a)[0]*(a)[0]) + ((a)[1]*(a)[1]) + ((a)[2]*(a)[2]) ) )

#define V3_DIST(a,b)            ( sqrt(   (((b)[0]-(a)[0])*((b)[0]-(a)[0])) \
                                        + (((b)[1]-(a)[1])*((b)[1]-(a)[1])) \
                                        + (((b)[2]-(a)[2])*((b)[2]-(a)[2])) ) )

#define V3_DIST2(a,b)           (         (((b)[0]-(a)[0])*((b)[0]-(a)[0])) \
                                        + (((b)[1]-(a)[1])*((b)[1]-(a)[1])) \
                                        + (((b)[2]-(a)[2])*((b)[2]-(a)[2]))   )

#define V3_CROSS(a,b,c)         { (a)[0]=(((b)[1]*(c)[2])-((b)[2]*(c)[1])); \
                                  (a)[1]=(((b)[2]*(c)[0])-((b)[0]*(c)[2])); \
                                  (a)[2]=(((b)[0]*(c)[1])-((b)[1]*(c)[0])); }

#define V3_ALIGN_SURFACE(a,b,n) { real8 _s = (   ((b)[0]*(n)[0])    \
                                               + ((b)[1]*(n)[1])    \
                                               + ((b)[2]*(n)[2]) ); \
                                  (a)[0] = (b)[0] - _s*(n)[0];      \
                                  (a)[1] = (b)[1] - _s*(n)[1];      \
                                  (a)[2] = (b)[2] - _s*(n)[2];      }

#define V3_NEG(a,b)             { (a)[0]=(-(b)[0]); \
                                  (a)[1]=(-(b)[1]); \
                                  (a)[2]=(-(b)[2]); }

#define V3_ABS(a,b)             { (a)[0]=fabs((b)[0]); \
                                  (a)[1]=fabs((b)[1]); \
                                  (a)[2]=fabs((b)[2]); }

#define V3_ABS_DIFF(a,b,c)      { (a)[0]=fabs((b)[0]-(c)[0]); \
                                  (a)[1]=fabs((b)[1]-(c)[1]); \
                                  (a)[2]=fabs((b)[2]-(c)[2]); }

#define V3_MIDPOINT(a,b,c)      { (a)[0]=((b)[0]+(c)[0])/2.0; \
                                  (a)[1]=((b)[1]+(c)[1])/2.0; \
                                  (a)[2]=((b)[2]+(c)[2])/2.0; }

#define V3_NORMALIZE(a,b)       { real8 _s = sqrt( ((b)[0]*(b)[0]) + ((b)[1]*(b)[1]) + ((b)[2]*(b)[2]) ); \
                                        _s = ( (_s>0.0) ? 1.0/_s : 0.0 ); \
                                  (a)[0]=_s*(b)[0]; \
                                  (a)[1]=_s*(b)[1]; \
                                  (a)[2]=_s*(b)[2]; }

#define V3_UNIT_VECTOR(a,b,c)   { (a)[0]=((c)[0]-(b)[0]); \
                                  (a)[1]=((c)[1]-(b)[1]); \
                                  (a)[2]=((c)[2]-(b)[2]); \
                                  real8 _s = sqrt( ((a)[0]*(a)[0]) + ((a)[1]*(a)[1]) + ((a)[2]*(a)[2]) ); \
                                        _s = ( (_s>0.0) ? 1.0/_s : 0.0 ); \
                                  (a)[0]=_s*(a)[0]; \
                                  (a)[1]=_s*(a)[1]; \
                                  (a)[2]=_s*(a)[2]; }

#define V3_RINT(a,b,s,l)        { (a)[0]=(rint((b)[0]*(s)[0])*l[0]); \
                                  (a)[1]=(rint((b)[1]*(s)[1])*l[1]); \
                                  (a)[2]=(rint((b)[2]*(s)[2])*l[2]); }

#define V3_SWAP(a,b)            { real8 _s = (a)[0]; (a)[0]=(b)[0]; (b)[0]=_s; \
                                        _s = (a)[1]; (a)[1]=(b)[1]; (b)[1]=_s; \
                                        _s = (a)[2]; (a)[2]=(b)[2]; (b)[2]=_s; }

#define V3_LOOP(a)              for (int i=0; (i<3); ++i) { a; }

#define V3_MIN(x,y,z)           ( (x)<(y) ? ( (x)<(z) ? (x) : (z)) : ( (y)<(z) ? (y) : (z)) )
#define V3_MAX(x,y,z)           ( (x)>(y) ? ( (x)>(z) ? (x) : (z)) : ( (y)>(z) ? (y) : (z)) )

#define V3_MIN_INDEX(x,y,z)     ( ((x)<(y)) ? ( ((x)<(z)) ? 0 : 2 ) \
                                            : ( ((y)<(z)) ? 1 : 2 ) );

#define V3_MAX_INDEX(x,y,z)     ( ((x)>(y)) ? ( ((x)>(z)) ? 0 : 2 ) \
                                            : ( ((y)>(z)) ? 1 : 2 ) );

#define V3_NEAR(a,b,tol)        (     fabs((a)[0]-(b)[0])>(tol)            \
                                   || fabs((a)[1]-(b)[1])>(tol)            \
                                   || fabs((a)[2]-(b)[2])>(tol)  ? 0 : 1 )

//-----------------------------------------------------------------------------------------------

extern int   V3_Near                (const V3     va, const V3    vb);
extern int   V3_Near                (const V3     va, const V3    vb, const real8 eps);
extern int   V3_Near                (const V3     v , const real8 x , const real8 y, const real8 z, const real8 eps);
extern int   V3_Near                (const real8 *v , const real8 x , const real8 y, const real8 z, const real8 eps);

extern void  V3_Normalize           (      real8 *v , const int vn);
extern void  V3_Scale               (      real8 *v , const int vn, const real8 s);
extern void  V3_Scale               (      real8 *v , const int vn, const real8 sx, const real8 sy, const real8 sz );

extern int   V3_Aligned             (const real8 *v , const real8  x , const real8 y, const real8 z, const real8 eps=1.0e-8);
extern int   V3_Aligned             (const real8 *v0, const real8 *v1,                               const real8 eps=1.0e-8);

extern void  V3_Ray_Plane_Intersect (const V3 pa, const V3 pb, const V3 p0, const V3 p1, const V3 p2, V3 pt, V3 tuv);

//-----------------------------------------------------------------------------------------------

#endif   // (_PDS_V3_H)
