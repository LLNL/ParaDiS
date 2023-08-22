#pragma once

#ifndef _PDS_M22_H
#define _PDS_M22_H

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "Typedefs.h"

typedef real8 M22[4];

//-----------------------------------------------------------------------------------------------
// Inline 2x2 Matrix Operations
//-----------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------

#define M22_SET1(m,s)      \
{                          \
   (m)[0]=(s); (m)[1]=(s); \
   (m)[2]=(s); (m)[3]=(s); \
}

//-----------------------------------------------------------------------------------------------

#define M22_SET(m,m00,m01,m10,m11)  \
{                                   \
   (m)[0]=(m00); (m)[1]=(m01);      \
   (m)[2]=(m10); (m)[3]=(m11);      \
}

//-----------------------------------------------------------------------------------------------

#define M22_COPY(q,p)               \
{                                   \
   (q)[0]=(p)[0]; (q)[1]=(p)[1];    \
   (q)[2]=(p)[2]; (q)[3]=(p)[3];    \
}

//-----------------------------------------------------------------------------------------------

#define M22_ADD(a,b,c)      \
{                           \
   (a)[0]=((b)[0]+(c)[0]);  \
   (a)[1]=((b)[1]+(c)[1]);  \
   (a)[2]=((b)[2]+(c)[2]);  \
   (a)[3]=((b)[3]+(c)[3]);  \
}

//-----------------------------------------------------------------------------------------------

#define M22_ADD_MS(a,b,s)   \
{                           \
   (a)[0]=((b)[0]+(s));     \
   (a)[1]=((b)[1]+(s));     \
   (a)[2]=((b)[2]+(s));     \
   (a)[3]=((b)[3]+(s));     \
}

//-----------------------------------------------------------------------------------------------

#define M22_ADD_SM(a,s,b)   \
{                           \
   (a)[0]=((s)+(b)[0]);     \
   (a)[1]=((s)+(b)[1]);     \
   (a)[2]=((s)+(b)[2]);     \
   (a)[3]=((s)+(b)[3]);     \
}

//-----------------------------------------------------------------------------------------------

#define M22_SUB(a,b,c)      \
{                           \
   (a)[0]=((b)[0]-(c)[0]);  \
   (a)[1]=((b)[1]-(c)[1]);  \
   (a)[2]=((b)[2]-(c)[2]);  \
   (a)[3]=((b)[3]-(c)[3]);  \
}

//-----------------------------------------------------------------------------------------------

#define M22_SUB_MS(a,b,s)   \
{                           \
   (a)[0]=((b)[0]-(s));     \
   (a)[1]=((b)[1]-(s));     \
   (a)[2]=((b)[2]-(s));     \
   (a)[3]=((b)[3]-(s));     \
}

//-----------------------------------------------------------------------------------------------

#define M22_SUB_SM(a,s,b)   \
{                           \
   (a)[0]=((s)-(b)[0]);     \
   (a)[1]=((s)-(b)[1]);     \
   (a)[2]=((s)-(b)[2]);     \
   (a)[3]=((s)-(b)[3]);     \
}

//-----------------------------------------------------------------------------------------------

#define M22_DOT(a,b,c)      \
{                           \
   (a)  = (b)[0]*(c)[0];    \
   (a) += (b)[1]*(c)[1];    \
   (a) += (b)[2]*(c)[2];    \
   (a) += (b)[3]*(c)[3];    \
}

//-----------------------------------------------------------------------------------------------

#define M22_IDENTITY(m)     \
{                           \
   M22_SET(m,1.0,0.0,       \
             0.0,1.0);      \
}

//-----------------------------------------------------------------------------------------------

#define M22_ZERO(m)         \
{                           \
   M22_SET(m,0.0,0.0,       \
             0.0,0.0);      \
}

//-----------------------------------------------------------------------------------------------

#define M22_ONE(m)          \
{                           \
   M22_SET(m,1.0,1.0,       \
             1.0,1.0);      \
}

//-----------------------------------------------------------------------------------------------

#define M22_RAND(m,min,max)                    \
{                                              \
   real8  s = ((real8)((max)-(min)))/RAND_MAX; \
                                               \
   (m)[0] = (min)+(s*rand());                  \
   (m)[1] = (min)+(s*rand());                  \
   (m)[2] = (min)+(s*rand());                  \
   (m)[3] = (min)+(s*rand());                  \
}

//-----------------------------------------------------------------------------------------------

#define M22_SCALE(m,x,y)  \
{                         \
   M22_SET(m,(x), 0.0,    \
             0.0, (y) );  \
}

//-----------------------------------------------------------------------------------------------

#define M22_DIAG(m,x,y)   \
{                         \
   M22_SET(m,(x), 0.0,    \
             0.0, (y) );  \
}

//-----------------------------------------------------------------------------------------------

#define M22_ROT_Z(m,rz)           \
{                                 \
   M22_SET(m,cos(rz),-sin(rz),    \
             sin(rz), cos(rz) );  \
}

//-----------------------------------------------------------------------------------------------

#define M22_FLIP_H(a,b)     \
{                           \
   M22_SET(a,b[1],b[0],     \
             b[3],b[2]);    \
}

//-----------------------------------------------------------------------------------------------

#define M22_FLIP_V(a,b)     \
{                           \
   M22_SET(a,b[2],b[3],     \
             b[0],b[1]);    \
}

//-----------------------------------------------------------------------------------------------

#define M22_RIGHT_TO_LEFT(m) \
{                            \
   M22_SET(m, 1.0, 0.0,      \
              0.0,-1.0 );    \
}

//-----------------------------------------------------------------------------------------------

#define M22_REFLECT_X(m)     \
{                            \
   M22_SET(m,-1.0, 0.0,      \
              0.0, 1.0 );    \
}

//-----------------------------------------------------------------------------------------------

#define M22_REFLECT_Y(m)     \
{                            \
   M22_SET(m, 1.0, 0.0,      \
              0.0,-1.0 );    \
}

//-----------------------------------------------------------------------------------------------

#define M22_MUL2(a,b)                           \
{                                               \
   real8 _t[4]; memcpy(_t,(a),sizeof(_t));      \
   (a)[0] = (_t[0]*(b)[0]) + (_t[1]*(b)[2]) ;   \
   (a)[1] = (_t[0]*(b)[1]) + (_t[1]*(b)[3]) ;   \
   (a)[2] = (_t[2]*(b)[0]) + (_t[3]*(b)[2]) ;   \
   (a)[3] = (_t[2]*(b)[1]) + (_t[3]*(b)[3]) ;   \
}

//-----------------------------------------------------------------------------------------------

#define M22_MUL(a,b,c)                            \
{                                                 \
   (a)[0] = ((b)[0]*(c)[0]) + ((b)[1]*(c)[2]) ;   \
   (a)[1] = ((b)[0]*(c)[1]) + ((b)[1]*(c)[3]) ;   \
   (a)[2] = ((b)[2]*(c)[0]) + ((b)[3]*(c)[2]) ;   \
   (a)[3] = ((b)[2]*(c)[1]) + ((b)[3]*(c)[3]) ;   \
}

//-----------------------------------------------------------------------------------------------

#define M22_MUL_MS(a,b,s) \
{                         \
   (a)[0]=((b)[0]*(s));   \
   (a)[1]=((b)[1]*(s));   \
   (a)[2]=((b)[2]*(s));   \
   (a)[3]=((b)[3]*(s));   \
}

//-----------------------------------------------------------------------------------------------

#define M22_MUL_SM(a,s,b) \
{                         \
   (a)[0]=((s)*(b)[0]);   \
   (a)[1]=((s)*(b)[1]);   \
   (a)[2]=((s)*(b)[2]);   \
   (a)[3]=((s)*(b)[3]);   \
}

//-----------------------------------------------------------------------------------------------

#define M22_DET(m)     \
(                      \
     ((m)[0]*(m)[3])   \
   - ((m)[2]*(m)[1])   \
)

//-----------------------------------------------------------------------------------------------

#define M22_INV(a,b)                                 \
{                                                    \
   real8 det = ((b)[0]*(b)[3]) - ((b)[2]*(b)[1]);    \
         det = ( (fabs(det)>0.0) ? 1.0/det : 0.0 );  \
                                                     \
   (a)[0] =  (b)[3] * det;                           \
   (a)[1] = -(b)[1] * det;                           \
   (a)[2] = -(b)[2] * det;                           \
   (a)[3] =  (b)[0] * det;                           \
}

#define M22_ADJOINT(a,b) \
{                        \
   (a)[0] =  (b)[3];     \
   (a)[1] = -(b)[1];     \
   (a)[2] = -(b)[2];     \
   (a)[3] =  (b)[0];     \
}

//-----------------------------------------------------------------------------------------------

#define M22_TRANSPOSE(q,p)      \
{                               \
   M22_SET(q, (p)[0], (p)[2],   \
              (p)[1], (p)[3] ); \
}

//-----------------------------------------------------------------------------------------------

#define M22_NEG(a,b)     \
{                        \
   (a)[0]= -((b)[0]);    \
   (a)[1]= -((b)[1]);    \
   (a)[2]= -((b)[2]);    \
   (a)[3]= -((b)[3]);    \
}

//-----------------------------------------------------------------------------------------------

#define M22_ABS(a,b)     \
{                        \
   (a)[0]=fabs((b)[0]);  \
   (a)[1]=fabs((b)[1]);  \
   (a)[2]=fabs((b)[2]);  \
   (a)[3]=fabs((b)[3]);  \
}

//-----------------------------------------------------------------------------------------------

#define M22_RINT(a,b)    \
{                        \
   (a)[0]=rint((b)[0]);  \
   (a)[1]=rint((b)[1]);  \
   (a)[2]=rint((b)[2]);  \
   (a)[3]=rint((b)[3]);  \
}

//-----------------------------------------------------------------------------------------------

#define M22_SIGN(a,b)    \
{                        \
   (a)[0]=( ((b)[0]<0.0) ? -1.0 : ( (((b)[0]>0.0)) ? +1.0 : 0.0) ); \
   (a)[1]=( ((b)[1]<0.0) ? -1.0 : ( (((b)[1]>0.0)) ? +1.0 : 0.0) ); \
   (a)[2]=( ((b)[2]<0.0) ? -1.0 : ( (((b)[2]>0.0)) ? +1.0 : 0.0) ); \
   (a)[3]=( ((b)[3]<0.0) ? -1.0 : ( (((b)[3]>0.0)) ? +1.0 : 0.0) ); \
}

//-----------------------------------------------------------------------------------------------

#define M22_NORMALIZE(a,b)                               \
{                                                        \
   real8 s =      fabs((b)[0]);                          \
         s = ( (s>fabs((b)[1])) ? s : fabs((b)[1]) );    \
         s = ( (s>fabs((b)[2])) ? s : fabs((b)[2]) );    \
         s = ( (s>fabs((b)[3])) ? s : fabs((b)[3]) );    \
         s = ( (s>0.0) ? (1.0/s) : 0.0 );                \
                                                         \
   (a)[0]=s*((b)[0]);                                    \
   (a)[1]=s*((b)[1]);                                    \
   (a)[2]=s*((b)[2]);                                    \
   (a)[3]=s*((b)[3]);                                    \
}

//-----------------------------------------------------------------------------------------------

#define V2_M22_V2_MUL(q,m,p)          \
{                                     \
   real8 _x = (p)[0];                 \
   real8 _y = (p)[1];                 \
                                      \
   (q)[0] = ((m)[0]*_x)+((m)[1]*_y);  \
   (q)[1] = ((m)[2]*_x)+((m)[3]*_y);  \
}

//-----------------------------------------------------------------------------------------------

#define V2_V2_M22_MUL(q,p,m)          \
{                                     \
   real8 _x = (p)[0];                 \
   real8 _y = (p)[1];                 \
                                      \
   (q)[0] = (_x*(m)[0])+(_y*(m)[2]);  \
   (q)[1] = (_x*(m)[1])+(_y*(m)[3]);  \
}

//-----------------------------------------------------------------------------------------------

#define MUL_M22_V2N(m,v,n)                \
{                                         \
   if (m && v && ((n)>0))                 \
   {                                      \
      for (int i=0,k=0; (i<n); ++i,k+=2)  \
      {                                   \
         real8 x = v[k  ];                \
         real8 y = v[k+1];                \
                                          \
         v[k  ] = (m[0]*x)+(m[1]*y);      \
         v[k+1] = (m[2]*x)+(m[3]*y);      \
      }                                   \
   }                                      \
}

//-----------------------------------------------------------------------------------------------

#define MUL_V2N_M22(m,v,n)                \
{                                         \
   if (m && v && ((n)>0))                 \
   {                                      \
      for (int i=0,k=0; (i<n); ++i,k+=2)  \
      {                                   \
         real8 x = v[k  ];                \
         real8 y = v[k+1];                \
                                          \
         v[k  ] = (x*m[0])+(y*m[2]);      \
         v[k+1] = (x*m[1])+(y*m[3]);      \
      }                                   \
   }                                      \
}

//-----------------------------------------------------------------------------------------------

#define MUL_V2N_M22_V2N(a,m,b,n)          \
{                                         \
   if (a && m && b && ((n)>0))            \
   {                                      \
      for (int i=0,k=0; (i<n); ++i,k+=2)  \
      {                                   \
         real8 x = b[k  ];                \
         real8 y = b[k+1];                \
                                          \
         a[k  ] = (m[0]*x)+(m[1]*y);      \
         a[k+1] = (m[2]*x)+(m[3]*y);      \
      }                                   \
   }                                      \
}

//-----------------------------------------------------------------------------------------------

#define MUL_V2N_V2N_M22(a,b,m,n)          \
{                                         \
   if (a && b && m && ((n)>0))            \
   {                                      \
      for (int i=0,k=0; (i<n); ++i,k+=3)  \
      {                                   \
         real8 x = b[k  ];                \
         real8 y = b[k+1];                \
                                          \
         a[k  ] = (x*m[0])+(y*m[2]);      \
         a[k+1] = (x*m[1])+(y*m[3]);      \
      }                                   \
   }                                      \
}

//-----------------------------------------------------------------------------------------------

#define M22_FOR_ALL(fn) { for (int i=0; (i<4); ++i) { fn; } }

//-----------------------------------------------------------------------------------------------

#define  M22_SIN(a,b)               { M22_FOR_ALL( (a)[i]=sin((b)[i]) ) }
#define  M22_COS(a,b)               { M22_FOR_ALL( (a)[i]=cos((b)[i]) ) }
#define  M22_TAN(a,b)               { M22_FOR_ALL( (a)[i]=tan((b)[i]) ) }
#define  M22_LOG(a,b)               { M22_FOR_ALL( (a)[i]=(((b)[i]>0.0) ? log  ((b)[i]) : 0.0) ) }
#define  M22_LOG2(a,b)              { M22_FOR_ALL( (a)[i]=(((b)[i]>0.0) ? log2 ((b)[i]) : 0.0) ) }
#define  M22_LOG10(a,b)             { M22_FOR_ALL( (a)[i]=(((b)[i]>0.0) ? log10((b)[i]) : 0.0) ) }
#define  M22_SQR(a,b)               { M22_FOR_ALL( (a)[i]=(b)[i]*(b)[i] ) }
#define  M22_SQRT(a,b)              { M22_FOR_ALL( (a)[i]=(((b)[i]>0.0) ? sqrt ((b)[i]) : 0.0) ) }

#define  M22_THRESHOLD(a,b,s)       { M22_FOR_ALL( (a)[i]=( (fabs((b)[i])<(s)) ? 0.0    : (b)[i] ) ) }
#define  M22_CAP(a,b,max)           { M22_FOR_ALL( (a)[i]=( ((b)[i]<(max))     ? (b)[i] : (max)  ) ) }
#define  M22_CONSTRAIN(a,b,min,max) { M22_FOR_ALL( (a)[i]=( ((b)[i]<(min))     ? (min)  : ( ((b)[i]<(max)) ? (b)[i] : (max) ) ) ) }

//-----------------------------------------------------------------------------------------------

class Matrix_2x2
{
   public :
      real8  mtx[4];

   public :
       Matrix_2x2(void)                     { M22_ZERO(mtx);       }
       Matrix_2x2(const Matrix_2x2     & m) { M22_COPY(mtx,m.mtx); }

       Matrix_2x2(const unsigned char    s) {          M22_SET1(mtx,s);   }
       Matrix_2x2(const unsigned short   s) {          M22_SET1(mtx,s);   }
       Matrix_2x2(const          short   s) {          M22_SET1(mtx,s);   }
       Matrix_2x2(const unsigned int     s) {          M22_SET1(mtx,s);   }
       Matrix_2x2(const          int     s) {          M22_SET1(mtx,s);   }
       Matrix_2x2(const          real8   s) {          M22_SET1(mtx,s);   }

       Matrix_2x2(const unsigned char   *s) { if (s) { M22_COPY(mtx,s); } }
       Matrix_2x2(const unsigned short  *s) { if (s) { M22_COPY(mtx,s); } }
       Matrix_2x2(const          short  *s) { if (s) { M22_COPY(mtx,s); } }
       Matrix_2x2(const unsigned int    *s) { if (s) { M22_COPY(mtx,s); } }
       Matrix_2x2(const          int    *s) { if (s) { M22_COPY(mtx,s); } }
       Matrix_2x2(const          real8  *s) { if (s) { M22_COPY(mtx,s); } }

       Matrix_2x2(const real8 m00, const real8 m01,
                  const real8 m10, const real8 m11  ) { M22_SET(mtx,m00,m01,m10,m11); }

      ~Matrix_2x2() {}

      const Matrix_2x2 & operator =  (const Matrix_2x2     & m) { M22_COPY(mtx,m.mtx); return(*this); }

      const Matrix_2x2 & operator =  (const unsigned char    s) { M22_SET1(mtx,s); return(*this); }
      const Matrix_2x2 & operator =  (const unsigned short   s) { M22_SET1(mtx,s); return(*this); }
      const Matrix_2x2 & operator =  (const          short   s) { M22_SET1(mtx,s); return(*this); }
      const Matrix_2x2 & operator =  (const unsigned int     s) { M22_SET1(mtx,s); return(*this); }
      const Matrix_2x2 & operator =  (const          int     s) { M22_SET1(mtx,s); return(*this); }
      const Matrix_2x2 & operator =  (const          real8   s) { M22_SET1(mtx,s); return(*this); }

      const Matrix_2x2 & operator =  (const unsigned char   *s) { if (s) { M22_COPY(mtx,s); } return(*this); }
      const Matrix_2x2 & operator =  (const unsigned short  *s) { if (s) { M22_COPY(mtx,s); } return(*this); }
      const Matrix_2x2 & operator =  (const          short  *s) { if (s) { M22_COPY(mtx,s); } return(*this); }
      const Matrix_2x2 & operator =  (const unsigned int    *s) { if (s) { M22_COPY(mtx,s); } return(*this); }
      const Matrix_2x2 & operator =  (const          int    *s) { if (s) { M22_COPY(mtx,s); } return(*this); }
      const Matrix_2x2 & operator =  (const          real8  *s) { if (s) { M22_COPY(mtx,s); } return(*this); }

            void         operator += (const Matrix_2x2     & m) { M22_ADD   (mtx,mtx,m.mtx); }
            void         operator += (const real8            s) { M22_ADD_MS(mtx,mtx,s    ); }
            void         operator -= (const Matrix_2x2     & m) { M22_SUB   (mtx,mtx,m.mtx); }
            void         operator -= (const real8            s) { M22_SUB_MS(mtx,mtx,s    ); }
            void         operator *= (const Matrix_2x2     & m) { Matrix_2x2 tmp(*this); M22_MUL(mtx,tmp.mtx,m.mtx); }
            void         operator *= (const real8            s) { M22_MUL_MS(mtx,mtx,s    ); }
            void         operator /= (const real8            s) { real8 x = ( (fabs(s)>0.0) ? 1.0/s : 0.0 ); M22_MUL_MS(mtx,mtx,x); }

            Matrix_2x2   operator +  (const Matrix_2x2     & m) { Matrix_2x2 tmp(*this); tmp+=m; return(tmp); }
            Matrix_2x2   operator +  (const real8            s) { Matrix_2x2 tmp(*this); tmp+=s; return(tmp); }
            Matrix_2x2   operator -  (const Matrix_2x2     & m) { Matrix_2x2 tmp(*this); tmp-=m; return(tmp); }
            Matrix_2x2   operator -  (const real8            s) { Matrix_2x2 tmp(*this); tmp-=s; return(tmp); }
            Matrix_2x2   operator *  (const Matrix_2x2     & m) { Matrix_2x2 tmp, tmp2(*this); M22_MUL(tmp.mtx,tmp2.mtx,m.mtx); return(tmp); }
            Matrix_2x2   operator *  (const real8            s) { Matrix_2x2 tmp(*this); tmp*=s; return(tmp); }

            Matrix_2x2   operator - ()                          { Matrix_2x2 tmp; M22_NEG(tmp.mtx,mtx); return(tmp); }

                         operator  real8  *() const { return ( (real8  *) mtx ); }

      real8 Dot           (const Matrix_2x2 & m) const { real8 s=0.0; M22_DOT(s,mtx,m.mtx); return(s); }

      void  Transpose     (void)                       { Matrix_2x2 tmp(*this); M22_TRANSPOSE(mtx,tmp.mtx); }
      void  Transpose     (const Matrix_2x2 & m )      { Matrix_2x2 tmp(m);     M22_TRANSPOSE(mtx,tmp.mtx); }

      void  Normalize     (void)                       { M22_NORMALIZE(mtx,  mtx); }
      void  Normalize     (const Matrix_2x2 & m )      { M22_NORMALIZE(mtx,m.mtx); }

      real8 Determinant   (void)                       { return(M22_DET(mtx)); }

      void  Inverse       (void)                       { Matrix_2x2 tmp(*this); M22_INV(mtx,tmp.mtx); }
      void  Inverse       (const Matrix_2x2 & m )      { Matrix_2x2 tmp(m);     M22_INV(mtx,tmp.mtx); }

      int   Compare       (const Matrix_2x2 & m, const real8 eps=4.0e-8) const;
      int   Compare       (const real8        s, const real8 eps=4.0e-8) const;
      int   Compare       (const real8       *s, const real8 eps=4.0e-8) const;

      int   Near          (const Matrix_2x2 & m, const real8 eps=4.0e-8) const { return( Compare(m,eps) ? 1 : 0 ); }
      int   Near          (const real8        s, const real8 eps=4.0e-8) const { return( Compare(s,eps) ? 1 : 0 ); }
      int   Near          (const real8       *s, const real8 eps=4.0e-8) const { return( Compare(s,eps) ? 1 : 0 ); }

      int   operator ==   (const Matrix_2x2 & m) const { return( Compare(m) ? 1 : 0 ); }
      int   operator ==   (const real8        s) const { return( Compare(s) ? 1 : 0 ); }
      int   operator ==   (const real8       *s) const { return( Compare(s) ? 1 : 0 ); }

      int   operator !=   (const Matrix_2x2 & m) const { return( Compare(m) ? 0 : 1 ); }
      int   operator !=   (const real8        s) const { return( Compare(s) ? 0 : 1 ); }
      int   operator !=   (const real8       *s) const { return( Compare(s) ? 0 : 1 ); }

      int   operator <    (const Matrix_2x2 & m) const { M22_FOR_ALL( if(mtx[i]>=m.mtx[i]) return(0) )  return(1); }
      int   operator <    (const real8        s) const { M22_FOR_ALL( if(mtx[i]>=s       ) return(0) )  return(1); }
      int   operator <    (const real8       *s) const { if (!s) return(0);
                                                         M22_FOR_ALL( if(mtx[i]>=s[i]    ) return(0) ) return(1); }

      int   operator <=   (const Matrix_2x2 & m) const { M22_FOR_ALL( if(mtx[i]> m.mtx[i]) return(0) ) return(1); }
      int   operator <=   (const real8        s) const { M22_FOR_ALL( if(mtx[i]> s       ) return(0) ) return(1); }
      int   operator <=   (const real8       *s) const { if (!s) return(0);
                                                         M22_FOR_ALL( if(mtx[i]> s[i]    ) return(0) ) return(1); }

      int   operator >    (const Matrix_2x2 & m) const { M22_FOR_ALL( if(mtx[i]<=m.mtx[i]) return(0) ) return(1); }
      int   operator >    (const real8        s) const { M22_FOR_ALL( if(mtx[i]<=s       ) return(0) ) return(1); }
      int   operator >    (const real8       *s) const { if (!s) return(0);
                                                         M22_FOR_ALL( if(mtx[i]<=s[i]    ) return(0) ) return(1); }

      int   operator >=   (const Matrix_2x2 & m) const { M22_FOR_ALL( if(mtx[i]< m.mtx[i]) return(0) ) return(1); }
      int   operator >=   (const real8        s) const { M22_FOR_ALL( if(mtx[i]< s       ) return(0) ) return(1); }
      int   operator >=   (const real8       *s) const { if (!s) return(0);
                                                         M22_FOR_ALL( if(mtx[i]< s[i]    ) return(0) ) return(1); }

      void  Identity      (void) { M22_IDENTITY(mtx); }
      void  Zero          (void) { M22_ZERO    (mtx); }
      void  One           (void) { M22_ONE     (mtx); }
      void  Rand          (const real8 min=0.0, const real8 max=1.0) { M22_RAND(mtx,min,max); }

      void  Abs           (void)                 { M22_ABS  (mtx,  mtx); }
      void  Abs           (const Matrix_2x2 & m) { M22_ABS  (mtx,m.mtx); }
      void  Rint          (void)                 { M22_RINT (mtx,  mtx); }
      void  Rint          (const Matrix_2x2 & m) { M22_RINT (mtx,m.mtx); }
      void  Sign          (void)                 { M22_SIGN (mtx,  mtx); }
      void  Sign          (const Matrix_2x2 & m) { M22_SIGN (mtx,m.mtx); }

      void  Sin           (void)                 { M22_SIN  (mtx,  mtx); }
      void  Sin           (const Matrix_2x2 & m) { M22_SIN  (mtx,m.mtx); }
      void  Cos           (void)                 { M22_COS  (mtx,  mtx); }
      void  Cos           (const Matrix_2x2 & m) { M22_COS  (mtx,m.mtx); }
      void  Tan           (void)                 { M22_TAN  (mtx,  mtx); }
      void  Tan           (const Matrix_2x2 & m) { M22_TAN  (mtx,m.mtx); }
      void  Log           (void)                 { M22_LOG  (mtx,  mtx); }
      void  Log           (const Matrix_2x2 & m) { M22_LOG  (mtx,m.mtx); }
      void  Log2          (void)                 { M22_LOG2 (mtx,  mtx); }
      void  Log2          (const Matrix_2x2 & m) { M22_LOG2 (mtx,m.mtx); }
      void  Log10         (void)                 { M22_LOG10(mtx,  mtx); }
      void  Log10         (const Matrix_2x2 & m) { M22_LOG10(mtx,m.mtx); }
      void  Sqr           (void)                 { M22_SQR  (mtx,  mtx); }
      void  Sqr           (const Matrix_2x2 & m) { M22_SQR  (mtx,m.mtx); }
      void  Sqrt          (void)                 { M22_SQRT (mtx,  mtx); }
      void  Sqrt          (const Matrix_2x2 & m) { M22_SQRT (mtx,m.mtx); }

      void  Threshold     (const real8 eps)                                        { M22_THRESHOLD(mtx,  mtx,eps); }
      void  Threshold     (const Matrix_2x2 & m, const real8 eps)                  { M22_THRESHOLD(mtx,m.mtx,eps); }

      void  Cap           (const real8 max)                                        { M22_CAP      (mtx,  mtx,max); }
      void  Cap           (const Matrix_2x2 & m, const real8 max)                  { M22_CAP      (mtx,m.mtx,max); }

      void  Constrain     (const real8 min, const real8 max)                       { M22_CONSTRAIN(mtx,  mtx,min,max); }
      void  Constrain     (const Matrix_2x2 & m, const real8 min, const real8 max) { M22_CONSTRAIN(mtx,m.mtx,min,max); }

      void  Scale         (const real8 sx, const real8 sy)                         { M22_SCALE   (mtx,sx,sy); }
      void  Diag          (const real8 sx, const real8 sy)                         { M22_DIAG    (mtx,sx,sy); }

      void  Rotate_Z      (const real8 rz) { M22_ROT_Z(mtx,rz); }

      void  Flip_H        (void)                  { real8 tmp[9];  M22_COPY(tmp,  mtx); M22_FLIP_H(mtx,tmp); }
      void  Flip_H        (const Matrix_2x2 & m)  { real8 tmp[9];  M22_COPY(tmp,m.mtx); M22_FLIP_H(mtx,tmp); }
      void  Flip_V        (void)                  { real8 tmp[9];  M22_COPY(tmp,  mtx); M22_FLIP_V(mtx,tmp); }
      void  Flip_V        (const Matrix_2x2 & m)  { real8 tmp[9];  M22_COPY(tmp,m.mtx); M22_FLIP_V(mtx,tmp); }

      void  Right_To_Left (void)           { M22_RIGHT_TO_LEFT(mtx); }
      void  Reflect_X     (void)           { M22_REFLECT_X    (mtx); }
      void  Reflect_Y     (void)           { M22_REFLECT_Y    (mtx); }

      void  Transform_M22_V2N (real8 *v ,                  const int n) const { MUL_M22_V2N    (mtx,v,n);     }
      void  Transform_V2N_M22 (real8 *v ,                  const int n) const { MUL_V2N_M22    (mtx,v,n);     }
      void  Transform_M22_V2N (real8 *va, const real8 *vb, const int n) const { MUL_V2N_M22_V2N(va,mtx,vb,n); }
      void  Transform_V2N_M22 (real8 *va, const real8 *vb, const int n) const { MUL_V2N_V2N_M22(va,vb,mtx,n); }

      void  Print         (FILE *fd=stdout, const char *fmt="%8.4lf ") const;
};

//-----------------------------------------------------------------------------------------------

inline Matrix_2x2 M22_Identity      (void) { Matrix_2x2 m; M22_IDENTITY     (m.mtx); return(m); }
inline Matrix_2x2 M22_Zero          (void) { Matrix_2x2 m; M22_ZERO         (m.mtx); return(m); }
inline Matrix_2x2 M22_One           (void) { Matrix_2x2 m; M22_ONE          (m.mtx); return(m); }

inline Matrix_2x2 M22_Right_To_Left (void) { Matrix_2x2 m; M22_RIGHT_TO_LEFT(m.mtx); return(m); }
inline Matrix_2x2 M22_Reflect_X     (void) { Matrix_2x2 m; M22_REFLECT_X    (m.mtx); return(m); }
inline Matrix_2x2 M22_Reflect_Y     (void) { Matrix_2x2 m; M22_REFLECT_Y    (m.mtx); return(m); }

extern int  M22_Near  (M22 ma, const M22   mb);
extern int  M22_Near  (M22 ma, const M22   mb, const real8 eps);
extern int  M22_Near  (M22 m , const real8 b , const real8 eps);

extern int  M22_Near  (M22 m        , const real8 m00, const real8 m01,
                                      const real8 m10, const real8 m11, const real8 eps);

extern int  M22_Near  (real8 m[2][2], const real8 m00, const real8 m01,
                                      const real8 m10, const real8 m11, const real8 eps);

extern void M22_Print (const M22 m, const char *title);

#endif
