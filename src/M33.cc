//-----------------------------------------------------------------------------------------------
// Module:  M33 - implements a simple 3x3 transform matrix
//-----------------------------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "V3.h"
#include "M33.h"

int Matrix_3x3::Compare (const Matrix_3x3 & m, const real8 eps) const
{
   for (int i=0; (i<9); i+=3)
   {  
      if ( fabs(mtx[i  ]-m.mtx[i  ])>eps ) return(0);
      if ( fabs(mtx[i+1]-m.mtx[i+1])>eps ) return(0);
      if ( fabs(mtx[i+2]-m.mtx[i+2])>eps ) return(0);
   }

   return(1);
}

int Matrix_3x3::Compare (const real8        s, const real8 eps) const
{
   for (int i=0; (i<9); i+=3)
   {  
      if ( fabs(mtx[i  ]-s)>eps ) return(0);
      if ( fabs(mtx[i+1]-s)>eps ) return(0);
      if ( fabs(mtx[i+2]-s)>eps ) return(0);
   }

   return(1);
}

int Matrix_3x3::Compare (const real8       *s, const real8 eps) const
{
   if (!s) return(0);

   for (int i=0; (i<9); i+=3)
   {  
      if ( fabs(mtx[i  ]-s[i  ])>eps ) return(0);
      if ( fabs(mtx[i+1]-s[i+1])>eps ) return(0);
      if ( fabs(mtx[i+2]-s[i+2])>eps ) return(0);
   }

   return(1);
}

//-----------------------------------------------------------------------------------------------

void Matrix_3x3::Print (FILE *fd, const char *xfmt ) const
{
   const char *fmt = ( xfmt ? xfmt : "%8.4lf " );

   for (int i=0; (i<9); i+=3)
   {
      fprintf(fd,fmt,mtx[i  ]);
      fprintf(fd,fmt,mtx[i+1]);
      fprintf(fd,fmt,mtx[i+2]); fprintf(fd,"\n");
   }
}

//-----------------------------------------------------------------------------------------------

int M33_Near (real8 *ma, const real8 *mb)
{
   return(M33_Near(ma,mb,1e-12));
}

//-----------------------------------------------------------------------------------------------

int M33_Near (real8 *ma, const real8 *mb, const real8 eps)
{
   if (ma && mb)
   {
      if (fabs(ma[0]-mb[0])>eps) return(0);
      if (fabs(ma[1]-mb[1])>eps) return(0);
      if (fabs(ma[2]-mb[2])>eps) return(0);
      if (fabs(ma[3]-mb[3])>eps) return(0);
      if (fabs(ma[4]-mb[4])>eps) return(0);
      if (fabs(ma[5]-mb[5])>eps) return(0);
      if (fabs(ma[6]-mb[6])>eps) return(0);
      if (fabs(ma[7]-mb[7])>eps) return(0);
      if (fabs(ma[8]-mb[8])>eps) return(0);

      return(1);
   }

   return(0);
}

//-----------------------------------------------------------------------------------------------

int M33_Near (real8 *m, const real8 b, const real8 eps )
{
   if (m)
   {
      if (fabs(m[0]-b)>eps) return(0);
      if (fabs(m[1]-b)>eps) return(0);
      if (fabs(m[2]-b)>eps) return(0);
      if (fabs(m[3]-b)>eps) return(0);
      if (fabs(m[4]-b)>eps) return(0);
      if (fabs(m[5]-b)>eps) return(0);
      if (fabs(m[6]-b)>eps) return(0);
      if (fabs(m[7]-b)>eps) return(0);
      if (fabs(m[8]-b)>eps) return(0);

      return(1);
   }

   return(0);
}

//-----------------------------------------------------------------------------------------------

int M33_Near (
   real8 *m, 
   const real8 m00, const real8 m01, const real8 m02,
   const real8 m10, const real8 m11, const real8 m12,
   const real8 m20, const real8 m21, const real8 m22, const real8 eps )
{
   if (m)
   {
      if (fabs(m[0]-m00)>eps) return(0);
      if (fabs(m[1]-m01)>eps) return(0);
      if (fabs(m[2]-m02)>eps) return(0);
      if (fabs(m[3]-m10)>eps) return(0);
      if (fabs(m[4]-m11)>eps) return(0);
      if (fabs(m[5]-m12)>eps) return(0);
      if (fabs(m[6]-m20)>eps) return(0);
      if (fabs(m[7]-m21)>eps) return(0);
      if (fabs(m[8]-m22)>eps) return(0);

      return(1);
   }

   return(0);
}

//-----------------------------------------------------------------------------------------------

int M33_Near (
   real8 m[3][3], 
   const real8 m00, const real8 m01, const real8 m02,
   const real8 m10, const real8 m11, const real8 m12,
   const real8 m20, const real8 m21, const real8 m22, const real8 eps)
{
   if (m)
   {
      if (fabs(m[0][0]-m00)>eps) return(0);
      if (fabs(m[0][1]-m01)>eps) return(0);
      if (fabs(m[0][2]-m02)>eps) return(0);
      if (fabs(m[1][0]-m10)>eps) return(0);
      if (fabs(m[1][1]-m11)>eps) return(0);
      if (fabs(m[1][2]-m12)>eps) return(0);
      if (fabs(m[2][0]-m20)>eps) return(0);
      if (fabs(m[2][1]-m21)>eps) return(0);
      if (fabs(m[2][2]-m22)>eps) return(0);

      return(1);
   }

   return(0);
}

//-----------------------------------------------------------------------------------------------

void M33_Print (const real8 *m, const char *title)
{
   if (m)
   {
      if (title) printf("<m3x3 title=\"%s\">\n",title);
      else       printf("<m3x3>\n");

      printf("%7.4lf,%7.4lf,%7.4lf,"  , m[0], m[1], m[2] );
      printf("%7.4lf,%7.4lf,%7.4lf,"  , m[3], m[4], m[5] );
      printf("%7.4lf,%7.4lf,%7.4lf \n", m[6], m[7], m[8] );

      printf("</m3x3>\n");
   }
}

