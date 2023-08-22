//-----------------------------------------------------------------------------------------------
// Module:  M22 - implements a simple 3x3 transform matrix
//-----------------------------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "M22.h"

int Matrix_2x2::Compare (const Matrix_2x2 & m, const real8 eps) const
{
   if ( fabs(mtx[0]-m.mtx[0])>eps ) return(0);
   if ( fabs(mtx[1]-m.mtx[1])>eps ) return(0);
   if ( fabs(mtx[2]-m.mtx[2])>eps ) return(0);
   if ( fabs(mtx[3]-m.mtx[3])>eps ) return(0);

   return(1);
}

int Matrix_2x2::Compare (const real8  s, const real8 eps) const
{
   if ( fabs(mtx[0]-s)>eps ) return(0);
   if ( fabs(mtx[1]-s)>eps ) return(0);
   if ( fabs(mtx[2]-s)>eps ) return(0);
   if ( fabs(mtx[3]-s)>eps ) return(0);

   return(1);
}

int Matrix_2x2::Compare (const real8  *s, const real8 eps) const
{
   if (!s) return(0);

   if ( fabs(mtx[0]-s[0])>eps ) return(0);
   if ( fabs(mtx[1]-s[1])>eps ) return(0);
   if ( fabs(mtx[2]-s[2])>eps ) return(0);
   if ( fabs(mtx[3]-s[3])>eps ) return(0);

   return(1);
}

//-----------------------------------------------------------------------------------------------

void Matrix_2x2::Print (FILE *fd, const char *xfmt ) const
{
   const char *fmt = ( xfmt ? xfmt : "%8.4lf " );

   fprintf(fd,fmt,mtx[0]);
   fprintf(fd,fmt,mtx[1]); fprintf(fd,"\n");
   fprintf(fd,fmt,mtx[2]);
   fprintf(fd,fmt,mtx[3]); fprintf(fd,"\n");
}

//-----------------------------------------------------------------------------------------------

int M22_Near (real8 *ma, const real8 *mb)
{
   return(M22_Near(ma,mb,1e-12));
}

//-----------------------------------------------------------------------------------------------

int M22_Near (real8 *ma, const real8 *mb, const real8 eps)
{
   if (ma && mb)
   {
      if (fabs(ma[0]-mb[0])>eps) return(0);
      if (fabs(ma[1]-mb[1])>eps) return(0);
      if (fabs(ma[2]-mb[2])>eps) return(0);
      if (fabs(ma[3]-mb[3])>eps) return(0);

      return(1);
   }

   return(0);
}

//-----------------------------------------------------------------------------------------------

int M22_Near (real8 *m, const real8 b, const real8 eps )
{
   if (m)
   {
      if (fabs(m[0]-b)>eps) return(0);
      if (fabs(m[1]-b)>eps) return(0);
      if (fabs(m[2]-b)>eps) return(0);
      if (fabs(m[3]-b)>eps) return(0);

      return(1);
   }

   return(0);
}

//-----------------------------------------------------------------------------------------------

int M22_Near (
   real8 *m, 
   const real8 m00, const real8 m01, 
   const real8 m10, const real8 m11, const real8 eps )
{
   if (m)
   {
      if (fabs(m[0]-m00)>eps) return(0);
      if (fabs(m[1]-m01)>eps) return(0);
      if (fabs(m[2]-m10)>eps) return(0);
      if (fabs(m[3]-m11)>eps) return(0);

      return(1);
   }

   return(0);
}

//-----------------------------------------------------------------------------------------------

int M22_Near (
   real8 m[3][3], 
   const real8 m00, const real8 m01, 
   const real8 m10, const real8 m11, const real8 eps)
{
   if (m)
   {
      if (fabs(m[0][0]-m00)>eps) return(0);
      if (fabs(m[0][1]-m01)>eps) return(0);
      if (fabs(m[1][0]-m10)>eps) return(0);
      if (fabs(m[1][1]-m11)>eps) return(0);

      return(1);
   }

   return(0);
}

//-----------------------------------------------------------------------------------------------

void M22_Print (const real8 *m, const char *title)
{
   if (m)
   {
      if (title) printf("<m2x2 title=\"%s\">\n",title);
      else       printf("<m2x2>\n");

      printf("%7.4lf,%7.4lf,"  , m[0], m[1] );
      printf("%7.4lf,%7.4lf,\n", m[2], m[3] );

      printf("</m2x2>\n");
   }
}

