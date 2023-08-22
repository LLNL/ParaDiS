#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include "mpi_portability.h"

#include "Typedefs.h"
#include "FFT.h"
#include "Matrix_t.h"
#include "Vector_t.h"

//---------------------------------------------------------------------------------------------------------

#define MATRIX_2x2(m00,m01,m10,m11)                        \
{                                                          \
   Allocate(2,2);                                          \
                                                           \
   real8 *m = mtx;                                         \
                                                           \
   if (m)                                                  \
   {                                                       \
      m[0]=(m00); m[1]=(m01);                              \
      m[2]=(m10); m[3]=(m11);                              \
   }                                                       \
}

#define MATRIX_3x3(m00,m01,m02,m10,m11,m12,m20,m21,m22)    \
{                                                          \
   Allocate(3,3);                                          \
                                                           \
   real8 *m = mtx;                                         \
                                                           \
   if (m)                                                  \
   {                                                       \
      m[0]=(m00); m[1]=(m01); m[2]=(m02);                  \
      m[3]=(m10); m[4]=(m11); m[5]=(m12);                  \
      m[6]=(m20); m[7]=(m21); m[8]=(m22);                  \
   }                                                       \
}

#define MATRIX_4x4(m00,m01,m02,m03,                        \
                   m10,m11,m12,m13,                        \
                   m20,m21,m22,m23,                        \
                   m30,m31,m32,m33 )                       \
{                                                          \
   Allocate(4,4);                                          \
                                                           \
   real8 *m = mtx;                                         \
                                                           \
   if (m)                                                  \
   {                                                       \
      m[ 0]=(m00); m[ 1]=(m01); m[ 2]=(m02); m[ 3]=(m03);  \
      m[ 4]=(m10); m[ 5]=(m11); m[ 6]=(m12); m[ 7]=(m13);  \
      m[ 8]=(m20); m[ 9]=(m21); m[10]=(m22); m[11]=(m23);  \
      m[12]=(m30); m[13]=(m31); m[14]=(m32); m[15]=(m33);  \
   }                                                       \
}

#define MATRIX_FOR_ALL(fn)                       \
{                                                \
   if (mtx && (mtx_w>0) && (mtx_h>0))            \
   {                                             \
      for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) \
      { fn; }                                    \
   }                                             \
}

//---------------------------------------------------------------------------------------------------------

static void _init(Matrix_t &m)
{
   m.mtx_w         = 0;
   m.mtx_h         = 0;
   m.mtx           = 0;
   m.mtx_rows      = 0;
   m.mtx_piv       = 0;

#ifdef PARALLEL
   m.mtx_mpi_cnt   = 0;
   m.mtx_mpi_bytes = 0;
   m.mtx_mpi       = 0;
#endif
}

Matrix_t::Matrix_t(void)                                               { _init(*this); }
Matrix_t::Matrix_t(const int w, const int h)                           { _init(*this); Allocate(w,h);          }
Matrix_t::Matrix_t(const int w, const int h, const int              s) { _init(*this); Allocate(w,h); *this=s; }
Matrix_t::Matrix_t(const int w, const int h, const real8            s) { _init(*this); Allocate(w,h); *this=s; }
Matrix_t::Matrix_t(const int w, const int h, const unsigned char   *s) { _init(*this); Allocate(w,h); *this=s; }
Matrix_t::Matrix_t(const int w, const int h, const unsigned short  *s) { _init(*this); Allocate(w,h); *this=s; }
Matrix_t::Matrix_t(const int w, const int h, const          short  *s) { _init(*this); Allocate(w,h); *this=s; }
Matrix_t::Matrix_t(const int w, const int h, const unsigned int    *s) { _init(*this); Allocate(w,h); *this=s; }
Matrix_t::Matrix_t(const int w, const int h, const          int    *s) { _init(*this); Allocate(w,h); *this=s; }
Matrix_t::Matrix_t(const int w, const int h, const          real8  *s) { _init(*this); Allocate(w,h); *this=s; }
Matrix_t::Matrix_t(const Matrix_t & m)                                 { _init(*this);                *this=m; }

Matrix_t::Matrix_t(const int w, const int h, real8 (*fn)(const int i, const int j))
{
   _init(*this); Allocate(w,h);

   real8 **m = mtx_rows;

   if (m)
   {
      for (int j=0; (j<h); ++j)
      for (int i=0; (i<w); ++i)
         m[j][i]=(*fn)(j,i);
   }
}

Matrix_t::Matrix_t
(
   const real8 m00, const real8 m01,
   const real8 m10, const real8 m11
)
{
   _init(*this); Allocate(2,2);

   real8 *m = mtx;

   if (m)
   {
      m[0]=m00; m[1]=m01;
      m[2]=m10; m[3]=m11;
   }
}

Matrix_t::Matrix_t
(
   const real8 m00, const real8 m01, const real8 m02,
   const real8 m10, const real8 m11, const real8 m12
)
{
   _init(*this); Allocate(3,2);

   real8 *m = mtx;

   if (m)
   {
      m[0]=m00; m[1]=m01; m[2]=m02;
      m[3]=m10; m[4]=m11; m[5]=m12;
   }
}

Matrix_t::Matrix_t
(
   const real8 m00, const real8 m01, const real8 m02,
   const real8 m10, const real8 m11, const real8 m12,
   const real8 m20, const real8 m21, const real8 m22
)
{
   _init(*this); Allocate(3,3);

   real8 *m = mtx;

   if (m)
   {
      m[0]=m00; m[1]=m01; m[2]=m02;
      m[3]=m10; m[4]=m11; m[5]=m12;
      m[6]=m20; m[7]=m21; m[8]=m22;
   }
}

Matrix_t::Matrix_t
(
   const real8 m00, const real8 m01, const real8 m02, const real8 m03,
   const real8 m10, const real8 m11, const real8 m12, const real8 m13,
   const real8 m20, const real8 m21, const real8 m22, const real8 m23
)
{
   _init(*this); Allocate(4,3);

   real8 *m = mtx;

   if (m)
   {
      m[0]=m00; m[1]=m01; m[ 2]=m02; m[ 3]=m03;
      m[4]=m10; m[5]=m11; m[ 6]=m12; m[ 7]=m13;
      m[8]=m20; m[9]=m21; m[10]=m22; m[11]=m23;
   }
}

Matrix_t::Matrix_t
(
   const real8 m00, const real8 m01, const real8 m02, const real8 m03,
   const real8 m10, const real8 m11, const real8 m12, const real8 m13,
   const real8 m20, const real8 m21, const real8 m22, const real8 m23,
   const real8 m30, const real8 m31, const real8 m32, const real8 m33
)
{
   _init(*this); Allocate(4,4);

   real8 *m = mtx;

   if (m)
   {
      m[ 0]=m00; m[ 1]=m01; m[ 2]=m02; m[ 3]=m03;
      m[ 4]=m10; m[ 5]=m11; m[ 6]=m12; m[ 7]=m13;
      m[ 8]=m20; m[ 9]=m21; m[10]=m22; m[11]=m23;
      m[12]=m30; m[13]=m31; m[14]=m32; m[15]=m33;
   }
}

Matrix_t::Matrix_t
(
   const real8 m00, const real8 m01, const real8 m02, const real8 m03, const real8 m04, const real8 m05, const real8 m06, const real8 m07,
   const real8 m10, const real8 m11, const real8 m12, const real8 m13, const real8 m14, const real8 m15, const real8 m16, const real8 m17,
   const real8 m20, const real8 m21, const real8 m22, const real8 m23, const real8 m24, const real8 m25, const real8 m26, const real8 m27,
   const real8 m30, const real8 m31, const real8 m32, const real8 m33, const real8 m34, const real8 m35, const real8 m36, const real8 m37,
   const real8 m40, const real8 m41, const real8 m42, const real8 m43, const real8 m44, const real8 m45, const real8 m46, const real8 m47,
   const real8 m50, const real8 m51, const real8 m52, const real8 m53, const real8 m54, const real8 m55, const real8 m56, const real8 m57,
   const real8 m60, const real8 m61, const real8 m62, const real8 m63, const real8 m64, const real8 m65, const real8 m66, const real8 m67,
   const real8 m70, const real8 m71, const real8 m72, const real8 m73, const real8 m74, const real8 m75, const real8 m76, const real8 m77
)
{
   _init(*this); Allocate(8,8);

   real8 **m = mtx_rows;

   if (m)
   {
      m[0][0]=m00; m[0][1]=m01; m[0][2]=m02; m[0][3]=m03; m[0][4]=m04; m[0][5]=m05; m[0][6]=m06; m[0][7]=m07;
      m[1][0]=m10; m[1][1]=m11; m[1][2]=m12; m[1][3]=m13; m[1][4]=m14; m[1][5]=m15; m[1][6]=m16; m[1][7]=m17;
      m[2][0]=m20; m[2][1]=m21; m[2][2]=m22; m[2][3]=m23; m[2][4]=m24; m[2][5]=m25; m[2][6]=m26; m[2][7]=m27;
      m[3][0]=m30; m[3][1]=m31; m[3][2]=m32; m[3][3]=m33; m[3][4]=m34; m[3][5]=m35; m[3][6]=m36; m[3][7]=m37;
      m[4][0]=m40; m[4][1]=m41; m[4][2]=m42; m[4][3]=m43; m[4][4]=m44; m[4][5]=m45; m[4][6]=m46; m[4][7]=m47;
      m[5][0]=m50; m[5][1]=m51; m[5][2]=m52; m[5][3]=m53; m[5][4]=m54; m[5][5]=m55; m[5][6]=m56; m[5][7]=m57;
      m[6][0]=m60; m[6][1]=m61; m[6][2]=m62; m[6][3]=m63; m[6][4]=m64; m[6][5]=m65; m[6][6]=m66; m[6][7]=m67;
      m[7][0]=m70; m[7][1]=m71; m[7][2]=m72; m[7][3]=m73; m[7][4]=m74; m[7][5]=m75; m[7][6]=m76; m[7][7]=m77;
   }
}

Matrix_t::~Matrix_t() { Recycle(); }

void Matrix_t::Recycle (void)
{
   if (mtx     ) { delete [] mtx     ; mtx     =0; }
   if (mtx_rows) { delete [] mtx_rows; mtx_rows=0; }
   if (mtx_piv ) { delete [] mtx_piv ; mtx_piv =0; }
#ifdef PARALLEL
   if (mtx_mpi ) { delete [] mtx_mpi ; mtx_mpi =0; mtx_mpi_cnt=0; mtx_mpi_bytes=0; }
#endif

   mtx_w = 0;
   mtx_h = 0;
}

void Matrix_t::Allocate (const int w, const int h)
{
   if (mtx_piv) { delete [] mtx_piv; mtx_piv=0; }

   if ( (w>0) && (h>0) )
   {
      if ( ((w*h)>(mtx_w*mtx_h)) && mtx     ) { delete [] mtx     ; mtx     =0; }
      if ( ((  h)>(      mtx_h)) && mtx_rows) { delete [] mtx_rows; mtx_rows=0; }

      if (!mtx     ) { mtx      = new real8 [w*h]; }
      if (!mtx_rows) { mtx_rows = new real8*[  h]; }

      if (mtx)             { memset(mtx,0,w*h*sizeof(real8) ); }
      if (mtx && mtx_rows) { for (int i=0,k=0; (i<h); ++i, k+=w)  { mtx_rows[i]=mtx+k; } }
   }

   mtx_w = (mtx ? w : 0 );
   mtx_h = (mtx ? h : 0 );
}

real8 **Matrix_t::Offset_Matrix(const int ioff, const int joff)
{
   real8 **r = ( (mtx && (mtx_h>0)) ? new real8 *[mtx_h] : 0 );

   if (r)
   { for (int i=0,k=-joff; (i<mtx_h); i++, k+=mtx_w) { r[i]=mtx+k; } }

   return(r-ioff);
}

const Matrix_t & Matrix_t::operator = (const Matrix_t &m)
{
   Allocate(m.mtx_w,m.mtx_h);

   if (mtx && mtx_rows && m.mtx && m.mtx_rows)
   {
      size_t bytes = mtx_w * sizeof(real8);

      for (int j=0; (j<mtx_h); ++j)
         memcpy(mtx_rows[j],m.mtx_rows[j],bytes);
   }

   return(*this);
}

const Matrix_t & Matrix_t::operator = (const unsigned char   *m) { if (mtx && m) { for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) mtx[i]= (real8) m[i]; } return(*this); }
const Matrix_t & Matrix_t::operator = (const unsigned short  *m) { if (mtx && m) { for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) mtx[i]= (real8) m[i]; } return(*this); }
const Matrix_t & Matrix_t::operator = (const          short  *m) { if (mtx && m) { for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) mtx[i]= (real8) m[i]; } return(*this); }
const Matrix_t & Matrix_t::operator = (const unsigned int    *m) { if (mtx && m) { for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) mtx[i]= (real8) m[i]; } return(*this); }
const Matrix_t & Matrix_t::operator = (const          int    *m) { if (mtx && m) { for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) mtx[i]= (real8) m[i]; } return(*this); }
const Matrix_t & Matrix_t::operator = (const          real8  *m) { if (mtx && m) { for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) mtx[i]= (real8) m[i]; } return(*this); }

const Matrix_t & Matrix_t::operator = (const          real8   m) { if (mtx     ) { for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) mtx[i]= (real8) m;    } return(*this); }

int Matrix_t::operator == (const Matrix_t & m) const { return( Compare(m) ? 1 : 0 ); }
int Matrix_t::operator == (const real8      s) const { return( Compare(s) ? 1 : 0 ); }
int Matrix_t::operator == (const real8     *s) const { return( Compare(s) ? 1 : 0 ); }
int Matrix_t::operator != (const Matrix_t & m) const { return( Compare(m) ? 0 : 1 ); }
int Matrix_t::operator != (const real8      s) const { return( Compare(s) ? 0 : 1 ); }
int Matrix_t::operator != (const real8     *s) const { return( Compare(s) ? 0 : 1 ); }

int Matrix_t::operator < (const Matrix_t & m) const
{
   if (mtx && m.mtx && (mtx_w>0) && (mtx_h>0) && (mtx_w==m.mtx_w) && (mtx_h==m.mtx_h))
   {
      for (int i=0,n=(mtx_w*mtx_h); (i<n); ++i)
         if (mtx[i]>=m.mtx[i]) return(0);

      return(1);
   }

   return(0);
}

int Matrix_t::operator < (const real8 s) const
{
   if (mtx && (mtx_w>0) && (mtx_h>0))
   {
      for (int i=0,n=(mtx_w*mtx_h); (i<n); ++i)
         if (mtx[i]>=s) return(0);

      return(1);
   }

   return(0);
}

int Matrix_t::operator < (const real8 *s) const
{
   if (mtx && s && (mtx_w>0) && (mtx_h>0))
   {
      for (int i=0,n=(mtx_w*mtx_h); (i<n); ++i)
         if (mtx[i]>=s[i]) return(0);

      return(1);
   }

   return(0);
}

int Matrix_t::operator <= (const Matrix_t & m) const
{
   if (mtx && m.mtx && (mtx_w>0) && (mtx_h>0) && (mtx_w==m.mtx_w) && (mtx_h==m.mtx_h))
   {
      for (int i=0,n=(mtx_w*mtx_h); (i<n); ++i)
         if (mtx[i]>m.mtx[i]) return(0);

      return(1);
   }

   return(0);
}

int Matrix_t::operator <= (const real8 s) const
{
   if (mtx && (mtx_w>0) && (mtx_h>0))
   {
      for (int i=0,n=(mtx_w*mtx_h); (i<n); ++i)
         if (mtx[i]>s) return(0);

      return(1);
   }

   return(0);
}

int Matrix_t::operator <= (const real8 *s) const
{
   if (mtx && s && (mtx_w>0) && (mtx_h>0))
   {
      for (int i=0,n=(mtx_w*mtx_h); (i<n); ++i)
         if (mtx[i]>s[i]) return(0);

      return(1);
   }

   return(0);
}

int Matrix_t::operator > (const Matrix_t & m) const
{
   if (mtx && m.mtx && (mtx_w>0) && (mtx_h>0) && (mtx_w==m.mtx_w) && (mtx_h==m.mtx_h))
   {
      for (int i=0,n=(mtx_w*mtx_h); (i<n); ++i)
         if (mtx[i]<=m.mtx[i]) return(0);

      return(1);
   }

   return(0);
}

int Matrix_t::operator > (const real8 s) const
{
   if (mtx && (mtx_w>0) && (mtx_h>0))
   {
      for (int i=0,n=(mtx_w*mtx_h); (i<n); ++i)
         if (mtx[i]<=s) return(0);

      return(1);
   }

   return(0);
}

int Matrix_t::operator > (const real8 *s) const
{
   if (mtx && s && (mtx_w>0) && (mtx_h>0))
   {
      for (int i=0,n=(mtx_w*mtx_h); (i<n); ++i)
         if (mtx[i]<=s[i]) return(0);

      return(1);
   }

   return(0);
}

int Matrix_t::operator >= (const Matrix_t & m) const
{
   if (mtx && m.mtx && (mtx_w>0) && (mtx_h>0) && (mtx_w==m.mtx_w) && (mtx_h==m.mtx_h))
   {
      for (int i=0,n=(mtx_w*mtx_h); (i<n); ++i)
         if (mtx[i]<m.mtx[i]) return(0);

      return(1);
   }

   return(0);
}

int Matrix_t::operator >= (const real8 s) const
{
   if (mtx && (mtx_w>0) && (mtx_h>0))
   {
      for (int i=0,n=(mtx_w*mtx_h); (i<n); ++i)
         if (mtx[i]<s) return(0);

      return(1);
   }

   return(0);
}

int Matrix_t::operator >= (const real8 *s) const
{
   if (mtx && s && (mtx_w>0) && (mtx_h>0))
   {
      for (int i=0,n=(mtx_w*mtx_h); (i<n); ++i)
         if (mtx[i]<s[i]) return(0);

      return(1);
   }

   return(0);
}

Matrix_t Matrix_t::operator + (const Matrix_t & m) { Matrix_t tmp(*this); tmp+=m; return(tmp); }
Matrix_t Matrix_t::operator + (const real8      s) { Matrix_t tmp(*this); tmp+=s; return(tmp); }

Matrix_t Matrix_t::operator - (const Matrix_t & m) { Matrix_t tmp(*this); tmp-=m; return(tmp); }
Matrix_t Matrix_t::operator - (const real8      s) { Matrix_t tmp(*this); tmp-=s; return(tmp); }

Matrix_t Matrix_t::operator * (const Matrix_t & m) { Matrix_t tmp; tmp.Multiply(*this,m); return(tmp); }
Matrix_t Matrix_t::operator * (const real8      s) { Matrix_t tmp(*this); tmp*=s; return(tmp); }

Matrix_t Matrix_t::operator / (const real8      s) { real8 sx = ( (fabs(s)>1.0e-14) ? (1.0/s) : 0.0 );
                                                     Matrix_t tmp(*this); tmp*=sx; return(tmp); }

void  Matrix_t::operator += (const Matrix_t  & m)
{
   int w = ( (mtx_w==m.mtx_w) ? mtx_w : 0 );
   int h = ( (mtx_h==m.mtx_h) ? mtx_h : 0 );

   real8 *ma =   mtx;
   real8 *mb = m.mtx;

   if ( (w>0) && (h>0) && ma && mb)
   {
      for (int i=0; (i<(mtx_w*mtx_h)); ++i)
         ma[i]+=mb[i];
   }
}

void  Matrix_t::operator -= (const Matrix_t  & m)
{
   int w = ( (mtx_w==m.mtx_w) ? mtx_w : 0 );
   int h = ( (mtx_h==m.mtx_h) ? mtx_h : 0 );

   real8 *ma =   mtx;
   real8 *mb = m.mtx;

   if ( (w>0) && (h>0) && ma && mb)
   {
      for (int i=0; (i<(mtx_w*mtx_h)); ++i)
         ma[i]-=mb[i];
   }
}

void Matrix_t::operator *= (const Matrix_t  & m)
{
   Matrix_t tmp(*this);

   (*this).Multiply(tmp,m);
}

void Matrix_t::operator += (const real8  s) { if (mtx     ) { for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) mtx[i]+=s   ;     }}
void Matrix_t::operator += (const real8 *s) { if (mtx && s) { for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) mtx[i]+=s[i];     }}
void Matrix_t::operator -= (const real8  s) { if (mtx     ) { for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) mtx[i]-=s   ;     }}
void Matrix_t::operator -= (const real8 *s) { if (mtx && s) { for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) mtx[i]-=s[i];     }}

void Matrix_t::operator *= (const real8  s) { if (mtx     ) { for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) mtx[i]*=s   ;     }}
void Matrix_t::operator *= (const real8 *s) { if (mtx && s) { for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) mtx[i]*=s[i];     }}

void Matrix_t::operator /= (const real8  s) { if (mtx     ) { real8 sx = ( (fabs(s)>1.0e-12) ? 1.0/s : 0 );
                                                                for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) mtx[i]*=sx;       }}

Matrix_t Matrix_t::operator - ()
{
   Matrix_t tmp(*this);

   if (mtx) { for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) mtx[i]=-mtx[i];  }

   return(tmp);
}

void Matrix_t::Add      (const Matrix_t & m)                         { *this+=m; }
void Matrix_t::Add      (const Matrix_t & ma, const Matrix_t & mb) { *this=ma; *this+=mb; }

void Matrix_t::Subtract (const Matrix_t & m)                         { *this-=m; }
void Matrix_t::Subtract (const Matrix_t & ma, const Matrix_t & mb) { *this=ma; *this-=mb; }

real8 *Matrix_t::Get_Row (const int j) const
{
   real8 *row = ( (mtx_w>0) ? new real8[mtx_w] : 0 );

   if (row) { memcpy(row,mtx_rows[j],mtx_w*sizeof(real8)); }

   return(row);
}

real8 *Matrix_t::Get_Column (const int i) const
{
   real8 *col = ( (mtx_h>0) ? new real8[mtx_h] : 0 );

   if (col)
   {
      for (int j=0; (j<mtx_h); ++j)
      { col[j] = mtx_rows[j][i]; }
   }

   return(col);
}

real8 *Matrix_t::Get_Diagonal (void) const
{
   real8 *diag = ( (mtx && (mtx_h>0)) ? new real8[mtx_h] : 0 );

   if (mtx && diag)
   { for (int i=0,j=0; (i<mtx_h); i++, j+=(mtx_w+1)) { diag[i]=mtx[j]; } }

   return(diag);
}

void Matrix_t::Set_Diagonal (const real8 *diag)
{
   if (mtx && diag)
   { for (int i=0,j=0; (i<mtx_h); i++, j+=(mtx_w+1)) { mtx[j]=diag[i]; } }
}

void Matrix_t::Swap_Rows (const int i, const int j)
{
   int    w  = mtx_w;
   int    h  = mtx_h;
   real8 *r0 = ( mtx && (0<=i) && (i<h) ? mtx+(i*w) : 0 );
   real8 *r1 = ( mtx && (0<=j) && (j<h) ? mtx+(j*w) : 0 );

   if (r0 && r1 && (r0!=r1))
      for (int k=0; (k<w); ++k) { real8 s=r0[k]; r0[k]=r1[k]; r1[k]=s; }
}

void Matrix_t::Swap_Columns (const int i, const int j)
{
   int    w  = mtx_w;
   int    h  = mtx_h;
   real8 *c0 = ( mtx && (0<=i) && (i<w) ? mtx+i : 0 );
   real8 *c1 = ( mtx && (0<=j) && (j<w) ? mtx+j : 0 );

   if (c0 && c1 && (c0!=c1))
      for (int k=0,m=0; (k<h); ++k,m+=w) { real8 s=c0[m]; c0[m]=c1[m]; c1[m]=s; }
}

void Matrix_t::Sum_Rows (const int i, const int j)
{
   int    w  = mtx_w;
   int    h  = mtx_h;
   real8 *r0 = ( mtx && (0<=i) && (i<h) ? mtx+(i*w) : 0 );
   real8 *r1 = ( mtx && (0<=j) && (j<h) ? mtx+(j*w) : 0 );

   if (r0 && r1 && (r0!=r1))
      for (int k=0; (k<w); ++k) { r0[k]+=r1[k]; }
}

void Matrix_t::Sum_Rows (const int i, const int j, const real8 s)
{
   int    w  = mtx_w;
   int    h  = mtx_h;
   real8 *r0 = ( mtx && (0<=i) && (i<h) ? mtx+(i*w) : 0 );
   real8 *r1 = ( mtx && (0<=j) && (j<h) ? mtx+(j*w) : 0 );

   if (r0 && r1 && (r0!=r1))
      for (int k=0; (k<w); ++k) { r0[k]+=r1[k]*s; }
}

real8 Matrix_t::Row_Min (const int j, int & indx)
{
   int    w    = mtx_w;
   int    h    = mtx_h;
   real8 *row  = ( (mtx_rows && (0<=j) && (j<h)) ? mtx_rows[j] : 0 );
   real8  min  = 0.0;
          indx = 0;

   if (row)
   {
      min = row[0];
      for (int i=1; (i<w); ++i)
         if (row[i]<min) { min=row[i]; indx=i; }
   }

   return(min);
}

real8 Matrix_t::Row_Max (const int j, int & indx)
{
   int    w    = mtx_w;
   int    h    = mtx_h;
   real8 *row  = ( (mtx_rows && (0<=j) && (j<h)) ? mtx_rows[j] : 0 );
   real8  max  = 0.0;
          indx = 0;

   if (row)
   {
      max = row[0];
      for (int i=1; (i<w); ++i)
         if (row[i]>max) { max=row[i]; indx=i; }
   }

   return(max);
}

real8 Matrix_t::Column_Min (const int i, int & indx)
{
   int    w    = mtx_w;
   int    h    = mtx_h;
   real8  min  = 0.0;
          indx = 0;

   if (mtx_rows && (0<=i) && (i<w) )
   {
      real8 **m = mtx_rows;

      min = m[0][i];
      for (int j=1; (j<h); ++j)
         if (m[j][i]<min) { min=m[j][i]; indx=j; }
   }

   return(min);
}

real8 Matrix_t::Column_Max (const int i, int & indx)
{
   int    w    = mtx_w;
   int    h    = mtx_h;
   real8  max  = 0.0;
          indx = 0;

   if (mtx_rows && (0<=i) && (i<w) )
   {
      real8 **m = mtx_rows;

      max = m[0][i];
      for (int j=1; (j<h); ++j)
         if (m[j][i]>max) { max=m[j][i]; indx=j; }
   }

   return(max);
}

real8 Matrix_t::Dot (const Matrix_t & m)
{
   real8 dot = 0.0;

   int w = ( (mtx_w==m.mtx_w) ? mtx_w : 0 );
   int h = ( (mtx_h==m.mtx_h) ? mtx_h : 0 );

   real8 *ma =   mtx;
   real8 *mb = m.mtx;

   if ( (w>0) && (h>0) && ma && mb)
   {
      for (int i=0; (i<mtx_w*mtx_h); ++i)
         dot += ma[i]*mb[i];
   }

   return(dot);
}

void Matrix_t::Tensor (const Matrix_t & ma, const Matrix_t & mb)
{
   int w = ma.mtx_w*mb.mtx_w;
   int h = ma.mtx_h*mb.mtx_h;

   Allocate(w,h);

   if (   mtx      && ma.mtx      && mb.mtx
       && mtx_rows && ma.mtx_rows && mb.mtx_rows )
   {
      int p = mb.mtx_h;
      int q = mb.mtx_w;

      real8 **c =    mtx_rows;
      real8 **a = ma.mtx_rows;
      real8 **b = mb.mtx_rows;

      for (int j=0; (j<h); j++)
      for (int i=0; (i<w); i++)
         c[j][i] = a[j/p][i/q] * b[j%p][i%q];
   }
}

void Matrix_t::Multiply (const Matrix_t & m)
{
   Matrix_t tmp(*this);

   Multiply(tmp,m);
}

static real8 rowcol_mult(real8 *r, real8 *c, int n, int m)
{
   real8 s = 0.0;

   for (int i=0,j=0; (i<n); ++i, j+=m)
      s += r[i]*c[j];

   return(s);
}

void Matrix_t::Multiply (const Matrix_t & ma, const Matrix_t & mb)
{
   int nw = mb.mtx_w;
   int nh = ma.mtx_h;

   Allocate(nw,nh);

   real8 **m = mtx_rows;

   for (int j=0; (j<nh); ++j)
   {
      real8 *row = ma.mtx_rows[j];
      real8 *col = mb.mtx_rows[0];

      for (int i=0; (i<nw); ++i, ++col)
         m[j][i] = rowcol_mult(row,col,ma.mtx_w,mb.mtx_w);
   }
}

void Matrix_t::Linear_Sum
(
   const Matrix_t & ma, const Matrix_t & mx,
   const Matrix_t & mb, const Matrix_t & my
)
{
   Matrix_t tmpa; tmpa.Multiply(ma,mx);
   Matrix_t tmpb; tmpb.Multiply(mb,my);

   *this  = tmpa;
   *this += tmpb;
}

void Matrix_t::Transpose (void)
{
   Matrix_t tmp(*this);

   Transpose(tmp);
}

void Matrix_t::Transpose (const Matrix_t & m)
{
   int h = m.mtx_w;
   int w = m.mtx_h;

   if (this!=&m)
   {
      Allocate(w,h);

      real8 **p = m.mtx_rows;
      real8 **q =   mtx_rows;

      for (int j=0; (j<h); ++j)
      for (int i=0; (i<w); ++i)
         q[j][i] = p[i][j];
   }
   else { Transpose(); }
}

void Matrix_t::Inverse (void)
{
   Matrix_t tmp(*this);

   int     h   = mtx_h;
   real8 **a   = tmp.mtx_rows;
   real8 **ai  =     mtx_rows;
   int    *p   = (int *) ( (h>0) ? new int[h+1] : 0 );
   real8   tol = 1.0e-8;

   ::LUP_Decompose(   a,p,h,tol);
   ::LUP_Invert   (ai,a,p,h);

   if (p) { delete [] p; }
}

void Matrix_t::Inverse (const Matrix_t & m)
{
   *this=m;

   Inverse();
}

static int nearly_eq(const double a, const double b, const double eps)
{
   if (a==b) return(1);                     // (bit-for-bit equal)

   // if either arg is near zero, can't use absolute tolerance

   if (   (fabs(a)<1.0e-12)
       || (fabs(b)<1.0e-12) )
      return(  (fabs(a-b)<eps) ? 1 : 0 );

   // otherwise compute relative error by scaling by largest arg...

   double rerr = ( (fabs(a)>fabs(b)) ?  fabs((a-b)/a) : fabs((a-b)/b) );

   return( (rerr<eps) ? 1 : 0 );
}

int Matrix_t::Compare (const real8 s, const real8 eps) const
{
   int w = mtx_w;
   int h = mtx_h;

   if (mtx && (w>0) && (h>0) )
   {
      for (int i=0, n=(w*h); (i<n); ++i)
         if (!nearly_eq(mtx[i],s,eps)) return(0);

      return(1);
   }

   return(0);
}

int Matrix_t::Compare (const real8 *s, const real8 eps) const
{
   int w = mtx_w;
   int h = mtx_h;

   if (mtx && mtx_rows && s && (w>0) && (h>0))
   {
      real8 **m = mtx_rows;

      for (int j=0,k=0; (j<h); ++j)
      for (int i=0    ; (i<w); ++i,++k)
         if (!nearly_eq(m[j][i],s[k],eps)) return(0);

      return(1);
   }

   return(0);
}

int Matrix_t::Compare (const Matrix_t & m, const real8 eps) const
{
   if (mtx_w!=m.mtx_w) return(0);
   if (mtx_h!=m.mtx_h) return(0);

   if (mtx && mtx_rows && m.mtx && m.mtx_rows)
   {
      real8 **ma =   mtx_rows;
      real8 **mb = m.mtx_rows;

      for (int j=0; (j<mtx_h); ++j)
      for (int i=0; (i<mtx_w); ++i)
         if (!nearly_eq(ma[j][i],mb[j][i],eps)) return(0);

      return(1);
   }

   return(0);
}

void Matrix_t::Identity (void)
{
   if (mtx && mtx_rows)
   {
      for (int j=0,k=0; (j<mtx_h); ++j, k+=mtx_w) { mtx_rows[j]=mtx+k; }

      memset(mtx,0,mtx_w*mtx_h*sizeof(real8) );

      int n = ( (mtx_w<mtx_h) ? mtx_w : mtx_h );

      for (int i=0,k=0; (i<n); ++i, k+=(mtx_w+1) ) { mtx[k]=1.0; }
   }
}

void Matrix_t::Identity (const int w, const int h)
{
   Allocate(w,h);
   Identity();
}

void Matrix_t::Pascal (void)
{
   if (mtx)
   {
      int     w = mtx_w;
      int     h = mtx_h;
      real8 **m = mtx_rows;

      for (int i=0; (i<w); ++i)  { m[0][i]=1.0; }
      for (int j=0; (j<h); ++j)  { m[j][0]=1.0; }

      for (int j=1; (j<h); ++j)
      for (int i=1; (i<w); ++i)
         m[j][i]=m[j][i-1]+m[j-1][i];
   }
}

void Matrix_t::Pascal (const int w, const int h)
{
   Allocate(w,h);
   Pascal();
}

void Matrix_t::Hilbert (void)
{
   if (mtx)
   {
      int     w = mtx_w;
      int     h = mtx_h;
      real8 **m = mtx_rows;

      for (int j=0; (j<h); ++j)
      for (int i=0; (i<w); ++i)
         m[j][i]=1.0/(i+j+1);
   }
}

void Matrix_t::Hilbert (const int w, const int h)
{
   Allocate(w,h);
   Hilbert();
}

void Matrix_t::DCT (void)
{
   int w = mtx_w;
   int h = mtx_h;

   if (mtx && (w>0) && (h>0))
   {
      int k=0;
      for (int i=0; (i<w); ++i, ++k) { mtx[k]=(1.0/sqrt(h)); }

      for (int j=1; (j<h); ++j)
         for (int i=0; (i<w); ++i,++k)
            mtx[k] = (sqrt(2.0/h)*cos(((2*i+1)*j*M_PI)/(2*8)));
   }
}

void Matrix_t::Walsh (const int n)
{
   if (n>0)
   {
      Matrix_t m2(1,1,1,-1);

      *this = m2;

      for (int i=2; (i<=n); i++)
      {
         Matrix_t tmp(*this);

         Tensor(m2,tmp);
      }
   }
}

void Matrix_t::DCT (const int w, const int h)
{
   Allocate(w,h);
   DCT();
}

void Matrix_t::Rand (const real8 min, const real8 max)
{
   if (mtx)
   {
      real8 *m = mtx;
      int    n = mtx_w * mtx_h;
      real8  s = ((real8)(max-min))/RAND_MAX;

      for (int i=0; (i<n); ++i)
            m[i]= min+(s*rand());
   }
}

void Matrix_t::Rand (const int w, const int h, const real8 min, const real8 max)
{
   Allocate(w,h);
   Rand(min,max);
}

void Matrix_t::Abs  (void) { if (mtx) { for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) mtx[i] = fabs(mtx[i]); } }
void Matrix_t::Rint (void) { if (mtx) { for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i) mtx[i] = rint(mtx[i]); } }

void Matrix_t::Ramp (                          const real8 x0, const real8 dx) { if (mtx) { real8 x=x0; for (int i=0, n=(mtx_w*mtx_h); (i<n); ++i ) { mtx[i]=x; x+=dx; } } }
void Matrix_t::Ramp (const int w, const int h, const real8 x0, const real8 dx) { Allocate(w,h); Ramp(x0,dx); }

real8 Matrix_t::Min (void) const
{
   real8   min=0.0;

   if (mtx)
   {
      min = mtx[0];

      for (int i=0; (i<(mtx_w*mtx_h)); ++i)
         min = ( (mtx[i]<min) ? mtx[i] : min );
   }

   return(min);
}

real8 Matrix_t::Min (int & iy, int & ix) const
{
   ix=iy=0;

   real8 min=0.0;

   if (mtx)
   {
      min = mtx[0];

      for (int i=0; (i<(mtx_w*mtx_h)); ++i)
         if (mtx[i]<min) { min=mtx[i]; iy=(i/mtx_w); ix=(i%mtx_w); }
   }

   return(min);
}

real8 Matrix_t::Abs_Min (void) const
{
   real8   min=0.0;

   if (mtx)
   {
      min = fabs(mtx[0]);

      for (int i=0; (i<(mtx_w*mtx_h)); ++i)
         min = ( (fabs(mtx[i])<min) ? fabs(mtx[i]) : min );
   }

   return(min);
}

real8 Matrix_t::Abs_Min (int & iy, int & ix) const
{
   ix=iy=0;

   real8 min=0.0;

   if (mtx)
   {
      min = fabs(mtx[0]);

      for (int i=0; (i<(mtx_w*mtx_h)); ++i)
         if (fabs(mtx[i])<min) { min=fabs(mtx[i]); iy=(i/mtx_w); ix=(i%mtx_w); }
   }

   return(min);
}

real8 Matrix_t::Max (void) const
{
   real8 max=0.0;

   if (mtx)
   {
      max = mtx[0];

      for (int i=0; (i<(mtx_w*mtx_h)); ++i)
         max = ( (mtx[i]>max) ? mtx[i] : max );
   }

   return(max);
}

real8 Matrix_t::Max (int & iy, int & ix) const
{
   ix=iy=0;

   real8 max=0.0;

   if (mtx)
   {
      max = mtx[0];

      for (int i=0; (i<(mtx_w*mtx_h)); ++i)
         if (mtx[i]>max) { max=mtx[i]; iy=(i/mtx_w); ix=(i%mtx_w); }
   }

   return(max);
}

real8 Matrix_t::Abs_Max (void) const
{
   real8   max=0.0;

   if (mtx)
   {
      max = fabs(mtx[0]);

      for (int i=0; (i<(mtx_w*mtx_h)); ++i)
         max = ( (fabs(mtx[i])>max) ? fabs(mtx[i]) : max );
   }

   return(max);
}

real8 Matrix_t::Abs_Max (int & iy, int & ix) const
{
   ix=iy=0;

   real8 max=0.0;

   if (mtx)
   {
      max = fabs(mtx[0]);

      for (int i=0; (i<(mtx_w*mtx_h)); ++i)
         if (fabs(mtx[i])>max) { max=fabs(mtx[i]); iy=(i/mtx_w); ix=(i%mtx_w); }
   }

   return(max);
}

real8 Matrix_t::Mean (void) const
{
   real8  mean=0.0;

   if (mtx)
   {
      for (int i=0; (i<(mtx_w*mtx_h)); ++i)
         mean+=mtx[i];

      mean/=(mtx_w*mtx_h);
   }

   return(mean);
}

real8 Matrix_t::Range (void) const
{
   real8  min=0.0;
   real8  max=0.0;

   if (mtx)
   {
      min=max=mtx[0];

      for (int i=0; (i<(mtx_w*mtx_h)); ++i)
      {
         min = ( (mtx[i]<min) ? mtx[i] : min );
         max = ( (mtx[i]>max) ? mtx[i] : max );
      }
   }

   return(max-min);
}

void Matrix_t::Moment (real8 & min, real8 & max, real8 & mean, real8 & range) const
{
   min=max=mean=range=0.0;

   if (mtx)
   {
      min=max=mtx[0];

      for (int i=0; (i<(mtx_w*mtx_h)); ++i)
      {
         min = ( (mtx[i]<min) ? mtx[i] : min );
         max = ( (mtx[i]>max) ? mtx[i] : max );
         mean += mtx[i];
      }

      mean  /= (mtx_w*mtx_h);
      range  = (max-min);
   }
}

void Matrix_t::Moment (real8 & min, real8 & max, real8 & mean, real8 & range, real8 & sdev) const
{
   min=max=mean=range=sdev=0.0;

   Moment(min,max,mean,range);

   if (mtx)
   {
      real8 v=0.0,var=0.0;

      for (int i=0; (i<(mtx_w*mtx_h)); ++i)
      {
         v    = (mean-mtx[i]);
         var += (v*v);
      }

      var  /= (mtx_w*mtx_h);
      sdev  = sqrt(var);
   }
}

// Norm()
//
// Returns the row-sum norm of a matrix.  The norm is computed by summing the absolute
// values of each row member. The norm is the maximum of the row sums.
//---------------------------------------------------------------------------------------------------------

real8 Matrix_t::Norm (void)
{
   real8   norm=0.0;
   real8 **m   = mtx_rows;

   if (m)
   {
      for (int j=0; (j<mtx_h); ++j)
      {
         real8 sum=0.0;

         for (int i=0; (i<mtx_w); ++i)
            sum += fabs(m[j][i]);

         norm = ((sum>norm) ? sum : norm );
      }
   }

   return(norm);
}

// Normalize()
//
// Divides all members of the matrix by the absolute value of the largest member of the matrix.
// All members of the matrix will then be in the range [-1.0:+1.0].
//---------------------------------------------------------------------------------------------------------

void Matrix_t::Normalize (void)
{
   if (mtx)
   {
      real8 max = fabs(mtx[0]);

      for(int i=0; (i<(mtx_w*mtx_h)); ++i)
         max = ( (fabs(mtx[i])>max) ? fabs(mtx[i]) : max );

      if (max>0.0)
      {
         max = (1.0/max);

         for(int i=0; (i<(mtx_w*mtx_h)); ++i)
            mtx[i]*=max;
      }
   }
}

// Max_Norm()
//
// Returns the absolute value of the largest element of the matrix
//---------------------------------------------------------------------------------------------------------

real8 Matrix_t::Max_Norm (void)
{
   real8 max=0.0;

   if (mtx)
   {
      max = fabs(mtx[0]);
      for(int i=0; (i<(mtx_w*mtx_h)); ++i)
         max = ( (fabs(mtx[i])>max) ? fabs(mtx[i]) : max );
   }

   return(max);
}

//---------------------------------------------------------------------------------------------------------

real8 Matrix_t::Condition (void)
{
   Matrix_t tmp;  tmp.Inverse(*this);

   real8 c = tmp.Norm() * Norm();

   return(c);
}

//---------------------------------------------------------------------------------------------------------

void Matrix_t::LUP_Decompose (real8 tol)
{
   if (mtx_piv) { delete [] mtx_piv; }

   mtx_piv = ( (mtx_h>0) ? new int[mtx_h+1] : 0 );

   if (mtx && mtx_rows && mtx_piv)
      ::LUP_Decompose(mtx_rows, mtx_piv, mtx_h, tol);
}

void Matrix_t::LUP_Solve (real8 *x, real8 *b)
{
   if (mtx && mtx_rows && mtx_piv && x && b)
      ::LUP_Solve (mtx_rows, mtx_piv, x, b, mtx_h);
}

void Matrix_t::LUP_Invert (real8 tol)
{
   Matrix_t tmp(*this);

   int     h   = mtx_h;
   real8 **a   = tmp.mtx_rows;
   real8 **ai  =     mtx_rows;
   int    *p   = (int *) ( (h>0) ? new int[h+1] : 0 );

   ::LUP_Decompose(   a,p,h,tol);
   ::LUP_Invert   (ai,a,p,h);

   if (p) { delete [] p; }
}

real8 Matrix_t::LUP_Determinant (void)
{
   real8 det=0.0;

   if (mtx && mtx_rows && mtx_piv)
      det = ::LUP_Determinant (mtx_rows, mtx_piv, mtx_h);

   return(det);
}

//---------------------------------------------------------------------------------------------------------

real8 Matrix_t::Residual (const Matrix_t & m)
{
   int   w = ( (mtx_w==m.mtx_w) ? mtx_w : 0 );
   int   h = ( (mtx_h==m.mtx_h) ? mtx_h : 0 );
   real8 r = 0.0;

   if (mtx && mtx_rows && m.mtx && m.mtx_rows && (w>0) && (h>0))
   {
      real8 **a =   mtx_rows;
      real8 **b = m.mtx_rows;

      for (int j=0; (j<h); ++j)
      for (int i=0; (i<w); ++i)
      { a[j][i] = fabs(a[j][i]-b[j][i]); r+=a[j][i]; }
   }

   return(r);
}

real8 Matrix_t::Residual (const Matrix_t & ma, const Matrix_t & mb)
{
   int   w = ( (ma.mtx_w==mb.mtx_w) ? ma.mtx_w : 0 );
   int   h = ( (ma.mtx_h==mb.mtx_h) ? ma.mtx_h : 0 );
   real8 r = 0.0;

   if ( (w>0) && (h>0) && (this!= &ma) && (this!= &mb) )
      Allocate(w,h);

   real8 **a =    mtx_rows;
   real8 **b = ma.mtx_rows;
   real8 **c = mb.mtx_rows;

   if ( (w>0) && (h>0) && a && b && c )
   {
      for (int j=0; (j<h); ++j)
      for (int i=0; (i<w); ++i)
      { a[j][i] = fabs(b[j][i]-c[j][i]); r+=a[j][i]; }
   }

   return(r);
}

real8 Matrix_t::Residual_Sum (const Matrix_t & m)
{
   int     w    = ( (mtx_w==m.mtx_w) ? mtx_w : 0 );
   int     h    = ( (mtx_h==m.mtx_h) ? mtx_h : 0 );
   real8   rsum = 0.0;
   real8 **a    = ( (  mtx &&   mtx_rows) ?   mtx_rows : 0 );
   real8 **b    = ( (m.mtx && m.mtx_rows) ? m.mtx_rows : 0 );

   if ( a && b && (w>0) && (h>0) )
   {
      for (int j=0; (j<h); ++j)
      for (int i=0; (i<w); ++i)
      { rsum += fabs(a[j][i]-b[j][i]); }
   }

   return(rsum);
}

real8 Matrix_t::Residual_Sum (const Matrix_t & ma, const Matrix_t & mb)
{
   int     w    = ( (ma.mtx_w==mb.mtx_w) ? ma.mtx_w : 0 );
   int     h    = ( (ma.mtx_h==mb.mtx_h) ? ma.mtx_h : 0 );
   real8   rsum = 0.0;
   real8 **a    = ( (ma.mtx && ma.mtx_rows) ? ma.mtx_rows : 0 );
   real8 **b    = ( (mb.mtx && mb.mtx_rows) ? mb.mtx_rows : 0 );

   if ( a && b && (w>0) && (h>0) )
   {
      for (int j=0; (j<h); ++j)
      for (int i=0; (i<w); ++i)
      { rsum += fabs(a[j][i]-b[j][i]); }
   }

   return(rsum);
}

real8 Matrix_t::RMS_Error (const Matrix_t & m) const
{
   real8 rms=0.0;

   int w = ( (mtx_w==m.mtx_w) ? mtx_w : 0 );
   int h = ( (mtx_h==m.mtx_h) ? mtx_h : 0 );

   if (mtx && m.mtx && (w>0) && (h>0))
   {
      for (int i=0; (i<(w*h)); ++i)
      {
         real8 s = (mtx[i]-m.mtx[i]);
         rms += (s*s);
      }

      rms=sqrt(rms)/(w*h);
   }

   return(rms);
}

real8 Matrix_t::RMS_Error (const Matrix_t & ma, const Matrix_t & mb) const { return ( ma.RMS_Error(mb) ); }

void Matrix_t::Threshold (const real8 eps)                  { MATRIX_FOR_ALL( mtx[i] = ( (fabs(mtx[i])<eps) ? 0.0 : mtx[i] )) }
void Matrix_t::Cap       (const real8 max)                  { MATRIX_FOR_ALL( mtx[i] = ( (mtx[i]<max) ? mtx[i] : max )) }
void Matrix_t::Constrain (const real8 min, const real8 max) { MATRIX_FOR_ALL( mtx[i] = ( (mtx[i]<min) ? min : ( (mtx[i]<max) ? mtx[i] : max ))) }

void Matrix_t::Sin   (void) { MATRIX_FOR_ALL( mtx[i]=sin(mtx[i])) }
void Matrix_t::Cos   (void) { MATRIX_FOR_ALL( mtx[i]=cos(mtx[i])) }
void Matrix_t::Tan   (void) { MATRIX_FOR_ALL( mtx[i]=tan(mtx[i])) }

void Matrix_t::Log   (void) { MATRIX_FOR_ALL( mtx[i]=((mtx[i]>0.0) ? log  (mtx[i]) : 0.0 )) }
void Matrix_t::Log2  (void) { MATRIX_FOR_ALL( mtx[i]=((mtx[i]>0.0) ? log2 (mtx[i]) : 0.0 )) }
void Matrix_t::Log10 (void) { MATRIX_FOR_ALL( mtx[i]=((mtx[i]>0.0) ? log10(mtx[i]) : 0.0 )) }
void Matrix_t::Sqrt  (void) { MATRIX_FOR_ALL( mtx[i]=((mtx[i]>0.0) ? sqrt (mtx[i]) : 0.0 )) }

void Matrix_t::Sign  (void) { MATRIX_FOR_ALL( mtx[i]=((mtx[i]<0.0) ? -1.0 : ((mtx[i]>0.0) ? 1.0 : 0.0 )) ) }


void Matrix_t::Identity_2x2 (void)
{
   MATRIX_2x2(1.0,0.0,
              0.0,1.0);
}

void Matrix_t::Identity_3x3 (void)
{
   MATRIX_3x3(1.0,0.0,0.0,
              0.0,1.0,0.0,
              0.0,0.0,1.0);
}

void Matrix_t::Identity_4x4 (void)
{
   MATRIX_4x4(1.0,0.0,0.0,0.0,
              0.0,1.0,0.0,0.0,
              0.0,0.0,1.0,0.0,
              0.0,0.0,0.0,1.0);
}

void Matrix_t::Rotate_3x3   (const real8 rx, const real8 ry, const real8 rz)
{
   Allocate(3,3);

   Matrix_t mx; mx.Rotate_3x3_X(rx);
   Matrix_t my; my.Rotate_3x3_X(ry);
   Matrix_t mz; mz.Rotate_3x3_X(rz);
   Matrix_t mt;

   mt.Multiply(mx,my);
      Multiply(mt,mz);
}

// (initialize from a 6-component rotation tensor)

void Matrix_t::Rotate_3x3 (const real8 m0, const real8 m1, const real8 m2, const real8 m3, const real8 m4, const real8 m5)
{
   MATRIX_3x3(m0,m5,m4,
              m5,m1,m3,
              m4,m3,m2);
}

void Matrix_t::Rotate_3x3_X (const real8 rx)
{
   MATRIX_3x3(1.0,  0.0  ,   0.0  ,
              0.0,cos(rx),-sin(rx),
              0.0,sin(rx), cos(rx) );
}

void Matrix_t::Rotate_3x3_Y (const real8 ry)
{
   MATRIX_3x3( cos(ry),0.0,sin(ry),
                 0.0  ,1.0,  0.0  ,
              -sin(ry),0.0,cos(ry) );
}

void Matrix_t::Rotate_3x3_Z (const real8 rz)
{
   MATRIX_3x3(cos(rz),-sin(rz), 0.0,
              sin(rz), cos(rz), 0.0,
                0.0  ,   0.0  , 1.0 );
}

void Matrix_t::Rotate_4x4   (const real8 rx, const real8 ry, const real8 rz)
{
   Allocate(4,4);

   Matrix_t mx; mx.Rotate_4x4_X(rx);
   Matrix_t my; my.Rotate_4x4_X(ry);
   Matrix_t mz; mz.Rotate_4x4_X(rz);
   Matrix_t mt;

   mt.Multiply(mx,my);
      Multiply(mt,mz);
}

void Matrix_t::Rotate_4x4_X (const real8 rx)
{
   MATRIX_4x4(1.0,  0.0  ,   0.0  , 0.0,
              0.0,cos(rx),-sin(rx), 0.0,
              0.0,sin(rx), cos(rx), 0.0,
              0.0,  0.0  ,   0.0  , 1.0 );
}

void Matrix_t::Rotate_4x4_Y (const real8 ry)
{
   MATRIX_4x4( cos(ry),0.0,sin(ry), 0.0,
                 0.0  ,1.0,  0.0  , 0.0,
              -sin(ry),0.0,cos(ry), 0.0,
                 0.0  ,0.0,  0.0  , 1.0 );
}

void Matrix_t::Rotate_4x4_Z (const real8 rz)
{
   MATRIX_4x4(cos(rz),-sin(rz), 0.0, 0.0,
              sin(rz), cos(rz), 0.0, 0.0,
                0.0  ,   0.0  , 1.0, 0.0,
                0.0  ,   0.0  , 0.0, 1.0 );
}

void Matrix_t::Scale_2x2 (const real8 sx, const real8 sy)
{
   MATRIX_2x2(sx , 0.0,
              0.0, sy  );
}

void Matrix_t::Scale_3x3 (const real8 sx, const real8 sy, const real8 sz)
{
   MATRIX_3x3(sx , 0.0, 0.0,
              0.0, sy , 0.0,
              0.0, 0.0, sz );
}

void Matrix_t::Scale_4x4 (const real8 sx, const real8 sy, const real8 sz)
{
   MATRIX_4x4(sx , 0.0, 0.0, 0.0,
              0.0, sy , 0.0, 0.0,
              0.0, 0.0, sz , 0.0,
              0.0, 0.0, 0.0, 1.0 );
}

void Matrix_t::Translate_3x3 (const real8 tx, const real8 ty)
{
   MATRIX_3x3(1.0, 0.0, tx,
              0.0, 1.0, ty,
              0.0, 0.0, 1.0);
}

void Matrix_t::Translate_4x4 (const real8 tx, const real8 ty, const real8 tz)
{
   MATRIX_4x4(1.0, 0.0, 0.0, tx,
              0.0, 1.0, 0.0, ty,
              0.0, 0.0, 1.0, tz,
              0.0, 0.0, 0.0, 1.0 );
}

void Matrix_t::Flip_H (void)
{
   int w = mtx_w;
   int h = mtx_h;

   if (mtx && (w>0) && (h>0) )
   {
      for (int j=0; (j<h); ++j)
      {
         real8 *r = mtx+(j*w);

         for (int i=0,k=(w-1); (i<k); ++i,--k)
         { real8 s=r[i]; r[i]=r[k]; r[k]=s; }
      }
   }
}

void Matrix_t::Flip_V (void)
{
   int w = mtx_w;
   int h = mtx_h;

   if (mtx && mtx_rows && (w>0) && (h>0) )
   {
      real8 *tmp   = new real8[w];
      size_t bytes = w * sizeof(real8);

      if (tmp)
      {
         for (int i=0,j=(h-1); (i<j); ++i, --j)
         {
            real8 *p0 = mtx_rows[i];
            real8 *p1 = mtx_rows[j];

            memcpy(tmp,p0 ,bytes);
            memcpy(p0 ,p1 ,bytes);
            memcpy(p1 ,tmp,bytes);
         }

         delete [] tmp;
      }
   }
}

void Matrix_t::Flip_H (const Matrix_t & m)
{
   int w = m.mtx_w;
   int h = m.mtx_h;

   Allocate(w,h);

   if (mtx && m.mtx && (w>0) && (h>0) )
   {
      for (int j=0; (j<h); ++j)
      {
         real8 *p = m.mtx+(j*w);
         real8 *q =   mtx+(j*w);

         for (int i=0,k=(w-1); (i<w); ++i,--k)
            q[i]=p[k];
      }
   }
}

void Matrix_t::Flip_V (const Matrix_t & m)
{
   int w = m.mtx_w;
   int h = m.mtx_h;

   Allocate(w,h);

   if (mtx && mtx_rows && m.mtx && m.mtx_rows && (w>0) && (h>0) )
   {
      for (int i=0,j=(h-1); (i<h); ++i,--j)
         memcpy(mtx_rows[i], m.mtx_rows[j], w*sizeof(real8) );
   }
}

void Matrix_t::Transform (real8 (*f)(const real8 s))
{
   real8 **m = mtx_rows;

   if (m && f)
   {
      for (int j=0; (j<mtx_h); ++j)
         for (int i=0; (i<mtx_w); ++i)
            m[j][i]= (*f)(m[j][i]);
   }
}

void Matrix_t::Transform (real8 (*f)(const int iy, const int ix, const real8 s))
{
   real8 **m = mtx_rows;

   if (m && f)
   {
      for (int j=0; (j<mtx_h); ++j)
         for (int i=0; (i<mtx_w); ++i)
            m[j][i]= (*f)(j,i,m[j][i]);
   }
}

void Matrix_t::Transform (const Matrix_t & mtx, real8 (*f)(const int iy, const int ix, const real8 s))
{
   real8 **m0 =     mtx_rows;
   real8 **m1 = mtx.mtx_rows;

   if (m0 && m1 && f)
   {
      for (int j=0; (j<mtx_h); ++j)
         for (int i=0; (i<mtx_w); ++i)
            m0[j][i]= (*f)(j,i,m1[j][i]);
   }
}

//---------------------------------------------------------------------------------------------------------

void Matrix_t::Real_To_Complex (const Matrix_t & m)
{
   int  w = m.mtx_w;
   int  h = m.mtx_h;

   Allocate(2*w,h);

   if (mtx && mtx_rows && m.mtx && m.mtx_rows && (w>0) && (h>0))
   {
      for (int j=0; (j<h); ++j)
      {
         real8 *p = m.mtx_rows[j];
         real8 *q =   mtx_rows[j];

         for (int i=0,k=0; (i<w); ++i)
         {
            q[k++]=p[i];
            q[k++]=0.0;
         }
      }
   }
}

void Matrix_t::Complex_Abs (const Matrix_t & m)
{
   int  w = m.mtx_w;
   int  h = m.mtx_h;

   Allocate(w/2,h);

   if (mtx && mtx_rows && m.mtx && m.mtx_rows && (w>0) && (h>0))
   {
      for (int j=0; (j<h); ++j)
      {
         real8 *p = m.mtx_rows[j];
         real8 *q =   mtx_rows[j];

         for (int i=0,k=0; (i<(w/2)); ++i)
         {
            real8 a = p[k++];
            real8 b = p[k++];

            q[i] = sqrt(a*a+b*b);
         }
      }
   }
}

// FFT_Transform()
//
// Loads a complex matrix from the source matrix provided and performs a two-dimensional
// discrete fourier transform (DFT) using FFT's.
//---------------------------------------------------------------------------------------------------------

void Matrix_t::FFT_Transform (const Matrix_t & m, const int dir, const int norm)
{
   int  w = m.mtx_w;
   int  h = m.mtx_h;
   int  n = m.mtx_w;

   if (this==&m) { printf("%s::%s(ln=%d) - cannot perform FFT transform in place\n"              , __FILE__, __func__, __LINE__); return; }
   if (w!=h)     { printf("%s::%s(ln=%d) - cannot perform FFT transform on non-square matrices\n", __FILE__, __func__, __LINE__); return; }

   Allocate(2*n,n);

   int   *bfly = FFT_Butterfly(n);
   real8 *col  = new real8[2*n];

   if (mtx && mtx_rows && m.mtx && m.mtx_rows && (n>0) && bfly && col)
   {
      // perform row-wise transforms...

      for (int i=0; (i<n); ++i)
      {
         real8 *p = m.mtx_rows[i];
         real8 *q =   mtx_rows[i];

         ::FFT_Load     (q,p,bfly,n);
         ::FFT_Transform(q,n,dir,norm);
      }

      // perform column-wise transforms...

      for (int i=0; (i<n); ++i)
      {
         real8 *p = mtx+(2*i);
         real8 *q = col;

         ::FFT_Load_Complex (q,p,bfly,n,2*n);
         ::FFT_Transform    (q,n,dir,norm);

         for (int j=0,k=0,l=0; (j<n); ++j)
         {
            p[l  ] = q[k++];
            p[l+1] = q[k++]; l+=(2*n);
         }
      }
   }

   if (bfly) { delete [] bfly; bfly=0; }
   if (col ) { delete [] col ; col =0; }
}

void Matrix_t::FFT_Power_Spectrum (const Matrix_t & m)
{
   int  w = m.mtx_w;
   int  h = m.mtx_h;

   if (this==&m) { printf("%s::%s(ln=%d) - cannot perform FFT power spectrum in place\n", __FILE__, __func__, __LINE__); return; }

   Allocate(w/2,h);

   if (mtx && mtx_rows && m.mtx && m.mtx_rows && (w>0) && (h>0))
   {
      for (int i=0; (i<h); ++i)
      {
         real8 *p = m.mtx_rows[i];
         real8 *q =   mtx_rows[i];

         ::FFT_Power_Spectrum(q,p,w);
      }
   }
}

void Matrix_t::Extract (const Matrix_t & m, const int x0, const int y0, const int w, const int h)
{
   Allocate(w,h);

   real8 *p = (    (0<=x0) && (x0<m.mtx_w)
                && (0<=y0) && (y0<m.mtx_h) ? m.mtx+(y0*m.mtx_w)+x0 : 0 );

   real8 *q = ( mtx );

   if (p && q)
   {
      size_t bytes = mtx_w * sizeof(real8);
      int    pw    = m.mtx_w;
      int    qw    =   mtx_w;

      for (int i=0; (i<mtx_h); ++i, p+=pw, q+=qw)
         memcpy(q,p,bytes);
   }
}

void Matrix_t::Overlay (const Matrix_t & m, const int x0, const int y0)
{
   real8 *p = ( m.mtx );

   real8 *q = (    (0<=x0) && (x0<mtx_w)
                && (0<=y0) && (y0<mtx_h) ? mtx+(y0*mtx_w)+x0 : 0 );

   if (p && q)
   {
      size_t bytes = m.mtx_w * sizeof(real8);
      int    pw    = m.mtx_w;
      int    qw    =   mtx_w;

      for (int i=0; (i<m.mtx_h); ++i, p+=pw, q+=qw)
         memcpy(q,p,bytes);
   }
}

//-----------------------------------------------------------------------------------------------

void Matrix_t::Print (FILE *fd, const char *xfmt ) const
{
   const char *fmt = ( xfmt ? xfmt : "%8.4lf " );

   if (mtx && mtx_rows)
   {
      for (int j=0; (j<mtx_h); ++j)
      {
         for (int i=0; (i<mtx_w); ++i)
            fprintf(fd,fmt,mtx_rows[j][i]);

         fprintf(fd,"\n");
      }
   }
}

void Matrix_t::Save (const char *path) const
{
   FILE *fd = ( path ? fopen(path,"w") : 0 );

   if (fd)
   {
      if (mtx)
      {
         fwrite(&mtx_w, sizeof(mtx_w), 1, fd);
         fwrite(&mtx_h, sizeof(mtx_h), 1, fd);
         fwrite( mtx, mtx_w * mtx_h * sizeof(real8), 1, fd);
      }

      fclose(fd);
   }
}

void Matrix_t::Load (const char *path)
{
   FILE *fd = ( path ? fopen(path,"r") : 0 );

   if (fd)
   {
      int w; fread(&w, sizeof(w), 1, fd);
      int h; fread(&h, sizeof(h), 1, fd);

      Allocate(w,h);

      if (mtx)
        fread(mtx, w*h*sizeof(real8), 1, fd);

      fclose(fd);
   }
}

// MPI communication support....
//-----------------------------------------------------------------------------------------------

int    Matrix_t::MPI_Count (void) { return( mtx && (mtx_w>0) && (mtx_h>0) ? (mtx_w*mtx_h)+2 : 0 ); }
size_t Matrix_t::MPI_Bytes (void) { return( MPI_Count()*sizeof(real8) ); }

real8 *Matrix_t::MPI_Pack (real8 *q)
{
   if (!q) { q = ( mtx && (mtx_w>0) && (mtx_h>0) ? new real8[mtx_w*mtx_h+2] : 0 ); }

   if (q)
   {
      *q++ = (real8) mtx_w;
      *q++ = (real8) mtx_h;

      if (mtx && (mtx_w>0) && (mtx_h>0))
      {
         memcpy(q,mtx,(mtx_w*mtx_h)*sizeof(real8));

         q+=(mtx_w*mtx_h);
      }
   }

   return(q);
}

real8 *Matrix_t::MPI_Unpack (real8 *p)
{
   if (p)
   {
      int w = (int) *p++;
      int h = (int) *p++;

      Allocate(w,h);

      if (mtx && (mtx_w>0) && (mtx_h>0))
      {
         memcpy(mtx,p,(mtx_w*mtx_h)*sizeof(real8));

         p+=(mtx_w*mtx_h);
      }
   }

   return(p);
}

#ifdef PARALLEL
void Matrix_t::MPI_Send (const int dst, const int mtag, MPI_Comm comm, const int sync)
{
   if (mtx_mpi) { delete [] mtx_mpi; mtx_mpi=0; mtx_mpi_cnt=0; }

   mtx_mpi_cnt = MPI_Count();
   mtx_mpi     = MPI_Pack ();

   if (mtx_mpi && (mtx_mpi_cnt>0) )
   {
      if (sync) ::MPI_Send ( mtx_mpi, mtx_mpi_cnt, MPI_DOUBLE, dst, mtag, comm );
      else      ::MPI_Isend( mtx_mpi, mtx_mpi_cnt, MPI_DOUBLE, dst, mtag, comm, &mtx_mpi_request);
   }
}

void Matrix_t::MPI_Receive (const int src, const int mtag, MPI_Comm comm, const int sync)
{
   if (mtx_mpi)
   {
      if (sync) ::MPI_Recv ( mtx_mpi, mtx_mpi_cnt, MPI_DOUBLE, src, mtag, comm, &mtx_mpi_status  );
      else      ::MPI_Irecv( mtx_mpi, mtx_mpi_cnt, MPI_DOUBLE, src, mtag, comm, &mtx_mpi_request );
   }
}

#endif

// LUP_Decompose()
//
// Performs an in-place, LU decomposition of a matrix.
//
// On exit, the source matrix is deomposed into two matrices :  L-E and U as
// A=(L-E)+U such that P*A=L*U.  Each element of of the pivot vector contains the
// column index where the permutation matrix has "1". The last element of the pivot
// vector, P[N]=S+N, where S is the number of row exchanges needed for determinant
// computation, det(P)=(-1)^S
//-----------------------------------------------------------------------------------------------

int LUP_Decompose
(
   real8 **a  ,  /// source matrix (NxN)
   int    *p  ,  /// pivot vector  (N+1)
   int     n  ,  /// size of the source matrix (NxN)
   real8   tol   /// error tolerance for degeneracy test
)
{
   if (a && p && (n>0) )
   {
      for(int i=0; (i<=n); ++i) { p[i]=i; }       // initialize the pivot vector

      for(int i=0; (i< n); ++i)
      {
         int   imax = i;
         real8 maxa = 0.0;
         real8 absa = 0.0;

         for(int k=i; (k<n); ++k)
            if((absa=fabs(a[k][i]))>maxa){ maxa=absa; imax=k; }

         if (maxa<tol) return(-1); //failure, matrix is degenerate

         if(imax!=i)
         {
            { int    tmp=p[i]; p[i]=p[imax]; p[imax]=tmp; } // (swap pivot)
            { real8 *tmp=a[i]; a[i]=a[imax]; a[imax]=tmp; } // (swap rows )

            p[n]++; //increment pivot count (for determinant)
         }

         for(int j=(i+1); (j<n); ++j)
         {
            a[j][i]/=a[i][i];

            for(int k=(i+1); (k<n); ++k)
               a[j][k] -= a[j][i]*a[i][k];
         }
      }
   }

   return(0);
}

// LUP_Solve()
//
// Will solve the system Ax=b using a matrix that was decomposed using LUP_Decompose().
// Inputs are the source decomposed matrix, pivot vector, and the right-hand-side (RHS)
// vector (b). The result is left in x.
//-----------------------------------------------------------------------------------------------

int LUP_Solve
(
   real8 **a,  ///< source matrix      (NxN, LUP decomposed)
   int    *p,  ///< pivot vector       (N)
   real8  *x,  ///< computed solution  (N)
   real8  *b,  ///< rhs vector         (N)
   int     n   ///< size of the system (NxN)
)
{
   if (a && p && x && b && (n>0))
   {
      for(int i=0; (i<n); ++i)
      {
         x[i] = b[p[i]];

         for(int k=0; (k<i); ++k)
            x[i] -= a[i][k]*x[k];
      }

      for(int i=(n-1); (i>=0); i--)
      {
         for(int k=(i+1); (k<n); ++k)
            x[i] -= a[i][k]*x[k];

         x[i] = x[i]/a[i][i];
      }
   }

   return(0);
}

// LUP_Invert()
//
// Will compute the inverse of a matrix that was decomposed using LUP_Decompose().
// Inputs are the source decomposed matrix and pivot vector. Note that only the
// first N elements of the pivot vector are used. Also note that the source
// and destination matrices must be different.
//-----------------------------------------------------------------------------------------------

void LUP_Invert
(
   real8 **ai,   ///< inverse matrix     (NxN, result)
   real8 **a ,   ///< source matrix      (NxN, LUP decomposed)
   int    *p ,   ///< pivot vector       (N)
   int     n     ///< size of the matrix (NxN)
)
{
   if (ai && a && (ai!=a) && p && (n>0))
   {
      for(int j=0; (j<n); ++j)
      {
         for(int i=0; (i<n); ++i)
         {
            ai[i][j] = ( (p[i]==j) ? 1.0 : 0.0 );

            for(int k=0; (k<i); ++k)
               ai[i][j] -= a[i][k]*ai[k][j];
         }

         for(int i=(n-1); (i>=0); i--)
         {
            for(int k=(i+1); (k<n); ++k)
               ai[i][j] -= a[i][k]*ai[k][j];

            ai[i][j] = ai[i][j]/a[i][i];
         }
      }
   }
}

// LUP_Determinant()
//
// Will return the determinant of a matrix that was decomposed using LUP_Decompose().
// Inputs are the decomposed matrix and pivot vector.  Note that the pivot vector is of
// length N+1 (not N) where the final entry of the pivot vector contains the number of
// row exchanges that occurred in the decomposition.
//-----------------------------------------------------------------------------------------------

real8 LUP_Determinant
(
   real8 **a,  ///< source matrix      (NxN, LUP decomposed)
   int    *p,  ///< pivot vector       (N+1)
   int     n   ///< size of the matrix (N)
)
{
   real8 det=0.0;

   if (a && p && (n>0))
   {
      det = a[0][0];

      for(int i=1; (i<n); ++i)
         det *= a[i][i];

      det = ( ((p[n]-n)%2==0) ? det : -det );
   }

   return(det);
}

//---------------------------------------------------------------------------------------------------------------------------------

Matrix_t   M2x2_Identity   (void) { Matrix_t m; m.Identity_2x2(); return(m); }
Matrix_t   M3x3_Identity   (void) { Matrix_t m; m.Identity_3x3(); return(m); }
Matrix_t   M4x4_Identity   (void) { Matrix_t m; m.Identity_4x4(); return(m); }

Matrix_t   M3x3_Rotate     (const real8 rx, const real8 ry, const real8 rz) { Matrix_t m; m.Rotate_3x3  (rx,ry,rz); return(m); }
Matrix_t   M3x3_Rotate_X   (const real8 rx)                                 { Matrix_t m; m.Rotate_3x3_X(rx);       return(m); }
Matrix_t   M3x3_Rotate_Y   (const real8 ry)                                 { Matrix_t m; m.Rotate_3x3_Y(ry);       return(m); }
Matrix_t   M3x3_Rotate_Z   (const real8 rz)                                 { Matrix_t m; m.Rotate_3x3_Z(rz);       return(m); }

Matrix_t   M4x4_Rotate     (const real8 rx, const real8 ry, const real8 rz) { Matrix_t m; m.Rotate_4x4  (rx,ry,rz); return(m); }
Matrix_t   M4x4_Rotate_X   (const real8 rx)                                 { Matrix_t m; m.Rotate_4x4_X(rx);       return(m); }
Matrix_t   M4x4_Rotate_Y   (const real8 ry)                                 { Matrix_t m; m.Rotate_4x4_Y(ry);       return(m); }
Matrix_t   M4x4_Rotate_Z   (const real8 rz)                                 { Matrix_t m; m.Rotate_4x4_Z(rz);       return(m); }

Matrix_t   M2x2_Scale      (const real8 sx, const real8 sy)                 { Matrix_t m; m.Scale_2x2     (sx,sy);    return(m); }
Matrix_t   M3x3_Scale      (const real8 sx, const real8 sy, const real8 sz) { Matrix_t m; m.Scale_3x3     (sx,sy,sz); return(m); }
Matrix_t   M4x4_Scale      (const real8 sx, const real8 sy, const real8 sz) { Matrix_t m; m.Scale_4x4     (sx,sy,sz); return(m); }
Matrix_t   M3x3_Translate  (const real8 tx, const real8 ty)                 { Matrix_t m; m.Translate_3x3 (tx,ty   ); return(m); }
Matrix_t   M4x4_Translate  (const real8 tx, const real8 ty, const real8 tz) { Matrix_t m; m.Translate_4x4 (tx,ty,tz); return(m); }

//---------------------------------------------------------------------------------------------------------------------------------

void       M2x2_Identity   (Matrix_t & m) { m.Identity_2x2(); }
void       M3x3_Identity   (Matrix_t & m) { m.Identity_3x3(); }
void       M4x4_Identity   (Matrix_t & m) { m.Identity_4x4(); }

void       M3x3_Rotate     (Matrix_t & m, const real8 rx, const real8 ry, const real8 rz) { m.Rotate_3x3  (rx,ry,rz); }
void       M3x3_Rotate_X   (Matrix_t & m, const real8 rx)                                 { m.Rotate_3x3_X(rx);       }
void       M3x3_Rotate_Y   (Matrix_t & m, const real8 ry)                                 { m.Rotate_3x3_Y(ry);       }
void       M3x3_Rotate_Z   (Matrix_t & m, const real8 rz)                                 { m.Rotate_3x3_Z(rz);       }

void       M4x4_Rotate     (Matrix_t & m, const real8 rx, const real8 ry, const real8 rz) { m.Rotate_4x4  (rx,ry,rz); }
void       M4x4_Rotate_X   (Matrix_t & m, const real8 rx)                                 { m.Rotate_4x4_X(rx);       }
void       M4x4_Rotate_Y   (Matrix_t & m, const real8 ry)                                 { m.Rotate_4x4_Y(ry);       }
void       M4x4_Rotate_Z   (Matrix_t & m, const real8 rz)                                 { m.Rotate_4x4_Z(rz);       }

void       M3x3_Scale      (Matrix_t & m, const real8 sx, const real8 sy, const real8 sz) { m.Scale_3x3     (sx,sy,sz); }
void       M4x4_Scale      (Matrix_t & m, const real8 sx, const real8 sy, const real8 sz) { m.Scale_4x4     (sx,sy,sz); }
void       M3x3_Translate  (Matrix_t & m, const real8 tx, const real8 ty)                 { m.Translate_3x3 (tx,ty   ); }
void       M4x4_Translate  (Matrix_t & m, const real8 tx, const real8 ty, const real8 tz) { m.Translate_4x4 (tx,ty,tz); }

//---------------------------------------------------------------------------------------------------------------------------------
// SVD Support...
//---------------------------------------------------------------------------------------------------------------------------------

static inline int    Min  (const int     a, const int    b ) { return( (a<b) ? a : b ); }
static inline real8  Max  (const real8   a, const real8  b ) { return( (a>b) ? a : b ); }
static inline real8  Sign (const real8   a, const real8  b ) { return( (b>=0.0)? fabs(a) : -fabs(a) ); }

// Pythag()
//
// Pythagorean function implementation.
//
// This function will return the result of sqrt(x*x + y*y)
// without the sensitivities to scale...
//
// Based on routine from Numerical Recipes in C.
//--------------------------------------------------------------

static real8  Pythag(const real8 a, const real8 b)
{
	real8 absa = fabs(a);
	real8 absb = fabs(b);

	     if (absa>absb) { real8 x = (absb/absa); return( absa*sqrt(1.0+(x*x)) ); }
	else if (absb>0.0 ) { real8 x = (absa/absb); return( absb*sqrt(1.0+(x*x)) ); }

    return(0.0);
}

// MatrixMN_SVD()
//
// Computes the Singular Value Decomposition of the source matrix.
// Upon exit, populates three matrix and vector references that
// satisfy
//    SVD(A) = U*W*V'
//
// On entry:
//   A : MxN zero-indexed Iliffe matrix to be decomposed (input,updated)
//   V : NxN zero-indexed Iliffe matrix
//   W : 1xN zero-indexed vector
//   R : 1xN zero-indexed scratch vector
//
//   The source matrix A is provided as an Iliffe matrix, where
//   A[i] points to each row of the source matrix.  The source matrix
//   is passed by reference and modified by this routine.
//
// On successful exit:
//   a : MxN zero-indexed Iliffe matrix contains the decomposed result
//   v : NxN zero-indexed Iliffe matrix containing the transpose of V (e.g V')
//   w : 1xN vector containing the diagonal elements of the W matrix
//
// Returns...
//   0 : (success)
//  -1 : (failure - iteration count exceeded).
//
// Adapted from Numerical Recipes in C, Second Edition
//--------------------------------------------------------------------------------------

static int MxN_SVD
(
          real8 **a,  ///< MxN zero-indexed Iliffe matrix to be decomposed (input,updated)
          real8 **v,  ///< NxN zero-indexed Iliffe matrix
          real8  *w,  ///< 1xN allocated vector
          real8  *r,  ///< 1xN allocated scratch vector
   const  int     m,  ///< size of the source matrix (M-rows)
   const  int     n   ///< size of the source matrix (N-columns)
)
{
   int   i,its,j,jj,k,l,nm;
   real8 anorm,c,f,g,h,s,scale,x,y,z,*rv1;

   g=scale=anorm=0.0;

   for (i=0; (i<n); i++)
   {
      l=i+1;
      r[i]=scale*g;
      g=s=scale=0.0;
      if (i<m)
      {
         for (k=i; (k<m); k++) { scale += fabs(a[k][i]); }

         if (scale)
         {
            for (k=i; (k<m); k++)
            {
               a[k][i] /= scale;
               s += a[k][i]*a[k][i];
            }

            f=a[i][i];
            g = -Sign(sqrt(s),f);
            h=f*g-s;
            a[i][i]=f-g;

            for (j=l; (j<n); j++)
            {
               for (s=0.0,k=i; (k<m); k++) { s += a[k][i]*a[k][j]; }
               f=s/h;
               for (      k=i; (k<m); k++) { a[k][j] += f*a[k][i]; }
            }
            for (k=i; (k<m); k++) { a[k][i] *= scale; }
         }
      }

      w[i]=scale*g;

      g=s=scale=0.0;

      if ( (i<m) && i != n-1)
      {
         for (k=l; (k<n); k++) scale += fabs(a[i][k]);
         if (scale)
         {
            for (k=l; (k<n); k++)
            {
               a[i][k] /= scale;
               s += a[i][k]*a[i][k];
            }
            f=a[i][l];
            g = -Sign(sqrt(s),f);
            h=f*g-s;
            a[i][l]=f-g;
            for (k=l; (k<n); k++) { r[k]=a[i][k]/h; }
            for (j=l; (j<m); j++)
            {
               for (s=0.0,k=l; k<n; k++) { s += a[j][k]*a[i][k]; }
               for (      k=l; k<n; k++) { a[j][k] += s*r[k];  }
            }
            for (k=l; k<n; k++) { a[i][k] *= scale; }
         }
      }
      anorm = Max( anorm, (fabs(w[i])+fabs(r[i])) );
   }

   for (i=n-1; (i>=0); i--)
   {
      if (i<(n-1))
      {
         if (g)
         {
            for (j=l; (j<n); j++) { v[j][i]=(a[i][j]/a[i][l])/g; }
            for (j=l; (j<n); j++)
            {
               for (s=0.0,k=l; (k<n); k++) { s += a[i][k]*v[k][j]; }
               for (      k=l; (k<n); k++) { v[k][j] += s*v[k][i]; }
            }
         }
         for (j=l; (j<n); j++) { v[i][j]=v[j][i]=0.0; }
      }
      v[i][i]=1.0;
      g=r[i];
      l=i;
   }

   for (i=Min(m,n)-1; (i>=0); i--)
   {
      l=i+1;
      g=w[i];

      for (j=l; (j<n); j++) { a[i][j]=0.0; }

      if (g)
      {
         g=1.0/g;

         for (j=l; (j<n); j++)
         {
            for (s=0.0,k=l; (k<m); k++) { s += a[k][i]*a[k][j]; }
            f=(s/a[i][i])*g;
            for (      k=i; (k<m); k++) { a[k][j] += f*a[k][i]; }
         }
         for (j=i; (j<m); j++) { a[j][i] *= g; }
      }
      else
         for (j=i; (j<m); j++) { a[j][i] = 0.0; }

      ++a[i][i];
   }

   for (k=(n-1); (k>=0); k--)
   {
      for (its=0; (its<=30); its++)
      {
         int flag=1;

         for (l=k; (l>=0); l--)
         {
            nm=l-1;
            if ( (fabs(r[l] )+anorm) == anorm) { flag=0; break; }
            if ( (fabs(w[nm])+anorm) == anorm) {         break; }
         }

         if (flag)
         {
            c=0.0;
            s=1.0;

            for (i=l; (i<k); i++)
            {
               f=s*r[i];
               r[i]=c*r[i];
               if ( (fabs(f)+anorm) == anorm) break;
               g=w[i];
               h=Pythag(f,g);
               w[i]=h;
               h=1.0/h;
               c=g*h;
               s = -f*h;
               for (j=0; (j<m); j++)
               {
                  y=a[j][nm];
                  z=a[j][i ];
                  a[j][nm]=y*c+z*s;
                  a[j][i ]=z*c-y*s;
               }
            }
         }

         z=w[k];

         if (l==k)
         {
            if (z<0.0)
            {
               w[k] = -z;
               for (j=0; (j<n); j++) { v[j][k] = -v[j][k]; }
            }
            break;
         }

         if (its>=30)
         {
            printf("%s::%s(%d) - maximum iteration count reached, solution invalid\n", __FILE__, __func__, __LINE__ );
            return(-1);
         }

         x=w[l];
         nm=k-1;
         y=w[nm];
         g=r[nm];
         h=r[k];
         f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
         g=Pythag(f,1.0);
         f=((x-z)*(x+z)+h*((y/(f+Sign(g,f)))-h))/x;
         c=s=1.0;

         for (j=l; (j<=nm); j++)
         {
            i=j+1;
            g=r[i];
            y=w[i];
            h=s*g;
            g=c*g;
            z=Pythag(f,h);
            r[j]=z;
            c=f/z;
            s=h/z;
            f=x*c+g*s;
            g = g*c-x*s;
            h=y*s;
            y *= c;

            for (jj=0; (jj<n); jj++)
            {
               x=v[jj][j];
               z=v[jj][i];
               v[jj][j]=x*c+z*s;
               v[jj][i]=z*c-x*s;
            }
            z=Pythag(f,h);
            w[j]=z;

            if (z)
            {
               z=1.0/z;
               c=f*z;
               s=h*z;
            }
            f=c*g+s*y;
            x=c*y-s*g;

            for (jj=0; (jj<m); jj++)
            {
               y=a[jj][j];
               z=a[jj][i];
               a[jj][j]=y*c+z*s;
               a[jj][i]=z*c-y*s;
            }
         }

         r[l]=0.0;
         r[k]=f;
         w[k]=x;
      }
   }

   return(0);
}

// SVD()
//
// Computes the Singular Value Decomposition of the source matrix.
// Upon exit, populates three matrix and vector references that
// satisfy
//    A = U*W*V'
//
// For the most part, I haven't modified the code beyond altering the
// interface to accommodate zero-indexed rather than 1-indexed matrices.
// I've also added the U matrix as a pass by reference. This enables the
// method to leave the original source matrix unmodified.
//
// On entry:
//   The source matrix A (the source instance)
//
// On successful exit:
//   umat - contains the U matrix (0..m-1,0...n-1);
//   wvec - contains the elements of the diagonals of the W matrix (0...n-1)
//   vmat - contains the transpose of the V matrix (0..n-1,0...n-1);
//
// Other notes...
//
// Returns zero on successful completion, otherwise posts an
// error message using the Set_Error() method and returns the following
// error codes below.
//
//   0 : (success)
//  -1 : iteration count exceeded.
//
// Adapted from Numerical Recipes in C, Second Edition
/////////////////////////////////////////////////////////////////////

int Matrix_t::SVD(Matrix_t & umat, Vector_t & wvec, Matrix_t & vmat)
{
   int    flag,i,its,j,jj,k,l,nm;
   real8  anorm,c,f,g,h,s,scale,x,y,z;

   umat = *this; // (leave the original untouched)

   int m = mtx_h;
   int n = mtx_w;

   Matrix_t wmat; wmat.Allocate(n,n);  wmat=0.0;
                  vmat.Allocate(n,n);  vmat=0.0;

                  wvec.Allocate(n);    wvec=0.0;
   Vector_t rvec; rvec.Allocate(n);    rvec=0.0;

   real8 **a = umat.mtx_rows;
   real8 **v = vmat.mtx_rows;
   real8  *r = rvec.vec;
   real8  *w = wvec.vec;

   int     err = MxN_SVD(a,v,w,r,m,n);

   wmat = 0.0;
   wmat.Set_Diagonal(wvec);

   return (err);
}

// SVD_Inverse()
//
// Will develop a pseudo-inverse using the U-V-W matrices computed using
// singular value decompostion.
//
// Essentially, performs the logic...
//
//       [A]  = [U]*[W]*[V]'
//   Inv([A]) = [V]*[diag(1/w)][U]'
//
// On entry:
//   The source instance contains the matrix to be inverted
//
// On successful exit:
//   m - contains the pseudo inverse of the source instance
//
// Error exit:
//   simply passes on the interation error from the SVD call
//-------------------------------------------------------------------------------

int Matrix_t::SVD_Inverse (Matrix_t & m)
{
   *this = m;

   Matrix_t  u, v;
   Vector_t  w;

   if ( SVD(u,w,v) < 0 ) return (-1);                           // u,w,v = svd(m)

   w.Reciprocol();                                              // w = diag(1/w)
   u.Transpose();                                               // u = u'

   *this = v;                                                   //  a = v
                                                                //
   for (int i=0,k=0; (i<mtx_h); i++)                            //
   for (int j=0;     (j<mtx_w); j++,k++) { mtx[k] *= w[j]; }    //  a = v x diag(1/w)
                                                                //
   *this *= u;                                                  //  a = v x diag(1/w) x u'

   return (0);
}

