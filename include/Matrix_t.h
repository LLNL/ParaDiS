#pragma once

#ifndef __PDS_MATRIX_MXN_H
#define __PDS_MATRIX_MXN_H

#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Typedefs.h"

class Vector_t;

class Matrix_t
{
   public :
      int          mtx_w          ;  ///< matrix width
      int          mtx_h          ;  ///< matrix height

      real8       *mtx            ;  ///< contiguous memory buffer allocated for this matrix
      real8      **mtx_rows       ;  ///< pointers to each row (allows for access via [i][j] operators
      int         *mtx_piv        ;  ///< lu decomposition pivot array

#ifdef PARALLEL
      int          mtx_mpi_cnt    ;  ///< count of packed MPI elements
      size_t       mtx_mpi_bytes  ;  ///< bytes of packed MPI elements
      real8       *mtx_mpi        ;  ///< buffer used for MPI communications
      MPI_Status   mtx_mpi_status ;  ///< most recent MPI send/receive status
      MPI_Request  mtx_mpi_request;  ///< most recent MPI send/receive request
#endif

   public :
       Matrix_t(void);
       Matrix_t(const int w, const int h);
       Matrix_t(const int w, const int h, const          int     s);
       Matrix_t(const int w, const int h, const          real8   s);
       Matrix_t(const int w, const int h, const unsigned char   *s);
       Matrix_t(const int w, const int h, const unsigned short  *s);
       Matrix_t(const int w, const int h, const          short  *s);
       Matrix_t(const int w, const int h, const unsigned int    *s);
       Matrix_t(const int w, const int h, const          int    *s);
       Matrix_t(const int w, const int h, const          real8  *s);
       Matrix_t(const int w, const int h, real8 (*fn)(const int i, const int j));
       Matrix_t(const Matrix_t & s);

       Matrix_t(const real8 m00, const real8 m01,
                const real8 m10, const real8 m11 );

       Matrix_t(const real8 m00, const real8 m01, const real8 m02,
                const real8 m10, const real8 m11, const real8 m12 );

       Matrix_t(const real8 m00, const real8 m01, const real8 m02,
                const real8 m10, const real8 m11, const real8 m12,
                const real8 m20, const real8 m21, const real8 m22 );

       Matrix_t(const real8 m00, const real8 m01, const real8 m02, const real8 m03,
                const real8 m10, const real8 m11, const real8 m12, const real8 m13,
                const real8 m20, const real8 m21, const real8 m22, const real8 m23 );

       Matrix_t(const real8 m00, const real8 m01, const real8 m02, const real8 m03,
                const real8 m10, const real8 m11, const real8 m12, const real8 m13,
                const real8 m20, const real8 m21, const real8 m22, const real8 m23,
                const real8 m30, const real8 m31, const real8 m32, const real8 m33 );

       Matrix_t(const real8 m00, const real8 m01, const real8 m02, const real8 m03, const real8 m04, const real8 m05, const real8 m06, const real8 m07,
                const real8 m10, const real8 m11, const real8 m12, const real8 m13, const real8 m14, const real8 m15, const real8 m16, const real8 m17,
                const real8 m20, const real8 m21, const real8 m22, const real8 m23, const real8 m24, const real8 m25, const real8 m26, const real8 m27,
                const real8 m30, const real8 m31, const real8 m32, const real8 m33, const real8 m34, const real8 m35, const real8 m36, const real8 m37,
                const real8 m40, const real8 m41, const real8 m42, const real8 m43, const real8 m44, const real8 m45, const real8 m46, const real8 m47,
                const real8 m50, const real8 m51, const real8 m52, const real8 m53, const real8 m54, const real8 m55, const real8 m56, const real8 m57,
                const real8 m60, const real8 m61, const real8 m62, const real8 m63, const real8 m64, const real8 m65, const real8 m66, const real8 m67,
                const real8 m70, const real8 m71, const real8 m72, const real8 m73, const real8 m74, const real8 m75, const real8 m76, const real8 m77 );

      ~Matrix_t();

       void  Recycle                (void);
       void  Allocate               (const int w, const int h);

       real8 **Offset_Matrix        (const int ioff, const int joff);

       const Matrix_t & operator =  (const Matrix_t       & m);
       const Matrix_t & operator =  (const real8            s);

       const Matrix_t & operator =  (const unsigned char   *s);
       const Matrix_t & operator =  (const unsigned short  *s);
       const Matrix_t & operator =  (const          short  *s);
       const Matrix_t & operator =  (const unsigned int    *s);
       const Matrix_t & operator =  (const          int    *s);
       const Matrix_t & operator =  (const          real8  *s);

       int      operator == (const Matrix_t  & m) const;
       int      operator == (const real8       s) const;
       int      operator == (const real8      *s) const;
       int      operator != (const Matrix_t  & m) const;
       int      operator != (const real8       s) const;
       int      operator != (const real8      *s) const;
       int      operator <  (const Matrix_t  & m) const;
       int      operator <  (const real8       s) const;
       int      operator <  (const real8      *s) const;
       int      operator <= (const Matrix_t  & m) const;
       int      operator <= (const real8       s) const;
       int      operator <= (const real8      *s) const;
       int      operator >  (const Matrix_t  & m) const;
       int      operator >  (const real8       s) const;
       int      operator >  (const real8      *s) const;
       int      operator >= (const Matrix_t  & m) const;
       int      operator >= (const real8       s) const;
       int      operator >= (const real8      *s) const;

       Matrix_t operator +  (const Matrix_t &  m);
       Matrix_t operator +  (const real8       s);

       Matrix_t operator -  (const Matrix_t &  m);
       Matrix_t operator -  (const real8       s);

       Matrix_t operator *  (const Matrix_t &  m);
       Matrix_t operator *  (const real8       s);

       Matrix_t operator /  (const real8       s);

       void     operator += (const Matrix_t  & m);
       void     operator += (const real8       s);
       void     operator += (const real8      *s);

       void     operator -= (const Matrix_t  & m);
       void     operator -= (const real8       s);
       void     operator -= (const real8      *s);

       void     operator *= (const Matrix_t  & m);
       void     operator *= (const real8       s);
       void     operator *= (const real8      *s);

       void     operator /= (const real8       s);

       Matrix_t operator -  ();

       operator  real8  *() const { return ( (real8  *) mtx     ); }
       operator  real8 **() const { return ( (real8 **) mtx_rows); }

       real8 *Get_Row         (const int j) const;
       real8 *Get_Column      (const int i) const;
       void   Swap_Rows       (const int i, const int j);
       void   Swap_Columns    (const int i, const int j);

       real8 *Get_Diagonal    (void) const;
       void   Set_Diagonal    (const real8 *diag);

       real8  Row_Min         (const int i, int & indx);
       real8  Row_Max         (const int i, int & indx);
       real8  Column_Min      (const int i, int & indx);
       real8  Column_Max      (const int i, int & indx);

       void   Sum_Rows        (const int i, const int j);
       void   Sum_Rows        (const int i, const int j, const real8 s);

       void   Add             (const Matrix_t & m );
       void   Add             (const Matrix_t & ma, const Matrix_t & mb);

       void   Subtract        (const Matrix_t & m );
       void   Subtract        (const Matrix_t & ma, const Matrix_t & mb);

       void   Multiply        (const Matrix_t & m );
       void   Multiply        (const Matrix_t & ma, const Matrix_t & mb);

       void   Linear_Sum      (const Matrix_t & a, const Matrix_t & x,
                               const Matrix_t & b, const Matrix_t & y);

       real8  Dot             (const Matrix_t & m );

       void   Tensor          (const Matrix_t & ma, const Matrix_t & mb);

       void   Transpose       (void);
       void   Transpose       (const Matrix_t & m );

       void   Inverse         (void);
       void   Inverse         (const Matrix_t & m );

       int    Compare         (const Matrix_t & m, const real8  eps=4.0e-8) const;
       int    Compare         (const real8      s, const real8  eps=4.0e-8) const;
       int    Compare         (const real8     *s, const real8  eps=4.0e-8) const;

       void   Identity        (void);
       void   Identity        (const int w, const int h);

       void   Pascal          (void);
       void   Pascal          (const int w, const int h);

       void   Hilbert         (void);
       void   Hilbert         (const int w, const int h);

       void   Walsh           (const int n);

       void   DCT             (void);
       void   DCT             (const int w, const int h);

       void   Rand            (                          const real8 min=0.0, const real8 max=1.0);
       void   Rand            (const int w, const int h, const real8 min=0.0, const real8 max=1.0);

       void   Abs             (void);
       void   Rint            (void);

       real8  Min             (int & iy, int & ij) const;
       real8  Max             (int & iy, int & ij) const;
       real8  Min             (void) const;
       real8  Max             (void) const;
       real8  Mean            (void) const;
       real8  Range           (void) const;
       void   Moment          (real8 & min, real8 & max, real8 & mean, real8 & range              ) const;
       void   Moment          (real8 & min, real8 & max, real8 & mean, real8 & range, real8 & sdev) const;

       real8  Abs_Min         (int & iy, int & ij) const;
       real8  Abs_Max         (int & iy, int & ij) const;
       real8  Abs_Min         (void) const;
       real8  Abs_Max         (void) const;

       void   Sign            (void);

       real8  Norm            (void);
       void   Normalize       (void);
       real8  Max_Norm        (void);

       real8  Condition       (void);

       void   LUP_Decompose   (real8 tol);
       void   LUP_Solve       (real8 *x, real8  *b);
       void   LUP_Invert      (real8 tol);
       real8  LUP_Determinant (void);

       int    SVD             (Matrix_t & umat, Vector_t & wvec, Matrix_t & vmat);
       int    SVD_Inverse     (Matrix_t & m);

       void   Ramp            (                          const real8 x0=0.0, const real8 dx=1.0);
       void   Ramp            (const int w, const int h, const real8 x0=0.0, const real8 dx=1.0);

       real8  Residual        (const Matrix_t & ma);
       real8  Residual        (const Matrix_t & ma, const Matrix_t & mb);
       real8  Residual_Sum    (const Matrix_t & ma);
       real8  Residual_Sum    (const Matrix_t & ma, const Matrix_t & mb);

       real8  RMS_Error       (const Matrix_t & ma) const;
       real8  RMS_Error       (const Matrix_t & ma, const Matrix_t & mb) const;

       void   Threshold       (const real8 eps);
       void   Cap             (const real8 max);
       void   Constrain       (const real8 min, const real8 max);

       void   Sin             (void);
       void   Cos             (void);
       void   Tan             (void);
       void   Log             (void);
       void   Log2            (void);
       void   Log10           (void);
       void   Sqrt            (void);

       void   Identity_2x2    (void);
       void   Identity_3x3    (void);
       void   Identity_4x4    (void);

       void   Rotate_3x3      (const real8 rx, const real8 ry, const real8 rz);
       void   Rotate_3x3      (const real8 m0, const real8 m1, const real8 m2, const real8 m3, const real8 m4, const real8 m5);
       void   Rotate_3x3_X    (const real8 rx);
       void   Rotate_3x3_Y    (const real8 ry);
       void   Rotate_3x3_Z    (const real8 rz);

       void   Rotate_4x4      (const real8 rx, const real8 ry, const real8 rz);
       void   Rotate_4x4_X    (const real8 rx);
       void   Rotate_4x4_Y    (const real8 ry);
       void   Rotate_4x4_Z    (const real8 rz);

       void   Scale_2x2       (const real8 sx, const real8 sy);
       void   Scale_3x3       (const real8 sx, const real8 sy, const real8 sz);
       void   Scale_4x4       (const real8 sx, const real8 sy, const real8 sz);
       void   Translate_3x3   (const real8 tx, const real8 ty);
       void   Translate_4x4   (const real8 tx, const real8 ty, const real8 tz);

       void   Flip_H          (void);
       void   Flip_V          (void);
       void   Flip_H          (const Matrix_t & m);
       void   Flip_V          (const Matrix_t & m);

       void   Transform       (                    real8   (*f)  (                            const real8 s));
       void   Transform       (                    real8   (*f)  (const int iy, const int ix, const real8 s));
       void   Transform       (const Matrix_t & m, real8   (*f)  (const int iy, const int ix, const real8 s));

       void   Real_To_Complex    (const Matrix_t & m);
       void   Complex_Abs        (const Matrix_t & m);

       void   FFT_Transform      (const Matrix_t & m, const int dir, const int norm);
       void   FFT_Power_Spectrum (const Matrix_t & m);

       void   Extract         (const Matrix_t & m, const int x0, const int y0, const int w, const int h);
       void   Overlay         (const Matrix_t & m, const int x0, const int y0);

       void   Print           (FILE *fd=stdout, const char *fmt="%8.4lf ") const;
       void   Save            (const char *path) const;
       void   Load            (const char *path);

       int     MPI_Count      (void);
       size_t  MPI_Bytes      (void);
       real8  *MPI_Pack       (real8 *q=0);
       real8  *MPI_Unpack     (real8 *p=0);

#ifdef PARALLEL
       void    MPI_Send       (const int dst, const int mtag, MPI_Comm comm, const int sync);
       void    MPI_Receive    (const int src, const int mtag, MPI_Comm comm, const int sync);
#endif
};

extern int        LUP_Decompose   (real8 **a , int *p, int n, real8 tol );
extern int        LUP_Solve       (real8 **a , int *p, real8  *x, real8  *b, int n );
extern void       LUP_Invert      (real8 **ai, real8 **a, int *p, int n );
extern real8      LUP_Determinant (real8 **a , int *p, int n );

extern Matrix_t   M2x2_Identity   (void);
extern Matrix_t   M3x3_Identity   (void);
extern Matrix_t   M4x4_Identity   (void);

extern Matrix_t   M3x3_Rotate     (const real8 rx, const real8 ry, const real8 rz);
extern Matrix_t   M3x3_Rotate_X   (const real8 rx);
extern Matrix_t   M3x3_Rotate_Y   (const real8 ry);
extern Matrix_t   M3x3_Rotate_Z   (const real8 rz);

extern Matrix_t   M4x4_Rotate     (const real8 rx, const real8 ry, const real8 rz);
extern Matrix_t   M4x4_Rotate_X   (const real8 rx);
extern Matrix_t   M4x4_Rotate_Y   (const real8 ry);
extern Matrix_t   M4x4_Rotate_Z   (const real8 rz);

extern Matrix_t   M2x2_Scale      (const real8 sx, const real8 sy);
extern Matrix_t   M3x3_Scale      (const real8 sx, const real8 sy, const real8 sz);
extern Matrix_t   M4x4_Scale      (const real8 sx, const real8 sy, const real8 sz);
extern Matrix_t   M3x3_Translate  (const real8 tx, const real8 ty);
extern Matrix_t   M4x4_Translate  (const real8 tx, const real8 ty, const real8 tz);

extern void       M2x2_Identity   (Matrix_t & m);
extern void       M3x3_Identity   (Matrix_t & m);
extern void       M4x4_Identity   (Matrix_t & m);

extern void       M3x3_Rotate     (Matrix_t & m, const real8 rx, const real8 ry, const real8 rz);
extern void       M3x3_Rotate_X   (Matrix_t & m, const real8 rx);
extern void       M3x3_Rotate_Y   (Matrix_t & m, const real8 ry);
extern void       M3x3_Rotate_Z   (Matrix_t & m, const real8 rz);

extern void       M4x4_Rotate     (Matrix_t & m, const real8 rx, const real8 ry, const real8 rz);
extern void       M4x4_Rotate_X   (Matrix_t & m, const real8 rx);
extern void       M4x4_Rotate_Y   (Matrix_t & m, const real8 ry);
extern void       M4x4_Rotate_Z   (Matrix_t & m, const real8 rz);

extern void       M2x2_Scale      (Matrix_t & m, const real8 sx, const real8 sy);
extern void       M3x3_Scale      (Matrix_t & m, const real8 sx, const real8 sy, const real8 sz);
extern void       M4x4_Scale      (Matrix_t & m, const real8 sx, const real8 sy, const real8 sz);
extern void       M3x3_Translate  (Matrix_t & m, const real8 tx, const real8 ty);
extern void       M4x4_Translate  (Matrix_t & m, const real8 tx, const real8 ty, const real8 tz);

#endif
