#pragma once

#ifndef __PDS_VECTOR_N_H
#define __PDS_VECTOR_N_H

#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "Typedefs.h"

class Vector_t
{
   public : 
      int          vec_n          ;  ///< vector length
      int          vec_max        ;  ///< vector length (max)
      real8       *vec            ;  ///< contiguous memory buffer allocated for this vector

#ifdef PARALLEL
      int          vec_mpi_cnt    ;  ///< count of packed MPI elements
      size_t       vec_mpi_bytes  ;  ///< bytes of packed MPI elements
      real8       *vec_mpi        ;  ///< buffer used for MPI communications
      MPI_Status   vec_mpi_status ;  ///< most recent MPI send/receive status
      MPI_Request  vec_mpi_request;  ///< most recent MPI send/receive request
#endif

   public : 
       Vector_t(void);
       Vector_t(const int n);
       Vector_t(const int n, const          int     s);
       Vector_t(const int n, const          real8   s);
       Vector_t(const int n, const unsigned char   *s);
       Vector_t(const int n, const unsigned short  *s);
       Vector_t(const int n, const          short  *s);
       Vector_t(const int n, const unsigned int    *s);
       Vector_t(const int n, const          int    *s);
       Vector_t(const int n, const          real8  *s);
       Vector_t(const int n, real8 (*fn)(const int i));

       Vector_t(const Vector_t & s);

       Vector_t(const real8 x, const real8 y, const real8 z );
       Vector_t(const real8 x, const real8 y, const real8 z, const real8 s );

      ~Vector_t();

       void  Recycle  (void);
       void  Allocate (const int   n);
       void  Append   (const real8 s);
       void  Append   (const real8 x, const real8 y, const real8 z);
       void  Append   (const real8 x, const real8 y, const real8 z, const real8 s);

       const Vector_t & operator = (const Vector_t       & m);
       const Vector_t & operator = (const int              s);
       const Vector_t & operator = (const real8            s);

       const Vector_t & operator = (const unsigned char   *s);
       const Vector_t & operator = (const unsigned short  *s);
       const Vector_t & operator = (const          short  *s);
       const Vector_t & operator = (const unsigned int    *s);
       const Vector_t & operator = (const          int    *s);
       const Vector_t & operator = (const          real8  *s);

       int      operator == (const Vector_t & m) const;
       int      operator == (const real8      s) const;
       int      operator == (const real8     *s) const;
       int      operator != (const Vector_t & m) const;
       int      operator != (const real8      s) const;
       int      operator != (const real8     *s) const;
       int      operator <  (const Vector_t & m) const;
       int      operator <  (const real8      s) const;
       int      operator <  (const real8     *s) const;
       int      operator <= (const Vector_t & m) const;
       int      operator <= (const real8      s) const;
       int      operator <= (const real8     *s) const;
       int      operator >  (const Vector_t & m) const;
       int      operator >  (const real8      s) const;
       int      operator >  (const real8     *s) const;
       int      operator >= (const Vector_t & m) const;
       int      operator >= (const real8      s) const;
       int      operator >= (const real8     *s) const;

       Vector_t operator +  (const Vector_t &  m);
       Vector_t operator +  (const real8       s);

       Vector_t operator -  (const Vector_t &  m);
       Vector_t operator -  (const real8       s);

       Vector_t operator *  (const Vector_t &  m);
       Vector_t operator *  (const real8       s);

       Vector_t operator /  (const real8       s);

       void     operator += (const Vector_t  & m);
       void     operator += (const real8       s);
       void     operator += (const real8      *s);
 
       void     operator -= (const Vector_t  & m);
       void     operator -= (const real8       s);
       void     operator -= (const real8      *s);
 
       void     operator *= (const Vector_t  & m);
       void     operator *= (const real8       s);
       void     operator *= (const real8      *s);
 
       void     operator /= (const real8       s);

       Vector_t operator - ();

       operator  real8  *() const { return ( (real8  *) vec ); }

       real8   Min                (const int i, int & indx);
       real8   Max                (const int i, int & indx);

       void    Add                (const Vector_t & m );
       void    Add                (const Vector_t & ma, const Vector_t & mb);

       void    Subtract           (const Vector_t & m );
       void    Subtract           (const Vector_t & ma, const Vector_t & mb);

       void    Multiply           (const Vector_t & m );
       void    Multiply           (const Vector_t & ma, const Vector_t & mb);

       real8   Dot                (const Vector_t & m );

       void    Inverse            (void);
       void    Inverse            (const Vector_t & m );

       void    Reciprocol         (void);
       void    Reciprocol         (const Vector_t & m );

       int     Compare            (const Vector_t & m, const real8  eps=4.0e-8) const;
       int     Compare            (const real8      s, const real8  eps=4.0e-8) const;
       int     Compare            (const real8     *s, const real8  eps=4.0e-8) const;

       void    Ramp               (             const real8 x0, const real8 dx);
       void    Ramp               (const int n, const real8 x0, const real8 dx);

       void    Rand               (             const real8 min, const real8 max);
       void    Rand               (const int n, const real8 min, const real8 max);

       void    Abs                (void);
       void    Rint               (void);
       void    Round              (const real8 eps);
       void    Sign               (void);

       real8   Min                (int & i) const;
       real8   Max                (int & i) const;
       real8   Min                (void) const;
       real8   Max                (void) const;
       real8   Mean               (void) const;
       real8   Range              (void) const;
       void    Moment             (real8 & min, real8 & max                                           ) const;
       void    Moment             (real8 & min, real8 & max, real8 & mean                             ) const;
       void    Moment             (real8 & min, real8 & max, real8 & mean, real8 & range              ) const;
       void    Moment             (real8 & min, real8 & max, real8 & mean, real8 & range, real8 & sdev) const;

       real8   Residual           (const Vector_t & va);
       real8   Residual           (const Vector_t & va, const Vector_t & vb);
       real8   Residual_Sum       (const Vector_t & va);
       real8   Residual_Sum       (const Vector_t & va, const Vector_t & vb);

       void    Threshold          (const real8 eps);
       void    Cap                (const real8 max);
       void    Constrain          (const real8 min, const real8 max);

       void    Sin                (void);
       void    Cos                (void);
       void    Tan                (void);
       void    Log                (void);
       void    Log2               (void);
       void    Log10              (void);
       void    Sqrt               (void);

       void    Flip               (void);
       void    Flip               (const Vector_t & m);

       void    Transform          (real8   (*f)  (const int i, const real8 s));

       void    Real_To_Complex    (const Vector_t & m);
       void    Complex_Abs        (const Vector_t & m);

       void    FFT_Transform      (const Vector_t & m, const int dir, const int norm);
       void    FFT_Power_Spectrum (const Vector_t & m);

       void    Extract            (const Vector_t & m, const int i0, const int w );
       void    Overlay            (const Vector_t & m, const int i0 );

       void    Print              (FILE *fd=stdout, const char *fmt="%8.4lf ") const;
       void    Save               (const char *path) const;
       void    Load               (const char *path);

       int     MPI_Count          (void);
       size_t  MPI_Bytes          (void);
       real8  *MPI_Pack           (real8 *q=0);
       real8  *MPI_Unpack         (real8 *p=0);

#ifdef PARALLEL
       void    MPI_Send           (const int dst, const int mtag, MPI_Comm comm, const int sync);
       void    MPI_Receive        (const int src, const int mtag, MPI_Comm comm, const int sync);
#endif
};

#endif  // __PDS_VECTOR_N_H
