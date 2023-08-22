#pragma once

#ifndef __PDS_COMPLEX_VECTOR_H
#define __PDS_COMPLEX_VECTOR_H

#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"
#include "cuda_portability.h"

#include "Typedefs.h"
#include "Complex_t.h"

//-----------------------------------------------------------------------------------------------
// Complex_Vector_t - implements a basic complex vector class
//-----------------------------------------------------------------------------------------------

class Complex_Vector_t
{
   public :
      real8   *vc;  ///< array of interleaved complex values
      int      vn;  ///< length of the array

   public :
       Complex_Vector_t(void);                                                        // constructors...
       Complex_Vector_t(const int n);                                                 //
       Complex_Vector_t(const int n, const real8 re                );                 //
       Complex_Vector_t(const int n, const real8 re, const real8 im);                 //
       Complex_Vector_t(const int n, const Complex_t            & c);                 //
       Complex_Vector_t(             const Complex_Vector_t     & v);                 //
       Complex_Vector_t(const int n, Complex_t (*f) (const int i, const int imax) );  //

      ~Complex_Vector_t();

             int                operator == (const Complex_Vector_t  &  s) const;      // comparison operators
             int                operator == (const Complex_t         &  s) const;      //
             int                operator != (const Complex_Vector_t  &  s) const;      //
             int                operator != (const Complex_t         &  s) const;      //

       const Complex_Vector_t & operator =  (const Complex_Vector_t  &  s);            // assignment from basic data types
       const Complex_Vector_t & operator =  (const unsigned char        s);            //
       const Complex_Vector_t & operator =  (const unsigned short       s);            //
       const Complex_Vector_t & operator =  (const          short       s);            //
       const Complex_Vector_t & operator =  (const unsigned int         s);            //
       const Complex_Vector_t & operator =  (const          int         s);            //
       const Complex_Vector_t & operator =  (const          real8       s);            //
       const Complex_Vector_t & operator =  (const          Complex_t & s);            //

       const Complex_Vector_t & operator =  (const unsigned char       *s);            // assignment from arrays of basic data types
       const Complex_Vector_t & operator =  (const unsigned short      *s);            // note - assumes the source array is at least as large
       const Complex_Vector_t & operator =  (const          short      *s);            //
       const Complex_Vector_t & operator =  (const unsigned int        *s);            //
       const Complex_Vector_t & operator =  (const          int        *s);            //
       const Complex_Vector_t & operator =  (const          real8      *s);            //
       const Complex_Vector_t & operator =  (const          Complex_t  *s);            //

       // complex addition, subtraction, multiplication...

             Complex_Vector_t   operator +  (const Complex_Vector_t &  v);
             Complex_Vector_t   operator +  (const Complex_t        &  c);
             Complex_Vector_t   operator +  (const int                 c);
             Complex_Vector_t   operator +  (const real8               c);

             Complex_Vector_t   operator -  (const Complex_Vector_t &  v);
             Complex_Vector_t   operator -  (const Complex_t        &  c);
             Complex_Vector_t   operator -  (const int                 c);
             Complex_Vector_t   operator -  (const real8               c);

             Complex_Vector_t   operator *  (const Complex_Vector_t &  v);
             Complex_Vector_t   operator *  (const Complex_t        &  c);
             Complex_Vector_t   operator *  (const int                 c);
             Complex_Vector_t   operator *  (const real8               c);

             Complex_Vector_t   operator /  (const Complex_Vector_t &  v);
             Complex_Vector_t   operator /  (const Complex_t        &  c);
             Complex_Vector_t   operator /  (const int                 c);
             Complex_Vector_t   operator /  (const real8               c);

             void               operator += (const Complex_Vector_t &  v);
             void               operator += (const Complex_t        &  c);
             void               operator += (const int                 c);
             void               operator += (const real8               c);

             void               operator -= (const Complex_Vector_t &  v);
             void               operator -= (const Complex_t        &  c);
             void               operator -= (const int                 c);
             void               operator -= (const real8               c);

             void               operator *= (const Complex_Vector_t &  v);
             void               operator *= (const Complex_t        &  c);
             void               operator *= (const int                 c);
             void               operator *= (const real8               c);

             void               operator /= (const Complex_Vector_t &  v);
             void               operator /= (const Complex_t        &  c);
             void               operator /= (const int                 c);
             void               operator /= (const real8               c);

             real8            *Cmplx_Abs  (void);   // returns vector of complex absolute values (allocated)
             Complex_Vector_t  Cmplx_Conj (void);   // returns vector of complex conjugates

             void              Print     (FILE *fd=stdout);

      operator  real8  *() const { return ( (real8  *) vc ); }

      Complex_t operator[] (int i) const { return ( vc && (0<=i) && (i<vn) ? Complex_t(vc[2*i],vc[2*i+1]) : Complex_t() ); }   // (note: RHS only)
};

//-----------------------------------------------------------------------------------------------

extern real8            *Cmplx_Abs  (const Complex_Vector_t s);
extern Complex_Vector_t  Cmplx_Conj (const Complex_Vector_t s);

//-----------------------------------------------------------------------------------------------

#endif   // __PDS_COMPLEX_VECTOR_H
