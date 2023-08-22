#pragma once

#ifndef __PDS_COMPLEX_TYPE_H
#define __PDS_COMPLEX_TYPE_H

#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"
#include "cuda_portability.h"

#include "Typedefs.h"

//-----------------------------------------------------------------------------------------------
// Complex_t - implements a basic complex class
//-----------------------------------------------------------------------------------------------

class Complex_t
{
   public :
       real8   re;  ///< real      component of complex number
       real8   im;  ///< imaginary component of complex number

   public :
       __cuda_hdev__ Complex_t(void)                         { re=0.0 ; im=0.0 ; }  // constructors...
       __cuda_hdev__ Complex_t(const real8 r )               { re=r   ; im=0.0 ; }  //    |
       __cuda_hdev__ Complex_t(const real8 r, const real8 i) { re=r   ; im=i   ; }  //    |
       __cuda_hdev__ Complex_t(const Complex_t & c)          { re=c.re; im=c.im; }  //    |

       __cuda_hdev__ ~Complex_t() {}  // destructor (nothing allocated)

       __cuda_hdev__ real8 Real(void) const { return(re); }  // accessors
       __cuda_hdev__ real8 Imag(void) const { return(im); }  //    |

       __cuda_hdev__ const Complex_t & operator =  (const Complex_t     & s) { re=s.re; im=s.im; return(*this); }  // assignment from basic data types
       __cuda_hdev__ const Complex_t & operator =  (const unsigned char   s) { re=s   ; im=0.0;  return(*this); }  //    |
       __cuda_hdev__ const Complex_t & operator =  (const unsigned short  s) { re=s   ; im=0.0;  return(*this); }  //    |
       __cuda_hdev__ const Complex_t & operator =  (const          short  s) { re=s   ; im=0.0;  return(*this); }  //    |
       __cuda_hdev__ const Complex_t & operator =  (const unsigned int    s) { re=s   ; im=0.0;  return(*this); }  //    |
       __cuda_hdev__ const Complex_t & operator =  (const          int    s) { re=s   ; im=0.0;  return(*this); }  //    |
       __cuda_hdev__ const Complex_t & operator =  (const          real8  s) { re=s   ; im=0.0;  return(*this); }  //    |

       __cuda_hdev__       int       operator == (const Complex_t  & s) const { return( (fabs(re-s.re)>1e-10) || (fabs(im-s.im)>1e-10) ? 0 : 1 ); }  // comparison operators
       __cuda_hdev__       int       operator == (const real8        s) const { return( (fabs(re-s   )>1e-10) || (fabs(im-s   )>1e-10) ? 0 : 1 ); }  //     |
       __cuda_hdev__       int       operator != (const Complex_t  & s) const { return( (fabs(re-s.re)>1e-10) || (fabs(im-s.im)>1e-10) ? 1 : 0 ); }  //     |
       __cuda_hdev__       int       operator != (const real8        s) const { return( (fabs(re-s   )>1e-10) || (fabs(im-s   )>1e-10) ? 1 : 0 ); }  //     |

       // complex addition, subtraction, multiplication...

       __cuda_hdev__ Complex_t operator +  (const Complex_t &  s) { return( Complex_t(re+s.re, im+s.im) ); }
       __cuda_hdev__ Complex_t operator +  (const real8        s) { return( Complex_t(re+s   , im     ) ); }

       __cuda_hdev__ Complex_t operator -  (const Complex_t &  s) { return( Complex_t(re-s.re, im-s.im) ); }
       __cuda_hdev__ Complex_t operator -  (const real8        s) { return( Complex_t(re-s   , im     ) ); }

       __cuda_hdev__ Complex_t operator *  (const Complex_t &  s) { real8 a=re, b=im, c=s.re, d=s.im; return( Complex_t( (a*c)-(b*d), (a*d)+(b*c) ) ); }
       __cuda_hdev__ Complex_t operator *  (const real8        s) { real8 a=re, b=im, c=s;            return( Complex_t( (a*c)      ,       (b*c) ) ); }

       __cuda_hdev__ Complex_t operator /  (const Complex_t &  s) { real8 a=re, b=im, c=s.re, d=s.im; real8 t=(c*c)+(d*d);  t=( (t>0.0) ? (1.0/t) : 0.0 ); return( Complex_t( (a*c+b*d)*t, (b*c-a*d)*t ) ); }
       __cuda_hdev__ Complex_t operator /  (const real8        s) { real8 a=re, b=im, c=s;            real8 t=(c*c);        t=( (t>0.0) ? (1.0/t) : 0.0 ); return( Complex_t( (a*c    )*t, (b*c    )*t ) ); }

       __cuda_hdev__ void      operator += (const Complex_t  & s) { re+=s.re; im+=s.im; }
       __cuda_hdev__ void      operator += (const real8        s) { re+=s   ;           }

       __cuda_hdev__ void      operator -= (const Complex_t  & s) { re-=s.re; im-=s.im; }
       __cuda_hdev__ void      operator -= (const real8        s) { re-=s   ;           }

       __cuda_hdev__ void      operator *= (const Complex_t  & s) { real8 a=re, b=im, c=s.re, d=s.im;  re=(a*c)-(b*d); im=(a*d)+(b*c); }
       __cuda_hdev__ void      operator *= (const real8        s) { real8 a=re, b=im, c=s;             re=(a*c)      ; im=      (b*c); }

       __cuda_hdev__ void      operator /= (const Complex_t &  s) { real8 a=re, b=im, c=s.re, d=s.im; real8 t=(c*c)+(d*d);  t=( (t>0.0) ? (1.0/t) : 0.0 ); re=(a*c+b*d)*t; im=(b*c-a*d)*t; }
       __cuda_hdev__ void      operator /= (const real8        s) { real8 a=re, b=im, c=s;            real8 t=(c*c)      ;  t=( (t>0.0) ? (1.0/t) : 0.0 ); re=(a*c    )*t; im=(b*c    )*t; }

       __cuda_hdev__ Complex_t operator -  ()   { return( Complex_t( -re, -im ) ); } // (negation)
};

//-----------------------------------------------------------------------------------------------

__cuda_hdev__ static inline Complex_t operator + (const Complex_t & ca, const Complex_t & cb) { return( Complex_t( (ca.re+cb.re), ( ca.im+cb.im) ) ); }
__cuda_hdev__ static inline Complex_t operator + (const Complex_t & ca, const real8       s ) { return( Complex_t( (ca.re+s    ), ( ca.im      ) ) ); }
__cuda_hdev__ static inline Complex_t operator + (const real8       s , const Complex_t & ca) { return( Complex_t( (ca.re+s    ), ( ca.im      ) ) ); }

__cuda_hdev__ static inline Complex_t operator - (const Complex_t & ca, const Complex_t & cb) { return( Complex_t( (ca.re-cb.re), ( ca.im-cb.im) ) ); }
__cuda_hdev__ static inline Complex_t operator - (const Complex_t & ca, const real8       s ) { return( Complex_t( (ca.re-s    ), ( ca.im      ) ) ); }
__cuda_hdev__ static inline Complex_t operator - (const real8       s , const Complex_t & ca) { return( Complex_t( (s-ca.re    ), (-ca.im      ) ) ); }

__cuda_hdev__ static inline Complex_t operator * (const Complex_t & ca, const Complex_t & cb) { real8 a=ca.re, b=ca.im, c=cb.re, d=cb.im; return( Complex_t( (a*c)-(b*d), (a*d)+(b*c) ) ); }
__cuda_hdev__ static inline Complex_t operator * (const Complex_t & ca, const real8       s ) {                                           return( Complex_t( (ca.re*s)  , (ca.im*s)   ) ); }
__cuda_hdev__ static inline Complex_t operator * (const real8       s , const Complex_t & ca) {                                           return( Complex_t( (ca.re*s)  , (ca.im*s)   ) ); }

//-----------------------------------------------------------------------------------------------

__cuda_hdev__ static inline real8     creal      (const Complex_t & s) { return( s.re ); }
__cuda_hdev__ static inline real8     cimag      (const Complex_t & s) { return( s.im ); }

__cuda_hdev__ static inline real8     Cmplx_Abs  (const Complex_t & s) { return( sqrt(s.re*s.re+s.im*s.im) ); }
__cuda_hdev__ static inline Complex_t Cmplx_Conj (const Complex_t & s) { return( Complex_t(s.re,-s.im)     ); }
__cuda_hdev__ static inline Complex_t Cmplx_Exp  (const Complex_t & s) { return( Complex_t(cos(s.im), sin(s.im)) * exp(s.re) ); }

__cuda_hdev__ static inline Complex_t Cmplx_Sqrt (const Complex_t & s)
{
    real8 r  = sqrt( (s.re*s.re) + (s.im*s.im) );
    real8 rc = sqrt((r + s.re)/2.0);
    real8 ic = sqrt((r - s.re)/2.0);
          ic = ( (s.im >= 0.0) ? ic : -ic );

    return( Complex_t(rc,ic) );
}

__cuda_hdev__ static inline Complex_t Cmplx_Dot (const real8 *va, const Complex_t *vb)
{
   real8 a,c,d,re,im;

   a=va[0]; c=vb[0].re; d=vb[0].im;  re =(a*c); im =(a*d);
   a=va[1]; c=vb[1].re; d=vb[1].im;  re+=(a*c); im+=(a*d);
   a=va[2]; c=vb[2].re; d=vb[2].im;  re+=(a*c); im+=(a*d);

   return( Complex_t(re,im) );
}

__cuda_hdev__ static inline Complex_t  Cmplx_Dot (const Complex_t  *va, const Complex_t  *vb)
{
   real8 a,b,c,d,re,im;

   a=va[0].re; b=va[0].im; c=vb[0].re; d=vb[0].im;  re =(a*c)-(b*d); im =(a*d)+(b*c);
   a=va[1].re; b=va[1].im; c=vb[1].re; d=vb[1].im;  re+=(a*c)-(b*d); im+=(a*d)+(b*c);
   a=va[2].re; b=va[2].im; c=vb[2].re; d=vb[2].im;  re+=(a*c)-(b*d); im+=(a*d)+(b*c);

   return( Complex_t(re,im) );
}

#endif   // __PDS_COMPLEX_TYPE_H
