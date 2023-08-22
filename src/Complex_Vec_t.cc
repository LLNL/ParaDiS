#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "mpi_portability.h"

#include "Typedefs.h"
#include "Complex_t.h"
#include "Complex_Vec_t.h"

//------------------------------------------------------------------------------------------------------------

#define VDELETE(a) if (a) { delete [] a; a=0; }

//------------------------------------------------------------------------------------------------------------

Complex_Vector_t::Complex_Vector_t(void)
{
   vc=0; vn=0;
}

//------------------------------------------------------------------------------------------------------------

Complex_Vector_t::Complex_Vector_t(const int n, const real8 re, const real8 im)
{
   vc = ( (n>0) ? new real8[2*n] : 0 );
   vn = (  vc   ? n : 0 );

   if (vc) { for (int i=0; (i<2*vn); i+=2) { vc[i]=re; vc[i+1]=im; } }
}

//------------------------------------------------------------------------------------------------------------

Complex_Vector_t::Complex_Vector_t(const int n)
{
   vc = ( (n>0) ? new real8[2*n] : 0 );
   vn = (  vc   ? n : 0 );

   if (vc) { for (int i=0; (i<2*vn); i+=2) { vc[i]=0.0; vc[i+1]=0.0; } }
}

//------------------------------------------------------------------------------------------------------------

Complex_Vector_t::Complex_Vector_t(const int n, const real8 re )
{
   vc = ( (n>0) ? new real8[2*n] : 0 );
   vn = (  vc   ? n : 0 );

   if (vc) { for (int i=0; (i<2*vn); i+=2) { vc[i]=re; vc[i+1]=0.0; } }
}

//------------------------------------------------------------------------------------------------------------

Complex_Vector_t::Complex_Vector_t(const int n, const Complex_t  & c)
{
   vc = ( (n>0) ? new real8[2*n] : 0 );
   vn = (  vc   ? n : 0 );

   const real8 re=c.re;
   const real8 im=c.im;

   if (vc) { for (int i=0; (i<2*vn); i+=2) { vc[i]=re; vc[i+1]=im; } }
}

//------------------------------------------------------------------------------------------------------------

Complex_Vector_t::Complex_Vector_t(const Complex_Vector_t & v)
{
   vc=0; vn=0;

   if (v.vc && (v.vn>0))
   {
      vc = new real8[2*v.vn];
      vn = v.vn;

      if (vc) { for (int i=0; (i<2*vn); i++) { vc[i]=v.vc[i]; } }
   }
}

Complex_Vector_t::Complex_Vector_t(const int n, Complex_t (*f) (const int i, const int imax) )
{
   vc = ( (n>0) ? new real8[2*n] : 0 );
   vn = (  vc   ? n : 0 );

   if (vc) 
   { 
      for (int i=0,j=0; (i<vn); i++, j+=2)
      { 
         Complex_t c = (*f)(i,n);
   
         vc[j]=c.re; vc[j+1]=c.im; 
      }
   } 
}

//------------------------------------------------------------------------------------------------------------

Complex_Vector_t::~Complex_Vector_t() { VDELETE(vc); vn=0; }

//------------------------------------------------------------------------------------------------------------

int Complex_Vector_t::operator == (const Complex_Vector_t  & s) const 
{
   if (vc && (vn>0) && s.vc && (s.vn>0) && (vn==s.vn) )
   {
      for (int i=0; (i<2*vn); i+=2)
      {
         if ( fabs(vc[i  ]-s.vc[i  ])>1.0e-12 ) return(0);
         if ( fabs(vc[i+1]-s.vc[i+1])>1.0e-12 ) return(0);
      }
 
      return(1);
   }

   return(0);
}

int Complex_Vector_t::operator == (const Complex_t & s) const 
{
   if (vc && (vn>0) )
   {
      for (int i=0; (i<2*vn); i+=2)
      {
         if ( fabs(vc[i  ]-s.re)>1.0e-12 ) return(0);
         if ( fabs(vc[i+1]-s.im)>1.0e-12 ) return(0);
      }
 
      return(1);
   }

   return(0);
}

int Complex_Vector_t::operator != (const Complex_Vector_t  & s) const { return( (*this==s) ? 0 : 1 ); } 
int Complex_Vector_t::operator != (const Complex_t         & s) const { return( (*this==s) ? 0 : 1 ); } 

//------------------------------------------------------------------------------------------------------------

const Complex_Vector_t & Complex_Vector_t::operator = (const Complex_Vector_t & v)
{
   if (this!=&v)
   {
      if (v.vn>0)
      {
         if (vc && (vn!=v.vn)) { delete [] vc; vc=0; }

         if (!vc) { vc = new real8[2*v.vn]; vn=v.vn; }
      }

      if ( vc && (vn>0) && v.vc && (v.vn>0) && (vn==v.vn) )
      { for (int i=0; (i<2*vn); i++) { vc[i]=v.vc[i]; } }
   }

   return(*this);
}

//------------------------------------------------------------------------------------------------------------

const Complex_Vector_t & Complex_Vector_t::operator = (const unsigned char        s) { if (vc && (vn>0)     ) {                         for (int i=0; (i<2*vn); i+=2) { vc[i]=s;    vc[i+1]=0.0;  } } return(*this); }
const Complex_Vector_t & Complex_Vector_t::operator = (const unsigned short       s) { if (vc && (vn>0)     ) {                         for (int i=0; (i<2*vn); i+=2) { vc[i]=s;    vc[i+1]=0.0;  } } return(*this); }
const Complex_Vector_t & Complex_Vector_t::operator = (const          short       s) { if (vc && (vn>0)     ) {                         for (int i=0; (i<2*vn); i+=2) { vc[i]=s;    vc[i+1]=0.0;  } } return(*this); }
const Complex_Vector_t & Complex_Vector_t::operator = (const unsigned int         s) { if (vc && (vn>0)     ) {                         for (int i=0; (i<2*vn); i+=2) { vc[i]=s;    vc[i+1]=0.0;  } } return(*this); }
const Complex_Vector_t & Complex_Vector_t::operator = (const          int         s) { if (vc && (vn>0)     ) {                         for (int i=0; (i<2*vn); i+=2) { vc[i]=s;    vc[i+1]=0.0;  } } return(*this); }
const Complex_Vector_t & Complex_Vector_t::operator = (const          real8       s) { if (vc && (vn>0)     ) {                         for (int i=0; (i<2*vn); i+=2) { vc[i]=s;    vc[i+1]=0.0;  } } return(*this); }
const Complex_Vector_t & Complex_Vector_t::operator = (const          Complex_t & s) { if (vc && (vn>0)     ) { real8 re=s.re, im=s.im; for (int i=0; (i<2*vn); i+=2) { vc[i]=re;   vc[i+1]=im;   } } return(*this); }

//------------------------------------------------------------------------------------------------------------

const Complex_Vector_t & Complex_Vector_t::operator = (const unsigned char       *s) { if (vc && (vn>0) && s) {                         for (int i=0; (i<2*vn); i+=2) { vc[i]=s[i];    vc[i+1]=0.0;     } } return(*this); }
const Complex_Vector_t & Complex_Vector_t::operator = (const unsigned short      *s) { if (vc && (vn>0) && s) {                         for (int i=0; (i<2*vn); i+=2) { vc[i]=s[i];    vc[i+1]=0.0;     } } return(*this); }
const Complex_Vector_t & Complex_Vector_t::operator = (const          short      *s) { if (vc && (vn>0) && s) {                         for (int i=0; (i<2*vn); i+=2) { vc[i]=s[i];    vc[i+1]=0.0;     } } return(*this); }
const Complex_Vector_t & Complex_Vector_t::operator = (const unsigned int        *s) { if (vc && (vn>0) && s) {                         for (int i=0; (i<2*vn); i+=2) { vc[i]=s[i];    vc[i+1]=0.0;     } } return(*this); }
const Complex_Vector_t & Complex_Vector_t::operator = (const          int        *s) { if (vc && (vn>0) && s) {                         for (int i=0; (i<2*vn); i+=2) { vc[i]=s[i];    vc[i+1]=0.0;     } } return(*this); }
const Complex_Vector_t & Complex_Vector_t::operator = (const          real8      *s) { if (vc && (vn>0) && s) {                         for (int i=0; (i<2*vn); i+=2) { vc[i]=s[i];    vc[i+1]=0.0;     } } return(*this); }
const Complex_Vector_t & Complex_Vector_t::operator = (const          Complex_t  *s) { if (vc && (vn>0) && s) {                         for (int i=0; (i<2*vn); i+=2) { vc[i]=s[i].re; vc[i+1]=s[i].im; } } return(*this); }

//------------------------------------------------------------------------------------------------------------

Complex_Vector_t Complex_Vector_t::operator +  (const Complex_Vector_t &  v) { Complex_Vector_t tmp(*this); tmp+=v;         return(tmp); }
Complex_Vector_t Complex_Vector_t::operator +  (const Complex_t        &  c) { Complex_Vector_t tmp(*this); tmp+=c;         return(tmp); }
Complex_Vector_t Complex_Vector_t::operator +  (const int                 s) { Complex_Vector_t tmp(*this); tmp+=(real8) s; return(tmp); }
Complex_Vector_t Complex_Vector_t::operator +  (const real8               s) { Complex_Vector_t tmp(*this); tmp+=s;         return(tmp); }

Complex_Vector_t Complex_Vector_t::operator -  (const Complex_Vector_t &  v) { Complex_Vector_t tmp(*this); tmp-=v;         return(tmp); }
Complex_Vector_t Complex_Vector_t::operator -  (const Complex_t        &  c) { Complex_Vector_t tmp(*this); tmp-=c;         return(tmp); }
Complex_Vector_t Complex_Vector_t::operator -  (const int                 s) { Complex_Vector_t tmp(*this); tmp-=(real8) s; return(tmp); }
Complex_Vector_t Complex_Vector_t::operator -  (const real8               s) { Complex_Vector_t tmp(*this); tmp-=s;         return(tmp); }

Complex_Vector_t Complex_Vector_t::operator *  (const Complex_Vector_t &  v) { Complex_Vector_t tmp(*this); tmp*=v;         return(tmp); }
Complex_Vector_t Complex_Vector_t::operator *  (const Complex_t        &  c) { Complex_Vector_t tmp(*this); tmp*=c;         return(tmp); }
Complex_Vector_t Complex_Vector_t::operator *  (const int                 s) { Complex_Vector_t tmp(*this); tmp*=(real8) s; return(tmp); }
Complex_Vector_t Complex_Vector_t::operator *  (const real8               s) { Complex_Vector_t tmp(*this); tmp*=s;         return(tmp); }

Complex_Vector_t Complex_Vector_t::operator /  (const Complex_Vector_t &  v) { Complex_Vector_t tmp(*this); tmp/=v;         return(tmp); }
Complex_Vector_t Complex_Vector_t::operator /  (const Complex_t        &  c) { Complex_Vector_t tmp(*this); tmp/=c;         return(tmp); }
Complex_Vector_t Complex_Vector_t::operator /  (const real8               s) { Complex_Vector_t tmp(*this); tmp/=s;         return(tmp); }

void             Complex_Vector_t::operator += (const Complex_Vector_t &  v) { if (vc && v.vc) {                         for (int i=0; (i<2*vn); i+=2) { vc[i]+=v.vc[i]; vc[i+1]+=v.vc[i+1];  } } }
void             Complex_Vector_t::operator += (const Complex_t        &  c) { if (vc)         { real8 re=c.re, im=c.im; for (int i=0; (i<2*vn); i+=2) { vc[i]+=re;      vc[i+1]+=im;         } } }
void             Complex_Vector_t::operator += (const real8               s) { if (vc)         {                         for (int i=0; (i<2*vn); i+=2) { vc[i]+=s;                            } } }

void             Complex_Vector_t::operator -= (const Complex_Vector_t &  v) { if (vc && v.vc) {                         for (int i=0; (i<2*vn); i+=2) { vc[i]-=v.vc[i]; vc[i+1]-=v.vc[i+1];  } } }
void             Complex_Vector_t::operator -= (const Complex_t        &  c) { if (vc)         { real8 re=c.re, im=c.im; for (int i=0; (i<2*vn); i+=2) { vc[i]-=re;      vc[i+1]-=im;         } } }
void             Complex_Vector_t::operator -= (const real8               s) { if (vc)         {                         for (int i=0; (i<2*vn); i+=2) { vc[i]-=s;                            } } }

void             Complex_Vector_t::operator += (const int                 s) { *this += (real8) s; }
void             Complex_Vector_t::operator -= (const int                 s) { *this -= (real8) s; }
void             Complex_Vector_t::operator *= (const int                 s) { *this *= (real8) s; }
void             Complex_Vector_t::operator /= (const int                 s) { *this /= (real8) s; }

void Complex_Vector_t::operator *= (const Complex_Vector_t &  v)
{
   if (vc && v.vc)
   {
      for (int i=0; (i<2*vn); i+=2)
      {
         real8 a=  vc[i  ];
         real8 b=  vc[i+1];
         real8 c=v.vc[i  ];
         real8 d=v.vc[i+1];

         vc[i  ] = (a*c)-(b*d);
         vc[i+1] = (a*d)+(b*c);
      }
   }
}

void Complex_Vector_t::operator *= (const Complex_t & s)
{
   if (vc)
   {
      real8 c=s.re;
      real8 d=s.im;

      for (int i=0; (i<2*vn); i+=2)
      {
         real8 a=vc[i  ];
         real8 b=vc[i+1];

         vc[i  ] = (a*c)-(b*d);
         vc[i+1] = (a*d)+(b*c);
      }
   }
}

void Complex_Vector_t::operator *= (const real8 s)
{
   if (vc)
   {
      for (int i=0; (i<2*vn); i+=2)
      {
         vc[i  ]*=s;
         vc[i+1]*=s;
      }
   }
}

void Complex_Vector_t::operator /= (const Complex_Vector_t &  v)
{
   if (vc && v.vc)
   {
      for (int i=0; (i<2*vn); i+=2)
      {
         real8 a=  vc[i  ];
         real8 b=  vc[i+1];
         real8 c=v.vc[i  ];
         real8 d=v.vc[i+1];
         real8 t=(c*c)+(d*d);
               t=( (t>0.0) ? (1.0/t) : 0.0 );

         vc[i  ] = (a*c+b*d)*t;
         vc[i+1] = (b*c-a*d)*t;
      }
   }
}

void Complex_Vector_t::operator /= (const Complex_t & s)
{
   if (vc)
   {
      real8 c=s.re;
      real8 d=s.im;
      real8 t=(c*c)+(d*d);
            t=( (t>0.0) ? (1.0/t) : 0.0 );

      for (int i=0; (i<2*vn); i+=2)
      {
         real8 a=vc[i  ];
         real8 b=vc[i+1];

         vc[i  ] = (a*c+b*d)*t;
         vc[i+1] = (b*c-a*d)*t;
      }
   }
}

void Complex_Vector_t::operator /= (const real8 s)
{
   if (vc && (fabs(s)>0.0) )
   {
      real8 t = (1.0/s);

      for (int i=0; (i<2*vn); i+=2)
      {
         vc[i  ]*=t;
         vc[i+1]*=t;
      }
   }
}

real8 *Complex_Vector_t::Cmplx_Abs (void)
{
   real8 *r = ( (vc && (vn>0)) ? new real8[vn] : 0 );

   if (r) 
   { 
      for (int i=0,j=0; (i<vn); i++, j+=2)  
      {
         real8 a = vc[j  ];
         real8 b = vc[j+1];

         r[i] = sqrt(a*a+b*b);
      }
   }
  
   return(r); 
}

Complex_Vector_t  Complex_Vector_t::Cmplx_Conj (void)
{
   Complex_Vector_t tmp(*this); 

   if (tmp.vc && vc)
   {
      for (int i=0; (i<2*vn); i+=2)
         tmp.vc[i+1] = -tmp.vc[i+1];   // im = -im
   }

   return(tmp); 
}

real8 *Cmplx_Abs (const Complex_Vector_t s)
{
   real8 *r = ( (s.vc && (s.vn>0)) ? new real8[s.vn] : 0 );

   if (r) 
   { 
      for (int i=0,j=0; (i<s.vn); i++, j+=2)  
      {
         real8 a = s.vc[j  ];
         real8 b = s.vc[j+1];

         r[i] = sqrt(a*a+b*b);
      }
   }
  
   return(r); 
}

Complex_Vector_t  Cmplx_Conj(const Complex_Vector_t s)
{
   Complex_Vector_t tmp(s); 

   if (tmp.vc)
   {
      for (int i=0; (i<2*s.vn); i+=2)
         tmp.vc[i+1] = -tmp.vc[i+1];   // im = -im
   }

   return(tmp); 
}

void Complex_Vector_t::Print (FILE *fd)
{
   if (fd && vc && (vn>0))
   {
      fprintf(fd,"#cmplx vn=%d\n", vn);
      fprintf(fd,"#cmplx vc\n"  );

      for (int i=0,j=0; (i<vn); i++, j+=2)
      {
         fprintf(fd,"%10.6lf %10.6lf%c", vc[j], vc[j+1], ( ((i+1)%8==0) ? '\n' : ' ') );
      }

      if (vn%8) fprintf(fd,"\n");
   }
}
