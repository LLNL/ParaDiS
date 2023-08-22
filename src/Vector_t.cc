#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "mpi_portability.h"

#include "Typedefs.h"
#include "Vector_t.h"

#ifdef PARALLEL
#define VINIT     { vec_n=0; vec_max=0; vec=0; vec_mpi_cnt=0; vec_mpi_bytes=0; vec_mpi=0; }
#define VRECYCLE  { if (vec    ) delete [] vec    ;     \
                    if (vec_mpi) delete [] vec_mpi; VINIT; }
#else
#define VINIT     { vec_n=0; vec_max=0; vec=0; }
#define VRECYCLE  { if (vec    ) delete [] vec    ; VINIT; }
#endif

#define VALLOC(n) { vec = ( ((n)>0) ? new real8[(n)] : 0 ); vec_n = vec_max = ( vec ? (n) : 0 ); }

Vector_t::Vector_t(void)                                  { VINIT; }
Vector_t::Vector_t(const int n)                           { VINIT; VALLOC(n); if (vec      ) { memset(vec,0,n*sizeof(real8));                 } }
Vector_t::Vector_t(const int n, const          int     s) { VINIT; VALLOC(n); if (vec      ) { for (int i=0; (i<n); ++i) vec[i]=(real8) s   ; } }
Vector_t::Vector_t(const int n, const          real8   s) { VINIT; VALLOC(n); if (vec      ) { for (int i=0; (i<n); ++i) vec[i]=        s   ; } }
Vector_t::Vector_t(const int n, const unsigned char   *s) { VINIT; VALLOC(n); if (vec && s ) { for (int i=0; (i<n); ++i) vec[i]=(real8) s[i]; } }
Vector_t::Vector_t(const int n, const unsigned short  *s) { VINIT; VALLOC(n); if (vec && s ) { for (int i=0; (i<n); ++i) vec[i]=(real8) s[i]; } }
Vector_t::Vector_t(const int n, const          short  *s) { VINIT; VALLOC(n); if (vec && s ) { for (int i=0; (i<n); ++i) vec[i]=(real8) s[i]; } }
Vector_t::Vector_t(const int n, const unsigned int    *s) { VINIT; VALLOC(n); if (vec && s ) { for (int i=0; (i<n); ++i) vec[i]=(real8) s[i]; } }
Vector_t::Vector_t(const int n, const          int    *s) { VINIT; VALLOC(n); if (vec && s ) { for (int i=0; (i<n); ++i) vec[i]=(real8) s[i]; } }
Vector_t::Vector_t(const int n, const          real8  *s) { VINIT; VALLOC(n); if (vec && s ) { for (int i=0; (i<n); ++i) vec[i]=        s[i]; } }

Vector_t::Vector_t(const int n, real8 (*fn)(const int i)) { VINIT; VALLOC(n); if (vec && fn) { for (int i=0; (i<n); ++i) { vec[i]=(*fn)(i); } } }

Vector_t::Vector_t(const Vector_t & s)                    { VINIT; VALLOC(s.vec_n); if (vec && s.vec ) for (int i=0; (i<s.vec_n); i++) vec[i]=s.vec[i]; }

Vector_t::Vector_t(const real8 x, const real8 y, const real8 z )                { VINIT; VALLOC(3); if (vec) { vec[0]=x; vec[1]=y; vec[2]=z; } }
Vector_t::Vector_t(const real8 x, const real8 y, const real8 z, const real8 s ) { VINIT; VALLOC(4); if (vec) { vec[0]=x; vec[1]=y; vec[2]=z; vec[3]=s; } }

Vector_t::~Vector_t() { VRECYCLE; }

void Vector_t::Allocate (const int n)
{
   if (vec_n!=n)
   {
     if (vec) { delete [] vec; vec_n=0; vec_max=0; vec=0; }

     vec     = ( (n>0) ? new real8[n] : 0 );
     vec_n   = ( vec ? n : 0 );
     vec_max = vec_n;
   }

   if (vec) { memset(vec,0,vec_n*sizeof(real8)); }
}

void Vector_t::Append (const real8 s)
{
   if ( !vec || (vec_n==vec_max) )
   {
      vec_max = ( (vec_max==0) ? 64 : 2*vec_max );

      real8 *tmp = ( (vec_max>0) ? new real8[vec_max] : 0 );

      if (       tmp) { memset(tmp,  0,vec_max*sizeof(real8)); }
      if (vec && tmp) { memcpy(tmp,vec,vec_n  *sizeof(real8)); }

      if (vec) { delete [] vec; }

      vec = tmp;
   }

   if (vec && (vec_n<vec_max) ) { vec[vec_n++]=s; }
}

void Vector_t::Append (const real8 x, const real8 y, const real8 z)
{
   Append(x);
   Append(y);
   Append(z);
}

void Vector_t::Append (const real8 x, const real8 y, const real8 z, const real8 s)
{
   Append(x);
   Append(y);
   Append(z);
   Append(s);
}

const Vector_t & Vector_t::operator = (const Vector_t & v)
{
   if (this!=&v)
   {
      Allocate(v.vec_n);

      if (vec && v.vec) { memcpy(vec,v.vec,vec_n*sizeof(real8)); }
   }

   return(*this);
}

const Vector_t & Vector_t::operator = (const unsigned char   *s) { if (vec && s) { for (int i=0; (i<vec_n); ++i) vec[i]=(real8) s[i]; }  return(*this); }
const Vector_t & Vector_t::operator = (const unsigned short  *s) { if (vec && s) { for (int i=0; (i<vec_n); ++i) vec[i]=(real8) s[i]; }  return(*this); }
const Vector_t & Vector_t::operator = (const          short  *s) { if (vec && s) { for (int i=0; (i<vec_n); ++i) vec[i]=(real8) s[i]; }  return(*this); }
const Vector_t & Vector_t::operator = (const unsigned int    *s) { if (vec && s) { for (int i=0; (i<vec_n); ++i) vec[i]=(real8) s[i]; }  return(*this); }
const Vector_t & Vector_t::operator = (const          int    *s) { if (vec && s) { for (int i=0; (i<vec_n); ++i) vec[i]=(real8) s[i]; }  return(*this); }
const Vector_t & Vector_t::operator = (const          real8  *s) { if (vec && s) { for (int i=0; (i<vec_n); ++i) vec[i]=        s[i]; }  return(*this); }

const Vector_t & Vector_t::operator = (const          int     s) { if (vec     ) { for (int i=0; (i<vec_n); ++i) vec[i]=(real8) s;    }  return(*this); }
const Vector_t & Vector_t::operator = (const          real8   s) { if (vec     ) { for (int i=0; (i<vec_n); ++i) vec[i]=        s;    }  return(*this); }

int Vector_t::operator == (const Vector_t & v) const { return( Compare(v) ? 1 : 0 ); }
int Vector_t::operator == (const real8      s) const { return( Compare(s) ? 1 : 0 ); }
int Vector_t::operator == (const real8     *s) const { return( Compare(s) ? 1 : 0 ); }
int Vector_t::operator != (const Vector_t & v) const { return( Compare(v) ? 0 : 1 ); }
int Vector_t::operator != (const real8      s) const { return( Compare(s) ? 0 : 1 ); }
int Vector_t::operator != (const real8     *s) const { return( Compare(s) ? 0 : 1 ); }

#define VLOOP(fn) for (int i=0; (i<vec_n); ++i) { fn; }
 
int Vector_t::operator <  (const Vector_t & v) const { if (vec && v.vec && (vec_n>0) && (vec_n==v.vec_n)) { VLOOP( if (vec[i]>=v.vec[i]) return(0) ) return(1); } return(0); }
int Vector_t::operator <  (const real8      s) const { if (vec && (vec_n>0)     )                         { VLOOP( if (vec[i]>=s       ) return(0) ) return(1); } return(0); }
int Vector_t::operator <  (const real8     *s) const { if (vec && (vec_n>0) && s)                         { VLOOP( if (vec[i]>=s[i]    ) return(0) ) return(1); } return(0); }

int Vector_t::operator <= (const Vector_t & v) const { if (vec && v.vec && (vec_n>0) && (vec_n==v.vec_n)) { VLOOP( if (vec[i]> v.vec[i]) return(0) ) return(1); } return(0); }
int Vector_t::operator <= (const real8      s) const { if (vec && (vec_n>0)     )                         { VLOOP( if (vec[i]> s       ) return(0) ) return(1); } return(0); }
int Vector_t::operator <= (const real8     *s) const { if (vec && (vec_n>0) && s)                         { VLOOP( if (vec[i]> s[i]    ) return(0) ) return(1); } return(0); }

int Vector_t::operator >  (const Vector_t & v) const { if (vec && v.vec && (vec_n>0) && (vec_n==v.vec_n)) { VLOOP( if (vec[i]<=v.vec[i]) return(0) ) return(1); } return(0); }
int Vector_t::operator >  (const real8      s) const { if (vec && (vec_n>0)     )                         { VLOOP( if (vec[i]<=s       ) return(0) ) return(1); } return(0); }
int Vector_t::operator >  (const real8     *s) const { if (vec && (vec_n>0) && s)                         { VLOOP( if (vec[i]<=s[i]    ) return(0) ) return(1); } return(0); }

int Vector_t::operator >= (const Vector_t & v) const { if (vec && v.vec && (vec_n>0) && (vec_n==v.vec_n)) { VLOOP( if (vec[i]< v.vec[i]) return(0) ) return(1); } return(0); }
int Vector_t::operator >= (const real8      s) const { if (vec && (vec_n>0)     )                         { VLOOP( if (vec[i]< s       ) return(0) ) return(1); } return(0); }
int Vector_t::operator >= (const real8     *s) const { if (vec && (vec_n>0) && s)                         { VLOOP( if (vec[i]< s[i]    ) return(0) ) return(1); } return(0); }

Vector_t Vector_t::operator + (const Vector_t &  v) { Vector_t tmp(*this); tmp+=v; return(tmp); }
Vector_t Vector_t::operator + (const real8       s) { Vector_t tmp(*this); tmp+=s; return(tmp); }

Vector_t Vector_t::operator - (const Vector_t &  v) { Vector_t tmp(*this); tmp-=v; return(tmp); }
Vector_t Vector_t::operator - (const real8       s) { Vector_t tmp(*this); tmp-=s; return(tmp); }

Vector_t Vector_t::operator * (const Vector_t &  v) { Vector_t tmp(*this); tmp*=v; return(tmp); }
Vector_t Vector_t::operator * (const real8       s) { Vector_t tmp(*this); tmp*=s; return(tmp); }

// Vector_t Vector_t::operator / (const real8       s) {}

void Vector_t::operator += (const Vector_t & v) { if ( vec && v.vec && (vec_n>0) && (vec_n==v.vec_n) ) VLOOP( vec[i]+=v.vec[i] ) }
void Vector_t::operator += (const real8      s) { if ( vec && (vec_n>0)                              ) VLOOP( vec[i]+=s        ) }
void Vector_t::operator += (const real8     *s) { if ( vec && (vec_n>0) && s                         ) VLOOP( vec[i]+=s[i]     ) }

void Vector_t::operator -= (const Vector_t & v) { if ( vec && v.vec && (vec_n>0) && (vec_n==v.vec_n) ) VLOOP( vec[i]-=v.vec[i] ) }
void Vector_t::operator -= (const real8      s) { if ( vec && (vec_n>0)                              ) VLOOP( vec[i]-=s        ) }
void Vector_t::operator -= (const real8     *s) { if ( vec && (vec_n>0) && s                         ) VLOOP( vec[i]-=s[i]     ) }

void Vector_t::operator *= (const Vector_t & v) { if ( vec && v.vec && (vec_n>0) && (vec_n==v.vec_n) ) VLOOP( vec[i]*=v.vec[i] ) }
void Vector_t::operator *= (const real8      s) { if ( vec && (vec_n>0)                              ) VLOOP( vec[i]*=s        ) }
void Vector_t::operator *= (const real8     *s) { if ( vec && (vec_n>0) && s                         ) VLOOP( vec[i]*=s[i]     ) }

void Vector_t::operator /= (const real8      s) { if ( vec && (vec_n>0)                              ) { real8 t = ( (s!=0.0) ? (1.0/s) : 0.0); VLOOP(vec[i]*=t) } }

Vector_t Vector_t::operator - ()
{
   Vector_t tmp(*this);

   if(tmp.vec && (tmp.vec_n>0) )
      for (int i=0; (i<tmp.vec_n); ++i)
         tmp.vec[i]=-tmp.vec[i];

   return(tmp);
}

// void Vector_t::Add                (const Vector_t & v ) {}
// void Vector_t::Add                (const Vector_t & va, const Vector_t & vb) {}

// void Vector_t::Subtract           (const Vector_t & v ) {}
// void Vector_t::Subtract           (const Vector_t & va, const Vector_t & vb) {}

// void Vector_t::Multiply           (const Vector_t & v ) {}
// void Vector_t::Multiply           (const Vector_t & va, const Vector_t & vb) {}

real8 Vector_t::Dot (const Vector_t & v )
{
   real8 dot=0.0;

   if (vec && v.vec)
   {
      int n = ( (vec_n<v.vec_n) ? vec_n : v.vec_n );

      for (int i=0; (i<n); ++i) { dot+=(vec[i]*v.vec[i]); }
   }

   return(dot);
}

void Vector_t::Inverse (void) 
{
   if (vec && (vec_n>0))
      for (int i=0; (i<vec_n); ++i) { vec[i] = ( (fabs(vec[i])>0.0) ? 1.0/vec[i] : 0.0 ); }
}

void Vector_t::Inverse (const Vector_t & v ) 
{ 
   if (this==&v)  {          Inverse(); }
   else           { *this=v; Inverse(); }
}

void Vector_t::Reciprocol (void) 
{
   if (vec && (vec_n>0))
      for (int i=0; (i<vec_n); ++i) { vec[i] = ( (fabs(vec[i])>0.0) ? 1.0/vec[i] : 0.0 ); }
}

void Vector_t::Reciprocol (const Vector_t & v ) 
{ 
   if (this==&v)  {          Reciprocol(); }
   else           { *this=v; Reciprocol(); }
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

int Vector_t::Compare (const real8 s, const real8 eps) const
{
   int n = vec_n;

   if (vec && (n>0))
   {
      for (int i=0; (i<n); ++i)
         if (!nearly_eq(vec[i],s,eps)) return(0);

      return(1);
   }

   return(0);
}

int Vector_t::Compare (const real8 *s, const real8  eps) const
{
   int n = vec_n;

   if (vec && (n>0) && s)
   {
      for (int i=0; (i<n); ++i)
         if (!nearly_eq(vec[i],s[i],eps)) return(0);

      return(1);
   }

   return(0);
}

int Vector_t::Compare (const Vector_t & v, const real8 eps) const
{
   if (vec_n!=v.vec_n) return(0);

   if (vec && v.vec)
   {
      real8 *va =   vec;
      real8 *vb = v.vec;

      for (int i=0; (i<vec_n); ++i)
         if (!nearly_eq(va[i],vb[i],eps)) return(0);

      return(1);
   }

   return(0);
}


void Vector_t::Ramp (             const real8 x0, const real8 dx) { if (vec) { real8 x=x0; for (int i=0; (i<vec_n); ++i) { vec[i]=x; x+=dx; } } }
void Vector_t::Ramp (const int n, const real8 x0, const real8 dx) { Allocate(n); Ramp(x0,dx); }

void Vector_t::Rand (const real8 min, const real8 max)
{
   if (vec && (vec_n>0))
   {
      real8 s = ((real8)(max-min))/RAND_MAX;

      for (int i=0; (i<vec_n); ++i)
            vec[i]= min+(s*rand());
   }
}

void Vector_t::Rand (const int n, const real8 min, const real8 max) { Allocate(n); Rand(min,max); }

void Vector_t::Abs   (void) { if (vec && (vec_n>0)) { for (int i=0; (i<vec_n); ++i) vec[i]=fabs(vec[i]); } }
void Vector_t::Rint  (void) { if (vec && (vec_n>0)) { for (int i=0; (i<vec_n); ++i) vec[i]=rint(vec[i]); } }

void Vector_t::Round (const real8 dec)
{
   if (vec && (vec_n>0) && (fabs(dec)>0.0) )
      for (int i=0; (i<vec_n); ++i) vec[i]=round(vec[i]*dec)/dec;
}

void Vector_t::Sign  (void) { if (vec && (vec_n>0)) { for (int i=0; (i<vec_n); ++i) vec[i]=( (vec[i]>0.0) ? 1.0 : ( (vec[i]<0.0) ? -1.0 : 0.0 ) ); } }

real8 Vector_t::Min (void) const
{
   real8 min  = ( vec ? vec[0] : 0.0 );

   if (vec)
      for (int i=0; (i<vec_n); ++i) 
         if (vec[i]<min) { min=vec[i]; }

   return(min);
}

real8 Vector_t::Max (void) const
{
   real8 max  = ( vec ? vec[0] : 0.0 );

   if (vec)
      for (int i=0; (i<vec_n); ++i) 
         if (vec[i]>max) { max=vec[i]; }

   return(max);
}

real8 Vector_t::Min (int & indx) const
{
   real8 min  = ( vec ? vec[0] : 0.0 );
         indx = 0;

   if (vec)
      for (int i=0; (i<vec_n); ++i) 
         if (vec[i]<min) { min=vec[i]; indx=i; }

   return(min);
}

real8 Vector_t::Max (int & indx) const
{
   real8 max  = ( vec ? vec[0] : 0.0 );
         indx = 0;

   if (vec)
      for (int i=0; (i<vec_n); ++i) 
         if (vec[i]>max) { max=vec[i]; indx=i; }

   return(max);
}

real8 Vector_t::Mean (void) const
{
   real8 mean = 0.0;

   if (vec && (vec_n>0))
   {
      for (int i=0; (i<vec_n); ++i) 
         mean+=vec[i];

      mean /= (real8) vec_n;
   }

   return(mean);
}

real8 Vector_t::Range (void) const
{
   real8 min = ( vec ? vec[0] : 0.0 );
   real8 max = min;

   if (vec && (vec_n>0))
   {
      for (int i=0; (i<vec_n); ++i) 
      {
         min = ( (vec[i]<min) ? vec[i] : min );
         max = ( (vec[i]>max) ? vec[i] : max );
      }
   }

   return(max-min);
}

void Vector_t::Moment (real8 & min, real8 & max ) const 
{
   min = ( vec ? vec[0] : 0.0 );
   max = min;

   if (vec && (vec_n>0))
   {
      for (int i=0; (i<vec_n); ++i) 
      {
         min = ( (vec[i]<min) ? vec[i] : min );
         max = ( (vec[i]>max) ? vec[i] : max );
      }
   }
}

void Vector_t::Moment (real8 & min, real8 & max, real8 & mean ) const 
{
   min  = ( vec ? vec[0] : 0.0 );
   max  = min;
   mean = 0.0;

   if (vec && (vec_n>0))
   {
      for (int i=0; (i<vec_n); ++i) 
      {
         real8 s = vec[i];

         min   = ( (s<min) ? s : min );
         max   = ( (s>max) ? s : max );
         mean += s;
      }

      mean /= vec_n;
   }
}

void Vector_t::Moment (real8 & min, real8 & max, real8 & mean, real8 & range ) const 
{
   Moment(min,max,mean);  range=(max-min);
}

void Vector_t::Moment (real8 & min, real8 & max, real8 & mean, real8 & range, real8 & sdev) const 
{
   sdev = 0.0;

   Moment(min,max,mean);  range=(max-min);

}

real8 Vector_t::Residual (const Vector_t & v)
{
   int    n = ( (vec_n==v.vec_n) ? vec_n : 0 );
   real8  r = 0.0;
   real8 *a =   vec;
   real8 *b = v.vec;

   if (vec && v.vec && (n>0))
   {
      for (int i=0; (i<n); ++i)
      { a[i] = fabs(a[i]-b[i]); r+=a[i]; }
   }

   return(r);
}

real8 Vector_t::Residual (const Vector_t & va, const Vector_t & vb)
{
   int   n = ( (va.vec_n==vb.vec_n) ? va.vec_n : 0 );
   real8 r = 0.0;

   if ( (n>0) && (this!= &va) && (this!= &vb) )
      Allocate(n);

   real8 *a = ( (   vec && (   vec_n>0)) ?    vec : 0 );
   real8 *b = ( (va.vec && (va.vec_n>0)) ? va.vec : 0 );
   real8 *c = ( (vb.vec && (vb.vec_n>0)) ? vb.vec : 0 );

   if ( (n>0) && a && b && c )
   {
      for (int i=0; (i<n); ++i)
      { a[i] = fabs(b[i]-c[i]); r+=a[i]; }
   }

   return(r);
}

real8 Vector_t::Residual_Sum (const Vector_t & v)
{
   int    n    = ( (vec_n==v.vec_n) ? vec_n : 0 );
   real8  rsum = 0.0;
   real8 *a    = ( (  vec && (  vec_n>0)) ?   vec : 0 );
   real8 *b    = ( (v.vec && (v.vec_n>0)) ? v.vec : 0 );

   if ( a && b && (n>0) )
   {
      for (int i=0; (i<n); ++i)
      { rsum += fabs(a[i]-b[i]); }
   }

   return(rsum);
}

real8 Vector_t::Residual_Sum (const Vector_t & va, const Vector_t & vb)
{
   int    n    = ( (va.vec_n==vb.vec_n) ? va.vec_n : 0 );
   real8  rsum = 0.0;
   real8 *a    = ( (va.vec && (va.vec_n>0)) ? va.vec : 0 );
   real8 *b    = ( (vb.vec && (va.vec_n>0)) ? vb.vec : 0 );

   if ( a && b && (n>0) )
   {
      for (int i=0; (i<n); ++i)
      { rsum += fabs(a[i]-b[i]); }
   }

   return(rsum);
}

void Vector_t::Threshold (const real8 eps)                  { for (int i=0; (i<vec_n); ++i) { vec[i] = ( (fabs(vec[i])<eps) ? 0.0 : vec[i] ); } }
void Vector_t::Cap       (const real8 max)                  { for (int i=0; (i<vec_n); ++i) { vec[i] = ( (vec[i]<max) ? vec[i] : max );       } }
void Vector_t::Constrain (const real8 min, const real8 max) { for (int i=0; (i<vec_n); ++i) { vec[i] = ( (vec[i]<min) ? min : ( (vec[i]<max) ? vec[i] : max )); } }

void Vector_t::Sin   (void) { if (vec && (vec_n>0)) { for (int i=0; (i<vec_n); ++i) vec[i]=sin  (vec[i]); } }
void Vector_t::Cos   (void) { if (vec && (vec_n>0)) { for (int i=0; (i<vec_n); ++i) vec[i]=cos  (vec[i]); } }
void Vector_t::Tan   (void) { if (vec && (vec_n>0)) { for (int i=0; (i<vec_n); ++i) vec[i]=tan  (vec[i]); } }
void Vector_t::Log   (void) { if (vec && (vec_n>0)) { for (int i=0; (i<vec_n); ++i) vec[i]=log  (vec[i]); } }
void Vector_t::Log2  (void) { if (vec && (vec_n>0)) { for (int i=0; (i<vec_n); ++i) vec[i]=log2 (vec[i]); } }
void Vector_t::Log10 (void) { if (vec && (vec_n>0)) { for (int i=0; (i<vec_n); ++i) vec[i]=log10(vec[i]); } }
void Vector_t::Sqrt  (void) { if (vec && (vec_n>0)) { for (int i=0; (i<vec_n); ++i) vec[i]=sqrt (vec[i]); } }

void Vector_t::Flip (void) 
{
   if (vec && (vec_n>0) )
   {
      for (int i=0, j=(vec_n-1); (i<j); ++i,--j) 
      { real8 s=vec[i]; vec[i]=vec[j]; vec[j]=s; }
   }
}

void Vector_t::Flip (const Vector_t & v)
{
   if (this==&v) { Flip(); }
   else
   {
      Allocate(v.vec_n);

      if (vec && v.vec && (vec_n>0) )
         for (int i=0,j=(vec_n-1); (i<vec_n); ++i,--j) 
            vec[i]=v.vec[j];
   }
}

void Vector_t::Transform (real8 (*fn)(const int i, const real8 s)) { if (vec && (vec_n>0)) { for (int i=0; (i<vec_n); ++i) vec[i]=(*fn)(i,vec[i]); } }

void Vector_t::Real_To_Complex (const Vector_t & v)
{
   int n = v.vec_n;

   Allocate(2*n);

   if (vec && v.vec)
   {
      for (int i=0,k=0; (i<n); ++i, k+=2)
      {
         vec[k  ] = v.vec[i];
         vec[k+1] = 0.0;
      }
   }
}

void Vector_t::Complex_Abs (const Vector_t & v)
{
   int n = v.vec_n;

   Allocate(n/2);

   if (vec && v.vec && (n>0))
   {
      real8 *p = v.vec;
      real8 *q =   vec;

      for (int i=0,k=0; (i<(n/2)); ++i)
      {
         real8 a = p[k++];
         real8 b = p[k++];

         q[i] = sqrt(a*a+b*b);
      }
   }
}

// void Vector_t::FFT_Transform      (const Vector_t & v, const int dir, const int norm) {}
// void Vector_t::FFT_Power_Spectrum (const Vector_t & v) {}

// void Vector_t::Extract            (const Vector_t & m, const int i0, const int w ) {}
// void Vector_t::Overlay            (const Vector_t & m, const int i0 ) {}

void Vector_t::Print (FILE *fd, const char *fmt) const
{
   const char *fmts = ( fmt ? (char *) fmt : "%8.4lf " );

   if (fd && vec && (vec_n>0) )
   {
      for (int i=0; (i<vec_n); ++i)
      {
         fprintf(fd,fmts,vec[i]);
         if (((i+1)%8)==0) fprintf(fd,"\n");
      }

      if (vec_n%8) fprintf(fd,"\n");
   }
}

// void Vector_t::Save (const char *path) const {}
// void Vector_t::Load (const char *path) {}

int    Vector_t::MPI_Count (void) { return( (vec && (vec_n>0)) ? (vec_n)+1 : 0 ); }
size_t Vector_t::MPI_Bytes (void) { return(((vec && (vec_n>0)) ? (vec_n)+1 : 0 ) * sizeof(real8)); }

real8 *Vector_t::MPI_Pack (real8 *q)
{
   if (!q) { q = ( vec ? new real8[vec_n+1] : 0 ); }

   if (q)
   {
      *q++ = (real8) vec_n;

      if (vec && (vec_n>0)) { memcpy(q,vec,(vec_n)*sizeof(real8)); q+=vec_n; }
   }

   return(q);
}

real8 *Vector_t::MPI_Unpack (real8 *p)
{
   if (p)
   {
      int n = (int) *p++;

      Allocate(n);

      if (vec && (vec_n>0) ) { memcpy(vec,p,vec_n*sizeof(real8)); p+=vec_n; }
   }

   return(p);
}

#ifdef PARALLEL

void Vector_t::MPI_Send (const int dst, const int mtag, MPI_Comm comm, const int sync)
{
   if (vec_mpi) { delete [] vec_mpi; vec_mpi=0; vec_mpi_cnt=0; }

   vec_mpi_cnt = MPI_Count();
   vec_mpi     = MPI_Pack ();

   if (vec_mpi && (vec_mpi_cnt>0) )
   {
      if (sync) ::MPI_Send ( vec_mpi, vec_mpi_cnt, MPI_DOUBLE, dst, mtag, comm );
      else      ::MPI_Isend( vec_mpi, vec_mpi_cnt, MPI_DOUBLE, dst, mtag, comm, &vec_mpi_request);
   }
}

void Vector_t::MPI_Receive (const int src, const int mtag, MPI_Comm comm, const int sync)
{
   if (vec_mpi)
   {
      if (sync) ::MPI_Recv ( vec_mpi, vec_mpi_cnt, MPI_DOUBLE, src, mtag, comm, &vec_mpi_status  );
      else      ::MPI_Irecv( vec_mpi, vec_mpi_cnt, MPI_DOUBLE, src, mtag, comm, &vec_mpi_request );
   }
}

#endif

