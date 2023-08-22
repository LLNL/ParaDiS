#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "mpi_portability.h"
#include "cuda_portability.h"

#include "Append_v.h"

//------------------------------------------------------------------------------------------------------------
// Scattered throughout the ParaDiS code are multiple instances of dynamically growing vectors
// of various data entities.  The implementations included in this module cover implementations
// that use the basic data types within ParaDiS.
//
// The basic premise is given an instance of a variable, the variable is added (appended) to an
// existing array.  The existing array may (or may not) contain existing elements.
//
// Whenever the existing vector is full, we must allocate a new vector with additional space, copy the 
// old vector into the new vector, and delete the old vector. We then append the data item to the end 
// of the newly allocated space.  Obviously - this can be horribly inefficient. To reduce the frequency
// of the occurrences of the allocation and copies, more space is allocated than necessary.  This enables
// subsequent calls to simply append data values until the space is exhausted.
//
// Note - three tuning parameters are available and defaulted in the API
//    tmin : the initial vector allocation size  (default=32)
//    t0   : growth factor numerator             (default=1 )
//    t1   : growth factor denominator           (default=1 )
//
// If you know in advance the approximate size of the vector - use the tmin parameter to set the 
// size of the initial allocation.  
// 
// If you know in advance the growth behavior of the overall vector - use the t0 and t1 parameters 
// to control the growth of the buffers.  
//
// Growth is computed using the following formula,
//    int n = { current vector size };
//    { new vector size } = n + (t0*n)/t1;
//
// By default, t0 and t1 are set to 1, which means that each time the vector is extended - it is
// doubled in size (n+n). If you know that the vector will grow quickly, set (t0>t1). If you know
// that the vector will grow slowly, set (t0<t1).
//
// Arguments...
//    vtype  *v   ,    ///< points to existing vector
//    int &   vcnt,    ///< number of elements in existing vector
//    int &   vmax,    ///< allocated size of the existing vector
//    vtype   val ,    ///< element to append
//    int     tmin,    ///< initial vector allocation size (default=32)
//    int     t0  ,    ///< growth factor numerator        (default=1)
//    int     t1       ///< growth factor denominator      (default=1)
//   
// Note - the various versions of Append(...) assume the source vector is assigned to the result.
//
// example...
//
//   float *v=0;
//   int    vn=0, vmax=0;
//
//   for (int i=0; (i<1000); i++)
//   { v = Append(v,vn,vmax, (float) i); }
// 
// ...will result in a vector of 1000 floats in an allocation of 1024 floats. 
// Furthermore, the vector will get allocated, copied, and deallocated 6 times 
// through the progression (32,64,128,256,512,and 1024).
//  
// Judicious use fo the tuning parameters can significantly cut down on the reallocations.
// Poor use of the tuning parameters can result in wasted space or repeated and unnecessary 
// reallocations (or both).
//
// Lastly - obviously - this code should be templated rather than implemented as expanded macros.
// If implemented as templates, some of the more complex data structures in ParaDiS such as
// segments or nodes could use this logic.  That implementation is left for later.
//------------------------------------------------------------------------------------------------------------

#define APPEND(vtype) \
vtype *Append (vtype *v, int & vcnt, int & vmax, vtype val, int tmin, int t0, int t1) \
{                                                            \
   if ( !v || (vcnt==vmax) )                                 \
   {                                                         \
      tmin = ( (tmin<1) ? 1 : tmin );                        \
      t0   = ( (t0  <1) ? 1 : t0   );                        \
      t1   = ( (t1  <1) ? 1 : t1   );                        \
                                                             \
      vmax += ( (vmax==0) ? tmin : (t0*vmax)/t1 );           \
                                                             \
      vtype *tmp = ( (vmax>0) ? new vtype[vmax] : 0 );       \
                                                             \
      if (tmp)                                               \
      {                                                      \
         if (v)         { memcpy(tmp     ,v,(     vcnt)*sizeof(vtype));  } \
         if (vcnt<vmax) { memset(tmp+vcnt,0,(vmax-vcnt)*sizeof(vtype));  } \
      }                                                      \
                                                             \
      if (v) { delete [] v; }                                \
                                                             \
      v = tmp;                                               \
   }                                                         \
                                                             \
   if (v && (vcnt<vmax) ) { v[vcnt]=val; vcnt++; }           \
                                                             \
   return(v);                                                \
}

APPEND(         char )
APPEND(unsigned char )
APPEND(         short)
APPEND(unsigned short)
APPEND(         int  )
APPEND(unsigned int  )
APPEND(         long )
APPEND(unsigned long )
APPEND(         long long   )
APPEND(unsigned long long   )
APPEND(         float       )
APPEND(         double      )
APPEND(         long double )
APPEND(         void *      )   // (vector of pointers)

#undef APPEND

// Append (block form)
//
// This version of Append allows the caller to append multiple values to the end 
// of an existing array.
//
// Arguments...
//    vtype  *v   ,    ///< points to existing vector
//    int &   vcnt,    ///< number of elements in existing vector
//    int &   vmax,    ///< allocated size of the existing vector
//    vtype  *vsrc,    ///< source block of elements to append
//    int     vn  ,    ///< number of elements to append
//    int     tmin,    ///< initial vector allocation size
//    int     t0  ,    ///< growth factor numerator 
//    int     t1       ///< growth factor denominator
//------------------------------------------------------------------------------------------------------------

#define APPEND(vtype) \
vtype *Append (vtype *v, int & vcnt, int & vmax, vtype *vsrc, int vn, int tmin, int t0, int t1) \
{                                                               \
   if (vsrc && (vn>0))                                          \
   {                                                            \
      if ( !v || ((vcnt+vn)>vmax) )                             \
      {                                                         \
         tmin = ( (tmin<1) ? 1 : tmin );                        \
         t0   = ( (t0  <1) ? 1 : t0   );                        \
         t1   = ( (t1  <1) ? 1 : t1   );                        \
                                                                \
         vmax += ( (vmax==0) ? tmin : (t0*vmax)/t1 );           \
         vmax  = ( ((vcnt+vn)<vmax) ? vmax : (vcnt+vn) );       \
                                                                \
         vtype *tmp = ( (vmax>0) ? new vtype[vmax] : 0 );       \
                                                                \
         if (tmp)                                               \
         {                                                      \
            if (v)         { memcpy(tmp     ,v,(     vcnt)*sizeof(vtype));  }  \
            if (vcnt<vmax) { memset(tmp+vcnt,0,(vmax-vcnt)*sizeof(vtype));  }  \
         }                                                      \
                                                                \
         if (v) { delete [] v; }                                \
                                                                \
         v = tmp;                                               \
      }                                                         \
                                                                \
      if (v && ((vcnt+vn)<vmax) )                               \
      {                                                         \
         memcpy(v+vcnt,vsrc,vn*sizeof(vtype));  vcnt+=vn;       \
      }                                                         \
   }                                                            \
                                                                \
   return(v);                                                   \
}

APPEND(         char )
APPEND(unsigned char )
APPEND(         short)
APPEND(unsigned short)
APPEND(         int  )
APPEND(unsigned int  )
APPEND(         long )
APPEND(unsigned long )
APPEND(         long long   )
APPEND(unsigned long long   )
APPEND(         float       )
APPEND(         double      )
APPEND(         long double )
APPEND(         void *      )   // (vector of pointers)

#undef APPEND

//------------------------------------------------------------------------------------------------------------
// Append()
//
// Generic version for arbitraty sized structs.
//
// A shallow copy of the structure will be appended to the existing array. 
// If the structure contains pointers to memory, only the pointers will be copied.
//------------------------------------------------------------------------------------------------------------

void *Append 
(
    void    *v    ,   ///< points to existing vector              (may be deleted)
    int    & vcnt ,   ///< number of elements in existing vector  (may be updated)
    int    & vmax ,   ///< allocated size of the existing vector  (may be updated)
    void    *val  ,   ///< points to an instance of a structure to be appended (shallow copied)
    size_t   vsize,   ///< size of each structure (bytes)
    int      tmin ,   ///< initial vector allocation size (default=32)
    int      t0   ,   ///< growth factor numerator        (default=1 )
    int      t1       ///< growth factor denominator      (default=1 )
                      ///     note - vector growth n = n+(t0*n)/t1;  (default n=(n+n)=(2n) )
)
{
   if (val && (vsize>0))
   {
      if ( !v || (vcnt==vmax) )
      {
         tmin = ( (tmin<1) ? 1 : tmin );
         t0   = ( (t0  <1) ? 1 : t0   );
         t1   = ( (t1  <1) ? 1 : t1   );

         vmax += ( (vmax==0) ? tmin : (t0*vmax)/t1 );
   
         unsigned char *tmp = ( (vmax>0) ? new unsigned char[vmax*vsize] : 0 );

         if (tmp)
         {
            unsigned char *p   = tmp;
            unsigned char *pe  = tmp + vmax*vsize;
   
            if (v)    { memcpy(p,v,vcnt*vsize);  p+=vcnt*vsize; }  // copy existing members to new vector
            if (p<pe) { memset(p,0,(pe-p));                     }  // set remaining buffer to zero
         }

         if (v) { delete [] ((unsigned char *) v); }
   
         v = (void *) tmp;
      }
   
      if (v && (vcnt<vmax) ) 
      { 
         unsigned char *p = ((unsigned char *) val);
         unsigned char *q = ((unsigned char *) v  ) + (vcnt*vsize);  vcnt++;
         memcpy(q,p,vsize);
      }
   }

   return(v);
}

