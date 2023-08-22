#pragma once

#ifndef _PDS_SSF_FV_H
#define _PDS_SSF_FV_H

#include "cuda_portability.h"

#include "Typedefs.h"
#include "V3.h"

#define __inline__ inline __attribute__((__always_inline__))

//------------------------------------------------------------------------------------------------------------
// SSF_FV_t class support
//
// This class was constructed to facilitate host-device data transfers 
//------------------------------------------------------------------------------------------------------------

class SSF_FV_t
{
   public:
      real8 f1[3];  ///< force on node 1 <xyz>
      real8 f2[3];  ///< force on node 2 <xyz>
      real8 f3[3];  ///< force on node 3 <xyz>
      real8 f4[3];  ///< force on node 4 <xyz>

   public:
      __cuda_hdev__  SSF_FV_t(void);
      __cuda_hdev__  SSF_FV_t(const SSF_FV_t & fv);

      __cuda_hdev__ ~SSF_FV_t() {}

      __cuda_hdev__  const SSF_FV_t & operator =  (const SSF_FV_t & fv);
      __cuda_hdev__        int        operator == (const SSF_FV_t & fv) const;
};

__cuda_hdev__  
__inline__ 
SSF_FV_t::SSF_FV_t(void) 
{ 
   V3_ZERO(f1);
   V3_ZERO(f2);
   V3_ZERO(f3);
   V3_ZERO(f4); 
}

__cuda_hdev__  
__inline__ 
SSF_FV_t::SSF_FV_t(const SSF_FV_t & fv) 
{ 
   V3_COPY(f1,fv.f1);
   V3_COPY(f2,fv.f2);
   V3_COPY(f3,fv.f3);
   V3_COPY(f4,fv.f4); 
}

__cuda_hdev__  
__inline__ 
const SSF_FV_t & SSF_FV_t::operator = (const SSF_FV_t & fv) 
{ 
   V3_COPY(f1,fv.f1);
   V3_COPY(f2,fv.f2);
   V3_COPY(f3,fv.f3);
   V3_COPY(f4,fv.f4); 
   return(*this); 
}

__cuda_hdev__
__inline__ 
int SSF_FV_t::operator == (const SSF_FV_t& fv) const
{
   const real8 eps=1.0e-6;

   if (    ( fabs(f1[0]-fv.f1[0]) > eps )
        || ( fabs(f1[1]-fv.f1[1]) > eps )
        || ( fabs(f1[2]-fv.f1[2]) > eps ) ) return(0);

   if (    ( fabs(f2[0]-fv.f2[0]) > eps )
        || ( fabs(f2[1]-fv.f2[1]) > eps )
        || ( fabs(f2[2]-fv.f2[2]) > eps ) ) return(0);

   if (    ( fabs(f3[0]-fv.f3[0]) > eps )
        || ( fabs(f3[1]-fv.f3[1]) > eps )
        || ( fabs(f3[2]-fv.f3[2]) > eps ) ) return(0);

   if (    ( fabs(f4[0]-fv.f4[0]) > eps )
        || ( fabs(f4[1]-fv.f4[1]) > eps )
        || ( fabs(f4[2]-fv.f4[2]) > eps ) ) return(0);

   return(1);
}

#undef __inline__ 

#endif  //  _PDS_SSF_FV_H
