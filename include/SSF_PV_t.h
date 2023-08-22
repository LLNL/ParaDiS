#pragma once

#ifndef _PDS_SSF_PV_H
#define _PDS_SSF_PV_H

#include "cuda_portability.h"

#include "Typedefs.h"
#include "V3.h"

#define __inline__ inline __attribute__((__always_inline__))

//------------------------------------------------------------------------------------------------------------
// SSF_PV_t class support
//
// This class was added as an optimization for the GPU host-device memory transfers.  A seg-seg force 
// pair vector (SSF_PV) includes four node indices that can be used to lookup node positions and
// two burgers vectors (one for each segment).
//------------------------------------------------------------------------------------------------------------

class SSF_PV_t
{
   public:
      unsigned short n1;     ///< index of node 1
      unsigned short n2;     ///< index of node 2
      unsigned short n3;     ///< index of node 3
      unsigned short n4;     ///< index of node 4
               real8 b1[3];  ///< burgers vector (node 1 to node 2) <xyz>
               real8 b3[3];  ///< burgers vector (node 3 to node 4) <xyz>

   public:
      __cuda_hdev__  SSF_PV_t(void);
      __cuda_hdev__  SSF_PV_t(const SSF_PV_t & pv);

      __cuda_hdev__ ~SSF_PV_t() {}

      __cuda_hdev__  const SSF_PV_t & operator =  (const SSF_PV_t & pv);
      __cuda_hdev__        int        operator == (const SSF_PV_t & pv) const;
};

//------------------------------------------------------------------------------------------------------------

__cuda_hdev__  
__inline__
SSF_PV_t::SSF_PV_t(void)
{ 
   n1=0;
   n2=0;
   n3=0;
   n4=0;
   V3_ZERO(b1);
   V3_ZERO(b3); 
}

__cuda_hdev__  
__inline__
SSF_PV_t::SSF_PV_t(const SSF_PV_t & pv) 
{ 
   n1=pv.n1; 
   n2=pv.n2; 
   n3=pv.n3; 
   n4=pv.n4; 
   V3_COPY(b1,pv.b1); 
   V3_COPY(b3,pv.b3); 
}

__cuda_hdev__  
__inline__
const SSF_PV_t & SSF_PV_t::operator = (const SSF_PV_t & pv) 
{ 
   n1=pv.n1; 
   n2=pv.n2; 
   n3=pv.n3; 
   n4=pv.n4; 
   V3_COPY(b1,pv.b1); 
   V3_COPY(b3,pv.b3); 
   return(*this); 
}

__cuda_hdev__
__inline__
int SSF_PV_t::operator == (const SSF_PV_t& pv) const
{
   const real8 eps=1.0e-8; 

   if (    (n1!=pv.n1)
        || (n2!=pv.n2)
        || (n3!=pv.n3)
        || (n4!=pv.n4) ) return(0);

   if (    ( fabs(b1[0]-pv.b1[0]) > eps )
        || ( fabs(b1[1]-pv.b1[1]) > eps )
        || ( fabs(b1[2]-pv.b1[2]) > eps ) ) return(0);

   if (    ( fabs(b3[0]-pv.b3[0]) > eps )
        || ( fabs(b3[1]-pv.b3[1]) > eps )
        || ( fabs(b3[2]-pv.b3[2]) > eps ) ) return(0);

   return(1);
}

#undef __inline__ 

#endif  //  _PDS_SSF_PV_H
