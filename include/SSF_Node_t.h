#pragma once

#ifndef _PDS_SSF_NODE_H
#define _PDS_SSF_NODE_H

#include <stdlib.h>

#include "cuda_portability.h"

#include "Typedefs.h"
#include "Node.h"

//------------------------------------------------------------------------------------------------------------

class SSF_Node_t
{
   public:
      unsigned long    indx      ;   ///< node index
      unsigned long    tag       ;   ///< node tag
               real8   pos  [3]  ;   ///< node position <xyz>
               real8   force[3]  ;   ///< node force    <xyz>

      unsigned long    narms     ;   ///< arm count
      unsigned long   *arm_tags  ;   ///< array of arm node tags
               real8  *arm_burg  ;   ///< array of arm burgers vectors
               real8  *arm_forces;   ///< array of arm forces

   public:
      __cuda_hdev__  SSF_Node_t(void);
      __cuda_host__  SSF_Node_t(const SSF_Node_t & n);
      __cuda_host__ ~SSF_Node_t();

      __cuda_hdev__  const SSF_Node_t & operator =  (const SSF_Node_t & n);
      __cuda_hdev__  const SSF_Node_t & operator =  (const Node_t     & n);

      __cuda_hdev__        int          operator == (const SSF_Node_t & n) const;

      __cuda_hdev__           void     Reset_Forces     (void);
      __cuda_hdev__           int      Arm_Index        (const unsigned long t);
      __cuda_hdev__           void     Add_Arm_Force    (const unsigned long t, const real8 *f);

      __cuda_hdev__           size_t   Serialized_Bytes (void) const;
      __cuda_host__  unsigned char    *Serialize        (const unsigned char *p) const;
};

//------------------------------------------------------------------------------------------------------------

__cuda_host__ size_t SSF_Nodes_Serialized_Bytes (const SSF_Node_t  *nodes, const int  ncnt);
__cuda_host__ size_t SSF_Nodes_Serialized_Bytes (const Node_t      *nodes, const int  ncnt);
__cuda_host__ size_t SSF_Nodes_Serialized_Bytes (      Node_t     **narr , const int  ncnt);

__cuda_host__ void   SSF_Nodes_Serialize        (unsigned char *buf, const SSF_Node_t *nodes, const int ncnt);

__cuda_host__ void   SSF_Nodes_Remap            (unsigned char *cpu,                     const int  ncnt);
__cuda_host__ void   SSF_Nodes_Remap            (unsigned char *cpu, unsigned char *gpu, const int  ncnt);

__cuda_host__
void SSF_Nodes_Send
(
         unsigned char   *cpu  ,   ///< points to a serialized node buffer on host   (cpu)
         unsigned char   *gpu  ,   ///< points to a serialized node buffer on device (gpu)
   const          size_t  bytes,   ///< size of the serialized node buffers
   const          int     ncnt     ///< number of serialized nodes
);

__cuda_host__
void SSF_Nodes_Recv
(
         unsigned char   *cpu  ,   ///< points to a serialized node buffer on host   (cpu)
         unsigned char   *gpu  ,   ///< points to a serialized node buffer on device (gpu)
   const          size_t  bytes,   ///< size of the serialized node buffers
   const          int     ncnt     ///< number of serialized nodes
);

//------------------------------------------------------------------------------------------------------------

#endif  //  _PDS_SSF_NODE_H
