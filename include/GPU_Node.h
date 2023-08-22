#pragma once

#ifndef _PDS_GPU_NODE_H
#define _PDS_GPU_NODE_H

#include <stdlib.h>

#include "cuda_portability.h"

#include "Typedefs.h"
#include "Home.h"
#include "Node.h"

//------------------------------------------------------------------------------------------------------------

class GPU_Node
{
   public :
      unsigned long        tag       ;   ///< node tag
      unsigned long        indx      ;   ///< node index
      unsigned long        narms     ;   ///< arm count
   
               real8       pos  [3]  ;   ///< node position <xyz>
               real8       force[3]  ;   ///< node force    <xyz>
   
      unsigned long       *arm_tags  ;   ///< array of arm tags (packed)
               real8      *arm_burg  ;   ///< array of arm burgers vectors
               real8      *arm_forces;   ///< array of arm forces
               GPU_Node  **arm_nodes ;   ///< array of arm node pointers

   public :
      __cuda_host__    GPU_Node(void);
      __cuda_host__    GPU_Node(const GPU_Node & n);
      __cuda_host__   ~GPU_Node();

      __cuda_host__    const GPU_Node & operator =      (const GPU_Node & n);
      __cuda_host__    const GPU_Node & operator =      (const Node_t   & n);

      __cuda_host__    int              operator ==     (const Node_t   & n) const;

      __cuda_hdev__    size_t           Bytes           (void) const;
      __cuda_hdev__    int              Arm_Index       (const unsigned long tag) const;
      __cuda_hdev__    GPU_Node        *Arm_Node        (const unsigned long tag) const;
      __cuda_device__  void             Arm_Force_Accum (const unsigned long tag, const real8 fx, const real8 fy, const real8 fz);
      __cuda_device__  void             Arm_Force_Accum (const unsigned long tag, const real8 *fv);

      __cuda_hdev__    void             Reset_Forces    (void);
      __cuda_hdev__    void             Reset_Forces    (const real8 fx, const real8 fy, const real8 fz);

      __cuda_host__    unsigned char   *Serialize       (unsigned char *p) const;
      __cuda_host__    size_t           Remap           (unsigned char *pc , unsigned char *pg, unsigned char *pc0, unsigned char *pg0 ) const;
   
      __cuda_host__    void             Print           (FILE *fd=stdout) const;
};

//------------------------------------------------------------------------------------------------------------

__cuda_host__
size_t GPU_Nodes_Bytes
(
   Node_t  *nodes,  ///< array of nodes
   int      ncnt    ///< number of nodes
);

__cuda_host__
size_t GPU_Nodes_Bytes
(
   Node_t **narr ,  ///< dense array of node pointers
   int      ncnt    ///< number of nodes
);

//------------------------------------------------------------------------------------------------------------

__cuda_host__
void GPU_Node_Array 
(
   GPU_Node     **narr,   ///< array to receive the gpu node pointers (host)
   unsigned char *p   ,   ///< points to a buffer containing serialized gpu nodes (host)
            int   pn      ///< number of serialized nodes
);

__cuda_host__
GPU_Node **GPU_Node_Array 
(
   unsigned char *p   ,   ///< points to a buffer containing serialized gpu nodes (host)
            int   pn      ///< number of serialized nodes
);

__cuda_host__
GPU_Node **GPU_Nodes_Sort 
(
   GPU_Node **nodes,      ///< points to a dense array of GPU node pointers
   const int  ncnt        ///< number of nodes
);

__cuda_host__
GPU_Node *GPU_Nodes_Find 
(
   GPU_Node      **nodes,  ///< sorted array of GPU node pointers
   const int       ncnt ,  ///< number of GPU nodes
   unsigned long   ntag    ///< tag to find
);

__cuda_host__
GPU_Node **GPU_Nodes_Serialize 
(
   unsigned char *pbuf ,  ///< points to a buffer to receive the serialized nodes (host)
   const Node_t  *nodes,  ///< points to the source nodes
   const int      ncnt    ///< number of nodes to serialize
);

__cuda_host__
void GPU_Nodes_Remap 
(
   unsigned char *gbuf ,  ///< points to the base pointer on the gpu
   unsigned char *pbuf ,  ///< points to the base pointer on the host 
   const int      ncnt    ///< number of nodes to remap
);

//------------------------------------------------------------------------------------------------------------

__cuda_host__
void GPU_Nodes_Reset
(
         GPU_Node **narr     , ///< dense array of gpu node pointers (gpu)
   const int        ncnt     , ///< number of nodes
   const real8      fx  = 0.0, ///< force (x)
   const real8      fy  = 0.0, ///< force (y)
   const real8      fz  = 0.0  ///< force (z)
);

//------------------------------------------------------------------------------------------------------------

#endif   // ifndef _PDS_GPU_NODE_H
