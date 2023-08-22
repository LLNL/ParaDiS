#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "cuda_portability.h"

#include "Typedefs.h"
#include "Home.h"
#include "Node.h"
#include "GPU_Node.h"
#include "V3.h"

//------------------------------------------------------------------------------------------------------------

#define VDELETE(a) if (a) { delete [] a; a=0; }

//------------------------------------------------------------------------------------------------------------

__cuda_host__ 
GPU_Node::GPU_Node(void) 
{ 
   tag        = 0; 
   indx       = 0;
   narms      = 0; 

   V3_ZERO(pos  );
   V3_ZERO(force);

   arm_tags   = 0;
   arm_burg   = 0;
   arm_forces = 0;
   arm_nodes  = 0;
}

//------------------------------------------------------------------------------------------------------------

__cuda_host__ 
GPU_Node::GPU_Node(const GPU_Node & n)
{
   tag        = n.tag  ;
   indx       = n.indx ;
   narms      = n.narms;
  
   V3_COPY(pos  ,n.pos  ); 
   V3_COPY(force,n.force); 

   arm_tags   = 0;
   arm_burg   = 0;
   arm_forces = 0;
   arm_nodes  = 0;

   if (narms>0)
   {
      arm_tags   = new unsigned long[  narms];
      arm_burg   = new real8        [3*narms];
      arm_forces = new real8        [3*narms];
      arm_nodes  = new GPU_Node    *[  narms];

      memcpy(arm_tags  , n.arm_tags  ,   narms*sizeof(unsigned long) );
      memcpy(arm_burg  , n.arm_burg  , 3*narms*sizeof(real8)         );
      memcpy(arm_forces, n.arm_forces, 3*narms*sizeof(real8)         );
      memcpy(arm_nodes , n.arm_nodes ,   narms*sizeof(GPU_Node    *) );
   }
}

//------------------------------------------------------------------------------------------------------------

__cuda_host__ 
GPU_Node::~GPU_Node()
{
   VDELETE(arm_tags  );
   VDELETE(arm_burg  );
   VDELETE(arm_forces);
   VDELETE(arm_nodes ); narms=0; 
}

//------------------------------------------------------------------------------------------------------------

__cuda_host__
const GPU_Node & GPU_Node::operator = (const GPU_Node & n)
{
   if (this!=&n)
   {
      VDELETE(arm_tags  );
      VDELETE(arm_burg  );
      VDELETE(arm_forces);
      VDELETE(arm_nodes ); narms=0;

      tag   = n.tag  ;
      indx  = n.indx ;
      narms = n.narms;

      V3_COPY(pos  ,n.pos  );
      V3_COPY(force,n.force);
      
      if (narms>0)
      {
         arm_tags   = new unsigned long[  narms];
         arm_burg   = new real8        [3*narms];
         arm_forces = new real8        [3*narms];
         arm_nodes  = new GPU_Node    *[  narms];
   
         memcpy(arm_tags  , n.arm_tags  ,   narms*sizeof(unsigned long) );
         memcpy(arm_burg  , n.arm_burg  , 3*narms*sizeof(real8)         );
         memcpy(arm_forces, n.arm_forces, 3*narms*sizeof(real8)         );
         memcpy(arm_nodes , n.arm_nodes ,   narms*sizeof(GPU_Node *)    );
      }
   }

   return(*this);
}

// Copy operator GPU_Node = Node_t
//
// Copies the elements of a source node into a GPU node.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
const GPU_Node & GPU_Node::operator = (const Node_t & n)
{
   VDELETE(arm_tags  );
   VDELETE(arm_burg  );
   VDELETE(arm_forces);
   VDELETE(arm_nodes ); narms=0;

   tag    = (unsigned long) n.myTag;
   indx   = n.myIndex;
   narms  = n.numNbrs;

   pos  [0] = n.x;
   pos  [1] = n.y;
   pos  [2] = n.z;

   force[0] = n.fX;
   force[1] = n.fY;
   force[2] = n.fZ;

   if (narms>0)
   {
      arm_tags   = new unsigned long[  narms];
      arm_burg   = new real8        [3*narms];
      arm_forces = new real8        [3*narms];
      arm_nodes  = new GPU_Node    *[  narms];

      for (unsigned int j=0,k=0; (j<narms); j++, k+=3)
      {
         arm_tags  [j  ] = (unsigned long) n.nbrTag  [j];
         arm_burg  [k  ] = n.burgX[j];
         arm_burg  [k+1] = n.burgY[j];
         arm_burg  [k+2] = n.burgZ[j];
         arm_forces[k  ] = n.armfx[j];
         arm_forces[k+1] = n.armfy[j];
         arm_forces[k+2] = n.armfz[j];
         arm_nodes [j  ] = 0;
      }
   }

   return(*this);
}

//------------------------------------------------------------------------------------------------------------

__cuda_host__
int GPU_Node::operator == (const Node_t & n) const
{
   const real8 eps = 1.0e-6;

   if ( tag   != (unsigned long) n.myTag   ) return(0);
   if ( indx  != (unsigned int ) n.myIndex ) return(0);
   if ( narms != (unsigned int ) n.numNbrs ) return(0);

   if ( fabs(pos[0]-n.x) > eps ) return(0);
   if ( fabs(pos[1]-n.y) > eps ) return(0);
   if ( fabs(pos[2]-n.z) > eps ) return(0);

   if (narms>0)
   {
      for (unsigned int j=0,k=0; (j<narms); j++, k+=3)
      {
         if ( arm_tags[j] != (unsigned long) n.nbrTag  [j] ) return(0);

         if ( fabs(arm_burg[k  ]-n.burgX[j]) > eps ) return(0);
         if ( fabs(arm_burg[k+1]-n.burgY[j]) > eps ) return(0);
         if ( fabs(arm_burg[k+2]-n.burgZ[j]) > eps ) return(0);
      }
   }

   return(1);
}

// Bytes()
//
// Will return the size (in bytes) of the space needed to pack a gpu node. Includes the space
// required for the main data structure, plus the additional space for the dynamic arm info.
//------------------------------------------------------------------------------------------------------------

__cuda_hdev__ 
size_t GPU_Node::Bytes (void) const
{ 
   return(           sizeof(GPU_Node     )      // node struct
           +   narms*sizeof(unsigned long)      // arm_tags
           + 3*narms*sizeof(real8        )      // arm_burg
           + 3*narms*sizeof(real8        )      // arm_forces
           +   narms*sizeof(GPU_Node *   ) );   // arm_nodes
}
   
//------------------------------------------------------------------------------------------------------------

__cuda_hdev__ 
int GPU_Node::Arm_Index (const unsigned long tag) const
{
   if (arm_tags)
   {
      for (unsigned int i=0; (i<narms); i++)
      { if (tag==arm_tags[i]) return(i); }
   }

   return(-1);
}

//------------------------------------------------------------------------------------------------------------

__cuda_hdev__ 
GPU_Node *GPU_Node::Arm_Node (const unsigned long tag) const
{
   if (arm_tags && arm_nodes)
   {
      for (unsigned int i=0; (i<narms); i++)
      { if (tag==arm_tags[i]) return(arm_nodes[i]); }
   }

   return(0);
}

//------------------------------------------------------------------------------------------------------------

__cuda_device__
void GPU_Node::Arm_Force_Accum (const unsigned long tag, const real8 *fv)
{
   int j=3*Arm_Index(tag);

   if (arm_forces && (j>=0))
   {
#ifdef __CUDACC__
      atomicAdd(arm_forces+j  , fv[0]);
      atomicAdd(arm_forces+j+1, fv[1]);
      atomicAdd(arm_forces+j+2, fv[2]);
#else
      arm_forces[j  ] += fv[0];    // (NOT THREAD SAFE!)
      arm_forces[j+1] += fv[1];    // (NOT THREAD SAFE!)
      arm_forces[j+2] += fv[2];    // (NOT THREAD SAFE!)
#endif
   }
}

__cuda_device__  
void GPU_Node::Arm_Force_Accum (const unsigned long tag, const real8 fx, const real8 fy, const real8 fz)
{
   int j=3*Arm_Index(tag);

   if (arm_forces && (j>=0))
   {
#ifdef __CUDACC__
      atomicAdd(arm_forces+j  , fx);
      atomicAdd(arm_forces+j+1, fy);
      atomicAdd(arm_forces+j+2, fz);
#else
      arm_forces[j  ] += fx;       // (NOT THREAD SAFE!)
      arm_forces[j+1] += fy;       // (NOT THREAD SAFE!)
      arm_forces[j+2] += fz;       // (NOT THREAD SAFE!)
#endif
   }
}

//------------------------------------------------------------------------------------------------------------

__cuda_hdev__
void GPU_Node::Reset_Forces (void)
{
   memset(force,0,sizeof(force));

   if (arm_forces && (narms>0))
      memset(arm_forces,0,3*narms*sizeof(real8));
}

__cuda_hdev__
void GPU_Node::Reset_Forces 
(
   const real8 fx,   ///< force (x)
   const real8 fy,   ///< force (y)
   const real8 fz    ///< force (z)
)
{
   memset(force,0,sizeof(force));

   if (arm_forces && (narms>0))
   {
      for (unsigned int i=0,k=0; (i<narms); i++,k+=3)
      {
         arm_forces[k  ] = fx;
         arm_forces[k+1] = fy;
         arm_forces[k+2] = fz;
      }
   }
}

// Serialize()
//
// Will serialize the current instance into a packing buffer. 
//------------------------------------------------------------------------------------------------------------

__cuda_host__
unsigned char *GPU_Node::Serialize (unsigned char *p) const
{
   if (p)
   {
      GPU_Node *n = (GPU_Node *) p;   p += sizeof(GPU_Node);

      memcpy(n,this,sizeof(GPU_Node));

      if (narms>0)
      {
         n->arm_tags   = (unsigned long *) p;  p +=   narms*sizeof(unsigned long);
         n->arm_burg   = (real8         *) p;  p += 3*narms*sizeof(real8        );
         n->arm_forces = (real8         *) p;  p += 3*narms*sizeof(real8        );
         n->arm_nodes  = (GPU_Node     **) p;  p +=   narms*sizeof(GPU_Node *   );

         memcpy(n->arm_tags  , arm_tags  ,   narms*sizeof(unsigned long) );
         memcpy(n->arm_burg  , arm_burg  , 3*narms*sizeof(real8)         );
         memcpy(n->arm_forces, arm_forces, 3*narms*sizeof(real8)         );
         memcpy(n->arm_nodes , arm_nodes ,   narms*sizeof(GPU_Node *)    );
      }
   }

   return(p);
}

// Remap()
//
// Before copying a buffer containing a collection of serialized gpu nodes, the
// pointers within the serialized nodes need to be relative to the base address
// of the device memory.  
//------------------------------------------------------------------------------------------------------------

__cuda_host__  
size_t GPU_Node::Remap 
(
   unsigned char *pc ,  ///< points to a serialized gpu node in host   memory
   unsigned char *pg ,  ///< points to a serialized gpu node in device memory
   unsigned char *pc0,  ///< points to the base address of the cpu packing buffer
   unsigned char *pg0   ///< points to the base address of the gpu packing buffer
) const
{
   unsigned char *p0 = pc;

   if (pc && pg && pc0 && pg0)
   {
      GPU_Node *n = (GPU_Node *) pc;   pc += sizeof(GPU_Node); 
                                       pg += sizeof(GPU_Node);

      if (n->narms>0)
      {
         GPU_Node **nptrs = n->arm_nodes;  // (save)

         size_t bytes =    narms*sizeof(unsigned long);  n->arm_tags   = (unsigned long *) pg;  pg += bytes; pc += bytes;
                bytes =  3*narms*sizeof(real8        );  n->arm_burg   = (real8         *) pg;  pg += bytes; pc += bytes;
                bytes =  3*narms*sizeof(real8        );  n->arm_forces = (real8         *) pg;  pg += bytes; pc += bytes;
                bytes =    narms*sizeof(GPU_Node *   );  n->arm_nodes  = (GPU_Node     **) pg;  pg += bytes; pc += bytes;

         for (int i=0; (i<narms); i++)
            nptrs[i] = (GPU_Node *) (pg0+(((unsigned char *) nptrs[i]) - pc0));
      }
   }

   return((size_t)(pc-p0));
}

//------------------------------------------------------------------------------------------------------------

__cuda_host__ 
void GPU_Node::Print (FILE *fd) const
{
   if (fd)
   {
      char tstr[32];   sprintf(tstr,"(%d,%d)", (int) (tag>>32), (int) (tag & 0xffffffff) );

      fprintf(fd,"%5lu %-10s %10.4lf %10.4lf %10.4lf  %4.1lf %4.1lf %4.1lf\n", indx, tstr, pos[0], pos[1], pos[2], force[0], force[1], force[2] );

      for (unsigned int i=0,j=0; (i<narms); i++,j+=3)
      {
         sprintf(tstr,"(%d,%d)", (int) (arm_tags[i]>>32), (int) (arm_tags[i] & 0xffffffff) );
         
         fprintf(fd,"   %-10s %10.6lf %10.6lf %10.6lf %4.1lf %4.1lf %4.1lf\n", 
                 tstr, 
                 arm_burg  [j], arm_burg  [j+1], arm_burg  [j+2], 
                 arm_forces[j], arm_forces[j+1], arm_forces[j+2] );
      }
   }
}

// GPU_Nodes_Bytes()
//
// Given an array of node pointers, will return the total space required to pack those nodes
// into a buffer that can be sent and used on the GPU. Note that the space includes an array
// of node pointers at the head of the serialized result.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
size_t GPU_Nodes_Bytes
(
   Node_t *nodes, ///< array of nodes
   int     ncnt   ///< number of nodes
)
{
   size_t bytes=0;

   if (nodes && (ncnt>0))
   {
      int narms=0;
      for (int i=0; (i<ncnt); i++) { narms += nodes[i].numNbrs; }

      bytes  =   ncnt *sizeof(GPU_Node *   );   // array of node pointers
      bytes +=   ncnt *sizeof(GPU_Node     );   // array of node structs
      bytes +=   narms*sizeof(unsigned long);   // arm_tags
      bytes += 3*narms*sizeof(real8        );   // arm_burg
      bytes += 3*narms*sizeof(real8        );   // arm_forces
      bytes +=   narms*sizeof(GPU_Node *   );   // arm_nodes
   }

   return(bytes);
}

__cuda_host__
size_t GPU_Nodes_Bytes
(
   Node_t **narr, ///< array of node pointers 
   int      ncnt  ///< number of nodes
)
{
   size_t bytes=0;

   if (narr && (ncnt>0))
   {
      int narms=0;
      for (int i=0; (i<ncnt); i++) { narms += ( narr[i] ? narr[i]->numNbrs : 0 ); }

      bytes  =   ncnt *sizeof(GPU_Node *   );   // array of node pointers
      bytes +=   ncnt *sizeof(GPU_Node     );   // array of node structs
      bytes +=   narms*sizeof(unsigned long);   // arm_tags
      bytes += 3*narms*sizeof(real8        );   // arm_burg
      bytes += 3*narms*sizeof(real8        );   // arm_forces
      bytes +=   narms*sizeof(GPU_Node *   );   // arm_nodes
   }

   return(bytes);
}

//------------------------------------------------------------------------------------------------------------

__cuda_device__
void GPU_Node_Add_Forces
(
         GPU_Node  **nodes,  ///< array of packed node pointers
   const int         i1   ,  ///< index of node 1
   const int         i2   ,  ///< index of node 2
   const real8      *f1   ,  ///< points to force contribution on node 1 <xyz>
   const real8      *f2      ///< points to force contribution on node 2 <xyz>
)
{
   GPU_Node      *n1 = nodes[i1];
   GPU_Node      *n2 = nodes[i2];

   unsigned long  t1 = (unsigned long) n1->tag;
   unsigned long  t2 = (unsigned long) n2->tag;

   n1->Arm_Force_Accum(t2,f1);
   n2->Arm_Force_Accum(t1,f2);
}

// GPU_Node_Array()
//
// Given a buffer of serialized gpu nodes, will populate and return an array 
// of gpu node pointers within the buffer.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void GPU_Node_Array (GPU_Node **narr, unsigned char *p, const int pn)
{
   if (narr && p && (pn>0))
   {
      for (int i=0; (i<pn); ++i)
      { narr[i] = (GPU_Node *) p;  p += narr[i]->Bytes(); }
   }
}

__cuda_host__
GPU_Node **GPU_Node_Array (unsigned char *p, const int pn)
{
   GPU_Node **narr = ( (p && (pn>0)) ? new GPU_Node *[pn] : 0 );

   if (narr && p && (pn>0))
   {
      for (int i=0; (i<pn); ++i)
      { narr[i] = (GPU_Node *) p;  p += narr[i]->Bytes(); }
   }

   return(narr);
}

// GPU_Nodes_Sort()
//
// Will sort a dense array of node pointers by tag.
//------------------------------------------------------------------------------------------------------------

static int ncmp (const void *p0, const void *p1)
{
   GPU_Node *n0 = ( p0 ? *((GPU_Node **) p0) : 0 );
   GPU_Node *n1 = ( p1 ? *((GPU_Node **) p1) : 0 );

   if (n0 && n1)
   {
      if (n0->tag < n1->tag) return(-1);  // (n0<n1)
      if (n0->tag > n1->tag) return(+1);  // (n0>n1)
   }

   return(0); // (fall through n0==n1)
}

__cuda_host__
GPU_Node **GPU_Nodes_Sort 
(
   GPU_Node **nodes,
   const int  ncnt
)
{
   if (nodes && (ncnt>0)) 
   { qsort( (void *) nodes, ncnt, sizeof(GPU_Node *), ncmp); }

   return(nodes);
}

// GPU_Nodes_Find()
//
// Given a dense array of sorted node pointers, will perform a binary search for the given tag.
//---------------------------------------------------------------------------------------------

__cuda_host__
GPU_Node *GPU_Nodes_Find 
(
   GPU_Node      **narr,  ///< sorted array of GPU node pointers
   const int       ncnt,  ///< number of GPU nodes
   unsigned long   tag    ///< tag to find
)
{
   if (narr && (ncnt>0))
   {
      int i0=0;
      int i1=(ncnt-1);
   
      while (i0<=i1)
      {   
         int i = (i0+i1)/2;

         GPU_Node *p = narr[i];

         int  cmp = ( (p->tag < tag) ? -1 : 
                      (p->tag > tag) ? +1 : 0 );

         if (cmp==0) { return(p); }
         if (cmp< 0) { i0=i+1; }
         else        { i1=i-1; }
      }   
   }

   return(0);
}

// GPU_Nodes_Serialize()
//
// Will serialize a block of GPU nodes into a packing buffer suitable for sending to the GPU.
// Note that the serialized nodes includes an array of sorted node pointers preceeding the 
// serialized nodes.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
GPU_Node **GPU_Nodes_Serialize 
(
   unsigned char *pc   ,  ///< points to a buffer to receive the serialized nodes
   const Node_t  *nodes,  ///< points to the source nodes
   const int      ncnt    ///< number of nodes to serialize
)
{
   // reserve space at the head of the buffer for an array of gpu node pointers.

   GPU_Node **narr = (GPU_Node **) pc;

   if (pc && nodes && (ncnt>0))
   {
      pc += ncnt*sizeof(GPU_Node *);  // skip to first serialized node

      for (int i=0; (i<ncnt); i++)
      {
         GPU_Node gpu_node;  gpu_node=nodes[i];

         narr[i] = (GPU_Node *) pc;

         pc = gpu_node.Serialize(pc);
      }
   } 

   // once the nodes are serialized in memory, we need to wire up
   // the node pointers relative to the packed buffer.

   // copy the node array (we need them sorted for the lookups)

   GPU_Node **narr_tmp = new GPU_Node *[ncnt];

   if (narr_tmp && narr) 
   { 
      memcpy(narr_tmp,narr, ncnt*sizeof(GPU_Node *) ); 

      GPU_Nodes_Sort(narr,ncnt); 

      for (int i=0; (i<ncnt); i++)
      {
         GPU_Node *pn    =   narr[i];
         int       narms = ( pn ? pn->narms : 0 );

         if (pn && pn->arm_tags && pn->arm_nodes && (narms>0) )
         {
            for (int j=0; (j<narms); ++j)
            { pn->arm_nodes[j] = GPU_Nodes_Find(narr,ncnt,pn->arm_tags[j]); }
         } 
      }
   }  

   VDELETE(narr_tmp);

   return(narr);
}

// GPU_Nodes_Remap()
//
// After an array of nodes are serialized into a buffer, all the node pointers will be
// relative to their location within the buffer on the host.  In order to use the serialized
// nodes on the GPU, we need to reset those pointers relative the the buffer on the device.
//---------------------------------------------------------------------------------------------

__cuda_host__
void GPU_Nodes_Remap 
(
   unsigned char *gbuf ,  ///< points to the base pointer on the gpu
   unsigned char *pbuf ,  ///< points to the base pointer on the host 
   const int      ncnt    ///< number of nodes to remap
)
{
   if (gbuf && pbuf && (ncnt>0))
   {
      GPU_Node **narr = (GPU_Node **) pbuf;

      for (int i=0; (i<ncnt); i++)
      {
         GPU_Node *node = narr[i];

         if (node)
         {
            node->arm_tags   = (unsigned long *) ( gbuf + (((unsigned char *) node->arm_tags    ) - pbuf) );
            node->arm_burg   = (real8         *) ( gbuf + (((unsigned char *) node->arm_burg    ) - pbuf) );
            node->arm_forces = (real8         *) ( gbuf + (((unsigned char *) node->arm_forces  ) - pbuf) );

            int narms = narr[i]->narms;

            for (int j=0; (j<narms); j++)
               node->arm_nodes[j] = (GPU_Node *) ( gbuf + (((unsigned char *) node->arm_nodes[j]) - pbuf) );

            node->arm_nodes = (GPU_Node **)      ( gbuf + (((unsigned char *) node->arm_nodes   ) - pbuf) );

            narr[i]         = (GPU_Node  *)      ( gbuf + (((unsigned char *) narr[i]           ) - pbuf) );
         }
      }
   }
}

#ifdef GPU_ENABLED
// GPU_Nodes_Reset_K
//
// This kernel function will reset the force vectors within an array of GPU node pointers.
// This was mainly created to test the GPU node unpacking routine, but may be generally
// useful if we need to reset the forces within a step.
//---------------------------------------------------------------------------------------------

__cuda_global__
void GPU_Nodes_Reset_K
(
         GPU_Node **nodes, ///< array of node pointers (on device)
   const int        ncnt , ///< number of nodes
   const real8      fx   , ///< force (x) (default=0.0)
   const real8      fy   , ///< force (y) (default=0.0)
   const real8      fz     ///< force (z) (default=0.0)
)
{
   int bw = (blockDim .x);   // bw  = cuda block  width
   int bx = (blockIdx .x);   // bx  = cuda block  index
   int tx = (threadIdx.x);   // tx  = cuda thread index
   int i  = bx*bw+tx;        // i   = thread index for this kernel

   GPU_Node *node = ( (nodes && (i<ncnt)) ? nodes[i] : 0 );

   if (node) { node->Reset_Forces(fx,fy,fz); }
}

__cuda_host__
void GPU_Nodes_Reset
(
         GPU_Node **nodes, ///< array of node pointers (on device)
   const int        ncnt , ///< number of nodes
   const real8      fx   , ///< force (x) (default=0.0)
   const real8      fy   , ///< force (y) (default=0.0)
   const real8      fz     ///< force (z) (default=0.0)
)
{
   if (nodes && (ncnt>0))
   {
      int nthrd = 64;
      int nblk  = (ncnt+(nthrd-1))/nthrd;

      GPU_Nodes_Reset_K<<<nblk,nthrd>>>(nodes,ncnt, fx,fy,fz);
      CUDART_CHECK(cudaGetLastError());
   }
}

#endif
