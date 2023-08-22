#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "mpi_portability.h"
#include "cuda_portability.h"

#include "Typedefs.h"
#include "SSF_Node_t.h"

//------------------------------------------------------------------------------------------------------------
// SSF_Node_t class support
//------------------------------------------------------------------------------------------------------------

// void constructor...

__cuda_hdev__
SSF_Node_t::SSF_Node_t(void)
{
   indx       = 0;
   tag        = 0;
   pos  [0]   = pos  [1] = pos  [2] = 0.0;
   force[0]   = force[1] = force[2] = 0.0;
   narms      = 0;
   arm_tags   = 0;
   arm_burg   = 0;
   arm_forces = 0;
}

// copy constructor...

__cuda_host__
SSF_Node_t::SSF_Node_t(const SSF_Node_t & n)
{
   indx       = n.indx;
   tag        = n.tag ;

   pos  [0]   = n.pos[0];
   pos  [1]   = n.pos[1];
   pos  [2]   = n.pos[2];

   force[0]   = n.force[0];
   force[1]   = n.force[1];
   force[2]   = n.force[2];

   narms      = n.narms;

   arm_tags   = ( (narms>0) ? new unsigned long [  narms] : 0 );
   arm_burg   = ( (narms>0) ? new          real8[3*narms] : 0 );
   arm_forces = ( (narms>0) ? new          real8[3*narms] : 0 );

   for (int i=0,j=0; (i<narms); i++,j+=3)
   {
      arm_tags  [i  ] = n.arm_tags  [i  ];

      arm_burg  [j  ] = n.arm_burg  [j  ];
      arm_burg  [j+1] = n.arm_burg  [j+1];
      arm_burg  [j+2] = n.arm_burg  [j+2];

      arm_forces[j  ] = n.arm_forces[j  ];
      arm_forces[j+1] = n.arm_forces[j+1];
      arm_forces[j+2] = n.arm_forces[j+2];
   }
}

// destructor...

__cuda_host__
SSF_Node_t::~SSF_Node_t()
{
   narms=0;
   if (arm_tags  ) { delete [] arm_tags  ; arm_tags  =0; }
   if (arm_burg  ) { delete [] arm_burg  ; arm_burg  =0; }
   if (arm_forces) { delete [] arm_forces; arm_forces=0; }
}

// LHS assignment from SSF_Node_t...

__cuda_hdev__
const SSF_Node_t & SSF_Node_t::operator = (const SSF_Node_t &n)
{
   if (this!=&n)
   {
      narms=0;
      if (arm_tags  ) { delete [] arm_tags  ; arm_tags  =0; }
      if (arm_burg  ) { delete [] arm_burg  ; arm_burg  =0; }
      if (arm_forces) { delete [] arm_forces; arm_forces=0; }

      indx       = n.indx;
      tag        = n.tag;

      pos  [0]   = n.pos[0];
      pos  [1]   = n.pos[1];
      pos  [2]   = n.pos[2];

      force[0]   = n.force[0];
      force[1]   = n.force[1];
      force[2]   = n.force[2];

      narms      = n.narms;

      arm_tags   = ( (narms>0) ? new unsigned long [  narms] : 0 );
      arm_burg   = ( (narms>0) ? new          real8[3*narms] : 0 );
      arm_forces = ( (narms>0) ? new          real8[3*narms] : 0 );

      for (int i=0,j=0; (i<narms); i++,j+=3)
      {
         arm_tags  [i  ] = n.arm_tags  [i  ];

         arm_burg  [j  ] = n.arm_burg  [j  ];
         arm_burg  [j+1] = n.arm_burg  [j+1];
         arm_burg  [j+2] = n.arm_burg  [j+2];

         arm_forces[j  ] = n.arm_forces[j  ];
         arm_forces[j+1] = n.arm_forces[j+1];
         arm_forces[j+2] = n.arm_forces[j+2];
      }
   }

   return(*this);
}

// LHS assignment from Node_t class...

__cuda_hdev__
const SSF_Node_t & SSF_Node_t::operator = (const Node_t & n)
{
   narms=0;
   if (arm_tags  ) { delete [] arm_tags  ; arm_tags  =0; }
   if (arm_burg  ) { delete [] arm_burg  ; arm_burg  =0; }
   if (arm_forces) { delete [] arm_forces; arm_forces=0; }

   indx       = n.myIndex;
   tag        = (unsigned long) n.myTag;

   pos  [0]   = n.x;
   pos  [1]   = n.y;
   pos  [2]   = n.z;

   force[0]   = n.fX;
   force[1]   = n.fY;
   force[2]   = n.fZ;

   narms      = n.numNbrs;

   arm_tags   = ( (narms>0) ? new unsigned long [  narms] : 0 );
   arm_burg   = ( (narms>0) ? new          real8[3*narms] : 0 );
   arm_forces = ( (narms>0) ? new          real8[3*narms] : 0 );

   for (int i=0,j=0; (i<narms); i++,j+=3)
   {
      arm_tags  [i  ] = (unsigned long) n.nbrTag[i];

      arm_burg  [j  ] = n.burgX[i];
      arm_burg  [j+1] = n.burgY[i];
      arm_burg  [j+2] = n.burgZ[i];

      arm_forces[j  ] = n.armfx[i];
      arm_forces[j+1] = n.armfy[i];
      arm_forces[j+2] = n.armfz[i];
   }

   return(*this);
}

// comparison with another node...

__cuda_hdev__
int  SSF_Node_t::operator == (const SSF_Node_t & n) const
{
   if ( indx != n.indx ) return(0);
   if ( tag  != n.tag  ) return(0);

   if (    ( fabs(pos  [0]-n.pos  [0])>1.0e-4 )
        || ( fabs(pos  [1]-n.pos  [1])>1.0e-4 )
        || ( fabs(pos  [2]-n.pos  [2])>1.0e-4 ) ) return(0);

   if (    ( fabs(force[0]-n.force[0])>1.0e-6 )
        || ( fabs(force[1]-n.force[1])>1.0e-6 )
        || ( fabs(force[2]-n.force[2])>1.0e-6 ) ) return(0);

   if ( narms != n.narms ) return(0);

   if ( arm_tags && n.arm_tags )
   {
      for (int i=0; (i<narms); i++)
      { if ( arm_tags[i] != n.arm_tags[i] ) return(0); }
   }

   if ( arm_burg && n.arm_burg )
   {
      for (int i=0; (i<3*narms); i++)
      { if ( fabs(arm_burg[i]-n.arm_burg[i])>1.0e-6 ) return(0); }
   }

   if ( arm_forces && n.arm_forces )
   {
      for (int i=0; (i<3*narms); i++)
      { if ( fabs(arm_forces[i]-n.arm_forces[i])>1.0e-6 ) return(0); }
   }

   return(1);  // (fall through)
}

// Reset_Forces()
//
// Will reset all the forces within an instance of a node.
//------------------------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Node_t::Reset_Forces (void)
{
   force[0]=force[1]=force[2]=0.0;

   if (arm_forces && (narms>0))
   {
      for (int i=0; (i<3*narms); i++)
      { arm_forces[i]=0.0; }
   }
}

// Arm_Index()
//
// Returns the index of the arm associated with a given tag.
// Returns -1 if the tag is not found.
//------------------------------------------------------------------------------------------------------------

__cuda_hdev__
int SSF_Node_t::Arm_Index (const unsigned long t)
{
   if (arm_tags && (narms>0))
   {
      for (int i=0; (i<narms); i++)
      { if (t==arm_tags[i]) return(i); }
   }

   return(-1);  // (not found)
}

// Add_Arm_Force()
//
// Adds the force to the arm associated with a given tag.
//------------------------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Node_t::Add_Arm_Force (const unsigned long t, const real8 *f)
{
   int j = 3*Arm_Index(t);

   if (j>=0)
   {
#ifdef __CUDA_ARCH__
      atomicAdd( arm_forces+j  , f[0] );
      atomicAdd( arm_forces+j+1, f[1] );
      atomicAdd( arm_forces+j+2, f[2] );
#else
      arm_forces[j  ] += f[0];  ///  NOT THREAD SAFE!
      arm_forces[j+1] += f[1];  ///  NOT THREAD SAFE!
      arm_forces[j+2] += f[2];  ///  NOT THREAD SAFE!
#endif
   }
}

//------------------------------------------------------------------------------------------------------------
// Node Serialization
//
// When exchanging data between the host and gpu, we need to serialize the elements of each node into
// a buffer that can be copied between the host and gpu.  Given the arm count for a particular node
// is variable, the serialized size of a partocular node is dependent on the arm count.
//
// When serialized, the node structure containing the non-variable elements are copied to the buffer
// and immediately followed by the serialized arm data...
//
// +------------------------------------------------------+
// | node structure                 sizeof(SSF_Node_t)    |
// + - - - - - - - - - - - - - - - - - - - - - - - - - - -+
// | node arm tags                n*sizeof(unsigned long) |
// | node arm burgers vectors   3*n*sizeof(real8)         |
// | node arm forces            3*n*sizeof(real8)         |
// +------------------------------------------------------+
//
// Note that the internal pointers within the serialized node structure are with respect to the
// memory space provided when serialized.  This means the pointers may be with respect to either
// the host or global memory on the device.
//------------------------------------------------------------------------------------------------------------

// Serialized_Bytes()
//
// Returns the space required to serialize a node and arms.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
size_t SSF_Node_t::Serialized_Bytes (void) const
{
   size_t bytes  =         sizeof(SSF_Node_t   );   // reserve space for structure
          bytes +=   narms*sizeof(unsigned long);   // reserve space for arm tags
          bytes += 3*narms*sizeof(real8        );   // reserve space for burgers vectors
          bytes += 3*narms*sizeof(real8        );   // reserve space for forces

   return(bytes);
};

// Serialize()
//
// Serializes the elements of a node into a buffer. Note this is a host-only method and the memory
// is assumed to be allocated on the host.
//
// Returns the next available position to serialize.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
unsigned char *SSF_Node_t::Serialize (const unsigned char *buf) const
{
   SSF_Node_t    *p  = (SSF_Node_t    *) buf;
   unsigned char *pc = (unsigned char *) buf;

   if (buf)
   {
      size_t bytes = sizeof(SSF_Node_t); memcpy(pc, (unsigned char *) this, bytes); pc+=bytes;

      if (narms>0)
      {
         bytes =   narms*sizeof(unsigned long ); memcpy(pc, arm_tags  , bytes); p->arm_tags   = (unsigned long  *) pc; pc+=bytes;
         bytes = 3*narms*sizeof(         real8); memcpy(pc, arm_burg  , bytes); p->arm_burg   = (         real8 *) pc; pc+=bytes;
         bytes = 3*narms*sizeof(         real8); memcpy(pc, arm_forces, bytes); p->arm_forces = (         real8 *) pc; pc+=bytes;
      }
   }

   return(pc);  // (return next position)
}

// SSF_Node_Serialized_Bytes()
//
// Returns the space needed to serialize an array of nodes into a buffer.
// Several APIs provided.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
size_t SSF_Nodes_Serialized_Bytes
(
   const SSF_Node_t *nodes,    ///< array  of nodes
   const int         ncnt      ///< number of nodes
)
{
   size_t bytes=0;

   if (nodes && (ncnt>0))
   {
      bytes = ncnt*sizeof(SSF_Node_t *);          // reserve space for node pointers

      for (int i=0; (i<ncnt); i++)
      { bytes += nodes[i].Serialized_Bytes(); }   // ass space for each node
   }

   return(bytes);
}

// SSF_Nodes_Serialized_Bytes()
//
// This version returns the space required based on an array of general node pointers.
// The space computed is based on the actual number of nodes and arms within the non-null
// nodes pointed to in the array (e.g. the array can be sparse).
//------------------------------------------------------------------------------------------------------------

__cuda_host__
size_t SSF_Nodes_Serialized_Bytes
(
   const Node_t *nodes,    ///< array  of nodes
   const int     ncnt      ///< number of nodes
)
{
   int acnt = 0;  // ncnt = total arm  count

   if (nodes)
   {
      for (int i=0; (i<ncnt); i++)
      { acnt+=nodes[i].numNbrs; }
   }

   // reserve space for an array of node pointers and serialized nodes...

   size_t bytes  =   ncnt*sizeof(SSF_Node_t  *);   // reserve space for node pointers
          bytes +=   ncnt*sizeof(SSF_Node_t   );   // reserve space for node structures
          bytes +=   acnt*sizeof(unsigned long);   // reserve space for arm tags
          bytes += 3*acnt*sizeof(real8        );   // reserve space for arm burgers vectors
          bytes += 3*acnt*sizeof(real8        );   // reserve space for arm forces

   return(bytes);
}

__cuda_host__
size_t SSF_Nodes_Serialized_Bytes
(
         Node_t **narr     ,    ///< array  of node pointers (may be sparse)
   const int      narr_cnt      ///< number of node pointers
)
{
   int ncnt = 0;  // ncnt = total node count
   int acnt = 0;  // ncnt = total arm  count

   if (narr)
   {
      for (int i=0; (i<narr_cnt); i++)
      { if (narr[i]) { ncnt++; acnt+=narr[i]->numNbrs; } }
   }

   // reserve space for an array of node pointers and serialized nodes...

   size_t bytes  =   ncnt*sizeof(SSF_Node_t  *);   // reserve space for node pointers
          bytes +=   ncnt*sizeof(SSF_Node_t   );   // reserve space for node structures
          bytes +=   acnt*sizeof(unsigned long);   // reserve space for arm tags
          bytes += 3*acnt*sizeof(real8        );   // reserve space for arm burgers vectors
          bytes += 3*acnt*sizeof(real8        );   // reserve space for arm forces

   return(bytes);
}

// SSF_Nodes_Serialize()
//
// Will serialize an array of nodes into a serialization buffer.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Nodes_Serialize
(
         unsigned char *buf  ,    ///< buffer to receive the nodes
   const SSF_Node_t    *nodes,    ///< array of nodes
   const int            ncnt      ///< number of nodes
)
{
   if (buf && nodes && (ncnt>0) )
   {
      unsigned char  *p  = buf;
      SSF_Node_t    **pn = (SSF_Node_t **) p;  p += ncnt*sizeof(SSF_Node_t *);  // reserve space for node pointers

      for (int i=0; (i<ncnt); i++)
      {
         pn[i] = (SSF_Node_t *) p;

         p = nodes[i].Serialize(p);
      }
   }
}

// SSF_Nodes_Remap
//
// Serialized nodes can be moved from host and gpu memory spaces fairly seamlessly.
// The key to being able to move serialized nodes between host and gpu is manipulating the pointers
// within the structures within the serialized node buffers to point to valid memory locations
// either on the host or within the device.  The following two routines will remap those pointers
// as necessary for each memory space.
//
// Note - those are host-only routines.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Nodes_Remap
(
         unsigned char *cpu ,   ///< points to a buffer containing serialized nodes (host)
   const          int   ncnt    ///< number of serialized nodes
)
{
   if (cpu && (ncnt>0))
   {
      unsigned char *pc = (unsigned char *) cpu;

      SSF_Node_t **nodes = (SSF_Node_t **) pc;  pc+=ncnt*sizeof(SSF_Node_t *);

      for (int i=0; (i<ncnt); i++)
      {
         SSF_Node_t *p = (SSF_Node_t *) pc;  pc+=sizeof(SSF_Node_t);

         nodes[i] = p;

         int narms = p->narms;

         p->arm_tags   = ( (narms>0) ? (unsigned long  *) pc : 0 ); pc +=   narms*sizeof(unsigned long);
         p->arm_burg   = ( (narms>0) ? (         real8 *) pc : 0 ); pc += 3*narms*sizeof(real8);
         p->arm_forces = ( (narms>0) ? (         real8 *) pc : 0 ); pc += 3*narms*sizeof(real8);
      }
   }
}

// SSF_Nodes_Remap
//
// Remaps the serialize node pointers for deployment on the device (gpu).
//------------------------------------------------------------------------------------------------------------
__cuda_host__
void SSF_Nodes_Remap
(
         unsigned char *cpu ,   ///< points to the serialized node buffer on host   (cpu)
         unsigned char *gpu ,   ///< points to the serialized node buffer on device (gpu)
   const          int   ncnt    ///< number of serialized nodes
)
{
   if (cpu && gpu && (ncnt>0))
   {
      unsigned char *pc = (unsigned char *) cpu;
      unsigned char *pg = (unsigned char *) gpu;

      SSF_Node_t **nodes = (SSF_Node_t **) pc;

      size_t bytes = ncnt*sizeof(SSF_Node_t *); pc+=bytes; pg+=bytes;

      for (int i=0; (i<ncnt); i++)
      {
         nodes[i] = (SSF_Node_t *) pg;

         SSF_Node_t *p = (SSF_Node_t *) pc;  pc+=sizeof(SSF_Node_t);
                                             pg+=sizeof(SSF_Node_t);

         int narms = p->narms;

         p->arm_tags   = ( (narms>0) ? (unsigned long  *) pg : 0 ); bytes =   narms*sizeof(unsigned long);  pc+=bytes; pg+=bytes;
         p->arm_burg   = ( (narms>0) ? (         real8 *) pg : 0 ); bytes = 3*narms*sizeof(real8);          pc+=bytes; pg+=bytes;
         p->arm_forces = ( (narms>0) ? (         real8 *) pg : 0 ); bytes = 3*narms*sizeof(real8);          pc+=bytes; pg+=bytes;
      }
   }
}

// SSF_Send_Nodes()
//
// Will remap the serialized nodes and send them to the device (gpu).
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Nodes_Send
(
         unsigned char   *cpu  ,   ///< points to a serialized node buffer on host   (cpu)
         unsigned char   *gpu  ,   ///< points to a serialized node buffer on device (gpu)
   const          size_t  bytes,   ///< size of the serialized node buffers
   const          int     ncnt     ///< number of serialized nodes
)
{
   SSF_Nodes_Remap(cpu,gpu,ncnt);   // remap nodes/pointers for device (gpu) memory

   // copy the serialized/remapped nodes to the gpu...

   if (cpu && gpu && (bytes>0))
   { CUDART_CHECK( cudaMemcpy( (void *) gpu, (void *) cpu, bytes, cudaMemcpyHostToDevice )); }
}

// SSF_Recv_Nodes()
//
// Will pull the serialized nodes from the device (gpu) and remap them for host (cpu) memory.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Nodes_Recv
(
         unsigned char   *cpu  ,   ///< points to a serialized node buffer on host   (cpu)
         unsigned char   *gpu  ,   ///< points to a serialized node buffer on device (gpu)
   const          size_t  bytes,   ///< size of the serialized node buffers
   const          int     ncnt     ///< number of serialized nodes
)
{
   // retrieve the serialized nodes from the gpu...

   if (cpu && gpu && (bytes>0))
   { CUDART_CHECK( cudaMemcpy( (void *) cpu, (void *) gpu, bytes, cudaMemcpyDeviceToHost )); }

   SSF_Nodes_Remap(cpu,ncnt);   // remap nodes/pointers for host (cpu) memory
}

