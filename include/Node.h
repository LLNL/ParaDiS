#pragma once

#ifndef _PDS_NODE_H
#define _PDS_NODE_H

/*--------------------------------------------------------------------------
 *
 *	Node.h	Define the struct that holds all relevant data for a single
 *		node, either native or ghost
 *
 *		Notice: make depend if new files are added, otherwise always
 *			make clean after changing this file  Wei Cai 04/09/2002
 *
 *------------------------------------------------------------------------*/

#include "mpi_portability.h"

#include "Typedefs.h"
#include "ParadisThread.h"
#include "Tag.h"
#include "Thread.h"
#include "Home.h"

/*
 *      Define the various permitted nodal constraints.  The nodal constraint
 *      is a single integer value in which each possible constraint is
 *      represented by a single bit in the integer.   This allows a node
 *      to be constrained in multiple ways.
 *
 *      WARNING: Currently, certain portions of ParaDiS (primarily the
 *               code dealing with topological changes) do not support
 *               pinning a node in only specific dimensions.  As a result,
 *               those portions of the code will treat a node that is
 *               pinned in ANY dimension as if it is pinned in ALL dimensions.
 */
#define CONSTRAINT_ANY      0xffffffff   /* -1 */
#define UNCONSTRAINED       0x00000000   /*  0 */
#define PINNED_X            0x00000001   /*  1 ONLY PARTIALLY SUPPORTED */
#define PINNED_Y            0x00000002   /*  2 ONLY PARTIALLY SUPPORTED */
#define PINNED_Z            0x00000004   /*  4 ONLY PARTIALLY SUPPORTED */
#define SURFACE_NODE        0x00000008   /*  8 */

#define PINNED_NODE         7            /* Special constraint which */
                                         /* represents all of:       */
                                         /*   PINNED_X               */
                                         /*   PINNED_Y               */
                                         /*   PINNED_Z               */

/*
 *      The following constraints are used in post-processing tools only
 */
#define MULTI_JUNCTION_NODE 0x00001000  /*  4096: Node has all 4 111 burgers */
                                        /*        vectors associated with it */

#define LOOP_NODE           0x00002000  /*  8192: Node is in a loop with a */
                                        /*        unique burgers vector    */

#define CROSS_NODE          0x00004000  /* 16384: Node is a non multi-    */
                                        /*        junction node with more */
                                        /*        than 2 neighbors        */


/*************************************************************************
 *
 *      Now define some macros for examining/manipulating node
 *      constraint values...
 *
 ************************************************************************/

/*
 *      If the constraint value <cVal> has any constraint corresponding
 *      to the bits set in the constraint list <cList>, return a value
 *      of 1, otherwise return zero.
 */
#define HAS_ANY_OF_CONSTRAINTS(cVal, cList) (((cVal) & (cList)) != 0)

/*
 *      If the constraint value <cVal> has all the constraints corresponding
 *      to the bits set in the constraint list <cList>, return a value
 *      of 1, otherwise return zero.
 */
#define HAS_ALL_OF_CONSTRAINTS(cVal, cList) (((cVal) & (cList)) == (cList))

/*
 *      If the constraint value <cVal> has all the constraints corresponding
 *      to the bits set in the constraint list <cList> and no other
 *      constraints, return a value of 1, otherwise return zero.
 */
#define HAS_ONLY_CONSTRAINTS(cVal, cList) ((cVal) == (cList))

/*
 *      If the constraint value <cVal> has none of the constraints
 *      corresponding to the bits set in the constraint list <cList>
 *      return a value of 1, otherwise return zero.
 */
#define HAS_NONE_OF_CONSTRAINTS(cVal, cList) (((cVal) & (cList)) == 0)

/*
 *      If the constraint value <cVal> represents an unconstrained node
 *      return a value of 1, otherwise return zero.
 */
#define HAS_NO_CONSTRAINTS(cVal) (cVal == UNCONSTRAINED)

/*
 *      Add all constraints specified in <cList> to <cVal>.
 */
#define ADD_CONSTRAINTS(cVal, cList) (cVal |= (cList))

/*
 *      Remove any constraints specified in <cList> from <cVal>.
 */
#define REMOVE_CONSTRAINTS(cVal, cList) (cVal &= (~(cList)))

/*
 *      Set the constraint value <cVal> to the constraint list
 *      specified in <cList>. 
 */
#define SET_CONSTRAINTS(cVal, cList) (cVal = (cList))


/*
 *      Define the bit flags that can be set for the node.  Note: these
 *      flags may be OR'ed together.
 */
#define NODE_RESET_FORCES    0x01
#define NODE_OUTSIDE_SURF    0x02
#define NO_COLLISIONS        0x04
#define NO_MESH_COARSEN      0x08
#define NODE_CHK_DBL_LINK    0x10

/*
 *      Used as bit flags to indicate the type of nodal data items
 *      preserved by calls to PreserveNodalData().
 */
#define NODE_POSITION 0x01
#define NODE_CURR_VEL 0x02
#define NODE_OLD_VEL  0x04

class Node_t 
{
   public:
      int      flags;
           
      real8    x , y , z ;                  //  nodal position 
      real8    fX, fY, fZ;                  //  nodal force: units=Pa*b^2) 
      real8    vX, vY, vZ;                  //  nodal velocity: units=burgMag/sec 
      
           
      real8    oldx  , oldy  , oldz  ;      //  for strain increment, Moono.Rhee 
      real8    oldvX , oldvY , oldvZ ;      //  nodal velocity at previous step 
      real8    currvX, currvY, currvZ;      //  nodal velocity at beginning of the 
                                            //  current step.  Only used during    
                                            //  timestep integration while vX, vY  
                                            //  and vZ are being determined.       
    
      Tag_t    myTag;                       // tag associated with this node

      int      myIndex;                     // index associated with this node (for gpu)
   
      int      numNbrs;                     // number of neighbors
      Tag_t   *nbrTag;                      // array of neighbor tags
      int     *nbrIndex;                    // array of node indices if neighbors (for gpu)
   
      real8   *armfx, *armfy, *armfz;       //  arm specific force contribution 
      real8   *burgX, *burgY, *burgZ;       //  burgers vector 
      real8   *nx   , *ny   , *nz   ;       //  glide plane 
   
      real8   *sigbLoc;                     //  sig.b on arms (numNbr*3) 
      real8   *sigbRem;                     //  sig.b on arms (numNbr*3) 
   
      int     *armCoordIndex;               //  Array of indices (1 per arm) into 
                                            //  the mirror domain's arrays of     
                                            //  coordinates (armX, armY, armZ)    
                                            //  of nodes' neighbors.  This array  
                                            //  is only present and useful while  
                                            //  task zero is downloading the data 
                                            //  from the remote domain for        
                                            //  generating output                 

      int     *armNodeIndex;                //  Array of node/arm indices into the 
                                            //  current node list.
   
      int      constraint;                  //  constraint =  1 : any surface node   
                                            //  constraint =  7 : Frank-Read end   
                                            //                    points, fixed   
                                            //  constraint =  9 : Jog      
                                            //  constraint = 10 : Junction      
                                            //  constraint = 20 : cross slip(default)
                                            //                    2-nodes non-deletable

   
      int   cellIdx;                        //  cell node is currently sorted into 
      int   cell2Idx;                       //  cell2 node is currently sorted into 
      int   cell2QentIdx;                   //  Index of this node's entry in the 
                                            //  home->cell2QentArray.             
   
      int   native;                         //  1 = native node, 0 = ghost node 
   
      Node_t   *next;                       //  pointer to the next node in the queue 
                                            //  (ghost or free)          
   
      Node_t   *nextInCell;                 //  used to queue node onto the current 
                                            //  containing cell             
   
      int   tmark;                          //   1: if marked as visited during topology traversal
                                            //   0: otherwise

#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
           int     multiNodeLife;
#endif
   
#ifdef _ARLFEM
           real8   surfaceNorm[3];
#endif
   
#ifdef _OPENMP
           omp_lock_t nodeLock;             // mutex lock for OpenMP threads
#endif
#ifdef PTHREADS
      pthread_mutex_t   node_plock;         // mutex lock for pthread threads
#endif
   
   public :
      Node_t(void); 
      Node_t(const Tag_t & tag); 
      Node_t(const Tag_t & tag, const int narms); 
      Node_t(const Node_t & n); 
     ~Node_t();

      const Node_t & operator =  (const Node_t & n);

            int operator <  (const Node_t & node) const;
            int operator >  (const Node_t & node) const;
            int operator == (const Node_t & node) const;
            int operator != (const Node_t & node) const;

               int      IsA_Discrete     (void) const { return ( (numNbrs==2) ? 1 : 0 ); }
               int      IsA_Junction     (void) const { return ( (numNbrs> 2) ? 1 : 0 ); }
               int      IsA_End          (void) const { return ( (numNbrs==1) ? 1 : 0 ); }
               int      Visited          (void) const { return ( (tmark  ==1) ? 1 : 0 ); }

               int      IsA_M_Node       (void) const;
               int      IsA_N_Node       (void) const;

               void     Position         (real8 & px, real8 & py, real8 & pz) const { px=x; py=y; pz=z; }
               void     Position         (real8  *p )                         const { if (p) { p[0]=x;  p[1]=y;  p[2]=z;  } }

               void     Force            (real8 & fx, real8 & fy, real8 & fz) const { fx=fX; fy=fY; fz=fZ; }
               void     Force            (real8  *f )                         const { if (f) { f[0]=fX; f[1]=fY; f[2]=fZ; } }

               void     Velocity         (real8 & vx, real8 & vy, real8 & vz) const { vx=vX; vy=vY; vz=vZ; }
               void     Velocity         (real8  *v )                         const { if (v) { v[0]=vX; v[1]=vY; v[2]=vZ; } }

               void     Allocate_Arms    (const int narms);
               void     Free_Arms        (void);
               void     Recycle          (void);

               size_t   Bytes            (void) const;

               int      Native           (void) const { return( native ? 1 : 0 ); }
               int      Ghost            (void) const { return( native ? 0 : 1 ); }
               int      Validate         (void) const;

               int      Inside           (const real8 x0, const real8 x1, 
                                          const real8 y0, const real8 y1, 
                                          const real8 z0, const real8 z1 ) const;

               int      Arm_Index        (const Tag_t & tag) const;
               int      Arm_Index        (const real8 bx, const real8 by, const real8 bz) const;

               void     Delete_Arm       (const Tag_t & tag);

               void     Forces_Reset     (void);
               void     Add_Arm_Forces   (void);
               void     Add_Arm_Force    (const Tag_t & tag, const real8 *fv); 

               void     Arm_Sum          (real8 *f) const;

               void     Advance          (const real8  dt,
                                          const real8  cx, const real8 cy, const real8 cz,
                                          const real8  lx, const real8 ly, const real8 lz,
                                          const real8  sx, const real8 sy, const real8 sz );

               void     Preserve         (const int items);

               size_t   GPU_Bytes        (void);
      unsigned char    *GPU_Pack         (unsigned char *q);
      unsigned char    *GPU_Unpack       (unsigned char *p);

               int      MPI_Packed_Count (void) const;
               real8   *MPI_Pack         (real8 *pbuf) const;
               real8   *MPI_Unpack       (real8 *pbuf);

               int      MPI_Pack_Forces_Count  (void) const;
               real8   *MPI_Pack_Forces        (real8 *pbuf) const;
               real8   *MPI_Unpack_Forces      (real8 *pbuf);

               char    *Parse            (char *p , const char *pe);
               void     Load             (FILE *fd, const int fmt=0);

               void     Print            (FILE *fd) const;
               void     Print            (FILE *fd, const int info, const int pos, const int vel, const int arm, const int bvec, const int glide) const;
};

// Forward declarations... 

extern  int             Node_Count         (Node_t  *node );
extern  int             Node_Count         (Node_t **narr , int   cnt );
extern  Node_t        **Node_Array         (Node_t  *node , int & cnt );
extern  Node_t        **Nodes_Vector       (Node_t  *nodes, int   cnt );

extern  Node_t         *Nodes_Sort         (Node_t  *narr , int   cnt, int (*cmp)(const void *, const void *)=0 );
extern  Node_t        **Nodes_Sort         (Node_t **narr , int   cnt, int (*cmp)(const void *, const void *)=0 );

extern  Node_t         *Nodes_Find         (Node_t  *narr , int   cnt, const Tag_t & tag);
extern  Node_t         *Nodes_Find         (Node_t **narr , int   cnt, const Tag_t & tag);

extern  Node_t         *Node_Append        (Node_t  *nvec , int & cnt, int & max, Node_t  &node );
extern  Node_t        **Node_Append        (Node_t **nvec , int & cnt, int & max, Node_t  *node );
extern  Node_t        **Nodes_Append       (Node_t **nvec , int & cnt, int & max, Node_t **nodes, const int ncnt );

extern  void            Nodes_Set_Native   (Node_t  *nodes, int   cnt, 
                                            const real8 x0, const real8 x1, 
                                            const real8 y0, const real8 y1, 
                                            const real8 z0, const real8 z1 );

extern  void            Nodes_Set_Native   (Node_t **narr , int   cnt, 
                                            const real8 x0, const real8 x1, 
                                            const real8 y0, const real8 y1, 
                                            const real8 z0, const real8 z1 );

extern  void            Node_Forces_Reset  (Node_t  *narr , int ncnt);
extern  void            Node_Forces_Reset  (Node_t **narr , int ncnt);

extern  void            Node_Forces_Sum    (Node_t  *narr , int ncnt);
extern  void            Node_Forces_Sum    (Node_t **narr , int ncnt);

extern  void            Nodes_Preserve     (Node_t  *narr , const int ncnt, const int items);
extern  void            Nodes_Preserve     (Node_t **narr , const int ncnt, const int items);

extern  void            Nodes_Advance      (Node_t  *narr, const int ncnt, const real8 dt,
                                            const real8    cx, const real8 cy, const real8 cz,
                                            const real8    lx, const real8 ly, const real8 lz,
                                            const real8    sx, const real8 sy, const real8 sz );

extern  void            Nodes_Advance      (Node_t **narr, const int ncnt, const real8 dt,
                                            const real8    cx, const real8 cy, const real8 cz,
                                            const real8    lx, const real8 ly, const real8 lz,
                                            const real8    sx, const real8 sy, const real8 sz );

extern  int             Nodes_Validate     (Node_t  *narr , int ncnt);
extern  int             Nodes_Validate     (Node_t **narr , int ncnt);

extern  Node_t         *Nodes_Load         (const char *path, int & cnt);
extern  void            Nodes_Save         (const char *path, Node_t *narr, int cnt);
extern  void            Nodes_Save_GNU     (const char *path, Node_t *narr, int cnt, const real8 lx, const real8 ly, const real8 lz );

extern  int             MPI_Packed_Count   (      Node_t **narr , const int ncnt);
extern  int             MPI_Packed_Count   (      Node_t **narr , const int ncnt, const int rank);
extern  int             MPI_Packed_Count   (const Node_t  *nodes, const int ncnt);
extern  int             MPI_Packed_Count   (const Node_t  *nodes, const int ncnt, const int rank);

extern  real8          *MPI_Pack_Nodes     (      Node_t **narr , const int ncnt);
extern  real8          *MPI_Pack_Nodes     (const Node_t  *nodes, const int ncnt);
extern  void            MPI_Unpack_Nodes   (      Node_t  *nodes, const int ncnt, const real8 *mbuf);

extern  unsigned char  *Nodes_CPU_Mbuf     (void);
extern  unsigned char  *Nodes_GPU_Mbuf     (void);
extern  size_t          Nodes_GPU_Bytes    (void);

extern  size_t          Nodes_GPU_Bytes    (Node_t *nodes, int cnt);
extern  void            Nodes_GPU_Allocate (Node_t *nodes, int cnt);
extern  void            Nodes_GPU_Allocate (size_t  bytes);
extern  void            Nodes_GPU_Free     (void);
extern  void            Nodes_GPU_Pack     (Node_t *nodes , int cnt, float *msec=0);
extern  void            Nodes_GPU_Unpack   (Node_t *nodes , int cnt, float *msec=0);
extern  void            Nodes_GPU_Send     (float  *msec=0);
extern  void            Nodes_GPU_Recv     (float  *msec=0);

#endif
