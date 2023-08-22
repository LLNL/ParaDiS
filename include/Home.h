#pragma once

#ifndef _PDS_HOME_H
#define _PDS_HOME_H

/****************************************************************************
 *
 *  Home.h  Define the Home struct which holds or points to all relevant
 *          data for this domain, plus define miscellaneous things
 *          that are used throughout the code.
 *
 *  Internal Units:
 *
 *    length in b (burgMag, read in meters, e.g. 2.725e-10)
 *    stress in Pa
 *    time   in second
 *
 *    force  in Pa*b^2
 *    force per unit length in Pa*b
 *    velocity in b/second
 *
 *    mobility in 1/Pa/second
 *
 ***************************************************************************/

#if defined(_ARLFEM) && defined(_ARLFEM_FULL_IMAGE_FORCES)
#include <H5FDdsm.h>
#include <H5FDdsmManager.h>
#include <hdf5.h>
#endif

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#ifdef USE_HDF
#include <hdf5.h>
#endif

#ifdef USE_SCR
#include <scr.h>
#endif

#include "mpi_portability.h"
#include "cuda_portability.h"

#if defined(USE_KINSOL) || defined(USE_ARKODE)
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_nvector.h>
#endif

#if defined(USE_ARKODE)
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_ls.h>
#endif

#include "Constants.h"
#include "ParadisThread.h"
#include "Typedefs.h"
#include "ParadisProto.h"
#include "FM.h"
#include "Node.h"
#include "NodeBlock.h"
#include "Param.h"
#include "Cell.h"
#include "Tag.h"
#include "TagMap.h"
#include "RemoteDomain.h"
#include "Tag.h"
#include "MirrorDomain.h"
#include "Topology.h"
#include "OpList.h"
#include "Timer.h"
#include "Util.h"
#include "Init.h"
#include "Parse.h"
#include "Matrix.h"
#include "DebugFunctions.h"
#include "Force.h"
#include "Segment.h"
#include "SegmentPair.h"

#ifdef ANISOTROPIC
#include "Anisotropic.h"
#include "AnisotropicVars_t.h"
#endif

#ifdef SPECTRAL
#include "Spectral.h"
#endif

/*
 *      Don't let the code compile with both PARALLEL and FULL_N2_FORCES
 *      defined...
 */
#if defined(PARALLEL) && defined(FULL_N2_FORCES)
#error "Cannot define FULL_N2_FORCES with PARALLEL (See makefile.setup)"
#endif

/*
 *      Don't let the user compile more than one variant of the FMM
 *      at a time...
 */
#if defined (TAYLORFMM) && defined (BBFMM)
#error "Cannot define TAYLORFMM and BBFMM together"
#endif

#if defined (TAYLORFMM) && defined (UNIFORMFMM)
#error "Cannot define TAYLORFMM and UNIFORMFMM together"
#endif

#if defined (BBFMM) && defined (UNIFORMFMM)
#error "Cannot define BBFMM and UNIFORMFMM together"
#endif


/*
 *      Anisotropy cannot be used with the 'rational' force, so don't
 *      let the compilation continue if both were selected.
 */
#if defined(ANISOTROPIC) && defined(USE_RATIONAL_SEG_FORCES)
#error "Cannot define ANISOTROPIC with USE_RATIONAL_SEG_FORCES. See makefile.setup file"
#endif

/*
 *      Add wrappers around the memory managment functions.  If
 *      memory debugging is enabled, these wrapper functions will
 *      allow us to do some tracking of allocated memory and some
 *      checks to help identify memory corruption and memory leaks.
 */
#define free(a)       ParadisFree   (__FILE__, __func__, __LINE__, (a))
#define malloc(a)     ParadisMalloc (__FILE__, __func__, __LINE__, (a))
#define calloc(a,b)   ParadisCalloc (__FILE__, __func__, __LINE__, (a), (b))
#define realloc(a,b)  ParadisRealloc(__FILE__, __func__, __LINE__, (a), (b))

void  ParadisFree   (const char *fname, const char *func, const int ln, void   *ptr );
void *ParadisMalloc (const char *fname, const char *func, const int ln, const size_t  size);
void *ParadisCalloc (const char *fname, const char *func, const int ln, const size_t  numElem, const size_t size);
void *ParadisRealloc(const char *fname, const char *func, const int ln, void   *ptr, const size_t size);

#ifdef DEBUG_MEM
void ParadisMemCheck         (void);
void ParadisMemCheckInit     (void);
void ParadisMemCheckTerminate(void);
#endif /* ifdef DEBUG_MEM */

/*
 *      On BGL, rint was *really* slow, so here's a kludge from one
 *      of the IBM guru's
 */
#ifdef _BGL
#define rint(x) (6.7553994410557440e+15 + (x) - 6.7553994410557440e+15)
#endif

#ifdef ESHELBY
/*
 *      Define a structure of information needed for tracking segments
 *      that intersect Eshelby inclusions.  This is only needed when using
 *      a mobility function that gives different mobilities to segments
 *      inside eshelby inclusions from those that are outside
 */
#define MAX_SEGPART_INTERSECTIONS 10

struct _SegPartIntersect_t {
        Tag_t tag1;  /* Tag of the segment's first endpoint. NOTE: must be */
                     /* "less than" tag2 when evaluated by OrderNodes()!   */

        Tag_t tag2;  /* Tag of the segment's secodn endpoint */

        int   numIntersections; /* Number of inclusions intersected by */
                                /* this segment.  Will never be greater*/
                                /* than MAX_SEGPART_INTERSECTIONS      */

        int   inclusionIndex[MAX_SEGPART_INTERSECTIONS]; /* Indices in the */
                                /* inclusion array of the inclusions that  */
                                /* intersect the segment                   */

        int   inclusionID[MAX_SEGPART_INTERSECTIONS]; /* ID of the inclusion */
                                /* intersecting the segment.  Curently only  */
                                /* needed when sending the intersect info    */
                                /* from one process to another.              */
};
#endif  // ESHELBY

/*
 *     Define a structure to hold the burgers vectors and associated
 *     glide planes for material for which we dynamically calculate
 *     those values.
 */
typedef struct {
        int   numBurgVectors;       /* Total number of burgers vectors */
                                    /* in <burgList> below             */

        int   numPlanes;            /* Total number of planes in */
                                    /* <planeList> below         */

        int   numSplinterableBurgs; /* Number of junction Burgers vectors  */
                                    /* that can be splintered              */

        int   *numPlanesPerBurg;    /* Number of planes associated with */
                                    /* each burgers vector in the list  */

        int   *burgFirstPlaneIndex; /* Index (0 offset) of the first plane */
                                    /* in <planeList> associated with each */
                                    /* burgers vector                      */

        int   *numGlissilePlanesPerBurg; /* Number of glissile planes    */
                                         /* associated with each burgers */
                                         /* vector in the <planeList>    */

        int   *planeType;          /* Array of <numPlanes> ints specifying   */
                                   /* the type for each plane in <planeList> */
                                   /* Types are:                             */
                                   /*     0 : basal                          */
                                   /*     1 : prismatic                      */
                                   /*     2 : pyramidal                      */
                                   /*     3 : second pyramidal               */
                                   /*     4 : sessile                        */

        real8 (*burgList)[3];       /* Array of burgers vectors */

        real8 (*planeList)[3];      /* Array of plane vectors associated with*/
                                    /* the burgers vectors in <burgList>     */

        real8 angleList[8];        /* Array of character angles*/
        int   lPlaneIndex[17];
        real8 (*lineList)[3];      /* Array of line vectors associated with*/
                                   /* the burgers vectors in <burgList>,   */
                                   /* planes vectors in <planeList> and    */
                                   /* angle in <angleList>                 */

        SplinterableBurg_t *splinterableBurgList;  /* Array of structs       */
                                    /* containing the junction Burgers       */
                                    /* vectors that can be splintered and    */
                                    /* to yield a lower energy configuration */
                                    /* and the burgers vectors into which    */
                                    /* they can splinter                     */
} BurgInfo_t;




struct _home {

        int         myDomain         ;  //< encoded MPI domain index for this domain
        int         numDomains       ;  //< number of MPI domains

        int         mpi_rank         ;  //< MPI rank index of this process
        int         mpi_size         ;  //< total number of MPI processes
        int         mpi_node_rank    ;  //< MPI rank index of this process on the compute node
        int         mpi_node_size    ;  //< number of MPI processes on the compute node
        char        mpi_hostname[256];  //< hostname of the node in which this process resides
                                        //   - note that if we are executing serially, the MPI variables
                                        //     will be initialized to rank=0, size=1, node_rank=1, etc.

        int         cycle            ;  //< current cycle
        int         lastCycle        ;  //< cycle to quit on

        long long   utc_cycle_start  ;  //< wall clock time at start of simulation (utc microsecs)
        long long   utc_cycle_prev   ;  //< wall clock time of previous cycle      (utc microsecs)
        long long   utc_cycle_curr   ;  //< wall clock time of current  cycle      (utc microsecs)

        Param_t   *param;

#ifdef ANISOTROPIC
        AnisotropicVars_t  anisoVars;
#ifdef TAYLORFMM
        FMAnisoTable_t    *fmAnisoTable;
#endif
#endif

#ifdef SPECTRAL
        Spectral_t *spectral;
#endif

/*
 *      The following two pointers are used for registering control
 *      and data file parameters and when reading/writing restart
 *      files.  They will only be set and used in task zero.
 */
        ParamList_t *ctrlParamList;
        ParamList_t *dataParamList;

/*
 *      The following three pointers are queue-heads of singly-linked lists of
 *      nodes. Nodes are pushed onto the head of a queue, but may be removed
 *      from anywhere in the queue. Blocks of nodes are allocated at one time
 *      for efficiency. These blocks are themselves kept in the nodeBlockQ,
 *      so they can be found if they ever need to be.
 */
        Node_t    *nativeNodeQ;
        Node_t    *ghostNodeQ;
        Node_t    *freeNodeQ;

        Node_t    *lastFreeNode;
        Node_t    *lastGhostNode;

        NodeBlock_t  *nodeBlockQ;

/*
 *      Using a queue for the ghost nodes makes for some fast insertions
 *      to the queue and moving the ghost nodes back onto the free node
 *      queue, but is not conducive to threading loops over the ghost nodes.
 *      Therefore, we also keep an array of ghost node pointers that can
 *      be used where the code must loop through all the ghost nodes.
 *      We'll also need to keep a count of the ghost nodes and the size
 *      of the allocated array (ie. max number of elements that can be placed
 *      on the current array).
 */
        Node_t **ghostNodeList;
        int    ghostNodeListSize;
        int    ghostNodeCount;

/*
 *      the nodeKeys array contains pointers to Node_t structs. For a node
 *      with a given tag, nodeKeys[tag.index] points to the node's struct.
 *      The recycle node heap contains the indices of tags that were in
 *      use, but have been freed. When a node is created or moves into the
 *      domain, its tag is assigned from the recycle node heap if possible.
 *      Note: node tags retrieved from the heap will always be such that
 *      the tag retrieved is the lowest available tag.
 */
        Node_t    **nodeKeys;
        int       newNodeKeyPtr;
        int       newNodeKeyMax;

        int       *recycledNodeHeap;
        int       recycledNodeHeapSize;
        int       recycledNodeHeapEnts;

#ifdef ESHELBY
/*
 *      Some variables for handling the list of Eshelby inclusions
 *      in the simulation
 */
        int          totInclusionCount;  /* Total number of inclusions   */
                                         /* on <eshelbyInclusions> array */
                                         /* includes ghost inclusions. */
                                         /* NOTE: all locally owned      */
                                         /* inclusions will preceed the  */
                                         /* ghost inclusions on the aaray*/

        int          locInclusionCount;  /* Number of locally owned in-     */
                                         /* clusions on <eshelbyInclusions> */

        EInclusion_t *eshelbyInclusions;
#endif  // ESHELBY

/*
 *      cellList keeps a list of all base cells of interest to this domain,
 *      whether native, ghost, or the base cell of a periodic ghost cell.
 *      The first nativeCellCount items in the list are native cells.
 *
 *      cellTable is a hash table for storing the standard cell structures
 *      required by the process.  Once allocated, the cell structure will
 *      remain on the hash table for the duration of the simulation.  This
 *      may result in some cells being unnecessarily saved, but eliminates
 *      the code to free and reallocate the needed cell structures every
 *      timestep.
 */
        int       *cellList;
        int       cellCount;
        int       nativeCellCount;

        Cell_t    *cellTable[CELL_HASH_TABLE_SIZE];

/*
 *      There are two classes of remote domains. The primary remote domains
 *      are those associated with any primary ghost nodes (i.e. any
 *      domain that intersects any cell that intersects (or is an immediate
 *      neighbor of a cell native to the current domain).  The second class
 *      of remote domains are those from which the current domain only
 *      requires secondary ghost nodes.
 */
        int       remoteDomainCount; /* Number of primary remote domains */

        int       secondaryRemoteDomainCount;  /* Number of secondary */
                                               /* remote domains only */

        int       *remoteDomains;    /* encoded indices of the neighbor */
                                     /* domains.  Includes primary and  */
                                     /* secondary remote domains        */

        RemoteDomain_t  **remoteDomainKeys; /* pointers to RemoteDomain_t    */
                                            /* structs of neighboring remote */
                                            /* domains                       */
/*
 *      To allow for multiple types of domain decomposition, we
 *      maintain a generic pointer to decomposition data here.
 *      The pointer will be cast to the proper type of pointer
 *      based on the type of decomposition specified in the control
 *      file.
 */
        void      *decomp;

/*
 *      Also keep the domain boundaries for this domain, in more
 *      accessible form... saves us from looking these values up
 *      in the domain decomposition multiple times during each step.
 */
        real8     domXmin;
        real8     domXmax;
        real8     domYmin;
        real8     domYmax;
        real8     domZmin;
        real8     domZmax;

/*
 *      For the RCB domain decomposition, the number of levels in the
 *      decomposition hierarchy is dependent on the number of domains
 *      in each dimension.  So, at startup determine the maximum levels
 *      and store the values for future reference.
 */
        int       xMaxLevel;
        int       yMaxLevel;
        int       zMaxLevel;

/*
 *      Some MPI stuff. Might be better if it was in a comm sub-structure
 */
#ifdef PARALLEL
        int       maxPackSiz;     /* byte length required to accomodate */
                                  /* the largest communication buffer  */
                                  /* being sent.  Nearly obsolete: is  */
                                  /* now only used in mirror node comm */

        MPI_Request  *inRequests; /* used for asynchronous I/O */
        MPI_Request  *outRequests;
        MPI_Status   *inStatus;
        MPI_Status   *outStatus;

#endif

/*
 *      array of Burgers vectors -- almost obsolete
 */
        real8     *burgX;
        real8     *burgY;
        real8     *burgZ;

        int       nburg;

/*
 *      Used for various types of output generation when passing nodal
 *      information from ALL domains to domain 0
 */
        MirrorDomain_t  **mirrorDomainKeys;     /* image of entire domains */
        int             currentMirrors;         /* 1 if mirrors already sent */
                                                /* this cycle; else 0        */
        char      *inBuf;
        char      *outBuf;

/*
 *      Operation list for topological changes across domain
 */
        char     *opListBuf;
        int       opListAllocLen;
        int       opListUsedLen;


        Timer_t   *timers;

/*
 *      Define values related to cell2 grid overlaid on standard
 *      cell structure during collision handling.
 */
        int       *cell2;
        C2Qent_t  *cell2QentArray;

        int       cell2nx;  /* number of 'cell2' elements in X dimension */
        int       cell2ny;  /* number of 'cell2' elements in Y dimension */
        int       cell2nz;  /* number of 'cell2' elements in Z dimension */
/*
 *      array to hold net charge tensor of each cell
 */
        real8     *cellCharge;

/*
 *      During initialization and node migration, nodes may move
 *      between domains.  Nodes moved between domains receive new
 *      ID tags.  The tagMap pointer will be used to hold an array
 *      of mappings between old and new tags for these moving
 *      nodes with tagMapSize and tagMapEnts indicating the number
 *      of elements currently allocated in the tagMap and the
 *      number of those elements that are actually in use.
 */
        TagMap_t  *tagMap;
        int       tagMapSize;
        int       tagMapEnts;

/*
 *      Remote force calculations use a couple arrays for doing
 *      Gauss-Legendre integration, so we'll keep them here where
 *      they'll be available whenever needed.
 */
        real8     *glPositions;
        real8     *glWeights;

/*
 *      Define things needed for the fast-multipole stuff
 */
        FMLayer_t *fmLayer;         /* Pointer to array of structures */
                                    /* (1 per layer of the FM hier-   */
                                    /* archy).                        */

        int       fmNumMPCoeff;     /* Number of coefficients required*/
                                    /* for multipole expansions of the*/
                                    /* order specified in the control */
                                    /* file.                          */

        int       fmNumExpansionCoeff; /* Number of coefficients required*/
                                       /* for FMM expansions of the      */
                                       /* order specified in the control */
                                       /* file.                          */

         int       eshelbyfmNumMPCoeff; /* Number of coefficients required*/
                                        /* for multipole expansions of the*/
                                        /* order specified in the control */
                                        /* file.                          */

        int       eshelbyfmNumTaylorCoeff;/* Number of coefficients required*/
                                    /* for talyor expansions of the   */
                                    /* order specified in the control */
                                    /* file.                          */

        int       eshelbyFMInitialized;  /* Flag indicating if the    */
                                    /* eshelby FM coefficients have   */
                                    /* been calculated -- since the   */
                                    /* inclusions do not move through */
                                    /* the simulation, we only need   */
                                    /* calculate the coefficients once*/

/*
 * Need some padded size for FFTs.
 */

  int fmNumMPCoeffsm;
  int fmNumExpansionCoefflg;


/*
 *      When using a mobility function that supports a different
 *      mobility inside eshelby inclusions than outside, we need
 *      a temporary list of inclusion-intersecting segments for
 *      the mobility function.
 */
        int segpartIntersectCnt;      /* number of segments in list */
        int segpartIntersectListSize; /* allocated size of the list */
        SegPartIntersect_t *segpartIntersectList;
#ifdef _OPENMP
        omp_lock_t segpartIntersectListLock;
#endif

/*
 *      Used only when load-balancing is done based on the number of
 *      segment to segment force calculations rather than the actual
 *      wallclock time spent.
 */
        int       cycleForceCalcCount;

/*
 *      For some of the material types, we dynamically calculate the list
 *      of possible burgers vectors (and associated glide planes) during
 *      initialization.  For those material types, we'll store all that
 *      information in the BurgListInfo_t structure for use throughout the
 *      code as needed.
 */
        BurgInfo_t burgData;

/*
 *      For parallel IO, each process needs to know things like
 *      which IO group it's in, the first and last processes in
 *      the group and so on.  Since those values are static, we
 *      can just calculate them once and store them in the Home_t
 *      structure for later use.
 */
        int   ioGroupNum;    /* IO group containing this process */

        int   firstInIOGroup, lastInIOGroup; /* first and last processes in */
                                             /* this IO group. */

        int   prevInIOGroup, nextInIOGroup;  /* Previous and next processes */
                                             /* in this IO group */

        int   isFirstInIOGroup;   /* Set to 1 if this process is the  */
                                  /* first in its IO group */
        int   isLastInIOGroup;    /* Set to 1 if this process is the  */
                                  /* last in its IO group */
#ifdef PARALLEL
        MPI_Comm commLastInIOGroup;  /* Communicator encompassing only   */
                                     /* the last process in each IO group*/
#endif


#ifdef _ARLFEM
        H5FDdsmManager * dsmManager; /* Distributed shared memory manager */

        hid_t dsmFAPL;               /* File access property list for */
                                     /* distributed shared memory     */
/*
 *      When calculating Yoffe stress, we need to build a global list
 *      of segments that intersect the free surfaces.  These two
 *      variables are used for dealing with that list.  See code
 *      for description of contents of the segment list array.
 */
        int   surfaceSegCount;
        real8 *surfaceSegList;
#endif

/*
 *      Numerous portions of ParaDiS perform MPI reduction operations
 *      that act on buffers with 1 element per process.  Rather than
 *      allocate and deallocate such buffers every time they are needed
 *      we'll allocate them once during initialization and save pointers
 *      to them here.
 */
        int *localReduceBuf;
        int *globalReduceBuf;

#ifdef USE_ARKODE
        void *arkodeLinearSolver;
        void *arkodeNonlinearSolver;
        void *arkodeMemBlock;
        int   arkodeInitFlag;
#endif

        int isFirstTaskOnNode;  /* Flag indicating of the current MPI task */
                                /* is the first task on a compute-node.    */
                                /* If the code is compiled without parallel*/
                                /* or shard memory support, this flag will */
                                /* always be set to 1.                     */

#if defined(PARALLEL) && defined(SHARED_MEM)
/*
 *      When support for both shared memory buffers and parallel execution
 *      is enabled, the code needs to identify MPI tasks co-located on
 *      compute nodes and set up MPI communicators specific to those
 *      tasks.  Here we declare some stuff needed for that.
 */
        MPI_Group  groupOnNode;
        MPI_Comm   commNode;

        MPI_Group  groupFirstOnNode;
        MPI_Comm   commFirstOnNode;
#endif
};

/* prototypes */
/* FIX ME!  Most of these prototypes should be moved elsewhere */

extern void     AddNode              (Home_t *home, Node_t *nodeA, Node_t *newNode, Node_t *nodeB);
extern void     CommSendMirrorNodes  (Home_t *home, int stage);
extern int      Connected            (Node_t *node1, Node_t *node2, int *armidx);
extern Node_t  *GetNewNativeNode     (Home_t *home);
extern Node_t  *GetNewGhostNode      (Home_t *home, int domain, int index);
extern void     Gnuplot              (Home_t *home, char *baseFileName, int ioGroup, int firstInGroup, int writePrologue, int writeEpilogue);

#endif
