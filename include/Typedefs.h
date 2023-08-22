#pragma once

#ifndef _PDS_TYPEDEFS_H
#define _PDS_TYPEDEFS_H

/***************************************************************************
 *
 *  Typedefs.h - This include file defines all the typedefs used for
 *               data structures in the code. It allows the typedefs
 *               to be used in various header files before the header
 *               file where the typedef'd structure is defined. It should
 *               be included at the top of any header file that references
 *               structures defined in other header files
 *
 *************************************************************************/

#include "Constants.h"

class Complex_t;

typedef float     real4;
typedef double    real8;
typedef Complex_t complex8;

typedef struct _binnode BinNode_t;
typedef struct _binseg BinSeg_t;
typedef struct _dcompdomain DcompDomain_t;
typedef struct _home Home_t;
typedef struct _indata InData_t;
typedef struct _innode InNode_t;
typedef struct _param Param_t;
typedef struct _remotedomain RemoteDomain_t;
typedef struct _sortnode SortNode_t;
typedef struct _unmappedarm_t UnMappedArm_t;
typedef struct _SegPartIntersect_t SegPartIntersect_t;

class Cell_t;
class Node_t;
class NodeBlock_t;
class Segment_t;
class SegmentPair_t;

typedef enum {
	Periodic=0,
	Free=1,
	Reflecting=2
} BoundType_t;


/*
 *      Define a structure used when creating cell2 queues for
 *      collision handling.  An array of these structures will
 *      be allocated and initialized in SortNodesForCollision.c
 */
typedef struct {
	Node_t	*node;
	int	next;
} C2Qent_t;


/*
 *      Define a structure of data containing arrays and
 *      miscellaneous data needed for writing binary restart
 *      files.
 */
typedef struct {
        int   firstInGroup;
        int   lastInGroup;
        int   nodeCount;
        int   segCount;
        int   *nodeIndex;
        int   *nodeConstraint;
        int   *nodeNumSegs;
        int   *segTags;
        real8 *nodePos;
        real8 *burgersVec;
        real8 *glidePlane;
} BinFileData_t;

#ifdef ESHELBY
/*
 *      Define a structure of the information needed for handling
 *      eshelby inclusions
 */
typedef struct {
   int   id;                 /* particle id mainly for debugging */
   int   cellID;             /* ID of cell owning particle's center */
   int   nextInCell;         /* array index of next particle in cell */
   real8 position[3];        /* Coordinates of center of particle */
   real8 radius[3];          /* particle semi-principal axes in units of b */
   real8 rotation[2][3];     /* rotation matrix : two vectors only. The third one is the cross product. */
   real8 strainField[6];     /* strain field for the particle. Used when all particles are spherical for now*/
} EInclusion_t;
#endif  // ESHELBY

/*
 *      Define a couple structures used for managing the control
 *      and data file parameter lists
 */
typedef struct {
        char varName[MAX_STRING_LEN];
        int  valType;
        int  valCnt;
        int  flags;
        void *valList;
} VarData_t;


typedef struct {
        int       paramCnt;
        VarData_t *varList;
} ParamList_t;

/*
 *      Define a structure with some info to keep track of how much
 *      of which types of data a remote domain will be receiving
 *      during the force calcs.
 */
typedef struct {
        int domID;  /* ID of domain to which data will be sent */

        int segCnt; /* number segments for which force data is being sent */

        int intersectDataCnt; /* number of seg/particle intersection structs */
                              /* being sent                                  */
} SegCommInfo_t;

/*
 *      For some materials, segments of certain burgers vectors can be
 *      splintered into two segments of a different burgers vector to
 *      yield a lower energy configuration.  This structure maps one of
 *      the splinterable reference burgers vectors to the new burgers
 *      vector into which it may be splintered.
 *
 *      Note: These burgers vectors are not normalized.
 */
typedef struct {
        real8 refBurg[3];  /* reference burgers vector */
        real8 newBurg[3];  /* burgers vector into which it can be split */
} SplinterableBurg_t;


/*
 *      All the mobility functions must have the same calling sequence, but
 *      some of the mobility functions may return information not required
 *      by others (for instance, non-planar mobilities will not raeturn any
 *      glide constraints).  So, we define a structure containing various
 *      mobility function parameters that may not be required by all the
 *      mobility functions.
 */
typedef struct {
        real8 invDragMatrix[3][3];    /* Array in which to return the  */
                                      /* inverse of the drag matrix    */
                                      /* used by the mobility function */
                                      /*                               */
                                      /* NOTE: ALL mobility functions  */
                                      /* Need to return this data      */

        real8 glideConstraints[3][3]; /* Array in which to return the    */
                                      /* normals to the glide planes     */
                                      /* to which the node is restricted */
                                      /*                               */
                                      /* NOTE: Only mobility functions   */
                                      /* that enforce glide planes need  */
                                      /* return this data.               */

        int   numGlideConstraints;    /* Number of glide plane normals  */
                                      /* returned in <glideConstraints> */
                                      /*                                */
                                      /* NOTE: Mobility functions that  */
                                      /* do not enforce glide planes    */
                                      /* must set this value to zero,   */
                                      /* all other mobility functions   */
                                      /* must set this to the number of */
                                      /* constraints returned in the    */
                                      /* <glideConstraints> array.      */
} MobArgs_t;

#endif  // _PDS_TYPEDEFS_H
