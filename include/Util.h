#pragma once

#ifndef _PDS_UTIL_H
#define _PDS_UTIL_H

/****************************************************************************
 *
 *      Util.h  Contains miscellaneous definitions, and prototypes for
 *              functions defined in Util.c and used throughout the
 *              rest of the code.
 *
 ***************************************************************************/

#include "Home.h"


#define DotProduct(vec1,vec2)       \
            ((vec1[0])*(vec2[0]) +  \
             (vec1[1])*(vec2[1]) +  \
             (vec1[2])*(vec2[2]))


/***************************************************************************
 *
 *                              Prototypes
 *
 **************************************************************************/
/*
 *      Structure in which we preserve the state of the segments from
 *      the original node.
 */
typedef struct {
        Tag_t nbrTag;
        real8 burg[3];
        real8 planeNorm[3];
        real8 fNode[3]; /* force on the local end of the segment */
        real8 fNbr[3];  /* force on the end of segment at the neighbor */
} SegState_t;

/*
 *      Structure in which we preserve some of the nodal data items from
 *      the original node.
 */
typedef struct {
        int    numNbrs;
        Tag_t  origNodeTag;
        real8  x, y, z;
        real8  fx, fy, fz;
        real8  vx, vy, vz;
        SegState_t *segState;
} NodeState_t;

/*
 *      3D to linear index mapping functions
 */
extern void    DecodeCellIdx                (Home_t *home, int idx, int *iX, int *jY, int *kZ);
extern void    DecodeCell2Idx               (Home_t *home, int idx, int *iX, int *jY, int *kZ);
extern void    DecodeDomainIdx              (Home_t *home, int idx, int *iXdom, int *jYdom, int *kZdom);
extern int     EncodeCellIdx                (Home_t *home, int iX, int jY, int kZ);
extern int     EncodeCell2Idx               (Home_t *home, int iX, int jY, int kZ);
extern int     EncodeDomainIdx              (Home_t *home, int iXdom, int jYdom, int kZdom);

/*
 *      Various numerical operations
 */
extern void    cross                        (real8 a[3], real8 b[3], real8 c[3]);
extern void    CSpline                      (real8 *x, real8 *y, real8 *y2, int numPoints);
extern void    CSplint                      (real8 *xa, real8 *ya, real8 *y2, int numPoints, real8 x, real8 *y);
extern void    DecompVec                    (real8 inVec[3], real8 vec1[3], real8 vec2[3], real8 ansVec[2]);
extern void    FindAbsMax                   (real8 *array, int numElements, real8 *maxVal, int *indexOfMax);
extern void    FindAbsMin                   (real8 *array, int numElements, real8 *minVal, int *indexOfMin);
extern void    FindMax                      (real8 *array, int numElements, real8 *maxVal, int *indexOfMax);
extern void    FindMaxExcluding             (real8 *array, int numElements, real8 *maxVal, int *indexOfMax, int indexOfMax1=-1, int indexOfMax2=-1);
extern void    FindMin                      (real8 *array, int numElements, real8 *minVal, int *indexOfMin);
extern void    GetPlaneNormFromPoints       (real8 p1[3], real8 p2[3], real8 p3[3], real8 normalVec[3]);
extern void    GetUnitVector                (int unitFlag, real8 vx0, real8 vy0, real8 vz0, real8 vx1, real8 vy1, real8 vz1, real8 *ux, real8 *uy, real8 *uz, real8 *disMag);
extern void    InterpolateL                 (real8 *xa, real8 *ya, int numPoints, real8 x, real8 *y);
extern int     InterpolateBL                (real8 *x, real8 *y, real8 *val, int numElem, real8 targetX, real8 targetY, real8 *result);
extern real8   Normal                       (real8 a[3]);
extern void    Normalize                    (real8 *ax, real8 *ay, real8 *az);
extern void    NormalizeVec                 (real8 vec[3]);
extern void    CrossVector                  (real8 a[3], real8 b[3], real8 c[3]);
extern void    NormalizedCrossVector        (real8 a[3], real8 b[3], real8 c[3]);
extern void    Orthogonalize                (real8 *ax, real8 *ay, real8 *az, real8 bx, real8 by, real8 bz);
extern void    xvector                      (real8 ax, real8 ay, real8 az, real8 bx, real8 by, real8 bz, real8 *cx, real8 *cy, real8 *cz);


/*
 *      Nodal manipulations and operations
 */
extern void    AllocNodeArms                (Node_t *node, int n);
extern void    FreeNode                     (Home_t *home, int index);
extern void    FreeNodeArms                 (Node_t *node);
extern Node_t *GetNeighborNode              (Home_t *home, Node_t *node, int n);
extern Node_t *GetNodeFromIndex             (Home_t *home, int domID, int index);
extern Node_t *GetNodeFromTag               (Home_t *home, Tag_t tag);
extern void    InsertArm                    (Home_t *home, Node_t *nodeA, Tag_t *nodeBtag, real8 bx, real8 by, real8 bz, real8 nx, real8 ny, real8 nz, int Log);
extern void    MarkNodeForceObsolete        (Home_t *home, Node_t *node);
extern void    PrintNode                    (Node_t *node);
extern void    ReallocNodeArms              (Node_t *node, int n);
extern void    RemoveNode                   (Home_t *home, Node_t *node, int Log);
extern void    RepositionNode               (Home_t *home, real8 newPos[3], Tag_t *tag, int globalOp);
extern void    ResetNodeArmForce            (Home_t *home, Node_t *node);
extern void    SubtractSegForce             (              Node_t *node1, Node_t *node2);


/*
 *      Recycled node handling functions
 */
extern int     GetFreeNodeTag               (Home_t *home);
extern int     GetRecycledNodeTag           (Home_t *home);
extern void    RecycleNodeTag               (Home_t *home, int tagIdx);


/*
 *      PBC Image reconciliation functions
 */

extern void    FoldBox                      (const Param_t *param, real8 *x, real8 *y, real8 *z);
extern void    FoldBox                      (const Param_t *param, real8 &x, real8 &y, real8 &z);
extern void    FoldBox                      (const Param_t *param, real8 *p);

extern void    PBCPOSITION                  (const Param_t *param, const real8 x0, const real8 y0, const real8 z0, real8 *x, real8 *y, real8 *z);
extern void    PBCPOSITION                  (const Param_t *param, const real8 x0, const real8 y0, const real8 z0, real8 &x, real8 &y, real8 &z);
extern void    PBCPOSITION                  (const Param_t *param, const real8 *p0, real8 *p);

extern void    ZImage                       (const Param_t *param, real8 *x, real8 *y, real8 *z);
extern void    ZImage                       (const Param_t *param, real8 &x, real8 &y, real8 &z);
extern void    ZImage                       (const Param_t *param, real8 *p);


/*
 *      Topological operation and modification support functions
 */
extern void    ChangeArmBurg                (Home_t *home, Node_t *node1, Tag_t *tag2,
                                            real8 bx, real8 by, real8 bz, real8 nx,
                                            real8 ny, real8 nz, int Log, real8 del_seg_factor);
extern int     ChangeConnection             (Home_t *home, Node_t *node1, Tag_t *tag2, Tag_t *tag3, int Log);
extern void    CompressArmLists             (Node_t *node);
extern int     GetArmID                     (              Node_t *node1, Node_t *node2);
extern void    GetBurgersVectorNormal       (              Node_t *node1, Node_t *node2, real8 *bx, real8 *by, real8 *bz, real8 *nx, real8 *ny, real8 *nz);
extern void    RecalcSegGlidePlane          (Home_t *home, Node_t *node1, Node_t *node2, int ignoreIfScrew);
extern void    ResetGlidePlane              (Home_t *home, real8 newPlane[3], Tag_t *tag1, Tag_t *tag2, int globalOp);
extern void    ResetSegForces               (Home_t *home, Node_t *nodeA, Tag_t *nodeBtag, real8 fx, real8 fy, real8 fz, int globalOp);
extern void    ResetSegForces2              (Home_t *home, Node_t *nodeA, Tag_t *nodeBtag, real8 f1x, real8 f1y, real8 f1z, real8 f2x, real8 f2y, real8 f2z, int globalOp);
extern void    ResetNodalVelocity           (Home_t *home, Node_t *nodeA, real8 vx, real8 vy, real8 vz, int globalOp);

/*
 *      Node ordering and sorting functions
 */
extern int     CollisionNodeOrder           (Home_t *home, Tag_t *tagA, Tag_t *tagB);
extern int     IntCompare                   (const void *a, const void *b);
extern int     DomainOwnsSeg                (Home_t *home, int opClass, int thisDomain, Tag_t *endTag);
extern int     NodeCmpByTag                 (const void *a, const void *b);
extern int     OrderNodes                   (const void *a, const void *b);
extern int     OrderTags                    (const void *a, const void *b);
extern void    SortNativeNodes              (Home_t *home);

/*
 *      Prototypes for output related functions
 */
extern void    Fatal                        (const char *format, ...);
extern void    Warning                      (const char *format, ...);

extern void    Plot                         (Home_t *home, int domIndex, int blkFlag);
extern void    WriteArms                    (Home_t *home, char *baseFileName, int ioGroup, int firstInGroup, int writePrologue, int writeEpilogue);
extern void    WriteDensityField            (Home_t *home, char *fileName);
extern void    WritePoleFig                 (Home_t *home, char *baseFileName, int ioGroup, int firstInGroup, int writePrologue, int writeEpilogue);
extern void    WritePovray                  (Home_t *home, char *baseFileName, int ioGroup, int firstInGroup, int writePrologue, int writeEpilogue);

extern void    FindCellCenter               (Param_t *param, real8 x, real8 y, real8 z, int type, real8 *xCenter, real8 *yCenter, real8 *zCenter);
extern void    LocateCell                   (Home_t *home, int *cellID, real8 coord[3]);
extern void    Meminfo                      (int *wss);
extern void    _CheckMemUsage               (Home_t *home, char *msg);
extern int     NodeHasSessileBurg           (Home_t *home, Node_t *node);
extern int     NodePinned                   (Home_t *home, Node_t *node, int planeIndex, real8 (*glidedir)[3], int numGlideDirs);
extern real8   randm                        (int *seed);
extern void    ReadTabulatedData            (char *fileName, int numCols, real8 ***colData, int *numRows);
extern void    ReadRijm                     (Home_t *home);
extern void    ReadRijmPBC                  (Home_t *home);
extern real8   Sign                         (real8 x);
extern int     StrEquiv                     (const char *s1, const char *s2);
extern void    testdeWitStress2             ();
extern void    Uniq                         (int*, int*);

extern void    GetGlideConstraintList       (Node_t *node, int *numGlideConstraints, real8 (*glideConstraints)[3]);

extern void    ApplyConstraintsToVelocity   (int numGlideConstraints, real8(*glideConstraints)[3], real8 velocity[3]);
extern void    ApplyShearVelocityCap        (Home_t *home, real8 velocity[3]);

extern int     Get2ndMaxIndex               (real8 *A, int numElements, int indexOfMax);
extern int     Get2ndMaxIndex               (real8 *A, int numElements, int indexOfMax1, int indexOfMax2, int indexOfMax3);

extern void    PrintAllNodes                (Home_t *home);

extern void ThreeExtremaIndex(real8 x, real8 y, real8 z, int *imax,int *imin);
extern void Sort3(real8 values[3], int  ind[3]);
extern void Print3(const char *msg, real8 A[3]   );
extern void Print3x3(const char *msg, real8 A[3][3]);
extern void Print2x3(const char *msg, real8 A[3][3]);

//
// Useful topological routines
//

extern void InitState(NodeState_t *origNodeState);

extern void PreserveState(Home_t *home, Node_t *node,
                          NodeState_t *origNodeState);

extern void RestoreState(Home_t *home, Node_t *node,
                         NodeState_t *origNodeState);

extern void DiscardState(NodeState_t *origNodeState);


extern void PrintVec(const char *msg, real8 A0, real8 A1, real8 A2, int i);

#ifdef ESHELBY
void PrintInclusion(EInclusion_t       *inclusion);
#endif

#endif
