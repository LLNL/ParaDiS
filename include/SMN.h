#pragma once

#ifndef _PDS_SMN_H
#define _PDS_SMN_H

/*************************************************************************
 *
 *      Module:      SMN.h  (Split Multi Nodes)
 *
 *      Description: Define various structures for data generated/needed
 *                   during the multi-node splitting process and prototype
 *                   the related functions.
 *
 *                   NOTE: Some functions are only for use with the
 *                         parallel multi-node splitting process, some
 *                         for the serial, and some for both.
 *
 *************************************************************************/

/*
 *      Define a structure in which to store pointers to both
 *      nodal endpoints of a segment.
 */
typedef struct {
        Node_t *node1;
        Node_t *node2;
}  SegNodes_t;

/*
 *      Structure to hold components of force at both endpoints
 *      of a dislocation segment.
 */
typedef struct {
        real8 f1[3];
        real8 f2[3];
}  SegForce_t;

/*
 *      Struct to hold the forces associated with the segments
 *      of a node.  Space provided for maximum number of segments
 *      although some may not be used.
 */
typedef struct {
        SegForce_t segForce[MAX_NBRS];
} SMN_Force_t;

/*
 *      Define a structure to hold forces calculated for the segments
 *      of two nodes resulting from the split of a multinode.
 *
 *      node1* values correspond to the multinode that was split
 *      node2* values correspond to the new node created during the split
 */
typedef struct {
        int         node1Segs;  /* # of segments associated with node1 */
        int         node2Segs;  /* # of segments associated with node2 */
        SMN_Force_t node1Force; /* segment forces for node1 */
        SMN_Force_t node2Force; /* segment forces for node2 */
} SMN_Split_t;

/*
 *      Define a structure to hold partial forces (and other info)
 *      calculated for a multinode for which we are evaluating
 *      splits in parallel on different domains.
 */
typedef struct {
        int status;  /* 0 == OK, 1 = failed to split, 2 == mobility error */
        int numSplits;
        Node_t *multiNode;
        SMN_Split_t *splitInfo; /* array: <numSplits> elements */
} SMN_Info_t;

/*
 *      Prototype the functions needed in support of the parallel
 *      multi-node splitting
 */
void BuildGhostMultiNodeList(Home_t *home, int *nodeCount,
        SMN_Info_t **nodeInfo);

void BuildLocalMultiNodeList(Home_t *home, int *nodeCount,
        SMN_Info_t **nodeInfo);

void BuildRecvDomList(Home_t *home, int lmnCount, SMN_Info_t *lmnInfo,
        int *numRecvDoms, int **recvDomList);

void BuildSMNLocalSegList(Home_t *home, SegNodes_t **segList, int *numSegs);

int BuildSplitList(int totalArms, int splitCnt, int level, int countOnly,
        int *currentSplit, int **splitList);

void CopyForcesToNode(Home_t *home, Node_t *node, SMN_Force_t *forces);

void GetSavedSegForces(Home_t *home, Node_t *node, SegForce_t *segForces,
        int *segIndex);

int  HasShortSeg(Home_t *home, Node_t *node, real8 shortSegLenSq);

void MergeGhostNode(Home_t *home, int opClass, Node_t *node1, Node_t *node2,
        real8 *position, Node_t **mergedNode, int *status);

void SaveSegForces(Home_t *home, Node_t *node, SegForce_t *segForces,
        int *segIndex);

void SMNEvalMultiNodeSplits(Home_t *home, int nodeCount, SMN_Info_t *nodeInfo);

void SMNEvalSplitForces(Home_t *home, int nodeCount, SMN_Info_t *nodeInfo);

void SMNForceComm(Home_t *home, int sendDomCount, int *sendDoms,
        int *sendBufLen, real8 **sendBuf, int recvDomCount, int *recvDoms,
        int lmnCount, SMN_Info_t *lmnInfo);

void SMNGetArmSets(int numNbrs, int *setCnt, int ***armSetList);

void SMNPackForceBuffers(Home_t *home, int gmnCount, SMN_Info_t *gmnInfo,
        int *numSendBufs, int **sendDoms, int **sendBufLens,
        real8 ***sendBufs);

void SMNPartialForce(Home_t *home, Node_t *node, int numSegs,
        SegNodes_t *segList, SMN_Force_t *nodeForces);

void SMNUnpackBuf(Home_t *home, int bufLen, real8 *buf, int lmnCount,
        SMN_Info_t *lmnInfo);

void SplitGhostMultiNodes(Home_t *home, int nodeCount, SMN_Info_t *nodeInfo);

int  SplitGhostNode(Home_t *home, int opClass, Node_t *node, real8 *pos1,
        real8 *pos2, real8 *vel1, real8 *vel2, int armCount, int *armList,
        Node_t **splitNode1, Node_t **splitNode2, int flags);

void SplitLocalMultiNodes(Home_t *home, int nodeCount, SMN_Info_t *nodeInfo);

void SplitMultiNodesParallel(Home_t *home);

#endif
