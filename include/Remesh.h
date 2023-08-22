#pragma once

#ifndef _PDS_REMESH_H
#define _PDS_REMESH_H

/***************************************************************************
 *
 *      Remesh.h   Declare function prototypes and other relevant
 *                 related data needed by modules handling the
 *                 mesh rediscretization
 *
 **************************************************************************/

#include "Home.h"

/*
 *      Define the set of statuses an entry in the coarsen or refine
 *      lists may be in.
 */
enum {
        REMESH_PENDING = 0,
        REMESH_COMMITTED,
        REMESH_REJECTED
};


/*
 *      Structure to hold information regarding nodes selected to be
 *      coarsened out during the first phase of MeshCoarsen().
 */
typedef struct {
        int   status;
        Tag_t nodeTag;
        Tag_t nbr1Tag;
        Tag_t nbr2Tag;
} CoarsenList_t;


/*
 *      Structure to hold information regarding segments selected to be
 *      refined out during the first phase of MeshRefine().
 */
typedef struct {
        int   status;
        int   haveNewPos; /* only for remesh version 3 */
        Tag_t node1Tag;
        Tag_t node2Tag;
        Tag_t newNodeTag;
        real8 vector[3];
        real8 f0Seg1[3];
        real8 f1Seg1[3];
        real8 f0Seg2[3];
        real8 f1Seg2[3];
        real8 newPos[3]; /* only for remesh version 3 */
} RefineList_t;


/*
 *      Prototype the global functions related to rediscretization
 */
void EstCoarsenForces(Home_t *home, Node_t *node1, Node_t *node2,
        Node_t *node3, real8 f0Seg[3], real8 f1Seg[3]);
void EstRefinementForces(Home_t *home, Node_t *node1, Node_t *node2,
        real8 newPos[3], real8 vec[3],
        real8 f0Seg1[3], real8 f1Seg1[3],
        real8 f0Seg2[3], real8 f1Seg2[3]);
void MeshCoarsen(Home_t *home);
int  ProceedWithCoarsen(CoarsenList_t *list, int listEnts, Tag_t *nodeTag,
        Tag_t *nbr1Tag, Tag_t *nbr2Tag);
void Remesh(Home_t *home);
void RemeshRule_2(Home_t *home);
void RemeshRule_3(Home_t *home);

#endif /* _REMESH_H */
