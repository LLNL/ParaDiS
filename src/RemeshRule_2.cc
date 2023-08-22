/*****************************************************************************
 *
 *      Module:         RemeshRule_2.c
 *      Description:    This module contains functions specific to the
 *                      version 3 remesh for rediscretizing dislocation
 *                      segments.
 *                      
 *      Included functions:
 *              MeshRefine()
 *              RemeshRule_2()
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "QueueOps.h"
#include "Mobility.h"
#include "Remesh.h"

static int dbgDom;


/*-------------------------------------------------------------------------
 *
 *      Function:    MeshRefine
 *      Description: 
 *
 *------------------------------------------------------------------------*/
static void MeshRefine(Home_t *home)
{

/*
 *      Initialize some of the constants we'll need during this call.
 */
        int      thisDomain = home->myDomain;
        Param_t *param      = home->param;

        real8    areaMax  = param->remeshAreaMax;
        real8    areaMax2 = areaMax * areaMax;
        real8    maxSeg2  = param->maxSeg * param->maxSeg;
        real8    delta    = 1.0e-16;

        int localRefineCnt = 0;
        
/*
 *      Allocate a block of structures to hold info about rediscretization
 *      events.  We allocate enough structures to refine twice as many segments
 *      as there are nodes which is probably overkill, but it means we
 *      shouldn't run out of entries and have to reallocate the block,
 *      plus minimizes the number of operations we'll have to do within
 *      the thread-critical sections of code.
 */
        int           refineListSize = 2 * home->newNodeKeyPtr;
        int           refineListEnts = 0;
        RefineList_t *refineList     = (RefineList_t *)NULL;

        if (refineListSize > 0) {
            refineList = (RefineList_t *)malloc(sizeof(RefineList_t) * refineListSize);
        }

/*
 *      First phase is to loop through all the nodes, making a list
 *      of segments to be refined.
 */
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            int numNodes, nIndex;
            int threadID, threadIterStart, threadIterEnd;

            numNodes = home->newNodeKeyPtr;

            GetThreadIterationIndices(numNodes, &threadID, &threadIterStart,
                                      &threadIterEnd);

            for (nIndex = threadIterStart; nIndex < threadIterEnd; nIndex++) {
                int    segID;
                int    splitOK[2], splitSegList[2];
                real8  s, area2;
                real8  dvec1xdt, dvec1ydt, dvec1zdt;
                real8  dvec2xdt, dvec2ydt, dvec2zdt;
                real8  dvec3xdt, dvec3ydt, dvec3zdt;
                real8  dr1dt, dr2dt, dr3dt, dsdt, darea2dt;
                real8  r1, r2, r3;
                real8  vec1[3], vec2[3], vec3[3];
                Node_t *node, *nbr1, *nbr2, *nbr;
        
                if ((node = home->nodeKeys[nIndex]) == (Node_t *)NULL) {
                    continue;
                }
        
/*
 *              It can happen that a node is pinned and has two neighbors.
 *              Prevent remesh refinement of such a node.
 */
                if (HAS_ANY_OF_CONSTRAINTS(node->constraint, PINNED_NODE)) 
                    continue;                
/*
 *              If the node has only two arms, use both arm lengths and
 *              the area in the triangle defined by the node and its 
 *              neighbors to decide whether one of its segments should
 *              be refined.
 */
                if (node->numNbrs == 2) {
        
/*
 *              Calculate the lengths of the node's 2 arms plus
 *              the distance between the two neighbor nodes.
 */
                    nbr1 = GetNeighborNode(home, node, 0);
                    nbr2 = GetNeighborNode(home, node, 1);
        
                    if ((nbr1 == (Node_t *)NULL) || (nbr2 == (Node_t *)NULL)) {
                        printf("WARNING: Neighbor not found at %s line %d\n",
                               __FILE__, __LINE__);
                        continue;
                    }

                    vec1[X] = nbr1->x - node->x;
                    vec1[Y] = nbr1->y - node->y;
                    vec1[Z] = nbr1->z - node->z;

                    vec2[X] = nbr2->x - node->x;
                    vec2[Y] = nbr2->y - node->y;
                    vec2[Z] = nbr2->z - node->z;

                    vec3[X] = vec2[X] - vec1[X];
                    vec3[Y] = vec2[Y] - vec1[Y];
                    vec3[Z] = vec2[Z] - vec1[Z];

                    ZImage(param, &vec1[X], &vec1[Y], &vec1[Z]);
                    ZImage(param, &vec2[X], &vec2[Y], &vec2[Z]);
                    ZImage(param, &vec3[X], &vec3[Y], &vec3[Z]);

                    r1 = sqrt(DotProduct(vec1, vec1));
                    r2 = sqrt(DotProduct(vec2, vec2));
                    r3 = sqrt(DotProduct(vec3, vec3));
        
                    s = 0.5 * (r1 + r2 + r3);
                    area2 = (s * (s - r1) * (s - r2) * (s - r3));
        
/*
 *                  Determine if the area of the triangle defined by the node
 *                  and its two neighbors is increasing or decreasing.
 */
                    dvec1xdt = nbr1->vX - node->vX;
                    dvec1ydt = nbr1->vY - node->vY;
                    dvec1zdt = nbr1->vZ - node->vZ;

                    dvec2xdt = nbr2->vX - node->vX;
                    dvec2ydt = nbr2->vY - node->vY;
                    dvec2zdt = nbr2->vZ - node->vZ;

                    dvec3xdt = dvec2xdt - dvec1xdt;
                    dvec3ydt = dvec2ydt - dvec1ydt;
                    dvec3zdt = dvec2zdt - dvec1zdt;

                    dr1dt = ((vec1[X] * dvec1xdt) + (vec1[Y] * dvec1ydt) +
                             (vec1[Z] * dvec1zdt)) / (r1 + delta);

                    dr2dt = ((vec2[X] * dvec2xdt) + (vec2[Y] * dvec2ydt) +
                             (vec2[Z] * dvec2zdt)) / (r2 + delta);

                    dr3dt = ((vec3[X] * dvec3xdt) + (vec3[Y] * dvec3ydt) +
                             (vec3[Z] * dvec3zdt)) / (r3 + delta);

                    dsdt = 0.5 * (dr1dt + dr2dt + dr3dt);

                    darea2dt = (dsdt * (s-r1) * (s-r2) * (s-r3));
                    darea2dt += s * (dsdt-dr1dt) * (s-r2) * (s-r3);
                    darea2dt += s * (s-r1) * (dsdt-dr2dt) * (s-r3);
                    darea2dt += s * (s-r1) * (s-r2) * (dsdt-dr3dt);

/*
 *                  Start with the assumption that both segments are
 *                  permitted to be refined, then do some checks to see
 *                  if indded we are permitted to refine the segments.
 */
                    splitOK[0] = 1;
                    splitOK[1] = 1;

/*
 *                  If both nodes of the segment are flagged to have forces
 *                  updated, forces and velocities on the node may not be
 *                  good enough to accurately position the new node, so
 *                  don't refine this segment this cycle.
 */
                    if (((node->flags & NODE_RESET_FORCES) != 0) &&
                        ((nbr1->flags & NODE_RESET_FORCES) != 0)) {
                        splitOK[0] = 0;
                    }

                    if (((node->flags & NODE_RESET_FORCES) != 0) &&
                        ((nbr2->flags & NODE_RESET_FORCES) != 0)) {
                        splitOK[1] = 0;
                    }

/*
 *                  Check if the current domain owns the segments.
 *                  It may only refine a segment it owns...
 */
                    if (splitOK[0]) {
                        splitOK[0] = DomainOwnsSeg(home, OPCLASS_REMESH,
                                                   thisDomain, &nbr1->myTag);
                    }

                    if (splitOK[1]) {
                        splitOK[1] = DomainOwnsSeg(home, OPCLASS_REMESH,
                                                   thisDomain, &nbr2->myTag);
                    }

                    splitSegList[0] = 1;
                    splitSegList[1] = 2;

/*
 *                  Check both segments (in the preferred order) and 
 *                  determine if we should bisect one.
 */
                    for (segID = 0; segID < 2; segID++) {
                        int   segLen;
                        real8 *vec;

                        if (splitSegList[segID] == 1) {
                            if (splitOK[0] == 0) {
                                continue;
                            }
                            nbr = nbr1;
                            vec = vec1;
                            segLen = r1;
                        } else {
                            if (splitOK[1] == 0) {
                                continue;
                            }
                            nbr = nbr2;
                            vec = vec2;
                            segLen = r2;
                        }


/*
 *                      If the segment length is less than the maxseg length,
 *                      AND (the area of the triangle defined by the node and
 *                      its two neighbors is less than the limit OR the area
 *                      is decreasing in size OR the segment is too small to
 *                      bisect) we don't refine it.
 */
                       if ((segLen < param->maxSeg) &&
                           ((area2 < areaMax2) ||
                            (segLen < param->minSeg * 2.0) ||
                            (darea2dt < 0.0))) {
                          continue;
                       }
 
/*
 *                      Add the segment to list of segments to be refined.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_REFINE_LIST)
#endif
                        {
                            int i, duplicateEntry = 0;

/*
 *                          Just a quick sanity check, but make sure there
 *                          is still space on the list...
 */
                            if (refineListEnts < refineListSize) {
/*
 *                              We're not requiring the segment to be owned by
 *                              the current node, so we need to veryify the
 *                              segment has not already been placed on the
 *                              list by another thread.
 */
                                for (i = 0; i < refineListEnts; i++) {
                                    RefineList_t *entry;

                                    entry = &refineList[i];

                                    if ((entry->node1Tag.domainID ==
                                         nbr->myTag.domainID)        &&
                                        (entry->node1Tag.index    ==
                                         nbr->myTag.index   )        &&
                                        (entry->node2Tag.domainID ==
                                         node->myTag.domainID)       &&
                                        (entry->node2Tag.index    ==
                                         node->myTag.index   )) {
                                        duplicateEntry = 1;
                                        break;
                                    }
                                }

                                if (!duplicateEntry) {

                                    refineList[refineListEnts].node1Tag =
                                            node->myTag;
                                    refineList[refineListEnts].node2Tag =
                                            nbr->myTag;
                                    refineList[refineListEnts].status  =
                                            REMESH_PENDING;
                                    VECTOR_COPY(
                                            refineList[refineListEnts].vector,
                                            vec);

                                    refineListEnts += 1;
                                }

                            }  /* end if (refineListEnts < refineListSize) */

                        }  /* end "omp critical" section */

                    }  /* for (seg = 0; seg < 2; ...) */

                } else {
                    int nbrIndex;
/*
 *                  For nodes with other than exactly two arms, we
 *                  just refine any arm exceeding the max allowable
 *                  segment length, but only do the refinement from the
 *                  domain "owning" the segment.
 *
 *                  A "split" is only considered a global operation
 *                  to be sent to remote domains if the segment being
 *                  split spans domains.
 */
                    for (nbrIndex = 0; nbrIndex < node->numNbrs; nbrIndex++) {
        
                        nbr1 = GetNeighborNode(home, node, nbrIndex);
        
                        if (nbr1 == (Node_t *)NULL) {
                            printf("WARNING: Neighbor not found at %s line %d\n",
                                   __FILE__, __LINE__);
                            continue;
                        }

                        if (OrderNodes(node, nbr1) >= 0) {
                            continue;
                        }

                        if (!DomainOwnsSeg(home, OPCLASS_REMESH,
                                           thisDomain, &nbr1->myTag)) {
                            continue;
                        }
        
                        vec1[X] = nbr1->x - node->x;
                        vec1[Y] = nbr1->y - node->y;
                        vec1[Z] = nbr1->z - node->z;
        
                        ZImage(param, &vec1[X], &vec1[Y], &vec1[Z]);
        
                        r1 = DotProduct(vec1, vec1);
        
/*
 *                      If the segment is too long, add it to the list of
 *                      segments to be refined
 */
                        if (r1 > maxSeg2) {
/*
 *                          If there is no more room on the list, reallocate
 *                          the buffer with some more entries.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_REFINE_LIST)
#endif
                            {
                                int i, duplicateEntry = 0;
/*
 *                              Just a quick sanity check, but make sure there
 *                              is still space on the list...
 */
                                if (refineListEnts < refineListSize) {
/*
 *                                  We're not requiring the segment to be
 *                                  owned by the current node, so we need to
 *                                  veryify the segment has not already been
 *                                  placed on the list by another thread.
 */
                                    for (i = 0; i < refineListEnts; i++) {
                                        RefineList_t *entry;

                                        entry = &refineList[i];

                                        if ((entry->node1Tag.domainID ==
                                             nbr1->myTag.domainID)        &&
                                            (entry->node1Tag.index    ==
                                             nbr1->myTag.index)           &&
                                            (entry->node2Tag.domainID ==
                                             node->myTag.domainID)        &&
                                            (entry->node2Tag.index    ==
                                             node->myTag.index)) {
                                            duplicateEntry = 1;
                                            break;
                                        }
                                    }

                                    if (!duplicateEntry) {

                                        refineList[refineListEnts].node1Tag =
                                                node->myTag;
                                        refineList[refineListEnts].node2Tag =
                                                nbr1->myTag;
                                        refineList[refineListEnts].status   =
                                                REMESH_PENDING;

                                        VECTOR_COPY(refineList[refineListEnts].vector, vec1);

                                        refineListEnts += 1;
                                    }
                                }

                            }  /* end "omp critical" section */
                        }
                    }
                }

            }  /* end loop over all nodes */

        }  /* end "omp parallel" section */

/*
 *      Next phase is to loop through the segments selected for
 *      refinement and bisect them.
 *
 *      Note: The restrictions implemented when adding segments to the list
 *            should be sufficient that we can refine every segment that
 *            made it onto the list.
 */
#ifdef _OPENMP
#pragma omp parallel for if (refineListEnts > 1)
#endif
        for (int segIndex=0; segIndex < refineListEnts; segIndex++) {
            int    splitStatus=0, globalOp;
            real8  vec1[3];
            real8  newPos[3], newVel[3];
            real8  nodePos[3], nodeVel[3];
            real8  f0Seg1[3], f1Seg1[3], f0Seg2[3], f1Seg2[3];
            Node_t *node1, *node2, *splitNode1, *splitNode2;

            node1 = GetNodeFromTag(home, refineList[segIndex].node1Tag);
            node2 = GetNodeFromTag(home, refineList[segIndex].node2Tag);
            VECTOR_COPY(vec1, refineList[segIndex].vector);

/*
 *          Splitting a segment screws with the associated nodes'
 *          segment lists, velocities, etc.  In order to prevent multiple
 *          threads from messing with the same nodes simultaneously
 *          lock the two nodes. 
 */
            LOCK(&node1->nodeLock);
            LOCK(&node2->nodeLock);

/*
 *          Get an estimate of forces that would change due to refining
 *          the segment.
 */
            newPos[X] = node1->x + (vec1[X] * 0.5);
            newPos[Y] = node1->y + (vec1[Y] * 0.5);
            newPos[Z] = node1->z + (vec1[Z] * 0.5);
        
            newVel[X] = (node1->vX + node2->vX) * 0.5;
            newVel[Y] = (node1->vY + node2->vY) * 0.5;
            newVel[Z] = (node1->vZ + node2->vZ) * 0.5;

/*
 *          Because the EstRefinementForces() call looks at the segment
 *          lists of <node1> and <node2> to find the connections between
 *          them, this call must be within the critical section or we
 *          risk another thread bisecting a segment attached to one of
 *          the nodes and modifying the node's segment list while
 *          EstRefinementForces is looking at the list.
 */
            EstRefinementForces(home, node1, node2, newPos, vec1,
                                f0Seg1, f1Seg1, f0Seg2, f1Seg2);

            FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);

/*
 *          This should be a global operation distributed out to remote
 *          domains only if the segment spans domains.
 */
            globalOp = (node2->myTag.domainID != node1->myTag.domainID);
        
            nodePos[X] = node1->x;
            nodePos[Y] = node1->y;
            nodePos[Z] = node1->z;
        
            nodeVel[X] = node1->vX;
            nodeVel[Y] = node1->vY;
            nodeVel[Z] = node1->vZ;
        
/*
 *          Need a critical section now because we'll be inserting a
 *          node and changing connections between nodes...
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_REFINE_SEG)
#endif
            {
/*
 *              When refining a segment, we always move exactly one arm
 *              from the original node to the node being created...
 */
                int armIndex = GetArmID(node1, node2);
                int armCount = 1;

                if (armIndex >= 0) {

                    splitStatus = SplitNode(home, OPCLASS_REMESH,
                                            node1, nodePos, newPos,
                                            nodeVel, newVel,
                                            armCount, &armIndex,
                                            globalOp, &splitNode1,
                                            &splitNode2, 0);

                    if (splitStatus == SPLIT_SUCCESS) {
                        localRefineCnt++;
                    }
                }

            }  /* end "omp critical" section */

/*
 *          It's safe to allow another thread to mess with the two
 *          nodes again.
 */
            UNLOCK(&node1->nodeLock);
            UNLOCK(&node2->nodeLock);

/*
 *          We refined the segment, so we need to save some more info related
 *          to this refine event.  Don't need to lock anything here since 
 *          the refine list will not be realloced here nor will any other
 *          thread be updating this entry.
 */
            if (splitStatus == SPLIT_SUCCESS) {
                refineList[segIndex].status = REMESH_COMMITTED;
                refineList[segIndex].newNodeTag = splitNode2->myTag;
                VECTOR_COPY(refineList[segIndex].f0Seg1, f0Seg1);
                VECTOR_COPY(refineList[segIndex].f1Seg1, f1Seg1);
                VECTOR_COPY(refineList[segIndex].f0Seg2, f0Seg2);
                VECTOR_COPY(refineList[segIndex].f1Seg2, f1Seg2);
            }
        }

/*
 *      Last step is to update nodal forces and velocities for
 *      nodes involved in the the refine events
 */
#ifdef _OPENMP
#pragma omp parallel for if (refineListEnts > 1)
#endif
        for (int segIndex=0; segIndex < refineListEnts; segIndex++) 
        {
            MobArgs_t mobArgs;
 
/*
 *          I don't think any of the refine events should have failed
 *          in the loop above, but just in case...
 */
            if (refineList[segIndex].status != REMESH_COMMITTED) {
                continue;
            }

            Node_t *node1   = GetNodeFromTag(home, refineList[segIndex].node1Tag  );
            Node_t *node2   = GetNodeFromTag(home, refineList[segIndex].node2Tag  );
            Node_t *newNode = GetNodeFromTag(home, refineList[segIndex].newNodeTag);

/*
 *          Mark the force and velocity data for some nodes as obsolete so
 *          that more accurate forces will be recalculated at the beginning
 *          of the next timestep.
 */
            MarkNodeForceObsolete(home, node1  );
            MarkNodeForceObsolete(home, node2  );
            MarkNodeForceObsolete(home, newNode);
        
/*
 *          Using the force estimates calculated in the previous loop, 
 *          reset nodal forces on all nodes involved in the refine event.
 *          Need to lock each node as we update it so no other thread
 *          updates the same node simultaneously
 */
            real8  f0Seg1[3];  VECTOR_COPY(f0Seg1, refineList[segIndex].f0Seg1);
            real8  f1Seg1[3];  VECTOR_COPY(f1Seg1, refineList[segIndex].f1Seg1);
            real8  f0Seg2[3];  VECTOR_COPY(f0Seg2, refineList[segIndex].f0Seg2);
            real8  f1Seg2[3];  VECTOR_COPY(f1Seg2, refineList[segIndex].f1Seg2);

/*
 *          Node 1
 */
            LOCK(&node1->nodeLock);
            ResetSegForces(home, node1, &newNode->myTag, f0Seg1[X], f0Seg1[Y], f0Seg1[Z], 1);
            (void)EvaluateMobility(home, node1, &mobArgs);
            UNLOCK(&node1->nodeLock);

/*
 *          New node
 */
            LOCK(&newNode->nodeLock);
            ResetSegForces(home, newNode, &node1->myTag, f1Seg1[X], f1Seg1[Y], f1Seg1[Z], 1); 
            ResetSegForces(home, newNode, &node2->myTag, f0Seg2[X], f0Seg2[Y], f0Seg2[Z], 1);
            (void)EvaluateMobility(home, newNode, &mobArgs);
            UNLOCK(&newNode->nodeLock);

/*
 *          Node 2
 */
            LOCK(&node2->nodeLock);
            ResetSegForces(home, node2, &newNode->myTag, f1Seg2[X], f1Seg2[Y], f1Seg2[Z], 1);
            (void)EvaluateMobility(home, node2, &mobArgs);
            UNLOCK(&node2->nodeLock);
        
/*
 *          When debugging, dump some info on refine events taking
 *          place and the nodes involved
 */
#ifdef DEBUG_TOPOLOGY_CHANGES

            if ((dbgDom < 0) || (dbgDom == home->myDomain)) 
            {
                printf("Remesh/refine:  (%d,%d)--(%d,%d) ==> "
                       "(%d,%d)--(%d,%d)--(%d,%d)\n",
                       node1->myTag.domainID, node1->myTag.index,
                       node2->myTag.domainID, node2->myTag.index,
                       node1->myTag.domainID, node1->myTag.index,
                       newNode->myTag.domainID, newNode->myTag.index,
                       node2->myTag.domainID, node2->myTag.index);
                PrintNode(node1);
                PrintNode(newNode);
                PrintNode(node2);
            }
#endif
        }

        
#ifdef DEBUG_LOG_MESH_REFINE
#ifdef PARALLEL
        int globalRefineCnt = 0;
        MPI_Reduce(&localRefineCnt, &globalRefineCnt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        int globalRefineCnt = localRefineCnt;
#endif
        if (home->myDomain == 0) {
            printf("  Remesh: refine count = %d\n", globalRefineCnt);
        }
#endif

/*
 *      Free the refine list now that we no longer need it.
 */
        if (refineList != (RefineList_t *)NULL) {
            free(refineList);
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    RemeshRule_2
 *      Description: Base function which invokes all subroutines needed
 *                   for handling dislocation segment rediscretiztation
 *                   using the second version of remesh.
 *
 *------------------------------------------------------------------------*/
void RemeshRule_2(Home_t *home)
{

#ifdef DEBUG_TOPOLOGY_DOMAIN
        dbgDom = DEBUG_TOPOLOGY_DOMAIN;
#else
        dbgDom = -1;
#endif

        MeshCoarsen(home);
        MeshRefine(home);

        return;
}
