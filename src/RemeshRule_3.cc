/*****************************************************************************
 *
 *      Module:         RemeshRule_3.c
 *      Description:    This module contains functions specific to
 *                      version 3 remesh for rediscretizting the dislocation
 *                      segments.  This version is very similar to
 *                      the version 2 remesh with some code added to
 *                      MeshRefine() to better calculate the position
 *                      at which to place new nodes when bisecting
 *                      the arms of a 2-node. (New nodes will be placed
 *                      on the arc containing the original node and
 *                      its two neighbors rather than just at the
 *                      midpoints of the original segments.)
 *                      
 *      Included functions:
 *              MeshRefine()
 *              RemeshRule_3()
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
 *
 *      Description: Check all the local dislocation segments looking
 *                   for any which should be discretized more finely.
 *                   All such segments found are added to a 'refine list'
 *                   after which the identified segments will be re-
 *                   discretized by bisecting the segment with a new node.
 *
 *------------------------------------------------------------------------*/
static void MeshRefine(Home_t *home)
{

/*
 *      Initialize some of the constants we'll need during this call.
 */
        int      thisDomain = home->myDomain;
        Param_t *param      = home->param;

        real8    areaMax    = param->remeshAreaMax;
        real8    areaMaxSq  = areaMax * areaMax;
        real8    maxSegSq   = param->maxSeg * param->maxSeg;
        real8    delta      = 1.0e-16;
        real8    eps        = 1.0e-12;

        int      localRefineCnt = 0;
        
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
        RefineList_t *refineList     = (RefineList_t *) NULL;

        if (refineListSize > 0) {
            refineList = (RefineList_t *)malloc(sizeof(RefineList_t) *
                                                refineListSize);
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
                    real8 p1, p2, p3;
                    real8 mid1[3], mid2[3];
                    real8 mv3[3], rhs[3];
                    real8 mat[3][3], invMat[3][3];
        
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
 *                  However, only calculate the effects of the velocity of
 *                  the primary node on that area (i.e. treat the neighboring
 *                  nodes as if they had zero velocity).  If we don't do
 *                  this, there are situations where although the area is
 *                  increasing the primary node is moving *toward* the
 *                  other two nodes and placing new nodes on an arc
 *                  containing the 3 nodes turns out to be the incorrect
 *                  placement.
 */
                    dr1dt = ((vec1[X] * -node->vX) + (vec1[Y] * -node->vY) +
                             (vec1[Z] * -node->vZ)) / (r1 + delta);

                    dr2dt = ((vec2[X] * -node->vX) + (vec2[Y] * -node->vY) +
                             (vec2[Z] * -node->vZ)) / (r2 + delta);

                    dr3dt = 0.0;

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
 *                  Set up some stuff for calculating the arc containing
 *                  the node and its two neighbors so we can position the
 *                  new nodes on the arc.
 */
                    mid1[X] = node->x + (vec1[X] * 0.5);
                    mid1[Y] = node->y + (vec1[Y] * 0.5);
                    mid1[Z] = node->z + (vec1[Z] * 0.5);

                    mid2[X] = node->x + (vec2[X] * 0.5);
                    mid2[Y] = node->y + (vec2[Y] * 0.5);
                    mid2[Z] = node->z + (vec2[Z] * 0.5);

                    p1 = (vec1[X]*mid1[X] + vec1[Y]*mid1[Y] + vec1[Z]*mid1[Z]);
                    p2 = (vec2[X]*mid2[X] + vec2[Y]*mid2[Y] + vec2[Z]*mid2[Z]);

                    cross(vec1, vec2, mv3);

                    p3 = (mv3[X]*node->x + mv3[Y]*node->y + mv3[Z]*node->z);

                    mat[0][X] = vec1[X];
                    mat[0][Y] = vec1[Y];
                    mat[0][Z] = vec1[Z];

                    mat[1][X] = vec2[X];
                    mat[1][Y] = vec2[Y];
                    mat[1][Z] = vec2[Z];

                    mat[2][X] = mv3[X];
                    mat[2][Y] = mv3[Y];
                    mat[2][Z] = mv3[Z];

                    rhs[X] = p1;
                    rhs[Y] = p2;
                    rhs[Z] = p3;

/*
 *                  Check both segments (in the preferred order) and 
 *                  determine if we should bisect one.
 */
                    for (segID = 0; segID < 2; segID++) {
                        int   segLen;
                        int   burgIndex;
                        real8 *vec, *mid;
                        real8 radius, norm, np, tmp;
                        real8 center[3], dir[3], burg[3], nplane[3];
                        real8 diff[3], newPos[3], tmp3[3];

                        if (splitSegList[segID] == 1) {
                            if (splitOK[0] == 0) {
                                continue;
                            }
                            nbr = nbr1;
                            vec = vec1;
                            segLen = r1;
                            mid = mid1;
                            burgIndex = 0;
                        } else {
                            if (splitOK[1] == 0) {
                                continue;
                            }
                            nbr = nbr2;
                            vec = vec2;
                            segLen = r2;
                            mid = mid2;
                            burgIndex = 1;
                        }

/*
 *                      If the segment length is less than the maxseg length,
 *                      AND (the area of the triangle defined by the node and
 *                      its two neighbors is less than the limit OR the area
 *                      is decreasing in size OR the segment is too small to
 *                      bisect) we don't refine it.
 */
                        if ((segLen < param->maxSeg) &&
                            ((area2 < areaMaxSq) ||
                             (segLen < param->minSeg * 2.0) ||
                             (darea2dt < 0.0))) {
                            continue;
                        }

/*
 *                      Now need to calculate the position along the arc at
 *                      which the new node should be added.
 */
                        if (area2 > areaMaxSq) {
                            if (Matrix33_Inverse(invMat, mat) < 0 ) 
                            { Fatal("%s::%s(%d) : Cannot invert 3X3 matrix!", __FILE__, __func__, __LINE__ ); }

                            Matrix33Vector3Multiply(invMat, rhs, center);

                            tmp3[X] = node->x - center[X];
                            tmp3[Y] = node->y - center[Y];
                            tmp3[Z] = node->z - center[Z];

                            radius = Normal(tmp3);

                            tmp3[X] = mid[X] - center[X];
                            tmp3[Y] = mid[Y] - center[Y];
                            tmp3[Z] = mid[Z] - center[Z];

                            norm = Normal(tmp3);

                            dir[X] = tmp3[X] / norm;
                            dir[Y] = tmp3[Y] / norm;
                            dir[Z] = tmp3[Z] / norm;

                            burg[X] = node->burgX[burgIndex];
                            burg[Y] = node->burgY[burgIndex];
                            burg[Z] = node->burgZ[burgIndex];

                            cross(burg, vec, nplane);
                            np = Normal(nplane);

                            if (np > eps) {
                                nplane[X] /= np;
                                nplane[Y] /= np;
                                nplane[Z] /= np;
                            } else {
                                nplane[X] = 0.0;
                                nplane[Y] = 0.0;
                                nplane[Z] = 0.0;
                            }

                            diff[X] = (radius*dir[X]) - (mid[X]-center[X]);
                            diff[Y] = (radius*dir[Y]) - (mid[Y]-center[Y]);
                            diff[Z] = (radius*dir[Z]) - (mid[Z]-center[Z]);

                            tmp = (diff[X]*nplane[X] +
                                   diff[Y]*nplane[Y] +
                                   diff[Z]*nplane[Z]);

                            diff[X] = diff[X] - (tmp * nplane[X]);
                            diff[Y] = diff[Y] - (tmp * nplane[Y]);
                            diff[Z] = diff[Z] - (tmp * nplane[Z]);

                            newPos[X] = mid[X] + diff[X];
                            newPos[Y] = mid[Y] + diff[Y];
                            newPos[Z] = mid[Z] + diff[Z];
                        } else {
                            newPos[X] = node->x + (vec[X] * 0.5);
                            newPos[Y] = node->y + (vec[Y] * 0.5);
                            newPos[Z] = node->z + (vec[Z] * 0.5);
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
                                    refineList[refineListEnts].status =
                                            REMESH_PENDING;
                                    refineList[refineListEnts].haveNewPos = 1;

                                    VECTOR_COPY(
                                            refineList[refineListEnts].vector,
                                            vec);
                                    VECTOR_COPY(
                                            refineList[refineListEnts].newPos,
                                            newPos);
    
                                    refineListEnts += 1;
                                }
                            }

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
                        if (r1 > maxSegSq) {
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
                                        refineList[refineListEnts].haveNewPos = 0;

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
        for (int segIndex = 0; segIndex < refineListEnts; segIndex++) {
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
 *          For some of the segments we'll be placing the new node
 *          on the arc containing the original 3 nodes.  In such cases
 *          we calculated the new position earlier, so use the known
 *          position.  Otherwise calculate it here.
 */
            if (refineList[segIndex].haveNewPos) {
                VECTOR_COPY(newPos, refineList[segIndex].newPos);
            } else {
                newPos[X] = node1->x + (vec1[X] * 0.5);
                newPos[Y] = node1->y + (vec1[Y] * 0.5);
                newPos[Z] = node1->z + (vec1[Z] * 0.5);
            }
        
            newVel[X] = (node1->vX + node2->vX) * 0.5;
            newVel[Y] = (node1->vY + node2->vY) * 0.5;
            newVel[Z] = (node1->vZ + node2->vZ) * 0.5;
        
/*
 *          Get an estimate of forces that would change due to refining
 *          the segment.
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
        for (int segIndex = 0; segIndex < refineListEnts; segIndex++) 
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
            MarkNodeForceObsolete(home, node1);
            MarkNodeForceObsolete(home, node2);
            MarkNodeForceObsolete(home, newNode);
        
/*
 *          Using the force estimates calculated in the previous loop, 
 *          reset nodal forces on all nodes involved in the refine event.
 *          Need to lock each node as we update it so no other thread
 *          updates the same node simultaneously
 */
            real8  f0Seg1[3]; VECTOR_COPY(f0Seg1, refineList[segIndex].f0Seg1);
            real8  f1Seg1[3]; VECTOR_COPY(f1Seg1, refineList[segIndex].f1Seg1);
            real8  f0Seg2[3]; VECTOR_COPY(f0Seg2, refineList[segIndex].f0Seg2);
            real8  f1Seg2[3]; VECTOR_COPY(f1Seg2, refineList[segIndex].f1Seg2);

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

            if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
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
        int globalRefineCnt=0;
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
 *      Function:    RemeshRule_3
 *      Description: Base function which invokes all subroutines needed
 *                   for handling dislocation segment rediscretization
 *                   using the third version of remesh.
 *
 *------------------------------------------------------------------------*/
void RemeshRule_3(Home_t *home)
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
