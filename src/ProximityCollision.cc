/*****************************************************************************
 *
 *      Module:         ProximityCollision.c
 *      Description:    This module contains various functions used
 *                      for detecting various types of collisions
 *                      (segment/segment, node/node, zipping) and
 *                      dealing with those collisions.  These functions
 *                      are specific to the type 1 collision handling
 *                      which is based primarily on proximity criteria.
 *
 *                      Each domain handles local collisions and
 *                      distributes the necessary topological changes
 *                      to remote domains via the same mechanism
 *                      employed by remesh.
 *
 *                      NOTE: Certain strict rules govern what topological
 *                      changes are permitted based on noda and segment
 *                      ownership.  See comments at the beginning of
 *                      the module Topology.c for details of the
 *                      rule.  Additional restrictions may be implemented
 *                      for the collision handling; see code below
 *                      for details of any such restrictions.
 *
 *
 *      Included functions:
 *
 *          FindIdealCollisionPoint()
 *          ProximityCollisions()
 *
 *****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "Comm.h"
#include "Mobility.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif /* ifdef _ARLFEM */

static int dbgDom;


/*---------------------------------------------------------------------------
 *
 *      Function:       FindIdealCollisionPoint
 *      Description:    Given two points and the glide plane containing
 *                      each, select the best location at which the two
 *                      points could be merged.
 *
 *      Arguments:
 *           pos1, pos2    Positions of the two points which will be collided
 *           gp1, gp2      Glide planes containing the respective points.
 *           cPoint        Array in which to return to the caller the
 *                         selected collision point.
 *
 *-------------------------------------------------------------------------*/
static void FindIdealCollisionPoint(Home_t *home, real8 pos1[3],
                                    real8 pos2[3], real8 gp1[3], real8 gp2[3],
                                    real8 cPoint[3])
{
        int   i, j, matOrder;
        real8 pos1Dotgp1, pos2Dotgp2;
        real8 newPlaneCond = 0.875;
        real8 vec[5], result[5];
        real8 mat[5][5], invmat[5][5];


        memset(mat, 0, sizeof(mat));

        for (i = 0; i < 3; i++) {
            vec[i] = 0.5 * (pos1[i] + pos2[i]);
        }

/*
 *      If strict glide planes are not being enforced, we
 *      can use 'fuzzy' glide planes so that we don't
 *      throw nodes long distances if the planes are
 *      close to parallel but don't quite meet the criteria
 *      for being parallel.
 *
 *      <factor> determines how 'fuzzy' planes are.  The smaller the
 *      value, the fuzzier the plane constraints.
 */
        if (home->param->enforceGlidePlanes == 0) {
            real8 factor = 5.0;
            real8 gp1TM[3][3], gp2TM[3][3];

            matOrder = 3;

            pos1Dotgp1 = DotProduct(pos1, gp1);
            pos2Dotgp2 = DotProduct(pos2, gp2);

            for (i = 0; i < 3; i++) {
                vec[i] += factor * (pos1Dotgp1 * gp1[i] + pos2Dotgp2 * gp2[i]);
            }

            Vec3TransposeAndMult(gp1, gp1TM);
            Vec3TransposeAndMult(gp2, gp2TM);

            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    mat[i][j] = (real8)(i==j) +
                                factor * (gp1TM[i][j] + gp2TM[i][j]);
                }
            }
        } else {
/*
 *          Strict glide planes *are* being enforced, so we have to find
 *          a point that satisfies both glide plane constraints.
 */
            real8 planeTest;

            planeTest = DotProduct(gp1, gp2);
            planeTest = planeTest * planeTest / 
                        (DotProduct(gp1,gp1) * DotProduct(gp2,gp2));


            if (planeTest < newPlaneCond) {
/*
 *              There are two separable planes
 */
                matOrder = 5;

                vec[3] = DotProduct(pos1, gp1);
                vec[4] = DotProduct(pos2, gp2);

                for (i = 0; i < 3; i++) {

                    mat[i][i] = 1.0;

                    mat[i][3] = gp1[i];
                    mat[i][4] = gp2[i];

                    mat[3][i] = gp1[i];
                    mat[4][i] = gp2[i];
                }

            } else {
/*
 *              There is only 1 unique plane
 */
                matOrder = 4;

                vec[3] = DotProduct(pos1, gp1);

                for (i = 0; i < 3; i++) {

                    mat[i][i] = 1.0;

                    mat[i][3] = gp1[i];
                    mat[3][i] = gp1[i];
                }
            }
        }

/*
 *      Calculate the collision point
 */

        if (!MatrixInvert((real8 *)mat, (real8 *)invmat, matOrder, 5)) {
            Fatal("FindIdealCollisionPoint(): Unable to invert matrix!");
        }

        MatrixMult((real8 *)invmat, matOrder, matOrder, 5,
                   (real8 *)vec, 1, 1, 
                   (real8 *)result, 1);

        VECTOR_COPY(cPoint, result);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       ProximityCollisions
 *      Description:    Loop though all native nodes, identify segments
 *                      or nodes that should be collided (or zipped) by
 *                      the current domain and handle all such collisions.
 *
 *                      Note the following restrictions on collisions:
 *
 *                      - A 'pinned' node be deleted during a collision
 *                      - A node which has been the target of a merge
 *                        (this timestep) may not be deleted during any 
 *                        other type of collision in this timestep.
 *
 *      All nodes with arms owned by another domain have previously
 *      been marked as non-deletable nodes for this cycle.
 *
 *-------------------------------------------------------------------------*/
void ProximityCollisions(Home_t *home)
{
        int     i, j, k, q, arm12, arm21, arm34, arm43;
        int     thisDomain, splitStatus, mergeStatus;
        int     armAB, armBA;
        int     globalOp = 1, didCollision;
        int     close2node1, close2node2, close2node3, close2node4;
        int     splitSeg1, splitSeg2;
        int     cell2Index, nbrCell2Index, nextIndex;
        int     cell2X, cell2Y, cell2Z, cx, cy, cz;
        int     retCode;
        real8   mindist2, dist2, ddist2dt, L1, L2, eps, half;
        real8   x1, y1, z1, vx1, vy1, vz1;
        real8   x2, y2, z2, vx2, vy2, vz2;
        real8   x3, y3, z3, vx3, vy3, vz3;
        real8   x4, y4, z4, vx4, vy4, vz4;
        real8   seg1Lx, seg1Ly, seg1Lz;
        real8   seg2Lx, seg2Ly, seg2Lz;
        real8   newx, newy, newz;
        real8   newvx, newvy, newvz;
        real8   p0[3], p1[3], pnew[3], burg1[3], burg2[3];
        real8   oldfp0[3], oldfp1[3];
        real8   f0seg1[3], f1seg1[3], f0seg2[3], f1seg2[3];
        real8   nodeVel[3], newNodeVel[3];
        real8   newPos[3], newVel[3];
        real8   vec1[3], vec2[3];
        Tag_t   oldTag1, oldTag2;
        Node_t  *node1, *node2, *node3, *node4, *tmpNbr;
        Node_t  *mergenode1, *mergenode2, *targetNode;
        Node_t  *splitNode1, *splitNode2;
        Param_t *param;
        MobArgs_t mobArgs;
#ifdef _ARLFEM
        int     resetSurfaceProperties;
        real8   femSurfaceNorm[3];
#endif

        thisDomain = home->myDomain;
        param      = home->param;

        mindist2 = param->rann * param->rann;
        eps      = 1.0e-12;
        half     = 0.5;

        int localCollisionCnt = 0;

#ifdef DEBUG_TOPOLOGY_DOMAIN
        dbgDom = DEBUG_TOPOLOGY_DOMAIN;
#else
        dbgDom = -1;
#endif

        TimerStart(home, COLLISION_HANDLING);

/*
 *      Start looping through native nodes looking for segments to collide...
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
            if (node1->flags & NO_COLLISIONS) continue;

            didCollision = 0;

/*
 *          Loop through all cell2s neighboring the node.  Only
 *          nodes in these neighboring cell2s are candidates for
 *          collisions.
 *
 *          NOTE: All nodes are assigned cell2 membership prior to
 *          entering collision handling.  However, nodes added during
 *          node separation and collision handling are not given
 *          cell2 membership for this cycle.  In most cases this
 *          is not a problem since those nodes have usually been
 *          exempted from subsequent collisions this cycle.  There
 *          are certain circumstances, though, in which the exemptions
 *          have not been set.  If we encounter one such node that
 *          has no cell2 membership, skip it for this cycle, or bad
 *          things happen.
 */
            cell2Index = node1->cell2Idx;
            if (cell2Index < 0) {
                continue;
            }

            DecodeCell2Idx(home, cell2Index, &cell2X, &cell2Y, &cell2Z);

            for (cx = cell2X - 1; cx <= cell2X + 1; cx++) {
             for (cy = cell2Y - 1; cy <= cell2Y + 1; cy++) {
              for (cz = cell2Z - 1; cz <= cell2Z + 1; cz++) {

                nbrCell2Index = EncodeCell2Idx(home, cx, cy, cz);

/*
 *              Loop though all nodes in the neighbor cell2
 */
                nextIndex = home->cell2[nbrCell2Index];

                while (nextIndex >= 0) {

                    if (didCollision) break;

                    node3 = home->cell2QentArray[nextIndex].node;
                    nextIndex = home->cell2QentArray[nextIndex].next;

                    if (node3 == (Node_t *)NULL) continue;
                    if (node3->flags & NO_COLLISIONS) continue;

                    if (CollisionNodeOrder(home, &node1->myTag,
                                           &node3->myTag) >= 0) {
                        continue;
                    }

/*
 *                  Loop over all arms of node1.  Skip any arms that
 *                  terminate at node3 (those hinge arms will be dealt
 *                  with later)
 */
                    for (arm12 = 0; arm12 < node1->numNbrs; arm12++) {

                        if (didCollision) break;

                        node2 = GetNodeFromTag(home, node1->nbrTag[arm12]);

                        if (node2 == (Node_t *)NULL) continue;
                        if (node2->flags & NO_COLLISIONS) continue;

                        if (CollisionNodeOrder(home, &node1->myTag,
                                               &node2->myTag) > 0) {
                            continue;
                        }

                        if ((node2->myTag.domainID == node3->myTag.domainID) &&
                            (node2->myTag.index    == node3->myTag.index   )) {
                            continue;
                        }

/*
 *                      Segment node1/node2 may only be used in a collision
 *                      if the segment is owned by the current domain.
 */
                        if (!DomainOwnsSeg(home, OPCLASS_COLLISION,
                                           thisDomain, &node2->myTag)) {
                            continue;
                        }

/*
 *                      Loop over all arms of node3.  Skip any segments that
 *                      terminate at node2 (those hinge arms will be dealt
 *                      with later), and skip any segments not owned by
 *                      node3 -- we'll deal with those segments when we
 *                      hit the node owning the segment
 *                      
 */
                        for (arm34 = 0; arm34 < node3->numNbrs; arm34++) {
                            real8 dx, dy, dz;

                            if (didCollision) break;

                            node4 = GetNodeFromTag(home, node3->nbrTag[arm34]);

                            if (node4 == (Node_t *)NULL) continue;
                            if (node4->flags & NO_COLLISIONS) continue;

                            if ((node4->myTag.domainID==node2->myTag.domainID)&&
                                (node4->myTag.index   ==node2->myTag.index)) {
                                continue;
                            }
            
                            if (CollisionNodeOrder(home, &node3->myTag,
                                                   &node4->myTag) > 0) {
                                continue;
                            }

/*
 *                          At this point, segment node3/node4 is owned by
 *                          node3.  If node3 is not native to this domain,
 *                          the segment may not be used in a collision since
 *                          the domain doing to collision must own both
 *                          segments.
 */
                            if (node3->myTag.domainID != thisDomain) {
                                continue;
                            }

                            x1 = node1->x; y1 = node1->y; z1 = node1->z;
                            x2 = node2->x; y2 = node2->y; z2 = node2->z;
                            x3 = node3->x; y3 = node3->y; z3 = node3->z;
                            x4 = node4->x; y4 = node4->y; z4 = node4->z;
    
                            vx1 = node1->vX; vy1 = node1->vY; vz1 = node1->vZ;
                            vx2 = node2->vX; vy2 = node2->vY; vz2 = node2->vZ;
                            vx3 = node3->vX; vy3 = node3->vY; vz3 = node3->vZ;
                            vx4 = node4->vX; vy4 = node4->vY; vz4 = node4->vZ;
    
                            PBCPOSITION(param, x1, y1, z1, &x2, &y2, &z2);
                            PBCPOSITION(param, x1, y1, z1, &x3, &y3, &z3);
                            PBCPOSITION(param, x3, y3, z3, &x4, &y4, &z4);
    
/*
 *                          It is possible to have a zero-length segment
 *                          (created by a previous collision).  If we find
 *                          such a segment, do not try to use it in any
 *                          subsequent collisions.
 */
                            dx = x1 - x2;
                            dy = y1 - y2;
                            dz = z1 - z2;

                            if ((dx*dx + dy*dy + dz*dz) < 1.0e-20) {
                                continue;
                            }

                            dx = x3 - x4;
                            dy = y3 - y4;
                            dz = z3 - z4;

                            if ((dx*dx + dy*dy + dz*dz) < 1.0e-20) {
                                continue;
                            }

/*
 *                          Find the minimum distance between the two segments
 *                          and determine if they should be collided.
 */
                            GetMinDist(x1, y1, z1, vx1, vy1, vz1,
                                       x2, y2, z2, vx2, vy2, vz2,
                                       x3, y3, z3, vx3, vy3, vz3,
                                       x4, y4, z4, vx4, vy4, vz4,
                                       &dist2, &ddist2dt, &L1, &L2);
    
/*
 *                          Requiring the time rate of change of dist2 to be 
 *                          larger than zero excludes connected links.
 *                          Non-endpoint condition looks for co-planar line
 *                          segments that are already intersecting whose
 *                          collision was not anticipated due to a large
 *                          timestep
 */
                            if (((dist2<mindist2) && (ddist2dt<-eps)) ||
                                (dist2<eps)) {
                                real8 p1[3], p2[3];
                                real8 cPoint[3], mPoint1[3], mPoint2[3];

                                burg1[X] = node1->burgX[arm12];
                                burg1[Y] = node1->burgY[arm12];
                                burg1[Z] = node1->burgZ[arm12];
    
                                burg2[X] = node3->burgX[arm34];
                                burg2[Y] = node3->burgY[arm34];
                                burg2[Z] = node3->burgZ[arm34];
    
/*
 *                              Calculate the ideal collision point between
 *                              the two segments using the 'closest' points
 *                              identified above.
 */
                                mPoint1[X] = x1 * (1.0-L1) + x2 * L1;
                                mPoint1[Y] = y1 * (1.0-L1) + y2 * L1;
                                mPoint1[Z] = z1 * (1.0-L1) + z2 * L1;

                                mPoint2[X] = x3 * (1.0-L2) + x4 * L2;
                                mPoint2[Y] = y3 * (1.0-L2) + y4 * L2;
                                mPoint2[Z] = z3 * (1.0-L2) + z4 * L2;

                                if (param->enforceGlidePlanes) {
                                    p1[X] = node1->nx[arm12];
                                    p1[Y] = node1->ny[arm12];
                                    p1[Z] = node1->nz[arm12];
                                    p2[X] = node3->nx[arm34];
                                    p2[Y] = node3->ny[arm34];
                                    p2[Z] = node3->nz[arm34];
                                } else {
                                    real8 vec[3], vel[3];
                                    vec[X] = x1 - x2;
                                    vec[Y] = y1 - y2;
                                    vec[Z] = z1 - z2;
                                    cross(vec, burg1, p1);
                                    if (Normal(p1) > eps) {
                                        NormalizeVec(p1);
                                    } else {
                                        vel[X] = node1->vX * (1.0-L1) +
                                                 node2->vX * L1;
                                        vel[Y] = node1->vY * (1.0-L1) +
                                                 node2->vY * L1;
                                        vel[Z] = node1->vZ * (1.0-L1) +
                                                 node2->vZ * L1;
                                        cross(vec, vel, p1);
                                        if (Normal(p1) > eps) {
                                            NormalizeVec(p1);
                                        } else {
                                            VECTOR_ZERO(p1);
                                        }
                                    }

                                    vec[X] = x3 - x4;
                                    vec[Y] = y3 - y4;
                                    vec[Z] = z3 - z4;
                                    cross(vec, burg2, p2);
                                    if (Normal(p2) > eps) {
                                        NormalizeVec(p2);
                                    } else {
                                        vel[X] = node3->vX * (1.0-L2) +
                                                 node4->vX * L2;
                                        vel[Y] = node3->vY * (1.0-L2) +
                                                 node4->vY * L2;
                                        vel[Z] = node3->vZ * (1.0-L2) +
                                                 node4->vZ * L2;
                                        cross(vec, vel, p2);
                                        if (Normal(p2) > eps) {
                                            NormalizeVec(p2);
                                        } else {
                                            VECTOR_ZERO(p2);
                                        }
                                    }
                                }

                                FindIdealCollisionPoint(home, mPoint1, mPoint2,
                                                        p1, p2, cPoint);
/*
 *                              Redefine L1 based on the 'ideal' collision
 *                              point and determine if we can use one of the
 *                              segment endpoints or if we have to bisect
 *                              segment 1 and add a new node
 *
 *                              Note: this calculation will be based solely
 *                              on nodal position; hence all velocity
 *                              components are set to zero.
 */
                                GetMinDist(x1, y1, z1, 0.0, 0.0, 0.0,
                                           x2, y2, z2, 0.0, 0.0, 0.0,
                                           cPoint[X], cPoint[Y], cPoint[Z],
                                           0.0, 0.0, 0.0,
                                           cPoint[X], cPoint[Y], cPoint[Z],
                                           0.0, 0.0, 0.0,
                                           &dist2, &ddist2dt, &L1, &L2);

/*
 *                              Identify the first node to be merged.  If the
 *                              collision point is close to one of the nodal
 *                              endpoints, use that node, otherwise insert a
 *                              new node in the segment.
 *
 *                              NOTE: The current domain owns node1 but may
 *                              not own node2.  If it does not own node2, we 
 *                              cannot allow the collision to use node2
 *                              even if the collision point is close to 
 *                              that node.
 */
                                seg1Lx = x1 - x2;
                                seg1Ly = y1 - y2;
                                seg1Lz = z1 - z2;
    
                                close2node1 = (((seg1Lx*seg1Lx +
                                                 seg1Ly*seg1Ly +
                                                 seg1Lz*seg1Lz) *
                                                (L1 * L1)) < mindist2);
    
                                close2node2 = (((seg1Lx*seg1Lx +
                                                 seg1Ly*seg1Ly +
                                                 seg1Lz*seg1Lz) *
                                                ((1-L1) * (1-L1))) < mindist2);

                                if ((node2->myTag.domainID != thisDomain) &&
                                    close2node2) {
                                    continue;
                                }

                                if (close2node1) {
                                    mergenode1 = node1;
                                    splitSeg1 = 0;
                                } else if (close2node2) {
                                    mergenode1 = node2;
                                    splitSeg1 = 0;
                                } else {
                                    splitSeg1 = 1;
                                }

/*
 *                              If we need to add a new node to the first
 *                              segment, do it now.
 */
                                if (splitSeg1) {

                                     newx = x1 * (1.0-L1) + x2*L1;
                                     newy = y1 * (1.0-L1) + y2*L1;
                                     newz = z1 * (1.0-L1) + z2*L1;
    
                                     newvx = vx1 * (1.0-L1) + vx2*L1;
                                     newvy = vy1 * (1.0-L1) + vy2*L1;
                                     newvz = vz1 * (1.0-L1) + vz2*L1;
    
/*
 *                                   Estimate resulting forces on all segments
 *                                   involved in the split.
 */
                                     arm21 = GetArmID(node2, node1);
    
                                     oldfp0[X] = node1->armfx[arm12];
                                     oldfp0[Y] = node1->armfy[arm12];
                                     oldfp0[Z] = node1->armfz[arm12];
    
                                     oldfp1[X] = node2->armfx[arm21];
                                     oldfp1[Y] = node2->armfy[arm21];
                                     oldfp1[Z] = node2->armfz[arm21];
    
                                     p0[X] = x1;   p1[X] = x2;
                                     p0[Y] = y1;   p1[Y] = y2;
                                     p0[Z] = z1;   p1[Z] = z2;
    
                                     pnew[X] = newx;
                                     pnew[Y] = newy;
                                     pnew[Z] = newz;
    
                                     newNodeVel[X] = newvx;
                                     newNodeVel[Y] = newvy;
                                     newNodeVel[Z] = newvz;

                                     nodeVel[X] = vx1;
                                     nodeVel[Y] = vy1;
                                     nodeVel[Z] = vz1;

                                     FindSubFSeg(home, p0, p1, burg1, oldfp0,
                                                 oldfp1, pnew, f0seg1, f1seg1,
                                                 f0seg2, f1seg2);
    
                                     oldTag1 = node1->myTag;
                                     oldTag2 = node2->myTag;

                                     FoldBox(param, &pnew[X], &pnew[Y], &pnew[Z]);

                                     splitStatus = SplitNode(home,
                                                             OPCLASS_COLLISION,
                                                             node1, p0, pnew,
                                                             nodeVel,
                                                             newNodeVel, 1,
                                                             &arm12, globalOp,
                                                             &splitNode1,
                                                             &splitNode2, 0);
/*
 *                                   If we were unable to split the node
 *                                   go back to looking for more collision
 *                                   candidates.
 */
                                     if (splitStatus == SPLIT_FAILED) {
                                         continue;
                                     }

/*
 *                                   The force estimates above are good enough
 *                                   for the remainder of this timestep, but
 *                                   mark the force and velocity data for
 *                                   some nodes as obsolete so that more
 *                                   accurate forces will be recalculated
 *                                   either at the end of this timestep, or
 *                                   the beginning of the next.
 */
                                     mergenode1 = splitNode2;

                                     MarkNodeForceObsolete(home, splitNode2);
    
                                     for (q = 0; q < splitNode2->numNbrs; q++) {
                                         tmpNbr = GetNodeFromTag(home, splitNode2->nbrTag[q]);
                                         if (tmpNbr != (Node_t *)NULL) {
                                             tmpNbr->flags |= NODE_RESET_FORCES;
                                         }
                                     }
    
/*
 *                                   Reset nodal forces on nodes involved in the
 *                                   split.
 */
                                     ResetSegForces(home, splitNode1,
                                                    &splitNode2->myTag,
                                                    f0seg1[X], f0seg1[Y],
                                                    f0seg1[Z], 1);
    
                                     ResetSegForces(home, splitNode2,
                                                    &splitNode1->myTag,
                                                    f1seg1[X], f1seg1[Y],
                                                    f1seg1[Z], 1);
    
                                     ResetSegForces(home, splitNode2,
                                                    &node2->myTag,
                                                    f0seg2[X], f0seg2[Y],
                                                    f0seg2[Z], 1);
    
                                     ResetSegForces(home, node2,
                                                    &splitNode2->myTag,
                                                    f1seg2[X], f1seg2[Y],
                                                    f1seg2[Z], 1);
    
                                     (void)EvaluateMobility(home, splitNode1,
                                                            &mobArgs);

                                     (void)EvaluateMobility(home, splitNode2,
                                                            &mobArgs);

                                     (void)EvaluateMobility(home, node2,
                                                            &mobArgs);
    
/*
 *                                  When debugging, dump some info on
 *                                  topological changes taking place and
 *                                  the nodes involved
 */
#ifdef DEBUG_TOPOLOGY_CHANGES
                                    if ((dbgDom < 0)||(dbgDom == home->myDomain)) {
                                        printf("  Split-1(SegCollision): "
                                               "(%d,%d)--(%d,%d) ==> "
                                               "(%d,%d)--(%d,%d)--(%d,%d)\n",
                                               oldTag1.domainID, oldTag1.index,
                                               oldTag2.domainID, oldTag2.index,
                                               splitNode1->myTag.domainID,
                                               splitNode1->myTag.index,
                                               splitNode2->myTag.domainID,
                                               splitNode2->myTag.index,
                                               node2->myTag.domainID,
                                               node2->myTag.index);
                                        PrintNode(splitNode1);
                                        PrintNode(splitNode2);
                                        PrintNode(node2);
                                     }
#endif
/*
 *                                   When we actually do a collision, we flag
 *                                   it so we don't attempt additional changes
 *                                   on either segment.  It's possible, however
 *                                   that we will split the segment but later
 *                                   determine a collision will not be done.
 *                                   Since the original node1/node2 segment
 *                                   has been bisected, we must treat it as
 *                                   if a collision had occurred so we do not
 *                                   attempt to collide the now non-existent
 *                                   node1/node2 segment.
 */
                                     didCollision = 1;

                                }  /* if (splitSeg1) */

/*
 *                              Redefine L2 based on the 'ideal' collision
 *                              point and determine if we can use one of the
 *                              segment endpoints or if we have to bisect
 *                              segment 2 and add a new node
 *
 *                              Note: this calculation will be based solely
 *                              on nodal position; hence all velocity
 *                              components are set to zero.
 */
                                GetMinDist(cPoint[X], cPoint[Y], cPoint[Z],
                                           0.0, 0.0, 0.0,
                                           cPoint[X], cPoint[Y], cPoint[Z],
                                           0.0, 0.0, 0.0,
                                           x3, y3, z3, 0.0, 0.0, 0.0,
                                           x4, y4, z4, 0.0, 0.0, 0.0,
                                           &dist2, &ddist2dt, &L1, &L2);

/*
 *                              Identify the second node to be merged
 *
 *                              Note: The current domain owns node3 but may not
 *                              own node4.  If it does not own node4, we 
 *                              cannot allow the collision to use node4
 *                              even if the collision point is close to 
 *                              that node.
 */
                                seg2Lx = x3 - x4;
                                seg2Ly = y3 - y4;
                                seg2Lz = z3 - z4;
    
                                close2node3 = (((seg2Lx*seg2Lx +
                                                 seg2Ly*seg2Ly +
                                                 seg2Lz*seg2Lz) *
                                                (L2 * L2)) < mindist2);
    
                                close2node4 = (((seg2Lx*seg2Lx +
                                                 seg2Ly*seg2Ly +
                                                 seg2Lz*seg2Lz) *
                                                ((1-L2) * (1-L2))) < mindist2);
    
                                if ((node4->myTag.domainID != thisDomain) &&
                                    close2node4) {
                                    continue;
                                }

                                if (close2node3) {
                                    mergenode2 = node3;
                                    splitSeg2 = 0;
                                } else if (close2node4) {
                                    mergenode2 = node4;
                                    splitSeg2 = 0;
                                } else {
                                    splitSeg2 = 1;
                                }

/*
 *                              If we need to add a new node to the second
 *                              segment, do it now.
 */
                                if (splitSeg2) {

                                     newx = x3 * (1.0-L2) + x4*L2;
                                     newy = y3 * (1.0-L2) + y4*L2;
                                     newz = z3 * (1.0-L2) + z4*L2;
    
                                     newvx = vx3 * (1.0-L2) + vx4*L2;
                                     newvy = vy3 * (1.0-L2) + vy4*L2;
                                     newvz = vz3 * (1.0-L2) + vz4*L2;
    
/*
 *                                   Estimate resulting forces on all segments
 *                                   involved in the split.
 */
                                     arm43 = GetArmID(node4, node3);
    
                                     burg1[X] = node3->burgX[arm34];
                                     burg1[Y] = node3->burgY[arm34];
                                     burg1[Z] = node3->burgZ[arm34];
    
                                     oldfp0[X] = node3->armfx[arm34];
                                     oldfp0[Y] = node3->armfy[arm34];
                                     oldfp0[Z] = node3->armfz[arm34];
    
                                     oldfp1[X] = node4->armfx[arm43];
                                     oldfp1[Y] = node4->armfy[arm43];
                                     oldfp1[Z] = node4->armfz[arm43];
    
                                     p0[X] = x3;   p1[X] = x4;
                                     p0[Y] = y3;   p1[Y] = y4;
                                     p0[Z] = z3;   p1[Z] = z4;
    
                                     pnew[X] = newx;
                                     pnew[Y] = newy;
                                     pnew[Z] = newz;
    
                                     newNodeVel[X] = newvx;
                                     newNodeVel[Y] = newvy;
                                     newNodeVel[Z] = newvz;

                                     nodeVel[X] = vx3;
                                     nodeVel[Y] = vy3;
                                     nodeVel[Z] = vz3;

                                     FindSubFSeg(home, p0, p1, burg1, oldfp0,
                                                 oldfp1, pnew, f0seg1, f1seg1,
                                                 f0seg2, f1seg2);
    
                                     oldTag1 = node3->myTag;
                                     oldTag2 = node4->myTag;

                                     FoldBox(param, &pnew[X], &pnew[Y], &pnew[Z]);

                                     splitStatus = SplitNode(home,
                                                            OPCLASS_COLLISION,
                                                            node3, p0,
                                                            pnew, nodeVel,
                                                            newNodeVel, 1,
                                                            &arm34, globalOp,
                                                            &splitNode1,
                                                            &splitNode2, 0);
/*
 *                                   If we were unable to split the node
 *                                   go back to looking for more collision
 *                                   candidates.
 */
                                     if (splitStatus == SPLIT_FAILED) {
                                         continue;
                                     }
/*
 *                                   The force estimates above are good enough
 *                                   for the remainder of this timestep, but
 *                                   mark the force and velocity data for some
 *                                   nodes as obsolete so that more accurate
 *                                   forces will be recalculated either at the
 *                                   end of this timestep, or the beginning of
 *                                   the next.
 */
                                     mergenode2 = splitNode2;

                                     MarkNodeForceObsolete(home, splitNode2);
    
                                     for (q = 0; q < splitNode2->numNbrs; q++) {
                                         tmpNbr = GetNodeFromTag(home, splitNode2->nbrTag[q]);
                                         if (tmpNbr != (Node_t *)NULL) {
                                             tmpNbr->flags |= NODE_RESET_FORCES;
                                         }
                                     }
    
/*
 *                                   Reset nodal forces on nodes involved in the
 *                                   split.
 */
                                     ResetSegForces(home, splitNode1,
                                                    &splitNode2->myTag,
                                                    f0seg1[X], f0seg1[Y],
                                                    f0seg1[Z], 1);
    
                                     ResetSegForces(home, splitNode2,
                                                    &splitNode1->myTag,
                                                    f1seg1[X], f1seg1[Y],
                                                    f1seg1[Z], 1);
    
                                     ResetSegForces(home, splitNode2,
                                                    &node4->myTag,
                                                    f0seg2[X], f0seg2[Y],
                                                    f0seg2[Z], 1);
    
                                     ResetSegForces(home, node4,
                                                    &splitNode2->myTag,
                                                    f1seg2[X], f1seg2[Y],
                                                    f1seg2[Z], 1);
 
                                     (void)EvaluateMobility(home, splitNode1,
                                                            &mobArgs);

                                     (void)EvaluateMobility(home, splitNode2,
                                                            &mobArgs);

                                     (void)EvaluateMobility(home, node4,
                                                            &mobArgs);

/*
 *                                   When debugging, dump some info on
 *                                   topological changes taking place and
 *                                   the nodes involved
 */
#ifdef DEBUG_TOPOLOGY_CHANGES
                                     if ((dbgDom < 0) ||
                                         (dbgDom == home->myDomain)) {
                                         printf("  Split-2(SegCollision): "
                                               "(%d,%d)--(%d,%d) ==> "
                                               "(%d,%d)--(%d,%d)--(%d,%d)\n",
                                               oldTag1.domainID,
                                               oldTag1.index,
                                               oldTag2.domainID,
                                               oldTag2.index,
                                               splitNode1->myTag.domainID,
                                               splitNode1->myTag.index,
                                               splitNode2->myTag.domainID,
                                               splitNode2->myTag.index,
                                               node4->myTag.domainID,
                                               node4->myTag.index);
                                         PrintNode(splitNode1);
                                         PrintNode(splitNode2);
                                         PrintNode(node4);
                                    }
#endif
                                }  /* if (splitSeg2) */
    
    
/*
 *                             Select the physical position at which to place
 *                             the node resulting from the merge of the
 *                             two nodes.  However, if it looks like the
 *                             node will be thrown a significant distance
 *                             during the collision (due to glide con-
 *                             straints), don't do the collision.
 */
                               retCode = FindCollisionPoint(home, mergenode1,
                                                            mergenode2,
                                                            &newPos[X],
                                                            &newPos[Y],
                                                            &newPos[Z]);

/*
 *                             If we were unable to determine the collision
 *                             point, don't collide the segments.
 */
                               if (retCode != 1) {
                                   continue;
                               }

                               vec1[X] = newPos[X] - mergenode1->x;
                               vec1[Y] = newPos[Y] - mergenode1->y;
                               vec1[Z] = newPos[Z] - mergenode1->z;

                               vec2[X] = newPos[X] - mergenode2->x;
                               vec2[Y] = newPos[Y] - mergenode2->y;
                               vec2[Z] = newPos[Z] - mergenode2->z;

                               if ((DotProduct(vec1, vec1) > 16 * mindist2) &&
                                   (DotProduct(vec2, vec2) > 16 * mindist2)) {
                                   continue;
                               }

                               FoldBox(param, &newPos[X],&newPos[Y],&newPos[Z]);
    
#ifdef _ARLFEM
/*
 *                             If colliding 2 surface nodes, we may have to
 *                             adjust the collision point so it too is on the
 *                             surface.
 */
                               resetSurfaceProperties = 0;

                               if (HAS_ANY_OF_CONSTRAINTS(
                                       mergenode1->constraint, SURFACE_NODE) &&
                                   HAS_ANY_OF_CONSTRAINTS(
                                       mergenode2->constraint, SURFACE_NODE)) {

                                   Node_t *seg1Node2, *seg2Node2;

                                   seg1Node2 = (mergenode1 == node1) ?
                                               node2 : node1;
                                   seg2Node2 = (mergenode2 == node3) ?
                                               node4 : node3;

                                   FEM_AdjustCollisionPoint(mergenode1,
                                                            seg1Node2,
                                                            mergenode2,
                                                            seg2Node2,
                                                            newPos,
                                                            femSurfaceNorm);
                                   resetSurfaceProperties = 1;
                               }
#endif  /* ifdef _ARLFEM */

                               newVel[X] = half * (mergenode1->vX +
                                                   mergenode2->vX);
                               newVel[Y] = half * (mergenode1->vY +
                                                   mergenode2->vY);
                               newVel[Z] = half * (mergenode1->vZ +
                                                   mergenode2->vZ);
    
    
                               oldTag1 = mergenode1->myTag;
                               oldTag2 = mergenode2->myTag;

                               MergeNode(home, OPCLASS_COLLISION, mergenode1,
                                         mergenode2, newPos, &targetNode,
                                         &mergeStatus, globalOp);
    
/*
 *                             If the merge did not succeed, go back and
 *                             continue looking for collision candidates.
 */
                               if ((mergeStatus & MERGE_SUCCESS) == 0) {
                                   continue;
                               }
#ifdef _ARLFEM
/*
 *                             Need to explicitly reset surface properties
 *                             after colliding 2 surface nodes.
 */
                               if (resetSurfaceProperties) {
                                   if (targetNode != (Node_t *)NULL) {
                                       VECTOR_COPY(targetNode->surfaceNorm,
                                                   femSurfaceNorm);
                                   }
                               }
#endif

/*
 *                             When debugging, dump some info on topological
 *                             changes taking place and the nodes involved
 */
#ifdef DEBUG_TOPOLOGY_CHANGES
                               if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                                   if (targetNode == (Node_t *)NULL) {
                                       printf("  Merge(SegCollision): "
                                              "(%d,%d) and (%d,%d)\n",
                                              oldTag1.domainID, oldTag1.index,
                                              oldTag2.domainID, oldTag2.index);
                                   } else {
                                       printf("  Merge(SegCollision): "
                                              "(%d,%d) and (%d,%d) at "
                                              "(%d,%d)\n",
                                              oldTag1.domainID, oldTag1.index,
                                              oldTag2.domainID, oldTag2.index,
                                              targetNode->myTag.domainID,
                                              targetNode->myTag.index);
                                       PrintNode(targetNode);
                                   }
                               }
#endif

#if 0
/*
 *                             Under rare circumstances a collision occurs
 *                             in which the selected collision point is 
 *                             at the exact same location as another node
 *                             which is already connected to targetNode.
 *                             (Typically could happen when glide planes
 *                             are being enforced)  If this does happen
 *                             it leaves a zero length segment between
 *                             two nodes which wreaks havoc with the code.
 *                             This block of code tries to identify such
 *                             situations and handle it by merging the
 *                             two nodes that are at the same position.
 *
 *                             The drawback is that it may be possible
 *                             to encounter such a situation in which
 *                             the two nodes cannot be merged due to
 *                             constraints imposed by domain boundaries
 *                             and node ownership rules.  If that happens
 *                             we're still screwed.
 */
                               if (targetNode != (Node_t *)NULL) {
                                   int a, zeroLenSeg = 0;
                                   real8 lenX, lenY, lenZ, len2;
                                   Node_t *nbr;
                                   for (a = 0; a < targetNode->numNbrs; a++) {
                                       nbr = GetNeighborNode(home,targetNode,a);
                                       if (nbr == (Node_t *)NULL) continue;
                                       lenX = nbr->x - targetNode->x;
                                       lenY = nbr->y - targetNode->y;
                                       lenZ = nbr->z - targetNode->z;
                                       ZImage(param, &lenX, &lenY, &lenZ);
                                       len2 = lenX*lenX + lenY*lenY + lenZ*lenZ;
/*
 *                                     Treat a segment length less than 1e-10
 *                                     as zero length 
 */
                                       if (len2 < 1.0e-20) {
                                           mergenode1 = targetNode;
                                           mergenode2 = nbr;
                                           newPos[X] = targetNode->x;
                                           newPos[Y] = targetNode->y;
                                           newPos[Z] = targetNode->z;
                                           zeroLenSeg = 1;
                                           break;
                                       }
                                   }

/*
 *                                 Found a zero-length connection between
 *                                 the target node and one of its neighbors
 *                                 so try to remove it
 */
                                   if (zeroLenSeg) {
                                       MergeNode(home, OPCLASS_COLLISION,
                                                 mergenode1, mergenode2,
                                                 newPos, &targetNode,
                                                 &mergeStatus, globalOp);
                                   }
                               }
#endif

/*
 *                             If the target node exists after the merge,
 *                             reset its velocity, and reset the topological
 *                             change exemptions for the target node; use
 *                             NodeTopologyExemptions() to get the basic
 *                             exemptions, then exempt all the node's arms
 *                             from additional segment collisions this cycle.
 */
                               if (targetNode != (Node_t *)NULL) {
/*
 *                                 If we are enforcing glide planes but
 *                                 allowing some fuzziness in the planes, we
 *                                 also need to recalculate the glide
 *                                 planes for the segments attched to the
 *                                 collision node.
 */
                                   if (param->enforceGlidePlanes &&
                                       param->allowFuzzyGlidePlanes) {
                                       int n;
                                       for (n=0; n<targetNode->numNbrs; n++) {
                                           tmpNbr = GetNodeFromTag(home,
                                                   targetNode->nbrTag[n]);
                                           RecalcSegGlidePlane(home,
                                                               targetNode,
                                                               tmpNbr, 1);
                                       }
                                   }
/*
 *                                 Estimate velocity so mobility function
 *                                 has a reasonable starting point
 */
                                   targetNode->vX = newVel[X];
                                   targetNode->vY = newVel[Y];
                                   targetNode->vZ = newVel[Z];
    
                                   (void)EvaluateMobility(home, targetNode,
                                                          &mobArgs);

                                   targetNode->flags |= NO_COLLISIONS;
#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
                                   targetNode->multiNodeLife = 0;
#endif
                               }
    
                               didCollision = 1;
                               localCollisionCnt++;

                            }  /* If conditions for collision are met */
                        }  /* Loop over node3 arms */
                    }  /* Loop over node1 arms */
                }  /* while(nextIndex >= 0) */
            }  /* Loop over neighboring cell2s */
           }
          }
        }  /* for (i = 0 ...) */

/*
 *      Now we have to loop for collisions on hinge joints (i.e zipping)
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
            if (node1->flags & NO_COLLISIONS) continue;

            for (j = 0; j < node1->numNbrs; j++) {

                if (node1->myTag.domainID != node1->nbrTag[j].domainID)
                    continue;

                for (k = j + 1; k < node1->numNbrs; k++) {
                    real8 dx, dy, dz;

                    if (node1->myTag.domainID != node1->nbrTag[k].domainID)
                        continue;

                    node3 = GetNodeFromTag(home, node1->nbrTag[j]);
                    node4 = GetNodeFromTag(home, node1->nbrTag[k]);

                    x1 = node1->x; y1 = node1->y; z1 = node1->z;
                    x3 = node3->x; y3 = node3->y; z3 = node3->z;
                    x4 = node4->x; y4 = node4->y; z4 = node4->z;

                    vx1 = node1->vX; vy1 = node1->vY; vz1 = node1->vZ;
                    vx3 = node3->vX; vy3 = node3->vY; vz3 = node3->vZ;
                    vx4 = node4->vX; vy4 = node4->vY; vz4 = node4->vZ;

                    PBCPOSITION(param, x1, y1, z1, &x3, &y3, &z3);
                    PBCPOSITION(param, x3, y3, z3, &x4, &y4, &z4);

                    x2  = x1;  y2  = y1;  z2  = z1;
                    vx2 = vx1; vy2 = vy1; vz2 = vz1;

/*
 *                  It is possible to have a zero-length segment
 *                  (created by a previous collision).  If we find
 *                  such a segment, do not try to use it in any
 *                  subsequent collisions.
 */
                    dx = x1 - x3;
                    dy = y1 - y3;
                    dz = z1 - z3;

                    if ((dx*dx + dy*dy + dz*dz) < 1.0e-20) {
                        continue;
                    }

                    dx = x1 - x4;
                    dy = y1 - y4;
                    dz = z1 - z4;

                    if ((dx*dx + dy*dy + dz*dz) < 1.0e-20) {
                        continue;
                    }

/*
 *                  Find the minimum distance between the node1/node4
 *                  segment and the point at node3 to determine if they
 *                  should be collided.
 */
                    GetMinDist(x1, y1, z1, vx1, vy1, vz1,
                               x4, y4, z4, vx4, vy4, vz4,
                               x3, y3, z3, vx3, vy3, vz3,
                               x3, y3, z3, vx3, vy3, vz3,
                               &dist2, &ddist2dt, &L1, &L2);

                    if ((dist2 < mindist2) && (ddist2dt < -eps)) {
                        real8 mPoint1[3], mPoint2[3], cPoint[3];
                        real8 vec[3], vel[3], p2[3];

/*
 *                      A collision should occur, so calculate the ideal
 *                      collision point between the node1/node4 segment and
 *                      point node3, using the 'closest' point identified above
 */
                        mPoint1[X] = x1 * (1.0-L1) + x4 * L1;
                        mPoint1[Y] = y1 * (1.0-L1) + y4 * L1;
                        mPoint1[Z] = z1 * (1.0-L1) + z4 * L1;

                        mPoint2[X] = x3;
                        mPoint2[Y] = y3;
                        mPoint2[Z] = z3;

                        if (param->enforceGlidePlanes) {
                            p1[X] = node1->nx[k];
                            p1[Y] = node1->ny[k];
                            p1[Z] = node1->nz[k];
                            p2[X] = node1->nx[j];
                            p2[Y] = node1->ny[j];
                            p2[Z] = node1->nz[j];
                        } else {
                            vec[X] = x1 - x4;
                            vec[Y] = y1 - y4;
                            vec[Z] = z1 - z4;
                            burg1[X] = node1->burgX[k];
                            burg1[Y] = node1->burgY[k];
                            burg1[Z] = node1->burgZ[k];
                            cross(vec, burg1, p1);
                            if (Normal(p1) > eps) {
                                NormalizeVec(p1);
                            } else {
                                vel[X] = node1->vX * (1.0-L1) + node4->vX * L1;
                                vel[Y] = node1->vY * (1.0-L1) + node4->vY * L1;
                                vel[Z] = node1->vZ * (1.0-L1) + node4->vZ * L1;
                                cross(vec, vel, p1);
                                if (Normal(p1) > eps) {
                                    NormalizeVec(p1);
                                } else {
                                    VECTOR_ZERO(p1);
                                }
                            }

                            vec[X] = x1 - x3;
                            vec[Y] = y1 - y3;
                            vec[Z] = z1 - z3;
                            burg2[X] = node1->burgX[j];
                            burg2[Y] = node1->burgY[j];
                            burg2[Z] = node1->burgZ[j];
                            cross(vec, burg2, p2);
                            if (Normal(p2) > eps) {
                                NormalizeVec(p2);
                            } else {
                                vel[X] = node3->vX;
                                vel[Y] = node3->vY;
                                vel[Z] = node3->vZ;
                                cross(vec, vel, p2);
                                if (Normal(p2) > eps) {
                                    NormalizeVec(p2);
                                } else {
                                    VECTOR_ZERO(p2);
                                }
                            }
                        }

                        FindIdealCollisionPoint(home, mPoint1, mPoint2,
                                                p1, p2, cPoint);

/*
 *                      Redefine L1 (the relative distance along the
 *                      node1/node4 segment at which the collision
 *                      should occur) based on the 'ideal' collision point
 *                      and determine if we can use one of the node1/node4
 *                      segment endpoints or if we have to bisect the 
 *                      segment and add a new node
 *
 *                      Note: this calculation will be based solely
 *                      on nodal position; hence all velocity
 *                      components are set to zero.
 */
                        GetMinDist(x1, y1, z1, 0.0, 0.0, 0.0,
                                   x4, y4, z4, 0.0, 0.0, 0.0,
                                   cPoint[X], cPoint[Y], cPoint[Z],
                                   0.0, 0.0, 0.0,
                                   cPoint[X], cPoint[Y], cPoint[Z],
                                   0.0, 0.0, 0.0,
                                   &dist2, &ddist2dt, &L1, &L2);

                        mergenode1 = node3;


                        seg1Lx = x1 - x4;
                        seg1Ly = y1 - y4;
                        seg1Lz = z1 - z4;
                        
/*
 *                      Determine if the collision point is close to one of
 *                      the existing nodes on the node1/node4 segment
 */
                        close2node1 = (((seg1Lx*seg1Lx +
                                         seg1Ly*seg1Ly +
                                         seg1Lz*seg1Lz) *
                                        (L1 * L1)) < mindist2);

                        close2node4 = (((seg1Lx*seg1Lx +
                                         seg1Ly*seg1Ly +
                                         seg1Lz*seg1Lz) *
                                        ((1.0-L1) * (1.0-L1))) < mindist2);

/*
 *                      If the collision point is close to one of the
 *                      endpoints of the node1/node4 segment, but the
 *                      node in question is not owned by the current
 *                      domain, skip the collision this time around.
 */
                        if ((close2node1&&(node1->myTag.domainID!=thisDomain))||
                            (close2node4&&(node4->myTag.domainID!=thisDomain))){
                            continue;
                        }

                        if (close2node1) {
                             mergenode2 = node1;
                        } else if (close2node4) {
                             mergenode2 = node4;
                        } else { 
/*
 *                           Collision point is not close enough to one of
 *                           the segment endpoints, so bisect the segment
 *                           and add a new node.
 */
                             newPos[X] = x1*(1.0-L1) + x4*L1;
                             newPos[Y] = y1*(1.0-L1) + y4*L1;
                             newPos[Z] = z1*(1.0-L1) + z4*L1;

                             newVel[X] = vx1*(1.0-L1) + vx4*L1;
                             newVel[Y] = vy1*(1.0-L1) + vy4*L1;
                             newVel[Z] = vz1*(1.0-L1) + vz4*L1;

                             nodeVel[X] = node3->vX;
                             nodeVel[Y] = node3->vY;
                             nodeVel[Z] = node3->vZ;

/*
 *                           Estimate resulting forces on all segments
 *                           involved in the split.
 */
                             armAB = GetArmID(node1, node4);
                             armBA = GetArmID(node4, node1);

                             burg1[X] = node1->burgX[armAB];
                             burg1[Y] = node1->burgY[armAB];
                             burg1[Z] = node1->burgZ[armAB];

                             oldfp0[X] = node1->armfx[armAB];
                             oldfp0[Y] = node1->armfy[armAB];
                             oldfp0[Z] = node1->armfz[armAB];

                             oldfp1[X] = node4->armfx[armBA];
                             oldfp1[Y] = node4->armfy[armBA];
                             oldfp1[Z] = node4->armfz[armBA];

                             p0[X] = x1;   p1[X] = x4;
                             p0[Y] = y1;   p1[Y] = y4;
                             p0[Z] = z1;   p1[Z] = z4;

                             FindSubFSeg(home, p0, p1, burg1, oldfp0,
                                         oldfp1, newPos, f0seg1, f1seg1,
                                         f0seg2, f1seg2);

                             FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);

/*
 *                           Attempt to split the segment.  If the split
 *                           fails, we can't continue with this collision.
 */
                             oldTag1 = node1->myTag;
                             oldTag2 = node4->myTag;

                             splitStatus = SplitNode(home, OPCLASS_COLLISION,
                                                     node1, p0, newPos,
                                                     nodeVel, newVel, 1,
                                                     &k, globalOp,
                                                     &splitNode1,
                                                     &splitNode2, 0);

                             if (splitStatus != SPLIT_SUCCESS) {
                                 continue;
                             }

/*
 *                           The force estimates above are good enough for
 *                           the remainder of this timestep, but mark the
 *                           force and velocity data for some nodes as obsolete
 *                           so that more accurate forces will be recalculated
 *                           either at the end of this timestep, or the
 *                           beginning of the next.
 */
                             mergenode2 = splitNode2;

                             MarkNodeForceObsolete(home, splitNode1);
                             MarkNodeForceObsolete(home, splitNode2);
                             MarkNodeForceObsolete(home, node4);

/*
 *                           Reset nodal forces for nodes involved in the
 *                           split.
 */
                             ResetSegForces(home,splitNode1,&splitNode2->myTag,
                                            f0seg1[X], f0seg1[Y], f0seg1[Z], 1);

                             ResetSegForces(home,splitNode2,&splitNode1->myTag,
                                            f1seg1[X], f1seg1[Y], f1seg1[Z], 1);

                             ResetSegForces(home, splitNode2, &node4->myTag,
                                            f0seg2[X], f0seg2[Y], f0seg2[Z], 1);

                             ResetSegForces(home, node4, &splitNode2->myTag,
                                            f1seg2[X], f1seg2[Y], f1seg2[Z], 1);

                             (void)EvaluateMobility(home, splitNode1,
                                                    &mobArgs);

                             (void)EvaluateMobility(home, splitNode2,
                                                    &mobArgs);

                             (void)EvaluateMobility(home, node4,
                                                    &mobArgs);

/*
 *                           When debugging, dump some info on topological
 *                           changes taking place and the nodes involved
 */
#ifdef DEBUG_TOPOLOGY_CHANGES
                             if ((dbgDom < 0)||(dbgDom == home->myDomain)) {
                                 printf("  Split-1(Hinge): "
                                        "(%d,%d)--(%d,%d) ==> "
                                        "(%d,%d)--(%d,%d)--(%d,%d)\n",
                                        oldTag1.domainID, oldTag1.index,
                                        oldTag2.domainID, oldTag2.index,
                                        splitNode1->myTag.domainID,
                                        splitNode1->myTag.index,
                                        splitNode2->myTag.domainID,
                                        splitNode2->myTag.index,
                                        node4->myTag.domainID,
                                        node4->myTag.index);
                                 PrintNode(splitNode1);
                                 PrintNode(splitNode2);
                                 PrintNode(node4);
                             }
#endif
                        }

/*
 *                      Select the physical position at which to place
 *                      the node resulting from the merge of the two nodes.
 *                      However, if it looks like the node will be thrown a
 *                      significant distance during the collision (due to
 *                      glide constraints), don't do the collision.
 */

                        retCode = FindCollisionPoint(home, mergenode1, mergenode2,
                                           &newPos[X], &newPos[Y], &newPos[Z]);

/*
 *                      If we were unable to determine the collision
 *                      point, don't collide the segments.
 */
                        if (retCode != 1) {
                            continue;
                        }

                        vec1[X] = newPos[X] - mergenode1->x;
                        vec1[Y] = newPos[Y] - mergenode1->y;
                        vec1[Z] = newPos[Z] - mergenode1->z;

                        vec2[X] = newPos[X] - mergenode2->x;
                        vec2[Y] = newPos[Y] - mergenode2->y;
                        vec2[Z] = newPos[Z] - mergenode2->z;

                        if ((DotProduct(vec1, vec1) > 16 * mindist2) &&
                            (DotProduct(vec2, vec2) > 16 * mindist2)) {
                            continue;
                        }

                        FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);

/*
 *                      The merge is going to happen, so there's a few
 *                      more nodes we'll need to mark for force/velocity
 *                      re-evaluation.
 */
                        MarkNodeForceObsolete(home, node1);
                        MarkNodeForceObsolete(home, node3);
                        MarkNodeForceObsolete(home, node4);

                        for (q = 0; q < node1->numNbrs; q++) {
                            tmpNbr = GetNodeFromTag(home, node1->nbrTag[q]);
                            if (tmpNbr == (Node_t *)NULL) continue;
                            tmpNbr->flags |= NODE_RESET_FORCES;
                        }

                        for (q = 0; q < node3->numNbrs; q++) {
                            tmpNbr = GetNodeFromTag(home, node3->nbrTag[q]);
                            if (tmpNbr == (Node_t *)NULL) continue;
                            tmpNbr->flags |= NODE_RESET_FORCES;
                        }

                        for (q = 0; q < node4->numNbrs; q++) {
                            tmpNbr = GetNodeFromTag(home, node4->nbrTag[q]);
                            if (tmpNbr == (Node_t *)NULL) continue;
                            tmpNbr->flags |= NODE_RESET_FORCES;
                        }

/*
 *                      Try the collision
 */
                        oldTag1 = mergenode1->myTag;
                        oldTag2 = mergenode2->myTag;

                        MergeNode(home, OPCLASS_COLLISION, mergenode1,
                                  mergenode2, newPos, &targetNode,
                                  &mergeStatus, globalOp);

                        localCollisionCnt++;

/*
 *                      If the target node exists after the merge, reevaulate
 *                      the node velocity, mark the forces as obsolete,
 *                      and reset the topological change exemptions for
 *                      the target node and exempt all the node's arms
 *                      from subsequent segment collisions this cycle.
 */
                        if (targetNode != (Node_t *)NULL) {

/*
 *                          If we are enforcing glide planes but
 *                          allowing some fuzziness in the planes, we
 *                          also need to recalculate the glide
 *                          planes for the segments attched to the
 *                          collision node.
 */
                            if (param->enforceGlidePlanes &&
                                param->allowFuzzyGlidePlanes) {
                                int n;
                                for (n = 0; n < targetNode->numNbrs; n++) {
                                    tmpNbr = GetNodeFromTag(home,
                                            targetNode->nbrTag[n]);
                                    RecalcSegGlidePlane(home, targetNode,
                                                        tmpNbr, 1);
                                }
                            }

                            (void)EvaluateMobility(home, targetNode, &mobArgs);

                            targetNode->flags |= NODE_RESET_FORCES;
                            targetNode->flags |= NO_COLLISIONS;

#ifdef DEBUG_TOPOLOGY_CHANGES
                            if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                                printf("  Merge(HingeCollision): "
                                       "(%d,%d) and (%d,%d) at (%d,%d)\n",
                                       oldTag1.domainID, oldTag1.index,
                                       oldTag2.domainID, oldTag2.index,
                                       targetNode->myTag.domainID,
                                       targetNode->myTag.index);
                                PrintNode(targetNode);
                            }
#endif
#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
                            targetNode->multiNodeLife = 0;
#endif
                        }

                    }  /* conditions for zip were met */
                }  /* for (k = 0...) */
            }  /* for (j = 0...) */
        }  /* for (i = 0...) */


#ifdef DEBUG_LOG_COLLISIONS
#ifdef PARALLEL
        int globalCollisionCnt = 0;
        MPI_Reduce(&localCollisionCnt, &globalCollisionCnt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        int globalCollisionCnt = localCollisionCnt;
#endif
        if (home->myDomain == 0) {
            printf("  Collision count = %d\n", globalCollisionCnt);
        }
#endif

        TimerStop(home, COLLISION_HANDLING);

        return;
}
