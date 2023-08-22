/*****************************************************************************
 *
 *      Module:         PredictiveCollision.c
 *      Description:    This module contains various functions used
 *                      for detecting various types of collisions
 *                      (segment/segment, node/node, zipping) and
 *                      dealing with those collisions.  These functions
 *                      are specific to the type 2 collision handling
 *                      which attempts to predict if there is a future
 *                      time at which segments will physicaly intersect,
 *                      and uses this data as collision criteria rather
 *                      than the simple proximity criteria of the old
 *                      mechanism
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
 *      Included functions:
 *
 *          FindCollisionPointAndTime()
 *          PredictiveCollisions()
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "Comm.h"
#include "Mobility.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif  /* ifdef _ARLFEM */

static int dbgDom;

/*---------------------------------------------------------------------------
 *
 *      Function:       SplitSegment
 *
 *-------------------------------------------------------------------------*/
static bool SplitSegment(Home_t *home, Node_t *node1, Node_t *node2,
			 int arm12, real8 L1, 
			 Node_t **mergenode1, int Opn)
{
        int     splitStatus;
        real8   newx, newy, newz;
        real8   newvx, newvy, newvz;
        real8   x1, y1, z1, vx1, vy1, vz1;
        real8   x2, y2, z2, vx2, vy2, vz2;

	int     arm21;
        real8   oldfp0[3], oldfp1[3];
	real8   pos0[3], pos1[3];
        real8   pnew[3], burg1[3];
	real8   nodeVel[3], newNodeVel[3];
        real8   f0seg1[3], f1seg1[3], f0seg2[3], f1seg2[3];
        Tag_t   oldTag1, oldTag2;
	Node_t  *tmpNbr;
        MobArgs_t mobArgs;	
	Param_t *param = home->param;
        Node_t  *splitNode1, *splitNode2;
		
/*
 *  Save node 1 and node 2 for later
 */
	vx1 = node1->vX; vy1 = node1->vY; vz1 = node1->vZ;
	vx2 = node2->vX; vy2 = node2->vY; vz2 = node2->vZ;
	
	x1 = node1->x; y1 = node1->y; z1 = node1->z;
	x2 = node2->x; y2 = node2->y; z2 = node2->z;
	PBCPOSITION(param, x1, y1, z1, &x2, &y2, &z2);
	
	newx = x1 * (1.0-L1) + x2*L1;
	newy = y1 * (1.0-L1) + y2*L1;
	newz = z1 * (1.0-L1) + z2*L1;
    
	newvx = vx1 * (1.0-L1) + vx2*L1;
	newvy = vy1 * (1.0-L1) + vy2*L1;
	newvz = vz1 * (1.0-L1) + vz2*L1;    
/*
 *      Estimate resulting forces on all segments
 *      involved in the split.
 */
	arm21 = GetArmID(node2, node1);
  
	oldfp0[X] = node1->armfx[arm12];
	oldfp0[Y] = node1->armfy[arm12];
	oldfp0[Z] = node1->armfz[arm12];
	
	oldfp1[X] = node2->armfx[arm21];
	oldfp1[Y] = node2->armfy[arm21];
	oldfp1[Z] = node2->armfz[arm21];
	
	pos0[X] = x1;   pos1[X] = x2;
	pos0[Y] = y1;   pos1[Y] = y2;
	pos0[Z] = z1;   pos1[Z] = z2;
	
	pnew[X] = newx;
	pnew[Y] = newy;
	pnew[Z] = newz;
	
	newNodeVel[X] = newvx;
	newNodeVel[Y] = newvy;
	newNodeVel[Z] = newvz;
	
	nodeVel[X] = vx1;
	nodeVel[Y] = vy1;
	nodeVel[Z] = vz1;
	
	burg1[X] = node1->burgX[arm12];
	burg1[Y] = node1->burgY[arm12];
	burg1[Z] = node1->burgZ[arm12];
	
	FindSubFSeg(home, pos0, pos1, burg1, oldfp0,
		    oldfp1, pnew, f0seg1, f1seg1,
		    f0seg2, f1seg2);
	
	FoldBox(param, &pnew[X], &pnew[Y], &pnew[Z]);
	
	oldTag1 = node1->myTag;
	oldTag2 = node2->myTag;
	
	splitStatus = SplitNode(home,
				OPCLASS_COLLISION,
				node1, pos0, pnew,
				nodeVel,
				newNodeVel, 1,
				&arm12, Opn,
				&splitNode1,
				&splitNode2, 0);
	*mergenode1 = splitNode2;
	MarkNodeForceObsolete(home, splitNode2);
	for (int q = 0; q < splitNode2->numNbrs; q++) {
	  tmpNbr = GetNodeFromTag(home, splitNode2->nbrTag[q]);
	  if (tmpNbr != (Node_t *)NULL) {
	    tmpNbr->flags |= NODE_RESET_FORCES;
	  }
	}				     

	if(splitStatus == SPLIT_FAILED)
	  return splitStatus;

/*
 * Reset nodal forces on nodes involved in the
 * split.
 */
	ResetSegForces(home, splitNode1,
		       &splitNode2->myTag,
		       f0seg1[X], f0seg1[Y],
		       f0seg1[Z], Opn);
	
	ResetSegForces(home, splitNode2,
		       &splitNode1->myTag,
		       f1seg1[X], f1seg1[Y],
		       f1seg1[Z], Opn);
	
	ResetSegForces(home, splitNode2,
		       &node2->myTag,
		       f0seg2[X], f0seg2[Y],
		       f0seg2[Z], Opn);
	
	ResetSegForces(home, node2,
		       &splitNode2->myTag,
		       f1seg2[X], f1seg2[Y],
		       f1seg2[Z], Opn);
	
	(void)EvaluateMobility(home, splitNode1,
			       &mobArgs);
	
	(void)EvaluateMobility(home, splitNode2,
			       &mobArgs);
	
	(void)EvaluateMobility(home, node2,
			       &mobArgs);
	
/*
 * When debugging, dump some info on
 * topological changes taking place and
 * the nodes involved
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
	return true;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       SplitHingeSegment
 *
 *-------------------------------------------------------------------------*/
static bool SplitHingeSegment(Home_t *home, Node_t *node1, Node_t *node4,
			      int armAB, real8 L1, 
			      Node_t **mergenode2, int Opn)
{
	int   splitStatus;  
	real8 x1,  y1, z1, x4, y4, z4;
	real8 vx1,vy1,vz1,vx4,vy4,vz4;
	
	
	
	real8 pos0[3], pos1[3];
	real8 newPos[3], newVel[3], nodeVel[3];
	int   armBA;
	real8 burg1[3], oldfp0[3], oldfp1[3];
	Tag_t oldTag1, oldTag2;
	
	real8 f0seg1[3], f1seg1[3], f0seg2[3], f1seg2[3];	
        MobArgs_t mobArgs;
	Param_t *param = home->param;
	Node_t  *splitNode1, *splitNode2;
	
	x1 = node1->x; y1 = node1->y; z1 = node1->z;
	x4 = node4->x; y4 = node4->y; z4 = node4->z;
	PBCPOSITION(param, x1, y1, z1, &x4, &y4, &z4);
	
	newPos[X] = x1*(1.0-L1) + x4*L1;
	newPos[Y] = y1*(1.0-L1) + y4*L1;
	newPos[Z] = z1*(1.0-L1) + z4*L1;
	
	newVel[X] = vx1*(1.0-L1) + vx4*L1;
	newVel[Y] = vy1*(1.0-L1) + vy4*L1;
	newVel[Z] = vz1*(1.0-L1) + vz4*L1;
	
	nodeVel[X] = node1->vX;
	nodeVel[Y] = node1->vY;
	nodeVel[Z] = node1->vZ;

/*
 *      Estimate resulting forces on all segments
 *      involved in the split.
 */
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
	
	pos0[X] = x1;   pos1[X] = x4;
	pos0[Y] = y1;   pos1[Y] = y4;
	pos0[Z] = z1;   pos1[Z] = z4;
			     
	FindSubFSeg(home, pos0, pos1, burg1, oldfp0,
		    oldfp1, newPos, f0seg1, f1seg1,
		    f0seg2, f1seg2);
	
	FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);

/*
 *      Attempt to split the segment.  If the split
 *      fails, we can't continue with this collision.
 */
	oldTag1 = node1->myTag;
	oldTag2 = node4->myTag;

	splitStatus = SplitNode(home, OPCLASS_COLLISION,
				node1, pos0, newPos,
				nodeVel, newVel, 1,
				&armAB, Opn,
				&splitNode1,
				&splitNode2, 0);
	if (splitStatus != SPLIT_SUCCESS) {
	  return false;
	}
/*
 *      The force estimates above are good enough for
 *      the remainder of this timestep, but mark the
 *      force and velocity data for some nodes as obsolete
 *      so that more accurate forces will be recalculated
 *      either at the end of this timestep, or the
 *      beginning of the next.
 */
			     
	*mergenode2 = splitNode2;

	MarkNodeForceObsolete(home, splitNode1);
	MarkNodeForceObsolete(home, splitNode2);
	MarkNodeForceObsolete(home, node4);			     

/*
 *      Reset nodal forces for nodes involved in the
 *      split.
 */
	ResetSegForces(home,splitNode1,&splitNode2->myTag,
		       f0seg1[X], f0seg1[Y], f0seg1[Z], Opn);

	ResetSegForces(home,splitNode2,&splitNode1->myTag,
		       f1seg1[X], f1seg1[Y], f1seg1[Z], Opn);
	
	ResetSegForces(home, splitNode2, &node4->myTag,
		       f0seg2[X], f0seg2[Y], f0seg2[Z], Opn);
	
	ResetSegForces(home, node4, &splitNode2->myTag,
		       f1seg2[X], f1seg2[Y], f1seg2[Z], Opn);
	
	(void)EvaluateMobility(home, splitNode1,
			       &mobArgs);
	
	(void)EvaluateMobility(home, splitNode2,
			       &mobArgs);
	
	(void)EvaluateMobility(home, node4,
			       &mobArgs);

/*
 *      When debugging, dump some info on topological
 *      changes taking place and the nodes involved
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
	return true;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       FindCollisionPointAndTime
 *      Description:    Given two segments and their velocities, determine
 *                      if the segments will physically intersect at any
 *                      given time in the future.  If so, return to the
 *                      caller the position at which the collision will 
 *                      occur, the time (in units of current deltaTT)
 *                      at which the collision would occur, etc.
 *                      
 *
 *      Arguments:
 *          p1,p2       Endpoint nodes of segment 1
 *          p3,p4       Endpoint nodes of segment 2
 *          v1,v2 v3,v4 Velocity of nodes p1 thru p4 respectively
 *          cPoint      Will be set to the coordinates of the collision
 *                      point if the two segments will collide.
 *          cTime       Set to the time (in units of current timestep
 *                      length (deltaTT)) at which the two segments will
 *                      collide, or negative if not colliding.
 *          minDist2    Min distance between segments, squared.
 *          L1, L2
 *
 *-------------------------------------------------------------------------*/
static void FindCollisionPointAndTime(Home_t *home, real8 p1[3], real8 p2[3],
                                  real8 p3[3], real8 p4[3], real8 v1[3],
                                  real8 v2[3], real8 v3[3], real8 v4[3],
                                  real8 cPoint[3], real8 *cTime,
                                  real8 *minDist2, real8 *L1, real8 *L2)
{
        int   i;
        real8 A, B, C, D, E, F, G, H, J;
        real8 dt, matDet, eps;
        real8 L1a, L2a, L3, L3a;
        real8 dist2, ddist2dt, minDist2Tol;
        real8 tmpMat22Inv[2][2];
        real8 tmpVec2a[2], tmpVec2b[2];
        real8 tmpVec3[3];
        real8 R0[3], seg1[3], seg2[3], vDiff[3];
        real8 uBar10[3], vBar10[3];
        real8 rhs[3];
        real8 mat[3][3], matInv[3][3];


        eps = 1.0e-12;
        dt = home->param->deltaTT;

        for (i = 0; i < 3; i++) {
            R0[i] = p1[i] - p3[i];
            seg1[i] = p2[i] - p1[i];
            seg2[i] = p4[i] - p3[i];
            vBar10[i] = (v1[i] + v2[i]) * 5.0 * dt;
            uBar10[i] = (v3[i] + v4[i]) * 5.0 * dt;
            vDiff[i] = vBar10[i] - uBar10[i];
        }

        A =   DotProduct(seg1, seg1);
        B = -(DotProduct(seg1, seg2));
        C =   DotProduct(seg1, vDiff);
        D =   DotProduct(seg2, seg2);
        E = -(DotProduct(seg2, vDiff));
        F =   DotProduct(vDiff, vDiff);
        G =   DotProduct(R0, seg1);
        H = -(DotProduct(R0, seg2));
        J =   DotProduct(R0, vDiff);

        mat[0][0] = A;  mat[0][1] = C;  mat[0][2] = B;
        mat[1][0] = C;  mat[1][1] = F;  mat[1][2] = E;
        mat[2][0] = B;  mat[2][1] = E;  mat[2][2] = D;

        rhs[0] = -G;
        rhs[1] = -J;
        rhs[2] = -H;

/*
 *      Check the solution of all the faces, edges and corners.
 */
        L3 = 0.0;  /* Start with zero time */

        GetMinDist(p1[0], p1[1], p1[2], v1[0], v1[1], v1[2],
                   p2[0], p2[1], p2[2], v2[0], v2[1], v2[2],
                   p3[0], p3[1], p3[2], v3[0], v3[1], v3[2],
                   p4[0], p4[1], p4[2], v4[0], v4[1], v4[2],
                   minDist2, &ddist2dt, L1, L2);

/*
 *      Check the end node pairs of the segments
 */
        if (F > eps) {

            L1a = 0.0;  /* Check for this node 1 on seg1 */
            L2a = 0.0;  /* Check for this node 1 on seg2 */
            L3a = J / -F;

            if (L3a > 0.0) {

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] + vDiff[i] * L3a;
                }

                dist2 = DotProduct(tmpVec3, tmpVec3);
              
                if (dist2 < *minDist2) {
                    *L1 = L1a;
                    *L2 = L2a;
                    L3 = L3a;
                    *minDist2 = dist2;
                }
            }

            L1a = 1.0;  /* Check for this node 2 on seg1 */
            L2a = 0.0;  /* Check for this node 1 on seg2 */
            L3a = (J + C) / -F;

            if (L3a > 0.0) {

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] + seg1[i] + vDiff[i] * L3a;
                }

                dist2 = DotProduct(tmpVec3, tmpVec3);

                if (dist2 < *minDist2) {
                    *L1 = L1a;
                    *L2 = L2a;
                    L3 = L3a;
                    *minDist2 = dist2;
                }
            }

            L1a = 0.0;  /* Check for this node 1 on seg1 */
            L2a = 1.0;  /* Check for this node 2 on seg2 */
            L3a = (J + E) / -F;

            if (L3a > 0.0) {

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] - seg2[i] + vDiff[i] * L3a;
                }

                dist2 = DotProduct(tmpVec3, tmpVec3);

                if (dist2 < *minDist2) {
                    *L1 = L1a;
                    *L2 = L2a;
                    L3 = L3a;
                    *minDist2 = dist2;
                }
            }

            L1a = 1.0;  /* Check for this node 2 on seg1 */
            L2a = 1.0;  /* Check for this node 2 on seg2 */
            L3a = (J + C + E) / -F;

            if (L3a > 0.0) {

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] + seg1[i] - seg2[i] + vDiff[i] * L3a;
                }

                dist2 = DotProduct(tmpVec3, tmpVec3);

                if (dist2 < *minDist2) {
                    *L1 = L1a;
                    *L2 = L2a;
                    L3 = L3a;
                    *minDist2 = dist2;
                }
            }

        } /* end if (F > eps) */

/*
 *      Check single end nodes on seg1
 */
        matDet = (mat[1][1]*mat[2][2])-(mat[2][1]*mat[1][2]);

        if (fabs(matDet)>eps) {
            L1a = 0.0;

            tmpMat22Inv[0][0] =  mat[2][2]/matDet;    // tmpMat22Inv = Adj(lower_right(mat))/Det(lower_right(mat))
            tmpMat22Inv[0][1] = -mat[1][2]/matDet;    //
            tmpMat22Inv[1][0] = -mat[2][1]/matDet;    //
            tmpMat22Inv[1][1] =  mat[1][1]/matDet;    //

            Matrix22_Vmul(tmpVec2a, tmpMat22Inv, &rhs[1] );

            if ((tmpVec2a[1] >= 0.0) &&
                (tmpVec2a[1] <= 1.0) &&
                (tmpVec2a[0] >= 0.0)) {

                L2a = tmpVec2a[1];
                L3a = tmpVec2a[0];

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] - (seg2[i] * L2a) + (vDiff[i] * L3a);
                }

                dist2 = DotProduct(tmpVec3, tmpVec3);

                if (dist2 < *minDist2) {
                    *L1 = L1a;
                    *L2 = L2a;
                    L3 = L3a;
                    *minDist2 = dist2;
                }
            }

            L1a = 1.0;

            tmpVec2a[0] = rhs[1] - mat[1][0];
            tmpVec2a[1] = rhs[2] - mat[2][0];
 
            /* note: tmpMat22Inv calculated above */
            Matrix22_Vmul(tmpVec2b, tmpMat22Inv, tmpVec2a );

            if ((tmpVec2b[1] >= 0.0) &&
                (tmpVec2b[1] <= 1.0) &&
                (tmpVec2b[0] >= 0.0)) {

                L2a = tmpVec2b[1];
                L3a = tmpVec2b[0];

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] - (seg2[i] * L2a) + (vDiff[i] * L3a);
                }

                dist2 = DotProduct(tmpVec3, tmpVec3);

                if (dist2 < *minDist2) {
                    *L1 = L1a;
                    *L2 = L2a;
                    L3 = L3a;
                    *minDist2 = dist2;
                }
            }

        } /* end if (fabs(D*F - E*E) > eps) */

/*
 *      Check single end nodes on seg2
 */
        matDet = (mat[0][0]*mat[1][1])-(mat[1][0]*mat[0][1]);

        if (fabs(matDet)>eps) {

            L2a = 0.0;

            tmpMat22Inv[0][0] =  mat[1][1]/matDet;    // tmpMat22Inv = Adj(upper_left(mat))/Det(upper_left(mat))
            tmpMat22Inv[0][1] = -mat[0][1]/matDet;    //
            tmpMat22Inv[1][0] = -mat[1][0]/matDet;    //
            tmpMat22Inv[1][1] =  mat[0][0]/matDet;    //

            Matrix22_Vmul(tmpVec2a, tmpMat22Inv, rhs );

            if ((tmpVec2a[0] >= 0.0) &&
                (tmpVec2a[0] <= 1.0) &&
                (tmpVec2a[1] >= 0.0)) {

                L1a = tmpVec2a[0];
                L3a = tmpVec2a[1];

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] + (seg1[i] * L1a) + (vDiff[i] * L3a);
                }

                dist2 = DotProduct(tmpVec3, tmpVec3);

                if (dist2 < *minDist2) {
                    *L1 = L1a;
                    *L2 = L2a;
                     L3 = L3a;
                    *minDist2 = dist2;
                }
            }
            
            L2a = 1.0;

            tmpVec2a[0] = rhs[0] - mat[0][2];
            tmpVec2a[1] = rhs[1] - mat[1][2];
 
            /* note: tmpMat22Inv calculated above */
            Matrix22_Vmul(tmpVec2b,tmpMat22Inv, tmpVec2a );

            if ((tmpVec2b[0] >= 0.0) &&
                (tmpVec2b[0] <= 1.0) &&
                (tmpVec2b[1] >= 0.0)) {

                L1a = tmpVec2b[0];
                L3a = tmpVec2b[1];

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] + (seg1[i] * L1a) - seg2[i] +
                                 (vDiff[i] * L3a);
                }

                dist2 = DotProduct(tmpVec3, tmpVec3);

                if (dist2 < *minDist2) {
                    *L1 = L1a;
                    *L2 = L2a;
                    L3 = L3a;
                    *minDist2 = dist2;
                }
            }

        }  /* end if (fabs(A*F - C*C) > eps) */

        matDet = Matrix33_Det(mat);

        if (fabs(matDet) > eps) {

            if ( Matrix33_Inverse(matInv, mat) < 0 ) 
            {
#if 0
                Fatal("Unable to invert matrix in FindCollisionPointAndTime");
#endif
                *L1 = 0.0;
                *L2 = 0.0;
                *cTime = -1;

                return;
            }

            Matrix33Vector3Multiply(matInv, rhs, tmpVec3);

            if ((tmpVec3[0] >= 0.0) && (tmpVec3[0] <= 1.0) &&
                (tmpVec3[2] >= 0.0) && (tmpVec3[2] <= 1.0) &&
                (tmpVec3[1] >= 0.0)) {
/*
 *              An interior point has been found
 */
                *L1 = tmpVec3[0];
                *L2 = tmpVec3[2];
                L3 = tmpVec3[1];

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] + (seg1[i] * *L1) - (seg2[i] * *L2) +
                                 (vDiff[i] * L3);
                }

                *minDist2 = DotProduct(tmpVec3, tmpVec3);
            }

        } /* end if (fabs(matDet) > eps) */


/*
 *      In general, if the dislocation segments on their current
 *      trajectories will not eventually intersect we don't want
 *      to do a collision and set the collision time negative to
 *      indicate no collision.  First, determine the tolerance for
 *      deciding if the segments will be within collision distance in
 *      the future.  Use a tight tolerance if strict glide planes
 *      are enforced, loosen up the tolerance a little in all
 *      other cases.
 */
        if (home->param->enforceGlidePlanes &&
            (home->param->allowFuzzyGlidePlanes == 0)) {
            minDist2Tol = eps;
        } else {
            minDist2Tol = 1.0;
        }


/*
 *      If collision criteria is not met, set collision time to a
 *      negative value.  Otherwise calculate the collision point and
 *      time.
 */
        if (*minDist2 > minDist2Tol) {

            *cTime = -1.0;

        } else {

            *cTime = 10.0 * L3;

            for (i = 0; i < 3; i++) {
                cPoint[i] = 0.5 * ((p1[i] + p3[i]) +
                                   (*L1 * seg1[i]) +
                                   (*L2 * seg2[i]) +
                                   (vBar10[i] + uBar10[i]) * L3);
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       PredictiveCollisions
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
void PredictiveCollisions(Home_t *home)
{
        int     i, j, k, q, arm12, arm21, arm34, arm43;
        int     thisDomain, splitStatus, mergeStatus;
        int     armAB, armBA;
        int     globalOp = 1, didCollision;
        int     close2node1, close2node2, close2node3, close2node4;
        int     splitSeg1, splitSeg2;
        int     cell2Index, nbrCell2Index, nextIndex;
        int     cell2X, cell2Y, cell2Z, cx, cy, cz;
        int     collisionConditionIsMet, adjustCollisionPoint=0;
        real8   mindist2, dist2, ddist2dt, L1, L2, half;
        real8   cTime, cPoint[3];
        real8   x1, y1, z1, vx1, vy1, vz1;
        real8   x2, y2, z2, vx2, vy2, vz2;
        real8   x3, y3, z3, vx3, vy3, vz3;
        real8   x4, y4, z4, vx4, vy4, vz4;
        real8   newx, newy, newz;
        real8   newvx, newvy, newvz;
        real8   pnew[3], burg1[3];
        real8   oldfp0[3], oldfp1[3];
        real8   f0seg1[3], f1seg1[3], f0seg2[3], f1seg2[3];
        real8   nodeVel[3], newNodeVel[3];
        real8   newPos[3], newVel[3];
        real8   vec1[3], vec2[3];
        real8   p1[3], p2[3], p3[3], p4[3];
        real8   v1[3], v2[3], v3[3], v4[3];
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
        for (i = 0; i < home->newNodeKeyPtr; i++) 
        {
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
    
                            p1[X] = x1;  p1[Y] = y1;  p1[Z] = z1;
                            p2[X] = x2;  p2[Y] = y2;  p2[Z] = z2;
                            p3[X] = x3;  p3[Y] = y3;  p3[Z] = z3;
                            p4[X] = x4;  p4[Y] = y4;  p4[Z] = z4;

                            v1[X] = vx1;  v1[Y] = vy1;  v1[Z] = vz1;
                            v2[X] = vx2;  v2[Y] = vy2;  v2[Z] = vz2;
                            v3[X] = vx3;  v3[Y] = vy3;  v3[Z] = vz3;
                            v4[X] = vx4;  v4[Y] = vy4;  v4[Z] = vz4;

/*
 *                          It is possible to have a zero-length segment
 *                          (created by a previous collision).  If we find
 *                          such a segment, do not try to use it in any
 *                          subsequent collisions.
 */
                            vec1[X] = x2 - x1;
                            vec1[Y] = y2 - y1;
                            vec1[Z] = z2 - z1;

                            if (DotProduct(vec1, vec1) < 1.0e-20) {
                                continue;
                            }

                            vec2[X] = x4 - x3;
                            vec2[Y] = y4 - y3;
                            vec2[Z] = z4 - z3;

                            if (DotProduct(vec2, vec2) < 1.0e-20) {
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
 *                          First check of segments already intersect.  If
 *                          not, find out if they will collide in the future.
 *
 *                          Note: If the separation between the segments is
 *                          less than 1b, treat them as if they are already
 *                          intersecting
 */
                            if (dist2 < 1.0) {
                                cPoint[X] = x1 + vec1[X] * L1;
                                cPoint[Y] = y1 + vec1[Y] * L1;
                                cPoint[Z] = z1 + vec1[Z] * L1;
                                collisionConditionIsMet = 1;
                            } else if (dist2 < mindist2) {
/*
 *                              FIX ME!  If segments are to far away but
 *                              moving fast, they can pass right through
 *                              each other with no collision
 *
 *                              Only do a rigorous treatment of points 
 *                              within the distance filter.  Find the 
 *                              collision point, collision time and the
 *                              points on the segments where the two segments
 *                              will be colliding.
 */
                                FindCollisionPointAndTime(home, p1, p2, p3, p4,
                                                      v1, v2, v3, v4, cPoint,
                                                      &cTime, &dist2, &L1, &L2);

                                collisionConditionIsMet = ((cTime > 0.0) &&
                                                           (cTime < 10.0));
                            } else {
                                collisionConditionIsMet = 0;
                            }

                            if (collisionConditionIsMet) {
/*
 *                              Segments are unconnected and colliding.
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
                                adjustCollisionPoint = 0;

                                vec1[X] = cPoint[X] - p1[X];
                                vec1[Y] = cPoint[Y] - p1[Y];
                                vec1[Z] = cPoint[Z] - p1[Z];
    
                                vec2[X] = cPoint[X] - p2[X];
                                vec2[Y] = cPoint[Y] - p2[Y];
                                vec2[Z] = cPoint[Z] - p2[Z];

                                close2node1 = (DotProduct(vec1,vec1)<mindist2);
                                close2node2 = (DotProduct(vec2,vec2)<mindist2);
    
                                if ((node2->myTag.domainID != thisDomain) &&
                                    close2node2) {
                                    continue;
                                }

                                if (close2node1) {
                                    mergenode1 = node1;
                                    splitSeg1 = 0;
                                    adjustCollisionPoint = 1;
                                } else if (close2node2) {
                                    mergenode1 = node2;
                                    splitSeg1 = 0;
                                    adjustCollisionPoint = 1;
                                } else {
                                    splitSeg1 = 1;
                                }

/*
 *                              If we need to add a new node to the first
 *                              segment, do it now.
 */
                                if (splitSeg1) {
                                     real8 pos0[3], pos1[3];

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
    
                                     pos0[X] = x1;   pos1[X] = x2;
                                     pos0[Y] = y1;   pos1[Y] = y2;
                                     pos0[Z] = z1;   pos1[Z] = z2;
    
                                     pnew[X] = newx;
                                     pnew[Y] = newy;
                                     pnew[Z] = newz;
    
                                     newNodeVel[X] = newvx;
                                     newNodeVel[Y] = newvy;
                                     newNodeVel[Z] = newvz;

                                     nodeVel[X] = vx1;
                                     nodeVel[Y] = vy1;
                                     nodeVel[Z] = vz1;

                                     burg1[X] = node1->burgX[arm12];
                                     burg1[Y] = node1->burgY[arm12];
                                     burg1[Z] = node1->burgZ[arm12];

                                     FindSubFSeg(home, pos0, pos1, burg1, oldfp0,
                                                 oldfp1, pnew, f0seg1, f1seg1,
                                                 f0seg2, f1seg2);
    
                                     oldTag1 = node1->myTag;
                                     oldTag2 = node2->myTag;

                                     FoldBox(param, &pnew[X], &pnew[Y], &pnew[Z]);

                                     splitStatus = SplitNode(home,
                                                             OPCLASS_COLLISION,
                                                             node1, pos0, pnew,
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
 *                              Identify the second node to be merged
 *
 *                              Note: The current domain owns node3 but may not
 *                              own node4.  If it does not own node4, we 
 *                              cannot allow the collision to use node4
 *                              even if the collision point is close to 
 *                              that node.
 */
                                vec1[X] = cPoint[X] - p3[X];
                                vec1[Y] = cPoint[Y] - p3[Y];
                                vec1[Z] = cPoint[Z] - p3[Z];
    
                                vec2[X] = cPoint[X] - p4[X];
                                vec2[Y] = cPoint[Y] - p4[Y];
                                vec2[Z] = cPoint[Z] - p4[Z];

                                close2node3 = (DotProduct(vec1,vec1)<mindist2);
                                close2node4 = (DotProduct(vec2,vec2)<mindist2);
    
                                if ((node4->myTag.domainID != thisDomain) &&
                                    close2node4) {
                                    continue;
                                }

                                if (close2node3) {
                                    mergenode2 = node3;
                                    splitSeg2 = 0;
                                    adjustCollisionPoint = 1;
                                } else if (close2node4) {
                                    mergenode2 = node4;
                                    splitSeg2 = 0;
                                    adjustCollisionPoint = 1;
                                } else {
                                    splitSeg2 = 1;
                                }

/*
 *                              If we need to add a new node to the second
 *                              segment, do it now.
 */
                                if (splitSeg2) {
                                     real8 pos0[3], pos1[3];

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
    
                                     pos0[X] = x3;   pos1[X] = x4;
                                     pos0[Y] = y3;   pos1[Y] = y4;
                                     pos0[Z] = z3;   pos1[Z] = z4;
    
                                     pnew[X] = newx;
                                     pnew[Y] = newy;
                                     pnew[Z] = newz;
    
                                     newNodeVel[X] = newvx;
                                     newNodeVel[Y] = newvy;
                                     newNodeVel[Z] = newvz;

                                     nodeVel[X] = vx3;
                                     nodeVel[Y] = vy3;
                                     nodeVel[Z] = vz3;

                                     FindSubFSeg(home, pos0, pos1, burg1, oldfp0,
                                                 oldfp1, pnew, f0seg1, f1seg1,
                                                 f0seg2, f1seg2);
    
                                     oldTag1 = node3->myTag;
                                     oldTag2 = node4->myTag;

                                     FoldBox(param, &pnew[X], &pnew[Y], &pnew[Z]);

                                     splitStatus = SplitNode(home,
                                                            OPCLASS_COLLISION,
                                                            node3, pos0,
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
 *                             If the initially selected collision point
 *                             was close enough to one of the segment
 *                             endpoints that we decided to use the
 *                             endpoint node for the merge, we need to
 *                             recalculate the final collision point.
 */
                               if (adjustCollisionPoint) {
                                   FindCollisionPoint(home, mergenode1,
                                                      mergenode2, &newPos[X],
                                                      &newPos[Y], &newPos[Z]);
                               } else {
                                   newPos[X] = cPoint[X];
                                   newPos[Y] = cPoint[Y];
                                   newPos[Z] = cPoint[Z];
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
#endif  /* ifdef _ARLFEM */

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
 *                  Find the minimum distance between the the node1/node4
 *                  segment and the point at node3 to determine if they
 *                  should be collided.
 */
                    GetMinDist(x1, y1, z1, vx1, vy1, vz1,
                               x4, y4, z4, vx4, vy4, vz4,
                               x3, y3, z3, vx3, vy3, vz3,
                               x3, y3, z3, vx3, vy3, vz3,
                               &dist2, &ddist2dt, &L1, &L2);

                    if (dist2 < mindist2) {

                        p1[X] = x1;  p1[Y] = y1;  p1[Z] = z1;
                        p3[X] = x3;  p3[Y] = y3;  p3[Z] = z3;
                        p4[X] = x4;  p4[Y] = y4;  p4[Z] = z4;

                        v1[X] = vx1;  v1[Y] = vy1;  v1[Z] = vz1;
                        v3[X] = vx3;  v3[Y] = vy3;  v3[Z] = vz3;
                        v4[X] = vx4;  v4[Y] = vy4;  v4[Z] = vz4;

                        FindCollisionPointAndTime(home, p1, p4, p3, p3,
                                              v1, v4, v3, v3, cPoint, &cTime,
                                              &dist2, &L1, &L2);

                        collisionConditionIsMet = ((cTime > 0.0) &&
                                                   (cTime < 10.0));
                    } else {
                        collisionConditionIsMet = 0;
                    }

                    if (collisionConditionIsMet) {

/*
 *                      A collision should occur, so calculate the ideal
 *                      collision point between the node1/node4 segment and
 *                      point node3, using the 'closest' point identified above
 */

/*
 *                      Node3 is used as one of the collision points, but see
 *                      if the collision point is close to one of the existing
 *                      nodes on the node1/node4 segment
 */
                        mergenode1 = node3;

                        vec1[X] = cPoint[X] - x1;
                        vec1[Y] = cPoint[Y] - y1;
                        vec1[Z] = cPoint[Z] - z1;

                        vec2[X] = cPoint[X] - x4;
                        vec2[Y] = cPoint[Y] - y4;
                        vec2[Z] = cPoint[Z] - z4;

                        close2node1 = (DotProduct(vec1,vec1)<mindist2);
                        close2node4 = (DotProduct(vec2,vec2)<mindist2);

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
                             adjustCollisionPoint = 1;
                        } else if (close2node4) {
                             mergenode2 = node4;
                             adjustCollisionPoint = 1;
                        } else { 
                             real8 pos0[3], pos1[3];

                             adjustCollisionPoint = 0;
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

                             pos0[X] = x1;   pos1[X] = x4;
                             pos0[Y] = y1;   pos1[Y] = y4;
                             pos0[Z] = z1;   pos1[Z] = z4;

                             FindSubFSeg(home, pos0, pos1, burg1, oldfp0,
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
                                                     node1, pos0, newPos,
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
 *                      If the initially selected collision point
 *                      was close enough to one of the segment
 *                      endpoints that we decided to use the
 *                      endpoint node for the merge, we need to
 *                      recalculate the final collision point.
 */
                        if (adjustCollisionPoint) {
                            FindCollisionPoint(home, mergenode1,
                                               mergenode2, &newPos[X],
                                               &newPos[Y], &newPos[Z]);
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

                            (void)EvaluateMobility(home, targetNode,
                                                   &mobArgs);

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
        int globalCollisionCnt=0;
        MPI_Reduce(&localCollisionCnt, &globalCollisionCnt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        int globalCollisionCnt = localCollisionCnt;
#endif
        if (home->myDomain == 0) {
            printf("Predictive Collision count = %d\n", globalCollisionCnt);
        }
#endif

        TimerStop(home, COLLISION_HANDLING);

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       PredictiveCollisionsLG
 *      
 *-------------------------------------------------------------------------*/
void PredictiveCollisionsLG(Home_t *home)
{
        int     i, j, k, q, arm12, arm21, arm34, arm43, arm14;
        int     thisDomain, splitStatus, mergeStatus, mergeStatus_tmp;
        int     armAB, armBA;
        int     globalOp = 1, localOp = 0, didCollision;
        int     close2node1, close2node2, close2node3, close2node4;
        int     splitSeg1, splitSeg2, splitSeg3;
        int     cell2Index, nbrCell2Index, nextIndex;
        int     cell2X, cell2Y, cell2Z, cx, cy, cz;
        int     collisionConditionIsMet, adjustCollisionPoint=0;
        real8   mindist2, dist2, ddist2dt, L1, L2, half;
        real8   cTime, cPoint[3];
        real8   x1, y1, z1, vx1, vy1, vz1;
        real8   x2, y2, z2, vx2, vy2, vz2;
        real8   x3, y3, z3, vx3, vy3, vz3;
        real8   x4, y4, z4, vx4, vy4, vz4;
        real8   newx, newy, newz;
        real8   newvx, newvy, newvz;
        real8   pnew[3], burg1[3];
        real8   oldfp0[3], oldfp1[3];
        real8   f0seg1[3], f1seg1[3], f0seg2[3], f1seg2[3];
        real8   nodeVel[3], newNodeVel[3];
        real8   newPos[3], newVel[3];
        real8   vec1[3], vec2[3];
        real8   p1[3], p2[3], p3[3], p4[3];
        real8   v1[3], v2[3], v3[3], v4[3];
        Tag_t   oldTag1, oldTag2;
        Node_t  *node1, *node2, *node3, *node4, *tmpNbr;
        Node_t  *mergenode1, *mergenode2, *targetNode;
        Node_t  *splitNode1, *splitNode2;
        Param_t *param;
        MobArgs_t mobArgs;
	int     globalOpRestore;
#ifdef _ARLFEM
        int     resetSurfaceProperties;
        real8   femSurfaceNorm[3];
#endif
	NodeState_t origNode1State, origNode2State;
	NodeState_t origNode3State, origNode4State;


	InitState(&origNode1State);
	InitState(&origNode2State);
	InitState(&origNode3State);
	InitState(&origNode4State);

        thisDomain = home->myDomain;
        param      = home->param;

        mindist2 = param->rann * param->rann;
        half     = 0.5;

        int localCollisionCnt  = 0;
	int globalCollisionCnt = 0;
	real8 pos0[3], pos2[3];
	Node_t *tmpNode;

#ifdef DEBUG_TOPOLOGY_DOMAIN
        dbgDom = DEBUG_TOPOLOGY_DOMAIN;
#else
        dbgDom = -1;
#endif

        TimerStart(home, COLLISION_HANDLING);

/*
 *      Start looping through native nodes looking for segments to collide...
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) 
        {
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
    
                            p1[X] = x1;  p1[Y] = y1;  p1[Z] = z1;
                            p2[X] = x2;  p2[Y] = y2;  p2[Z] = z2;
                            p3[X] = x3;  p3[Y] = y3;  p3[Z] = z3;
                            p4[X] = x4;  p4[Y] = y4;  p4[Z] = z4;

                            v1[X] = vx1;  v1[Y] = vy1;  v1[Z] = vz1;
                            v2[X] = vx2;  v2[Y] = vy2;  v2[Z] = vz2;
                            v3[X] = vx3;  v3[Y] = vy3;  v3[Z] = vz3;
                            v4[X] = vx4;  v4[Y] = vy4;  v4[Z] = vz4;

/*
 *                          It is possible to have a zero-length segment
 *                          (created by a previous collision).  If we find
 *                          such a segment, do not try to use it in any
 *                          subsequent collisions.
 */
                            vec1[X] = x2 - x1;
                            vec1[Y] = y2 - y1;
                            vec1[Z] = z2 - z1;

                            if (DotProduct(vec1, vec1) < 1.0e-20) {
                                continue;
                            }

                            vec2[X] = x4 - x3;
                            vec2[Y] = y4 - y3;
                            vec2[Z] = z4 - z3;

                            if (DotProduct(vec2, vec2) < 1.0e-20) {
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
 *                          First check of segments already intersect.  If
 *                          not, find out if they will collide in the future.
 *
 *                          Note: If the separation between the segments is
 *                          less than 1b, treat them as if they are already
 *                          intersecting
 */
                            if (dist2 < 1.0) {
                                cPoint[X] = x1 + vec1[X] * L1;
                                cPoint[Y] = y1 + vec1[Y] * L1;
                                cPoint[Z] = z1 + vec1[Z] * L1;
                                collisionConditionIsMet = 1;
                            } else if (dist2 < mindist2) {
/*
 *                              FIX ME!  If segments are to far away but
 *                              moving fast, they can pass right through
 *                              each other with no collision
 *
 *                              Only do a rigorous treatment of points 
 *                              within the distance filter.  Find the 
 *                              collision point, collision time and the
 *                              points on the segments where the two segments
 *                              will be colliding.
 */
                                FindCollisionPointAndTime(home, p1, p2, p3, p4,
                                                      v1, v2, v3, v4, cPoint,
                                                      &cTime, &dist2, &L1, &L2);

                                collisionConditionIsMet = ((cTime > 0.0) &&
                                                           (cTime < 10.0));
                            } else {
                                collisionConditionIsMet = 0;
                            }

                            if (collisionConditionIsMet) {
/*
 *                              Segments are unconnected and colliding.
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
                                adjustCollisionPoint = 0;

                                vec1[X] = cPoint[X] - p1[X];
                                vec1[Y] = cPoint[Y] - p1[Y];
                                vec1[Z] = cPoint[Z] - p1[Z];
    
                                vec2[X] = cPoint[X] - p2[X];
                                vec2[Y] = cPoint[Y] - p2[Y];
                                vec2[Z] = cPoint[Z] - p2[Z];

                                close2node1 = (DotProduct(vec1,vec1)<mindist2);
                                close2node2 = (DotProduct(vec2,vec2)<mindist2);
    
                                if ((node2->myTag.domainID != thisDomain) &&
                                    close2node2) {
                                    continue;
                                }

                                if (close2node1) {
                                    mergenode1 = node1;
                                    splitSeg1 = 0;
                                    adjustCollisionPoint = 1;
                                } else if (close2node2) {
                                    mergenode1 = node2;
                                    splitSeg1 = 0;
                                    adjustCollisionPoint = 1;
                                } else {
                                    splitSeg1 = 1;
                                }

				if( (close2node1) && (node1->myTag.domainID != thisDomain) ) continue;
				if( (close2node2) && (node2->myTag.domainID != thisDomain) ) continue;

/*
 *                              Identify the second node to be merged
 *
 *                              Note: The current domain owns node3 but may not
 *                              own node4.  If it does not own node4, we 
 *                              cannot allow the collision to use node4
 *                              even if the collision point is close to 
 *                              that node.
 */
                                vec1[X] = cPoint[X] - p3[X];
                                vec1[Y] = cPoint[Y] - p3[Y];
                                vec1[Z] = cPoint[Z] - p3[Z];
    
                                vec2[X] = cPoint[X] - p4[X];
                                vec2[Y] = cPoint[Y] - p4[Y];
                                vec2[Z] = cPoint[Z] - p4[Z];

                                close2node3 = (DotProduct(vec1,vec1)<mindist2);
                                close2node4 = (DotProduct(vec2,vec2)<mindist2);
    
                                if ((node4->myTag.domainID != thisDomain) &&
                                    close2node4) {
                                    continue;
                                }

                                if (close2node3) {
                                    mergenode2 = node3;
                                    splitSeg2 = 0;
                                    adjustCollisionPoint = 1;
                                } else if (close2node4) {
                                    mergenode2 = node4;
                                    splitSeg2 = 0;
                                    adjustCollisionPoint = 1;
                                } else {
                                    splitSeg2 = 1;
                                }

				if( (close2node3) && (node3->myTag.domainID != thisDomain) ) continue;
				if( (close2node4) && (node4->myTag.domainID != thisDomain) ) continue;

				
/*
 *                              If we need to add a new node to the first
 *                              segment, do it now.
 */
                                if (splitSeg1) {
				  PreserveState(home, node1, &origNode1State);
				  PreserveState(home, node2, &origNode2State);

				  splitStatus = SplitSegment(home, node1, node2, 
							     arm12, L1, 
							     &mergenode1, localOp);
/*
 *                                If we were unable to split the node
 *                                go back to looking for more collision
 *                                candidates.
 */				  
				  if (splitStatus == SPLIT_FAILED) {
				    continue;
				  }
				}

/*
 *                              If we need to add a new node to the second
 *                              segment, do it now.
 */
                                if (splitSeg2) {

				  PreserveState(home, node3, &origNode3State);
				  PreserveState(home, node4, &origNode4State);

				  splitStatus = SplitSegment(home, node3, node4,
							     arm34, L2, 
							     &mergenode2, localOp);
/*
 *                                If we were unable to split the node
 *                                go back to looking for more collision
 *                                candidates.
 */				  
				  if (splitStatus == SPLIT_FAILED) {
				    continue;
				  }
				}
			        
    
/*
 *                             If the initially selected collision point
 *                             was close enough to one of the segment
 *                             endpoints that we decided to use the
 *                             endpoint node for the merge, we need to
 *                             recalculate the final collision point.
 */
				if (adjustCollisionPoint) {
				  FindCollisionPoint(home, mergenode1,
						     mergenode2, &newPos[X],
						     &newPos[Y], &newPos[Z]);
				} else {
				  newPos[X] = cPoint[X];
				  newPos[Y] = cPoint[Y];
				  newPos[Z] = cPoint[Z];
				}

				FoldBox(param, &newPos[X],&newPos[Y],&newPos[Z]);

				mergeStatus = 0;
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


				pos0[0] = node1->x; pos0[1] = node1->y; pos0[2] = node1->z;
				pos2[0] = node3->x; pos2[1] = node3->y; pos2[2] = node3->z;
				
				globalOpRestore = 0;
				mergeStatus = 0;

				if(MergeNodePreCondition(home, OPCLASS_COLLISION, mergenode1, mergenode2) != -1){
				  globalOpRestore = 1;
				  
				  // switch localOp to globalOp
				  if(splitSeg1) {
				    MergeNode(home, OPCLASS_COLLISION, mergenode1, node1, pos0, &tmpNode, &mergeStatus, localOp);
				    if(mergeStatus == 0)
				      Fatal("[%s:%d] Fail to merge two nodes (%d,%d) and (%d,%d) ",
					    __func__,__LINE__, node1->myTag.domainID, node1->myTag.index,
					    mergenode1->myTag.domainID, mergenode1->myTag.index);
				    if(node1->myTag != tmpNode->myTag)
				      Fatal("[%s:%d] Node (%d,%d) was removed! ",
					    __func__,__LINE__, node1->myTag.domainID, node1->myTag.index);

				    
				    RestoreState(home, node1, &origNode1State);
				    RestoreState(home, node2, &origNode2State);

				    splitStatus = SplitSegment(home, node1, node2,
							       arm12, L1, 
							       &mergenode1, globalOp);
				    
				    if (splitStatus == SPLIT_FAILED)
				      Fatal("[%s:%d] Fail to split node (%d,%d) ",
					    __func__,__LINE__, node1->myTag.domainID, node1->myTag.index);
				    didCollision = 1;				    
				  }
				

				  if(splitSeg2) {
				    
				    MergeNode(home, OPCLASS_COLLISION, mergenode2, node3, pos2, &tmpNode, &mergeStatus, localOp);
				    if(mergeStatus == 0)
				      Fatal("[%s:%d] Fail to merge two nodes (%d,%d) and (%d,%d) ",
					    __func__,__LINE__, node2->myTag.domainID, node2->myTag.index,
					    mergenode2->myTag.domainID, mergenode2->myTag.index);
				    
				    if(node3->myTag != tmpNode->myTag)
				      Fatal("[%s:%d] Node (%d,%d) was removed! ",
					    __func__,__LINE__, node1->myTag.domainID, node1->myTag.index);
				    
				    RestoreState(home, node3, &origNode3State);
				    RestoreState(home, node4, &origNode4State);
				    
				    splitStatus = SplitSegment(home, node3, node4, 
							       arm34, L2, 
							       &mergenode2, globalOp);
				    
				    if (splitStatus == SPLIT_FAILED)
				      Fatal("[%s:%d] Fail to split node (%d,%d)  ",
					    __func__,__LINE__, node3->myTag.domainID, node3->myTag.index);
				    
				  }
				  
				  real8 mergenode1Pos[3] = {mergenode1->x, mergenode1->y, mergenode1->z};
				  real8 mergenode2Pos[3] = {mergenode2->x, mergenode2->y, mergenode2->z};
				  MergeNode(home, OPCLASS_COLLISION,
					    mergenode1, mergenode2, newPos,
					    &targetNode,
					    &mergeStatus, globalOp);

				  if ((mergeStatus & MERGE_SUCCESS) == 0) {
				    continue;
				  }

				  if (targetNode != (Node_t *)NULL) {
				    targetNode->flags |= NO_COLLISIONS;
				  }
				} else {
				  mergeStatus = 0;
				}

/*
 *                              If the merge did not succeed, go back and
 *                              continue looking for collision candidates.
 */
				if ((mergeStatus & MERGE_SUCCESS) == 0) {

				  if(globalOpRestore){
				    if(splitSeg1) {
				      MergeNode(home, OPCLASS_COLLISION, mergenode1, node1, pos0, &tmpNode, &mergeStatus, globalOp);
				      if(mergeStatus == 0)
					Fatal("[%s:%d] Fail to merge two nodes (%d,%d) and (%d,%d) ",
					      __func__,__LINE__, node1->myTag.domainID, node1->myTag.index,
					      mergenode1->myTag.domainID, mergenode1->myTag.index);
				      if(node1->myTag != tmpNode->myTag)
					Fatal("[%s:%d] Node (%d,%d) was removed! ",
					      __func__,__LINE__, node1->myTag.domainID, node1->myTag.index);
				      
				    }
				    if(splitSeg2) {
				      MergeNode(home, OPCLASS_COLLISION, mergenode2, node3, pos2, &tmpNode, &mergeStatus, globalOp);
				      if(mergeStatus == 0)
					Fatal("[%s:%d] Fail to merge two nodes (%d,%d) and (%d,%d) ",
					      __func__,__LINE__, node2->myTag.domainID, node2->myTag.index,
					      mergenode2->myTag.domainID, mergenode2->myTag.index);
				      if(node3->myTag != tmpNode->myTag)
					Fatal("[%s:%d] Node (%d,%d) was removed! ",
					      __func__,__LINE__, node3->myTag.domainID, node3->myTag.index);
				      
				    }
				  } else {
				    if(splitSeg1){
				      MergeNode(home, OPCLASS_COLLISION, mergenode1, node1, pos0, &tmpNode, &mergeStatus, localOp);
				      if(mergeStatus == 0)
					Fatal("[%s:%d] Fail to merge two nodes (%d,%d) and (%d,%d) ",
					      __func__,__LINE__, node1->myTag.domainID, node1->myTag.index,
					      mergenode1->myTag.domainID, mergenode1->myTag.index);
				      if(node1->myTag != tmpNode->myTag){
					Fatal("[%s:%d] Node (%d,%d) was removed! ",
					      __func__,__LINE__, node1->myTag.domainID, node1->myTag.index);
				      }
				    }
				    if(splitSeg2) {
				      MergeNode(home, OPCLASS_COLLISION, mergenode2, node3, pos2, &tmpNode, &mergeStatus, localOp);
				      if(mergeStatus == 0)
					Fatal("[%s:%d] Fail to merge two nodes (%d,%d) and (%d,%d) ",
					      __func__,__LINE__, node2->myTag.domainID, node2->myTag.index,
					      mergenode2->myTag.domainID, mergenode2->myTag.index);
				    }
				  }
				  if (splitSeg1) {				  
				    RestoreState(home, node1, &origNode1State);
				    RestoreState(home, node2, &origNode2State);
				  }
				  if (splitSeg2) {
				    RestoreState(home, node3, &origNode3State);
				    RestoreState(home, node4, &origNode4State);
				  }
				  continue;
				}

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
 *                  Find the minimum distance between the the node1/node4
 *                  segment and the point at node3 to determine if they
 *                  should be collided.
 */
                    GetMinDist(x1, y1, z1, vx1, vy1, vz1,
                               x4, y4, z4, vx4, vy4, vz4,
                               x3, y3, z3, vx3, vy3, vz3,
                               x3, y3, z3, vx3, vy3, vz3,
                               &dist2, &ddist2dt, &L1, &L2);

                    if (dist2 < mindist2) {

                        p1[X] = x1;  p1[Y] = y1;  p1[Z] = z1;
                        p3[X] = x3;  p3[Y] = y3;  p3[Z] = z3;
                        p4[X] = x4;  p4[Y] = y4;  p4[Z] = z4;

                        v1[X] = vx1;  v1[Y] = vy1;  v1[Z] = vz1;
                        v3[X] = vx3;  v3[Y] = vy3;  v3[Z] = vz3;
                        v4[X] = vx4;  v4[Y] = vy4;  v4[Z] = vz4;

                        FindCollisionPointAndTime(home, p1, p4, p3, p3,
                                              v1, v4, v3, v3, cPoint, &cTime,
                                              &dist2, &L1, &L2);

                        collisionConditionIsMet = ((cTime > 0.0) &&
                                                   (cTime < 10.0));
                    } else {
                        collisionConditionIsMet = 0;
                    }

                    if (collisionConditionIsMet) {

		        DiscardState(&origNode1State);
		        DiscardState(&origNode2State);
		        DiscardState(&origNode3State);
		        DiscardState(&origNode4State);

/*
 *                      A collision should occur, so calculate the ideal
 *                      collision point between the node1/node4 segment and
 *                      point node3, using the 'closest' point identified above
 */

/*
 *                      Node3 is used as one of the collision points, but see
 *                      if the collision point is close to one of the existing
 *                      nodes on the node1/node4 segment
 */
                        mergenode1 = node3;

                        vec1[X] = cPoint[X] - x1;
                        vec1[Y] = cPoint[Y] - y1;
                        vec1[Z] = cPoint[Z] - z1;

                        vec2[X] = cPoint[X] - x4;
                        vec2[Y] = cPoint[Y] - y4;
                        vec2[Z] = cPoint[Z] - z4;

                        close2node1 = (DotProduct(vec1,vec1)<mindist2);
                        close2node4 = (DotProduct(vec2,vec2)<mindist2);

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

			splitSeg3 = 0;
                        if (close2node1) {
                             mergenode2 = node1;
                             adjustCollisionPoint = 1;
                        } else if (close2node4) {
                             mergenode2 = node4;
                             adjustCollisionPoint = 1;
                        } else {
			  splitSeg3 = 1;
			}

/*
 *                      If the collision point is close to one of the
 *                      endpoints of the node1/node4 segment, but the
 *                      node in question is not owned by the current
 *                      domain, skip the collision this time around.
 */
			if((close2node1)||(close2node4)){
			  if( ( (mergenode2->myTag==node1->myTag)  &&
			        (node1->myTag.domainID!=thisDomain)) ||
			      ( (mergenode2->myTag==node4->myTag)  &&
				(node4->myTag.domainID!=thisDomain)) ){
                            continue;
			  }
			}

			if(splitSeg3) {

			  PreserveState(home, node1, &origNode1State);
			  PreserveState(home, node3, &origNode3State);
			  PreserveState(home, node4, &origNode4State);

			  splitStatus = SplitHingeSegment(home, node1, node4,
							  arm14, L1,
							  &mergenode2,
							  localOp);
			  if(splitStatus == false)
			    continue;

			  newPos[0] = mergenode2->x;
			  newPos[1] = mergenode2->y;
			  newPos[2] = mergenode2->z;

			}

/*
 *                      If the initially selected collision point
 *                      was close enough to one of the segment
 *                      endpoints that we decided to use the
 *                      endpoint node for the merge, we need to
 *                      recalculate the final collision point.
 */
                        if (adjustCollisionPoint) {
                            FindCollisionPoint(home, mergenode1,
                                               mergenode2, &newPos[X],
                                               &newPos[Y], &newPos[Z]);
                        }

                        FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);

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

			globalOpRestore = 0;
			
			real8 pos0[3];
			pos0[0] = node1->x;  pos0[1] = node1->y;  pos0[2] = node1->z;
			mergeStatus = 0;

			if(MergeNodePreCondition(home, OPCLASS_COLLISION, mergenode1, mergenode2) != -1 ){

			  int   nd1dom, nd1id, nd4dom, nd4id;
			  nd1dom = node1->myTag.domainID;
			  nd1id = node1->myTag.index;
			  nd4dom = node4->myTag.domainID;
			  nd4id = node4->myTag.index;			  
			  globalOpRestore = 1;

			  if(splitSeg3){
			    MergeNode(home, OPCLASS_COLLISION,  mergenode2, node1, 
				      pos0, &tmpNode, &mergeStatus, localOp);
			    
			    if(mergeStatus == 0)
			      Fatal("[%s:%d] Fail to merge two nodes (%d,%d) and (%d,%d) ",
				    __func__,__LINE__, node1->myTag.domainID, node1->myTag.index,
				    mergenode2->myTag.domainID, mergenode2->myTag.index);
			    if(node1->myTag != tmpNode->myTag)
			      Fatal("[%s:%d] Node (%d,%d) was removed! ",
				    __func__,__LINE__, node1->myTag.domainID, node1->myTag.index);

			    RestoreState(home, node1, &origNode1State);
			    RestoreState(home, node3, &origNode3State);			    
			    RestoreState(home, node4, &origNode4State);

			    splitStatus =  SplitHingeSegment(home, node1, node4,
							     arm14, L1, 
							     &mergenode2,
							     globalOp);
			  
			    if(splitStatus == SPLIT_FAILED)
			      Fatal("[%s:%d] Fail to split node (%d,%d)  ",
				    __func__,__LINE__, node1->myTag.domainID, node1->myTag.index);
			    
			    newPos[0] = mergenode2->x;
			    newPos[1] = mergenode2->y;
			    newPos[2] = mergenode2->z;
			    
/*
 *                          If the initially selected collision point
 *                          was close enough to one of the segment
 *                          endpoints that we decided to use the
 *                          endpoint node for the merge, we need to
 *                          recalculate the final collision point.
 */
			    if (adjustCollisionPoint) {
			      FindCollisionPoint(home, mergenode1,
						 mergenode2, &newPos[X],
						 &newPos[Y], &newPos[Z]);
			    }

			    FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);

			  }

/*
 *                      The merge is going to happen, so there's a few
 *                      more nodes we'll need to mark for force/velocity
 *                      re-evaluation.
 */
			  MarkNodeForceObsolete(home, node1);
			  MarkNodeForceObsolete(home, node3);
			  MarkNodeForceObsolete(home, node4);
			  
			  MergeNode(home, OPCLASS_COLLISION,
				    mergenode1, mergenode2, newPos,
				    &targetNode,
				    &mergeStatus, globalOp);

			  if ( (mergeStatus & MERGE_SUCCESS) && (targetNode != (Node_t *)NULL) ){
			    targetNode->flags |= NO_COLLISIONS;

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



                            (void)EvaluateMobility(home, targetNode,
                                                   &mobArgs);

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
			} else {
			  mergeStatus = 0;
			}

			if(mergeStatus == 0){
			  if(globalOpRestore){
			    if(splitSeg3){
			      MergeNode(home, OPCLASS_COLLISION, mergenode2, node1, pos0, &tmpNode, &mergeStatus, globalOp);
			      if(mergeStatus == 0)
				Fatal("[%s:%d] Fail to merge two nodes (%d,%d) and (%d,%d) ",
				      __func__,__LINE__, node1->myTag.domainID, node1->myTag.index,
				      mergenode2->myTag.domainID, mergenode2->myTag.index);
			      if(node1->myTag != tmpNode->myTag)
				Fatal("[%s:%d] Node (%d,%d) was removed! ",
				      __func__,__LINE__, node1->myTag.domainID, node1->myTag.index);
			      
			    }
			  } else {

			    if(splitSeg3){ 
			      MergeNode(home, OPCLASS_COLLISION, mergenode2, node1, pos0, &tmpNode, &mergeStatus, localOp);
			      if(mergeStatus == 0)
				Fatal("[%s:%d] Fail to merge two nodes (%d,%d) and (%d,%d) ",
				      __func__,__LINE__, node1->myTag.domainID, node1->myTag.index,
				      mergenode2->myTag.domainID, mergenode2->myTag.index);
			      if(node1->myTag != tmpNode->myTag)
				Fatal("[%s:%d] Node (%d,%d) was removed! ",
				      __func__,__LINE__, node1->myTag.domainID, node1->myTag.index);			      
			    }

			  }

			  if(splitSeg3){

			    RestoreState(home, node1, &origNode1State);
			    RestoreState(home, node3, &origNode3State);
			    RestoreState(home, node4, &origNode4State);
			    
			  }
			}
			localCollisionCnt++;


                    }  /* conditions for zip were met */
                }  /* for (k = 0...) */
            }  /* for (j = 0...) */
        }  /* for (i = 0...) */

/*
 *      Now we have to loop for two node collisions.
 *      The first target node is the node with 3 neighbors
 *      The second merged node is its one of neighbor
 *      The distance between these two nodes should be small
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
            if (node1->flags & NO_COLLISIONS) continue;

            for (j = 0; j < node1->numNbrs; j++) {

                if (node1->myTag.domainID != node1->nbrTag[j].domainID)
                    continue;
		
		x1 = node1->x;  y1 = node1->y;  z1 = node1->z;
		
		node2 = GetNodeFromTag(home, node1->nbrTag[j]);
		if (node2->flags & NO_COLLISIONS) continue;
	    
		x2 = node2->x;  y2 = node2->y;  z2 = node2->z;
		
		PBCPOSITION(param, x1, y1, z1, &x2, &y2, &z2);
		
                real8 vec1x = x2 - x1;
		real8 vec1y = y2 - y1;
		real8 vec1z = z2 - z1;
		ZImage(param, &vec1x, &vec1y, &vec1z);
		
                real8 r1 = sqrt(vec1x*vec1x + vec1y*vec1y + vec1z*vec1z);

		if(r1 > param->rann)
		  continue;

		real8 v1x = node1->vX;
		real8 v1y = node1->vY;
		real8 v1z = node1->vZ;
		real8 v2x = node2->vX;
		real8 v2y = node2->vY;
		real8 v2z = node2->vZ;

		
		newPos[0] = x2;
		newPos[1] = y2;
		newPos[2] = z2;

		if (param->enforceGlidePlanes){
		  printf("[WARNING] Two nodes (%d,%d) and (%d,%d) are collide, but it is not gurantee its merged node satisfies the glide plane constraints\n",
			 node1->myTag.domainID, node1->myTag.index, node2->myTag.domainID, node2->myTag.index);
		}
		
		MergeNode(home, OPCLASS_COLLISION, node1,
			  node2, newPos, &targetNode,
			  &mergeStatus, globalOp);

		if (targetNode != (Node_t *)NULL) {
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

		  (void)EvaluateMobility(home, targetNode,
					 &mobArgs);

		  targetNode->flags |= NODE_RESET_FORCES;
		  targetNode->flags |= NO_COLLISIONS;
		}
	    }
	}

#ifdef DEBUG_LOG_COLLISIONS
#ifdef PARALLEL
        int globalCollisionCnt=0;
        MPI_Reduce(&localCollisionCnt, &globalCollisionCnt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        int globalCollisionCnt = localCollisionCnt;
#endif
        if (home->myDomain == 0) {
            printf("Predictive Collision count = %d\n", globalCollisionCnt);
        }
#endif

        TimerStop(home, COLLISION_HANDLING);

        return;
}
